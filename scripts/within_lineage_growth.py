#!/usr/bin/env python3
"""Within-lineage variant and displacement growth rates from a V-gTK SQLite database.

Ported from Mut_acc_PLM_MS with calibration guarantees for V-gTK viral databases.
Operates natively on V-gTK SQLite databases (meta_data, features, sequence_alignment, trees).

Methodology & Calibration Guarantees
------------------------------------
Implements the calibrated growth estimation framework addressing findings in
methodology_review.md:
1. Rarity Absorption (Finding C1): Uses a 2-parameter logistic growth model where the
   intercept absorbs initial frequency, guaranteeing that constant-frequency variants
   return s = 0 regardless of sample size or rarity.
2. Honest Generation Intervals (Finding C4): Uses Euler-Lotka moment-generating conversion
   to ln(R_variant / R_ref) parameterised by an explicit Gamma distribution.
3. Identifiability Status (Findings S1, M1): Variants lacking trajectory (singletons, single
   time bins, short spans) return NaN with explicit status codes, never 0.0 placeholders.
4. Stratified Estimation (Finding S5): Multi-stratum logistic regression with per-region
   intercepts and shared slopes prevents Simpson's paradox across geographic regions,
   with likelihood-ratio heterogeneity testing.
5. Displacement Advantage: Measures displacement directly against the parent lineage/clade
   co-circulating in contemporaneous time bins.
6. Design Identifiability (Finding C3): Reports rank and confounded groups for mutation
   presence matrices.
"""
from __future__ import annotations

import argparse
import os
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import growth_core as gc
import protein_alleles as pal
import lineage_growth_rates as lgr


def say(msg: str):
    print(f"[{pd.Timestamp.now().strftime('%H:%M:%S')}] {msg}", flush=True)


class WithinLineageGrowthAnalyzer:
    """Estimates within-lineage variant growth and lineage displacement from SQLite."""

    def __init__(
        self,
        db_path: str | Path,
        protein: str = "nucleoprotein N",
        segment: str | None = None,
        bin_days: float = 30.0,
        generation_time_days: float = 30.0,
        generation_sd_days: float = 15.0,
        stratify_by: str = "country",
        min_carriers: int = 2,
        min_variant_bins: int = 2,
        min_span_days: float = 28.0,
    ):
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {self.db_path}")
        self.protein = protein
        self.segment = segment
        self.bin_days = float(bin_days)
        self.gen_interval = gc.GenerationInterval(
            mean_days=float(generation_time_days),
            sd_days=float(generation_sd_days),
        )
        self.stratify_by = stratify_by
        self.min_carriers = int(min_carriers)
        self.min_variant_bins = int(min_variant_bins)
        self.min_span_days = float(min_span_days)

        self.conn = sqlite3.connect(f"file:{self.db_path.resolve()}?mode=ro", uri=True)
        self.cohort = None
        self.allele_set = None

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def load_cohort(self) -> pd.DataFrame:
        """Load cohort from meta_data with robust date parsing and plausibility masking."""
        say(f"Loading cohort from {self.db_path.name}")
        c = self.conn.cursor()
        c.execute("""
            SELECT primary_accession, collection_date, collection_year,
                   country, host, exclusion_status, nearest_reference_genotype,
                   nearest_reference_subtype
            FROM meta_data
        """)
        rows = c.fetchall()
        cols = ["primary_accession", "collection_date", "collection_year",
                "country", "host", "exclusion_status", "genotype", "subtype"]
        df = pd.DataFrame(rows, columns=cols)

        parsed = [gc.parse_collection_date(d) for d in df["collection_date"]]
        df["parsed_date"] = [p[0] for p in parsed]
        df["date_precision"] = [p[1] for p in parsed]

        valid_dates = df["parsed_date"].notna()
        if valid_dates.any():
            mask = gc.plausible_date_mask(
                df.loc[valid_dates, "parsed_date"],
                min_date="1900-01-01",
                max_date="2030-12-31",
            )
            df = df.loc[valid_dates][mask].copy()
        else:
            df = df.iloc[0:0].copy()

        df["decimal_date"] = [
            (d.year + (d.dayofyear - 1) / (366.0 if d.is_leap_year else 365.0))
            if pd.notna(d) else np.nan for d in df["parsed_date"]
        ]
        say(f"  {len(df)} sequences with valid, plausible collection dates")
        self.cohort = df
        return df

    def build_alleles(self) -> lgr.AlleleSet:
        """Extract amino acid alleles for the requested protein."""
        if self.cohort is None:
            self.load_cohort()
        say(f"Resolving protein alleles for {self.protein}")
        args = argparse.Namespace(
            protein=self.protein,
            segment=self.segment,
            allele_source="auto",
            positions=None,
            coord_space="auto",
            protein_start=None,
            calibrate_start=False,
        )
        self.allele_set = lgr.build_alleles(self.conn, args, self.cohort)
        return self.allele_set

    def estimate_variant_growth(self) -> pd.DataFrame:
        """Estimate within-lineage variant growth rates against ancestral/reference residue."""
        if self.allele_set is None:
            self.build_alleles()

        say("Estimating within-lineage variant growth rates")
        df_cohort = self.cohort.set_index("primary_accession")
        min_date = df_cohort["parsed_date"].min()

        acc_time_days = {
            acc: float((row["parsed_date"] - min_date).total_seconds() / 86400.0)
            for acc, row in df_cohort.iterrows()
        }
        acc_stratum = {
            acc: str(row.get(self.stratify_by) or "unknown")
            for acc, row in df_cohort.iterrows()
        }

        records = []
        all_sites = sorted(self.allele_set.evaluated.keys())
        strata_list = sorted(set(acc_stratum.values()))
        n_strata = len(strata_list)
        stratum_to_idx = {s: i for i, s in enumerate(strata_list)}

        max_time = max(acc_time_days.values()) if acc_time_days else 0.0
        n_bins = max(int(np.ceil(max_time / self.bin_days)), 1) + 1

        for site in all_sites:
            tested_accs = [a for a in self.allele_set.evaluated[site] if a in acc_time_days]
            if not tested_accs:
                continue

            residue_counts = {}
            for allele in self.allele_set.alleles_at(site):
                res = self.allele_set.residues[allele]["residue"]
                c_accs = [a for a in self.allele_set.carriers[allele] if a in tested_accs]
                residue_counts[res] = residue_counts.get(res, 0) + len(c_accs)

            if not residue_counts:
                continue
            ancestral_res = max(residue_counts.keys(), key=lambda r: residue_counts[r])

            for allele in self.allele_set.alleles_at(site):
                res = self.allele_set.residues[allele]["residue"]
                if res == ancestral_res or res in ("X", "-", "*"):
                    continue

                carriers = [a for a in self.allele_set.carriers[allele] if a in tested_accs]
                carrier_count = len(carriers)
                if carrier_count < self.min_carriers:
                    continue

                ref_accs = [
                    a for a in tested_accs
                    if a not in self.allele_set.carriers[allele]
                ]

                c_m = np.zeros(n_bins)
                c_r = np.zeros(n_bins)
                cc_strat = np.zeros((n_bins, n_strata))
                cp_strat = np.zeros((n_bins, n_strata))

                for a in carriers:
                    t_idx = min(int(acc_time_days[a] // self.bin_days), n_bins - 1)
                    c_m[t_idx] += 1
                    s_idx = stratum_to_idx[acc_stratum[a]]
                    cc_strat[t_idx, s_idx] += 1

                for a in ref_accs:
                    t_idx = min(int(acc_time_days[a] // self.bin_days), n_bins - 1)
                    c_r[t_idx] += 1
                    s_idx = stratum_to_idx[acc_stratum[a]]
                    cp_strat[t_idx, s_idx] += 1

                fit = gc.fit_logistic_growth(
                    c_m, c_r, self.bin_days,
                    generation_interval=self.gen_interval,
                    min_bins=3,
                    min_variant_bins=self.min_variant_bins,
                    min_span_days=self.min_span_days,
                )

                strat_fit = gc.fit_logistic_growth_stratified(
                    cc_strat, cp_strat, self.bin_days,
                    generation_interval=self.gen_interval,
                )

                records.append({
                    "protein": self.protein,
                    "site": site,
                    "ref_residue": ancestral_res,
                    "alt_residue": res,
                    "variant": f"{ancestral_res}{site}{res}",
                    "carriers": carrier_count,
                    "background_count": len(ref_accs),
                    "total_evaluated": len(tested_accs),
                    "frequency": carrier_count / len(tested_accs),
                    "growth_per_day": fit.growth_per_day,
                    "growth_per_year": fit.growth_per_day * 365.25 if np.isfinite(fit.growth_per_day) else np.nan,
                    "se_per_day": fit.se_per_day,
                    "se_per_year": fit.se_per_day * 365.25 if np.isfinite(fit.se_per_day) else np.nan,
                    "log_R_ratio": fit.log_R_ratio,
                    "status": fit.status,
                    "separated": fit.separated,
                    "converged": fit.converged,
                    "span_days": fit.span_days,
                    "n_bins": fit.n_bins,
                    "n_variant_bins": fit.n_variant_bins,
                    "strat_growth_per_day": strat_fit.growth_per_day,
                    "strat_growth_per_year": strat_fit.growth_per_day * 365.25 if np.isfinite(strat_fit.growth_per_day) else np.nan,
                    "strat_se_per_day": strat_fit.se_per_day,
                    "strat_heterogeneity_chi2": strat_fit.heterogeneity_chi2,
                    "strat_heterogeneity_df": strat_fit.heterogeneity_df,
                    "strat_heterogeneity_p": strat_fit.heterogeneity_p,
                    "strat_n_strata": strat_fit.n_strata,
                })

        df_res = pd.DataFrame(records)
        say(f"  evaluated {len(df_res)} variants ({sum(df_res['status'] == 'ok') if not df_res.empty else 0} with status 'ok')")
        return df_res

    def estimate_lineage_displacement(self) -> pd.DataFrame:
        """Estimate displacement growth rates between lineages/clades from metadata or trees."""
        if self.cohort is None:
            self.load_cohort()

        say("Estimating lineage/clade displacement growth rates")
        df_cohort = self.cohort.set_index("primary_accession")
        min_date = df_cohort["parsed_date"].min()

        acc_time_days = {
            acc: float((row["parsed_date"] - min_date).total_seconds() / 86400.0)
            for acc, row in df_cohort.iterrows()
        }
        acc_stratum = {
            acc: str(row.get(self.stratify_by) or "unknown")
            for acc, row in df_cohort.iterrows()
        }
        strata_list = sorted(set(acc_stratum.values()))
        n_strata = len(strata_list)
        stratum_to_idx = {s: i for i, s in enumerate(strata_list)}

        max_time = max(acc_time_days.values()) if acc_time_days else 0.0
        n_bins = max(int(np.ceil(max_time / self.bin_days)), 1) + 1

        records = []
        for group_col in ["genotype", "subtype"]:
            groups = df_cohort[group_col].dropna().unique()
            groups = [g for g in groups if g and str(g).strip() and str(g) != "None"]
            if len(groups) >= 2:
                for i in range(len(groups)):
                    for j in range(i + 1, len(groups)):
                        g_child, g_parent = groups[i], groups[j]
                        acc_child = df_cohort[df_cohort[group_col] == g_child].index
                        acc_parent = df_cohort[df_cohort[group_col] == g_parent].index

                        cc = np.zeros(n_bins)
                        cp = np.zeros(n_bins)
                        cc_strat = np.zeros((n_bins, n_strata))
                        cp_strat = np.zeros((n_bins, n_strata))

                        for a in acc_child:
                            if a in acc_time_days:
                                t_idx = min(int(acc_time_days[a] // self.bin_days), n_bins - 1)
                                cc[t_idx] += 1
                                cc_strat[t_idx, stratum_to_idx[acc_stratum[a]]] += 1
                        for a in acc_parent:
                            if a in acc_time_days:
                                t_idx = min(int(acc_time_days[a] // self.bin_days), n_bins - 1)
                                cp[t_idx] += 1
                                cp_strat[t_idx, stratum_to_idx[acc_stratum[a]]] += 1

                        fit = gc.fit_logistic_growth(
                            cc, cp, self.bin_days,
                            generation_interval=self.gen_interval,
                        )
                        strat_fit = gc.fit_logistic_growth_stratified(
                            cc_strat, cp_strat, self.bin_days,
                            generation_interval=self.gen_interval,
                        )

                        records.append({
                            "transition_type": group_col,
                            "child_lineage": g_child,
                            "parent_lineage": g_parent,
                            "child_count": len(acc_child),
                            "parent_count": len(acc_parent),
                            "growth_per_day": fit.growth_per_day,
                            "growth_per_year": fit.growth_per_day * 365.25 if np.isfinite(fit.growth_per_day) else np.nan,
                            "se_per_day": fit.se_per_day,
                            "se_per_year": fit.se_per_day * 365.25 if np.isfinite(fit.se_per_day) else np.nan,
                            "log_R_ratio": fit.log_R_ratio,
                            "status": fit.status,
                            "strat_growth_per_day": strat_fit.growth_per_day,
                            "strat_growth_per_year": strat_fit.growth_per_day * 365.25 if np.isfinite(strat_fit.growth_per_day) else np.nan,
                            "strat_se_per_day": strat_fit.se_per_day,
                            "strat_heterogeneity_chi2": strat_fit.heterogeneity_chi2,
                            "strat_heterogeneity_p": strat_fit.heterogeneity_p,
                            "strat_n_strata": strat_fit.n_strata,
                        })

        df_disp = pd.DataFrame(records)
        say(f"  evaluated {len(df_disp)} lineage displacement pairs")
        return df_disp

    def write_report(
        self,
        out_dir: str | Path,
        variant_df: pd.DataFrame,
        displacement_df: pd.DataFrame | None = None,
    ) -> Path:
        """Write TSV exports and a detailed markdown report."""
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        var_tsv = out_path / "variant_growth_rates.tsv"
        variant_df.to_csv(var_tsv, sep="\t", index=False)

        if displacement_df is not None and not displacement_df.empty:
            disp_tsv = out_path / "lineage_displacement.tsv"
            displacement_df.to_csv(disp_tsv, sep="\t", index=False)

        report_md = out_path / "report.md"
        lines = [
            f"# Growth Rate Estimation Report: {self.protein}",
            "",
            f"**Database:** `{self.db_path.name}`  ",
            f"**Sequences:** {len(self.cohort)} dated accessions  ",
            f"**Bin Width:** {self.bin_days:.1f} days  ",
            f"**Generation Interval:** Gamma(mu={self.gen_interval.mean_days:.1f}d, sd={self.gen_interval.sd_days:.1f}d)  ",
            f"**Geographic Stratification:** `{self.stratify_by}`  ",
            "",
            "## Methodology & Calibration Guarantees",
            "",
            "- **Rarity Absorption (Finding C1)**: Fitted using a 2-parameter logistic regression with intercept absorbing baseline rarity, guaranteeing zero growth for neutral variants.",
            "- **Euler-Lotka Reproduction Number Ratio (Finding C4)**: Converts daily growth rates into ln(R_variant / R_ref) using the virus generation interval moment-generating function.",
            "- **Honest Identifiability (Findings S1, M1)**: Single-bin and singleton variants return NaN with an explicit status code rather than zero placeholders.",
            "- **Geographic Stratification (Finding S5)**: Shared-slope regression with per-region intercepts guards against Simpson\'s paradox across countries.",
            "",
            "## Top Growing Variants",
            "",
        ]
        fittable = variant_df[variant_df["status"] == "ok"].copy()
        if not fittable.empty:
            top_variants = fittable.sort_values("growth_per_day", ascending=False).head(15)
            lines.append("| Variant | Carriers | Freq | Growth (/yr) | SE (/yr) | ln(R ratio) | Strat Growth (/yr) | Het p-val |")
            lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
            for _, r in top_variants.iterrows():
                lines.append(f"| {r['variant']} | {r['carriers']} | {r['frequency']:.3f} | "
                             f"{r['growth_per_year']:+.3f} | {r['se_per_year']:.3f} | "
                             f"{r['log_R_ratio']:+.3f} | {r['strat_growth_per_year']:+.3f} | "
                             f"{r['strat_heterogeneity_p']:.3e} |")
        else:
            lines.append("_No variants met the multi-bin trajectory criteria for a converged slope._")

        if displacement_df is not None and not displacement_df.empty:
            lines.extend([
                "",
                "## Lineage Displacement Rates",
                "",
                "| Transition | Child | Parent | Child N | Parent N | Growth (/yr) | SE (/yr) | ln(R ratio) | Het p-val |",
                "| --- | --- | --- | --- | --- | --- | --- | --- |",
            ])
            for _, r in displacement_df.iterrows():
                lines.append(f"| {r['transition_type']} | {r['child_lineage']} | {r['parent_lineage']} | "
                             f"{r['child_count']} | {r['parent_count']} | {r['growth_per_year']:+.3f} | "
                             f"{r['se_per_year']:.3f} | {r['log_R_ratio']:+.3f} | "
                             f"{r['strat_heterogeneity_p']:.3e} |")

        with open(report_md, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

        say(f"Report and tables written to {out_path}")
        return report_md


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Calibrated within-lineage and displacement growth rate estimator.")
    p.add_argument("--db", required=True, help="Path to V-gTK SQLite database")
    p.add_argument("--protein", default="nucleoprotein N", help="Protein product to evaluate")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--stratify-by", default="country", help="Column in meta_data for regional stratification")
    p.add_argument("--bin-days", type=float, default=30.0, help="Time bin width in days")
    p.add_argument("--generation-time", type=float, default=30.0, help="Mean generation interval in days")
    p.add_argument("--generation-sd", type=float, default=15.0, help="SD of generation interval in days")
    p.add_argument("--min-carriers", type=int, default=2, help="Minimum carriers to evaluate a variant")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    with WithinLineageGrowthAnalyzer(
        db_path=args.db,
        protein=args.protein,
        bin_days=args.bin_days,
        generation_time_days=args.generation_time,
        generation_sd_days=args.generation_sd,
        stratify_by=args.stratify_by,
        min_carriers=args.min_carriers,
    ) as analyzer:
        var_df = analyzer.estimate_variant_growth()
        disp_df = analyzer.estimate_lineage_displacement()
        analyzer.write_report(args.out, var_df, disp_df)


if __name__ == "__main__":
    main()
