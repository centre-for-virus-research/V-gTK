import argparse
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_MIN_YEARLY_SEQUENCES = 25


def _load_yearly_cohort(conn):
    return pd.read_sql_query(
        """
        SELECT DISTINCT
            primary_accession,
            CAST(collection_year AS INTEGER) AS collection_year
        FROM meta_data
        WHERE TRIM(COALESCE(primary_accession, '')) != ''
          AND COALESCE(TRIM(exclusion_status), '') != '1'
          AND TRIM(COALESCE(collection_year, '')) GLOB '[0-9][0-9][0-9][0-9]'
        """,
        conn,
    )


def _load_relevant_mutation_counts(conn):
    return pd.read_sql_query(
        """
        SELECT
            primary_accession,
            total_relevant_mutation_count
        FROM sequence_relevant_mutation_summary
        """,
        conn,
    )


def _load_completed_signature_counts(conn):
    return pd.read_sql_query(
        """
        SELECT
            primary_accession,
            COUNT(DISTINCT signature_id) AS completed_signature_count
        FROM completed_signatures_only
        GROUP BY primary_accession
        """,
        conn,
    )


def _load_signature_drug_hits(conn):
    return pd.read_sql_query(
        """
        SELECT DISTINCT
            completed_signatures_only.primary_accession,
            TRIM(mutation_catalog.drug) AS drug
        FROM completed_signatures_only
        INNER JOIN mutation_catalog
            ON mutation_catalog.signature_id = completed_signatures_only.signature_id
        WHERE TRIM(COALESCE(mutation_catalog.drug, '')) != ''
        """,
        conn,
    )


def _smooth_weighted_proportions(years, proportions, sample_sizes, span_fraction=0.22):
    years = np.asarray(years, dtype=float)
    proportions = np.asarray(proportions, dtype=float)
    sample_sizes = np.asarray(sample_sizes, dtype=float)

    if len(years) == 0:
        return np.array([])
    if len(years) == 1:
        return proportions.copy()

    year_range = max(years.max() - years.min(), 1.0)
    bandwidth = max(year_range * span_fraction, 2.0)
    smoothed = np.empty_like(proportions, dtype=float)

    for index, target_year in enumerate(years):
        distances = np.abs(years - target_year) / bandwidth
        kernel_weights = np.where(distances < 1, (1 - distances**3) ** 3, 0.0)
        combined_weights = kernel_weights * sample_sizes
        if not np.any(combined_weights > 0):
            smoothed[index] = proportions[index]
            continue
        smoothed[index] = np.average(proportions, weights=combined_weights)

    return np.clip(smoothed, 0.0, 1.0)


def _wilson_95_ci(successes, totals):
    successes = pd.Series(successes, dtype=float)
    totals = pd.Series(totals, dtype=float).clip(lower=1)
    z = 1.96
    proportions = successes / totals
    denominator = 1 + (z**2 / totals)
    centre = (proportions + (z**2 / (2 * totals))) / denominator
    margin = (
        (z / denominator)
        * np.sqrt((proportions * (1 - proportions) / totals) + (z**2 / (4 * totals**2)))
    )
    lower = (centre - margin).clip(lower=0.0)
    upper = (centre + margin).clip(upper=1.0)
    return lower, upper


def build_yearly_summaries(db_path, min_yearly_sequences=DEFAULT_MIN_YEARLY_SEQUENCES):
    conn = sqlite3.connect(str(db_path))
    try:
        cohort = _load_yearly_cohort(conn)
        mutation_counts = _load_relevant_mutation_counts(conn)
        signature_counts = _load_completed_signature_counts(conn)
        signature_drug_hits = _load_signature_drug_hits(conn)
    finally:
        conn.close()

    if cohort.empty:
        raise ValueError("No non-excluded sequences with a valid four-digit collection_year were found in meta_data")

    cohort = cohort.copy()
    cohort["collection_year"] = cohort["collection_year"].astype(int)

    mutation_counts = mutation_counts.rename(columns={"total_relevant_mutation_count": "relevant_mutation_count"})
    mutation_counts["relevant_mutation_count"] = pd.to_numeric(
        mutation_counts["relevant_mutation_count"], errors="coerce"
    ).fillna(0)

    signature_counts["completed_signature_count"] = pd.to_numeric(
        signature_counts["completed_signature_count"], errors="coerce"
    ).fillna(0)

    mutation_cohort = cohort.merge(mutation_counts, on="primary_accession", how="left")
    mutation_cohort["relevant_mutation_count"] = mutation_cohort["relevant_mutation_count"].fillna(0)

    signature_cohort = cohort.merge(signature_counts, on="primary_accession", how="left")
    signature_cohort["completed_signature_count"] = signature_cohort["completed_signature_count"].fillna(0)

    yearly_denominator = (
        cohort.groupby("collection_year", as_index=False)
        .agg(total_sequences=("primary_accession", "nunique"))
        .sort_values("collection_year")
    )
    eligible_years = yearly_denominator.loc[
        yearly_denominator["total_sequences"] >= int(min_yearly_sequences),
        "collection_year",
    ]

    if eligible_years.empty:
        raise ValueError(
            f"No collection_year values met the minimum yearly sample size threshold of {int(min_yearly_sequences)}"
        )

    mutation_cohort = mutation_cohort[mutation_cohort["collection_year"].isin(eligible_years)]
    signature_cohort = signature_cohort[signature_cohort["collection_year"].isin(eligible_years)]
    cohort = cohort[cohort["collection_year"].isin(eligible_years)]

    yearly_mutation_summary = (
        mutation_cohort.groupby("collection_year", as_index=False)
        .agg(
            sequence_count=("primary_accession", "nunique"),
            mean_relevant_mutation_count=("relevant_mutation_count", "mean"),
        )
        .sort_values("collection_year")
    )

    yearly_signature_summary = (
        signature_cohort.groupby("collection_year", as_index=False)
        .agg(
            sequence_count=("primary_accession", "nunique"),
            mean_completed_signature_count=("completed_signature_count", "mean"),
        )
        .sort_values("collection_year")
    )
    yearly_denominator = yearly_denominator[yearly_denominator["collection_year"].isin(eligible_years)]

    if signature_drug_hits.empty:
        yearly_drug_proportion = pd.DataFrame(
            columns=[
                "collection_year",
                "drug",
                "sequences_with_signature",
                "total_sequences",
                "proportion",
                "smoothed_proportion",
            ]
        )
    else:
        signature_drug_hits = cohort.merge(signature_drug_hits, on="primary_accession", how="inner")
        yearly_drug_counts = (
            signature_drug_hits.groupby(["collection_year", "drug"], as_index=False)
            .agg(sequences_with_signature=("primary_accession", "nunique"))
            .sort_values(["drug", "collection_year"])
        )
        drug_year_grid = (
            yearly_denominator.assign(key=1)
            .merge(
                pd.DataFrame({"drug": sorted(yearly_drug_counts["drug"].unique()), "key": 1}),
                on="key",
                how="inner",
            )
            .drop(columns="key")
        )
        yearly_drug_proportion = drug_year_grid.merge(
            yearly_drug_counts,
            on=["collection_year", "drug"],
            how="left",
        )
        yearly_drug_proportion["sequences_with_signature"] = yearly_drug_proportion[
            "sequences_with_signature"
        ].fillna(0)
        yearly_drug_proportion["proportion"] = (
            yearly_drug_proportion["sequences_with_signature"]
            / yearly_drug_proportion["total_sequences"]
        )
        ci_lower, ci_upper = _wilson_95_ci(
            yearly_drug_proportion["sequences_with_signature"],
            yearly_drug_proportion["total_sequences"],
        )
        yearly_drug_proportion["ci_lower"] = ci_lower
        yearly_drug_proportion["ci_upper"] = ci_upper
        yearly_drug_proportion["smoothed_proportion"] = (
            yearly_drug_proportion.sort_values(["drug", "collection_year"])
            .groupby("drug", group_keys=False)
            .apply(
                lambda group: pd.Series(
                    _smooth_weighted_proportions(
                        group["collection_year"].to_numpy(),
                        group["proportion"].to_numpy(),
                        group["total_sequences"].to_numpy(),
                    ),
                    index=group.index,
                )
            )
        )

    yearly_sample_size = yearly_denominator.rename(columns={"total_sequences": "sequence_count"})

    return {
        "yearly_mutation_summary": yearly_mutation_summary,
        "yearly_signature_summary": yearly_signature_summary,
        "yearly_drug_proportion": yearly_drug_proportion,
        "yearly_sample_size": yearly_sample_size,
    }


def _plot_mean_series(ax, df, value_column, title, ylabel, color):
    ax.plot(df["collection_year"], df[value_column], color=color, linewidth=2.4, marker="o", markersize=4.2)
    ax.set_title(title, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.tick_params(axis="both", labelsize=15)


def _apply_year_axis(ax, yearly_sample_size):
    if yearly_sample_size.empty:
        return
    min_year = int(yearly_sample_size["collection_year"].min())
    max_year = int(yearly_sample_size["collection_year"].max())
    ax.set_xlim(min_year - 0.5, max_year + 0.5)


def _plot_drug_proportions(ax, df, yearly_sample_size, show_ci_cloud=False):
    if df.empty:
        ax.text(0.5, 0.5, "No completed signature drug mappings available", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    ax2 = ax.twinx()
    ax2.bar(
        yearly_sample_size["collection_year"],
        yearly_sample_size["sequence_count"],
        width=0.9,
        color="#d9d9d9",
        alpha=0.35,
        zorder=0,
    )
    ax2.set_ylabel("Yearly sample size", fontsize=17, color="#5a5a5a", labelpad=18)
    ax2.tick_params(axis="y", labelsize=15, colors="#5a5a5a", pad=3)
    ax2.spines["top"].set_visible(False)

    cmap = plt.get_cmap("tab20")
    for index, drug in enumerate(sorted(df["drug"].unique())):
        drug_df = df[df["drug"] == drug].sort_values("collection_year")
        color = cmap(index % 20)
        marker_sizes = 18 + 40 * np.sqrt(drug_df["total_sequences"] / drug_df["total_sequences"].max())
        if show_ci_cloud and {"ci_lower", "ci_upper"}.issubset(drug_df.columns):
            ax.fill_between(
                drug_df["collection_year"],
                drug_df["ci_lower"],
                drug_df["ci_upper"],
                color=color,
                alpha=0.12,
                linewidth=0,
                zorder=1,
            )
        ax.scatter(
            drug_df["collection_year"],
            drug_df["proportion"],
            color=color,
            s=marker_sizes,
            alpha=0.28,
            linewidths=0,
            zorder=3,
        )
        ax.plot(
            drug_df["collection_year"],
            drug_df["smoothed_proportion"],
            label=drug,
            color=color,
            linewidth=2.7,
            zorder=4,
        )
    ax.set_title("Completed-signature frequency by antiviral", fontsize=24)
    ax.set_ylabel("Proportion of sequences", fontsize=18)
    ax.set_xlabel("Collection year", fontsize=18)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.tick_params(axis="both", labelsize=16)
    _apply_year_axis(ax, yearly_sample_size)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.18, 1.02),
        frameon=False,
        title="Antiviral",
        fontsize=16,
        title_fontsize=18,
        ncol=1,
    )


def generate_plots(db_path, output_dir, min_yearly_sequences=DEFAULT_MIN_YEARLY_SEQUENCES):
    db_path = Path(db_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries = build_yearly_summaries(db_path, min_yearly_sequences=min_yearly_sequences)
    mutation_summary = summaries["yearly_mutation_summary"]
    signature_summary = summaries["yearly_signature_summary"]
    drug_proportion = summaries["yearly_drug_proportion"]
    yearly_sample_size = summaries["yearly_sample_size"]

    figure, axes = plt.subplots(3, 1, figsize=(16, 20), constrained_layout=False)

    _plot_mean_series(
        axes[0],
        mutation_summary,
        "mean_relevant_mutation_count",
        "Average relevant mutation count per year",
        "Mean relevant mutation count",
        "#1f4e79",
    )
    _apply_year_axis(axes[0], yearly_sample_size)

    _plot_mean_series(
        axes[1],
        signature_summary,
        "mean_completed_signature_count",
        "Average completed signature count per year",
        "Mean completed signature count",
        "#8c2d04",
    )
    axes[1].set_xlabel("Collection year", fontsize=16)
    _apply_year_axis(axes[1], yearly_sample_size)

    _plot_drug_proportions(axes[2], drug_proportion, yearly_sample_size)
    figure.subplots_adjust(right=0.76, hspace=0.34)

    combined_plot_path = output_dir / "hcv_mutation_trends.png"
    figure.savefig(combined_plot_path, dpi=260, bbox_inches="tight")
    plt.close(figure)

    mutation_plot_path = output_dir / "mean_relevant_mutation_count_by_year.png"
    signature_plot_path = output_dir / "mean_completed_signature_count_by_year.png"
    drug_plot_path = output_dir / "antiviral_signature_proportion_by_year.png"
    drug_plot_with_ci_path = output_dir / "antiviral_signature_proportion_by_year_with_ci.png"

    for data_frame, value_column, title, ylabel, color, plot_path in [
        (
            mutation_summary,
            "mean_relevant_mutation_count",
            "Average relevant mutation count per year",
            "Mean relevant mutation count",
            "#1f4e79",
            mutation_plot_path,
        ),
        (
            signature_summary,
            "mean_completed_signature_count",
            "Average completed signature count per year",
            "Mean completed signature count",
            "#8c2d04",
            signature_plot_path,
        ),
    ]:
        figure, ax = plt.subplots(figsize=(11, 5.5))
        _plot_mean_series(ax, data_frame, value_column, title, ylabel, color)
        ax.set_xlabel("Collection year", fontsize=16)
        _apply_year_axis(ax, yearly_sample_size)
        figure.savefig(plot_path, dpi=240, bbox_inches="tight")
        plt.close(figure)

    figure, ax = plt.subplots(figsize=(18, 11))
    _plot_drug_proportions(ax, drug_proportion, yearly_sample_size)
    figure.subplots_adjust(right=0.69)
    figure.savefig(drug_plot_path, dpi=260, bbox_inches="tight")
    plt.close(figure)

    figure, ax = plt.subplots(figsize=(18, 11))
    _plot_drug_proportions(ax, drug_proportion, yearly_sample_size, show_ci_cloud=True)
    figure.subplots_adjust(right=0.69)
    figure.savefig(drug_plot_with_ci_path, dpi=260, bbox_inches="tight")
    plt.close(figure)

    mutation_summary.to_csv(output_dir / "yearly_relevant_mutation_summary.tsv", sep="\t", index=False)
    signature_summary.to_csv(output_dir / "yearly_completed_signature_summary.tsv", sep="\t", index=False)
    drug_proportion.to_csv(output_dir / "yearly_antiviral_signature_proportion.tsv", sep="\t", index=False)
    yearly_sample_size.to_csv(output_dir / "yearly_sample_size.tsv", sep="\t", index=False)

    summaries["plot_paths"] = {
        "combined": combined_plot_path,
        "mutations": mutation_plot_path,
        "completed_signatures": signature_plot_path,
        "drug_proportions": drug_plot_path,
        "drug_proportions_with_ci": drug_plot_with_ci_path,
    }
    return summaries


def main():
    parser = argparse.ArgumentParser(
        description="Plot HCV mutation and completed-signature trends by collection year from a SQLite DB"
    )
    parser.add_argument("--db", required=True, help="Path to the SQLite DB")
    parser.add_argument("--output_dir", required=True, help="Directory to write plot files into")
    parser.add_argument(
        "--min_yearly_sequences",
        type=int,
        default=DEFAULT_MIN_YEARLY_SEQUENCES,
        help="Only plot years with at least this many non-excluded sequences",
    )
    args = parser.parse_args()

    summaries = generate_plots(args.db, args.output_dir, min_yearly_sequences=args.min_yearly_sequences)
    print(f"Saved plots to {Path(args.output_dir).resolve()}")
    print(f"Years plotted: {int(summaries['yearly_mutation_summary']['collection_year'].min())} - {int(summaries['yearly_mutation_summary']['collection_year'].max())}")
    print(
        "Sample size range: "
        f"{int(summaries['yearly_sample_size']['sequence_count'].min())} - "
        f"{int(summaries['yearly_sample_size']['sequence_count'].max())} sequences per year"
    )


if __name__ == "__main__":
    # Tiny test DB example:
    # python scripts/PlotHcvMutationTrends.py --db /tmp/hcv_mutation_plot_test.db --output_dir /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/plots_mutations --min_yearly_sequences 2
    # Full DB example:
    # python scripts/PlotHcvMutationTrends.py --db /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/HCV_full.db --output_dir /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/plots_mutations --min_yearly_sequences 25
    main()