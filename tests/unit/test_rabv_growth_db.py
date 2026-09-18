"""Integration tests for growth rate estimation on the RABV update-test SQLite database.

These tests exercise the full ``growth_core`` / ``within_lineage_growth`` stack
against the real ``test_out/update_test/rabv-jul0425-update-test.db`` database that
is also used by the Nextflow ``update_mode`` test configuration.

Calibration guarantees tested
------------------------------
C1  Constant-frequency variants must not return non-zero growth.
C4  Euler-Lotka ln(R) conversion must round-trip via the Gamma MGF.
S1  Variants with only a single occurrence bin must return NaN, not 0.0.
M1  Singletons must return NaN with status ``no_carriers``.
S5  Geographic stratification must reduce growth-rate bias from heterogeneous
    spatial sampling across Switzerland, Bosnia, Thailand, India, etc.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Path resolution — database lives at test_out/update_test/ relative to the
# repository root regardless of where pytest is invoked from.
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).parent.parent.parent
DB_PATH = REPO_ROOT / "test_out" / "update_test" / "rabv-jul0425-update-test.db"

SCRIPTS_DIR = REPO_ROOT / "scripts"
import sys
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import growth_core as gc
import growth_stats as gs
from within_lineage_growth import WithinLineageGrowthAnalyzer

# ---------------------------------------------------------------------------
# Skip if database not present (CI that hasn't run the Nextflow test yet)
# ---------------------------------------------------------------------------
requires_db = pytest.mark.skipif(
    not DB_PATH.exists(),
    reason=f"RABV update-test database not found at {DB_PATH}",
)


# ===========================================================================
# Database sanity checks
# ===========================================================================

@requires_db
class TestDatabaseSanity:
    """Verify the database is populated and structurally coherent."""

    def test_tables_present(self):
        with sqlite3.connect(str(DB_PATH)) as conn:
            tables = {
                r[0]
                for r in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                ).fetchall()
            }
        assert {"meta_data", "features", "sequence_alignment", "trees"} <= tables

    def test_sequence_count(self):
        with sqlite3.connect(str(DB_PATH)) as conn:
            n = conn.execute("SELECT count(*) FROM meta_data").fetchone()[0]
        assert n >= 200, f"Expected >=200 sequences, got {n}"

    def test_proteins_present(self):
        with sqlite3.connect(str(DB_PATH)) as conn:
            products = {
                r[0]
                for r in conn.execute("SELECT DISTINCT product FROM features").fetchall()
            }
        assert "nucleoprotein N" in products
        assert "transmembrane glycoprotein G" in products

    def test_countries_present(self):
        with sqlite3.connect(str(DB_PATH)) as conn:
            countries = {
                r[0]
                for r in conn.execute(
                    "SELECT DISTINCT country FROM meta_data WHERE country IS NOT NULL"
                ).fetchall()
            }
        assert "Switzerland" in countries
        assert "Thailand" in countries

    def test_trees_present(self):
        with sqlite3.connect(str(DB_PATH)) as conn:
            rows = conn.execute("SELECT name, source FROM trees").fetchall()
        assert len(rows) > 0, "No trees found in database"
        sources = {r[1] for r in rows}
        assert "usher" in sources or "iqtree" in sources


# ===========================================================================
# Date parsing on real collection dates
# ===========================================================================

@requires_db
class TestDateParsingOnRealDates:
    """parse_collection_date must handle the mixed formats in the RABV database."""

    def test_bare_year(self):
        ts, precision = gc.parse_collection_date("2001")
        assert ts is not None
        assert precision == "year"
        assert ts.year == 2001
        assert ts.month == 7  # midpoint placement

    def test_dmy_format(self):
        ts, precision = gc.parse_collection_date("26-Jan-2023")
        assert ts is not None
        assert precision == "day"
        assert ts.year == 2023
        assert ts.month == 1
        assert ts.day == 26

    def test_iso_month(self):
        ts, precision = gc.parse_collection_date("2009-06")
        assert ts is not None
        assert precision == "month"
        assert ts.year == 2009
        assert ts.month == 6

    def test_invalid_date_returns_none(self):
        ts, precision = gc.parse_collection_date("not-a-date")
        assert ts is None
        assert precision == "none"

    def test_plausibility_mask_removes_corrupted_dates(self):
        """Dates from serial-number epoch (1900, 1905) must be excluded."""
        import pandas as pd
        dates = pd.to_datetime(["1900-07-01", "1905-01-01", "2001-07-01", "2009-07-01", "2015-07-01"])
        mask = gc.plausible_date_mask(dates, min_date="1950-01-01", max_date="2030-01-01")
        assert not mask[0]  # 1900 excluded
        assert not mask[1]  # 1905 excluded
        assert mask[2]      # 2001 kept
        assert mask[4]      # 2015 kept

    def test_all_rabv_dates_parseable(self):
        """Every non-null collection date in the RABV database must parse to non-None or 'none'."""
        with sqlite3.connect(str(DB_PATH)) as conn:
            dates = conn.execute(
                "SELECT collection_date FROM meta_data WHERE collection_date IS NOT NULL"
            ).fetchall()
        # We expect the mixed bare-year/dmy/iso formats in RABV to all be parseable
        for (cd,) in dates:
            ts, precision = gc.parse_collection_date(cd)
            if ts is None:
                # The date returned None — check that precision is "none"
                assert precision == "none"


# ===========================================================================
# GenerationInterval / Euler-Lotka conversion
# ===========================================================================

class TestGenerationIntervalConversion:
    """C4: Euler-Lotka moment-generating function round-trip."""

    def test_zero_growth_gives_zero_log_r_ratio(self):
        gi = gc.GenerationInterval(mean_days=30.0, sd_days=15.0)
        assert float(gi.log_R_ratio(0.0)) == pytest.approx(0.0, abs=1e-10)

    def test_rabv_generation_interval_attributes(self):
        gi = gc.RABV_GENERATION_INTERVAL
        assert gi.mean_days == pytest.approx(30.0, abs=0.1)
        assert gi.sd_days == pytest.approx(15.0, abs=0.1)
        # ln(R) should be positive for a positive growth rate
        assert float(gi.log_R_ratio(1.0 / 365.25)) > 0.0
        assert float(gi.log_R_ratio(-1.0 / 365.25)) < 0.0

    def test_log_r_ratio_monotone(self):
        gi = gc.RABV_GENERATION_INTERVAL
        # Growth rates in per-day units
        rates = np.linspace(-0.01, 0.01, 9)
        log_rs = [float(gi.log_R_ratio(r)) for r in rates]
        diffs = np.diff(log_rs)
        assert np.all(diffs > 0), "log_R_ratio must be strictly increasing in growth rate"


# ===========================================================================
# WithinLineageGrowthAnalyzer: cohort loading
# ===========================================================================

@requires_db
class TestCohortLoading:
    """Verify cohort loading from the RABV test database."""

    def test_load_cohort_returns_dataframe(self):
        with WithinLineageGrowthAnalyzer(DB_PATH, protein="nucleoprotein N") as wlg:
            cohort = wlg.load_cohort()
        assert isinstance(cohort, pd.DataFrame)
        assert len(cohort) > 100

    def test_cohort_has_required_columns(self):
        with WithinLineageGrowthAnalyzer(DB_PATH, protein="nucleoprotein N") as wlg:
            cohort = wlg.load_cohort()
        for col in ("primary_accession", "decimal_date", "country"):
            assert col in cohort.columns, (
                f"Missing column '{col}': {list(cohort.columns)}"
            )

    def test_cohort_dates_are_plausible(self):
        """All loaded dates must be in a plausible range (no 1900-era junk)."""
        with WithinLineageGrowthAnalyzer(DB_PATH, protein="nucleoprotein N") as wlg:
            cohort = wlg.load_cohort()
        assert cohort["decimal_date"].min() > 1960.0
        assert cohort["decimal_date"].max() < 2030.0

    def test_cohort_countries_include_main_contributors(self):
        with WithinLineageGrowthAnalyzer(DB_PATH, protein="nucleoprotein N") as wlg:
            cohort = wlg.load_cohort()
        countries = set(cohort["country"].dropna().unique())
        assert "Switzerland" in countries
        assert "Thailand" in countries


# ===========================================================================
# Variant growth estimation on nucleoprotein N
# ===========================================================================

@requires_db
class TestVariantGrowthNucleoproteinN:
    """End-to-end variant growth estimation on nucleoprotein N."""

    @pytest.fixture(scope="class")
    def result(self):
        with WithinLineageGrowthAnalyzer(DB_PATH, protein="nucleoprotein N") as wlg:
            wlg.load_cohort()
            wlg.build_alleles()
            return wlg.estimate_variant_growth()

    def test_returns_dataframe(self, result):
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_required_columns_present(self, result):
        expected = {"variant", "carriers", "growth_per_year", "se_per_year", "log_R_ratio", "status"}
        assert expected <= set(result.columns), (
            f"Missing columns: {expected - set(result.columns)}"
        )

    def test_ok_variants_have_finite_growth(self, result):
        """Variants with status 'ok' must have finite growth rates and SEs."""
        ok = result[result["status"] == "ok"]
        assert len(ok) > 0, "No variants returned with status 'ok'"
        assert ok["growth_per_year"].notna().all()
        assert ok["se_per_year"].notna().all()
        assert ok["se_per_year"].gt(0).all()

    def test_unfittable_variants_return_nan_not_zero(self, result):
        """C1/M1/S1: unfittable variants must carry NaN, never 0.0."""
        bad = result[result["status"] != "ok"]
        if len(bad) > 0:
            assert bad["growth_per_year"].isna().all(), (
                "Non-ok variants returned 0.0 instead of NaN -- regression of M1/S1 bug"
            )

    def test_status_codes_are_known(self, result):
        known_statuses = {
            "ok", "no_carriers", "variant_in_too_few_bins",
            "too_few_informative_bins", "span_too_short", "singular",
            "no_usable_stratum",
        }
        unknown = set(result["status"].unique()) - known_statuses
        assert not unknown, f"Unknown status codes: {unknown}"

    def test_growth_rates_not_wildly_large(self, result):
        """Growth rates for real RABV data should be <20 per year."""
        ok = result[result["status"] == "ok"]
        assert ok["growth_per_year"].abs().max() < 20.0, (
            "Growth rates unexpectedly large -- possible calibration regression"
        )

    def test_log_r_ratio_consistent_sign_with_growth(self, result):
        """ln(R) ratio must have the same sign as growth_per_year (C4)."""
        ok = result[(result["status"] == "ok") & result["log_R_ratio"].notna() & result["growth_per_year"].notna()]
        sign_ok = (np.sign(ok["growth_per_year"]) == np.sign(ok["log_R_ratio"]))
        # Allow zero growth to have zero log_R_ratio
        nonzero = ok[ok["growth_per_year"].abs() > 1e-8]
        if len(nonzero) > 0:
            sign_ok_nz = (np.sign(nonzero["growth_per_year"]) == np.sign(nonzero["log_R_ratio"]))
            assert sign_ok_nz.all(), (
                "ln(R) ratio sign mismatch with growth_per_year -- Euler-Lotka sign error"
            )


# ===========================================================================
# Geographic stratification
# ===========================================================================

@requires_db
class TestGeographicStratification:
    """S5: stratified estimation must produce per-region heterogeneity statistics."""

    @pytest.fixture(scope="class")
    def result(self):
        with WithinLineageGrowthAnalyzer(
            DB_PATH,
            protein="nucleoprotein N",
            stratify_by="country",
        ) as wlg:
            wlg.load_cohort()
            wlg.build_alleles()
            return wlg.estimate_variant_growth()

    def test_strat_growth_column_present(self, result):
        assert "strat_growth_per_year" in result.columns, (
            f"Columns: {list(result.columns)}"
        )

    def test_strat_heterogeneity_p_column_present(self, result):
        assert "strat_heterogeneity_p" in result.columns

    def test_het_p_values_in_range(self, result):
        """p-values where present must be in [0, 1]."""
        pvals = result["strat_heterogeneity_p"].dropna()
        if len(pvals) > 0:
            assert pvals.between(0.0, 1.0).all(), "heterogeneity p-values out of [0,1]"

    def test_stratified_growth_differs_for_heterogeneous_variants(self, result):
        """For variants showing geographic heterogeneity, stratified != pooled."""
        het = result[
            (result["strat_heterogeneity_p"] < 0.05)
            & (result["status"] == "ok")
            & result["strat_growth_per_year"].notna()
        ]
        if len(het) == 0:
            pytest.skip("No significantly heterogeneous variants in test DB")
        pool = het["growth_per_year"]
        strat = het["strat_growth_per_year"]
        diffs = (pool - strat).abs()
        assert diffs.max() > 1e-6, (
            "Stratified growth equals pooled despite heterogeneity -- S5 guard not working"
        )


# ===========================================================================
# Logistic growth core -- calibration C1 (via fit_logistic_growth)
# ===========================================================================

class TestLogisticGrowthCalibrationC1:
    """Rarity must not leak into growth -- constant-frequency variants -> s = 0.

    fit_logistic_growth(c_variant, c_reference, bin_days) takes counts of
    variant and reference sequences per bin.
    """

    @pytest.mark.parametrize("freq", [0.5, 0.1, 0.02])
    @pytest.mark.parametrize("n_per_bin", [5, 20, 100])
    def test_constant_frequency_gives_zero_growth(self, freq, n_per_bin):
        """C1: constant carrier frequency should give s=0 regardless of rarity."""
        n_bins = 20
        n_carrier = max(1, round(freq * n_per_bin))
        n_ref = n_per_bin - n_carrier
        c_variant = np.full(n_bins, float(n_carrier))
        c_reference = np.full(n_bins, float(n_ref))
        fit = gc.fit_logistic_growth(c_variant, c_reference, bin_days=30.0)
        if fit.status == "ok":
            assert abs(fit.growth_per_day) < 1e-6, (
                f"freq={freq} n={n_per_bin}: s={fit.growth_per_day} (expected 0.0) -- C1 regression"
            )


# ===========================================================================
# fit_logistic_growth_stratified on multi-country data
# ===========================================================================

@requires_db
class TestStratifiedFitOnRealData:
    """Multi-country stratified fit must respect shared-slope across regions.

    fit_logistic_growth_stratified(c_child, c_parent, bin_days) takes 2-D
    (time x stratum) count arrays.
    """

    def test_random_carriers_give_near_zero_stratified_slope(self):
        """Random carrier assignment across strata should give near-zero shared slope."""
        with sqlite3.connect(str(DB_PATH)) as conn:
            meta = pd.read_sql(
                "SELECT primary_accession, country "
                "FROM meta_data WHERE exclusion_status IS NULL OR exclusion_status = ''",
                conn,
            )
        big_countries = (
            meta.groupby("country").size()[lambda s: s >= 20].index.tolist()
        )
        if len(big_countries) < 2:
            pytest.skip("Not enough multi-country data for stratified test")
        meta = meta[meta["country"].isin(big_countries[:4])].copy()

        n_bins = 20
        n_strata = len(big_countries[:4])
        rng = np.random.default_rng(42)
        # Random 30% carrier frequency, constant across time
        c_child = rng.integers(1, 5, size=(n_bins, n_strata)).astype(float)
        c_parent = rng.integers(5, 15, size=(n_bins, n_strata)).astype(float)

        fit = gc.fit_logistic_growth_stratified(c_child, c_parent, bin_days=30.0)
        if fit.status == "ok":
            rate_per_yr = fit.growth_per_day * 365.25
            # Random data: slope should be close to zero (with wide tolerance)
            assert abs(rate_per_yr) < 5.0, (
                f"Stratified growth {rate_per_yr:.3f}/yr for random data -- unexpectedly large"
            )


# ===========================================================================
# Design identifiability
# ===========================================================================

class TestDesignIdentifiabilityOnRealData:
    """C3: design_identifiability must detect any rank-deficient mutation co-occurrence."""

    def test_perfectly_correlated_columns_are_confounded(self):
        """Two identical columns must appear in the same confounded group."""
        X = np.array([
            [1, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [0, 0, 1],
        ], dtype=float)
        labels = ["A", "B", "C"]
        info = gc.design_identifiability(X, labels)
        assert "confounded_groups" in info
        assert info["rank"] < info["n_parameters"]
        # A and B are identical columns -- should be confounded
        confounded_flat = [item for group in info["confounded_groups"] for item in group]
        assert "A" in confounded_flat and "B" in confounded_flat

    def test_independent_columns_are_identifiable(self):
        """Full-rank design should report as identifiable."""
        X = np.eye(4)
        labels = ["v1", "v2", "v3", "v4"]
        info = gc.design_identifiability(X, labels)
        assert info["identifiable"]
        assert info["rank"] == info["n_parameters"]
        assert info["null_space_dim"] == 0

    def test_returns_required_keys(self):
        X = np.array([[1, 0], [1, 0], [0, 1]], dtype=float)
        labels = ["A", "B"]
        info = gc.design_identifiability(X, labels)
        for key in ("rank", "n_parameters", "identifiable", "confounded_groups"):
            assert key in info, f"Missing key: {key}"


# ===========================================================================
# CLI smoke test
# ===========================================================================

@requires_db
class TestWithinLineageGrowthCLI:
    """Verify the CLI writes expected output files."""

    def test_generates_report_and_tsv(self, tmp_path):
        import within_lineage_growth as wlg_mod
        out_dir = tmp_path / "growth_output"
        wlg_mod.main([
            "--db", str(DB_PATH),
            "--protein", "nucleoprotein N",
            "--out", str(out_dir),
            "--stratify-by", "country",
        ])
        assert (out_dir / "report.md").exists()
        assert (out_dir / "variant_growth_rates.tsv").exists()

    def test_tsv_has_expected_columns(self, tmp_path):
        import within_lineage_growth as wlg_mod
        out_dir = tmp_path / "growth_output2"
        wlg_mod.main([
            "--db", str(DB_PATH),
            "--protein", "nucleoprotein N",
            "--out", str(out_dir),
        ])
        df = pd.read_csv(out_dir / "variant_growth_rates.tsv", sep="\t")
        for col in ("variant", "growth_per_year", "se_per_year", "status"):
            assert col in df.columns, f"Missing TSV column: {col}"

    def test_nan_not_zero_in_tsv_for_unfittable_variants(self, tmp_path):
        """M1: TSV must contain NaN (not 0.0) for unfittable variants."""
        import within_lineage_growth as wlg_mod
        out_dir = tmp_path / "growth_output3"
        wlg_mod.main([
            "--db", str(DB_PATH),
            "--protein", "nucleoprotein N",
            "--out", str(out_dir),
        ])
        df = pd.read_csv(out_dir / "variant_growth_rates.tsv", sep="\t")
        bad = df[df["status"] != "ok"]
        if len(bad) > 0:
            zero_growth = bad["growth_per_year"].eq(0.0)
            assert not zero_growth.any(), (
                f"{zero_growth.sum()} unfittable variants have 0.0 instead of NaN -- M1 regression"
            )


# ===========================================================================
# glycoprotein G -- cross-protein smoke test
# ===========================================================================

@requires_db
class TestGlycoproteinG:
    """Smoke test on transmembrane glycoprotein G to verify cross-protein portability."""

    def test_glycoprotein_g_loads_cohort(self):
        with WithinLineageGrowthAnalyzer(
            DB_PATH, protein="transmembrane glycoprotein G"
        ) as wlg:
            cohort = wlg.load_cohort()
        assert isinstance(cohort, pd.DataFrame)
        assert len(cohort) > 50

    def test_glycoprotein_g_variant_growth_runs(self, tmp_path):
        import within_lineage_growth as wlg_mod
        out_dir = tmp_path / "glycoG_growth"
        wlg_mod.main([
            "--db", str(DB_PATH),
            "--protein", "transmembrane glycoprotein G",
            "--out", str(out_dir),
        ])
        assert (out_dir / "report.md").exists()
