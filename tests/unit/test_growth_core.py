"""Tests for scripts/growth_core.py.

These exist because the estimator they replace passed every smoke test and was
still wrong. ``fit_variant_renewal_fitness`` ran, converged, produced plausible
looking figures, and reported a 46-sigma transmission disadvantage for a variant
held at a constant 1% frequency. Nothing short of feeding it data with a known
answer would have caught that, so that is what the bulk of this module does.

The calibration tests are the point:

  constant frequency  a variant that is not changing must return growth 0, at
                      every frequency and every sample size. The old estimator
                      returned -1.40 at 1% because its first time bin asserted a
                      50:50 expectation, so rarity leaked into the slope.
  known sweep         a logistic sweep at a known rate must come back at that
                      rate, not a fraction of it.
  unfittable input    a singleton has no trajectory and must return NaN. The old
                      code returned 0.0, which ranked lookup failures above
                      every real variant once the rest had been dragged negative.

See tmp_explore/growth_rates/beta_testing/methodology_review.md for the full diagnosis.
"""
from __future__ import annotations

import numpy as np
import pytest

import growth_core as gc


def constant_frequency(freq: float, per_bin: int, n_bins: int = 26):
    """Counts for a variant pinned at ``freq`` for the whole window."""
    return (np.full(n_bins, per_bin * freq, dtype=float),
            np.full(n_bins, per_bin * (1.0 - freq), dtype=float))


def logistic_sweep(rate_per_bin: float, per_bin: int = 1000, n_bins: int = 26,
                   midpoint: float = 15.0):
    """Counts for a variant sweeping logistically at a known rate."""
    t = np.arange(n_bins, dtype=float)
    p = 1.0 / (1.0 + np.exp(-rate_per_bin * (t - midpoint)))
    return np.round(per_bin * p), np.round(per_bin * (1.0 - p))


# ---------------------------------------------------------------------------
# Calibration: the tests the previous estimator failed
# ---------------------------------------------------------------------------
class TestConstantFrequencyIsZeroGrowth:
    """A variant that is not changing frequency has growth rate zero."""

    @pytest.mark.parametrize("freq", [0.5, 0.2, 0.05, 0.01, 0.002])
    @pytest.mark.parametrize("per_bin", [50, 500, 5000])
    def test_growth_is_zero(self, freq, per_bin):
        c_m, c_r = constant_frequency(freq, per_bin)
        fit = gc.fit_logistic_growth(c_m, c_r, 14.0)
        assert fit.status == "ok"
        assert fit.growth_per_day == pytest.approx(0.0, abs=1e-9)
        assert fit.log_R_ratio == pytest.approx(0.0, abs=1e-9)

    def test_rarity_does_not_leak_into_growth(self):
        """The regression test for finding C1, stated directly.

        The old estimator's answer was a function of frequency: 50% -> 0.000,
        5% -> -0.357, 1% -> -1.398. Here all three must agree.
        """
        growths = [
            gc.fit_logistic_growth(*constant_frequency(f, 5000), 14.0).growth_per_day
            for f in (0.5, 0.05, 0.01)
        ]
        assert max(abs(g) for g in growths) < 1e-9

    def test_bias_does_not_survive_more_data(self):
        """Whatever bias remains must shrink with N, i.e. be noise not structure."""
        wide = [gc.fit_logistic_growth(*constant_frequency(0.01, n), 14.0)
                for n in (50, 500, 5000)]
        assert all(abs(f.growth_per_day) < 1e-9 for f in wide)
        # Standard errors shrink as sqrt(N), so a real bias would show up as a
        # growing z-score. The old code's did: -1.288 at N=50 to -1.398 at N=5000.
        z = [abs(f.growth_per_day) / f.se_per_day for f in wide]
        assert max(z) < 0.01


class TestKnownSweepIsRecovered:
    @pytest.mark.parametrize("rate_per_bin", [0.1, 0.3, 0.6])
    def test_recovers_rate_within_one_percent(self, rate_per_bin):
        c_m, c_r = logistic_sweep(rate_per_bin)
        fit = gc.fit_logistic_growth(c_m, c_r, 14.0)
        expected_per_day = rate_per_bin / 14.0
        assert fit.status == "ok"
        assert fit.growth_per_day == pytest.approx(expected_per_day, rel=0.01)

    def test_declining_variant_gives_negative_growth(self):
        c_m, c_r = logistic_sweep(-0.3)
        fit = gc.fit_logistic_growth(c_m, c_r, 14.0)
        assert fit.growth_per_day == pytest.approx(-0.3 / 14.0, rel=0.01)

    def test_growth_sign_survives_a_low_starting_frequency(self):
        """A rare *but growing* variant must score positive.

        Under the old estimator rarity dominated, so genuinely expanding rare
        variants came back negative. This is the case that mattered
        scientifically and it must now come out right.
        """
        c_m, c_r = logistic_sweep(0.3, per_bin=1000, midpoint=34.0)
        fit = gc.fit_logistic_growth(c_m, c_r, 14.0)
        assert fit.growth_per_day > 0
        assert fit.growth_per_day == pytest.approx(0.3 / 14.0, rel=0.05)
        assert c_m.sum() / (c_m.sum() + c_r.sum()) < 0.1   # genuinely rare overall


# ---------------------------------------------------------------------------
# Unfittable input must be NaN, never zero
# ---------------------------------------------------------------------------
class TestUnfittableReturnsNaN:
    def test_no_carriers(self):
        fit = gc.fit_logistic_growth(np.zeros(26), np.full(26, 700.0), 14.0)
        assert fit.status == "no_carriers"
        assert np.isnan(fit.growth_per_day)

    @pytest.mark.parametrize("carriers", [1, 2, 9])
    def test_single_bin_variant_regardless_of_count(self, carriers):
        c_m = np.zeros(26)
        c_m[20] = carriers
        fit = gc.fit_logistic_growth(c_m, np.full(26, 700.0), 14.0)
        assert fit.status == "variant_in_too_few_bins"
        assert np.isnan(fit.growth_per_day)
        assert fit.carriers == carriers

    def test_short_span_is_rejected(self):
        c_m = np.array([1.0, 2.0])
        c_r = np.array([100.0, 100.0])
        fit = gc.fit_logistic_growth(c_m, c_r, 14.0, min_bins=2)
        assert fit.status == "span_too_short"
        assert np.isnan(fit.growth_per_day)

    def test_nan_is_distinguishable_from_neutral(self):
        """The M1 regression test: unfittable must not look like neutral."""
        unfittable = gc.fit_logistic_growth(np.zeros(26), np.full(26, 700.0), 14.0)
        neutral = gc.fit_logistic_growth(*constant_frequency(0.05, 500), 14.0)
        assert np.isnan(unfittable.growth_per_day)
        assert neutral.growth_per_day == pytest.approx(0.0, abs=1e-9)
        assert not np.isnan(neutral.growth_per_day)


class TestFitValidation:
    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match="same shape"):
            gc.fit_logistic_growth(np.zeros(5), np.zeros(6), 14.0)

    def test_negative_counts_raise(self):
        with pytest.raises(ValueError, match="non-negative"):
            gc.fit_logistic_growth(np.array([-1.0, 2.0, 3.0]), np.ones(3), 14.0)

    def test_separation_is_flagged_not_silently_infinite(self):
        c_m = np.concatenate([np.zeros(13), np.full(13, 500.0)])
        c_r = np.concatenate([np.full(13, 500.0), np.zeros(13)])
        fit = gc.fit_logistic_growth(c_m, c_r, 14.0)
        assert np.isfinite(fit.growth_per_day)
        assert fit.separated
        assert fit.growth_per_day > 0

    def test_intercept_absorbs_level_not_slope(self):
        """Two variants with the same trajectory shape at different levels must
        get the same slope and different intercepts."""
        lo = gc.fit_logistic_growth(*logistic_sweep(0.3, per_bin=1000), 14.0)
        c_m, c_r = logistic_sweep(0.3, per_bin=1000)
        # Same shape, ten-fold rarer: scale the variant arm only.
        hi = gc.fit_logistic_growth(c_m / 10.0, c_r, 14.0)
        assert hi.growth_per_day == pytest.approx(lo.growth_per_day, rel=0.05)
        assert hi.intercept < lo.intercept


# ---------------------------------------------------------------------------
# Generation interval
# ---------------------------------------------------------------------------
class TestGenerationInterval:
    def test_moments_round_trip(self):
        gi = gc.GenerationInterval(3.2, 1.3)
        assert gi.shape * gi.scale == pytest.approx(3.2)
        assert np.sqrt(gi.shape) * gi.scale == pytest.approx(1.3)

    def test_fortnightly_bins_put_all_mass_in_lag_zero(self):
        """The C4 diagnosis, as an assertion.

        The old code spread weights 0.75/0.20/0.05 over lags 1-3 (14-56 days)
        and excluded lag 0. For this generation interval lag 0 holds all of the
        mass, so the kernel it used and the distribution it documented did not
        overlap at all. That also means a renewal model is degenerate at this
        bin width, which is why growth_core estimates a growth rate instead.
        """
        w = gc.H3N2_GENERATION_INTERVAL.discretise(14.0, n_lags=6)
        assert w[0] == pytest.approx(1.0, abs=1e-5)
        assert w[1:].sum() < 1e-5

    def test_short_bins_give_a_real_distribution(self):
        w = gc.H3N2_GENERATION_INTERVAL.discretise(1.0, n_lags=8)
        assert w.sum() == pytest.approx(1.0)
        assert w[0] < 0.1               # little same-day transmission
        assert w.argmax() in (2, 3)     # peaks near the 3.2-day mean

    def test_log_R_ratio_is_zero_at_zero_growth(self):
        assert gc.H3N2_GENERATION_INTERVAL.log_R_ratio(0.0) == pytest.approx(0.0)

    def test_log_R_ratio_approximates_growth_per_generation(self):
        gi = gc.H3N2_GENERATION_INTERVAL
        for r in (0.001, 0.005, 0.01):
            assert gi.log_R_ratio(r) == pytest.approx(gi.mean_days * r, rel=0.02)

    def test_log_R_ratio_is_monotone_and_signed(self):
        gi = gc.H3N2_GENERATION_INTERVAL
        vals = [gi.log_R_ratio(r) for r in (-0.05, -0.01, 0.0, 0.01, 0.05)]
        assert all(b > a for a, b in zip(vals, vals[1:]))
        assert vals[0] < 0 < vals[-1]

    def test_fortnightly_rate_is_not_a_reproduction_ratio(self):
        """Quantifies the ~4.4x overstatement in the old Figure 6B labels."""
        gi = gc.H3N2_GENERATION_INTERVAL
        per_fortnight = 1.0067          # the old fitted f_v for lineage K
        per_day = per_fortnight / 14.0
        assert np.exp(gi.log_R_ratio(per_day)) == pytest.approx(1.24, abs=0.05)
        assert np.exp(per_fortnight) == pytest.approx(2.74, abs=0.05)

    def test_bad_bin_width_raises(self):
        with pytest.raises(ValueError, match="positive"):
            gc.H3N2_GENERATION_INTERVAL.discretise(0.0)


# ---------------------------------------------------------------------------
# Multinomial lineage model
# ---------------------------------------------------------------------------
class TestMultinomialLogisticGrowth:
    @staticmethod
    def _two_lineage_tensor(rate_per_bin, n_bins=30, per_bin=800, n_regions=2):
        t = np.arange(n_bins, dtype=float)
        p = 1.0 / (1.0 + np.exp(-rate_per_bin * (t - n_bins / 2)))
        C = np.zeros((n_bins, n_regions, 2))
        for r in range(n_regions):
            C[:, r, 1] = np.round(per_bin * p)
            C[:, r, 0] = np.round(per_bin * (1.0 - p))
        return C

    def test_recovers_a_known_two_lineage_advantage(self):
        C = self._two_lineage_tensor(0.3)
        s, se, logR, n_cells = gc.fit_multinomial_logistic_growth(C, 14.0)
        assert s[0] == 0.0
        assert s[1] == pytest.approx(0.3 / 14.0, rel=0.05)
        assert se[1] > 0
        assert n_cells > 0

    def test_reference_lineage_is_pinned_at_zero(self):
        C = self._two_lineage_tensor(0.3)
        s, _, logR, _ = gc.fit_multinomial_logistic_growth(C, 14.0, reference_index=0)
        assert s[0] == 0.0
        assert logR[0] == pytest.approx(0.0)

    def test_region_specific_prevalence_is_not_read_as_growth(self):
        """A lineage that merely dominates one region must not score as growing.

        The old floored-Lambda likelihood could not tell these apart, because an
        absent lineage was given a non-zero expected share everywhere.
        """
        n_bins = 30
        C = np.zeros((n_bins, 2, 2))
        C[:, 0, 0] = 900.0      # region 0: lineage 0 dominant, flat
        C[:, 0, 1] = 100.0
        C[:, 1, 0] = 100.0      # region 1: lineage 1 dominant, also flat
        C[:, 1, 1] = 900.0
        s, _, _, _ = gc.fit_multinomial_logistic_growth(C, 14.0)
        assert abs(s[1]) < 1e-3

    def test_sparse_cells_are_excluded_not_floored(self):
        C = self._two_lineage_tensor(0.3)
        C[0, :, :] = 0.0        # an empty leading bin, as every real series has
        s, _, _, _ = gc.fit_multinomial_logistic_growth(C, 14.0, min_cell_total=3)
        assert s[1] == pytest.approx(0.3 / 14.0, rel=0.05)

    def test_rejects_bad_shape(self):
        with pytest.raises(ValueError, match="time, region, lineage"):
            gc.fit_multinomial_logistic_growth(np.zeros((4, 4)), 14.0)

    def test_rejects_when_no_cell_is_informative(self):
        with pytest.raises(ValueError, match="min_cell_total"):
            gc.fit_multinomial_logistic_growth(np.zeros((5, 2, 3)), 14.0)


# ---------------------------------------------------------------------------
# Design identifiability
# ---------------------------------------------------------------------------
class TestDesignIdentifiability:
    def test_flags_an_underdetermined_design(self):
        X = np.array([[1, 1, 0], [1, 1, 1]], dtype=float)
        rep = gc.design_identifiability(X, ["a", "b", "c"])
        assert rep["n_observations"] == 2
        assert rep["n_parameters"] == 3
        assert not rep["identifiable"]
        assert rep["null_space_dim"] == 1

    def test_identifies_perfectly_confounded_columns(self):
        # 'a' and 'b' occur in exactly the same rows, so they cannot be separated.
        X = np.array([[1, 1, 0], [1, 1, 1], [0, 0, 1]], dtype=float)
        rep = gc.design_identifiability(X, ["a", "b", "c"])
        assert rep["largest_confounded_group"] == 2
        assert sorted(rep["confounded_groups"][0]) == ["a", "b"]

    def test_full_rank_design_is_identifiable(self):
        rep = gc.design_identifiability(np.eye(4), list("abcd"))
        assert rep["identifiable"]
        assert rep["null_space_dim"] == 0
        assert rep["confounded_groups"] == []


# ---------------------------------------------------------------------------
# Correlations with uncertainty
# ---------------------------------------------------------------------------
class TestCorrelationWithCI:
    def test_recovers_a_strong_correlation_with_a_bracketing_ci(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=300)
        y = 0.8 * x + rng.normal(scale=0.6, size=300)
        out = gc.correlation_with_ci(x, y, n_boot=400)
        assert out["ci_low"] < out["estimate"] < out["ci_high"]
        assert out["ci_low"] > 0

    def test_ci_brackets_zero_for_noise(self):
        rng = np.random.default_rng(1)
        out = gc.correlation_with_ci(rng.normal(size=200), rng.normal(size=200),
                                     n_boot=400)
        assert out["ci_low"] < 0 < out["ci_high"]

    def test_nans_are_dropped_and_counted_not_propagated(self):
        """The M3 regression test: one NaN must not blank the whole result."""
        x = np.arange(50.0)
        y = 2 * x
        y[7] = np.nan
        out = gc.correlation_with_ci(x, y, n_boot=200)
        assert np.isfinite(out["estimate"])
        assert out["n"] == 49
        assert out["n_dropped"] == 1

    def test_cluster_bootstrap_widens_the_interval(self):
        """Nested lineages are not independent rows, and the CI must say so."""
        rng = np.random.default_rng(2)
        groups = np.repeat(np.arange(6), 60)
        offset = rng.normal(scale=2.0, size=6)[groups]
        x = rng.normal(size=360) + offset
        y = x + rng.normal(scale=0.5, size=360)
        naive = gc.correlation_with_ci(x, y, n_boot=500)
        clustered = gc.correlation_with_ci(x, y, groups=groups, n_boot=500)
        naive_w = naive["ci_high"] - naive["ci_low"]
        clustered_w = clustered["ci_high"] - clustered["ci_low"]
        assert clustered_w > naive_w

    def test_spearman_method_available(self):
        x = np.arange(60.0)
        out = gc.correlation_with_ci(x, np.exp(x / 20.0), method="spearman", n_boot=200)
        assert out["estimate"] == pytest.approx(1.0)

    def test_too_few_points_returns_nan_not_an_exception(self):
        out = gc.correlation_with_ci([1.0, 2.0], [1.0, 2.0], n_boot=50)
        assert np.isnan(out["estimate"])


class TestWeightedPearson:
    def test_matches_numpy_when_unweighted(self):
        rng = np.random.default_rng(3)
        x, y = rng.normal(size=100), rng.normal(size=100)
        assert gc.weighted_pearson(x, y) == pytest.approx(np.corrcoef(x, y)[0, 1])

    def test_weights_shift_the_estimate_toward_the_weighted_points(self):
        x = np.array([0.0, 1.0, 2.0, 3.0])
        y = np.array([0.0, 1.0, 2.0, -9.0])
        # Down-weighting the outlier should push the correlation up toward 1.
        w = np.array([1.0, 1.0, 1.0, 1e-6])
        assert gc.weighted_pearson(x, y, w) > gc.weighted_pearson(x, y)

    def test_zero_variance_returns_nan(self):
        assert np.isnan(gc.weighted_pearson(np.ones(10), np.arange(10.0)))


# ---------------------------------------------------------------------------
# Collection dates
# ---------------------------------------------------------------------------
class TestParseCollectionDate:
    """The V-gTK database stores dates in ten shapes; two of them are traps."""

    @pytest.mark.parametrize("raw,expected", [
        ("2024-09-04", (2024, 9, 4)),
        ("04-Sep-2024", (2024, 9, 4)),
        ("22-Aug-2024", (2024, 8, 22)),
        ("1-Jan-2020", (2020, 1, 1)),
    ])
    def test_day_precision(self, raw, expected):
        d, precision = gc.parse_collection_date(raw)
        assert precision == "day"
        assert (d.year, d.month, d.day) == expected

    def test_year_only_is_flagged_not_silently_january_first(self):
        """``pd.to_datetime("2025")`` gives 2025-01-01 without complaint.

        Lineage K's whole sampled span is 401 days, so a bare year placed on
        1 January puts a sequence six months from where it belongs. The
        precision label is what lets the caller refuse it.
        """
        d, precision = gc.parse_collection_date("2025")
        assert precision == "year"
        assert d.month == 7 and d.day == 1      # midpoint, not 1 January

    @pytest.mark.parametrize("raw", ["2022-06", "Feb-2022"])
    def test_month_precision_is_mid_month(self, raw):
        d, precision = gc.parse_collection_date(raw)
        assert precision == "month"
        assert 10 <= d.day <= 20

    @pytest.mark.parametrize("raw", [
        "unknown", "", None, "Jan-2020/Feb-2020", "01-Jan-2020/02-Jan-2020",
        "13/04/2024", "2024-13-01", "2024-02-30", "99-Zzz-2024",
    ])
    def test_unusable_values_are_rejected(self, raw):
        d, precision = gc.parse_collection_date(raw)
        assert d is None and precision == "none"

    def test_precision_rank_orders_day_before_month_before_year(self):
        assert (gc.PRECISION_RANK["day"] < gc.PRECISION_RANK["month"]
                < gc.PRECISION_RANK["year"] < gc.PRECISION_RANK["none"])


class TestPlausibleDateMask:
    def test_rejects_spreadsheet_serial_corruption(self):
        """1899-12-30 is Excel day zero; 1905-07-13 is serial 2021.

        These parse cleanly as day-precision dates, so only a plausibility
        window catches them. One of them shifted lineage J.2's bin grid by 116
        years.
        """
        # Passed as raw strings: 1480-10-07 is outside the nanosecond range,
        # so the mask has to survive input pandas cannot vectorise.
        dates = ["1899-12-30", "1905-07-13", "1480-10-07",
                 "2024-01-01", "2025-06-01"]
        mask = gc.plausible_date_mask(dates, min_date="2000-01-01",
                                      max_date="2026-12-31")
        assert list(mask) == [False, False, False, True, True]

    def test_out_of_range_date_parses_but_the_window_rejects_it(self):
        """Division of labour: the parser judges shape, the window judges sense.

        ``1480-10-07`` is a well-formed day-precision date, so the parser
        returns it. Only the plausibility window knows it cannot be an influenza
        collection date.
        """
        d, precision = gc.parse_collection_date("1480-10-07")
        assert precision == "day" and d is not None
        mask = gc.plausible_date_mask([d], min_date="2000-01-01",
                                      max_date="2026-12-31")
        assert not mask[0]

    def test_group_window_rejects_a_date_wrong_for_its_lineage(self):
        """2002 is plausible for H3N2 but not for a lineage that emerged in 2022."""
        import pandas as pd
        dates = list(pd.date_range("2022-06-01", periods=40, freq="14D")) \
            + [pd.Timestamp("2002-12-06")]
        groups = ["J"] * 41
        mask = gc.plausible_date_mask(dates, groups, min_date="2000-01-01",
                                      max_date="2026-12-31",
                                      group_margin_days=730.0)
        assert mask[:40].all()
        assert not mask[40]

    def test_group_window_keeps_a_genuine_early_tail(self):
        import pandas as pd
        dates = list(pd.date_range("2016-06-01", periods=60, freq="14D"))
        mask = gc.plausible_date_mask(dates, ["C.1"] * 60, min_date="2000-01-01",
                                      max_date="2026-12-31")
        assert mask.all()

    def test_small_groups_are_not_trimmed(self):
        import pandas as pd
        dates = pd.to_datetime(["2020-01-01", "2024-01-01", "2025-01-01"])
        mask = gc.plausible_date_mask(dates, ["X"] * 3, min_date="2000-01-01",
                                      max_date="2026-12-31")
        assert mask.all()


# ---------------------------------------------------------------------------
# Meta-analysis
# ---------------------------------------------------------------------------
class TestMetaAnalyse:
    def test_agreeing_strata_give_the_inverse_variance_mean(self):
        out = gc.meta_analyse([0.10, 0.10, 0.10], [0.01, 0.01, 0.01])
        assert out["estimate"] == pytest.approx(0.10)
        assert out["tau2"] == pytest.approx(0.0)
        assert out["se"] == pytest.approx(0.01 / np.sqrt(3), rel=1e-6)

    def test_disagreeing_strata_widen_the_interval(self):
        """The fix for the over-precise lineage CIs.

        Cochran's Q rejected homogeneity for 6 of 10 lineages, some at
        p < 1e-13, yet the fixed-effect pool reported intervals as narrow as
        +/-0.0014/day.
        """
        est, se = [0.05, 0.20, 0.35], [0.01, 0.01, 0.01]
        fixed = gc.meta_analyse(est, se, method="fixed")
        random = gc.meta_analyse(est, se, method="random")
        assert random["se"] > 5 * fixed["se"]
        assert random["tau2"] > 0
        assert random["i2"] > 0.9

    def test_i_squared_is_zero_when_strata_agree(self):
        out = gc.meta_analyse([0.1, 0.1, 0.1], [0.02, 0.02, 0.02])
        assert out["i2"] == pytest.approx(0.0)

    def test_single_stratum_passes_through(self):
        out = gc.meta_analyse([0.2], [0.03])
        assert out["estimate"] == pytest.approx(0.2)
        assert out["se"] == pytest.approx(0.03)
        assert out["k"] == 1

    def test_empty_input_is_nan_not_an_exception(self):
        out = gc.meta_analyse([], [])
        assert np.isnan(out["estimate"])
        assert out["k"] == 0

    def test_non_finite_strata_are_dropped(self):
        out = gc.meta_analyse([0.1, np.nan, 0.1], [0.01, 0.01, np.nan])
        assert out["k"] == 1
        assert out["estimate"] == pytest.approx(0.1)


# ---------------------------------------------------------------------------
# Stratified fit
# ---------------------------------------------------------------------------
class TestStratifiedGrowthFit:
    """One shared slope, per-stratum intercepts. Replaces meta-analysing
    per-stratum point estimates, which with two strata produced intervals driven
    by an unestimable tau^2 rather than by the data."""

    @staticmethod
    def _panel(rate_per_bin, baselines, per_bin=1000, n_bins=26, mid=13.0):
        t = np.arange(n_bins, dtype=float)
        cc, cp = [], []
        for b in baselines:
            p = 1.0 / (1.0 + np.exp(-(rate_per_bin * (t - mid) + b)))
            cc.append(np.round(per_bin * p))
            cp.append(np.round(per_bin * (1.0 - p)))
        return np.array(cc).T, np.array(cp).T

    def test_recovers_a_shared_slope_across_different_baselines(self):
        cc, cp = self._panel(0.3, [-2.0, 0.0, 2.0])
        fit = gc.fit_logistic_growth_stratified(cc, cp, 14.0)
        assert fit.status == "ok"
        assert fit.growth_per_day == pytest.approx(0.3 / 14.0, rel=0.01)
        assert fit.n_strata == 3

    def test_baseline_differences_do_not_leak_into_the_slope(self):
        flat = gc.fit_logistic_growth_stratified(*self._panel(0.3, [0.0, 0.0, 0.0]), 14.0)
        spread = gc.fit_logistic_growth_stratified(*self._panel(0.3, [-3.0, 0.0, 3.0]), 14.0)
        assert spread.growth_per_day == pytest.approx(flat.growth_per_day, rel=0.02)

    @staticmethod
    def _two_strata_differing_slopes(r0=0.20, r1=0.40, per_bin=1000, n_bins=26):
        t = np.arange(n_bins, dtype=float)
        cc, cp = [], []
        for r in (r0, r1):
            p = 1.0 / (1.0 + np.exp(-r * (t - n_bins / 2)))
            cc.append(np.round(per_bin * p))
            cp.append(np.round(per_bin * (1.0 - p)))
        return np.array(cc).T, np.array(cp).T

    def test_agreeing_strata_give_the_same_answer_as_pooling(self):
        """When the strata agree, tau^2 is 0 and random effects reduces to fixed,
        so the two approaches must coincide. This is the case where nothing was
        wrong, and it is here to bound the claim made by the next test."""
        cc, cp = self._panel(0.3, [-2.0, 2.0])
        strat = gc.fit_logistic_growth_stratified(cc, cp, 14.0)
        per = [gc.fit_logistic_growth(cc[:, r], cp[:, r], 14.0) for r in range(2)]
        pooled = gc.meta_analyse([f.growth_per_day for f in per],
                                 [f.se_per_day for f in per], method="random")
        assert strat.heterogeneity_p > 0.5
        assert strat.se_per_day == pytest.approx(pooled["se"], rel=0.05)

    def test_disagreeing_strata_blow_up_random_effects_but_not_the_stratified_fit(self):
        """The defect this exists to fix, quantified.

        Two regions with genuinely different slopes is the real-data situation
        (fitted heterogeneity p < 0.001 for three of the seven transitions). With
        k=2, DerSimonian-Laird estimates tau^2 from one degree of freedom and the
        interval balloons; the reported symptom was J.2 -> K coming back at
        [-0.030, +0.214], spanning zero for the fastest lineage in the dataset.
        The stratified fit's interval stays proportionate to the data.
        """
        cc, cp = self._two_strata_differing_slopes()
        strat = gc.fit_logistic_growth_stratified(cc, cp, 14.0)
        per = [gc.fit_logistic_growth(cc[:, r], cp[:, r], 14.0) for r in range(2)]
        pooled = gc.meta_analyse([f.growth_per_day for f in per],
                                 [f.se_per_day for f in per], method="random")
        assert strat.heterogeneity_p < 1e-6, "the strata must genuinely disagree"
        assert pooled["se"] > 20 * strat.se_per_day
        # and the point estimate stays between the two stratum slopes
        assert min(f.growth_per_day for f in per) <= strat.growth_per_day \
            <= max(f.growth_per_day for f in per)

    def test_a_disagreeing_pool_can_span_zero_where_the_stratified_fit_does_not(self):
        """The exact shape of the reported bug: both strata clearly positive, yet
        the random-effects interval covers zero."""
        cc, cp = self._two_strata_differing_slopes(0.05, 0.60)
        strat = gc.fit_logistic_growth_stratified(cc, cp, 14.0)
        per = [gc.fit_logistic_growth(cc[:, r], cp[:, r], 14.0) for r in range(2)]
        pooled = gc.meta_analyse([f.growth_per_day for f in per],
                                 [f.se_per_day for f in per], method="random")
        assert all(f.growth_per_day > 0 for f in per)
        assert strat.growth_per_day - 1.96 * strat.se_per_day > 0, (
            "the stratified interval must stay off zero when every stratum is "
            "clearly positive")
        pooled_lo = pooled["estimate"] - 1.96 * pooled["se"]
        assert pooled_lo < strat.growth_per_day - 1.96 * strat.se_per_day

    def test_heterogeneity_fires_only_when_slopes_differ(self):
        same = gc.fit_logistic_growth_stratified(*self._panel(0.3, [-1.0, 1.0]), 14.0)
        assert same.heterogeneity_p > 0.01
        t = np.arange(26, dtype=float)
        cc, cp = [], []
        for r in (0.1, 0.5):
            p = 1.0 / (1.0 + np.exp(-r * (t - 13)))
            cc.append(np.round(1000 * p))
            cp.append(np.round(1000 * (1 - p)))
        diff = gc.fit_logistic_growth_stratified(np.array(cc).T, np.array(cp).T, 14.0)
        assert diff.heterogeneity_p < 1e-6
        assert diff.heterogeneity_df == 1

    def test_single_stratum_matches_the_two_parameter_fit(self):
        t = np.arange(26, dtype=float)
        p = 1.0 / (1.0 + np.exp(-0.3 * (t - 13)))
        cc, cp = np.round(1000 * p), np.round(1000 * (1 - p))
        a = gc.fit_logistic_growth(cc, cp, 14.0)
        b = gc.fit_logistic_growth_stratified(cc, cp, 14.0)
        assert b.growth_per_day == pytest.approx(a.growth_per_day, rel=1e-6)
        assert b.n_strata == 1
        assert b.heterogeneity_df == 0

    def test_strata_without_both_arms_are_dropped(self):
        cc, cp = self._panel(0.3, [0.0, 0.0])
        cc[:, 1] = 0.0            # second stratum never sees the child
        fit = gc.fit_logistic_growth_stratified(cc, cp, 14.0)
        assert fit.n_strata == 1

    def test_no_usable_stratum_is_nan_not_a_crash(self):
        fit = gc.fit_logistic_growth_stratified(np.zeros((5, 2)), np.zeros((5, 2)), 14.0)
        assert fit.status == "no_usable_stratum"
        assert np.isnan(fit.growth_per_day)

    def test_scale_invariance_across_strata(self):
        cc, cp = self._panel(0.3, [-1.0, 1.0], per_bin=200)
        a = gc.fit_logistic_growth_stratified(cc, cp, 14.0, slope_prior_sd=0.0)
        b = gc.fit_logistic_growth_stratified(cc * 7, cp * 7, 14.0, slope_prior_sd=0.0)
        assert b.growth_per_day == pytest.approx(a.growth_per_day, rel=1e-8)

    def test_non_finite_and_negative_counts_are_refused(self):
        cc, cp = self._panel(0.3, [0.0])
        bad = cc.copy(); bad[3, 0] = np.nan
        with pytest.raises(ValueError, match="finite"):
            gc.fit_logistic_growth_stratified(bad, cp, 14.0)
        bad = cc.copy(); bad[3, 0] = -1.0
        with pytest.raises(ValueError, match="non-negative"):
            gc.fit_logistic_growth_stratified(bad, cp, 14.0)

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="same shape"):
            gc.fit_logistic_growth_stratified(np.zeros((5, 2)), np.zeros((5, 3)), 14.0)

    def test_one_bin_stratum_is_lossless_to_drop(self):
        """Why min_bins_per_stratum is 2 and not 3.

        Once a stratum carries its own intercept, a stratum observed in a single
        bin contributes nothing to the shared slope: the intercept absorbs it
        exactly. So excluding it must leave both the estimate and its standard
        error untouched, which is what makes 2 the principled threshold. At 3 the
        exclusion discarded genuine slope information.
        """
        t = np.arange(26, dtype=float)
        p = 1.0 / (1.0 + np.exp(-0.3 * (t - 13)))
        a, b = np.round(1000 * p), np.round(1000 * (1 - p))
        alone = gc.fit_logistic_growth_stratified(a[:, None], b[:, None], 14.0)

        cc = np.zeros((26, 2)); cp = np.zeros((26, 2))
        cc[:, 0], cp[:, 0] = a, b
        cc[7, 1], cp[7, 1] = 500, 500        # one bin only
        withit = gc.fit_logistic_growth_stratified(cc, cp, 14.0)

        assert withit.growth_per_day == pytest.approx(alone.growth_per_day, rel=1e-9)
        assert withit.se_per_day == pytest.approx(alone.se_per_day, rel=1e-9)

    def test_two_bin_stratum_is_kept_and_moves_the_estimate(self):
        t = np.arange(26, dtype=float)
        p = 1.0 / (1.0 + np.exp(-0.3 * (t - 13)))
        a, b = np.round(1000 * p), np.round(1000 * (1 - p))
        alone = gc.fit_logistic_growth_stratified(a[:, None], b[:, None], 14.0)

        cc = np.zeros((26, 2)); cp = np.zeros((26, 2))
        cc[:, 0], cp[:, 0] = a, b
        cc[5, 1], cp[5, 1] = 50, 950
        cc[20, 1], cp[20, 1] = 900, 100
        withit = gc.fit_logistic_growth_stratified(cc, cp, 14.0)

        assert withit.n_strata == 2
        assert withit.growth_per_day != pytest.approx(alone.growth_per_day, rel=1e-6)

    def test_thin_strata_are_merged_rather_than_all_discarded(self):
        """Several one-bin regions can jointly earn one intercept."""
        cc = np.zeros((10, 3)); cp = np.zeros((10, 3))
        for j, tb in enumerate([(1,), (4,), (8,)]):
            for t in tb:
                cc[t, j], cp[t, j] = 40 + 20 * t, 60 - 2 * t
        fit = gc.fit_logistic_growth_stratified(cc, cp, 14.0)
        assert fit.n_strata == 1, "the three one-bin regions should merge into one"
        assert np.isfinite(fit.growth_per_day)
