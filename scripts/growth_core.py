"""Identifiable growth-rate estimators and calibration primitives for V-gTK.

Ported from Mut_acc_PLM_MS with calibration guarantees for within-lineage variant
and lineage-level displacement growth models.

Background & Calibration Guarantees
------------------------------------
The previous estimators in renewal fitness scripts exhibited critical structural issues:
1. Rarity Leakage (Finding C1): A model with no intercept forced variant level/rarity
   into the slope; constant 1% frequency variants returned severe false deleterious
   rates (e.g. s = -1.40 +/- 0.03). Fitting an explicit intercept absorbs initial frequency,
   ensuring the slope strictly measures change in frequency over time.
2. Generation Interval Mismatch (Finding C4): Fortnightly bins (14d) exceed generation
   intervals (~3.2d for flu, variable for RABV). This module provides an explicit
   Gamma-distributed GenerationInterval and Wallinga-Lipsitch Euler-Lotka moment-generating
   conversion to ln(R_variant / R_ref).
3. Identifiability Gates & Non-zero Placeholders (Findings S1, M1): Variants without sufficient
   trajectory (singletons, single time bins, insufficient span) return NaN with a specific
   identifiability status (never 0.0, which made lookup failures look fit). Standard errors
   from Fisher information are retained for inverse-variance weighting.
4. Stratified Estimation (Finding S5): Stratified logistic regression models shared slopes
   across geographic strata with per-stratum intercepts, avoiding Simpson's paradox when
   sampling depth varies across regions, coupled with likelihood-ratio heterogeneity testing.
5. Design Identifiability (Finding C3): Directly diagnoses rank deficiency and confounded
   groups in mutation-presence regressions.
6. Robust Date Parsing (Finding M9): Distinguishes day, month, and bare-year precisions
   with midpoint imputation, and applies two-stage plausibility masking to discard
   spreadsheet serial corruptions (1899, 1905) and lineage temporal outliers.

Designed to be lightweight and executable in the `vgtk` conda environment (pure numpy/math,
no mandatory scipy or statsmodels).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math
import re
from typing import Sequence

import numpy as np

__all__ = [
    "GenerationInterval",
    "H3N2_GENERATION_INTERVAL",
    "RABV_GENERATION_INTERVAL",
    "LogisticGrowthFit",
    "fit_logistic_growth",
    "StratifiedGrowthFit",
    "fit_logistic_growth_stratified",
    "fit_multinomial_logistic_growth",
    "design_identifiability",
    "correlation_with_ci",
    "weighted_pearson",
    "parse_collection_date",
    "plausible_date_mask",
    "meta_analyse",
    "PRECISION_RANK",
    "gammainc_upper",
    "chi2_sf",
]

# ---------------------------------------------------------------------------
# Pure-Python / NumPy statistical primitives (no scipy dependency required)
# ---------------------------------------------------------------------------

def _gamma_series(a: float, x: float) -> float:
    """Lower incomplete gamma via Taylor series P(a, x) = gamma(a, x)/Gamma(a)."""
    if x <= 0:
        return 0.0
    term = 1.0 / a
    total = term
    for n in range(1, 200):
        term *= x / (a + n)
        total += term
        if abs(term) < abs(total) * 1e-15:
            break
    log_gamma_a = math.lgamma(a)
    val = total * math.exp(-x + a * math.log(x) - log_gamma_a)
    return float(np.clip(val, 0.0, 1.0))


def _gamma_continued_fraction(a: float, x: float) -> float:
    """Upper incomplete gamma Q(a, x) = Gamma(a, x)/Gamma(a) via Legendre continued fraction."""
    tiny = 1e-30
    b0 = 0.0
    c0 = 1.0 / tiny
    d0 = 0.0

    b1 = x + 1.0 - a
    c1 = b1
    d1 = 1.0 / b1 if b1 != 0 else 1.0 / tiny

    f = c1
    for m in range(1, 200):
        a_m = -m * (m - a)
        b_m = b1 + 2.0
        d1 = b_m + a_m * d1
        if abs(d1) < tiny:
            d1 = tiny
        d1 = 1.0 / d1

        c1 = b_m + a_m / c1
        if abs(c1) < tiny:
            c1 = tiny

        del_f = c1 * d1
        f *= del_f
        b1 = b_m
        if abs(del_f - 1.0) < 1e-15:
            break

    log_gamma_a = math.lgamma(a)
    val = math.exp(-x + a * math.log(x) - log_gamma_a) * (1.0 / f)
    return float(np.clip(val, 0.0, 1.0))


def gammainc_upper(a: float, x: float) -> float:
    """Normalized upper incomplete gamma Q(a, x) = Gamma(a, x) / Gamma(a)."""
    if x <= 0:
        return 1.0
    if a <= 0:
        return 0.0
    if x < a + 1.0:
        return float(1.0 - _gamma_series(a, x))
    return float(_gamma_continued_fraction(a, x))


def chi2_sf(x: float, df: int | float) -> float:
    """Chi-squared survival function P(X >= x) for X ~ Chi2(df)."""
    if not math.isfinite(x) or x <= 0 or df <= 0:
        return 1.0 if x <= 0 else 0.0
    return gammainc_upper(df / 2.0, x / 2.0)


def _rankdata(a: np.ndarray) -> np.ndarray:
    """Compute fractional ranks (handling ties by average) in pure NumPy."""
    a = np.asarray(a)
    n = a.size
    if n == 0:
        return np.array([])
    order = a.argsort()
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1, dtype=float)

    sorted_a = a[order]
    dup = np.r_[True, sorted_a[1:] != sorted_a[:-1], True]
    dup_idx = np.flatnonzero(dup)
    for i in range(len(dup_idx) - 1):
        start, end = dup_idx[i], dup_idx[i + 1]
        if end - start > 1:
            ranks[order[start:end]] = (start + 1 + end) / 2.0
    return ranks


# ---------------------------------------------------------------------------
# Combining estimates across strata
# ---------------------------------------------------------------------------
def meta_analyse(estimates, standard_errors, method: str = "random") -> dict:
    """Combine per-stratum estimates, defaulting to random effects.

    DerSimonian-Laird random effects adds the between-region variance ``tau^2``
    to each weight, so genuine geographic disagreement widens the interval
    instead of being averaged away (Finding S5).
    """
    g = np.asarray(estimates, dtype=float)
    s = np.asarray(standard_errors, dtype=float)
    ok = np.isfinite(g) & np.isfinite(s) & (s > 0)
    g, s = g[ok], s[ok]
    k = g.size
    out = {"k": int(k), "method": method, "tau2": 0.0, "q": np.nan,
           "q_df": max(k - 1, 0), "i2": np.nan}
    if k == 0:
        return {**out, "estimate": np.nan, "se": np.nan}
    if k == 1:
        return {**out, "estimate": float(g[0]), "se": float(s[0])}

    w = 1.0 / s ** 2
    fixed = float((w * g).sum() / w.sum())
    q = float((w * (g - fixed) ** 2).sum())
    out["q"] = q
    out["i2"] = float(max(0.0, (q - (k - 1)) / q)) if q > 0 else 0.0

    if method == "fixed":
        return {**out, "estimate": fixed, "se": float(np.sqrt(1.0 / w.sum()))}

    denom = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = 0.0 if denom <= 0 else max(0.0, (q - (k - 1)) / denom)
    out["tau2"] = float(tau2)
    w_re = 1.0 / (s ** 2 + tau2)
    return {**out,
            "estimate": float((w_re * g).sum() / w_re.sum()),
            "se": float(np.sqrt(1.0 / w_re.sum()))}


# ---------------------------------------------------------------------------
# Collection dates & plausibility masks
# ---------------------------------------------------------------------------
_MONTH_ABBR = {m: i + 1 for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun",
     "jul", "aug", "sep", "oct", "nov", "dec"])}

_RE_ISO_DAY = re.compile(r"^(\d{4})-(\d{2})-(\d{2})$")
_RE_DMY_DAY = re.compile(r"^(\d{1,2})-([A-Za-z]{3})-(\d{4})$")
_RE_ISO_MONTH = re.compile(r"^(\d{4})-(\d{2})$")
_RE_MY_MONTH = re.compile(r"^([A-Za-z]{3})-(\d{4})$")
_RE_YEAR = re.compile(r"^(\d{4})$")

PRECISION_RANK = {"day": 0, "month": 1, "year": 2, "none": 3}


def parse_collection_date(raw):
    """Parse one collection_date string into ``(pandas.Timestamp | None, precision)``.

    ``precision`` is one of ``"day"``, ``"month"``, ``"year"`` or ``"none"``.
    Month- and year-precision dates are placed at the midpoint of their period,
    not the first day, so that imputation error is centred rather than
    systematically early (Finding M9). Ambiguous formats return ``"none"``.
    """
    import pandas as pd

    if raw is None:
        return None, "none"
    s = str(raw).strip()
    if not s:
        return None, "none"

    m = _RE_ISO_DAY.match(s)
    if m:
        y, mo, d = (int(g) for g in m.groups())
        try:
            return pd.Timestamp(year=y, month=mo, day=d), "day"
        except ValueError:
            return None, "none"
    m = _RE_DMY_DAY.match(s)
    if m:
        d, mon, y = m.group(1), m.group(2).lower(), m.group(3)
        if mon not in _MONTH_ABBR:
            return None, "none"
        try:
            return pd.Timestamp(year=int(y), month=_MONTH_ABBR[mon], day=int(d)), "day"
        except ValueError:
            return None, "none"
    m = _RE_ISO_MONTH.match(s)
    if m:
        y, mo = int(m.group(1)), int(m.group(2))
        if not 1 <= mo <= 12:
            return None, "none"
        return pd.Timestamp(year=y, month=mo, day=1) + pd.Timedelta(days=14), "month"
    m = _RE_MY_MONTH.match(s)
    if m:
        mon, y = m.group(1).lower(), int(m.group(2))
        if mon not in _MONTH_ABBR:
            return None, "none"
        return pd.Timestamp(year=y, month=_MONTH_ABBR[mon], day=1) + pd.Timedelta(days=14), "month"
    m = _RE_YEAR.match(s)
    if m:
        return pd.Timestamp(year=int(m.group(1)), month=7, day=1), "year"
    return None, "none"


def plausible_date_mask(dates, groups=None, *, min_date, max_date,
                        group_margin_days: float = 730.0,
                        group_quantile: float = 0.01):
    """Boolean mask over ``dates`` keeping only plausible collection dates.

    Applies an absolute date window (removing spreadsheet-serial corruptions such
    as 1899-12-30 or 1905) followed by a per-group robust quantile window to prevent
    isolated outlier sequences from distorting temporal regression slopes.
    """
    import pandas as pd

    raw = pd.Series(list(dates))
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            d = pd.to_datetime(raw, errors="coerce", format="mixed")
    except Exception:
        def _one(v):
            try:
                return pd.Timestamp(v)
            except Exception:
                return pd.NaT
        d = pd.Series([_one(v) for v in raw])
    ok = d.notna() & (d >= pd.Timestamp(min_date)) & (d <= pd.Timestamp(max_date))
    if groups is None:
        return ok.to_numpy()

    g = pd.Series(list(groups), index=d.index)
    margin = pd.Timedelta(days=float(group_margin_days))
    for key, idx in g.groupby(g, observed=True).groups.items():
        sel = d.loc[idx][ok.loc[idx]]
        if sel.size < 10:
            continue
        lo = sel.quantile(group_quantile) - margin
        hi = sel.quantile(1.0 - group_quantile) + margin
        bad = idx[(d.loc[idx] < lo) | (d.loc[idx] > hi)]
        ok.loc[bad] = False
    return ok.to_numpy()


# ---------------------------------------------------------------------------
# Generation interval & Euler-Lotka reproduction number conversion
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class GenerationInterval:
    """A Gamma generation-interval distribution parameterised by mean and SD.

    Enables computing probability masses over discrete time bins and converting
    logistic growth rates into reproduction number ratios ln(R_variant / R_ref)
    via Euler-Lotka (Finding C4).
    """

    mean_days: float = 3.2
    sd_days: float = 1.3

    @property
    def shape(self) -> float:
        return (self.mean_days / self.sd_days) ** 2

    @property
    def scale(self) -> float:
        return self.sd_days ** 2 / self.mean_days

    def discretise(self, bin_days: float, n_lags: int = 6) -> np.ndarray:
        """Probability mass of the generation interval in each lag bin."""
        if bin_days <= 0:
            raise ValueError("bin_days must be positive")
        k, theta = self.shape, self.scale
        edges = np.arange(n_lags + 1, dtype=float) * bin_days
        cdfs = np.array([1.0 - gammainc_upper(k, e / theta) for e in edges])
        mass = np.diff(cdfs)
        total = mass.sum()
        if total <= 0:
            raise ValueError("generation interval has no mass in the requested lags")
        return mass / total

    def log_R_ratio(self, growth_per_day, ref_growth_per_day: float = 0.0):
        """Convert a logistic growth-rate advantage to ``ln(R_variant / R_ref)``.

        Euler-Lotka gives ``1 = R * M_w(-r)`` for exponential growth at rate ``r``.
        For Gamma(shape k, scale theta), M_w(-r) = (1 + theta*r)^-k, so
        ``ln R = k * ln(1 + theta*r)`` (Wallinga & Lipsitch 2007).
        """
        k, theta = self.shape, self.scale
        r_ref = float(ref_growth_per_day)
        r_var = r_ref + np.asarray(growth_per_day, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = k * (np.log1p(theta * r_var) - np.log1p(theta * r_ref))
            out = np.where(1.0 + theta * r_var > 0, out, np.nan)
        return out if getattr(out, "ndim", 0) else float(out)


H3N2_GENERATION_INTERVAL = GenerationInterval(mean_days=3.2, sd_days=1.3)
RABV_GENERATION_INTERVAL = GenerationInterval(mean_days=30.0, sd_days=15.0)


# ---------------------------------------------------------------------------
# Two-parameter binomial logistic growth
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LogisticGrowthFit:
    """Result of one variant-versus-reference growth fit."""

    growth_per_day: float
    se_per_day: float
    log_R_ratio: float
    intercept: float
    n_bins: int
    n_variant_bins: int
    carriers: int
    span_days: float
    converged: bool
    separated: bool
    status: str

    def as_dict(self) -> dict:
        return asdict(self)


_UNFIT = dict(
    growth_per_day=np.nan, se_per_day=np.nan, log_R_ratio=np.nan,
    intercept=np.nan, converged=False, separated=False,
)


def fit_logistic_growth(
    c_variant: Sequence[float],
    c_reference: Sequence[float],
    bin_days: float = 14.0,
    *,
    generation_interval: GenerationInterval = H3N2_GENERATION_INTERVAL,
    slope_prior_sd: float = 0.25,
    min_bins: int = 3,
    min_variant_bins: int = 2,
    min_span_days: float = 28.0,
    max_iter: int = 100,
    tol: float = 1e-9,
) -> LogisticGrowthFit:
    """Fit ``logit(p_t) = a + s*t`` for variant frequency against a reference.

    The intercept ``a`` absorbs the variant's starting frequency, so the slope
    ``s`` measures change in frequency and nothing else (fixes Finding C1).

    Identifiability gates return NaN with a specific status (Finding S1, M1).
    A weak Gaussian prior is placed on the slope only to keep separated series finite.
    """
    c_m = np.asarray(c_variant, dtype=float)
    c_r = np.asarray(c_reference, dtype=float)
    if c_m.shape != c_r.shape:
        raise ValueError("variant and reference count arrays must be the same shape")
    if c_m.ndim != 1:
        raise ValueError("count arrays must be one-dimensional")
    if np.any(c_m < 0) or np.any(c_r < 0):
        raise ValueError("counts must be non-negative")
    if not np.all(np.isfinite(c_m)) or not np.all(np.isfinite(c_r)):
        raise ValueError("counts must be finite (got NaN or inf)")

    carriers = int(round(c_m.sum()))
    n_tot = c_m + c_r
    informative = n_tot > 0
    variant_bins = int(np.count_nonzero(c_m > 0))
    n_inf = int(np.count_nonzero(informative))

    def _fail(status: str) -> LogisticGrowthFit:
        return LogisticGrowthFit(
            n_bins=n_inf, n_variant_bins=variant_bins, carriers=carriers,
            span_days=0.0, status=status, **_UNFIT,
        )

    if carriers == 0:
        return _fail("no_carriers")
    if variant_bins < min_variant_bins:
        return _fail("variant_in_too_few_bins")
    if n_inf < min_bins:
        return _fail("too_few_informative_bins")

    idx = np.flatnonzero(informative)
    span_days = float((idx[-1] - idx[0]) * bin_days)
    if span_days < min_span_days:
        return _fail("span_too_short")

    x = (idx - idx.mean()) * float(bin_days)
    cm, ct = c_m[idx], n_tot[idx]

    prior_prec = 0.0 if not np.isfinite(slope_prior_sd) or slope_prior_sd <= 0         else 1.0 / float(slope_prior_sd) ** 2
    beta = np.zeros(2)
    converged = False
    for _ in range(max_iter):
        eta = beta[0] + beta[1] * x
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30.0, 30.0)))
        resid = cm - ct * p
        grad = np.array([resid.sum(), float(resid @ x)])
        grad[1] -= prior_prec * beta[1]

        w = ct * p * (1.0 - p)
        h00 = -w.sum()
        h01 = -float(w @ x)
        h11 = -float(w @ (x * x)) - prior_prec
        hess = np.array([[h00, h01], [h01, h11]])

        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.solve(hess - 1e-8 * np.eye(2), grad)
        beta = beta - step
        if np.max(np.abs(step)) < tol:
            converged = True
            break

    eta = beta[0] + beta[1] * x
    p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30.0, 30.0)))
    w = ct * p * (1.0 - p)
    fisher = np.array([
        [w.sum(), float(w @ x)],
        [float(w @ x), float(w @ (x * x)) + prior_prec],
    ])
    try:
        cov = np.linalg.inv(fisher)
        se = float(np.sqrt(max(cov[1, 1], 0.0)))
    except np.linalg.LinAlgError:
        se = np.nan

    x_var = x[cm > 0]
    x_ref = x[(ct - cm) > 0]
    separated = bool(
        x_var.size and x_ref.size
        and (x_var.min() > x_ref.max() or x_ref.min() > x_var.max())
    )

    growth = float(beta[1])
    return LogisticGrowthFit(
        growth_per_day=growth,
        se_per_day=se,
        log_R_ratio=float(generation_interval.log_R_ratio(growth)),
        intercept=float(beta[0]),
        n_bins=n_inf,
        n_variant_bins=variant_bins,
        carriers=carriers,
        span_days=span_days,
        converged=converged,
        separated=separated,
        status="ok" if converged else "not_converged",
    )


# ---------------------------------------------------------------------------
# Stratified fit: one shared slope, per-stratum intercepts
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class StratifiedGrowthFit:
    growth_per_day: float
    se_per_day: float
    log_R_ratio: float
    n_strata: int
    n_cells: int
    heterogeneity_chi2: float
    heterogeneity_df: int
    heterogeneity_p: float
    converged: bool
    status: str

    def as_dict(self) -> dict:
        return asdict(self)


def fit_logistic_growth_stratified(
    c_child, c_parent, bin_days: float = 14.0, *,
    generation_interval: GenerationInterval = None,
    slope_prior_sd: float = 0.25,
    min_bins_per_stratum: int = 2,
    max_iter: int = 200,
    tol: float = 1e-10,
) -> StratifiedGrowthFit:
    """Fit ``logit p[t, r] = a_r + s*t`` over a (time x stratum) count pair.

    One shared growth rate, one intercept per stratum. Controls for geography
    without collapsing strata to fragile point estimates (Finding S5).
    Includes a likelihood-ratio test of heterogeneity across strata.
    """
    gi = generation_interval or H3N2_GENERATION_INTERVAL
    cc = np.asarray(c_child, dtype=float)
    cp = np.asarray(c_parent, dtype=float)
    if cc.shape != cp.shape:
        raise ValueError("child and parent count arrays must be the same shape")
    if cc.ndim == 1:
        cc, cp = cc[:, None], cp[:, None]
    if cc.ndim != 2:
        raise ValueError("counts must be (time, stratum)")
    if np.any(cc < 0) or np.any(cp < 0):
        raise ValueError("counts must be non-negative")
    if not np.all(np.isfinite(cc)) or not np.all(np.isfinite(cp)):
        raise ValueError("counts must be finite (got NaN or inf)")

    n_t = cc.shape[0]
    tot = cc + cp

    def _viable(a, b):
        return (np.count_nonzero((a + b) > 0) >= min_bins_per_stratum
                and a.sum() > 0 and b.sum() > 0)

    own = [r for r in range(cc.shape[1]) if _viable(cc[:, r], cp[:, r])]
    thin = [r for r in range(cc.shape[1]) if r not in own and tot[:, r].sum() > 0]
    cols = [(cc[:, r], cp[:, r]) for r in own]
    if thin:
        m_c = cc[:, thin].sum(axis=1)
        m_p = cp[:, thin].sum(axis=1)
        if _viable(m_c, m_p):
            cols.append((m_c, m_p))
    if not cols:
        m_c, m_p = cc.sum(axis=1), cp.sum(axis=1)
        if _viable(m_c, m_p):
            cols = [(m_c, m_p)]
    keep = list(range(len(cols)))
    if not cols:
        return StratifiedGrowthFit(
            n_strata=0, n_cells=0, status="no_usable_stratum",
            growth_per_day=np.nan, se_per_day=np.nan, log_R_ratio=np.nan,
            heterogeneity_chi2=np.nan, heterogeneity_df=0,
            heterogeneity_p=np.nan, converged=False)

    x_all = np.arange(n_t, dtype=float) * float(bin_days)
    rows = []
    for slot, (a, b) in enumerate(cols):
        m = (a + b) > 0
        rows.append((x_all[m], a[m], (a + b)[m], np.full(int(m.sum()), slot)))
    x = np.concatenate([a for a, _, _, _ in rows])
    y = np.concatenate([b for _, b, _, _ in rows])
    n = np.concatenate([c for _, _, c, _ in rows])
    g = np.concatenate([d for _, _, _, d in rows]).astype(int)
    x = x - x.mean()
    n_s = len(keep)

    prior_prec = 0.0 if not np.isfinite(slope_prior_sd) or slope_prior_sd <= 0         else 1.0 / float(slope_prior_sd) ** 2

    def _loglik(beta):
        eta = beta[g] + beta[-1] * x
        eta = np.clip(eta, -30.0, 30.0)
        return float(np.sum(y * eta - n * np.log1p(np.exp(eta))))

    beta = np.zeros(n_s + 1)
    converged = False
    for _ in range(max_iter):
        eta = np.clip(beta[g] + beta[-1] * x, -30.0, 30.0)
        p = 1.0 / (1.0 + np.exp(-eta))
        resid = y - n * p
        w = n * p * (1.0 - p)

        grad = np.zeros(n_s + 1)
        np.add.at(grad, g, resid)
        grad[-1] = float(resid @ x) - prior_prec * beta[-1]

        d_aa = np.zeros(n_s)
        np.add.at(d_aa, g, w)
        d_as = np.zeros(n_s)
        np.add.at(d_as, g, w * x)
        d_ss = float(w @ (x * x)) + prior_prec
        H = np.zeros((n_s + 1, n_s + 1))
        H[np.arange(n_s), np.arange(n_s)] = d_aa
        H[np.arange(n_s), n_s] = d_as
        H[n_s, np.arange(n_s)] = d_as
        H[n_s, n_s] = d_ss
        try:
            step = np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.solve(H + 1e-9 * np.eye(n_s + 1), grad)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            converged = True
            break

    eta = np.clip(beta[g] + beta[-1] * x, -30.0, 30.0)
    p = 1.0 / (1.0 + np.exp(-eta))
    w = n * p * (1.0 - p)
    d_aa = np.zeros(n_s)
    np.add.at(d_aa, g, w)
    d_as = np.zeros(n_s)
    np.add.at(d_as, g, w * x)
    H = np.zeros((n_s + 1, n_s + 1))
    H[np.arange(n_s), np.arange(n_s)] = d_aa
    H[np.arange(n_s), n_s] = d_as
    H[n_s, np.arange(n_s)] = d_as
    H[n_s, n_s] = float(w @ (x * x)) + prior_prec
    try:
        se = float(np.sqrt(max(np.linalg.inv(H)[n_s, n_s], 0.0)))
    except np.linalg.LinAlgError:
        se = np.nan

    chi2 = np.nan
    df = max(n_s - 1, 0)
    p_val = np.nan
    if n_s > 1:
        ll_shared = _loglik(beta)
        ll_sep = 0.0
        for slot in range(n_s):
            m = g == slot
            b = np.zeros(2)
            for _ in range(max_iter):
                e = np.clip(b[0] + b[1] * x[m], -30.0, 30.0)
                pp = 1.0 / (1.0 + np.exp(-e))
                r_ = y[m] - n[m] * pp
                ww = n[m] * pp * (1.0 - pp)
                gr = np.array([r_.sum(), float(r_ @ x[m]) - prior_prec * b[1]])
                hh = np.array([[ww.sum(), float(ww @ x[m])],
                               [float(ww @ x[m]), float(ww @ (x[m] ** 2)) + prior_prec]])
                try:
                    st = np.linalg.solve(hh, gr)
                except np.linalg.LinAlgError:
                    break
                b = b + st
                if np.max(np.abs(st)) < tol:
                    break
            e = np.clip(b[0] + b[1] * x[m], -30.0, 30.0)
            ll_sep += float(np.sum(y[m] * e - n[m] * np.log1p(np.exp(e))))
        chi2 = float(max(2.0 * (ll_sep - ll_shared), 0.0))
        p_val = float(chi2_sf(chi2, df)) if df > 0 else np.nan

    growth = float(beta[-1])
    return StratifiedGrowthFit(
        growth_per_day=growth,
        se_per_day=se,
        log_R_ratio=float(gi.log_R_ratio(growth)),
        n_strata=n_s,
        n_cells=int(x.size),
        heterogeneity_chi2=chi2,
        heterogeneity_df=df,
        heterogeneity_p=p_val,
        converged=converged,
        status="ok" if converged else "not_converged",
    )


# ---------------------------------------------------------------------------
# Multinomial logistic growth across lineages (pure numpy L-BFGS solver)
# ---------------------------------------------------------------------------
def _lbfgs_minimize(func_grad, x0, m=10, max_iter=500, tol=1e-5):
    """Pure NumPy L-BFGS optimizer with Armijo backtracking line search."""
    x = np.array(x0, dtype=float)
    f, g = func_grad(x)
    s_list, y_list, rho_list = [], [], []
    for it in range(max_iter):
        if np.max(np.abs(g)) < tol:
            break
        q = g.copy()
        alpha_list = []
        for s, y, rho in reversed(list(zip(s_list, y_list, rho_list))):
            alpha = rho * float(np.dot(s, q))
            alpha_list.append(alpha)
            q -= alpha * y
        gamma = (float(np.dot(s_list[-1], y_list[-1])) / float(np.dot(y_list[-1], y_list[-1]))) if s_list else 1.0
        r = gamma * q
        for (s, y, rho), alpha in zip(zip(s_list, y_list, rho_list), reversed(alpha_list)):
            beta = rho * float(np.dot(y, r))
            r += s * (alpha - beta)
        p = -r

        step_size = 1.0
        c1 = 1e-4
        step_factor = 0.5
        dg_init = float(np.dot(g, p))
        if dg_init >= 0:
            p = -g
            dg_init = float(np.dot(g, p))

        x_new, f_new, g_new = x, f, g
        for _ in range(35):
            x_candidate = x + step_size * p
            f_cand, g_cand = func_grad(x_candidate)
            if f_cand <= f + c1 * step_size * dg_init:
                x_new, f_new, g_new = x_candidate, f_cand, g_cand
                break
            step_size *= step_factor
        else:
            x_new, f_new, g_new = x + step_size * p, func_grad(x + step_size * p)[0], func_grad(x + step_size * p)[1]

        s_k = x_new - x
        y_k = g_new - g
        sy = float(np.dot(s_k, y_k))
        if sy > 1e-10:
            if len(s_list) >= m:
                s_list.pop(0)
                y_list.pop(0)
                rho_list.pop(0)
            s_list.append(s_k)
            y_list.append(y_k)
            rho_list.append(1.0 / sy)
        x, f, g = x_new, f_new, g_new
    return x


def fit_multinomial_logistic_growth(
    counts: np.ndarray,
    bin_days: float = 14.0,
    *,
    generation_interval: GenerationInterval = H3N2_GENERATION_INTERVAL,
    reference_index: int = 0,
    slope_prior_sd: float = 0.25,
    min_cell_total: int = 3,
    max_iter: int = 500,
):
    """Fit per-lineage growth advantages from a (time, region, lineage) count tensor.

    Model: ``log p[t, r, k] ~ a[r, k] + s[k] * t``, with ``a[.,ref] = s[ref] = 0``.
    Returns ``(growth_per_day, se_per_day, log_R_ratio, n_cells)``.
    """
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 3:
        raise ValueError("counts must be a (time, region, lineage) tensor")
    n_t, n_r, n_k = counts.shape
    if not 0 <= reference_index < n_k:
        raise ValueError("reference_index out of range")

    t_centre = (np.arange(n_t, dtype=float) - (n_t - 1) / 2.0) * float(bin_days)
    cell_tot = counts.sum(axis=2)
    keep_t, keep_r = np.nonzero(cell_tot >= min_cell_total)
    if keep_t.size == 0:
        raise ValueError("no spatiotemporal cell reaches min_cell_total")
    C = counts[keep_t, keep_r, :]
    x = t_centre[keep_t]
    reg = keep_r

    free_k = [k for k in range(n_k) if k != reference_index]
    n_free = len(free_k)
    regions_present = sorted(set(reg.tolist()))
    r_slot = {r: j for j, r in enumerate(regions_present)}
    n_int = len(regions_present) * n_free

    prior_prec = 0.0 if not np.isfinite(slope_prior_sd) or slope_prior_sd <= 0         else 1.0 / float(slope_prior_sd) ** 2

    int_index = np.full((len(regions_present), n_free), -1, dtype=int)
    for r_j in range(len(regions_present)):
        for k_j in range(n_free):
            int_index[r_j, k_j] = r_j * n_free + k_j
    cell_r_slot = np.array([r_slot[r] for r in reg])

    def unpack(theta):
        a = np.zeros((len(regions_present), n_k))
        s = np.zeros(n_k)
        flat_int = theta[:n_int].reshape(len(regions_present), n_free)
        for k_j, k in enumerate(free_k):
            a[:, k] = flat_int[:, k_j]
            s[k] = theta[n_int + k_j]
        return a, s

    def nll_grad(theta):
        a, s = unpack(theta)
        eta = a[cell_r_slot, :] + s[None, :] * x[:, None]
        eta -= eta.max(axis=1, keepdims=True)
        expo = np.exp(eta)
        denom = expo.sum(axis=1, keepdims=True)
        logp = eta - np.log(denom)
        p = expo / denom
        n_cell = C.sum(axis=1, keepdims=True)
        nll = -float((C * logp).sum())
        resid = n_cell * p - C

        g = np.zeros_like(theta)
        for k_j, k in enumerate(free_k):
            np.add.at(g, int_index[cell_r_slot, k_j], resid[:, k])
            g[n_int + k_j] = float(resid[:, k] @ x)
        slopes = theta[n_int:]
        nll += 0.5 * prior_prec * float(slopes @ slopes)
        g[n_int:] += prior_prec * slopes
        return nll, g

    theta0 = np.zeros(n_int + n_free)
    try:
        from scipy import optimize
        res = optimize.minimize(nll_grad, theta0, jac=True, method="L-BFGS-B",
                                options={"maxiter": max_iter})
        theta_opt = res.x
    except ImportError:
        theta_opt = _lbfgs_minimize(nll_grad, theta0, max_iter=max_iter)

    _, s_hat = unpack(theta_opt)

    eps = 1e-6
    n_par = theta_opt.size
    H = np.zeros((n_par, n_par))
    for i in range(n_par):
        xp, xm = theta_opt.copy(), theta_opt.copy()
        xp[i] += eps
        xm[i] -= eps
        H[:, i] = (nll_grad(xp)[1] - nll_grad(xm)[1]) / (2 * eps)
    H = 0.5 * (H + H.T)
    cov = np.linalg.pinv(H)
    se_hat = np.zeros(n_k)
    for k_j, k in enumerate(free_k):
        se_hat[k] = float(np.sqrt(max(cov[n_int + k_j, n_int + k_j], 0.0)))

    log_R = np.zeros(n_k)
    for k in range(n_k):
        log_R[k] = float(generation_interval.log_R_ratio(s_hat[k]))
    return s_hat, se_hat, log_R, int(C.shape[0])


# ---------------------------------------------------------------------------
# Design identifiability
# ---------------------------------------------------------------------------
def design_identifiability(X: np.ndarray, labels: Sequence) -> dict:
    """Report what a lineage-by-mutation design can and cannot identify (Finding C3)."""
    X = np.asarray(X, dtype=float)
    n_obs, n_par = X.shape
    rank = int(np.linalg.matrix_rank(X))
    groups: dict[tuple, list] = {}
    for j in range(n_par):
        groups.setdefault(tuple(X[:, j]), []).append(labels[j])
    confounded = sorted((v for v in groups.values() if len(v) > 1), key=len, reverse=True)
    return {
        "n_observations": n_obs,
        "n_parameters": n_par,
        "rank": rank,
        "null_space_dim": n_par - rank,
        "n_distinct_patterns": len(groups),
        "identifiable": n_par <= rank,
        "confounded_groups": confounded,
        "largest_confounded_group": max((len(v) for v in groups.values()), default=0),
    }


# ---------------------------------------------------------------------------
# Correlations with honest uncertainty
# ---------------------------------------------------------------------------
def weighted_pearson(x, y, w=None) -> float:
    """Pearson correlation, optionally weighted (e.g. by inverse variance)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if w is not None:
        w = np.asarray(w, dtype=float)
        ok &= np.isfinite(w) & (w > 0)
    if ok.sum() < 3:
        return np.nan
    x, y = x[ok], y[ok]
    ww = np.ones_like(x) if w is None else w[ok]
    ww = ww / ww.sum()
    mx, my = float(ww @ x), float(ww @ y)
    cxy = float(ww @ ((x - mx) * (y - my)))
    cxx = float(ww @ (x - mx) ** 2)
    cyy = float(ww @ (y - my) ** 2)
    mss_x = float(ww @ (x * x))
    mss_y = float(ww @ (y * y))
    if cxx <= 1e-14 * max(mss_x, 1e-300) or cyy <= 1e-14 * max(mss_y, 1e-300):
        return np.nan
    return cxy / np.sqrt(cxx * cyy)


def correlation_with_ci(
    x, y, *, weights=None, method: str = "pearson", n_boot: int = 2000,
    groups=None, seed: int = 20260917, ci: float = 95.0,
) -> dict:
    """Correlation plus a bootstrap confidence interval."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    w = None
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        ok &= np.isfinite(w) & (w > 0)
    g = None if groups is None else np.asarray(groups)[ok]
    x, y = x[ok], y[ok]
    if w is not None:
        w = w[ok]
    n = x.size
    out = {"n": int(n), "n_dropped": int(np.size(ok) - n), "method": method}
    if n < 3:
        return {**out, "estimate": np.nan, "ci_low": np.nan, "ci_high": np.nan}

    def _stat(xi, yi, wi):
        if method == "spearman":
            if np.unique(xi).size < 2 or np.unique(yi).size < 2:
                return np.nan
            rx, ry = _rankdata(xi), _rankdata(yi)
            return weighted_pearson(rx, ry, wi)
        return weighted_pearson(xi, yi, wi)

    est = _stat(x, y, w)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    if g is None:
        for b in range(n_boot):
            i = rng.integers(0, n, n)
            boots[b] = _stat(x[i], y[i], None if w is None else w[i])
    else:
        uniq = np.unique(g)
        index_by_group = {u: np.flatnonzero(g == u) for u in uniq}
        for b in range(n_boot):
            pick = rng.integers(0, uniq.size, uniq.size)
            i = np.concatenate([index_by_group[uniq[j]] for j in pick])
            boots[b] = _stat(x[i], y[i], None if w is None else w[i])
    boots = boots[np.isfinite(boots)]
    lo, hi = (np.nan, np.nan) if boots.size < 20 else np.percentile(
        boots, [(100 - ci) / 2, 100 - (100 - ci) / 2])
    return {**out, "estimate": est, "ci_low": float(lo), "ci_high": float(hi),
            "n_boot": int(boots.size)}
