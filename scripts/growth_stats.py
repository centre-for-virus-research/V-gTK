#!/usr/bin/env python3
"""Statistical primitives for lineage growth-rate estimation.

Everything here is deliberately dependency-light: the `vgtk` environment ships
numpy, pandas, biopython and matplotlib and **no scipy or statsmodels**, so the
distributions, the GLM solver and the hypothesis tests are implemented directly
rather than imported. Adding scipy to `environment.yml` for a bolt-on analysis
would change the environment every pipeline run has to solve, which is a much
bigger change than the ~300 lines below.

The module is split into three layers:

1. **Distributions** - normal, chi-squared and Student-t tails, plus the normal
   quantile. These are what turn a test statistic into a p-value or a
   confidence interval.
2. **A generalised linear model** fitted by iteratively reweighted least
   squares, with binomial (logit) and Poisson (log) families, optional
   quasi-likelihood dispersion, and optional Firth penalisation for the
   separated fits that small lineages produce constantly.
3. **Growth-rate estimators** built on top of those: the logistic /
   multinomial-logistic frequency models, the Poisson count model, the
   Frequency Increment Test, and the non-parametric trend tests.

Nothing here knows about trees, databases or viruses.
"""

import math

import numpy as np

#: Probabilities are clipped this far away from 0 and 1 before any log or
#: division. Without it a perfectly separated fit (every late sequence a
#: carrier) divides by zero in the IRLS weights instead of simply running the
#: coefficient off to infinity, which the caller can at least detect.
_EPS = 1e-10

#: IRLS is declared converged when the largest absolute coefficient change
#: falls below this. Deviance-based stopping was tried first and stalls on
#: separated fits, where the deviance is flat while the coefficient is still
#: doubling every iteration.
_TOL = 1e-8

_MAX_ITER = 100


# ---------------------------------------------------------------------------
# Distributions
# ---------------------------------------------------------------------------

def norm_cdf(z):
	"""Standard normal CDF, via the error function in the standard library."""
	return 0.5 * math.erfc(-float(z) / math.sqrt(2.0))


def norm_sf(z):
	"""Upper tail of the standard normal."""
	return 0.5 * math.erfc(float(z) / math.sqrt(2.0))


def two_sided_normal_p(z):
	"""Two-sided p-value for a z statistic, safe for NaN and infinite z."""
	z = float(z)
	if not math.isfinite(z):
		return 0.0 if math.isinf(z) else float('nan')
	return math.erfc(abs(z) / math.sqrt(2.0))


#: Coefficients of Acklam's rational approximation to the normal quantile. It is
#: accurate to ~1.15e-9 in relative terms, and the Halley refinement below takes
#: it to full double precision - which matters because the 95% CI multiplier is
#: computed from it once and then applied to every estimate in the run.
_ACKLAM_A = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
			 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
_ACKLAM_B = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
			 6.680131188771972e+01, -1.328068155288572e+01)
_ACKLAM_C = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
			 -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
_ACKLAM_D = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
			 3.754408661907416e+00)
_ACKLAM_PLOW = 0.02425


def norm_ppf(p):
	"""Standard normal quantile."""
	p = float(p)
	if not 0.0 < p < 1.0:
		if p == 0.0:
			return float('-inf')
		if p == 1.0:
			return float('inf')
		return float('nan')
	if p < _ACKLAM_PLOW:
		q = math.sqrt(-2 * math.log(p))
		x = (((((_ACKLAM_C[0] * q + _ACKLAM_C[1]) * q + _ACKLAM_C[2]) * q + _ACKLAM_C[3]) * q
			  + _ACKLAM_C[4]) * q + _ACKLAM_C[5]) / ((((_ACKLAM_D[0] * q + _ACKLAM_D[1]) * q
													  + _ACKLAM_D[2]) * q + _ACKLAM_D[3]) * q + 1)
	elif p <= 1 - _ACKLAM_PLOW:
		q = p - 0.5
		r = q * q
		x = (((((_ACKLAM_A[0] * r + _ACKLAM_A[1]) * r + _ACKLAM_A[2]) * r + _ACKLAM_A[3]) * r
			  + _ACKLAM_A[4]) * r + _ACKLAM_A[5]) * q / (((((_ACKLAM_B[0] * r + _ACKLAM_B[1]) * r
															+ _ACKLAM_B[2]) * r + _ACKLAM_B[3]) * r
														  + _ACKLAM_B[4]) * r + 1)
	else:
		q = math.sqrt(-2 * math.log(1 - p))
		x = -(((((_ACKLAM_C[0] * q + _ACKLAM_C[1]) * q + _ACKLAM_C[2]) * q + _ACKLAM_C[3]) * q
			   + _ACKLAM_C[4]) * q + _ACKLAM_C[5]) / ((((_ACKLAM_D[0] * q + _ACKLAM_D[1]) * q
													   + _ACKLAM_D[2]) * q + _ACKLAM_D[3]) * q + 1)
	# One Halley step against the exact CDF.
	err = norm_cdf(x) - p
	density = math.exp(-x * x / 2) / math.sqrt(2 * math.pi)
	if density > 0:
		u = err / density
		x = x - u / (1 + x * u / 2)
	return x


def _gamma_series(a, x):
	"""Lower regularised incomplete gamma P(a, x) by series expansion."""
	if x <= 0:
		return 0.0
	ap = a
	total = 1.0 / a
	delta = total
	for _ in range(1000):
		ap += 1
		delta *= x / ap
		total += delta
		if abs(delta) < abs(total) * 1e-15:
			break
	return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gamma_continued_fraction(a, x):
	"""Upper regularised incomplete gamma Q(a, x) by continued fraction."""
	tiny = 1e-300
	b = x + 1.0 - a
	c = 1.0 / tiny
	d = 1.0 / b if b != 0 else 1.0 / tiny
	h = d
	for i in range(1, 1000):
		an = -i * (i - a)
		b += 2.0
		d = an * d + b
		if abs(d) < tiny:
			d = tiny
		c = b + an / c
		if abs(c) < tiny:
			c = tiny
		d = 1.0 / d
		delta = d * c
		h *= delta
		if abs(delta - 1.0) < 1e-15:
			break
	return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def gammainc_upper(a, x):
	"""Regularised upper incomplete gamma Q(a, x) = 1 - P(a, x)."""
	a = float(a)
	x = float(x)
	if x < 0 or a <= 0:
		return float('nan')
	if x == 0:
		return 1.0
	if x < a + 1.0:
		return 1.0 - _gamma_series(a, x)
	return _gamma_continued_fraction(a, x)


def chi2_sf(x, df):
	"""Upper tail of the chi-squared distribution."""
	x = float(x)
	if not math.isfinite(x):
		return 0.0 if x > 0 else float('nan')
	if x <= 0:
		return 1.0
	return gammainc_upper(df / 2.0, x / 2.0)


def _betacf(a, b, x):
	"""Continued fraction for the incomplete beta function (Lentz's method)."""
	tiny = 1e-300
	qab, qap, qam = a + b, a + 1.0, a - 1.0
	c = 1.0
	d = 1.0 - qab * x / qap
	if abs(d) < tiny:
		d = tiny
	d = 1.0 / d
	h = d
	for m in range(1, 300):
		m2 = 2 * m
		aa = m * (b - m) * x / ((qam + m2) * (a + m2))
		d = 1.0 + aa * d
		if abs(d) < tiny:
			d = tiny
		c = 1.0 + aa / c
		if abs(c) < tiny:
			c = tiny
		d = 1.0 / d
		h *= d * c
		aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
		d = 1.0 + aa * d
		if abs(d) < tiny:
			d = tiny
		c = 1.0 + aa / c
		if abs(c) < tiny:
			c = tiny
		d = 1.0 / d
		delta = d * c
		h *= delta
		if abs(delta - 1.0) < 1e-15:
			break
	return h


def betainc(a, b, x):
	"""Regularised incomplete beta I_x(a, b)."""
	a, b, x = float(a), float(b), float(x)
	if x <= 0:
		return 0.0
	if x >= 1:
		return 1.0
	front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
					 + a * math.log(x) + b * math.log1p(-x))
	if x < (a + 1.0) / (a + b + 2.0):
		return front * _betacf(a, b, x) / a
	return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_sf(t, df):
	"""Upper tail of Student's t with `df` degrees of freedom."""
	t = float(t)
	df = float(df)
	if df <= 0 or not math.isfinite(t):
		return float('nan')
	tail = 0.5 * betainc(df / 2.0, 0.5, df / (df + t * t))
	return tail if t >= 0 else 1.0 - tail


def two_sided_t_p(t, df):
	"""Two-sided p-value for a t statistic."""
	t = float(t)
	if not math.isfinite(t) or df <= 0:
		return float('nan')
	return betainc(df / 2.0, 0.5, df / (df + t * t))


# ---------------------------------------------------------------------------
# Hypothesis tests that do not need a model
# ---------------------------------------------------------------------------

def fisher_exact_2x2(a, b, c, d):
	"""Two-sided Fisher exact test on [[a, b], [c, d]].

	Returns ``(odds_ratio, p_value)``. The odds ratio uses the Haldane-Anscombe
	0.5 correction when a cell is empty, so an infinite ratio is reported as a
	large finite number rather than an inf that every downstream format string
	then has to special-case.
	"""
	a, b, c, d = int(a), int(b), int(c), int(d)
	n = a + b + c + d
	if n == 0:
		return float('nan'), float('nan')
	row1, row2 = a + b, c + d
	col1, col2 = a + c, b + d
	if min(row1, row2, col1, col2) == 0:
		return float('nan'), 1.0

	def log_prob(x):
		return (math.lgamma(row1 + 1) + math.lgamma(row2 + 1)
				+ math.lgamma(col1 + 1) + math.lgamma(col2 + 1)
				- math.lgamma(n + 1) - math.lgamma(x + 1)
				- math.lgamma(row1 - x + 1) - math.lgamma(col1 - x + 1)
				- math.lgamma(row2 - col1 + x + 1))

	low = max(0, col1 - row2)
	high = min(row1, col1)
	observed = log_prob(a)
	# The 1e-7 relative slack is the conventional guard against a table whose
	# probability is equal to the observed one only up to floating-point noise
	# being excluded from the tail.
	total = 0.0
	for x in range(low, high + 1):
		p = math.exp(log_prob(x))
		if log_prob(x) <= observed + 1e-7:
			total += p
	p_value = min(1.0, total)
	if b * c == 0 or a * d == 0:
		odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
	else:
		odds = (a * d) / (b * c)
	return odds, p_value


def mann_kendall(times, values):
	"""Non-parametric monotonic trend test with Sen's slope.

	Returns a dict with the S statistic, Kendall's tau, the tie-corrected
	normal approximation and its p-value, and Sen's slope (the median pairwise
	slope, in units of `values` per unit of `times`).

	This is the estimator that survives a trajectory the logistic model cannot
	describe - a frequency that rises and then plateaus, or one driven by a
	handful of very large sampling batches.
	"""
	times = np.asarray(times, dtype=float)
	values = np.asarray(values, dtype=float)
	keep = np.isfinite(times) & np.isfinite(values)
	times, values = times[keep], values[keep]
	order = np.argsort(times, kind='mergesort')
	times, values = times[order], values[order]
	n = len(values)
	out = {'n': n, 'S': float('nan'), 'tau': float('nan'), 'z': float('nan'),
		   'p_value': float('nan'), 'sen_slope': float('nan')}
	if n < 3:
		return out

	diff = np.sign(values[None, :] - values[:, None])
	s = float(np.sum(np.triu(diff, 1)))
	_, counts = np.unique(values, return_counts=True)
	tie_term = float(np.sum(counts * (counts - 1) * (2 * counts + 5)))
	variance = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0
	if variance <= 0:
		return out
	if s > 0:
		z = (s - 1) / math.sqrt(variance)
	elif s < 0:
		z = (s + 1) / math.sqrt(variance)
	else:
		z = 0.0

	slopes = []
	for i in range(n - 1):
		dt = times[i + 1:] - times[i]
		dv = values[i + 1:] - values[i]
		good = dt != 0
		if np.any(good):
			slopes.append(dv[good] / dt[good])
	sen = float(np.median(np.concatenate(slopes))) if slopes else float('nan')

	out.update({'S': s, 'tau': 2.0 * s / (n * (n - 1)), 'z': z,
				'p_value': two_sided_normal_p(z), 'sen_slope': sen})
	return out


def benjamini_hochberg(p_values):
	"""Benjamini-Hochberg q-values, NaN-preserving.

	Every run here tests hundreds of alleles or thousands of tree nodes, so an
	uncorrected p-value column would be actively misleading. NaNs (a fit that
	did not converge) are carried through as NaN rather than being treated as
	p = 1, which would inflate the number of tests and deflate everyone else's
	q-value.
	"""
	p = np.asarray(p_values, dtype=float)
	q = np.full(p.shape, np.nan)
	finite = np.isfinite(p)
	if not np.any(finite):
		return q
	idx = np.flatnonzero(finite)
	order = idx[np.argsort(p[idx], kind='mergesort')]
	m = len(order)
	ranks = np.arange(1, m + 1)
	adjusted = p[order] * m / ranks
	adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
	q[order] = np.minimum(adjusted, 1.0)
	return q


def wilson_interval(successes, trials, confidence=0.95):
	"""Wilson score interval for a binomial proportion.

	Preferred over the Wald interval because most per-time-bin counts here are
	small or at a boundary, where Wald intervals leave the unit interval.
	"""
	if trials <= 0:
		return float('nan'), float('nan')
	z = norm_ppf(0.5 + confidence / 2.0)
	phat = successes / trials
	denom = 1 + z * z / trials
	centre = (phat + z * z / (2 * trials)) / denom
	half = z * math.sqrt(phat * (1 - phat) / trials + z * z / (4 * trials * trials)) / denom
	return max(0.0, centre - half), min(1.0, centre + half)


# ---------------------------------------------------------------------------
# Generalised linear models
# ---------------------------------------------------------------------------

class GlmFit(object):
	"""Result of one IRLS fit. Plain attributes, so it survives `vars()`."""

	def __init__(self, **kwargs):
		self.params = kwargs.get('params')
		self.cov = kwargs.get('cov')
		self.se = kwargs.get('se')
		self.loglik = kwargs.get('loglik', float('nan'))
		self.penalised_loglik = kwargs.get('penalised_loglik', float('nan'))
		self.deviance = kwargs.get('deviance', float('nan'))
		self.dispersion = kwargs.get('dispersion', 1.0)
		self.n = kwargs.get('n', 0)
		self.df_resid = kwargs.get('df_resid', 0)
		self.converged = kwargs.get('converged', False)
		self.iterations = kwargs.get('iterations', 0)
		self.family = kwargs.get('family', '')
		self.separated = kwargs.get('separated', False)

	def z_values(self):
		with np.errstate(divide='ignore', invalid='ignore'):
			return np.where(self.se > 0, self.params / self.se, np.nan)

	def wald_p(self):
		return np.array([two_sided_normal_p(z) for z in self.z_values()])

	def conf_int(self, confidence=0.95):
		z = norm_ppf(0.5 + confidence / 2.0)
		return self.params - z * self.se, self.params + z * self.se


def _link_inverse(eta, family):
	if family == 'binomial':
		# Clip the linear predictor, not the probability: exp(800) overflows to
		# inf and then inf/inf gives NaN, which poisons the whole fit rather
		# than saturating it.
		return 1.0 / (1.0 + np.exp(-np.clip(eta, -500, 500)))
	return np.exp(np.clip(eta, -500, 500))


def _solve_spd(matrix, rhs):
	"""Solve a symmetric positive-definite system, ridging if it is singular.

	Returns ``(solution, inverse, ridged)``. Separation and collinear designs
	both make X'WX singular; refusing to fit at all would drop exactly the
	fast-growing lineages the analysis exists to find, so a tiny ridge is added
	and the caller is told it happened.
	"""
	scale = float(np.mean(np.abs(np.diag(matrix)))) or 1.0
	for ridge in (0.0, 1e-10, 1e-6, 1e-3):
		try:
			damped = matrix + ridge * scale * np.eye(matrix.shape[0])
			inverse = np.linalg.inv(damped)
			solution = inverse.dot(rhs)
			if np.all(np.isfinite(solution)):
				return solution, inverse, ridge > 0
		except np.linalg.LinAlgError:
			continue
	inverse = np.linalg.pinv(matrix)
	return inverse.dot(rhs), inverse, True


def fit_glm(X, y, family='binomial', trials=None, offset=None, firth=False,
			max_iter=_MAX_ITER, tol=_TOL):
	"""Fit a GLM by iteratively reweighted least squares.

	`family` is 'binomial' (logit link; `y` is a proportion and `trials` the
	number of observations behind it) or 'poisson' (log link; `y` is a count and
	`offset` the log of the exposure).

	`firth=True` applies Firth's penalised likelihood to a binomial fit. This is
	not a refinement - it is what makes the small-lineage results usable at all.
	A clade whose carriers are *all* in the last two years is perfectly
	separated: the maximum-likelihood slope is infinite, the standard error is
	infinite with it, and the Wald p-value comes back at 1.0 for what is the
	strongest signal in the dataset. Firth's penalty removes the first-order
	bias and always produces a finite estimate.
	"""
	X = np.asarray(X, dtype=float)
	y = np.asarray(y, dtype=float)
	n, p = X.shape
	trials = np.ones(n) if trials is None else np.asarray(trials, dtype=float)
	offset = np.zeros(n) if offset is None else np.asarray(offset, dtype=float)
	if family not in ('binomial', 'poisson'):
		raise ValueError("family must be 'binomial' or 'poisson'")
	if firth and family != 'binomial':
		raise ValueError('Firth penalisation is only implemented for the binomial family')

	beta = np.zeros(p)
	converged = False
	ridged = False
	iterations = 0
	inverse = np.eye(p)

	for iterations in range(1, max_iter + 1):
		eta = X.dot(beta) + offset
		mu = _link_inverse(eta, family)
		if family == 'binomial':
			mu = np.clip(mu, _EPS, 1 - _EPS)
			weights = trials * mu * (1 - mu)
			residual = trials * (y - mu)
		else:
			mu = np.maximum(mu, _EPS)
			weights = trials * mu
			residual = trials * (y - mu)

		weights = np.maximum(weights, _EPS)
		information = X.T.dot(X * weights[:, None])

		if firth:
			# Hat values of the weighted design; the penalty adds h/2 "half
			# observations" at each end, which is exactly what stops a
			# separated fit running away.
			_, info_inv, _ = _solve_spd(information, np.zeros(p))
			sqrt_w = np.sqrt(weights)
			scaled = X * sqrt_w[:, None]
			hat = np.einsum('ij,jk,ik->i', scaled, info_inv, scaled)
			residual = residual + hat * (0.5 - mu)

		step, inverse, step_ridged = _solve_spd(information, X.T.dot(residual))
		ridged = ridged or step_ridged
		# Step halving keeps a wild first step (common when every observation in
		# one time bin has the same outcome) from throwing the fit somewhere the
		# link function has already saturated.
		scale = 1.0
		for _ in range(20):
			candidate = beta + scale * step
			if np.all(np.isfinite(candidate)) and np.max(np.abs(candidate)) < 1e4:
				break
			scale /= 2.0
		else:
			candidate = beta
		delta = np.max(np.abs(candidate - beta))
		beta = candidate
		if delta < tol:
			converged = True
			break

	eta = X.dot(beta) + offset
	mu = _link_inverse(eta, family)
	if family == 'binomial':
		mu = np.clip(mu, _EPS, 1 - _EPS)
		weights = trials * mu * (1 - mu)
		successes = trials * y
		failures = trials * (1 - y)
		loglik = float(np.sum(successes * np.log(mu) + failures * np.log1p(-mu)))
		saturated = np.clip(y, _EPS, 1 - _EPS)
		deviance = 2.0 * float(np.sum(successes * np.log(saturated / mu)
									  + failures * np.log((1 - saturated) / (1 - mu))))
		pearson = (successes - trials * mu) / np.sqrt(np.maximum(trials * mu * (1 - mu), _EPS))
	else:
		mu = np.maximum(mu, _EPS)
		weights = trials * mu
		loglik = float(np.sum(trials * (y * np.log(mu) - mu)))
		with np.errstate(divide='ignore', invalid='ignore'):
			ratio = np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / mu), 0.0)
		deviance = 2.0 * float(np.sum(trials * (ratio - (y - mu))))
		pearson = (y - mu) * np.sqrt(trials) / np.sqrt(np.maximum(mu, _EPS))

	information = X.T.dot(X * np.maximum(weights, _EPS)[:, None])
	_, inverse, inv_ridged = _solve_spd(information, np.zeros(p))
	ridged = ridged or inv_ridged
	sign, logdet = np.linalg.slogdet(information)
	penalised = loglik + 0.5 * logdet if sign > 0 else loglik

	df_resid = max(n - p, 1)
	dispersion = float(np.sum(pearson ** 2)) / df_resid
	se = np.sqrt(np.maximum(np.diag(inverse), 0.0))

	return GlmFit(params=beta, cov=inverse, se=se, loglik=loglik,
				  penalised_loglik=penalised, deviance=deviance,
				  dispersion=dispersion, n=n, df_resid=df_resid,
				  converged=converged, iterations=iterations, family=family,
				  separated=ridged or not converged)


def likelihood_ratio_p(X, y, column, family='binomial', trials=None, offset=None, firth=False):
	"""p-value for dropping one design column, by (penalised) likelihood ratio.

	The Wald test is unreliable exactly where this analysis lives - large
	coefficients with large standard errors - so the reported significance for
	every growth rate comes from here, and Wald is kept alongside only for
	comparison.
	"""
	X = np.asarray(X, dtype=float)
	keep = [j for j in range(X.shape[1]) if j != column]
	full = fit_glm(X, y, family=family, trials=trials, offset=offset, firth=firth)
	if not keep:
		return float('nan'), full
	reduced = fit_glm(X[:, keep], y, family=family, trials=trials, offset=offset, firth=firth)
	if firth:
		statistic = 2.0 * (full.penalised_loglik - reduced.penalised_loglik)
	else:
		statistic = 2.0 * (full.loglik - reduced.loglik)
	if not math.isfinite(statistic) or statistic < 0:
		statistic = 0.0
	return chi2_sf(statistic, 1), full


# ---------------------------------------------------------------------------
# Growth-rate estimators
# ---------------------------------------------------------------------------

#: Below this many observations in either outcome class, the maximum-likelihood
#: logistic fit is biased away from zero badly enough to matter, so Firth's
#: penalty is used even when the fit did converge.
_FIRTH_MIN_CLASS = 25

#: Above this many observations, an individual-level logistic fit is collapsed
#: onto its distinct sampling times first. See logistic_growth for why that is
#: exact rather than an approximation.
_AGGREGATE_ABOVE = 2000

LOGISTIC_KEYS = (
	'n', 'n_carrier', 'n_background', 'time_min', 'time_max', 'time_span',
	'growth_advantage_per_year', 'std_error', 'ci_low', 'ci_high',
	'wald_p', 'lrt_p', 'intercept', 'fitted_freq_start', 'fitted_freq_end',
	'observed_freq', 'doubling_time_years', 'method', 'converged',
	'dispersion', 'overdispersed',
)


def _empty_logistic():
	result = {key: float('nan') for key in LOGISTIC_KEYS}
	result.update({'n': 0, 'n_carrier': 0, 'n_background': 0, 'method': 'none',
				   'converged': False, 'overdispersed': False})
	return result


def doubling_time(rate):
	"""Time for the quantity to double at exponential rate `rate`.

	A negative rate returns a negative number: the magnitude is the halving
	time. Returning NaN for declining lineages, as a doubling time strictly
	should, throws away the only number that says how fast they are going.
	"""
	rate = float(rate)
	if not math.isfinite(rate) or rate == 0:
		return float('nan')
	return math.log(2.0) / rate


def reproduction_number(rate, generation_time, generation_sd=0.0):
	"""Convert an exponential growth rate into a reproduction number.

	With a fixed generation interval this is the textbook R = exp(rT). With a
	gamma-distributed interval of mean T and standard deviation s it is the
	Wallinga-Lipsitch moment-generating result R = (1 + r s^2 / T)^(T^2 / s^2),
	which is materially smaller than exp(rT) for a realistically dispersed
	interval - assuming a fixed interval systematically overstates R.
	"""
	rate = float(rate)
	generation_time = float(generation_time)
	if not math.isfinite(rate) or generation_time <= 0:
		return float('nan')
	if generation_sd and generation_sd > 0:
		shape = (generation_time / generation_sd) ** 2
		base = 1.0 + rate * generation_sd ** 2 / generation_time
		if base <= 0:
			return float('nan')
		return base ** shape
	return math.exp(rate * generation_time)


def logistic_growth(times, carrier, trials=None, confidence=0.95, firth='auto',
					quasi=False):
	"""Logistic (frequency) growth rate of a lineage against its background.

	Fits ``logit(f(t)) = a + s.t`` where f is the frequency of the lineage among
	the sampled sequences. Under a two-type exponential model the slope `s` is
	exactly the difference in exponential growth rates between the lineage and
	everything it is competing with, which is why it - and not the raw growth in
	counts - is the standard measure of a variant's fitness advantage. It is
	also the reason this is robust to sampling effort: multiplying the number of
	sequences collected in a year by ten changes the counts and leaves the
	frequency alone.

	`trials` turns the fit from individual-level into binned counts (`carrier`
	is then the proportion in each bin). `quasi=True` estimates a dispersion
	parameter from the binned fit and inflates the standard errors by it, which
	is the honest thing to do with genomic surveillance data: sequences arrive
	in correlated batches - one outbreak, one submitting lab - so the binomial
	variance is an underestimate and the uncorrected confidence interval is too
	narrow.
	"""
	times = np.asarray(times, dtype=float)
	carrier = np.asarray(carrier, dtype=float)
	weights = np.ones_like(times) if trials is None else np.asarray(trials, dtype=float)
	keep = np.isfinite(times) & np.isfinite(carrier) & (weights > 0)
	times, carrier, weights = times[keep], carrier[keep], weights[keep]

	result = _empty_logistic()
	total = float(np.sum(weights))
	carriers = float(np.sum(carrier * weights))
	result.update({'n': int(round(total)), 'n_carrier': int(round(carriers)),
				   'n_background': int(round(total - carriers))})
	if len(times) < 3 or total <= 0:
		return result
	result['observed_freq'] = carriers / total
	result.update({'time_min': float(np.min(times)), 'time_max': float(np.max(times))})
	result['time_span'] = result['time_max'] - result['time_min']
	if result['time_span'] <= 0 or carriers <= 0 or carriers >= total:
		return result

	use_firth = firth is True
	if firth == 'auto':
		use_firth = min(carriers, total - carriers) < _FIRTH_MIN_CLASS

	# Sequences sharing a sampling date carry identical information, and the
	# binomial logit likelihood for grouped counts is *the same function* as the
	# one for the individual rows behind them. Firth's penalty is a function of
	# the Fisher information, which grouping likewise leaves unchanged, so the
	# penalised fit is invariant too - there are tests pinning both to 1e-8.
	#
	# Collapsing is therefore exact rather than an approximation, and it is what
	# makes a full influenza segment fit in minutes instead of hours. Excluding
	# the Firth path from it was much the worse half of the mistake: the rare
	# alleles are the ones that need the penalty *and* the large majority of any
	# protein, so they were each doing a hat-matrix fit over all 81,000 rows,
	# twice, for an answer identical to the one over a few thousand.
	fit_times, fit_response, fit_weights = times, carrier, weights
	if trials is None and len(times) > _AGGREGATE_ABOVE:
		unique_times, inverse = np.unique(times, return_inverse=True)
		counts = np.bincount(inverse).astype(float)
		fit_times = unique_times
		fit_response = np.bincount(inverse, weights=carrier) / counts
		fit_weights = counts

	centre = float(np.mean(fit_times))
	design = np.column_stack([np.ones_like(fit_times), fit_times - centre])
	fit = fit_glm(design, fit_response, 'binomial', trials=fit_weights, firth=use_firth)
	if not use_firth and firth == 'auto' and fit.separated:
		use_firth = True
		fit = fit_glm(design, fit_response, 'binomial', trials=fit_weights, firth=True)

	slope = float(fit.params[1])
	se = float(fit.se[1])
	dispersion = float(fit.dispersion)
	overdispersed = bool(quasi and trials is not None and dispersion > 1.0)
	if overdispersed:
		se *= math.sqrt(dispersion)

	if overdispersed:
		# With a quasi-likelihood dispersion the reference distribution is t,
		# not normal, and the likelihood ratio no longer has its chi-squared
		# calibration - so the interval and the p-value both come from the
		# scaled Wald statistic here.
		critical = -norm_ppf((1 - confidence) / 2.0) if fit.df_resid > 30 else 2.0
		p_value = two_sided_t_p(slope / se if se > 0 else float('nan'), fit.df_resid)
		lrt_p = float('nan')
	else:
		critical = -norm_ppf((1 - confidence) / 2.0)
		p_value = two_sided_normal_p(slope / se) if se > 0 else float('nan')
		lrt_p, _ = likelihood_ratio_p(design, fit_response, 1, 'binomial',
									  trials=fit_weights, firth=use_firth)

	intercept = float(fit.params[0]) - slope * centre
	result.update({
		'growth_advantage_per_year': slope,
		'std_error': se,
		'ci_low': slope - critical * se,
		'ci_high': slope + critical * se,
		'wald_p': p_value,
		'lrt_p': lrt_p,
		'intercept': intercept,
		'fitted_freq_start': float(_link_inverse(np.array([intercept + slope * result['time_min']]), 'binomial')[0]),
		'fitted_freq_end': float(_link_inverse(np.array([intercept + slope * result['time_max']]), 'binomial')[0]),
		'doubling_time_years': doubling_time(slope),
		'method': 'firth' if use_firth else 'ml',
		'converged': bool(fit.converged),
		'dispersion': dispersion,
		'overdispersed': overdispersed,
	})
	return result


def poisson_growth(times, counts, exposure=None, confidence=0.95, quasi=True):
	"""Exponential growth rate of a lineage's absolute count over time.

	``log(E[count]) = a + r.t`` (+ log exposure). Without an exposure this is the
	growth of the lineage *as sampled*, which mixes the virus's growth with
	sequencing effort; with `exposure` set to the total number of sequences per
	bin it becomes a relative rate again and should agree with the logistic
	slope while frequencies are low.

	Both are reported because they disagree informatively: a lineage growing in
	count while flat in frequency is riding an expanding epidemic, not
	outcompeting anything.
	"""
	times = np.asarray(times, dtype=float)
	counts = np.asarray(counts, dtype=float)
	keep = np.isfinite(times) & np.isfinite(counts)
	if exposure is not None:
		exposure = np.asarray(exposure, dtype=float)
		keep = keep & np.isfinite(exposure) & (exposure > 0)
	times, counts = times[keep], counts[keep]
	offset = np.log(exposure[keep]) if exposure is not None else None

	out = {'n_bins': len(times), 'growth_rate_per_year': float('nan'),
		   'std_error': float('nan'), 'ci_low': float('nan'), 'ci_high': float('nan'),
		   'p_value': float('nan'), 'doubling_time_years': float('nan'),
		   'dispersion': float('nan'), 'converged': False,
		   'relative_to_sampling': exposure is not None}
	if len(times) < 3 or np.sum(counts) <= 0 or np.ptp(times) <= 0:
		return out

	centre = float(np.mean(times))
	design = np.column_stack([np.ones_like(times), times - centre])
	fit = fit_glm(design, counts, 'poisson', offset=offset)
	slope = float(fit.params[1])
	se = float(fit.se[1])
	dispersion = float(fit.dispersion)
	if quasi and dispersion > 1.0:
		se *= math.sqrt(dispersion)
		p_value = two_sided_t_p(slope / se if se > 0 else float('nan'), fit.df_resid)
	else:
		p_value = two_sided_normal_p(slope / se) if se > 0 else float('nan')
	critical = -norm_ppf((1 - confidence) / 2.0)
	out.update({'growth_rate_per_year': slope, 'std_error': se,
				'ci_low': slope - critical * se, 'ci_high': slope + critical * se,
				'p_value': p_value, 'doubling_time_years': doubling_time(slope),
				'dispersion': dispersion, 'converged': bool(fit.converged)})
	return out


def multinomial_growth(times, labels, pivot=None, confidence=0.95, max_iter=200, tol=1e-9,
					   prior_sd=5.0):
	"""Multinomial logistic growth rates for several competing variants at once.

	This is the model behind the variant growth-advantage estimates published
	through the pandemic (and `evofr`'s MLR): every variant's frequency is
	modelled jointly against one pivot, so the estimates are constrained to sum
	to a valid set of frequencies. Fitting each variant separately against "all
	the others" cannot do that - the others are themselves growing, so each
	pairwise fit silently uses a different, moving baseline, and a set of such
	estimates need not be mutually consistent.

	A weak Gaussian prior of standard deviation `prior_sd` (on the log-odds per
	year scale) is placed on the slopes. It is not cosmetic: a residue seen
	twice, both times late, is perfectly separated, the unpenalised maximum
	likelihood for it is infinite, and without the prior the fit returns
	coefficients of 1e31 with a p-value of exactly zero - a number that looks
	like the strongest result in the table and means nothing. A prior sd of 5
	per year is far wider than any real growth advantage, so it changes a
	well-determined estimate by a negligible amount and is decisive only where
	the data determine nothing at all.

	Returns ``(records, info)``: one record per category, with the pivot's
	growth advantage fixed at zero by construction.
	"""
	times = np.asarray(times, dtype=float)
	labels = np.asarray(labels, dtype=object)
	keep = np.isfinite(times) & np.array([lab is not None and str(lab) != '' for lab in labels])
	times, labels = times[keep], labels[keep]

	categories, counts = np.unique(labels.astype(str), return_counts=True)
	info = {'n': int(len(times)), 'n_categories': int(len(categories)),
			'pivot': None, 'converged': False, 'iterations': 0}
	if len(categories) < 2 or len(times) < 4 or np.ptp(times) <= 0:
		return [], info

	if pivot is None or pivot not in set(categories):
		pivot = str(categories[int(np.argmax(counts))])
	info['pivot'] = pivot
	others = [c for c in categories if c != pivot]

	centre = float(np.mean(times))
	design = np.column_stack([np.ones_like(times), times - centre])
	n, p = design.shape
	k = len(others)
	indicator = np.column_stack([(labels.astype(str) == c).astype(float) for c in others])

	beta = np.zeros((k, p))
	converged = False
	covariance = np.eye(k * p)
	# The prior applies to the slopes only; an intercept is a base frequency and
	# has no scale a prior could sensibly be set on.
	penalty = np.zeros(k * p)
	if prior_sd and prior_sd > 0:
		penalty[1::p] = 1.0 / (prior_sd ** 2)
	info['prior_sd'] = float(prior_sd) if prior_sd else None
	for iteration in range(1, max_iter + 1):
		eta = design.dot(beta.T)
		eta = np.clip(eta, -500, 500)
		exp_eta = np.exp(eta)
		denom = 1.0 + np.sum(exp_eta, axis=1)
		probs = exp_eta / denom[:, None]
		flat = beta.reshape(-1)
		gradient = (design.T.dot(indicator - probs)).T.reshape(-1) - penalty * flat
		hessian = np.diag(penalty).astype(float)
		for a in range(k):
			for b in range(k):
				weight = probs[:, a] * ((1.0 if a == b else 0.0) - probs[:, b])
				block = design.T.dot(design * weight[:, None])
				hessian[a * p:(a + 1) * p, b * p:(b + 1) * p] += block
		step, covariance, _ = _solve_spd(hessian, gradient)
		step = step.reshape(k, p)
		# Reject a step that leaves the region the prior makes plausible rather
		# than taking a scaled version of it: a halved infinite step is still
		# infinite, and the loop used to accept one once it ran out of halvings.
		scale = 1.0
		while scale > 1e-6 and np.max(np.abs(beta + scale * step)) > 50.0:
			scale /= 2.0
		if np.max(np.abs(beta + scale * step)) > 50.0:
			break
		beta = beta + scale * step
		info['iterations'] = iteration
		if np.max(np.abs(scale * step)) < tol:
			converged = True
			break
	info['converged'] = converged

	critical = -norm_ppf((1 - confidence) / 2.0)
	standard_errors = np.sqrt(np.maximum(np.diag(covariance), 0.0)).reshape(k, p)
	records = [{'category': pivot, 'is_pivot': True, 'n': int(counts[list(categories).index(pivot)]),
				'growth_advantage_per_year': 0.0, 'std_error': float('nan'),
				'ci_low': float('nan'), 'ci_high': float('nan'), 'p_value': float('nan')}]
	for index, category in enumerate(others):
		slope = float(beta[index, 1])
		se = float(standard_errors[index, 1])
		records.append({
			'category': str(category), 'is_pivot': False,
			'n': int(counts[list(categories).index(category)]),
			'growth_advantage_per_year': slope, 'std_error': se,
			'ci_low': slope - critical * se, 'ci_high': slope + critical * se,
			'p_value': two_sided_normal_p(slope / se) if se > 0 else float('nan'),
		})
	return records, info


def frequency_increment_test(times, frequencies, sample_sizes=None):
	"""Frequency Increment Test (Feder, Kryazhimskiy & Plotkin 2014).

	Rescales successive frequency increments by the drift variance expected
	between those two time points and asks whether the rescaled increments have
	mean zero. Under neutrality they do, whatever the demography; a consistent
	upward drift in a trajectory is evidence of selection that does not depend
	on the trajectory being logistic.

	It is the right companion to the logistic fit precisely because it assumes
	so much less: a trajectory that rises, plateaus and falls has no meaningful
	logistic slope but still gives an interpretable FIT statistic over its
	rising phase.
	"""
	times = np.asarray(times, dtype=float)
	frequencies = np.asarray(frequencies, dtype=float)
	order = np.argsort(times, kind='mergesort')
	times, frequencies = times[order], frequencies[order]
	if sample_sizes is not None:
		sample_sizes = np.asarray(sample_sizes, dtype=float)[order]
		# Agresti-Coull style shrinkage keeps a bin that happens to be all-carrier
		# in the series instead of deleting it, which on short series is most of
		# the data.
		frequencies = (frequencies * sample_sizes + 0.5) / (sample_sizes + 1.0)

	out = {'n_time_points': 0, 'fit_statistic': float('nan'), 'p_value': float('nan'),
		   'mean_increment': float('nan'), 'df': 0}
	usable = np.isfinite(times) & np.isfinite(frequencies) & (frequencies > 0) & (frequencies < 1)
	times, frequencies = times[usable], frequencies[usable]
	if len(times) < 3:
		return out

	delta_t = np.diff(times)
	good = delta_t > 0
	if np.sum(good) < 2:
		return out
	numerator = np.diff(frequencies)[good]
	denominator = np.sqrt(2.0 * frequencies[:-1][good] * (1.0 - frequencies[:-1][good]) * delta_t[good])
	increments = numerator / denominator
	increments = increments[np.isfinite(increments)]
	if len(increments) < 2:
		return out

	mean = float(np.mean(increments))
	sd = float(np.std(increments, ddof=1))
	df = len(increments) - 1
	if sd <= 0:
		return out
	statistic = mean / (sd / math.sqrt(len(increments)))
	out.update({'n_time_points': int(len(times)), 'fit_statistic': statistic,
				'p_value': two_sided_t_p(statistic, df), 'mean_increment': mean, 'df': df})
	return out


def bin_trajectory(times, carrier, bin_width=1.0, origin=None, confidence=0.95):
	"""Bin a carrier/background series into a frequency trajectory.

	Returns a list of per-bin dicts (midpoint, counts, frequency, Wilson
	interval). The bins are what the binned logistic fit, the Poisson count fit,
	the FIT and the trend tests all consume, so they are produced once here.
	"""
	times = np.asarray(times, dtype=float)
	carrier = np.asarray(carrier, dtype=float)
	keep = np.isfinite(times) & np.isfinite(carrier)
	times, carrier = times[keep], carrier[keep]
	if len(times) == 0:
		return []
	if bin_width <= 0:
		raise ValueError('bin_width must be positive')
	origin = float(np.min(times)) if origin is None else float(origin)

	index = np.floor((times - origin) / bin_width).astype(int)
	rows = []
	for value in np.unique(index):
		mask = index == value
		total = int(np.sum(mask))
		successes = int(np.sum(carrier[mask] > 0.5))
		low, high = wilson_interval(successes, total, confidence)
		rows.append({
			'bin_start': origin + value * bin_width,
			'bin_mid': origin + (value + 0.5) * bin_width,
			'n': total,
			'n_carrier': successes,
			'frequency': successes / total if total else float('nan'),
			'ci_low': low,
			'ci_high': high,
			'mean_time': float(np.mean(times[mask])),
		})
	return rows
