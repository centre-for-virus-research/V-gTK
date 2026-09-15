import math

import numpy as np

import growth_stats as gs


# --- distributions -------------------------------------------------------
# Checked against published table values rather than against another
# implementation, because there is no other implementation in this environment.

def test_normal_quantile_matches_the_textbook_975_point():
	assert abs(gs.norm_ppf(0.975) - 1.959963985) < 1e-8
	assert abs(gs.norm_ppf(0.5)) < 1e-12
	assert gs.norm_ppf(0.0) == float('-inf')


def test_chi_squared_upper_tail():
	# The 5% point of chi-squared on 1 df is 3.841459.
	assert abs(gs.chi2_sf(3.841459, 1) - 0.05) < 1e-6
	assert abs(gs.chi2_sf(10.0, 3) - 0.018566) < 1e-5
	assert gs.chi2_sf(0.0, 1) == 1.0


def test_student_t_two_sided():
	# t = 2.228 on 10 df is the 5% two-sided point.
	assert abs(gs.two_sided_t_p(2.228, 10) - 0.05) < 1e-4


def test_fisher_exact_reproduces_the_tea_tasting_table():
	odds, p_value = gs.fisher_exact_2x2(3, 1, 1, 3)
	assert abs(p_value - 0.4857142857) < 1e-9
	assert abs(odds - 9.0) < 1e-9


def test_fisher_exact_corrects_an_empty_cell_instead_of_returning_infinity():
	odds, p_value = gs.fisher_exact_2x2(10, 0, 0, 10)
	assert math.isfinite(odds) and odds > 1
	assert p_value < 1e-4


def test_benjamini_hochberg_is_monotone_and_keeps_nan():
	q = gs.benjamini_hochberg([0.01, 0.02, 0.03, 0.9, float('nan')])
	assert np.all(np.diff(q[:4]) >= -1e-12)
	assert np.isnan(q[4])
	# A NaN must not count towards the number of tests.
	assert abs(q[0] - 0.04) < 1e-12


def test_wilson_interval_stays_inside_the_unit_interval_at_the_boundary():
	low, high = gs.wilson_interval(0, 10)
	# Floating point leaves a value at the 1e-17 level rather than a clean zero.
	assert low < 1e-12 and 0 < high < 1


# --- GLM -----------------------------------------------------------------

def test_binomial_glm_recovers_known_coefficients():
	rng = np.random.default_rng(7)
	times = rng.uniform(0, 10, 4000)
	probability = 1 / (1 + np.exp(-(-3 + 0.8 * times)))
	outcome = (rng.random(4000) < probability).astype(float)
	design = np.column_stack([np.ones_like(times), times])
	fit = gs.fit_glm(design, outcome, 'binomial')
	assert fit.converged
	assert abs(fit.params[0] + 3) < 0.25
	assert abs(fit.params[1] - 0.8) < 0.05


def test_poisson_glm_recovers_known_coefficients_and_unit_dispersion():
	rng = np.random.default_rng(11)
	times = rng.uniform(0, 8, 3000)
	counts = rng.poisson(np.exp(1 + 0.3 * times)).astype(float)
	design = np.column_stack([np.ones_like(times), times])
	fit = gs.fit_glm(design, counts, 'poisson')
	assert abs(fit.params[1] - 0.3) < 0.02
	assert 0.8 < fit.dispersion < 1.2


def test_firth_rescues_a_perfectly_separated_fit():
	# Every late observation is a carrier and every early one is not: the
	# maximum-likelihood slope is infinite and its Wald test is powerless.
	times = np.array([1.0, 2, 3, 4, 5, 6])
	outcome = np.array([0.0, 0, 0, 1, 1, 1])
	design = np.column_stack([np.ones(6), times])

	plain = gs.fit_glm(design, outcome, 'binomial')
	assert plain.params[1] > 10 and plain.se[1] > 100
	assert gs.two_sided_normal_p(plain.params[1] / plain.se[1]) > 0.9

	penalised = gs.fit_glm(design, outcome, 'binomial', firth=True)
	assert 0 < penalised.params[1] < 5
	assert penalised.se[1] < 5
	assert gs.likelihood_ratio_p(design, outcome, 1, firth=True)[0] < 0.05


# --- growth estimators ---------------------------------------------------

def test_logistic_growth_recovers_a_known_selection_coefficient():
	rng = np.random.default_rng(3)
	times = rng.uniform(2010, 2020, 5000)
	logit = -0.9 * 2015 - 1 + 0.9 * times
	carrier = (rng.random(5000) < 1 / (1 + np.exp(-logit))).astype(float)
	result = gs.logistic_growth(times, carrier)
	assert abs(result['growth_advantage_per_year'] - 0.9) < 0.06
	assert result['ci_low'] < 0.9 < result['ci_high']
	assert result['lrt_p'] < 1e-10
	assert result['method'] == 'ml'
	assert abs(result['doubling_time_years'] - math.log(2) / result['growth_advantage_per_year']) < 1e-9


def test_logistic_growth_uses_firth_for_a_rare_allele():
	times = np.concatenate([np.linspace(2000, 2020, 200), np.linspace(2018, 2020, 4)])
	carrier = np.concatenate([np.zeros(200), np.ones(4)])
	result = gs.logistic_growth(times, carrier)
	assert result['method'] == 'firth'
	assert math.isfinite(result['growth_advantage_per_year'])


def test_quasi_binomial_widens_the_interval_when_bins_are_overdispersed():
	times = np.arange(2000.0, 2020.0)
	trials = np.full(20, 100.0)
	# Frequencies scattered far more than binomial sampling allows.
	rng = np.random.default_rng(5)
	base = 1 / (1 + np.exp(-(-0.2 * 2010 + 0.2 * times)))
	noisy = np.clip(base + rng.normal(0, 0.15, 20), 0.02, 0.98)
	plain = gs.logistic_growth(times, noisy, trials=trials, quasi=False)
	quasi = gs.logistic_growth(times, noisy, trials=trials, quasi=True)
	assert quasi['dispersion'] > 2
	assert quasi['std_error'] > plain['std_error']


def test_multinomial_growth_is_consistent_and_survives_a_separated_category():
	rng = np.random.default_rng(3)
	times = list(rng.uniform(2010, 2020, 300))
	labels = ['L'] * 300
	times += [2019.9, 2019.95]          # a rare, perfectly separated residue
	labels += ['E', 'E']
	times += list(rng.uniform(2010, 2020, 80))
	labels += ['M'] * 80
	records, info = gs.multinomial_growth(np.array(times), np.array(labels, dtype=object))
	assert info['converged'] and info['pivot'] == 'L'
	by_category = {record['category']: record for record in records}
	assert by_category['L']['growth_advantage_per_year'] == 0.0
	assert abs(by_category['M']['growth_advantage_per_year']) < 0.3
	# Without the prior this comes back as 1e31 with a p-value of exactly zero.
	assert 0 < by_category['E']['growth_advantage_per_year'] < 20


def test_frequency_increment_test_separates_a_driven_trajectory_from_a_flat_one():
	times = np.arange(0.0, 12.0)
	rising = 1 / (1 + np.exp(-(-6 + 1.0 * times)))
	sizes = np.full(12, 200.0)
	driven = gs.frequency_increment_test(times, rising, sizes)
	flat = gs.frequency_increment_test(times, np.full(12, 0.4), sizes)
	assert driven['p_value'] < 0.05 and driven['mean_increment'] > 0
	assert not (flat['p_value'] < 0.05)


def test_mann_kendall_detects_a_monotone_rise_and_reports_sens_slope():
	times = np.arange(10.0)
	result = gs.mann_kendall(times, 0.05 * times + 0.1)
	assert result['tau'] == 1.0
	assert result['p_value'] < 0.01
	assert abs(result['sen_slope'] - 0.05) < 1e-9


def test_reproduction_number_is_smaller_for_a_dispersed_generation_interval():
	fixed = gs.reproduction_number(0.5, 1.0, 0.0)
	dispersed = gs.reproduction_number(0.5, 1.0, 0.5)
	assert abs(fixed - math.exp(0.5)) < 1e-12
	assert dispersed < fixed


def test_bin_trajectory_counts_and_bounds_each_bin():
	times = np.array([2000.1, 2000.9, 2001.5, 2001.6, 2002.2])
	carrier = np.array([1.0, 0, 1, 1, 0])
	bins = gs.bin_trajectory(times, carrier, 1.0, origin=2000.0)
	assert [entry['n'] for entry in bins] == [2, 2, 1]
	assert [entry['n_carrier'] for entry in bins] == [1, 2, 0]
	assert all(entry['ci_low'] <= entry['frequency'] <= entry['ci_high'] for entry in bins)


def test_aggregating_identical_dates_is_exact_not_approximate():
	# Above the aggregation threshold the individual rows are collapsed onto
	# their distinct sampling dates. The binomial logit likelihood for grouped
	# counts is the same function as the one for the rows behind them, so the
	# estimate must match the ungrouped fit to numerical precision. That is what
	# licenses the speedup; an approximation here would quietly bias every rate.
	rng = np.random.default_rng(17)
	dates = np.repeat(np.arange(2000.0, 2020.0), 400)      # 8000 rows, 20 dates
	probability = 1 / (1 + np.exp(-(-0.4 * 2010 + 0.4 * dates)))
	carrier = (rng.random(dates.size) < probability).astype(float)

	grouped = gs.logistic_growth(dates, carrier)
	assert grouped['n'] == dates.size

	saved = gs._AGGREGATE_ABOVE
	try:
		gs._AGGREGATE_ABOVE = 10 ** 9                       # force the ungrouped path
		individual = gs.logistic_growth(dates, carrier)
	finally:
		gs._AGGREGATE_ABOVE = saved

	assert grouped['method'] == individual['method'] == 'ml'
	assert abs(grouped['growth_advantage_per_year']
			   - individual['growth_advantage_per_year']) < 1e-8
	assert abs(grouped['std_error'] - individual['std_error']) < 1e-8
	assert abs(grouped['lrt_p'] - individual['lrt_p']) < 1e-10


def test_firth_is_also_invariant_to_grouping():
	# Firth's penalty is a function of the Fisher information, which grouping
	# does not change, so the penalised fit must match the ungrouped one. This is
	# what licenses aggregating the rare alleles - the ones that need the penalty
	# and are also the large majority of any protein.
	rng = np.random.default_rng(4)
	dates = np.repeat(np.arange(2000.0, 2025.0), 400)
	carrier = np.zeros(dates.size)
	late = np.flatnonzero(dates >= 2022)
	carrier[rng.choice(late, 12, replace=False)] = 1.0

	design = np.column_stack([np.ones_like(dates), dates - dates.mean()])
	individual = gs.fit_glm(design, carrier, 'binomial', firth=True)

	unique, inverse = np.unique(dates, return_inverse=True)
	counts = np.bincount(inverse).astype(float)
	grouped = gs.fit_glm(np.column_stack([np.ones_like(unique), unique - dates.mean()]),
						 np.bincount(inverse, weights=carrier) / counts,
						 'binomial', trials=counts, firth=True)

	assert abs(individual.params[1] - grouped.params[1]) < 1e-8
	assert abs(individual.se[1] - grouped.se[1]) < 1e-8
	# And the estimate is finite and positive, which is the point of Firth here.
	assert 0 < grouped.params[1] < 10
