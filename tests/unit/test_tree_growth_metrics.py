import math

import numpy as np
import pytest

import tree_growth_metrics as tgm


SMALL = "((A:1,B:1)n2:2,(C:1,(D:0.5,E:0.5)n4:0.5)n3:2)n1:0;"


# --- parsing -------------------------------------------------------------

def test_parses_tips_internal_labels_and_branch_lengths():
	tree = tgm.parse_newick(SMALL)
	assert sorted(tree.name[tip] for tip in tree.tips) == ['A', 'B', 'C', 'D', 'E']
	assert tree.name[tree.root] == 'n1'
	assert tree.length[tree.name_to_node['D']] == 0.5
	assert dict(zip(tree.name, tree.tip_counts()))['n3'] == 3


def test_parses_a_tree_with_no_branch_lengths_and_quoted_labels():
	tree = tgm.parse_newick("(('a b','c,d'),e);")
	assert sorted(tree.name[tip] for tip in tree.tips) == ['a b', 'c,d', 'e']
	assert np.all(np.isnan(tree.length[tree.tips]))


def test_postorder_places_every_node_after_its_descendants():
	tree = tgm.parse_newick(SMALL)
	position = {node: index for index, node in enumerate(tree.postorder)}
	for node in range(tree.n_nodes):
		for kid in tree.children[node]:
			assert position[kid] < position[node]


def test_unbalanced_parentheses_are_rejected():
	with pytest.raises(ValueError):
		tgm.parse_newick("((A,B);")


# --- restriction ---------------------------------------------------------

def test_induced_subtree_preserves_path_lengths_and_suppresses_unifurcations():
	tree = tgm.parse_newick(SMALL)
	keep = [tree.name_to_node[name] for name in ('A', 'D', 'E')]
	sub, source = tgm.induced_subtree(tree, keep)
	lengths = dict(zip(sub.name, sub.length))
	# A was 1 below n2 which was 2 below the root; n2 has one kept child and goes.
	assert lengths['A'] == 3.0
	assert lengths['n4'] == 2.5
	assert sorted(sub.name[tip] for tip in sub.tips) == ['A', 'D', 'E']
	assert tree.name[source[0]] == sub.name[0]


def test_induced_subtree_reroots_on_the_mrca_of_what_is_kept():
	tree = tgm.parse_newick(SMALL)
	sub, _ = tgm.induced_subtree(tree, [tree.name_to_node[name] for name in ('D', 'E')])
	assert sub.name[sub.root] == 'n4'


# --- dating --------------------------------------------------------------

DATES = {'A': 2010.0, 'B': 2012.0, 'C': 2015.0, 'D': 2018.0, 'E': 2019.0}


def test_mrca_bound_dates_each_node_to_its_earliest_descendant():
	tree = tgm.parse_newick(SMALL)
	times, info = tgm.date_nodes(tree, DATES, 'mrca-bound')
	by_name = dict(zip(tree.name, times))
	assert by_name['n4'] == 2018.0
	assert by_name['n1'] == 2010.0
	assert info['method'] == 'mrca-bound'


def test_root_to_tip_recovers_a_planted_clock_rate():
	# Divergence built as exactly 0.01 per year since a root at 2000.
	newick = "((A:0.10,B:0.12):0.05,(C:0.14,D:0.16):0.05);"
	dates = {'A': 2015.0, 'B': 2017.0, 'C': 2019.0, 'D': 2021.0}
	tree = tgm.parse_newick(newick)
	times, info = tgm.date_nodes(tree, dates, 'root-to-tip')
	# Too few dated tips for the clock to be trusted, so it must fall back.
	assert info['method'] == 'mrca-bound'

	# With enough tips it is used, and recovers the rate and the root date.
	tips = []
	for index in range(40):
		year = 2000.0 + index * 0.5
		tips.append('T%d:%0.4f' % (index, 0.01 * (year - 1990.0)))
		dates['T%d' % index] = year
	tree = tgm.parse_newick('(' + ','.join(tips) + ');')
	times, info = tgm.date_nodes(tree, dates, 'root-to-tip')
	assert info['method'] == 'root-to-tip'
	assert abs(info['clock_rate'] - 0.01) < 1e-6
	assert abs(info['root_time'] - 1990.0) < 1e-6


def test_branch_length_dating_shifts_the_root_onto_the_tip_dates():
	tree = tgm.parse_newick("((A:5,B:7):3,(C:9,D:11):3);")
	dates = {'A': 2008.0, 'B': 2010.0, 'C': 2012.0, 'D': 2014.0}
	times, info = tgm.date_nodes(tree, dates, 'branch-lengths')
	assert info['method'] == 'branch-lengths'
	# Branch lengths are taken as years; the root shift is the least-squares fit
	# of the tips onto their dates, which here is exact.
	assert abs(times[tree.name_to_node['A']] - 2008.0) < 1e-6
	assert abs(times[tree.name_to_node['D']] - 2014.0) < 1e-6
	assert times[tree.root] < times[tree.name_to_node['A']]


def test_time_branch_lengths_are_never_negative():
	tree = tgm.parse_newick(SMALL)
	times, _ = tgm.date_nodes(tree, DATES, 'mrca-bound')
	lengths = tgm.time_branch_lengths(tree, times)
	assert np.all(lengths >= 0)


# --- independent origins -------------------------------------------------

def test_a_derived_trait_in_one_clade_is_one_origin():
	tree = tgm.parse_newick("((((A,B),(C,D)))x,E)r;")
	result = tgm.reconstruct_origins(tree, {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 0})
	assert result['root_state'] == 'absent'
	assert (result['origins_min'], result['origins_max']) == (1, 1)
	assert len(result['clusters']) == 1


def test_a_wild_type_residue_is_reconstructed_at_the_root_not_as_many_origins():
	# On a polytomy, forcing the root to absent makes each carrier its own
	# origin; letting parsimony choose recognises the ancestral state instead.
	tree = tgm.parse_newick("(A,B,C,D,E)r;")
	states = {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 0}
	auto = tgm.reconstruct_origins(tree, states)
	forced = tgm.reconstruct_origins(tree, states, root_state='absent')
	assert auto['root_state'] == 'present'
	assert auto['origins_min'] == 1 and auto['losses_min'] == 1
	assert forced['origins_min'] == 4


def test_tips_with_no_state_constrain_nothing_but_can_be_spanned():
	tree = tgm.parse_newick("(A,B,C,U1,U2,U3)r;")
	result = tgm.reconstruct_origins(tree, {'A': 1, 'B': 1, 'C': 0})
	assert result['parsimony_cost'] == 1.0


def test_the_origin_range_widens_when_parsimony_is_ambiguous():
	# Either the ancestor was present and lost twice, or gained twice; both
	# cost 2, so the number of gains is not identified.
	tree = tgm.parse_newick("((A,B)x,(C,D)y)r;")
	result = tgm.reconstruct_origins(tree, {'A': 1, 'B': 0, 'C': 1, 'D': 0})
	assert result['parsimony_cost'] == 2.0
	assert (result['origins_min'], result['origins_max']) == (1, 2)
	assert (result['losses_min'], result['losses_max']) == (0, 2)


def test_clusters_collect_the_tips_of_each_origin():
	tree = tgm.parse_newick("(((A,B)x,C)p,(D,E)y,F)r;")
	result = tgm.reconstruct_origins(tree, {'A': 1, 'B': 1, 'C': 0, 'D': 0, 'E': 0, 'F': 0})
	members = [sorted(tree.name[tip] for tip in cluster['tips'])
			   for cluster in result['clusters']]
	assert members == [['A', 'B']]


# --- rates ---------------------------------------------------------------

def test_local_branching_index_is_largest_where_most_branching_happens():
	# Same branch lengths either side, six tips against two: the LBI measures
	# how much branch length sits within tau of a node, so the busier clade wins.
	dense = ','.join('D%d:0.3' % index for index in range(6))
	tree = tgm.parse_newick("((%s)dense:1,(E:0.3,F:0.3)sparse:1)r:0;" % dense)
	raw, normalised = tgm.local_branching_index(tree, np.nan_to_num(tree.length), tau=1.0)
	by_name = dict(zip(tree.name, normalised))
	assert by_name['dense'] > by_name['sparse']
	assert abs(np.max(normalised) - 1.0) < 1e-12
	assert np.all(raw >= 0)


def test_magallon_sanderson_reduces_to_its_closed_forms_when_extinction_is_zero():
	assert abs(tgm.magallon_sanderson(100, 10, 0.0, crown=True) - math.log(50) / 10) < 1e-9
	assert abs(tgm.magallon_sanderson(100, 10, 0.0, crown=False) - math.log(100) / 10) < 1e-9
	# A higher assumed extinction fraction always implies a slower net rate.
	assert tgm.magallon_sanderson(100, 10, 0.9) < tgm.magallon_sanderson(100, 10, 0.0)


def test_yule_rate_and_gamma_on_a_simulated_pure_birth_tree():
	rng = np.random.default_rng(11)
	lineages, elapsed, total = 2, 0.0, 0.0
	branching = []
	while lineages < 400:
		elapsed += rng.exponential(1.0 / (lineages * 0.7))
		total += lineages * (elapsed - (branching[-1] if branching else 0.0))
		branching.append(elapsed)
		lineages += 1
	total = 0.0
	previous = 0.0
	for index, time in enumerate(branching):
		total += (index + 2) * (time - previous)
		previous = time
	assert abs(tgm.yule_birth_rate(lineages, total) - 0.7) < 0.1
	# Gamma is standard normal under pure birth, so a Yule tree sits near zero.
	assert abs(tgm.gamma_statistic([elapsed - time for time in branching])) < 3


def test_coalescent_recovers_a_simulated_exponential_growth_rate():
	rng = np.random.default_rng(5)
	n0, rate, lineages, tau = 1000.0, 1.5, 60, 0.0
	intervals, coalescences = [], []
	while lineages > 1:
		pairs = lineages * (lineages - 1) / 2
		draw = -math.log(rng.random()) * n0 * rate / pairs + math.exp(rate * tau)
		nxt = math.log(draw) / rate
		intervals.append((tau, nxt, lineages))
		coalescences.append(nxt)
		tau, lineages = nxt, lineages - 1
	result = tgm.coalescent_exponential_growth(intervals, coalescences)
	assert result['converged']
	assert abs(result['growth_rate_per_year'] - 1.5) < 0.3
	assert result['ci_low'] < 1.5 < result['ci_high']


def test_lineages_through_time_removes_a_lineage_at_each_serial_sample():
	tree = tgm.parse_newick("((A,B)x,(C,D)y)r;")
	times = np.array([2000.0, 2002.0, 2004.0, 2010.0, 2012.0, 2011.0, 2013.0])
	# root, x, A, B, y, C, D in node order from the parser.
	order = {name: index for index, name in enumerate(tree.name)}
	node_times = np.zeros(tree.n_nodes)
	for name, value in zip(['r', 'x', 'A', 'B', 'y', 'C', 'D'], times):
		node_times[order[name]] = value
	stamps, counts = tgm.lineages_through_time(tree, node_times)
	assert counts[-1] == 0
	assert max(counts) <= 4


def test_clade_scan_finds_a_planted_fast_growing_clade():
	# One clade sampled only late, one sampled throughout.
	fast = ['F%d' % i for i in range(40)]
	slow = ['S%d' % i for i in range(60)]
	tree = tgm.parse_newick('((%s)fast,(%s)slow)r;' % (','.join(fast), ','.join(slow)))
	tip_bin = np.full(tree.n_nodes, -1, dtype=np.int64)
	for index, name in enumerate(fast):
		tip_bin[tree.name_to_node[name]] = 6 + index % 4
	for index, name in enumerate(slow):
		tip_bin[tree.name_to_node[name]] = index % 10
	counts = tgm.clade_bin_counts(tree, tip_bin, 10)
	assert counts[tree.root].sum() == 100
	results = tgm.scan_clade_growth(tree, counts, np.arange(10.0) + 0.5,
									min_clade_size=10, max_clade_fraction=0.95)
	best = results[0]
	assert best['node_label'] == 'fast'
	assert best['growth_advantage_per_year'] > 0


def test_bin_counting_refuses_an_allocation_that_would_not_fit():
	tree = tgm.parse_newick("(A,B);")
	with pytest.raises(MemoryError):
		tgm.clade_bin_counts(tree, np.zeros(tree.n_nodes, dtype=np.int64),
							 tgm.MAX_BIN_CELLS)


def test_expansion_score_is_high_for_a_recently_sampled_clade():
	tree = tgm.parse_newick("((A,B)recent,(C,D)old)r;")
	tip_time = np.full(tree.n_nodes, np.nan)
	for name, value in (('A', 2020.0), ('B', 2021.0), ('C', 2012.0), ('D', 2013.0)):
		tip_time[tree.name_to_node[name]] = value
	scores = tgm.expansion_score(tree, tip_time, window=5.0)
	by_name = dict(zip(tree.name, scores))
	assert by_name['recent'] > by_name['old']


def test_skyline_and_the_coalescent_ml_agree_on_the_same_genealogy():
	# Two independent readings of the same growth rate: the maximum-likelihood
	# fit of one exponential, and the slope of log(Ne) across the skyline. They
	# must agree in sign and magnitude on data that really is exponential -
	# which is what makes a disagreement on real data informative rather than a
	# bug.
	rng = np.random.default_rng(5)
	n0, rate, lineages, tau = 1000.0, 1.5, 80, 0.0
	intervals, coalescences = [], []
	while lineages > 1:
		pairs = lineages * (lineages - 1) / 2
		draw = -math.log(rng.random()) * n0 * rate / pairs + math.exp(rate * tau)
		nxt = math.log(draw) / rate
		intervals.append((tau, nxt, lineages))
		coalescences.append(nxt)
		tau, lineages = nxt, lineages - 1

	maximum_likelihood = tgm.coalescent_exponential_growth(
		intervals, coalescences)['growth_rate_per_year']
	points = tgm.skyline(intervals, coalescences, group=10)
	taus = np.array([point['tau'] for point in points])
	sizes = np.log(np.array([point['ne'] for point in points]))
	# tau runs backwards, so a population growing forwards shrinks into the past.
	from_skyline = -float(np.polyfit(taus, sizes, 1)[0])

	assert from_skyline > 0 and maximum_likelihood > 0
	assert abs(from_skyline - maximum_likelihood) < 0.5
	assert abs(from_skyline - 1.5) < 0.5


def test_skyline_grouping_reduces_the_number_of_points():
	intervals = [(float(i), float(i + 1), 10) for i in range(30)]
	coalescences = [float(i + 1) for i in range(30)]
	assert len(tgm.skyline(intervals, coalescences, group=1)) == 30
	assert len(tgm.skyline(intervals, coalescences, group=10)) == 3
