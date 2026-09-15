#!/usr/bin/env python3
"""Tree-side machinery for lineage growth-rate estimation.

Three things live here, in order:

1. **A newick reader** that produces flat arrays rather than an object graph.
   The UShER trees in a V-gTK database routinely carry >130,000 tips, and a
   recursive Bio.Phylo clade tree of that size costs a large multiple of the
   memory and cannot be traversed without raising the recursion limit. Arrays
   also make every later pass a single loop with no attribute lookups.

2. **Node dating.** Almost every phylogenetic growth estimator needs node
   *times*, and an UShER newick does not have them - its branch lengths are
   mutation counts, and in a database built from a de-duplicated alignment most
   of them are zero. Three sources of node times are supported and the one used
   is always recorded in the output, because they are not interchangeable.

3. **The estimators themselves**: parsimony reconstruction of independent
   origins, the local branching index, lineages-through-time, the
   method-of-moments diversification rates, the coalescent skyline and the
   exponential-growth coalescent.

Nothing here touches the database.
"""

import math
import re

import numpy as np

import growth_stats as gs

# A quoted label may legally contain commas and parentheses, so the quoted
# form has to be matched before the structural characters are - otherwise a
# tip called 'Congo, DR/2013' silently becomes two tips.
_TOKEN = re.compile(r"'(?:[^']|'')*'[^(),;]*|[(),;]|[^(),;]+")

#: Newick node labels may be single-quoted, in which case '' is a literal quote.
_QUOTED = re.compile(r"^'(.*)'$", re.DOTALL)


class Tree(object):
	"""A rooted tree held as parallel arrays.

	`parent[i]` is -1 for the root. `children[i]` is a list. `postorder` lists
	every node after all of its descendants; `preorder` is its reverse ordering
	by construction, which is all the traversals below need.
	"""

	def __init__(self, parent, children, name, length):
		self.parent = np.asarray(parent, dtype=np.int64)
		self.children = children
		self.name = name
		self.length = np.asarray(length, dtype=float)
		self.n_nodes = len(parent)
		self.root = int(np.flatnonzero(self.parent < 0)[0]) if self.n_nodes else -1
		self.is_tip = np.array([len(kids) == 0 for kids in children], dtype=bool)
		self.postorder = self._build_postorder()
		self.tips = np.flatnonzero(self.is_tip)
		self.name_to_node = {}
		for index, label in enumerate(name):
			if label and self.is_tip[index]:
				self.name_to_node.setdefault(label, index)

	def _build_postorder(self):
		"""Iterative post-order. Recursion overflows on ladder-shaped trees."""
		order = []
		if self.n_nodes == 0:
			return np.zeros(0, dtype=np.int64)
		stack = [self.root]
		while stack:
			node = stack.pop()
			order.append(node)
			stack.extend(self.children[node])
		order.reverse()
		return np.asarray(order, dtype=np.int64)

	@property
	def preorder(self):
		return self.postorder[::-1]

	def tip_counts(self):
		"""Number of tips at or below every node."""
		counts = self.is_tip.astype(np.int64)
		for node in self.postorder:
			if not self.is_tip[node]:
				counts[node] = sum(counts[kid] for kid in self.children[node])
		return counts

	def depths(self):
		"""Cumulative branch length from the root to every node."""
		depth = np.zeros(self.n_nodes, dtype=float)
		lengths = np.nan_to_num(self.length, nan=0.0)
		for node in self.preorder:
			if self.parent[node] >= 0:
				depth[node] = depth[self.parent[node]] + lengths[node]
		return depth


def _clean_label(text):
	label = text.strip()
	match = _QUOTED.match(label)
	if match:
		label = match.group(1).replace("''", "'")
	return label


def _split_label_length(token):
	"""Split ``label:length`` while leaving a quoted label containing ':' alone."""
	text = token.strip()
	if not text:
		return '', float('nan')
	if text.startswith("'"):
		closing = text.find("'", 1)
		while closing != -1 and closing + 1 < len(text) and text[closing + 1] == "'":
			closing = text.find("'", closing + 2)
		if closing != -1:
			label = _clean_label(text[:closing + 1])
			rest = text[closing + 1:].strip()
			if rest.startswith(':'):
				try:
					return label, float(rest[1:].split('[')[0])
				except ValueError:
					return label, float('nan')
			return label, float('nan')
	if ':' in text:
		label, _, number = text.partition(':')
		try:
			return _clean_label(label), float(number.split('[')[0])
		except ValueError:
			return _clean_label(label), float('nan')
	return _clean_label(text), float('nan')


def parse_newick(text):
	"""Parse a newick string into a :class:`Tree`."""
	if text is None:
		raise ValueError('no newick string supplied')
	text = str(text).strip()
	if not text:
		raise ValueError('empty newick string')

	parent, children, name, length = [], [], [], []

	def new_node(par):
		index = len(parent)
		parent.append(par if par is not None else -1)
		children.append([])
		name.append('')
		length.append(float('nan'))
		if par is not None:
			children[par].append(index)
		return index

	stack = []
	root = None
	just_closed = None
	for match in _TOKEN.finditer(text):
		token = match.group(0)
		if token == '(':
			node = new_node(stack[-1] if stack else None)
			if root is None:
				root = node
			stack.append(node)
			just_closed = None
		elif token == ',':
			just_closed = None
		elif token == ')':
			if not stack:
				raise ValueError('unbalanced parentheses in newick string')
			just_closed = stack.pop()
		elif token == ';':
			break
		else:
			if token.strip() == '':
				continue
			label, branch = _split_label_length(token)
			if just_closed is not None:
				name[just_closed] = label
				length[just_closed] = branch
				just_closed = None
			else:
				node = new_node(stack[-1] if stack else None)
				if root is None:
					root = node
				name[node] = label
				length[node] = branch
	if stack:
		raise ValueError('unbalanced parentheses in newick string')
	if root is None:
		raise ValueError('newick string contained no nodes')
	return Tree(parent, children, name, length)


# ---------------------------------------------------------------------------
# Node dating
# ---------------------------------------------------------------------------

DATING_METHODS = ('mrca-bound', 'root-to-tip', 'branch-lengths')

#: A root-to-tip regression on fewer dated tips than this is not worth
#: believing, and neither is one that explains less variance than
#: MIN_CLOCK_R_SQUARED. Both thresholds are conservative: the fallback is a
#: perfectly usable dating method, so there is no reason to accept a marginal
#: molecular clock in preference to it.
MIN_CLOCK_TIPS = 20
MIN_CLOCK_R_SQUARED = 0.1


def date_nodes(tree, tip_times, method='mrca-bound'):
	"""Assign a time to every node. Returns ``(times, info)``.

	``mrca-bound`` sets each internal node to the earliest sampling date among
	its descendants. That is a *bound*, not an estimate - the real common
	ancestor is older - and it is deliberately the default, because it is the
	only one of the three that needs nothing from the branch lengths. On an
	UShER tree written out of a de-duplicated alignment the branch lengths are
	mostly zero, and the two clock-based methods below would silently return a
	tree in which every node has the same age. Rates computed against it are
	biased towards being *too fast* (the clade looks younger than it is), and
	that bias is reported rather than hidden.

	``root-to-tip`` is the TempEst regression: fit divergence from the root
	against sampling date over the dated tips, read the substitution rate off
	the slope, and convert each node's divergence into a time. It needs
	informative branch lengths and a clock-like tree, and refuses when it does
	not have them.

	``branch-lengths`` takes the branch lengths to be times already, and shifts
	the root so the tips land on their sampling dates as closely as possible -
	for a tree that has come out of a dating program.
	"""
	if method not in DATING_METHODS:
		raise ValueError('unknown dating method: %s' % method)

	times = np.full(tree.n_nodes, np.nan)
	for node, label in enumerate(tree.name):
		if tree.is_tip[node] and label in tip_times:
			times[node] = tip_times[label]
	dated_tips = np.flatnonzero(np.isfinite(times) & tree.is_tip)
	info = {'method': method, 'n_dated_tips': int(len(dated_tips)), 'fallback': None,
			'clock_rate': float('nan'), 'r_squared': float('nan'),
			'root_time': float('nan')}
	if len(dated_tips) == 0:
		raise ValueError('no tip in this tree carries a usable date')

	if method in ('root-to-tip', 'branch-lengths'):
		depth = tree.depths()
		finite_lengths = np.isfinite(tree.length) & (tree.length > 0)
		if method == 'root-to-tip':
			if len(dated_tips) < MIN_CLOCK_TIPS or np.sum(finite_lengths) < 2:
				info['fallback'] = 'too few dated tips or uninformative branch lengths'
			else:
				x = times[dated_tips]
				y = depth[dated_tips]
				if np.ptp(x) <= 0 or np.ptp(y) <= 0:
					info['fallback'] = 'no variation in dates or in root-to-tip divergence'
				else:
					slope, intercept = np.polyfit(x, y, 1)
					correlation = float(np.corrcoef(x, y)[0, 1])
					info['clock_rate'] = float(slope)
					info['r_squared'] = correlation ** 2
					if slope <= 0 or info['r_squared'] < MIN_CLOCK_R_SQUARED:
						info['fallback'] = ('molecular clock too weak '
											'(rate=%.3g, R2=%.3f)' % (slope, info['r_squared']))
					else:
						root_time = -intercept / slope
						info['root_time'] = float(root_time)
						return root_time + depth / slope, info
		else:
			if np.sum(finite_lengths) < 2:
				info['fallback'] = 'branch lengths are all zero or missing'
			else:
				shift = float(np.mean(times[dated_tips] - depth[dated_tips]))
				info['root_time'] = shift
				return shift + depth, info
		info['method'] = 'mrca-bound'

	node_times = np.where(tree.is_tip, times, np.inf)
	for node in tree.postorder:
		if not tree.is_tip[node]:
			child_times = [node_times[kid] for kid in tree.children[node]]
			node_times[node] = min(child_times) if child_times else np.inf
	node_times[~np.isfinite(node_times)] = np.nan
	info['root_time'] = float(node_times[tree.root]) if np.isfinite(node_times[tree.root]) else float('nan')
	return node_times, info


def time_branch_lengths(tree, node_times):
	"""Branch lengths in time units, from dated nodes. Negatives are clamped."""
	lengths = np.zeros(tree.n_nodes, dtype=float)
	for node in range(tree.n_nodes):
		par = tree.parent[node]
		if par < 0:
			continue
		if np.isfinite(node_times[node]) and np.isfinite(node_times[par]):
			lengths[node] = max(0.0, node_times[node] - node_times[par])
	return lengths


# ---------------------------------------------------------------------------
# Independent origins by parsimony
# ---------------------------------------------------------------------------

#: State codes for the binary trait reconstruction.
_ABSENT, _PRESENT = 0, 1


def reconstruct_origins(tree, tip_states, gain_cost=1.0, loss_cost=1.0,
						root_state='auto', resolution='deltran'):
	"""Reconstruct where a binary trait arose, and how often independently.

	`tip_states` maps a tip label to 1 (carrier), 0 (non-carrier) or None (not
	determined - the tip is left free and costs nothing either way).

	Returns a dict with:

	* ``gains_min`` / ``gains_max`` and ``losses_min`` / ``losses_max`` - the
	  smallest and largest number of independent acquisitions (and reversions)
	  over **all** most-parsimonious reconstructions. This range, not a single
	  number, is the honest answer: parsimony very often does not identify a
	  unique history, and a claim that a mutation arose "17 times
	  independently" made from one arbitrarily-chosen reconstruction is a claim
	  about the tie-breaking rule. The two minima may be attained by different
	  reconstructions.
	* ``root_state`` - the state at the root. With ``root_state='auto'`` it is
	  whichever is cheaper, and that choice matters more than anything else
	  here: forcing the root to "absent" makes an *ancestral* residue look like
	  a mutation that arose independently in every sequence carrying it, which
	  on a polytomous tree is once per tip. An allele reconstructed as present
	  at the root is the wild type, and its interesting count is the losses.
	* ``node_state`` - one concrete most-parsimonious reconstruction, used to
	  cut the tree into clusters.
	* ``clusters`` - one entry per independent origin: the node where the trait
	  appears and every tip descended from it.

	Independent origins are what let a frequency increase be attributed to the
	mutation rather than to the lineage carrying it: one origin that happens to
	sit in an expanding clade is one observation, however many sequences it
	contains, while the same allele rising from several separate origins is
	repeated evidence.
	"""
	costs = np.array([[0.0, gain_cost], [loss_cost, 0.0]])
	big = float('inf')
	n = tree.n_nodes
	cost = np.zeros((n, 2))
	# Counted transitions: index 0 counts gains (absent -> present), index 1
	# losses. Each is minimised and maximised independently over the optimal
	# set, which is what makes them a range rather than one history's tally.
	counted = np.array([[[0, 0], [1, 0]], [[0, 1], [0, 0]]])
	low = np.zeros((n, 2, 2))
	high = np.zeros((n, 2, 2))

	for node in tree.postorder:
		if tree.is_tip[node]:
			state = tip_states.get(tree.name[node])
			if state is None:
				cost[node] = (0.0, 0.0)
			elif int(state) == _PRESENT:
				cost[node] = (big, 0.0)
			else:
				cost[node] = (0.0, big)
			continue
		for state in (_ABSENT, _PRESENT):
			total = 0.0
			low_counts = np.zeros(2)
			high_counts = np.zeros(2)
			for kid in tree.children[node]:
				options = [cost[kid][s] + costs[state][s] for s in (_ABSENT, _PRESENT)]
				best = min(options)
				if not math.isfinite(best):
					total = big
					break
				total += best
				candidates = [s for s in (_ABSENT, _PRESENT) if options[s] <= best + 1e-12]
				for which in (0, 1):
					low_counts[which] += min(low[kid][s][which] + counted[state][s][which]
											 for s in candidates)
					high_counts[which] += max(high[kid][s][which] + counted[state][s][which]
											  for s in candidates)
			cost[node][state] = total
			low[node][state] = low_counts
			high[node][state] = high_counts

	root = tree.root
	# Both root states can be equally parsimonious - "gained twice" and "present
	# at the root and lost twice" cost the same - and when they are, the origin
	# count genuinely ranges over both. Reporting the range from whichever state
	# argmin happened to return would hide exactly the ambiguity the range exists
	# to expose.
	best_cost = float(np.min(cost[root]))
	optimal = [state for state in (_ABSENT, _PRESENT)
			   if math.isfinite(cost[root][state]) and cost[root][state] <= best_cost + 1e-12]
	if root_state in (None, 'auto'):
		# On a tie, reconstruct with the state that implies fewer independent
		# origins: the conservative reading for any claim about independence.
		start = min(optimal, key=lambda state: low[root][state][0] + (1 if state == _PRESENT else 0))
	else:
		start = _PRESENT if root_state in (_PRESENT, 'present', True) else _ABSENT
		if not math.isfinite(cost[root][start]):
			start = int(np.argmin(cost[root]))
		optimal = [start]

	node_state = np.full(n, -1, dtype=np.int8)
	node_state[root] = start
	prefer_parent = resolution != 'acctran'
	for node in tree.preorder:
		if node_state[node] < 0:
			continue
		for kid in tree.children[node]:
			options = [cost[kid][s] + costs[node_state[node]][s] for s in (_ABSENT, _PRESENT)]
			best = min(options)
			candidates = [s for s in (_ABSENT, _PRESENT) if options[s] <= best + 1e-12]
			if len(candidates) == 1:
				node_state[kid] = candidates[0]
			elif prefer_parent:
				node_state[kid] = node_state[node] if node_state[node] in candidates else candidates[0]
			else:
				other = [s for s in candidates if s != node_state[node]]
				node_state[kid] = other[0] if other else candidates[0]

	clusters = []
	for node in tree.preorder:
		if node_state[node] != _PRESENT:
			continue
		par = tree.parent[node]
		if par >= 0 and node_state[par] == _PRESENT:
			continue
		members = []
		stack = [node]
		while stack:
			current = stack.pop()
			if tree.is_tip[current]:
				members.append(current)
				continue
			for kid in tree.children[current]:
				if node_state[kid] == _PRESENT:
					stack.append(kid)
		clusters.append({'origin_node': int(node), 'tips': members})

	realised = int(sum(1 for node in range(n)
					   if node_state[node] == _PRESENT
					   and (tree.parent[node] < 0 or node_state[tree.parent[node]] == _ABSENT)))
	# A trait present at the root was not gained anywhere in the tree, but it is
	# still one independent present lineage. Keeping the two counts apart is the
	# difference between "this arose once" and "this is the ancestral residue".
	def at_root(state):
		return 1 if state == _PRESENT else 0

	return {
		'gains_min': int(min(low[root][state][0] for state in optimal)),
		'gains_max': int(max(high[root][state][0] for state in optimal)),
		'origins_min': int(min(low[root][state][0] + at_root(state) for state in optimal)),
		'origins_max': int(max(high[root][state][0] + at_root(state) for state in optimal)),
		'losses_min': int(min(low[root][state][1] for state in optimal)),
		'losses_max': int(max(high[root][state][1] for state in optimal)),
		'origins_realised': realised,
		'parsimony_cost': float(cost[root][start]),
		'root_state': 'present' if start == _PRESENT else 'absent',
		'node_state': node_state,
		'clusters': clusters,
		'resolution': resolution,
	}


# ---------------------------------------------------------------------------
# Local branching index
# ---------------------------------------------------------------------------

def local_branching_index(tree, branch_lengths, tau):
	"""Local branching index (Neher, Russell & Shraiman 2014).

	The total branch length in a node's neighbourhood, discounted exponentially
	with distance on a scale tau. It is the standard tree-shape proxy for
	fitness and it is *not* a rate: it says which parts of the tree are
	branching densely right now, which is what a growing lineage looks like
	before there is enough time series to fit a growth rate to it.

	Computed by the usual two message-passing sweeps, so it costs one pass up
	and one pass down whatever the size of the tree.
	"""
	if tau is None or tau <= 0:
		raise ValueError('tau must be positive')
	lengths = np.nan_to_num(np.asarray(branch_lengths, dtype=float), nan=0.0)
	up = np.zeros(tree.n_nodes)
	down = np.zeros(tree.n_nodes)

	for node in tree.postorder:
		total = sum(up[kid] for kid in tree.children[node])
		decay = math.exp(-lengths[node] / tau)
		up[node] = total * decay + tau * (1.0 - decay)

	for node in tree.preorder:
		kids = tree.children[node]
		if not kids:
			continue
		sibling_total = sum(up[kid] for kid in kids)
		for kid in kids:
			message = down[node] + (sibling_total - up[kid])
			decay = math.exp(-lengths[kid] / tau)
			down[kid] = message * decay + tau * (1.0 - decay)

	lbi = np.array([down[node] + sum(up[kid] for kid in tree.children[node])
					for node in range(tree.n_nodes)])
	peak = float(np.max(lbi)) if tree.n_nodes else 0.0
	return lbi, (lbi / peak if peak > 0 else lbi)


# ---------------------------------------------------------------------------
# Lineages through time and the method-of-moments rates
# ---------------------------------------------------------------------------

def lineages_through_time(tree, node_times, nodes=None):
	"""Lineage count as a function of time. Returns ``(times, counts)``.

	Serial sampling is handled properly: a tip *removes* a lineage at its
	sampling date, where the textbook LTT plot for contemporaneous tips never
	has to. Ignoring that on virus data - where sampling is spread over the
	whole history - inflates the apparent lineage count towards the present and
	so inflates every rate read off the plot.
	"""
	if nodes is None:
		nodes = range(tree.n_nodes)
	events = []
	for node in nodes:
		time = node_times[node]
		if not np.isfinite(time):
			continue
		if tree.is_tip[node]:
			events.append((float(time), -1))
		else:
			kids = len(tree.children[node])
			if kids > 1:
				events.append((float(time), kids - 1))
	if not events:
		return np.zeros(0), np.zeros(0)
	events.sort()
	times, counts = [], []
	current = 1
	for time, change in events:
		current += change
		times.append(time)
		counts.append(max(current, 0))
	return np.asarray(times, dtype=float), np.asarray(counts, dtype=float)


def ltt_slope(times, counts, min_lineages=2):
	"""Net diversification rate from the slope of log(lineages) against time.

	The oldest and most direct phylogenetic growth estimate there is: under a
	birth-death process the reconstructed lineage count grows as exp((lambda -
	mu) t), so the slope of a semi-log LTT plot estimates the net rate. It says
	nothing about birth and death separately, which is why the
	method-of-moments estimators below are reported next to it.
	"""
	times = np.asarray(times, dtype=float)
	counts = np.asarray(counts, dtype=float)
	keep = np.isfinite(times) & (counts >= min_lineages)
	out = {'n_points': int(np.sum(keep)), 'slope_per_year': float('nan'),
		   'std_error': float('nan'), 'r_squared': float('nan'), 'p_value': float('nan')}
	if np.sum(keep) < 3:
		return out
	x = times[keep]
	y = np.log(counts[keep])
	if np.ptp(x) <= 0:
		return out
	slope, intercept = np.polyfit(x, y, 1)
	fitted = slope * x + intercept
	residual = y - fitted
	dof = len(x) - 2
	if dof <= 0:
		return out
	sigma2 = float(np.sum(residual ** 2)) / dof
	sxx = float(np.sum((x - np.mean(x)) ** 2))
	se = math.sqrt(sigma2 / sxx) if sxx > 0 else float('nan')
	total = float(np.sum((y - np.mean(y)) ** 2))
	out.update({'slope_per_year': float(slope), 'std_error': se,
				'r_squared': 1.0 - float(np.sum(residual ** 2)) / total if total > 0 else float('nan'),
				'p_value': gs.two_sided_t_p(slope / se, dof) if se and se > 0 else float('nan')})
	return out


def magallon_sanderson(n_tips, age, epsilon=0.0, crown=True):
	"""Method-of-moments net diversification rate (Magallon & Sanderson 2001).

	Uses only the size of a clade and its age, so it works on a lineage far too
	small or too poorly sampled to fit anything to. `epsilon` is the relative
	extinction rate mu/lambda, which cannot be estimated from a clade's size and
	age alone - it has to be assumed, and the usual practice, followed here, is
	to report the rate over a range of assumed values rather than pick one.
	"""
	n_tips = float(n_tips)
	age = float(age)
	if n_tips < 2 or age <= 0 or not 0.0 <= epsilon < 1.0:
		return float('nan')
	if not crown:
		return math.log(n_tips * (1.0 - epsilon) + epsilon) / age
	inner = ((n_tips / 2.0) * (1.0 - epsilon ** 2) + 2.0 * epsilon
			 + ((1.0 - epsilon) / 2.0) * math.sqrt(n_tips ** 2 * (1.0 - epsilon) ** 2
												   + 4.0 * n_tips * epsilon))
	if inner <= 0:
		return float('nan')
	return (math.log(inner) - math.log(2.0)) / age


def yule_birth_rate(n_tips, tree_length, crown=True):
	"""Maximum-likelihood pure-birth rate: (events) / (total tree length in time).

	The pure-birth special case has a closed form, so unlike the birth-death
	rate it needs no assumption about extinction and no numerical optimisation.
	It is an *upper* bound on the net diversification rate of the same clade.
	"""
	n_tips = float(n_tips)
	tree_length = float(tree_length)
	if tree_length <= 0 or n_tips < 3:
		return float('nan')
	events = n_tips - 2.0 if crown else n_tips - 1.0
	if events <= 0:
		return float('nan')
	return events / tree_length


def gamma_statistic(branching_times):
	"""Pybus & Harvey (2000) gamma: where in the tree the branching happened.

	Standard-normal under a constant-rate pure-birth process. Negative means
	branching is concentrated towards the root - the deceleration expected of
	any lineage approaching the limit of its niche or its susceptible pool -
	and positive means a recent burst. It is a shape statistic, so it is the
	check on whether a single exponential growth rate is a fair summary at all.

	`branching_times` are node ages measured backwards from the present.
	"""
	ages = np.sort(np.asarray([t for t in branching_times if np.isfinite(t)], dtype=float))[::-1]
	n = len(ages) + 1
	if n < 4:
		return float('nan')
	# Interval g_j is the time during which exactly j lineages existed.
	boundaries = np.concatenate([ages, [0.0]])
	intervals = boundaries[:-1] - boundaries[1:]
	j = np.arange(2, n + 1, dtype=float)
	weighted = j * intervals
	total = float(np.sum(weighted))
	if total <= 0:
		return float('nan')
	partial = np.cumsum(weighted)[:-1]
	numerator = float(np.sum(partial)) / (n - 2) - total / 2.0
	denominator = total * math.sqrt(1.0 / (12.0 * (n - 2)))
	if denominator <= 0:
		return float('nan')
	return numerator / denominator


# ---------------------------------------------------------------------------
# Coalescent estimators
# ---------------------------------------------------------------------------

def coalescent_intervals(tree, node_times, nodes=None):
	"""Intervals of constant lineage count, backwards from the most recent tip.

	Returns ``(intervals, coalescent_taus, info)`` where each interval is
	``(tau_start, tau_end, k)``. A polytomy is expanded into ``children - 1``
	simultaneous coalescences, which no coalescent model actually allows -
	`info['polytomy_fraction']` reports how much of the tree was handled that
	way so a reader can see when the coalescent numbers should not be trusted.
	UShER trees built from de-duplicated alignments are heavily polytomous and
	this is usually the reason to prefer the frequency-based estimates.
	"""
	if nodes is None:
		nodes = list(range(tree.n_nodes))
	events = []
	polytomies = 0
	internal = 0
	for node in nodes:
		time = node_times[node]
		if not np.isfinite(time):
			continue
		if tree.is_tip[node]:
			events.append((float(time), 'sample', 1))
		else:
			kids = len(tree.children[node])
			if kids < 2:
				continue
			internal += 1
			if kids > 2:
				polytomies += 1
			events.append((float(time), 'coalescence', kids - 1))
	info = {'n_samples': sum(c for _, kind, c in events if kind == 'sample'),
			'n_coalescences': sum(c for _, kind, c in events if kind == 'coalescence'),
			'polytomy_fraction': (polytomies / internal) if internal else float('nan')}
	if not events:
		return [], [], info

	most_recent = max(time for time, _, _ in events)
	# tau runs backwards from the most recent sample.
	timeline = sorted(((most_recent - time, kind, count) for time, kind, count in events),
					  key=lambda row: row[0])
	intervals = []
	coalescent_taus = []
	lineages = 0
	previous = 0.0
	for tau, kind, count in timeline:
		if tau > previous and lineages >= 2:
			intervals.append((previous, tau, lineages))
		previous = tau
		if kind == 'sample':
			lineages += count
		else:
			for _ in range(count):
				if lineages >= 2:
					coalescent_taus.append(tau)
					lineages -= 1
	return intervals, coalescent_taus, info


def skyline(intervals, coalescent_taus, group=1):
	"""Classic (and grouped) coalescent skyline estimate of Ne through time.

	Each coalescent event gets Ne = (the coalescent-weighted waiting time that
	preceded it), which is its maximum-likelihood estimate under a locally
	constant population size. With serial sampling the waiting time accumulates
	over every sub-interval since the previous coalescence, which is what makes
	this valid on heterochronous virus data.

	`group` averages over consecutive events - Strimmer & Pybus's generalised
	skyline - because the classic estimator's variance per point is enormous
	and an ungrouped plot is mostly noise.
	"""
	if not coalescent_taus:
		return []
	waiting = {}
	for start, end, k in intervals:
		weight = k * (k - 1) / 2.0
		if weight <= 0:
			continue
		waiting.setdefault(round(end, 12), 0.0)
		waiting[round(end, 12)] += weight * (end - start)

	ordered = sorted(set(round(t, 12) for t in coalescent_taus))
	points = []
	for tau in ordered:
		value = waiting.get(tau, 0.0)
		if value > 0:
			points.append({'tau': tau, 'ne': value})
	if group <= 1 or not points:
		return points
	grouped = []
	for start in range(0, len(points), group):
		chunk = points[start:start + group]
		grouped.append({'tau': float(np.mean([p['tau'] for p in chunk])),
						'ne': float(np.mean([p['ne'] for p in chunk])),
						'n_events': len(chunk)})
	return grouped


def _exponential_coalescent_profile(rate, intervals, coalescent_taus):
	"""Profile log-likelihood of the exponential-growth coalescent at `rate`.

	Ne(tau) = N0 exp(-rate.tau) with tau measured into the past, so `rate` is
	the growth rate forwards in time. N0 is profiled out analytically, which
	turns a two-parameter optimisation into a one-dimensional one and removes
	the usual ridge between N0 and the rate.
	"""
	m = len(coalescent_taus)
	if m == 0:
		return float('-inf')
	total = 0.0
	for start, end, k in intervals:
		weight = k * (k - 1) / 2.0
		if weight <= 0:
			continue
		if abs(rate) < 1e-12:
			total += weight * (end - start)
		else:
			try:
				total += weight * (math.exp(rate * end) - math.exp(rate * start)) / rate
			except OverflowError:
				return float('-inf')
	if total <= 0 or not math.isfinite(total):
		return float('-inf')
	linear = sum(rate * tau for tau in coalescent_taus)
	if not math.isfinite(linear):
		return float('-inf')
	return linear - m * math.log(total / m) - m


def coalescent_exponential_growth(intervals, coalescent_taus, bounds=(-10.0, 10.0),
								  grid=201, tolerance=1e-6):
	"""Maximum-likelihood exponential growth rate under the coalescent.

	This is the phylodynamic counterpart of the frequency-based estimate: it
	reads the growth rate off the *shape* of the genealogy - how the coalescent
	events bunch towards the present - rather than off how common the lineage
	is in the sample. The two are independent lines of evidence, and are worth
	comparing for exactly that reason.

	The 95% interval comes from the profile likelihood (a drop of 1.92 units),
	not from a curvature approximation, because the profile is markedly
	asymmetric for small clades.
	"""
	out = {'growth_rate_per_year': float('nan'), 'ci_low': float('nan'),
		   'ci_high': float('nan'), 'n0': float('nan'), 'loglik': float('nan'),
		   'n_coalescences': len(coalescent_taus), 'converged': False}
	if len(coalescent_taus) < 3 or not intervals:
		return out

	low, high = bounds
	candidates = np.linspace(low, high, grid)
	values = [_exponential_coalescent_profile(r, intervals, coalescent_taus) for r in candidates]
	best = int(np.argmax(values))
	if not math.isfinite(values[best]):
		return out

	# Golden-section refinement inside the best grid cell.
	left = candidates[max(best - 1, 0)]
	right = candidates[min(best + 1, grid - 1)]
	phi = (math.sqrt(5.0) - 1.0) / 2.0
	a, b = left, right
	c, d = b - phi * (b - a), a + phi * (b - a)
	fc, fd = (_exponential_coalescent_profile(c, intervals, coalescent_taus),
			  _exponential_coalescent_profile(d, intervals, coalescent_taus))
	for _ in range(200):
		if abs(b - a) < tolerance:
			break
		if fc > fd:
			b, d, fd = d, c, fc
			c = b - phi * (b - a)
			fc = _exponential_coalescent_profile(c, intervals, coalescent_taus)
		else:
			a, c, fc = c, d, fd
			d = a + phi * (b - a)
			fd = _exponential_coalescent_profile(d, intervals, coalescent_taus)
	rate = (a + b) / 2.0
	peak = _exponential_coalescent_profile(rate, intervals, coalescent_taus)

	def crosses(direction):
		step = (high - low) / grid
		position = rate
		for _ in range(2000):
			position += direction * step
			if not low - 1 <= position <= high + 1:
				return float('nan')
			if _exponential_coalescent_profile(position, intervals, coalescent_taus) < peak - 1.920729:
				return position
		return float('nan')

	# Recover N0 at the optimum from the same profiling identity.
	total = 0.0
	for start, end, k in intervals:
		weight = k * (k - 1) / 2.0
		if weight <= 0:
			continue
		if abs(rate) < 1e-12:
			total += weight * (end - start)
		else:
			total += weight * (math.exp(rate * end) - math.exp(rate * start)) / rate

	out.update({'growth_rate_per_year': float(rate), 'ci_low': crosses(-1),
				'ci_high': crosses(1), 'n0': total / len(coalescent_taus),
				'loglik': float(peak), 'converged': True})
	return out


# ---------------------------------------------------------------------------
# Scanning the tree for fast-growing clades
# ---------------------------------------------------------------------------

#: Binning every clade against every time bin is the one array in this module
#: whose size is (nodes x bins). A 130k-tip tree with monthly bins over a
#: century would be several gigabytes, so the caller is stopped before it
#: allocates rather than after.
MAX_BIN_CELLS = 60_000_000


def subtree_nodes(tree, node):
	"""Every node at or below `node`, root first."""
	collected = []
	stack = [node]
	while stack:
		current = stack.pop()
		collected.append(current)
		stack.extend(tree.children[current])
	return collected


def clade_bin_counts(tree, tip_bin, n_bins):
	"""``(nodes x bins)`` counts of cohort tips, accumulated up the tree.

	`tip_bin` is an array over nodes giving each tip's time bin, or -1 for a tip
	that is not in the cohort. Building this once and summing children into
	parents is what makes a whole-tree scan affordable: the alternative - taking
	the descendant set of each node in turn - is quadratic on a ladder-shaped
	tree, and virus trees are frequently ladder-shaped.
	"""
	if tree.n_nodes * n_bins > MAX_BIN_CELLS:
		raise MemoryError(
			'a %d-node tree against %d time bins needs %d cells; widen --bin-width '
			'or restrict the cohort' % (tree.n_nodes, n_bins, tree.n_nodes * n_bins))
	counts = np.zeros((tree.n_nodes, n_bins), dtype=np.int32)
	for node in tree.postorder:
		if tree.is_tip[node]:
			index = tip_bin[node]
			if index >= 0:
				counts[node, index] = 1
		else:
			for kid in tree.children[node]:
				counts[node] += counts[kid]
	return counts


def clade_time_bounds(tree, tip_time):
	"""Earliest and latest cohort sampling date under every node."""
	earliest = np.full(tree.n_nodes, np.inf)
	latest = np.full(tree.n_nodes, -np.inf)
	for node in tree.postorder:
		if tree.is_tip[node]:
			value = tip_time[node]
			if np.isfinite(value):
				earliest[node] = value
				latest[node] = value
		else:
			for kid in tree.children[node]:
				if earliest[kid] < earliest[node]:
					earliest[node] = earliest[kid]
				if latest[kid] > latest[node]:
					latest[node] = latest[kid]
	earliest[~np.isfinite(earliest)] = np.nan
	latest[~np.isfinite(latest)] = np.nan
	return earliest, latest


def scan_clade_growth(tree, bin_counts, bin_mid, min_clade_size=20,
					  max_clade_fraction=0.9, max_nodes=20000, quasi=True,
					  min_time_bins=3, min_bin_total=1):
	"""Logistic growth rate of every clade's frequency against the whole cohort.

	This is "lineage growth rates through the tree" in its most direct form: a
	clade is a lineage, its frequency among sampled sequences is what a fitness
	advantage moves, and the logistic slope is the advantage. Running it over
	every clade rather than a chosen few is what turns it from a test of one
	hypothesis into a scan, so the results are multiple-testing corrected by the
	caller.

	Clades spanning more than `max_clade_fraction` of the cohort are skipped:
	their "frequency" is most of the tree, so the fit is measuring the
	background against itself.
	"""
	totals = bin_counts[tree.root].astype(float)
	cohort_size = float(totals.sum())
	clade_sizes = bin_counts.sum(axis=1)
	# Bins holding almost nothing are not evidence about a frequency, and a long
	# tail of them - a handful of sequences per decade before surveillance
	# started - otherwise dominates the time axis every slope is fitted against.
	usable_bins = totals >= max(1, min_bin_total)
	if cohort_size <= 0 or int(np.sum(usable_bins)) < min_time_bins:
		return []

	times = np.asarray(bin_mid, dtype=float)[usable_bins]
	totals_used = totals[usable_bins]

	candidates = [node for node in range(tree.n_nodes)
				  if not tree.is_tip[node]
				  and node != tree.root
				  and clade_sizes[node] >= min_clade_size
				  and clade_sizes[node] <= max_clade_fraction * cohort_size]
	candidates.sort(key=lambda node: -clade_sizes[node])
	truncated = len(candidates) > max_nodes
	candidates = candidates[:max_nodes]

	results = []
	for node in candidates:
		carriers = bin_counts[node][usable_bins].astype(float)
		occupied = int(np.sum(carriers > 0))
		if occupied < min_time_bins:
			continue
		fit = gs.logistic_growth(times, carriers / totals_used, trials=totals_used,
								 quasi=quasi)
		if not np.isfinite(fit['growth_advantage_per_year']):
			continue
		fit.update({'node': int(node), 'node_label': tree.name[node],
					'clade_size': int(clade_sizes[node]),
					'clade_fraction': float(clade_sizes[node]) / cohort_size,
					'n_time_bins_occupied': occupied})
		results.append(fit)
	results.sort(key=lambda row: -row['growth_advantage_per_year'])
	if truncated and results:
		results[0]['scan_truncated'] = True
	return results


def expansion_score(tree, tip_time, window, reference_time=None):
	"""Tips sampled in the last `window` years over those in the window before.

	The crude, assumption-free companion to the model fits: a lineage that has
	doubled its share of recent sequencing is visible here with no time series,
	no clock and no model, which makes it the sanity check when the logistic
	fit says something surprising.
	"""
	finite = tip_time[np.isfinite(tip_time)]
	if finite.size == 0:
		return np.full(tree.n_nodes, np.nan)
	if reference_time is None:
		reference_time = float(np.max(finite))
	recent_cut = reference_time - window
	earlier_cut = reference_time - 2 * window

	recent = np.zeros(tree.n_nodes)
	earlier = np.zeros(tree.n_nodes)
	for node in tree.postorder:
		if tree.is_tip[node]:
			value = tip_time[node]
			if np.isfinite(value):
				if value > recent_cut:
					recent[node] = 1
				elif value > earlier_cut:
					earlier[node] = 1
		else:
			for kid in tree.children[node]:
				recent[node] += recent[kid]
				earlier[node] += earlier[kid]
	# Haldane correction: a lineage with no earlier sequences is the most
	# striking case there is, and must not come back as a division by zero.
	return (recent + 0.5) / (earlier + 0.5)


def mean_terminal_branch_length(tree, branch_lengths, nodes=None):
	"""Mean length of the terminal branches under a set of nodes.

	Short terminal branches relative to the rest of the tree are the signature
	of a lineage that has expanded recently: many tips have had little time to
	accumulate private changes since their common ancestor.
	"""
	if nodes is None:
		nodes = tree.tips
	values = [branch_lengths[node] for node in nodes
			  if tree.is_tip[node] and np.isfinite(branch_lengths[node])]
	return float(np.mean(values)) if values else float('nan')


def induced_subtree(tree, keep_tips):
	"""The tree restricted to `keep_tips`, with unifurcations suppressed.

	Returns ``(subtree, source_index)`` where ``source_index[i]`` is the node in
	the original tree that new node `i` came from. Branch lengths are summed
	along every suppressed path, so distances between surviving nodes are
	preserved.

	Restricting matters for the phylodynamic estimators, not just for speed: a
	lineages-through-time plot or a coalescent fit computed on a tree still
	containing tips that are not in the cohort is describing a different
	population from the one the frequencies were measured in.
	"""
	keep = np.zeros(tree.n_nodes, dtype=bool)
	for tip in keep_tips:
		keep[tip] = True
	kept_below = keep.astype(np.int64)
	for node in tree.postorder:
		if not tree.is_tip[node]:
			kept_below[node] = sum(kept_below[kid] for kid in tree.children[node])

	if kept_below[tree.root] == 0:
		raise ValueError('none of the requested tips is in this tree')

	# The new root is the deepest node that still has more than one kept child.
	root = tree.root
	while True:
		branching = [kid for kid in tree.children[root] if kept_below[kid] > 0]
		if len(branching) == 1 and not tree.is_tip[branching[0]]:
			root = branching[0]
			continue
		if len(branching) == 1 and tree.is_tip[branching[0]]:
			root = branching[0]
		break

	parent, children, name, length, source = [], [], [], [], []

	def add(old, new_parent, branch):
		index = len(parent)
		parent.append(new_parent if new_parent is not None else -1)
		children.append([])
		name.append(tree.name[old])
		length.append(branch)
		source.append(old)
		if new_parent is not None:
			children[new_parent].append(index)
		return index

	lengths = np.nan_to_num(tree.length, nan=0.0)
	new_root = add(root, None, 0.0)
	stack = [(root, new_root)]
	while stack:
		old, new = stack.pop()
		for kid in tree.children[old]:
			if kept_below[kid] == 0:
				continue
			# Walk down through any chain of single-child-of-interest nodes.
			current = kid
			distance = lengths[kid]
			while not tree.is_tip[current]:
				branching = [grandchild for grandchild in tree.children[current]
							 if kept_below[grandchild] > 0]
				if len(branching) != 1:
					break
				current = branching[0]
				distance += lengths[current]
			child_index = add(current, new, distance)
			if not tree.is_tip[current]:
				stack.append((current, child_index))

	return Tree(parent, children, name, length), np.asarray(source, dtype=np.int64)
