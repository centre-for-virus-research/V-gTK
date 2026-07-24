#!/usr/bin/env python3
"""Assign genotype/subtype (a.k.a. major/minor clade) to query tips using the
phylogenetic tree instead of a single best BLAST hit.

The pipeline already builds a full-resolution tree (IQ-TREE backbone + UShER
placement, or a resume/update tree) whose tips include both the queries and the
labelled reference sequences. Deriving a query's genotype from the clade its
tree neighbours belong to is far more robust than the best-BLAST-hit heuristic:
a tip cannot disagree with all of its own neighbours, which is exactly the
artefact best-identity assignment produces (a partial genome or an off-clade but
highly-similar reference can win the top hit).

The core routine walks, for every query tip, up towards the root and returns the
consensus label of the nearest ancestor that encloses at least one labelled
reference tip. Genotype and subtype are resolved independently because some
references carry a genotype but no subtype.
"""

from collections import Counter
from io import StringIO

from Bio import Phylo


def _clean_label(value):
	text = "" if value is None else str(value).strip()
	if text.lower() in {"", "na", "n/a", "nan", "none", "null"}:
		return ""
	return text


def _consensus(counter):
	"""Most common label in a Counter, ties broken alphabetically for a stable
	result regardless of tree traversal / dict ordering."""
	if not counter:
		return ""
	return sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


DEFAULT_MIN_PURITY = 0.6
DEFAULT_MAX_CLADE_FRACTION = 0.5
# Below this many tips, "fraction of the whole tree" is not a meaningful locality
# signal (in a 4-tip tree any useful clade spans most of it), so the fraction
# guard only starts biting once a tree is big enough for it to mean something.
MIN_CLADE_TIP_ALLOWANCE = 25


def assign_labels_from_tree(newick, reference_labels,
							min_purity=DEFAULT_MIN_PURITY,
							max_clade_fraction=DEFAULT_MAX_CLADE_FRACTION):
	"""Assign labels to the query tips of a tree from their nearest labelled refs.

	A label is only accepted when the enclosing clade is genuinely informative
	about the query's neighbourhood. Without that guard, a query sitting in a
	reference-free region of the tree climbs until it hits a near-root node
	spanning most of the tree, and the "consensus" degenerates into "whichever
	genotype happens to have the most reference sequences overall" — which is
	noise, not phylogeny. Two guards prevent that:

	  * max_clade_fraction: the deciding clade must not span more than this
	    fraction of the tree's tips.
	  * min_purity: the winning label must account for at least this share of
	    the reference tips inside the deciding clade.

	When either guard fails the label is returned blank so the caller can fall
	back to whatever it had (e.g. the BLAST-derived assignment).

	Args:
		newick: newick string, or a path to a newick file.
		reference_labels: {accession: {"genotype": str, "subtype": str}} for the
			reference tips whose clade is known. Missing/blank values are ignored.
		min_purity: minimum winning-label share of references in the deciding clade.
		max_clade_fraction: maximum share of the tree the deciding clade may span.

	Returns:
		{accession: {"genotype": str, "subtype": str}} for every tip present in
		the tree that is NOT a reference (i.e. the queries). References are left
		out — the caller already knows their labels. Values are "" when no
		sufficiently local/pure reference neighbourhood supports a call.
	"""
	if not newick:
		return {}
	text = str(newick).strip()
	if not text:
		return {}
	# Accept either a newick string or a path to a newick file.
	if "(" not in text and ";" not in text:
		try:
			with open(text, "r", encoding="utf-8") as handle:
				text = handle.read().strip()
		except (OSError, ValueError):
			return {}
	if not text:
		return {}

	ref_geno = {}
	ref_sub = {}
	for acc, labels in (reference_labels or {}).items():
		acc = _clean_label(acc)
		if not acc:
			continue
		geno = _clean_label((labels or {}).get("genotype"))
		sub = _clean_label((labels or {}).get("subtype"))
		if geno:
			ref_geno[acc] = geno
		if sub:
			ref_sub[acc] = sub

	try:
		tree = Phylo.read(StringIO(text), "newick")
	except Exception:
		return {}

	# Parent pointers so each query tip can climb toward the root.
	parents = {}
	for clade in tree.find_clades(order="level"):
		for child in clade.clades:
			parents[id(child)] = clade

	# Post-order aggregation: each node holds a Counter of the reference genotype
	# and subtype labels among its descendant tips. Children precede parents in
	# post-order, so a parent can merge counters its children already produced.
	geno_counts = {}
	sub_counts = {}
	tip_counts = {}
	query_tips = []
	total_tips = 0
	for clade in tree.find_clades(order="postorder"):
		g_counter = Counter()
		s_counter = Counter()
		if clade.is_terminal():
			tips_here = 1
			total_tips += 1
			name = _clean_label(clade.name)
			if name and (name in ref_geno or name in ref_sub):
				if name in ref_geno:
					g_counter[ref_geno[name]] += 1
				if name in ref_sub:
					s_counter[ref_sub[name]] += 1
			elif name:
				query_tips.append((name, clade))
		else:
			tips_here = 0
			for child in clade.clades:
				g_counter.update(geno_counts.get(id(child), Counter()))
				s_counter.update(sub_counts.get(id(child), Counter()))
				tips_here += tip_counts.get(id(child), 0)
		geno_counts[id(clade)] = g_counter
		sub_counts[id(clade)] = s_counter
		tip_counts[id(clade)] = tips_here

	tip_cap = max(MIN_CLADE_TIP_ALLOWANCE, total_tips * max_clade_fraction) if total_tips else 0

	def _nearest(clade, counts):
		"""Nearest enclosing clade holding references, subject to the locality and
		purity guards. Climbing past a failing clade only makes both worse, so a
		failure ends the search with a blank label."""
		node = clade
		while node is not None:
			counter = counts.get(id(node))
			if counter:
				if tip_cap and tip_counts.get(id(node), 0) > tip_cap:
					return ""
				total = sum(counter.values())
				winner = _consensus(counter)
				if total and (counter[winner] / total) < min_purity:
					return ""
				return winner
			node = parents.get(id(node))
		return ""

	assignments = {}
	for name, clade in query_tips:
		# First query tip with a given name wins (duplicate tip names are rare
		# but harmless); keep it deterministic.
		if name in assignments:
			continue
		assignments[name] = {
			"genotype": _nearest(clade, geno_counts),
			"subtype": _nearest(clade, sub_counts),
		}
	return assignments
