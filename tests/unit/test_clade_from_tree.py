from pathlib import Path

from clade_from_tree import assign_labels_from_tree


REFS = {
	"R1a": {"genotype": "1", "subtype": "a"},
	"R1b": {"genotype": "1", "subtype": "a"},
	"R2a": {"genotype": "2", "subtype": "b"},
	"R2b": {"genotype": "2", "subtype": "b"},
}


def test_query_inherits_clade_of_nearest_reference_neighbourhood():
	# Q1 sits inside the genotype-2 clade, Q2 inside genotype-1.
	newick = "((R2a:0.1,(Q1:0.1,R2b:0.1):0.1):0.2,(R1a:0.1,(Q2:0.1,R1b:0.1):0.1):0.2);"
	result = assign_labels_from_tree(newick, REFS)
	assert result["Q1"] == {"genotype": "2", "subtype": "b"}
	assert result["Q2"] == {"genotype": "1", "subtype": "a"}
	# References are not re-reported; the caller already knows them.
	assert "R1a" not in result


def test_query_uses_majority_when_neighbourhood_is_mixed():
	# Q sits with two genotype-1 refs and one genotype-2 ref -> majority is 1.
	newick = "(((R1a:0.1,R1b:0.1):0.1,Q:0.1):0.1,R2a:0.2);"
	result = assign_labels_from_tree(newick, REFS)
	assert result["Q"]["genotype"] == "1"


def test_subtype_resolved_independently_from_genotype():
	# The nearest ref has a genotype but no subtype; the next one supplies subtype.
	refs = {
		"RG": {"genotype": "3", "subtype": ""},
		"RS": {"genotype": "3", "subtype": "c"},
	}
	newick = "((Q:0.1,RG:0.1):0.1,RS:0.2);"
	result = assign_labels_from_tree(newick, refs)
	assert result["Q"]["genotype"] == "3"
	assert result["Q"]["subtype"] == "c"


def test_returns_blank_when_tree_has_no_labelled_reference():
	newick = "(Q1:0.1,Q2:0.1);"
	result = assign_labels_from_tree(newick, REFS)
	assert result["Q1"] == {"genotype": "", "subtype": ""}
	assert result["Q2"] == {"genotype": "", "subtype": ""}


def test_accepts_a_path_to_a_newick_file(tmp_path: Path):
	tree_path = tmp_path / "tree.nwk"
	tree_path.write_text("((R2a:0.1,Q1:0.1):0.2,R1a:0.2);\n", encoding="utf-8")
	result = assign_labels_from_tree(str(tree_path), REFS)
	assert result["Q1"]["genotype"] == "2"


def test_empty_and_malformed_inputs_return_empty_dict():
	assert assign_labels_from_tree("", REFS) == {}
	assert assign_labels_from_tree(None, REFS) == {}
	assert assign_labels_from_tree("not a tree at all", REFS) == {}


def test_blank_when_only_reference_neighbourhood_spans_most_of_the_tree():
	"""Regression: a query in a reference-free region used to climb to a near-root
	node and inherit whichever genotype had the most references overall. That is
	noise, not phylogeny, so it must come back blank and let the caller fall back."""
	# Q hangs off a large reference-free cluster; the only clade containing any
	# reference is essentially the whole tree.
	queries = ",".join(f"F{i}:0.01" for i in range(60))
	newick = f"((Q:0.01,({queries}):0.01):0.01,(R1a:0.1,(R2a:0.1,R2b:0.1):0.1):0.1);"
	result = assign_labels_from_tree(newick, REFS)
	assert result["Q"] == {"genotype": "", "subtype": ""}


def test_blank_when_deciding_clade_is_too_impure():
	# Deciding clade holds one ref of each genotype -> 50% purity, below the
	# 0.6 default, so no call is made.
	newick = "((Q:0.1,(R1a:0.1,R2a:0.1):0.1):0.1,R1b:0.1);"
	result = assign_labels_from_tree(newick, REFS, max_clade_fraction=1.0)
	assert result["Q"]["genotype"] == ""


def test_purity_and_locality_guards_are_tunable():
	newick = "((Q:0.1,(R1a:0.1,R2a:0.1):0.1):0.1,R1b:0.1);"
	# Relaxing purity lets the (alphabetically stable) plurality through.
	result = assign_labels_from_tree(newick, REFS, min_purity=0.5, max_clade_fraction=1.0)
	assert result["Q"]["genotype"] == "1"


def test_local_pure_neighbourhood_still_assigns():
	# A tight, unambiguous neighbourhood is unaffected by the guards.
	newick = "(((Q:0.01,R2a:0.01):0.01,R2b:0.01):0.5,(R1a:0.1,R1b:0.1):0.5);"
	result = assign_labels_from_tree(newick, REFS)
	assert result["Q"] == {"genotype": "2", "subtype": "b"}


def test_reference_labels_may_be_empty():
	assert assign_labels_from_tree("(Q1:0.1,Q2:0.1);", {}) == {
		"Q1": {"genotype": "", "subtype": ""},
		"Q2": {"genotype": "", "subtype": ""},
	}
