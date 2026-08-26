"""Edge cases in tree-free mode, clade/genotype assignment and tree handling.

Lens: clade_from_tree.assign_labels_from_tree, CladeAssignment, BuildTreeManifest,
PrepareUsherResume, ReplaceUsherTreeInDb, and the CreateSqliteDB tier stack
(direct ref-list label > tree neighbourhood > EPA-ng > BLAST top hit).

Every case below is anchored to something the shipped profiles actually produce.
The two reference lists wired into nextflow.config are the recurring trigger:

  * generic/rabv/ref_list_clades.txt  - 437 of its 497 rows carry genotype "NA"
    and subtype "NA", i.e. the majority of RABV references are unlabelled.
  * generic/hcv/ref_list_subtype_genotype.txt - 111 of 437 rows carry a real
    genotype but subtype "NA".

"NA" there is a missing-value sentinel, and clade_from_tree._clean_label maps it
(plus "n/a", "nan", "none", "null") to "" for tip names, genotypes and subtypes
alike - one sentinel vocabulary applied to three columns with different
semantics. That is the same shape as the influenza segment/"NA" collision.
"""

import csv
import io
import contextlib
import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

import BuildTreeManifest
import PrepareUsherResume
from clade_from_tree import assign_labels_from_tree
from CreateSqliteDB import CreateSqliteDB
from ReplaceUsherTreeInDb import replace_tree


REPO_ROOT = Path(__file__).resolve().parents[2]
HCV_DB = REPO_ROOT / "test_out" / "HCV_OM_test" / "HCV_OM_test.db"
TREEFREE_DB = REPO_ROOT / "test_out" / "basic_test_treefree" / "rabv-jul0425.db"
RABV_CLADE_REF_LIST = REPO_ROOT / "generic" / "rabv" / "ref_list_clades.txt"
HCV_CLADE_REF_LIST = REPO_ROOT / "generic" / "hcv" / "ref_list_subtype_genotype.txt"


REFS = {
	"R1a": {"genotype": "1", "subtype": "1a"},
	"R1b": {"genotype": "1", "subtype": "1a"},
	"R2a": {"genotype": "2", "subtype": "2b"},
	"R2b": {"genotype": "2", "subtype": "2b"},
}


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _write_ref_list(path: Path, rows):
	"""Reference list in the shipped 5-column shape:
	primary_accession / accession_type / segment / genotype / subtype."""
	with path.open("w", newline="", encoding="utf-8") as handle:
		writer = csv.writer(handle, delimiter="\t")
		for row in rows:
			writer.writerow(row)
	return path


def _db(tmp_path: Path, **kwargs) -> CreateSqliteDB:
	"""CreateSqliteDB carrying only the clade-relevant inputs.

	_add_reference_columns / _tree_based_reference_labels never touch the other
	input files, so placeholders are enough and no build has to be run. Keyword
	arguments only, so that new constructor parameters cannot silently shift a
	positional here.
	"""
	base = dict(
		meta_data="unused", features="unused", pad_aln="unused", gene_info=None,
		m49_countries="unused", m49_interm_region="unused", m49_regions="unused",
		m49_sub_regions="unused", proj_settings="unused", fasta_sequence_file="unused",
		insertions="unused", host_taxa_file="unused", base_dir=str(tmp_path),
		output_dir="SqliteDB", db_name="clade_lens", db_status="new db",
	)
	base.update(kwargs)
	return CreateSqliteDB(**base)


def _meta(accessions):
	return pd.DataFrame({"primary_accession": list(accessions)})


def _aln(pairs):
	"""sequence_alignment shape: which reference each accession was aligned to.

	alignment_name IS the BLAST top hit for a query, and is what tree-free mode
	falls back to for nearest_reference_genotype.
	"""
	if not pairs:
		return pd.DataFrame(columns=["primary_accession", "alignment_name"])
	return pd.DataFrame(pairs, columns=["primary_accession", "alignment_name"])


def _row(df, accession):
	hit = df[df["primary_accession"] == accession].iloc[0]
	return hit["nearest_reference_genotype"], hit["nearest_reference_subtype"]


# ==========================================================================
# clade_from_tree: post-order aggregation
# ==========================================================================

def test_two_genotype_tie_at_deciding_node_is_a_no_call():
	"""A query sitting between one genotype-1 and one genotype-2 reference must
	not be given either label.

	Real trigger: UShER is run with -C, which collapses zero-mutation internal
	nodes into polytomies, so identical-distance queries routinely land in a node
	holding an even mix of references. Breaking that tie alphabetically would
	hand every such query genotype "1" for HCV.
	"""
	result = assign_labels_from_tree("(Q1:0.1,R1a:0.1,R2a:0.1);", REFS)
	assert result["Q1"] == {"genotype": "", "subtype": ""}


def test_polytomy_majority_beats_the_purity_floor():
	"""A 2:1 polytomy (0.667 purity) clears DEFAULT_MIN_PURITY=0.6 and is called.

	The same -C collapse that creates ties also creates large legitimate
	polytomies; those must still resolve or a collapsed UShER tree would assign
	nothing at all.
	"""
	result = assign_labels_from_tree("(Q1:0.1,R1a:0.1,R1b:0.1,R2a:0.1);", REFS)
	assert result["Q1"]["genotype"] == "1"


def test_negative_and_scientific_notation_branch_lengths_are_tolerated():
	"""Assignment is topological, so it must survive the branch lengths real
	tools emit: IQ-TREE writes 1e-06 for near-identical sequences and can emit
	small negative lengths under some models.
	"""
	assert assign_labels_from_tree("((Q1:-0.001,R1a:1e-06):2E-3,R2a:0.5);", REFS)["Q1"]["genotype"] == "1"


def test_missing_semicolon_and_trailing_whitespace_still_parse():
	"""A tree copied out of a log, or truncated by a trailing-newline strip, can
	arrive without its terminating semicolon; Biopython accepts it and so must we.
	"""
	assert assign_labels_from_tree("  ((Q1:0.1,R1a:0.1):0.1,R2a:0.5)  \n", REFS)["Q1"]["genotype"] == "1"


def test_single_tip_and_reference_free_trees_return_blanks_not_crashes():
	"""A segment with one surviving sequence produces a one-tip tree, and an
	influenza ref_list carries no genotype column at all, so a whole tree can be
	reference-free. Both must degrade to "no call", never to an exception.
	"""
	assert assign_labels_from_tree("Q1;", REFS) == {"Q1": {"genotype": "", "subtype": ""}}
	assert assign_labels_from_tree("(Q1:0.1,Q2:0.1);", {}) == {
		"Q1": {"genotype": "", "subtype": ""},
		"Q2": {"genotype": "", "subtype": ""},
	}


def test_quoted_labels_with_commas_parentheses_and_colons_keep_their_full_name():
	"""GISAID-style strain names ("A/swine/Iowa/1/2020 (H1N1)") reach the tree as
	single-quoted newick labels containing the very characters newick uses as
	delimiters. If the label were split on them, the accession key would no
	longer match meta_data and the assignment would be dropped.
	"""
	result = assign_labels_from_tree("(('Q,1(a):x':0.1,R1a:0.1):0.1,R2a:0.5);", REFS)
	assert "Q,1(a):x" in result
	assert result["Q,1(a):x"]["genotype"] == "1"


def test_very_long_tip_label_is_not_truncated():
	"""Concatenated GenBank definition lines can exceed 250 characters; the key
	must survive intact or it will not join back to primary_accession."""
	long_name = "Q" + "x" * 400
	result = assign_labels_from_tree(f"(({long_name}:0.1,R1a:0.1):0.1,R2a:0.5);", REFS)
	assert long_name in result


@pytest.mark.xfail(
	reason="A tip whose label is literally 'NA' (or n/a, nan, none, null) is dropped: "
		   "_clean_label applies the missing-value sentinel vocabulary to tip NAMES as "
		   "well as to clade values, so the tip is neither counted as a reference nor "
		   "emitted as a query and silently receives no assignment.",
	strict=False,
)
def test_tip_named_na_is_still_assignable():
	"""Same shape as the influenza segment/"NA" collision: a legitimate label that
	looks like a null sentinel. Influenza tip labels can be gene/segment names
	("NA" = neuraminidase) once GISAID strain identifiers are merged in, and a
	tree tip called NA then vanishes from clade assignment without a warning.
	"""
	result = assign_labels_from_tree("((NA:0.1,R1a:0.1):0.1,R2a:0.5);", REFS)
	assert "NA" in result


@pytest.mark.xfail(
	reason="Duplicate tip labels are counted once per occurrence, so a reference that "
		   "appears twice in the alignment inflates its own share of the deciding "
		   "clade and converts a deliberate no-call into a confident call.",
	strict=False,
)
def test_duplicate_reference_tip_does_not_inflate_purity():
	"""Real trigger: --gisaid_dir merges bring the same accession in from GenBank
	and GISAID, and UsherPlacement.expand_identical_sequence_tree_outputs
	deliberately re-inserts identical sequences as extra polytomy members. Either
	route can put one reference into the tree twice.

	(Q1,R1a,R2a) is a 1:1 tie and correctly yields no call. Duplicating R1a must
	not turn that into "genotype 1" - the evidence is still one reference each.
	"""
	tie = assign_labels_from_tree("(Q1:0.1,R1a:0.1,R2a:0.1);", REFS)
	assert tie["Q1"]["genotype"] == ""
	duplicated = assign_labels_from_tree("(Q1:0.1,R1a:0.1,R1a:0.1,R2a:0.1);", REFS)
	assert duplicated["Q1"]["genotype"] == ""


@pytest.mark.xfail(
	reason="genotype and subtype are resolved by two independent climbs, so a query "
		   "can be handed a genotype from one clade and a subtype from a different, "
		   "contradictory clade (e.g. genotype 2 + subtype 1b).",
	strict=False,
)
def test_genotype_and_subtype_come_from_the_same_neighbourhood():
	"""The live trigger is the shipped HCV reference list: 111 of its 437 rows
	carry a real genotype but subtype "NA". When a query's nearest labelled
	neighbour is one of those, the genotype climb stops there while the subtype
	climb keeps going until it reaches a reference of a *different* genotype -
	and the pair written to meta_data is one no reference actually has.

	Here R_NEAR is genotype 2 / subtype NA and R_FAR is genotype 1 / subtype 1b.
	"""
	refs = {"R_NEAR": {"genotype": "2", "subtype": ""}, "R_FAR": {"genotype": "1", "subtype": "1b"}}
	result = assign_labels_from_tree("((Q1:0.1,R_NEAR:0.1):0.1,R_FAR:0.5);", refs)
	genotype, subtype = result["Q1"]["genotype"], result["Q1"]["subtype"]
	assert not (genotype and subtype) or subtype.startswith(genotype)


@pytest.mark.xfail(
	reason="A reference whose ref-list genotype AND subtype are both the 'NA' sentinel "
		   "is not recognised as a reference at all: it is emitted as a query and given "
		   "an invented clade, which CladeAssignment.process()/CreateSqliteDB then write "
		   "into the same column as curated labels with no provenance.",
	strict=False,
)
def test_unlabelled_reference_is_not_handed_an_invented_clade():
	"""generic/rabv/ref_list_clades.txt (used by two nextflow profiles) has 437 of
	497 rows at genotype=NA, subtype=NA. Under assign_labels_from_tree the
	majority of RABV *references* are therefore treated as queries and each one
	is assigned whatever its neighbours happen to be. Nothing downstream can tell
	that label apart from the ~60 curated ones.
	"""
	refs = {"R1a": {"genotype": "1", "subtype": "1a"}, "R_UNLABELLED": {"genotype": "", "subtype": ""}}
	result = assign_labels_from_tree("((R_UNLABELLED:0.1,R1a:0.1):0.1,Q1:0.5);", refs)
	assert "R_UNLABELLED" not in result


def test_malformed_newick_returns_empty_without_raising():
	"""IQ-TREE and UShER can leave a half-written tree behind when a run is killed.
	The call must not explode - but note it also cannot be distinguished from
	"tree had no usable references", which is what
	test_truncated_treefile_is_not_stored_or_is_warned_about covers.
	"""
	assert assign_labels_from_tree("((Q1:0.1,R1a:0.1)", REFS) == {}
	assert assign_labels_from_tree("", REFS) == {}
	assert assign_labels_from_tree(None, REFS) == {}


@pytest.mark.xfail(
	reason="Any string containing neither '(' nor ';' is treated as a file path, so a "
		   "one-tip newick written without its semicolon is opened as a filename, fails, "
		   "and silently yields no assignments.",
	strict=False,
)
def test_bare_single_tip_newick_is_not_mistaken_for_a_path():
	"""A segment reduced to a single sequence can be written as just its label.
	Being read as a path makes it indistinguishable from a missing tree.
	"""
	assert assign_labels_from_tree("Q1", REFS) == {"Q1": {"genotype": "", "subtype": ""}}


# ==========================================================================
# CreateSqliteDB: the tier stack (direct > tree > EPA-ng > BLAST top hit)
# ==========================================================================

def test_direct_reference_label_wins_over_every_inferred_source(tmp_path: Path):
	"""A reference's own curated ref_list label must never be replaced by a tree
	or BLAST inference, or a curated dataset would drift on every rebuild.
	"""
	ref_list = _write_ref_list(tmp_path / "ref.txt", [
		["R1a", "reference", "1", "1", "1a"],
		["R2a", "reference", "1", "2", "2b"],
	])
	tree = tmp_path / "usher.nh"
	tree.write_text("((R1a:0.1,R2a:0.1):0.1,Q1:0.5);\n", encoding="utf-8")
	db = _db(tmp_path, reference_tsv=str(ref_list), usher_tree=str(tree))
	out = db._add_reference_columns(_meta(["R1a", "R2a"]), _aln([("R1a", "R2a")]))
	assert _row(out, "R1a") == ("1", "1a")
	assert _row(out, "R2a") == ("2", "2b")


def test_tree_neighbourhood_wins_over_blast_top_hit(tmp_path: Path):
	"""The whole point of the tree tier: a query whose BLAST top hit is an
	off-clade but highly similar reference must still take its tree neighbours'
	genotype. A partial genome hitting a conserved region is the classic case.
	"""
	ref_list = _write_ref_list(tmp_path / "ref.txt", [
		["R1a", "reference", "1", "1", "1a"],
		["R2a", "reference", "1", "2", "2b"],
	])
	tree = tmp_path / "usher.nh"
	tree.write_text("((Q1:0.1,R1a:0.1):0.1,R2a:0.5);\n", encoding="utf-8")
	db = _db(tmp_path, reference_tsv=str(ref_list), usher_tree=str(tree))
	# BLAST says R2a (genotype 2); the tree says genotype 1.
	out = db._add_reference_columns(_meta(["Q1"]), _aln([("Q1", "R2a")]))
	assert _row(out, "Q1") == ("1", "1a")


def test_epa_ng_assignment_wins_over_blast_but_loses_to_the_tree(tmp_path: Path):
	"""Precedence must be strict: EPA-ng only exists to cover runs with no tree,
	so it must never displace a tree call, and must always displace BLAST.
	"""
	ref_list = _write_ref_list(tmp_path / "ref.txt", [
		["R1a", "reference", "1", "1", "1a"],
		["R2a", "reference", "1", "2", "2b"],
	])
	clade_tsv = tmp_path / "clade_assignments.tsv"
	clade_tsv.write_text("primary_accession\tgenotype\tsubtype\nQ1\t3\t3c\nQ2\t3\t3c\n", encoding="utf-8")

	# No tree: EPA-ng beats the BLAST hit R2a.
	db = _db(tmp_path, reference_tsv=str(ref_list), clade_assignments=str(clade_tsv))
	out = db._add_reference_columns(_meta(["Q1"]), _aln([("Q1", "R2a")]))
	assert _row(out, "Q1") == ("3", "3c")

	# With a tree: the tree beats EPA-ng.
	tree = tmp_path / "usher.nh"
	tree.write_text("((Q1:0.1,R1a:0.1):0.1,R2a:0.5);\n", encoding="utf-8")
	db_tree = _db(tmp_path, reference_tsv=str(ref_list), clade_assignments=str(clade_tsv), usher_tree=str(tree))
	out_tree = db_tree._add_reference_columns(_meta(["Q1"]), _aln([("Q1", "R2a")]))
	assert _row(out_tree, "Q1") == ("1", "1a")


def test_tree_free_mode_falls_back_to_the_blast_top_hit(tmp_path: Path):
	"""--tree_free skips phylogenetics and (per vgtk-init.nf) also skips the
	EPA-ng fallback, so alignment_name - the BLAST top hit - is the only source
	left. This pins that the fallback fires at all.
	"""
	ref_list = _write_ref_list(tmp_path / "ref.txt", [["R2a", "reference", "1", "2", "2b"]])
	db = _db(tmp_path, reference_tsv=str(ref_list))
	out = db._add_reference_columns(_meta(["Q1"]), _aln([("Q1", "R2a")]))
	assert _row(out, "Q1") == ("2", "2b")


def test_tree_free_fallback_from_an_unlabelled_reference_is_blank_not_wrong(tmp_path: Path):
	"""When a query's BLAST top hit is one of the many unlabelled references
	(genotype "NA" in generic/rabv/ref_list_clades.txt), the fallback must return
	nothing rather than inventing a label. Absent is recoverable; wrong is not.
	"""
	ref_list = _write_ref_list(tmp_path / "ref.txt", [
		["R_NA", "reference", "1", "NA", "NA"],
		["R2a", "reference", "1", "2", "2b"],
	])
	db = _db(tmp_path, reference_tsv=str(ref_list))
	out = db._add_reference_columns(_meta(["Q1"]), _aln([("Q1", "R_NA")]))
	assert _row(out, "Q1") == ("", "")


def test_query_absent_from_the_tree_keeps_the_blast_answer(tmp_path: Path):
	"""IQ-TREE backbones hold only MMseqs cluster representatives, so most queries
	are simply not in the tree. Those must fall through to BLAST rather than being
	blanked by the tree tier.
	"""
	ref_list = _write_ref_list(tmp_path / "ref.txt", [
		["R1a", "reference", "1", "1", "1a"],
		["R2a", "reference", "1", "2", "2b"],
	])
	tree = tmp_path / "iqtree.treefile"
	tree.write_text("((Q_IN_TREE:0.1,R1a:0.1):0.1,R2a:0.5);\n", encoding="utf-8")
	db = _db(tmp_path, reference_tsv=str(ref_list), iqtree_file=str(tree))
	out = db._add_reference_columns(_meta(["Q_ABSENT"]), _aln([("Q_ABSENT", "R2a")]))
	assert _row(out, "Q_ABSENT") == ("2", "2b")


def test_tree_tip_absent_from_meta_data_does_not_create_a_row(tmp_path: Path):
	"""UShER trees can retain tips for sequences that were later filtered out of
	the matrix. Those must not be resurrected as meta_data rows.
	"""
	ref_list = _write_ref_list(tmp_path / "ref.txt", [["R1a", "reference", "1", "1", "1a"]])
	tree = tmp_path / "usher.nh"
	tree.write_text("((GHOST:0.1,R1a:0.1):0.1,Q1:0.5);\n", encoding="utf-8")
	db = _db(tmp_path, reference_tsv=str(ref_list), usher_tree=str(tree))
	out = db._add_reference_columns(_meta(["Q1"]), _aln([]))
	assert list(out["primary_accession"]) == ["Q1"]


@pytest.mark.xfail(
	reason="_tree_based_reference_labels merges the genotype and subtype fields "
		   "independently across candidate trees, so a query can end up with a genotype "
		   "from the UShER tree and a subtype from the IQ-TREE tree even when the two "
		   "trees place it in different genotype clades.",
	strict=False,
)
def test_genotype_and_subtype_come_from_the_same_tree(tmp_path: Path):
	"""Segmented builds stage one UShER tree per segment plus an IQ-TREE backbone,
	so several candidate trees routinely contain the same accession. Here the
	UShER tree only knows a genotype-2 reference (subtype "NA") and the IQ-TREE
	backbone only knows a genotype-1 reference - the merged row (2, 1b) is
	supported by neither tree.
	"""
	ref_list = _write_ref_list(tmp_path / "ref.txt", [
		["R_GENO_ONLY", "reference", "1", "2", "NA"],
		["R_SUBTYPED", "reference", "1", "1", "1b"],
	])
	usher = tmp_path / "usher.nh"
	usher.write_text("(Q1:0.1,R_GENO_ONLY:0.1);\n", encoding="utf-8")
	iqtree = tmp_path / "iqtree.treefile"
	iqtree.write_text("(Q1:0.1,R_SUBTYPED:0.1);\n", encoding="utf-8")
	db = _db(tmp_path, reference_tsv=str(ref_list), usher_tree=str(usher), iqtree_file=str(iqtree))
	genotype, subtype = _row(db._add_reference_columns(_meta(["Q1"]), _aln([])), "Q1")
	assert not (genotype and subtype) or subtype.startswith(genotype)


# ==========================================================================
# BuildTreeManifest
# ==========================================================================

def test_manifest_with_no_trees_is_header_only(tmp_path: Path):
	"""vgtk-init.nf gates --tree_manifest on `wc -l > 1`, so a tree-free run must
	produce a manifest with a header and nothing else.
	"""
	rows = BuildTreeManifest.build_tree_manifest(tmp_path / "m.tsv", iqtree_dir=tmp_path / "absent", usher_dir=None)
	assert rows == []
	assert (tmp_path / "m.tsv").read_text(encoding="utf-8").strip() == "source\tname\tsegment_key\tpath"


def test_manifest_strips_the_staged_directory_prefixes(tmp_path: Path):
	"""Nextflow stages per-segment output dirs as IQTree_MMseqClusters_<key> /
	Usher_MMseqClusters_<key>; the manifest's segment_key must be the bare key so
	CreateSqliteDB can map it onto a segment number.
	"""
	seg = tmp_path / "iq" / "IQTree_MMseqClusters_4"
	seg.mkdir(parents=True)
	(seg / "seg4.treefile").write_text("(A:0.1,B:0.1);\n", encoding="utf-8")
	rows = BuildTreeManifest.build_tree_manifest(tmp_path / "m.tsv", iqtree_dir=tmp_path / "iq")
	assert [(r["source"], r["segment_key"]) for r in rows] == [("iqtree", "4")]


@pytest.mark.xfail(
	reason="A zero-byte .treefile is emitted into the manifest as if it were a finished "
		   "tree; CreateSqliteDB then reads it, gets '', and drops the row, so the "
		   "segment silently ends up with no tree in the DB and nothing is logged.",
	strict=False,
)
def test_empty_treefile_is_not_advertised_as_a_tree(tmp_path: Path):
	"""IQ-TREE creates its output files before the search finishes, so a killed or
	still-running job leaves a 0-byte .treefile behind. The manifest is what tells
	the DB builder a segment has a tree, so an empty file must not appear in it.
	"""
	seg = tmp_path / "iq" / "IQTree_MMseqClusters_2"
	seg.mkdir(parents=True)
	(seg / "seg2.treefile").write_text("", encoding="utf-8")
	rows = BuildTreeManifest.build_tree_manifest(tmp_path / "m.tsv", iqtree_dir=tmp_path / "iq")
	assert rows == []


@pytest.mark.xfail(
	reason="Two staged directories that reduce to the same segment_key (e.g. "
		   "'IQTree_MMseqClusters_PB1' and a bare 'PB1') both emit a row with the same "
		   "source/name/segment_key, so the trees table gains two rows claiming to be "
		   "the same segment's tree and update mode's DELETE-then-append duplicates them.",
	strict=False,
)
def test_manifest_segment_keys_are_unique(tmp_path: Path):
	"""Resume/update runs stage both the fresh per-segment directory and the
	carried-over one, which is exactly how two directories collapse onto one key.
	"""
	root = tmp_path / "iq"
	for name, newick in [("IQTree_MMseqClusters_PB1", "(A:0.1,B:0.1);"), ("PB1", "(C:0.1,D:0.1);")]:
		d = root / name
		d.mkdir(parents=True)
		(d / "t.treefile").write_text(newick + "\n", encoding="utf-8")
	rows = BuildTreeManifest.build_tree_manifest(tmp_path / "m.tsv", iqtree_dir=root)
	keys = [r["segment_key"] for r in rows]
	assert len(keys) == len(set(keys)), f"duplicate segment keys in manifest: {keys}"


@pytest.mark.xfail(
	reason="A truncated newick is copied verbatim into trees.newick with no warning: "
		   "_read_tree_file only checks for emptiness, assign_labels_from_tree swallows "
		   "the parse error internally so CreateSqliteDB's own except-and-warn guard "
		   "never fires, and the corrupt tree only surfaces if TEST_DB_VALIDATION runs.",
	strict=False,
)
def test_truncated_treefile_is_not_stored_or_is_warned_about(tmp_path: Path):
	"""A tree file that is written but not finished (killed IQ-TREE/UShER job, a
	full disk) parses as neither empty nor valid. Storing it means every consumer
	of the trees table - and ValidateDbTree itself, which calls Phylo.read
	unguarded - hits the error much later, in a run that reported success.
	"""
	truncated = tmp_path / "seg.treefile"
	truncated.write_text("((A:0.1,B:0.1)", encoding="utf-8")
	manifest = tmp_path / "tree_manifest.tsv"
	manifest.write_text(
		"source\tname\tsegment_key\tpath\n" f"iqtree\tiqtree_1\t1\t{truncated}\n", encoding="utf-8"
	)
	db = _db(tmp_path, tree_manifest=str(manifest))
	buffer = io.StringIO()
	with contextlib.redirect_stdout(buffer):
		stored = [db._read_tree_file(entry["path"]) for entry in db._load_tree_manifest(str(manifest))]
	assert stored == [None] or "warn" in buffer.getvalue().lower(), (
		f"truncated newick stored verbatim with no warning: {stored!r}"
	)


# ==========================================================================
# PrepareUsherResume
# ==========================================================================

def test_resume_prefers_the_uncondensed_tree(tmp_path: Path):
	"""UShER's final-tree.nh is the collapsed output; uncondensed-final-tree.nh is
	the one guaranteed to carry every placed sample. Resuming from the collapsed
	tree would lose tips, so the uncondensed file must win when both exist.
	"""
	chunk = tmp_path / "chunk_0001"
	chunk.mkdir()
	(chunk / "final-tree.nh").write_text("(R1:0.1);\n", encoding="utf-8")
	(chunk / "uncondensed-final-tree.nh").write_text("(R1:0.1,P1:0.2);\n", encoding="utf-8")
	_, tree_path = PrepareUsherResume.find_latest_chunk_with_tree(tmp_path)
	assert tree_path.name == "uncondensed-final-tree.nh"


def test_resume_tree_round_trip_keeps_awkward_tip_labels(tmp_path: Path):
	"""normalise_tree re-writes the tree through Biopython. Influenza/GISAID strain
	names carry '/', spaces and '|', and any mangling here makes
	subset_alignment fail to find the sequence - or worse, hand UShER a tip name
	that no longer matches the alignment.
	"""
	src = tmp_path / "in.nh"
	src.write_text("('A/swine/Iowa/1/2020 (H1N1)':0.1,gi|123|ref|NC_001542:0.2);\n", encoding="utf-8")
	out = tmp_path / "out.nh"
	PrepareUsherResume.normalise_tree(src, out)
	assert set(PrepareUsherResume.read_tree_terminals(out)) == {
		"A/swine/Iowa/1/2020 (H1N1)",
		"gi|123|ref|NC_001542",
	}


def test_empty_resume_tree_raises_rather_than_producing_a_stub(tmp_path: Path):
	"""An empty chunk tree must fail loudly - resuming from a stub would silently
	re-place every sample from chunk 0.
	"""
	src = tmp_path / "in.nh"
	src.write_text("\n", encoding="utf-8")
	with pytest.raises(ValueError):
		PrepareUsherResume.normalise_tree(src, tmp_path / "out.nh")


@pytest.mark.xfail(
	reason="_chunk_index regex-searches the whole path string, so a run directory whose "
		   "own path contains 'chunk_<n>' gives every child the parent's index; the "
		   "latest completed chunk is then picked by filesystem iteration order, not by "
		   "chunk number, and 'chunks/' is misread as a chunk directory too.",
	strict=False,
)
def test_chunk_index_reads_the_chunk_directory_not_its_parents(tmp_path: Path):
	"""Iterative UShER runs are commonly resumed by pointing --run_dir at a
	previous attempt's directory, and those get named after the chunk they died
	on ("chunk_1_rerun", "usher_chunk_3"). Resuming from an older chunk than the
	one that finished silently discards completed placement work.
	"""
	run_dir = tmp_path / "chunk_1_rerun"
	for name, newick in [("chunk_0001", "(R1:0.1,P1:0.2);"), ("chunk_0002", "(R1:0.1,P1:0.2,P2:0.3);")]:
		d = run_dir / name
		d.mkdir(parents=True)
		(d / "uncondensed-final-tree.nh").write_text(newick + "\n", encoding="utf-8")
	assert PrepareUsherResume._chunk_index(run_dir / "chunk_0002") == 2
	latest_dir, _ = PrepareUsherResume.find_latest_chunk_with_tree(run_dir)
	assert latest_dir.name == "chunk_0002"


@pytest.mark.xfail(
	reason="normalise_tree writes branch lengths with Biopython's 5-decimal default, so "
		   "UShER-scale lengths (mutations/site, ~1e-6) are all flattened to 0.00000 in "
		   "resume_tree.nwk.",
	strict=False,
)
def test_resume_tree_keeps_usher_scale_branch_lengths(tmp_path: Path):
	"""UShER trees for closely related viral genomes carry branch lengths around
	1e-6; a resume tree in which every one of them is exactly zero has lost all
	its distance information.
	"""
	src = tmp_path / "in.nh"
	src.write_text("(A:0.000001,B:0.0000015,C:0.1);\n", encoding="utf-8")
	out = tmp_path / "out.nh"
	PrepareUsherResume.normalise_tree(src, out)
	from io import StringIO
	from Bio import Phylo
	lengths = [t.branch_length for t in Phylo.read(StringIO(out.read_text(encoding="utf-8")), "newick").get_terminals()]
	assert all(length for length in lengths), f"branch lengths flattened to zero: {lengths}"


# ==========================================================================
# ReplaceUsherTreeInDb
# ==========================================================================

def _tree_db(path: Path, rows):
	conn = sqlite3.connect(str(path))
	conn.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, segment_key TEXT, newick TEXT)")
	conn.executemany("INSERT INTO trees(name, source, segment, segment_key, newick) VALUES (?,?,?,?,?)", rows)
	conn.commit()
	conn.close()
	return path


def test_replace_tree_refuses_an_ambiguous_or_missing_selector(tmp_path: Path):
	"""Replacing the wrong row, or silently replacing none, would leave the DB
	claiming a tree it does not have. Both must raise.
	"""
	db_path = _tree_db(tmp_path / "t.db", [
		("usher", "usher", "1", "", "(A:0.1);"),
		("usher", "usher", "1", "", "(B:0.1);"),
	])
	replacement = tmp_path / "new.nh"
	replacement.write_text("(A:0.1,B:0.1);\n", encoding="utf-8")
	with pytest.raises(ValueError):
		replace_tree(db_path, replacement, source="usher", name="usher", segment="1")
	with pytest.raises(ValueError):
		replace_tree(db_path, replacement, source="usher", name="nope", segment="1")


def test_replace_tree_refuses_an_empty_replacement(tmp_path: Path):
	"""An empty replacement file (a failed rerun) must not blank out a good tree."""
	db_path = _tree_db(tmp_path / "t.db", [("usher", "usher", "1", "", "(A:0.1);")])
	empty = tmp_path / "empty.nh"
	empty.write_text("\n", encoding="utf-8")
	with pytest.raises(ValueError):
		replace_tree(db_path, empty, source="usher", name="usher", segment="1")
	conn = sqlite3.connect(str(db_path))
	assert conn.execute("SELECT newick FROM trees").fetchone()[0] == "(A:0.1);"
	conn.close()


@pytest.mark.xfail(
	reason="replace_tree only checks that the replacement file is non-empty; a truncated "
		   "newick overwrites a good stored tree and the corruption is not detected until "
		   "something downstream tries to parse it.",
	strict=False,
)
def test_replace_tree_rejects_an_unparseable_replacement(tmp_path: Path):
	"""This script exists to hand-patch a tree into a finished DB, so it is run
	against real published databases - exactly the place where a truncated file
	must not be allowed to destroy the only copy of the tree.
	"""
	db_path = _tree_db(tmp_path / "t.db", [("usher", "usher", "1", "", "(A:0.1,B:0.1);")])
	truncated = tmp_path / "bad.nh"
	truncated.write_text("((A:0.1,B:0.1)", encoding="utf-8")
	with pytest.raises(ValueError):
		replace_tree(db_path, truncated, source="usher", name="usher", segment="1")


# ==========================================================================
# Real databases (read-only)
# ==========================================================================

@pytest.mark.skipif(not HCV_DB.exists(), reason=f"HCV test DB not present: {HCV_DB}")
def test_hcv_query_genotype_subtype_pairs_are_supported_by_a_reference():
	"""Guards the impossible-pair bug on a real build: every (genotype, subtype)
	pair assigned to a query must be a pair some reference actually carries. A
	query labelled genotype 2 / subtype 1b would mean the two fields were
	resolved from different clades.
	"""
	conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
	try:
		rows = conn.execute(
			"SELECT nearest_reference_genotype, nearest_reference_subtype, accession_type "
			"FROM meta_data "
			"WHERE COALESCE(nearest_reference_genotype,'') <> '' "
			"  AND COALESCE(nearest_reference_subtype,'') <> ''"
		).fetchall()
	finally:
		conn.close()
	reference_pairs = {(g, s) for g, s, kind in rows if kind in ("reference", "master")}
	assert reference_pairs, "no labelled reference rows to validate against"
	unsupported = sorted({(g, s) for g, s, kind in rows if kind not in ("reference", "master")} - reference_pairs)
	assert not unsupported, f"query clade pairs no reference carries: {unsupported}"


@pytest.mark.skipif(not HCV_DB.exists(), reason=f"HCV test DB not present: {HCV_DB}")
def test_hcv_stored_trees_all_parse():
	"""Every newick in the trees table must be parseable. This is the check that a
	truncated .treefile written into the DB would fail - ValidateDbTree only runs
	under test=1, so nothing else catches it on a production build.
	"""
	from io import StringIO
	from Bio import Phylo

	conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
	try:
		rows = conn.execute("SELECT name, newick FROM trees").fetchall()
	finally:
		conn.close()
	assert rows, "HCV DB has no trees"
	for name, newick in rows:
		assert newick and newick.strip(), f"tree {name!r} stored empty"
		tips = [t.name for t in Phylo.read(StringIO(newick), "newick").get_terminals() if t.name]
		assert tips, f"tree {name!r} has no named tips"


@pytest.mark.skipif(not TREEFREE_DB.exists(), reason=f"tree-free test DB not present: {TREEFREE_DB}")
def test_tree_free_db_has_no_trees_but_still_carries_the_clade_columns():
	"""--tree_free stores no trees at all, yet nearest_reference_genotype /
	nearest_reference_subtype must still exist so the schema does not change
	between run modes and downstream queries keep working.
	"""
	conn = sqlite3.connect(f"file:{TREEFREE_DB}?mode=ro", uri=True)
	try:
		tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
		if "trees" in tables:
			assert conn.execute("SELECT COUNT(*) FROM trees").fetchone()[0] == 0
		columns = {r[1] for r in conn.execute("PRAGMA table_info(meta_data)")}
	finally:
		conn.close()
	assert {"nearest_reference_genotype", "nearest_reference_subtype"} <= columns


@pytest.mark.skipif(not RABV_CLADE_REF_LIST.exists(), reason=f"missing {RABV_CLADE_REF_LIST}")
def test_shipped_rabv_clade_list_is_mostly_unlabelled():
	"""Documents the live trigger for the unlabelled-reference finding rather than
	asserting a fix: two nextflow profiles point ref_list at this file, and the
	overwhelming majority of its rows use "NA" for both clade columns, so most
	RABV references reach assign_labels_from_tree looking like queries.
	"""
	rows = [line.rstrip("\n").split("\t") for line in RABV_CLADE_REF_LIST.read_text(encoding="utf-8").splitlines() if line.strip()]
	sentinel = {"", "na", "n/a", "nan", "none", "null"}
	unlabelled = [r for r in rows if len(r) >= 5 and r[3].strip().lower() in sentinel and r[4].strip().lower() in sentinel]
	assert len(unlabelled) > len(rows) / 2, (
		"expected the shipped RABV clade list to be mostly NA/NA; "
		f"got {len(unlabelled)}/{len(rows)}"
	)


@pytest.mark.skipif(not HCV_CLADE_REF_LIST.exists(), reason=f"missing {HCV_CLADE_REF_LIST}")
def test_shipped_hcv_clade_list_has_genotype_only_references():
	"""Documents the live trigger for the split genotype/subtype climb: a quarter
	of the HCV references carry a genotype but subtype "NA", which is precisely
	the configuration that lets the subtype search escape into another genotype's
	clade.
	"""
	rows = [line.rstrip("\n").split("\t") for line in HCV_CLADE_REF_LIST.read_text(encoding="utf-8").splitlines() if line.strip()]
	sentinel = {"", "na", "n/a", "nan", "none", "null"}
	genotype_only = [r for r in rows if len(r) >= 5 and r[3].strip().lower() not in sentinel and r[4].strip().lower() in sentinel]
	assert genotype_only, "expected genotype-only (subtype NA) references in the shipped HCV list"
