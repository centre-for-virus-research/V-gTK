"""Tests for CladeAssignment: the tree-based method (preferred) and the EPA-ng
fallback (used only when no tree exists). Both paths are covered end to end,
including how reference labels are sourced, how gappa output is parsed, and how
the external commands are constructed."""

import csv
import subprocess
from pathlib import Path

import pytest

from CladeAssignment import CladeAssignment


# --------------------------------------------------------------------------- helpers


def write_gb_matrix(path, rows, fieldnames):
	with path.open("w", newline="", encoding="utf-8") as handle:
		writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
		writer.writeheader()
		for row in rows:
			writer.writerow(row)


def read_gb_matrix(path):
	with path.open(newline="", encoding="utf-8") as handle:
		return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv(path):
	with Path(path).open(newline="", encoding="utf-8") as handle:
		return list(csv.DictReader(handle, delimiter="\t"))


def make_taxon_files(tmp_path, major_rows, minor_rows):
	major = tmp_path / "major_clades.tsv"
	minor = tmp_path / "minor_clades.tsv"
	major.write_text("".join(f"{a}\t{c}\n" for a, c in major_rows), encoding="utf-8")
	minor.write_text("".join(f"{a}\t{c}\n" for a, c in minor_rows), encoding="utf-8")
	return major, minor


def make_reference_tsv(tmp_path, rows):
	"""HCV-style reference list: accession, status, segment, genotype, subtype."""
	path = tmp_path / "ref_list.tsv"
	path.write_text(
		"primary_accession\tstatus\tsegment\tgenotype\tsubtype\n"
		+ "".join(f"{a}\treference\t1\t{g}\t{s}\n" for a, g, s in rows),
		encoding="utf-8",
	)
	return path


DEFAULT_MATRIX_FIELDS = ["primary_accession", "gi_number", "accession_type", "exclusion_status"]


def make_matrix(tmp_path, entries, fieldnames=None):
	"""entries: list of (accession, accession_type)."""
	path = tmp_path / "gB_matrix.tsv"
	write_gb_matrix(
		path,
		[
			{"primary_accession": a, "gi_number": a, "accession_type": t, "exclusion_status": "0"}
			for a, t in entries
		],
		fieldnames or DEFAULT_MATRIX_FIELDS,
	)
	return path


def processor(tmp_path, gb_matrix=None, major=None, minor=None, **kwargs):
	return CladeAssignment(
		major_clade=str(major) if major else None,
		minor_clade=str(minor) if minor else None,
		padded_alignment=str(tmp_path / "padded.fasta"),
		base_dir=str(tmp_path),
		output_dir="CladeAssignment",
		gb_matrix=str(gb_matrix) if gb_matrix else None,
		threads="2",
		iqtree_model="GTR+G",
		**kwargs,
	)


# =========================================================================
# Reference label sourcing (taxon files, ref_list, both, neither)
# =========================================================================


def test_labels_from_taxon_files(tmp_path: Path):
	major, minor = make_taxon_files(tmp_path, [("R1", "1"), ("R2", "2")], [("R1", "a")])
	proc = processor(tmp_path, major=major, minor=minor)
	assert proc._load_reference_labels() == ({"R1": "1", "R2": "2"}, {"R1": "a"})


def test_labels_from_reference_tsv_when_no_taxon_files(tmp_path: Path):
	"""HCV-style datasets have no taxon files; genotype/subtype come from ref_list."""
	ref_tsv = make_reference_tsv(tmp_path, [("R1", "1", "a"), ("R2", "2", "b")])
	proc = processor(tmp_path, reference_tsv=str(ref_tsv))
	major, minor = proc._load_reference_labels()
	assert major == {"R1": "1", "R2": "2"}
	assert minor == {"R1": "a", "R2": "b"}


def test_taxon_files_take_precedence_over_reference_tsv(tmp_path: Path):
	ref_tsv = make_reference_tsv(tmp_path, [("R1", "1", "a"), ("R2", "2", "b")])
	# Curated taxon file overrides R1's genotype and adds R3.
	major, minor = make_taxon_files(tmp_path, [("R1", "9"), ("R3", "3")], [("R1", "z")])
	proc = processor(tmp_path, major=major, minor=minor, reference_tsv=str(ref_tsv))
	got_major, got_minor = proc._load_reference_labels()
	assert got_major == {"R1": "9", "R2": "2", "R3": "3"}
	assert got_minor == {"R1": "z", "R2": "b"}


def test_reference_tsv_na_subtype_is_dropped(tmp_path: Path):
	ref_tsv = make_reference_tsv(tmp_path, [("R1", "1", "NA")])
	proc = processor(tmp_path, reference_tsv=str(ref_tsv))
	major, minor = proc._load_reference_labels()
	assert major == {"R1": "1"}
	assert minor == {}


def test_process_raises_when_no_reference_labels_available(tmp_path: Path):
	matrix = make_matrix(tmp_path, [("Q1", "query")])
	proc = processor(tmp_path, gb_matrix=matrix)
	with pytest.raises(ValueError, match="No reference clade labels available"):
		proc.process()


def test_read_taxon_file_tolerates_header_and_blank_rows(tmp_path: Path):
	taxon = tmp_path / "taxon.tsv"
	taxon.write_text(
		"accession\tclade\n"   # header, skipped
		"REF1\tArctic\n"
		"REF2\t\n"              # blank clade, skipped
		"REF3\tAsian\n",
		encoding="utf-8",
	)
	assert CladeAssignment._read_taxon_file(str(taxon)) == {"REF1": "Arctic", "REF3": "Asian"}


def test_read_taxon_file_missing_returns_empty():
	assert CladeAssignment._read_taxon_file("/nonexistent/taxon.tsv") == {}


# =========================================================================
# TREE METHOD
# =========================================================================


def test_tree_method_assigns_queries_from_nearest_reference(tmp_path: Path, monkeypatch):
	major, minor = make_taxon_files(tmp_path, [("REF1", "1"), ("REF2", "2")], [("REF1", "a"), ("REF2", "b")])
	matrix = make_matrix(tmp_path, [("REF1", "reference"), ("REF2", "reference"), ("Q1", "query"), ("Q2", "query")])
	tree = tmp_path / "final-tree.nh"
	tree.write_text("((REF1:0.1,Q1:0.1):0.2,(REF2:0.1,Q2:0.1):0.2);\n", encoding="utf-8")

	proc = processor(tmp_path, gb_matrix=matrix, major=major, minor=minor, usher_tree=str(tree))
	monkeypatch.setattr(proc, "assign_with_epa_ng", lambda *a, **k: pytest.fail("EPA-ng must not run with a tree"))

	proc.process()

	rows = {r["primary_accession"]: r for r in read_gb_matrix(matrix)}
	assert (rows["Q1"]["major_clade"], rows["Q1"]["minor_clade"]) == ("1", "a")
	assert (rows["Q2"]["major_clade"], rows["Q2"]["minor_clade"]) == ("2", "b")
	assert (rows["REF1"]["major_clade"], rows["REF1"]["minor_clade"]) == ("1", "a")


def test_tree_method_works_from_reference_tsv_labels(tmp_path: Path):
	"""The tree path must work for ref_list-only datasets too."""
	ref_tsv = make_reference_tsv(tmp_path, [("REF1", "1", "a"), ("REF2", "2", "b")])
	matrix = make_matrix(tmp_path, [("REF1", "reference"), ("REF2", "reference"), ("Q1", "query")])
	tree = tmp_path / "final-tree.nh"
	tree.write_text("((REF1:0.1,Q1:0.1):0.2,REF2:0.2);\n", encoding="utf-8")

	proc = processor(tmp_path, gb_matrix=matrix, reference_tsv=str(ref_tsv), usher_tree=str(tree))
	proc.process()

	rows = {r["primary_accession"]: r for r in read_gb_matrix(matrix)}
	assert (rows["Q1"]["major_clade"], rows["Q1"]["minor_clade"]) == ("1", "a")


def test_usher_tree_preferred_over_iqtree(tmp_path: Path):
	major, minor = make_taxon_files(tmp_path, [("REF1", "1"), ("REF2", "2")], [("REF1", "a"), ("REF2", "b")])
	usher = tmp_path / "final-tree.nh"
	iq = tmp_path / "x.treefile"
	usher.write_text("((REF1:0.1,Q1:0.1):0.2,REF2:0.2);\n", encoding="utf-8")
	iq.write_text("((REF2:0.1,Q1:0.1):0.2,REF1:0.2);\n", encoding="utf-8")

	proc = processor(tmp_path, major=major, minor=minor, usher_tree=str(usher), iqtree_tree=str(iq))
	path, source = proc._select_tree()
	assert source == "usher" and path == str(usher)

	# Falls back to IQ-TREE when no UShER tree is present.
	proc2 = processor(tmp_path, major=major, minor=minor, iqtree_tree=str(iq))
	path2, source2 = proc2._select_tree()
	assert source2 == "iqtree" and path2 == str(iq)


def test_select_tree_ignores_sentinels_and_missing_files(tmp_path: Path):
	proc = processor(tmp_path, usher_tree="null", iqtree_tree="UNSET")
	assert proc._select_tree() == (None, None)
	proc2 = processor(tmp_path, usher_tree=str(tmp_path / "does_not_exist.nh"))
	assert proc2._select_tree() == (None, None)


def test_assign_from_tree_returns_blank_for_unsupported_neighbourhood(tmp_path: Path):
	major, minor = make_taxon_files(tmp_path, [("REF1", "1")], [("REF1", "a")])
	tree = tmp_path / "t.nh"
	tree.write_text("(Q1:0.1,Q2:0.1);\n", encoding="utf-8")  # no references at all
	proc = processor(tmp_path, major=major, minor=minor)
	got_major, got_minor = proc.assign_from_tree(str(tree), {"REF1": "1"}, {"REF1": "a"})
	assert got_major == {"Q1": "", "Q2": ""}
	assert got_minor == {"Q1": "", "Q2": ""}


# =========================================================================
# EPA-ng FALLBACK
# =========================================================================


def fake_epa_environment(proc, tmp_path, monkeypatch, major_rows, minor_rows):
	"""Stub iqtree/epa-ng/gappa. Records commands and writes the gappa per_query
	outputs the real tools would produce."""
	calls = []
	out_dir = tmp_path / "CladeAssignment"

	monkeypatch.setattr("shutil.which", lambda name: f"/env/bin/{name}" if name == "iqtree3" else None)

	def fake_run(cmd, check=True, **kwargs):
		calls.append(list(cmd))
		tool = Path(cmd[0]).name
		if tool == "iqtree3":
			prefix = cmd[cmd.index("-pre") + 1]
			Path(prefix + ".treefile").write_text("(R1:0.1,R2:0.1);\n", encoding="utf-8")
		elif tool == "epa-ng":
			outdir = Path(cmd[cmd.index("--outdir") + 1])
			outdir.mkdir(parents=True, exist_ok=True)
			(outdir / "epa_result.jplace").write_text("{}", encoding="utf-8")
		elif tool == "gappa":
			gappa_out = Path(cmd[cmd.index("--out-dir") + 1])
			gappa_out.mkdir(parents=True, exist_ok=True)
			rows = major_rows if gappa_out.name == "gappa_major_clades" else minor_rows
			with (gappa_out / "per_query.tsv").open("w", encoding="utf-8") as handle:
				handle.write("name\tLWR\tfract\taLWR\tafract\ttaxopath\n")
				for name, alwr, taxopath in rows:
					handle.write(f"{name}\t{alwr}\t{alwr}\t{alwr}\t{alwr}\t{taxopath}\n")
		return None

	monkeypatch.setattr("subprocess.run", fake_run)
	return calls, out_dir


def test_epa_ng_fallback_end_to_end(tmp_path: Path, monkeypatch):
	ref_tsv = make_reference_tsv(tmp_path, [("R1", "1", "a"), ("R2", "2", "b")])
	matrix = make_matrix(tmp_path, [("R1", "reference"), ("R2", "reference"), ("Q1", "query"), ("Q2", "query")])
	padded = tmp_path / "padded.fasta"
	padded.write_text(">R1\nATGC\n>R2\nATGA\n>Q1\nATGT\n>Q2\nATGG\n", encoding="utf-8")

	proc = processor(tmp_path, gb_matrix=matrix, reference_tsv=str(ref_tsv))  # no tree
	proc.padded_alignment = str(padded)
	calls, out_dir = fake_epa_environment(
		proc, tmp_path, monkeypatch,
		major_rows=[("Q1", "0.9", "1"), ("Q2", "0.8", "2")],
		minor_rows=[("Q1", "0.9", "a"), ("Q2", "0.8", "b")],
	)

	proc.process()

	tools = [Path(c[0]).name for c in calls]
	assert tools == ["iqtree3", "epa-ng", "gappa", "gappa"]

	rows = {r["primary_accession"]: r for r in read_gb_matrix(matrix)}
	assert (rows["Q1"]["major_clade"], rows["Q1"]["minor_clade"]) == ("1", "a")
	assert (rows["Q2"]["major_clade"], rows["Q2"]["minor_clade"]) == ("2", "b")


def test_epa_ng_derives_taxon_file_when_dataset_has_none(tmp_path: Path, monkeypatch):
	"""ref_list-only datasets have no curated taxon file; gappa still needs one."""
	ref_tsv = make_reference_tsv(tmp_path, [("R1", "1", "a"), ("R2", "2", "b")])
	matrix = make_matrix(tmp_path, [("R1", "reference"), ("Q1", "query")])
	padded = tmp_path / "padded.fasta"
	padded.write_text(">R1\nATGC\n>Q1\nATGT\n", encoding="utf-8")

	proc = processor(tmp_path, gb_matrix=matrix, reference_tsv=str(ref_tsv))
	proc.padded_alignment = str(padded)
	calls, out_dir = fake_epa_environment(
		proc, tmp_path, monkeypatch,
		major_rows=[("Q1", "0.9", "1")], minor_rows=[("Q1", "0.9", "a")],
	)

	proc.process()

	gappa_calls = [c for c in calls if Path(c[0]).name == "gappa"]
	major_taxon = Path(gappa_calls[0][gappa_calls[0].index("--taxon-file") + 1])
	assert major_taxon.name == "derived_major_clades.tsv"
	assert major_taxon.read_text(encoding="utf-8").strip() == "R1\t1\nR2\t2"


def test_epa_ng_uses_curated_taxon_file_when_present(tmp_path: Path, monkeypatch):
	major, minor = make_taxon_files(tmp_path, [("R1", "Arctic")], [("R1", "AL1a")])
	matrix = make_matrix(tmp_path, [("R1", "reference"), ("Q1", "query")])
	padded = tmp_path / "padded.fasta"
	padded.write_text(">R1\nATGC\n>Q1\nATGT\n", encoding="utf-8")

	proc = processor(tmp_path, gb_matrix=matrix, major=major, minor=minor)
	proc.padded_alignment = str(padded)
	calls, _ = fake_epa_environment(
		proc, tmp_path, monkeypatch,
		major_rows=[("Q1", "0.9", "Arctic")], minor_rows=[("Q1", "0.9", "AL1a")],
	)

	proc.process()

	gappa_calls = [c for c in calls if Path(c[0]).name == "gappa"]
	assert gappa_calls[0][gappa_calls[0].index("--taxon-file") + 1] == str(major)
	assert gappa_calls[1][gappa_calls[1].index("--taxon-file") + 1] == str(minor)


def test_epa_ng_command_construction(tmp_path: Path, monkeypatch):
	ref_tsv = make_reference_tsv(tmp_path, [("R1", "1", "a")])
	matrix = make_matrix(tmp_path, [("R1", "reference"), ("Q1", "query")])
	padded = tmp_path / "padded.fasta"
	padded.write_text(">R1\nATGC\n>Q1\nATGT\n", encoding="utf-8")

	proc = processor(tmp_path, gb_matrix=matrix, reference_tsv=str(ref_tsv))
	proc.padded_alignment = str(padded)
	calls, _ = fake_epa_environment(
		proc, tmp_path, monkeypatch,
		major_rows=[("Q1", "0.9", "1")], minor_rows=[("Q1", "0.9", "a")],
	)

	proc.process()

	iq = next(c for c in calls if Path(c[0]).name == "iqtree3")
	assert "-m" in iq and iq[iq.index("-m") + 1] == "GTR+G"
	assert iq[iq.index("-nt") + 1] == "2"
	assert iq[iq.index("-s") + 1].endswith("reference_aln.fa")

	epa = next(c for c in calls if Path(c[0]).name == "epa-ng")
	assert "--redo" in epa
	assert epa[epa.index("-s") + 1].endswith("reference_aln.fa")
	assert epa[epa.index("-q") + 1].endswith("query_aln.fa")
	assert epa[epa.index("-t") + 1].endswith("ref_tree.treefile")

	gappa = next(c for c in calls if Path(c[0]).name == "gappa")
	assert gappa[1:3] == ["examine", "assign"]
	assert "--per-query-results" in gappa
	assert gappa[gappa.index("--jplace-path") + 1].endswith("epa_result.jplace")


def test_alignment_splits_query_and_reference_by_matrix_type(tmp_path: Path):
	ref_tsv = make_reference_tsv(tmp_path, [("REF1", "1", "a")])
	matrix = make_matrix(tmp_path, [("REF1", "reference"), ("Q1", "query")])
	padded = tmp_path / "padded.fasta"
	padded.write_text(">REF1\nATGC\n>Q1\nATGT\n", encoding="utf-8")

	proc = processor(tmp_path, gb_matrix=matrix, reference_tsv=str(ref_tsv))
	proc.padded_alignment = str(padded)

	query_path, reference_path = proc.alignment()
	query_text = Path(query_path).read_text(encoding="utf-8")
	reference_text = Path(reference_path).read_text(encoding="utf-8")
	assert ">Q1" in query_text and ">REF1" not in query_text
	assert ">REF1" in reference_text and ">Q1" not in reference_text


def test_alignment_skips_excluded_rows(tmp_path: Path):
	"""exclusion_status=1 rows are not treated as queries."""
	matrix = tmp_path / "gB_matrix.tsv"
	write_gb_matrix(
		matrix,
		[
			{"primary_accession": "REF1", "gi_number": "REF1", "accession_type": "reference", "exclusion_status": "0"},
			{"primary_accession": "QBAD", "gi_number": "QBAD", "accession_type": "query", "exclusion_status": "1"},
		],
		DEFAULT_MATRIX_FIELDS,
	)
	padded = tmp_path / "padded.fasta"
	padded.write_text(">REF1\nATGC\n>QBAD\nATGT\n", encoding="utf-8")

	proc = processor(tmp_path, gb_matrix=matrix)
	proc.padded_alignment = str(padded)
	query_path, reference_path = proc.alignment()
	# Excluded query is not in the query set, so it lands on the reference side.
	assert ">QBAD" not in Path(query_path).read_text(encoding="utf-8")


# ---------------------------------------------------- gappa output parsing


def test_parse_gappa_per_query_keeps_best_hit_not_last_row(tmp_path: Path):
	"""Regression: the original code kept whichever row appeared LAST, which for a
	boundary query is the lowest-confidence alternative clade."""
	per_query = tmp_path / "per_query.tsv"
	per_query.write_text(
		"name\tLWR\tfract\taLWR\tafract\ttaxopath\n"
		"MK548366\t0.80\t0.80\t0.80\t0.80\tgenotype_1\n"
		"MK548366\t0.20\t0.20\t0.20\t0.20\tgenotype_2\n",
		encoding="utf-8",
	)
	assert CladeAssignment.parse_gappa_per_query(str(per_query)) == {"MK548366": "genotype_1"}


def test_parse_gappa_per_query_best_hit_regardless_of_row_order(tmp_path: Path):
	per_query = tmp_path / "per_query.tsv"
	per_query.write_text(
		"name\tLWR\tfract\taLWR\tafract\ttaxopath\n"
		"Q1\t0.10\t0.10\t0.10\t0.10\tlow\n"
		"Q1\t0.95\t0.95\t0.95\t0.95\thigh\n"
		"Q1\t0.50\t0.50\t0.50\t0.50\tmid\n",
		encoding="utf-8",
	)
	assert CladeAssignment.parse_gappa_per_query(str(per_query)) == {"Q1": "high"}


def test_parse_gappa_per_query_handles_multiple_queries(tmp_path: Path):
	per_query = tmp_path / "per_query.tsv"
	per_query.write_text(
		"name\tLWR\tfract\taLWR\tafract\ttaxopath\n"
		"Q1\t0.9\t0.9\t0.9\t0.9\tA\n"
		"Q2\t0.4\t0.4\t0.4\t0.4\tB\n"
		"Q2\t0.6\t0.6\t0.6\t0.6\tC\n",
		encoding="utf-8",
	)
	assert CladeAssignment.parse_gappa_per_query(str(per_query)) == {"Q1": "A", "Q2": "C"}


def test_parse_gappa_per_query_falls_back_to_lwr_without_alwr(tmp_path: Path):
	per_query = tmp_path / "per_query.tsv"
	per_query.write_text(
		"name\tLWR\ttaxopath\n"
		"Q1\t0.3\tlow\n"
		"Q1\t0.7\thigh\n",
		encoding="utf-8",
	)
	assert CladeAssignment.parse_gappa_per_query(str(per_query)) == {"Q1": "high"}


def test_parse_gappa_per_query_skips_blank_and_malformed_rows(tmp_path: Path):
	per_query = tmp_path / "per_query.tsv"
	per_query.write_text(
		"name\tLWR\taLWR\ttaxopath\n"
		"\t0.9\t0.9\torphan\n"        # no name
		"Q1\t0.9\t0.9\t\n"             # no taxopath
		"Q1\tnotanumber\tzzz\tvalid\n"  # unparseable score still usable
		,
		encoding="utf-8",
	)
	assert CladeAssignment.parse_gappa_per_query(str(per_query)) == {"Q1": "valid"}


def test_parse_gappa_per_query_missing_file_returns_empty():
	assert CladeAssignment.parse_gappa_per_query("/nonexistent/per_query.tsv") == {}
	assert CladeAssignment.parse_gappa_per_query(None) == {}


@pytest.mark.parametrize(
	"taxopath,expected",
	[
		("Bats;EF-W1", "EF-W1"),
		("Cosmopolitan", "Cosmopolitan"),
		("A;B;C", "C"),
		("", ""),
		("Bats;", "Bats"),
	],
)
def test_leaf_taxon(taxopath, expected):
	assert CladeAssignment._leaf_taxon(taxopath) == expected


# ---------------------------------------------------- iqtree binary resolution


def test_resolve_iqtree_binary_prefers_iqtree3(monkeypatch):
	monkeypatch.setattr("shutil.which", lambda name: "/env/bin/iqtree3" if name == "iqtree3" else None)
	assert CladeAssignment._resolve_iqtree_binary() == "iqtree3"


def test_resolve_iqtree_binary_falls_back_to_iqtree2(monkeypatch):
	monkeypatch.setattr("shutil.which", lambda name: "/env/bin/iqtree2" if name == "iqtree2" else None)
	assert CladeAssignment._resolve_iqtree_binary() == "iqtree2"


def test_resolve_iqtree_binary_raises_when_none_available(monkeypatch):
	monkeypatch.setattr("shutil.which", lambda name: None)
	with pytest.raises(FileNotFoundError, match="No iqtree binary"):
		CladeAssignment._resolve_iqtree_binary()


# =========================================================================
# Output writing
# =========================================================================


def test_writes_assignments_tsv_for_downstream_consumption(tmp_path: Path):
	major, minor = make_taxon_files(tmp_path, [("REF1", "1")], [("REF1", "a")])
	tree = tmp_path / "final-tree.nh"
	tree.write_text("((REF1:0.1,Q1:0.1):0.2,Q2:0.5);\n", encoding="utf-8")
	out = tmp_path / "clade_assignments.tsv"

	proc = processor(tmp_path, major=major, minor=minor, usher_tree=str(tree), assignments_out=str(out))
	proc.process()

	rows = {r["primary_accession"]: r for r in read_tsv(out)}
	assert rows["REF1"]["genotype"] == "1"
	assert rows["Q1"]["genotype"] == "1"
	assert set(rows["REF1"].keys()) == {"primary_accession", "genotype", "subtype"}


def test_gb_matrix_write_keys_on_primary_accession_and_overwrites_existing(tmp_path: Path):
	matrix = tmp_path / "gB_matrix.tsv"
	write_gb_matrix(
		matrix,
		[{"primary_accession": "Q1", "gi_number": "12345", "major_clade": "OLD", "minor_clade": "OLD"}],
		["primary_accession", "gi_number", "major_clade", "minor_clade"],
	)
	proc = processor(tmp_path, gb_matrix=matrix)
	proc.write_clades_to_gb_matrix({"Q1": "2"}, {"Q1": "b"})

	rows = read_gb_matrix(matrix)
	assert rows[0]["major_clade"] == "2"
	assert rows[0]["minor_clade"] == "b"
	# No duplicated clade columns on rewrite.
	assert rows[0].keys() == {"primary_accession", "gi_number", "major_clade", "minor_clade"}


def test_gb_matrix_write_is_optional(tmp_path: Path):
	"""Driver runs only need the assignments TSV; gb_matrix stays untouched."""
	major, minor = make_taxon_files(tmp_path, [("REF1", "1")], [("REF1", "a")])
	tree = tmp_path / "final-tree.nh"
	tree.write_text("(REF1:0.1,Q1:0.1);\n", encoding="utf-8")
	out = tmp_path / "assignments.tsv"

	proc = processor(tmp_path, major=major, minor=minor, usher_tree=str(tree), assignments_out=str(out))
	assert proc.gb_matrix is None
	proc.process()
	assert out.exists()
