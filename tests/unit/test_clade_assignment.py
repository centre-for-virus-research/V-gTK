import csv
from pathlib import Path

import pytest

from CladeAssignment import CladeAssignment


def _write_gb_matrix(path, rows, fieldnames):
	with path.open("w", newline="", encoding="utf-8") as handle:
		writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
		writer.writeheader()
		for row in rows:
			writer.writerow(row)


def _read_gb_matrix(path):
	with path.open(newline="", encoding="utf-8") as handle:
		return list(csv.DictReader(handle, delimiter="\t"))


def _make_taxon_files(tmp_path, major_rows, minor_rows):
	major = tmp_path / "major_clades.tsv"
	minor = tmp_path / "minor_clades.tsv"
	major.write_text("".join(f"{acc}\t{clade}\n" for acc, clade in major_rows), encoding="utf-8")
	minor.write_text("".join(f"{acc}\t{clade}\n" for acc, clade in minor_rows), encoding="utf-8")
	return major, minor


def _processor(tmp_path, major, minor, gb_matrix, **kwargs):
	return CladeAssignment(
		major_clade=str(major),
		minor_clade=str(minor),
		padded_alignment=str(tmp_path / "padded.fasta"),
		base_dir=str(tmp_path),
		output_dir="CladeAssignment",
		gb_matrix=str(gb_matrix),
		threads="2",
		iqtree_model="GTR+G",
		**kwargs,
	)


def test_process_uses_tree_and_writes_query_clades(tmp_path: Path, monkeypatch):
	major, minor = _make_taxon_files(
		tmp_path,
		[("REF1", "1"), ("REF2", "2")],
		[("REF1", "a"), ("REF2", "b")],
	)
	gb_matrix = tmp_path / "gB_matrix.tsv"
	fieldnames = ["primary_accession", "gi_number", "accession_type", "exclusion_status"]
	_write_gb_matrix(
		gb_matrix,
		[
			{"primary_accession": "REF1", "gi_number": "REF1", "accession_type": "reference", "exclusion_status": "0"},
			{"primary_accession": "REF2", "gi_number": "REF2", "accession_type": "reference", "exclusion_status": "0"},
			{"primary_accession": "Q1", "gi_number": "Q1", "accession_type": "query", "exclusion_status": "0"},
			{"primary_accession": "Q2", "gi_number": "Q2", "accession_type": "query", "exclusion_status": "0"},
		],
		fieldnames,
	)
	usher_tree = tmp_path / "final-tree.nh"
	# Q1 clusters with REF1 (clade 1a), Q2 with REF2 (clade 2b).
	usher_tree.write_text("((REF1:0.1,Q1:0.1):0.2,(REF2:0.1,Q2:0.1):0.2);\n", encoding="utf-8")

	processor = _processor(tmp_path, major, minor, gb_matrix, usher_tree=str(usher_tree))

	# EPA-ng must not run when a tree is available.
	def fail_epa(*args, **kwargs):
		raise AssertionError("EPA-ng fallback should not run when a tree is supplied")
	monkeypatch.setattr(processor, "assign_with_epa_ng", fail_epa)

	processor.process()

	rows = {r["primary_accession"]: r for r in _read_gb_matrix(gb_matrix)}
	assert (rows["Q1"]["major_clade"], rows["Q1"]["minor_clade"]) == ("1", "a")
	assert (rows["Q2"]["major_clade"], rows["Q2"]["minor_clade"]) == ("2", "b")
	# References keep their curated clade.
	assert (rows["REF1"]["major_clade"], rows["REF1"]["minor_clade"]) == ("1", "a")
	assert (rows["REF2"]["major_clade"], rows["REF2"]["minor_clade"]) == ("2", "b")


def test_process_falls_back_to_epa_ng_when_no_tree(tmp_path: Path, monkeypatch):
	major, minor = _make_taxon_files(tmp_path, [("REF1", "1")], [("REF1", "a")])
	gb_matrix = tmp_path / "gB_matrix.tsv"
	fieldnames = ["primary_accession", "gi_number", "accession_type", "exclusion_status"]
	_write_gb_matrix(
		gb_matrix,
		[
			{"primary_accession": "REF1", "gi_number": "REF1", "accession_type": "reference", "exclusion_status": "0"},
			{"primary_accession": "Q1", "gi_number": "Q1", "accession_type": "query", "exclusion_status": "0"},
		],
		fieldnames,
	)

	processor = _processor(tmp_path, major, minor, gb_matrix)  # no tree

	called = {"epa": False}
	def fake_epa():
		called["epa"] = True
		return {"Q1": "5"}, {"Q1": "z"}
	monkeypatch.setattr(processor, "assign_with_epa_ng", fake_epa)

	processor.process()

	assert called["epa"] is True
	rows = {r["primary_accession"]: r for r in _read_gb_matrix(gb_matrix)}
	assert (rows["Q1"]["major_clade"], rows["Q1"]["minor_clade"]) == ("5", "z")


def test_parse_gappa_per_query_keeps_best_hit_not_last_row(tmp_path: Path):
	per_query = tmp_path / "per_query.tsv"
	# gappa lists several candidate taxopaths per query. The genotype-1 row has the
	# highest aLWR; a spurious genotype-2 row appears LAST in the file.
	per_query.write_text(
		"name\tLWR\tfract\taLWR\tafract\ttaxopath\n"
		"MK548366\t0.80\t0.80\t0.80\t0.80\tgenotype_1\n"
		"MK548366\t0.20\t0.20\t0.20\t0.20\tgenotype_2\n",
		encoding="utf-8",
	)
	best = CladeAssignment.parse_gappa_per_query(str(per_query))
	assert best == {"MK548366": "genotype_1"}


def test_parse_gappa_per_query_takes_leaf_of_taxopath():
	# End-to-end leaf extraction from a semicolon taxopath.
	assert CladeAssignment._leaf_taxon("Bats;EF-W1") == "EF-W1"
	assert CladeAssignment._leaf_taxon("Cosmopolitan") == "Cosmopolitan"
	assert CladeAssignment._leaf_taxon("") == ""


def test_gb_matrix_write_keys_on_primary_accession_and_overwrites_existing(tmp_path: Path):
	major, minor = _make_taxon_files(tmp_path, [], [])
	gb_matrix = tmp_path / "gB_matrix.tsv"
	# gi_number differs from primary_accession, and stale clade columns exist.
	fieldnames = ["primary_accession", "gi_number", "major_clade", "minor_clade"]
	_write_gb_matrix(
		gb_matrix,
		[
			{"primary_accession": "Q1", "gi_number": "12345", "major_clade": "OLD", "minor_clade": "OLD"},
		],
		fieldnames,
	)

	processor = _processor(tmp_path, major, minor, gb_matrix)
	processor.write_clades_to_gb_matrix({"Q1": "2"}, {"Q1": "b"})

	rows = _read_gb_matrix(gb_matrix)
	assert rows[0]["major_clade"] == "2"
	assert rows[0]["minor_clade"] == "b"
	# Columns are not duplicated on rewrite.
	assert rows[0].keys() == {"primary_accession", "gi_number", "major_clade", "minor_clade"}


def test_read_taxon_file_tolerates_header_and_blank_rows(tmp_path: Path):
	taxon = tmp_path / "taxon.tsv"
	taxon.write_text(
		"accession\tclade\n"      # header row, skipped
		"REF1\tArctic\n"
		"REF2\t\n"                  # blank clade, skipped
		"REF3\tAsian\n",
		encoding="utf-8",
	)
	mapping = CladeAssignment._read_taxon_file(str(taxon))
	assert mapping == {"REF1": "Arctic", "REF3": "Asian"}


def test_alignment_splits_query_and_reference_by_matrix_type(tmp_path: Path):
	major, minor = _make_taxon_files(tmp_path, [("REF1", "1")], [("REF1", "a")])
	gb_matrix = tmp_path / "gB_matrix.tsv"
	fieldnames = ["primary_accession", "gi_number", "accession_type", "exclusion_status"]
	_write_gb_matrix(
		gb_matrix,
		[
			{"primary_accession": "REF1", "gi_number": "REF1", "accession_type": "reference", "exclusion_status": "0"},
			{"primary_accession": "Q1", "gi_number": "Q1", "accession_type": "query", "exclusion_status": "0"},
		],
		fieldnames,
	)
	padded = tmp_path / "padded.fasta"
	padded.write_text(">REF1\nATGC\n>Q1\nATGT\n", encoding="utf-8")

	processor = _processor(tmp_path, major, minor, gb_matrix)
	processor.padded_alignment = str(padded)

	query_path, reference_path = processor.alignment()
	assert ">Q1" in Path(query_path).read_text(encoding="utf-8")
	assert ">REF1" not in Path(query_path).read_text(encoding="utf-8")
	assert ">REF1" in Path(reference_path).read_text(encoding="utf-8")
