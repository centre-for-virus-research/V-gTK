"""BuildReferenceAlignment: the backbone every reference is projected through.

The guarantees under test are the ones the pipeline relies on downstream:

* the master row degaps to the master record (CalcAlignmentCord's 5' trim is 0);
* gaps inside a master gene are whole codons (codon reads stay in frame);
* every reference base is either in its row or in dropped_insertions.tsv;
* the output file resolves through projectability exactly as PadAlignment opens it.

Most tests drive the pure functions. The end-to-end tests need MAFFT and are
skipped without it.
"""

import csv
import shutil
from pathlib import Path

import pytest

import BuildReferenceAlignment as bra
import projectability

HAS_MAFFT = shutil.which("mafft") is not None
needs_mafft = pytest.mark.skipif(not HAS_MAFFT, reason="mafft not on PATH")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def test_translate_marks_stops_and_ambiguity():
	assert bra.translate("ATGTAANNNTGG") == "M*XW"
	assert bra.translate("ATGA") == "M"  # trailing partial codon ignored


def test_best_frame_prefers_the_frame_without_stops():
	coding = "ATGGCTAAAGGTCTGGAACGTTAA"  # MAKGLER*
	offset, n_codons, protein, stops = bra.best_frame("C" + coding)
	assert offset == 1
	assert protein == "MAKGLER*"
	assert stops == 0
	assert n_codons == 8


def test_back_translate_expands_gaps_to_codons():
	assert bra.back_translate("M-K", "ATGAAA") == "ATG---AAA"
	with pytest.raises(ValueError):
		bra.back_translate("MK", "ATGAAACCC")


def test_choose_codon_regions_keeps_longest_non_overlapping_clean_cds():
	master = "ATG" + "GCT" * 30 + "TAA" + "CCCCCC"
	cds = [
		{"start": 1, "end": 96, "strand": "+", "phase": "0", "product": "long"},
		{"start": 1, "end": 30, "strand": "+", "phase": "0", "product": "nested"},
		{"start": 1, "end": 20, "strand": "+", "phase": "0", "product": "not codons"},
		{"start": 1, "end": 96, "strand": "-", "phase": "0", "product": "minus"},
	]
	chosen, skipped = bra.choose_codon_regions(cds, master)
	assert [c["product"] for c in chosen] == ["long"]
	reasons = {entry["product"]: reason for entry, reason in skipped}
	assert reasons["not codons"] == "not a whole number of codons"
	assert reasons["minus"] == "not plus strand"


def test_choose_codon_regions_skips_a_master_cds_with_internal_stops():
	master = "ATGTAAGCTGCTTAA"
	chosen, skipped = bra.choose_codon_regions(
		[{"start": 1, "end": 15, "strand": "+", "phase": "0", "product": "broken"}], master)
	assert chosen == []
	assert skipped[0][1] == "internal stop in the master"


def test_mafft_strategy_uses_linsi_only_within_budget():
	assert bra.mafft_strategy(10, 1000, protein=False)[0] == "L-INS-i"
	assert bra.mafft_strategy(238, 9646, protein=False)[0] == "FFT-NS-i"
	assert bra.mafft_strategy(238, 3011, protein=True)[0] == "L-INS-i"
	assert bra.mafft_strategy(5000, 1000, protein=False)[0] == "FFT-NS-2"
	assert bra.mafft_strategy(10, 1000, protein=False, override="fftnsi")[0] == "FFT-NS-i"


def test_removed_bases_reports_runs_after_the_preceding_column():
	record = "AAACCCGGGTTT"
	row = "AAA---GGG-TT"
	# CCC missing after column 2, one T missing after column 8 (greedy left match
	# puts the dropped T at the end)
	runs = bra.removed_bases(record, row)
	assert runs == [(2, "CCC"), (11, "T")]


def test_removed_bases_rejects_a_row_the_record_cannot_produce():
	with pytest.raises(ValueError):
		bra.removed_bases("AAAA", "AAGA")


def test_master_gap_runs_out_of_frame_flags_non_codon_runs():
	master_row = "ATG---AAA--CCC"
	regions = [{"start": 1, "end": 9, "product": "gene"}]
	assert bra.master_gap_runs_out_of_frame(master_row, regions) == [("gene", 2)]


# ---------------------------------------------------------------------------
# Assembly and filtering
# ---------------------------------------------------------------------------

def test_assemble_and_filter_keep_codon_units_whole():
	skeleton = {
		"M": "CCATGAAACCCTAAGG",
		"R1": "CCATGAAACCCTAAGG",
	}
	region = {"start": 3, "end": 14, "product": "gene"}
	pieces = [{
		"M": ("", "ATG---AAACCCTAA", ""),
		"R1": ("T", "ATGGGGAAACCCTAA", ""),
	}]
	rows, units = bra.assemble(skeleton, "M", [region], pieces)
	assert rows["M"] == "CC-ATG---AAACCCTAAGG"
	assert rows["R1"] == "CCTATGGGGAAACCCTAAGG"
	codon_units = [u for u in units if u[2] == "codon"]
	assert all(width == 3 for _, width, _ in codon_units)

	# min support 2: the single-reference codon insertion and the single leftover
	# base both go, and they go whole.
	filtered, all_gap, low = bra.filter_columns(rows, "M", units, min_support=2)
	assert filtered["M"] == "CCATGAAACCCTAAGG"
	assert filtered["R1"] == "CCATGAAACCCTAAGG"
	assert low == 2 and all_gap == 0

	kept, _, low = bra.filter_columns(rows, "M", units, min_support=1)
	assert kept == rows and low == 0


# ---------------------------------------------------------------------------
# Reference grouping
# ---------------------------------------------------------------------------

def _write(path: Path, text: str) -> str:
	path.write_text(text, encoding="utf-8")
	return str(path)


def test_unsegmented_build_uses_one_group_with_the_master_first(tmp_path):
	ref_list = _write(tmp_path / "refs.tsv", "R1\treference\t1\nM\tmaster\t1\nR2\treference\t1\n")
	records = [("R1", "ACGT" * 10), ("M", "ACGT" * 10), ("R2", "ACGT" * 10)]
	groups = bra.group_references(records, ref_list, "N")
	assert groups == [("1", "M", ["M", "R1", "R2"])]


def test_unsegmented_master_without_segment_is_labelled_zero(tmp_path):
	ref_list = _write(tmp_path / "refs.tsv", "M\tmaster\nR1\treference\n")
	groups = bra.group_references([("M", "ACGT" * 5), ("R1", "ACGT" * 5)], ref_list, "N")
	assert groups[0][0] == "0"


def test_segmented_references_without_a_segment_go_to_the_closest_master(tmp_path):
	seg1 = "ACGTTGCAAGGCTTACCGATCGATCGGATCCATGCAAGTCC" * 4
	seg2 = "TTTTGGGGCCCCAAAATGTGTGCACACAGAGATCTCTC" * 4
	ref_list = _write(tmp_path / "refs.tsv", "M1\tmaster\t1\nM2\tmaster\t2\nR1\treference\t1\n")
	records = [("M1", seg1), ("M2", seg2), ("R1", seg1), ("X", seg2[:-4] + "AAAA")]
	groups = {label: members for label, _, members in bra.group_references(records, ref_list, "Y")}
	assert groups["1"] == ["M1", "R1"]
	assert groups["2"] == ["M2", "X"]


# ---------------------------------------------------------------------------
# End to end with MAFFT
# ---------------------------------------------------------------------------

GENE = ("ATGGCTAGCAAAGGTGAAGAACTGTTTACCGGTGTTGTTCCGATTCTGGTTGAACTGGATGGT"
		"GATGTTAACGGTCATAAATTTAGCGTTAGCGGTGAAGGTGAAGGTGATGCGACCTATGGTAAA"
		"CTGACCCTGAAATTTATTTGCACCACCGGTAAACTGCCGGTTCCGTGGCCGACCCTGGTTACC"
		"ACCCTGACCTATGGTGTTCAGTGCTTTAGCCGTTATCCGGATCATATGAAACGTCATGATTAA")
UTR5 = "GGACTTCAGTCGAACCTTGA"
UTR3 = "CCTGGATCACTAGGCTTTAA"


def _mutate(seq, every, base_cycle="ACGT"):
	out = list(seq)
	for i in range(5, len(out), every):
		if i % 3 == 2:  # third codon positions only, so the protein is unchanged
			out[i] = base_cycle[(base_cycle.index(out[i]) + 1) % 4] if out[i] in base_cycle else out[i]
	return "".join(out)


def _build_inputs(tmp_path):
	codon_insert = "GCTGCTGCT"  # three codons, shared by R1 and R2
	private_insert = "AAACCC"    # two codons, R3 only
	insert_at = 60               # after codon 20
	master = UTR5 + GENE + UTR3
	r1 = UTR5 + _mutate(GENE[:insert_at], 7) + codon_insert + GENE[insert_at:] + UTR3
	r2 = UTR5 + GENE[:insert_at] + codon_insert + _mutate(GENE[insert_at:], 11) + UTR3
	r3 = UTR5 + "T" + GENE[:120] + private_insert + GENE[120:] + UTR3
	r4 = UTR5[:-2] + GENE + UTR3  # shorter 5' UTR
	fasta = "".join(f">{acc}\n{seq}\n" for acc, seq in
					[("MASTER", master), ("R1", r1), ("R2", r2), ("R3", r3), ("R4", r4)])
	ref_fasta = _write(tmp_path / "refs.fa", fasta)
	ref_list = _write(tmp_path / "ref_list.txt",
					  "MASTER\tmaster\t1\nR1\treference\t1\nR2\treference\t1\nR3\treference\t1\nR4\treference\t1\n")
	gene_start = len(UTR5) + 1
	gene_end = len(UTR5) + len(GENE)
	gff = _write(tmp_path / "MASTER.gff3",
				 "##gff-version 3\n"
				 f"MASTER\tRefSeq\tregion\t1\t{len(master)}\t.\t+\t.\tID=MASTER\n"
				 f"MASTER\tRefSeq\tCDS\t{gene_start}\t{gene_end}\t.\t+\t0\tID=cds1;product=test protein\n")
	return ref_fasta, ref_list, gff, {"MASTER": master, "R1": r1, "R2": r2, "R3": r3, "R4": r4}


def _read_rows(path):
	return dict(bra.read_aligned_fasta(path))


@needs_mafft
def test_end_to_end_backbone_keeps_shared_codon_insertion_and_records_private_one(tmp_path):
	ref_fasta, ref_list, gff, raw = _build_inputs(tmp_path)
	out = tmp_path / "ref_set_aligned"
	reports = bra.build(ref_fasta, ref_list, [gff], str(out), threads=2, min_insertion_support=2)

	backbone = out / "refset_1_aln.fasta"
	assert backbone.exists()
	rows = _read_rows(backbone)
	assert list(rows)[0] == "MASTER"
	assert len({len(r) for r in rows.values()}) == 1

	# master row is the whole record, and its gene has only whole-codon gaps
	assert rows["MASTER"].replace("-", "") == raw["MASTER"]
	region = {"start": len(UTR5) + 1, "end": len(UTR5) + len(GENE), "product": "test protein"}
	assert bra.master_gap_runs_out_of_frame(rows["MASTER"], [region]) == []
	assert rows["MASTER"].count("-") >= 9  # the shared three-codon insertion has columns

	# R1 and R2 keep every base; R3's private codons are dropped and recorded
	assert rows["R1"].replace("-", "") == raw["R1"]
	assert rows["R2"].replace("-", "") == raw["R2"]
	with open(out / "dropped_insertions.tsv", newline="") as handle:
		dropped = {row["primary_accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
	assert "R3" in dropped
	recorded = sum(len(part.split(":", 1)[1]) for part in dropped["R3"]["insertion"].split(";"))
	assert len(rows["R3"].replace("-", "")) + recorded == len(raw["R3"])
	assert dropped["R3"]["reference"] == "MASTER"
	assert dropped["R3"]["segment"] == ""  # unsegmented build

	report = reports[0]
	assert report["codon_regions"] == 1
	assert report["references"] == 4
	assert (out / "build_report.tsv").exists()


@needs_mafft
def test_end_to_end_backbone_resolves_the_way_padalignment_opens_it(tmp_path):
	ref_fasta, ref_list, gff, _ = _build_inputs(tmp_path)
	out = tmp_path / "ref_set_aligned"
	bra.build(ref_fasta, ref_list, [gff], str(out), threads=2)
	resolved = projectability.find_precomputed_reference_alignment(str(out), "1")
	assert resolved == str(out / "refset_1_aln.fasta")
	ids = projectability.read_fasta_ids(resolved)
	assert ids == {"MASTER", "R1", "R2", "R3", "R4"}


@needs_mafft
def test_min_support_one_keeps_every_reference_base(tmp_path):
	ref_fasta, ref_list, gff, raw = _build_inputs(tmp_path)
	out = tmp_path / "ref_set_aligned"
	bra.build(ref_fasta, ref_list, [gff], str(out), threads=2, min_insertion_support=1)
	rows = _read_rows(out / "refset_1_aln.fasta")
	for acc, seq in raw.items():
		assert rows[acc].replace("-", "") == seq


@needs_mafft
def test_without_a_gff_the_backbone_is_nucleotide_only_but_still_valid(tmp_path):
	ref_fasta, ref_list, _, raw = _build_inputs(tmp_path)
	out = tmp_path / "ref_set_aligned"
	reports = bra.build(ref_fasta, ref_list, [], str(out), threads=2)
	rows = _read_rows(out / "refset_1_aln.fasta")
	assert rows["MASTER"].replace("-", "") == raw["MASTER"]
	assert reports[0]["codon_regions"] == 0


def test_main_reports_a_missing_reference_fasta_as_an_error(tmp_path, capsys):
	ref_list = _write(tmp_path / "refs.tsv", "M\tmaster\t1\n")
	code = bra.main(["--ref_fasta", str(tmp_path / "absent.fa"), "--ref_list", ref_list,
					 "--output_dir", str(tmp_path / "out")])
	assert code == 2
	assert "ERROR" in capsys.readouterr().err
