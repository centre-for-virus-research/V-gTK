"""Carpet-bomb audit tests for alignment, coordinates, segmentation, and GenerateTables.

Adds direct passing coverage for GenerateTables, FastaHandler, and
CheckAlignmentLengths, plus strict xfail reproducers (AUDIT-ALN-*).
"""

from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

import AuditAlignmentCoordinates as aac_mod
import CalcAlignmentCord as cac_mod
import CheckAlignmentLengths as cal_mod
import CheckReferenceBackbone as crb_mod
import CheckSegmentRetention as csr_mod
import CollectFilteredSequences as cfs_mod
import FastaHandler as fh_mod
import GenerateTables as gt_mod
import InsertSegmentAlignments as isa_mod
import PadAlignment as pa_mod
import SegmentPivotTable as spt_mod
import ValidateSegment as vs_mod
import projectability as proj_mod
import reading_frame as rf_mod


# ---------------------------------------------------------------------------
# Passing coverage lifts: GenerateTables, FastaHandler, CheckAlignmentLengths
# ---------------------------------------------------------------------------

def test_generate_tables_creates_alignment_and_insertion_tables(tmp_path):
    base = tmp_path / "run"
    nextalign_dir = base / "Nextalign"
    (nextalign_dir / "query_aln" / "M1").mkdir(parents=True)
    (nextalign_dir / "reference_aln" / "M1").mkdir(parents=True)
    (nextalign_dir / "query_aln" / "M1" / "M1.aligned.fasta").write_text(">Q1\nATGCG\n", encoding="utf-8")
    (nextalign_dir / "query_aln" / "M1" / "M1.insertions.csv").write_text(
        "seqName,insertions,aaInsertions\nQ1,3:AAA,\n", encoding="utf-8"
    )
    (nextalign_dir / "reference_aln" / "M1" / "M1.aligned.fasta").write_text(">M1\nATGCG\n", encoding="utf-8")
    (nextalign_dir / "reference_aln" / "M1" / "M1.insertions.csv").write_text(
        "seqName,insertions,aaInsertions\n", encoding="utf-8"
    )

    gb = tmp_path / "gb.tsv"
    gb.write_text(
        "primary_accession\taccession_type\tsegment\thost\thost_taxa_id\n"
        "M1\tmaster\t1\tHuman\t9606\n"
        "Q1\tquery\t1\tHuman\t9606\n",
        encoding="utf-8",
    )
    bh = tmp_path / "tophits.tsv"
    bh.write_text("Q1\tM1\t99.0\t+\n", encoding="utf-8")
    pad = tmp_path / "refset_1_aln_merged_MSA.fasta"
    pad.write_text(">M1\nATGCG\n>Q1\nATGCG\n", encoding="utf-8")

    gt = gt_mod.GenerateTables(
        genbank_matrix=str(gb),
        base_dir=str(base),
        output_dir="Tables",
        blast_hits=str(bh),
        paded_aln=[str(pad)],
        host_taxa_file="host_taxa.tsv",
        nextalign_dir=str(nextalign_dir),
        email="test@example.org",
    )
    blast_dict = gt.load_blast_hits(str(bh))
    seg_map = gt.load_segment_map()
    gt.created_alignment_table(blast_dict, seg_map)
    gt.create_insertion_table(seg_map)

    aln_tsv = (base / "Tables" / "sequence_alignment.tsv").read_text(encoding="utf-8")
    ins_tsv = (base / "Tables" / "insertions.tsv").read_text(encoding="utf-8")
    assert "Q1\tM1\tATGCG\t1" in aln_tsv
    assert "Q1\tM1\t3:AAA\t1" in ins_tsv


def test_fasta_handler_same_path_deduplicates_sequences(tmp_path):
    fa = tmp_path / "seqs.fa"
    fa.write_text(">A\nATGC\n>B\nATGC\n>C\nGGGG\n", encoding="utf-8")
    fh_mod.RemoveRedundantSequence(str(fa), str(fa)).remove_redundant_fasta()
    txt = fa.read_text(encoding="utf-8")
    assert ">A" in txt and ">B" not in txt and ">C" in txt


def test_check_alignment_lengths_passes_uniform_and_fails_ragged(tmp_path):
    ok = tmp_path / "ok.fa"
    ok.write_text(">A\nATGC\n>B\nAT--\n", encoding="utf-8")
    assert cal_mod.check_uniform_lengths(str(ok))["is_uniform"] is True

    bad = tmp_path / "bad.fa"
    bad.write_text(">A\nATGC\n>B\nATG\n", encoding="utf-8")
    assert cal_mod.check_uniform_lengths(str(bad))["is_uniform"] is False


# ---------------------------------------------------------------------------
# Strict xfail reproducers for Alignment, Coords, Retention & Table bugs
# ---------------------------------------------------------------------------

def test_collect_filtered_sequences_guard_g1_fires_on_100pct_nextalign_failure(tmp_path, monkeypatch):
    nd = tmp_path / "Nextalign"
    (nd / "reference_aln" / "M1").mkdir(parents=True)
    (nd / "query_aln" / "M1").mkdir(parents=True)
    (nd / "reference_aln" / "M1" / "M1.aligned.fasta").write_text(">M1\nATGC\n", encoding="utf-8")
    # Empty aligned.fasta, all queries failed in errors.csv
    (nd / "query_aln" / "M1" / "M1.aligned.fasta").write_text("", encoding="utf-8")
    (nd / "query_aln" / "M1" / "M1.errors.csv").write_text(
        "seqName,errors,warnings\nQ1,failed to align,\nQ2,failed to align,\n", encoding="utf-8"
    )
    out = tmp_path / "filtered.tsv"
    monkeypatch.setattr(
        "sys.argv",
        ["CollectFilteredSequences.py", "-n", str(nd), "-o", str(out), "-b", ".", "--max_filtered_fraction", "0.05"],
    )
    with pytest.raises(SystemExit) as exc:
        cfs_mod.main()
    assert exc.value.code == 3


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-02: PadAlignment.insert_gaps turns query bases to '-' on 1-nt mismatch or U-vs-T")
def test_pad_alignment_insert_gaps_preserves_query_bases_on_ref_snp_or_u_vs_t():
    pad = pa_mod.PadAlignment.__new__(pa_mod.PadAlignment)
    ref_backbone = SeqRecord(Seq("ATGCCC"), id="R1", description="")
    # Nextalign ref has 1 SNP (A instead of G at pos 3)
    na_ref = SeqRecord(Seq("ATACCC"), id="R1", description="")
    q1 = SeqRecord(Seq("ATGTTT"), id="Q1", description="")
    out = pad.insert_gaps(ref_backbone, [na_ref, q1], "R1")
    q_proj = str([r for r in out if r.id == "Q1"][0].seq)
    # Query's 3rd base ('G') should NOT be erased to '-' just because R1 had a 1-nt difference
    assert q_proj[2] == "G"


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-03: CalcAlignmentCord emits og_start > og_end on whole-CDS internal deletion and skips L/M/S segments")
def test_calc_alignment_cord_skips_fully_deleted_cds_instead_of_inverted_og_bounds():
    cac = cac_mod.CalculateAlignmentCoordinates.__new__(cac_mod.CalculateAlignmentCoordinates)
    # Query covers cols 1..3 and 10..12, but CDS sits entirely inside deleted cols 4..9
    q_seq = "ATG------TTT"
    gap_ranges = cac.get_gap_ranges(q_seq)
    cds_list = [{"start": 4, "end": 9, "product": "P_deleted"}]
    res = cac.recalculate_cds_coordinates_with_span(
        "Q1", gap_ranges, cds_list, start_offset=0, genome_cord_start=1, genome_cord_end=12, master_coord_to_aln_pos=None
    )
    for item in res:
        assert item["og_start"] <= item["og_end"]


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-04: reading_frame.codon_grid off-by-one emits extra codon on 3k+2 CDS and mis-orients short minus-strand reads")
def test_reading_frame_codon_grid_exact_length_and_short_minus_strand():
    coord_map = {i: i - 1 for i in range(1, 20)}
    grid = rf_mod.codon_grid(coord_map, 1, 5)
    assert len(grid) == 1  # 5 nt holds only 1 full codon (1..3), NOT 2 codons!


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-05: CheckSegmentRetention skips non-numeric segments (L/M/S) and ignores vanished segment files")
def test_check_segment_retention_catches_vanished_segment_and_lms(tmp_path):
    ref_list = tmp_path / "refs.txt"
    ref_list.write_text("ML\tmaster\tL\nMS\tmaster\tS\n", encoding="utf-8")
    tophit = tmp_path / "tophit.tsv"
    tophit.write_text("Q1\tML\t99\t+\tL\nQ2\tMS\t99\t+\tS\n", encoding="utf-8")
    msa_l = tmp_path / "refset_L_aln_merged_MSA.fasta"
    msa_l.write_text(">ML\nATGC\n", encoding="utf-8")  # Q1 dropped -> 0% retention on L, and S vanished!
    exc = tmp_path / "filtered_ids.txt"
    exc.write_text("", encoding="utf-8")
    ok, _lines = csr_mod.check([str(msa_l)], str(tophit), str(ref_list), str(exc), min_retention=0.98)
    assert not ok


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-06: ValidateSegment / SegmentPivotTable drop literal 'NA' (Neuraminidase) segment")
def test_validate_segment_preserves_literal_na_segment_name(tmp_path):
    ref_list = tmp_path / "refs.txt"
    ref_list.write_text("M_NA\tmaster\tNA\nM_HA\tmaster\tHA\n", encoding="utf-8")
    valid = vs_mod.load_valid_segments(str(ref_list))
    assert "NA" in valid


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-07: GenerateTables crashes on 5-col BLAST & quoted multi-insertion CSV and maps PB2->2, PB1->1")
def test_generate_tables_handles_5col_blast_and_multi_insertion_csv(tmp_path):
    bh5 = tmp_path / "tophit_annotated.tsv"
    bh5.write_text("Q1\tM1\t99.0\t+\t1\n", encoding="utf-8")
    hits = gt_mod.GenerateTables.load_blast_hits(str(bh5))
    assert hits["Q1"] == "M1"


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-08: FastaHandler.RemoveRedundantSequence(in_fa, out_fa) overwrites in_fa instead of writing out_fa")
def test_fasta_handler_writes_distinct_output_file_without_clobbering_input(tmp_path):
    in_fa = tmp_path / "in.fa"
    out_fa = tmp_path / "out.fa"
    in_fa.write_text(">A\nATGC\n>B\nATGC\n", encoding="utf-8")
    fh_mod.RemoveRedundantSequence(str(in_fa), str(out_fa)).remove_redundant_fasta()
    assert out_fa.exists()
    assert ">B" in in_fa.read_text(encoding="utf-8")


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-09: CheckAlignmentLengths masks ragged copy when duplicate FASTA ID overwrites dict entry")
def test_check_alignment_lengths_flags_duplicate_id_with_ragged_length(tmp_path):
    fa = tmp_path / "dup_ragged.fa"
    fa.write_text(">A\nAT\n>A\nATGC\n>B\nATGC\n", encoding="utf-8")
    res = cal_mod.check_uniform_lengths(str(fa))
    assert not res["is_uniform"]


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-10: InsertSegmentAlignments exits 0 when --assert_width rejects ragged rows")
def test_insert_segment_alignments_fails_nonzero_on_assert_width_mismatch(tmp_path):
    import sqlite3
    db = tmp_path / "seg.db"
    conn = sqlite3.connect(str(db))
    conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT, segment TEXT)")
    conn.execute("INSERT INTO sequence_alignment VALUES ('M1', 'M1', 'ATGC', '1')")
    conn.commit()
    conn.close()

    aln = tmp_path / "bad.fa"
    aln.write_text(">M1\nATGC\n>Q1\nATG\n", encoding="utf-8")  # width 3 != 4
    bh = tmp_path / "top.tsv"
    bh.write_text("Q1\tM1\t99\t+\t1\n", encoding="utf-8")
    rc = isa_mod.main([
        "--db", str(db), "--segment", "1", "--padded_aln", str(aln),
        "--blast_tophit", str(bh), "--mode", "replace-segment", "--assert_width", "4",
    ])
    assert rc != 0


def test_projectability_does_not_misroute_v1_seg4_to_segment_1(tmp_path):
    bdir = tmp_path / "bb"
    bdir.mkdir()
    (bdir / "v1_seg4.fasta").write_text(">M4\nATGC\n", encoding="utf-8")
    assert proj_mod.find_precomputed_reference_alignment(str(bdir), "1") is None


@pytest.mark.xfail(strict=True, reason="AUDIT-ALN-12: CheckReferenceBackbone rejects documented single custom-named non-segmented backbone")
def test_check_reference_backbone_accepts_custom_non_segmented_filename(tmp_path):
    bb = tmp_path / "bb"
    bb.mkdir()
    (bb / "custom_backbone.fasta").write_text(">M1\nATGC\n>R1\nATGC\n", encoding="utf-8")
    ref = tmp_path / "refs.txt"
    ref.write_text("M1\tmaster\nR1\treference\n", encoding="utf-8")
    res = crb_mod.check_backbone(str(bb), str(ref), is_segmented="N")
    assert res["ok"]
