import shutil
from pathlib import Path
import types
import sqlite3

import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from PadAlignment import PadAlignment


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "pad_alignment"
# The DB this pointed at (rabv-jul0425.db) has never existed in the repo, so seven
# tests - four of them update-mode - skipped silently in CI and locally. This one
# is tracked and its features schema matches the fixtures below.
REAL_UPDATE_DB = REPO_ROOT / "test_data" / "RABV_test" / "rabv-11Jun26.db"


def _pick_master_and_segment_from_real_db():
    if not REAL_UPDATE_DB.exists():
        pytest.skip(f"Real update DB not found at {REAL_UPDATE_DB}")

    conn = sqlite3.connect(str(REAL_UPDATE_DB))
    try:
        cur = conn.cursor()
        cur.execute(
            "SELECT primary_accession, segment FROM meta_data "
            "WHERE lower(coalesce(accession_type,''))='master' "
            "AND trim(coalesce(segment,'')) != '' LIMIT 1"
        )
        row = cur.fetchone()
    finally:
        conn.close()

    if not row:
        pytest.skip("No master accession found in real update DB")
    master = str(row[0]).strip()
    segment = str(row[1]).strip() if row[1] is not None else "1"
    return master, segment


def _copy_fixture_tree(tmp_path: Path):
    input_dir = tmp_path / "query_aln"
    (input_dir / "REF_OK").mkdir(parents=True, exist_ok=True)
    (input_dir / "REF_ORPHAN").mkdir(parents=True, exist_ok=True)

    shutil.copyfile(DATA_DIR / "REF_OK.aligned.fasta", input_dir / "REF_OK" / "REF_OK.aligned.fasta")
    shutil.copyfile(
        DATA_DIR / "REF_ORPHAN.aligned.fasta",
        input_dir / "REF_ORPHAN" / "REF_ORPHAN.aligned.fasta",
    )

    return input_dir


def _write_fasta(path: Path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")


def _write_update_alignment_db(path: Path, meta_rows, alignment_rows):
    conn = sqlite3.connect(str(path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)"
        )
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
            meta_rows,
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, alignment_name, alignment, segment) VALUES (?, ?, ?, ?)",
            alignment_rows,
        )
        conn.commit()
    finally:
        conn.close()


def _make_processor(tmp_path: Path):
    return PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        new_outputfile=False,
    )


def test_find_orphan_query_references_identifies_unprojectable_refs(tmp_path: Path):
    input_dir = _copy_fixture_tree(tmp_path)
    processor = _make_processor(tmp_path)

    orphan_refs = processor.find_orphan_query_references(
        str(DATA_DIR / "master_ref.aligned.fasta"),
        str(input_dir),
    )

    assert set(orphan_refs.keys()) == {"REF_ORPHAN"}
    assert sorted(orphan_refs["REF_ORPHAN"]) == ["Q_ORPHAN_1", "Q_ORPHAN_2"]


def test_process_master_alignment_warns_and_skips_orphan_refs(tmp_path: Path, capsys):
    input_dir = _copy_fixture_tree(tmp_path)
    processor = _make_processor(tmp_path)

    processor.process_master_alignment(
        reference_alignment_file=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(input_dir),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
    )

    out = capsys.readouterr().out
    assert "[warn]" in out
    assert "REF_ORPHAN" in out

    merged_path = tmp_path / "pad_out" / "master_ref.aligned_merged_MSA.fasta"
    assert merged_path.exists()

    merged_ids = {r.id for r in SeqIO.parse(str(merged_path), "fasta")}
    assert "Q_OK_1" in merged_ids
    assert "Q_OK_2" in merged_ids
    assert "Q_ORPHAN_1" not in merged_ids
    assert "Q_ORPHAN_2" not in merged_ids


def test_process_all_masters_uses_precomputed_by_segment_then_falls_back(tmp_path: Path):
    processor = _make_processor(tmp_path)

    precomputed_dir = tmp_path / "ref_set_aligned"
    precomputed_dir.mkdir(parents=True, exist_ok=True)
    (precomputed_dir / "refset_1_aln.fasta").write_text(">REF_OK\nACGT\n", encoding="utf-8")

    nextalign_dir = tmp_path / "Nextalign"
    (nextalign_dir / "reference_aln" / "MASTER2").mkdir(parents=True, exist_ok=True)
    (nextalign_dir / "reference_aln" / "MASTER2" / "MASTER2.aligned.fasta").write_text(">MASTER2\nACGT\n", encoding="utf-8")

    called_refs = []

    def _capture_process(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False, segment_value=None):
        called_refs.append(Path(reference_alignment_file).name)

    processor.process_master_alignment = types.MethodType(_capture_process, processor)

    processor.process_all_masters(
        master_list=["MASTER1", "MASTER2"],
        nextalign_dir=str(nextalign_dir),
        master_segment_map={"MASTER1": "1", "MASTER2": "2"},
        precomputed_ref_dir=str(precomputed_dir),
    )

    assert "refset_1_aln.fasta" in called_refs
    assert "MASTER2.aligned.fasta" in called_refs


def test_insert_gaps_pads_when_sequence_shorter_than_reference(tmp_path: Path):
    processor = _make_processor(tmp_path)

    reference_aligned = "ACGT--ACGT"
    subalignment = [SeqRecord(Seq("ACGTAC"), id="Q1")]

    updated = processor.insert_gaps(reference_aligned, subalignment, "Q1")
    assert len(updated) == 1
    assert str(updated[0].seq) == "ACGT--AC--"
    assert len(str(updated[0].seq)) == len(reference_aligned)


def test_insert_gaps_with_insertion_coordinate_shift(tmp_path: Path):
    processor = _make_processor(tmp_path)

    # reference_aligned has length 10, ungapped length 8: 'ACGTACGT'
    reference_aligned = "ACGT--ACGT"
    
    # Subalignment has:
    # 1. The genotype-specific reference containing insertions relative to master:
    #    'AA-CGTACGT' -> ungapped sequence 'AACGTACGT' (has extra 'AA' at start)
    # 2. A query sequence aligned to it, starting after the extra 'AA':
    #    '---CGTAC--' -> representing the query aligned
    subalignment = [
        SeqRecord(Seq("AA-CGTACGT"), id="REF1"),
        SeqRecord(Seq("---CGTAC--"), id="Q1"),
    ]

    updated = processor.insert_gaps(reference_aligned, subalignment, "REF1")
    assert len(updated) == 2
    
    # We expect Q1 to align correctly relative to the master-aligned reference coordinates.
    # Master-aligned reference starts with 'A' (of 'ACGT...').
    # But Q1's first matching nucleotide ('C' at index 3 of '---CGT...') aligns to the 'C' of 'REF1'
    # which is the third nucleotide of 'REF1' (non-gap coordinate 2, since first is 'A' (0), second 'A' (1), third 'C' (2)).
    # So 'C' should map to the third nucleotide of master_dq_raw ('C' of 'ACGT...'),
    # which is at index 1 of reference_aligned ('A'(0), 'C'(1)).
    # So gapped Q1 should start at index 1 of the projected sequence!
    # Expected alignment:
    # reference_aligned: A C G T - - A C G T
    # projected Q1:      - C G T - - A C - -
    assert str(updated[1].seq) == "-CGT--AC--"
    assert len(str(updated[1].seq)) == len(reference_aligned)


def test_insert_gaps_with_subalignment_fallback(tmp_path: Path):
    processor = _make_processor(tmp_path)

    reference_aligned = "ACGT--ACGT"
    # Q1 is in subalignment but ref_id is passed as 'NON_EXISTENT'
    subalignment = [SeqRecord(Seq("ACGTAC"), id="Q1")]

    # Should fall back to subalignment[0] ('Q1') as the nextalign reference and run cleanly
    updated = processor.insert_gaps(reference_aligned, subalignment, "NON_EXISTENT")
    assert len(updated) == 1
    assert str(updated[0].seq) == "ACGT--AC--"


def test_insert_gaps_exact_match(tmp_path: Path):
    processor = _make_processor(tmp_path)
    reference_aligned = "ACGTACGT"
    subalignment = [
        SeqRecord(Seq("ACGTACGT"), id="REF1"),
        SeqRecord(Seq("ACGTACGT"), id="Q1")
    ]
    updated = processor.insert_gaps(reference_aligned, subalignment, "REF1")
    assert len(updated) == 2
    assert str(updated[0].seq) == "ACGTACGT"
    assert str(updated[1].seq) == "ACGTACGT"


def test_insert_gaps_with_deletion_in_subalignment_ref(tmp_path: Path):
    processor = _make_processor(tmp_path)
    reference_aligned = "ACGT--ACGT"
    subalignment = [
        SeqRecord(Seq("AC--ACGT"), id="REF1"),
        SeqRecord(Seq("AC--ACGT"), id="Q1")
    ]
    updated = processor.insert_gaps(reference_aligned, subalignment, "REF1")
    assert len(updated) == 2
    assert str(updated[1].seq) == "AC----ACGT"


def test_insert_gaps_multiple_queries(tmp_path: Path):
    processor = _make_processor(tmp_path)
    reference_aligned = "ACGT--ACGT"
    subalignment = [
        SeqRecord(Seq("AA-CGTACGT"), id="REF1"),
        SeqRecord(Seq("---CGTAC--"), id="Q1"),
        SeqRecord(Seq("AA-CGT----"), id="Q2"),
        SeqRecord(Seq("------ACGT"), id="Q3"),
    ]
    updated = processor.insert_gaps(reference_aligned, subalignment, "REF1")
    assert len(updated) == 4
    assert str(updated[1].seq) == "-CGT--AC--"
    assert str(updated[2].seq) == "ACGT------"
    assert str(updated[3].seq) == "------ACGT"


def test_insert_gaps_query_ends_early(tmp_path: Path):
    processor = _make_processor(tmp_path)
    reference_aligned = "ACGT--ACGT"
    subalignment = [
        SeqRecord(Seq("AA-CGTACGT"), id="REF1"),
        SeqRecord(Seq("AA-CGT----"), id="Q1"),
    ]
    updated = processor.insert_gaps(reference_aligned, subalignment, "REF1")
    assert len(updated) == 2
    assert str(updated[1].seq) == "ACGT------"


def test_insert_gaps_query_starts_late(tmp_path: Path):
    processor = _make_processor(tmp_path)
    reference_aligned = "ACGT--ACGT"
    subalignment = [
        SeqRecord(Seq("AA-CGTACGT"), id="REF1"),
        SeqRecord(Seq("------ACGT"), id="Q1"),
    ]
    updated = processor.insert_gaps(reference_aligned, subalignment, "REF1")
    assert len(updated) == 2
    assert str(updated[1].seq) == "------ACGT"


def test_insert_gaps_with_multiple_complex_shifts(tmp_path: Path):
    processor = _make_processor(tmp_path)
    reference_aligned = "ACGT--ACGTACGT"
    subalignment = [
        SeqRecord(Seq("AA-CGT--AC--ACGT"), id="REF1"),
        SeqRecord(Seq("---CGT--AC--ACGT"), id="Q1"),
    ]
    updated = processor.insert_gaps(reference_aligned, subalignment, "REF1")
    assert len(updated) == 2
    assert str(updated[1].seq) == "-CGT--AC--ACGT"


def test_process_all_masters_strict_mode_fails_when_db_segment_backbone_missing(tmp_path: Path):
    master, segment = _pick_master_and_segment_from_real_db()

    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        strict_segment_backbone=True,
    )

    precomputed_dir = tmp_path / "ref_backbones"
    precomputed_dir.mkdir(parents=True, exist_ok=True)

    wrong_segment = "999" if segment != "999" else "998"
    (precomputed_dir / f"refset_{wrong_segment}_aln.fasta").write_text(">DUMMY\nATGC\n", encoding="utf-8")

    nextalign_dir = tmp_path / "Nextalign"
    nextalign_dir.mkdir(parents=True, exist_ok=True)

    with pytest.raises(ValueError, match="Missing DB/precomputed backbone"):
        processor.process_all_masters(
            master_list=[master],
            nextalign_dir=str(nextalign_dir),
            master_segment_map={master: segment},
            precomputed_ref_dir=str(precomputed_dir),
        )


def test_process_all_masters_ambiguous_db_segment_backbones_pick_deterministic_file(tmp_path: Path):
    master, segment = _pick_master_and_segment_from_real_db()

    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        strict_segment_backbone=False,
    )

    precomputed_dir = tmp_path / "ref_backbones"
    precomputed_dir.mkdir(parents=True, exist_ok=True)

    # Create ambiguous numeric matches for the same segment.
    file_b = precomputed_dir / f"zzz_segment_{segment}.fasta"
    file_a = precomputed_dir / f"aaa_segment_{segment}.fasta"
    file_b.write_text(">REF_B\nATGC\n", encoding="utf-8")
    file_a.write_text(">REF_A\nATGC\n", encoding="utf-8")

    chosen = {}

    def _capture_process(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False, segment_value=None):
        chosen["ref"] = Path(reference_alignment_file).name
        chosen["segment"] = str(segment_value)

    processor.process_master_alignment = types.MethodType(_capture_process, processor)

    processor.process_all_masters(
        master_list=[master],
        nextalign_dir=str(tmp_path / "Nextalign"),
        master_segment_map={master: segment},
        precomputed_ref_dir=str(precomputed_dir),
    )

    assert chosen["segment"] == segment
    assert chosen["ref"] == "aaa_segment_" + segment + ".fasta"


def test_export_update_backbones_writes_segment_fastas(tmp_path: Path, basic_update_db: Path):
    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        update_db=str(basic_update_db),
    )

    out_dir = tmp_path / "db_ref_backbones"
    resolved = processor.export_update_backbones(str(out_dir))

    assert resolved == str(out_dir)
    assert (out_dir / "refset_1_aln.fasta").exists()
    assert (out_dir / "refset_2_aln.fasta").exists()
    assert ">REF1" in (out_dir / "refset_1_aln.fasta").read_text(encoding="utf-8")


def test_get_master_segment_map_uses_update_db_and_normalizes_blank_segment(tmp_path: Path):
    update_db = tmp_path / "update.db"
    conn = sqlite3.connect(str(update_db))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
            [("NC_001542", "master", ""), ("REF2", "reference", "")],
        )
        conn.commit()
    finally:
        conn.close()

    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        update_db=str(update_db),
    )

    assert processor.get_master_segment_map("ignored.tsv") == {"NC_001542": "0"}


def test_process_all_masters_update_db_blank_segment_uses_refset_zero(tmp_path: Path):
    update_db = tmp_path / "update.db"
    conn = sqlite3.connect(str(update_db))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
            [("NC_001542", "master", ""), ("REF2", "reference", "")],
        )
        conn.commit()
    finally:
        conn.close()

    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        strict_segment_backbone=True,
        update_db=str(update_db),
    )

    precomputed_dir = tmp_path / "ref_backbones"
    precomputed_dir.mkdir(parents=True, exist_ok=True)
    (precomputed_dir / "refset_0_aln.fasta").write_text(">NC_001542\nATGC\n", encoding="utf-8")

    chosen = {}

    def _capture_process(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False, segment_value=None):
        chosen["ref"] = Path(reference_alignment_file).name
        chosen["segment"] = str(segment_value)

    processor.process_master_alignment = types.MethodType(_capture_process, processor)
    processor.process_all_masters(
        master_list=processor.get_master_list("ignored.tsv"),
        nextalign_dir=str(tmp_path / "Nextalign"),
        master_segment_map=processor.get_master_segment_map("ignored.tsv"),
        precomputed_ref_dir=str(precomputed_dir),
    )

    assert chosen["ref"] == "refset_0_aln.fasta"
    assert chosen["segment"] == "0"


def test_export_update_backbones_deduplicates_same_accession(tmp_path: Path):
    update_db = tmp_path / "update.db"
    conn = sqlite3.connect(str(update_db))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT)"
        )
        cur.execute(
            "INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
            ("REF1", "master", ""),
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, alignment_name, alignment) VALUES (?, ?, ?)",
            [
                ("REF1", "OTHER_BACKBONE", "AAAA"),
                ("REF1", "REF1", "TTTT"),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        update_db=str(update_db),
    )

    out_dir = tmp_path / "db_ref_backbones"
    processor.export_update_backbones(str(out_dir))
    content = (out_dir / "refset_0_aln.fasta").read_text(encoding="utf-8")
    assert content.count(">REF1") == 1
    assert "TTTT" in content


def test_export_update_backbones_normalizes_text_segment_labels(tmp_path: Path):
    update_db = tmp_path / "update.db"
    _write_update_alignment_db(
        update_db,
        meta_rows=[
            ("MASTER2", "master", "segment 02"),
            ("REF2A", "reference", "segment-02"),
        ],
        alignment_rows=[
            ("MASTER2", "MASTER2", "A-CGT", "segment 02"),
            ("REF2A", "REF2A", "ATCGT", "segment-02"),
        ],
    )

    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        update_db=str(update_db),
    )

    out_dir = tmp_path / "db_ref_backbones"
    processor.export_update_backbones(str(out_dir))

    assert processor.get_master_segment_map("ignored.tsv") == {"MASTER2": "2"}
    assert (out_dir / "refset_2_aln.fasta").exists()
    content = (out_dir / "refset_2_aln.fasta").read_text(encoding="utf-8")
    assert ">MASTER2" in content
    assert ">REF2A" in content


def test_resolve_precomputed_ref_dir_prefers_explicit_dir_over_db_export(tmp_path: Path):
    update_db = tmp_path / "update.db"
    _write_update_alignment_db(
        update_db,
        meta_rows=[("MASTER0", "master", "")],
        alignment_rows=[("MASTER0", "MASTER0", "ATGC", "")],
    )

    explicit_dir = tmp_path / "explicit_backbones"
    explicit_dir.mkdir(parents=True, exist_ok=True)
    (explicit_dir / "refset_0_aln.fasta").write_text(">MASTER0\nTTTT\n", encoding="utf-8")

    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        update_db=str(update_db),
    )

    resolved = processor.resolve_precomputed_ref_dir(str(explicit_dir))

    assert resolved == str(explicit_dir)
    assert not (tmp_path / "pad_out" / "db_ref_backbones").exists()


def test_process_all_masters_update_mode_end_to_end_handles_orphans_and_reference_only_segments(tmp_path: Path, capsys):
    update_db = tmp_path / "update.db"
    _write_update_alignment_db(
        update_db,
        meta_rows=[
            ("MASTER0", "master", ""),
            ("REF_A", "reference", ""),
            ("REF_B", "reference", ""),
            ("REF_C", "reference", ""),
        ],
        alignment_rows=[
            ("MASTER0", "MASTER0", "A-CGT", ""),
            ("REF_A", "REF_A", "ATCGT", ""),
            ("REF_B", "REF_B", "A-GGT", ""),
            ("REF_C", "REF_C", "AC-GT", ""),
            ("REF_A", "LEGACY_CLUSTER", "TTTTT", ""),
        ],
    )

    query_dir = tmp_path / "query_aln"
    _write_fasta(query_dir / "REF_A" / "REF_A.aligned.fasta", [("Q_A_1", "ATGGT")])
    _write_fasta(query_dir / "REF_C" / "REF_C.aligned.fasta", [("Q_C_1", "ACAGT")])
    _write_fasta(query_dir / "ORPHAN_REF" / "ORPHAN_REF.aligned.fasta", [("Q_ORPHAN", "TTTT")])

    processor = PadAlignment(
        reference_alignment=None,
        input_dir=str(query_dir),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        segment_manifest_out=str(tmp_path / "pad_manifest.tsv"),
        strict_segment_backbone=True,
        update_db=str(update_db),
    )

    resolved = processor.resolve_precomputed_ref_dir(None)
    processor.process_all_masters(
        master_list=processor.get_master_list("ignored.tsv"),
        nextalign_dir=str(tmp_path / "Nextalign"),
        master_segment_map=processor.get_master_segment_map("ignored.tsv"),
        precomputed_ref_dir=resolved,
    )

    out = capsys.readouterr().out
    assert "Using DB-derived precomputed segment alignments" in out
    assert "ORPHAN_REF" in out

    merged = tmp_path / "pad_out" / "refset_0_aln_merged_MSA.fasta"
    assert merged.exists()
    merged_records = list(SeqIO.parse(str(merged), "fasta"))
    merged_ids = [record.id for record in merged_records]
    assert merged_ids.count("REF_A") == 1
    assert "MASTER0" in merged_ids
    assert "REF_B" in merged_ids
    assert "Q_A_1" in merged_ids
    assert "Q_C_1" in merged_ids
    assert "Q_ORPHAN" not in merged_ids

    by_id = {record.id: str(record.seq) for record in merged_records}
    assert by_id["Q_A_1"] == "ATGGT"
    assert by_id["Q_C_1"] == "AC-AG"
    assert all(row["segment"] == "0" for row in processor.segment_manifest_rows)


def test_process_all_masters_uses_only_backbone_file_when_segment_unknown_and_single_candidate_exists(tmp_path: Path):
    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
        strict_segment_backbone=True,
    )

    precomputed_dir = tmp_path / "ref_backbones"
    precomputed_dir.mkdir(parents=True, exist_ok=True)
    (precomputed_dir / "sole_backbone.fasta").write_text(">MASTERX\nATGC\n", encoding="utf-8")

    chosen = {}

    def _capture_process(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False, segment_value=None):
        chosen["ref"] = Path(reference_alignment_file).name
        chosen["segment"] = segment_value

    processor.process_master_alignment = types.MethodType(_capture_process, processor)
    processor.process_all_masters(
        master_list=["MASTERX"],
        nextalign_dir=str(tmp_path / "Nextalign"),
        master_segment_map={"MASTERX": None},
        precomputed_ref_dir=str(precomputed_dir),
    )

    assert chosen["ref"] == "sole_backbone.fasta"
    assert chosen["segment"] is None


def _guide_processor(tmp_path: Path, guide, skip_ids=None):
    """A processor wired to a caller-written guide alignment under ``tmp_path``."""
    return PadAlignment(
        reference_alignment=str(guide),
        input_dir=str(tmp_path / "in"),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "out"),
        keep_intermediate_files=True,
        skip_ids=skip_ids,
    )


def _merged_records(tmp_path: Path, guide: Path):
    merged = tmp_path / "out" / guide.name.replace(".fasta", "_merged_MSA.fasta")
    assert merged.exists(), f"no merged MSA written at {merged}"
    return list(SeqIO.parse(str(merged), "fasta"))


class TestAccessionIdentityIsCanonical:
    """One genome, one spelling of its name.

    Everything this script joins on is the bare accession: the ``query_aln/<ref>/``
    directory NextalignAlignment created, the padded filename, the skip list, the
    merged headers, the segment manifest. The versioned spelling is legitimate only
    in ``meta_data.accession_version``, where GenBankFetcher uses it to notice a
    revised record.

    Each test below pins a place where one side of a comparison was normalised and
    the other was not - a comparison that cannot match, and fails silently rather
    than loudly. Every input the repo ships is bare already, so none of this changes
    a shipped run; it changes what happens the first time a versioned value arrives,
    which an NCBI FASTA download or a hand-built ``--precomputed_ref_dir`` refset
    does by default.
    """

    def test_skip_ids_are_loaded_bare_and_labels_survive(self, tmp_path: Path):
        """``filtered_sequences_ids.txt`` was read verbatim and matched against
        version-stripped record ids, so a versioned line excluded nothing."""
        ids_file = tmp_path / "filtered_sequences_ids.txt"
        ids_file.write_text(
            "NC_001542.1\nKX148218\ncluster.1\nA/swine/Iowa/4.1/1976\n\n",
            encoding="utf-8",
        )

        assert PadAlignment._load_skip_ids(str(ids_file)) == {
            "NC_001542",
            "KX148218",
            "cluster.1",
            "A/swine/Iowa/4.1/1976",
        }

    def test_a_versioned_skip_id_excludes_the_sequence_it_names(self, tmp_path: Path):
        """The failure was *open*: a sequence CollectFilteredSequences excluded was
        padded in anyway and reached the tree and the DB as though it had aligned."""
        guide = tmp_path / "refset_1_aln.fasta"
        _write_fasta(guide, [("NC_001542", "ACGTACGT")])
        _write_fasta(
            tmp_path / "in" / "NC_001542" / "NC_001542.aligned.fasta",
            [("NC_001542", "ACGTACGT"), ("KX148218", "ACGTACGA"), ("MK123456", "ACGTACGC")],
        )
        ids_file = tmp_path / "filtered_sequences_ids.txt"
        ids_file.write_text("KX148218.1\n", encoding="utf-8")

        processor = _guide_processor(tmp_path, guide, skip_ids=str(ids_file))
        processor.process_master_alignment(
            str(guide), str(tmp_path / "in"), str(tmp_path), str(tmp_path / "out"), True
        )

        merged_ids = [record.id for record in _merged_records(tmp_path, guide)]
        assert "KX148218" not in merged_ids
        assert "MK123456" in merged_ids

    def test_a_versioned_backbone_header_is_emitted_bare(self, tmp_path: Path):
        """The backbone row used to keep its original header while every key derived
        from it was canonical, putting both spellings of one genome into the merged
        MSA - one as the backbone row, one as the nextalign reference row."""
        guide = tmp_path / "refset_1_aln.fasta"
        guide.write_text(
            ">NC_001542.1 Rabies lyssavirus, complete genome\nACGTACGT\n", encoding="utf-8"
        )
        _write_fasta(
            tmp_path / "in" / "NC_001542" / "NC_001542.aligned.fasta",
            [("NC_001542", "ACGTACGT"), ("KX148218", "ACGTACGA")],
        )

        processor = _guide_processor(tmp_path, guide)
        processor.process_master_alignment(
            str(guide), str(tmp_path / "in"), str(tmp_path), str(tmp_path / "out"), True
        )

        records = _merged_records(tmp_path, guide)
        assert records[0].id == "NC_001542"
        # Free text after the accession is annotation, not identity - it stays.
        assert records[0].description == "NC_001542 Rabies lyssavirus, complete genome"
        merged_text = (
            tmp_path / "out" / "refset_1_aln_merged_MSA.fasta"
        ).read_text(encoding="utf-8")
        assert "NC_001542.1" not in merged_text

    def test_get_master_list_returns_bare_accessions(self, tmp_path: Path):
        """The other side of this lookup is a directory NextalignAlignment named
        with ``accession_from_filename``, which is always bare."""
        processor = _guide_processor(tmp_path, tmp_path / "unused.fasta")

        assert processor.get_master_list("NC_001542.1,KX148218") == ["NC_001542", "KX148218"]

    def test_a_versioned_master_still_finds_its_reference_aln_directory(self, tmp_path: Path):
        """Missing the directory is not a soft failure: in strict update mode it
        aborts the run, and otherwise the master's whole segment leaves the MSA."""
        nextalign_dir = tmp_path / "Nextalign"
        ref_aln = nextalign_dir / "reference_aln" / "NC_001542"
        ref_aln.mkdir(parents=True)
        _write_fasta(ref_aln / "NC_001542.aligned.fasta", [("NC_001542", "ACGT")])

        processor = _guide_processor(tmp_path, tmp_path / "unused.fasta")
        used = []

        def _capture_process(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False, segment_value=None):
            used.append(Path(reference_alignment_file).name)

        processor.process_master_alignment = types.MethodType(_capture_process, processor)
        processor.process_all_masters(
            master_list=processor.get_master_list("NC_001542.1"),
            nextalign_dir=str(nextalign_dir),
            master_segment_map={},
            precomputed_ref_dir=None,
        )

        assert used == ["NC_001542.aligned.fasta"]

    def test_update_db_master_list_and_segment_map_share_one_key(self, tmp_path: Path):
        """``process_all_masters`` looks the segment up with an entry from the master
        list, so the two have to be keyed the same way or the segment reads None."""
        update_db = tmp_path / "update.db"
        _write_update_alignment_db(
            update_db,
            meta_rows=[("NC_001542.1", "master", "1")],
            alignment_rows=[("NC_001542.1", "NC_001542.1", "ACGT", "1")],
        )

        processor = PadAlignment(
            reference_alignment=None,
            input_dir=str(tmp_path / "in"),
            base_dir=str(tmp_path),
            output_dir=str(tmp_path / "out"),
            keep_intermediate_files=True,
            update_db=str(update_db),
        )

        assert processor.get_master_list("ignored.tsv") == ["NC_001542"]
        assert processor.get_master_segment_map("ignored.tsv") == {"NC_001542": "1"}

    def test_read_fasta_headers_drop_versions_without_truncating_labels(self, tmp_path: Path):
        """``split('.')[0]`` truncated a strain name mid-name, so the orphan warning
        named an accession that does not exist."""
        fasta = tmp_path / "headers.fasta"
        _write_fasta(
            fasta,
            [("NC_001542.1", "ACGT"), ("A/swine/Iowa/4.1/1976", "ACGT"), ("EPI_ISL_402124", "ACGT")],
        )

        assert PadAlignment._read_fasta_headers(str(fasta)) == [
            "NC_001542",
            "A/swine/Iowa/4.1/1976",
            "EPI_ISL_402124",
        ]

    def test_orphan_queries_are_reported_under_their_real_names(self, tmp_path: Path):
        guide = tmp_path / "refset_1_aln.fasta"
        _write_fasta(guide, [("NC_001542", "ACGT")])
        _write_fasta(
            tmp_path / "in" / "REF_ORPHAN" / "REF_ORPHAN.aligned.fasta",
            [("REF_ORPHAN", "ACGT"), ("A/swine/Iowa/4.1/1976", "ACGT"), ("KX148218.1", "ACGT")],
        )

        processor = _guide_processor(tmp_path, guide)
        orphans = processor.find_orphan_query_references(str(guide), str(tmp_path / "in"))

        assert orphans == {"REF_ORPHAN": ["A/swine/Iowa/4.1/1976", "KX148218"]}

    def test_a_versioned_query_directory_is_not_reported_as_an_orphan(self, tmp_path: Path):
        """The directory name and the backbone id name the same genome; comparing the
        two spellings raised a warning telling the operator to investigate sequences
        that were in fact projected normally."""
        guide = tmp_path / "refset_1_aln.fasta"
        _write_fasta(guide, [("NC_001542", "ACGT")])
        _write_fasta(
            tmp_path / "in" / "NC_001542.1" / "NC_001542.1.aligned.fasta",
            [("NC_001542.1", "ACGT"), ("KX148218", "ACGT")],
        )

        processor = _guide_processor(tmp_path, guide)

        assert processor.find_orphan_query_references(str(guide), str(tmp_path / "in")) == {}

    def test_insert_gaps_matches_a_reference_whose_name_contains_a_dot(self, tmp_path: Path):
        """Truncating the label meant no nextalign reference row was found and the
        projection fell through to index matching, which assumes the query and the
        master share coordinates - here they do not, and the bases slide."""
        processor = _guide_processor(tmp_path, tmp_path / "unused.fasta")

        subalignment = [
            SeqRecord(Seq("AA-CGTACGT"), id="A/swine/Iowa/4.1/1976"),
            SeqRecord(Seq("---CGTAC--"), id="Q1"),
        ]
        updated = processor.insert_gaps("ACGT--ACGT", subalignment, "A/swine/Iowa/4.1/1976")

        assert str(updated[1].seq) == "-CGT--AC--"
