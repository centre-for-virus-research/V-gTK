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
REAL_UPDATE_DB = REPO_ROOT / "test_data" / "RABV_test" / "rabv-jul0425.db"


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
        output_dir="pad_out",
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

    updated = processor.insert_gaps(reference_aligned, subalignment)
    assert len(updated) == 1
    assert str(updated[0].seq) == "ACGT--AC--"
    assert len(str(updated[0].seq)) == len(reference_aligned)


def test_process_all_masters_strict_mode_fails_when_db_segment_backbone_missing(tmp_path: Path):
    master, segment = _pick_master_and_segment_from_real_db()

    processor = PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir="pad_out",
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
        output_dir="pad_out",
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
        output_dir="pad_out",
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
        output_dir="pad_out",
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
        output_dir="pad_out",
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
        output_dir="pad_out",
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
        output_dir="pad_out",
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
        output_dir="pad_out",
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
        output_dir="pad_out",
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
        output_dir="pad_out",
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
