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
