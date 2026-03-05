import shutil
from pathlib import Path
import types

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from PadAlignment import PadAlignment


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "pad_alignment"


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
