import csv
from pathlib import Path
import shutil
import subprocess

import pytest

from BlastAlignment import BlastAlignment


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "blast_alignment"


def read_tsv_as_dicts(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def make_processor(tmp_path: Path, gb_matrix_path: Path, is_segmented: str = "N"):
    return BlastAlignment(
        query_fasta=str(DATA_DIR / "query.fa"),
        db_fasta=str(DATA_DIR / "ref.fa"),
        base_dir=str(tmp_path),
        output_dir="blast_out",
        output_file="query_tophits.tsv",
        is_segmented_virus=is_segmented,
        master_acc=str(DATA_DIR / "master_refs.tsv"),
        is_update="N",
        keep_blast_tmp_dir="N",
        gb_matrix=str(gb_matrix_path),
        segment_file=str(DATA_DIR / "master_refs.tsv"),
    )


def test_get_exclusion_list_refs(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    refs = processor.get_exclusion_list_refs()
    assert refs == {"REF_EXCL"}


def test_get_master_list_from_file(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    assert processor.get_master_list() == ["REF_MASTER"]


def test_get_master_list_from_comma_string(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    processor.master_acc = "REF_MASTER,REF_OTHER"
    assert processor.get_master_list() == ["REF_MASTER", "REF_OTHER"]


def test_ref_segments_static():
    segment_map = BlastAlignment.ref_segments(str(DATA_DIR / "query_uniq_tophit_annotated.tsv"))
    assert segment_map == {
        "1": {"REF_MASTER": ["Q1"]},
        "2": {"REF_OTHER": ["Q4"]},
    }


def test_write_filtered_ref_fasta_excludes_exclusion_list(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    processor.write_filtered_ref_fasta(str(out_dir), {"REF_EXCL"})

    actual = (out_dir / "ref_seq_filtered.fa").read_text(encoding="utf-8")
    expected = (DATA_DIR / "expected_ref_seq_filtered.fa").read_text(encoding="utf-8")
    assert actual == expected


def test_update_gb_matrix_marks_no_hit_and_exclusion_refs(tmp_path: Path):
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    out_dir = tmp_path / "blast"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "query_tophits.tsv").write_text(
        (DATA_DIR / "query_tophits.tsv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    (out_dir / "query_uniq_tophits.tsv").write_text(
        (DATA_DIR / "query_uniq_tophits.tsv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    processor = make_processor(tmp_path, gb_copy)
    processor.update_gB_matrix(
        query_fasta=str(DATA_DIR / "query.fa"),
        query_tophit_uniq=str(out_dir / "query_uniq_tophits.tsv"),
        gB_matrix_file=str(gb_copy),
    )

    actual = read_tsv_as_dicts(gb_copy)
    expected = read_tsv_as_dicts(DATA_DIR / "expected_gb_matrix_after_update.tsv")
    assert actual == expected


def test_process_non_segmented_orchestration_without_external_tools(tmp_path: Path, monkeypatch):
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    processor = make_processor(tmp_path, gb_copy, is_segmented="N")
    out_dir = tmp_path / "blast_out"
    out_dir.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(processor, "check_blast_exists", lambda _cmd: True)
    monkeypatch.setattr(processor, "check_tool_exists", lambda _cmd: True)
    monkeypatch.setattr(processor, "run_makeblastdb", lambda _tmp: None)

    def fake_run_blastn(output_dir, _query_file):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        (Path(output_dir) / "query_tophits.tsv").write_text(
            (DATA_DIR / "query_tophits.tsv").read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    def fake_process_non_segmented_virus(output_dir, _query_fasta):
        (Path(output_dir) / "query_uniq_tophits.tsv").write_text(
            (DATA_DIR / "query_uniq_tophits.tsv").read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    monkeypatch.setattr(processor, "run_blastn", fake_run_blastn)
    monkeypatch.setattr(processor, "process_non_segmented_virus", fake_process_non_segmented_virus)

    processor.process()

    actual_matrix = read_tsv_as_dicts(gb_copy)
    expected_matrix = read_tsv_as_dicts(DATA_DIR / "expected_gb_matrix_after_update.tsv")
    assert actual_matrix == expected_matrix

    filtered_fa = out_dir / "ref_seq_filtered.fa"
    assert filtered_fa.exists()
    assert "REF_EXCL" not in filtered_fa.read_text(encoding="utf-8")


def test_update_gb_matrix_raises_on_malformed_unique_hits_row(tmp_path: Path):
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    out_dir = tmp_path / "blast"
    out_dir.mkdir(parents=True, exist_ok=True)
    bad_uniq = out_dir / "query_uniq_tophits.tsv"
    bad_uniq.write_text("Q1\tREF_MASTER\t99.9\n", encoding="utf-8")

    processor = make_processor(tmp_path, gb_copy)
    with pytest.raises(ValueError, match="Malformed unique BLAST hit row"):
        processor.update_gB_matrix(
            query_fasta=str(DATA_DIR / "query.fa"),
            query_tophit_uniq=str(bad_uniq),
            gB_matrix_file=str(gb_copy),
        )


def test_process_raises_when_gb_matrix_missing_gi_number(tmp_path: Path):
    gb_bad = tmp_path / "gb_bad.tsv"
    gb_bad.write_text("primary_accession\tsequence\nQ1\tATGC\n", encoding="utf-8")

    processor = make_processor(tmp_path, gb_bad)
    with pytest.raises(ValueError, match="missing required columns: gi_number"):
        processor.process()


def test_process_raises_when_reference_fasta_empty(tmp_path: Path):
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    empty_ref = tmp_path / "empty_ref.fa"
    empty_ref.write_text("", encoding="utf-8")

    processor = make_processor(tmp_path, gb_copy)
    processor.db_fasta = str(empty_ref)
    processor.master_acc = ""

    with pytest.raises(ValueError, match="Reference FASTA has no sequences"):
        processor.process()


def test_process_raises_when_master_missing_in_reference_fasta(tmp_path: Path):
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    processor = make_processor(tmp_path, gb_copy)
    processor.master_acc = "NC_001542"

    with pytest.raises(ValueError, match="Master reference accession\(s\) not found"):
        processor.process()


def test_run_makeblastdb_raises_on_subprocess_failure(tmp_path: Path, monkeypatch):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")

    def fail_run(*args, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=args[0])

    monkeypatch.setattr(subprocess, "run", fail_run)
    with pytest.raises(RuntimeError, match="makeblastdb failed"):
        processor.run_makeblastdb(str(tmp_path))


def test_hydrate_update_reference_assets_uses_db_as_reference_source(tmp_path: Path, basic_update_db: Path):
    query_fa = tmp_path / "query.fa"
    query_fa.write_text(">Q_NEW\nATGC\n", encoding="utf-8")
    gb_matrix = tmp_path / "gb_matrix.tsv"
    gb_matrix.write_text(
        "division\tgi_number\tsequence\nVRL\tQ_NEW\tATGC\n",
        encoding="utf-8",
    )

    processor = BlastAlignment(
        query_fasta=str(query_fa),
        db_fasta=str(tmp_path / "unused.fa"),
        base_dir=str(tmp_path),
        output_dir="blast_out",
        output_file="query_tophits.tsv",
        is_segmented_virus="Y",
        master_acc=None,
        is_update="N",
        keep_blast_tmp_dir="N",
        gb_matrix=str(gb_matrix),
        segment_file=None,
        update_db=str(basic_update_db),
    )

    processor.hydrate_update_reference_assets()

    ref_fasta = Path(processor.db_fasta)
    assert ref_fasta.exists()
    ref_text = ref_fasta.read_text(encoding="utf-8")
    assert ">REF1" in ref_text
    assert ">REF2" in ref_text
    assert processor.segment_file is not None
    assert Path(processor.segment_file).exists()
    assert processor.master_acc == processor.segment_file
    ref_list_text = Path(processor.segment_file).read_text(encoding="utf-8")
    assert "REF1\tmaster\t1" in ref_list_text
    assert "REF2\treference\t2" in ref_list_text


# --- seqkit reverse-complement: failures must not pass silently -------------
# Regression cover for minus-strand queries vanishing from a run. The call went
# through os.system(" ".join(cmd)) inside an "except subprocess.CalledProcessError"
# block; os.system cannot raise that, so a seqkit failure printed success, left
# an empty redirect target, and dropped every reverse-complement submission.


def test_reverse_complement_raises_when_seqkit_exits_nonzero(tmp_path: Path, monkeypatch):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    src = tmp_path / "minus.fa"
    src.write_text(">Q1\nAAAA\n", encoding="utf-8")

    def fail_run(*args, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=args[0], stderr=b"seqkit: boom")

    monkeypatch.setattr(subprocess, "run", fail_run)
    with pytest.raises(RuntimeError, match="seqkit failed to reverse-complement"):
        processor._reverse_complement_fasta(str(src), str(tmp_path / "out.fa"))


def test_reverse_complement_raises_when_seqkit_not_on_path(tmp_path: Path, monkeypatch):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    src = tmp_path / "minus.fa"
    src.write_text(">Q1\nAAAA\n", encoding="utf-8")

    def missing_tool(*args, **kwargs):
        raise FileNotFoundError("seqkit")

    monkeypatch.setattr(subprocess, "run", missing_tool)
    with pytest.raises(RuntimeError, match="seqkit is not on PATH"):
        processor._reverse_complement_fasta(str(src), str(tmp_path / "out.fa"))


def test_reverse_complement_raises_when_sequences_are_dropped(tmp_path: Path, monkeypatch):
    """seqkit exiting 0 after writing short output is still silent data loss."""
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    src = tmp_path / "minus.fa"
    src.write_text(">Q1\nAAAA\n>Q2\nCCCC\n>Q3\nGGGG\n", encoding="utf-8")
    out = tmp_path / "out.fa"

    def truncating_run(cmd, check=False, stdout=None, stderr=None, **kwargs):
        stdout.write(b">Q1\nTTTT\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(subprocess, "run", truncating_run)
    with pytest.raises(RuntimeError, match=r"lost sequences .*3 in, 1 out"):
        processor._reverse_complement_fasta(str(src), str(out))


@pytest.mark.skipif(shutil.which("seqkit") is None, reason="seqkit not installed")
def test_reverse_complement_roundtrip_with_real_seqkit(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    src = tmp_path / "minus.fa"
    src.write_text(">Q1\nAAAAGGGG\n", encoding="utf-8")
    out = tmp_path / "out.fa"

    processor._reverse_complement_fasta(str(src), str(out))

    text = out.read_text(encoding="utf-8")
    assert ">Q1" in text
    assert "CCCCTTTT" in text.upper()


# --- concatenation ----------------------------------------------------------


def test_concat_fasta_files_separates_records_without_trailing_newline(tmp_path: Path):
    first = tmp_path / "a.fa"
    second = tmp_path / "b.fa"
    first.write_text(">Q1\nAAAA", encoding="utf-8")  # no trailing newline
    second.write_text(">Q2\nCCCC\n", encoding="utf-8")
    out = tmp_path / "all.fa"

    BlastAlignment._concat_fasta_files([str(first), str(second)], str(out))

    assert out.read_text(encoding="utf-8") == ">Q1\nAAAA\n>Q2\nCCCC\n"


def test_concat_fasta_files_handles_paths_containing_spaces(tmp_path: Path):
    spaced_dir = tmp_path / "dir with space"
    spaced_dir.mkdir()
    src = spaced_dir / "seg 1.fa"
    src.write_text(">Q1\nAAAA\n", encoding="utf-8")
    out = tmp_path / "all.fa"

    BlastAlignment._concat_fasta_files([str(src)], str(out))

    assert out.read_text(encoding="utf-8") == ">Q1\nAAAA\n"


# --- missing queries become real exclusions ---------------------------------


def test_update_gb_matrix_marks_missing_query_sequences_as_excluded(tmp_path: Path):
    """A query with a hit but no sequence must not stay an included row.

    These were previously skipped with a [warn] only, leaving DB rows with no
    alignment and no tree tip.
    """
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    out_dir = tmp_path / "blast"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "query_uniq_tophits.tsv").write_text(
        (DATA_DIR / "query_uniq_tophits.tsv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    processor = make_processor(tmp_path, gb_copy)
    # Q1 has a good plus-strand hit, so it would otherwise be included.
    processor.missing_query_accs.add("Q1")
    processor.update_gB_matrix(
        query_fasta=str(DATA_DIR / "query.fa"),
        query_tophit_uniq=str(out_dir / "query_uniq_tophits.tsv"),
        gB_matrix_file=str(gb_copy),
    )

    rows = {row["gi_number"]: row for row in read_tsv_as_dicts(gb_copy)}
    assert rows["Q1"]["exclusion_status"] == "1"
    assert rows["Q1"]["exclusion_criteria"] == "excluded: no sequence in strand-normalised query FASTA"


def test_missing_query_accs_starts_empty(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    assert processor.missing_query_accs == set()
