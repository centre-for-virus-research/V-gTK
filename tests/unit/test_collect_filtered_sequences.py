import csv
import subprocess
import sys
from pathlib import Path

import pytest

from CollectFilteredSequences import (
    collect_filtered_sequences,
    write_filtered_ids_only,
    write_filtered_list,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "CollectFilteredSequences.py"
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "collect_filtered_sequences_edge"


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["seqName", "errors", "warnings"])
        writer.writerows(rows)


def test_collect_filtered_sequences_only_keeps_query_alignment_rows_with_errors(tmp_path: Path):
    nextalign_dir = tmp_path / "Nextalign"

    write_csv(
        nextalign_dir / "query_aln" / "REF_A" / "chunk.errors.csv",
        [
            ["SEQ_1", "frame shift", ""],
            ["SEQ_2", "", "minor issue"],
        ],
    )
    write_csv(
        nextalign_dir / "reference_aln" / "REF_B" / "part.errors.csv",
        [["SEQ_3", "too many Ns", "warning text"]],
    )

    filtered = collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "out.tsv"))

    assert sorted(filtered.keys()) == ["SEQ_1"]
    assert filtered["SEQ_1"]["reference"] == "REF_A"
    assert "SEQ_3" not in filtered


def test_write_outputs(tmp_path: Path):
    filtered = {
        "SEQ_B": {"reference": "R2", "error": "error B", "warnings": ""},
        "SEQ_A": {"reference": "R1", "error": "error A", "warnings": "warn"},
    }

    out_tsv = tmp_path / "filtered_sequences.tsv"
    write_filtered_list(filtered, str(out_tsv))
    ids_file = write_filtered_ids_only(filtered, str(out_tsv))

    lines = out_tsv.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "seq_name\treference\terror\twarnings"
    # sorted by seq_name
    assert lines[1].startswith("SEQ_A\t")
    assert lines[2].startswith("SEQ_B\t")

    assert Path(ids_file).read_text(encoding="utf-8").splitlines() == ["SEQ_A", "SEQ_B"]


def test_cli_creates_empty_outputs_when_no_errors(tmp_path: Path):
    nextalign_dir = tmp_path / "Nextalign"
    write_csv(
        nextalign_dir / "query_aln" / "REF_X" / "empty.errors.csv",
        [["SEQ_OK", "", "only warning"]],
    )

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-n",
            str(nextalign_dir),
            "-o",
            "filtered.tsv",
            "-b",
            str(tmp_path),
        ],
        check=True,
    )

    out_tsv = tmp_path / "filtered.tsv"
    out_ids = tmp_path / "filtered_ids.txt"

    assert out_tsv.exists()
    assert out_ids.exists()
    assert out_tsv.read_text(encoding="utf-8").strip() == "seq_name\treference\terror\twarnings"
    assert out_ids.read_text(encoding="utf-8") == ""


def test_collect_filtered_sequences_with_fixture_dataset(tmp_path: Path):
    nextalign_dir = DATA_DIR / "Nextalign"
    out_tsv = tmp_path / "fixture_filtered.tsv"

    filtered = collect_filtered_sequences(str(nextalign_dir), str(out_tsv))
    write_filtered_list(filtered, str(out_tsv))
    out_ids = write_filtered_ids_only(filtered, str(out_tsv))

    assert out_tsv.read_text(encoding="utf-8") == (DATA_DIR / "expected_filtered.tsv").read_text(encoding="utf-8")
    assert Path(out_ids).read_text(encoding="utf-8") == (DATA_DIR / "expected_ids.txt").read_text(encoding="utf-8")


def test_collect_filtered_sequences_marks_unprojectable_queries(tmp_path: Path):
    nextalign_dir = tmp_path / "Nextalign"

    # reference_aln contains only REF_OK as projectable into master coordinates
    ref_dir = nextalign_dir / "reference_aln" / "MASTER1"
    ref_dir.mkdir(parents=True, exist_ok=True)
    (ref_dir / "MASTER1.aligned.fasta").write_text(
        ">MASTER1\nATGC\n>REF_OK\nATGC\n",
        encoding="utf-8",
    )

    # query_aln for REF_OK should remain unflagged
    q_ok_dir = nextalign_dir / "query_aln" / "REF_OK"
    q_ok_dir.mkdir(parents=True, exist_ok=True)
    (q_ok_dir / "REF_OK.aligned.fasta").write_text(
        ">REF_OK\nATGC\n>Q_OK\nATGC\n",
        encoding="utf-8",
    )

    # query_aln for REF_ORPHAN cannot be projected because REF_ORPHAN absent in reference_aln
    q_orphan_dir = nextalign_dir / "query_aln" / "REF_ORPHAN"
    q_orphan_dir.mkdir(parents=True, exist_ok=True)
    (q_orphan_dir / "REF_ORPHAN.aligned.fasta").write_text(
        ">REF_ORPHAN\nATGC\n>Q_ORPHAN\nATGC\n",
        encoding="utf-8",
    )

    filtered = collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "filtered.tsv"))

    assert "Q_ORPHAN" in filtered
    assert filtered["Q_ORPHAN"]["reference"] == "REF_ORPHAN"
    assert "cannot be projected" in filtered["Q_ORPHAN"]["error"]
    assert "Q_OK" not in filtered


def test_collect_filtered_sequences_raises_when_nextalign_dir_missing(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="Nextalign directory not found"):
        collect_filtered_sequences(str(tmp_path / "missing"), str(tmp_path / "out.tsv"))


def test_collect_filtered_sequences_raises_on_malformed_errors_csv(tmp_path: Path):
    nextalign_dir = tmp_path / "Nextalign"
    bad_csv = nextalign_dir / "query_aln" / "REF_A" / "bad.errors.csv"
    bad_csv.parent.mkdir(parents=True, exist_ok=True)
    bad_csv.write_text("seq\terr\nQ1\toops\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Malformed nextalign errors file"):
        collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "out.tsv"))
