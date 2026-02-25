from pathlib import Path

import pytest

from CheckAlignmentLengths import check_uniform_lengths


def write_fasta(path: Path, records):
    with path.open("w", encoding="utf-8") as handle:
        for rec_id, seq in records:
            handle.write(f">{rec_id}\n{seq}\n")


def test_check_uniform_lengths_pass(tmp_path: Path):
    fasta = tmp_path / "uniform.fasta"
    write_fasta(fasta, [("a", "ACGT"), ("b", "TGCA"), ("c", "NNNN")])

    result = check_uniform_lengths(str(fasta))
    assert result["is_uniform"] is True
    assert result["count"] == 3
    assert result["min_len"] == 4
    assert result["max_len"] == 4
    assert result["unique_lengths"] == [4]
    assert result["non_majority_ids"] == []


def test_check_uniform_lengths_fail(tmp_path: Path):
    fasta = tmp_path / "mixed.fasta"
    write_fasta(fasta, [("a", "ACGT"), ("b", "TGCAA"), ("c", "NNNN")])

    result = check_uniform_lengths(str(fasta))
    assert result["is_uniform"] is False
    assert result["count"] == 3
    assert result["min_len"] == 4
    assert result["max_len"] == 5
    assert result["unique_lengths"] == [4, 5]
    assert result["non_majority_ids"] == ["b"]


def test_check_uniform_lengths_empty_raises(tmp_path: Path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="No FASTA records found"):
        check_uniform_lengths(str(fasta))
