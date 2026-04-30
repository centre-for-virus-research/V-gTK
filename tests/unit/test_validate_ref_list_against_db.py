from pathlib import Path

import pytest

from ValidateRefListAgainstDb import main as validate_ref_main


class _Args:
    def __init__(self, ref_list: str, db: str):
        self.ref_list = ref_list
        self.db = db


def test_validate_ref_list_against_db_accepts_matching_refs(tmp_path: Path, basic_update_db: Path):
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\nREF2\treference\n", encoding="utf-8")

    validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))


def test_validate_ref_list_against_db_fails_on_missing_reference(tmp_path: Path, basic_update_db: Path):
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\nREF_MISSING\treference\n", encoding="utf-8")

    with pytest.raises(ValueError, match="not present in DB"):
        validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))


def test_validate_ref_list_against_db_ignores_non_reference_rows(tmp_path: Path, basic_update_db: Path):
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\nSOMETHING\texclusion_list\n", encoding="utf-8")

    validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))


def test_validate_ref_list_against_db_accepts_headered_refs(tmp_path: Path, basic_update_db: Path):
    ref = tmp_path / "refs.tsv"
    ref.write_text(
        "primary_accession\tstatus\tsegment\nREF1\tmaster\t1\nREF2\treference\t1\n",
        encoding="utf-8",
    )

    validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))
