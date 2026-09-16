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


# --------------------------------------------------------------------------
# References must be tips in the stored UShER tree
# --------------------------------------------------------------------------

import sqlite3

from ValidateRefListAgainstDb import newick_tips

REPO_ROOT = Path(__file__).resolve().parents[2]
RABV_UPDATE_DB = REPO_ROOT / "test_data" / "RABV_test" / "rabv-7Jul26.db"
RABV_TEST_REFS = REPO_ROOT / "test_data" / "rabv_test_ref_list.txt"


def _store_trees(db: Path, *rows):
    """Replace the fixture's stored trees (it ships one UShER tree holding REF1 and REF2)."""
    conn = sqlite3.connect(str(db))
    conn.execute("CREATE TABLE IF NOT EXISTS trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT)")
    conn.execute("DELETE FROM trees")
    conn.executemany("INSERT INTO trees (name, source, newick) VALUES (?, ?, ?)", rows)
    conn.commit()
    conn.close()


def test_fixture_tree_holds_the_fixture_references(tmp_path: Path, basic_update_db: Path, capsys):
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\nREF2\treference\n", encoding="utf-8")
    validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))
    assert "tips in the stored UShER tree" in capsys.readouterr().out


def test_db_without_usher_tree_is_not_checked(tmp_path: Path, basic_update_db: Path, capsys):
    _store_trees(basic_update_db)
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\nREF2\treference\n", encoding="utf-8")
    validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))
    assert "No UShER tree" in capsys.readouterr().out


def test_newick_tips_reads_labels_and_drops_versions():
    assert newick_tips("((NC_001542.1:0.1,'Q 2':0.2)0.9:0.3,REF2:0.4);") == {"NC_001542", "Q 2", "REF2"}


def test_reference_present_in_usher_tree_passes(tmp_path: Path, basic_update_db: Path):
    _store_trees(basic_update_db, ("usher", "usher", "((REF1:1,REF2:1):1,Q1:1);"))
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\nREF2\treference\n", encoding="utf-8")
    validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))


def test_reference_missing_from_usher_tree_fails(tmp_path: Path, basic_update_db: Path):
    _store_trees(basic_update_db, ("usher", "usher", "(REF1:1,Q1:1);"))
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\nREF2\treference\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not tips in the update DB's UShER tree.*REF2"):
        validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))


def test_iqtree_representatives_are_not_required_to_hold_references(tmp_path: Path, basic_update_db: Path):
    _store_trees(basic_update_db, ("iqtree", "iqtree", "(REF1:1,Q1:1);"),
                 ("usher", "usher", "((REF1:1,REF2:1):1,Q1:1);"))
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\nREF2\treference\n", encoding="utf-8")
    validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))


def test_segmented_references_may_sit_in_different_segment_trees(tmp_path: Path, basic_update_db: Path):
    _store_trees(basic_update_db, ("usher_seg1", "usher", "(REF1:1,Q1:1);"),
                 ("usher_seg2", "usher", "(REF2:1,Q2:1);"))
    ref = tmp_path / "refs.tsv"
    ref.write_text("REF1\tmaster\t1\nREF2\treference\t2\n", encoding="utf-8")
    validate_ref_main(_Args(ref_list=str(ref), db=str(basic_update_db)))


@pytest.mark.skipif(not RABV_UPDATE_DB.exists(), reason="RABV update test DB not present")
def test_update_test_profile_references_are_all_in_its_usher_tree(capsys):
    validate_ref_main(_Args(ref_list=str(RABV_TEST_REFS), db=str(RABV_UPDATE_DB)))
    assert "tips in the stored UShER tree" in capsys.readouterr().out
