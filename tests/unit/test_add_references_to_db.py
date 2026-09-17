"""The implemented parts of the surgical reference-update tool.

scripts/AddReferencesToDb.py is work in progress. These tests cover what it
already promises: a plan accepts only pure additions, the backbone export matches
update mode, and verify-unchanged catches any edit to an original row.
"""

import shutil
import sqlite3
from pathlib import Path

import pytest

import AddReferencesToDb as ARD

REPO_ROOT = Path(__file__).resolve().parents[2]


def _list(tmp_path, text):
    path = tmp_path / "new_refs.tsv"
    path.write_text(text, encoding="utf-8")
    return path


def test_plan_lists_only_the_new_references(tmp_path, basic_update_db):
    ref_list = _list(tmp_path, "REF1\tmaster\t1\nREF2\treference\t2\nNEW_LINEAGE\treference\t1\n")
    rows = ARD.plan(str(basic_update_db), str(ref_list))
    assert [r["primary_accession"] for r in rows] == ["NEW_LINEAGE"]
    assert rows[0]["accession_type"] == "reference" and rows[0]["segment"] == "1"


def test_plan_refuses_a_removed_reference(tmp_path, basic_update_db):
    ref_list = _list(tmp_path, "REF1\tmaster\t1\nNEW_LINEAGE\treference\t1\n")
    with pytest.raises(ARD.PlanError, match="not in the new list"):
        ARD.plan(str(basic_update_db), str(ref_list))


def test_plan_refuses_a_changed_type(tmp_path, basic_update_db):
    ref_list = _list(tmp_path, "REF1\treference\t1\nREF2\tmaster\t1\nNEW\treference\t1\n")
    with pytest.raises(ARD.PlanError, match="REF1: master in the database, reference in the list"):
        ARD.plan(str(basic_update_db), str(ref_list))


def test_plan_refuses_a_new_master(tmp_path, basic_update_db):
    ref_list = _list(tmp_path, "REF1\tmaster\t1\nREF2\treference\t1\nNEW\tmaster\t2\n")
    with pytest.raises(ARD.PlanError, match="new master"):
        ARD.plan(str(basic_update_db), str(ref_list))


def test_plan_with_nothing_new_exits_non_zero(tmp_path, basic_update_db):
    ref_list = _list(tmp_path, "REF1\tmaster\t1\nREF2\treference\t2\n")
    assert ARD.main(["plan", "--db", str(basic_update_db), "--ref_list", str(ref_list),
                     "--output", str(tmp_path / "plan.tsv")]) == 1


def test_export_backbone_holds_the_database_references(tmp_path, basic_update_db):
    out = ARD.export_backbone(str(basic_update_db), str(tmp_path / "backbone"))
    headers = {line[1:].strip() for f in Path(out).glob("refset_*_aln.fasta")
               for line in f.read_text().splitlines() if line.startswith(">")}
    assert {"REF1", "REF2"} <= headers


def test_verify_unchanged_accepts_pure_additions(tmp_path, basic_update_db):
    copy = tmp_path / "copy.db"
    shutil.copyfile(basic_update_db, copy)
    conn = sqlite3.connect(str(copy))
    columns = [r[1] for r in conn.execute("PRAGMA table_info(meta_data)")]
    conn.execute(f"INSERT INTO meta_data (primary_accession, accession_type) VALUES ('NEW', 'reference')")
    conn.execute(f"CREATE TABLE {ARD.LOG_TABLE} (primary_accession TEXT)")
    conn.commit()
    conn.close()
    assert ARD.verify_unchanged(str(basic_update_db), str(copy)) == []


def test_verify_unchanged_catches_an_edited_row(tmp_path, basic_update_db):
    copy = tmp_path / "copy.db"
    shutil.copyfile(basic_update_db, copy)
    conn = sqlite3.connect(str(copy))
    conn.execute("UPDATE meta_data SET accession_type = 'query' WHERE primary_accession = 'REF2'")
    conn.commit()
    conn.close()
    problems = ARD.verify_unchanged(str(basic_update_db), str(copy))
    assert problems and problems[0].startswith("meta_data: 1 original row")


def test_verify_unchanged_catches_a_new_column(tmp_path, basic_update_db):
    copy = tmp_path / "copy.db"
    shutil.copyfile(basic_update_db, copy)
    conn = sqlite3.connect(str(copy))
    conn.execute("ALTER TABLE meta_data ADD COLUMN added TEXT")
    conn.commit()
    conn.close()
    assert any("columns changed" in p for p in ARD.verify_unchanged(str(basic_update_db), str(copy)))


def test_insert_never_writes_over_the_original(tmp_path, basic_update_db):
    with pytest.raises(ValueError, match="original"):
        ARD.insert(str(basic_update_db), str(basic_update_db), [], "unused.fasta")


def test_entry_point_is_separate_and_marked_work_in_progress():
    text = (REPO_ROOT / "vgtk-reference-update.nf").read_text()
    assert "NOT EXPECTED TO WORK YET" in text
    assert "AddReferencesToDb.py" not in (REPO_ROOT / "vgtk-init.nf").read_text()


def test_plan_refuses_a_moved_segment(tmp_path, basic_update_db):
    ref_list = _list(tmp_path, "REF1\tmaster\t1\nREF2\treference\t1\nNEW\treference\t1\n")
    with pytest.raises(ARD.PlanError, match="REF2: segment 2 in the database, 1 in the list"):
        ARD.plan(str(basic_update_db), str(ref_list))
