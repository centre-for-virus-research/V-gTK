import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "ValidateDbTree.py"


def _add_required_alignment_tables(db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS insertions (primary_accession TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS host_taxa (primary_accession TEXT)")
        cur.execute("INSERT INTO insertions(primary_accession) VALUES ('REF1')")
        cur.execute("INSERT INTO host_taxa(primary_accession) VALUES ('REF1')")
        conn.commit()
    finally:
        conn.close()


def test_validate_db_tree_update_integrity_passes_on_seed_db(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    outdir = tmp_path / "out"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(outdir),
            "--check-update-integrity",
                "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "Update integrity checks:" in report
    assert "last_batch_id:" in report


def test_validate_db_tree_update_integrity_fails_without_audit_tables(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DROP TABLE update_batches")
        conn.execute("DROP TABLE update_table_deltas")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "update audit tables are missing" in result.stderr


def test_validate_db_tree_update_integrity_detects_segment_contamination(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("UPDATE features SET segment='9' WHERE accession='Q_OLD'")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "segment contamination detected" in result.stderr
