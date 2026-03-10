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
        cur.execute("CREATE TABLE IF NOT EXISTS host_taxa (taxonomy_id TEXT)")
        cur.execute("INSERT INTO insertions(primary_accession) VALUES ('REF1')")
        cur.execute("INSERT INTO host_taxa(taxonomy_id) VALUES ('111')")
        cur.execute("INSERT INTO host_taxa(taxonomy_id) VALUES ('222')")
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
    assert "[sequence_alignment vs meta_data]" in report
    assert "[features vs meta_data]" in report
    assert "[host_taxa vs meta_data]" in report
    assert "Update integrity checks:" in report
    assert "last_batch_id:" in report


def test_validate_db_tree_fails_when_features_do_not_match_expected_accessions(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DELETE FROM features WHERE accession='REF2'")
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
    assert "DB consistency checks failed for features vs meta_data" in result.stderr


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


def test_validate_db_tree_update_integrity_ignores_placeholder_cluster_assignments_for_excluded_rows(tmp_path: Path):
    db_path = tmp_path / "placeholder_cluster.db"
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE meta_data (primary_accession TEXT, cluster_98pct TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, cluster_98pct, segment) VALUES (?, ?, ?)",
            [
                ("REF1", "REF1", ""),
                ("Q_PLACEABLE", "REF1", ""),
                ("Q_UNPLACED", "NA- see tree", ""),
            ],
        )
        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.executemany(
            "INSERT INTO sequences(header, sequence) VALUES (?, ?)",
            [("REF1", "ATGC"), ("Q_PLACEABLE", "ATGT"), ("Q_UNPLACED", "ATGA")],
        )
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, segment) VALUES (?, ?, ?, ?)",
            [("REF1", "REF1", "REF1", ""), ("Q_PLACEABLE", "Q_PLACEABLE", "REF1", "")],
        )
        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
        cur.execute(
            "INSERT INTO trees(name, source, newick, segment, segment_key) VALUES (?, ?, ?, ?, ?)",
            ("usher", "usher", "(REF1:0.1,Q_PLACEABLE:0.2);", "", None),
        )
        cur.execute("CREATE TABLE excluded_accessions (primary_accession TEXT, reason TEXT)")
        cur.execute(
            "INSERT INTO excluded_accessions(primary_accession, reason) VALUES (?, ?)",
            ("Q_UNPLACED", "alignment_filtering"),
        )
        cur.execute("CREATE TABLE update_batches (batch_id TEXT, mode TEXT, started_at TEXT, finished_at TEXT)")
        cur.execute("INSERT INTO update_batches(batch_id, mode, started_at, finished_at) VALUES ('b1', 'update', '2026-01-01', '2026-01-01')")
        cur.execute("CREATE TABLE update_table_deltas (batch_id TEXT, table_name TEXT, before_count INTEGER, after_count INTEGER, delta INTEGER)")
        cur.execute("INSERT INTO update_table_deltas(batch_id, table_name, before_count, after_count, delta) VALUES ('b1', 'meta_data', 2, 3, 1)")
        cur.execute("CREATE TABLE features (accession TEXT, segment TEXT)")
        cur.executemany("INSERT INTO features(accession, segment) VALUES (?, ?)", [("REF1", ""), ("Q_PLACEABLE", "")])
        cur.execute("CREATE TABLE insertions (primary_accession TEXT)")
        cur.execute("CREATE TABLE host_taxa (primary_accession TEXT)")
        conn.commit()
    finally:
        conn.close()

    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(outdir),
            "--check-update-integrity",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "Missing centroids in tree: 0" in report
    assert "Missing in tree (meta_data -> tree): 0" in report
