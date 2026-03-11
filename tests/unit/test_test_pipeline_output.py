import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "TestPipelineOutput.py"


def _write_tsv(path: Path, header: str, rows: List[str]) -> None:
    body = "\n".join([header] + rows) + "\n"
    path.write_text(body, encoding="utf-8")


def _create_valid_db(path: Path) -> None:
    con = sqlite3.connect(path)
    cur = con.cursor()
    cur.execute("CREATE TABLE meta_data (primary_accession TEXT, segment TEXT)")
    cur.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, alignment_name TEXT)")
    cur.execute("CREATE TABLE excluded_accessions (primary_accession TEXT, reason TEXT)")
    cur.execute("INSERT INTO meta_data VALUES ('A1', '1')")
    cur.execute("INSERT INTO sequence_alignment VALUES ('A1', 'REF1')")
    cur.execute("INSERT INTO meta_data VALUES ('REF1', '1')")
    con.commit()
    con.close()


def _run_non_segmented(tmp_path: Path, sqlite_db: Path, blast_hits_content: Optional[str] = None):
    blast_hits = tmp_path / "blast.tsv"
    if blast_hits_content is None:
        blast_hits_content = "q\tr\tscore\tstrand\nq1\tr1\t100\t+\n"
    blast_hits.write_text(blast_hits_content, encoding="utf-8")

    gb_matrix = tmp_path / "gb.tsv"
    _write_tsv(gb_matrix, "primary_accession\thost", ["A1\thuman"])

    output = tmp_path / "out.txt"
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--mode",
            "non_segmented",
            "--blast_hits",
            str(blast_hits),
            "--gb_matrix",
            str(gb_matrix),
            "--sqlite_db",
            str(sqlite_db),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
    )


def test_non_segmented_missing_db_file_fails_informatively(tmp_path: Path):
    missing_db = tmp_path / "does_not_exist.db"
    result = _run_non_segmented(tmp_path, missing_db)

    assert result.returncode == 2
    assert "input validation error" in result.stderr.lower()
    assert "sqlite_db does not exist" in result.stderr
    assert "traceback" not in result.stderr.lower()


def test_non_segmented_missing_required_table_fails_informatively(tmp_path: Path):
    db = tmp_path / "bad.db"
    con = sqlite3.connect(db)
    cur = con.cursor()
    cur.execute("CREATE TABLE meta_data (primary_accession TEXT)")
    con.commit()
    con.close()

    result = _run_non_segmented(tmp_path, db)

    assert result.returncode == 2
    assert "missing required table" in result.stderr.lower()
    assert "sequence_alignment" in result.stderr
    assert "excluded_accessions" in result.stderr
    assert "traceback" not in result.stderr.lower()


def test_non_segmented_missing_required_column_fails_informatively(tmp_path: Path):
    db = tmp_path / "bad_columns.db"
    con = sqlite3.connect(db)
    cur = con.cursor()
    cur.execute("CREATE TABLE meta_data (primary_accession TEXT)")
    cur.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, alignment_name TEXT)")
    cur.execute("CREATE TABLE excluded_accessions (accession TEXT, reason TEXT)")
    con.commit()
    con.close()

    result = _run_non_segmented(tmp_path, db)

    assert result.returncode == 2
    assert "missing required column" in result.stderr.lower()
    assert "excluded_accessions: primary_accession" in result.stderr
    assert "traceback" not in result.stderr.lower()


def test_non_segmented_bad_blast_format_fails_cleanly(tmp_path: Path):
    db = tmp_path / "ok.db"
    _create_valid_db(db)

    result = _run_non_segmented(tmp_path, db, blast_hits_content="q,r,score,strand\nq1,r1,100,+\n")

    assert result.returncode == 1
    assert "BLAST file has 1 columns, expected 4" in result.stdout
    assert "traceback" not in result.stderr.lower()


def test_non_segmented_empty_blast_file_warns_but_does_not_fail(tmp_path: Path):
    db = tmp_path / "ok.db"
    con = sqlite3.connect(db)
    cur = con.cursor()
    cur.execute("CREATE TABLE meta_data (primary_accession TEXT, segment TEXT)")
    cur.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, alignment_name TEXT)")
    cur.execute("CREATE TABLE excluded_accessions (primary_accession TEXT, reason TEXT)")
    cur.execute("INSERT INTO meta_data VALUES ('A1', '1')")
    cur.execute("INSERT INTO meta_data VALUES ('REF1', '1')")
    cur.execute("INSERT INTO sequence_alignment VALUES ('A1', 'REF1')")
    cur.execute("INSERT INTO sequence_alignment VALUES ('REF1', 'REF1')")
    con.commit()
    con.close()

    gb_matrix = tmp_path / "gb.tsv"
    _write_tsv(gb_matrix, "primary_accession\thost\texclusion_status", ["A1\thuman\t1"])

    output = tmp_path / "out.txt"
    blast_hits = tmp_path / "blast.tsv"
    blast_hits.write_text("", encoding="utf-8")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--mode",
            "non_segmented",
            "--blast_hits",
            str(blast_hits),
            "--gb_matrix",
            str(gb_matrix),
            "--sqlite_db",
            str(db),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "all matrix records are excluded" in result.stdout
    assert "traceback" not in result.stderr.lower()


def test_non_segmented_missing_alignment_reports_examples_and_duplicates(tmp_path: Path):
    db = tmp_path / "diagnostic.db"
    con = sqlite3.connect(db)
    cur = con.cursor()
    cur.execute("CREATE TABLE meta_data (primary_accession TEXT, segment TEXT)")
    cur.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, alignment_name TEXT)")
    cur.execute("CREATE TABLE excluded_accessions (primary_accession TEXT, reason TEXT)")
    cur.execute("INSERT INTO meta_data VALUES ('A1', '1')")
    cur.execute("INSERT INTO meta_data VALUES ('A2', '1')")
    cur.execute("INSERT INTO sequence_alignment VALUES ('A1', 'A1')")
    cur.execute("INSERT INTO sequence_alignment VALUES ('A1', 'A1')")
    con.commit()
    con.close()

    gb_matrix = tmp_path / "gb.tsv"
    _write_tsv(gb_matrix, "primary_accession\thost", ["A1\thuman", "A2\tdog"])

    output = tmp_path / "out.txt"
    blast_hits = tmp_path / "blast.tsv"
    blast_hits.write_text("q\tr\tscore\tstrand\nq1\tr1\t100\t+\n", encoding="utf-8")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--mode",
            "non_segmented",
            "--blast_hits",
            str(blast_hits),
            "--gb_matrix",
            str(gb_matrix),
            "--sqlite_db",
            str(db),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 1
    assert "duplicated aligned sequence_id values: 1" in result.stdout
    assert "example duplicated sequence_alignment entries" in result.stdout
    assert "2 | A1 | alignment_name(s): A1" in result.stdout
    assert "example non-excluded meta_data accessions missing from sequence_alignment" in result.stdout
    assert "A2" in result.stdout
    assert "1 meta_data sequences missing alignment (examples: A2)" in result.stdout
    assert "traceback" not in result.stderr.lower()
