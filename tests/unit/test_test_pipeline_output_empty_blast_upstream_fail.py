import sqlite3
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "TestPipelineOutput.py"


def _write_tsv(path: Path, header: str, rows):
    path.write_text("\n".join([header] + list(rows)) + "\n", encoding="utf-8")


def test_non_segmented_empty_blast_with_non_excluded_rows_fails(tmp_path: Path):
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

    blast_hits = tmp_path / "blast.tsv"
    blast_hits.write_text("", encoding="utf-8")

    gb_matrix = tmp_path / "gb.tsv"
    _write_tsv(gb_matrix, "primary_accession\thost\texclusion_status", ["A1\thuman\t0"])

    output = tmp_path / "out.txt"
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
    assert "upstream BLAST/alignment failure" in result.stdout
