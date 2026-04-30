import sqlite3
import subprocess
from pathlib import Path

import pandas as pd

from DownloadGFF import NCBI_GFF_Downloader


def test_get_accession_list_uses_update_db_master_rows(tmp_path: Path):
    update_db = tmp_path / "update.db"
    with sqlite3.connect(str(update_db)) as conn:
        conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
        conn.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
            [("REF1", "master", "1"), ("REF2", "reference", "2"), ("Q1", "query", "1")],
        )
        conn.commit()

    downloader = NCBI_GFF_Downloader("ignored.tsv", str(tmp_path), "Gff", update_db=str(update_db))
    assert downloader.get_accession_list() == ["REF1"]


def test_get_accession_list_from_master_marked_file(tmp_path: Path):
    ref_file = tmp_path / "refs.tsv"
    pd.DataFrame(
        [
            ["ACC_MASTER_1", "master"],
            ["ACC_REF_1", "reference"],
            ["ACC_MASTER_2", "MASTER"],
        ]
    ).to_csv(ref_file, sep="\t", header=False, index=False)

    d = NCBI_GFF_Downloader(str(ref_file), str(tmp_path), "Gff")
    assert d.get_accession_list() == ["ACC_MASTER_1", "ACC_MASTER_2"]


def test_get_accession_list_falls_back_to_first_column(tmp_path: Path):
    ref_file = tmp_path / "refs.tsv"
    pd.DataFrame([["A1"], ["A2"]]).to_csv(ref_file, sep="\t", header=False, index=False)

    d = NCBI_GFF_Downloader(str(ref_file), str(tmp_path), "Gff")
    assert d.get_accession_list() == ["A1", "A2"]


def test_get_accession_list_from_headered_master_marked_file(tmp_path: Path):
    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text(
        "primary_accession\tstatus\tsegment\nACC_MASTER_1\tmaster\t1\nACC_REF_1\treference\t1\n",
        encoding="utf-8",
    )

    d = NCBI_GFF_Downloader(str(ref_file), str(tmp_path), "Gff")
    assert d.get_accession_list() == ["ACC_MASTER_1"]


def test_get_accession_list_from_comma_separated_string(tmp_path: Path):
    d = NCBI_GFF_Downloader("A1, A2 ,,A3", str(tmp_path), "Gff")
    assert d.get_accession_list() == ["A1", "A2", "A3"]


def test_download_gff_skips_existing_and_handles_errors(tmp_path: Path, monkeypatch):
    out_dir = tmp_path / "Gff"
    out_dir.mkdir(parents=True, exist_ok=True)

    existing = out_dir / "ACC_SKIP.gff3"
    existing.write_text("already here\n", encoding="utf-8")

    d = NCBI_GFF_Downloader("ACC_SKIP,ACC_OK,ACC_FAIL", str(tmp_path), "Gff")

    calls = []

    def fake_run(command, check, stdout, stderr):
        acc = command[4]
        calls.append(acc)
        if acc == "ACC_FAIL":
            raise subprocess.CalledProcessError(returncode=1, cmd=command, stderr=b"boom")
        stdout.write("##gff-version 3\n")

    monkeypatch.setattr("DownloadGFF.subprocess.run", fake_run)

    d.download_gff()

    # Existing file is skipped; subprocess called only for ACC_OK and ACC_FAIL
    assert calls == ["ACC_OK", "ACC_FAIL"]
    assert (out_dir / "ACC_OK.gff3").read_text(encoding="utf-8").startswith("##gff-version 3")
    assert (out_dir / "ACC_FAIL.gff3").exists()
