import csv
import subprocess
import sys
from pathlib import Path

import pytest

from CalcAlignmentCord import CalculateAlignmentCoordinates


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "calc_alignment_cord"
SCRIPT_PATH = REPO_ROOT / "scripts" / "CalcAlignmentCord.py"


def read_tsv_as_dicts(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def make_processor(tmp_path: Path):
    return CalculateAlignmentCoordinates(
        paded_alignment=str(DATA_DIR / "padded_alignment"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(DATA_DIR / "query_uniq_tophits.tsv"),
    )


def test_get_master_list_and_gff_resolution(tmp_path: Path):
    processor = make_processor(tmp_path)
    assert processor.get_master_list() == ["MASTER1"]
    assert processor.get_gff_for_master("MASTER1").endswith("MASTER1.gff3")


def test_gap_helpers_and_cds_recalculation(tmp_path: Path):
    processor = make_processor(tmp_path)

    gaps = processor.get_gap_ranges("ATG-CGTAA")
    assert gaps == [[4, 4]]
    assert processor.count_gaps_before_position(gaps, 3) == 0
    assert processor.count_gaps_before_position(gaps, 4) == 1
    assert processor.count_gaps_before_position(gaps, 9) == 1

    cds_list = [
        {"start": "1", "end": "4", "product": "P1"},
        {"start": "5", "end": "9", "product": "P2"},
    ]
    adjusted = processor.recalculate_cds_coordinates("Q_A", gaps, cds_list, start_offset=1)
    assert adjusted == [[1, 3], [4, 8]]


def test_get_products_for_range(tmp_path: Path):
    processor = make_processor(tmp_path)
    cds_list = [
        {"start": "1", "end": "4", "product": "P1"},
        {"start": "5", "end": "9", "product": "P2"},
    ]
    products = processor.get_products_for_range(cds_list, [4, 8])
    assert [p["product"] for p in products] == ["P1", "P2"]


def test_load_blast_hits(tmp_path: Path):
    processor = make_processor(tmp_path)
    assert processor.load_blast_hits() == {"Q_A": "REF_A"}


def test_find_gaps_in_fasta_matches_expected(tmp_path: Path):
    processor = make_processor(tmp_path)
    processor.find_gaps_in_fasta()

    actual = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    expected = read_tsv_as_dicts(DATA_DIR / "expected_features.tsv")
    assert actual == expected


def test_cli_generates_expected_features(tmp_path: Path):
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-i",
            str(DATA_DIR / "padded_alignment"),
            "-g",
            str(DATA_DIR / "MASTER1.gff3"),
            "-b",
            str(tmp_path),
            "-d",
            "Tables",
            "-o",
            "features.tsv",
            "-m",
            str(DATA_DIR / "master_list.tsv"),
            "-bh",
            str(DATA_DIR / "query_uniq_tophits.tsv"),
        ],
        check=True,
    )

    actual = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    expected = read_tsv_as_dicts(DATA_DIR / "expected_features.tsv")
    assert actual == expected


def test_find_gaps_in_fasta_raises_when_alignment_dir_missing(tmp_path: Path):
    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(tmp_path / "missing_dir"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(DATA_DIR / "query_uniq_tophits.tsv"),
    )
    with pytest.raises(FileNotFoundError, match="Padded alignment directory not found"):
        processor.find_gaps_in_fasta()


def test_load_blast_hits_raises_on_malformed_row(tmp_path: Path):
    bad_hits = tmp_path / "bad_hits.tsv"
    bad_hits.write_text("Q_A\tREF_A\t99.9\n", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(DATA_DIR / "padded_alignment"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(bad_hits),
    )

    with pytest.raises(ValueError, match="Malformed BLAST hits row"):
        processor.load_blast_hits()


def test_find_gaps_in_fasta_update_scope_emits_only_scoped_accessions(tmp_path: Path):
    scope = tmp_path / "scope.tsv"
    scope.write_text("primary_accession\nNON_EXISTENT\n", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(DATA_DIR / "padded_alignment"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(DATA_DIR / "query_uniq_tophits.tsv"),
        update_scope_tsv=str(scope),
    )
    processor.find_gaps_in_fasta()

    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    assert rows == []


def test_find_gaps_in_fasta_update_db_skips_existing_feature_accessions(tmp_path: Path):
    update_db = tmp_path / "update.db"
    import sqlite3

    conn = sqlite3.connect(str(update_db))
    try:
        conn.execute("CREATE TABLE features (accession TEXT)")
        conn.execute("INSERT INTO features(accession) VALUES ('Q_A')")
        conn.commit()
    finally:
        conn.close()

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(DATA_DIR / "padded_alignment"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(DATA_DIR / "query_uniq_tophits.tsv"),
        update_db=str(update_db),
    )
    processor.find_gaps_in_fasta()

    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    assert rows
    assert all(row["accession"] != "Q_A" for row in rows)


def test_find_gaps_in_fasta_raises_on_segment_mismatch(tmp_path: Path):
    segment_map = tmp_path / "segment_map.tsv"
    segment_map.write_text(
        "primary_accession\tsegment\n"
        "Q_A\t2\n"
        "MASTER1\t1\n",
        encoding="utf-8",
    )

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(DATA_DIR / "padded_alignment"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(DATA_DIR / "query_uniq_tophits.tsv"),
        segment_map_tsv=str(segment_map),
    )

    with pytest.raises(ValueError, match="Segment mismatch"):
        processor.find_gaps_in_fasta()
