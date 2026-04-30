import csv
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from CalcAlignmentCord import CalculateAlignmentCoordinates
from GffToDictionary import GffDictionary


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "calc_alignment_cord"
SCRIPT_PATH = REPO_ROOT / "scripts" / "CalcAlignmentCord.py"
REAL_UPDATE_DB = REPO_ROOT / "test_data" / "RABV_test" / "rabv-jul0425.db"


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


@pytest.fixture()
def real_update_db_copy(tmp_path: Path) -> Path:
    if not REAL_UPDATE_DB.exists():
        pytest.skip(f"Real update DB not found at {REAL_UPDATE_DB}")
    dst = tmp_path / "rabv-jul0425.copy.db"
    shutil.copyfile(REAL_UPDATE_DB, dst)
    return dst


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
    assert adjusted == [
        {"start": 1, "end": 3, "product": "P1"},
        {"start": 4, "end": 8, "product": "P2"},
    ]


def test_get_products_for_range(tmp_path: Path):
    processor = make_processor(tmp_path)
    cds_list = [
        {"start": "1", "end": "4", "product": "P1"},
        {"start": "5", "end": "9", "product": "P2"},
    ]
    products = processor.get_products_for_range(cds_list, [4, 8])
    assert [p["product"] for p in products] == ["P1", "P2"]


def test_gff_dictionary_prefers_mature_protein_regions_for_effective_cds(tmp_path: Path):
    gff_path = tmp_path / "hcv.gff3"
    gff_path.write_text(
        "##gff-version 3\n"
        "NC_004102.1\tRefSeq\tgene\t342\t9377\t.\t+\t.\tID=gene-HCVgp1;gene=POLY\n"
        "NC_004102.1\tRefSeq\tCDS\t342\t9377\t.\t+\t0\tID=cds-NP_671491.1;product=polyprotein\n"
        "NC_004102.1\tRefSeq\tmature_protein_region_of_CDS\t3420\t5312\t.\t+\t.\tID=id-NS3;product=protease/helicase protein NS3\n"
        "NC_004102.1\tRefSeq\tmature_protein_region_of_CDS\t6258\t7601\t.\t+\t.\tID=id-NS5A;product=nonstructural protein NS5A\n"
        "NC_004102.1\tRefSeq\tCDS\t342\t369\t.\t+\t0\tID=cds-NP_803170.1;product=protein F\n",
        encoding="utf-8",
    )

    gff_dict = GffDictionary(str(gff_path)).gff_dict

    assert [entry["product"] for entry in gff_dict["CDS"]] == [
        "protease/helicase protein NS3",
        "nonstructural protein NS5A",
    ]


def test_load_blast_hits(tmp_path: Path):
    processor = make_processor(tmp_path)
    assert processor.load_blast_hits() == {"Q_A": "REF_A"}


def test_find_gaps_in_fasta_matches_expected(tmp_path: Path):
    processor = make_processor(tmp_path)
    processor.find_gaps_in_fasta()

    actual = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    expected = read_tsv_as_dicts(DATA_DIR / "expected_features.tsv")
    assert actual == expected


def test_find_gaps_in_fasta_preserves_projected_product_identity(tmp_path: Path):
    processor = make_processor(tmp_path)
    processor.find_gaps_in_fasta()

    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    q_a_rows = [row for row in rows if row["accession"] == "Q_A"]

    assert q_a_rows == [
        {
            "accession": "Q_A",
            "master_ref_accession": "MASTER1",
            "reference_accession": "REF_A",
            "aln_start": "1",
            "aln_end": "9",
            "cds_start": "1",
            "cds_end": "3",
            "product": "P1",
        },
        {
            "accession": "Q_A",
            "master_ref_accession": "MASTER1",
            "reference_accession": "REF_A",
            "aln_start": "1",
            "aln_end": "9",
            "cds_start": "4",
            "cds_end": "8",
            "product": "P2",
        },
    ]


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


def test_find_gaps_in_fasta_raises_when_update_scope_missing_primary_accession(tmp_path: Path):
    bad_scope = tmp_path / "bad_scope.tsv"
    bad_scope.write_text("wrong_col\nQ_A\n", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(DATA_DIR / "padded_alignment"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(DATA_DIR / "query_uniq_tophits.tsv"),
        update_scope_tsv=str(bad_scope),
    )

    with pytest.raises(ValueError, match="Update scope TSV is missing required column"):
        processor.find_gaps_in_fasta()


def test_find_gaps_in_fasta_raises_when_segment_map_missing_columns(tmp_path: Path):
    bad_segment_map = tmp_path / "bad_segment_map.tsv"
    bad_segment_map.write_text("primary_accession\nQ_A\n", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(DATA_DIR / "padded_alignment"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(DATA_DIR / "query_uniq_tophits.tsv"),
        segment_map_tsv=str(bad_segment_map),
    )

    with pytest.raises(ValueError, match="Segment map TSV is missing required columns"):
        processor.find_gaps_in_fasta()


def test_find_gaps_in_fasta_raises_when_update_db_features_missing_accession(tmp_path: Path):
    bad_db = tmp_path / "bad_update.db"
    conn = sqlite3.connect(str(bad_db))
    try:
        conn.execute("CREATE TABLE features (foo TEXT)")
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
        update_db=str(bad_db),
    )

    with pytest.raises(ValueError, match="features table is missing required column: accession"):
        processor.find_gaps_in_fasta()


def test_find_gaps_in_fasta_update_mode_real_db_skips_existing_scoped_accession(tmp_path: Path, real_update_db_copy: Path):
    conn = sqlite3.connect(str(real_update_db_copy))
    try:
        cols = [row[1] for row in conn.execute("PRAGMA table_info(features)").fetchall()]
        if "accession" not in cols:
            pytest.skip("Real DB features table lacks accession column")

        values = {c: None for c in cols}
        values["accession"] = "Q_A"
        insert_cols = list(values.keys())
        placeholders = ",".join(["?"] * len(insert_cols))
        conn.execute(
            f"INSERT INTO features ({','.join(insert_cols)}) VALUES ({placeholders})",
            [values[c] for c in insert_cols],
        )
        conn.commit()
    finally:
        conn.close()

    scope = tmp_path / "scope.tsv"
    scope.write_text("primary_accession\nQ_A\n", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(DATA_DIR / "padded_alignment"),
        master_gff=[str(DATA_DIR / "MASTER1.gff3")],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(DATA_DIR / "master_list.tsv"),
        blast_uniq_hits=str(DATA_DIR / "query_uniq_tophits.tsv"),
        update_db=str(real_update_db_copy),
        update_scope_tsv=str(scope),
    )

    processor.find_gaps_in_fasta()
    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    assert rows == []
