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
        {"start": 1, "end": 4, "og_start": 1, "og_end": 3, "product": "P1"},
        {"start": 5, "end": 9, "og_start": 4, "og_end": 8, "product": "P2"},
    ]


def test_recalculate_cds_coordinates_keeps_reference_feature_span_for_partial_sequences(tmp_path: Path):
    processor = make_processor(tmp_path)

    gaps = processor.get_gap_ranges("--GCCGT--")
    cds_list = [
        {"start": "1", "end": "4", "product": "P1"},
        {"start": "5", "end": "9", "product": "P2"},
    ]

    adjusted = processor.recalculate_cds_coordinates_with_span(
        "Q_B",
        gaps,
        cds_list,
        start_offset=3,
        genome_cord_start=3,
        genome_cord_end=7,
    )

    assert adjusted == [
        {"start": 3, "end": 4, "og_start": 1, "og_end": 2, "product": "P1"},
        {"start": 5, "end": 7, "og_start": 3, "og_end": 5, "product": "P2"},
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


def test_find_gaps_in_fasta_update_mode_recalculates_existing_accessions(tmp_path: Path):
    update_db = tmp_path / "update.db"
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
            "cds_end": "4",
            "cds_start_OG_seq": "1",
            "cds_end_OG_seq": "3",
            "product": "P1",
        },
        {
            "accession": "Q_A",
            "master_ref_accession": "MASTER1",
            "reference_accession": "REF_A",
            "aln_start": "1",
            "aln_end": "9",
            "cds_start": "5",
            "cds_end": "9",
            "cds_start_OG_seq": "4",
            "cds_end_OG_seq": "8",
            "product": "P2",
        },
    ]


def test_find_gaps_in_fasta_partial_sequences_keep_overlapping_reference_feature_coordinates(tmp_path: Path):
    processor = make_processor(tmp_path)
    processor.find_gaps_in_fasta()

    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    q_b_rows = [row for row in rows if row["accession"] == "Q_B"]

    assert q_b_rows == [
        {
            "accession": "Q_B",
            "master_ref_accession": "MASTER1",
            "reference_accession": "MASTER1",
            "aln_start": "3",
            "aln_end": "7",
            "cds_start": "3",
            "cds_end": "4",
            "cds_start_OG_seq": "1",
            "cds_end_OG_seq": "2",
            "product": "P1",
        },
        {
            "accession": "Q_B",
            "master_ref_accession": "MASTER1",
            "reference_accession": "MASTER1",
            "aln_start": "3",
            "aln_end": "7",
            "cds_start": "5",
            "cds_end": "7",
            "cds_start_OG_seq": "3",
            "cds_end_OG_seq": "5",
            "product": "P2",
        },
    ]


def test_find_gaps_in_fasta_skips_features_outside_partial_sequence_span(tmp_path: Path):
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir()
    (alignment_dir / "MASTERX.aligned_merged_MSA.fasta").write_text(
        ">MASTERX\n"
        "ATGCCGTAA\n"
        ">Q_PARTIAL\n"
        "ATG------\n",
        encoding="utf-8",
    )

    gff_path = tmp_path / "MASTERX.gff3"
    gff_path.write_text(
        "##gff-version 3\n"
        "MASTERX\tRefSeq\tCDS\t1\t4\t.\t+\t0\tID=cds1;product=P1\n"
        "MASTERX\tRefSeq\tCDS\t5\t9\t.\t+\t0\tID=cds2;product=P2\n",
        encoding="utf-8",
    )

    master_list = tmp_path / "master_list.tsv"
    master_list.write_text("MASTERX\n", encoding="utf-8")

    blast_hits = tmp_path / "query_uniq_tophits.tsv"
    blast_hits.write_text("Q_PARTIAL\tREF_PARTIAL\t99.0\tplus\n", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir),
        master_gff=[str(gff_path)],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(master_list),
        blast_uniq_hits=str(blast_hits),
    )

    processor.find_gaps_in_fasta()

    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    q_partial_rows = [row for row in rows if row["accession"] == "Q_PARTIAL"]

    assert q_partial_rows == [
        {
            "accession": "Q_PARTIAL",
            "master_ref_accession": "MASTERX",
            "reference_accession": "REF_PARTIAL",
            "aln_start": "1",
            "aln_end": "3",
            "cds_start": "1",
            "cds_end": "3",
            "cds_start_OG_seq": "1",
            "cds_end_OG_seq": "3",
            "product": "P1",
        }
    ]


def test_find_gaps_in_fasta_resolves_segmented_refset_files_to_matching_master(tmp_path: Path):
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir()
    (alignment_dir / "refset_2_aln_merged_MSA.fasta").write_text(
        ">MASTER2\n"
        "ATGCCGTAA\n"
        ">Q_SEG\n"
        "ATG-CGTAA\n",
        encoding="utf-8",
    )

    gff1 = tmp_path / "MASTER1.gff3"
    gff1.write_text(
        "##gff-version 3\n"
        "MASTER1\tRefSeq\tCDS\t1\t4\t.\t+\t0\tID=cds1;product=P1\n",
        encoding="utf-8",
    )
    gff2 = tmp_path / "MASTER2.gff3"
    gff2.write_text(
        "##gff-version 3\n"
        "MASTER2\tRefSeq\tCDS\t1\t4\t.\t+\t0\tID=cds2;product=P2\n",
        encoding="utf-8",
    )

    ref_list = tmp_path / "ref_list.tsv"
    ref_list.write_text(
        "MASTER1\tmaster\t1\n"
        "MASTER2\tmaster\t2\n",
        encoding="utf-8",
    )

    blast_hits = tmp_path / "query_uniq_tophits.tsv"
    blast_hits.write_text("Q_SEG\tREF_SEG\t99.0\tplus\n", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir),
        master_gff=[str(gff1), str(gff2)],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(ref_list),
        blast_uniq_hits=str(blast_hits),
    )

    processor.find_gaps_in_fasta()

    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    q_seg_rows = [row for row in rows if row["accession"] == "Q_SEG"]

    assert q_seg_rows == [
        {
            "accession": "Q_SEG",
            "master_ref_accession": "MASTER2",
            "reference_accession": "REF_SEG",
            "aln_start": "1",
            "aln_end": "9",
            "cds_start": "1",
            "cds_end": "4",
            "cds_start_OG_seq": "1",
            "cds_end_OG_seq": "3",
            "product": "P2",
            "genome_coverage": "75.00",
        }
    ]


def test_find_gaps_in_fasta_segmented_refset_emits_segment_column_from_segment_map(tmp_path: Path):
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir()
    (alignment_dir / "refset_2_aln_merged_MSA.fasta").write_text(
        ">MASTER2\n"
        "ATGCCGTAA\n"
        ">Q_SEG\n"
        "ATG-CGTAA\n",
        encoding="utf-8",
    )

    gff1 = tmp_path / "MASTER1.gff3"
    gff1.write_text(
        "##gff-version 3\n"
        "MASTER1\tRefSeq\tCDS\t1\t4\t.\t+\t0\tID=cds1;product=P1\n",
        encoding="utf-8",
    )
    gff2 = tmp_path / "MASTER2.gff3"
    gff2.write_text(
        "##gff-version 3\n"
        "MASTER2\tRefSeq\tCDS\t1\t4\t.\t+\t0\tID=cds2;product=P2\n",
        encoding="utf-8",
    )

    ref_list = tmp_path / "ref_list.tsv"
    ref_list.write_text(
        "MASTER1\tmaster\t1\n"
        "MASTER2\tmaster\t2\n",
        encoding="utf-8",
    )

    blast_hits = tmp_path / "query_uniq_tophits.tsv"
    blast_hits.write_text("Q_SEG\tREF_SEG\t99.0\tplus\n", encoding="utf-8")

    segment_map = tmp_path / "segment_map.tsv"
    segment_map.write_text(
        "primary_accession\tsegment\n"
        "MASTER2\t2\n"
        "REF_SEG\t2\n"
        "Q_SEG\t2\n",
        encoding="utf-8",
    )

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir),
        master_gff=[str(gff1), str(gff2)],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(ref_list),
        blast_uniq_hits=str(blast_hits),
        segment_map_tsv=str(segment_map),
    )

    processor.find_gaps_in_fasta()

    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    q_seg_rows = [row for row in rows if row["accession"] == "Q_SEG"]

    assert q_seg_rows == [
        {
            "accession": "Q_SEG",
            "master_ref_accession": "MASTER2",
            "reference_accession": "REF_SEG",
            "aln_start": "1",
            "aln_end": "9",
            "cds_start": "1",
            "cds_end": "4",
            "cds_start_OG_seq": "1",
            "cds_end_OG_seq": "3",
            "product": "P2",
            "genome_coverage": "75.00",
            "segment": "2",
        }
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
        update_db=str(update_db),
        update_scope_tsv=str(scope),
    )
    processor.find_gaps_in_fasta()

    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    assert rows == []


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


def test_stress_recalculate_cds_coordinates_extreme_gaps(tmp_path: Path):
    processor = make_processor(tmp_path)
    # 1. Alternating gaps and bases
    gaps = processor.get_gap_ranges("-A-C-G-T-")
    assert gaps == [[1, 1], [3, 3], [5, 5], [7, 7], [9, 9]]
    
    # 2. Query has 100% gaps (only deletion)
    gaps_all = processor.get_gap_ranges("---------")
    assert gaps_all == [[1, 9]]
    
    # 3. Gap offsets flanking the CDS
    cds_list = [{"start": "4", "end": "6", "product": "P1"}]
    # gaps cover first 3 positions and last 3 positions
    gaps_flank = processor.get_gap_ranges("---ACG---")
    # gap ranges: [1, 3] and [7, 9]
    assert gaps_flank == [[1, 3], [7, 9]]
    
    # Recalculate query covering overlap. Without span clamping.
    adjusted = processor.recalculate_cds_coordinates_with_span(
        "Q_STRESS", gaps_flank, cds_list, start_offset=4, genome_cord_start=None, genome_cord_end=None
    )
    # Overlap start=4, gaps before 4 is 3. So og_start = 4 - 3 = 1
    # Overlap end=6, gaps before 6 is 3. So og_end = 6 - 3 = 3
    assert adjusted == [
        {"start": 4, "end": 6, "og_start": 1, "og_end": 3, "product": "P1"}
    ]


def test_stress_recalculate_cds_coordinates_partial_overlap_boundaries(tmp_path: Path):
    processor = make_processor(tmp_path)
    cds_list = [{"start": "100", "end": "200", "product": "P_MID"}]
    gaps = [] # no gaps
    
    # 1. Span completely before CDS
    adjusted = processor.recalculate_cds_coordinates_with_span(
        "Q", gaps, cds_list, start_offset=1, genome_cord_start=10, genome_cord_end=50
    )
    assert adjusted == []
    
    # 2. Span completely after CDS
    adjusted = processor.recalculate_cds_coordinates_with_span(
        "Q", gaps, cds_list, start_offset=1, genome_cord_start=250, genome_cord_end=300
    )
    assert adjusted == []
    
    # 3. Span overlaps CDS by exactly 1 base at start
    adjusted = processor.recalculate_cds_coordinates_with_span(
        "Q", gaps, cds_list, start_offset=1, genome_cord_start=50, genome_cord_end=100
    )
    assert adjusted == [{"start": 100, "end": 100, "og_start": 100, "og_end": 100, "product": "P_MID"}]
    
    # 4. Span overlaps CDS by exactly 1 base at end
    adjusted = processor.recalculate_cds_coordinates_with_span(
        "Q", gaps, cds_list, start_offset=1, genome_cord_start=200, genome_cord_end=250
    )
    assert adjusted == [{"start": 200, "end": 200, "og_start": 200, "og_end": 200, "product": "P_MID"}]


def test_stress_recalculate_cds_coordinates_out_of_bounds_span(tmp_path: Path):
    processor = make_processor(tmp_path)
    cds_list = [{"start": "10", "end": "20", "product": "P"}]
    
    # 1. Invalid coordinates: start > end (should be skipped)
    adjusted = processor.recalculate_cds_coordinates_with_span(
        "Q", [], cds_list, start_offset=1, genome_cord_start=25, genome_cord_end=5
    )
    assert adjusted == []
    
    # 2. Huge/out-of-bounds coordinates (clamped correctly)
    adjusted = processor.recalculate_cds_coordinates_with_span(
        "Q", [], cds_list, start_offset=1, genome_cord_start=1, genome_cord_end=999999
    )
    assert adjusted == [{"start": 10, "end": 20, "og_start": 10, "og_end": 20, "product": "P"}]


def test_stress_find_gaps_in_fasta_empty_and_corrupt_files(tmp_path: Path):
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir()
    
    # FASTA containing lowercase, invalid characters (X, N, ?, etc.)
    (alignment_dir / "MASTERX.aligned_merged_MSA.fasta").write_text(
        ">MASTERX\n"
        "atgccgtaA\n"
        ">Q_CORRUPT\n"
        "aTgXXnN??\n",
        encoding="utf-8",
    )
    
    gff_path = tmp_path / "MASTERX.gff3"
    gff_path.write_text(
        "##gff-version 3\n"
        "MASTERX\tRefSeq\tCDS\t1\t9\t.\t+\t0\tID=cds1;product=P_CORRUPT\n",
        encoding="utf-8",
    )
    
    master_list = tmp_path / "master_list.tsv"
    master_list.write_text("MASTERX\n", encoding="utf-8")
    
    blast_hits = tmp_path / "query_uniq_tophits.tsv"
    blast_hits.write_text("Q_CORRUPT\tREF_CORRUPT\t99.0\tplus\n", encoding="utf-8")
    
    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir),
        master_gff=[str(gff_path)],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(master_list),
        blast_uniq_hits=str(blast_hits),
    )
    processor.find_gaps_in_fasta()
    
    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    assert len(rows) == 2  # MASTERX and Q_CORRUPT
    assert rows[1]["accession"] == "Q_CORRUPT"
    assert rows[1]["cds_start_OG_seq"] == "1"
    assert rows[1]["cds_end_OG_seq"] == "9"


def test_stress_find_gaps_in_fasta_segment_checking_edge_cases(tmp_path: Path):
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir()
    (alignment_dir / "refset_1_aln_merged_MSA.fasta").write_text(
        ">MASTER1\n"
        "ATGC\n"
        ">Q_A\n"
        "ATGC\n",
        encoding="utf-8",
    )
    gff1 = tmp_path / "MASTER1.gff3"
    gff1.write_text(
        "##gff-version 3\n"
        "MASTER1\tRefSeq\tCDS\t1\t4\t.\t+\t0\tID=cds1;product=P1\n",
        encoding="utf-8",
    )
    ref_list = tmp_path / "ref_list.tsv"
    ref_list.write_text("MASTER1\tmaster\t1\n", encoding="utf-8")
    blast_hits = tmp_path / "query_uniq_tophits.tsv"
    blast_hits.write_text("Q_A\tREF_A\t99.0\tplus\n", encoding="utf-8")
    
    # 1. Segment map with blank/whitespace entries
    segment_map = tmp_path / "segment_map.tsv"
    segment_map.write_text(
        "primary_accession\tsegment\n"
        "MASTER1\t 1 \n"   # spaces should be stripped
        "REF_A\t1\n"
        "Q_A\t\n",         # blank segment is handled safely
        encoding="utf-8",
    )
    
    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir),
        master_gff=[str(gff1)],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(ref_list),
        blast_uniq_hits=str(blast_hits),
        segment_map_tsv=str(segment_map),
    )
    processor.find_gaps_in_fasta()
    rows = read_tsv_as_dicts(tmp_path / "Tables" / "features.tsv")
    assert rows[1]["segment"] == ""  # Q_A segment maps to empty

    # 2. Conflicting segment mappings raise ValueError
    bad_segment_map = tmp_path / "bad_segment_map.tsv"
    bad_segment_map.write_text(
        "primary_accession\tsegment\n"
        "MASTER1\t1\n"
        "REF_A\t2\n"  # conflicts with segment of MASTER1 (1) because both are aligned together
        "Q_A\t1\n",
        encoding="utf-8",
    )
    processor_conflict = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir),
        master_gff=[str(gff1)],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(ref_list),
        blast_uniq_hits=str(blast_hits),
        segment_map_tsv=str(bad_segment_map),
    )
    with pytest.raises(ValueError, match="Segment mismatch for Q_A"):
        processor_conflict.find_gaps_in_fasta()

