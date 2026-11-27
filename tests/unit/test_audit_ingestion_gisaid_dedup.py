"""Carpet-bomb audit tests for ingestion, XML parsing, GISAID, dedup, and metadata.

Adds passing coverage for gisaid_tidy, HostTaxaTable (all 4 outputs), and
ValidateRefList, plus strict xfail reproducers (AUDIT-ING-*).
"""

import runpy
import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

import DedupStrainAlignment as dsa_mod
import FilterAndExtractSequences as fae_mod
import GenBankFetcher as gbf_mod
import GenBankParser as gbp_mod
import GffToDictionary as gff_mod
import HostTaxaTable as htt_mod
import ValidateMatrix as vm_mod
import ValidateStrain as vs_mod
import accession_utils as au_mod
import cross_db_dedup as cdd_mod
import date_utils as du_mod
import gisaid_tidy as gt_mod
import merge_into_gB_matrix as mgm_mod

REPO_ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------
# Passing coverage lifts: gisaid_tidy, HostTaxaTable, ValidateRefList
# ---------------------------------------------------------------------------

def test_gisaid_tidy_long_and_wide_format_with_db_exclusion(tmp_path):
    data_dir = tmp_path / "raw"
    out_dir = tmp_path / "out"
    data_dir.mkdir()
    db_path = tmp_path / "prev.db"
    conn = sqlite3.connect(str(db_path))
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, exclusion_status TEXT)")
    conn.execute("INSERT INTO meta_data VALUES ('EPI100001', '0')")
    conn.commit()
    conn.close()

    # Wide-format TSV + FASTA
    (data_dir / "meta.tsv").write_text(
        "Isolate_Id\tIsolate_Name\tHA Segment_Id\tNA Segment_Id\tLocation\n"
        "EPI_ISL_1\tA/duck/UK/1/2024\tEPI100001\tEPI100002\tEurope / United Kingdom\n",
        encoding="utf-8",
    )
    (data_dir / "seqs.fas").write_text(
        ">EPI100002|NA|EPI_ISL_1\nATGCGATGCGATGCG\n"
        ">EPI100001|HA|EPI_ISL_1\nATGCGATGCGATGCG\n",
        encoding="utf-8",
    )

    tidy = gt_mod.GISAIDTidy(
        data_dir=str(data_dir),
        output_dir=str(out_dir),
        filetype="tsv",
        log_file=str(out_dir / "run.log"),
        db_file=str(db_path),
    )
    tidy.process()

    meta_out = out_dir / "metadata.tsv"
    fas_out = out_dir / "all_nuc.fas"
    assert meta_out.exists() and fas_out.exists()
    df = pd.read_csv(meta_out, sep="\t", dtype=str, keep_default_na=False)
    # EPI100001 already in DB, only EPI100002 should remain
    assert list(df["Segment_Id"]) == ["EPI100002"]
    assert df.loc[0, "Segment"] == "NA"
    assert ">EPI100002" in fas_out.read_text(encoding="utf-8")


def test_host_taxa_table_writes_all_four_tables(tmp_path):
    gb = tmp_path / "gb.tsv"
    gb.write_text("primary_accession\thost_taxa_id\nACC1\t9606\nACC2\tNA\n", encoding="utf-8")
    names = tmp_path / "names.dmp"
    names.write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|\n"
        "9606\t|\thuman\t|\t\t|\tgenbank common name\t|\n"
        "9605\t|\tHomo\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "9605\t|\t1\t|\tgenus\t|\n"
        "9606\t|\t9605\t|\tspecies\t|\n",
        encoding="utf-8",
    )

    ht = htt_mod.HostTaxaTable(
        gb_matrix=str(gb),
        output_dir="HostTaxa",
        names_file=str(names),
        nodes_file=str(nodes),
        base_dir=str(tmp_path),
        host_output_file="Host_taxa.tsv",
        child_output_file="Host_taxa_children.tsv",
        lineage_output_file="Host_taxa_lineage.tsv",
        lineage_lookup_output_file="Host_taxa_lineage_lookup.tsv",
    )
    ht.write_tables()
    for fname in (
        "Host_taxa.tsv",
        "Host_taxa_children.tsv",
        "Host_taxa_lineage.tsv",
        "Host_taxa_lineage_lookup.tsv",
    ):
        assert (tmp_path / "HostTaxa" / fname).exists()


def test_validate_ref_list_passes_valid_segmented_and_unsegmented(tmp_path, monkeypatch):
    ref_n = tmp_path / "ref_n.txt"
    ref_n.write_text("NC_001542\tmaster\nKF154996\treference\n", encoding="utf-8")
    out_n = tmp_path / "out_n.txt"
    monkeypatch.setattr(sys, "argv", ["ValidateRefList.py", "-r", str(ref_n), "-s", "N", "-o", str(out_n)])
    runpy.run_path(str(REPO_ROOT / "scripts/ValidateRefList.py"), run_name="__main__")
    assert "passed" in out_n.read_text(encoding="utf-8").lower()

    ref_y = tmp_path / "ref_y.txt"
    ref_y.write_text("M1\tmaster\t1\nR1\treference\t1\nEX1\texclusion_list\tB\n", encoding="utf-8")
    out_y = tmp_path / "out_y.txt"
    monkeypatch.setattr(sys, "argv", ["ValidateRefList.py", "-r", str(ref_y), "-s", "Y", "-o", str(out_y)])
    runpy.run_path(str(REPO_ROOT / "scripts/ValidateRefList.py"), run_name="__main__")
    assert out_y.exists()


# ---------------------------------------------------------------------------
# Promoted regression tests for fixed Ingestion & GFF bugs (WP-04) & remaining xfails
# ---------------------------------------------------------------------------

def test_genbank_fetcher_preserves_meta_exclusions_when_excluded_table_exists(tmp_path):
    db = tmp_path / "u.db"
    conn = sqlite3.connect(str(db))
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, exclusion_status TEXT)")
    conn.execute("INSERT INTO meta_data VALUES ('BAD001', '1')")
    conn.execute("CREATE TABLE excluded_accessions (primary_accession TEXT)")
    conn.commit()
    conn.close()

    fetcher = gbf_mod.GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="t@example.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir=str(tmp_path),
        update_file=None,
    )
    _, excluded = fetcher._load_db_accession_context(str(db))
    assert "BAD001" in excluded


def test_genbank_parser_handles_missing_xml_subtags_and_empty_sequence(tmp_path):
    xml = tmp_path / "minimal.xml"
    xml.write_text(
        "<GBSet><GBSeq>"
        "<GBSeq_primary-accession>ACC1</GBSeq_primary-accession>"
        "<GBSeq_accession-version>ACC1.1</GBSeq_accession-version>"
        "<GBSeq_length>0</GBSeq_length>"
        "<GBSeq_division>VRL</GBSeq_division>"
        "<GBSeq_feature-table><GBFeature><GBFeature_key>source</GBFeature_key><GBFeature_quals>"
        "<GBQualifier><GBQualifier_name>country</GBQualifier_name></GBQualifier>"
        "</GBFeature_quals></GBFeature></GBSeq_feature-table>"
        "<GBSeq_sequence/>"
        "</GBSeq></GBSet>",
        encoding="utf-8",
    )
    parser = gbp_mod.GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=None,
        exclusion_list=None,
        is_segmented_virus="N",
    )
    rows = parser.xml_to_tsv(str(xml), {}, [])
    assert len(rows) == 1
    assert rows[0]["sequence"] == ""


def test_genbank_parser_does_not_overwrite_exclusion_list_reason(tmp_path):
    xml_dir = tmp_path / "xml"
    xml_dir.mkdir()
    (xml_dir / "batch.xml").write_text(
        "<GBSet>"
        "<GBSeq>"
        "<GBSeq_locus>M1</GBSeq_locus><GBSeq_length>1000</GBSeq_length><GBSeq_strandedness>ss</GBSeq_strandedness>"
        "<GBSeq_moltype>RNA</GBSeq_moltype><GBSeq_topology>linear</GBSeq_topology><GBSeq_division>VRL</GBSeq_division>"
        "<GBSeq_update-date>01-JAN-2024</GBSeq_update-date><GBSeq_create-date>01-JAN-2024</GBSeq_create-date>"
        "<GBSeq_definition>Master</GBSeq_definition><GBSeq_primary-accession>M1</GBSeq_primary-accession>"
        "<GBSeq_accession-version>M1.1</GBSeq_accession-version><GBSeq_other-seqids/>"
        "<GBSeq_source>Virus</GBSeq_source><GBSeq_organism>Virus</GBSeq_organism><GBSeq_taxonomy>Viruses</GBSeq_taxonomy>"
        "<GBSeq_feature-table/><GBSeq_sequence>" + ("A" * 1000) + "</GBSeq_sequence>"
        "</GBSeq>"
        "<GBSeq>"
        "<GBSeq_locus>EX1</GBSeq_locus><GBSeq_length>10</GBSeq_length><GBSeq_strandedness>ss</GBSeq_strandedness>"
        "<GBSeq_moltype>RNA</GBSeq_moltype><GBSeq_topology>linear</GBSeq_topology><GBSeq_division>VRL</GBSeq_division>"
        "<GBSeq_update-date>01-JAN-2024</GBSeq_update-date><GBSeq_create-date>01-JAN-2024</GBSeq_create-date>"
        "<GBSeq_definition>Short outgroup</GBSeq_definition><GBSeq_primary-accession>EX1</GBSeq_primary-accession>"
        "<GBSeq_accession-version>EX1.1</GBSeq_accession-version><GBSeq_other-seqids/>"
        "<GBSeq_source>Virus</GBSeq_source><GBSeq_organism>Virus</GBSeq_organism><GBSeq_taxonomy>Viruses</GBSeq_taxonomy>"
        "<GBSeq_feature-table/><GBSeq_sequence>ATGCATGCAT</GBSeq_sequence>"
        "</GBSeq>"
        "</GBSet>",
        encoding="utf-8",
    )
    ref_list = tmp_path / "refs.txt"
    ref_list.write_text("M1\tmaster\nEX1\texclusion_list\n", encoding="utf-8")
    parser = gbp_mod.GenBankParser(
        input_dir=str(xml_dir),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path),
        ref_list=str(ref_list),
        exclusion_list=None,
        is_segmented_virus="N",
        min_length_ratio=0.5,
    )
    parser.process()
    df = pd.read_csv(tmp_path / "gB_matrix_raw.tsv", sep="\t", dtype=str, keep_default_na=False)
    ex_row = df[df["primary_accession"] == "EX1"].iloc[0]
    assert "exclusion_list" in ex_row["exclusion_criteria"] or "exclusion flag" in ex_row["exclusion_criteria"]


def test_filter_and_extract_normalises_versions_and_skips_missing_fasta_seq(tmp_path):
    ref_list = tmp_path / "refs.txt"
    ref_list.write_text("NC_002204\tmaster\t1\n", encoding="utf-8")
    fa = tmp_path / "seqs.fa"
    fa.write_text(">NC_002204.1\n" + ("A" * 500) + "\n", encoding="utf-8")
    mat = tmp_path / "gb.tsv"
    mat.write_text(
        "gi_number\tprimary_accession\taccession_type\texclusion_status\texclusion_criteria\tdivision\tlength\tn\n"
        "NC_002204.1\tNC_002204.1\t\t0\t\tVRL\t500\t0\n"
        "Q_MISSING\tQ_MISSING\t\t0\t\tVRL\t500\t0\n",
        encoding="utf-8",
    )
    fae = fae_mod.FilterAndExtractSequences(
        genbank_matrix=str(mat),
        sequence_file=str(fa),
        genbank_matrix_filtered=str(tmp_path / "gb_filt.tsv"),
        ref_file=str(ref_list),
        base_dir=str(tmp_path),
        output_dir="out",
        total_length=10,
        real_length=10,
        prop_ambigious=None,
        segmented_virus="Y",
        gb_division=None,
        valid_divisions=["VRL"],
        seq_type=None,
    )
    fae.filter_columns()
    ref_out = (tmp_path / "out" / "ref_seq.fa").read_text(encoding="utf-8")
    q_out = (tmp_path / "out" / "query_seq.fa").read_text(encoding="utf-8")
    assert "NC_002204" in ref_out
    assert ">Q_MISSING\n\n" not in q_out


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-05: merge_into_gB_matrix leaves missing-FASTA rows with exclusion_status != '1' and breaks on FASTA descriptions")
def test_merge_into_gb_matrix_handles_fasta_descriptions_and_sets_exclusion_status(tmp_path):
    fa = tmp_path / "ext.fas"
    fa.write_text(">EPI999 A/duck/2024\nATGCATGCAT\n", encoding="utf-8")
    merger = mgm_mod.NormalizeAndMerge.__new__(mgm_mod.NormalizeAndMerge)
    merger.input_fasta = str(fa)
    fmap = merger.parse_fasta()
    assert "EPI999" in fmap


def test_gff_to_dictionary_preserves_standalone_cds_and_handles_fasta_section(tmp_path):
    gff = tmp_path / "mixed.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "\n"
        "NC_1\tRefSeq\tCDS\t10\t100\t.\t+\t0\tID=cds1;product=ORF1ab\n"
        "NC_1\tRefSeq\tmature_protein_region_of_CDS\t10\t50\t.\t+\t0\tParent=cds1;product=nsp1\n"
        "NC_1\tRefSeq\tCDS\t200\t400\t.\t+\t0\tID=cds2;product=Spike\n"
        "##FASTA\n"
        ">NC_1\nATGC\n",
        encoding="utf-8",
    )
    d = gff_mod.GffDictionary(str(gff)).create_gff_dict()
    products = {c["product"] for c in d.get("CDS", [])}
    assert "Spike" in products and "nsp1" in products


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-07: ValidateMatrix rejects ISO 8601 dates ('2020-03-15') and corrupts 'Nasal swab' via unanchored NA regex")
def test_validate_matrix_accepts_iso_dates_and_preserves_nasal_comment():
    assert vm_mod.ValidateMatrix._validate_date_str("2020-03-15") == "2020-03-15"


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-08: ValidateStrain fails to strip (H10N8) and silently drops ragged rows with None fields")
def test_validate_strain_cleans_two_digit_subtypes_and_preserves_ragged_rows(tmp_path):
    assert vs_mod.clean_strain_name("A/duck/Jiangxi/1/2013(H10N8)") == "A/duck/Jiangxi/1/2013"


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-09: cross_db_dedup extracts 'subgroup_a' from non-flu definition and maps bunyavirus 'M' -> '7'")
def test_cross_db_dedup_does_not_mangle_non_flu_definitions_or_bunyavirus_m_segment():
    iso = cdd_mod.extract_canonical_isolate(
        {"definition": "Human orthopneumovirus (subgroup A), complete genome", "strain": "", "isolate": "", "Parsed_strain": ""}
    )
    assert iso == ""
    assert cdd_mod.resolve_row_segment_hint({"segment": "M"}) == "M"


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-10: DedupStrainAlignment._load_protected_accessions protects exclusion_list entries")
def test_dedup_strain_alignment_does_not_protect_exclusion_list_entries(tmp_path):
    ref_list = tmp_path / "refs.txt"
    ref_list.write_text("M1\tmaster\t1\nEX1\texclusion_list\tB\n", encoding="utf-8")
    protected = dsa_mod._load_protected_accessions(str(ref_list))
    assert "M1" in protected and "EX1" not in protected


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-11: gisaid_tidy.extract_segment_epi rejects EPI_ISL_ only and drops Influenza C/D HE/P3 columns")
def test_gisaid_tidy_supports_epi_isl_for_unsegmented_gisaid():
    epi = gt_mod.extract_segment_epi(">hRSV/A/2024|EPI_ISL_123456|2024-01-01", {"EPI_ISL_123456"}, {})
    assert epi == "EPI_ISL_123456"


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-12: date_utils maps '05/12/65' to 2065 and drops full month '15-March-2020' and INSDC ranges")
def test_date_utils_handles_1900s_two_digit_years_and_full_month_names():
    c1 = du_mod.split_date_components("05/12/65")
    assert c1["year"] == 1965
    c2 = du_mod.split_date_components("15-March-2020")
    assert c2 == {"day": 15, "month": 3, "year": 2020}


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-13: ValidateRefList accepts 0 masters (is_segmented=N) and duplicate accessions")
def test_validate_ref_list_rejects_zero_masters(tmp_path, monkeypatch):
    bad = tmp_path / "no_master.txt"
    bad.write_text("R1\treference\nR2\treference\n", encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["ValidateRefList.py", "-r", str(bad), "-s", "N"])
    with pytest.raises((ValueError, SystemExit)):
        runpy.run_path(str(REPO_ROOT / "scripts/ValidateRefList.py"), run_name="__main__")


@pytest.mark.xfail(strict=True, reason="AUDIT-ING-14: accession_utils fails to strip version on WGS accessions (JAAABC010000001.1) and turns NaN into 'nan'")
def test_accession_utils_handles_wgs_and_nan():
    assert au_mod.normalise_accession("JAAABC010000001.1") == "JAAABC010000001"
    assert au_mod.normalise_accession(float("nan")) is None
