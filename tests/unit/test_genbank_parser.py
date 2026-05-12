from pathlib import Path
from typing import Optional
import sqlite3

import pandas as pd
import pytest

from GenBankParser import GenBankParser


def _gbseq_xml(
    primary_accession: str,
    sequence: str,
    country: str = "UK:London",
    collection_date: str = "12-Jan-2020",
    accession_version: Optional[str] = None,
    segment: Optional[str] = None,
) -> str:
    accession_version = accession_version or f"{primary_accession}.1"
    segment_tag = (
        f"<GBQualifier><GBQualifier_name>segment</GBQualifier_name><GBQualifier_value>{segment}</GBQualifier_value></GBQualifier>"
        if segment is not None
        else ""
    )
    return f"""
<GBSeq>
  <GBSeq_locus>{primary_accession}</GBSeq_locus>
  <GBSeq_length>{len(sequence)}</GBSeq_length>
  <GBSeq_moltype>genomic RNA</GBSeq_moltype>
  <GBSeq_topology>linear</GBSeq_topology>
  <GBSeq_division>VRL</GBSeq_division>
  <GBSeq_update-date>01-JAN-2024</GBSeq_update-date>
  <GBSeq_create-date>01-JAN-2020</GBSeq_create-date>
  <GBSeq_definition>Example definition {primary_accession}</GBSeq_definition>
  <GBSeq_primary-accession>{primary_accession}</GBSeq_primary-accession>
  <GBSeq_accession-version>{accession_version}</GBSeq_accession-version>
  <GBSeq_source>Virus</GBSeq_source>
  <GBSeq_organism>Virus organism</GBSeq_organism>
  <GBSeq_taxonomy>Viruses; Example</GBSeq_taxonomy>
  <GBSeq_feature-table>
    <GBFeature>
      <GBFeature_key>source</GBFeature_key>
      <GBFeature_quals>
        <GBQualifier><GBQualifier_name>country</GBQualifier_name><GBQualifier_value>{country}</GBQualifier_value></GBQualifier>
        <GBQualifier><GBQualifier_name>host</GBQualifier_name><GBQualifier_value>bat</GBQualifier_value></GBQualifier>
        <GBQualifier><GBQualifier_name>collection_date</GBQualifier_name><GBQualifier_value>{collection_date}</GBQualifier_value></GBQualifier>
                {segment_tag}
      </GBFeature_quals>
    </GBFeature>
    <GBFeature>
      <GBFeature_key>gene</GBFeature_key>
      <GBFeature_location>1..3</GBFeature_location>
      <GBFeature_quals>
        <GBQualifier><GBQualifier_name>gene</GBQualifier_name><GBQualifier_value>N</GBQualifier_value></GBQualifier>
      </GBFeature_quals>
    </GBFeature>
  </GBSeq_feature-table>
  <GBSeq_references>
    <GBReference>
      <GBReference_reference>1</GBReference_reference>
      <GBReference_position>1..{len(sequence)}</GBReference_position>
      <GBReference_authors><GBAuthor>A. Author</GBAuthor></GBReference_authors>
      <GBReference_title>Test title</GBReference_title>
      <GBReference_journal>Test journal</GBReference_journal>
      <GBReference_pubmed>12345</GBReference_pubmed>
    </GBReference>
  </GBSeq_references>
  <GBSeq_sequence>{sequence.lower()}</GBSeq_sequence>
</GBSeq>
"""


def _write_xml(path: Path, gbseq_blocks):
    xml = "<GBSet>\n" + "\n".join(gbseq_blocks) + "\n</GBSet>\n"
    path.write_text(xml, encoding="utf-8")


def _write_update_db(path: Path, existing_accessions=None, excluded_accessions=None):
    existing_accessions = existing_accessions or []
    excluded_accessions = excluded_accessions or []
    conn = sqlite3.connect(str(path))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE meta_data (primary_accession TEXT, exclusion_status TEXT, exclusion_criteria TEXT)")
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, exclusion_status, exclusion_criteria) VALUES (?, '', '')",
            [(acc,) for acc in existing_accessions],
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, exclusion_status, exclusion_criteria) VALUES (?, '1', 'excluded')",
            [(acc,) for acc in excluded_accessions],
        )
        conn.commit()
    finally:
        conn.close()


def test_load_ref_list_supports_two_and_three_column_lines(tmp_path: Path):
    ref_file = tmp_path / "ref_list.tsv"
    ref_file.write_text("REF1\tmaster\t1\nREF2\treference\nBADLINE\n", encoding="utf-8")

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="N",
    )

    refs = parser.load_ref_list(str(ref_file))
    assert refs == {"REF1": "master", "REF2": "reference"}


def test_load_ref_list_supports_headered_lines(tmp_path: Path):
    ref_file = tmp_path / "ref_list.tsv"
    ref_file.write_text(
        "primary_accession\tstatus\tsegment\tgenotype\tsubtype\nREF1\tmaster\t1\t1\ta\nREF2\treference\t1\t1\tb\n",
        encoding="utf-8",
    )

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="Y",
    )

    refs = parser.load_ref_list(str(ref_file))
    assert refs == {"REF1": "master", "REF2": "reference"}


def test_xml_to_tsv_parses_and_applies_exclusion(tmp_path: Path):
    xml_file = tmp_path / "batch-1.xml"
    _write_xml(
        xml_file,
        [
            _gbseq_xml("REF1", "ATGCNN"),
            _gbseq_xml("Q1", "AATTCG", country="France:Paris", collection_date="2021"),
        ],
    )

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=None,
        exclusion_list=None,
        is_segmented_virus="N",
    )

    rows = parser.xml_to_tsv(str(xml_file), {"REF1": "master"}, ["Q1"])
    by_acc = {r["primary_accession"]: r for r in rows}

    assert by_acc["REF1"]["accession_type"] == "master"
    assert by_acc["Q1"]["accession_type"] == "excluded"
    assert by_acc["Q1"]["exclusion_status"] == "1"
    assert by_acc["Q1"]["country"] == "France"
    assert by_acc["Q1"]["geo_loc"] == "Paris"
    assert by_acc["REF1"]["n"] == 2


def test_xml_to_tsv_non_segmented_forces_segment_one(tmp_path: Path):
    xml_file = tmp_path / "batch-1.xml"
    _write_xml(
        xml_file,
        [
            _gbseq_xml("Q1", "AATTCG", segment="NS3"),
        ],
    )

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=None,
        exclusion_list=None,
        is_segmented_virus="N",
    )

    rows = parser.xml_to_tsv(str(xml_file), {}, [])
    assert rows[0]["segment"] == "1"


def test_select_xml_files_raises_on_missing_critical_refs(tmp_path: Path):
    xml_file = tmp_path / "batch-1.xml"
    _write_xml(xml_file, [_gbseq_xml("REF1", "ATGC")])

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=None,
        exclusion_list=None,
        is_segmented_virus="N",
        test_run=True,
        test_xml_limit=2,
        require_refs=True,
    )

    with pytest.raises(ValueError, match="Missing CRITICAL reference accessions"):
        parser.select_xml_files({"REF1": "master", "REF_MISSING": "reference"}, ["batch-1.xml"])


def test_select_xml_files_update_mode_allows_missing_critical_refs_if_in_db(tmp_path: Path):
    xml_file = tmp_path / "batch-1.xml"
    _write_xml(xml_file, [_gbseq_xml("REF1", "ATGC")])

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=None,
        exclusion_list=None,
        is_segmented_virus="N",
        test_run=True,
        test_xml_limit=2,
        require_refs=True,
        update="dummy.db",
    )
    parser._existing_accessions = {"REF_MISSING"}

    selected = parser.select_xml_files({"REF1": "master", "REF_MISSING": "reference"}, ["batch-1.xml"])
    assert selected == ["batch-1.xml"]


def test_process_writes_matrix_and_fasta(tmp_path: Path):
    in_dir = tmp_path / "xml"
    in_dir.mkdir(parents=True, exist_ok=True)

    _write_xml(in_dir / "batch-1.xml", [_gbseq_xml("REF1", "ATGC")])
    _write_xml(in_dir / "batch-2.xml", [_gbseq_xml("Q1", "AATT")])

    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\n", encoding="utf-8")

    parser = GenBankParser(
        input_dir=str(in_dir),
        base_dir=str(tmp_path),
        output_dir="GenBank-matrix",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="N",
    )
    parser.process()

    matrix_path = tmp_path / "GenBank-matrix" / "gB_matrix_raw.tsv"
    fasta_path = tmp_path / "GenBank-matrix" / "sequences.fa"

    assert matrix_path.exists()
    assert fasta_path.exists()

    df = pd.read_csv(matrix_path, sep="\t")
    assert set(df["primary_accession"].tolist()) == {"REF1", "Q1"}
    assert "sequence" not in df.columns

    fasta_text = fasta_path.read_text(encoding="utf-8")
    assert ">REF1" in fasta_text
    assert ">Q1" in fasta_text


def test_load_existing_accessions_from_db_requires_meta_data_table(tmp_path: Path):
    bad_db = tmp_path / "bad.db"
    conn = sqlite3.connect(str(bad_db))
    try:
        conn.execute("CREATE TABLE other_table (id TEXT)")
        conn.commit()
    finally:
        conn.close()

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=None,
        exclusion_list=None,
        is_segmented_virus="N",
        update=str(bad_db),
    )

    with pytest.raises(ValueError, match="meta_data"):
        parser.load_existing_accessions_from_db(str(bad_db))


def test_process_update_mode_filters_by_existing_and_excluded_accessions(tmp_path: Path):
    xml_dir = tmp_path / "GenBank-XML"
    xml_dir.mkdir(parents=True, exist_ok=True)

    _write_xml(
        xml_dir / "batch-1.xml",
        [
            _gbseq_xml("REF1", "ATGC"),
            _gbseq_xml("Q1", "AATT"),
            _gbseq_xml("NEW1", "ATAT"),
        ],
    )

    update_db = tmp_path / "previous.db"
    _write_update_db(update_db, existing_accessions=["REF1"], excluded_accessions=["Q1"])

    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\n", encoding="utf-8")

    parser = GenBankParser(
        input_dir=None,
        base_dir=str(tmp_path),
        output_dir="GenBank-matrix",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="N",
        update=str(update_db),
    )
    parser.process()

    matrix_path = tmp_path / "GenBank-matrix" / "gB_matrix_raw.tsv"
    fasta_path = tmp_path / "GenBank-matrix" / "sequences.fa"

    assert matrix_path.exists()
    assert fasta_path.exists()

    df = pd.read_csv(matrix_path, sep="\t")
    assert df["primary_accession"].tolist() == ["REF1", "NEW1"]

    fasta_text = fasta_path.read_text(encoding="utf-8")
    assert ">NEW1" in fasta_text
    assert ">Q1" not in fasta_text
    assert ">REF1" in fasta_text


def test_process_update_mode_keeps_reference_rows(tmp_path: Path):
    xml_dir = tmp_path / "GenBank-XML"
    xml_dir.mkdir(parents=True, exist_ok=True)

    _write_xml(
        xml_dir / "batch-1.xml",
        [
            _gbseq_xml("REF1", "ATGC"),
            _gbseq_xml("Q1", "AATT"),
            _gbseq_xml("NEW1", "ATAT"),
        ],
    )

    update_db = tmp_path / "previous.db"
    _write_update_db(update_db, existing_accessions=["REF1", "Q1"], excluded_accessions=[])

    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\n", encoding="utf-8")

    parser = GenBankParser(
        input_dir=None,
        base_dir=str(tmp_path),
        output_dir="GenBank-matrix",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="N",
        update=str(update_db),
    )
    parser.process()

    df = pd.read_csv(tmp_path / "GenBank-matrix" / "gB_matrix_raw.tsv", sep="\t")
    assert set(df["primary_accession"].tolist()) == {"REF1", "NEW1"}


def test_process_update_mode_keeps_same_accession_across_segments(tmp_path: Path):
    xml_dir = tmp_path / "GenBank-XML"
    xml_dir.mkdir(parents=True, exist_ok=True)

    _write_xml(
        xml_dir / "batch-1.xml",
        [
            _gbseq_xml("SEG_ACC", "ATGC", segment="1"),
            _gbseq_xml("SEG_ACC", "ATGA", segment="2"),
        ],
    )

    update_db = tmp_path / "previous.db"
    _write_update_db(update_db, existing_accessions=[], excluded_accessions=[])

    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\n", encoding="utf-8")

    parser = GenBankParser(
        input_dir=None,
        base_dir=str(tmp_path),
        output_dir="GenBank-matrix",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="Y",
        update=str(update_db),
    )
    parser.process()

    df = pd.read_csv(tmp_path / "GenBank-matrix" / "gB_matrix_raw.tsv", sep="\t")
    seg_rows = df[df["primary_accession"] == "SEG_ACC"]
    assert sorted(seg_rows["segment"].astype(str).tolist()) == ["1", "2"]


def test_process_update_mode_require_refs_uses_db_for_missing_reference(tmp_path: Path):
    xml_dir = tmp_path / "GenBank-XML"
    xml_dir.mkdir(parents=True, exist_ok=True)

    _write_xml(
        xml_dir / "batch-1.xml",
        [
            _gbseq_xml("NEW1", "ATAT"),
        ],
    )

    update_db = tmp_path / "previous.db"
    _write_update_db(update_db, existing_accessions=["REF1", "REF_MISSING"], excluded_accessions=[])

    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\nREF_MISSING\treference\n", encoding="utf-8")

    parser = GenBankParser(
        input_dir=None,
        base_dir=str(tmp_path),
        output_dir="GenBank-matrix",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="N",
        require_refs=True,
        update=str(update_db),
    )
    parser.process()

    df = pd.read_csv(tmp_path / "GenBank-matrix" / "gB_matrix_raw.tsv", sep="\t")
    assert df["primary_accession"].tolist() == ["NEW1"]
