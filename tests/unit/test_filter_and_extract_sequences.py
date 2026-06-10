import csv
import sqlite3
from pathlib import Path
import io
from contextlib import redirect_stdout

from FilterAndExtractSequences import FilterAndExtractSequences


def _write_tsv(path: Path, rows, header):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def _read_tsv(path: Path):
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_filter_columns_writes_query_ref_and_exclusions(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [
            ["VRL", "REF1", "4", "0", "0", ""],
            ["VRL", "Q1", "4", "0", "0", ""],
            ["VRL", "Q2", "3", "0", "0", ""],
        ],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )

    seqs.write_text(
        ">REF1\nATGC\n>Q1\nAATT\n>Q2\nAAA\n",
        encoding="utf-8",
    )
    refs.write_text("REF1\treference\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=4,
        real_length=4,
        prop_ambigious=None,
        segmented_virus="N",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )
    processor.process()

    query_fa = (tmp_path / "Sequences" / "query_seq.fa").read_text(encoding="utf-8")
    ref_fa = (tmp_path / "Sequences" / "ref_seq.fa").read_text(encoding="utf-8")

    assert ">Q1" in query_fa
    assert ">Q2" not in query_fa
    assert ">REF1" in ref_fa

    rows = _read_tsv(matrix)
    by_acc = {r["gi_number"]: r for r in rows}
    assert by_acc["Q2"]["exclusion_status"] == "1"


def test_segmented_mode_writes_exclusion_refs_file(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [["VRL", "REF1", "4", "0", "0", ""]],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )
    seqs.write_text(">REF1\nATGC\n", encoding="utf-8")
    refs.write_text("REF1\tmaster\t1\nEXCL1\texclusion_list\t1\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=1,
        real_length=1,
        prop_ambigious=None,
        segmented_virus="Y",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )
    processor.process()

    exclusion_refs = (tmp_path / "Sequences" / "exclusion_refs.txt").read_text(encoding="utf-8").splitlines()
    assert exclusion_refs == ["EXCL1"]


def test_exclusion_list_reference_is_marked_excluded_in_matrix(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [
            ["VRL", "REF1", "4", "0", "0", ""],
            ["VRL", "EXCL1", "4", "0", "0", ""],
            ["VRL", "Q1", "4", "0", "0", ""],
        ],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )
    seqs.write_text(">REF1\nATGC\n>EXCL1\nATGC\n>Q1\nAATT\n", encoding="utf-8")
    refs.write_text("REF1\tmaster\t1\nEXCL1\texclusion_list\t1\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=1,
        real_length=1,
        prop_ambigious=None,
        segmented_virus="Y",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )
    processor.process()

    rows = _read_tsv(matrix)
    by_acc = {r["gi_number"]: r for r in rows}
    assert by_acc["EXCL1"]["exclusion_status"] == "1"
    assert by_acc["EXCL1"]["exclusion_criteria"] == "excluded by reference list exclusion flag"

    ref_fa = (tmp_path / "Sequences" / "ref_seq.fa").read_text(encoding="utf-8")
    assert ">REF1" in ref_fa
    assert ">EXCL1" not in ref_fa


def test_existing_exclusion_status_is_preserved(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [["VRL", "QX", "4", "0", "1", "already excluded"]],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )
    seqs.write_text(">QX\nATGC\n", encoding="utf-8")
    refs.write_text("REF1\treference\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=1,
        real_length=1,
        prop_ambigious=None,
        segmented_virus="N",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )
    processor.process()

    rows = _read_tsv(matrix)
    assert rows[0]["exclusion_status"] == "1"
    assert rows[0]["exclusion_criteria"] == "already excluded"


def test_update_db_reference_source_overrides_ref_file(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"
    update_db = tmp_path / "update.db"

    with sqlite3.connect(str(update_db)) as conn:
        conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
        conn.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
            [("REF_DB", "master", "1"), ("Q1", "query", "1")],
        )
        conn.commit()

    _write_tsv(
        matrix,
        [["VRL", "REF_DB", "4", "0", "0", ""], ["VRL", "Q1", "4", "0", "0", ""]],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )
    seqs.write_text(">REF_DB\nATGC\n>Q1\nAATT\n", encoding="utf-8")
    refs.write_text("REF_FILE\treference\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=1,
        real_length=1,
        prop_ambigious=None,
        segmented_virus="N",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
        update_db=str(update_db),
    )
    processor.process()

    ref_fa = (tmp_path / "Sequences" / "ref_seq.fa").read_text(encoding="utf-8")
    query_fa = (tmp_path / "Sequences" / "query_seq.fa").read_text(encoding="utf-8")
    assert ">REF_DB" in ref_fa
    assert ">Q1" in query_fa


def test_filter_columns_accepts_headered_ref_list(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [["VRL", "REF1", "4", "0", "0", ""], ["VRL", "Q1", "4", "0", "0", ""]],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )
    seqs.write_text(">REF1\nATGC\n>Q1\nAATT\n", encoding="utf-8")
    refs.write_text(
        "primary_accession\tstatus\tsegment\tgenotype\nREF1\treference\t1\t1\n",
        encoding="utf-8",
    )

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=1,
        real_length=1,
        prop_ambigious=None,
        segmented_virus="N",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )
    processor.process()

    ref_fa = (tmp_path / "Sequences" / "ref_seq.fa").read_text(encoding="utf-8")
    assert ">REF1" in ref_fa


def test_filter_columns_warns_when_all_queries_filtered(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [
            ["VRL", "REF1", "1000", "0", "0", "", "reference"],
            ["VRL", "Q1", "100", "0", "0", "sequence_too_short: real_length=100 < 500", "query"],
            ["PAT", "Q2", "1000", "0", "0", "GenBank division PAT not in valid division list", "query"],
        ],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria", "accession_type"],
    )
    seqs.write_text(">REF1\nATGC\n>Q1\nATGC\n>Q2\nATGC\n", encoding="utf-8")
    refs.write_text("REF1\treference\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=500,
        real_length=500,
        prop_ambigious=None,
        segmented_virus="N",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )

    stdout = io.StringIO()
    with redirect_stdout(stdout):
        processor.process()

    output = stdout.getvalue()
    assert "WARNING: all 2 query rows were filtered before BLAST; query_seq.fa is empty" in output
    assert "[FilterAndExtractSequences] Query Exclusion Reasons:" in output
    assert "GenBank division PAT not in valid division list" in output
    assert "sequence_too_short: real_length=100 < 500" in output

    summary = (tmp_path / "Sequences" / "filter_summary.txt").read_text(encoding="utf-8")
    assert "query rows written to query_seq.fa: 0" in summary
    assert "query rows excluded before BLAST: 2" in summary
