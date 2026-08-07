# pyright: reportMissingImports=false
import csv
import subprocess
import sys
from pathlib import Path

import pandas as pd

from merge_into_gB_matrix import NormalizeAndMerge


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "merge_into_gB_matrix.py"

BASE_COLUMNS = [
    "primary_accession",
    "accession_version",
    "locus",
    "gi_number",
    "sequence",
    "real_length",
    "exclusion",
    "country",
    "host",
]


def read_tsv_as_dicts(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, header, rows):
    lines = ["\t".join(header)]
    lines += ["\t".join(row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def make_base_matrix(tmp_path: Path) -> Path:
    return write_tsv(
        tmp_path / "gB_matrix.tsv",
        BASE_COLUMNS,
        [
            [
                "NC_001542",
                "NC_001542.1",
                "NC_001542",
                "12345",
                "ATGCATGC",
                "8",
                "",
                "United Kingdom",
                "Canis lupus",
            ]
        ],
    )


def make_mapping(tmp_path: Path) -> Path:
    return write_tsv(
        tmp_path / "column_mapping.tsv",
        ["Collection_Date", "collection_date"],
        [
            ["Location", "country"],
            ["Host_Species", "host"],
            # Blank and short rows must be tolerated.
            [""],
            ["   ", "ignored"],
        ],
    )


def make_input_tsv(tmp_path: Path, rows=None) -> Path:
    header = ["Segment_Id", "Collection_Date", "Location", "Host_Species", "Extra_Junk"]
    rows = rows if rows is not None else [
        ["EPI001", "2020-01-01", "Kenya", "Bos taurus", "junk1"],
        ["EPI002", "2021-02-02", "Tanzania", "Canis lupus", "junk2"],
    ]
    return write_tsv(tmp_path / "metadata.tsv", header, rows)


def make_fasta(tmp_path: Path, entries=None) -> Path:
    entries = entries if entries is not None else {"EPI001": ["ATGCNN", "atgc"]}
    lines = []
    for name, chunks in entries.items():
        lines.append(f">{name}")
        lines.extend(chunks)
    (tmp_path / "all_nuc.fas").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return tmp_path / "all_nuc.fas"


def make_job(tmp_path: Path, **overrides) -> NormalizeAndMerge:
    """Build a job, only generating the default inputs that were not overridden."""
    defaults = {
        "gb_matrix": lambda: str(make_base_matrix(tmp_path)),
        "input_tsv": lambda: str(make_input_tsv(tmp_path)),
        "input_fasta": lambda: str(make_fasta(tmp_path)),
        "mapping_file": lambda: str(make_mapping(tmp_path)),
        "output": lambda: str(tmp_path / "out" / "gB_matrix_merged.tsv"),
        "key": lambda: "Segment_Id",
        "dataset_source": lambda: "gisaid",
        "log_file": lambda: str(tmp_path / "gisaid_processing.log"),
    }
    kwargs = {
        name: overrides[name] if name in overrides else build()
        for name, build in defaults.items()
    }
    kwargs.update(overrides)
    return NormalizeAndMerge(**kwargs)


def test_read_mapping_skips_blank_and_incomplete_rows(tmp_path: Path):
    job = make_job(tmp_path)

    assert job.read_mapping() == {
        "Collection_Date": "collection_date",
        "Location": "country",
        "Host_Species": "host",
    }


def test_parse_fasta_joins_wrapped_lines_and_ignores_blanks(tmp_path: Path):
    fasta = tmp_path / "seqs.fa"
    fasta.write_text(
        ">EPI001 \nATGC\n\natgc\n>EPI002\nTTTT\n",
        encoding="utf-8",
    )
    job = make_job(tmp_path, input_fasta=str(fasta))

    assert job.parse_fasta() == {"EPI001": "ATGCatgc", "EPI002": "TTTT"}


def test_parse_fasta_returns_empty_without_fasta(tmp_path: Path):
    job = make_job(tmp_path, input_fasta=None)

    assert job.parse_fasta() == {}


def test_parse_fasta_ignores_a_file_with_no_records(tmp_path: Path):
    fasta = tmp_path / "headerless.fa"
    fasta.write_text("ATGC\nATGC\n", encoding="utf-8")
    job = make_job(tmp_path, input_fasta=str(fasta))

    assert job.parse_fasta() == {}


def test_write_duplicates_report_skips_empty_frames(tmp_path: Path):
    job = make_job(tmp_path)

    assert job.write_duplicates_report_df(None) is None
    assert job.write_duplicates_report_df(pd.DataFrame()) is None
    assert not (tmp_path / "out" / "merge_matrix_duplicates.tsv").exists()


def test_count_real_length_counts_acgt_only(tmp_path: Path):
    job = make_job(tmp_path)

    assert job.count_real_length("ATGCNNatgc") == 8
    assert job.count_real_length("") == 0
    assert job.count_real_length(None) == 0
    assert job.count_real_length(float("nan")) == 0


def test_collapse_duplicate_columns_is_a_no_op_without_duplicates(tmp_path: Path):
    job = make_job(tmp_path)
    df = pd.DataFrame({"a": ["1"], "b": ["2"]})

    assert job.collapse_duplicate_columns(df) is df


def test_collapse_duplicate_columns_coalesces_left_to_right(tmp_path: Path):
    job = make_job(tmp_path)
    df = pd.DataFrame(
        [["", "second", "keep"], ["first", "ignored", "keep2"]],
        columns=["country", "country", "host"],
    )

    collapsed = job.collapse_duplicate_columns(df)

    assert list(collapsed.columns) == ["country", "host"]
    assert collapsed["country"].tolist() == ["second", "first"]
    assert collapsed["host"].tolist() == ["keep", "keep2"]


def test_process_appends_normalized_rows(tmp_path: Path):
    out = tmp_path / "out" / "gB_matrix_merged.tsv"
    make_job(tmp_path, output=str(out)).process()

    rows = read_tsv_as_dicts(out)
    assert len(rows) == 3

    base_row, first, second = rows
    assert base_row["primary_accession"] == "NC_001542"
    assert base_row["dataset_source"] == "genbank"

    # Key is copied into every gB accession field.
    for field in ("primary_accession", "accession_version", "locus", "gi_number"):
        assert first[field] == "EPI001"
    assert first["dataset_source"] == "gisaid"
    assert first["Segment_Id"] == "EPI001"

    # Mapped columns are renamed onto the base schema.
    assert first["country"] == "Kenya"
    assert first["host"] == "Bos taurus"
    assert first["collection_date"] == "2020-01-01"

    # FASTA sequence is attached and ambiguous bases excluded from real_length.
    assert first["sequence"] == "ATGCNNatgc"
    assert first["real_length"] == "8"
    assert first["exclusion"] == ""

    # No FASTA match for EPI002.
    assert second["sequence"] == ""
    assert second["real_length"] == "0"
    assert second["exclusion"] == "missing fasta sequence"

    # Unmapped input columns are retained after the base columns by default.
    header = list(pd.read_csv(out, sep="\t", dtype=str, nrows=0).columns)
    assert header[: len(BASE_COLUMNS)] == BASE_COLUMNS
    assert "Extra_Junk" in header
    assert first["Extra_Junk"] == "junk1"


def test_process_reports_duplicates_and_missing_keys(tmp_path: Path):
    input_tsv = make_input_tsv(
        tmp_path,
        rows=[
            ["EPI001", "2020-01-01", "Kenya", "Bos taurus", "first"],
            ["EPI001", "2020-01-01", "Kenya", "Bos taurus", "dupe"],
            ["", "2022-03-03", "Uganda", "Canis lupus", "nokey"],
        ],
    )
    out = tmp_path / "out" / "gB_matrix_merged.tsv"
    make_job(tmp_path, input_tsv=str(input_tsv), output=str(out)).process()

    rows = read_tsv_as_dicts(out)
    assert [r["primary_accession"] for r in rows] == ["NC_001542", "EPI001"]
    assert rows[1]["Extra_Junk"] == "first"

    dup_rows = read_tsv_as_dicts(out.parent / "merge_matrix_duplicates.tsv")
    assert [(r["reason"], r["row_key"], r["Extra_Junk"]) for r in dup_rows] == [
        ("missing-key", "", "nokey"),
        ("duplicate", "EPI001", "dupe"),
    ]


def test_process_without_duplicates_writes_no_report(tmp_path: Path):
    out = tmp_path / "out" / "gB_matrix_merged.tsv"
    make_job(tmp_path, output=str(out)).process()

    assert not (out.parent / "merge_matrix_duplicates.tsv").exists()


def test_process_passes_base_through_when_key_column_missing(tmp_path: Path):
    out = tmp_path / "out" / "gB_matrix_merged.tsv"
    log_file = tmp_path / "run.log"
    make_job(tmp_path, output=str(out), key="Missing_Id", log_file=str(log_file)).process()

    rows = read_tsv_as_dicts(out)
    assert len(rows) == 1
    assert rows[0]["primary_accession"] == "NC_001542"
    assert rows[0]["dataset_source"] == "genbank"
    assert "Key column 'Missing_Id' not found" in log_file.read_text(encoding="utf-8")


def test_process_drop_unmapped_removes_extra_input_columns(tmp_path: Path):
    out = tmp_path / "out" / "gB_matrix_merged.tsv"
    make_job(tmp_path, output=str(out), drop_unmapped=True).process()

    header = list(pd.read_csv(out, sep="\t", dtype=str, nrows=0).columns)
    assert "Extra_Junk" not in header
    assert header[: len(BASE_COLUMNS)] == BASE_COLUMNS
    # Mapping targets, required columns and the key survive the drop.
    for column in ("country", "host", "collection_date", "dataset_source", "Segment_Id"):
        assert column in header


def test_process_collapses_columns_colliding_with_mapping_targets(tmp_path: Path):
    # `Location` maps onto `country`, which already exists in the input TSV.
    input_tsv = write_tsv(
        tmp_path / "metadata.tsv",
        ["Segment_Id", "country", "Location"],
        [
            ["EPI001", "", "Kenya"],
            ["EPI002", "Malawi", "Tanzania"],
        ],
    )
    out = tmp_path / "out" / "gB_matrix_merged.tsv"
    log_file = tmp_path / "run.log"
    make_job(
        tmp_path, input_tsv=str(input_tsv), output=str(out), log_file=str(log_file)
    ).process()

    rows = read_tsv_as_dicts(out)
    assert [r["country"] for r in rows[1:]] == ["Kenya", "Malawi"]
    assert "Duplicate column labels after mapping" in log_file.read_text(encoding="utf-8")


def test_process_without_mapping_keeps_source_column_names(tmp_path: Path):
    empty_mapping = tmp_path / "empty_mapping.tsv"
    empty_mapping.write_text("", encoding="utf-8")
    out = tmp_path / "out" / "gB_matrix_merged.tsv"
    make_job(tmp_path, mapping_file=str(empty_mapping), output=str(out)).process()

    header = list(pd.read_csv(out, sep="\t", dtype=str, nrows=0).columns)
    assert "Location" in header and "Host_Species" in header
    assert "collection_date" not in header

    rows = read_tsv_as_dicts(out)
    assert rows[1]["Location"] == "Kenya"
    # The unmapped base column stays empty for the appended rows.
    assert rows[1]["country"] == ""


def test_process_creates_missing_output_directory(tmp_path: Path):
    out = tmp_path / "nested" / "deeper" / "gB_matrix_merged.tsv"
    make_job(tmp_path, output=str(out)).process()

    assert out.exists()


def test_cli_matches_nextflow_invocation(tmp_path: Path):
    base = make_base_matrix(tmp_path)
    input_tsv = make_input_tsv(tmp_path)
    fasta = make_fasta(tmp_path)
    mapping = make_mapping(tmp_path)
    out = tmp_path / "gB_matrix_merged.tsv"

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-g", str(base),
            "-t", str(input_tsv),
            "-f", str(fasta),
            "-o", str(out),
            "-k", "Segment_Id",
            "--dataset_source", "gisaid",
            "-m", str(mapping),
            "--log_file", str(tmp_path / "gisaid_processing.log"),
        ],
        check=True,
        cwd=tmp_path,
    )

    rows = read_tsv_as_dicts(out)
    assert [r["primary_accession"] for r in rows] == ["NC_001542", "EPI001", "EPI002"]
    assert [r["dataset_source"] for r in rows] == ["genbank", "gisaid", "gisaid"]
