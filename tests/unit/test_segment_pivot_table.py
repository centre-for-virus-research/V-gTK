import csv
from pathlib import Path

import pytest

from SegmentPivotTable import (
    DEFAULT_REQUIRED_SEGMENTS as required_segments,
    elect_isolate_key_column,
    natural_segment_sort,
    normalise_segment_label,
    pivot_data,
    required_segments_from_ref_list,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURES = REPO_ROOT / "test_data" / "unit" / "segment_pivot"


def write_tsv(path: Path, rows, header):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv_dicts(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_summary(path: Path):
    return {row["key"]: row["value"] for row in read_tsv_dicts(path)}


def test_pivot_data_complete_and_incomplete(tmp_path: Path):
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"

    rows = []
    for seg in range(1, 9):
        rows.append([f"A{seg}", "StrainComplete", str(seg), ""])

    rows.extend(
        [
            ["B1", "StrainIncomplete", "1", ""],
            ["B2", "StrainIncomplete", "2.0", ""],
            ["EXCL", "StrainExcluded", "3", "manual exclusion"],
        ]
    )

    write_tsv(
        input_tsv,
        rows,
        header=["primary_accession", "Parsed_strain", "segment_validated", "exclusion"],
    )

    pivot_data(str(input_tsv), str(output_tsv), required_segments)

    out = read_tsv_dicts(output_tsv)
    by_strain = {r["Parsed_strain"]: r for r in out}

    assert by_strain["StrainComplete"]["Complete_status"] == "Complete"
    assert by_strain["StrainComplete"]["1"] == "A1"
    assert by_strain["StrainComplete"]["8"] == "A8"

    assert by_strain["StrainIncomplete"]["Complete_status"] == "Incomplete"
    assert by_strain["StrainIncomplete"]["1"] == "B1"
    assert by_strain["StrainIncomplete"]["2"] == "B2"

    assert "StrainExcluded" not in by_strain


def test_pivot_data_missing_required_columns_raises(tmp_path: Path):
    input_tsv = tmp_path / "bad.tsv"
    output_tsv = tmp_path / "pivot.tsv"

    write_tsv(
        input_tsv,
        [["A1", "StrainX", ""]],
        header=["primary_accession", "Parsed_strain", "exclusion"],
    )

    with pytest.raises(ValueError, match="missing one or more required columns"):
        pivot_data(str(input_tsv), str(output_tsv), required_segments)


def test_pivot_data_missing_every_isolate_key_raises(tmp_path: Path):
    input_tsv = tmp_path / "bad.tsv"
    output_tsv = tmp_path / "pivot.tsv"

    write_tsv(
        input_tsv,
        [["A1", "1"]],
        header=["primary_accession", "segment_validated"],
    )

    with pytest.raises(ValueError, match="missing one or more required columns"):
        pivot_data(str(input_tsv), str(output_tsv), required_segments)


def test_exclusion_column_is_optional(tmp_path: Path):
    """A matrix with no `exclusion` column pivots instead of raising."""
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"

    write_tsv(
        input_tsv,
        [["A1", "StrainX", "1"], ["A2", "StrainX", "2"]],
        header=["primary_accession", "Parsed_strain", "segment_validated"],
    )

    pivot_data(str(input_tsv), str(output_tsv), required_segments)

    out = read_tsv_dicts(output_tsv)
    assert out[0]["Parsed_strain"] == "StrainX"
    assert out[0]["Complete_status"] == "Incomplete"


def test_flu_golden_file(tmp_path: Path):
    """Bare defaults reproduce the pre-generalisation output byte for byte."""
    output_tsv = tmp_path / "pivot.tsv"

    pivot_data(str(FIXTURES / "flu_matrix.tsv"), str(output_tsv))

    assert output_tsv.read_bytes() == (FIXTURES / "flu_expected.tsv").read_bytes()


def test_required_segments_from_influenza_ref_list():
    """Master rows give 1-8; the B/C/D exclusion_list decoys stay out."""
    pytest.importorskip("pandas")
    ref_list = REPO_ROOT / "generic" / "influenza" / "ref_list_refmast.txt"

    segments, source = required_segments_from_ref_list(str(ref_list))

    assert segments == ["1", "2", "3", "4", "5", "6", "7", "8"]
    assert source == "ref_list:master"


def test_lasv_letter_segments_golden(tmp_path: Path):
    """A 2-segment virus with no Parsed_strain column pivots on `isolate`."""
    pytest.importorskip("pandas")
    output_tsv = tmp_path / "pivot.tsv"

    pivot_data(
        str(FIXTURES / "lasv_matrix.tsv"),
        str(output_tsv),
        ref_list=str(FIXTURES / "lasv_ref_list.tsv"),
    )

    assert output_tsv.read_bytes() == (FIXTURES / "lasv_expected.tsv").read_bytes()

    out = read_tsv_dicts(output_tsv)
    assert list(out[0]) == ["isolate", "L", "S", "Complete_status"]


def test_letter_segment_does_not_raise(tmp_path: Path):
    """Regression guard: `str(int(float('L')))` used to kill the process."""
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"

    write_tsv(
        input_tsv,
        [["A1", "iso1", "L", ""]],
        header=["primary_accession", "isolate", "segment_validated", "exclusion"],
    )

    summary = pivot_data(str(input_tsv), str(output_tsv), ["L", "S"])

    assert summary["n_isolates"] == 1
    assert read_tsv_dicts(output_tsv)[0]["L"] == "A1"


def test_key_election(tmp_path: Path):
    header = ["primary_accession", "Parsed_strain", "isolate", "strain", "segment_validated"]

    with_parsed = tmp_path / "with_parsed.tsv"
    write_tsv(with_parsed, [["A1", "PS", "iso", "st", "1"]], header)
    rows = read_tsv_dicts(with_parsed)
    column, coverage, coverages = elect_isolate_key_column(
        rows, ("Parsed_strain", "isolate", "strain")
    )
    assert column == "Parsed_strain"
    assert coverage == 1.0
    assert coverages["isolate"] == 1.0

    no_parsed = tmp_path / "no_parsed.tsv"
    write_tsv(
        no_parsed,
        [["A1", "iso", "st", "1"]],
        ["primary_accession", "isolate", "strain", "segment_validated"],
    )
    column, _, coverages = elect_isolate_key_column(
        read_tsv_dicts(no_parsed), ("Parsed_strain", "isolate", "strain")
    )
    assert column == "isolate"
    assert "Parsed_strain" not in coverages

    # Nothing clears min_coverage -> best-covered candidate still wins.
    sparse = tmp_path / "sparse.tsv"
    write_tsv(
        sparse,
        [["A1", "", "st", "1"], ["A2", "", "", "2"], ["A3", "", "", "3"]],
        ["primary_accession", "isolate", "strain", "segment_validated"],
    )
    column, coverage, _ = elect_isolate_key_column(
        read_tsv_dicts(sparse), ("isolate", "strain")
    )
    assert column == "strain"
    assert coverage == pytest.approx(1 / 3)

    # No candidate populated at all -> no elected column.
    empty = tmp_path / "empty.tsv"
    write_tsv(
        empty,
        [["A1", "", "", "1"]],
        ["primary_accession", "isolate", "strain", "segment_validated"],
    )
    column, coverage, _ = elect_isolate_key_column(
        read_tsv_dicts(empty), ("isolate", "strain")
    )
    assert column is None
    assert coverage == 0.0


def test_blank_key_becomes_singleton(tmp_path: Path):
    input_tsv = tmp_path / "matrix.tsv"
    header = ["primary_accession", "isolate", "segment_validated", "exclusion"]
    write_tsv(
        input_tsv,
        [
            ["A1", "iso1", "L", ""],
            ["A2", "iso1", "S", ""],
            ["B1", "", "L", ""],
            ["B2", "", "S", ""],
        ],
        header,
    )

    singleton = tmp_path / "singleton.tsv"
    pivot_data(str(input_tsv), str(singleton), ["L", "S"], extended_columns=True)
    rows = {r["isolate"]: r for r in read_tsv_dicts(singleton)}
    assert rows["iso1"]["Complete_status"] == "Complete"
    # The two blank-key records stay apart instead of fusing into a fake genome.
    assert rows["B1"]["isolate_key_source"] == "accession_fallback"
    assert rows["B1"]["Complete_status"] == "Incomplete"
    assert rows["B2"]["Complete_status"] == "Incomplete"

    legacy = tmp_path / "legacy.tsv"
    pivot_data(str(input_tsv), str(legacy), ["L", "S"], blank_isolate="group")
    legacy_rows = {r["isolate"]: r for r in read_tsv_dicts(legacy)}
    # Legacy behaviour: every blank key collapses into one bucket, which can be
    # reported Complete from unrelated records.
    assert legacy_rows[""]["Complete_status"] == "Complete"


def test_natural_segment_order(tmp_path: Path):
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"
    segments = [str(n) for n in range(1, 12)]
    write_tsv(
        input_tsv,
        [[f"A{n}", "iso1", n, ""] for n in segments],
        ["primary_accession", "isolate", "segment_validated", "exclusion"],
    )

    pivot_data(str(input_tsv), str(output_tsv), segments)

    header = list(read_tsv_dicts(output_tsv)[0])
    assert header == ["isolate"] + segments + ["Complete_status"]
    assert natural_segment_sort(["10", "2", "1", "L"]) == ["1", "2", "10", "L"]


def test_duplicate_accessions_are_joined_and_counted(tmp_path: Path):
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"
    summary_tsv = tmp_path / "pivot.summary.tsv"
    write_tsv(
        input_tsv,
        [["A1", "iso1", "L", ""], ["A2", "iso1", "S", ""], ["A6", "iso1", "S", ""]],
        ["primary_accession", "isolate", "segment_validated", "exclusion"],
    )

    summary = pivot_data(
        str(input_tsv),
        str(output_tsv),
        ["L", "S"],
        extended_columns=True,
        summary_file=str(summary_tsv),
    )

    row = read_tsv_dicts(output_tsv)[0]
    assert row["S"] == "A2,A6"
    assert row["Duplicate_segments"] == "S:2"
    assert row["Complete_status"] == "Complete"
    assert summary["n_duplicate_isolate_segment_cells"] == 1
    assert read_summary(summary_tsv)["n_duplicate_isolate_segment_cells"] == "1"


def test_unexpected_segment_is_dropped_and_counted(tmp_path: Path):
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"
    write_tsv(
        input_tsv,
        [["A1", "iso1", "L", ""], ["A2", "iso1", "X", ""]],
        ["primary_accession", "isolate", "segment_validated", "exclusion"],
    )

    summary = pivot_data(str(input_tsv), str(output_tsv), ["L", "S"])

    assert list(read_tsv_dicts(output_tsv)[0]) == ["isolate", "L", "S", "Complete_status"]
    assert summary["n_rows_unexpected_segment"] == 1
    assert summary["unexpected_segment_labels"] == "X"


def test_exclusion_zero_does_not_drop_rows(tmp_path: Path):
    """`exclusion=0` means not excluded; the old bare-truthiness test dropped it."""
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"
    write_tsv(
        input_tsv,
        [["A1", "iso1", "L", "0"], ["A2", "iso1", "S", ""]],
        ["primary_accession", "isolate", "segment_validated", "exclusion"],
    )

    summary = pivot_data(str(input_tsv), str(output_tsv), ["L", "S"])

    assert summary["n_rows_excluded"] == 0
    assert read_tsv_dicts(output_tsv)[0]["Complete_status"] == "Complete"


def test_summary_file_always_written(tmp_path: Path):
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"
    summary_tsv = tmp_path / "pivot.summary.tsv"
    write_tsv(
        input_tsv,
        [["A1", "iso1", "not found", ""]],
        ["primary_accession", "isolate", "segment_validated", "exclusion"],
    )

    pivot_data(
        str(input_tsv), str(output_tsv), ["L", "S"], summary_file=str(summary_tsv)
    )

    summary = read_summary(summary_tsv)
    assert summary["required_segments_source"] == "argument"
    assert summary["n_isolates"] == "0"
    assert read_tsv_dicts(output_tsv) == []


def test_normalise_segment_label():
    assert normalise_segment_label("2.0") == "2"
    assert normalise_segment_label(" 08 ") == "8"
    assert normalise_segment_label("L") == "L"
    assert normalise_segment_label("not found") is None
    assert normalise_segment_label("") is None
