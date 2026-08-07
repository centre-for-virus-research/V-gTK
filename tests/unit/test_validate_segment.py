import csv
from pathlib import Path

from ValidateSegment import annotate_matrix, build_reference_map, build_segment_map


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "validate_segment_edge"


def write_tsv(path: Path, rows, header=None):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        if header:
            writer.writerow(header)
        writer.writerows(rows)


def read_tsv_dicts(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_build_maps(tmp_path: Path):
    ann = tmp_path / "annotated.tsv"
    write_tsv(
        ann,
        [
            ["Q1", "REF_A", "99.9", "+", "1"],
            ["Q2", "REF_B", "99.0", "-", "9"],
        ],
    )

    assert build_segment_map(str(ann)) == {"Q1": "1", "Q2": "9"}
    assert build_reference_map(str(ann)) == {"Q1": "REF_A", "Q2": "REF_B"}


def test_annotate_matrix_preserves_existing_exclusion_by_default(tmp_path: Path):
    matrix = tmp_path / "gB_matrix.tsv"
    write_tsv(
        matrix,
        [
            ["Q1", ""],
            ["Q2", "manual exclusion"],
            ["Q3", ""],
        ],
        header=["primary_accession", "exclusion"],
    )

    ann = tmp_path / "annotated.tsv"
    write_tsv(
        ann,
        [
            ["Q1", "REF_A", "99.9", "+", "1"],
            ["Q2", "REF_B", "88.1", "-", "9"],
        ],
    )

    out_file = tmp_path / "out.tsv"
    annotate_matrix(
        str(matrix),
        build_segment_map(str(ann)),
        build_reference_map(str(ann)),
        str(out_file),
        overwrite=False,
        overwrite_exclusions=False,
    )

    rows = read_tsv_dicts(out_file)
    assert rows[0]["segment_validated"] == "1"
    assert rows[0]["closest_reference"] == "REF_A"
    assert rows[0]["exclusion"] == ""

    # Existing exclusion must be preserved when overwrite_exclusions=False
    assert rows[1]["segment_validated"] == "9"
    assert rows[1]["closest_reference"] == "REF_B"
    assert rows[1]["exclusion"] == "manual exclusion"

    # No hit -> not significant BLAST hit
    assert rows[2]["segment_validated"] == "not found"
    assert rows[2]["closest_reference"] == "not found"
    assert rows[2]["exclusion"] == "not significant BLAST hit"


def test_annotate_matrix_can_overwrite_exclusions(tmp_path: Path):
    matrix = tmp_path / "gB_matrix.tsv"
    write_tsv(
        matrix,
        [["Q1", "manual exclusion"]],
        header=["primary_accession", "exclusion"],
    )

    ann = tmp_path / "annotated.tsv"
    write_tsv(ann, [["Q1", "REF_X", "90.0", "+", "10"]])

    out_file = tmp_path / "out.tsv"
    annotate_matrix(
        str(matrix),
        build_segment_map(str(ann)),
        build_reference_map(str(ann)),
        str(out_file),
        overwrite=False,
        overwrite_exclusions=True,
    )

    rows = read_tsv_dicts(out_file)
    assert rows[0]["segment_validated"] == "10"
    assert rows[0]["closest_reference"] == "REF_X"
    assert rows[0]["exclusion"] == "non IAV genomic sequence"


def test_annotate_matrix_does_not_mark_reference_rows_as_no_hit(tmp_path: Path):
    matrix = tmp_path / "gB_matrix.tsv"
    write_tsv(
        matrix,
        [
            ["REF1", "reference", "4", ""],
            ["Q1", "query", "", ""],
        ],
        header=["primary_accession", "accession_type", "segment", "exclusion"],
    )

    ann = tmp_path / "annotated.tsv"
    write_tsv(ann, [["Q1", "REF1", "99.9", "+", "4"]])

    out_file = tmp_path / "out.tsv"
    annotate_matrix(
        str(matrix),
        build_segment_map(str(ann)),
        build_reference_map(str(ann)),
        str(out_file),
        overwrite=False,
        overwrite_exclusions=False,
    )

    rows = read_tsv_dicts(out_file)
    assert rows[0]["segment_validated"] == "4"
    assert rows[0]["closest_reference"] == "REF1"
    assert rows[0]["exclusion"] == ""

    assert rows[1]["segment_validated"] == "4"
    assert rows[1]["closest_reference"] == "REF1"
    assert rows[1]["exclusion"] == ""


def test_edge_dataset_preserve_existing_exclusions(tmp_path: Path):
    matrix = tmp_path / "input_matrix.tsv"
    matrix.write_text((DATA_DIR / "input_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    ann = DATA_DIR / "annotated_hits.tsv"
    out_file = tmp_path / "out_preserve.tsv"

    annotate_matrix(
        str(matrix),
        build_segment_map(str(ann)),
        build_reference_map(str(ann)),
        str(out_file),
        overwrite=False,
        overwrite_exclusions=False,
    )

    assert read_tsv_dicts(out_file) == read_tsv_dicts(DATA_DIR / "expected_preserve.tsv")


def test_edge_dataset_force_overwrite_exclusions(tmp_path: Path):
    matrix = tmp_path / "input_matrix.tsv"
    matrix.write_text((DATA_DIR / "input_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    ann = DATA_DIR / "annotated_hits.tsv"
    out_file = tmp_path / "out_overwrite.tsv"

    annotate_matrix(
        str(matrix),
        build_segment_map(str(ann)),
        build_reference_map(str(ann)),
        str(out_file),
        overwrite=False,
        overwrite_exclusions=True,
    )

    assert read_tsv_dicts(out_file) == read_tsv_dicts(DATA_DIR / "expected_overwrite.tsv")


def test_valid_segments_from_ref_list(tmp_path: Path):
    """With a reference list, letter segments survive instead of being excluded."""
    matrix = tmp_path / "gB_matrix.tsv"
    write_tsv(
        matrix,
        [["Q1", "", "", ""], ["Q2", "", "", ""]],
        header=["primary_accession", "segment_validated", "closest_reference", "exclusion"],
    )
    ann = tmp_path / "annotated.tsv"
    write_tsv(ann, [["Q1", "REF_L", "99.9", "+", "L"], ["Q2", "REF_X", "99.0", "+", "X"]])
    out_file = tmp_path / "out.tsv"

    annotate_matrix(
        str(matrix),
        build_segment_map(str(ann)),
        build_reference_map(str(ann)),
        str(out_file),
        overwrite=False,
        overwrite_exclusions=False,
        valid_segments={"l", "s"},
    )

    rows = {r["primary_accession"]: r for r in read_tsv_dicts(out_file)}
    assert rows["Q1"]["segment_validated"] == "L"
    assert rows["Q1"]["exclusion"] == ""
    assert rows["Q2"]["exclusion"] == "segment not in reference segment set"


def test_valid_segments_flu_equivalence(tmp_path: Path):
    """The influenza reference list partitions exactly as the hardcoded 1-8 range."""
    import pytest

    pytest.importorskip("pandas")
    from ValidateSegment import load_valid_segments

    valid = load_valid_segments(str(REPO_ROOT / "generic" / "influenza" / "ref_list_refmast.txt"))
    assert valid == {"1", "2", "3", "4", "5", "6", "7", "8"}

    matrix = tmp_path / "gB_matrix.tsv"
    write_tsv(
        matrix,
        [["Q1", "", "", ""], ["Q2", "", "", ""], ["Q3", "", "", ""], ["Q4", "", "", ""]],
        header=["primary_accession", "segment_validated", "closest_reference", "exclusion"],
    )
    ann = tmp_path / "annotated.tsv"
    write_tsv(
        ann,
        [
            ["Q1", "R", "99", "+", "1"],
            ["Q2", "R", "99", "+", "08"],
            ["Q3", "R", "99", "+", "10"],
            ["Q4", "R", "99", "+", "2.5"],
        ],
    )
    out_file = tmp_path / "out.tsv"

    annotate_matrix(
        str(matrix),
        build_segment_map(str(ann)),
        build_reference_map(str(ann)),
        str(out_file),
        overwrite=False,
        overwrite_exclusions=False,
        valid_segments=valid,
    )

    rows = {r["primary_accession"]: r for r in read_tsv_dicts(out_file)}
    assert rows["Q1"]["exclusion"] == ""
    assert rows["Q2"]["exclusion"] == ""
    assert rows["Q3"]["exclusion"] == "segment not in reference segment set"
    assert rows["Q4"]["exclusion"] == "segment not in reference segment set"
