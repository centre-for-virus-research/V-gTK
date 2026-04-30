import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "ValidateRefList.py"
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "validate_ref_list_edge"


def test_validate_ref_list_non_segmented_ok(tmp_path: Path):
    ref = tmp_path / "ref.tsv"
    ref.write_text("ACC1\tmaster\nACC2\treference\nACC3\texclusion_list\n", encoding="utf-8")

    out = tmp_path / "validated.txt"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-r",
            str(ref),
            "-s",
            "N",
            "-o",
            str(out),
        ],
        check=True,
    )

    assert out.read_text(encoding="utf-8").strip() == "Reference list validation passed."


def test_validate_ref_list_segmented_missing_column_fails(tmp_path: Path):
    ref = tmp_path / "bad.tsv"
    ref.write_text("ACC1\tmaster\nACC2\treference\n", encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-r", str(ref), "-s", "Y"],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "reference list must contain three columns" in result.stderr.lower()


def test_validate_ref_list_invalid_indicator_fails(tmp_path: Path):
    ref = tmp_path / "bad_indicator.tsv"
    ref.write_text("ACC1\tmaster\t1\nACC2\tbadvalue\t1\n", encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-r", str(ref), "-s", "Y"],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "invalid values" in result.stderr.lower()


def test_validate_ref_list_exclusion_only_segment_is_allowed(tmp_path: Path):
    ref = tmp_path / "ref.tsv"
    ref.write_text(
        "ACC1\tmaster\t1\n"
        "ACC2\treference\t1\n"
        "ACC3\texclusion_list\t9\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-r", str(ref), "-s", "Y"],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "validation passed" in result.stdout.lower()


def test_validate_ref_list_edge_missing_master_in_segment_fails():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-r",
            str(DATA_DIR / "segmented_missing_master.tsv"),
            "-s",
            "Y",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "no master found for segment" in result.stderr.lower()


def test_validate_ref_list_edge_non_segmented_allows_extra_columns(tmp_path: Path):
    out_file = tmp_path / "ok.txt"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-r",
            str(DATA_DIR / "non_segmented_extra_columns.tsv"),
            "-s",
            "N",
            "-o",
            str(out_file),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert out_file.read_text(encoding="utf-8").strip() == "Reference list validation passed."


def test_validate_ref_list_edge_segmented_exclusion_only_passes():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-r",
            str(DATA_DIR / "segmented_exclusion_only_ok.tsv"),
            "-s",
            "Y",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "validation passed" in result.stdout.lower()


def test_validate_ref_list_accepts_headered_segmented_file(tmp_path: Path):
    ref = tmp_path / "headered.tsv"
    ref.write_text(
        "primary_accession\tstatus\tsegment\tgenotype\tsubtype\n"
        "NC_004102\tmaster\t1\t1\ta\n"
        "EU781827\treference\t1\t1\tb\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-r", str(ref), "-s", "Y"],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "validation passed" in result.stdout.lower()
