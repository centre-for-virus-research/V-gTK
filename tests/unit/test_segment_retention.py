"""Tests for the per-segment retention guard.

This is the check that would have caught the HA/NA subtype-loss defect on its
first run. The real numbers it is modelled on (db_summary.txt, 26 Aug 2026):

    segment 1: 462,975 passed /    925 failed   ratio 0.998
    segment 4: 205,397 passed / 367,894 failed  ratio 0.358   <- HA
    segment 6: 248,412 passed / 268,067 failed  ratio 0.481   <- NA
    segment 7: 492,067 passed /  1,044 failed   ratio 0.998

Every healthy segment is above 0.998; the two broken ones are nowhere near.
A 0.98 threshold separates them with a wide margin.
"""

import subprocess
import sys
from pathlib import Path

import pytest

from CheckSegmentRetention import check, segment_of_merged_file

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "CheckSegmentRetention.py"


def build_inputs(tmp_path: Path, retained_per_segment, assigned_per_segment,
                 excused=(), references=(("REF4", "4"), ("REF6", "6"))):
    """Write a merged MSA per segment plus the BLAST tophit and ref-segment files."""
    ref_seg = tmp_path / "ref_seq_seg.tsv"
    ref_seg.write_text("".join(f"{acc}\t{seg}\n" for acc, seg in references), encoding="utf-8")

    merged_files = []
    for segment, accs in sorted(retained_per_segment.items()):
        path = tmp_path / f"refset_{segment}_aln_merged_MSA.fasta"
        # PadAlignment prepends the backbone reference row to every block, so the
        # reference legitimately appears here and must not be counted as a query.
        ref_acc = next(a for a, s in references if s == segment)
        body = f">{ref_acc}\nACGT\n" + "".join(f">{a}\nACGT\n" for a in accs)
        path.write_text(body, encoding="utf-8")
        merged_files.append(str(path))

    tophit = tmp_path / "query_uniq_tophit_annotated.tsv"
    rows = []
    for segment, accs in sorted(assigned_per_segment.items()):
        ref_acc = next(a for a, s in references if s == segment)
        for a in accs:
            rows.append(f"{a}\t{ref_acc}\t99.0\tplus\t{segment}\n")
    tophit.write_text("".join(rows), encoding="utf-8")

    filtered_ids = tmp_path / "filtered_sequences_ids.txt"
    filtered_ids.write_text("".join(f"{a}\n" for a in excused), encoding="utf-8")

    return merged_files, str(tophit), str(ref_seg), str(filtered_ids)


def test_healthy_segments_pass(tmp_path):
    accs = [f"Q{i}" for i in range(100)]
    merged, tophit, ref_seg, filtered = build_inputs(
        tmp_path, {"4": accs, "6": accs}, {"4": accs, "6": accs})
    assert check(merged, tophit, ref_seg, filtered, min_retention=0.98) == []


def test_the_real_defect_is_caught(tmp_path):
    """Segment 4 keeps 36% of what BLAST assigned to it - the actual failure."""
    assigned = [f"Q{i}" for i in range(1000)]
    retained = assigned[:358]
    merged, tophit, ref_seg, filtered = build_inputs(
        tmp_path, {"4": retained}, {"4": assigned})

    failures = check(merged, tophit, ref_seg, filtered, min_retention=0.98)

    assert len(failures) == 1
    segment, got, expected, ratio, _examples = failures[0]
    assert segment == "4"
    assert (got, expected) == (358, 1000)
    assert ratio == pytest.approx(0.358)


def test_genuine_nextalign_failures_are_excused(tmp_path):
    """The ~6,300 sequences nextalign truly could not align are legitimate losses
    and must not trip the guard."""
    assigned = [f"Q{i}" for i in range(100)]
    excused = assigned[:5]
    retained = assigned[5:]
    merged, tophit, ref_seg, filtered = build_inputs(
        tmp_path, {"4": retained}, {"4": assigned}, excused=excused)
    assert check(merged, tophit, ref_seg, filtered, min_retention=0.98) == []


def test_one_bad_segment_does_not_mask_the_others(tmp_path):
    """Segments are reported independently - the real build had 6 good, 2 bad."""
    good = [f"G{i}" for i in range(100)]
    assigned6 = [f"B{i}" for i in range(100)]
    merged, tophit, ref_seg, filtered = build_inputs(
        tmp_path, {"4": good, "6": assigned6[:48]}, {"4": good, "6": assigned6})

    failures = check(merged, tophit, ref_seg, filtered, min_retention=0.98)
    assert [f[0] for f in failures] == ["6"]


def test_failure_lists_missing_accessions(tmp_path):
    """The operator needs to know which sequences vanished, not just that some did."""
    assigned = [f"Q{i}" for i in range(10)]
    merged, tophit, ref_seg, filtered = build_inputs(
        tmp_path, {"4": assigned[:2]}, {"4": assigned})
    failures = check(merged, tophit, ref_seg, filtered, min_retention=0.98)
    assert failures[0][4], "expected example missing accessions"
    assert all(a.startswith("Q") for a in failures[0][4])


@pytest.mark.parametrize("name,expected", [
    ("refset_4_aln_merged_MSA.fasta", "4"),
    ("/a/b/refset_12_aln_merged_MSA.fasta", "12"),
    ("NC_001542.aligned_merged_MSA.fasta", None),
])
def test_segment_inferred_from_filename(name, expected):
    assert segment_of_merged_file(name) == expected


def test_cli_exits_nonzero_on_failure(tmp_path):
    assigned = [f"Q{i}" for i in range(1000)]
    merged, tophit, ref_seg, filtered = build_inputs(
        tmp_path, {"4": assigned[:358]}, {"4": assigned})

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-m", *merged, "-t", tophit,
         "-s", ref_seg, "-f", filtered, "--min_retention", "0.98"],
        capture_output=True, text=True,
    )
    assert result.returncode == 4, result.stdout + result.stderr
    assert "retained 358/1000" in result.stderr


def test_cli_exits_zero_when_healthy(tmp_path):
    assigned = [f"Q{i}" for i in range(100)]
    merged, tophit, ref_seg, filtered = build_inputs(
        tmp_path, {"4": assigned}, {"4": assigned})

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-m", *merged, "-t", tophit,
         "-s", ref_seg, "-f", filtered],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "All segments passed" in result.stdout
