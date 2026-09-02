import csv
import subprocess
import sys
from pathlib import Path

import pytest

from CollectFilteredSequences import (
    collect_filtered_sequences,
    write_filtered_ids_only,
    write_filtered_list,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "CollectFilteredSequences.py"
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "collect_filtered_sequences_edge"


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["seqName", "errors", "warnings"])
        writer.writerows(rows)


def test_collect_filtered_sequences_only_keeps_query_alignment_rows_with_errors(tmp_path: Path):
    nextalign_dir = tmp_path / "Nextalign"

    write_csv(
        nextalign_dir / "query_aln" / "REF_A" / "chunk.errors.csv",
        [
            ["SEQ_1", "frame shift", ""],
            ["SEQ_2", "", "minor issue"],
        ],
    )
    write_csv(
        nextalign_dir / "reference_aln" / "REF_B" / "part.errors.csv",
        [["SEQ_3", "too many Ns", "warning text"]],
    )

    filtered = collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "out.tsv"))

    assert sorted(filtered.keys()) == ["SEQ_1"]
    assert filtered["SEQ_1"]["reference"] == "REF_A"
    assert "SEQ_3" not in filtered


def test_write_outputs(tmp_path: Path):
    filtered = {
        "SEQ_B": {"reference": "R2", "error": "error B", "warnings": ""},
        "SEQ_A": {"reference": "R1", "error": "error A", "warnings": "warn"},
    }

    out_tsv = tmp_path / "filtered_sequences.tsv"
    write_filtered_list(filtered, str(out_tsv))
    ids_file = write_filtered_ids_only(filtered, str(out_tsv))

    lines = out_tsv.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "seq_name\treference\terror\twarnings"
    # sorted by seq_name
    assert lines[1].startswith("SEQ_A\t")
    assert lines[2].startswith("SEQ_B\t")

    assert Path(ids_file).read_text(encoding="utf-8").splitlines() == ["SEQ_A", "SEQ_B"]


def test_cli_creates_empty_outputs_when_no_errors(tmp_path: Path):
    nextalign_dir = tmp_path / "Nextalign"
    write_csv(
        nextalign_dir / "query_aln" / "REF_X" / "empty.errors.csv",
        [["SEQ_OK", "", "only warning"]],
    )

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-n",
            str(nextalign_dir),
            "-o",
            "filtered.tsv",
            "-b",
            str(tmp_path),
        ],
        check=True,
    )

    out_tsv = tmp_path / "filtered.tsv"
    out_ids = tmp_path / "filtered_ids.txt"

    assert out_tsv.exists()
    assert out_ids.exists()
    assert out_tsv.read_text(encoding="utf-8").strip() == "seq_name\treference\terror\twarnings"
    assert out_ids.read_text(encoding="utf-8") == ""


def test_collect_filtered_sequences_with_fixture_dataset(tmp_path: Path):
    nextalign_dir = DATA_DIR / "Nextalign"
    out_tsv = tmp_path / "fixture_filtered.tsv"

    filtered = collect_filtered_sequences(str(nextalign_dir), str(out_tsv))
    write_filtered_list(filtered, str(out_tsv))
    out_ids = write_filtered_ids_only(filtered, str(out_tsv))

    assert out_tsv.read_text(encoding="utf-8") == (DATA_DIR / "expected_filtered.tsv").read_text(encoding="utf-8")
    assert Path(out_ids).read_text(encoding="utf-8") == (DATA_DIR / "expected_ids.txt").read_text(encoding="utf-8")


def test_collect_filtered_sequences_marks_unprojectable_queries(tmp_path: Path):
    nextalign_dir = tmp_path / "Nextalign"

    # reference_aln contains only REF_OK as projectable into master coordinates
    ref_dir = nextalign_dir / "reference_aln" / "MASTER1"
    ref_dir.mkdir(parents=True, exist_ok=True)
    (ref_dir / "MASTER1.aligned.fasta").write_text(
        ">MASTER1\nATGC\n>REF_OK\nATGC\n",
        encoding="utf-8",
    )

    # query_aln for REF_OK should remain unflagged
    q_ok_dir = nextalign_dir / "query_aln" / "REF_OK"
    q_ok_dir.mkdir(parents=True, exist_ok=True)
    (q_ok_dir / "REF_OK.aligned.fasta").write_text(
        ">REF_OK\nATGC\n>Q_OK\nATGC\n",
        encoding="utf-8",
    )

    # query_aln for REF_ORPHAN cannot be projected because REF_ORPHAN absent in reference_aln
    q_orphan_dir = nextalign_dir / "query_aln" / "REF_ORPHAN"
    q_orphan_dir.mkdir(parents=True, exist_ok=True)
    (q_orphan_dir / "REF_ORPHAN.aligned.fasta").write_text(
        ">REF_ORPHAN\nATGC\n>Q_ORPHAN\nATGC\n",
        encoding="utf-8",
    )

    filtered = collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "filtered.tsv"))

    assert "Q_ORPHAN" in filtered
    assert filtered["Q_ORPHAN"]["reference"] == "REF_ORPHAN"
    assert "cannot be projected" in filtered["Q_ORPHAN"]["error"]
    assert "Q_OK" not in filtered


def test_collect_filtered_sequences_raises_when_nextalign_dir_missing(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="Nextalign directory not found"):
        collect_filtered_sequences(str(tmp_path / "missing"), str(tmp_path / "out.tsv"))


def test_collect_filtered_sequences_raises_on_malformed_errors_csv(tmp_path: Path):
    nextalign_dir = tmp_path / "Nextalign"
    bad_csv = nextalign_dir / "query_aln" / "REF_A" / "bad.errors.csv"
    bad_csv.parent.mkdir(parents=True, exist_ok=True)
    bad_csv.write_text("seq\terr\nQ1\toops\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Malformed nextalign errors file"):
        collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "out.tsv"))


# ---------------------------------------------------------------------------
# Regression tests for the HA/NA subtype-loss defect (MISSING_H3_report.md).
#
# collect_unprojectable_queries used to decide projectability from
# Nextalign/reference_aln/, i.e. "could nextalign align this reference to a
# master?". For segmented influenza that is the wrong question: segment 4's
# master AB573800 is H1N1, so 87 of the 117 HA references failed to align to it
# and every query in those groups - 366,623 sequences, 99.98% of all H3 - was
# marked unprojectable and deleted by PadAlignment --skip_ids. The right question
# is "is this reference a row of the backbone PadAlignment will open?", and
# refset_4_aln.fasta contains all 117.
# ---------------------------------------------------------------------------


def _divergent_subtype_fixture(tmp_path: Path):
    """Segment 4 in miniature: a master plus two references too divergent for
    nextalign to align against it, all three present in the segment backbone."""
    nextalign_dir = tmp_path / "Nextalign"

    # nextalign managed only the master itself against the master
    ref_aln = nextalign_dir / "reference_aln" / "MASTER_H1"
    ref_aln.mkdir(parents=True, exist_ok=True)
    (ref_aln / "MASTER_H1.aligned.fasta").write_text(">MASTER_H1\nACGTACGT\n", encoding="utf-8")

    # queries exist for the master and for both divergent references
    for ref, queries in (("MASTER_H1", ["Q_H1"]), ("REF_H3", ["Q_H3_1", "Q_H3_2"]), ("REF_H7", ["Q_H7"])):
        d = nextalign_dir / "query_aln" / ref
        d.mkdir(parents=True, exist_ok=True)
        body = f">{ref}\nACGTACGT\n" + "".join(f">{q}\nACGTACGT\n" for q in queries)
        (d / f"{ref}.aligned.fasta").write_text(body, encoding="utf-8")

    # the per-segment backbone holds all three references
    ref_dir = tmp_path / "ref_set_aligned"
    ref_dir.mkdir(parents=True, exist_ok=True)
    (ref_dir / "refset_4_aln.fasta").write_text(
        ">MASTER_H1\nACGTACGT\n>REF_H3\nACGTACGT\n>REF_H7\nACGTACGT\n", encoding="utf-8")

    ref_list = tmp_path / "ref_list_refmast.txt"
    ref_list.write_text(
        "MASTER_H1\tmaster\t4\nREF_H3\treference\t4\nREF_H7\treference\t4\n", encoding="utf-8")

    return nextalign_dir, ref_dir, ref_list


def test_precomputed_ref_dir_overrides_reference_aln(tmp_path: Path):
    """The fix: divergent-subtype queries survive because the backbone holds their
    reference, even though reference_aln does not. Fails before the fix."""
    nextalign_dir, ref_dir, ref_list = _divergent_subtype_fixture(tmp_path)

    filtered = collect_filtered_sequences(
        str(nextalign_dir), str(tmp_path / "out.tsv"),
        precomputed_ref_dir=str(ref_dir), ref_list=str(ref_list),
    )

    assert filtered == {}, f"nothing should be filtered, got {sorted(filtered)}"


def test_reference_aln_alone_would_have_deleted_them(tmp_path: Path):
    """Documents the old behaviour, which is still correct for non-segmented
    builds: without a backbone directory, reference_aln is the only authority."""
    nextalign_dir, _ref_dir, _ref_list = _divergent_subtype_fixture(tmp_path)

    filtered = collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "out.tsv"))

    assert set(filtered) == {"Q_H3_1", "Q_H3_2", "Q_H7"}
    assert all("cannot be projected" in v["error"] for v in filtered.values())
    assert "Q_H1" not in filtered


def test_ref_list_alone_is_enough_when_backbones_are_absent(tmp_path: Path):
    """--ref_list without --precomputed_ref_dir still routes through the shared
    predicate, which falls back to reference_aln per master."""
    nextalign_dir, _ref_dir, ref_list = _divergent_subtype_fixture(tmp_path)

    filtered = collect_filtered_sequences(
        str(nextalign_dir), str(tmp_path / "out.tsv"), ref_list=str(ref_list))

    assert set(filtered) == {"Q_H3_1", "Q_H3_2", "Q_H7"}


def test_empty_projectable_set_raises_rather_than_deleting_everything(tmp_path: Path):
    """An empty projectable set marks every query unprojectable. That is the
    shape of the original incident, so it must be a hard error."""
    nextalign_dir, _ref_dir, _ref_list = _divergent_subtype_fixture(tmp_path)
    empty_dir = tmp_path / "empty_backbones"
    empty_dir.mkdir()
    (empty_dir / "refset_4_aln.fasta").write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="zero projectable reference accessions"):
        collect_filtered_sequences(
            str(nextalign_dir), str(tmp_path / "out.tsv"),
            precomputed_ref_dir=str(empty_dir), ref_list=str(_ref_list))


def test_count_query_sequences_excludes_reference_rows(tmp_path: Path):
    """The G1 denominator counts queries, not the reference row PadAlignment
    re-inserts from the backbone."""
    from CollectFilteredSequences import count_query_sequences

    nextalign_dir, _, _ = _divergent_subtype_fixture(tmp_path)
    assert count_query_sequences(str(nextalign_dir)) == 4  # Q_H1, Q_H3_1, Q_H3_2, Q_H7


def test_guard_g1_aborts_when_too_much_is_filtered(tmp_path: Path):
    """CLI-level: filtering 3 of 4 queries must exit non-zero, not warn."""
    nextalign_dir, _ref_dir, _ref_list = _divergent_subtype_fixture(tmp_path)

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH),
         "-n", str(nextalign_dir), "-o", "filtered.tsv", "-b", str(tmp_path),
         "--max_filtered_fraction", "0.02"],
        capture_output=True, text=True,
    )
    assert result.returncode == 3, result.stdout + result.stderr
    assert "exceeds --max_filtered_fraction" in result.stderr


def test_guard_g1_passes_when_the_backbone_is_supplied(tmp_path: Path):
    """Same inputs, same threshold, but asking the right question - exit 0."""
    nextalign_dir, ref_dir, ref_list = _divergent_subtype_fixture(tmp_path)

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH),
         "-n", str(nextalign_dir), "-o", "filtered.tsv", "-b", str(tmp_path),
         "--precomputed_ref_dir", str(ref_dir), "--ref_list", str(ref_list),
         "--max_filtered_fraction", "0.02"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
