import csv
from pathlib import Path

import pytest
from Bio import SeqIO

from MMseqsClustering import (
    apply_trimming,
    parse_ranges,
    run_mmseqs_clustering,
    strip_alignment_gaps,
    trim_alignment,
    write_aligned_representatives,
)


def write_fasta(path: Path, records):
    with path.open("w", encoding="utf-8") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


def test_parse_ranges_with_open_bounds():
    # alignment length is used for open-ended ranges
    assert parse_ranges(["1:3", "5:"], alignment_length=8) == [(0, 3), (4, 8)]
    assert parse_ranges([":2"], alignment_length=8) == [(0, 2)]


def test_trim_alignment_and_apply_trimming(tmp_path: Path):
    input_dir = tmp_path / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    in_fasta = input_dir / "segA.fasta"
    write_fasta(in_fasta, [("S1", "ABCDEFGH"), ("S2", "12345678")])

    out_fasta = tmp_path / "trimmed.fasta"
    trim_alignment(str(in_fasta), str(out_fasta), ["1:2", "5:6"])

    recs = list(SeqIO.parse(str(out_fasta), "fasta"))
    assert str(recs[0].seq) == "ABEF"
    assert str(recs[1].seq) == "1256"

    cds = tmp_path / "cds.tsv"
    with cds.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["input_fasta", "ranges"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"input_fasta": "segA.fasta", "ranges": "2:3,7:8"})

    trimmed_map = apply_trimming(str(cds), str(input_dir), str(tmp_path / "trimmed_dir"))
    assert "segA.fasta" in trimmed_map

    trimmed_records = list(SeqIO.parse(trimmed_map["segA.fasta"], "fasta"))
    assert str(trimmed_records[0].seq) == "BCGH"


def test_run_mmseqs_clustering_builds_expected_command_chain(tmp_path: Path, monkeypatch):
    """--keep-gaps preserves the original command chain exactly."""
    in_fasta = tmp_path / "input.fasta"
    write_fasta(in_fasta, [("S1", "AAAA"), ("S2", "AAAT")])

    calls = []

    def fake_run(cmd, check=True):
        calls.append(cmd)
        return 0

    monkeypatch.setattr("MMseqsClustering.subprocess.run", fake_run)

    run_mmseqs_clustering(str(in_fasta), str(tmp_path / "out"), min_seq_id=0.9, threads=3, strip_gaps=False)

    assert len(calls) == 7
    assert calls[0][:2] == ["mmseqs", "createdb"]
    assert calls[1][:2] == ["mmseqs", "cluster"]
    assert "--min-seq-id" in calls[1]
    assert "--threads" in calls[1]
    assert calls[2][:2] == ["mmseqs", "createtsv"]
    assert calls[3][:2] == ["mmseqs", "createseqfiledb"]
    assert calls[4][:2] == ["mmseqs", "result2flat"]
    assert calls[5][:2] == ["mmseqs", "createsubdb"]
    assert calls[6][:2] == ["mmseqs", "convert2fasta"]


def test_strip_alignment_gaps_removes_gaps_and_keeps_ids(tmp_path: Path):
    aligned = tmp_path / "aligned.fasta"
    write_fasta(aligned, [("S1", "AC--GT"), ("S2", "A.CGT-"), ("S3", "------")])
    out = tmp_path / "ungapped.seq"

    written = strip_alignment_gaps(str(aligned), str(out))

    assert written == 2  # the all-gap record is dropped
    records = {r.id: str(r.seq) for r in SeqIO.parse(str(out), "fasta")}
    assert records == {"S1": "ACGT", "S2": "ACGT"}


def test_write_aligned_representatives_uses_the_padded_alignment(tmp_path: Path):
    """Representatives must come back aligned and equal length: IQ_TREE aborts on
    a ragged alignment, and VERY_FAST_TREE has no such guard at all."""
    aligned = tmp_path / "aligned.fasta"
    write_fasta(aligned, [("S1", "AC--GT"), ("S2", "ACGGGT"), ("S3", "AC-CGT")])
    tsv = tmp_path / "clusters.tsv"
    tsv.write_text("S1\tS1\nS1\tS3\nS2\tS2\n", encoding="utf-8")
    out = tmp_path / "reps.fasta"

    count = write_aligned_representatives(str(aligned), str(tsv), str(out))

    assert count == 2
    records = {r.id: str(r.seq) for r in SeqIO.parse(str(out), "fasta")}
    assert records == {"S1": "AC--GT", "S2": "ACGGGT"}
    assert len({len(s) for s in records.values()}) == 1


def test_write_aligned_representatives_rejects_missing_representative(tmp_path: Path):
    aligned = tmp_path / "aligned.fasta"
    write_fasta(aligned, [("S1", "ACGT")])
    tsv = tmp_path / "clusters.tsv"
    tsv.write_text("S1\tS1\nGHOST\tGHOST\n", encoding="utf-8")

    with pytest.raises(ValueError, match="GHOST"):
        write_aligned_representatives(str(aligned), str(tsv), str(tmp_path / "reps.fasta"))


def test_run_mmseqs_clustering_strips_gaps_and_rebuilds_aligned_reps(tmp_path: Path, monkeypatch):
    in_fasta = tmp_path / "input.fasta"
    write_fasta(in_fasta, [("S1", "AC--GT"), ("S2", "ACGGGT")])
    out_dir = tmp_path / "out"

    calls = []

    def fake_run(cmd, check=True):
        calls.append(cmd)
        # Stand in for `mmseqs createtsv`, which write_aligned_representatives reads.
        if cmd[:2] == ["mmseqs", "createtsv"]:
            Path(cmd[5]).write_text("S1\tS1\nS1\tS2\n", encoding="utf-8")
        return 0

    monkeypatch.setattr("MMseqsClustering.subprocess.run", fake_run)

    run_mmseqs_clustering(str(in_fasta), str(out_dir), min_seq_id=0.9, threads=3)

    # createseqfiledb/result2flat/createsubdb/convert2fasta are all skipped.
    assert [c[1] for c in calls] == ["createdb", "cluster", "createtsv"]

    base = out_dir / "input"
    # The scratch unaligned copy must not survive: four call sites fall back to a
    # generic `*.fasta` glob over this directory.
    assert not list(base.glob("*_ungapped*"))
    assert not list(base.rglob("*_ungapped*"))
    assert [p.name for p in base.glob("*.fasta")] == ["input_cluster_rep.fasta"]

    reps = {r.id: str(r.seq) for r in SeqIO.parse(str(base / "input_cluster_rep.fasta"), "fasta")}
    assert reps == {"S1": "AC--GT"}


def test_run_mmseqs_clustering_passes_max_seqs_through(tmp_path: Path, monkeypatch):
    in_fasta = tmp_path / "input.fasta"
    write_fasta(in_fasta, [("S1", "ACGT")])

    calls = []

    def fake_run(cmd, check=True):
        calls.append(cmd)
        if cmd[:2] == ["mmseqs", "createtsv"]:
            Path(cmd[5]).write_text("S1\tS1\n", encoding="utf-8")
        return 0

    monkeypatch.setattr("MMseqsClustering.subprocess.run", fake_run)

    run_mmseqs_clustering(str(in_fasta), str(tmp_path / "out"), min_seq_id=0.9, threads=3, max_seqs=2000)

    cluster_cmd = next(c for c in calls if c[:2] == ["mmseqs", "cluster"])
    assert "--max-seqs" in cluster_cmd
    assert cluster_cmd[cluster_cmd.index("--max-seqs") + 1] == "2000"

    default_calls = []
    monkeypatch.setattr("MMseqsClustering.subprocess.run", lambda cmd, check=True: (
        default_calls.append(cmd),
        Path(cmd[5]).write_text("S1\tS1\n", encoding="utf-8") if cmd[:2] == ["mmseqs", "createtsv"] else None,
        0,
    )[-1])
    run_mmseqs_clustering(str(in_fasta), str(tmp_path / "out2"), min_seq_id=0.9, threads=3)
    cluster_cmd = next(c for c in default_calls if c[:2] == ["mmseqs", "cluster"])
    assert "--max-seqs" not in cluster_cmd
