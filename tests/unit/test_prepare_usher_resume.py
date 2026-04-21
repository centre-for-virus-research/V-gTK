from pathlib import Path

from PrepareUsherResume import prepare_resume


def write_fasta(path: Path, records):
	with path.open("w", encoding="utf-8") as handle:
		for seq_id, sequence in records:
			handle.write(f">{seq_id}\n{sequence}\n")


def test_prepare_resume_uses_latest_completed_chunk_and_merges_remaining_ids(tmp_path: Path):
	full_fasta = tmp_path / "full.fasta"
	write_fasta(
		full_fasta,
		[
			("REF1", "AAAA"),
			("PLACED1", "AAAT"),
			("PLACED2", "AACC"),
			("REMAIN1", "CCCC"),
			("REMAIN2", "GGGG"),
			("REMAIN3", "TTTT"),
		],
	)

	run_dir = tmp_path / "out_seg4_ushertest2"
	run_dir.mkdir()
	chunk_0001 = run_dir / "chunk_0001"
	chunk_0001.mkdir()
	(chunk_0001 / "uncondensed-final-tree.nh").write_text("(REF1:0.1,PLACED1:0.2);\n", encoding="utf-8")
	chunk_0002 = run_dir / "chunk_0002"
	chunk_0002.mkdir()
	(chunk_0002 / "uncondensed-final-tree.nh").write_text("(REF1:0.1,PLACED1:0.2,PLACED2:0.3);\n", encoding="utf-8")
	chunk_0003 = run_dir / "chunk_0003"
	chunk_0003.mkdir()

	chunks_dir = run_dir / "chunks"
	chunks_dir.mkdir()
	write_fasta(chunks_dir / "chunk_0001.fasta", [("REF1", "AAAA"), ("PLACED1", "AAAT")])
	write_fasta(chunks_dir / "chunk_0002.fasta", [("REF1", "AAAA"), ("PLACED2", "AACC")])
	write_fasta(chunks_dir / "chunk_0003.fasta", [("REF1", "AAAA"), ("REMAIN1", "CCCC"), ("REMAIN2", "GGGG")])
	write_fasta(chunks_dir / "chunk_0004.fasta", [("REF1", "AAAA"), ("REMAIN3", "TTTT")])

	resume_dir = tmp_path / "resume"
	outputs = prepare_resume(str(run_dir), str(full_fasta), str(resume_dir))

	assert Path(outputs["resume_tree"]).read_text(encoding="utf-8").strip().endswith(";")
	assert Path(outputs["existing_ids"]).read_text(encoding="utf-8").strip().splitlines() == ["PLACED1", "PLACED2", "REF1"]
	resume_ids = [line[1:] for line in Path(outputs["resume_alignment"]).read_text(encoding="utf-8").splitlines() if line.startswith(">")]
	assert resume_ids == ["REF1", "PLACED1", "PLACED2", "REMAIN1", "REMAIN2", "REMAIN3"]
	manifest_text = Path(outputs["manifest"]).read_text(encoding="utf-8")
	assert "latest_completed_chunk\tchunk_0002" in manifest_text
	assert "remaining_chunk_count\t2" in manifest_text
