import os
import sqlite3
from pathlib import Path

import pytest

from UsherPlacement import UsherPlacement


def test_prepare_update_assets_exports_tree_and_existing_ids(tmp_path: Path, basic_update_db: Path):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(">REF1\nATGC\n>Q_NEW\nATGT\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		update_db=str(basic_update_db),
	)

	tree_path, ids_path = processor.prepare_update_assets()

	assert Path(tree_path).exists()
	assert Path(ids_path).exists()
	assert "REF1" in Path(ids_path).read_text(encoding="utf-8")
	assert "Q_EXCL" not in Path(ids_path).read_text(encoding="utf-8")
	assert "Q_OLD" in Path(tree_path).read_text(encoding="utf-8")


def test_update_mode_uses_tree_terminals_as_existing_ids_and_rehydrates_missing_db_alignments(tmp_path: Path):
	db_path = tmp_path / "update.db"
	conn = sqlite3.connect(str(db_path))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
		cur.executemany(
			"INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
			[("REF1", "master", "1"), ("Q_OLD", "query", "1"), ("Q_NEW", "query", "1")],
		)
		cur.execute(
			"CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, alignment TEXT, segment TEXT)"
		)
		cur.executemany(
			"INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, alignment, segment) VALUES (?, ?, ?, ?, ?)",
			[("REF1", "REF1", "REF1", "ATGC", "1"), ("Q_OLD", "Q_OLD", "REF1", "ATGT", "1")],
		)
		cur.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, segment_key TEXT, newick TEXT)")
		cur.execute(
			"INSERT INTO trees(name, source, segment, segment_key, newick) VALUES (?, ?, ?, ?, ?)",
			("usher_seg1", "usher", "1", "seg_1", "(REF1:0.1);"),
		)
		conn.commit()
	finally:
		conn.close()

	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(">REF1\nATGC\n>Q_NEW\nATGA\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		update_db=str(db_path),
	)

	tree_path, ids_path = processor.prepare_update_assets()
	assert Path(ids_path).read_text(encoding="utf-8").strip().splitlines() == ["REF1"]

	merged_path = processor.build_update_alignment_input(tree_path)
	merged_text = Path(merged_path).read_text(encoding="utf-8")
	assert ">Q_NEW" in merged_text
	assert ">Q_OLD" in merged_text


def test_split_alignment_into_chunks_writes_ref_into_each_chunk(tmp_path: Path):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(
		">REF1\nATGC\n>Q1\nATGA\n>Q2\nATGT\n>Q3\nATGG\n>Q4\nATCC\n>Q5\nATCA\n",
		encoding="utf-8",
	)

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		chunk_size=2,
		chunk_threshold=3,
	)

	chunk_paths = processor.split_alignment_into_chunks(str(msa), "REF1", ["Q1", "Q2", "Q3", "Q4", "Q5"])

	assert len(chunk_paths) == 3
	for chunk_path in chunk_paths:
		text = Path(chunk_path).read_text(encoding="utf-8")
		assert text.startswith(">REF1\nATGC\n")

	assert Path(chunk_paths[0]).read_text(encoding="utf-8").count(">") == 3
	assert Path(chunk_paths[1]).read_text(encoding="utf-8").count(">") == 3
	assert Path(chunk_paths[2]).read_text(encoding="utf-8").count(">") == 2


def test_split_alignment_into_chunks_falls_back_to_python_when_seqkit_missing(tmp_path: Path, monkeypatch):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(
		">REF1\nATGC\n>Q1\nATGA\n>Q2\nATGT\n>Q3\nATGG\n",
		encoding="utf-8",
	)

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		chunk_size=2,
	)

	monkeypatch.setattr("shutil.which", lambda name: None)

	chunk_paths = processor.split_alignment_into_chunks(str(msa), "REF1", ["Q1", "Q2", "Q3"])

	assert len(chunk_paths) == 2
	assert Path(chunk_paths[0]).read_text(encoding="utf-8").startswith(">REF1\nATGC\n")


def test_split_alignment_into_chunks_uses_seqkit_outputs_when_available(tmp_path: Path, monkeypatch):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(
		">REF1\nATGC\n>Q1\nATGA\n>Q2\nATGT\n>Q3\nATGG\n",
		encoding="utf-8",
	)

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		chunk_size=2,
	)

	monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/seqkit")

	run_calls = []
	def fake_run(cmd, check=True):
		run_calls.append(cmd)
		chunk_dir = Path(processor.output_dir) / "chunks"
		raw_dir = chunk_dir / "raw"
		raw_dir.mkdir(parents=True, exist_ok=True)
		if cmd[1] == "grep":
			if str(cmd[-1]).endswith("placeable_only.fasta"):
				Path(cmd[-1]).write_text(">Q1\nATGA\n>Q2\nATGT\n>Q3\nATGG\n", encoding="utf-8")
		elif cmd[1] == "split2":
			(raw_dir / "placeable_only.part_001.fasta").write_text(">Q1\nATGA\n>Q2\nATGT\n", encoding="utf-8")
			(raw_dir / "placeable_only.part_002.fasta").write_text(">Q3\nATGG\n", encoding="utf-8")
		return None

	monkeypatch.setattr("subprocess.run", fake_run)

	chunk_paths = processor.split_alignment_into_chunks(str(msa), "REF1", ["Q1", "Q2", "Q3"])

	assert len(chunk_paths) == 2
	assert any(cmd[1] == "split2" for cmd in run_calls)
	for chunk_path in chunk_paths:
		assert Path(chunk_path).read_text(encoding="utf-8").startswith(">REF1\nATGC\n")


def test_run_update_mode_chunks_iteratively_for_large_alignment(tmp_path: Path, monkeypatch):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(
		">REF1\nATGC\n>Q1\nATGA\n>Q2\nATGT\n>Q3\nATGG\n>Q4\nATCC\n>Q5\nATCA\n",
		encoding="utf-8",
	)
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	seed_tree = output_dir / "seed_tree.nwk"
	seed_tree.write_text("(REF1:0.1);\n", encoding="utf-8")
	ids_file = output_dir / "existing_ids_segment_1.txt"
	ids_file.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
		update_db=str(tmp_path / "dummy.db"),
		chunk_size=2,
		chunk_threshold=3,
	)

	monkeypatch.setattr(processor, "prepare_update_assets", lambda: (str(seed_tree), str(ids_file)))
	monkeypatch.setattr(processor, "build_update_alignment_input", lambda tree_file: str(msa))

	vcf_calls = []
	def fake_build_vcf(ref_id, exclude_ids_file=None, alignment_fasta=None):
		assert alignment_fasta is not None
		vcf_calls.append(Path(alignment_fasta).name)
		vcf_path = output_dir / f"{Path(alignment_fasta).stem}.vcf"
		vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
		return str(vcf_path)
	monkeypatch.setattr(processor, "build_vcf", fake_build_vcf)

	run_calls = []
	def fake_run_usher(tree_file, vcf_path, chunk_output_dir):
		run_calls.append((os.path.basename(tree_file), os.path.basename(chunk_output_dir)))
		Path(chunk_output_dir).mkdir(parents=True, exist_ok=True)
		chunk_idx = len(run_calls)
		(Path(chunk_output_dir) / "uncondensed-final-tree.nh").write_text(
			f"(REF1:0.1,Q{chunk_idx}:0.2);\n",
			encoding="utf-8",
		)
		(Path(chunk_output_dir) / "usher.pb").write_text("pb", encoding="utf-8")
		return str(Path(chunk_output_dir) / "uncondensed-final-tree.nh")
	monkeypatch.setattr(processor, "run_usher", fake_run_usher)

	processor.run()

	assert len(vcf_calls) == 3
	assert run_calls[0][0] == "seed_tree.nwk"
	assert run_calls[1][0] == "uncondensed-final-tree.nh"
	assert run_calls[2][0] == "uncondensed-final-tree.nh"
	assert (output_dir / "uncondensed-final-tree.nh").exists()