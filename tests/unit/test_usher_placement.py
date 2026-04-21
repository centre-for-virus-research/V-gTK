import os
import sqlite3
import subprocess
from pathlib import Path

import pytest

from UsherPlacement import UsherPlacement


def write_fasta(path: Path, records):
	with path.open("w", encoding="utf-8") as handle:
		for seq_id, sequence in records:
			handle.write(f">{seq_id}\n{sequence}\n")


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


def test_collapse_identical_sequences_prefers_existing_tree_anchor(tmp_path: Path):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(
		">REF1\nATGC\n>C1\nATGA\n>Q1\nATGA\n>Q2\nATTT\n>Q3\nATTT\n",
		encoding="utf-8",
	)

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
	)

	plan = processor.collapse_identical_sequences(str(msa), "REF1", {"C1"})

	assert plan["placeable_ids"] == ["Q2"]
	assert plan["anchor_to_members"] == {"C1": ["Q1"], "Q2": ["Q3"]}
	assert plan["member_to_anchor"] == {"Q1": "C1", "Q3": "Q2"}


def test_expand_identical_sequence_tree_outputs_adds_zero_length_polytomy(tmp_path: Path):
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	for name in ["uncondensed-final-tree.nh", "final-tree.nh"]:
		(output_dir / name).write_text("(REF1:0.1,Q2:0.2);\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(tmp_path / "segment_1_dedup.fasta"),
		output_dir=str(output_dir),
	)

	processor.expand_identical_sequence_tree_outputs({"Q2": ["Q3", "Q4"]})

	for name in ["uncondensed-final-tree.nh", "final-tree.nh"]:
		text = (output_dir / name).read_text(encoding="utf-8")
		assert "Q2:0.00000" in text
		assert "Q3:0.00000" in text
		assert "Q4:0.00000" in text


def test_run_skips_usher_when_only_identical_duplicates_need_attaching(tmp_path: Path, monkeypatch):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(
		">REF1\nATGC\n>C1\nATGA\n>Q1\nATGA\n",
		encoding="utf-8",
	)
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	cluster_rep = tmp_path / "cluster_rep.fasta"
	cluster_rep.write_text(">REF1\nATGC\n>C1\nATGA\n", encoding="utf-8")
	seed_tree = tmp_path / "seed.treefile"
	seed_tree.write_text("(REF1:0.1,C1:0.2);\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
		mmseq_cluster_dir=str(tmp_path / "mmseq"),
		iqtree_dir=str(tmp_path / "iqtree"),
	)

	monkeypatch.setattr(processor, "resolve_non_update_assets", lambda: (str(cluster_rep), str(seed_tree)))

	def fail_run_usher(tree_file, vcf_path, chunk_output_dir):
		raise AssertionError("run_usher should not be called when only duplicate attachment is required")

	monkeypatch.setattr(processor, "run_usher", fail_run_usher)

	processor.run()

	tree_text = (output_dir / "uncondensed-final-tree.nh").read_text(encoding="utf-8")
	assert "Q1:0.00000" in tree_text
	report_lines = (output_dir / "identical_sequence_groups.tsv").read_text(encoding="utf-8").strip().splitlines()
	assert report_lines == [
		"anchor_id\tmember_id\tanchor_requires_placement",
		"C1\tQ1\t0",
	]
	assert (output_dir / "exclude_ids.txt").read_text(encoding="utf-8").strip().splitlines() == ["C1"]


def test_run_chunked_mode_only_places_unique_anchors_from_duplicate_swarms(tmp_path: Path, monkeypatch):
	msa = tmp_path / "segment_1_dedup.fasta"
	write_fasta(
		msa,
		[
			("REF1", "AAAA"),
			("C1", "AAAT"),
			("SWARM_A1", "CCCC"),
			("SWARM_A2", "CCCC"),
			("SWARM_A3", "CCCC"),
			("SWARM_B1", "GGGG"),
			("SWARM_B2", "GGGG"),
			("SWARM_C1", "TTTA"),
			("SWARM_C2", "TTTA"),
			("SWARM_C3", "TTTA"),
			("SWARM_C4", "TTTA"),
			("UNIQ1", "ACTG"),
		],
	)
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	cluster_rep = tmp_path / "cluster_rep.fasta"
	write_fasta(cluster_rep, [("REF1", "AAAA"), ("C1", "AAAT")])
	seed_tree = tmp_path / "seed.treefile"
	seed_tree.write_text("(REF1:0.1,C1:0.2);\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
		mmseq_cluster_dir=str(tmp_path / "mmseq"),
		iqtree_dir=str(tmp_path / "iqtree"),
		chunk_size=2,
		chunk_threshold=3,
	)

	monkeypatch.setattr(processor, "resolve_non_update_assets", lambda: (str(cluster_rep), str(seed_tree)))

	vcf_inputs = {}
	def fake_build_vcf(ref_id, exclude_ids_file=None, alignment_fasta=None):
		chunk_ids = processor._read_ids_from_fasta(alignment_fasta)
		vcf_path = output_dir / f"{Path(alignment_fasta).stem}.vcf"
		vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
		vcf_inputs[str(vcf_path)] = chunk_ids
		return str(vcf_path)
	monkeypatch.setattr(processor, "build_vcf", fake_build_vcf)

	placed_ids = set()
	def fake_run_usher(tree_file, vcf_path, chunk_output_dir):
		for seq_id in vcf_inputs[vcf_path]:
			if seq_id != "REF1":
				placed_ids.add(seq_id)
			assert seq_id not in {"SWARM_A2", "SWARM_A3", "SWARM_B2", "SWARM_C2", "SWARM_C3", "SWARM_C4"}
		Path(chunk_output_dir).mkdir(parents=True, exist_ok=True)
		all_tips = ["REF1", "C1"] + sorted(placed_ids)
		(Path(chunk_output_dir) / "uncondensed-final-tree.nh").write_text(
			"(" + ",".join(f"{tip}:0.2" for tip in all_tips) + ");\n",
			encoding="utf-8",
		)
		(Path(chunk_output_dir) / "usher.pb").write_text("pb", encoding="utf-8")
		return str(Path(chunk_output_dir) / "uncondensed-final-tree.nh")
	monkeypatch.setattr(processor, "run_usher", fake_run_usher)

	processor.run()

	assert placed_ids == {"SWARM_A1", "SWARM_B1", "SWARM_C1", "UNIQ1"}
	assert len(vcf_inputs) == 2
	all_chunk_ids = [seq_id for chunk_ids in vcf_inputs.values() for seq_id in chunk_ids if seq_id != "REF1"]
	assert sorted(all_chunk_ids) == ["SWARM_A1", "SWARM_B1", "SWARM_C1", "UNIQ1"]

	tree_text = (output_dir / "uncondensed-final-tree.nh").read_text(encoding="utf-8")
	for seq_id in ["SWARM_A2", "SWARM_A3", "SWARM_B2", "SWARM_C2", "SWARM_C3", "SWARM_C4"]:
		assert f"{seq_id}:0.00000" in tree_text

	report_rows = (output_dir / "identical_sequence_groups.tsv").read_text(encoding="utf-8").strip().splitlines()
	assert report_rows == [
		"anchor_id\tmember_id\tanchor_requires_placement",
		"SWARM_A1\tSWARM_A2\t1",
		"SWARM_A1\tSWARM_A3\t1",
		"SWARM_B1\tSWARM_B2\t1",
		"SWARM_C1\tSWARM_C2\t1",
		"SWARM_C1\tSWARM_C3\t1",
		"SWARM_C1\tSWARM_C4\t1",
	]


def test_run_update_mode_collapses_against_existing_tree_and_new_duplicate_anchor(tmp_path: Path, monkeypatch):
	msa = tmp_path / "segment_1_dedup.fasta"
	write_fasta(
		msa,
		[
			("REF1", "AAAA"),
			("OLD1", "CCCC"),
			("OLD_DUP1", "CCCC"),
			("OLD_DUP2", "CCCC"),
			("NEW1", "GGGG"),
			("NEW2", "GGGG"),
			("NEW3", "TTTT"),
		],
	)
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	seed_tree = output_dir / "seed_tree.nwk"
	seed_tree.write_text("(REF1:0.1,OLD1:0.2);\n", encoding="utf-8")
	ids_file = output_dir / "existing_ids_segment_1.txt"
	ids_file.write_text("REF1\nOLD1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
		update_db=str(tmp_path / "dummy.db"),
	)

	monkeypatch.setattr(processor, "prepare_update_assets", lambda: (str(seed_tree), str(ids_file)))
	monkeypatch.setattr(processor, "build_update_alignment_input", lambda tree_file: str(msa))

	vcf_calls = []
	def fake_build_vcf(ref_id, exclude_ids_file=None, alignment_fasta=None):
		vcf_calls.append({
			"exclude_ids": Path(exclude_ids_file).read_text(encoding="utf-8").strip().splitlines(),
			"alignment_ids": processor._read_ids_from_fasta(alignment_fasta),
		})
		vcf_path = output_dir / "all_samples.vcf"
		vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
		return str(vcf_path)
	monkeypatch.setattr(processor, "build_vcf", fake_build_vcf)

	def fake_run_usher(tree_file, vcf_path, chunk_output_dir):
		Path(chunk_output_dir).mkdir(parents=True, exist_ok=True)
		(Path(chunk_output_dir) / "uncondensed-final-tree.nh").write_text(
			"(REF1:0.1,OLD1:0.2,NEW1:0.2,NEW3:0.2);\n",
			encoding="utf-8",
		)
		(Path(chunk_output_dir) / "final-tree.nh").write_text(
			"(REF1:0.1,OLD1:0.2,NEW1:0.2,NEW3:0.2);\n",
			encoding="utf-8",
		)
		return str(Path(chunk_output_dir) / "uncondensed-final-tree.nh")
	monkeypatch.setattr(processor, "run_usher", fake_run_usher)

	processor.run()

	assert len(vcf_calls) == 1
	assert vcf_calls[0]["alignment_ids"] == ["REF1", "OLD1", "OLD_DUP1", "OLD_DUP2", "NEW1", "NEW2", "NEW3"]
	assert sorted(vcf_calls[0]["exclude_ids"]) == sorted(["OLD1", "OLD_DUP1", "OLD_DUP2", "NEW2"])

	tree_text = (output_dir / "uncondensed-final-tree.nh").read_text(encoding="utf-8")
	assert "OLD_DUP1:0.00000" in tree_text
	assert "OLD_DUP2:0.00000" in tree_text
	assert "NEW2:0.00000" in tree_text

	report_rows = (output_dir / "identical_sequence_groups.tsv").read_text(encoding="utf-8").strip().splitlines()
	assert report_rows == [
		"anchor_id\tmember_id\tanchor_requires_placement",
		"OLD1\tOLD_DUP1\t0",
		"OLD1\tOLD_DUP2\t0",
		"NEW1\tNEW2\t1",
	]


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


def test_run_chunked_mode_reports_sequence_progress(tmp_path: Path, monkeypatch, capsys):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(
		">REF1\nATGC\n>Q1\nATGA\n>Q2\nATGT\n>Q3\nATGG\n>Q4\nATCC\n>Q5\nATCA\n",
		encoding="utf-8",
	)
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	cluster_rep = tmp_path / "cluster_rep.fasta"
	cluster_rep.write_text(">REF1\nATGC\n", encoding="utf-8")
	seed_tree = tmp_path / "seed.treefile"
	seed_tree.write_text("(REF1:0.1);\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
		mmseq_cluster_dir=str(tmp_path / "mmseq"),
		iqtree_dir=str(tmp_path / "iqtree"),
		chunk_size=2,
		chunk_threshold=3,
	)

	monkeypatch.setattr(processor, "resolve_non_update_assets", lambda: (str(cluster_rep), str(seed_tree)))

	def fake_build_vcf(ref_id, exclude_ids_file=None, alignment_fasta=None):
		vcf_path = output_dir / f"{Path(alignment_fasta).stem}.vcf"
		vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
		return str(vcf_path)
	monkeypatch.setattr(processor, "build_vcf", fake_build_vcf)

	def fake_run_usher(tree_file, vcf_path, chunk_output_dir):
		Path(chunk_output_dir).mkdir(parents=True, exist_ok=True)
		(Path(chunk_output_dir) / "uncondensed-final-tree.nh").write_text("(REF1:0.1,Q1:0.2);\n", encoding="utf-8")
		(Path(chunk_output_dir) / "usher.pb").write_text("pb", encoding="utf-8")
		return str(Path(chunk_output_dir) / "uncondensed-final-tree.nh")
	monkeypatch.setattr(processor, "run_usher", fake_run_usher)

	processor.run()

	stdout = capsys.readouterr().out
	assert "batches complete: 0/3; remaining: 3." in stdout
	assert "batches complete: 1/3; remaining: 2." in stdout
	assert "batches complete: 2/3; remaining: 1." in stdout
	assert "batches complete: 3/3; remaining: 0." in stdout


def test_run_usher_writes_verbose_output_to_chunk_log(tmp_path: Path, monkeypatch):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "input.fasta"),
		output_dir=str(tmp_path / "usher_out"),
		threads=4,
	)
	chunk_dir = tmp_path / "usher_out" / "chunk_0001"
	tree_file = tmp_path / "seed.treefile"
	vcf_path = tmp_path / "all_samples.vcf"
	tree_file.write_text("(REF1:0.1);\n", encoding="utf-8")
	vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")

	run_calls = []

	class FakeCompleted:
		def __init__(self, stdout="", stderr=""):
			self.stdout = stdout
			self.stderr = stderr

	def fake_run(cmd, **kwargs):
		run_calls.append((cmd, kwargs))
		if cmd == ["usher", "--help"]:
			return FakeCompleted(stdout="usher usage -T threads", stderr="")
		log_handle = kwargs["stdout"]
		log_handle.write("verbose usher output\n")
		(chunk_dir / "uncondensed-final-tree.nh").write_text("(REF1:0.1,Q1:0.2);\n", encoding="utf-8")
		return FakeCompleted()

	monkeypatch.setattr("subprocess.run", fake_run)

	result = processor.run_usher(str(tree_file), str(vcf_path), str(chunk_dir))

	assert result == str(chunk_dir / "uncondensed-final-tree.nh")
	assert (chunk_dir / "usher.verbose.log").read_text(encoding="utf-8") == "verbose usher output\n"
	assert run_calls[1][1]["stderr"] == subprocess.STDOUT


def test_run_resume_mode_uses_supplied_tree_and_existing_ids(tmp_path: Path, monkeypatch):
	msa = tmp_path / "resume_alignment.fasta"
	msa.write_text(
		">REF1\nAAAA\n>PLACED1\nAAAT\n>Q1\nCCCC\n>Q2\nCCCC\n",
		encoding="utf-8",
	)
	output_dir = tmp_path / "resume_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	starter_tree = tmp_path / "resume_tree.nwk"
	starter_tree.write_text("(REF1:0.1,PLACED1:0.2);\n", encoding="utf-8")
	existing_ids = tmp_path / "resume_existing_ids.txt"
	existing_ids.write_text("REF1\nPLACED1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
		starter_tree=str(starter_tree),
		existing_ids_file=str(existing_ids),
	)

	vcf_calls = []
	def fake_build_vcf(ref_id, exclude_ids_file=None, alignment_fasta=None):
		vcf_calls.append({
			"alignment_ids": processor._read_ids_from_fasta(alignment_fasta),
			"exclude_ids": Path(exclude_ids_file).read_text(encoding="utf-8").strip().splitlines(),
		})
		vcf_path = output_dir / "resume.vcf"
		vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
		return str(vcf_path)
	monkeypatch.setattr(processor, "build_vcf", fake_build_vcf)

	def fake_run_usher(tree_file, vcf_path, chunk_output_dir):
		assert tree_file == str(starter_tree)
		Path(chunk_output_dir).mkdir(parents=True, exist_ok=True)
		(Path(chunk_output_dir) / "uncondensed-final-tree.nh").write_text(
			"(REF1:0.1,PLACED1:0.2,Q1:0.3);\n",
			encoding="utf-8",
		)
		(Path(chunk_output_dir) / "final-tree.nh").write_text(
			"(REF1:0.1,PLACED1:0.2,Q1:0.3);\n",
			encoding="utf-8",
		)
		return str(Path(chunk_output_dir) / "uncondensed-final-tree.nh")
	monkeypatch.setattr(processor, "run_usher", fake_run_usher)

	processor.run()

	assert len(vcf_calls) == 1
	assert vcf_calls[0]["alignment_ids"] == ["REF1", "PLACED1", "Q1", "Q2"]
	assert sorted(vcf_calls[0]["exclude_ids"]) == ["PLACED1", "Q2"]
	tree_text = (output_dir / "uncondensed-final-tree.nh").read_text(encoding="utf-8")
	assert "Q2:0.00000" in tree_text


def test_run_non_update_mode_chunks_iteratively_for_large_alignment(tmp_path: Path, monkeypatch):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(
		">REF1\nATGC\n>C1\nATGA\n>C2\nATGT\n>Q1\nATGG\n>Q2\nATCC\n>Q3\nATCA\n>Q4\nATCG\n",
		encoding="utf-8",
	)
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	cluster_rep = tmp_path / "cluster_rep.fasta"
	cluster_rep.write_text(">REF1\nATGC\n>C1\nATGA\n>C2\nATGT\n", encoding="utf-8")
	seed_tree = tmp_path / "seed.treefile"
	seed_tree.write_text("(REF1:0.1,C1:0.2,C2:0.2);\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
		mmseq_cluster_dir=str(tmp_path / "mmseq"),
		iqtree_dir=str(tmp_path / "iqtree"),
		chunk_size=2,
		chunk_threshold=3,
	)

	monkeypatch.setattr(processor, "resolve_non_update_assets", lambda: (str(cluster_rep), str(seed_tree)))

	vcf_calls = []
	def fake_build_vcf(ref_id, exclude_ids_file=None, alignment_fasta=None):
		assert alignment_fasta is not None
		assert exclude_ids_file is None
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

	assert len(vcf_calls) == 2
	assert run_calls[0][0] == "seed.treefile"
	assert run_calls[1][0] == "uncondensed-final-tree.nh"
	assert (output_dir / "uncondensed-final-tree.nh").exists()