"""Unit tests for scripts/UsherPlacement.py.

The first half covers the ordinary paths: update/resume/standard asset
resolution, identical-sequence collapsing, chunked placement and the protobuf
helpers.

The second half, from "Adversarial coverage" onward, covers the boundaries:
malformed FASTA, hostile accession identifiers, duplicate records, stale
working directories, absent external binaries, degenerate trees, and the
placement-order guarantees the module's own docstrings promise. Every test
there began as a reproduced defect or an unasserted guarantee, and its
docstring says what went wrong and why it mattered, so a regression explains
itself rather than merely going red.
"""

import os
import re
import sqlite3
import subprocess
from io import StringIO
from pathlib import Path

import pytest
from Bio import Phylo

from UsherPlacement import UsherPlacement


def write_fasta(path: Path, records):
	with path.open("w", encoding="utf-8") as handle:
		for seq_id, sequence in records:
			handle.write(f">{seq_id}\n{sequence}\n")


def make(tmp_path: Path, name="segment_1_dedup.fasta", **kwargs):
	"""A processor with its output directory already created."""
	kwargs.setdefault("padded_aln", str(tmp_path / name))
	kwargs.setdefault("output_dir", str(tmp_path / "out"))
	processor = UsherPlacement(**kwargs)
	os.makedirs(processor.output_dir, exist_ok=True)
	return processor


def ladder_newick(depth, tip_name="ANCHOR"):
	"""A maximally unbalanced (caterpillar) tree - what iterative placement of an
	increasingly divergent lineage actually produces."""
	newick = f"{tip_name}:0.1"
	for idx in range(depth):
		newick = f"(x{idx}:0.1,{newick})"
	return newick + ";"


def install_run_fakes(processor, monkeypatch, output_dir, ref_id, existing_tips, force_python_split=True):
	"""Stub out the external-tool boundary of run() (VCF construction and the two
	USHER protobuf steps) and record how they were driven. The mutation-annotated
	backbone is built once via build_backbone_pb, then every chunk is placed via
	place_onto_pb with the protobuf carried forward. Returns a records dict."""
	output_dir = Path(output_dir)
	records = {
		"vcf": [],
		"vcf_to_ids": {},
		"backbone": [],
		"place": [],
		"placed_ids": [],
	}

	if force_python_split:
		# Force the deterministic pure-python chunk splitter regardless of whether
		# seqkit happens to be installed on the machine running the tests.
		monkeypatch.setattr("shutil.which", lambda name: None)

	def fake_build_vcf(ref, exclude_ids_file=None, alignment_fasta=None, vcf_path=None):
		ids = processor._read_ids_from_fasta(alignment_fasta)
		if vcf_path is None:
			vcf_path = str(output_dir / f"{Path(alignment_fasta).stem}.vcf")
		Path(vcf_path).parent.mkdir(parents=True, exist_ok=True)
		Path(vcf_path).write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
		records["vcf"].append({
			"ids": ids,
			"vcf_path": vcf_path,
			"exclude_ids_file": exclude_ids_file,
			"alignment_fasta": alignment_fasta,
		})
		records["vcf_to_ids"][vcf_path] = ids
		return vcf_path
	monkeypatch.setattr(processor, "build_vcf", fake_build_vcf)

	def fake_build_backbone_pb(tree_file, backbone_ids, ref, alignment_fasta):
		records["backbone"].append({
			"tree_file": tree_file,
			"backbone_ids": list(backbone_ids),
			"ref_id": ref,
			"alignment_fasta": alignment_fasta,
		})
		pb = output_dir / "backbone" / "backbone.pb"
		pb.parent.mkdir(parents=True, exist_ok=True)
		pb.write_text("pb", encoding="utf-8")
		return str(pb)
	monkeypatch.setattr(processor, "build_backbone_pb", fake_build_backbone_pb)

	def fake_place_onto_pb(input_pb, vcf_path, chunk_output_dir):
		chunk_ids = records["vcf_to_ids"].get(vcf_path, [])
		for seq_id in chunk_ids:
			if seq_id != ref_id:
				records["placed_ids"].append(seq_id)
		Path(chunk_output_dir).mkdir(parents=True, exist_ok=True)
		tips = list(existing_tips) + sorted(set(records["placed_ids"]))
		newick = "(" + ",".join(f"{tip}:0.2" for tip in tips) + ");\n"
		(Path(chunk_output_dir) / "uncondensed-final-tree.nh").write_text(newick, encoding="utf-8")
		(Path(chunk_output_dir) / "final-tree.nh").write_text(newick, encoding="utf-8")
		output_pb = Path(chunk_output_dir) / "usher.pb"
		output_pb.write_text("pb", encoding="utf-8")
		records["place"].append({
			"input_pb": input_pb,
			"vcf_path": vcf_path,
			"chunk_output_dir": str(chunk_output_dir),
			"output_pb": str(output_pb),
		})
		return str(output_pb)
	monkeypatch.setattr(processor, "place_onto_pb", fake_place_onto_pb)

	return records


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


def test_prepare_update_assets_prefers_usher_tree_over_iqtree(tmp_path: Path):
	db_path = tmp_path / "update.db"
	conn = sqlite3.connect(str(db_path))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, segment_key TEXT, newick TEXT)")
		cur.executemany(
			"INSERT INTO trees(name, source, segment, segment_key, newick) VALUES (?, ?, ?, ?, ?)",
			[
				("iqtree_seg1", "iqtree", "1", "seg_1", "(REF1:0.1,IQ_ONLY:0.2);"),
				("usher_seg1", "usher", "1", "seg_1", "(REF1:0.1,USHER_ONLY:0.2);"),
			],
		)
		conn.commit()
	finally:
		conn.close()

	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(">REF1\nATGC\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		update_db=str(db_path),
	)

	tree_path, _ = processor.prepare_update_assets()
	tree_text = Path(tree_path).read_text(encoding="utf-8")
	assert "USHER_ONLY" in tree_text
	assert "IQ_ONLY" not in tree_text


def test_prepare_update_assets_raises_when_segment_tree_missing(tmp_path: Path):
	db_path = tmp_path / "update.db"
	conn = sqlite3.connect(str(db_path))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, segment_key TEXT, newick TEXT)")
		cur.execute(
			"INSERT INTO trees(name, source, segment, segment_key, newick) VALUES (?, ?, ?, ?, ?)",
			("usher_seg2", "usher", "2", "seg_2", "(REF2:0.1);"),
		)
		conn.commit()
	finally:
		conn.close()

	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(">REF1\nATGC\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		update_db=str(db_path),
	)

	with pytest.raises(ValueError, match="Missing tree for segment 1"):
		processor.prepare_update_assets()


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


def test_split_alignment_into_chunks_falls_back_when_seqkit_errors(tmp_path: Path, monkeypatch):
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

	def fake_run(cmd, check=True):
		raise subprocess.CalledProcessError(1, cmd)

	monkeypatch.setattr("subprocess.run", fake_run)

	chunk_paths = processor.split_alignment_into_chunks(str(msa), "REF1", ["Q1", "Q2", "Q3"])

	assert len(chunk_paths) == 2
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

	def fail_backbone(*args, **kwargs):
		raise AssertionError("build_backbone_pb should not be called when only duplicate attachment is required")

	def fail_place(*args, **kwargs):
		raise AssertionError("place_onto_pb should not be called when only duplicate attachment is required")

	monkeypatch.setattr(processor, "build_backbone_pb", fail_backbone)
	monkeypatch.setattr(processor, "place_onto_pb", fail_place)

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
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1", "C1"])

	processor.run()

	placed = set(records["placed_ids"])
	assert placed == {"SWARM_A1", "SWARM_B1", "SWARM_C1", "UNIQ1"}

	# Only unique anchors ever entered a VCF; duplicate swarm members never did.
	all_chunk_ids = [seq_id for call in records["vcf"] for seq_id in call["ids"]]
	duplicate_members = {"SWARM_A2", "SWARM_A3", "SWARM_B2", "SWARM_C2", "SWARM_C3", "SWARM_C4"}
	assert not (duplicate_members & set(all_chunk_ids))
	# 4 anchors at chunk_size=2 -> 2 chunks/placements.
	assert len(records["vcf"]) == 2
	assert len(records["place"]) == 2

	tree_text = (output_dir / "uncondensed-final-tree.nh").read_text(encoding="utf-8")
	for seq_id in duplicate_members:
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
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1", "OLD1"])

	processor.run()

	# Backbone protobuf is built once, from the existing tree and its tips.
	assert len(records["backbone"]) == 1
	assert records["backbone"][0]["tree_file"] == str(seed_tree)
	assert records["backbone"][0]["backbone_ids"] == ["REF1", "OLD1"]

	# Only the two brand-new unique anchors are placed; existing tips and all
	# duplicates (old and new) are excluded from the VCF.
	assert len(records["vcf"]) == 1
	assert records["vcf"][0]["ids"] == ["REF1", "NEW1", "NEW3"]
	assert set(records["placed_ids"]) == {"NEW1", "NEW3"}

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


def test_run_update_mode_chunks_iteratively_and_carries_protobuf_forward(tmp_path: Path, monkeypatch):
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
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])

	processor.run()

	# 5 placeable at chunk_size=2 -> 3 chunks.
	assert len(records["vcf"]) == 3
	assert len(records["place"]) == 3

	backbone_pb = str(output_dir / "backbone" / "backbone.pb")
	# First placement runs against the freshly built backbone protobuf; each
	# subsequent placement runs against the protobuf produced by the prior chunk.
	assert records["place"][0]["input_pb"] == backbone_pb
	assert records["place"][1]["input_pb"] == records["place"][0]["output_pb"]
	assert records["place"][2]["input_pb"] == records["place"][1]["output_pb"]
	assert set(records["placed_ids"]) == {"Q1", "Q2", "Q3", "Q4", "Q5"}
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
	install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])

	processor.run()

	stdout = capsys.readouterr().out
	assert "batches complete: 0/3; remaining: 3." in stdout
	assert "batches complete: 1/3; remaining: 2." in stdout
	assert "batches complete: 2/3; remaining: 1." in stdout
	assert "batches complete: 3/3; remaining: 0." in stdout


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

	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1", "PLACED1"])

	processor.run()

	# The supplied starter tree and its tips seed the backbone protobuf.
	assert len(records["backbone"]) == 1
	assert records["backbone"][0]["tree_file"] == str(starter_tree)
	assert records["backbone"][0]["backbone_ids"] == ["REF1", "PLACED1"]

	# Only Q1 needs placing; Q2 is an identical duplicate attached afterwards.
	assert len(records["vcf"]) == 1
	assert records["vcf"][0]["ids"] == ["REF1", "Q1"]
	assert set(records["placed_ids"]) == {"Q1"}

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
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1", "C1", "C2"])

	processor.run()

	# Cluster reps C1/C2 are backbone tips (excluded); Q1..Q4 are placed in 2 chunks.
	assert len(records["vcf"]) == 2
	assert len(records["place"]) == 2
	assert records["backbone"][0]["tree_file"] == str(seed_tree)
	assert records["backbone"][0]["backbone_ids"] == ["REF1", "C1", "C2"]
	assert set(records["placed_ids"]) == {"Q1", "Q2", "Q3", "Q4"}
	# Backbone tips are never re-placed.
	all_chunk_ids = [seq_id for call in records["vcf"] for seq_id in call["ids"]]
	assert "C1" not in all_chunk_ids
	assert "C2" not in all_chunk_ids
	assert (output_dir / "uncondensed-final-tree.nh").exists()


# ---------------------------------------------------------------------------
# Focused tests for the mutation-annotated-backbone / protobuf helpers.
# ---------------------------------------------------------------------------


def test_place_onto_pb_writes_verbose_log_and_returns_updated_pb(tmp_path: Path, monkeypatch):
	output_dir = tmp_path / "usher_out"
	chunk_dir = output_dir / "chunk_0001"
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "input.fasta"),
		output_dir=str(output_dir),
		threads=4,
	)
	input_pb = tmp_path / "backbone.pb"
	input_pb.write_text("pb", encoding="utf-8")
	vcf_path = tmp_path / "all_samples.vcf"
	vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")

	run_calls = []

	class FakeCompleted:
		def __init__(self, stdout="", stderr=""):
			self.stdout = stdout
			self.stderr = stderr

	def fake_run(cmd, **kwargs):
		run_calls.append((cmd, kwargs))
		if cmd == ["usher", "--help"]:
			return FakeCompleted(stdout="usage: usher -T <threads>")
		kwargs["stdout"].write("verbose usher output\n")
		Path(chunk_dir).mkdir(parents=True, exist_ok=True)
		(chunk_dir / "uncondensed-final-tree.nh").write_text("(REF1:0.1,Q1:0.2);\n", encoding="utf-8")
		(chunk_dir / "usher.pb").write_text("pb", encoding="utf-8")
		return FakeCompleted()

	monkeypatch.setattr("subprocess.run", fake_run)

	result = processor.place_onto_pb(str(input_pb), str(vcf_path), str(chunk_dir))

	assert result == str(chunk_dir / "usher.pb")
	assert (chunk_dir / "usher.verbose.log").read_text(encoding="utf-8") == "verbose usher output\n"
	usher_cmd = run_calls[1][0]
	assert usher_cmd[:3] == ["usher", "-i", str(input_pb)]
	assert usher_cmd[-2:] == ["-T", "4"]
	assert run_calls[1][1]["stderr"] == subprocess.STDOUT


def test_place_onto_pb_raises_when_usher_produces_no_tree(tmp_path: Path, monkeypatch):
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "input.fasta"),
		output_dir=str(output_dir),
	)
	input_pb = tmp_path / "backbone.pb"
	input_pb.write_text("pb", encoding="utf-8")
	vcf_path = tmp_path / "chunk.vcf"
	vcf_path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")

	monkeypatch.setattr(processor, "_append_threads", lambda cmd: cmd)
	# USHER produces neither a tree nor a protobuf.
	monkeypatch.setattr(processor, "_run_usher_cmd", lambda *args, **kwargs: None)

	with pytest.raises(FileNotFoundError):
		processor.place_onto_pb(str(input_pb), str(vcf_path), str(output_dir / "chunk_0001"))


def test_build_backbone_pb_raises_when_usher_makes_no_protobuf(tmp_path: Path, monkeypatch):
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "input.fasta"),
		output_dir=str(output_dir),
	)
	tree_file = tmp_path / "seed.nwk"
	tree_file.write_text("(REF1:0.1);\n", encoding="utf-8")

	monkeypatch.setattr(processor, "collect_backbone_fasta", lambda *a, **k: str(tmp_path / "bb.fasta"))
	monkeypatch.setattr(processor, "build_vcf", lambda *a, **k: str(tmp_path / "bb.vcf"))
	monkeypatch.setattr(processor, "_append_threads", lambda cmd: cmd)
	# USHER runs but never writes the protobuf.
	monkeypatch.setattr(processor, "_run_usher_cmd", lambda *a, **k: None)

	with pytest.raises(FileNotFoundError):
		processor.build_backbone_pb(str(tree_file), ["REF1"], "REF1", str(tmp_path / "input.fasta"))


def test_resolve_usher_tree_output_prefers_final_tree(tmp_path: Path):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	(output_dir / "uncondensed-final-tree.nh").write_text("(A);\n", encoding="utf-8")
	assert UsherPlacement._resolve_usher_tree_output(str(output_dir)).endswith("uncondensed-final-tree.nh")

	(output_dir / "final-tree.nh").write_text("(A);\n", encoding="utf-8")
	assert UsherPlacement._resolve_usher_tree_output(str(output_dir)).endswith("final-tree.nh")


def test_resolve_usher_tree_output_raises_when_missing(tmp_path: Path):
	output_dir = tmp_path / "empty"
	output_dir.mkdir(parents=True, exist_ok=True)
	with pytest.raises(FileNotFoundError):
		UsherPlacement._resolve_usher_tree_output(str(output_dir))


def test_append_threads_probes_usher_once_and_caches(tmp_path: Path, monkeypatch):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "a.fasta"),
		output_dir=str(tmp_path / "out"),
		threads=8,
	)

	help_calls = []

	class FakeCompleted:
		def __init__(self, stdout="", stderr=""):
			self.stdout = stdout
			self.stderr = stderr

	def fake_run(cmd, **kwargs):
		help_calls.append(cmd)
		return FakeCompleted(stdout="usage: usher -T <n>")

	monkeypatch.setattr("subprocess.run", fake_run)

	cmd1 = processor._append_threads(["usher", "-i", "x"])
	cmd2 = processor._append_threads(["usher", "-i", "y"])

	assert cmd1[-2:] == ["-T", "8"]
	assert cmd2[-2:] == ["-T", "8"]
	# The --help probe is cached across calls.
	assert len(help_calls) == 1


def test_append_threads_omits_flag_when_unsupported(tmp_path: Path, monkeypatch):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "a.fasta"),
		output_dir=str(tmp_path / "out"),
		threads=8,
	)

	class FakeCompleted:
		def __init__(self, stdout="", stderr=""):
			self.stdout = stdout
			self.stderr = stderr

	monkeypatch.setattr("subprocess.run", lambda cmd, **kwargs: FakeCompleted(stdout="usage: usher [options]"))

	cmd = processor._append_threads(["usher", "-i", "x"])
	assert "-T" not in cmd


def test_run_usher_cmd_sets_thread_env_and_captures_output(tmp_path: Path, monkeypatch):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "a.fasta"),
		output_dir=str(output_dir),
		threads=3,
	)

	captured = {}

	def fake_run(cmd, **kwargs):
		captured["env"] = kwargs["env"]
		captured["stderr"] = kwargs["stderr"]
		kwargs["stdout"].write("hello\n")
		return None

	monkeypatch.setattr("subprocess.run", fake_run)

	processor._run_usher_cmd(["usher"], str(output_dir), log_name="probe.log")

	assert (output_dir / "probe.log").read_text(encoding="utf-8") == "hello\n"
	assert captured["stderr"] == subprocess.STDOUT
	for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
		assert captured["env"][var] == "3"


def test_collect_backbone_fasta_pulls_backbone_tips_from_db(tmp_path: Path):
	db_path = tmp_path / "update.db"
	conn = sqlite3.connect(str(db_path))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment TEXT)")
		cur.executemany(
			"INSERT INTO sequence_alignment(primary_accession, alignment) VALUES (?, ?)",
			[("OLD1", "CCCC"), ("OLD2", "GGGG")],
		)
		conn.commit()
	finally:
		conn.close()

	# The working alignment only holds the reference; backbone tips must be
	# rehydrated from the update DB or the backbone branches carry no mutations.
	alignment = tmp_path / "aln.fasta"
	alignment.write_text(">REF1\nATGC\n", encoding="utf-8")
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)

	processor = UsherPlacement(
		padded_aln=str(alignment),
		output_dir=str(output_dir),
		update_db=str(db_path),
	)

	out_path = processor.collect_backbone_fasta(["OLD1", "OLD2"], "REF1", str(alignment))
	ids = processor._read_ids_from_fasta(out_path)
	assert set(ids) == {"REF1", "OLD1", "OLD2"}


def test_collect_backbone_fasta_raises_when_reference_absent(tmp_path: Path):
	alignment = tmp_path / "aln.fasta"
	alignment.write_text(">OLD1\nCCCC\n", encoding="utf-8")
	output_dir = tmp_path / "usher_out"
	output_dir.mkdir(parents=True, exist_ok=True)

	processor = UsherPlacement(
		padded_aln=str(alignment),
		output_dir=str(output_dir),
	)

	with pytest.raises(ValueError, match="Reference sequence 'REF1' not found"):
		processor.collect_backbone_fasta(["OLD1"], "REF1", str(alignment))


def test_fetch_db_alignments_skips_missing_empty_and_nan(tmp_path: Path):
	db_path = tmp_path / "update.db"
	conn = sqlite3.connect(str(db_path))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment TEXT)")
		cur.executemany(
			"INSERT INTO sequence_alignment(primary_accession, alignment) VALUES (?, ?)",
			[("A", "AAAA"), ("B", ""), ("C", "nan"), ("D", None), ("E", "CCCC")],
		)
		conn.commit()
	finally:
		conn.close()

	processor = UsherPlacement(
		padded_aln=str(tmp_path / "aln.fasta"),
		output_dir=str(tmp_path / "out"),
		update_db=str(db_path),
	)

	result = processor._fetch_db_alignments(["A", "B", "C", "D", "E", "MISSING"])
	assert result == {"A": "AAAA", "E": "CCCC"}


def test_fetch_db_alignments_returns_empty_without_db(tmp_path: Path):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "aln.fasta"),
		output_dir=str(tmp_path / "out"),
	)
	assert processor._fetch_db_alignments(["A", "B"]) == {}


def test_build_vcf_refuses_to_drop_the_exclude_filter_when_faToVcf_fails(tmp_path: Path, monkeypatch):
	"""This used to retry without -excludeFile and return a VCF. That was wrong:
	-excludeFile is what keeps backbone tips out of the placement VCF, so the
	silent retry re-submitted every existing tree tip for placement and
	duplicated them in the tree. A failure to exclude must abort."""
	msa = tmp_path / "aln.fasta"
	msa.write_text(">REF1\nATGC\n>Q1\nATGA\n", encoding="utf-8")
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	exclude = tmp_path / "exclude.txt"
	exclude.write_text("Q2\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
	)

	cmds = []

	def fake_run(cmd, check=True):
		cmds.append(cmd)
		if any(str(arg).startswith("-excludeFile=") for arg in cmd):
			raise subprocess.CalledProcessError(1, cmd)
		return None

	monkeypatch.setattr("subprocess.run", fake_run)

	vcf_path = output_dir / "out.vcf"
	with pytest.raises(RuntimeError, match="Refusing to retry"):
		processor.build_vcf(
			"REF1",
			exclude_ids_file=str(exclude),
			alignment_fasta=str(msa),
			vcf_path=str(vcf_path),
		)

	# Exactly one attempt: no unfiltered second call was made.
	assert len(cmds) == 1
	assert any(str(arg).startswith("-excludeFile=") for arg in cmds[0])


def test_build_vcf_uses_exclude_file_when_successful(tmp_path: Path, monkeypatch):
	msa = tmp_path / "aln.fasta"
	msa.write_text(">REF1\nATGC\n>Q1\nATGA\n", encoding="utf-8")
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	exclude = tmp_path / "exclude.txt"
	exclude.write_text("Q2\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(output_dir),
	)

	cmds = []
	monkeypatch.setattr("subprocess.run", lambda cmd, check=True: cmds.append(cmd))

	vcf_path = output_dir / "out.vcf"
	result = processor.build_vcf(
		"REF1",
		exclude_ids_file=str(exclude),
		alignment_fasta=str(msa),
		vcf_path=str(vcf_path),
	)

	assert result == str(vcf_path)
	# Only one call: the exclude-file attempt succeeded, no retry.
	assert len(cmds) == 1
	assert any(str(arg).startswith("-excludeFile=") for arg in cmds[0])


def test_promote_final_usher_outputs_copies_present_files(tmp_path: Path):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	final_chunk = output_dir / "chunk_0002"
	final_chunk.mkdir(parents=True, exist_ok=True)
	(final_chunk / "final-tree.nh").write_text("(A);\n", encoding="utf-8")
	(final_chunk / "uncondensed-final-tree.nh").write_text("(A);\n", encoding="utf-8")
	(final_chunk / "usher.pb").write_text("pb", encoding="utf-8")
	# mutation-paths.txt / placement_stats.tsv / all_samples.vcf intentionally absent.

	processor = UsherPlacement(
		padded_aln=str(tmp_path / "aln.fasta"),
		output_dir=str(output_dir),
	)

	processor.promote_final_usher_outputs(str(final_chunk))

	assert (output_dir / "final-tree.nh").exists()
	assert (output_dir / "uncondensed-final-tree.nh").exists()
	assert (output_dir / "usher.pb").exists()
	assert not (output_dir / "mutation-paths.txt").exists()
	assert not (output_dir / "all_samples.vcf").exists()


# ---------------------------------------------------------------------------
# Path / value normalization helpers.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value", [None, "", "   ", "null", "NULL", "Null", "UNSET"])
def test_normalize_optional_path_returns_none_for_sentinels(value):
	assert UsherPlacement._normalize_optional_path(value) is None


def test_normalize_optional_path_strips_real_paths():
	assert UsherPlacement._normalize_optional_path("  /data/update.db  ") == "/data/update.db"


def test_constructor_treats_sentinel_paths_as_unset(tmp_path: Path):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "aln.fasta"),
		output_dir=str(tmp_path / "out"),
		update_db="null",
		mmseq_cluster_dir="UNSET",
		iqtree_dir="   ",
		starter_tree=None,
	)
	assert processor.update_db is None
	assert processor.mmseq_cluster_dir is None
	assert processor.iqtree_dir is None
	assert processor.starter_tree is None


def test_constructor_clamps_threads_and_chunk_sizes(tmp_path: Path):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "aln.fasta"),
		output_dir=str(tmp_path / "out"),
		threads=0,
		chunk_size=0,
		chunk_threshold=-5,
	)
	assert processor.threads == 1
	assert processor.chunk_size == 1
	assert processor.chunk_threshold == 1


@pytest.mark.parametrize(
	"value,expected",
	[
		(None, "0"),
		("", "0"),
		("  ", "0"),
		("1", "1"),
		("segment_2", "2"),
		("seg", "seg"),
	],
)
def test_normalize_segment(value, expected):
	assert UsherPlacement._normalize_segment(value) == expected


@pytest.mark.parametrize(
	"filename,expected",
	[
		("segment_1_dedup.fasta", "1"),
		("segment_12_dedup.fasta", "12"),
		("no_digits.fasta", "0"),
	],
)
def test_segment_from_padded_alignment(tmp_path: Path, filename, expected):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / filename),
		output_dir=str(tmp_path / "out"),
	)
	assert processor._segment_from_padded_alignment() == expected


def test_read_tree_terminals_handles_empty_file(tmp_path: Path):
	empty_tree = tmp_path / "empty.nwk"
	empty_tree.write_text("   \n", encoding="utf-8")
	assert UsherPlacement._read_tree_terminals(str(empty_tree)) == set()


def test_run_raises_when_alignment_has_single_sequence(tmp_path: Path, monkeypatch):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(">REF1\nATGC\n", encoding="utf-8")
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
	)

	monkeypatch.setattr(processor, "resolve_non_update_assets", lambda: (str(cluster_rep), str(seed_tree)))

	with pytest.raises(ValueError, match="Need at least 2 sequences"):
		processor.run()


def test_resolve_reference_id_falls_back_when_cluster_rep_absent_from_alignment(tmp_path: Path):
	cluster_rep = tmp_path / "cluster_rep.fasta"
	cluster_rep.write_text(">GHOST\nATGC\n", encoding="utf-8")
	msa = tmp_path / "aln.fasta"
	msa.write_text(">REF1\nATGC\n>Q1\nATGA\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "out"),
	)

	# GHOST is not in the alignment, so the resolver falls back to the first MSA id.
	ref_id = processor.resolve_reference_id(cluster_rep=str(cluster_rep), alignment_fasta=str(msa))
	assert ref_id == "REF1"


def test_resolve_reference_id_raises_on_empty_alignment(tmp_path: Path):
	msa = tmp_path / "empty.fasta"
	msa.write_text("", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "out"),
	)

	with pytest.raises(ValueError, match="Could not resolve reference ID"):
		processor.resolve_reference_id(alignment_fasta=str(msa))


def test_resolve_resume_assets_requires_both_files(tmp_path: Path):
	starter_tree = tmp_path / "tree.nwk"
	starter_tree.write_text("(REF1:0.1);\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(tmp_path / "aln.fasta"),
		output_dir=str(tmp_path / "out"),
		starter_tree=str(starter_tree),
		existing_ids_file=None,
	)

	with pytest.raises(ValueError, match="starter_tree and existing_ids_file"):
		processor.resolve_resume_assets()


# ===========================================================================
# Adversarial coverage
#
# Everything below attacks a boundary rather than a happy path. Each test began
# as a defect reproduced against the code as it stood, or as a guarantee the
# module claimed but never asserted.
# ===========================================================================


# ---------------------------------------------------------------------------
# 1. Malformed / hostile FASTA input
# ---------------------------------------------------------------------------


def test_read_ids_from_fasta_survives_a_bare_header(tmp_path: Path):
	"""A header that is a bare '>' - an empty description, which NCBI downloads
	and several aligners emit - used to raise IndexError out of .split()[0] and
	take down resolve_reference_id before any placement work started."""
	fasta = tmp_path / "bare.fasta"
	fasta.write_text(">\nATGC\n>Q1\nATGA\n", encoding="utf-8")
	assert UsherPlacement._read_ids_from_fasta(str(fasta)) == ["Q1"]


def test_read_ids_from_fasta_tolerates_non_utf8_bytes(tmp_path: Path):
	"""One non-UTF-8 byte in one public-database header used to abort the whole
	run with a UnicodeDecodeError naming a byte offset rather than a record."""
	fasta = tmp_path / "latin1.fasta"
	fasta.write_bytes(b">Q\xe9_ISOLATE\nATGC\n>Q2\nATGA\n")
	assert UsherPlacement._read_ids_from_fasta(str(fasta))[-1] == "Q2"


def test_read_ids_from_fasta_takes_only_the_first_header_token(tmp_path: Path):
	"""Locks the id convention shared with Bio.SeqIO's record.id: everything after
	the first whitespace is a description and must not become part of the id."""
	fasta = tmp_path / "desc.fasta"
	fasta.write_text(">Q1 Rabies lyssavirus strain X\nATGC\n", encoding="utf-8")
	assert UsherPlacement._read_ids_from_fasta(str(fasta)) == ["Q1"]


def test_build_update_alignment_input_handles_alignment_without_trailing_newline(tmp_path: Path):
	"""The working alignment was copied with a raw read() and rehydrated records
	appended after it. With no trailing newline the first rehydrated header was
	glued onto the last sequence line, so that record silently disappeared AND
	the sequence before it was corrupted with header text."""
	aln = tmp_path / "segment_1_dedup.fasta"
	aln.write_text(">REF1\nATGC\n>Q_NEW\nATGA", encoding="utf-8")  # no trailing newline
	tree = tmp_path / "seed.nwk"
	tree.write_text("(REF1:0.1);\n", encoding="utf-8")

	processor = make(tmp_path, padded_aln=str(aln))
	processor._load_update_missing_alignment_rows = lambda *a, **k: [("Q_OLD", "ATGT")]

	merged = processor.build_update_alignment_input(str(tree))
	assert UsherPlacement._read_ids_from_fasta(merged) == ["REF1", "Q_NEW", "Q_OLD"]


def test_build_update_alignment_input_rejects_fasta_injection_from_db_content(tmp_path: Path):
	"""Rehydrated DB alignments are written with an f-string. A stored alignment
	containing a newline followed by '>' - trivially produced by a bad import into
	sequence_alignment - injected extra FASTA records, putting accessions that
	exist in no metadata table into the published tree."""
	aln = tmp_path / "segment_1_dedup.fasta"
	aln.write_text(">REF1\nATGC\n", encoding="utf-8")
	tree = tmp_path / "seed.nwk"
	tree.write_text("(REF1:0.1);\n", encoding="utf-8")

	processor = make(tmp_path, padded_aln=str(aln))
	processor._load_update_missing_alignment_rows = lambda *a, **k: [
		("Q_OLD", "ATGT\n>PHANTOM\nGGGG")
	]

	with pytest.raises(ValueError, match="not a sequence"):
		processor.build_update_alignment_input(str(tree))


def test_collapse_identical_sequences_is_case_insensitive(tmp_path: Path):
	"""Bucketing on the raw str(record.seq) treated a soft-masked copy as a
	distinct genotype. _informative_length upper()s before counting, so the two
	halves of the module disagreed about what case means."""
	aln = tmp_path / "case.fasta"
	write_fasta(aln, [("REF1", "AAAA"), ("Q1", "CCCC"), ("Q2", "cccc"), ("Q3", "CCCC")])
	processor = make(tmp_path, padded_aln=str(aln))

	plan = processor.collapse_identical_sequences(str(aln), "REF1", set())

	assert plan["placeable_ids"] == ["Q1"]
	assert sorted(plan["anchor_to_members"]["Q1"]) == ["Q2", "Q3"]


def test_collapse_identical_sequences_deduplicates_repeated_accessions(tmp_path: Path):
	"""A repeated accession entered placeable_ids twice and put two VCF sample
	columns with one name in front of usher; len(placeable_ids) also over-counted
	the batch total reported to the user."""
	aln = tmp_path / "dup.fasta"
	write_fasta(aln, [("REF1", "AAAA"), ("Q1", "CCCC"), ("Q1", "GGGG"), ("Q2", "TTTT")])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=2)

	plan = processor.collapse_identical_sequences(str(aln), "REF1", set())
	assert plan["placeable_ids"] == ["Q1", "Q2"]

	chunks = processor._split_alignment_into_chunks_python(str(aln), "REF1", plan["placeable_ids"])
	for chunk in chunks:
		ids = UsherPlacement._read_ids_from_fasta(chunk)
		assert len(ids) == len(set(ids)), f"{chunk} has duplicate sample names: {ids}"


def test_uninformative_sequences_are_not_submitted_for_placement(tmp_path: Path):
	"""A sequence with zero informative bases says nothing about where it belongs,
	yet it still became a placeable anchor and usher attached it arbitrarily.
	All-N and all-gap records also formed two separate groups."""
	aln = tmp_path / "uninformative.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("REAL", "ACGA"), ("ALL_N", "NNNN"), ("ALL_GAP", "----")])
	processor = make(tmp_path, padded_aln=str(aln))

	plan = processor.collapse_identical_sequences(str(aln), "REF1", set())
	assert plan["placeable_ids"] == ["REAL"]


def test_informative_length_ignores_ambiguity_codes_gaps_and_case():
	"""The quality score must count only A/C/G/T/U, in either case."""
	assert UsherPlacement._informative_length("ACGT") == 4
	assert UsherPlacement._informative_length("acgt") == 4
	assert UsherPlacement._informative_length("ACGU") == 4
	assert UsherPlacement._informative_length("NNNN") == 0
	assert UsherPlacement._informative_length("----") == 0
	# Every IUPAC ambiguity code is uninformative for placement.
	assert UsherPlacement._informative_length("RYSWKMBDHVN") == 0
	assert UsherPlacement._informative_length("ACGT-N?ACGT") == 8


# ---------------------------------------------------------------------------
# 2. Segment inference
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
	"filename,expected",
	[
		("h3n2_segment_4_dedup.fasta", "4"),
		("IAV_H1N1_segment_2_dedup.fasta", "2"),
		("rabv_2024_segment_1_dedup.fasta", "1"),
	],
)
def test_segment_from_padded_alignment_reads_the_segment_field(tmp_path: Path, filename, expected):
	"""This was re.search(r'(\\d+)') over the basename, taking the FIRST digit run,
	so a subtype ('h3n2' -> 3), a serotype ('H1N1' -> 1) or a year ('2024' -> 2024)
	won over the real segment. In update mode that segment selects the starter
	tree, so a segment-4 alignment was silently placed onto the segment-3 tree -
	the same digit-scraping mistake segment_utils exists to end."""
	path = tmp_path / filename
	path.parent.mkdir(parents=True, exist_ok=True)
	processor = UsherPlacement(padded_aln=str(path), output_dir=str(tmp_path / "out"))
	assert processor._segment_from_padded_alignment() == expected


@pytest.mark.parametrize(
	"filename,expected",
	[
		("segment_1_dedup.fasta", "1"),
		("segment_12_dedup.fasta", "12"),
		("run3/segment_5_dedup.fasta", "5"),
		("no_digits_here.fasta", "0"),
	],
)
def test_segment_from_padded_alignment_is_correct_without_confounding_digits(tmp_path: Path, filename, expected):
	"""The parse is only safe when the basename carries no digits before the
	segment field - note the enclosing directory is ignored, so 'run3/' is fine."""
	path = tmp_path / filename
	path.parent.mkdir(parents=True, exist_ok=True)
	processor = UsherPlacement(padded_aln=str(path), output_dir=str(tmp_path / "out"))
	assert processor._segment_from_padded_alignment() == expected


# ---------------------------------------------------------------------------
# 3. Update-DB robustness
# ---------------------------------------------------------------------------


def _db_with_trees(path: Path, rows):
	conn = sqlite3.connect(str(path))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, newick TEXT)")
		cur.executemany("INSERT INTO trees(name, source, segment, newick) VALUES (?, ?, ?, ?)", rows)
		conn.commit()
	finally:
		conn.close()
	return path


def test_update_mode_reports_a_missing_metadata_table_cleanly(tmp_path: Path):
	"""prepare_update_assets guards the trees table with a sqlite_master lookup,
	but the alignment rehydration went straight to pd.read_sql_query, so a DB
	built by an older pipeline version died mid-run with a pandas DatabaseError
	instead of the error type every other missing asset here raises."""
	db = _db_with_trees(tmp_path / "notables.db", [("t", "usher", "1", "(REF1:0.1);")])
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ATGC"), ("Q1", "ATGA")])

	processor = make(tmp_path, padded_aln=str(aln), update_db=str(db))
	tree_path, _ = processor.prepare_update_assets()

	with pytest.raises((ValueError, FileNotFoundError)):
		processor.build_update_alignment_input(tree_path)


def test_prepare_update_assets_breaks_source_ties_deterministically(tmp_path: Path):
	"""Two rows sharing the winning source priority were resolved by keeping
	whichever sqlite returned first, and the query had no ORDER BY. Row order
	after a VACUUM or a differently-ordered insert is not stable, so two runs
	against equivalent DBs could start from different backbone trees."""
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ATGC")])

	first = _db_with_trees(tmp_path / "a.db", [
		("older", "usher", "1", "(REF1:0.1,OLDER:0.2);"),
		("newer", "usher", "1", "(REF1:0.1,NEWER:0.2);"),
	])
	second = _db_with_trees(tmp_path / "b.db", [
		("newer", "usher", "1", "(REF1:0.1,NEWER:0.2);"),
		("older", "usher", "1", "(REF1:0.1,OLDER:0.2);"),
	])

	picked = []
	for idx, db in enumerate((first, second)):
		processor = make(tmp_path, padded_aln=str(aln),
						 output_dir=str(tmp_path / f"out{idx}"), update_db=str(db))
		tree_path, _ = processor.prepare_update_assets()
		picked.append("NEWER" in Path(tree_path).read_text(encoding="utf-8"))

	assert picked[0] == picked[1], "tree choice depends on sqlite row order"


def test_prepare_update_assets_ignores_trees_from_other_segments(tmp_path: Path):
	"""A segment_key/segment mismatch must never let a foreign tree through."""
	db = _db_with_trees(tmp_path / "mixed.db", [
		("seg2", "usher", "2", "(REF2:0.1,OTHER:0.2);"),
		("seg1", "iqtree", "1", "(REF1:0.1,MINE:0.2);"),
	])
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ATGC")])

	processor = make(tmp_path, padded_aln=str(aln), update_db=str(db))
	tree_path, _ = processor.prepare_update_assets()
	text = Path(tree_path).read_text(encoding="utf-8")
	assert "MINE" in text and "OTHER" not in text


def test_fetch_db_alignments_batches_beyond_the_sqlite_variable_limit(tmp_path: Path):
	"""The 900-per-statement batching must survive an accession list far larger
	than SQLITE_MAX_VARIABLE_NUMBER (999 on stock builds)."""
	db = tmp_path / "many.db"
	conn = sqlite3.connect(str(db))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment TEXT)")
		cur.executemany(
			"INSERT INTO sequence_alignment(primary_accession, alignment) VALUES (?, ?)",
			[(f"ACC{i:05d}", "ACGT") for i in range(2500)],
		)
		conn.commit()
	finally:
		conn.close()

	processor = make(tmp_path, update_db=str(db))
	result = processor._fetch_db_alignments([f"ACC{i:05d}" for i in range(2500)])
	assert len(result) == 2500


# ---------------------------------------------------------------------------
# 4. Placement ordering - the module's central correctness claim
# ---------------------------------------------------------------------------


def test_order_ids_by_quality_ranks_by_informative_bases(tmp_path: Path):
	aln = tmp_path / "q.fasta"
	write_fasta(aln, [
		("REF1", "ACGTACGTAC"),
		("PARTIAL", "AC" + "N" * 8),
		("FULL", "ACGTACGTAC"),
		("HALF", "ACGTAC" + "N" * 4),
	])
	processor = make(tmp_path, padded_aln=str(aln))

	ordered, scores = processor.order_ids_by_quality(str(aln), ["PARTIAL", "FULL", "HALF"])

	assert ordered == ["FULL", "HALF", "PARTIAL"]
	assert scores == {"FULL": 10, "HALF": 6, "PARTIAL": 2}


def test_order_ids_by_quality_breaks_ties_on_accession_for_reproducibility(tmp_path: Path):
	aln = tmp_path / "ties.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("ZZZ", "ACGT"), ("AAA", "ACGT"), ("MMM", "ACGT")])
	processor = make(tmp_path, padded_aln=str(aln))

	first, _ = processor.order_ids_by_quality(str(aln), ["ZZZ", "AAA", "MMM"])
	second, _ = processor.order_ids_by_quality(str(aln), ["MMM", "ZZZ", "AAA"])

	assert first == second == ["AAA", "MMM", "ZZZ"]


def test_order_ids_by_quality_sorts_ids_absent_from_the_alignment_last(tmp_path: Path):
	aln = tmp_path / "ghost.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("REAL", "ACGA"), ("EMPTY", "NNNN")])
	processor = make(tmp_path, padded_aln=str(aln))

	ordered, scores = processor.order_ids_by_quality(str(aln), ["REAL", "GHOST", "EMPTY"])

	assert ordered[-1] == "GHOST"
	assert "GHOST" not in scores
	# An id that is present but scores 0 must still outrank a wholly absent one.
	assert ordered == ["REAL", "EMPTY", "GHOST"]


def test_python_split_preserves_rank_order_within_and_across_batches(tmp_path: Path):
	"""The docstring promises rank order is the literal placement order 'end to
	end'. The splitter streams records in ALIGNMENT order, so the rewrite pass
	that restores rank order inside each batch is what makes that true - assert
	it against an alignment stored in the exact reverse of quality order."""
	aln = tmp_path / "reversed.fasta"
	write_fasta(aln, [
		("REF1", "ACGTACGTAC"),
		("WORST", "A" + "N" * 9),
		("BAD", "ACG" + "N" * 7),
		("GOOD", "ACGTACG" + "N" * 3),
		("BEST", "ACGTACGTAC"),
	])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=2)

	ordered, _ = processor.order_ids_by_quality(str(aln), ["WORST", "BAD", "GOOD", "BEST"])
	assert ordered == ["BEST", "GOOD", "BAD", "WORST"]

	chunks = processor._split_alignment_into_chunks_python(str(aln), "REF1", ordered)

	assert len(chunks) == 2
	assert UsherPlacement._read_ids_from_fasta(chunks[0]) == ["REF1", "BEST", "GOOD"]
	assert UsherPlacement._read_ids_from_fasta(chunks[1]) == ["REF1", "BAD", "WORST"]


def test_placement_order_report_batch_column_matches_the_real_batching(tmp_path: Path):
	aln = tmp_path / "report.fasta"
	write_fasta(aln, [("REF1", "ACGT")] + [(f"Q{i}", "ACG" + "ACGT"[i % 4]) for i in range(5)])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=2)

	ids = [f"Q{i}" for i in range(5)]
	ordered, scores = processor.order_ids_by_quality(str(aln), ids)
	report = Path(processor.write_placement_order_report(ordered, scores))
	chunks = processor._split_alignment_into_chunks_python(str(aln), "REF1", ordered)

	rows = [line.split("\t") for line in report.read_text(encoding="utf-8").strip().splitlines()[1:]]
	assert [r[0] for r in rows] == ["1", "2", "3", "4", "5"]
	assert [r[1] for r in rows] == ordered

	for acc, batch in ((r[1], int(r[3])) for r in rows):
		assert acc in UsherPlacement._read_ids_from_fasta(chunks[batch - 1]), \
			f"{acc} is reported in batch {batch} but is not in chunk {batch}"


def test_run_places_highest_quality_sequences_in_the_first_batch(tmp_path: Path, monkeypatch):
	"""End-to-end: the first VCF handed to usher must hold the best sequences,
	and every later batch places onto the protobuf the previous one produced."""
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [
		("REF1", "ACGTACGTAC"),
		("JUNK", "AC" + "N" * 8),
		("BEST", "ACGTACGTAA"),
		("MEH", "ACGTA" + "N" * 5),
		("GOOD", "ACGTACGT" + "NN"),
	])
	output_dir = tmp_path / "out"
	starter = tmp_path / "seed.nwk"
	starter.write_text("(REF1:0.1);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir),
		chunk_size=2, chunk_threshold=3,
		starter_tree=str(starter), existing_ids_file=str(existing),
	)
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])

	processor.run()

	assert [call["ids"] for call in records["vcf"]] == [
		["REF1", "BEST", "GOOD"],
		["REF1", "MEH", "JUNK"],
	]
	assert records["place"][1]["input_pb"] == records["place"][0]["output_pb"]


def test_run_with_input_order_skips_ranking_and_writes_no_order_report(tmp_path: Path, monkeypatch):
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("JUNK", "ANNN"), ("BEST", "ACGA")])
	output_dir = tmp_path / "out"
	starter = tmp_path / "seed.nwk"
	starter.write_text("(REF1:0.1);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir), chunk_size=5,
		starter_tree=str(starter), existing_ids_file=str(existing),
		placement_order="input",
	)
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])

	processor.run()

	assert records["vcf"][0]["ids"] == ["REF1", "JUNK", "BEST"]
	assert not (output_dir / "placement_order.tsv").exists()


@pytest.mark.parametrize("bad", ["QUALITY_DESC", "random", "1", "  "])
def test_constructor_rejects_an_unknown_placement_order(tmp_path: Path, bad):
	with pytest.raises(ValueError, match="placement_order"):
		UsherPlacement(
			padded_aln=str(tmp_path / "a.fasta"),
			output_dir=str(tmp_path / "out"),
			placement_order=bad,
		)


@pytest.mark.parametrize("value,expected", [("QUALITY", "quality"), ("  Input  ", "input"), (None, "quality")])
def test_constructor_normalises_placement_order(tmp_path: Path, value, expected):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "a.fasta"),
		output_dir=str(tmp_path / "out"),
		placement_order=value,
	)
	assert processor.placement_order == expected


# ---------------------------------------------------------------------------
# 5. Chunking: stale state, empty batches, ordering of parts
# ---------------------------------------------------------------------------


def test_python_split_does_not_emit_a_reference_only_chunk(tmp_path: Path):
	"""total_chunks was derived from len(placeable_ids) while only ids found in
	the alignment were streamed into a chunk file, so an id present in the tree
	bookkeeping but absent from the alignment produced a trailing chunk holding
	nothing but the reference - still handed to faToVcf and usher as a batch."""
	aln = tmp_path / "ghost.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("Q1", "ACGA"), ("Q2", "ACGG")])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=2)

	chunks = processor._split_alignment_into_chunks_python(str(aln), "REF1", ["Q1", "Q2", "GHOST"])

	for chunk in chunks:
		ids = UsherPlacement._read_ids_from_fasta(chunk)
		assert ids != ["REF1"], f"{chunk} contains the reference and nothing to place"


def test_seqkit_split_ignores_leftover_parts_from_a_previous_run(tmp_path: Path, monkeypatch):
	"""The seqkit path collected its parts with os.listdir without clearing the
	directory first. Re-running a segment into an existing output directory -
	exactly what resume mode does - made leftover part files from the previous run
	become chunk_0001, shifting the real batches down and feeding usher sequences
	the current run never selected."""
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("Q1", "ACGA"), ("Q2", "ACGG")])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=2)

	raw_dir = Path(processor.output_dir) / "chunks" / "raw"
	raw_dir.mkdir(parents=True, exist_ok=True)
	(raw_dir / "placeable_only.part_001.fasta").write_text(">STALE\nAAAA\n", encoding="utf-8")

	monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/seqkit")

	def fake_run(cmd, check=True):
		raw_dir.mkdir(parents=True, exist_ok=True)
		if cmd[1] == "grep":
			Path(cmd[-1]).write_text(">Q1\nACGA\n>Q2\nACGG\n", encoding="utf-8")
		elif cmd[1] == "split2":
			(raw_dir / "placeable_only.part_002.fasta").write_text(
				">Q1\nACGA\n>Q2\nACGG\n", encoding="utf-8")
		return None
	monkeypatch.setattr("subprocess.run", fake_run)

	chunks = processor.split_alignment_into_chunks(str(aln), "REF1", ["Q1", "Q2"], ordered=False)

	placed = {i for c in chunks for i in UsherPlacement._read_ids_from_fasta(c)}
	assert "STALE" not in placed


def test_seqkit_part_files_are_ordered_numerically_not_lexicographically(tmp_path: Path, monkeypatch):
	"""sorted(os.listdir(...)) sorts seqkit's part names as text. seqkit pads to
	three digits, so once a split produces 1000+ parts 'part_1000' sorts next to
	'part_100' and ahead of 'part_999'. Batches then reached usher out of order,
	destroying the placement ranking and desynchronising chunk_NNNN directories
	from their contents."""
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ACGT")] + [(f"Q{i:04d}", "ACGA") for i in range(3)])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=1)

	raw_dir = Path(processor.output_dir) / "chunks" / "raw"
	monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/seqkit")

	def fake_run(cmd, check=True):
		raw_dir.mkdir(parents=True, exist_ok=True)
		if cmd[1] == "grep":
			Path(cmd[-1]).write_text(">Q0000\nACGA\n", encoding="utf-8")
		elif cmd[1] == "split2":
			for part, seq_id in ((999, "FIRST"), (1000, "SECOND"), (1001, "THIRD")):
				(raw_dir / f"placeable_only.part_{part:03d}.fasta").write_text(
					f">{seq_id}\nACGA\n", encoding="utf-8")
		return None
	monkeypatch.setattr("subprocess.run", fake_run)

	chunks = processor.split_alignment_into_chunks(str(aln), "REF1", ["Q0000"], ordered=False)
	order = [UsherPlacement._read_ids_from_fasta(c)[1] for c in chunks]
	assert order == ["FIRST", "SECOND", "THIRD"]


def test_python_split_stays_correct_past_the_open_file_handle_cap(tmp_path: Path):
	"""max_open=256 evicts and reopens handles in append mode; every record must
	still land in exactly one batch, with the reference at the head of each."""
	count = 600
	aln = tmp_path / "many.fasta"
	write_fasta(aln, [("REF1", "ACGT")] + [(f"Q{i:04d}", "ACGA" if i % 2 else "ACGG") for i in range(count)])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=1)

	ids = [f"Q{i:04d}" for i in range(count)]
	chunks = processor._split_alignment_into_chunks_python(str(aln), "REF1", ids)

	assert len(chunks) == count
	seen = []
	for idx, chunk in enumerate(chunks):
		chunk_ids = UsherPlacement._read_ids_from_fasta(chunk)
		assert chunk_ids[0] == "REF1"
		assert chunk_ids[1:] == [ids[idx]]
		seen.extend(chunk_ids[1:])
	assert seen == ids


def test_split_alignment_into_chunks_raises_when_the_reference_is_absent(tmp_path: Path):
	aln = tmp_path / "noref.fasta"
	write_fasta(aln, [("Q1", "ACGA"), ("Q2", "ACGG")])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=1)

	with pytest.raises(ValueError, match="Reference ID REF1 not found"):
		processor._split_alignment_into_chunks_python(str(aln), "REF1", ["Q1", "Q2"])


def test_split_alignment_into_chunks_returns_nothing_for_an_empty_id_list(tmp_path: Path):
	aln = tmp_path / "a.fasta"
	write_fasta(aln, [("REF1", "ACGT")])
	processor = make(tmp_path, padded_aln=str(aln))
	assert processor.split_alignment_into_chunks(str(aln), "REF1", [], ordered=True) == []
	assert processor.split_alignment_into_chunks(str(aln), "REF1", [], ordered=False) == []


def test_chunk_threshold_suppresses_chunking_below_the_trigger(tmp_path: Path, monkeypatch):
	"""--chunk_threshold is documented as 'Trigger iterative chunked placement
	when sequences-to-place exceed this count' and was stored and clamped, but
	run() never read it: placement was always chunked at chunk_size, so the flag
	did nothing and small runs paid for an extra usher round-trip per batch."""
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ACGT")] + [(f"Q{i}", "ACG" + "ACGT"[i % 4]) for i in range(5)])
	output_dir = tmp_path / "out"
	starter = tmp_path / "seed.nwk"
	starter.write_text("(REF1:0.1);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir),
		chunk_size=2, chunk_threshold=10_000,
		starter_tree=str(starter), existing_ids_file=str(existing),
	)
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])

	processor.run()

	# Far below the 10,000 threshold, so one batch is expected.
	assert len(records["place"]) == 1


# ---------------------------------------------------------------------------
# 6. Tree rewriting: depth, degenerate shapes, hostile names
# ---------------------------------------------------------------------------


def test_expand_identical_sequence_tree_outputs_handles_a_deep_ladder_tree(tmp_path: Path):
	"""_replace_terminal_with_polytomy recursed once per tree level and Bio.Phylo's
	newick writer recurses too, so this raised RecursionError on a caterpillar tree
	only a few hundred levels deep - the shape iterative single-sample placement
	produces. It runs AFTER promote_final_usher_outputs, so the run died at the
	very end with all usher work done and every collapsed duplicate lost."""
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	newick = ladder_newick(500)
	for name in ("uncondensed-final-tree.nh", "final-tree.nh"):
		(output_dir / name).write_text(newick + "\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.expand_identical_sequence_tree_outputs({"ANCHOR": ["DUP1", "DUP2"]})

	text = (output_dir / "final-tree.nh").read_text(encoding="utf-8")
	assert "DUP1" in text and "DUP2" in text


def test_expand_identical_sequence_tree_outputs_handles_a_root_only_tree(tmp_path: Path):
	"""The search only ever walked clade.clades, so an anchor that was the root
	itself - a two-sequence segment, or any tree usher emits as a single named
	node - was never found. The duplicates were dropped from the tree behind a
	[warn] and a zero exit status, so the pipeline reported success while
	silently losing samples."""
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	for name in ("uncondensed-final-tree.nh", "final-tree.nh"):
		(output_dir / name).write_text("ANCHOR:0.1;\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.expand_identical_sequence_tree_outputs({"ANCHOR": ["DUP1"]})

	assert "DUP1" in (output_dir / "final-tree.nh").read_text(encoding="utf-8")


def test_expand_identical_sequence_tree_outputs_round_trips_newick_metacharacters(tmp_path: Path):
	"""Accessions containing ',', ':' or brackets must be quoted on the way out,
	or the published tree stops parsing."""
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	for name in ("uncondensed-final-tree.nh", "final-tree.nh"):
		(output_dir / name).write_text("(REF1:0.1,ANCHOR:0.2);\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.expand_identical_sequence_tree_outputs({"ANCHOR": ["A,B", "C:D", "E F"]})

	terminals = UsherPlacement._read_tree_terminals(str(output_dir / "final-tree.nh"))
	assert {"REF1", "ANCHOR", "A,B", "C:D"} <= terminals


def test_expand_identical_sequence_tree_outputs_preserves_the_anchor_branch_length(tmp_path: Path):
	"""The duplicate polytomy must hang off the anchor's original branch, not
	shorten it - otherwise every collapsed swarm shrinks the tree it sits in."""
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	for name in ("uncondensed-final-tree.nh", "final-tree.nh"):
		(output_dir / name).write_text("(REF1:0.1,ANCHOR:0.25);\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.expand_identical_sequence_tree_outputs({"ANCHOR": ["DUP1"]})

	tree = Phylo.read(StringIO((output_dir / "final-tree.nh").read_text(encoding="utf-8")), "newick")
	by_name = {t.name: t for t in tree.get_terminals()}
	assert pytest.approx(0.0) == by_name["ANCHOR"].branch_length
	assert pytest.approx(0.0) == by_name["DUP1"].branch_length
	# Total root-to-tip distance is unchanged by the expansion.
	assert pytest.approx(0.25, abs=1e-6) == tree.distance(tree.root, by_name["ANCHOR"])


def test_expand_identical_sequence_tree_outputs_warns_but_survives_a_missing_anchor(tmp_path: Path, capsys):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	for name in ("uncondensed-final-tree.nh", "final-tree.nh"):
		(output_dir / name).write_text("(REF1:0.1,OTHER:0.2);\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.expand_identical_sequence_tree_outputs({"ABSENT": ["DUP1"]})

	assert "Could not find identical-sequence anchor 'ABSENT'" in capsys.readouterr().out


def test_read_tree_terminals_rejects_a_malformed_newick(tmp_path: Path):
	tree = tmp_path / "broken.nwk"
	tree.write_text("(REF1:0.1,,,;\n", encoding="utf-8")
	with pytest.raises(Exception):
		UsherPlacement._read_tree_terminals(str(tree))


def test_copy_tree_outputs_tolerates_a_starter_tree_already_at_the_destination(tmp_path: Path):
	"""The natural resume invocation reuses the previous run's output directory
	('--starter_tree out/final-tree.nh --output_dir out'), which makes source and
	destination the same path; shutil.copyfile then raised SameFileError after the
	run had already decided there was nothing to place."""
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	starter = output_dir / "final-tree.nh"
	starter.write_text("(REF1:0.1,Q1:0.2);\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.copy_tree_outputs(str(starter))

	assert (output_dir / "uncondensed-final-tree.nh").exists()


# ---------------------------------------------------------------------------
# 7. External-tool boundary
# ---------------------------------------------------------------------------


def test_run_command_reports_a_missing_binary_instead_of_raising(tmp_path: Path):
	"""run_command caught only CalledProcessError. When faToVcf is not on PATH -
	a wrong conda env, a partially built container - subprocess.run raises
	FileNotFoundError, which escaped as a bare traceback instead of the handled
	failure path the caller was written for."""
	processor = make(tmp_path)
	assert processor.run_command(["definitely-not-a-real-binary-xyz"], "[warn] missing") is False


def test_append_threads_reports_a_missing_usher_binary_clearly(tmp_path: Path, monkeypatch):
	"""check=False only suppresses a non-zero exit; a missing usher binary still
	raises FileNotFoundError from the probe itself, which surfaced as an unhandled
	exception in the middle of thread detection."""
	processor = make(tmp_path, threads=4)

	def fake_run(cmd, **kwargs):
		raise FileNotFoundError(2, "No such file or directory: 'usher'")
	monkeypatch.setattr("subprocess.run", fake_run)

	with pytest.raises(RuntimeError, match="usher"):
		processor._append_threads(["usher", "-i", "x"])


def test_append_threads_does_not_duplicate_the_flag_across_commands(tmp_path: Path, monkeypatch):
	"""The probe is cached; each command must still get exactly one -T pair."""
	processor = make(tmp_path, threads=6)

	class FakeCompleted:
		stdout = "usage: usher -T [ --threads ] arg"
		stderr = ""

	monkeypatch.setattr("subprocess.run", lambda cmd, **kwargs: FakeCompleted())

	cmd = processor._append_threads(["usher", "-i", "x"])
	cmd = processor._append_threads(cmd)
	assert cmd.count("-T") == 2, "re-appending is the caller's bug, not the method's"

	fresh = processor._append_threads(["usher", "-i", "y"])
	assert fresh.count("-T") == 1


def test_append_threads_detects_support_from_stderr_only_help(tmp_path: Path, monkeypatch):
	"""Some usher builds print usage to stderr; the probe reads both streams."""
	processor = make(tmp_path, threads=2)

	class FakeCompleted:
		stdout = ""
		stderr = "Options:\n  -T [ --threads ] arg   number of threads\n"

	monkeypatch.setattr("subprocess.run", lambda cmd, **kwargs: FakeCompleted())
	assert processor._append_threads(["usher"])[-2:] == ["-T", "2"]


def test_build_vcf_does_not_silently_place_excluded_ids_on_retry(tmp_path: Path, monkeypatch):
	"""The fallback silently dropped -excludeFile when the first faToVcf call
	failed. -excludeFile is what keeps backbone tips out of the placement VCF, so
	the retry re-submitted every existing tree tip for placement, duplicating
	them. A failure to exclude must abort, not proceed unfiltered."""
	aln = tmp_path / "aln.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("BACKBONE", "ACGA"), ("Q1", "ACGG")])
	exclude = tmp_path / "exclude.txt"
	exclude.write_text("BACKBONE\n", encoding="utf-8")
	processor = make(tmp_path, padded_aln=str(aln))

	def fake_run(cmd, check=True):
		if any(str(a).startswith("-excludeFile=") for a in cmd):
			raise subprocess.CalledProcessError(1, cmd)
		return None
	monkeypatch.setattr("subprocess.run", fake_run)

	with pytest.raises(Exception):
		processor.build_vcf(
			"REF1", exclude_ids_file=str(exclude),
			alignment_fasta=str(aln), vcf_path=str(tmp_path / "out.vcf"),
		)


def test_build_vcf_refuses_option_like_identifiers(tmp_path: Path, monkeypatch):
	"""faToVcf takes its alignment and output positionally with no '--' separator,
	so a value beginning with '-' is parsed as an option. ref_id can fall back to
	whatever the first alignment record happens to be, and '-' is legal in a FASTA
	header, so this is reachable from input data."""
	aln = tmp_path / "aln.fasta"
	write_fasta(aln, [("-startTree", "ACGT"), ("Q1", "ACGA")])
	processor = make(tmp_path, padded_aln=str(aln))
	monkeypatch.setattr("subprocess.run", lambda cmd, check=True: None)

	with pytest.raises(ValueError):
		processor.build_vcf(
			"-startTree", alignment_fasta=str(aln), vcf_path=str(tmp_path / "out.vcf"),
		)


def test_run_usher_cmd_truncates_a_stale_log_rather_than_appending(tmp_path: Path, monkeypatch):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	(output_dir / "usher.verbose.log").write_text("PREVIOUS RUN\n", encoding="utf-8")
	processor = make(tmp_path, output_dir=str(output_dir), threads=2)

	def fake_run(cmd, **kwargs):
		kwargs["stdout"].write("fresh\n")
	monkeypatch.setattr("subprocess.run", fake_run)

	processor._run_usher_cmd(["usher"], str(output_dir))
	text = (output_dir / "usher.verbose.log").read_text(encoding="utf-8")
	assert text == "fresh\n" and "PREVIOUS RUN" not in text


def test_threads_are_capped_at_the_available_cpu_count(tmp_path: Path):
	"""threads was only clamped from below, so --threads 512 went straight to
	usher's -T and into OMP_NUM_THREADS, oversubscribing a shared node."""
	processor = make(tmp_path, threads=100_000)
	assert processor.threads <= (os.cpu_count() or 1)


def test_run_usher_cmd_pins_every_blas_thread_variable(tmp_path: Path, monkeypatch):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	processor = make(tmp_path, output_dir=str(output_dir), threads=3)
	monkeypatch.setenv("OMP_NUM_THREADS", "64")

	captured = {}

	def fake_run(cmd, **kwargs):
		captured.update(kwargs)
		kwargs["stdout"].write("")
	monkeypatch.setattr("subprocess.run", fake_run)

	processor._run_usher_cmd(["usher"], str(output_dir))

	# A hostile ambient value must be overridden, not inherited.
	for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
		assert captured["env"][var] == "3"


# ---------------------------------------------------------------------------
# 8. Mode resolution and asset selection
# ---------------------------------------------------------------------------


def test_resume_mode_takes_precedence_over_update_db(tmp_path: Path, monkeypatch):
	"""Passing both --starter_tree and --update_db is ambiguous; resume wins today
	and prepare_update_assets is never consulted. Pin it so the precedence cannot
	drift silently."""
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("Q1", "ACGA")])
	output_dir = tmp_path / "out"
	starter = tmp_path / "seed.nwk"
	starter.write_text("(REF1:0.1);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir),
		starter_tree=str(starter), existing_ids_file=str(existing),
		update_db=str(tmp_path / "never_opened.db"),
	)

	def boom():
		raise AssertionError("prepare_update_assets must not run in resume mode")
	monkeypatch.setattr(processor, "prepare_update_assets", boom)
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])

	processor.run()
	assert records["backbone"][0]["tree_file"] == str(starter)


def test_resolve_resume_assets_rejects_a_starter_tree_that_does_not_exist(tmp_path: Path):
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")
	processor = make(
		tmp_path,
		starter_tree=str(tmp_path / "absent.nwk"),
		existing_ids_file=str(existing),
	)
	with pytest.raises(FileNotFoundError, match="Resume starter tree not found"):
		processor.resolve_resume_assets()


def test_resolve_non_update_assets_falls_back_to_any_fasta_with_a_warning(tmp_path: Path, capsys):
	mmseq = tmp_path / "mmseq"
	mmseq.mkdir()
	(mmseq / "clusters.fasta").write_text(">REF1\nACGT\n", encoding="utf-8")
	iqtree = tmp_path / "iqtree"
	iqtree.mkdir()
	(iqtree / "run.treefile").write_text("(REF1:0.1);\n", encoding="utf-8")

	processor = make(tmp_path, mmseq_cluster_dir=str(mmseq), iqtree_dir=str(iqtree))
	cluster_rep, tree_file = processor.resolve_non_update_assets()

	assert cluster_rep.endswith("clusters.fasta")
	assert tree_file.endswith("run.treefile")
	assert "*_cluster_rep.fasta not found" in capsys.readouterr().out


def test_find_first_prefers_an_exact_name_over_a_suffix_match(tmp_path: Path):
	"""_find_first walked files in os.walk order and tested exact_name and
	pattern_suffix inside the same per-file loop, so whichever file the filesystem
	yielded first won even when a later file was the requested exact name."""
	root = tmp_path / "tree"
	(root / "sub").mkdir(parents=True)
	(root / "sub" / "other_cluster_rep.fasta").write_text(">A\nA\n", encoding="utf-8")
	(root / "sub" / "wanted.fasta").write_text(">B\nA\n", encoding="utf-8")

	found = UsherPlacement._find_first(str(root), pattern_suffix=".fasta", exact_name="wanted.fasta")
	assert Path(found).name == "wanted.fasta"


def test_collect_backbone_fasta_warns_when_a_backbone_tip_has_no_sequence(tmp_path: Path, capsys):
	"""A tip with no genotype leaves its branch mutation-less - the exact failure
	build_backbone_pb exists to prevent - so the warning must name it."""
	aln = tmp_path / "aln.fasta"
	write_fasta(aln, [("REF1", "ACGT")])
	processor = make(tmp_path, padded_aln=str(aln))

	processor.collect_backbone_fasta(["GHOST_TIP"], "REF1", str(aln))

	out = capsys.readouterr().out
	assert "1 backbone tip(s) have no aligned sequence" in out
	assert "GHOST_TIP" in out


def test_collect_backbone_fasta_prefers_the_working_alignment_over_the_db(tmp_path: Path):
	"""The working alignment is the current truth; a stale DB row must not win."""
	db = tmp_path / "update.db"
	conn = sqlite3.connect(str(db))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment TEXT)")
		cur.execute("INSERT INTO sequence_alignment VALUES ('OLD1', 'TTTT')")
		conn.commit()
	finally:
		conn.close()

	aln = tmp_path / "aln.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("OLD1", "ACGA")])
	processor = make(tmp_path, padded_aln=str(aln), update_db=str(db))

	out_path = processor.collect_backbone_fasta(["OLD1"], "REF1", str(aln))
	text = Path(out_path).read_text(encoding="utf-8")
	assert "ACGA" in text and "TTTT" not in text


def test_run_raises_before_touching_usher_when_the_alignment_is_empty(tmp_path: Path, monkeypatch):
	aln = tmp_path / "segment_1_dedup.fasta"
	aln.write_text("", encoding="utf-8")
	output_dir = tmp_path / "out"
	starter = tmp_path / "seed.nwk"
	starter.write_text("(REF1:0.1);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir),
		starter_tree=str(starter), existing_ids_file=str(existing),
	)
	monkeypatch.setattr(processor, "build_backbone_pb", lambda *a, **k: pytest.fail("usher reached"))

	with pytest.raises(ValueError):
		processor.run()


# ---------------------------------------------------------------------------
# 9. Progress reporting
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
	"completed,total,expected",
	[
		(0, 0, "[" + "-" * 30 + "]"),
		(1, 0, "[" + "-" * 30 + "]"),
		(0, 4, "[" + "-" * 30 + "]"),
		(4, 4, "[" + "#" * 30 + "]"),
		(8, 4, "[" + "#" * 30 + "]"),
		(-1, 4, "[" + "-" * 30 + "]"),
	],
)
def test_render_progress_bar_clamps_degenerate_counts(completed, total, expected):
	assert UsherPlacement._render_progress_bar(completed, total) == expected


def test_log_chunk_progress_never_reports_negative_remaining(capsys):
	UsherPlacement.log_chunk_progress(
		UsherPlacement(padded_aln="a", output_dir="b"), completed_chunks=9, total_chunks=4
	)
	assert "remaining: 0." in capsys.readouterr().out


# ---------------------------------------------------------------------------
# 10. Regression cover for the fixes themselves
#
# The sections above describe the defects. These assert the specific mechanism
# each fix relies on, so a later refactor that reintroduces one fails here with
# a pointed message rather than only through an end-to-end test.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
	"filename,expected",
	[
		# Real names this pipeline produces. Segmented builds carry refset_<N>.
		("refset_6_aln_merged_MSA_dedup.fasta", "6"),
		("refset_11_aln_merged_MSA_dedup.fasta", "11"),
		("refset_0_aln_merged_MSA_dedup.fasta", "0"),
		# Non-segmented builds are named for their reference accession. Those
		# digits are an accession, not a segment, so they name no segment - which
		# is what prepare_update_assets' COALESCE(segment, '0') stores for them.
		("NC_001542.aligned_merged_MSA_dedup.fasta", "0"),
		("AF009606.aligned_merged_MSA_dedup.fasta", "0"),
		("CY087822.aligned_merged_MSA_dedup.fasta", "0"),
		# Subtypes and years in the name must not win over the segment field.
		("h3n2_segment_4_dedup.fasta", "4"),
		("IAV_H1N1_segment_2_dedup.fasta", "2"),
		("rabv_2024_segment_1_dedup.fasta", "1"),
		("seg7.fasta", "7"),
		("RNA_4_dedup.fasta", "4"),
	],
)
def test_segment_resolution_over_real_pipeline_filenames(tmp_path: Path, filename, expected):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / filename), output_dir=str(tmp_path / "out"))
	assert processor._segment_from_padded_alignment() == expected


@pytest.mark.parametrize("supplied,expected", [("4", "4"), ("04", "4"), ("segment 7", "7"), (" 2 ", "2")])
def test_explicit_segment_overrides_the_filename(tmp_path: Path, supplied, expected):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "h3n2_segment_9_dedup.fasta"),
		output_dir=str(tmp_path / "out"),
		segment=supplied,
	)
	assert processor._segment_from_padded_alignment() == expected


@pytest.mark.parametrize("sentinel", [None, "", "   ", "null", "UNSET"])
def test_segment_sentinels_fall_back_to_the_filename(tmp_path: Path, sentinel):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "segment_3_dedup.fasta"),
		output_dir=str(tmp_path / "out"),
		segment=sentinel,
	)
	assert processor._segment_from_padded_alignment() == "3"


def test_a_filename_naming_two_different_segments_is_refused(tmp_path: Path):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "refset_6_aln_segment_2_dedup.fasta"),
		output_dir=str(tmp_path / "out"),
	)
	with pytest.raises(ValueError, match="more than one segment"):
		processor._segment_from_padded_alignment()


def test_a_filename_repeating_one_segment_is_not_ambiguous(tmp_path: Path):
	processor = UsherPlacement(
		padded_aln=str(tmp_path / "segment_4_refset_4_dedup.fasta"),
		output_dir=str(tmp_path / "out"),
	)
	assert processor._segment_from_padded_alignment() == "4"


def test_update_mode_finds_the_tree_for_a_non_segmented_build(tmp_path: Path):
	"""End-to-end consequence of the segment fix: a non-segmented build stores
	segment NULL, so the alignment side has to resolve to '0' for the two to meet.
	Under the old first-digit-run parse this asked for segment '009606' and update
	mode could never start."""
	db = tmp_path / "update.db"
	conn = sqlite3.connect(str(db))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE trees (source TEXT, segment TEXT, newick TEXT)")
		cur.execute("INSERT INTO trees VALUES ('usher', NULL, '(AF009606:0.1,OLD1:0.2);')")
		conn.commit()
	finally:
		conn.close()

	aln = tmp_path / "AF009606.aligned_merged_MSA_dedup.fasta"
	write_fasta(aln, [("AF009606", "ACGT"), ("Q1", "ACGA")])
	processor = make(tmp_path, padded_aln=str(aln), update_db=str(db))

	tree_path, ids_path = processor.prepare_update_assets()

	assert "OLD1" in Path(tree_path).read_text(encoding="utf-8")
	assert Path(ids_path).read_text(encoding="utf-8").split() == ["AF009606", "OLD1"]
	assert Path(tree_path).name == "tree_segment_0.nwk"


@pytest.mark.parametrize("depth", [500, 2000, 5000])
def test_expansion_survives_very_deep_ladder_trees(tmp_path: Path, depth):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	for name in ("uncondensed-final-tree.nh", "final-tree.nh"):
		(output_dir / name).write_text(ladder_newick(depth) + "\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.expand_identical_sequence_tree_outputs({"ANCHOR": ["DUP1", "DUP2"]})

	terminals = UsherPlacement._read_tree_terminals(str(output_dir / "final-tree.nh"))
	assert {"ANCHOR", "DUP1", "DUP2"} <= terminals
	assert len(terminals) == depth + 3


def test_every_anchor_is_expanded_in_a_single_pass(tmp_path: Path):
	"""One walk handles all anchors; the old code re-searched the whole tree per
	anchor, which was O(anchors x tree) on top of being recursive."""
	anchors = [f"A{i:03d}" for i in range(200)]
	newick = "(" + ",".join(f"{a}:0.1" for a in anchors) + ");"
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	for name in ("uncondensed-final-tree.nh", "final-tree.nh"):
		(output_dir / name).write_text(newick + "\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.expand_identical_sequence_tree_outputs({a: [f"{a}_dup"] for a in anchors})

	terminals = UsherPlacement._read_tree_terminals(str(output_dir / "final-tree.nh"))
	assert len(terminals) == 400
	for anchor in anchors:
		assert f"{anchor}_dup" in terminals


def test_root_only_tree_becomes_a_polytomy_keeping_its_branch_length(tmp_path: Path):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	for name in ("uncondensed-final-tree.nh", "final-tree.nh"):
		(output_dir / name).write_text("ANCHOR:0.3;\n", encoding="utf-8")

	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fasta"), output_dir=str(output_dir))
	processor.expand_identical_sequence_tree_outputs({"ANCHOR": ["DUP1", "DUP2"]})

	tree = Phylo.read(StringIO((output_dir / "final-tree.nh").read_text(encoding="utf-8")), "newick")
	assert {t.name for t in tree.get_terminals()} == {"ANCHOR", "DUP1", "DUP2"}
	for terminal in tree.get_terminals():
		assert pytest.approx(0.0) == terminal.branch_length


@pytest.mark.parametrize(
	"newick,expected",
	[
		("(A:0.1,B:0.2);", 1),
		("((A:0.1,B:0.2):0.1,C:0.3);", 2),
		("A:0.1;", 0),
		# Brackets inside a quoted label are label text, not nesting.
		("('weird(name)here':0.1,B:0.2);", 1),
		('("also,(quoted)":0.1,B:0.2);', 1),
	],
)
def test_estimate_newick_depth_ignores_brackets_inside_quoted_labels(newick, expected):
	assert UsherPlacement._estimate_newick_depth(newick) == expected


def test_read_tree_terminals_handles_a_deep_tree(tmp_path: Path):
	tree = tmp_path / "deep.nwk"
	tree.write_text(ladder_newick(3000) + "\n", encoding="utf-8")
	assert len(UsherPlacement._read_tree_terminals(str(tree))) == 3001


def test_run_uses_one_batch_below_the_chunk_threshold(tmp_path: Path, monkeypatch):
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ACGT")] + [(f"Q{i}", "ACG" + "ACGT"[i % 4]) for i in range(5)])
	output_dir = tmp_path / "out"
	starter = tmp_path / "seed.nwk"
	starter.write_text("(REF1:0.1);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir),
		chunk_size=2, chunk_threshold=10_000,
		starter_tree=str(starter), existing_ids_file=str(existing),
	)
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])

	processor.run()

	assert len(records["place"]) == 1
	# The order report must describe the batching that actually happened.
	rows = (output_dir / "placement_order.tsv").read_text(encoding="utf-8").strip().splitlines()[1:]
	assert {row.split("\t")[3] for row in rows} == {"1"}


def test_run_chunks_once_the_threshold_is_exceeded(tmp_path: Path, monkeypatch):
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ACGT")] + [(f"Q{i}", "ACG" + "ACGT"[i % 4]) for i in range(5)])
	output_dir = tmp_path / "out"
	starter = tmp_path / "seed.nwk"
	starter.write_text("(REF1:0.1);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir),
		chunk_size=2, chunk_threshold=2,
		starter_tree=str(starter), existing_ids_file=str(existing),
	)
	records = install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])

	processor.run()

	# Q0/Q4 share a genotype and Q3 matches the reference, leaving 3 placeable
	# anchors; 3 > threshold 2, so they are split into batches of chunk_size 2.
	assert len(records["place"]) == 2
	assert [len(call["ids"]) - 1 for call in records["vcf"]] == [2, 1]


def test_splitter_always_honours_an_explicit_chunk_size(tmp_path: Path):
	"""The threshold is run()-level policy. Asking the splitter for batches of N
	must still give batches of N regardless of chunk_threshold."""
	aln = tmp_path / "a.fasta"
	write_fasta(aln, [("REF1", "ACGT")] + [(f"Q{i}", "ACG" + "ACGT"[i % 4]) for i in range(4)])
	processor = make(tmp_path, padded_aln=str(aln), chunk_size=2, chunk_threshold=1_000_000)

	chunks = processor._split_alignment_into_chunks_python(
		str(aln), "REF1", [f"Q{i}" for i in range(4)])
	assert len(chunks) == 2


def test_prepare_chunk_dir_removes_everything_from_a_previous_run(tmp_path: Path):
	processor = make(tmp_path)
	chunk_dir = Path(processor.output_dir) / "chunks"
	(chunk_dir / "raw").mkdir(parents=True, exist_ok=True)
	(chunk_dir / "chunk_0009.fasta").write_text(">STALE\nAAAA\n", encoding="utf-8")
	(chunk_dir / "raw" / "placeable_only.part_001.fasta").write_text(">STALE\nAAAA\n", encoding="utf-8")

	processor._prepare_chunk_dir()

	assert list(chunk_dir.iterdir()) == []


@pytest.mark.parametrize(
	"names,expected_order",
	[
		(["p.part_002.fa", "p.part_001.fa", "p.part_010.fa"], [1, 2, 10]),
		(["p.part_999.fa", "p.part_1000.fa", "p.part_1001.fa", "p.part_100.fa"], [100, 999, 1000, 1001]),
	],
)
def test_natural_sort_key_orders_parts_numerically(names, expected_order):
	ordered = sorted(names, key=UsherPlacement._natural_sort_key)
	numbers = [int(re.search(r"part_(\d+)", n).group(1)) for n in ordered]
	assert numbers == expected_order


def test_present_placeable_ids_drops_absent_and_repeated_ids(tmp_path: Path, capsys):
	aln = tmp_path / "a.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("Q1", "ACGA"), ("Q2", "ACGG")])
	processor = make(tmp_path, padded_aln=str(aln))

	assert processor._present_placeable_ids(str(aln), ["Q1", "GHOST", "Q1", "Q2"]) == ["Q1", "Q2"]
	assert "GHOST" in capsys.readouterr().out


def test_collapse_reports_duplicates_and_uninformative_sequences(tmp_path: Path):
	aln = tmp_path / "messy.fasta"
	with aln.open("w", encoding="utf-8") as handle:
		handle.write(">REF1\nACGT\n>Q1\nACGA\n>Q1\nTTTT\n>ALLN\nNNNN\n>GAPS\n----\n")
	processor = make(tmp_path, padded_aln=str(aln))

	plan = processor.collapse_identical_sequences(str(aln), "REF1", set())

	assert plan["placeable_ids"] == ["Q1"]
	assert plan["duplicate_ids"] == ["Q1"]
	assert plan["uninformative_ids"] == ["ALLN", "GAPS"]

	report = Path(processor.write_excluded_sequence_report(plan))
	assert report.read_text(encoding="utf-8").strip().splitlines() == [
		"sequence_id\treason",
		"Q1\tduplicate_accession_in_alignment",
		"ALLN\tfewer_than_min_informative_bases",
		"GAPS\tfewer_than_min_informative_bases",
	]


def test_min_informative_bases_zero_restores_the_old_permissive_behaviour(tmp_path: Path):
	aln = tmp_path / "alln.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("ALLN", "NNNN")])
	processor = make(tmp_path, padded_aln=str(aln), min_informative_bases=0)

	plan = processor.collapse_identical_sequences(str(aln), "REF1", set())
	assert plan["placeable_ids"] == ["ALLN"]
	assert plan["uninformative_ids"] == []


def test_min_informative_bases_can_demand_more_coverage(tmp_path: Path):
	aln = tmp_path / "partial.fasta"
	write_fasta(aln, [("REF1", "ACGTACGTAC"), ("THIN", "AC" + "N" * 8), ("FULL", "ACGTACGTAA")])
	processor = make(tmp_path, padded_aln=str(aln), min_informative_bases=5)

	plan = processor.collapse_identical_sequences(str(aln), "REF1", set())
	assert plan["placeable_ids"] == ["FULL"]
	assert plan["uninformative_ids"] == ["THIN"]


def test_an_existing_tree_tip_is_never_dropped_for_being_uninformative(tmp_path: Path):
	"""Dropping a backbone tip would delete it from the tree it is already in."""
	aln = tmp_path / "a.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("OLD_THIN", "NNNN"), ("Q1", "ACGA")])
	processor = make(tmp_path, padded_aln=str(aln), min_informative_bases=2)

	plan = processor.collapse_identical_sequences(str(aln), "REF1", {"OLD_THIN"})
	assert plan["uninformative_ids"] == []
	assert plan["placeable_ids"] == ["Q1"]


def test_collapse_groups_mixed_case_against_an_existing_tree_anchor(tmp_path: Path):
	aln = tmp_path / "case.fasta"
	write_fasta(aln, [("REF1", "AAAA"), ("OLD1", "CCCC"), ("NEW1", "cccc"), ("NEW2", "CcCc")])
	processor = make(tmp_path, padded_aln=str(aln))

	plan = processor.collapse_identical_sequences(str(aln), "REF1", {"OLD1"})

	assert plan["placeable_ids"] == []
	assert sorted(plan["anchor_to_members"]["OLD1"]) == ["NEW1", "NEW2"]


def test_hash_grouping_does_not_merge_distinct_sequences(tmp_path: Path):
	"""Bucketing on a digest must not collapse genuinely different genotypes."""
	records = [("REF1", "ACGT" * 10)]
	records += [(f"Q{i:03d}", ("ACGT" * 10)[:39] + "ACGT"[i % 4] if i % 4 else ("ACGT" * 10)[:38] + "AA")
				for i in range(200)]
	aln = tmp_path / "many.fasta"
	write_fasta(aln, records)
	processor = make(tmp_path, padded_aln=str(aln))

	plan = processor.collapse_identical_sequences(str(aln), "REF1", set())

	distinct = {str(r[1]).upper() for r in records[1:]}
	assert len(plan["placeable_ids"]) == len(distinct - {("ACGT" * 10).upper()})


def test_run_always_writes_the_excluded_sequences_report(tmp_path: Path, monkeypatch):
	"""Even when nothing was excluded - an absent file is not evidence."""
	aln = tmp_path / "segment_1_dedup.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("Q1", "ACGA")])
	output_dir = tmp_path / "out"
	starter = tmp_path / "seed.nwk"
	starter.write_text("(REF1:0.1);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir),
		starter_tree=str(starter), existing_ids_file=str(existing),
	)
	install_run_fakes(processor, monkeypatch, output_dir, "REF1", ["REF1"])
	processor.run()

	assert (output_dir / "excluded_sequences.tsv").read_text(encoding="utf-8").strip() == "sequence_id\treason"


@pytest.mark.parametrize(
	"kwargs",
	[
		{"ref_id": "-startTree"},
		{"alignment_fasta": "-oops.fasta"},
		{"vcf_path": "-out.vcf"},
		{"exclude_ids_file": "-ids.txt"},
	],
)
def test_build_vcf_rejects_every_option_like_argument(tmp_path: Path, monkeypatch, kwargs):
	aln = tmp_path / "aln.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("Q1", "ACGA")])
	processor = make(tmp_path, padded_aln=str(aln))
	monkeypatch.setattr("subprocess.run", lambda cmd, check=True: None)

	call = {
		"ref_id": "REF1",
		"alignment_fasta": str(aln),
		"vcf_path": str(tmp_path / "out.vcf"),
	}
	call.update(kwargs)
	ref_id = call.pop("ref_id")
	with pytest.raises(ValueError, match="command-line"):
		processor.build_vcf(ref_id, **call)


def test_build_vcf_raises_when_faToVcf_is_absent(tmp_path: Path, monkeypatch):
	aln = tmp_path / "aln.fasta"
	write_fasta(aln, [("REF1", "ACGT")])
	processor = make(tmp_path, padded_aln=str(aln))

	def missing(cmd, check=True):
		raise FileNotFoundError(2, "No such file or directory: 'faToVcf'")
	monkeypatch.setattr("subprocess.run", missing)

	with pytest.raises(RuntimeError, match="faToVcf failed"):
		processor.build_vcf("REF1", alignment_fasta=str(aln), vcf_path=str(tmp_path / "o.vcf"))


def test_run_usher_cmd_raises_a_clear_error_when_usher_is_absent(tmp_path: Path, monkeypatch):
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	processor = make(tmp_path, output_dir=str(output_dir))

	def missing(cmd, **kwargs):
		raise FileNotFoundError(2, "usher")
	monkeypatch.setattr("subprocess.run", missing)

	with pytest.raises(RuntimeError, match="not on PATH"):
		processor._run_usher_cmd(["usher", "-i", "x"], str(output_dir))


def test_place_onto_pb_does_not_accept_a_previous_runs_outputs(tmp_path: Path, monkeypatch):
	"""A usher invocation that silently produces nothing must not pass the output
	checks on files left behind by an earlier run into the same directory."""
	output_dir = tmp_path / "out"
	chunk_dir = output_dir / "chunk_0001"
	chunk_dir.mkdir(parents=True, exist_ok=True)
	(chunk_dir / "final-tree.nh").write_text("(STALE:0.1);\n", encoding="utf-8")
	(chunk_dir / "uncondensed-final-tree.nh").write_text("(STALE:0.1);\n", encoding="utf-8")
	(chunk_dir / "usher.pb").write_text("stale pb", encoding="utf-8")

	processor = make(tmp_path, output_dir=str(output_dir))
	monkeypatch.setattr(processor, "_append_threads", lambda cmd: cmd)
	monkeypatch.setattr(processor, "_run_usher_cmd", lambda *a, **k: None)

	vcf = tmp_path / "c.vcf"
	vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
	with pytest.raises(FileNotFoundError):
		processor.place_onto_pb(str(tmp_path / "in.pb"), str(vcf), str(chunk_dir))


@pytest.mark.parametrize("requested", [1, 2, 10**6])
def test_threads_stay_within_the_machine(tmp_path: Path, requested):
	processor = make(tmp_path, threads=requested)
	assert 1 <= processor.threads <= (os.cpu_count() or 1)
	assert processor.threads <= requested


def test_find_first_is_stable_across_calls(tmp_path: Path):
	"""os.walk yields entries in inode order on ext4, so without the sort which
	cluster representative gets picked is not reproducible across machines."""
	root = tmp_path / "many"
	for sub in ("c", "a", "b"):
		(root / sub).mkdir(parents=True)
		(root / sub / f"{sub}_cluster_rep.fasta").write_text(">X\nA\n", encoding="utf-8")

	answers = {UsherPlacement._find_first(str(root), pattern_suffix="_cluster_rep.fasta")
			   for _ in range(5)}
	assert len(answers) == 1
	assert Path(answers.pop()).name == "a_cluster_rep.fasta"


def test_find_first_exact_name_wins_over_a_suffix_match_listed_earlier(tmp_path: Path):
	root = tmp_path / "tree"
	(root / "sub").mkdir(parents=True)
	(root / "sub" / "aaa_cluster_rep.fasta").write_text(">A\nA\n", encoding="utf-8")
	(root / "sub" / "zzz_wanted.fasta").write_text(">B\nA\n", encoding="utf-8")

	found = UsherPlacement._find_first(
		str(root), pattern_suffix=".fasta", exact_name="zzz_wanted.fasta")
	assert Path(found).name == "zzz_wanted.fasta"


@pytest.mark.parametrize(
	"stored,expected",
	[
		("ACGT", "ACGT"),
		("  ACGT  ", "ACGT"),
		("ACG\nT", "ACGT"),
		("AC GT\r\n", "ACGT"),
	],
)
def test_sanitise_db_sequence_strips_whitespace(stored, expected):
	assert UsherPlacement._sanitise_db_sequence("ACC1", stored) == expected


def test_sanitise_db_sequence_refuses_a_value_containing_a_header():
	with pytest.raises(ValueError, match="not a sequence"):
		UsherPlacement._sanitise_db_sequence("ACC1", "ACGT\n>EVIL\nTTTT")


def test_collect_backbone_fasta_sanitises_db_rows(tmp_path: Path):
	db = tmp_path / "update.db"
	conn = sqlite3.connect(str(db))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment TEXT)")
		cur.execute("INSERT INTO sequence_alignment VALUES ('OLD1', 'AC GT\n')")
		conn.commit()
	finally:
		conn.close()

	aln = tmp_path / "aln.fasta"
	write_fasta(aln, [("REF1", "ACGT")])
	processor = make(tmp_path, padded_aln=str(aln), update_db=str(db))

	out_path = processor.collect_backbone_fasta(["OLD1"], "REF1", str(aln))
	assert UsherPlacement._read_ids_from_fasta(out_path) == ["REF1", "OLD1"]
	assert ">OLD1\nACGT\n" in Path(out_path).read_text(encoding="utf-8")


def test_a_non_utf8_alignment_survives_the_whole_collapse_path(tmp_path: Path):
	aln = tmp_path / "latin1.fasta"
	aln.write_bytes(b">REF1\nACGT\n>Q\xe9_ISOLATE\nACGA\n>Q2\nACGG\n")
	processor = make(tmp_path, padded_aln=str(aln))

	plan = processor.collapse_identical_sequences(str(aln), "REF1", set())

	assert len(plan["placeable_ids"]) == 2
	assert "Q2" in plan["placeable_ids"]
	# The bad byte is replaced, not fatal, and the id still round-trips.
	assert UsherPlacement._read_ids_from_fasta(str(aln))[1] in plan["placeable_ids"]


def test_a_non_utf8_tree_file_can_still_be_read(tmp_path: Path):
	tree = tmp_path / "latin1.nwk"
	tree.write_bytes(b"(REF1:0.1,Q\xe9:0.2);\n")
	assert len(UsherPlacement._read_tree_terminals(str(tree))) == 2


def test_resume_into_the_previous_output_directory_reuses_the_tree(tmp_path: Path, monkeypatch):
	"""The full SameFileError scenario: resume from out/final-tree.nh into out/,
	with nothing new to place."""
	output_dir = tmp_path / "out"
	output_dir.mkdir(parents=True, exist_ok=True)
	starter = output_dir / "final-tree.nh"
	starter.write_text("(REF1:0.1,PLACED1:0.2);\n", encoding="utf-8")
	(output_dir / "uncondensed-final-tree.nh").write_text(
		"(REF1:0.1,PLACED1:0.2);\n", encoding="utf-8")
	existing = tmp_path / "existing.txt"
	existing.write_text("REF1\nPLACED1\n", encoding="utf-8")
	aln = tmp_path / "aln.fasta"
	write_fasta(aln, [("REF1", "ACGT"), ("PLACED1", "ACGA"), ("DUP1", "ACGA")])

	processor = UsherPlacement(
		padded_aln=str(aln), output_dir=str(output_dir),
		starter_tree=str(starter), existing_ids_file=str(existing),
	)
	monkeypatch.setattr(processor, "build_backbone_pb",
						lambda *a, **k: pytest.fail("nothing needed placing"))

	processor.run()

	terminals = UsherPlacement._read_tree_terminals(str(output_dir / "final-tree.nh"))
	assert terminals == {"REF1", "PLACED1", "DUP1"}
