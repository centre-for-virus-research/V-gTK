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


def test_build_vcf_retries_without_exclude_when_exclude_fails(tmp_path: Path, monkeypatch):
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
	result = processor.build_vcf(
		"REF1",
		exclude_ids_file=str(exclude),
		alignment_fasta=str(msa),
		vcf_path=str(vcf_path),
	)

	assert result == str(vcf_path)
	# First attempt used -excludeFile and failed; the retry dropped it.
	assert len(cmds) == 2
	assert any(str(arg).startswith("-excludeFile=") for arg in cmds[0])
	assert not any(str(arg).startswith("-excludeFile=") for arg in cmds[1])


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
