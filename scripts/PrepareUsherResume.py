#!/usr/bin/env python3
import argparse
import os
import re
from io import StringIO
from pathlib import Path

from Bio import Phylo, SeqIO


def _chunk_index(path_value):
	match = re.search(r"chunk_(\d+)", str(path_value))
	if not match:
		return None
	return int(match.group(1))


def _tree_candidates(chunk_dir):
	for name in ["uncondensed-final-tree.nh", "final-tree.nh"]:
		path = Path(chunk_dir) / name
		if path.is_file():
			yield path


def find_latest_chunk_with_tree(run_dir):
	run_path = Path(run_dir)
	chunk_dirs = sorted(
		[path for path in run_path.iterdir() if path.is_dir() and _chunk_index(path) is not None],
		key=lambda path: _chunk_index(path),
	)
	for chunk_dir in reversed(chunk_dirs):
		for tree_path in _tree_candidates(chunk_dir):
			return chunk_dir, tree_path
	raise FileNotFoundError(f"No completed chunk tree found in {run_dir}")


def normalise_tree(tree_path, output_path):
	with open(tree_path, "r", encoding="utf-8") as handle:
		newick = handle.read().strip()
	if not newick:
		raise ValueError(f"Tree file is empty: {tree_path}")
	tree = Phylo.read(StringIO(newick), "newick")
	with open(output_path, "w", encoding="utf-8") as handle:
		Phylo.write(tree, handle, "newick")
	return output_path


def read_tree_terminals(tree_path):
	with open(tree_path, "r", encoding="utf-8") as handle:
		newick = handle.read().strip()
	if not newick:
		return []
	tree = Phylo.read(StringIO(newick), "newick")
	return sorted(term.name for term in tree.get_terminals() if term.name)


def read_ids_from_fasta(fasta_path):
	ids = []
	with open(fasta_path, "r", encoding="utf-8") as handle:
		for line in handle:
			if line.startswith(">"):
				ids.append(line[1:].strip().split()[0])
	return ids


def find_remaining_chunk_fastas(run_dir, latest_chunk_index):
	chunks_dir = Path(run_dir) / "chunks"
	if not chunks_dir.is_dir():
		raise FileNotFoundError(f"Chunk FASTA directory not found: {chunks_dir}")
	chunk_fastas = sorted(
		[path for path in chunks_dir.iterdir() if path.is_file() and _chunk_index(path) is not None],
		key=lambda path: _chunk_index(path),
	)
	remaining = [path for path in chunk_fastas if _chunk_index(path) > latest_chunk_index]
	if not remaining:
		raise ValueError(f"No remaining chunk FASTAs found after chunk {latest_chunk_index:04d}")
	return remaining


def subset_alignment(full_fasta, wanted_ids, output_fasta):
	wanted = set(wanted_ids)
	found = set()
	with open(output_fasta, "w", encoding="utf-8") as handle:
		for record in SeqIO.parse(full_fasta, "fasta"):
			seq_id = str(record.id).strip()
			if seq_id in wanted and seq_id not in found:
				SeqIO.write([record], handle, "fasta")
				found.add(seq_id)
	missing = sorted(wanted - found)
	if missing:
		raise ValueError(
			f"{len(missing)} resume sequence(s) were not found in {full_fasta}: {', '.join(missing[:10])}"
		)
	return output_fasta


def prepare_resume(run_dir, padded_aln, resume_dir):
	resume_path = Path(resume_dir)
	resume_path.mkdir(parents=True, exist_ok=True)
	latest_chunk_dir, latest_tree = find_latest_chunk_with_tree(run_dir)
	latest_index = _chunk_index(latest_chunk_dir)
	remaining_chunk_fastas = find_remaining_chunk_fastas(run_dir, latest_index)

	tree_out = resume_path / "resume_tree.nwk"
	ids_out = resume_path / "resume_existing_ids.txt"
	fasta_out = resume_path / "resume_alignment.fasta"
	manifest_out = resume_path / "resume_manifest.tsv"

	normalise_tree(latest_tree, tree_out)
	tree_ids = read_tree_terminals(tree_out)
	remaining_ids = []
	for chunk_fasta in remaining_chunk_fastas:
		remaining_ids.extend(read_ids_from_fasta(chunk_fasta))

	merged_ids = []
	seen = set()
	for seq_id in tree_ids + remaining_ids:
		if seq_id and seq_id not in seen:
			merged_ids.append(seq_id)
			seen.add(seq_id)

	subset_alignment(padded_aln, merged_ids, fasta_out)
	with open(ids_out, "w", encoding="utf-8") as handle:
		for seq_id in tree_ids:
			handle.write(seq_id + "\n")

	with open(manifest_out, "w", encoding="utf-8") as handle:
		handle.write("key\tvalue\n")
		handle.write(f"latest_completed_chunk\t{latest_chunk_dir.name}\n")
		handle.write(f"latest_tree_source\t{latest_tree}\n")
		handle.write(f"remaining_chunk_count\t{len(remaining_chunk_fastas)}\n")
		handle.write(f"tree_tip_count\t{len(tree_ids)}\n")
		handle.write(f"resume_alignment_count\t{len(merged_ids)}\n")

	return {
		"latest_chunk_dir": str(latest_chunk_dir),
		"latest_tree": str(latest_tree),
		"resume_tree": str(tree_out),
		"existing_ids": str(ids_out),
		"resume_alignment": str(fasta_out),
		"manifest": str(manifest_out),
		"remaining_chunks": [str(path) for path in remaining_chunk_fastas],
	}


def main():
	parser = argparse.ArgumentParser(description="Prepare resume inputs for a failed iterative USHER run")
	parser.add_argument("--run_dir", required=True, help="Existing iterative USHER output directory containing chunk_* subdirectories")
	parser.add_argument("--padded_aln", required=True, help="Full padded alignment FASTA used for the original run")
	parser.add_argument("--resume_dir", required=True, help="Directory where resume inputs should be written")
	args = parser.parse_args()

	outputs = prepare_resume(args.run_dir, args.padded_aln, args.resume_dir)
	print(f"[info] Latest completed chunk: {Path(outputs['latest_chunk_dir']).name}")
	print(f"[info] Normalised resume tree: {outputs['resume_tree']}")
	print(f"[info] Existing tree tip IDs: {outputs['existing_ids']}")
	print(f"[info] Resume alignment FASTA: {outputs['resume_alignment']}")
	print(f"[info] Manifest: {outputs['manifest']}")


if __name__ == "__main__":
	main()