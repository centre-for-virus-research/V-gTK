#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sqlite3
import subprocess
import shutil
from io import StringIO

import pandas as pd
from Bio import Phylo, SeqIO


class UsherPlacement:
	def __init__(self, padded_aln, output_dir, mmseq_cluster_dir=None, iqtree_dir=None, update_db=None, threads=1, test_mode=False, chunk_size=50000, chunk_threshold=100000):
		self.padded_aln = padded_aln
		self.output_dir = output_dir
		self.mmseq_cluster_dir = self._normalize_optional_path(mmseq_cluster_dir)
		self.iqtree_dir = self._normalize_optional_path(iqtree_dir)
		self.update_db = self._normalize_optional_path(update_db)
		self.threads = max(1, int(threads))
		self.test_mode = str(test_mode).strip() == "1" if not isinstance(test_mode, bool) else test_mode
		self.chunk_size = max(1, int(chunk_size))
		self.chunk_threshold = max(1, int(chunk_threshold))

	@staticmethod
	def _normalize_optional_path(path_value):
		if path_value is None:
			return None
		path_str = str(path_value).strip()
		if not path_str or path_str.lower() == "null" or path_str == "UNSET":
			return None
		return path_str

	@staticmethod
	def _normalize_segment(value):
		if value is None:
			return "0"
		s = str(value).strip()
		if not s:
			return "0"
		digits = ''.join(ch for ch in s if ch.isdigit())
		return digits if digits else s

	@staticmethod
	def _find_first(path_root, pattern_suffix=None, exact_name=None):
		for root, _, files in os.walk(path_root):
			for name in files:
				if exact_name and name == exact_name:
					return os.path.join(root, name)
				if pattern_suffix and name.endswith(pattern_suffix):
					return os.path.join(root, name)
		return None

	@staticmethod
	def _read_ids_from_fasta(fasta_path):
		ids = []
		with open(fasta_path, "r", encoding="utf-8") as handle:
			for line in handle:
				if line.startswith(">"):
					ids.append(line[1:].strip().split()[0])
		return ids

	@staticmethod
	def _read_tree_terminals(tree_path):
		with open(tree_path, "r", encoding="utf-8") as handle:
			newick = handle.read().strip()
		if not newick:
			return set()
		tree = Phylo.read(StringIO(newick), "newick")
		return {term.name for term in tree.get_terminals() if term.name}

	def _segment_from_padded_alignment(self):
		name = os.path.basename(self.padded_aln)
		name = re.sub(r"_dedup\.fasta$", "", name)
		match = re.search(r"(\d+)", name)
		return match.group(1) if match else "0"

	def prepare_update_assets(self):
		if not self.update_db:
			raise ValueError("update_db is required for update-mode USHER placement")
		if not os.path.isfile(self.update_db):
			raise FileNotFoundError(f"Update DB not found: {self.update_db}")
		os.makedirs(self.output_dir, exist_ok=True)

		segment = self._segment_from_padded_alignment()
		tree_out = os.path.join(self.output_dir, f"tree_segment_{segment}.nwk")
		ids_out = os.path.join(self.output_dir, f"existing_ids_segment_{segment}.txt")

		conn = sqlite3.connect(self.update_db)
		try:
			cur = conn.cursor()

			tree_written = False
			if cur.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='trees'").fetchone():
				rows = cur.execute(
					"SELECT COALESCE(segment, '0') AS segment, COALESCE(source, '') AS source, newick "
					"FROM trees WHERE newick IS NOT NULL AND TRIM(newick) != ''"
				).fetchall()
				priority = {"usher": 0, "iqtree": 1, "veryfasttree": 2}
				best = None
				for row in rows:
					row_seg = self._normalize_segment(row[0])
					if row_seg != segment:
						continue
					source = str(row[1]).strip().lower()
					score = priority.get(source, 9)
					if best is None or score < best[0]:
						best = (score, row[2])
				if best and best[1]:
					with open(tree_out, "w", encoding="utf-8") as handle:
						handle.write(str(best[1]).strip() + "\n")
					tree_written = True

			if not tree_written:
				raise ValueError(f"Missing tree for segment {segment} in update DB")

		finally:
			conn.close()

		tree_ids = sorted(self._read_tree_terminals(tree_out))

		with open(ids_out, "w", encoding="utf-8") as handle:
			for acc in tree_ids:
				handle.write(acc + "\n")

		return tree_out, ids_out

	def _load_update_missing_alignment_rows(self, segment, tree_ids, current_ids):
		conn = sqlite3.connect(self.update_db)
		try:
			df_meta = pd.read_sql_query("SELECT primary_accession, accession_type, segment FROM meta_data", conn)
			df_aln = pd.read_sql_query("SELECT * FROM sequence_alignment", conn)
			excluded = set()
			cur = conn.cursor()
			if cur.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='excluded_accessions'").fetchone():
				for (acc,) in cur.execute("SELECT primary_accession FROM excluded_accessions"):
					if acc:
						excluded.add(str(acc).strip())
		finally:
			conn.close()

		if "primary_accession" not in df_aln.columns and "sequence_id" in df_aln.columns:
			df_aln["primary_accession"] = df_aln["sequence_id"]

		if "alignment" not in df_aln.columns:
			for column in ["aligned_seq", "sequence", "aln", "alignment_seq"]:
				if column in df_aln.columns:
					df_aln["alignment"] = df_aln[column]
					break

		if "alignment" not in df_aln.columns or "primary_accession" not in df_aln.columns:
			return []

		if "alignment_name" not in df_aln.columns:
			df_aln["alignment_name"] = ""

		df_meta["segment"] = df_meta["segment"].apply(self._normalize_segment).fillna("0")
		target_meta = df_meta[df_meta["segment"] == segment].copy()
		if target_meta.empty:
			return []

		target_meta["primary_accession"] = target_meta["primary_accession"].fillna("").astype(str).str.strip()
		needed_ids = {
			acc for acc in target_meta["primary_accession"].tolist()
			if acc and acc not in excluded and acc not in tree_ids and acc not in current_ids
		}
		if not needed_ids:
			return []

		part = df_aln[df_aln["primary_accession"].fillna("").astype(str).str.strip().isin(needed_ids)].copy()
		if part.empty:
			return []

		part["primary_accession"] = part["primary_accession"].fillna("").astype(str).str.strip()
		part["alignment_name"] = part["alignment_name"].fillna("").astype(str).str.strip()
		part["alignment"] = part["alignment"].fillna("").astype(str).str.strip()
		part["_preferred_backbone"] = (
			part["alignment_name"].str.lower() == part["primary_accession"].str.lower()
		).astype(int)
		part = part.sort_values(
			by=["primary_accession", "_preferred_backbone", "alignment_name"],
			ascending=[True, False, True],
		)

		rows = []
		for _, row in part.drop_duplicates(subset=["primary_accession"]).iterrows():
			acc = str(row.get("primary_accession", "")).strip()
			seq = str(row.get("alignment", "")).strip()
			if acc and seq and seq.lower() != "nan":
				rows.append((acc, seq))
		return rows

	def build_update_alignment_input(self, tree_file):
		segment = self._segment_from_padded_alignment()
		tree_ids = self._read_tree_terminals(tree_file)
		current_ids = set(self._read_ids_from_fasta(self.padded_aln))
		missing_rows = self._load_update_missing_alignment_rows(segment, tree_ids, current_ids)
		if not missing_rows:
			return self.padded_aln

		merged_fasta = os.path.join(self.output_dir, "update_all_samples.fasta")
		with open(merged_fasta, "w", encoding="utf-8") as out_handle:
			with open(self.padded_aln, "r", encoding="utf-8") as in_handle:
				out_handle.write(in_handle.read())
			for acc, seq in missing_rows:
				out_handle.write(f">{acc}\n{seq}\n")

		print(
			f"[info] Added {len(missing_rows)} historical DB alignment(s) missing from the starter tree "
			f"for update placement on segment {segment}."
		)
		return merged_fasta

	def _count_sequences_to_place(self, alignment_fasta, ref_id, existing_ids):
		aln_ids = self._read_ids_from_fasta(alignment_fasta)
		placeable_ids = [acc for acc in aln_ids if acc != ref_id and acc not in existing_ids]
		return aln_ids, placeable_ids

	def collapse_identical_sequences(self, alignment_fasta, ref_id, existing_ids):
		existing_tree_ids = set(existing_ids)
		if ref_id:
			existing_tree_ids.add(ref_id)

		aln_ids = []
		seq_groups = {}
		for record in SeqIO.parse(alignment_fasta, "fasta"):
			seq_id = str(record.id).strip()
			if not seq_id:
				continue
			aln_ids.append(seq_id)
			seq_groups.setdefault(str(record.seq), []).append(seq_id)

		placeable_ids = []
		anchor_to_members = {}
		member_to_anchor = {}

		for group_ids in seq_groups.values():
			anchor_id = next((acc for acc in group_ids if acc in existing_tree_ids), group_ids[0])
			members_to_add = [
				acc for acc in group_ids
				if acc != anchor_id and acc not in existing_tree_ids
			]
			if members_to_add:
				anchor_to_members[anchor_id] = members_to_add
				for member_id in members_to_add:
					member_to_anchor[member_id] = anchor_id
			if anchor_id != ref_id and anchor_id not in existing_tree_ids:
				placeable_ids.append(anchor_id)

		return {
			"alignment_ids": aln_ids,
			"placeable_ids": placeable_ids,
			"anchor_to_members": anchor_to_members,
			"member_to_anchor": member_to_anchor,
		}

	def write_identical_sequence_report(self, collapse_plan):
		report_path = os.path.join(self.output_dir, "identical_sequence_groups.tsv")
		with open(report_path, "w", encoding="utf-8", newline="") as handle:
			writer = csv.writer(handle, delimiter="\t")
			writer.writerow(["anchor_id", "member_id", "anchor_requires_placement"])
			placeable_ids = set(collapse_plan["placeable_ids"])
			for anchor_id, member_ids in collapse_plan["anchor_to_members"].items():
				for member_id in member_ids:
					writer.writerow([anchor_id, member_id, int(anchor_id in placeable_ids)])
		return report_path

	@staticmethod
	def _render_progress_bar(completed, total, width=30):
		if total <= 0:
			return f"[{'-' * width}]"
		ratio = min(max(completed / total, 0), 1)
		filled = min(width, int(round(ratio * width)))
		return f"[{'#' * filled}{'-' * (width - filled)}]"

	def log_chunk_progress(self, completed_chunks, total_chunks):
		bar = self._render_progress_bar(completed_chunks, total_chunks)
		remaining_chunks = max(0, total_chunks - completed_chunks)
		print(
			f"[progress] {bar} batches complete: {completed_chunks}/{total_chunks}; "
			f"remaining: {remaining_chunks}.",
			flush=True,
		)

	@staticmethod
	def _usher_verbose_log_path(output_dir):
		return os.path.join(output_dir, "usher.verbose.log")

	def _split_alignment_into_chunks_python(self, alignment_fasta, ref_id, placeable_ids):
		if not placeable_ids:
			return []

		placeable_set = set(placeable_ids)
		ref_record = None
		for record in SeqIO.parse(alignment_fasta, "fasta"):
			if record.id == ref_id:
				ref_record = record
				break
		if ref_record is None:
			raise ValueError(f"Reference ID {ref_id} not found in {alignment_fasta}")

		chunk_dir = os.path.join(self.output_dir, "chunks")
		os.makedirs(chunk_dir, exist_ok=True)
		chunk_paths = []
		chunk_handle = None
		chunk_index = 0
		written_in_chunk = 0

		def open_chunk():
			nonlocal chunk_handle, chunk_index, written_in_chunk
			chunk_index += 1
			written_in_chunk = 0
			chunk_path = os.path.join(chunk_dir, f"chunk_{chunk_index:04d}.fasta")
			chunk_paths.append(chunk_path)
			chunk_handle = open(chunk_path, "w", encoding="utf-8")
			SeqIO.write([ref_record], chunk_handle, "fasta")

		def close_chunk():
			nonlocal chunk_handle
			if chunk_handle:
				chunk_handle.close()
				chunk_handle = None

		try:
			for record in SeqIO.parse(alignment_fasta, "fasta"):
				if record.id == ref_id or record.id not in placeable_set:
					continue
				if chunk_handle is None or written_in_chunk >= self.chunk_size:
					close_chunk()
					open_chunk()
				SeqIO.write([record], chunk_handle, "fasta")
				written_in_chunk += 1
		finally:
			close_chunk()

		return chunk_paths

	def split_alignment_into_chunks(self, alignment_fasta, ref_id, placeable_ids):
		if not placeable_ids:
			return []
		chunk_dir = os.path.join(self.output_dir, "chunks")
		os.makedirs(chunk_dir, exist_ok=True)

		seqkit = shutil.which("seqkit")
		if not seqkit:
			return self._split_alignment_into_chunks_python(alignment_fasta, ref_id, placeable_ids)

		ref_fasta = os.path.join(chunk_dir, "reference.fasta")
		ids_file = os.path.join(chunk_dir, "placeable_ids.txt")
		placeable_fasta = os.path.join(chunk_dir, "placeable_only.fasta")
		raw_chunk_dir = os.path.join(chunk_dir, "raw")
		os.makedirs(raw_chunk_dir, exist_ok=True)

		ref_record = None
		for record in SeqIO.parse(alignment_fasta, "fasta"):
			if record.id == ref_id:
				ref_record = record
				break
		if ref_record is None:
			raise ValueError(f"Reference ID {ref_id} not found in {alignment_fasta}")
		with open(ref_fasta, "w", encoding="utf-8") as handle:
			SeqIO.write([ref_record], handle, "fasta")

		with open(ids_file, "w", encoding="utf-8") as handle:
			for acc in placeable_ids:
				handle.write(acc + "\n")

		try:
			subprocess.run([seqkit, "grep", "-n", "-f", ids_file, alignment_fasta, "-o", placeable_fasta], check=True)
			subprocess.run([seqkit, "split2", "-s", str(self.chunk_size), "-O", raw_chunk_dir, placeable_fasta], check=True)
		except subprocess.CalledProcessError:
			return self._split_alignment_into_chunks_python(alignment_fasta, ref_id, placeable_ids)

		raw_chunks = sorted(
			os.path.join(raw_chunk_dir, name)
			for name in os.listdir(raw_chunk_dir)
			if name.lower().endswith((".fa", ".fasta", ".fna"))
		)
		if not raw_chunks:
			return self._split_alignment_into_chunks_python(alignment_fasta, ref_id, placeable_ids)

		chunk_paths = []
		for idx, raw_chunk in enumerate(raw_chunks, start=1):
			chunk_path = os.path.join(chunk_dir, f"chunk_{idx:04d}.fasta")
			chunk_paths.append(chunk_path)
			with open(chunk_path, "w", encoding="utf-8") as out_handle:
				with open(ref_fasta, "r", encoding="utf-8") as ref_handle:
					out_handle.write(ref_handle.read())
				with open(raw_chunk, "r", encoding="utf-8") as raw_handle:
					out_handle.write(raw_handle.read())

		return chunk_paths

	def resolve_non_update_assets(self):
		if not self.mmseq_cluster_dir or not os.path.isdir(self.mmseq_cluster_dir):
			raise FileNotFoundError(f"MMseqs cluster directory not found: {self.mmseq_cluster_dir}")
		if not self.iqtree_dir or not os.path.isdir(self.iqtree_dir):
			raise FileNotFoundError(f"IQ-TREE directory not found: {self.iqtree_dir}")

		cluster_rep = self._find_first(self.mmseq_cluster_dir, pattern_suffix="_cluster_rep.fasta")
		if not cluster_rep:
			cluster_rep = self._find_first(self.mmseq_cluster_dir, pattern_suffix=".fasta")
			if cluster_rep:
				print(f"[warn] *_cluster_rep.fasta not found; falling back to {cluster_rep}")
		if not cluster_rep:
			raise FileNotFoundError(f"No FASTA found in {self.mmseq_cluster_dir}")

		tree_file = self._find_first(self.iqtree_dir, pattern_suffix=".treefile")
		if not tree_file:
			raise FileNotFoundError(f"No treefile found in {self.iqtree_dir}")

		return cluster_rep, tree_file

	def resolve_reference_id(self, cluster_rep=None, alignment_fasta=None):
		if cluster_rep:
			cluster_ids = self._read_ids_from_fasta(cluster_rep)
			if cluster_ids:
				ref_id = cluster_ids[0]
				target_fasta = alignment_fasta or self.padded_aln
				if ref_id in set(self._read_ids_from_fasta(target_fasta)):
					return ref_id
				print(f"[warn] Reference ID '{ref_id}' is not present in {target_fasta}; using first MSA ID instead.")

		target_fasta = alignment_fasta or self.padded_aln
		aln_ids = self._read_ids_from_fasta(target_fasta)
		if not aln_ids:
			raise ValueError(f"Could not resolve reference ID from {target_fasta}")
		return aln_ids[0]

	def write_ids_file(self, rel_name, ids):
		path = os.path.join(self.output_dir, rel_name)
		with open(path, "w", encoding="utf-8") as handle:
			for item in ids:
				handle.write(str(item).strip() + "\n")
		return path

	def run_command(self, command, retry_message=None):
		try:
			subprocess.run(command, check=True)
			return True
		except subprocess.CalledProcessError:
			if retry_message:
				print(retry_message)
			return False

	def build_vcf(self, ref_id, exclude_ids_file=None, alignment_fasta=None):
		alignment_fasta = alignment_fasta or self.padded_aln
		vcf_path = os.path.join(self.output_dir, "all_samples.vcf")
		base_cmd = ["faToVcf", f"-ref={ref_id}"]
		if exclude_ids_file:
			cmd = base_cmd + [f"-excludeFile={exclude_ids_file}", alignment_fasta, vcf_path]
			if self.run_command(cmd, "[warn] faToVcf failed with -excludeFile; retrying without exclude filter."):
				return vcf_path
		cmd = base_cmd + [alignment_fasta, vcf_path]
		subprocess.run(cmd, check=True)
		return vcf_path

	@staticmethod
	def _resolve_usher_tree_output(output_dir):
		for name in ["final-tree.nh", "uncondensed-final-tree.nh"]:
			path = os.path.join(output_dir, name)
			if os.path.isfile(path):
				return path
		raise FileNotFoundError(f"USHER tree output not found in {output_dir}")

	def run_usher(self, tree_file, vcf_path, output_dir):
		os.makedirs(output_dir, exist_ok=True)
		usher_help = subprocess.run(["usher", "--help"], capture_output=True, text=True, check=False)
		usher_cmd = [
			"usher",
			"-v", vcf_path,
			"-t", tree_file,
			"-d", output_dir,
			"-o", os.path.join(output_dir, "usher.pb"),
			"-C", "-u",
		]
		if " -T " in (usher_help.stdout + usher_help.stderr):
			usher_cmd.extend(["-T", str(self.threads)])

		env = os.environ.copy()
		for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
			env[var] = str(self.threads)
		verbose_log = self._usher_verbose_log_path(output_dir)
		with open(verbose_log, "w", encoding="utf-8") as verbose_handle:
			subprocess.run(usher_cmd, check=True, env=env, stdout=verbose_handle, stderr=subprocess.STDOUT)
		return self._resolve_usher_tree_output(output_dir)

	def promote_final_usher_outputs(self, final_output_dir):
		for name in ["final-tree.nh", "uncondensed-final-tree.nh", "usher.pb", "mutation-paths.txt", "placement_stats.tsv", "all_samples.vcf"]:
			src = os.path.join(final_output_dir, name)
			if os.path.isfile(src):
				shutil.copyfile(src, os.path.join(self.output_dir, name))

	def copy_tree_outputs(self, tree_file):
		for name in ["uncondensed-final-tree.nh", "final-tree.nh"]:
			shutil.copyfile(tree_file, os.path.join(self.output_dir, name))

	@staticmethod
	def _replace_terminal_with_polytomy(clade, anchor_id, member_ids):
		for idx, child in enumerate(list(clade.clades)):
			if child.is_terminal() and child.name == anchor_id:
				replacement = child.__class__(branch_length=child.branch_length)
				replacement.clades.append(child.__class__(name=anchor_id, branch_length=0.0))
				for member_id in member_ids:
					replacement.clades.append(child.__class__(name=member_id, branch_length=0.0))
				clade.clades[idx] = replacement
				return True
			if UsherPlacement._replace_terminal_with_polytomy(child, anchor_id, member_ids):
				return True
		return False

	def expand_identical_sequence_tree_outputs(self, anchor_to_members):
		if not anchor_to_members:
			return
		for name in ["uncondensed-final-tree.nh", "final-tree.nh"]:
			tree_path = os.path.join(self.output_dir, name)
			if not os.path.isfile(tree_path):
				continue
			with open(tree_path, "r", encoding="utf-8") as handle:
				newick = handle.read().strip()
			if not newick:
				continue
			tree = Phylo.read(StringIO(newick), "newick")
			for anchor_id, member_ids in anchor_to_members.items():
				if not self._replace_terminal_with_polytomy(tree.root, anchor_id, member_ids):
					print(f"[warn] Could not find identical-sequence anchor '{anchor_id}' in {tree_path}; skipping expansion for {len(member_ids)} sequence(s).")
			with open(tree_path, "w", encoding="utf-8") as handle:
				Phylo.write(tree, handle, "newick")

	def run(self):
		os.makedirs(self.output_dir, exist_ok=True)

		cluster_rep = None
		alignment_fasta = self.padded_aln
		if self.update_db:
			tree_file, existing_ids_file = self.prepare_update_assets()
			alignment_fasta = self.build_update_alignment_input(tree_file)
		else:
			cluster_rep, tree_file = self.resolve_non_update_assets()
			centroid_ids = self._read_ids_from_fasta(cluster_rep)
			self.write_ids_file("centroid_ids.txt", centroid_ids)
			self.write_ids_file("aln_ids.txt", self._read_ids_from_fasta(self.padded_aln))
			existing_ids_file = self.write_ids_file("exclude_ids.txt", centroid_ids[1:])

		ref_id = self.resolve_reference_id(cluster_rep=cluster_rep, alignment_fasta=alignment_fasta)
		existing_ids = []
		if os.path.isfile(existing_ids_file):
			existing_ids = self._read_text_lines(existing_ids_file)
		collapse_plan = self.collapse_identical_sequences(alignment_fasta, ref_id, set(existing_ids))
		self.write_identical_sequence_report(collapse_plan)
		aln_ids = collapse_plan["alignment_ids"]
		placeable_ids = collapse_plan["placeable_ids"]
		collapsed_count = len(collapse_plan["member_to_anchor"])
		if collapsed_count:
			print(
				f"[info] Collapsed {collapsed_count} identical sequence(s) onto "
				f"{len(collapse_plan['anchor_to_members'])} anchor node(s) before UShER placement.",
				flush=True,
			)
		if len(aln_ids) <= 1:
			raise ValueError(f"Need at least 2 sequences in {alignment_fasta}, found {len(aln_ids)}")

		if not placeable_ids:
			print("[info] No sequences require UShER placement after collapsing identical sequences and checking the current tree; reusing existing tree.")
			self.copy_tree_outputs(tree_file)
			self.expand_identical_sequence_tree_outputs(collapse_plan["anchor_to_members"])
			return

		if len(placeable_ids) > self.chunk_threshold:
			chunk_fastas = self.split_alignment_into_chunks(alignment_fasta, ref_id, placeable_ids)
			print(
				f"[info] Splitting {len(placeable_ids)} placement sequence(s) into {len(chunk_fastas)} chunk(s) "
				f"of up to {self.chunk_size} sequence(s) each."
			)
			current_tree = tree_file
			final_chunk_dir = None
			print(
				f"[info] Detailed USHER output for each batch is written under {os.path.join(self.output_dir, 'chunk_*', 'usher.verbose.log')}",
				flush=True,
			)
			self.log_chunk_progress(0, len(chunk_fastas))
			for idx, chunk_fasta in enumerate(chunk_fastas, start=1):
				chunk_dir = os.path.join(self.output_dir, f"chunk_{idx:04d}")
				vcf_path = self.build_vcf(ref_id, alignment_fasta=chunk_fasta)
				current_tree = self.run_usher(current_tree, vcf_path, chunk_dir)
				final_chunk_dir = chunk_dir
				self.log_chunk_progress(idx, len(chunk_fastas))
			if final_chunk_dir:
				self.promote_final_usher_outputs(final_chunk_dir)
				self.expand_identical_sequence_tree_outputs(collapse_plan["anchor_to_members"])
			return

		exclude_ids_file = None
		exclude_ids = [
			acc for acc in aln_ids
			if acc != ref_id and acc not in set(placeable_ids)
		]
		exclude_count = len(exclude_ids)
		if exclude_count < (len(aln_ids) - 1):
			exclude_ids_file = self.write_ids_file("exclude_ids.txt", exclude_ids)
		else:
			mode = " in test mode" if self.test_mode else ""
			print(f"[warn] Exclude list would remove all non-reference sequences ({exclude_count}/{len(aln_ids)}){mode}; running faToVcf without -excludeFile.")

		vcf_path = self.build_vcf(ref_id, exclude_ids_file=exclude_ids_file, alignment_fasta=alignment_fasta)
		self.run_usher(tree_file, vcf_path, self.output_dir)
		self.expand_identical_sequence_tree_outputs(collapse_plan["anchor_to_members"])

	@staticmethod
	def _read_text_lines(path):
		with open(path, "r", encoding="utf-8") as handle:
			return [line.strip() for line in handle if line.strip()]


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Run USHER placement for either standard or update mode")
	parser.add_argument("--padded_aln", required=True, help="Padded alignment FASTA")
	parser.add_argument("--output_dir", required=True, help="Output directory")
	parser.add_argument("--mmseq_cluster_dir", default=None, help="MMseqs cluster directory for standard mode")
	parser.add_argument("--iqtree_dir", default=None, help="IQ-TREE directory for standard mode")
	parser.add_argument("--update_db", default=None, help="Existing update DB for update-mode tree and ID extraction")
	parser.add_argument("--threads", default=1, type=int, help="Thread count")
	parser.add_argument("--test_mode", default="0", help="Whether test mode is enabled (1/0)")
	parser.add_argument("--chunk_size", default=5000, type=int, help="Maximum sequences per iterative placement chunk when chunking is triggered")
	parser.add_argument("--chunk_threshold", default=10000, type=int, help="Trigger iterative chunked placement when sequences-to-place exceed this count")
	args = parser.parse_args()
	
	UsherPlacement(
		padded_aln=args.padded_aln,
		output_dir=args.output_dir,
		mmseq_cluster_dir=args.mmseq_cluster_dir,
		iqtree_dir=args.iqtree_dir,
		update_db=args.update_db,
		threads=args.threads,
		test_mode=args.test_mode,
		chunk_size=args.chunk_size,
		chunk_threshold=args.chunk_threshold,
	).run()