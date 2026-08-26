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
import segment_utils


class UsherPlacement:
	# Bases that carry phylogenetic information. N, gaps and the IUPAC ambiguity
	# codes tell usher nothing about which branch a sample belongs on, so they do
	# not count toward a sequence's quality score.
	INFORMATIVE_BASES = ("A", "C", "G", "T", "U")

	def __init__(self, padded_aln, output_dir, mmseq_cluster_dir=None, iqtree_dir=None, update_db=None, threads=1, test_mode=False, chunk_size=50000, chunk_threshold=100000, starter_tree=None, existing_ids_file=None, placement_order="quality"):
		self.padded_aln = padded_aln
		self.output_dir = output_dir
		self.mmseq_cluster_dir = self._normalize_optional_path(mmseq_cluster_dir)
		self.iqtree_dir = self._normalize_optional_path(iqtree_dir)
		self.update_db = self._normalize_optional_path(update_db)
		self.threads = max(1, int(threads))
		self.test_mode = str(test_mode).strip() == "1" if not isinstance(test_mode, bool) else test_mode
		self.chunk_size = max(1, int(chunk_size))
		self.chunk_threshold = max(1, int(chunk_threshold))
		self.starter_tree = self._normalize_optional_path(starter_tree)
		self.existing_ids_file = self._normalize_optional_path(existing_ids_file)
		self.placement_order = str(placement_order or "quality").strip().lower()
		if self.placement_order not in ("quality", "input"):
			raise ValueError(f"placement_order must be 'quality' or 'input', got {placement_order!r}")
		self._usher_supports_T = None

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
		"""Delegates to :mod:`segment_utils` - the single normalisation authority.

		This used to scrape every digit out of the string: ``4.0`` became segment 40,
		and the polymerase names inverted - ``PB2`` (segment 1) became ``2`` and ``PB1``
		(segment 2) became ``1``. Missing still maps to ``"0"`` for this call site.
		"""
		normalised = segment_utils.normalise_segment(value)
		if normalised is None or normalised.casefold() in segment_utils.PANDAS_NULL_TOKENS:
			return "0"
		return normalised

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
			meta_cols = {row[1] for row in cur.execute("PRAGMA table_info(meta_data)").fetchall()}
			if "exclusion_status" in meta_cols:
				for (acc,) in cur.execute(
					"SELECT primary_accession FROM meta_data WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) <> '' "
					"AND LOWER(TRIM(COALESCE(CAST(exclusion_status AS TEXT), ''))) NOT IN ('', '0', 'false', 'no', 'na', 'none', 'nan')"
				):
					if acc:
						excluded.add(str(acc).strip())
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

	@classmethod
	def _informative_length(cls, sequence):
		"""Unambiguous bases in a sequence: everything that is not N, a gap or an
		IUPAC ambiguity code. str.count runs in C, so this stays cheap over a
		full segment alignment."""
		seq = str(sequence).upper()
		return sum(seq.count(base) for base in cls.INFORMATIVE_BASES)

	def order_ids_by_quality(self, alignment_fasta, ids):
		"""Rank ids by informative base count, longest first.

		usher adds samples one at a time, in VCF column order, each onto the tree
		built so far, and each batch is placed onto the protobuf produced by the
		previous batch. So this ranking is the literal placement order, end to end.

		Placing well-covered sequences first means a partial sequence is later
		matched against a tree that already holds the full-length relatives it
		belongs next to, instead of being stranded on whichever branch happens to
		fit the handful of sites it covers - a misplacement that never gets
		revisited once made.
		"""
		wanted = set(ids)
		scores = {}
		for record in SeqIO.parse(alignment_fasta, "fasta"):
			seq_id = str(record.id).strip()
			if seq_id in wanted:
				scores[seq_id] = self._informative_length(record.seq)

		# Anything missing from the alignment sorts last rather than vanishing.
		# The id is a tie-break so the order is deterministic across runs.
		ordered = sorted(ids, key=lambda acc: (-scores.get(acc, -1), acc))
		return ordered, scores

	def write_placement_order_report(self, ordered_ids, scores):
		report_path = os.path.join(self.output_dir, "placement_order.tsv")
		with open(report_path, "w", encoding="utf-8", newline="") as handle:
			writer = csv.writer(handle, delimiter="\t")
			writer.writerow(["rank", "accession", "informative_bases", "batch"])
			for rank, acc in enumerate(ordered_ids, start=1):
				writer.writerow([rank, acc, scores.get(acc, ""), ((rank - 1) // self.chunk_size) + 1])
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

		ref_record = None
		for record in SeqIO.parse(alignment_fasta, "fasta"):
			if record.id == ref_id:
				ref_record = record
				break
		if ref_record is None:
			raise ValueError(f"Reference ID {ref_id} not found in {alignment_fasta}")

		chunk_dir = os.path.join(self.output_dir, "chunks")
		os.makedirs(chunk_dir, exist_ok=True)

		# placeable_ids carries the order sequences should be placed in, so a
		# sequence's rank decides its batch: batch 1 holds the highest-quality
		# sequences. Records stream straight to their batch file, so nothing
		# beyond the id -> batch map is held in memory. Order *within* each batch
		# is fixed up afterwards - see the rank-order rewrite below.
		id_to_chunk = {acc: rank // self.chunk_size for rank, acc in enumerate(placeable_ids)}
		total_chunks = (len(placeable_ids) + self.chunk_size - 1) // self.chunk_size
		chunk_paths = [
			os.path.join(chunk_dir, f"chunk_{idx + 1:04d}.fasta")
			for idx in range(total_chunks)
		]

		# Seed every batch with the reference up front, then append members.
		for chunk_path in chunk_paths:
			with open(chunk_path, "w", encoding="utf-8") as handle:
				SeqIO.write([ref_record], handle, "fasta")

		# Cap concurrent file handles so a small chunk_size over a large alignment
		# cannot exhaust the process limit; evicted batches reopen in append mode.
		max_open = 256
		handles = {}

		def handle_for(chunk_idx):
			handle = handles.get(chunk_idx)
			if handle is None:
				if len(handles) >= max_open:
					_, victim = handles.popitem()
					victim.close()
				handle = open(chunk_paths[chunk_idx], "a", encoding="utf-8")
				handles[chunk_idx] = handle
			return handle

		try:
			for record in SeqIO.parse(alignment_fasta, "fasta"):
				chunk_idx = id_to_chunk.get(str(record.id).strip())
				if chunk_idx is None:
					continue
				SeqIO.write([record], handle_for(chunk_idx), "fasta")
		finally:
			for handle in handles.values():
				handle.close()

		# The streaming pass above wrote each batch in alignment order. That is not
		# good enough: faToVcf emits VCF sample columns in FASTA order, and
		# place_onto_pb runs usher with no -s/-S/-A sort flag, so usher places
		# samples in VCF column order, one at a time, each onto the tree built so
		# far. Order within a batch therefore matters exactly as much as order
		# across batches, so rewrite each batch in rank order. Only one batch is
		# held in memory at a time.
		#
		# Do not add usher's -A/--sort-before-placement-3 here: it would re-sort
		# each batch by ambiguous-base count on its own terms and override this
		# ranking, and it can only see one batch, never the global order.
		rank_of = {acc: rank for rank, acc in enumerate(placeable_ids)}
		for chunk_path in chunk_paths:
			records = list(SeqIO.parse(chunk_path, "fasta"))
			ref_records = [r for r in records if str(r.id).strip() == ref_id]
			members = [r for r in records if str(r.id).strip() != ref_id]
			members.sort(key=lambda r: rank_of.get(str(r.id).strip(), len(rank_of)))
			with open(chunk_path, "w", encoding="utf-8") as handle:
				SeqIO.write(ref_records + members, handle, "fasta")

		return chunk_paths

	def split_alignment_into_chunks(self, alignment_fasta, ref_id, placeable_ids, ordered=False):
		if not placeable_ids:
			return []
		if ordered:
			# seqkit grep/split2 emit records in alignment order regardless of the
			# order of the id file, which would throw the quality ranking away, so
			# the ordered path stays in Python.
			return self._split_alignment_into_chunks_python(alignment_fasta, ref_id, placeable_ids)
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
			subprocess.run([seqkit, "grep", "-j", str(self.threads), "-n", "-f", ids_file, alignment_fasta, "-o", placeable_fasta], check=True)
			subprocess.run([seqkit, "split2", "-j", str(self.threads), "-s", str(self.chunk_size), "-O", raw_chunk_dir, placeable_fasta], check=True)
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

	def resolve_resume_assets(self):
		if not self.starter_tree or not self.existing_ids_file:
			raise ValueError("starter_tree and existing_ids_file are both required for resume mode")
		if not os.path.isfile(self.starter_tree):
			raise FileNotFoundError(f"Resume starter tree not found: {self.starter_tree}")
		if not os.path.isfile(self.existing_ids_file):
			raise FileNotFoundError(f"Resume existing IDs file not found: {self.existing_ids_file}")
		return self.starter_tree, self.existing_ids_file

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

	def build_vcf(self, ref_id, exclude_ids_file=None, alignment_fasta=None, vcf_path=None):
		alignment_fasta = alignment_fasta or self.padded_aln
		vcf_path = vcf_path or os.path.join(self.output_dir, "all_samples.vcf")
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

	def _append_threads(self, cmd):
		"""Add -T <threads> if this usher build accepts it (cached probe)."""
		if self._usher_supports_T is None:
			usher_help = subprocess.run(["usher", "--help"], capture_output=True, text=True, check=False)
			self._usher_supports_T = " -T " in (usher_help.stdout + usher_help.stderr)
		if self._usher_supports_T:
			cmd.extend(["-T", str(self.threads)])
		return cmd

	def _run_usher_cmd(self, cmd, output_dir, log_name="usher.verbose.log"):
		env = os.environ.copy()
		for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
			env[var] = str(self.threads)
		log_path = os.path.join(output_dir, log_name)
		with open(log_path, "w", encoding="utf-8") as log_handle:
			subprocess.run(cmd, check=True, env=env, stdout=log_handle, stderr=subprocess.STDOUT)

	def _fetch_db_alignments(self, accessions):
		"""Fetch aligned (reference-coordinate) sequences for the given accessions
		from the update DB's sequence_alignment table. Returns {accession: seq}."""
		result = {}
		accessions = [str(a).strip() for a in accessions if a and str(a).strip()]
		if not accessions or not self.update_db or not os.path.isfile(self.update_db):
			return result
		conn = sqlite3.connect(self.update_db)
		try:
			cur = conn.cursor()
			step = 900
			for start in range(0, len(accessions), step):
				chunk = accessions[start:start + step]
				placeholders = ",".join("?" * len(chunk))
				for acc, seq in cur.execute(
					f"SELECT primary_accession, alignment FROM sequence_alignment "
					f"WHERE primary_accession IN ({placeholders})", chunk
				):
					acc = str(acc).strip()
					seq = "" if seq is None else str(seq).strip()
					if acc and seq and seq.lower() != "nan" and acc not in result:
						result[acc] = seq
		finally:
			conn.close()
		return result

	def collect_backbone_fasta(self, backbone_ids, ref_id, alignment_fasta):
		"""Assemble an aligned FASTA holding the reference plus every backbone
		(existing-tree) tip's genotype, so faToVcf/usher can reconstruct the
		mutations on the backbone branches. Without this the backbone is a
		mutation-less scaffold and every new sample collapses onto the root."""
		needed = {str(a).strip() for a in backbone_ids if a and str(a).strip()}
		needed.add(ref_id)
		out_path = os.path.join(self.output_dir, "backbone_alignment.fasta")
		found = set()
		with open(out_path, "w", encoding="utf-8") as out_handle:
			# 1) whatever is already in the working alignment (always holds the ref)
			for record in SeqIO.parse(alignment_fasta, "fasta"):
				rid = str(record.id).strip()
				if rid in needed and rid not in found:
					SeqIO.write([record], out_handle, "fasta")
					found.add(rid)
			# 2) the backbone tips themselves, from the update DB alignment table
			remaining = needed - found
			if remaining:
				for acc, seq in self._fetch_db_alignments(remaining).items():
					out_handle.write(f">{acc}\n{seq}\n")
					found.add(acc)
			# 3) last resort: the on-disk padded alignment, if distinct
			remaining = needed - found
			if remaining and os.path.abspath(alignment_fasta) != os.path.abspath(self.padded_aln):
				for record in SeqIO.parse(self.padded_aln, "fasta"):
					rid = str(record.id).strip()
					if rid in remaining and rid not in found:
						SeqIO.write([record], out_handle, "fasta")
						found.add(rid)

		if ref_id not in found:
			raise ValueError(f"Reference sequence '{ref_id}' not found for backbone VCF construction")
		missing = needed - found - {ref_id}
		if missing:
			preview = ", ".join(sorted(missing)[:5])
			print(
				f"[warn] {len(missing)} backbone tip(s) have no aligned sequence available; "
				f"they will carry no mutations in the backbone (e.g. {preview}).",
				flush=True,
			)
		print(
			f"[info] Assembled backbone genotypes for {len(found) - 1} tip(s) + reference "
			f"for mutation-annotated-tree construction.",
			flush=True,
		)
		return out_path

	def build_backbone_pb(self, tree_file, backbone_ids, ref_id, alignment_fasta):
		"""Build a mutation-annotated tree (protobuf) from the backbone tree and
		its tips' genotypes. Returns the path to the saved .pb."""
		backbone_dir = os.path.join(self.output_dir, "backbone")
		os.makedirs(backbone_dir, exist_ok=True)
		backbone_fasta = self.collect_backbone_fasta(backbone_ids, ref_id, alignment_fasta)
		backbone_vcf = self.build_vcf(
			ref_id, alignment_fasta=backbone_fasta,
			vcf_path=os.path.join(backbone_dir, "backbone.vcf"),
		)
		backbone_pb = os.path.join(backbone_dir, "backbone.pb")
		cmd = ["usher", "-t", tree_file, "-v", backbone_vcf, "-d", backbone_dir, "-o", backbone_pb]
		self._append_threads(cmd)
		print("[info] Building mutation-annotated backbone protobuf ...", flush=True)
		self._run_usher_cmd(cmd, backbone_dir, log_name="usher.backbone.log")
		if not os.path.isfile(backbone_pb):
			raise FileNotFoundError(f"USHER did not produce a backbone protobuf at {backbone_pb}")
		return backbone_pb

	def place_onto_pb(self, input_pb, vcf_path, output_dir):
		"""Place the samples in vcf_path onto the mutation-annotated tree in
		input_pb, saving an updated protobuf. Returns the new .pb path so the
		next batch places against the backbone PLUS everything added here."""
		os.makedirs(output_dir, exist_ok=True)
		output_pb = os.path.join(output_dir, "usher.pb")
		cmd = [
			"usher",
			"-i", input_pb,
			"-v", vcf_path,
			"-d", output_dir,
			"-o", output_pb,
			"-C", "-u",
		]
		self._append_threads(cmd)
		self._run_usher_cmd(cmd, output_dir)
		self._resolve_usher_tree_output(output_dir)
		if not os.path.isfile(output_pb):
			raise FileNotFoundError(f"USHER did not produce an updated protobuf at {output_pb}")
		return output_pb

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
		# backbone_ids = tips of the starting tree; their genotypes must feed the
		# mutation-annotated-tree build so the backbone branches carry mutations.
		if self.starter_tree or self.existing_ids_file:
			tree_file, existing_ids_file = self.resolve_resume_assets()
			backbone_ids = self._read_text_lines(existing_ids_file)
		elif self.update_db:
			tree_file, existing_ids_file = self.prepare_update_assets()
			alignment_fasta = self.build_update_alignment_input(tree_file)
			backbone_ids = self._read_text_lines(existing_ids_file)
		else:
			cluster_rep, tree_file = self.resolve_non_update_assets()
			centroid_ids = self._read_ids_from_fasta(cluster_rep)
			self.write_ids_file("centroid_ids.txt", centroid_ids)
			self.write_ids_file("aln_ids.txt", self._read_ids_from_fasta(self.padded_aln))
			existing_ids_file = self.write_ids_file("exclude_ids.txt", centroid_ids[1:])
			# every cluster representative is a backbone (IQ-TREE) tip, ref included
			backbone_ids = centroid_ids

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

		# Build a mutation-annotated backbone protobuf ONCE. This is the fix:
		# usher places samples by the mutations they share with backbone branches,
		# so the backbone tips' genotypes must be present when the .pb is built.
		# Previously usher was handed a bare newick with a VCF that excluded the
		# backbone tips, leaving every backbone branch mutation-less -> every new
		# sample was equidistant from all branches and collapsed onto the root.
		backbone_pb = self.build_backbone_pb(tree_file, backbone_ids, ref_id, alignment_fasta)

		# Place samples in batches, carrying the protobuf forward. Because each
		# batch's placements are saved into the .pb handed to the next batch,
		# later batches place against the backbone PLUS everything already added,
		# not just the static reference set.
		ordered = self.placement_order == "quality"
		if ordered:
			placeable_ids, quality_scores = self.order_ids_by_quality(alignment_fasta, placeable_ids)
			report_path = self.write_placement_order_report(placeable_ids, quality_scores)
			if placeable_ids:
				best = quality_scores.get(placeable_ids[0], 0)
				worst = quality_scores.get(placeable_ids[-1], 0)
				print(
					f"[info] Placement ordered by sequence quality: {best} informative bases in the "
					f"first sequence down to {worst} in the last. Order written to {report_path}",
					flush=True,
				)

		chunk_fastas = self.split_alignment_into_chunks(alignment_fasta, ref_id, placeable_ids, ordered=ordered)
		print(
			f"[info] Placing {len(placeable_ids)} sequence(s) in {len(chunk_fastas)} batch(es) "
			f"of up to {self.chunk_size} sequence(s) each onto the backbone protobuf."
		)
		print(
			f"[info] Detailed USHER output for each batch is written under {os.path.join(self.output_dir, 'chunk_*', 'usher.verbose.log')}",
			flush=True,
		)
		current_pb = backbone_pb
		final_chunk_dir = None
		self.log_chunk_progress(0, len(chunk_fastas))
		for idx, chunk_fasta in enumerate(chunk_fastas, start=1):
			chunk_dir = os.path.join(self.output_dir, f"chunk_{idx:04d}")
			os.makedirs(chunk_dir, exist_ok=True)
			vcf_path = self.build_vcf(
				ref_id, alignment_fasta=chunk_fasta,
				vcf_path=os.path.join(chunk_dir, "all_samples.vcf"),
			)
			current_pb = self.place_onto_pb(current_pb, vcf_path, chunk_dir)
			final_chunk_dir = chunk_dir
			self.log_chunk_progress(idx, len(chunk_fastas))
		if final_chunk_dir:
			self.promote_final_usher_outputs(final_chunk_dir)
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
	parser.add_argument("--starter_tree", default=None, help="Existing tree to resume iterative placement from")
	parser.add_argument("--existing_ids_file", default=None, help="Text file listing tip IDs already present in --starter_tree")
	parser.add_argument("--threads", default=1, type=int, help="Thread count")
	parser.add_argument("--test_mode", default="0", help="Whether test mode is enabled (1/0)")
	parser.add_argument("--chunk_size", default=5000, type=int, help="Maximum sequences per iterative placement chunk when chunking is triggered")
	parser.add_argument("--chunk_threshold", default=10000, type=int, help="Trigger iterative chunked placement when sequences-to-place exceed this count")
	parser.add_argument(
		"--placement_order",
		default="quality",
		choices=["quality", "input"],
		help="Order sequences are placed in. 'quality' (default) places the sequences with the most "
			 "unambiguous bases first, so partial sequences are matched against a tree that already "
			 "contains their full-length relatives. 'input' keeps alignment order.",
	)
	args = parser.parse_args()
	
	UsherPlacement(
		padded_aln=args.padded_aln,
		output_dir=args.output_dir,
		mmseq_cluster_dir=args.mmseq_cluster_dir,
		iqtree_dir=args.iqtree_dir,
		update_db=args.update_db,
		starter_tree=args.starter_tree,
		existing_ids_file=args.existing_ids_file,
		threads=args.threads,
		test_mode=args.test_mode,
		chunk_size=args.chunk_size,
		chunk_threshold=args.chunk_threshold,
		placement_order=args.placement_order,
	).run()