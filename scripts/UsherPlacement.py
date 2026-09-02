#!/usr/bin/env python3
import argparse
import csv
import hashlib
import os
import re
import sqlite3
import subprocess
import shutil
import sys
import threading
from io import StringIO

import pandas as pd
from Bio import Phylo, SeqIO
import segment_utils


class UsherPlacement:
	# Bases that carry phylogenetic information. N, gaps and the IUPAC ambiguity
	# codes tell usher nothing about which branch a sample belongs on, so they do
	# not count toward a sequence's quality score.
	INFORMATIVE_BASES = ("A", "C", "G", "T", "U")

	# Filename forms that genuinely name a segment. Anchored on an explicit
	# segment token for exactly the reason segment_utils refuses to scrape digits
	# out of arbitrary text: the accession-style alignment names this pipeline
	# also produces ("AF009606.aligned_merged_MSA_dedup.fasta") are full of digits
	# that are not segments.
	_SEGMENT_IN_FILENAME_RE = re.compile(
		r"(?:^|[_.\-])(?:refset|segments|segment|sgt|seg|vrna|rna)[_.\s:#\-]*(\d+)(?:[_.\-]|$)",
		re.IGNORECASE,
	)

	# Bio.Phylo parses and writes newick recursively, and a caterpillar tree -
	# the shape iterative single-sample placement produces - is as deep as it is
	# wide. Tree I/O therefore runs on a worker thread with a large stack and a
	# raised recursion limit; 256 MB carries a 50,000-level ladder.
	_TREE_THREAD_STACK_BYTES = 256 * 1024 * 1024
	_TREE_RECURSION_PER_LEVEL = 30
	_TREE_RECURSION_BASE = 2000

	def __init__(self, padded_aln, output_dir, mmseq_cluster_dir=None, iqtree_dir=None, update_db=None, threads=1, test_mode=False, chunk_size=50000, chunk_threshold=100000, starter_tree=None, existing_ids_file=None, placement_order="quality", segment=None, min_informative_bases=1):
		self.padded_aln = padded_aln
		self.output_dir = output_dir
		self.mmseq_cluster_dir = self._normalize_optional_path(mmseq_cluster_dir)
		self.iqtree_dir = self._normalize_optional_path(iqtree_dir)
		self.update_db = self._normalize_optional_path(update_db)
		# Clamped from above as well as below: an unbounded -T is handed straight
		# to usher and to OMP_NUM_THREADS, which oversubscribes a shared node.
		self.threads = min(max(1, int(threads)), max(1, os.cpu_count() or 1))
		# Accepted for pipeline-call compatibility (vgtk-init.nf always passes
		# --test_mode); this script has no test-only behaviour of its own.
		self.test_mode = str(test_mode).strip() == "1" if not isinstance(test_mode, bool) else test_mode
		self.chunk_size = max(1, int(chunk_size))
		self.chunk_threshold = max(1, int(chunk_threshold))
		self.starter_tree = self._normalize_optional_path(starter_tree)
		self.existing_ids_file = self._normalize_optional_path(existing_ids_file)
		self.min_informative_bases = max(0, int(min_informative_bases))
		self.placement_order = str(placement_order or "quality").strip().lower()
		if self.placement_order not in ("quality", "input"):
			raise ValueError(f"placement_order must be 'quality' or 'input', got {placement_order!r}")
		segment = self._normalize_optional_path(segment)
		self.segment = self._normalize_segment(segment) if segment is not None else None
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
		"""Locate a file under ``path_root``, exact name first.

		Both passes walk in sorted order. os.walk yields directory entries in
		os.scandir order, which is inode/hash order on ext4, so without the sort
		which cluster representative or treefile gets picked is not reproducible
		across machines. The two passes also mean an exact_name match anywhere
		beats a suffix match that merely happened to be listed first.
		"""
		def walk_sorted():
			for root, dirs, files in os.walk(path_root):
				dirs.sort()
				for name in sorted(files):
					yield root, name

		if exact_name:
			for root, name in walk_sorted():
				if name == exact_name:
					return os.path.join(root, name)
		if pattern_suffix:
			for root, name in walk_sorted():
				if name.endswith(pattern_suffix):
					return os.path.join(root, name)
		return None

	@staticmethod
	def _open_text(path, mode="r"):
		"""Open a pipeline text file.

		errors="replace" on the way in: a single non-UTF-8 byte anywhere in a
		public-database header used to abort the whole run with a UnicodeDecodeError
		naming a byte offset rather than a record. One mangled character in one
		header is a far better outcome than losing the placement.
		"""
		if "r" in mode:
			return open(path, mode, encoding="utf-8", errors="replace")
		return open(path, mode, encoding="utf-8")

	@classmethod
	def _parse_fasta(cls, fasta_path):
		"""SeqIO.parse over a handle we control the decoding of.

		SeqIO.parse(path) opens the file itself with the locale default encoding
		and no error handler, so it has the same non-UTF-8 fragility as a bare
		open().
		"""
		with cls._open_text(fasta_path) as handle:
			for record in SeqIO.parse(handle, "fasta"):
				yield record

	@classmethod
	def _read_ids_from_fasta(cls, fasta_path):
		"""Sequence ids in file order.

		A header that is a bare ">" (an empty description - NCBI downloads and
		several aligners emit them) used to raise IndexError out of
		.split()[0] and take down resolve_reference_id before any placement
		work started. Such records have no id, so they are skipped.
		"""
		ids = []
		with cls._open_text(fasta_path) as handle:
			for line in handle:
				if line.startswith(">"):
					parts = line[1:].strip().split()
					if parts:
						ids.append(parts[0])
		return ids

	@staticmethod
	def _sanitise_db_sequence(accession, sequence):
		"""Make a DB-stored alignment safe to write into a FASTA.

		Rows are written with an f-string, so an alignment value containing a
		newline followed by ">" - trivially produced by a bad import into
		sequence_alignment - injected extra FASTA records, putting accessions
		that exist in no metadata table into the published tree. Alignment
		strings never legitimately contain whitespace, so all of it is dropped;
		a surviving ">" means the value is not a sequence at all.
		"""
		cleaned = "".join(str(sequence).split())
		if ">" in cleaned:
			raise ValueError(
				f"Stored alignment for {accession!r} contains a '>' and is not a sequence; "
				"refusing to write it into a FASTA"
			)
		return cleaned

	@classmethod
	def _run_with_deep_stack(cls, func, depth_hint):
		"""Run ``func`` on a worker thread sized for a ``depth_hint``-deep tree.

		Bio.Phylo's newick reader and writer both recurse once per tree level, so
		a caterpillar tree a few hundred levels deep - the shape usher produces
		when it places samples one at a time down a diverging lineage - raised
		RecursionError. Raising sys.setrecursionlimit alone is not enough: the C
		stack overflows and the interpreter dies, so the work runs on a thread
		with a large stack instead.
		"""
		limit = cls._TREE_RECURSION_BASE + max(0, int(depth_hint)) * cls._TREE_RECURSION_PER_LEVEL
		box = {}

		def worker():
			previous = sys.getrecursionlimit()
			if limit > previous:
				sys.setrecursionlimit(limit)
			try:
				box["result"] = func()
			except BaseException as exc:  # re-raised on the calling thread
				box["error"] = exc
			finally:
				sys.setrecursionlimit(previous)

		previous_stack = threading.stack_size(cls._TREE_THREAD_STACK_BYTES)
		try:
			thread = threading.Thread(target=worker)
			thread.start()
			thread.join()
		finally:
			threading.stack_size(previous_stack)
		if "error" in box:
			raise box["error"]
		return box.get("result")

	@staticmethod
	def _estimate_newick_depth(newick):
		"""Deepest bracket nesting, counted without parsing the tree."""
		depth = 0
		deepest = 0
		quote = None
		for char in newick:
			if quote:
				if char == quote:
					quote = None
			elif char in "'\"":
				quote = char
			elif char == "(":
				depth += 1
				if depth > deepest:
					deepest = depth
			elif char == ")":
				depth -= 1
		return deepest

	@classmethod
	def _read_newick(cls, newick):
		return cls._run_with_deep_stack(
			lambda: Phylo.read(StringIO(newick), "newick"),
			cls._estimate_newick_depth(newick),
		)

	@classmethod
	def _write_newick(cls, tree, path, depth_hint):
		def write():
			with cls._open_text(path, "w") as handle:
				Phylo.write(tree, handle, "newick")
		cls._run_with_deep_stack(write, depth_hint)

	@classmethod
	def _read_tree_terminals(cls, tree_path):
		with cls._open_text(tree_path) as handle:
			newick = handle.read().strip()
		if not newick:
			return set()

		# get_terminals() is a recursive depth-first walk in Bio.Phylo, so it needs
		# the same deep stack the parse does - collecting the names on the calling
		# thread would just move the RecursionError one line down.
		def read_and_collect():
			tree = Phylo.read(StringIO(newick), "newick")
			return {term.name for term in tree.get_terminals() if term.name}

		return cls._run_with_deep_stack(read_and_collect, cls._estimate_newick_depth(newick))

	def _segment_from_padded_alignment(self):
		r"""Resolve the segment this alignment belongs to.

		``--segment`` wins when supplied. Otherwise the basename is searched for
		an explicit segment token ("refset_6_aln_merged_MSA_dedup.fasta",
		"segment_4_dedup.fasta"). Anything else - including the accession-style
		"AF009606.aligned_merged_MSA_dedup.fasta" a non-segmented build produces -
		names no segment, so it resolves to "0", which is exactly what
		prepare_update_assets' COALESCE(segment, '0') stores for those builds.

		This used to be re.search(r"(\d+)") over the whole basename, taking the
		FIRST digit run. That made 'h3n2_segment_4' segment 3, 'H1N1_segment_2'
		segment 1, and 'AF009606.aligned_merged_MSA' segment 9606 - so a subtype
		in the filename silently selected another segment's starter tree, and
		update mode could never find a tree for a non-segmented virus at all. It
		is the same digit-scraping mistake segment_utils exists to end.
		"""
		if self.segment is not None:
			return self.segment
		name = os.path.basename(self.padded_aln)
		name = re.sub(r"_dedup\.fasta$", "", name)
		found = {
			segment_utils.normalise_segment(value)
			for value in self._SEGMENT_IN_FILENAME_RE.findall(name)
		}
		found.discard(None)
		if len(found) > 1:
			raise ValueError(
				f"Alignment filename {name!r} names more than one segment ({sorted(found)}); "
				"pass --segment to say which one is meant"
			)
		return found.pop() if found else "0"

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
					"FROM trees WHERE newick IS NOT NULL AND TRIM(newick) != '' "
					"ORDER BY segment, source, newick"
				).fetchall()
				priority = {"usher": 0, "iqtree": 1, "veryfasttree": 2}
				# Candidates are ranked on their own content, never on the order
				# sqlite happened to return them in. The old rule kept whichever
				# equal-priority row arrived first with no ORDER BY at all, so a
				# VACUUM or a differently-ordered insert could change which
				# backbone a rerun started from, silently.
				candidates = []
				for row in rows:
					if self._normalize_segment(row[0]) != segment:
						continue
					newick = str(row[2]).strip()
					if not newick:
						continue
					source = str(row[1]).strip().lower()
					candidates.append((priority.get(source, 9), source, newick))
				candidates.sort()
				if candidates:
					best = candidates[0]
					tied = [c for c in candidates if c[0] == best[0] and c[2] != best[2]]
					if tied:
						print(
							f"[warn] Update DB holds {len(tied) + 1} distinct '{best[1]}' trees for "
							f"segment {segment}; taking the lexicographically first so reruns agree. "
							"De-duplicate the trees table to remove the ambiguity.",
							flush=True,
						)
					with self._open_text(tree_out, "w") as handle:
						handle.write(best[2] + "\n")
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

	@staticmethod
	def _require_tables(cur, table_names):
		"""Fail with the module's own error type when a table is absent.

		prepare_update_assets guards the trees table with a sqlite_master lookup,
		but the alignment rehydration used to go straight to pd.read_sql_query, so
		a DB built by an older pipeline version died mid-run with a pandas
		DatabaseError instead of the ValueError/FileNotFoundError every other
		missing asset in this module raises.
		"""
		present = {
			row[0] for row in cur.execute(
				"SELECT name FROM sqlite_master WHERE type IN ('table', 'view')"
			)
		}
		missing = [name for name in table_names if name not in present]
		if missing:
			raise ValueError(
				f"Update DB is missing required table(s): {', '.join(sorted(missing))}"
			)

	def _load_update_missing_alignment_rows(self, segment, tree_ids, current_ids):
		conn = sqlite3.connect(self.update_db)
		try:
			cur = conn.cursor()
			self._require_tables(cur, ["meta_data", "sequence_alignment"])
			df_meta = pd.read_sql_query("SELECT primary_accession, accession_type, segment FROM meta_data", conn)
			df_aln = pd.read_sql_query("SELECT * FROM sequence_alignment", conn)
			excluded = set()
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
		with self._open_text(merged_fasta, "w") as out_handle:
			with self._open_text(self.padded_aln) as in_handle:
				existing = in_handle.read()
			out_handle.write(existing)
			# Without this the first rehydrated header is glued onto the last
			# sequence line of an alignment that lacks a trailing newline: that
			# record silently disappears AND the sequence before it is corrupted
			# with header text. Nothing downstream notices.
			if existing and not existing.endswith("\n"):
				out_handle.write("\n")
			for acc, seq in missing_rows:
				out_handle.write(f">{acc}\n{self._sanitise_db_sequence(acc, seq)}\n")

		print(
			f"[info] Added {len(missing_rows)} historical DB alignment(s) missing from the starter tree "
			f"for update placement on segment {segment}."
		)
		return merged_fasta

	def collapse_identical_sequences(self, alignment_fasta, ref_id, existing_ids):
		"""Group the alignment by genotype and pick one anchor per group.

		Three things this does that the naive version did not:

		* Sequences are bucketed on a hash of the UPPER-CASED sequence. Upper-case
		  because a soft-masked copy is the same genotype - _informative_length
		  already upper()s, so keying on the raw string had the two halves of this
		  module disagreeing about what case means. Hashed because the raw strings
		  were retained as dict keys, holding the entire alignment in memory
		  (measured at 1.04x the file) in a module whose own splitter streams.
		* A repeated accession is kept once. The same id appearing twice with
		  different sequences used to enter placeable_ids twice and put two VCF
		  sample columns with one name in front of usher.
		* A sequence with fewer than ``min_informative_bases`` unambiguous bases is
		  not placed. An all-N or all-gap record says nothing about where it
		  belongs, so usher attaches it arbitrarily; they are reported instead.
		"""
		existing_tree_ids = set(existing_ids)
		if ref_id:
			existing_tree_ids.add(ref_id)

		aln_ids = []
		seen_ids = set()
		duplicate_ids = []
		uninformative_ids = []
		seq_groups = {}
		for record in self._parse_fasta(alignment_fasta):
			seq_id = str(record.id).strip()
			if not seq_id:
				continue
			if seq_id in seen_ids:
				duplicate_ids.append(seq_id)
				continue
			seen_ids.add(seq_id)
			aln_ids.append(seq_id)
			sequence = str(record.seq).upper()
			if (seq_id not in existing_tree_ids
					and self._informative_length(sequence) < self.min_informative_bases):
				uninformative_ids.append(seq_id)
				continue
			digest = hashlib.sha256(sequence.encode("utf-8", "replace")).hexdigest()
			seq_groups.setdefault(digest, []).append(seq_id)

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
			"duplicate_ids": duplicate_ids,
			"uninformative_ids": uninformative_ids,
		}

	def write_excluded_sequence_report(self, collapse_plan):
		"""Record every sequence dropped before placement, with the reason.

		Dropping is only defensible if nothing disappears quietly, so the report
		is written whenever run() collapses, even when it is empty.
		"""
		report_path = os.path.join(self.output_dir, "excluded_sequences.tsv")
		with self._open_text(report_path, "w") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["sequence_id", "reason"])
			for seq_id in collapse_plan.get("duplicate_ids", []):
				writer.writerow([seq_id, "duplicate_accession_in_alignment"])
			for seq_id in collapse_plan.get("uninformative_ids", []):
				writer.writerow([seq_id, "fewer_than_min_informative_bases"])
		return report_path

	def write_identical_sequence_report(self, collapse_plan):
		report_path = os.path.join(self.output_dir, "identical_sequence_groups.tsv")
		with self._open_text(report_path, "w") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
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
		for record in self._parse_fasta(alignment_fasta):
			seq_id = str(record.id).strip()
			if seq_id in wanted and seq_id not in scores:
				scores[seq_id] = self._informative_length(record.seq)

		# Anything missing from the alignment sorts last rather than vanishing.
		# The id is a tie-break so the order is deterministic across runs.
		ordered = sorted(ids, key=lambda acc: (-scores.get(acc, -1), acc))
		return ordered, scores

	def write_placement_order_report(self, ordered_ids, scores, chunk_size=None):
		"""``chunk_size`` is the size batches were ACTUALLY split at, which is not
		self.chunk_size when the run stays under chunk_threshold. Passing the
		effective value keeps the report's batch column agreeing with the chunk
		files it describes."""
		chunk_size = max(1, int(chunk_size or self.chunk_size))
		report_path = os.path.join(self.output_dir, "placement_order.tsv")
		with self._open_text(report_path, "w") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["rank", "accession", "informative_bases", "batch"])
			for rank, acc in enumerate(ordered_ids, start=1):
				writer.writerow([rank, acc, scores.get(acc, ""), ((rank - 1) // chunk_size) + 1])
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

	def _resolve_chunk_size(self, placeable_count):
		"""The batch size run() should split at.

		--chunk_threshold is documented as "Trigger iterative chunked placement
		when sequences-to-place exceed this count". It was stored, clamped and
		then never read: placement was always chunked at chunk_size, so the flag
		did nothing and small runs paid for an extra usher round-trip per batch.

		This is a run()-level policy only. The splitters take an explicit
		chunk_size and always honour it, so "split this into batches of N" stays
		a thing a caller can ask for directly.
		"""
		if placeable_count <= self.chunk_threshold:
			return max(1, placeable_count)
		return self.chunk_size

	def _prepare_chunk_dir(self):
		"""A clean chunks/ directory.

		The seqkit path collects its parts with os.listdir, so leftover part files
		from a previous run into the same output directory - which is exactly what
		resume mode does - became chunk_0001 and shifted every real batch down,
		feeding usher sequences this run never selected.
		"""
		chunk_dir = os.path.join(self.output_dir, "chunks")
		if os.path.isdir(chunk_dir):
			shutil.rmtree(chunk_dir)
		os.makedirs(chunk_dir, exist_ok=True)
		return chunk_dir

	def _present_placeable_ids(self, alignment_fasta, placeable_ids):
		"""placeable_ids restricted to records the alignment actually holds, in
		the given order, each id once.

		total_chunks was derived from len(placeable_ids) while only ids found in
		the alignment were streamed into a chunk file, so an id present in the
		tree bookkeeping but absent from the alignment produced a trailing chunk
		holding nothing but the reference - still handed to faToVcf and usher as
		a real batch."""
		present = {str(r.id).strip() for r in self._parse_fasta(alignment_fasta)}
		missing = [acc for acc in dict.fromkeys(placeable_ids) if acc not in present]
		if missing:
			preview = ", ".join(missing[:5])
			print(
				f"[warn] {len(missing)} id(s) queued for placement are not in {alignment_fasta} "
				f"and will be skipped (e.g. {preview}).",
				flush=True,
			)
		return [acc for acc in dict.fromkeys(placeable_ids) if acc in present]

	def _split_alignment_into_chunks_python(self, alignment_fasta, ref_id, placeable_ids, chunk_size=None):
		if not placeable_ids:
			return []

		placeable_ids = self._present_placeable_ids(alignment_fasta, placeable_ids)
		if not placeable_ids:
			return []
		chunk_size = max(1, int(chunk_size or self.chunk_size))

		ref_record = None
		for record in self._parse_fasta(alignment_fasta):
			if str(record.id).strip() == ref_id:
				ref_record = record
				break
		if ref_record is None:
			raise ValueError(f"Reference ID {ref_id} not found in {alignment_fasta}")

		chunk_dir = self._prepare_chunk_dir()

		# placeable_ids carries the order sequences should be placed in, so a
		# sequence's rank decides its batch: batch 1 holds the highest-quality
		# sequences. Records stream straight to their batch file, so nothing
		# beyond the id -> batch map is held in memory. Order *within* each batch
		# is fixed up afterwards - see the rank-order rewrite below.
		id_to_chunk = {acc: rank // chunk_size for rank, acc in enumerate(placeable_ids)}
		total_chunks = (len(placeable_ids) + chunk_size - 1) // chunk_size
		chunk_paths = [
			os.path.join(chunk_dir, f"chunk_{idx + 1:04d}.fasta")
			for idx in range(total_chunks)
		]

		# Seed every batch with the reference up front, then append members.
		for chunk_path in chunk_paths:
			with self._open_text(chunk_path, "w") as handle:
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

		written = set()
		try:
			for record in self._parse_fasta(alignment_fasta):
				seq_id = str(record.id).strip()
				chunk_idx = id_to_chunk.get(seq_id)
				if chunk_idx is None or seq_id in written:
					continue
				written.add(seq_id)
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
			records = list(self._parse_fasta(chunk_path))
			ref_records = [r for r in records if str(r.id).strip() == ref_id][:1]
			members = [r for r in records if str(r.id).strip() != ref_id]
			members.sort(key=lambda r: rank_of.get(str(r.id).strip(), len(rank_of)))
			with self._open_text(chunk_path, "w") as handle:
				SeqIO.write(ref_records + members, handle, "fasta")

		return chunk_paths

	@staticmethod
	def _natural_sort_key(path):
		"""Sort seqkit part files by their number, not as text.

		seqkit pads part numbers to three digits, so once a split produces 1000+
		parts a plain sorted() puts 'part_1000' next to 'part_100' and ahead of
		'part_999'. Batches then reached usher out of order, destroying the
		placement ranking and desynchronising chunk_NNNN directories from their
		contents."""
		name = os.path.basename(path)
		return [int(token) if token.isdigit() else token.lower()
				for token in re.split(r"(\d+)", name)]

	def split_alignment_into_chunks(self, alignment_fasta, ref_id, placeable_ids, ordered=False, chunk_size=None):
		if not placeable_ids:
			return []
		if ordered:
			# seqkit grep/split2 emit records in alignment order regardless of the
			# order of the id file, which would throw the quality ranking away, so
			# the ordered path stays in Python.
			return self._split_alignment_into_chunks_python(alignment_fasta, ref_id, placeable_ids, chunk_size)

		seqkit = shutil.which("seqkit")
		if not seqkit:
			return self._split_alignment_into_chunks_python(alignment_fasta, ref_id, placeable_ids, chunk_size)

		placeable_ids = self._present_placeable_ids(alignment_fasta, placeable_ids)
		if not placeable_ids:
			return []
		chunk_size = max(1, int(chunk_size or self.chunk_size))
		chunk_dir = self._prepare_chunk_dir()

		ref_fasta = os.path.join(chunk_dir, "reference.fasta")
		ids_file = os.path.join(chunk_dir, "placeable_ids.txt")
		placeable_fasta = os.path.join(chunk_dir, "placeable_only.fasta")
		raw_chunk_dir = os.path.join(chunk_dir, "raw")
		os.makedirs(raw_chunk_dir, exist_ok=True)

		ref_record = None
		for record in self._parse_fasta(alignment_fasta):
			if str(record.id).strip() == ref_id:
				ref_record = record
				break
		if ref_record is None:
			raise ValueError(f"Reference ID {ref_id} not found in {alignment_fasta}")
		with self._open_text(ref_fasta, "w") as handle:
			SeqIO.write([ref_record], handle, "fasta")

		with self._open_text(ids_file, "w") as handle:
			for acc in placeable_ids:
				handle.write(acc + "\n")

		try:
			subprocess.run([seqkit, "grep", "-j", str(self.threads), "-n", "-f", ids_file, alignment_fasta, "-o", placeable_fasta], check=True)
			subprocess.run([seqkit, "split2", "-j", str(self.threads), "-s", str(chunk_size), "-O", raw_chunk_dir, placeable_fasta], check=True)
		except (subprocess.CalledProcessError, FileNotFoundError):
			return self._split_alignment_into_chunks_python(alignment_fasta, ref_id, placeable_ids, chunk_size)

		raw_chunks = sorted(
			(os.path.join(raw_chunk_dir, name)
			 for name in os.listdir(raw_chunk_dir)
			 if name.lower().endswith((".fa", ".fasta", ".fna"))),
			key=self._natural_sort_key,
		)
		if not raw_chunks:
			return self._split_alignment_into_chunks_python(alignment_fasta, ref_id, placeable_ids, chunk_size)

		chunk_paths = []
		for idx, raw_chunk in enumerate(raw_chunks, start=1):
			chunk_path = os.path.join(chunk_dir, f"chunk_{idx:04d}.fasta")
			with self._open_text(chunk_path, "w") as out_handle:
				with self._open_text(ref_fasta) as ref_handle:
					ref_text = ref_handle.read()
				out_handle.write(ref_text)
				if ref_text and not ref_text.endswith("\n"):
					out_handle.write("\n")
				with self._open_text(raw_chunk) as raw_handle:
					out_handle.write(raw_handle.read())
			# A batch holding only the reference has nothing to place; handing it
			# to faToVcf and usher is a wasted round-trip at best.
			if len(self._read_ids_from_fasta(chunk_path)) > 1:
				chunk_paths.append(chunk_path)
			else:
				os.remove(chunk_path)

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
		with self._open_text(path, "w") as handle:
			for item in ids:
				handle.write(str(item).strip() + "\n")
		return path

	def run_command(self, command, retry_message=None):
		"""Run a command, reporting failure rather than raising.

		FileNotFoundError is caught alongside CalledProcessError: check=True only
		suppresses a non-zero exit, so a tool that is not on PATH at all - a wrong
		conda env, a partially built container - used to escape this handler as a
		bare traceback instead of the failure path the caller was written for.
		"""
		try:
			subprocess.run(command, check=True)
			return True
		except FileNotFoundError:
			print(f"[warn] Command not found: {command[0]}")
			if retry_message:
				print(retry_message)
			return False
		except subprocess.CalledProcessError:
			if retry_message:
				print(retry_message)
			return False

	@staticmethod
	def _reject_option_like(**values):
		"""faToVcf takes its alignment and output positionally and has no '--'
		separator, so a value beginning with '-' is parsed as an option. ref_id
		can be whatever the first alignment record happens to be, and '-' is legal
		in a FASTA header, so this is reachable from input data."""
		for name, value in values.items():
			if value is not None and str(value).startswith("-"):
				raise ValueError(
					f"{name} {value!r} starts with '-' and would be read as a command-line "
					"option by faToVcf; rename it or pass an explicit path"
				)

	def build_vcf(self, ref_id, exclude_ids_file=None, alignment_fasta=None, vcf_path=None):
		alignment_fasta = alignment_fasta or self.padded_aln
		vcf_path = vcf_path or os.path.join(self.output_dir, "all_samples.vcf")
		self._reject_option_like(
			ref_id=ref_id, alignment_fasta=alignment_fasta,
			vcf_path=vcf_path, exclude_ids_file=exclude_ids_file,
		)
		base_cmd = ["faToVcf", f"-ref={ref_id}"]
		if exclude_ids_file:
			# No unfiltered fallback. -excludeFile is what keeps backbone tips out
			# of the placement VCF, so retrying without it silently re-submitted
			# every existing tree tip for placement and duplicated them in the tree.
			cmd = base_cmd + [f"-excludeFile={exclude_ids_file}", alignment_fasta, vcf_path]
			if not self.run_command(cmd):
				raise RuntimeError(
					f"faToVcf failed with -excludeFile={exclude_ids_file}. Refusing to retry "
					"without it: that would re-place every sequence the exclude file exists "
					"to hold back."
				)
			return vcf_path
		cmd = base_cmd + [alignment_fasta, vcf_path]
		if not self.run_command(cmd):
			raise RuntimeError(f"faToVcf failed building {vcf_path}")
		return vcf_path

	@staticmethod
	def _resolve_usher_tree_output(output_dir):
		for name in ["final-tree.nh", "uncondensed-final-tree.nh"]:
			path = os.path.join(output_dir, name)
			if os.path.isfile(path):
				return path
		raise FileNotFoundError(f"USHER tree output not found in {output_dir}")

	def _append_threads(self, cmd):
		"""Add -T <threads> if this usher build accepts it (cached probe).

		check=False only suppresses a non-zero exit; a usher that is not on PATH
		raises FileNotFoundError from the probe itself, which used to surface as an
		unhandled exception in the middle of thread detection.
		"""
		if self._usher_supports_T is None:
			try:
				usher_help = subprocess.run(["usher", "--help"], capture_output=True, text=True, check=False)
			except FileNotFoundError:
				raise RuntimeError(
					"usher is not on PATH. Activate the pipeline environment (environment.yml) "
					"before running USHER placement."
				)
			self._usher_supports_T = " -T " in (usher_help.stdout + usher_help.stderr)
		if self._usher_supports_T:
			cmd.extend(["-T", str(self.threads)])
		return cmd

	def _run_usher_cmd(self, cmd, output_dir, log_name="usher.verbose.log"):
		env = os.environ.copy()
		for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
			env[var] = str(self.threads)
		log_path = os.path.join(output_dir, log_name)
		with self._open_text(log_path, "w") as log_handle:
			try:
				subprocess.run(cmd, check=True, env=env, stdout=log_handle, stderr=subprocess.STDOUT)
			except FileNotFoundError:
				raise RuntimeError(
					f"{cmd[0]} is not on PATH. Activate the pipeline environment "
					"(environment.yml) before running USHER placement."
				)

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
		with self._open_text(out_path, "w") as out_handle:
			# 1) whatever is already in the working alignment (always holds the ref)
			for record in self._parse_fasta(alignment_fasta):
				rid = str(record.id).strip()
				if rid in needed and rid not in found:
					SeqIO.write([record], out_handle, "fasta")
					found.add(rid)
			# 2) the backbone tips themselves, from the update DB alignment table
			remaining = needed - found
			if remaining:
				for acc, seq in self._fetch_db_alignments(remaining).items():
					out_handle.write(f">{acc}\n{self._sanitise_db_sequence(acc, seq)}\n")
					found.add(acc)
			# 3) last resort: the on-disk padded alignment, if distinct
			remaining = needed - found
			if remaining and os.path.abspath(alignment_fasta) != os.path.abspath(self.padded_aln):
				for record in self._parse_fasta(self.padded_aln):
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
		# Clear any output left by a previous run into this directory: otherwise a
		# usher invocation that silently produces nothing passes the
		# _resolve_usher_tree_output / isfile checks below on stale files.
		for stale in ("final-tree.nh", "uncondensed-final-tree.nh", "usher.pb"):
			stale_path = os.path.join(output_dir, stale)
			if os.path.isfile(stale_path):
				os.remove(stale_path)
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
		"""Publish an unchanged starter tree as this run's result.

		The natural resume invocation reuses the previous run's output directory
		("--starter_tree out/final-tree.nh --output_dir out"), which makes source
		and destination the same path; shutil.copyfile then raised SameFileError
		after the run had already decided there was nothing to place.
		"""
		for name in ["uncondensed-final-tree.nh", "final-tree.nh"]:
			dest = os.path.join(self.output_dir, name)
			if os.path.abspath(tree_file) == os.path.abspath(dest):
				continue
			shutil.copyfile(tree_file, dest)

	@staticmethod
	def _make_polytomy(node, anchor_id, member_ids):
		"""An internal node on the anchor's original branch, carrying the anchor
		and its identical members as zero-length tips. Keeping the branch length
		on the new internal node is what leaves root-to-tip distances unchanged."""
		cls = node.__class__
		replacement = cls(branch_length=node.branch_length)
		replacement.clades.append(cls(name=anchor_id, branch_length=0.0))
		for member_id in member_ids:
			replacement.clades.append(cls(name=member_id, branch_length=0.0))
		return replacement

	@classmethod
	def _expand_anchors_in_tree(cls, root, anchor_to_members):
		"""Attach every anchor's members in ONE iterative pass. Returns the set of
		anchors that were found.

		Was a recursive search re-run per anchor: O(anchors x tree), and it raised
		RecursionError on a caterpillar tree only a few hundred levels deep - the
		shape iterative single-sample placement produces. It also only ever looked
		at clade.clades, so an anchor that was the root itself was never found and
		its duplicates were dropped from the tree behind a [warn] and a zero exit.
		"""
		remaining = {a: list(m) for a, m in anchor_to_members.items() if m}
		expanded = set()
		if not remaining:
			return expanded

		# The root has no parent to swap it out of, so it is converted in place.
		if root.is_terminal() and root.name in remaining:
			anchor_id = root.name
			polytomy = cls._make_polytomy(root, anchor_id, remaining.pop(anchor_id))
			root.name = None
			root.clades.extend(polytomy.clades)
			expanded.add(anchor_id)

		stack = [root]
		while stack and remaining:
			node = stack.pop()
			for idx, child in enumerate(node.clades):
				if child.is_terminal() and child.name in remaining:
					anchor_id = child.name
					node.clades[idx] = cls._make_polytomy(child, anchor_id, remaining.pop(anchor_id))
					expanded.add(anchor_id)
				else:
					stack.append(child)
		return expanded

	@classmethod
	def _replace_terminal_with_polytomy(cls, clade, anchor_id, member_ids):
		"""Single-anchor wrapper kept for callers outside this module."""
		return anchor_id in cls._expand_anchors_in_tree(clade, {anchor_id: member_ids})

	def expand_identical_sequence_tree_outputs(self, anchor_to_members):
		if not anchor_to_members:
			return
		for name in ["uncondensed-final-tree.nh", "final-tree.nh"]:
			tree_path = os.path.join(self.output_dir, name)
			if not os.path.isfile(tree_path):
				continue
			with self._open_text(tree_path) as handle:
				newick = handle.read().strip()
			if not newick:
				continue
			depth = self._estimate_newick_depth(newick)
			tree = self._read_newick(newick)
			expanded = self._expand_anchors_in_tree(tree.root, anchor_to_members)
			for anchor_id, member_ids in anchor_to_members.items():
				if member_ids and anchor_id not in expanded:
					print(
						f"[warn] Could not find identical-sequence anchor '{anchor_id}' in "
						f"{tree_path}; skipping expansion for {len(member_ids)} sequence(s)."
					)
			# +2 levels: every expanded anchor gains an internal node above it.
			self._write_newick(tree, tree_path, depth + 2)

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
		excluded_report = self.write_excluded_sequence_report(collapse_plan)
		aln_ids = collapse_plan["alignment_ids"]
		placeable_ids = collapse_plan["placeable_ids"]
		if collapse_plan["duplicate_ids"]:
			print(
				f"[warn] {len(collapse_plan['duplicate_ids'])} repeated accession(s) in "
				f"{alignment_fasta} were kept once each; see {excluded_report}",
				flush=True,
			)
		if collapse_plan["uninformative_ids"]:
			print(
				f"[warn] {len(collapse_plan['uninformative_ids'])} sequence(s) have fewer than "
				f"{self.min_informative_bases} unambiguous base(s) and carry no placement signal; "
				f"they are not placed. See {excluded_report}",
				flush=True,
			)
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
		# One batch unless the run is genuinely large: --chunk_threshold finally
		# means what its help text says.
		effective_chunk_size = self._resolve_chunk_size(len(placeable_ids))
		if ordered:
			placeable_ids, quality_scores = self.order_ids_by_quality(alignment_fasta, placeable_ids)
			report_path = self.write_placement_order_report(
				placeable_ids, quality_scores, chunk_size=effective_chunk_size)
			if placeable_ids:
				best = quality_scores.get(placeable_ids[0], 0)
				worst = quality_scores.get(placeable_ids[-1], 0)
				print(
					f"[info] Placement ordered by sequence quality: {best} informative bases in the "
					f"first sequence down to {worst} in the last. Order written to {report_path}",
					flush=True,
				)

		chunk_fastas = self.split_alignment_into_chunks(
			alignment_fasta, ref_id, placeable_ids, ordered=ordered, chunk_size=effective_chunk_size)
		if not chunk_fastas:
			print("[info] Nothing left to place after checking the alignment; reusing existing tree.")
			self.copy_tree_outputs(tree_file)
			self.expand_identical_sequence_tree_outputs(collapse_plan["anchor_to_members"])
			return
		print(
			f"[info] Placing {len(placeable_ids)} sequence(s) in {len(chunk_fastas)} batch(es) "
			f"of up to {effective_chunk_size} sequence(s) each onto the backbone protobuf."
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

	@classmethod
	def _read_text_lines(cls, path):
		with cls._open_text(path) as handle:
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
		"--segment",
		default=None,
		help="Segment this alignment belongs to. Authoritative when given. Without it the "
			 "segment is read from an explicit segment token in the alignment filename "
			 "(refset_6_aln..., segment_4_dedup.fasta); an accession-style name such as "
			 "AF009606.aligned_merged_MSA_dedup.fasta names no segment and resolves to '0'.",
	)
	parser.add_argument(
		"--min_informative_bases",
		default=1,
		type=int,
		help="Skip placement for sequences with fewer than this many unambiguous (A/C/G/T/U) "
			 "bases. They carry no signal about where they belong, so usher attaches them "
			 "arbitrarily. Skipped sequences are listed in excluded_sequences.tsv. 0 disables.",
	)
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
		segment=args.segment,
		min_informative_bases=args.min_informative_bases,
	).run()