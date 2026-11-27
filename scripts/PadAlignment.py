import os
import shutil
import argparse
import re
import sqlite3
import pandas as pd
from os.path import join
from Bio import SeqIO
from Bio.Seq import Seq
from ExportRefListFromUpdateDb import load_master_accessions, load_master_accessions_from_file, load_reference_file_table, load_reference_rows
import accession_utils
import segment_utils
import projectability

class PadAlignment:
	def __init__(self, reference_alignment, input_dir, base_dir, output_dir, keep_intermediate_files, new_outputfile=False, segment_manifest_out=None, strict_segment_backbone=False, update_db=None, skip_ids=None):
		self.reference_alignment = reference_alignment
		self.input_dir = input_dir
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.keep_intermediate_files = keep_intermediate_files
		self.new_outputfile = new_outputfile
		self.segment_manifest_out = segment_manifest_out
		self.strict_segment_backbone = strict_segment_backbone
		self.update_db = self._normalize_optional_path(update_db)
		self.segment_manifest_rows = []
		self.skip_ids: set = self._load_skip_ids(skip_ids)
		# Union of every backbone's accessions across all masters, filled in by
		# process_all_masters. Used only to scope the orphan warning run-wide
		# instead of per-segment; None means 'driven with -r, single backbone'.
		self.run_projectable_refs = None
		self._merged_already_deduped = False

	@staticmethod
	def _load_skip_ids(skip_ids_arg) -> set:
		"""Accessions to exclude from the subalignments, in canonical (bare) form.

		The ids come from ``CollectFilteredSequences``' ``filtered_sequences_ids.txt``
		and used to be kept exactly as written, while the record ids they are matched
		against in :meth:`process_master_alignment` had their version scraped off. A
		versioned line - ``NC_001542.1`` - could therefore never match anything and the
		sequence was padded in rather than skipped. That failure is *open*: a sequence
		excluded upstream reached the tree and the DB looking as though it had aligned
		cleanly. Both sides run through :mod:`accession_utils` now, so the comparison is
		between two canonical keys. A non-accession label - a cluster or group name that
		happens to contain a dot - survives verbatim and still matches itself.
		"""
		if not skip_ids_arg:
			return set()
		path_str = str(skip_ids_arg).strip()
		if not path_str or path_str.lower() in ('null', 'none', 'unset'):
			return set()
		if not os.path.isfile(path_str):
			print(f"[warn] skip_ids file not found: {path_str}; proceeding without skipping.")
			return set()
		ids = set()
		with open(path_str, 'r', encoding='utf-8') as fh:
			for line in fh:
				raw = line.strip().split()[0] if line.strip() else ''
				acc = accession_utils.normalise_accession(raw)
				if acc:
					ids.add(acc)
		print(f"[skip_ids] Loaded {len(ids)} accessions to skip from {path_str}")
		return ids

	@staticmethod
	def _normalize_optional_path(path_value):
		if path_value is None:
			return None
		path_str = str(path_value).strip()
		if not path_str or path_str.lower() == 'null' or path_str == 'UNSET':
			return None
		return path_str

	def get_master_list(self, master_acc):
		"""Master accessions in canonical (bare) form.

		The other side of this lookup is the ``reference_aln/<master>/`` tree, which
		NextalignAlignment names with :func:`accession_utils.accession_from_filename`
		and is therefore always bare. A reference list or update DB carrying
		``NC_001542.1`` used to be joined against a directory called ``NC_001542``,
		miss, and report "Reference alignment for ... not found" - which in strict
		update mode aborts the run and otherwise drops that master's whole segment
		out of the merged MSA. Both sides are canonical now.
		"""
		return [
			accession_utils.normalise_accession(master)
			for master in self._raw_master_list(master_acc)
			if accession_utils.normalise_accession(master)
		]

	def _raw_master_list(self, master_acc):
		if self.update_db:
			return load_master_accessions(self.update_db)
		if os.path.isfile(master_acc):
			try:
				return load_master_accessions_from_file(master_acc)
			except:
				return []
		else:
			return [x.strip() for x in master_acc.split(',') if x.strip()]

	def get_master_segment_map(self, master_acc):
		"""``{master accession: segment}``, keyed the way :meth:`get_master_list` reads.

		The keys have to be canonical for the same reason the master list is: this map
		is looked up in :meth:`process_all_masters` with an entry *from* that list, so
		one side carrying ``NC_001542.1`` and the other ``NC_001542`` yields segment
		``None`` - and a master with no segment picks the wrong backbone or, in strict
		mode, none at all.
		"""
		if self.update_db:
			master_segment = {}
			refs = load_reference_rows(self.update_db)
			masters = refs[refs["accession_type"].fillna("").astype(str).str.lower() == "master"]
			if masters.empty:
				masters = refs
			for _, row in masters.iterrows():
				acc = accession_utils.normalise_accession(row["primary_accession"])
				segment = self._normalize_segment(row.get("segment")) or "0"
				if acc:
					master_segment[acc] = segment
			return master_segment

		master_segment = {}
		if not os.path.isfile(master_acc):
			return master_segment
		try:
			df = load_reference_file_table(master_acc)
		except Exception:
			return master_segment

		if not df["segment"].astype(str).str.strip().ne("").any():
			return master_segment

		masters = df[df["accession_type"].astype(str).str.strip().str.lower() == 'master']
		for _, row in masters.iterrows():
			acc = accession_utils.normalise_accession(row["primary_accession"])
			if acc:
				master_segment[acc] = self._normalize_segment(row["segment"])
		return master_segment

	@staticmethod
	def _normalize_segment(segment_value):
		"""Delegates to :mod:`segment_utils` - the single normalisation authority.

		This used to scrape every digit out of the string, which turned ``4.0`` into
		segment 40 and, worse, inverted the polymerase segments: ``PB2`` (segment 1)
		became ``2`` and ``PB1`` (segment 2) became ``1``. Missing still maps to
		``None`` here, because callers of this particular helper depend on it.
		"""
		normalised = segment_utils.normalise_segment(segment_value)
		if normalised is None or normalised.casefold() in segment_utils.PANDAS_NULL_TOKENS:
			return None
		return normalised

	def find_precomputed_reference_alignment(self, precomputed_ref_dir, segment_value):
		"""
		Find a precomputed segment alignment file in precomputed_ref_dir.
		Expected canonical naming is refset_<segment>_aln.fasta, with a
		numeric-segment fallback for custom naming.

		Delegates to :func:`projectability.find_precomputed_reference_alignment`.
		CollectFilteredSequences has to resolve the *same* file to decide whether a
		query is projectable; when the two carried separate copies of this logic they
		drifted, and the stricter copy silently deleted 633,988 sequences. One
		implementation, two callers.
		"""
		return projectability.find_precomputed_reference_alignment(precomputed_ref_dir, segment_value)

	def export_update_backbones(self, output_dir):
		if not self.update_db:
			return None
		if not os.path.isfile(self.update_db):
			raise FileNotFoundError(f"Update DB not found: {self.update_db}")

		os.makedirs(output_dir, exist_ok=True)
		conn = sqlite3.connect(self.update_db)
		try:
			df_meta = pd.read_sql_query("SELECT primary_accession, accession_type, segment FROM meta_data", conn)
			df_aln = pd.read_sql_query("SELECT * FROM sequence_alignment", conn)
		finally:
			conn.close()

		if "primary_accession" not in df_aln.columns and "sequence_id" in df_aln.columns:
			df_aln["primary_accession"] = df_aln["sequence_id"]

		if "alignment" not in df_aln.columns:
			for column in ["aligned_seq", "sequence", "aln", "alignment_seq"]:
				if column in df_aln.columns:
					df_aln["alignment"] = df_aln[column]
					break

		if "alignment" not in df_aln.columns:
			raise ValueError("sequence_alignment must contain an alignment-like column")

		if "alignment_name" not in df_aln.columns:
			df_aln["alignment_name"] = "0"

		df_meta["segment"] = df_meta["segment"].apply(self._normalize_segment).fillna("0")
		df_meta["accession_type"] = df_meta["accession_type"].fillna("query").astype(str)

		df_join = df_aln.merge(
			df_meta[["primary_accession", "accession_type", "segment"]],
			on="primary_accession",
			how="left",
		)

		if "segment" not in df_join.columns:
			if "segment_y" in df_join.columns:
				df_join["segment"] = df_join["segment_y"]
			elif "segment_x" in df_join.columns:
				df_join["segment"] = df_join["segment_x"]
			else:
				df_join["segment"] = "0"

		df_join["segment"] = df_join["segment"].apply(self._normalize_segment).fillna("0")
		df_join = df_join[df_join["accession_type"].fillna("").str.lower().isin(["master", "reference"])]

		written = 0
		for segment, part in df_join.groupby("segment"):
			part = part.copy()
			part["_alignment_name"] = part["alignment_name"].fillna("").astype(str).str.strip()
			part["_primary_accession"] = part["primary_accession"].fillna("").astype(str).str.strip()
			part["_preferred_backbone"] = (
				part["_alignment_name"].str.lower() == part["_primary_accession"].str.lower()
			).astype(int)
			part = part.sort_values(
				by=["_primary_accession", "_preferred_backbone", "_alignment_name"],
				ascending=[True, False, True],
			)
			rows = []
			for _, row in part.drop_duplicates(subset=["primary_accession"]).iterrows():
				acc = str(row.get("primary_accession", "")).strip()
				seq = str(row.get("alignment", "")).strip()
				if acc and seq and seq.lower() != "nan":
					rows.append((acc, seq))
			if not rows:
				continue
			out_fa = os.path.join(output_dir, f"refset_{segment}_aln.fasta")
			with open(out_fa, "w", encoding="utf-8") as handle:
				for acc, seq in rows:
					handle.write(f">{acc}\n{seq}\n")
			written += 1

		return output_dir if written else None

	def resolve_precomputed_ref_dir(self, precomputed_ref_dir):
		precomputed_ref_dir = self._normalize_optional_path(precomputed_ref_dir)
		if precomputed_ref_dir:
			if not os.path.isdir(precomputed_ref_dir):
				raise FileNotFoundError(f"Precomputed reference alignment directory not found: {precomputed_ref_dir}")
			return precomputed_ref_dir

		if not self.update_db:
			return None

		db_ref_dir = os.path.join(self.base_dir, self.output_dir, "db_ref_backbones")
		resolved = self.export_update_backbones(db_ref_dir)
		if resolved:
			print(f"Using DB-derived precomputed segment alignments from {resolved}")
		return resolved

	def process_all_masters(self, master_list, nextalign_dir, master_segment_map=None, precomputed_ref_dir=None):
		# Resolve the run-wide projectable set once, from the same backbones the
		# loop below opens. Scoping the orphan report to a single segment made it
		# fire on every pass with a number three orders of magnitude too large.
		try:
			run_refs, _source = projectability.projectable_reference_ids(
				nextalign_dir=nextalign_dir,
				precomputed_ref_dir=precomputed_ref_dir,
				master_segment_map=master_segment_map,
			)
			self.run_projectable_refs = run_refs or None
		except Exception as exc:
			print(f"[warn] Could not resolve run-wide projectable references: {exc}")
			self.run_projectable_refs = None
		self._merged_already_deduped = False

		for master in master_list:
			ref_aln_file = None
			segment_val = (master_segment_map or {}).get(master)
			if precomputed_ref_dir and os.path.isdir(precomputed_ref_dir):
				precomputed_ref = self.find_precomputed_reference_alignment(precomputed_ref_dir, segment_val)
				if precomputed_ref:
					ref_aln_file = precomputed_ref
					print(
						f"Using precomputed segment alignment for master {master} "
						f"(segment {segment_val}): {ref_aln_file}"
					)
				else:
					print(
						f"[warn] No precomputed segment alignment found for master {master} "
						f"(segment {segment_val}). Falling back to nextalign reference output."
					)
					if self.strict_segment_backbone:
						raise ValueError(f"Missing DB/precomputed backbone for segment {segment_val} (master {master}) in strict update mode")

			if not ref_aln_file:
				ref_aln_file = join(nextalign_dir, "reference_aln", master, f"{master}.aligned.fasta")

			if os.path.exists(ref_aln_file):
				print(f"Processing master {master}...")
				self.process_master_alignment(ref_aln_file, self.input_dir, self.base_dir, self.output_dir, self.keep_intermediate_files, segment_val)
			else:
				print(f"Reference alignment for {master} not found at {ref_aln_file}")
				if self.strict_segment_backbone:
					raise FileNotFoundError(f"Reference alignment missing for master {master} segment {segment_val}: {ref_aln_file}")

	@staticmethod
	def _relabel_record(record, canonical_id):
		"""Re-label ``record`` with its canonical accession, annotation intact.

		Setting ``record.id`` on its own is not enough. Biopython writes the
		*description* as the FASTA title whenever the description's first token equals
		the id, and for a header with no spaces the description *is* the whole header -
		so a record read as ``NC_001542.1`` kept being written out versioned however
		the id was changed. Only the leading accession token is rewritten; any free
		text after it is annotation and is left exactly as it was.
		"""
		original = record.id
		record.id = canonical_id
		record.name = canonical_id
		description = record.description or ""
		if description == original:
			record.description = canonical_id
		elif description.startswith(original + " "):
			record.description = canonical_id + description[len(original):]
		return record

	@staticmethod
	def _build_projection_plan(ref_aligned_str, nextalign_ref_seq):
		"""Compile (is_slice, q_start, q_end, gap_str) runs once per reference backbone row."""
		master_dq_raw = ref_aligned_str.upper().replace('-', '')
		query_dq_raw = nextalign_ref_seq.upper().replace('-', '')

		n, m = len(master_dq_raw), len(query_dq_raw)
		width = max(100, abs(n - m) + 50)
		dp = [[-999999] * (2 * width + 1) for _ in range(n + 1)]
		dp[0][width] = 0
		for j in range(1, min(width + 1, m + 1)):
			dp[0][j + width] = -j
		for i in range(1, n + 1):
			for offset in range(2 * width + 1):
				j = i + offset - width
				if j < 0 or j > m:
					continue
				val1 = -999999
				if offset + 1 <= 2 * width:
					val1 = dp[i-1][offset+1] - 1
				val2 = -999999
				if offset - 1 >= 0:
					val2 = dp[i][offset-1] - 1
				val3 = -999999
				if j > 0:
					match_score = 1 if master_dq_raw[i-1] == query_dq_raw[j-1] else -1
					val3 = dp[i-1][offset] + match_score
				dp[i][offset] = max(val1, val2, val3)
		mapping = {x: None for x in range(n)}
		i = n
		offset = m - n + width
		while i > 0 and 0 <= offset <= 2 * width:
			j = i + offset - width
			val = dp[i][offset]
			if j > 0:
				match_score = 1 if master_dq_raw[i-1] == query_dq_raw[j-1] else -1
				if val == dp[i-1][offset] + match_score:
					if master_dq_raw[i-1] == query_dq_raw[j-1]:
						mapping[i-1] = j-1
					i -= 1
					continue
			if offset + 1 <= 2 * width and val == dp[i-1][offset+1] - 1:
				i -= 1
				offset += 1
				continue
			if offset - 1 >= 0 and val == dp[i][offset-1] - 1:
				offset -= 1
				continue
			break

		nextalign_raw_to_col = {}
		raw_idx = 0
		for col_idx, char in enumerate(nextalign_ref_seq):
			if char != '-':
				nextalign_raw_to_col[raw_idx] = col_idx
				raw_idx += 1

		col_targets = []
		master_raw_idx = 0
		max_col = -1
		for char in ref_aligned_str:
			if char == '-':
				col_targets.append(-1)
			else:
				query_raw_idx = mapping.get(master_raw_idx)
				col_nextalign = nextalign_raw_to_col.get(query_raw_idx) if query_raw_idx is not None else None
				if col_nextalign is not None:
					col_targets.append(col_nextalign)
					if col_nextalign > max_col:
						max_col = col_nextalign
				else:
					col_targets.append(-1)
				master_raw_idx += 1

		plan = []
		idx = 0
		total = len(col_targets)
		while idx < total:
			c = col_targets[idx]
			if c == -1:
				g_end = idx + 1
				while g_end < total and col_targets[g_end] == -1:
					g_end += 1
				plan.append((False, 0, 0, '-' * (g_end - idx)))
				idx = g_end
			else:
				s_start = c
				s_end = c + 1
				idx += 1
				while idx < total and col_targets[idx] == s_end:
					s_end += 1
					idx += 1
				plan.append((True, s_start, s_end, ''))
		return plan, col_targets, max_col

	@staticmethod
	def _project_seq_str(query_seq, plan, col_targets, max_col):
		if len(query_seq) > max_col:
			return ''.join(query_seq[s:e] if is_slice else gap_str for is_slice, s, e, gap_str in plan)
		qlen = len(query_seq)
		return ''.join(query_seq[c] if (c >= 0 and c < qlen) else '-' for c in col_targets)

	def insert_gaps(self, reference_aligned, subalignment_seqs, ref_id):
		ref_aligned_str = str(reference_aligned)
		ref_key = accession_utils.normalise_accession(ref_id)
		nextalign_ref_rec = None
		if ref_key is not None:
			for r in subalignment_seqs:
				if accession_utils.normalise_accession(r.id) == ref_key:
					nextalign_ref_rec = r
					break
		if not nextalign_ref_rec:
			ref_with_gaps_list = list(reference_aligned)
			updated_sequences = []
			for seq_record in subalignment_seqs:
				sequence = list(str(seq_record.seq))
				gapped_sequence = []
				seq_idx = 0
				for char in ref_with_gaps_list:
					if char == '-':
						gapped_sequence.append('-')
					else:
						if seq_idx < len(sequence):
							gapped_sequence.append(sequence[seq_idx])
							seq_idx += 1
						else:
							gapped_sequence.append('-')
				gapped_seq_str = ''.join(gapped_sequence)
				seq_record.seq = Seq(gapped_seq_str)
				updated_sequences.append(seq_record)
			return updated_sequences

		plan, col_targets, max_col = self._build_projection_plan(ref_aligned_str, str(nextalign_ref_rec.seq))
		updated_sequences = []
		for seq_record in subalignment_seqs:
			seq_record.seq = Seq(self._project_seq_str(str(seq_record.seq), plan, col_targets, max_col))
			updated_sequences.append(seq_record)
		return updated_sequences

	@staticmethod
	def _read_fasta_headers(fasta_path):
		"""Record ids of ``fasta_path``, canonicalised, in file order.

		The result is compared against backbone ids and against the ``query_aln/<ref>/``
		directory name, both of which are canonical, so the headers have to be too.
		``split('.')[0]`` got that right only for the versioned case: it also truncated
		every label with a dot in it, which turned one orphaned query into a different,
		non-existent accession in the warning the operator has to act on.
		"""
		if not os.path.exists(fasta_path):
			return []
		return accession_utils.normalise_accession_series(
			record.id for record in SeqIO.parse(fasta_path, "fasta")
		)

	def find_orphan_query_references(self, reference_alignment_file, input_dir):
		"""
		Find query_aln reference directories that cannot be projected because their
		reference is absent from *every* backbone this run opens.

		The comparison is run-wide, not per-segment. It used to be per-segment, and
		because this runs once per master that made it meaningless: on the segment-1
		pass every segment-2..8 reference is trivially "absent", so the warning
		reported 282-370 references and ~3.4M query sequences on all eight passes.
		A reference is only a genuine orphan if no segment's backbone contains it.

		``self.run_projectable_refs`` is populated by :meth:`process_all_masters`.
		When this object is driven directly with ``-r`` (single backbone, no master
		list) it stays None and the old single-file comparison is the right one.

		Returns a dict:
			{ref_id: [query_accession_1, query_accession_2, ...]}
		"""
		projectable_refs = getattr(self, "run_projectable_refs", None)
		if not projectable_refs:
			master_alignment = SeqIO.to_dict(SeqIO.parse(reference_alignment_file, "fasta"))
			projectable_refs = {
				accession_utils.normalise_accession(ref_id)
				for ref_id in master_alignment.keys()
			}

		orphans = {}
		if not os.path.isdir(input_dir):
			return orphans

		for dir_name in os.listdir(input_dir):
			ref_path = os.path.join(input_dir, dir_name)
			if not os.path.isdir(ref_path):
				continue
			# The directory is named by NextalignAlignment from
			# ``accession_from_filename``, so its canonical form is the key the backbone
			# ids and the record headers are held under; the name on disk still builds
			# the path.
			ref_id = accession_utils.accession_from_filename(dir_name)
			if ref_id in projectable_refs:
				continue
			aln_file = os.path.join(ref_path, f"{dir_name}.aligned.fasta")
			query_ids = [q for q in self._read_fasta_headers(aln_file) if q != ref_id]
			if query_ids:
				orphans[ref_id] = query_ids

		return orphans

	@staticmethod
	def _iter_fasta_raw(fasta_path, chunk_bytes=16 * 1024 * 1024):
		"""Yield (id_token, full_header_without_gt, seq_str) streaming in 16 MiB blocks."""
		with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
			carry = ""
			first_chunk = True
			while True:
				chunk = fh.read(chunk_bytes)
				if not chunk:
					break
				data = carry + chunk
				if first_chunk:
					gt = data.find(">")
					if gt == -1:
						carry = ""
						continue
					data = data[gt + 1:]
					first_chunk = False
				parts = data.split("\n>")
				carry = parts[-1]
				for blk in parts[:-1]:
					hdr, sep_nl, body = blk.partition("\n")
					if sep_nl or hdr:
						h = hdr.strip()
						rid = h.split(None, 1)[0] if h else ""
						yield rid, h, body.replace("\n", "").replace("\r", "")
			if carry and not first_chunk:
				hdr, sep_nl, body = carry.partition("\n")
				if sep_nl or hdr:
					h = hdr.strip()
					rid = h.split(None, 1)[0] if h else ""
					yield rid, h, body.replace("\n", "").replace("\r", "")

	@staticmethod
	def _relabel_header_str(rid, full_header, canonical_id):
		if not full_header or full_header == rid:
			return canonical_id
		if full_header.startswith(rid + " "):
			return canonical_id + full_header[len(rid):]
		return canonical_id

	def process_master_alignment(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False, segment_value=None):
		master_alignment = SeqIO.to_dict(SeqIO.parse(reference_alignment_file, "fasta"))
		merged_output_file = os.path.join(base_dir, output_dir, os.path.basename(reference_alignment_file).replace(".fasta", "_merged_MSA.fasta"))
		os.makedirs(join(base_dir, output_dir), exist_ok=True)
		os.makedirs(join(output_dir), exist_ok=True)

		seg_norm = self._normalize_segment(segment_value) or "0"
		seen_merged_keys = set()
		wrote_any_merged = False
		manifest_direct_file = None
		manifest_direct_handle = None
		if getattr(self, "_manifest_out_path", None):
			manifest_dir = os.path.dirname(self._manifest_out_path)
			if manifest_dir:
				os.makedirs(manifest_dir, exist_ok=True)
			need_hdr = not os.path.exists(self._manifest_out_path) or os.path.getsize(self._manifest_out_path) == 0
			manifest_direct_handle = open(self._manifest_out_path, "a", encoding="utf-8")
			if need_hdr:
				manifest_direct_handle.write("accession\tsegment\toutput_file\n")
			manifest_direct_file = self._manifest_out_path

		def _emit_to_merged(mh, rec_id, header_line, seq_str):
			nonlocal wrote_any_merged
			wrote_any_merged = True
			dedup_key = rec_id.split('|')[0] if '|' in rec_id else rec_id
			if dedup_key not in seen_merged_keys:
				seen_merged_keys.add(dedup_key)
				mh.write(f">{header_line}\n{seq_str}\n")
			if manifest_direct_handle is not None:
				manifest_direct_handle.write(f"{rec_id}\t{seg_norm}\t{merged_output_file}\n")
			elif len(self.segment_manifest_rows) < 250000:
				self.segment_manifest_rows.append({
					"accession": rec_id,
					"segment": seg_norm,
					"output_file": merged_output_file,
				})

		try:
			with open(merged_output_file, "w", encoding="utf-8") as merged_handle:
				for record_id, ref_record in master_alignment.items():
					ref_aligned_str = str(ref_record.seq)
					ref_id = accession_utils.normalise_accession(record_id) or record_id
					self._relabel_record(ref_record, ref_id)
					ref_header_line = ref_record.description if ref_record.description else ref_id
					subalignment_file = os.path.join(input_dir, f"{ref_id}/{ref_id}.aligned.fasta")
					output_file = os.path.join(output_dir, f"{ref_id}_aligned_padded.fasta")

					if os.path.exists(subalignment_file):
						print(f"Processing subalignment for {ref_id} using {subalignment_file}")
						ref_key = accession_utils.normalise_accession(ref_id)
						nextalign_ref_seq = None
						if ref_key is not None:
							for rid, _, rseq in self._iter_fasta_raw(subalignment_file):
								if accession_utils.normalise_accession(rid) == ref_key:
									nextalign_ref_seq = rseq
									break

						plan_tuple = self._build_projection_plan(ref_aligned_str, nextalign_ref_seq) if nextalign_ref_seq is not None else None
						ref_chars_fallback = list(ref_aligned_str) if plan_tuple is None else None
						skipped = 0

						sub_handle = open(output_file, "w", encoding="utf-8") if keep_intermediate_files else None
						try:
							if sub_handle is not None:
								sub_handle.write(f">{ref_header_line}\n{ref_aligned_str}\n")
							_emit_to_merged(merged_handle, ref_id, ref_header_line, ref_aligned_str)

							for rid, full_hdr, qseq in self._iter_fasta_raw(subalignment_file):
								norm_rid = accession_utils.normalise_accession(rid)
								if self.skip_ids and norm_rid in self.skip_ids:
									skipped += 1
									continue
								if plan_tuple is not None:
									gapped_str = self._project_seq_str(qseq, plan_tuple[0], plan_tuple[1], plan_tuple[2])
								else:
									q_idx = 0
									qlen = len(qseq)
									out_chars = []
									for ch in ref_chars_fallback:
										if ch == '-':
											out_chars.append('-')
										elif q_idx < qlen:
											out_chars.append(qseq[q_idx])
											q_idx += 1
										else:
											out_chars.append('-')
									gapped_str = ''.join(out_chars)
								if sub_handle is not None:
									sub_handle.write(f">{full_hdr}\n{gapped_str}\n")
								_emit_to_merged(merged_handle, rid, full_hdr, gapped_str)
						finally:
							if sub_handle is not None:
								sub_handle.close()

						if skipped:
							print(f"[skip_ids] Skipped {skipped} sequence(s) from {ref_id} subalignment")
						if keep_intermediate_files:
							print(f"Saved updated alignment to {output_file}")
					else:
						print(f"Subalignment file {subalignment_file} not found. Adding reference only for {ref_id}.")
						if keep_intermediate_files:
							with open(output_file, "w", encoding="utf-8") as output_handle:
								output_handle.write(f">{ref_header_line}\n{ref_aligned_str}\n")
							print(f"Saved reference-only alignment to {output_file}")
						_emit_to_merged(merged_handle, ref_id, ref_header_line, ref_aligned_str)
		finally:
			if manifest_direct_handle is not None:
				manifest_direct_handle.close()

		if wrote_any_merged:
			self._merged_already_deduped = True
			print(f"Saved merged alignment to {merged_output_file}")
		else:
			if os.path.exists(merged_output_file):
				os.remove(merged_output_file)
			print(f"No sequences found for {reference_alignment_file}, skipping merged file creation.")

		orphan_refs = self.find_orphan_query_references(reference_alignment_file, input_dir)
		if orphan_refs:
			orphan_query_count = sum(len(v) for v in orphan_refs.values())
			preview_refs = ', '.join(sorted(orphan_refs.keys())[:10])
			print(
				f"[warn] {len(orphan_refs)} query_aln reference(s) are not present in master-projected alignment "
				f"and were skipped ({orphan_query_count} query sequence(s))."
			)
			print(f"[warn] Skipped reference examples: {preview_refs}")

		if not keep_intermediate_files:
			for record_id in master_alignment:
				ref_id = accession_utils.normalise_accession(record_id) or record_id
				padded_file = os.path.join(output_dir, f"{ref_id}_aligned_padded.fasta")
				if os.path.exists(padded_file):
					os.remove(padded_file)
					print(f"Deleted intermediate file {padded_file}")
			if os.path.abspath(output_dir) != os.path.abspath(os.path.join(base_dir, output_dir)):
				shutil.rmtree(output_dir, ignore_errors=True)

	def find_fasta_file(self, input_dir,new_outputfile=False):
		directory = join(self.base_dir, self.output_dir)
		if not os.path.exists(directory):
			return None
		if new_outputfile:
			return os.path.join(directory, "new_output.fasta")
		else:
			for file in os.listdir(directory):
				if file.endswith(".fasta") or file.endswith(".fa"):
					return os.path.join(directory, file)
			return None

	def remove_redundant_sequences(self):
		if getattr(self, "_merged_already_deduped", False):
			return
		input_file = self.find_fasta_file(join(self.base_dir, self.output_dir), self.new_outputfile)
		if not input_file:
			print(f"No fasta file found in {join(self.base_dir, self.output_dir)} to remove redundant sequences from.")
			return

		seen = set()
		tmp_file = input_file + ".dedup.tmp"
		try:
			with open(tmp_file, "w", encoding="utf-8") as out_f:
				for rid, full_hdr, seq_str in self._iter_fasta_raw(input_file):
					accession = rid.split('|')[0] if '|' in rid else rid
					if accession not in seen:
						seen.add(accession)
						out_f.write(f">{full_hdr}\n{seq_str}\n")
			os.replace(tmp_file, input_file)
		except Exception as e:
			if os.path.exists(tmp_file):
				os.remove(tmp_file)
			print(f"Error removing redundant sequences: {e}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Insert gaps from master alignment into corresponding subalignments.")
	parser.add_argument("-r", "--reference_alignment", help="Path to master alignment file (FASTA format).")
	parser.add_argument("-i", "--input_dir", help="Directory containing subalignment files (Nextalign output).", default="tmp/Nextalign/query_aln")
	parser.add_argument("-d", "--base_dir", help="Base directory.", default="tmp")
	parser.add_argument("-o", "--output_dir", help="Directory to save padded subalignments and merged files.", default="Pad-alignment")
	parser.add_argument("--keep_intermediate_files", action="store_true", help="Keep intermediate files (padded subalignment). Default: disabled (files will be removed).")
	parser.add_argument("-n","--new_outputfile", action="store_true", help="New output file name for the final merged alignment.")
	parser.add_argument("-m", "--master_acc", help="Path to ref_list file (TSV with columns: accession, type, segment) OR comma-separated master accession IDs. For segmented viruses, the script extracts all 'master' entries to process each segment separately.")
	parser.add_argument("-nd", "--nextalign_dir", help="Path to Nextalign output directory containing reference_aln/ and query_aln/ subdirectories.")
	parser.add_argument("--precomputed_ref_dir", default=None, help="Optional directory containing precomputed segment alignments (e.g. refset_<segment>_aln.fasta). If absent or unmatched, falls back to nextalign reference_aln outputs.")
	parser.add_argument("--segment_manifest_out", default=None, help="Optional TSV path for accession-to-segment-to-output mapping")
	parser.add_argument("--strict_segment_backbone", action="store_true", help="Fail if segment-specific precomputed backbone is missing")
	parser.add_argument("--update_db", default=None, help="Existing update DB; when set, derive per-segment backbone alignments from sequence_alignment")
	parser.add_argument("--skip_ids", default=None, help="Path to a file listing accession IDs (one per line) to exclude from subalignments (e.g. filtered_sequences_ids.txt from CollectFilteredSequences)")
 
	args = parser.parse_args()

	processor = PadAlignment(
		args.reference_alignment,
		args.input_dir,
		args.base_dir,
		args.output_dir,
		args.keep_intermediate_files,
		args.new_outputfile,
		args.segment_manifest_out,
		args.strict_segment_backbone or bool(PadAlignment._normalize_optional_path(args.update_db)),
		args.update_db,
		args.skip_ids,
	)

	if args.segment_manifest_out and os.path.exists(args.segment_manifest_out):
		os.remove(args.segment_manifest_out)
	processor._manifest_out_path = args.segment_manifest_out

	if args.master_acc and args.nextalign_dir:
		masters = processor.get_master_list(args.master_acc)
		master_segment_map = processor.get_master_segment_map(args.master_acc)
		precomputed_ref_dir = processor.resolve_precomputed_ref_dir(args.precomputed_ref_dir)
		processor.process_all_masters(
			masters,
			args.nextalign_dir,
			master_segment_map=master_segment_map,
			precomputed_ref_dir=precomputed_ref_dir,
		)
	elif args.reference_alignment:
		processor.process_master_alignment(
			reference_alignment_file=args.reference_alignment,
			input_dir=args.input_dir,
			base_dir=args.base_dir,
			output_dir=args.output_dir,
			keep_intermediate_files=args.keep_intermediate_files
		)
	else:
		print("Error: Either -r (reference alignment) or both -m (master acc) and -nd (nextalign dir) must be provided.")

	processor.remove_redundant_sequences()

	if args.segment_manifest_out:
		if os.path.exists(args.segment_manifest_out):
			seen_m = set()
			tmp_m = args.segment_manifest_out + ".tmp"
			with open(args.segment_manifest_out, "r", encoding="utf-8") as mf_in, \
			     open(tmp_m, "w", encoding="utf-8") as mf_out:
				hdr = mf_in.readline()
				if hdr:
					mf_out.write(hdr)
				for line in mf_in:
					if line not in seen_m:
						seen_m.add(line)
						mf_out.write(line)
			os.replace(tmp_m, args.segment_manifest_out)
		else:
			manifest_dir = os.path.dirname(args.segment_manifest_out)
			if manifest_dir:
				os.makedirs(manifest_dir, exist_ok=True)
			df_manifest = pd.DataFrame(processor.segment_manifest_rows).drop_duplicates()
			df_manifest.to_csv(args.segment_manifest_out, sep='\t', index=False)

# -----------------------------------------------------------------------------
# UPDATE-MODE PLAN FOR PADALIGNMENT (COMMENT ONLY)
# Goal: project incoming partial/new sequences onto an existing DB-consistent
# alignment backbone, while preserving historical insertion structure.
#
# 1) Inputs required for robust update mode:
#    - Existing DB-derived per-segment backbone alignments (or equivalent files
#      exported from sequence_alignment) as primary projection targets.
#    - Incoming nextalign query_aln outputs for newly fetched/changed accessions.
#    - Master list + segment mapping from ref_list/metadata.
#
# 2) Segment-first processing contract:
#    - Partition all operations by segment before padding.
#    - Build one work unit per segment: {segment_key, masters, references,
#      incoming query alignments, existing backbone alignment}.
#    - Unsegmented viruses run as a single segment bucket.
#    - Never allow cross-segment projection (segment N query -> segment M
#      backbone is invalid and should be rejected).
#
# 3) Backbone selection order per segment:
#    - Preferred: existing DB backbone for that segment (preserves historic
#      insertion columns and coordinate compatibility).
#    - Fallback 1: precomputed segment reference alignment
#      (e.g., refset_<segment>_aln.fasta).
#    - Fallback 2: nextalign reference_aln output.
#    - If no segment-specific backbone exists, fail that segment explicitly
#      (do not silently reuse another segment's backbone).
#
# 4) Gap projection behavior per segment:
#    - For each reference subalignment in the segment, insert gaps according to
#      the segment backbone reference alignment.
#    - Preserve existing columns from historic backbone exactly.
#    - If incoming data introduces genuinely novel insertion columns, append them
#      in a deterministic segment-local manner (stable ordering), without
#      removing historical columns.
#
# 5) Merging outputs:
#    - Produce segment-scoped merged MSA outputs first.
#    - Optional global merged output is a concatenation of segment outputs only
#      when downstream expects a combined file; otherwise keep per-segment files
#      as canonical update artifacts.
#    - Keep mapping manifest: accession -> segment -> output file.
#
# 6) Orphans and unresolved references:
#    - Keep current orphan detection, but report per segment.
#    - Distinguish:
#      * orphan due to missing reference in segment backbone
#      * orphan due to unknown/missing segment assignment
#    - Route unresolved records to exclusion/quarantine list for DB audit tables.
#
# 7) Idempotency and determinism requirements:
#    - Rerunning same update batch with same inputs yields byte-stable segment
#      outputs (or equivalent sequence+coordinate content).
#    - Deduplication should be segment-aware if accession reuse across segments
#      is possible (key by accession+segment where required).
#
# 8) Hand-off to CalcAlignmentCord and DB update:
#    - Emit explicit segment metadata alongside padded outputs (filename
#      convention or manifest TSV) so coordinate calculation selects the correct
#      master GFF per segment.
#    - Ensure only update-scope accessions are forwarded for feature recalculation
#      and DB upsert.
# -----------------------------------------------------------------------------

