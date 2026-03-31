import os
import shutil
import argparse
import re
import sqlite3
import pandas as pd
from os.path import join
from Bio import SeqIO
from Bio.Seq import Seq
from ExportRefListFromUpdateDb import load_master_accessions, load_reference_rows

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

	@staticmethod
	def _load_skip_ids(skip_ids_arg) -> set:
		"""Load a set of accession IDs to skip from a file path or return empty set."""
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
				acc = line.strip().split()[0] if line.strip() else ''
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
		if self.update_db:
			return load_master_accessions(self.update_db)
		if os.path.isfile(master_acc):
			try:
				df = pd.read_csv(master_acc, sep='\t', header=None, dtype=str)
				if df.shape[1] >= 2:
					# Filter for master if available
					if df[1].str.lower().eq('master').any():
						masters = df[df[1].str.lower() == 'master']
						return masters[0].tolist()
						
					# If no explicit master, exclude exclusion_list entries
					df = df[df[1].str.lower() != 'exclusion_list']
						
				return df[0].tolist()
			except:
				return []
		else:
			return [x.strip() for x in master_acc.split(',') if x.strip()]

	def get_master_segment_map(self, master_acc):
		if self.update_db:
			master_segment = {}
			refs = load_reference_rows(self.update_db)
			masters = refs[refs["accession_type"].fillna("").astype(str).str.lower() == "master"]
			if masters.empty:
				masters = refs
			for _, row in masters.iterrows():
				acc = str(row["primary_accession"]).strip()
				segment = self._normalize_segment(row.get("segment")) or "0"
				if acc:
					master_segment[acc] = segment
			return master_segment

		master_segment = {}
		if not os.path.isfile(master_acc):
			return master_segment
		try:
			df = pd.read_csv(master_acc, sep='\t', header=None, dtype=str)
		except Exception:
			return master_segment

		if df.shape[1] < 3:
			return master_segment

		df[0] = df[0].astype(str).str.strip()
		df[1] = df[1].astype(str).str.strip().str.lower()
		df[2] = df[2].astype(str).str.strip()

		masters = df[df[1] == 'master']
		for _, row in masters.iterrows():
			master_segment[row[0]] = row[2]
		return master_segment

	@staticmethod
	def _normalize_segment(segment_value):
		if segment_value is None:
			return None
		segment_str = str(segment_value).strip()
		if not segment_str:
			return None
		match = re.search(r"(\d+)", segment_str)
		if not match:
			return None
		return str(int(match.group(1)))

	def find_precomputed_reference_alignment(self, precomputed_ref_dir, segment_value):
		"""
		Find a precomputed segment alignment file in precomputed_ref_dir.
		Expected canonical naming is refset_<segment>_aln.fasta, with a
		numeric-segment fallback for custom naming.
		"""
		if not precomputed_ref_dir or not os.path.isdir(precomputed_ref_dir):
			return None

		segment = self._normalize_segment(segment_value)
		if not segment:
			fasta_files = sorted(
				[
					os.path.join(precomputed_ref_dir, fname)
					for fname in os.listdir(precomputed_ref_dir)
					if fname.lower().endswith((".fasta", ".fa"))
				]
			)
			if len(fasta_files) == 1:
				print(
					f"[warn] No segment value supplied; using sole precomputed alignment {os.path.basename(fasta_files[0])}."
				)
				return fasta_files[0]
			return None

		preferred = os.path.join(precomputed_ref_dir, f"refset_{segment}_aln.fasta")
		if os.path.exists(preferred):
			return preferred

		fallback_matches = []
		for fname in os.listdir(precomputed_ref_dir):
			if not fname.lower().endswith((".fasta", ".fa")):
				continue
			seg_match = re.search(r"(\d+)", fname)
			if seg_match and seg_match.group(1) == segment:
				fallback_matches.append(os.path.join(precomputed_ref_dir, fname))

		if len(fallback_matches) == 1:
			return fallback_matches[0]
		if len(fallback_matches) > 1:
			print(
				f"[warn] Multiple precomputed alignment files matched segment {segment}: "
				f"{', '.join(os.path.basename(x) for x in sorted(fallback_matches))}. "
				f"Using {os.path.basename(sorted(fallback_matches)[0])}."
			)
			return sorted(fallback_matches)[0]

		return None

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

	def insert_gaps(self, reference_aligned, subalignment_seqs):
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

	@staticmethod
	def _read_fasta_headers(fasta_path):
		headers = []
		if not os.path.exists(fasta_path):
			return headers
		for record in SeqIO.parse(fasta_path, "fasta"):
			headers.append(record.id.split('.')[0])
		return headers

	def find_orphan_query_references(self, reference_alignment_file, input_dir):
		"""
		Find query_aln reference directories that cannot be projected because their
		reference is absent from the master-projected reference alignment file.

		Returns a dict:
			{ref_id: [query_accession_1, query_accession_2, ...]}
		"""
		master_alignment = SeqIO.to_dict(SeqIO.parse(reference_alignment_file, "fasta"))
		projectable_refs = {ref_id.split('.')[0] for ref_id in master_alignment.keys()}

		orphans = {}
		if not os.path.isdir(input_dir):
			return orphans

		for ref_id in os.listdir(input_dir):
			ref_path = os.path.join(input_dir, ref_id)
			if not os.path.isdir(ref_path):
				continue
			if ref_id in projectable_refs:
				continue
			aln_file = os.path.join(ref_path, f"{ref_id}.aligned.fasta")
			query_ids = [q for q in self._read_fasta_headers(aln_file) if q != ref_id]
			if query_ids:
				orphans[ref_id] = query_ids

		return orphans

	def process_master_alignment(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False, segment_value=None):
		master_alignment = SeqIO.to_dict(SeqIO.parse(reference_alignment_file, "fasta"))
		merged_sequences = []
		for ref_id, ref_record in master_alignment.items():
			ref_aligned = ref_record.seq
			ref_id = ref_id.split('.')[0]
			subalignment_file = os.path.join(input_dir, f"{ref_id}/{ref_id}.aligned.fasta")
			if os.path.exists(subalignment_file):
				print(f"Processing subalignment for {ref_id} using {subalignment_file}")
				subalignment_seqs = list(SeqIO.parse(subalignment_file, "fasta"))
				if self.skip_ids:
					before = len(subalignment_seqs)
					subalignment_seqs = [r for r in subalignment_seqs if r.id.split('.')[0] not in self.skip_ids]
					skipped = before - len(subalignment_seqs)
					if skipped:
						print(f"[skip_ids] Skipped {skipped} sequence(s) from {ref_id} subalignment")
				updated_seqs = self.insert_gaps(ref_aligned, subalignment_seqs)
				
				# Add the reference sequence to the list of sequences
				updated_seqs.insert(0, ref_record)

				os.makedirs(join(output_dir), exist_ok=True)
				output_file = os.path.join(output_dir, f"{ref_id}_aligned_padded.fasta")
				with open(join(output_file), "w") as output_handle:
					SeqIO.write(updated_seqs, output_handle, "fasta")
					print(f"Saved updated alignment to {output_file}")
				merged_sequences.extend(updated_seqs)
			else:
				print(f"Subalignment file {subalignment_file} not found. Adding reference only for {ref_id}.")
				# Even if no subalignment exists (no queries hit this ref), we must include the ref itself
				# so it appears in the final merged output.
				updated_seqs = [ref_record]
				
				os.makedirs(join(output_dir), exist_ok=True)
				output_file = os.path.join(output_dir, f"{ref_id}_aligned_padded.fasta")
				with open(join(output_file), "w") as output_handle:
					SeqIO.write(updated_seqs, output_handle, "fasta")
					print(f"Saved reference-only alignment to {output_file}")
				merged_sequences.extend(updated_seqs)

		if merged_sequences:
			merged_output_file = os.path.join(base_dir, output_dir, os.path.basename(reference_alignment_file).replace(".fasta", "_merged_MSA.fasta"))
			os.makedirs(join(base_dir, output_dir), exist_ok=True)
			with open(merged_output_file, "w") as merged_output_handle:
				SeqIO.write(merged_sequences, merged_output_handle, "fasta")
				print(f"Saved merged alignment to {merged_output_file}")
			for rec in merged_sequences:
				self.segment_manifest_rows.append({
					"accession": rec.id,
					"segment": self._normalize_segment(segment_value) or "0",
					"output_file": merged_output_file,
				})
		else:
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
			for ref_id in master_alignment:
				padded_file = os.path.join(output_dir, f"{ref_id.split('.')[0]}_aligned_padded.fasta")
				if os.path.exists(padded_file):
					os.remove(padded_file)
					print(f"Deleted intermediate file {padded_file}")
			shutil.rmtree(output_dir)

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
		
		input_file = self.find_fasta_file(join(self.base_dir, self.output_dir),self.new_outputfile) 
		if not input_file:
			print(f"No fasta file found in {join(self.base_dir, self.output_dir)} to remove redundant sequences from.")
			return

		unique_records = {}
		try:
			for record in SeqIO.parse(input_file, "fasta"):
				accession = record.id.split('|')[0] if '|' in record.id else record.id
				if accession not in unique_records:
					unique_records[accession] = record
			with open(input_file, "w") as output_handle:
				SeqIO.write(unique_records.values(), output_handle, "fasta")
		except Exception as e:
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

