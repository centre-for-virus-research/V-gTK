import os
import re
import sys
import sqlite3
import pandas as pd
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
from GffToDictionary import GffDictionary
from CalcGenomeCords import CalculateGenomeCoordinates 
from ExportRefListFromUpdateDb import load_master_accessions_from_file, load_reference_file_table

class CalculateAlignmentCoordinates:

	def __init__(self, paded_alignment, master_gff, tmp_dir, output_dir, output_file, master_accession, blast_uniq_hits, update_db=None, update_scope_tsv=None, segment_map_tsv=None):
		self.paded_alignment = paded_alignment
		self.master_gff = master_gff
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.output_file = output_file
		self.master_accession = master_accession
		self.blast_uniq_hits = blast_uniq_hits
		self.update_db = update_db
		self.update_scope_tsv = update_scope_tsv
		self.segment_map_tsv = segment_map_tsv

	@staticmethod
	def _normalize_segment_value(value):
		if value is None:
			return ""
		text = str(value).strip()
		if not text or text.lower() == "nan":
			return ""
		digits = ''.join(ch for ch in text if ch.isdigit())
		return digits if digits else text

	@staticmethod
	def _infer_segment_from_alignment_name(fasta_file):
		name = os.path.basename(fasta_file)
		patterns = [
			r"(?:^|[_-])refset[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])segment[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])seg[_-]?(\d+)(?:$|[_-])",
		]
		for pattern in patterns:
			match = re.search(pattern, name, flags=re.IGNORECASE)
			if match:
				return match.group(1)
		return ""

	def load_existing_feature_accessions(self):
		if not self.update_db:
			return set()
		if not os.path.isfile(self.update_db):
			raise FileNotFoundError(f"Update DB not found: {self.update_db}")
		conn = sqlite3.connect(self.update_db)
		try:
			row = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='features'").fetchone()
			if row is None:
				return set()
			cols = [r[1] for r in conn.execute("PRAGMA table_info(features)").fetchall()]
			if "accession" not in cols:
				raise ValueError("Update DB features table is missing required column: accession")
			df = pd.read_sql_query("SELECT accession FROM features WHERE accession IS NOT NULL", conn)
			return set(df["accession"].astype(str).str.strip().tolist())
		finally:
			conn.close()

	def load_update_scope_accessions(self):
		if not self.update_scope_tsv or not os.path.isfile(self.update_scope_tsv):
			return set()
		df = pd.read_csv(self.update_scope_tsv, sep='\t', dtype=str)
		if 'primary_accession' not in df.columns:
			raise ValueError(f"Update scope TSV is missing required column: primary_accession ({self.update_scope_tsv})")
		return set(df['primary_accession'].fillna('').astype(str).str.strip().tolist())

	def load_segment_map(self):
		if not self.segment_map_tsv or not os.path.isfile(self.segment_map_tsv):
			return {}
		df = pd.read_csv(self.segment_map_tsv, sep='\t', dtype=str)
		if 'primary_accession' not in df.columns or 'segment' not in df.columns:
			raise ValueError(f"Segment map TSV is missing required columns: primary_accession, segment ({self.segment_map_tsv})")
		df['primary_accession'] = df['primary_accession'].fillna('').astype(str).str.strip()
		df['segment'] = df['segment'].fillna('').astype(str).str.strip()
		return dict(zip(df['primary_accession'], df['segment']))

	def get_master_list(self):
		if os.path.isfile(self.master_accession):
			try:
				return load_master_accessions_from_file(self.master_accession)
			except:
				return []
		else:
			return [x.strip() for x in self.master_accession.split(',') if x.strip()]

	def load_master_segment_map(self):
		if not os.path.isfile(self.master_accession):
			return {}
		try:
			ref_df = load_reference_file_table(self.master_accession)
		except Exception:
			return {}
		if ref_df.empty or 'accession_type' not in ref_df.columns:
			return {}
		masters = ref_df[
			ref_df['accession_type'].fillna('').astype(str).str.strip().str.lower() == 'master'
		].copy()
		if masters.empty or 'segment' not in masters.columns:
			return {}
		masters['segment'] = masters['segment'].map(self._normalize_segment_value)
		masters['primary_accession'] = masters['primary_accession'].fillna('').astype(str).str.strip()
		masters = masters[(masters['segment'] != '') & (masters['primary_accession'] != '')]
		masters = masters.drop_duplicates(subset=['segment'], keep='first')
		return dict(zip(masters['segment'], masters['primary_accession']))

	def resolve_master_for_alignment(self, fasta_file, masters, master_segment_map=None):
		for master in masters:
			if fasta_file.startswith(master):
				return master

		segment = self._infer_segment_from_alignment_name(fasta_file)
		if segment and master_segment_map:
			master = master_segment_map.get(segment)
			if master:
				return master

		if len(masters) == 1:
			return masters[0]

		return None

	def get_gff_for_master(self, master):
		if isinstance(self.master_gff, list):
			for gff in self.master_gff:
				if master in os.path.basename(gff):
					return gff
		elif isinstance(self.master_gff, str):
			if master in os.path.basename(self.master_gff):
				return self.master_gff
		return None

	def get_gap_ranges(self, sequence):
		gap_ranges = []
		start = None

		for i, char in enumerate(sequence):
			if char == '-':
				if start is None:
					start = i + 1
			else:
				if start is not None:
					gap_ranges.append([start, i])
					start = None

		if start is not None:
			gap_ranges.append([start, len(sequence)])

		return gap_ranges

	def count_gaps_before_position(self, gap_ranges, position):
		count = 0
		for start, end in gap_ranges:
			if end < position:
				count += (end - start + 1)
			elif start <= position <= end:
				count += (position - start + 1)
		return count

	def format_genome_coverage(self, query_alignment, master_coord_to_aln_pos, feature_start, feature_end):
		try:
			feature_start = int(feature_start)
			feature_end = int(feature_end)
		except (TypeError, ValueError):
			return "NA"

		if feature_start > feature_end:
			feature_start, feature_end = feature_end, feature_start

		aln_indexes = []
		for master_pos in range(feature_start, feature_end + 1):
			if master_pos in master_coord_to_aln_pos:
				# Convert 1-based alignment position to 0-based string index
				aln_indexes.append(master_coord_to_aln_pos[master_pos] - 1)

		if not aln_indexes:
			return "NA"

		covered = 0
		for aln_index in aln_indexes:
			if aln_index >= len(query_alignment):
				continue

			base = query_alignment[aln_index]
			if base != '-' and base.upper() != 'N':
				covered += 1

		coverage = (covered / len(aln_indexes)) * 100
		return f"{coverage:.2f}"

	def recalculate_cds_coordinates_with_span(self, sequence_id, gap_ranges, cds_list, start_offset, genome_cord_start=None, genome_cord_end=None, master_coord_to_aln_pos=None):
		adjusted_coords = []
		clamp_to_span = genome_cord_start not in (None, "NA") and genome_cord_end not in (None, "NA")
		span_start = None
		span_end = None
		if clamp_to_span:
			span_start = int(str(genome_cord_start))
			span_end = int(str(genome_cord_end))

		for cds in cds_list:
			cds_start = int(cds['start'])
			cds_end = int(cds['end'])
			if clamp_to_span:
				assert span_start is not None and span_end is not None
				overlap_start = max(cds_start, span_start)
				overlap_end = min(cds_end, span_end)
				if overlap_start > overlap_end:
					continue
			else:
				overlap_start = cds_start
				overlap_end = cds_end

			gaps_before_start = self.count_gaps_before_position(gap_ranges, overlap_start)
			gaps_before_end = self.count_gaps_before_position(gap_ranges, overlap_end)

			adj_start = overlap_start - gaps_before_start
			adj_end = overlap_end - gaps_before_end

			if master_coord_to_aln_pos is not None:
				aln_start = master_coord_to_aln_pos.get(overlap_start, overlap_start)
				aln_end = master_coord_to_aln_pos.get(overlap_end, overlap_end)
			else:
				aln_start = overlap_start
				aln_end = overlap_end

			adjusted_entry = {
				'start': aln_start,
				'end': aln_end,
				'og_start': adj_start,
				'og_end': adj_end,
				'feature_start': cds_start,
				'feature_end': cds_end,
				'product': cds['product'],
			}
			if adjusted_entry not in adjusted_coords:
				adjusted_coords.append(adjusted_entry)

		return adjusted_coords

	def load_blast_hits(self):
		if not self.blast_uniq_hits or not os.path.isfile(self.blast_uniq_hits):
			raise FileNotFoundError(f"BLAST unique hits file not found: {self.blast_uniq_hits}")
		acc_dict = {}
		for i in open(self.blast_uniq_hits):
			parts = i.strip().split('\t')
			if len(parts) != 4:
				raise ValueError(f"Malformed BLAST hits row in {self.blast_uniq_hits}: {i.strip()}")
			query, ref, score, strand = parts
			acc_dict[query] = ref
		return acc_dict

	def load_historical_alignment_lengths(self):
		if not self.update_db or not os.path.isfile(self.update_db):
			return {}
		conn = sqlite3.connect(self.update_db)
		try:
			row = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='sequence_alignment'").fetchone()
			if row is None:
				return {}
			
			cols = [r[1] for r in conn.execute("PRAGMA table_info(sequence_alignment)").fetchall()]
			aln_col = None
			for candidate in ["alignment", "aligned_seq", "sequence", "aln", "alignment_seq"]:
				if candidate in cols:
					aln_col = candidate
					break
			
			if not aln_col:
				return {}
				
			id_col = "primary_accession" if "primary_accession" in cols else "sequence_id" if "sequence_id" in cols else None
			if not id_col:
				return {}

			df = pd.read_sql_query(f"SELECT {id_col} as acc, LENGTH({aln_col}) as aln_len FROM sequence_alignment", conn)
			return dict(zip(df["acc"].astype(str), df["aln_len"]))
		except Exception as e:
			print(f"[warn] Could not load historic alignment lengths: {e}", file=sys.stderr)
			return {}
		finally:
			conn.close()

	def find_gaps_in_fasta(self):
		os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)
		update_scope_accessions = self.load_update_scope_accessions()
		existing_features = self.load_existing_feature_accessions()
		segment_map = self.load_segment_map()
		historical_lengths = self.load_historical_alignment_lengths()

		fasta_file_dir = self.paded_alignment
		if not fasta_file_dir or not os.path.isdir(fasta_file_dir):
			raise FileNotFoundError(f"Padded alignment directory not found: {fasta_file_dir}")
		fasta_files = [f for f in os.listdir(fasta_file_dir) if os.path.isfile(join(fasta_file_dir, f))]
		if not fasta_files:
			raise ValueError(f"No alignment files found in directory: {fasta_file_dir}")

		blast_dict = self.load_blast_hits()
		masters = self.get_master_list()
		master_segment_map = self.load_master_segment_map()
		if not masters:
			raise ValueError("No master accession could be resolved from --master_accession")

		header = [
			"accession",
			"master_ref_accession",
			"reference_accession",
			"aln_start",
			"aln_end",
			"cds_start",
			"cds_end",
			"cds_start_OG_seq",
			"cds_end_OG_seq",
			"product",
			"genome_coverage"
		]
		if segment_map:
			header.append("segment")

		with open(join(self.tmp_dir, self.output_dir, self.output_file), "w") as out_f:
			out_f.write("\t".join(header) + "\n")

			for fasta_file in fasta_files:
				current_master = self.resolve_master_for_alignment(fasta_file, masters, master_segment_map)
				if not current_master:
					print(f"Could not determine master for {fasta_file}. Skipping.")
					continue

				gff_file = self.get_gff_for_master(current_master)
				if not gff_file:
					print(f"No GFF found for master {current_master}. Skipping.")
					continue

				gff_dict = GffDictionary(gff_file).gff_dict
				cds_list = gff_dict['CDS']

				calc = CalculateGenomeCoordinates(join(fasta_file_dir, fasta_file), current_master)
				genome_coords = calc.extract_alignment_coordinates()
				
				fasta_records = list(SeqIO.parse(join(fasta_file_dir, fasta_file), "fasta"))
				master_record = next((r for r in fasta_records if r.id == current_master), None)
				if not master_record:
					master_record = fasta_records[0]
				master_alignment = str(master_record.seq)
				master_coord_to_aln_pos = {}
				master_res_count = 0
				for align_pos, base in enumerate(master_alignment, start=1):
					if base != "-":
						master_res_count += 1
						master_coord_to_aln_pos[master_res_count] = align_pos

				for record in fasta_records:
					record_id = str(record.id).strip()
					sequence = str(record.seq)
					current_len = len(sequence)
					
					old_len = historical_lengths.get(record_id)
					backbone_expanded = old_len is not None and current_len > old_len

					if update_scope_accessions and record_id not in update_scope_accessions:
						if not backbone_expanded:
							continue
						else:
							print(f"[info] Forcing coordinate recalculation for historic ID {record_id} due to backbone expansion ({old_len} -> {current_len})")

					if update_scope_accessions and record_id in existing_features and not backbone_expanded:
						continue

					if record.id in genome_coords:
						master_acc, genome_cord_start, genome_cord_end = genome_coords[record.id]
					else:
						genome_cord_start, genome_cord_end = "NA", "NA"

					gaps = self.get_gap_ranges(sequence)

					if gaps and gaps[0][0] == 1:
						start_offset = gaps[0][1] + 1
					else:
						start_offset = 1

					adjusted = self.recalculate_cds_coordinates_with_span(
						record.id,
						gaps,
						cds_list,
						start_offset,
						genome_cord_start=genome_cord_start,
						genome_cord_end=genome_cord_end,
						master_coord_to_aln_pos=master_coord_to_aln_pos,
					)

					for each_cords in adjusted:
						reference_acc = blast_dict[record.id] if record.id in blast_dict else current_master
						if segment_map:
							record_segment = segment_map.get(record.id, "")
							master_segment = segment_map.get(current_master, record_segment)
							ref_segment = segment_map.get(reference_acc, record_segment)
							if record_segment and master_segment and record_segment != master_segment:
								raise ValueError(f"Segment mismatch for {record.id}: record={record_segment}, master={master_segment}")
							if record_segment and ref_segment and record_segment != ref_segment:
								raise ValueError(f"Segment mismatch for {record.id}: record={record_segment}, ref={ref_segment}")

						# Execute coverage evaluation using tracking maps
						genome_coverage = self.format_genome_coverage(
							sequence,
							master_coord_to_aln_pos,
							each_cords.get('feature_start'),
							each_cords.get('feature_end')
						)

						data = [
							record.id,
							current_master,
							reference_acc,
							str(genome_cord_start),
							str(genome_cord_end),
							str(each_cords['start']),
							str(each_cords['end']),
							str(each_cords['og_start']),
							str(each_cords['og_end']),
							each_cords['product'],
							str(genome_coverage)
						]
						if segment_map:
							data.append(segment_map.get(record.id, ""))
						out_f.write('\t'.join(data) + "\n")

if __name__ == "__main__":
	parser = ArgumentParser(description='Calculates the genome and cds coordinates for a given sequences')
	parser.add_argument('-i', '--paded_alignment', help='Sequence file directory, it can be single or multiple fasta sequence files.', required=True)
	parser.add_argument('-b', '--tmp_dir', help='Base directory', default="tmp")
	parser.add_argument('-d', '--output_dir', help='Output directory where processed data and results are stored', default='Tables')
	parser.add_argument('-o', '--output_file', help='Output file name', default='features.tsv')
	parser.add_argument('-m', '--master_accession', help='Master accession', required=True)
	parser.add_argument('-bh', '--blast_uniq_hits', help='Blast unique hits file', default='tmp/Blast/query_uniq_tophits.tsv')
	parser.add_argument('-g', '--master_gff', help='Master GFF3 file(s)', required=True, nargs='+')
	parser.add_argument('--update_db', help='Existing DB path; when set, only emit feature rows for accessions not already in DB features table', default=None)
	parser.add_argument('--update_scope_tsv', help='TSV with primary_accession column; when set, only recalculate coordinates for these accessions', default=None)
	parser.add_argument('--segment_map_tsv', help='TSV with primary_accession and segment columns for segment-consistency checks', default=None)
	args = parser.parse_args()

	processor = CalculateAlignmentCoordinates(args.paded_alignment, args.master_gff, args.tmp_dir, args.output_dir, args.output_file, args.master_accession, args.blast_uniq_hits, args.update_db, args.update_scope_tsv, args.segment_map_tsv)
	try:
		processor.find_gaps_in_fasta()
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		sys.exit(2)