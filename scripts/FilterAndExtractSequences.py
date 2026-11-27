import os
import csv
import shutil
import tempfile
from collections import Counter
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
from ExportRefListFromUpdateDb import load_reference_dict, load_reference_file_dict
import accession_utils

#genbank_divisions = ['VRL', 'PAT', 'SYN', 'ENV']

class FilterAndExtractSequences:
	def __init__(self, genbank_matrix, sequence_file, genbank_matrix_filtered, ref_file, base_dir, output_dir, total_length, real_length, prop_ambigious, segmented_virus, gb_division, valid_divisions, seq_type, update_db=None, max_n_count=None, max_ambig_fraction=None):
		self.genbank_matrix = genbank_matrix
		self.sequence_file = sequence_file
		self.genbank_matrix_filtered = genbank_matrix_filtered
		self.ref_file = ref_file
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.total_length = total_length
		self.real_length = real_length
		self.prop_ambigious = prop_ambigious
		self.segmented_virus = segmented_virus
		self.gb_division = gb_division
		self.valid_divisions = valid_divisions
		self.seq_type = seq_type
		self.update_db = update_db
		self.max_n_count = int(max_n_count) if max_n_count is not None and str(max_n_count).strip().lower() not in ("", "null", "none", "unset") else None
		self.max_ambig_fraction = float(max_ambig_fraction) if max_ambig_fraction is not None and str(max_ambig_fraction).strip().lower() not in ("", "null", "none", "unset") else None
		os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)

	def fasta_to_dict(self, allowed_ids=None):
		seq_dict = {}
		cur_id = None
		chunks = []
		def _store(rec_id, seq_chunks):
			seq_str = "".join(seq_chunks)
			seq_dict[rec_id] = seq_str
			norm = accession_utils.normalise_accession(rec_id)
			if norm and norm not in seq_dict:
				seq_dict[norm] = seq_str

		with open(self.sequence_file, "r", encoding="utf-8", errors="replace") as fh:
			for raw_line in fh:
				line = raw_line.strip()
				if not line:
					continue
				if line.startswith(">"):
					if cur_id is not None:
						_store(cur_id, chunks)
					parts = line[1:].strip().split(None, 1)
					cand_id = parts[0] if parts else ""
					cand_norm = accession_utils.normalise_accession(cand_id) if cand_id else ""
					if cand_id and (
						allowed_ids is None
						or cand_id in allowed_ids
						or (cand_norm and cand_norm in allowed_ids)
					):
						cur_id = cand_id
						chunks = []
					else:
						cur_id = None
						chunks = []
				elif cur_id is not None:
					chunks.append(line)
			if cur_id is not None:
				_store(cur_id, chunks)
		return seq_dict

	def read_ref_list(self):
		if self.update_db:
			ref_dict = load_reference_dict(self.update_db)
			if ref_dict:
				return ref_dict
		return load_reference_file_dict(self.ref_file)

	def read_ref_list_segmented_virus(self):
		"""Return a dict {accession: type} where type is master/reference/exclusion_list."""
		if self.update_db:
			ref_dict = load_reference_dict(self.update_db)
			if ref_dict:
				return ref_dict
		return load_reference_file_dict(self.ref_file)

	def check_gb_division(self):
		if self.gb_division is not None:
			for each_val in self.gb_division:
				if each_val not in self.genbank_division:
					return False
			return True
		return True

	@staticmethod
	def _classify_accession(row, accession, ref_list, raw_accession=None):
		candidates = [accession]
		if raw_accession and raw_accession != accession:
			candidates.append(raw_accession)
		norm = accession_utils.normalise_accession(accession)
		if norm and norm not in candidates:
			candidates.append(norm)
		for cand in candidates:
			ref_type = str(ref_list.get(cand, '')).strip().lower()
			if ref_type in ('master', 'reference', 'exclusion_list'):
				return ref_type
		row_type = str(row.get('accession_type', '')).strip().lower()
		if row_type:
			return row_type
		for cand in candidates:
			ref_type = str(ref_list.get(cand, '')).strip().lower()
			if ref_type:
				return ref_type
		return 'query'

	def _write_filter_summary(self, total_queries, written_queries, excluded_query_reasons, written_refs):
		summary_lines = [
			"[FilterAndExtractSequences] Summary:",
			f"  - reference/master rows written to ref_seq.fa: {written_refs}",
			f"  - query rows written to query_seq.fa: {written_queries}",
			f"  - query rows excluded before BLAST: {sum(excluded_query_reasons.values())}",
		]
		if excluded_query_reasons:
			summary_lines.append("[FilterAndExtractSequences] Query Exclusion Reasons:")
			for reason, count in excluded_query_reasons.most_common(10):
				summary_lines.append(f"  - {count:<4} : {reason}")
		if total_queries > 0 and written_queries == 0:
			summary_lines.append(
				f"[FilterAndExtractSequences] WARNING: all {total_queries} query rows were filtered before BLAST; query_seq.fa is empty"
			)

		summary_text = "\n".join(summary_lines)
		print(summary_text)
		summary_path = join(self.base_dir, self.output_dir, 'filter_summary.txt')
		with open(summary_path, 'w', encoding='utf-8') as handle:
			handle.write(summary_text + "\n")

	def _collect_needed_fasta_ids(self, ref_list):
		"""Pre-scan metadata so multi-million-row runs only load FASTA records that can survive QC."""
		needed = set()
		with open(self.genbank_matrix, encoding="utf-8", errors="replace") as file:
			csv_reader = csv.DictReader(file, delimiter='\t')
			for row in csv_reader:
				if str(row.get('exclusion_status', '')).strip() == '1':
					continue
				raw_acc = (row.get('gi_number') or row.get('primary_accession') or '').strip()
				if not raw_acc:
					continue
				accession = accession_utils.normalise_accession(raw_acc) or raw_acc
				acc_type = self._classify_accession(row, accession, ref_list, raw_accession=raw_acc)
				if acc_type == 'exclusion_list':
					continue
				if acc_type != 'query':
					needed.add(raw_acc)
					needed.add(accession)
					continue
				gb_division = (row.get('division') or '').strip()
				data_source = str(row.get('data_source') or row.get('dataset_source') or '').strip().lower()
				is_genbank_row = data_source in ('', 'ncbi', 'genbank')
				if is_genbank_row and gb_division:
					if self.gb_division is not None and gb_division not in self.gb_division:
						continue
					if self.valid_divisions is not None and gb_division not in self.valid_divisions:
						continue
				try:
					seq_len = int(row['length'])
					n_count = int(row['n'])
				except (KeyError, ValueError, TypeError):
					# Missing numeric length/n in matrix: inspect sequence directly
					needed.add(raw_acc)
					needed.add(accession)
					continue
				if seq_len < self.total_length or (seq_len - n_count) < self.real_length:
					continue
				if self.max_n_count is not None and n_count > self.max_n_count:
					continue
				if self.max_ambig_fraction is not None and seq_len > 0:
					raw_real = str(row.get('real_length', '')).strip()
					if raw_real != '':
						try:
							atgc_count = int(raw_real)
							if (max(0, seq_len - atgc_count) / seq_len) > self.max_ambig_fraction:
								continue
						except (ValueError, TypeError):
							pass
				needed.add(raw_acc)
				needed.add(accession)
		return needed

	def filter_columns(self):
		#if not self.check_gb_division():
		#	print("Error: Invalid GenBank division specified.")
		#	return
		exclusion_dict = {}
	
		# Both methods now return dict {accession: type}
		if self.segmented_virus == "Y":
			raw_ref_list = self.read_ref_list_segmented_virus()
		else:
			raw_ref_list = self.read_ref_list()

		ref_list = dict(raw_ref_list or {})
		for acc_key, acc_type in list(ref_list.items()):
			norm_k = accession_utils.normalise_accession(acc_key)
			if norm_k and norm_k not in ref_list:
				ref_list[norm_k] = acc_type

		needed_ids = self._collect_needed_fasta_ids(ref_list)
		sequence_dict = self.fasta_to_dict(allowed_ids=needed_ids)

		# Identify which refs are exclusion_list type
		exclusion_list_refs = {
			accession_utils.normalise_accession(acc) or acc
			for acc, acc_type in raw_ref_list.items()
			if str(acc_type).strip().lower() == 'exclusion_list'
		}
		if exclusion_list_refs:
			print(f"Found {len(exclusion_list_refs)} exclusion_list references: {', '.join(sorted(exclusion_list_refs))}")
			# Write exclusion ref list for downstream use by BlastAlignment
			exclusion_refs_file = join(self.base_dir, self.output_dir, 'exclusion_refs.txt')
			with open(exclusion_refs_file, 'w') as ef:
				for acc in sorted(exclusion_list_refs):
					ef.write(f"{acc}\n")

		query_seq_file = join(self.base_dir, self.output_dir, 'query_seq.fa')
		ref_seq_file = join(self.base_dir, self.output_dir, 'ref_seq.fa')
		excluded_query_reasons = Counter()
		total_queries = 0
		written_queries = 0
		written_refs = 0

		tmp_path = None
		try:
			with open(query_seq_file, 'w') as write_query_seq, open(ref_seq_file, 'w') as write_ref_seq:
				with open(self.genbank_matrix) as file, tempfile.NamedTemporaryFile('w', delete=False, newline='') as tmpfile:
					tmp_path = tmpfile.name
					csv_reader = csv.DictReader(file, delimiter='\t')
					fieldnames = list(csv_reader.fieldnames or [])

					if 'exclusion_status' not in fieldnames:
						fieldnames += ['exclusion_status']
					if 'exclusion_criteria' not in fieldnames:
						fieldnames += ['exclusion_criteria']

					writer = csv.DictWriter(tmpfile, fieldnames=fieldnames, delimiter='\t')
					writer.writeheader()

					for row in csv_reader:
						raw_acc = (row.get('gi_number') or row.get('primary_accession') or '').strip()
						accession = accession_utils.normalise_accession(raw_acc) or raw_acc
						acc_type = self._classify_accession(row, accession, ref_list, raw_accession=raw_acc)
						is_query = acc_type == 'query'
						if is_query:
							total_queries += 1

						if str(row.get('exclusion_status', '')).strip() == '1':
							if is_query:
								reason = (row.get('exclusion_criteria') or '').strip() or 'metadata_exclusion'
								excluded_query_reasons[reason] += 1
							writer.writerow(row)
							continue

						gb_division = (row.get('division') or '').strip()
						data_source = str(row.get('data_source') or row.get('dataset_source') or '').strip().lower()
						is_genbank_row = data_source in ('', 'ncbi', 'genbank')

						sequence = sequence_dict.get(raw_acc) or sequence_dict.get(accession) or ""

						try:
							seq_len = int(row['length'])
						except (ValueError, TypeError, KeyError):
							seq_len = len(sequence)

						try:
							n_count = int(row['n'])
						except (ValueError, TypeError, KeyError):
							n_count = sequence.upper().count('N') if sequence else 0

						seq_len_without_n = seq_len - n_count

						exclusion_status = 0
						exclusion_criteria = ""
						if acc_type == 'exclusion_list':
							exclusion_status = 1
							exclusion_criteria = "excluded by reference list exclusion flag"
						elif not sequence.strip():
							exclusion_status = 1
							exclusion_criteria = "missing_or_empty_sequence"

						if is_query and exclusion_status == 0:
							if is_genbank_row and gb_division:
								if self.gb_division is not None and gb_division not in self.gb_division:
									exclusion_status = 1
									exclusion_criteria = f"GenBank division {gb_division} not in user-defined exclusion list"

								elif self.valid_divisions is not None and gb_division not in self.valid_divisions:
									exclusion_status = 1
									exclusion_criteria = f"GenBank division {gb_division} not in valid division list"

							if seq_len < self.total_length:
								exclusion_status = 1
								exclusion_criteria = "Sequence length less than total length"

							if seq_len_without_n < self.real_length:
								exclusion_status = 1
								exclusion_criteria = "Real sequence length (length - Ns) less than the threshold"

							if self.max_n_count is not None and n_count > self.max_n_count:
								exclusion_status = 1
								exclusion_criteria = f"N count ({n_count}) exceeds max_n_count ({self.max_n_count})"

							if self.max_ambig_fraction is not None and seq_len > 0:
								try:
									atgc_count = int(row['real_length']) if 'real_length' in row and str(row['real_length']).strip() != '' else sum(sequence.upper().count(b) for b in 'ATGC')
								except (ValueError, TypeError):
									atgc_count = sum(sequence.upper().count(b) for b in 'ATGC') if sequence else seq_len_without_n
								ambig_count = max(0, seq_len - atgc_count)
								ambig_frac = ambig_count / seq_len
								if ambig_frac > self.max_ambig_fraction:
									exclusion_status = 1
									exclusion_criteria = f"Ambiguous/N fraction ({ambig_frac:.5f}) exceeds max_ambig_fraction ({self.max_ambig_fraction})"

						existing_criteria = (row.get('exclusion_criteria') or '').strip()
						if existing_criteria:
							row['exclusion_criteria'] = f"{existing_criteria}; {exclusion_criteria}" if exclusion_criteria else existing_criteria
						else:
							row['exclusion_criteria'] = exclusion_criteria

						row['exclusion_status'] = exclusion_status
						writer.writerow(row)

						if exclusion_status == 1:
							exclusion_dict[accession] = row['exclusion_criteria']
							if is_query:
								excluded_query_reasons[row['exclusion_criteria'] or 'metadata_exclusion'] += 1
							continue

						if acc_type != 'query':
							write_ref_seq.write(f">{accession}\n{sequence}\n")
							written_refs += 1
						else:
							write_query_seq.write(f">{accession}\n{sequence}\n")
							written_queries += 1

			# Replace original matrix file with updated temp file
			shutil.move(tmp_path, self.genbank_matrix)
			tmp_path = None
		finally:
			if tmp_path and os.path.exists(tmp_path):
				try:
					os.remove(tmp_path)
				except OSError:
					pass
		self._write_filter_summary(total_queries, written_queries, excluded_query_reasons, written_refs)

	def process(self):
		self.filter_columns()

 
if __name__ == "__main__":
	parser = ArgumentParser(description='Filter sequences and prepare sequences for BLAST alignment')
	parser.add_argument('-g', '--genbank_matrix', help='GenBank matrix file', required=True)
	parser.add_argument('-sf', '--sequence_file', help='fasta sequence file', default='tmp/GenBank-matrix/sequences.fa')
	parser.add_argument('-f', '--genbank_matrix_filtered', help='Filtered GenBank matrix file storage directory', default='tmp/GenBank-matrix')
	parser.add_argument('-r', '--ref_file', help='Text file containing list of reference sequence accessions', required=True)
	parser.add_argument('-b', '--base_dir', help='Base directory', default='tmp')
	parser.add_argument('-o', '--output_dir', help='Output directory to process the data', default="Sequences")
	parser.add_argument('-l', '--total_length', help='Total length of the sequence should be more or equal to the length provided', default=1, type=int)
	parser.add_argument('-n', '--real_length', help='Length of sequences without N', default=1, type=int)
	parser.add_argument('-a', '--prop_ambigious_data', help='Proportion of ambiguous data to be excluded.', nargs='+', default=None, type=str)
	parser.add_argument('-v', '--segmented_virus', help='Is segmented virus', default="N")
	parser.add_argument(
		'-d', '--genbank_division',
		help=(
				'GenBank division to be excluded. General divisions are:\n'
				'    "VRL" = Viral sequences\n'
				'    "ENV" = Environmental sequences\n'
				'    "PAT" = Patented sequences\n'
				'    "SYN" = Synthetic and chimeric sequences\n'
			),
				nargs='+',
				default=['VRL', 'PAT', 'SYN', 'ENV', 'PHG'],
				type=str
    )
	parser.add_argument('-vd', '--valid_divisions', help="Valid GenBank divisions to be considered for the analysis", nargs='+', default=['VRL', 'ENV', 'PHG'], type=str)
	parser.add_argument('-s', '--seq_type', help='Sequence type', default=None, type=str)
	parser.add_argument('--update_db', help='Existing update DB; when set, derive reference accessions from DB instead of ref_file', default=None)
	parser.add_argument('--max_n_count', help='Maximum allowed N bases in query sequences (default: unlimited)', default=None, type=int)
	parser.add_argument('--max_ambig_fraction', help='Maximum allowed non-ATGC fraction in query sequences (default: unlimited)', default=None, type=float)
	args = parser.parse_args()

	processor = FilterAndExtractSequences(
		genbank_matrix=args.genbank_matrix,
		sequence_file = args.sequence_file,
		genbank_matrix_filtered=args.genbank_matrix_filtered,
		ref_file=args.ref_file,
		base_dir = args.base_dir,
		output_dir=args.output_dir,
		total_length=args.total_length,
		real_length=args.real_length,
		prop_ambigious=args.prop_ambigious_data,
		segmented_virus=args.segmented_virus,
		gb_division=args.genbank_division,
		valid_divisions=args.valid_divisions,
		seq_type=args.seq_type,
		update_db=args.update_db,
		max_n_count=args.max_n_count,
		max_ambig_fraction=args.max_ambig_fraction,
	)
	processor.process()
