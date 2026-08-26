#python blastAlignment.py -s Y -f generic-influenza/ref_list.txt

import os
import csv
import sys
import sqlite3
import shutil
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
import read_file
from ExportRefListFromUpdateDb import load_master_accessions_from_file, load_reference_file_table, load_reference_rows
import segment_utils

class BlastAlignment:
	def __init__(self, query_fasta, db_fasta, base_dir, output_dir, output_file, is_segmented_virus, master_acc, is_update, keep_blast_tmp_dir, gb_matrix, segment_file=None, update_db=None, threads=1):
		self.query_fasta = query_fasta
		self.db_fasta = db_fasta
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.output_file = output_file
		self.is_segmented_virus = is_segmented_virus
		self.segment_file = segment_file
		self.master_acc = master_acc
		self.gb_matrix = gb_matrix
		self.is_update = is_update
		self.keep_blast_tmp_dir = keep_blast_tmp_dir
		self.update_db = update_db
		self.threads = max(1, int(threads))
		self.db_file_name = os.path.basename(db_fasta)

	@staticmethod
	def normalize_query_id(header):
		if header is None:
			return ""
		return header.strip().split()[0]

	@staticmethod
	def _require_file(path, label):
		if not path or not os.path.isfile(path):
			raise FileNotFoundError(f"{label} file not found: {path}")

	@staticmethod
	def _validate_tsv_columns(path, required_columns, label):
		with open(path, newline='', encoding='utf-8') as handle:
			reader = csv.DictReader(handle, delimiter='\t')
			fieldnames = reader.fieldnames or []
		missing = [c for c in required_columns if c not in fieldnames]
		if missing:
			raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")

	@staticmethod
	def _fasta_ids(path):
		return [record.id for record in SeqIO.parse(path, "fasta")]

	@staticmethod
	def _normalize_segment(value):
		"""Delegates to :mod:`segment_utils` - the single normalisation authority.

		This used to scrape every digit out of the string, which turned ``4.0`` into
		segment 40 and, worse, inverted the polymerase segments: ``PB2`` (segment 1)
		became ``2`` and ``PB1`` (segment 2) became ``1``. Missing still maps to
		``"0"`` here, because callers of this particular helper depend on it.
		"""
		normalised = segment_utils.normalise_segment(value)
		if normalised is None or normalised.casefold() in segment_utils.PANDAS_NULL_TOKENS:
			return "0"
		return normalised

	def hydrate_update_reference_assets(self):
		if not self.update_db:
			return
		if not os.path.isfile(self.update_db):
			raise FileNotFoundError(f"Update DB not found: {self.update_db}")

		output_dir = join(self.base_dir, self.output_dir)
		os.makedirs(output_dir, exist_ok=True)
		ref_list_path = join(output_dir, "ref_list_from_update_db.tsv")
		ref_fasta_path = join(output_dir, "ref_seq_from_update_db.fa")

		refs = load_reference_rows(self.update_db)

		conn = sqlite3.connect(self.update_db)
		try:
			df_seq = pd.read_sql_query("SELECT header, sequence FROM sequences", conn)
		finally:
			conn.close()
		if refs.empty:
			raise ValueError(f"Update DB does not contain any master/reference rows in meta_data: {self.update_db}")

		df_seq["header"] = df_seq["header"].fillna("").astype(str).str.strip()
		df_seq["sequence"] = df_seq["sequence"].fillna("").astype(str).str.strip()
		ref_seq_df = refs.merge(df_seq, left_on="primary_accession", right_on="header", how="left")
		missing_sequences = sorted(
			ref_seq_df.loc[ref_seq_df["sequence"].eq("") | ref_seq_df["sequence"].isna(), "primary_accession"].dropna().unique().tolist()
		)
		if missing_sequences:
			raise ValueError(
				f"Update DB is missing sequences for reference accession(s): {missing_sequences}"
			)

		with open(ref_list_path, "w", encoding="utf-8") as ref_handle:
			for _, row in refs.iterrows():
				ref_handle.write(
					f"{row['primary_accession']}\t{row['accession_type']}\t{row['segment']}\n"
				)

		with open(ref_fasta_path, "w", encoding="utf-8") as fasta_handle:
			for _, row in ref_seq_df.drop_duplicates(subset=["primary_accession"]).iterrows():
				fasta_handle.write(f">{row['primary_accession']}\n{row['sequence']}\n")

		self.db_fasta = ref_fasta_path
		self.db_file_name = os.path.basename(ref_fasta_path)
		self.segment_file = ref_list_path
		self.master_acc = ref_list_path
	
	def read_gb_matrix(self):
		self._require_file(self.gb_matrix, "GenBank matrix")
		self._validate_tsv_columns(self.gb_matrix, ["gi_number", "sequence"], "GenBank matrix")
		accessions = {}
		with open(self.gb_matrix, newline='', encoding='utf-8') as file:
			reader = csv.DictReader(file, delimiter='\t') 
			for row in reader:
				accessions[row["gi_number"]] = row["sequence"]

		return accessions

	def get_exclusion_list_refs(self):
		"""Read the ref_list/segment file and return a set of accessions marked as exclusion_list."""
		exclusion_refs = set()
		# Use segment_file if available (segmented), otherwise master_acc (which may be a file)
		ref_file = self.segment_file or self.master_acc
		if ref_file and os.path.isfile(ref_file):
			try:
				df = load_reference_file_table(ref_file)
				exclusion_refs = set(
					df[df["accession_type"].astype(str).str.strip().str.lower() == "exclusion_list"]["primary_accession"]
					.astype(str)
					.str.strip()
					.tolist()
				)
			except Exception:
				pass
		if exclusion_refs:
			print(f"[exclusion_list] Found {len(exclusion_refs)} exclusion_list references: {', '.join(sorted(exclusion_refs))}")
		return exclusion_refs

	def update(self, query_tmp_dir):
		if os.path.exists(self.query_fasta):
			query_accession = [record.id for record in SeqIO.parse(self.query_fasta, "fasta")]
			db_accession = [record.id for record in SeqIO.parse(self.db_fasta, "fasta")]
			all_accessions = self.read_gb_matrix()
			missing_accessions = [acc for acc in all_accessions.keys() if acc not in query_accession and acc not in db_accession]
			missing_accessions_count = len(missing_accessions)
			if missing_accessions_count > 0:
				print(f'{missing_accessions_count} new sequences to process')
				write_file = open(join(query_tmp_dir, "query.fa"), 'w')
				for each_missing_acc in missing_accessions:
					write_file.write(">" + each_missing_acc)
					write_file.write("\n")
					write_file.write(all_accessions[each_missing_acc])
					write_file.write("\n")
				# create temp query.fa wich should be used by either segmented or non-segmented blast analysis
			else:
				print("No new accession id's to process")

		else:
			print("No query file exists for update, you need to run it without update to perform blast on existing sequences")
		
	@staticmethod
	def check_blast_exists(command):
		try:
			subprocess.run([command, '-version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			print(f"{command} is working.")
			return True
		except subprocess.CalledProcessError:
			print(f"{command} is not installed correctly.")
			return False
		except FileNotFoundError:
			print(f"{command} is not found on the system.")
			return False
	
	@staticmethod
	def delete_directory(dir_path):
		if os.path.exists(dir_path):
			if os.path.isdir(dir_path):
				shutil.rmtree(dir_path) 
				print(f"Directory '{dir_path}' has been deleted. Ignore the message")
			else:
				print(f"'{dir_path}' exists but is not a directory. Ignore the message")
		else:
			print(f"Directory '{dir_path}' does not exist. Ignore the message")

	@staticmethod
	def ref_segments(query_tophits_annotated):
		segment_dict = {}
		for each_line in open(query_tophits_annotated):
			query, ref, score, strand, segment = each_line.strip().split('\t')
			if segment not in segment_dict:
				segment_dict[segment] = {}
        
			if ref not in segment_dict[segment]:
				segment_dict[segment][ref] = []
        
			segment_dict[segment][ref].append(query)

		return segment_dict
				
	def run_makeblastdb(self, tmp_dir):
		os.makedirs(join(tmp_dir, 'DB'), exist_ok=True)
		command = [
			'makeblastdb',
			'-in', self.db_fasta,
			'-out', join(tmp_dir, 'DB', self.db_file_name),
			'-title', "alignment",
			'-dbtype', 'nucl'
		]
		try:
			subprocess.run(command, check=True)
			print(f"makeblastdb ran successfully on {self.db_fasta}")
		except subprocess.CalledProcessError as e:
			raise RuntimeError(f"makeblastdb failed for reference FASTA '{self.db_fasta}': {e}") from e

	def run_blastn(self, output_dir, query_file):
		command = [
			'blastn',
			'-query', query_file,
			'-db', join(output_dir, "DB", self.db_file_name),
			'-task', 'blastn',
			'-max_target_seqs', '1',
			'-max_hsps', '1',
			'-out', join(output_dir, self.output_file),
			'-outfmt', "6 qacc sacc pident sstrand",
			'-num_threads', str(self.threads)
			]
		try:
			subprocess.run(command, check=True)
			print(f"blastn ran successfully. Results saved in {self.output_file}")
		except subprocess.CalledProcessError as e:
			raise RuntimeError(f"blastn failed for query '{query_file}' against db '{join(output_dir, 'DB', self.db_file_name)}': {e}") from e

	def get_master_list(self):
		if os.path.isfile(self.master_acc):
			try:
				return load_master_accessions_from_file(self.master_acc)
			except:
				return []
		else:
			return [x.strip() for x in self.master_acc.split(',') if x.strip()]

	def write_filtered_ref_fasta(self, output_dir, exclusion_refs):
		filtered_ref_path = join(output_dir, "ref_seq_filtered.fa")
		input_refs = read_file.fasta(self.db_fasta)
		count = 0
		with open(filtered_ref_path, 'w') as out_f:
			for header, seq in input_refs:
				# Use the first token of the header as accession for matching
				acc = header.strip().split()[0]
				if acc not in exclusion_refs:
					out_f.write(f">{header}\n{seq}\n")
				else:
					count += 1
		if count:
			print(f"[exclusion_list] Excluded {count} references from ref_seq_filtered.fa")

	def write_master_seq(self, output_dir):
		masters = self.get_master_list()
		with open(self.db_fasta, "r") as infile:
			records = SeqIO.to_dict(SeqIO.parse(infile, "fasta"))
			
		for acc in masters:
			if acc in records:
				with open(join(output_dir, acc + '.fasta'), "w") as outfile:
					SeqIO.write(records[acc], outfile, "fasta")
					print(f"Sequence '{acc}' has been saved to {join(output_dir, acc)}")
			else:
				print(f"Sequence ID '{acc}' not found in {self.db_fasta}")
	
	def process_non_segmented_virus(self, output_dir, query_fasta):
		input_file = join(output_dir, "query_tophits.tsv")
		query_tophit_uniq = join(output_dir, "query_uniq_tophits.tsv")
		grouped_fasta = join(output_dir, "grouped_fasta")
		sorted_fasta = join(output_dir, "sorted_fasta")
		merged_fasta = join(output_dir, "merged_fasta")
		sorted_all = join(output_dir, "sorted_all")		
		ref_seq_dir = join(output_dir, "ref_seqs")
		master_seq = join(output_dir, "master_seq")

		os.makedirs(grouped_fasta, exist_ok=True)
		os.makedirs(ref_seq_dir, exist_ok=True)
		os.makedirs(sorted_fasta, exist_ok=True)
		os.makedirs(merged_fasta, exist_ok=True)
		os.makedirs(sorted_all, exist_ok=True)
		os.makedirs(master_seq, exist_ok=True)

		exclusion_refs = self.get_exclusion_list_refs()
		records = {}
		values = {}

		self.write_master_seq(master_seq)
	
		with open(input_file, newline='') as file:
			reader = csv.reader(file, delimiter='\t')
			for row in reader:
				col1, col2, col3, col4 = row[0], row[1], float(row[2]), row[3]
				if col1 in values:
					existing_value = values[col1]
					if col3 > existing_value:
						records[col1] = [col1, col2, col3, col4]
						values[col1] = col3
				else:
					records[col1] = [col1, col2, col3, col4]
					values[col1] = col3

		with open(query_tophit_uniq, 'w', newline='') as file:
			writer = csv.writer(file, delimiter='\t')
			excluded_count = 0
			for record in records.values():
				# Skip queries whose best hit is an exclusion_list reference
				if record[1] in exclusion_refs:
					excluded_count += 1
					continue
				writer.writerow(record)
			if excluded_count:
				print(f"[exclusion_list] Excluded {excluded_count} queries matching exclusion_list references")
	
		seq_dicts = {}
		query_seqs = read_file.fasta(query_fasta)  
		for rows in query_seqs:
			raw_header = rows[0].strip()
			seq = rows[1].strip()
			seq_dicts[raw_header] = seq
			norm_header = self.normalize_query_id(raw_header)
			if norm_header and norm_header != raw_header:
				seq_dicts[norm_header] = seq

		#Seperate plus and minus strand sequences
		for each_line in open(query_tophit_uniq, 'r'):
			query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
			if query_acc not in seq_dicts:
				print(f"[warn] Missing query sequence for accession: {query_acc}")
				continue
			if strand == "plus":
				with open(join(sorted_fasta, "plus.fa"), 'a') as file_plus:
					file_plus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")
			else:
				with open(join(sorted_fasta, "minus.fa"), 'a') as file_minus:
					file_minus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")

		for each_file in os.listdir(sorted_fasta):
			if "minus" in each_file:
				command = ["seqkit", "seq", "-j", str(self.threads), "-r", "-p", "-v", "-t", "dna", join(sorted_fasta, each_file), ">", join(merged_fasta, each_file)]
			else:
				command = ["cp", join(sorted_fasta, each_file), join(merged_fasta, each_file)]

			try:
				print(' '.join(command))
				os.system(' '.join(command))
				print(f"seqkit ran successfully for {each_file}")
			except subprocess.CalledProcessError as e:
				print(f"Error running seqkit: {e}")

		file_list = []
		for each_file in os.listdir(merged_fasta):
			prefix = merged_fasta + "/"
			file_list.append(prefix + each_file)
		
		command = ["cat", " ".join(file_list), ">", join(sorted_all, "query_seq.fa")]
		print(' '.join(command))
		try:
			os.system(' '.join(command))
			print(f"concatenation successful")
		except subprocess.CalledProcessError as e:
			print(f"Error in concatenation: {e}")

		grouped_dict = {}
		for each_line in open(query_tophit_uniq, 'r'):
			query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
			if ref_acc not in grouped_dict:
				grouped_dict[ref_acc] = [query_acc]
			else:
				grouped_dict[ref_acc].append(query_acc)

		print('\n'.join(grouped_dict.keys()))
		seq_dicts = {}
		query_seqs = read_file.fasta(join(sorted_all, "query_seq.fa"))
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		for each_ref_acc, list_of_query_acc in grouped_dict.items():
			# Don't write grouped fasta or ref_seqs for exclusion_list references
			if each_ref_acc in exclusion_refs:
				continue
			with open(join(grouped_fasta, each_ref_acc + '.fasta'), 'a') as write_file:
				for each_query_acc in list_of_query_acc:
					seqs = seq_dicts[each_query_acc]
					write_file.write(">" + each_query_acc + '\n')
					for i in range(0, len(seqs), 80):
						write_file.write(seqs[i:i + 80] + '\n')
		
		#writing reference sequences into individual fasta files
		ref_seqs = read_file.fasta(self.db_fasta)
		for rows in ref_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()
			
		for ref_accs in grouped_dict.keys():
			# Don't write ref_seq files for exclusion_list references
			if ref_accs in exclusion_refs:
				continue
			seqs = seq_dicts[ref_accs]
			with open(join(ref_seq_dir, ref_accs + '.fasta'), 'w') as write_file:
				write_file.write(">" + ref_accs + '\n')
				for i in range(0, len(seqs), 80):
					write_file.write(seqs[i:i + 80] + '\n')
			
	def process_segment_virus(self, input_file, uniq_hit_output, segment_file, annotated_output):
		uniq_hits = {}
		segments = {}
		segment_assigned = {}
		seg_dict = {}

		exclusion_refs = self.get_exclusion_list_refs()

		os.makedirs(join(self.base_dir, self.output_dir, "grouped_fasta"), exist_ok=True)
		os.makedirs(join(self.base_dir, self.output_dir, "segment_sorted"), exist_ok=True)
		os.makedirs(join(self.base_dir, self.output_dir, "segment_merged_fasta"), exist_ok=True)
		os.makedirs(join(self.base_dir, self.output_dir, "segment_sorted_all"), exist_ok=True)
		os.makedirs(join(self.base_dir, self.output_dir, "ref_seqs"), exist_ok=True)
		os.makedirs(join(self.base_dir, self.output_dir, "master_seq"), exist_ok=True)

		segment_sorted = join(self.base_dir, self.output_dir, "segment_sorted")
		segment_merged = join(self.base_dir, self.output_dir, "segment_merged_fasta")
		segment_sorted_all =join(self.base_dir, self.output_dir, "segment_sorted_all")
		grouped_fasta = join(self.base_dir, self.output_dir, "grouped_fasta")
		ref_seqs = join(self.base_dir, self.output_dir, "ref_seqs")
		master_seq = join(self.base_dir, self.output_dir, "master_seq")
		
		self.write_master_seq(master_seq)

		with open(input_file, newline='') as file:
			reader = csv.reader(file, delimiter='\t')
			for row in reader:
				col1, col2, col3, col4 = row[0], row[1], row[2], row[3]
				if col1 in uniq_hits:
					if float(col3) > float(uniq_hits[col1][2]):
						uniq_hits[col1] = [col1, col2, col3, col4]
				else:
					uniq_hits[col1] = [col1, col2, col3, col4]

		# Filter out queries whose best hit is an exclusion_list reference
		excluded_query_count = 0
		with open(uniq_hit_output, 'w') as write_uniq_hits:
			for k, v in uniq_hits.items():
				if v[1] in exclusion_refs:
					excluded_query_count += 1
					continue
				write_uniq_hits.write('\t'.join(v) + '\n')
		if excluded_query_count:
			print(f"[exclusion_list] Excluded {excluded_query_count} queries matching exclusion_list references")
		
		# Build segment lookup, but also track which segments are exclusion-only
		exclusion_segments = set()
		for line in open(segment_file):
				parts = line.strip().split('\t')
				if len(parts) >= 2:
					accession = parts[0]
					acc_type = parts[1].strip().lower() if len(parts) >= 2 else ''
					segment = parts[-1] # Assume last column is segment
				else:
					continue

				accession = accession.split('|')[0]
				segments[accession] = segment
				if acc_type == 'exclusion_list':
					exclusion_segments.add(segment)

		# Log which segments are being excluded
		if exclusion_segments:
			print(f"[exclusion_list] Segments with exclusion_list refs: {', '.join(sorted(exclusion_segments))}")

		with open(annotated_output, 'w') as write_segment:
			with open(uniq_hit_output, newline='') as file:
				reader = csv.reader(file, delimiter='\t')
				for row in reader:
					col1, col2, col3, col4 = row[0], row[1], row[2], row[3]
					if col2 in segments:
						val = [col1, col2, col3, col4, segments[col2]]
						segment_assigned[col1] = val
						write_segment.write('\t'.join(val) + '\n')
					else:
						print("Could not find the segment for :" + col2)

		fasta_seqs = read_file.fasta(self.query_fasta)
		for each_seq in fasta_seqs:
			raw_header = each_seq[0].strip()
			header = self.normalize_query_id(raw_header)
			sequence = each_seq[1].strip()

			assigned_header = None
			if header in segment_assigned:
				assigned_header = header
			elif raw_header in segment_assigned:
				assigned_header = raw_header

			if assigned_header:
				accession, reference, score, strand, segment = segment_assigned[assigned_header]
				if strand == "plus":
					with open(join(segment_sorted, f"seg_{segment}_plus.fa"), 'a') as write_file:
						write_file.write(f">{assigned_header}\n{sequence}\n")
				else:
						with open(join(segment_sorted, f"seg_{segment}_minus.fa"), 'a') as write_file:
							write_file.write(f">{assigned_header}\n{sequence}\n")
					
		for each_seg in os.listdir(segment_sorted):
			name, seg_num, strand = each_seg.split('.')[0].split('_')
			print(f"Processing {name}-{strand}")
			if strand == "minus":
				command = ["seqkit", "seq", "-j", str(self.threads), "-r", "-p", "-v", "-t", "dna", join(segment_sorted, each_seg), ">", join(segment_merged, each_seg)]
			else:
				command = ["cp", join(segment_sorted, each_seg), join(segment_merged, each_seg)]

			try:
				print(' '.join(command))
				os.system(' '.join(command))
				print(f"seqkit ran successfully for {each_seg}")
			except subprocess.CalledProcessError as e:
				print(f"Error running seqkit: {e}")

		for each_seg in os.listdir(segment_merged):
			name, seg_num, strand = each_seg.split('.')[0].split('_')
			if name + "_" + seg_num not in seg_dict:
				seg_dict[name + "_" + seg_num] = [each_seg] 
			else:
				seg_dict[name + "_" + seg_num].append(each_seg)

		for each_segment, files in seg_dict.items():
			prefix = segment_merged + "/"
			output_file = join(segment_sorted_all, each_segment + ".fa")
			file_list = [prefix + item for item in files]
			command = ["cat", " ".join(file_list), ">", output_file]
			print(' '.join(command))
			try:
				os.system(' '.join(command))
				print(f"{each_segment} concatenated sucessfully")
			except subprocess.CalledProcessError as e:
				print(f"Error in concatenation: {e}")

		segment_dictionary = self.ref_segments(join(self.base_dir, self.output_dir, "query_uniq_tophit_annotated.tsv"))
		missing_queries = set()
		for segment, ref_acc in segment_dictionary.items():

			seq_dict = {}
			segment_fa = read_file.fasta(join(segment_sorted_all, f'seg_{segment}.fa'))
			for header, sequence in segment_fa:
				seq_dict[header] = sequence
			print(f"loaded seg_{segment}.fa")

			for each_ref_acc, query_acc in ref_acc.items():
				# Don't write grouped fasta for exclusion_list references
				if each_ref_acc in exclusion_refs:
					continue
				write_file = open(join(grouped_fasta, each_ref_acc + ".fa"), 'w')
				for each_query in query_acc:
						if each_query not in seq_dict:
							missing_queries.add(each_query)
							continue
						write_file.write(">" + each_query + "\n" + seq_dict[each_query] + "\n")
				write_file.close()

		if missing_queries:
			print(f"[warn] {len(missing_queries)} query sequences were missing after segment FASTA build; skipped while writing grouped FASTA")
			preview = ', '.join(sorted(missing_queries)[:10])
			print(f"[warn] Missing query accession examples: {preview}")

		reference_seqs = read_file.fasta(self.db_fasta)
	
		for header, sequence in reference_seqs:
			# Don't write ref_seq files for exclusion_list references
			if header in exclusion_refs:
				continue
			write_file = open(join(ref_seqs, header + ".fa"), "w")
			write_file.write(">" + header + "\n" + sequence + "\n")
			write_file.close()

	def update_gB_matrix(self, query_fasta, query_tophit_uniq, gB_matrix_file):
		self._require_file(query_fasta, "query fasta")
		self._require_file(query_tophit_uniq, "unique BLAST hits")
		self._require_file(gB_matrix_file, "GenBank matrix")
		self._validate_tsv_columns(gB_matrix_file, ["gi_number"], "GenBank matrix")

		exclusion_refs = self.get_exclusion_list_refs()

		uniq_blast_acc = {}
		with open(query_tophit_uniq) as f:
			for line in f:
				parts = line.strip().split("\t")
				if len(parts) != 4:
					raise ValueError(f"Malformed unique BLAST hit row in {query_tophit_uniq}: {line.strip()}")
				query_acc, ref_acc, score, strand = parts
				uniq_blast_acc[query_acc] = ref_acc

		# Also read the raw blast hits to identify which queries hit exclusion refs
		# (these were filtered from query_uniq_tophits but we need to mark them specifically)
		raw_blast_file = os.path.join(os.path.dirname(query_tophit_uniq), "query_tophits.tsv")
		exclusion_hit_queries = set()
		if os.path.exists(raw_blast_file) and exclusion_refs:
			raw_hits = {}
			with open(raw_blast_file) as f:
				for line in f:
					parts = line.strip().split("\t")
					if len(parts) >= 4:
						q, r, s, st = parts[0], parts[1], float(parts[2]), parts[3]
						if q not in raw_hits or s > raw_hits[q][1]:
							raw_hits[q] = (r, s)
			for q, (r, s) in raw_hits.items():
				if r in exclusion_refs:
					exclusion_hit_queries.add(q)

		query_acc_status = {}
		read_query_obj = read_file.fasta(query_fasta)
		for each_seq_obj in read_query_obj:
			header = each_seq_obj[0]
			if header not in uniq_blast_acc:
				query_acc_status[header] = 1  # excluded
			else:
				query_acc_status[header] = 0  # not excluded

		updated_rows = []
		with open(gB_matrix_file, newline='') as infile:
			reader = csv.DictReader(infile, delimiter='\t')
			fieldnames = reader.fieldnames or []

			if 'exclusion_status' not in fieldnames:
				fieldnames.append('exclusion_status')
			if 'exclusion_criteria' not in fieldnames:
				fieldnames.append('exclusion_criteria')

			for row in reader:
				gi = row.get('gi_number')
				status = query_acc_status.get(gi, 0)

				existing_criteria = row.get('exclusion_criteria', '')
				if gi in exclusion_hit_queries:
					new_criteria = 'excluded: best BLAST hit is an exclusion_list reference'
					if existing_criteria:
						row['exclusion_criteria'] = f"{existing_criteria}; {new_criteria}"
					else:
						row['exclusion_criteria'] = new_criteria
					row['exclusion_status'] = '1'
				elif status == 1 and gi in query_acc_status:
					new_criteria = 'excluded due to no hit'
					if existing_criteria:
						row['exclusion_criteria'] = f"{existing_criteria}; {new_criteria}"
					else:
							row['exclusion_criteria'] = new_criteria
					row['exclusion_status'] = '1'
				else:
					row['exclusion_criteria'] = existing_criteria

				updated_rows.append(row)

		with open(gB_matrix_file, 'w', newline='') as outfile:
			writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
			writer.writeheader()
			writer.writerows(updated_rows)

	'''
	def update_gB_matrix(self, query_fasta, query_tophit_uniq, gB_matrix_file):
		uniq_blast_acc = {}
		with open(query_tophit_uniq) as f:
			for line in f:
				query_acc, ref_acc, score, strand = line.strip().split("\t")
				uniq_blast_acc.setdefault(query_acc, ref_acc)	

		query_acc_status = {}
		read_query_obj = read_file.fasta(query_fasta)
		for each_seq_obj in read_query_obj:
			header = each_seq_obj[0]
			seq = each_seq_obj[1]
			if header not in uniq_blast_acc:
				query_acc_status[header] = 1 # excluded
			else:
				query_acc_status[header] = 0 # not excluded

		df = pd.read_csv(gB_matrix_file, sep='\t')

		df['exclusion_status'] = df['gi_number'].map(query_acc_status).fillna(0).astype(int)
		
		df['exclusion_criteria'] = np.where(
			(df['exclusion_status'] == 1) & (df['gi_number'].isin(query_acc_status.keys())),
			'excluded due to no hit', 
			''
			)

		df.to_csv(gB_matrix_file, sep='\t', index=False)

	'''

	def process(self):
		print(f'Using {self.gb_matrix} GenBank matrix file')
		self.hydrate_update_reference_assets()
		self._require_file(self.query_fasta, "Query FASTA")
		self._require_file(self.db_fasta, "Reference FASTA")
		self._require_file(self.gb_matrix, "GenBank matrix")
		self._validate_tsv_columns(self.gb_matrix, ["gi_number"], "GenBank matrix")
		query_ids = self._fasta_ids(self.query_fasta)
		if not query_ids:
			raise ValueError(f"Query FASTA has no sequences: {self.query_fasta}")
		db_ids = self._fasta_ids(self.db_fasta)
		if not db_ids:
			raise ValueError(f"Reference FASTA has no sequences: {self.db_fasta}")
		if self.is_segmented_virus == 'Y':
			self._require_file(self.segment_file, "Segment annotation")
		if os.path.isfile(self.master_acc):
			self._require_file(self.master_acc, "Master accession list")

		master_list = self.get_master_list()
		if master_list:
			missing_masters = [acc for acc in master_list if acc not in db_ids]
			if missing_masters:
				raise ValueError(
					f"Master reference accession(s) not found in reference FASTA '{self.db_fasta}': {missing_masters}"
				)

		if not self.check_blast_exists("blastn"):
			raise RuntimeError("blastn executable not available")
		if not self.check_blast_exists("makeblastdb"):
			raise RuntimeError("makeblastdb executable not available")
		if self.is_segmented_virus == 'Y' and not self.segment_file:
			raise ValueError("Missing segment file with accession and segment information")

		if self.is_segmented_virus == 'Y':
			self.run_makeblastdb(join(self.base_dir, self.output_dir))
			self.run_blastn(join(self.base_dir, self.output_dir), self.query_fasta)
			for each_segment_dir in [join(self.base_dir, self.output_dir, "segment_sorted"), join(self.base_dir, self.output_dir, "segment_sorted_all"), join(self.base_dir, self.output_dir, "segment_merged_fasta")]:
				self.delete_directory(each_segment_dir)

			self.process_segment_virus(
				join(self.base_dir, self.output_dir, "query_tophits.tsv"),
				join(self.base_dir, self.output_dir, "query_uniq_tophits.tsv"),
				self.segment_file,
				join(self.base_dir, self.output_dir, "query_uniq_tophit_annotated.tsv")
			)
			# Mark excluded sequences in the GenBank matrix (including exclusion_list hits)
			self.update_gB_matrix(self.query_fasta, join(self.base_dir, self.output_dir, "query_uniq_tophits.tsv"), self.gb_matrix)
			self.write_filtered_ref_fasta(join(self.base_dir, self.output_dir), self.get_exclusion_list_refs())
		else:
			if self.is_update == 'Y':	
				blast_tmp_dir = join(self.base_dir, self.output_dir, "tmp_dir")
				os.makedirs(blast_tmp_dir, exist_ok=True)
				self.update(blast_tmp_dir)
				self.run_makeblastdb(blast_tmp_dir)
				self.run_blastn(blast_tmp_dir, join(blast_tmp_dir, "query.fa"))
				self.process_non_segmented_virus(blast_tmp_dir, join(blast_tmp_dir, "query.fa"))
				self.update_gB_matrix(join(blast_tmp_dir, "query.fa"), join(blast_tmp_dir, "query_uniq_tophits.tsv"), self.gb_matrix)
				if self.keep_blast_tmp_dir == 'N':
					shutil.rmtree(blast_tmp_dir)
			else:
				self.run_makeblastdb(join(self.base_dir, self.output_dir))
				self.run_blastn(join(self.base_dir, self.output_dir), self.query_fasta)
				self.process_non_segmented_virus(join(self.base_dir, self.output_dir), self.query_fasta)
				self.update_gB_matrix(self.query_fasta, join(self.base_dir, self.output_dir, "query_uniq_tophits.tsv"), self.gb_matrix)
				self.write_filtered_ref_fasta(join(self.base_dir, self.output_dir), self.get_exclusion_list_refs())


if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the BLAST alignment of query sequences against the given reference sequences')
	parser.add_argument('-q', '--query_fa', help='query fasta file', default="tmp/Sequences/query_seq.fa")
	parser.add_argument('-r', '--ref_fa', help='Blast DB fasta file. (Note: program will consider this file as db file and create the blast index for the given file)', default="tmp/Sequences/ref_seq.fa")
	parser.add_argument('-b', '--base_dir', help='Base directory', default="tmp")
	parser.add_argument('-t', '--output_dir', help='Output directory', default="Blast")
	parser.add_argument('-o', '--output_file', help='output file', default='query_tophits.tsv')
	parser.add_argument('-s', '--is_segmented_virus', help='Type Y for segmented virus else N', default='N')
	parser.add_argument('-f', '--segment_file', help='File containing information about the segments')
	parser.add_argument('-m', '--master_acc', help='Master accession. Example Rabies Virus uses NC_001542 as master reference', default=None)
	parser.add_argument('-u', '--is_update', help='If you have new downloaded sequence to blast then use this option, it will avoid performing blast on existing sequences', default='N')
	parser.add_argument('-k', '--keep_blast_tmp_dir', help='Retains the blast temp directory for debug purpose', default='N')
	parser.add_argument('-g', '--gb_matrix', help='GenBank matrix file', default='tmp/GenBank-matrix/gB_matrix_raw.tsv')
	parser.add_argument('--update_db', help='Existing SQLite DB used as the source of truth for reference accessions and sequences in update mode', default=None)
	parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
	args = parser.parse_args()

	if not args.master_acc and not args.update_db:
		parser.error('--master_acc is required unless --update_db is provided')
	if args.is_segmented_virus == 'Y' and not args.segment_file and not args.update_db:
		parser.error('--segment_file is required for segmented viruses unless --update_db is provided')

	processor = BlastAlignment(
		args.query_fa,
		args.ref_fa,
		args.base_dir,
		args.output_dir,
		args.output_file,
		args.is_segmented_virus,
		args.master_acc,
		args.is_update,
		args.keep_blast_tmp_dir,
		args.gb_matrix,
		args.segment_file,
		args.update_db,
		threads=args.threads
		)
	try:
		processor.process()
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		sys.exit(2)
