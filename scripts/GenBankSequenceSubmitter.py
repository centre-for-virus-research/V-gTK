#python GenBankSequenceSubmitter-1.py -q ./../test_data/example-2/fasta-seq/ -m ./../test_data/example-2/metadata.tsv -n ./../test_data/example-2/template.sbt -gff NC_001542.gff3 -db ./../../test_version/TING/tmp-rabv/SqliteDB/rabv_apr0825.db

import os
import re
import csv
import uuid
import glob
import shutil
import sqlite3
import requests
import read_file
import subprocess
from Bio import SeqIO
from time import sleep
from Bio.Seq import Seq
from os.path import join
from argparse import ArgumentParser
#from PadAlignment import PadAlignment
from GffToDictionary import GffDictionary
#from FeatureCalculator  import FeatureCordCalculator 


class GenBankSequenceSubmitter:
	def __init__(self, sequence_dir, tmp_dir, output_dir, metadata, ncbi_submission_template, gff_file, db, gaps_to_ignore):
		self.sequence_dir = sequence_dir
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.metadata = metadata
		self.ncbi_submission_template = ncbi_submission_template
		self.min_seed = "44"
		self.seed_spacing = "50"
		self.min_match_rate = "0.1"
		self.gff_file = gff_file
		self.db = db
		self.gaps_to_ignore = gaps_to_ignore
		#self.feature_cords = FeatureCordCalculator()
		self.analysis_dir = str(uuid.uuid4())

	@staticmethod
	def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]

	def table2asn(self):
		os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)
		command = [
			'table2asn',
			'-indir', join(self.tmp_dir, self.output_dir, "tmp"),
			'-t', self.ncbi_template
		]
		try:
			subprocess.run(command, check=True)
			print(f"Table2asn ran successfully on {self.squence_dir}")
		except subprocess.CalledProcessError as e:
			print(f"Error running Table2asn: {e}")

	def run_makeblastdb(self, db_fasta):
		db_file_name = db_fasta
		command = [
			'makeblastdb',
			'-in', db_fasta,
			'-out', join(db_file_name),
			'-title', "alignment",
			'-dbtype', 'nucl'
		]
		try:
			subprocess.run(command, check=True)
			print(f"makeblastdb ran successfully on {db_file_name}")
		except subprocess.CalledProcessError as e:
			print(f"Error running makeblastdb: {e}")

	def blastn(self, query_path, db_path):
		db_file_name = os.path.basename(db_path)
		output_file = join(self.tmp_dir, "analysis", self.analysis_dir, "query_tophits.tsv")
		command = [
			'blastn',
			'-query', query_path,
			'-db', join(self.tmp_dir, "analysis", self.analysis_dir, "DB", db_file_name),
			'-task', 'blastn',
			'-max_target_seqs', '1',
			'-max_hsps', '1',
			'-out', output_file,
			'-outfmt', "6 qacc sacc pident sstrand"
		]
		try:
			subprocess.run(command, check=True)
			print(f"blastn ran successfully. Results saved in {output_file}")
		except subprocess.CalledProcessError as e:
			print(f"Error running blastn: {e}")

	def metadata_header(self):
		headers = ["SeqID", "organism", "genotype", "host", "isolate", "isolation_source", "collection-date", "note", "country", "note"]
		#headers = ["FastaID", "sample_id", "isolation_source", "country" ,"place_sampled" ,"host_species" ,"collection_year","notes", "collection_date"]
		return headers

	def load_metadata(self):
		headers = self.metadata_header() 
		df = {}
		with open(self.metadata, mode="r", encoding="utf-8") as file:
			reader = csv.DictReader(file, delimiter="\t")

			missing_headers = [col for col in headers if col not in reader.fieldnames]
			if missing_headers:
				raise ValueError(f"Missing columns in header: {', '.join(missing_headers)}")

			for row in reader:
				seq_id = row["SeqID"]
				if seq_id in df:
					raise ValueError(f"Duplicate SeqID found: {SeqID}")

				df[seq_id] = {key: row[key] for key in headers if key in row}

		return df

	def load_fasta(self):
		fasta_dict = {}

		for filename in os.listdir(self.sequence_dir):
			if filename.endswith(".fasta") or filename.endswith(".fa"):
				filepath = os.path.join(self.sequence_dir, filename)
				with open(filepath, "r") as file:
					for record in SeqIO.parse(file, "fasta"):
						if record.id in fasta_dict:
							raise ValueError(f"Duplicate accession found: {accession}. Ensure unique accessions.")
						fasta_dict[record.id] = str(record.seq)

		return fasta_dict

	def data_prep(self):
		os.makedirs(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp"), exist_ok=True)
		header = self.metadata_header()
		meta_data = self.load_metadata()
		sequence_obj = self.load_fasta()
		
		missing_in_fasta = [seq_id for seq_id in sequence_obj.keys() if seq_id not in meta_data.keys()]
		if len(missing_in_fasta) < 1:
			write_query_seq = open(join(self.tmp_dir, "analysis", self.analysis_dir, "query.fa"), 'w')
			#creating meta-data file
			for each_acc, data in meta_data.items():
				write_file = open(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp", each_acc + ".src"), "w")
				write_file.write("\t".join(header))
				write_file.write("\n")
				write_file.write("\t".join(data.values()))
				write_file.close()

			#creating fasta file
			for acc, seq in sequence_obj.items():
				write_query_seq.write(">" + acc + "\n" + seq + "\n")
				write_file = open(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp", acc + ".fsa"), "w")
				write_file.write(">" + acc)
				write_file.write("\n")
				write_file.write(seq)
				write_file.write("\n")
				write_file.close()
			write_query_seq.close()
		else:
			raise ValueError(f"Following sequences accssions are not found on meta data {missing_in_fasta}")		

	def extract_ref_seq(self):
		db_path = join(self.tmp_dir, "analysis", self.analysis_dir, "DB")
		os.makedirs(db_path, exist_ok=True)
		write_file = open(join(db_path, "db.fa"), 'w')
		conn = sqlite3.connect(self.db)
		cursor = conn.cursor()

		cursor.execute("SELECT primary_accession, accession_type FROM meta_data where accession_type='reference' OR accession_type = 'master'")
		ref_accs = cursor.fetchall()
		ref_accs_list = [(item[0], item[1]) for item in ref_accs]
		for each_acc, accession_type in ref_accs_list:
			cursor.execute("SELECT sequence FROM sequences WHERE header = ?", (each_acc,))
			if each_acc == "NC_001542": accession_type="master"
			result = cursor.fetchone()
			if result:
				sequence = result[0]
				acc_type = "|" + accession_type
				write_file.write(">" + each_acc + acc_type)
				write_file.write("\n")
				write_file.write(sequence)
				write_file.write("\n")
		write_file.close()

	def blast_analysis(self):

		values = {}
		records = {}

		query_path = join(self.tmp_dir, "analysis", self.analysis_dir, "query.fa")
		ref_path = join(self.tmp_dir, "analysis", self.analysis_dir, "ref.fa")
		db_path = join(self.tmp_dir, "analysis", self.analysis_dir, "DB", "db.fa")
		query_tophits = join(self.tmp_dir, "analysis", self.analysis_dir, "query_tophits.tsv")
		query_tophits_uniq = join(self.tmp_dir, "analysis", self.analysis_dir, "query_tophits_uniq.tsv")
		sorted_fasta = join(self.tmp_dir, "analysis", self.analysis_dir, "sorted_fasta")
		merged_fasta = join(self.tmp_dir, "analysis", self.analysis_dir, "merged_fasta")
		sorted_all = join(self.tmp_dir, "analysis", self.analysis_dir, "sorted_all")
		grouped_fasta = join(self.tmp_dir, "analysis", self.analysis_dir, "grouped_fasta")
		ref_seq_dir = join(self.tmp_dir, "analysis", self.analysis_dir, "reference_sequences")
		master_seq_dir = join(self.tmp_dir, 'analysis', self.analysis_dir, 'master_sequences')
		master_and_reference_merged = join(self.tmp_dir, 'analysis', self.analysis_dir, 'master_and_reference_merged')

		self.run_makeblastdb(db_path)
		self.blastn(query_path, db_path)

		with open(query_tophits, newline='') as file:
			reader = csv.reader(file, delimiter='\t')
			for row in reader:
				#query, ref, identity, strand,
				col1, col2, col3, col4 = row[0], row[1].split("|")[0], float(row[2]), row[3]
				if col1 in values:
					existing_value = values[col1]
					if col3 > existing_value:
						records[col1] = [col1, col2, col3, col4]
						values[col1] = col3
				else:
					records[col1] = [col1, col2, col3, col4]
					values[col1] = col3

		with open(query_tophits_uniq, 'w', newline='') as file:
			writer = csv.writer(file, delimiter='\t')
			for record in records.values():
				writer.writerow(record)

		seq_dicts = {}
		query_seqs = read_file.fasta(query_path)
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		#Seperate plus and minus strand sequences
		for each_line in open(query_tophits_uniq, 'r'):
			query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
			if strand == "plus":
				with open(join(sorted_fasta, "plus.fa"), 'a') as file_plus:
					file_plus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")
			else:
				with open(join(sorted_fasta, "minus.fa"), 'a') as file_minus:
					file_minus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")

		for each_file in os.listdir(sorted_fasta):
			if "minus" in each_file:
				command = ["seqkit", "seq", "-r", "-p", "-v", "-t", "dna", join(sorted_fasta, each_file), ">", join(merged_fasta, each_file)]
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
			print(f"concatenation sucessful")
		except subprocess.CalledProcessError as e:
			print(f"Error in concatenation: {e}")

		grouped_dict = {}
		for each_line in open(query_tophits_uniq, 'r'):
			query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
			if ref_acc not in grouped_dict:
				grouped_dict[ref_acc] = [query_acc]
			else:
				grouped_dict[ref_acc].append(query_acc)

		seq_dicts = {}
		master_seq_dict = {}
		query_seqs = read_file.fasta(join(sorted_all, "query_seq.fa"))
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		for each_ref_acc, list_of_query_acc in grouped_dict.items():
			with open(join(grouped_fasta, each_ref_acc + '.fasta'), 'a') as write_file:
				for each_query_acc in list_of_query_acc:
					seqs = seq_dicts[each_query_acc]
					write_file.write(">" + each_query_acc + '\n')
					for i in range(0, len(seqs), 80):
						write_file.write(seqs[i:i + 80] + '\n')
	
		ref_seqs = read_file.fasta(db_path)
		for rows in ref_seqs:
			accs, accs_type = rows[0].split("|")
			seq_dicts[accs] = rows[1].strip()
			if accs_type == "master":
				master_seq_dict[accs] = rows[1].strip()

		for ref_acc_info in grouped_dict.keys():
			if "|" in ref_acc_info:
				split_info = ref_acc_info.split("|")
				ref_accs = split_info[0]
				acc_type = split_info[1]
				seqs = seq_dicts[ref_accs]
			else:
				ref_accs = ref_acc_info
				acc_type = ""
				seqs = seq_dicts[ref_acc_info]

			with open(join(ref_seq_dir, ref_accs + '.fasta'), 'w') as write_file:
				write_file.write(">" + ref_accs + '\n')
				for i in range(0, len(seqs), 80):
					write_file.write(seqs[i:i + 80] + '\n')

		#writing master sequences
		master_fasta = ""
		for master_acc, seqs in master_seq_dict.items():
			master_fasta = join(master_seq_dir, master_acc + '.fasta')
			with open(join(master_seq_dir, master_acc + '.fasta'), 'w') as write_file:
				write_file.write(">" + master_acc + '\n')
				for i in range(0, len(seqs), 80):
					write_file.write(seqs[i:i + 80] + '\n')

		
		#merge reference sequence
		write_file = open(ref_path, 'w')
		for each_ref_fa in os.listdir(ref_seq_dir):
			df = read_file.fasta(join(ref_seq_dir, each_ref_fa))
			for rows in df:
				write_file.write(">" + rows[0].strip() + "\n" + rows[1].strip() + "\n")
		write_file.close()

	
		for each_ref in os.listdir(ref_seq_dir):
			os.system('cat ' + master_fasta + ' ' + f'{join(ref_seq_dir, each_ref)}' + ' >' f'{join(master_and_reference_merged, each_ref)}') 
			print(f"Merging complete {join(master_and_reference_merged, each_ref)}") 

	def mafft_ref_sequences(self, ref_acc_file, output_dir):
		output_file = self.path_to_basename(ref_acc_file) + "_aln.fasta"
		output_path = join(output_dir, output_file)
		try:
			with open(output_path, 'w') as out_f:
				subprocess.run(['mafft', ref_acc_file], stdout=out_f, check=True)
				print(f"mafft ran successfully on {output_file}")
		except subprocess.CalledProcessError as e:
			print(f"Error running MAFFT: {e}")

	def mafft_query_sequences(self, query_file, ref_alignment_file, output_dir):
		output_file = self.path_to_basename(ref_alignment_file) + "_aln.fasta"
		output_path = join(output_dir, output_file)
		ref_alignment_file = ref_alignment_file.split('.')[0] + "_aln.fasta"
		try:
			with open(output_path, 'w') as out_f:
				subprocess.run(['mafft', '--add', query_file, '--keeplength', ref_alignment_file], stdout=out_f, check=True)
				print(f"mafft ran successfully on {output_file}")
		except subprocess.CalledProcessError as e:
			print(f"Error running MAFFT: {e}")

	def extract_matching_sequences(self, fasta_dir, output_file=None):
		results = {}

		for file in os.listdir(fasta_dir):
			if file.endswith("_aln.fasta"):
				accession = file.split("_")[0]
				file_path = os.path.join(fasta_dir, file)

				for record in SeqIO.parse(file_path, "fasta"):
					if accession in record.id:
						results[accession] = str(record.seq)
						break

		if output_file:
			with open(output_file, "w") as out_f:
				for acc, seq in results.items():
					out_f.write(f">{acc}\n{seq}\n")

	def find_cds_for_coordinates(self, gff_dict, query_start, query_end):
		matching_cds = []
		for feature, cds_list in gff_dict.items():
			if feature=='CDS':
				for cds in cds_list:
					cds_start = int(cds['start'])
					cds_end = int(cds['end'])
					if query_end >= cds_start and query_start <= cds_end:
						overlap_start = max(query_start, cds_start)
						overlap_end = min(query_end, cds_end)
						matching_cds.append({
							'start': overlap_start,
							'end': overlap_end,
							'product': cds['product']
					})

		return matching_cds

	def table2asn_coordinates(self, gff_dict, start, end, tolerance=10):
		for cds_entry in gff_dict['CDS']:
			cds_start = int(cds_entry['start'])
			cds_end = int(cds_entry['end'])
			product = cds_entry['product']

			# Check if given coordinates fall within tolerance range of the actual CDS
			if (cds_start - tolerance) <= start <= (cds_start + tolerance) and \
				(cds_end - tolerance) <= end <= (cds_end + tolerance):

				# Determine 5' and 3' partial status
				if start < cds_start:
					start_str = f"<{start}"  # 5' partial
				else:
					start_str = str(start)

				if end > cds_end:
					end_str = f">{end}"  # 3' partial
				else:
					end_str = str(end)

				return [start_str, end_str, product]

		return None  # No match within tolerance

	def extract_gene_name(self, protein_name):
		words = protein_name.split()
		return words[-1] if words[-1].isalpha() else words[-2]

	def translate_six_frames(self, seq):
		frames = [seq[i:].translate(to_stop=False) for i in range(3)]
		rev_seq = seq.reverse_complement()
		frames += [rev_seq[i:].translate(to_stop=False) for i in range(3)]
		return frames

	def find_orfs(self, protein_seq):
		orfs = []
		seq_str = str(protein_seq)
		start = 0
		while start < len(seq_str):
			start = seq_str.find("M", start)
			if start == -1:
				break
			end = seq_str.find("*", start)
			if end != -1:
				orfs.append(seq_str[start:end+1])
			start += 1
		return orfs

	def select_best_orf(self, orfs):
		return max(orfs, key=len) if orfs else None


	def recalculate_cds_cords(self, dna_sequence: Seq, start: int, end: int):
		best_orf = None
		best_length = 0
		dna_seq = Seq(dna_sequence.replace("-", ""))
		for frame in self.translate_six_frames(Seq(dna_seq)):
			orfs = self.find_orfs(frame)
			longest_orf = self.select_best_orf(orfs)
			if longest_orf and len(longest_orf) > best_length:
				best_orf = longest_orf
				best_length = len(longest_orf)

		if best_orf:
			recalculated_end = start + (best_length * 3) - 1
			return start, recalculated_end
		else:
			return start, end

	def get_gap_ranges(self, sequence):
		gap_ranges = []
		start = None

		for i, char in enumerate(sequence):
			if char == '-':
				if start is None:
					start = i + 1  # Convert to 1-based indexing
			else:
				if start is not None:
					gap_ranges.append([start, i])
					start = None

		if start is not None:
			gap_ranges.append([start, len(sequence)])

		return gap_ranges

	def count_gaps_before_position(self, gap_ranges, position):
		"""Count how many positions are removed before a given alignment position."""
		count = 0
		for start, end in gap_ranges:
			if end < position:
				count += (end - start + 1)
			elif start <= position <= end:
				count += (position - start + 1)
		return count

	def recalculate_cds_coordinates(self, sequence_id, gap_ranges, cds_list, start_offset):
		adjusted_coords = []
		for cds in cds_list:
			cds_start = int(cds['start'])
			cds_end = int(cds['end'])

			gaps_before_start = self.count_gaps_before_position(gap_ranges, cds_start)
			gaps_before_end = self.count_gaps_before_position(gap_ranges, cds_end)

			adj_start = cds_start - gaps_before_start
			adj_end = cds_end - gaps_before_end

			# Apply start offset (e.g., 71 if first gap ends at 70)
			#adj_start += start_offset
			#adj_end += start_offset
			adj_start = cds_start - gaps_before_start + (start_offset - 1)
			adj_end = cds_end - gaps_before_end + (start_offset - 1)


			adjusted_coords.append([adj_start, adj_end])
		return adjusted_coords

	def find_gaps_in_fasta(self, fasta_file, gff_file):
		gff_dict = GffDictionary(gff_file).gff_dict
		cds_list = gff_dict['CDS']

		data = []

		unique_acc_list = []
		meta_data = self.load_metadata()

		for record in SeqIO.parse(fasta_file, "fasta"):
			if record.id in meta_data:
				if record.id not in unique_acc_list:
					unique_acc_list.append(record.id)

					sequence = str(record.seq)
					gaps = self.get_gap_ranges(sequence)

					# Calculate start offset: just after the first gap
					if gaps and gaps[0][0] == 1:
						start_offset = gaps[0][1] + 1
					else:
						start_offset = 1

					adjusted = self.recalculate_cds_coordinates(record.id, gaps, cds_list, start_offset)

					#print(f">{record.id}")
					#print(adjusted)
					data.append([record.id, adjusted, sequence])
		return data

	# Example usage:
	#find_gaps_in_fasta("query_kx_aligned.fa", "NC_001542.gff3")

	def table2asn(self, input_dir, ncbi_template):
		command = [
			'table2asn' ,
			'-indir', input_dir,
			'-t', ncbi_template
		]

		command_str = " ".join(command)
		print(f"Executing command: {command_str}")

		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{input_dir} completed successfully.")
		else:
			print(f"{input_dir} failed with return code {return_code}")

	def run_table2asn(self):
		self.table2asn(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp"), self.ncbi_submission_template)

	def move_sqn_to_output_dir(self):
		output_dir = join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "SQN_files")
		tmp_dir = join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp")
		os.makedirs(output_dir, exist_ok=True)
		sqn_files = glob.glob(os.path.join(tmp_dir, '*.sqn'))
		if sqn_files:
			for file_path in sqn_files:
				shutil.copy(file_path, output_dir)
			print(f"GenBank SQN files copied sucessfully to location : {output_dir}")
		else:
				print("No *.SQN file exists")

	def process(self):
		os.makedirs(join(self.tmp_dir, 'analysis'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'sorted_fasta'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'merged_fasta'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'sorted_all'), exist_ok=True)
		grouped_fasta = join(self.tmp_dir, 'analysis', self.analysis_dir, 'grouped_fasta')
		os.makedirs(grouped_fasta, exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'reference_sequences'), exist_ok=True)
		query_aln_output_dir = join(self.tmp_dir, self.analysis_dir, "query_aln")
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'master_sequences'), exist_ok=True)
		master_and_reference_merged = join(self.tmp_dir, "analysis", self.analysis_dir, "master_and_reference_merged")
		mafft_reference_alignment = join(self.tmp_dir, "analysis", self.analysis_dir, "mafft_reference_alignment")
		reference_alignments =  join(self.tmp_dir, "analysis", self.analysis_dir, "reference_alignments")
		query_ref_alignment = join(self.tmp_dir, "analysis", self.analysis_dir, "query_ref_alignment")
		table2asn_tmp = join(self.tmp_dir, "analysis", self.analysis_dir, "Table2asn", "tmp")
		create_tmp_dir = os.makedirs(table2asn_tmp, exist_ok=True)

		os.makedirs(master_and_reference_merged, exist_ok=True)
		os.makedirs(mafft_reference_alignment, exist_ok=True)
		os.makedirs(reference_alignments, exist_ok=True)
		os.makedirs(query_ref_alignment, exist_ok=True)
					
		self.data_prep()
		self.extract_ref_seq()
		self.blast_analysis() 

		for each_ref_merged_acc in os.listdir(master_and_reference_merged):
			self.mafft_ref_sequences(join(master_and_reference_merged, each_ref_merged_acc), mafft_reference_alignment)

		for each_ref_aln in os.listdir(mafft_reference_alignment):
			output_file = self.path_to_basename(each_ref_aln)
			output_file = join(reference_alignments, output_file + '.fasta')
			self.extract_matching_sequences(join(mafft_reference_alignment), output_file)

		for each_query_seq in os.listdir(grouped_fasta):
			reference_file = self.path_to_basename(each_query_seq)
			self.mafft_query_sequences(join(grouped_fasta, each_query_seq), join(reference_alignments, reference_file), query_ref_alignment)

		meta_data = self.load_metadata()
		gff_dict = GffDictionary(self.gff_file).gff_dict

		non_redundant_acc = []
		for each_query_ref_aln in os.listdir(query_ref_alignment):
			partial_gaps = self.find_gaps_in_fasta(join(query_ref_alignment, each_query_ref_aln), self.gff_file)
			for each_acc in partial_gaps:
			
		
				acc_num = each_acc[0]
				coordinates = each_acc[1]

				with open(join(table2asn_tmp, acc_num + ".tbl"), "w") as out:
					out.write(f">Feature {acc_num}\n")

					for each_cords in coordinates:
						cord_start = each_cords[0]
						cord_end = each_cords[1]
						table2asn_cords = self.table2asn_coordinates(gff_dict, cord_start, cord_end)

						start, end, product = table2asn_cords

						start_clean = start.lstrip("<>")
						end_clean = end.lstrip("<>")

						gene_symbol = product.split()[-1]

						product_name = product.replace(gene_symbol, "").strip(" ,")

						out.write(f"{start_clean}\t{end_clean}\tgene\n")
						out.write(f"\t\t\tgene\t{gene_symbol}\n")

						out.write(f"{start}\t{end}\tCDS\n")
						out.write(f"\t\t\tproduct\t{product_name}\n")
						out.write(f"\t\t\tgene\t{gene_symbol}\n")

		self.run_table2asn()
		self.move_sqn_to_output_dir()
	
if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the nextalign of each sequence')
	parser.add_argument('-q', '--sequence_dir', help='Sequence file directory, it can be single or multiple fasta sequencce files.', required=True)
	parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp")
	parser.add_argument('-o', '--output_dir', help='Output directory where processed data and results are stored', default='Table2asn')
	parser.add_argument('-m', '--metadata', help='A tab limited meta data file name', required=True)
	parser.add_argument('-n', '--ncbi_submission_template', help='NCBI submission template file generated from "https://submit.ncbi.nlm.nih.gov/genbank/template/submission/"', required=True)
	parser.add_argument('-gff', '--gff_file', help='GFF3 file containing CDS information. This file will be used for assigning the coding regions for the input sequence files', required=True)
	parser.add_argument('-db', '--vgtk-db', help='VGTK Database', required=True)
	parser.add_argument('-gp', '--gaps_to_ignore', help='Lenght of gaps to ignore', default=30)
	
	args = parser.parse_args()

	processor = GenBankSequenceSubmitter(args.sequence_dir, args.tmp_dir, args.output_dir, args.metadata, args.ncbi_submission_template, args.gff_file, args.vgtk_db, args.gaps_to_ignore)
	processor.process()
	
