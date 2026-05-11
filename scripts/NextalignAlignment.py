import os
import re
import csv
import read_file
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join
from argparse import ArgumentParser
from TextFileHandler import TextFileLoader
from FastaHandler import RemoveRedundantSequence

class NextalignAlignment:
	def __init__(self, gb_matrix, query_dir, ref_dir, ref_fa_file, master_seq_dir, tmp_dir, master_ref, nextalign_dir, reference_alignment):
		self.gb_matrix = gb_matrix
		self.query_dir = query_dir
		self.ref_dir = ref_dir
		self.ref_fa_file = ref_fa_file
		self.master_seq_dir = master_seq_dir
		self.master_ref = master_ref
		self.tmp_dir = tmp_dir
		self.min_seed = "44"
		self.seed_spacing = "50"
		self.min_match_rate = "0.1"
		self.reference_alignment = reference_alignment
		self.nextalign_dir = nextalign_dir
	
	@staticmethod
	def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]

	def get_master_list(self):
		if os.path.isfile(self.master_ref):
			try:
				df = pd.read_csv(self.master_ref, sep='\t', header=None, dtype=str)
				if df.shape[1] >= 2:
					if df[1].str.lower().eq('master').any():
						masters = df[df[1].str.lower() == 'master']
						return masters[0].tolist()
				return df[0].tolist()
			except:
				return []
		else:
			return [x.strip() for x in self.master_ref.split(',') if x.strip()]

	def nextalign_master(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(ref_acc_path)
		command = [
			'nextalign', 'run',
			'--min-seeds', f'{self.min_seed}',
			'--seed-spacing', f'{self.seed_spacing}',
			'--min-match-rate', f'{self.min_match_rate}',
			'--input-ref', ref_acc_path,
			'--output-all', join(query_aln_op, f'{accession}'),
			'--output-basename', f'{accession}',
			'--include-reference',
			query_acc_path
			]

		command_str = " ".join(command)
		print(f"Executing command: {command_str}")

		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{accession} completed successfully.")
		else:
			print(f"{accession} failed with return code {return_code}")

	def nextalign_query(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(query_acc_path)
		command = [
			'nextalign', 'run',
			'--min-seeds', f'{self.min_seed}',
			'--seed-spacing', f'{self.seed_spacing}',
			'--min-match-rate', f'{self.min_match_rate}',
			'--input-ref', ref_acc_path,
			'--output-all', join(query_aln_op, f'{accession}'),
			'--output-basename', f'{accession}',
			'--include-reference',
			query_acc_path
		]
        
		command_str = " ".join(command)
        
		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{accession} completed successfully.")
		else:
			print(f"{accession} failed with return code {return_code}")

	def update_gb_matrix(self, alignment_dir: list, gB_matrix_file):
		failed_accessions = {}

		for each_aln_type in alignment_dir:
			for each_aln in os.listdir(each_aln_type):
				error_file_path = join(each_aln_type, each_aln, each_aln + ".errors.csv")
				if os.path.exists(error_file_path):
					with open(error_file_path) as f:
						for line in f:
							if "In sequence" in line:
								acc = line.split(",")[0].strip()
								error = line.split(",")[1].strip()
								failed_accessions.setdefault(acc, []).append(error)

		updated_rows = []
		with open(gB_matrix_file, newline='') as csvfile:
			reader = csv.DictReader(csvfile, delimiter='\t')
			fieldnames = reader.fieldnames

			if 'exclusion_status' not in fieldnames:
				fieldnames.append('exclusion_status')
			if 'exclusion_criteria' not in fieldnames:
				fieldnames.append('exclusion_criteria')

			for row in reader:
				gi = row.get('gi_number')
				if gi in failed_accessions:
					row['exclusion_status'] = '1'
					existing_criteria = row.get('exclusion_criteria', '')
					new_criteria = '; '.join(failed_accessions[gi])
					if existing_criteria:
						row['exclusion_criteria'] = existing_criteria + '; ' + new_criteria
					else:
						row['exclusion_criteria'] = new_criteria
				else:
					row['exclusion_status'] = row.get('exclusion_status', '0')
					row['exclusion_criteria'] = row.get('exclusion_criteria', '')
				updated_rows.append(row)

		with open(gB_matrix_file, 'w', newline='') as csvfile:
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
			writer.writeheader()
			writer.writerows(updated_rows)

	def process(self):
		query_aln_output_dir = join(self.tmp_dir, self.nextalign_dir, "query_aln")
		ref_aln_output_dir = join(self.tmp_dir, self.nextalign_dir, "reference_aln")
		
		os.makedirs(query_aln_output_dir, exist_ok=True)
		os.makedirs(ref_aln_output_dir, exist_ok=True)

		#query_table = join(self.table_dir, "query_features.tsv")
		
		if self.reference_alignment:
			#self.master_feature_table(self.reference_alignment)
			for each_query_file in os.listdir(self.query_dir):
				ref_file = each_query_file
				self.nextalign_query(
					join(self.query_dir, each_query_file),
					join(self.ref_dir, ref_file), 
					query_aln_output_dir
				)
			self.update_gb_matrix([query_aln_output_dir], self.gb_matrix)
		
		else:
			#align query against reference sequence alignment
			for each_query_file in os.listdir(self.query_dir):
				ref_file = each_query_file
				self.nextalign_query(join(self.query_dir, each_query_file),join(self.ref_dir, ref_file), query_aln_output_dir)
		
			#align master and reference sequence alignment
			for each_ref in os.listdir(self.master_seq_dir):
				self.nextalign_master(self.ref_fa_file,join(self.master_seq_dir, each_ref),ref_aln_output_dir)

			masters = self.get_master_list()
			for master in masters:
				input_seq = join(ref_aln_output_dir, master, master + ".aligned.fasta")
				if os.path.exists(input_seq):
					output_seq = input_seq
					unique_seqs = RemoveRedundantSequence(input_seq, output_seq)
					unique_seqs.remove_redundant_fasta()

			self.update_gb_matrix([query_aln_output_dir, ref_aln_output_dir], self.gb_matrix)

if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the nextalign of each sequence')
	parser.add_argument('-g', '--gB_matrix', help='GenBank matrix (meta data) file.', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-q', '--query_dir', help='Query file directory.', default="tmp/Blast/grouped_fasta")
	parser.add_argument('-r', '--ref_dir', help='Reference fasta directory', default="tmp/Blast/ref_seqs")
	parser.add_argument('-f', '--ref_fa_file', help='Reference fasta sequences', default="tmp/Sequences/ref_seq.fa")
	parser.add_argument('-ms', '--master_seq_dir', help='Master sequence directory', default="tmp/Blast/master_seq")
	parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp")
	parser.add_argument('-m', '--master_ref', help='Master reference accession. Generally, the Ref Seq accession. In case of Rabies it is NC_001542', required=True)
	parser.add_argument('-n', '--nextalign_dir', help='Nextalign output to be saved', default="Nextalign")
	parser.add_argument('-ra', '--ref_alignment_file', help='Use your own reference alignment file instead of Nextalign perfoms the alignment of reference against the master reference sequence')
	args = parser.parse_args()

	processor = NextalignAlignment(args.gB_matrix, args.query_dir, args.ref_dir, args.ref_fa_file, args.master_seq_dir, args.tmp_dir, args.master_ref, args.nextalign_dir, args.ref_alignment_file)
	processor.process()
