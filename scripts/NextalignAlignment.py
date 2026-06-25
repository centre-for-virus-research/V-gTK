import os
import re
import csv
import sys
import read_file
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join
from argparse import ArgumentParser
from TextFileHandler import TextFileLoader
from FastaHandler import RemoveRedundantSequence
from ExportRefListFromUpdateDb import load_master_accessions, load_master_accessions_from_file

class NextalignAlignment:
	def __init__(self, gb_matrix, query_dir, ref_dir, ref_fa_file, master_seq_dir, tmp_dir, master_ref, nextalign_dir, reference_alignment, update_db=None):
		self.gb_matrix = gb_matrix
		self.query_dir = query_dir
		self.ref_dir = ref_dir
		self.ref_fa_file = ref_fa_file
		self.master_seq_dir = master_seq_dir
		self.master_ref = master_ref
		self.update_db = update_db
		self.tmp_dir = tmp_dir
		self.reference_alignment = reference_alignment
		self.nextalign_dir = nextalign_dir
		
		# Cascading parameter profiles for dynamic relaxation
		self.relaxation_profiles = [
			{"min_seeds": 44, "seed_spacing": 50, "min_match_rate": 0.1},
			{"min_seeds": 30, "seed_spacing": 40, "min_match_rate": 0.1},
			{"min_seeds": 20, "seed_spacing": 30, "min_match_rate": 0.05},
			{"min_seeds": 10, "seed_spacing": 20, "min_match_rate": 0.05},
			{"min_seeds": 5,  "seed_spacing": 15, "min_match_rate": 0.01}
		]
	
	@staticmethod
	def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]

	def get_master_list(self):
		if self.update_db:
			masters = load_master_accessions(self.update_db)
			if masters:
				return masters
		if os.path.isfile(self.master_ref):
			try:
				return load_master_accessions_from_file(self.master_ref)
			except:
				return []
		else:
			return [x.strip() for x in self.master_ref.split(',') if x.strip()]

	def _validate_alignment_integrity(self, aligned_fasta_path):
		"""
		Evaluates the structural validity of the alignment by checking nucleotide
		identity in the final 30% of the alignment window (conserved downstream genes)
		for all query sequences in the alignment against the reference sequence.
		Returns (bool_passed, min_calculated_identity).
		"""
		if not os.path.exists(aligned_fasta_path):
			return False, 0.0
			
		records = list(SeqIO.parse(aligned_fasta_path, "fasta"))
		if len(records) < 2:
			return False, 0.0
			
		ref_seq = str(records[0].seq).upper()
		total_columns = len(ref_seq)
		start_idx = int(total_columns * 0.70) # Target the 3' conserved regions
		
		min_identity = 1.0
		all_passed = True
		
		for r_idx in range(1, len(records)):
			query_seq = str(records[r_idx].seq).upper()
			
			matches = 0
			valid_positions = 0
			
			for i in range(start_idx, total_columns):
				n1 = ref_seq[i]
				n2 = query_seq[i]
				if n1 == '-' and n2 == '-':
					continue
				valid_positions += 1
				if n1 == n2 and n1 != '-':
					matches += 1
					
			if valid_positions == 0:
				identity = 0.0
			else:
				identity = matches / valid_positions
				
			if identity < min_identity:
				min_identity = identity
				
			if identity < 0.50:
				all_passed = False
				print(f"[validation] Sequence {records[r_idx].id} failed frame validation against reference {records[0].id}. Downstream Identity: {identity:.2%}")
				
		return all_passed, min_identity

	def nextalign_master(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(ref_acc_path)
		output_subdir = join(query_aln_op, f'{accession}')
		aligned_fasta = join(output_subdir, f'{accession}.aligned.fasta')
		
		# Loop dynamically through relaxation parameters
		for attempt, profile in enumerate(self.relaxation_profiles, 1):
			command = [
				'nextalign', 'run',
				'--min-seeds', str(profile["min_seeds"]),
				'--seed-spacing', str(profile["seed_spacing"]),
				'--min-match-rate', str(profile["min_match_rate"]),
				'--input-ref', ref_acc_path,
				'--output-all', output_subdir,
				'--output-basename', f'{accession}',
				'--include-reference',
				query_acc_path
			]

			command_str = " ".join(command)
			print(f"Executing Master Alignment (Attempt {attempt}): {command_str}")
			
			return_code = os.system(command_str)
			
			if return_code == 0:
				passed, identity = self._validate_alignment_integrity(aligned_fasta)
				if passed:
					print(f"SUCCESS: {accession} aligned cleanly (Downstream Identity: {identity:.2%}) using seeds={profile['min_seeds']}.")
					return
				else:
					print(f"WARNING: Nextalign completed for {accession} but failed frame validation. Downstream Identity: {identity:.2%}. Frame-shift suspected near hypervariable locus.")
			else:
				print(f"WARNING: Nextalign process execution failed for {accession} on attempt {attempt}.")
				
		# If all relaxation attempts fail, crash hard to preserve database parity
		print(f"\nCRITICAL ERROR: Fundamental alignment failure for reference sequence: {accession}")
		print(f"Reason: Out-of-frame insertions/deletions or extreme sequence divergence could not be resolved across all heuristic profiles.")
		print(f"Aborting execution to protect SQLite coordinate lookups.")
		raise RuntimeError(f"Alignment Frame Failure: {accession} aborted pipeline execution.")

	def nextalign_query(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(query_acc_path)
		# Use baseline profile for intra-group queries as they are tightly clustered
		base_profile = self.relaxation_profiles[0]
		command = [
			'nextalign', 'run',
			'--min-seeds', str(base_profile["min_seeds"]),
			'--seed-spacing', str(base_profile["seed_spacing"]),
			'--min-match-rate', str(base_profile["min_match_rate"]),
			'--input-ref', ref_acc_path,
			'--output-all', join(query_aln_op, f'{accession}'),
			'--output-basename', f'{accession}',
			'--include-reference',
			query_acc_path
		]
		
		command_str = " ".join(command)
		return_code = os.system(command_str)
		if return_code != 0:
			print(f"Query alignment error: {accession} exited with code {return_code}")

	@staticmethod
	def _read_nextalign_errors(error_file_path, reference_accession):
		failed = {}
		with open(error_file_path, newline='', encoding='utf-8') as handle:
			reader = csv.DictReader(handle)
			fieldnames = reader.fieldnames or []
			if {'seqName', 'errors'}.issubset(set(fieldnames)):
				for row in reader:
					acc = str(row.get('seqName', '')).strip()
					error = str(row.get('errors', '')).strip()
					if not acc or not error or acc == reference_accession:
						continue
					failed.setdefault(acc, []).append(error)
				return failed

		handle.seek(0)
		reader = csv.reader(handle)
		for row in reader:
			if len(row) < 2:
				continue
			acc = str(row[0]).strip()
			error = str(row[1]).strip()
			if not acc or not error or acc == reference_accession:
				continue
			failed.setdefault(acc, []).append(error)
		return failed

	def update_gb_matrix(self, alignment_dir: list, gB_matrix_file):
		failed_accessions = {}

		for each_aln_type in alignment_dir:
			if os.path.basename(os.path.normpath(each_aln_type)) != "query_aln":
				continue
			for each_aln in os.listdir(each_aln_type):
				error_file_path = join(each_aln_type, each_aln, each_aln + ".errors.csv")
				if os.path.exists(error_file_path):
					for acc, errors in self._read_nextalign_errors(error_file_path, each_aln).items():
						failed_accessions.setdefault(acc, []).extend(errors)

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
				acc_type = str(row.get('accession_type', '')).strip().lower()
				if gi in failed_accessions and acc_type not in {'reference', 'master', 'exclusion_list'}:
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
		
		if self.reference_alignment:
			for each_query_file in os.listdir(self.query_dir):
				ref_file = each_query_file
				self.nextalign_query(
					join(self.query_dir, each_query_file),
					join(self.ref_dir, ref_file), 
					query_aln_output_dir
				)
			self.update_gb_matrix([query_aln_output_dir], self.gb_matrix)
		
		else:
			for each_query_file in os.listdir(self.query_dir):
				ref_file = each_query_file
				self.nextalign_query(join(self.query_dir, each_query_file), join(self.ref_dir, ref_file), query_aln_output_dir)
		
			for each_ref in os.listdir(self.master_seq_dir):
				self.nextalign_master(self.ref_fa_file, join(self.master_seq_dir, each_ref), ref_aln_output_dir)

			masters = self.get_master_list()
			for master in masters:
				input_seq = join(ref_aln_output_dir, master, master + ".aligned.fasta")
				if os.path.exists(input_seq):
					output_seq = input_seq
					unique_seqs = RemoveRedundantSequence(input_seq, output_seq)
					unique_seqs.remove_redundant_fasta()

			self.update_gb_matrix([query_aln_output_dir, ref_aln_output_dir], self.gb_matrix)

if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the nextalign of each sequence with adaptive step-down error control.')
	parser.add_argument('-g', '--gB_matrix', help='GenBank matrix file.', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-q', '--query_dir', help='Query file directory.', default="tmp/Blast/grouped_fasta")
	parser.add_argument('-r', '--ref_dir', help='Reference fasta directory', default="tmp/Blast/ref_seqs")
	parser.add_argument('-f', '--ref_fa_file', help='Reference fasta sequences', default="tmp/Sequences/ref_seq.fa")
	parser.add_argument('-ms', '--master_seq_dir', help='Master sequence directory', default="tmp/Blast/master_seq")
	parser.add_argument('-t', '--tmp_dir', help='Temp directory', default="tmp")
	parser.add_argument('-m', '--master_ref', help='Master reference accession.', required=True)
	parser.add_argument('-n', '--nextalign_dir', help='Nextalign output directory', default="Nextalign")
	parser.add_argument('-ra', '--ref_alignment_file', help='Use custom reference alignment file')
	parser.add_argument('--update_db', help='Existing update DB', default=None)
	args = parser.parse_args()

	processor = NextalignAlignment(args.gB_matrix, args.query_dir, args.ref_dir, args.ref_fa_file, args.master_seq_dir, args.tmp_dir, args.master_ref, args.nextalign_dir, args.ref_alignment_file, args.update_db)
	processor.process()