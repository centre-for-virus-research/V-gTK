import os
import re
import csv
import time
import shutil
import read_file
import subprocess
import urllib.error
from Bio import Entrez
from Bio.Seq import Seq
from os.path import join
from itertools import islice
from os.path import join, normpath
from argparse import ArgumentParser
from collections import defaultdict

'''
Normal mode:
	python scripts/GenerateTables.py
	
Update mode:
	python scripts/GenerateTables.py --genbank_matrix tmp/Update/GenBank-matrix/gB_matrix_raw.tsv --blast_hits tmp/Update/Blast/query_uniq_tophits.tsv --nextalign_dir tmp/Update/Nextalign/ --paded_aln tmp/Update/Pad-alignment/alUnc509RefseqsMafftHandModified.fa --update
'''

class GenerateTables:
	def __init__(self, genbank_matrix, base_dir, output_dir, blast_hits, paded_aln, nextalign_dir, email):
		self.genbank_matrix = genbank_matrix
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.blast_hits = blast_hits
		self.paded_aln = paded_aln
		self.nextalign_dir = nextalign_dir
		self.email = email
		os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)

	def fetch_taxonomy_details(self, tax_id, max_retries=5, delay=2):
		Entrez.email = self.email
		for attempt in range(1, max_retries + 1):
			try:
				handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
				records = Entrez.read(handle)
				time.sleep(1)
				handle.close()

				if records:
					tax_record = records[0]
					taxonomy_info = {
						"Scientific Name": tax_record.get("ScientificName", "N/A"),
						"Taxonomy ID": tax_record.get("TaxId", "N/A"),
						"Rank": tax_record.get("Rank", "N/A"),
						"Lineage": tax_record.get("Lineage", "N/A"),
						"Other Names": tax_record.get("OtherNames", {}).get("Synonym", []),
					}
					return taxonomy_info
				else:
					return "No taxonomy details found for the given ID."

			except urllib.error.HTTPError as e:
				print(f"HTTPError on attempt {attempt} for TaxID {tax_id}: {e}")
				if attempt == max_retries:
					print("Max retries reached. Skipping this TaxID.")
					return {
						"Scientific Name": "N/A",
						"Taxonomy ID": tax_id,
						"Rank": "N/A",
						"Lineage": "N/A",
						"Other Names": [],
						}
				else:
					time.sleep(delay)

	def build_insertion_map(self):
		"""
		Returns dict: {accession: insertion_string}
		If an accession appears multiple times, insertions are joined with '; ' (unique, stable order).
		"""
		ins_map = defaultdict(list)

		# nextalign_dir layout: <nextalign_dir>/<query_aln|reference_aln>/<accession>/<accession>.insertions.csv
		for aln_dir in os.listdir(self.nextalign_dir):
			aln_root = join(self.nextalign_dir, aln_dir)
			if not os.path.isdir(aln_root):
				continue

			for acc in os.listdir(aln_root):
				acc_dir = join(aln_root, acc)
				if not os.path.isdir(acc_dir):
					continue

				ins_csv = join(acc_dir, f"{acc}.insertions.csv")
				if not os.path.exists(ins_csv):
					continue

				with open(ins_csv) as f:
					# skip header
					for line in islice(f, 1, None):
						parts = line.strip().split(",")
						if len(parts) < 2:
							continue
						accession = parts[0].strip()
						insertion = parts[1].strip()
						if insertion:
							ins_map[accession].append(insertion)

		# de-dup + join
		out = {}
		for acc, vals in ins_map.items():
			seen = set()
			uniq = []
			for v in vals:
				if v not in seen:
					seen.add(v)
					uniq.append(v)
			out[acc] = "; ".join(uniq)
		return out

	def host_taxa_file_check(self):
		host_tax_id_list = []
		host_file = join(self.base_dir, self.output_dir, self.host_taxa_file)
    
		if os.path.exists(host_file):
			with open(host_file, newline='') as f:
				reader = csv.DictReader(f, delimiter='\t', fieldnames=["taxonomy_id"])
				for row in reader:
					host_tax_id_list.append(row["taxonomy_id"])

		return host_tax_id_list

	def host_table(self):
		existing_host_taxa_id = self.host_taxa_file_check()
		taxa_file = join(self.base_dir, self.output_dir, self.host_taxa_file)
		file_exists = os.path.isfile(taxa_file)
		write_file = open(taxa_file, 'a')
		if not file_exists or os.path.getsize(taxa_file) == 0:
				header = ['taxonomy_id', 'scientific_name', 'rank', 'lineage']
				write_file.write('\t'.join(header) + '\n')

		host_taxa_list = []

		with open(self.genbank_matrix) as file:
			csv_reader = csv.DictReader(file, delimiter='\t')
			for row in csv_reader:
				if row['host_taxa_id'] not in host_taxa_list:
					if 'NA' not in row['host_taxa_id']:
						if row['host_taxa_id'] not in existing_host_taxa_id:
							host_taxa_list.append(row['host_taxa_id'])

		host_taxa_list = [item for item in host_taxa_list if item]
		host_taxa_list = [x for x in host_taxa_list if x != 'nan']
		print(f"\nTotal {len(host_taxa_list)} taxa id informations to fetch\n")	
		for idx, each_host_taxaid in enumerate(host_taxa_list, start=1):
			print(f" - Fetching taxa details for {each_host_taxaid} which is {idx} of {len(host_taxa_list)}")
			host_taxa = self.fetch_taxonomy_details(each_host_taxaid)
			tax_id = host_taxa['Taxonomy ID'] if host_taxa['Taxonomy ID'] is not None else 'NA'
			scientific_name = host_taxa['Scientific Name'] if host_taxa['Scientific Name'] is not None else 'NA'
			rank = host_taxa['Rank']	if host_taxa['Rank'] is not None else 'NA'
			lineage = host_taxa['Lineage'] if host_taxa['Lineage'] is not None else 'NA'
			write_file.write(tax_id + '\t' + scientific_name + '\t' + rank + '\t' + lineage + '\n')
		write_file.close() 

	@staticmethod
	def load_blast_hits(blast_hit_file):
		acc_dict = {}
		for i in open(blast_hit_file):
			query, ref, score, strand = i.strip().split('\t')
			acc_dict[query] = ref
		return acc_dict

	''''
	def created_alignment_table(self, blast_dict):
		accessions = {}
		missing_accs = []
		seqs = {}
		header = ["primary_accession", "alignment_name", "alignment"]
		write_file = open(join(self.base_dir, self.output_dir, "sequence_alignment.tsv"), 'w')
		write_file.write("\t".join(header) + "\n")
		
		# Handle single file or list of files
		paded_aln_files = self.paded_aln if isinstance(self.paded_aln, list) else [self.paded_aln]
		
		for paded_file in paded_aln_files:
			if not os.path.exists(paded_file):
				print(f"Warning: Padded alignment file not found: {paded_file}")
				continue
			
			rds = read_file.fasta(paded_file)
			for rows in rds:
				
				if rows[0] not in accessions:
					seqs[rows[0]] = rows[1]
					accessions[rows[0]] = 1
					if rows[0] in blast_dict:
						write_file.write(rows[0].strip() + '\t' + blast_dict[rows[0].strip()] + '\t' + rows[1] + '\n')
					else:
						missing_accs.append(rows[0].strip())
		
		for each_ref_aln in os.listdir(join(self.nextalign_dir)):
			for each_ref_aln_file in os.listdir(join(self.nextalign_dir, each_ref_aln)):
				rds = read_file.fasta(join(self.nextalign_dir, each_ref_aln, each_ref_aln_file, each_ref_aln_file + ".aligned.fasta"))
				for rows in rds:
					if rows[0].strip() in missing_accs:
						write_file.write(rows[0].strip() + '\t' + each_ref_aln_file + '\t' + seqs[rows[0]] + '\n')
						accessions[rows[0].strip()] = 1

			
		write_file.close()
		print("Removing the Sequence redundancy")
		self.remove_redundancy_from_alignment(join(self.base_dir, self.output_dir, "sequence_alignment.tsv"))
	'''
	def created_alignment_table(self, blast_dict):
		insertion_map = self.build_insertion_map()  # NEW

		accessions = {}
		missing_accs = []
		seqs = {}

		header = ["primary_accession", "alignment_name", "alignment", "insertion"]  # NEW COLUMN
		write_file = open(join(self.base_dir, self.output_dir, "sequence_alignment.tsv"), 'w')
		write_file.write("\t".join(header) + "\n")

		paded_aln_files = self.paded_aln if isinstance(self.paded_aln, list) else [self.paded_aln]

		for paded_file in paded_aln_files:
			if not os.path.exists(paded_file):
				print(f"Warning: Padded alignment file not found: {paded_file}")
				continue

			rds = read_file.fasta(paded_file)
			for rows in rds:
				if rows[0] not in accessions:
					acc = rows[0].strip()
					seqs[acc] = rows[1]
					accessions[acc] = 1

					ins_val = insertion_map.get(acc, "")  # NEW
					
					if acc in blast_dict:
						write_file.write(acc + '\t' + blast_dict[acc] + '\t' + rows[1] + '\t' + ins_val + '\n')
					else:
						missing_accs.append(acc)

		for each_ref_aln in os.listdir(join(self.nextalign_dir)):
			for each_ref_aln_file in os.listdir(join(self.nextalign_dir, each_ref_aln)):
				rds = read_file.fasta(join(self.nextalign_dir, each_ref_aln, each_ref_aln_file, each_ref_aln_file + ".aligned.fasta"))
				for rows in rds:
					acc = rows[0].strip()
					if acc in missing_accs:
						ins_val = insertion_map.get(acc, "")  # NEW
						write_file.write(acc + '\t' + each_ref_aln_file + '\t' + seqs[acc] + '\t' + ins_val + '\n')
						accessions[acc] = 1

		write_file.close()
		print("Removing the Sequence redundancy")
		self.remove_redundancy_from_alignment(join(self.base_dir, self.output_dir, "sequence_alignment.tsv"))
	
	def remove_redundancy_from_alignment(self, file_path, delimiter="\t"):
		seen = set()
		lines_to_write = []

		with open(file_path, "r") as infile:
			for line in infile:
				if not line.strip():
					continue  # skip empty lines
				key = line.strip().split(delimiter)[0]
				if key not in seen:
					seen.add(key)
					lines_to_write.append(line)

		with open(file_path, "w") as outfile:
			outfile.writelines(lines_to_write)
		print(f"Redundancy removed. File updated: {file_path}")

	def process(self):
		
		blast_dictionary = self.load_blast_hits(self.blast_hits)
		#self.load_gb_matrix()
		self.created_alignment_table(blast_dictionary)
		#self.host_table()

if __name__ == "__main__":
	parser = ArgumentParser(description='Generating tables sqlite DB')
	parser.add_argument('-g', '--genbank_matrix', help='Genbank matrix table', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-b', '--base_dir', help='base directory', default="tmp")
	parser.add_argument('-o', '--output_dir', help='output directory to store all the db-ready tsv files', default="Tables")
	#parser.add_argument('-f', '--host_taxa', help='Host taxa file name, default="host_taxa.tsv"', default="host_taxa.tsv")
	parser.add_argument('-bh', '--blast_hits', help='BLASTN unique hits', default="tmp/Blast/query_uniq_tophits.tsv")
	parser.add_argument('-p', '--paded_aln', help='Paded alignment file', nargs='+', default=["tmp/Pad-alignment/NC_001542.aligned_merged_MSA.fasta"])
	parser.add_argument('-n', '--nextalign_dir', help='Nextalign aligned directory', default="tmp/Nextalign/")
	parser.add_argument('-e', '--email', help='Email id', default='your-email@example.com')

	args = parser.parse_args()
	
	processor = GenerateTables(args.genbank_matrix, args.base_dir, args.output_dir, args.blast_hits, args.paded_aln, args.nextalign_dir, args.email)
	processor.process()