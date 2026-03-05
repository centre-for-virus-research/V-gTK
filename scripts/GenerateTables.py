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
from argparse import ArgumentParser
from collections import defaultdict

class GenerateTables:
	def __init__(self, genbank_matrix, base_dir, output_dir, blast_hits, paded_aln, host_taxa_file, nextalign_dir, email):
		self.genbank_matrix = genbank_matrix
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.blast_hits = blast_hits
		self.paded_aln = paded_aln
		self.host_taxa_file = host_taxa_file
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

	def load_segment_map(self):
		segment_map = {}
		if not self.genbank_matrix or not os.path.isfile(self.genbank_matrix):
			return segment_map
		with open(self.genbank_matrix, newline='') as handle:
			reader = csv.DictReader(handle, delimiter='\t')
			if not reader.fieldnames or 'primary_accession' not in reader.fieldnames or 'segment' not in reader.fieldnames:
				return segment_map
			for row in reader:
				acc = (row.get('primary_accession') or '').strip()
				seg = (row.get('segment') or '').strip()
				if acc:
					segment_map[acc] = seg
		return segment_map

	def created_alignment_table(self, blast_dict, segment_map=None):
		segment_map = segment_map or {}
		accessions = {}
		missing_accs = []
		seqs = {}
		header = ["primary_accession", "alignment_name", "alignment", "segment"]
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
					seg = segment_map.get(rows[0].strip(), "")
					if rows[0] in blast_dict:
						write_file.write(rows[0].strip() + '\t' + blast_dict[rows[0].strip()] + '\t' + rows[1] + '\t' + seg + '\n')
					else:
						missing_accs.append(rows[0].strip())

		for each_ref_aln in os.listdir(join(self.nextalign_dir)):
			for each_ref_aln_file in os.listdir(join(self.nextalign_dir, each_ref_aln)):
				rds = read_file.fasta(join(self.nextalign_dir, each_ref_aln, each_ref_aln_file, each_ref_aln_file + ".aligned.fasta"))
				for rows in rds:
					if rows[0].strip() in missing_accs:
						seg = segment_map.get(rows[0].strip(), "")
						write_file.write(rows[0].strip() + '\t' + each_ref_aln_file + '\t' + seqs[rows[0]] + '\t' + seg + '\n')
						accessions[rows[0].strip()] = 1

			
		write_file.close()
		print("Removing the Sequence redundancy")
		self.remove_redundancy_from_alignment(join(self.base_dir, self.output_dir, "sequence_alignment.tsv"))

	def remove_redundancy_from_alignment(self, file_path, delimiter="\t"):
		seen = set()
		lines_to_write = []

		with open(file_path, "r") as infile:
			header = infile.readline()
			if header:
				lines_to_write.append(header)
				cols = header.strip().split(delimiter)
				segment_idx = cols.index("segment") if "segment" in cols else None
			else:
				segment_idx = None
			for line in infile:
				if not line.strip():
					continue  # skip empty lines
				parts = line.strip().split(delimiter)
				acc = parts[0] if len(parts) > 0 else ""
				seg = parts[segment_idx] if segment_idx is not None and segment_idx < len(parts) else ""
				key = (acc, seg)
				if key not in seen:
					seen.add(key)
					lines_to_write.append(line)

		with open(file_path, "w") as outfile:
			outfile.writelines(lines_to_write)
		print(f"Redundancy removed. File updated: {file_path}")

	def create_insertion_table(self, segment_map=None):
		segment_map = segment_map or {}
		write_file = open(join(self.base_dir, self.output_dir, "insertions.tsv"), 'w')
		header = ["primary_accession", "reference", "insertion", "segment"]
		write_file.write("\t".join(header) + "\n")
		for aln_dir in os.listdir(self.nextalign_dir):
			for each_aln_dir in os.listdir(join(self.nextalign_dir, aln_dir)):
				with open(join(self.nextalign_dir, aln_dir, each_aln_dir, each_aln_dir + ".insertions.csv")) as f:
					for each_line in islice(f, 1, None):
						accession, insertion, aa_insertion = each_line.strip().split(",")
						if len(insertion) > 0:
							data = [accession, each_aln_dir, insertion, segment_map.get(accession, "")]
							write_file.write("\t".join(data) + "\n")

		write_file.close()
					
	def process(self):
		
		blast_dictionary = self.load_blast_hits(self.blast_hits)
		segment_map = self.load_segment_map()
		#self.load_gb_matrix()
		self.created_alignment_table(blast_dictionary, segment_map)
		self.host_table()
		self.create_insertion_table(segment_map)

if __name__ == "__main__":
	parser = ArgumentParser(description='Generating tables sqlite DB')
	parser.add_argument('-g', '--genbank_matrix', help='Genbank matrix table', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-b', '--base_dir', help='base directory', default="tmp")
	parser.add_argument('-o', '--output_dir', help='output directory to store all the db-ready tsv files', default="Tables")
	parser.add_argument('-f', '--host_taxa', help='Host taxa file name, default="host_taxa.tsv"', default="host_taxa.tsv")
	parser.add_argument('-bh', '--blast_hits', help='BLASTN unique hits', default="tmp/Blast/query_uniq_tophits.tsv")
	parser.add_argument('-p', '--paded_aln', help='Paded alignment file', nargs='+', default=["tmp/Pad-alignment/NC_001542.aligned_merged_MSA.fasta"])
	parser.add_argument('-n', '--nextalign_dir', help='Nextalign aligned directory', default="tmp/Nextalign/")
	parser.add_argument('-e', '--email', help='Email id', default='your-email@example.com')
	args = parser.parse_args()
	
	processor = GenerateTables(args.genbank_matrix, args.base_dir, args.output_dir, args.blast_hits, args.paded_aln, args.host_taxa, args.nextalign_dir, args.email)
	processor.process()
