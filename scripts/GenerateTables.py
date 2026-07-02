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
from collections import Counter

class GenerateTables:
	def __init__(self, genbank_matrix, base_dir, output_dir, blast_hits, paded_aln, host_taxa_file, nextalign_dir, email, reference_tsv=None):
		self.genbank_matrix = genbank_matrix
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.blast_hits = blast_hits
		self.paded_aln = paded_aln
		self.host_taxa_file = host_taxa_file
		self.nextalign_dir = nextalign_dir
		self.email = email
		self.reference_tsv = reference_tsv
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

	def _print_host_value_counts(self):
		host_counts = Counter()
		host_taxa_id_counts = Counter()
		with open(self.genbank_matrix) as file:
			csv_reader = csv.DictReader(file, delimiter='\t')
			for row in csv_reader:
				host_val = (row.get('host') or '').strip() or '<blank>'
				host_taxa_id_val = (row.get('host_taxa_id') or '').strip() or '<blank>'
				host_counts[host_val] += 1
				host_taxa_id_counts[host_taxa_id_val] += 1

		print("\n[debug] value_counts for column 'host' (top 20)")
		print("count\thost")
		for value, count in host_counts.most_common(20):
			print(f"{count}\t{value}")

		print("\n[debug] value_counts for column 'host_taxa_id' (top 20)")
		print("count\thost_taxa_id")
		for value, count in host_taxa_id_counts.most_common(20):
			print(f"{count}\t{value}")

	def host_table(self):
		existing_host_taxa_id = self.host_taxa_file_check()
		taxa_file = join(self.base_dir, self.output_dir, self.host_taxa_file)
		self._print_host_value_counts()
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
			if not isinstance(host_taxa, dict):
				print(f"[warn] Taxonomy lookup for {each_host_taxaid} returned non-dict response: {host_taxa}")
				host_taxa = {
					"Taxonomy ID": each_host_taxaid,
					"Scientific Name": "N/A",
					"Rank": "N/A",
					"Lineage": "N/A",
				}
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
		if self.genbank_matrix and os.path.isfile(self.genbank_matrix):
			with open(self.genbank_matrix, newline='') as handle:
				reader = csv.DictReader(handle, delimiter='\t')
				if reader.fieldnames and 'primary_accession' in reader.fieldnames and 'segment' in reader.fieldnames:
					for row in reader:
						acc = (row.get('primary_accession') or '').strip()
						seg = (row.get('segment') or '').strip()
						if acc:
							segment_map[acc] = seg

		# Override reference segments using reference_tsv
		if self.reference_tsv and os.path.exists(self.reference_tsv):
			try:
				from ExportRefListFromUpdateDb import load_reference_file_table
				ref_df = load_reference_file_table(self.reference_tsv)
				if "primary_accession" in ref_df.columns and "segment" in ref_df.columns:
					for _, row in ref_df.iterrows():
						acc = str(row["primary_accession"]).strip()
						seg = str(row["segment"]).strip()
						digits = ''.join(ch for ch in seg if ch.isdigit())
						normalized_seg = digits if digits else seg
						if acc and normalized_seg:
							segment_map[acc] = normalized_seg
			except Exception as e:
				print(f"[warn] Failed to load reference segments from reference_tsv: {e}")

		return segment_map

	def created_alignment_table(self, blast_dict, segment_map=None):
		segment_map = segment_map or {}
		accessions = {}
		missing_accs = []
		seqs = {}
		written_accessions = set()
		inferred_segs = {}
		header = ["primary_accession", "alignment_name", "alignment", "segment"]
		write_file = open(join(self.base_dir, self.output_dir, "sequence_alignment.tsv"), 'w')
		write_file.write("\t".join(header) + "\n")
		
		# Build segment_to_master map from genbank_matrix
		segment_to_master = {}
		with open(self.genbank_matrix) as file:
			csv_reader = csv.DictReader(file, delimiter='\t')
			for row in csv_reader:
				acc_type = str(row.get('accession_type', '')).strip().lower()
				if acc_type == 'master':
					acc = str(row.get('primary_accession', '')).strip()
					seg = str(row.get('segment', '')).strip()
					digits = ''.join(ch for ch in seg if ch.isdigit())
					normalized_seg = digits if digits else seg
					if acc and normalized_seg:
						segment_to_master[normalized_seg] = acc

		# Handle single file or list of files
		paded_aln_files = self.paded_aln if isinstance(self.paded_aln, list) else [self.paded_aln]
		
		for paded_file in paded_aln_files:
			if not os.path.exists(paded_file):
				print(f"Warning: Padded alignment file not found: {paded_file}")
				continue
				
			# Infer segment from filename (e.g. refset_4_aln_merged_MSA.fasta -> "4")
			match = re.search(r"refset_(\d+)_", os.path.basename(paded_file))
			inferred_seg = match.group(1) if match else ""
			
			rds = read_file.fasta(paded_file)
			for rows in rds:
				acc = rows[0].strip()
				if acc not in accessions:
					seqs[acc] = rows[1]
					accessions[acc] = 1
					inferred_segs[acc] = inferred_seg
					seg = segment_map.get(acc, "")
					if not seg:
						seg = inferred_seg
					if acc in blast_dict:
						write_file.write(acc + '\t' + blast_dict[acc] + '\t' + rows[1] + '\t' + seg + '\n')
						written_accessions.add(acc)
					else:
						missing_accs.append(acc)

		for each_ref_aln in os.listdir(join(self.nextalign_dir)):
			for each_ref_aln_file in os.listdir(join(self.nextalign_dir, each_ref_aln)):
				rds = read_file.fasta(join(self.nextalign_dir, each_ref_aln, each_ref_aln_file, each_ref_aln_file + ".aligned.fasta"))
				for rows in rds:
					acc = rows[0].strip()
					if acc in missing_accs:
						seg = segment_map.get(acc, "")
						if not seg:
							seg = inferred_segs.get(acc, "")
						write_file.write(acc + '\t' + each_ref_aln_file + '\t' + seqs[acc] + '\t' + seg + '\n')
						written_accessions.add(acc)

		# Resolve remaining missing reference/master accessions that failed Nextalign
		still_missing = []
		for acc in missing_accs:
			if acc not in written_accessions:
				seg = segment_map.get(acc, "")
				if not seg:
					seg = inferred_segs.get(acc, "")
				digits = ''.join(ch for ch in seg if ch.isdigit())
				normalized_seg = digits if digits else seg
				master_acc = segment_to_master.get(normalized_seg)
				if master_acc:
					write_file.write(acc + '\t' + master_acc + '\t' + seqs[acc] + '\t' + seg + '\n')
					written_accessions.add(acc)
				else:
					still_missing.append(acc)
		if still_missing:
			print(f"[warn] Could not resolve master sequences for {len(still_missing)} accessions: {', '.join(still_missing[:10])}")
			
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
	parser.add_argument('-r', '--reference_tsv', help='Optional reference list TSV/txt file', default=None)
	args = parser.parse_args()
	
	processor = GenerateTables(args.genbank_matrix, args.base_dir, args.output_dir, args.blast_hits, args.paded_aln, args.host_taxa, args.nextalign_dir, args.email, args.reference_tsv)
	processor.process()
