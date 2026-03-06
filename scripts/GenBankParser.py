import os
import sqlite3
import pandas as pd
from os.path import join
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from date_utils import split_date_components
from ExportRefListFromUpdateDb import load_reference_dict

# add source example NCBI or GISAID, or user define, add temp sequences 
# temp sequences which should available for temp purpose 
class GenBankParser:
	def __init__(
		self,
		input_dir,
		base_dir,
		output_dir,
		ref_list,
		exclusion_list,
		is_segmented_virus,
		test_run=False,
		test_xml_limit=10,
		require_refs=False,
		update=None,
	):
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.ref_list = ref_list
		self.exclusion_list = exclusion_list
		self.is_segmented_virus = is_segmented_virus
		self.test_run = test_run
		self.test_xml_limit = test_xml_limit
		self.require_refs = require_refs
		self.update_db_path = update
		self.update_mode = bool(update)

		if input_dir:
			self.input_dir = input_dir
		else:
			self.input_dir = join(self.base_dir, 'GenBank-XML')

		self.out_root = join(self.base_dir, self.output_dir)
		os.makedirs(self.out_root, exist_ok=True)
		self._existing_accessions = set()

	def count_ATGCN(self, sequence):
		nucl_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
		sequence = (sequence or '').upper()
		for each_nucl in sequence:
			if each_nucl in nucl_dict:
				nucl_dict[each_nucl] += 1
		return nucl_dict['A'], nucl_dict['T'], nucl_dict['G'], nucl_dict['C'], nucl_dict['N']

	def load_ref_list(self, acc_list_file):
		if self.update_mode and self.update_db_path:
			ref_dict = load_reference_dict(self.update_db_path)
			if ref_dict:
				return ref_dict
		if acc_list_file is None:
			return []
		ref_dict = {}
		try:
			with open(acc_list_file) as file:
				for line in file:
					parts = line.strip().split("\t")
					if len(parts) == 3:
						accession, accession_type, segment = parts
					elif len(parts) == 2:
						accession, accession_type = parts
					else:
						print(f"Warning: Skipping malformed line in ref list: {line.strip()}")
						continue
						
					if accession not in ref_dict:
						ref_dict[accession] = accession_type
		except FileNotFoundError:
			raise FileNotFoundError(f"Reference list not found: {acc_list_file}")
		return ref_dict

	
	def load_exclusion_list(self, acc_list_file):
		if acc_list_file is None:
			return []
		exclusion_list = set()
		try:
			with open (acc_list_file) as file:
				for line in file:
					accession = line.strip()
					if accession:
						exclusion_list.add(accession)
		except FileNotFoundError:
			print(f"Warning: File {acc_list_file} not found. Proceeding with all the available sequences")
		return list(exclusion_list)

	def load_existing_accessions_from_db(self, db_path):
		if not db_path:
			raise ValueError('Update DB path not provided')
		if not os.path.exists(db_path):
			raise FileNotFoundError(f"Update DB not found: {db_path}")

		existing = set()
		conn = sqlite3.connect(db_path)
		try:
			cursor = conn.cursor()

			cursor.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='meta_data' LIMIT 1")
			if cursor.fetchone() is None:
				raise ValueError('Update DB must contain table: meta_data')

			cursor.execute("PRAGMA table_info(meta_data)")
			meta_cols = {row[1] for row in cursor.fetchall()}
			if 'primary_accession' not in meta_cols:
				raise ValueError('Update DB table meta_data must contain column: primary_accession')

			cursor.execute("SELECT primary_accession FROM meta_data WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) != ''")
			existing.update([row[0].strip() for row in cursor.fetchall() if row[0]])

			cursor.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='excluded_accessions' LIMIT 1")
			if cursor.fetchone() is not None:
				cursor.execute("PRAGMA table_info(excluded_accessions)")
				excl_cols = {row[1] for row in cursor.fetchall()}
				if 'primary_accession' in excl_cols:
					cursor.execute("SELECT primary_accession FROM excluded_accessions WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) != ''")
					existing.update([row[0].strip() for row in cursor.fetchall() if row[0]])
		finally:
			conn.close()

		print(f"[update] Loaded {len(existing)} existing/excluded accessions from {db_path}")
		return existing

	def xml_to_tsv(self, xml_file, ref_seq_dict: dict, exclusion_acc_list: list):
		tree = ET.parse(xml_file)
		root = tree.getroot()
		data = []
		
		#ref_seq_dict = self.load_ref_list(self.ref_list)
		#exclusion_acc_list = self.load_exclusion_list(self.exclusion_list)
				
		for gbseq in root.findall('GBSeq'):
			content = {}
			content['locus'] = gbseq.find('GBSeq_locus').text
			content['length'] = gbseq.find('GBSeq_length').text
			content['data_source'] = 'ncbi'
			content['strandedness'] = gbseq.find('GBSeq_strandedness').text if gbseq.find('GBSeq_strandedness') is not None else None
			content['molecule_type'] = gbseq.find('GBSeq_moltype').text
			content['topology'] = gbseq.find('GBSeq_topology').text
			content['division'] = gbseq.find('GBSeq_division').text
			content['update_date'] = gbseq.find('GBSeq_update-date').text
			content['create_date'] = gbseq.find('GBSeq_create-date').text
			content['definition'] = gbseq.find('GBSeq_definition').text
			content['primary_accession'] = gbseq.find('GBSeq_primary-accession').text
			content['accession_version'] = gbseq.find('GBSeq_accession-version').text

			content['gi_number'] = content['primary_accession']
			content['source'] = gbseq.find('GBSeq_source').text
			content['organism'] = gbseq.find('GBSeq_organism').text
			content['taxonomy'] = gbseq.find('GBSeq_taxonomy').text
			
			if content['primary_accession'] in ref_seq_dict:
				content['accession_type'] = ref_seq_dict[content['primary_accession']]
			else:
				content['accession_type'] = 'query'

			if self.update_mode and self._existing_accessions:
				acc_type = str(content['accession_type']).strip().lower()
				if content['primary_accession'] in self._existing_accessions and acc_type not in {'master', 'reference'}:
					continue
	
			if content['primary_accession'] in exclusion_acc_list:
				content['accession_type'] = 'excluded'
				content['exclusion_criteria'] = 'excluded by the user'
				content['exclusion_status'] = '1'
			else:
				content['exclusion_criteria'] = ''
				content['exclusion_status'] = '0'

			mol_type = ''
			strain = ''
			isolate = ''
			isolation_source = ''
			db_xref = ''
			country = ''
			host = ''
			collection_date = ''
			segment = ''
			serotype = ''
			genes = []
			cds = []

			for gb_feature in gbseq.findall('GBSeq_feature-table/GBFeature'):
				if gb_feature.find('GBFeature_key').text == 'source':
					for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
						name = qualifier.find('GBQualifier_name').text if qualifier.find('GBQualifier_name') is not None else None
						value = qualifier.find('GBQualifier_value').text if qualifier.find('GBQualifier_value') is not None else None
						if name == 'mol_type':
							mol_type = value
						if name == 'strain':
							strain = value
						elif name == 'isolate':
							isolate = value
						elif name == 'isolation_source':
							isolation_source = value
						elif name == 'db_xref':
							db_xref = value
						elif name == 'country' or name == "geo_loc_name":
							country = value
						elif name == 'host':
							host = value
						elif name == 'collection_date':
							collection_date = value
						elif name == 'segment':
							segment = value
						elif name == 'serotype':
							serotype = value
				elif gb_feature.find('GBFeature_key').text == 'gene':
						gene_info = {}
						gene_info['gene_location'] = gb_feature.find('GBFeature_location').text
						for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
							name = qualifier.find('GBQualifier_name').text if qualifier.find('GBQualifier_name') is not None else None
							value = qualifier.find('GBQualifier_value').text if qualifier.find('GBQualifier_value') is not None else None
							if name == 'gene':
								gene_info['gene_name'] = value
								genes.append(gene_info)
				elif gb_feature.find('GBFeature_key').text == 'CDS':
						cds_info = {}
						cds_info['cds_location'] = gb_feature.find('GBFeature_location').text if gb_feature.find('GBFeature_location') is not None else None
						for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
							name = qualifier.find('GBQualifier_name').text if qualifier.find('GBQualifier_name') is not None else None
							value = qualifier.find('GBQualifier_value').text if qualifier.find('GBQualifier_value') is not None else None
							cds_info[name] = value
						cds.append(cds_info)

			pubmed_ids = []
			for reference in gbseq.findall('GBSeq_references/GBReference'):
				pubmed_tag = reference.find("GBReference_pubmed")
				if pubmed_tag is not None:
					pubmed_id = pubmed_tag.text
					if pubmed_id:
						pubmed_ids.append(pubmed_id)

			content['pubmed_id'] = '; '.join(pubmed_ids)
			content['mol_type'] = mol_type
			content['strain'] = strain
			content['isolate'] = isolate
			content['isolation_source'] = isolation_source
			content['db_xref'] = db_xref
			
			if ":" in country:
				tmp_country = country.split(":")
				content['country'] = tmp_country[0]
				content['geo_loc'] = tmp_country[1] if len(tmp_country) > 1 else ""
			else:
				content['country'] = country
				content['geo_loc'] = ""
		
			content['host'] = host
			content['collection_date'] = collection_date
			
			split_collection_date = split_date_components(content['collection_date'])
			content['collection_day'] = split_collection_date['day']
			content['collection_mon'] = split_collection_date['month']
			content['collection_year'] = split_collection_date['year']

			content['segment'] = segment
			content['serotype'] = serotype
			sequence = gbseq.find('GBSeq_sequence')
			content['sequence'] = sequence.text if sequence is not None else ''
			content['genes'] = '; '.join([f"{gene['gene_name']}({gene['gene_location']})" for gene in genes])
			content['cds_info'] = '; '.join([f"{k}: {v}" for each_cds in cds for k, v in each_cds.items()])
			nucl_count = self.count_ATGCN(content['sequence'])
			content['a'] = nucl_count[0]
			content['t'] = nucl_count[1]
			content['g'] = nucl_count[2]
			content['c'] = nucl_count[3]
			content['n'] = nucl_count[4]
			atgc_list = [int(content['a']), int(content['t']), int(content['g']), int(content['c'])]
			content['real_length'] = sum(atgc_list)
			content['comment'] = "NA"
			references = []
			for reference in gbseq.findall('GBSeq_references/GBReference'):
				authors = [author.text for author in reference.findall('GBReference_authors/GBAuthor')]
				reference_entry = {
					'reference_number': reference.find('GBReference_reference').text if reference.find('GBReference_reference') is not None else '',
					'position': reference.find('GBReference_position').text if reference.find('GBReference_position') is not None else '',
					'authors': ', '.join([a for a in authors if a]),
					'title': reference.find('GBReference_title').text if reference.find('GBReference_title') is not None else '',
					'journal': reference.find('GBReference_journal').text if reference.find('GBReference_journal') is not None else '',
				}
				references.append(reference_entry)

			content['reference_number'] = '; '.join([r['reference_number'] for r in references if r['reference_number']])
			content['position'] = '; '.join([r['position'] for r in references if r['position']])
			content['authors'] = '; '.join([r['authors'] for r in references if r['authors']])
			content['title'] = '; '.join([r['title'] for r in references if r['title']])
			content['journal'] = '; '.join([r['journal'] for r in references if r['journal']])
			data.append(content)
	
		return data

	def count_xml_files(self):
		directory = self.input_dir
		if not os.path.isdir(directory):
			print(f"Error: '{directory}' is not a valid directory.")
			return 0

		xml_files = [f for f in os.listdir(directory) if f.lower().endswith('.xml')]
		return len(xml_files)

	def list_xml_files(self):
		if not os.path.isdir(self.input_dir):
			print(f"Error: '{self.input_dir}' is not a valid directory.")
			return []
		return sorted([f for f in os.listdir(self.input_dir) if f.lower().endswith('.xml')])

	def find_ref_xml_files(self, ref_accessions, xml_files):
		if not ref_accessions:
			return set(), []

		ref_set = set(ref_accessions)
		found_refs = set()
		ref_xml_files = set()
		for xml_name in xml_files:
			xml_path = join(self.input_dir, xml_name)
			try:
				for _, elem in ET.iterparse(xml_path, events=("end",)):
					if elem.tag == "GBSeq":
						acc_elem = elem.find('GBSeq_primary-accession')
						if acc_elem is not None:
							acc = acc_elem.text
							if acc in ref_set:
								found_refs.add(acc)
								ref_xml_files.add(xml_name)
						elem.clear()
				if found_refs == ref_set:
					break
			except ET.ParseError as e:
				print(f"Warning: Skipping malformed XML {xml_name}: {e}")

		missing_refs = sorted(list(ref_set - found_refs))
		return ref_xml_files, missing_refs

	def select_xml_files(self, ref_seq_dict, xml_files):
		if not self.test_run:
			return xml_files

		if self.test_xml_limit <= 0:
			return []

		ref_accessions = list(ref_seq_dict.keys())
		ref_xml_files, missing_refs = self.find_ref_xml_files(ref_accessions, xml_files)
		
		# Check if any missing refs are CRITICAL (master/reference). 
		# We can tolerate missing exclusion_list refs.
		if self.require_refs and missing_refs:
			critical_missing = []
			for m in missing_refs:
				# Default to 'reference' if type not specified, so treat as critical
				rtype = ref_seq_dict.get(m, 'reference').strip().lower()
				if rtype != 'exclusion_list':
					if self.update_mode and m in self._existing_accessions:
						continue
					critical_missing.append(m)
			
			if critical_missing:
				raise ValueError(f"Missing CRITICAL reference accessions in XML input: {', '.join(critical_missing)}")
			else:
				print(f"Warning: {len(missing_refs)} exclusion_list references were not found in XML input. Proceeding as they are optional.")

		selected = list(ref_xml_files)
		remaining = [f for f in xml_files if f not in ref_xml_files]
		for f in remaining:
			if len(selected) >= self.test_xml_limit:
				break
			selected.append(f)

		return selected

	def process(self):

		ref_seq_dict = self.load_ref_list(self.ref_list)
		exclusion_acc_list = self.load_exclusion_list(self.exclusion_list)

		if self.update_mode:
			self._existing_accessions = self.load_existing_accessions_from_db(self.update_db_path)
			# Update-mode references are DB-backed source-of-truth; if missing in DB we still
			# allow parsing from XML, but emit deterministic diagnostics.
			ref_accs = sorted(ref_seq_dict.keys())
			missing_in_db = [r for r in ref_accs if r not in self._existing_accessions]
			if missing_in_db:
				print(f"[update] {len(missing_in_db)} references from ref_list are not present in DB and will be parsed from XML if available")

		xml_files = self.list_xml_files()
		if self.require_refs and self.ref_list is None and not self.update_mode:
			raise ValueError("Reference list is required when --require_refs is set")

		# Pass full dict, not just keys, to support type checking
		selected_xml_files = self.select_xml_files(ref_seq_dict, xml_files)

		if self.require_refs and not self.test_run:
			ref_accessions = list(ref_seq_dict.keys())
			_, missing_refs = self.find_ref_xml_files(ref_accessions, xml_files)
			if missing_refs:
				critical_missing = []
				for m in missing_refs:
					rtype = ref_seq_dict.get(m, 'reference').strip().lower()
					if rtype != 'exclusion_list':
						if self.update_mode and m in self._existing_accessions:
							continue
						critical_missing.append(m)
				
				if critical_missing:
					raise ValueError(f"Missing CRITICAL reference accessions in XML input: {', '.join(critical_missing)}")
				else:
					print(f"Warning: {len(missing_refs)} exclusion_list references were not found in XML input.")

		total_xml = len(selected_xml_files)
		if total_xml == 0:
			print(f"No XML files found in: {self.input_dir}")
			return

		count = 1
		merged_data = []
		for each_xml in selected_xml_files:
			print(f"Parsing: {count} of {total_xml}: " + each_xml)
			data = self.xml_to_tsv(join(self.input_dir, each_xml), ref_seq_dict, exclusion_acc_list)
			merged_data.extend(data)
			count+=1

		if not merged_data:
			if self.update_mode:
				print('[update] No new accessions found in XMLs (nothing to write).')
			else:
				print('No records parsed.')
			return

		df = pd.DataFrame(merged_data)
  		# this was subsetting on locus-??
		if 'segment' in df.columns:
			df = df.drop_duplicates(subset=['primary_accession', 'segment'], keep="last")
		else:
			df = df.drop_duplicates(subset='primary_accession', keep="last")
		with open(join(self.out_root, 'sequences.fa'), 'w') as fasta_file:
			for _, row in df.iterrows():
				fasta_file.write(f">{row['primary_accession']}\n{row['sequence']}\n")

		df.drop(columns=['sequence'], inplace=True)
		df.to_csv(join(self.out_root, 'gB_matrix_raw.tsv'), sep='\t', index=False)
		print(f"Input XML dir used: {self.input_dir}")
		print(f"Output dir used:    {self.out_root}")
		
# if running in update mode, need to check all the references are in the xml input or if they're in the DB
# if they're in the DB then pull them out of the DB and add to the output tsv
# # then pull the raw fastas out of the DB for the existing references and 

if __name__ == "__main__":
	parser = ArgumentParser(description='Extract GenBank XML files to a TSV table')
	parser.add_argument('-d', '--input_dir', help='Input directory (optional override)', default=None)
	parser.add_argument('-b', '--base_dir', help='Base directory', default='tmp')
	parser.add_argument('-o', '--output_dir', help='Output directory', default='GenBank-matrix')
	parser.add_argument('-r', '--ref_list', help='Set of reference accessions', required=True)
	parser.add_argument('-e', '--exclusion_list', help='Set of sequence accssions to be excluded')
	parser.add_argument('-s', '--is_segmented_virus', help='Is segmented virus (Y/N)', default='N')
	parser.add_argument('--test_run', action='store_true', help='Parse only a subset of XML files')
	parser.add_argument('--test_xml_limit', type=int, default=10, help='Max number of XML files to parse in test mode')
	parser.add_argument('--require_refs', action='store_true', help='Fail if any ref_list accession is missing in XML input')
	parser.add_argument('--update', default=None, help='Path to existing SQLite DB generated by CreateSqliteDB.py; when set, skip accessions already present or excluded in that DB')
	args = parser.parse_args()

	gb_parser = GenBankParser(
		args.input_dir,
		args.base_dir,
		args.output_dir,
		args.ref_list,
		args.exclusion_list,
		args.is_segmented_virus,
		test_run=args.test_run,
		test_xml_limit=args.test_xml_limit,
		require_refs=args.require_refs,
		update=args.update,
	)
	gb_parser.process()

