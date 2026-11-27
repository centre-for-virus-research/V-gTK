import os
import sqlite3
import pandas as pd
from os.path import join
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from date_utils import split_date_components
from ExportRefListFromUpdateDb import load_reference_dict, load_reference_file_dict
import accession_utils
import segment_utils

# add source example NCBI or GISAID, or user define, add temp sequences 
# temp sequences which should available for temp purpose 
class GenBankParser:
	_MATRIX_COLUMNS = [
		'locus', 'length', 'data_source', 'strandedness', 'molecule_type',
		'topology', 'division', 'update_date', 'create_date', 'definition',
		'primary_accession', 'accession_version', 'gi_number', 'source',
		'organism', 'taxonomy', 'accession_type', 'exclusion_criteria',
		'exclusion_status', 'pubmed_id', 'mol_type', 'strain', 'isolate',
		'isolation_source', 'db_xref', 'country', 'geo_loc', 'host',
		'collection_date', 'collection_day', 'collection_mon', 'collection_year',
		'segment', 'serotype', 'genes', 'cds_info', 'a', 't', 'g', 'c', 'n',
		'real_length', 'comment', 'reference_number', 'position', 'authors',
		'title', 'journal',
	]

	def __init__(
		self,
		input_dir=None,
		base_dir='tmp',
		output_dir='GenBank-matrix',
		ref_list=None,
		exclusion_list=None,
		is_segmented_virus='N',
		test_run=False,
		test_xml_limit=10,
		require_refs=False,
		update=None,
		min_length_ratio=0.9,
		allow_external_refs=False,
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
		self.min_length_ratio = float(min_length_ratio) if min_length_ratio is not None else None
		self.allow_external_refs = bool(allow_external_refs)

		if input_dir:
			self.input_dir = input_dir
		else:
			self.input_dir = join(self.base_dir, 'GenBank-XML')

		self.out_root = join(self.base_dir, self.output_dir)
		os.makedirs(self.out_root, exist_ok=True)
		self._existing_accessions = set()
		self._existing_versions = {}

	def _normalize_segment_value(self, raw_segment):
		"""Canonicalise the GenBank ``/segment=`` qualifier.

		This used to return the qualifier verbatim for segmented viruses, i.e. no
		normalisation at all, which is how ``'segment 7'``, ``'RNA 4'``, ``'MA'``
		and ``'M'`` reached the definitive ``segment`` column of a real influenza
		build - 62 rows of it. Submitters write that column freehand, so it has to
		be normalised at the point of entry like any other free text.
		"""
		if self.is_segmented_virus != 'Y':
			return '1'
		normalised = segment_utils.normalise_segment(raw_segment)
		if normalised is None or normalised.casefold() in segment_utils.PANDAS_NULL_TOKENS:
			return ''
		return normalised

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
			return {}
		try:
			ref_dict = load_reference_file_dict(acc_list_file)
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
		"""Which accessions does the update DB already hold, and at which version?

		Returns the set of BARE ``meta_data.primary_accession`` values, as it
		always has, and additionally populates :attr:`_existing_versions` - a
		``{bare_accession: highest_stored_version}`` map read from
		``meta_data.accession_version``, the one column that is allowed to keep
		its ``.<digits>`` suffix.

		The version map exists because the returned set is not enough to decide
		what to skip. ``xml_to_tsv`` skipped every accession in that set, and the
		set is keyed on the bare accession, so a record GenBank had *revised* -
		``PV547761.1`` becoming ``PV547761.2`` - looked identical to one that had
		not changed at all. GenBankFetcher downloads a revision deliberately and
		records it in update_accessions.tsv; the parser then dropped it and the
		run still exited 0, so ``--update_db`` could never apply a corrected
		sequence.

		Two shapes deliberately stay out of the map, because "no stored version"
		must keep meaning "skip, exactly as before":

		* a DB with no ``accession_version`` column at all (older builds, and the
		  shared test fixture) - nothing is knowable, so nothing is re-parsed;
		* an accession whose stored version is absent or unparseable.

		Excluded accessions are removed from the map as well. A revision must not
		resurrect a record an operator excluded by hand: ``exclusion_status`` is
		cleared-on-update by CreateSqliteDB, so re-parsing one would overwrite the
		stored exclusion with ``'0'``. GenBankFetcher never re-fetches an excluded
		base either, so the two agree.
		"""
		if not db_path:
			raise ValueError('Update DB path not provided')
		if not os.path.exists(db_path):
			raise FileNotFoundError(f"Update DB not found: {db_path}")

		existing = set()
		excluded = set()
		versions = {}
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
			if 'exclusion_status' in meta_cols:
				cursor.execute(
					"SELECT primary_accession FROM meta_data WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) != '' "
					"AND LOWER(TRIM(COALESCE(CAST(exclusion_status AS TEXT), ''))) NOT IN ('', '0', 'false', 'no', 'na', 'none', 'nan')"
				)
				excluded.update([row[0].strip() for row in cursor.fetchall() if row[0]])
				existing.update(excluded)

			if 'accession_version' in meta_cols:
				cursor.execute(
					"SELECT primary_accession, accession_version FROM meta_data "
					"WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) != '' "
					"AND accession_version IS NOT NULL AND TRIM(CAST(accession_version AS TEXT)) != ''"
				)
				for stored_accession, stored_version in cursor.fetchall():
					bare = accession_utils.normalise_accession(stored_accession)
					if bare is None:
						continue
					_, version = accession_utils.split_accession_version(stored_version)
					if version is None:
						continue
					# A segmented record has one row per segment; the newest
					# version any of them carries is what GenBank last gave us.
					if version > versions.get(bare, -1):
						versions[bare] = version

			cursor.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='excluded_accessions' LIMIT 1")
			if cursor.fetchone() is not None:
				cursor.execute("PRAGMA table_info(excluded_accessions)")
				excl_cols = {row[1] for row in cursor.fetchall()}
				if 'primary_accession' in excl_cols:
					cursor.execute("SELECT primary_accession FROM excluded_accessions WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) != ''")
					excluded_rows = [row[0].strip() for row in cursor.fetchall() if row[0]]
					existing.update(excluded_rows)
					excluded.update(excluded_rows)
		finally:
			conn.close()

		for excluded_accession in excluded:
			bare = accession_utils.normalise_accession(excluded_accession)
			if bare is not None:
				versions.pop(bare, None)

		self._existing_versions = versions
		print(f"[update] Loaded {len(existing)} existing/excluded accessions from {db_path}")
		if versions:
			print(f"[update] Loaded stored accession_version for {len(versions)} of them; revised records will be re-parsed")
		return existing

	def _is_revision_of_existing(self, primary_accession, accession_version):
		"""Is this XML record NEWER than the copy the update DB already holds?

		The only reason to re-parse an accession the DB already has. Both sides go
		through :func:`accession_utils.split_accession_version` so that the bare
		identity key and the version are separated by the same rule - the old
		``rpartition('.')`` idiom would read the ``.1`` of a strain name as a
		version number.

		Returns ``False`` whenever the comparison cannot be made - the DB has no
		stored version for this accession, or the incoming record carries no
		parseable one. Unprovable is not the same as newer, and the safe answer is
		the historical one: skip.
		"""
		bare = accession_utils.normalise_accession(primary_accession)
		if bare is None:
			return False
		stored_version = self._existing_versions.get(bare)
		if stored_version is None:
			return False
		_, incoming_version = accession_utils.split_accession_version(accession_version)
		if incoming_version is None:
			return False
		return incoming_version > stored_version

	@staticmethod
	def _xml_text(parent, path, default=''):
		if parent is None:
			return default
		elem = parent.find(path)
		if elem is None or elem.text is None:
			return default
		return elem.text

	@staticmethod
	def _lookup_norm(mapping, acc):
		if not mapping or not acc:
			return None
		if acc in mapping:
			return mapping[acc]
		norm = accession_utils.normalise_accession(acc)
		if norm and norm in mapping:
			return mapping[norm]
		return None

	def xml_to_tsv(self, xml_file, ref_seq_dict: dict, exclusion_acc_list: list):
		try:
			tree = ET.parse(xml_file)
		except ET.ParseError as exc:
			raise RuntimeError(f"Malformed GenBank XML in {xml_file}: {exc}") from exc
		root = tree.getroot()
		data = []
		norm_ref_dict = dict(ref_seq_dict or {})
		for k, v in list(norm_ref_dict.items()):
			nk = accession_utils.normalise_accession(k)
			if nk and nk not in norm_ref_dict:
				norm_ref_dict[nk] = v
		norm_exclusions = set(exclusion_acc_list or [])
		for ex in list(norm_exclusions):
			nx = accession_utils.normalise_accession(ex)
			if nx:
				norm_exclusions.add(nx)
		existing_norm = set(self._existing_accessions)
		for ex in list(self._existing_accessions):
			nx = accession_utils.normalise_accession(ex)
			if nx:
				existing_norm.add(nx)

		for gbseq in root.findall('GBSeq'):
			content = {}
			content['locus'] = self._xml_text(gbseq, 'GBSeq_locus')
			content['length'] = self._xml_text(gbseq, 'GBSeq_length', '0')
			content['data_source'] = 'ncbi'
			stranded = gbseq.find('GBSeq_strandedness')
			content['strandedness'] = stranded.text if (stranded is not None and stranded.text is not None) else None
			content['molecule_type'] = self._xml_text(gbseq, 'GBSeq_moltype')
			content['topology'] = self._xml_text(gbseq, 'GBSeq_topology')
			content['division'] = self._xml_text(gbseq, 'GBSeq_division')
			content['update_date'] = self._xml_text(gbseq, 'GBSeq_update-date')
			content['create_date'] = self._xml_text(gbseq, 'GBSeq_create-date')
			content['definition'] = self._xml_text(gbseq, 'GBSeq_definition')
			content['primary_accession'] = self._xml_text(gbseq, 'GBSeq_primary-accession')
			content['accession_version'] = self._xml_text(gbseq, 'GBSeq_accession-version')

			content['gi_number'] = content['primary_accession']
			content['source'] = self._xml_text(gbseq, 'GBSeq_source')
			content['organism'] = self._xml_text(gbseq, 'GBSeq_organism')
			content['taxonomy'] = self._xml_text(gbseq, 'GBSeq_taxonomy')

			ref_type = self._lookup_norm(norm_ref_dict, content['primary_accession'])
			content['accession_type'] = ref_type if ref_type else 'query'

			# Skip what the update DB already has - UNLESS GenBank has revised it.
			if self.update_mode and existing_norm:
				acc_type = str(content['accession_type']).strip().lower()
				norm_pa = accession_utils.normalise_accession(content['primary_accession']) or content['primary_accession']
				in_existing = content['primary_accession'] in existing_norm or norm_pa in existing_norm
				if in_existing and acc_type not in {'master', 'reference'}:
					if not self._is_revision_of_existing(content['primary_accession'], content['accession_version']):
						continue
					print(
						f"[update] Re-parsing revised record {content['accession_version']} "
						f"(stored version {self._existing_versions.get(norm_pa)})"
					)

			norm_pa = accession_utils.normalise_accession(content['primary_accession']) or content['primary_accession']
			if content['primary_accession'] in norm_exclusions or norm_pa in norm_exclusions:
				content['accession_type'] = 'excluded'
				content['exclusion_criteria'] = 'excluded by the user'
				content['exclusion_status'] = '1'
			elif str(content['accession_type']).strip().lower() == 'exclusion_list':
				content['exclusion_criteria'] = 'excluded by reference list exclusion flag'
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
				fkey = self._xml_text(gb_feature, 'GBFeature_key')
				if fkey == 'source':
					for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
						name = self._xml_text(qualifier, 'GBQualifier_name', None)
						value = self._xml_text(qualifier, 'GBQualifier_value', '')
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
				elif fkey == 'gene':
						gene_info = {}
						gene_info['gene_location'] = self._xml_text(gb_feature, 'GBFeature_location')
						for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
							name = self._xml_text(qualifier, 'GBQualifier_name', None)
							value = self._xml_text(qualifier, 'GBQualifier_value', '')
							if name == 'gene':
								gene_info['gene_name'] = value
								genes.append(gene_info)
				elif fkey == 'CDS':
						cds_info = {}
						cds_info['cds_location'] = self._xml_text(gb_feature, 'GBFeature_location', None)
						for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
							name = self._xml_text(qualifier, 'GBQualifier_name', None)
							value = self._xml_text(qualifier, 'GBQualifier_value', '')
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

			country = country or ''
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

			content['segment'] = self._normalize_segment_value(segment)
			content['serotype'] = serotype
			content['sequence'] = self._xml_text(gbseq, 'GBSeq_sequence', '')
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
					'reference_number': self._xml_text(reference, 'GBReference_reference'),
					'position': self._xml_text(reference, 'GBReference_position'),
					'authors': ', '.join([a for a in authors if a]),
					'title': self._xml_text(reference, 'GBReference_title'),
					'journal': self._xml_text(reference, 'GBReference_journal'),
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
		# We can tolerate missing exclusion_list refs, and when allow_external_refs
		# is enabled, non-INSDC IDs (e.g. GISAID EPI_ISL_*) are deferred to MERGE_GISAID.
		if self.require_refs and missing_refs:
			critical_missing = []
			deferred_external = []
			for m in missing_refs:
				# Default to 'reference' if type not specified, so treat as critical
				rtype = ref_seq_dict.get(m, 'reference').strip().lower()
				if rtype != 'exclusion_list':
					if self.update_mode and m in self._existing_accessions:
						continue
					if self.allow_external_refs and accession_utils.is_gisaid_accession(m):
						deferred_external.append(m)
						continue
					critical_missing.append(m)
			
			if deferred_external:
				print(f"[info] Deferring {len(deferred_external)} non-INSDC reference(s) to external/GISAID merge (e.g. {', '.join(deferred_external[:5])}).")
			if critical_missing:
				raise ValueError(f"Missing CRITICAL reference accessions in XML input: {', '.join(critical_missing)}")
			elif not deferred_external:
				print(f"Warning: {len(missing_refs)} exclusion_list references were not found in XML input. Proceeding as they are optional.")

		selected = list(ref_xml_files)
		remaining = [f for f in xml_files if f not in ref_xml_files]
		for f in remaining:
			if len(selected) >= self.test_xml_limit:
				break
			selected.append(f)

		return selected

	def _load_db_reference_lengths(self):
		"""Load minimum reference real_length values from update_db (global and per-segment)."""
		if not self.update_db_path or not os.path.isfile(self.update_db_path):
			return {}
		seg_lengths = {}
		conn = sqlite3.connect(self.update_db_path)
		try:
			cur = conn.cursor()
			meta_cols = {r[1] for r in cur.execute("PRAGMA table_info(meta_data)")}
			if "primary_accession" not in meta_cols or "accession_type" not in meta_cols:
				return {}
			seg_expr = "segment" if "segment" in meta_cols else "'1'"
			if "real_length" in meta_cols:
				rows = cur.execute(
					f"SELECT {seg_expr}, real_length FROM meta_data "
					"WHERE LOWER(TRIM(COALESCE(accession_type, ''))) IN ('master', 'reference') "
					"AND real_length IS NOT NULL"
				).fetchall()
				for seg_val, rl_val in rows:
					try:
						rl_int = int(float(rl_val))
					except (TypeError, ValueError):
						continue
					if rl_int > 0:
						seg_key = self._normalize_segment_value(seg_val)
						seg_lengths.setdefault(seg_key, []).append(rl_int)
		except sqlite3.Error:
			return {}
		finally:
			conn.close()
		return {k: min(v) for k, v in seg_lengths.items() if v}

	def _write_empty_outputs(self):
		with open(join(self.out_root, 'sequences.fa'), 'w', encoding='utf-8') as fasta_file:
			pass
		pd.DataFrame(columns=self._MATRIX_COLUMNS).to_csv(join(self.out_root, 'gB_matrix_raw.tsv'), sep='\t', index=False)

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
				deferred_external = []
				for m in missing_refs:
					rtype = ref_seq_dict.get(m, 'reference').strip().lower()
					if rtype != 'exclusion_list':
						if self.update_mode and m in self._existing_accessions:
							continue
						if self.allow_external_refs and not accession_utils.is_accession(m):
							deferred_external.append(m)
							continue
						critical_missing.append(m)
				
				if deferred_external:
					print(f"[info] Deferring {len(deferred_external)} non-INSDC reference(s) to external/GISAID merge (e.g. {', '.join(deferred_external[:5])}).")
				if critical_missing:
					raise ValueError(f"Missing CRITICAL reference accessions in XML input: {', '.join(critical_missing)}")
				elif not deferred_external:
					print(f"Warning: {len(missing_refs)} exclusion_list references were not found in XML input.")

		total_xml = len(selected_xml_files)
		if total_xml == 0:
			print(f"No XML files found in: {self.input_dir}")
			if self.update_mode:
				self._write_empty_outputs()
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
				print('[update] No new accessions found in XMLs (writing empty headered outputs).')
				self._write_empty_outputs()
			else:
				print('No records parsed.')
			return

		# --- length filter ---
		if self.min_length_ratio and self.min_length_ratio > 0:
			if self.is_segmented_virus == 'Y':
				ref_len_by_seg = {}
				for rec in merged_data:
					if str(rec.get('accession_type', '')).lower() in ('master', 'reference'):
						rl = int(rec.get('real_length', 0))
						seg_key = rec.get('segment', '') or ''
						if rl > 0:
							ref_len_by_seg.setdefault(seg_key, []).append(rl)
				seg_min_ref = {k: min(v) for k, v in ref_len_by_seg.items() if v}
				if self.update_mode:
					for seg_k, db_min in self._load_db_reference_lengths().items():
						if seg_k not in seg_min_ref:
							seg_min_ref[seg_k] = db_min
				if seg_min_ref:
					for rec in merged_data:
						acc_type = str(rec.get('accession_type', '')).lower()
						if acc_type in ('master', 'reference', 'excluded', 'exclusion_list'):
							continue
						seg_key = rec.get('segment', '') or ''
						min_ref_length = seg_min_ref.get(seg_key)
						if not min_ref_length:
							continue
						length_cutoff = min_ref_length * self.min_length_ratio
						rl = int(rec.get('real_length', 0))
						if rl < length_cutoff:
							print(f"[length_filter] Excluding {rec.get('primary_accession', '?')} (seg={seg_key}): real_length={rl} < {length_cutoff:.0f}")
							rec['exclusion_status'] = '1'
							rec['exclusion_criteria'] = (
								f'sequence_too_short: real_length={rl} < {length_cutoff:.0f} '
								f'(min_ref={min_ref_length} * ratio={self.min_length_ratio})'
							)
				else:
					print('[length_filter] No reference sequences found; skipping length filter.')
			else:
				ref_lengths = [
					int(rec.get('real_length', 0))
					for rec in merged_data
					if str(rec.get('accession_type', '')).lower() in ('master', 'reference')
					and int(rec.get('real_length', 0)) > 0
				]
				if not ref_lengths and self.update_mode:
					db_mins = list(self._load_db_reference_lengths().values())
					ref_lengths.extend(db_mins)
				if ref_lengths:
					min_ref_length = min(ref_lengths)
					length_cutoff = min_ref_length * self.min_length_ratio
					print(f'[length_filter] min ref real_length={min_ref_length}, cutoff={length_cutoff:.0f} (ratio={self.min_length_ratio})')
					for rec in merged_data:
						acc_type = str(rec.get('accession_type', '')).lower()
						if acc_type in ('master', 'reference', 'excluded', 'exclusion_list'):
							continue
						rl = int(rec.get('real_length', 0))
						if rl < length_cutoff:
							print(f"[length_filter] Excluding {rec.get('primary_accession', '?')}: real_length={rl} < {length_cutoff:.0f}")
							rec['exclusion_status'] = '1'
							rec['exclusion_criteria'] = (
								f'sequence_too_short: real_length={rl} < {length_cutoff:.0f} '
								f'(min_ref={min_ref_length} * ratio={self.min_length_ratio})'
							)
				else:
					print('[length_filter] No reference sequences found; skipping length filter.')

		df = pd.DataFrame(merged_data)
  		# this was subsetting on locus-??
		if 'segment' in df.columns:
			df = df.drop_duplicates(subset=['primary_accession', 'segment'], keep="last")
		else:
			df = df.drop_duplicates(subset='primary_accession', keep="last")
		with open(join(self.out_root, 'sequences.fa'), 'w') as fasta_file:
			for _, row in df.iterrows():
				seq_val = row.get('sequence') or ''
				fasta_file.write(f">{row['primary_accession']}\n{seq_val}\n")

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
	parser.add_argument('--min_length_ratio', type=float, default=0.9, help='Minimum sequence length as a fraction of the shortest reference real_length; sequences below this are marked excluded (default: 0.9)')
	parser.add_argument('--allow_external_refs', action='store_true', help='Allow non-INSDC reference IDs (e.g. GISAID EPI_ISL_*) to be satisfied by a downstream external merge instead of failing --require_refs')
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
		min_length_ratio=args.min_length_ratio,
		allow_external_refs=args.allow_external_refs,
	)
	gb_parser.process()

