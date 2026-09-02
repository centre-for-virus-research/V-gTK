import os
import csv
import sqlite3
import time
import random
from time import sleep
from os.path import join
from argparse import ArgumentParser
from ExportRefListFromUpdateDb import load_reference_file_accessions
import accession_utils

try:
	import requests
except ModuleNotFoundError:
	class _RequestsFallback:
		class exceptions:
			class JSONDecodeError(Exception):
				pass

			class ChunkedEncodingError(Exception):
				pass

			class ConnectionError(Exception):
				pass

			class HTTPError(Exception):
				def __init__(self, *args, response=None, **kwargs):
					super().__init__(*args)
					self.response = response

		@staticmethod
		def get(*args, **kwargs):
			raise ModuleNotFoundError("requests is required for network operations")

	requests = _RequestsFallback()

class GenBankFetcher:
	def __init__(self, taxid, base_url, email, output_dir, batch_size, sleep_time, base_dir, update_file, test_run=False, ref_list=None, test_fetch_count=100):
		self.taxid = taxid
		self.base_url = base_url
		self.email = email
		self.output_dir = output_dir
		self.batch_size = batch_size
		self.efetch_batch_size = 100
		self.sleep_time = sleep_time
		self.base_dir = base_dir
		self.update_file = update_file
		self.test_run = test_run
		self.ref_list = ref_list
		self.test_fetch_count = test_fetch_count
		self.update_log = "update_accessions.tsv"
		
		# Load NCBI API key if available and not in pytest
		import sys
		if "pytest" not in sys.modules:
			self.api_key = self._load_ncbi_key()
		else:
			self.api_key = None
		
		# Ensure output directories exist immediately
		os.makedirs(self.output_dir, exist_ok=True)
		os.makedirs(join(self.output_dir, self.base_dir), exist_ok=True)

	def _load_ncbi_key(self):
		key_path = "/home3/oml4h/RABV-gTK/dev/ncbi_key"
		if os.path.exists(key_path):
			try:
				with open(key_path, "r", encoding="utf-8") as f:
					key = f.read().strip()
				if key and len(key) == 36 and all(c.isalnum() or c in "-_" for c in key):
					return key
			except Exception as e:
				print(f"Warning: Could not read NCBI key from {key_path}: {e}")
		return None

	def _get(self, url):
		if self.api_key:
			if "?" in url:
				url += f"&api_key={self.api_key}"
			else:
				url += f"?api_key={self.api_key}"
		return requests.get(url)

	def get_record_count(self):
		search_url = f"{self.base_url}esearch.fcgi?db=nucleotide&term=txid{self.taxid}[Organism:exp]&retmode=json&email={self.email}"
		response = self._get(search_url)
		response.raise_for_status()
		data = response.json()
		return int(data['esearchresult']['count'])

	def fetch_ids(self):
		start_all = time.time()

		# 1) Initialize search history and retrieve WebEnv, QueryKey and count
		t0 = time.time()
		hist_url = (
			f"{self.base_url}esearch.fcgi?db=nucleotide"
			f"&term=txid{self.taxid}[Organism:exp]"
			f"&retmax=0&idtype=acc"
			f"&usehistory=y&email={self.email}&retmode=json"
		)
		resp = self._get(hist_url)
		resp.raise_for_status()
		hist = resp.json()["esearchresult"]
		
		count = int(hist["count"])
		print(f"[fetch_ids] ESearch history → count={count:,} took {time.time() - t0:.1f}s")
		
		if count == 0:
			return []

		webenv   = hist["webenv"]
		querykey = hist["querykey"]

		# 2) Paginate in chunks of self.batch_size
		all_accs = []
		for start in range(0, count, self.batch_size):
			t1 = time.time()
			page_url = (
				f"{self.base_url}esearch.fcgi?db=nucleotide"
				f"&WebEnv={webenv}&query_key={querykey}"
				f"&retstart={start}&retmax={self.batch_size}"
				f"&idtype=acc&retmode=json"
			)
			
			max_retries = 5
			for attempt in range(max_retries):
				try:
					page_resp = self._get(page_url)
					page_resp.raise_for_status()
					data = page_resp.json()
					if "esearchresult" not in data:
						raise ValueError("Incomplete JSON response: missing 'esearchresult'")
					page = data["esearchresult"]["idlist"]
					break
				except (requests.exceptions.JSONDecodeError, requests.exceptions.ChunkedEncodingError, requests.exceptions.ConnectionError, requests.exceptions.HTTPError, ValueError) as e:
					if attempt == max_retries - 1:
						print(f"Failed to fetch IDs for chunk starting at {start} after {max_retries} attempts.")
						raise e
					# Check if it's a 429 error and wait longer if so
					is_429 = isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 429
					wait_time = (self.sleep_time * (attempt + 1)) + (10 if is_429 else 0)
					
					print(f"Error fetching IDs ({type(e).__name__}) for chunk starting at {start}. Retrying in {wait_time}s... (Attempt {attempt+1}/{max_retries})")
					sleep(wait_time)

			page = [acc.split('.', 1)[0] for acc in page]
			all_accs.extend(page)
			print(f"[fetch_ids] chunk {start:,}-{min(start+self.batch_size, count):,} "
				f"({len(page):,} records) took {time.time() - t1:.1f}s")

		print(f"[fetch_ids] total accs fetched: {len(all_accs):,} in {time.time() - start_all:.1f}s")
		return all_accs

	def fetch_accs(self):
		start_all = time.time()

		t0 = time.time()
		hist_url = (
			f"{self.base_url}esearch.fcgi?db=nucleotide"
			f"&term=txid{self.taxid}[Organism:exp]"
			f"&retmax=0&idtype=acc"
			f"&usehistory=y&email={self.email}&retmode=json"
		)
		resp = self._get(hist_url)
		resp.raise_for_status()
		hist = resp.json()["esearchresult"]

		count = int(hist["count"])
		print(f"[fetch_accs] ESearch history → count={count:,} took {time.time() - t0:.1f}s")

		if count == 0:
			return []

		webenv = hist["webenv"]
		querykey = hist["querykey"]

		all_accs = []
		for start in range(0, count, self.batch_size):
			t1 = time.time()
			page_url = (
				f"{self.base_url}esearch.fcgi?db=nucleotide"
				f"&WebEnv={webenv}&query_key={querykey}"
				f"&retstart={start}&retmax={self.batch_size}"
				f"&idtype=acc&retmode=json"
			)

			max_retries = 5
			for attempt in range(max_retries):
				try:
					page_resp = self._get(page_url)
					page_resp.raise_for_status()
					data = page_resp.json()
					if "esearchresult" not in data:
						raise ValueError("Incomplete JSON response: missing 'esearchresult'")
					page = data["esearchresult"]["idlist"]
					break
				except (requests.exceptions.JSONDecodeError, requests.exceptions.ChunkedEncodingError, requests.exceptions.ConnectionError, requests.exceptions.HTTPError, ValueError) as e:
					if attempt == max_retries - 1:
						print(f"Failed to fetch accession versions for chunk starting at {start} after {max_retries} attempts.")
						raise e

					is_429 = isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 429
					wait_time = (self.sleep_time * (attempt + 1)) + (10 if is_429 else 0)
					print(f"Error fetching accession versions ({type(e).__name__}) for chunk starting at {start}. Retrying in {wait_time}s... (Attempt {attempt+1}/{max_retries})")
					sleep(wait_time)

			all_accs.extend(page)
			print(f"[fetch_accs] chunk {start:,}-{min(start+self.batch_size, count):,} ({len(page):,} records) took {time.time() - t1:.1f}s")

		print(f"[fetch_accs] total accs fetched: {len(all_accs):,} in {time.time() - start_all:.1f}s")
		return all_accs
	def iter_accs(self):
			t0 = time.time()
			hist_url = (
				f"{self.base_url}esearch.fcgi?db=nucleotide"
				f"&term=txid{self.taxid}[Organism:exp]"
				f"&retmax=0&idtype=acc"
				f"&usehistory=y&email={self.email}&retmode=json"
			)
			
			# 1) Wrap the initial history call in a retry loop
			max_retries = 5
			hist = None
			for attempt in range(max_retries):
				try:
					resp = self._get(hist_url)
					resp.raise_for_status()
					
					# Sanitise control characters if any exist before decoding JSON
					# This replaces unescaped control chars with an empty string
					sanitised_text = "".join(ch for ch in resp.text if ord(ch) >= 32 or ch in "\n\r\t")
					
					import json
					data = json.loads(sanitised_text)
					hist = data["esearchresult"]
					break
				except (requests.exceptions.JSONDecodeError, requests.exceptions.ChunkedEncodingError, 
						requests.exceptions.ConnectionError, requests.exceptions.HTTPError, KeyError, ValueError) as e:
					if attempt == max_retries - 1:
						print(f"Failed to initialize ESearch history after {max_retries} attempts.")
						raise e
					
					is_429 = isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 429
					wait_time = (self.sleep_time * (attempt + 1)) + (10 if is_429 else 0)
					print(f"Error initializing history ({type(e).__name__}). Retrying in {wait_time}s... (Attempt {attempt+1}/{max_retries})")
					sleep(wait_time)

			count = int(hist["count"])
			print(f"[fetch_accs] ESearch history → count={count:,} took {time.time() - t0:.1f}s")

			if count == 0:
				return

			webenv = hist["webenv"]
			querykey = hist["querykey"]

			# 2) Pagination loop continues unchanged...
			for start in range(0, count, self.batch_size):
				t1 = time.time()
				page_url = (
					f"{self.base_url}esearch.fcgi?db=nucleotide"
					f"&WebEnv={webenv}&query_key={querykey}"
					f"&retstart={start}&retmax={self.batch_size}"
					f"&idtype=acc&retmode=json"
				)

				max_retries = 5
				for attempt in range(max_retries):
					try:
						page_resp = self._get(page_url)
						page_resp.raise_for_status()
						
						# Apply identical sanitisation to the chunk results
						sanitised_page = "".join(ch for ch in page_resp.text if ord(ch) >= 32 or ch in "\n\r\t")
						import json
						data = json.loads(sanitised_page)
						
						if "esearchresult" not in data:
							raise ValueError("Incomplete JSON response: missing 'esearchresult'")
						page = data["esearchresult"]["idlist"]
						break
					except (requests.exceptions.JSONDecodeError, requests.exceptions.ChunkedEncodingError, requests.exceptions.ConnectionError, requests.exceptions.HTTPError, ValueError) as e:
						if attempt == max_retries - 1:
							print(f"Failed to fetch accession versions for chunk starting at {start} after {max_retries} attempts.")
							raise e

						is_429 = isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 429
						wait_time = (self.sleep_time * (attempt + 1)) + (10 if is_429 else 0)
						print(f"Error fetching accession versions ({type(e).__name__}) for chunk starting at {start}. Retrying in {wait_time}s... (Attempt {attempt+1}/{max_retries})")
						sleep(wait_time)

				print(f"[fetch_accs] chunk {start:,}-{min(start+self.batch_size, count):,} ({len(page):,} records) took {time.time() - t1:.1f}s")
				yield page


	def fetch_genbank_data(self, ids):
		# Use a separate internal batch size just for efetch
		batch_n = self.efetch_batch_size
		if self.test_run and not self.update_file:
			ids = random.sample(ids, k=min(self.test_fetch_count, len(ids)))
			batch_n=10

		if self.ref_list and not self.update_file:
			try:
				ref_list = load_reference_file_accessions(self.ref_list)
				ids.extend(ref_list)
				print(f"Added {len(ref_list)} IDs from reference list.")
			except Exception as e:
				print(f"Warning: Could not read reference list: {e}")
		
		# Remove duplicates
		ids = list(set(ids))
	
		max_ids=len(ids)
		for i in range(0, max_ids, batch_n):
			chunk = ids[i:i+batch_n]
			ids_str = ",".join(chunk)
			url = (
				f"{self.base_url}efetch.fcgi?db=nucleotide"
				f"&id={ids_str}"
				f"&retmode=xml&email={self.email}"
			)
			
			# Retry logic for network interruptions
			max_retries = 5
			for attempt in range(max_retries):
				try:
					resp = self._get(url)
					resp.raise_for_status()

					# Pass batch_n here to name the file as before
					self.save_data(resp.text, i + batch_n)
					break # Success, exit retry loop
				except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ConnectionError, requests.exceptions.HTTPError) as e:
					if attempt == max_retries - 1:
						print(f"Failed to download batch starting with index {i} after {max_retries} attempts.")
						raise e
					
					# Check if it's a 429 error and wait longer if so
					is_429 = isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 429
					wait_time = (self.sleep_time * (attempt + 1)) + (10 if is_429 else 0)

					print(f"Network error ({type(e).__name__}) for batch starting with index {i}. Retrying in {wait_time}s... (Attempt {attempt+1}/{max_retries})")
					sleep(wait_time)

			print(f"Downloaded XML {i+1:,}–{min(i+batch_n, len(ids)):,}")
			sleep(self.sleep_time)



	def save_data(self, data, batch_size):
		os.makedirs(self.output_dir, exist_ok=True)
		os.makedirs(join(self.output_dir, self.base_dir), exist_ok=True)

		base_filename = join(self.output_dir, self.base_dir, f"batch-{batch_size}")
		filename = f"{base_filename}.xml"
		counter = 1
		while os.path.exists(filename):
			filename = f"{base_filename}_{counter}.xml"
			counter += 1

		with open(filename, "w") as file:
				file.write(data)
		print(f"Data written to: {filename}")

	def download(self):
		if self.test_run:
			try:
				search_url = (
					f"{self.base_url}esearch.fcgi?db=nucleotide"
					f"&term=txid{self.taxid}[Organism:exp]"
					f"&retmax={self.test_fetch_count}&idtype=acc"
					f"&usehistory=y&email={self.email}&retmode=json"
				)
				response = self._get(search_url)
				response.raise_for_status()
				data = response.json()
				ids = data["esearchresult"]["idlist"]
			except Exception as e:
				print(f"Warning: Could not fetch IDs for test run: {e}")
				ids = []
		else:
			ids = self.fetch_ids()
			
		print(f"Found {len(ids)} IDs")
		self.fetch_genbank_data(ids)

	def _detect_meta_data_acc_col(self, conn, preferred='accession_version'):
		cursor = conn.cursor()
		cursor.execute("PRAGMA table_info(meta_data)")
		cols = [row[1] for row in cursor.fetchall()]

		if preferred in cols:
			return preferred

		fallbacks = ['primary_accession', 'locus', 'accession']
		for col in fallbacks:
			if col in cols:
				print(f"[db] preferred column {preferred!r} not found; using {col!r} instead")
				return col

		raise ValueError(
			"Could not find an accession column in meta_data. "
			f"Tried preferred={preferred!r} and fallbacks={fallbacks}. Found columns={cols}"
		)

	def _split_accession_version(self, accession_version):
		if accession_version is None:
			return None, None
		acc_ver = str(accession_version).strip()
		if not acc_ver:
			return None, None

		base, sep, version = acc_ver.rpartition('.')
		if not sep:
			return acc_ver, None

		try:
			return base, int(version)
		except ValueError:
			return base, None

	def _canonical_exclusion_keys(self, values):
		"""Bare identity keys for the curator's exclusion list.

		The exclusion set is matched against the *bare* accession parsed out of
		an NCBI id, but the rows feeding it come from whichever meta_data column
		:meth:`_detect_meta_data_acc_col` picked - and the default, and the
		column every shipped DB actually has, is ``accession_version``, which
		holds ``PV547728.1``. A bare key is never equal to a versioned string,
		so the guard was dead code: a record the curator had deliberately
		excluded got silently re-downloaded the moment GenBank revised it.

		It only ever bit on a revision because an *unrevised* excluded record is
		caught earlier by the ``local_exact`` test - which means it bit on
		exactly the run where honouring the exclusion matters, and looked fine
		on every other run.

		Normalising at load rather than at the comparison keeps this set
		canonical whichever column was used, and matches the
		``excluded_accessions`` table branch, which already stores bare
		``primary_accession`` values. :mod:`accession_utils` is the authority on
		what "bare" means: a strain name or a cluster label that happens to
		contain a dot is left whole rather than truncated by ``split('.')[0]``.
		"""
		keys = set()
		for value in values:
			key = accession_utils.normalise_accession(value)
			if key:
				keys.add(key)
		return keys

	def _load_db_accession_context(self, db_path, meta_acc_col='accession_version'):
		if not db_path:
			raise ValueError('Update DB path not provided')
		if not os.path.exists(db_path):
			raise FileNotFoundError(f"SQLite DB not found: {db_path}")

		conn = sqlite3.connect(db_path)
		try:
			cursor = conn.cursor()
			cursor.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='meta_data' LIMIT 1")
			if cursor.fetchone() is None:
				raise ValueError('SQLite DB must contain table: meta_data')

			acc_col = self._detect_meta_data_acc_col(conn, preferred=meta_acc_col)
			cursor.execute("PRAGMA table_info(meta_data)")
			meta_cols = {row[1] for row in cursor.fetchall()}
			cursor.execute(f"SELECT {acc_col} FROM meta_data WHERE {acc_col} IS NOT NULL")
			meta_accessions = [row[0] for row in cursor.fetchall() if row[0]]

			excluded_primary = set()
			if 'exclusion_status' in meta_cols:
				cursor.execute(
					f"SELECT {acc_col} FROM meta_data WHERE {acc_col} IS NOT NULL "
					"AND LOWER(TRIM(COALESCE(CAST(exclusion_status AS TEXT), ''))) NOT IN ('', '0', 'false', 'no', 'na', 'none', 'nan')"
				)
				excluded_primary.update(self._canonical_exclusion_keys(row[0] for row in cursor.fetchall()))
			cursor.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='excluded_accessions' LIMIT 1")
			if cursor.fetchone() is not None:
				cursor.execute("PRAGMA table_info(excluded_accessions)")
				excl_cols = {row[1] for row in cursor.fetchall()}
				if 'primary_accession' in excl_cols:
					cursor.execute("SELECT primary_accession FROM excluded_accessions WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) != ''")
					excluded_primary = self._canonical_exclusion_keys(row[0] for row in cursor.fetchall())
		finally:
			conn.close()

		print(f"[db] using meta_data column: {acc_col}")
		print("first 10 values in DB:", meta_accessions[:10])
		return meta_accessions, excluded_primary

	def _load_ref_accessions(self):
		if not self.ref_list or not os.path.exists(self.ref_list):
			return []
		return sorted(set(load_reference_file_accessions(self.ref_list, allowed_types={"master", "reference"})))

	def _write_update_log(self, updated_versions, new_accessions):
		if not self.update_log:
			return

		log_path = join(self.output_dir, self.update_log)
		with open(log_path, 'w', newline='', encoding='utf-8') as file:
			writer = csv.writer(file, delimiter='\t')
			writer.writerow(['old_accession_version', 'new_accession_version'])
			for old_acc, new_acc in updated_versions:
				writer.writerow([old_acc, new_acc])
			for new_acc in new_accessions:
				writer.writerow(['NA', new_acc])
		print(f"[update] Wrote update log: {log_path}")

	def _compute_missing_ids(self, ncbi_ids, meta_accessions, excluded_primary):
		# Both sides of the exclusion test are reduced to the bare identity key.
		# local_exact and local_versions below deliberately keep the versioned
		# spelling: that is what detects a revised record.
		excluded_bare = self._canonical_exclusion_keys(excluded_primary or ())
		local_exact = set(str(x).strip() for x in meta_accessions if x is not None)
		local_versions = {}
		for acc_ver in meta_accessions:
			base, version = self._split_accession_version(acc_ver)
			if not base:
				continue
			if version is None:
				local_versions.setdefault(base, None)
				continue
			if base not in local_versions:
				local_versions[base] = version
				continue
			prev = local_versions.get(base)
			if prev is None:
				continue
			if prev is not None and version > prev:
				local_versions[base] = version

		updated_versions = []
		new_accessions = []
		missing_ids = []

		for ncbi_acc in ncbi_ids:
			ncbi_acc = str(ncbi_acc).strip()
			if not ncbi_acc:
				continue

			if ncbi_acc in local_exact:
				continue

			base, version = self._split_accession_version(ncbi_acc)
			if not base:
				continue

			if accession_utils.normalise_accession(base) in excluded_bare:
				continue

			if base in local_versions:
				local_ver = local_versions[base]
				if local_ver is None:
					continue
				if version is not None and version > local_ver:
					old_acc = f"{base}.{local_ver}" if local_ver > 0 else base
					updated_versions.append((old_acc, ncbi_acc))
					missing_ids.append(ncbi_acc)
				continue

			if base not in local_versions:
				new_accessions.append(ncbi_acc)
				missing_ids.append(ncbi_acc)

		return missing_ids, updated_versions, new_accessions

	def update(self, db_path, meta_acc_col='accession_version'):
		meta_accessions, excluded_primary = self._load_db_accession_context(db_path, meta_acc_col=meta_acc_col)
		print(f"Found {len(meta_accessions)} local accessions (from DB)")

		if self.test_run:
			ids = []
			updated_versions = []
			new_accessions = []
			missing_ids = []
			limit = self.test_fetch_count
			for page in self.iter_accs():
				ids.extend(page)
				page_missing, page_updates, page_new = self._compute_missing_ids(page, meta_accessions, excluded_primary)
				updated_versions.extend(page_updates)
				new_accessions.extend(page_new)
				missing_ids.extend(page_missing)
				if len(missing_ids) >= limit:
					missing_ids = missing_ids[:limit]
					break
		else:
			ids = self.fetch_accs()
			missing_ids, updated_versions, new_accessions = self._compute_missing_ids(ids, meta_accessions, excluded_primary)

		print("first 10 IDs in NCBI:", ids[:10])

		# Update mode is DB-backed for references: do not pull reference/master accessions from GenBank.
		if self.ref_list:
			print("[update] Reference list supplied, but reference/master accessions are sourced from DB only in update mode")

		if updated_versions:
			print('Updated accessions found (old -> new):')
			for old_acc, new_acc in updated_versions[:50]:
				print(f"  {old_acc} -> {new_acc}")
			if len(updated_versions) > 50:
				print(f"  ... ({len(updated_versions)-50} more)")

		print(f"Found {len(missing_ids)} missing/new IDs total (updates={len(updated_versions)}, brand-new={len(new_accessions)})")

		if self.test_run and missing_ids:
			if new_accessions:
				missing_ids = list(dict.fromkeys(new_accessions))[:self.test_fetch_count]
				print(f"[update][test_run] selected {len(missing_ids)} brand-new accessions not present in DB")
			else:
				missing_ids = list(dict.fromkeys(missing_ids))[:self.test_fetch_count]
				print(f"[update][test_run] no brand-new accessions found; selected {len(missing_ids)} update candidates")

		self._write_update_log(updated_versions, new_accessions)
		if missing_ids:
			self.fetch_genbank_data(missing_ids)
		else:
			print('[update] Nothing to download.')

# still appending references to fetch list in update mode, which is probably stupid if there's loads?
# better method- see if references have been updated / any differences between refs nominated in db vs file
# then only add to fetch list if not already present in db or if updated version is newer than db version (if version info available)

if __name__ == "__main__":
	parser = ArgumentParser(description='This script downloads and updates GenBank XML files for a given taxonomic group identified by an NCBI taxonomic ID, usually a species')
	parser.add_argument('-t', '--taxid', help='TaxID example: 11292', required=True)
	parser.add_argument('-o', '--tmp_dir', help='Output directory where all the XML files are stored', default='tmp')
	parser.add_argument('-b', '--batch_size', help='Max number of XML files to pull and merge in a single file', default=100000, type=int)
	parser.add_argument('--update', help='Path to existing SQLite DB generated by CreateSqliteDB.py; downloads only updated/new accessions not present or excluded in DB')
	parser.add_argument('--meta_acc_col', default='accession_version', help='Column in meta_data used for version-aware update checks (default: accession_version)')
	parser.add_argument('--update_log', default='update_accessions.tsv', help='TSV file (under tmp_dir) to write old->new and new accession mappings')
	parser.add_argument('-u', '--base_url', help='Base URL to download the XML files', default='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/')
	parser.add_argument('-e', '--email', help='Email ID', default='your_email@example.com')
	parser.add_argument('-s', '--sleep_time', help='Delay after each set of information fetch', default=2, type=int)
	parser.add_argument('-d', '--base_dir', help='Directory where all the XML files are stored', default='GenBank-XML')
	parser.add_argument("--test_run", action="store_true", help="Run a test fetching only a few records for quick testing")
	parser.add_argument('--test_fetch_count', default=100, type=int, help='Number of accessions to fetch in test mode (default: 100)')
	parser.add_argument('--ref_list', help='Reference accession list for test run', default=None)
	args = parser.parse_args()

	fetcher = GenBankFetcher(
		taxid=args.taxid,
		base_url=args.base_url,
		email=args.email,
		output_dir=args.tmp_dir,
		batch_size=args.batch_size,
		sleep_time=args.sleep_time,
		base_dir = args.base_dir,
		update_file = args.update,
		test_run = args.test_run,
		ref_list = args.ref_list,
		test_fetch_count = args.test_fetch_count
	)
	fetcher.update_log = args.update_log

	if args.update:
		fetcher.update(args.update, meta_acc_col=args.meta_acc_col)
	else:
		fetcher.download()

