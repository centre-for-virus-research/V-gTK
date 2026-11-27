import os
import csv
import json
import sqlite3
import time
import random
import inspect
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
from os.path import join
from urllib.parse import quote_plus
from argparse import ArgumentParser
from ExportRefListFromUpdateDb import load_reference_file_accessions
import accession_utils

try:
	import requests
except ModuleNotFoundError:
	class _RequestsFallback:
		class exceptions:
			class RequestException(Exception):
				pass

			class JSONDecodeError(RequestException):
				pass

			class ChunkedEncodingError(RequestException):
				pass

			class ConnectionError(RequestException):
				pass

			class Timeout(RequestException):
				pass

			class HTTPError(RequestException):
				def __init__(self, *args, response=None, **kwargs):
					super().__init__(*args)
					self.response = response

		@staticmethod
		def get(*args, **kwargs):
			raise ModuleNotFoundError("requests is required for network operations")

	requests = _RequestsFallback()

class GenBankFetcher:
	def __init__(self, taxid, base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/', email='your_email@example.com', output_dir='tmp', batch_size=100000, sleep_time=2, base_dir='GenBank-XML', update_file=None, test_run=False, ref_list=None, test_fetch_count=100, entrez_filter=None, efetch_batch_size=100, efetch_workers=1):
		self.taxid = taxid
		self.base_url = base_url
		self.email = email
		self.output_dir = output_dir
		self.batch_size = batch_size
		self.efetch_batch_size = max(1, int(efetch_batch_size))
		self.efetch_workers = max(1, int(efetch_workers))
		self.sleep_time = sleep_time
		self.base_dir = base_dir
		self.update_file = update_file
		self.test_run = test_run
		self.ref_list = ref_list
		self.test_fetch_count = test_fetch_count
		self.entrez_filter = str(entrez_filter).strip() if entrez_filter and str(entrez_filter).strip().lower() not in ("null", "none", "unset") else None
		self.update_log = "update_accessions.tsv"
		self._last_request = 0.0
		self._throttle_lock = threading.Lock()
		self._save_lock = threading.Lock()
		
		# Load NCBI API key if available and not in pytest
		import sys
		if "pytest" not in sys.modules:
			self.api_key = self._load_ncbi_key()
		else:
			self.api_key = None

		# Pace outgoing requests under NCBI's limit: 10 req/s with an API key,
		# 3 req/s without. Disabled under pytest so unit tests stay fast.
		self.min_request_interval = 0.0 if "pytest" in sys.modules else (0.11 if self.api_key else 0.36)
		
		# Ensure output directories exist immediately
		os.makedirs(self.output_dir, exist_ok=True)
		os.makedirs(join(self.output_dir, self.base_dir), exist_ok=True)

	def _load_ncbi_key(self):
		# Resolved relative to this script so the key travels with the checkout
		# rather than being pinned to one server's home directory. Without a key
		# NCBI throttles by IP at 3 req/s instead of 10 req/s.
		env_key = os.environ.get("NCBI_API_KEY", "").strip()
		if env_key:
			if self._valid_ncbi_key(env_key):
				return env_key
			print("Warning: NCBI_API_KEY is set but does not look like a valid key; ignoring")

		repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		candidates = [
			join(repo_root, "dev", "ncbi_key"),
			join(repo_root, "ncbi_key"),
			os.path.expanduser("~/.ncbi_key"),
		]
		for key_path in candidates:
			if not os.path.exists(key_path):
				continue
			try:
				with open(key_path, "r", encoding="utf-8") as f:
					key = f.read().strip()
			except Exception as e:
				print(f"Warning: Could not read NCBI key from {key_path}: {e}")
				continue
			if self._valid_ncbi_key(key):
				print(f"[ncbi] using API key from {key_path} (10 req/s)")
				return key
			print(f"Warning: NCBI key in {key_path} is malformed; ignoring")

		print("Warning: no NCBI API key found (looked at $NCBI_API_KEY, "
			f"{candidates[0]}). Falling back to 3 req/s IP-based rate limiting.")
		return None

	@staticmethod
	def _valid_ncbi_key(key):
		return bool(key) and len(key) == 36 and all(c.isalnum() or c in "-_" for c in key)

	def _get(self, url):
		if self.api_key:
			if "?" in url:
				url += f"&api_key={self.api_key}"
			else:
				url += f"?api_key={self.api_key}"
		self._throttle()
		sig = inspect.signature(requests.get)
		accepts_kwargs = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values())
		if accepts_kwargs or "timeout" in sig.parameters:
			return requests.get(url, timeout=(20, 180))
		return requests.get(url)

	@staticmethod
	def _is_429(exc):
		resp = getattr(exc, "response", None)
		return resp is not None and getattr(resp, "status_code", None) == 429

	def _esearch_json_with_retry(self, url, context, max_retries=5):
		timeout_exc = getattr(requests.exceptions, "Timeout", TimeoutError)
		retryable = (
			requests.exceptions.JSONDecodeError,
			requests.exceptions.ChunkedEncodingError,
			requests.exceptions.ConnectionError,
			requests.exceptions.HTTPError,
			timeout_exc,
			KeyError,
			ValueError,
		)
		for attempt in range(max_retries):
			try:
				resp = self._get(url)
				resp.raise_for_status()
				data = self._parse_json(resp, context)
				if "esearchresult" not in data:
					raise ValueError("Incomplete JSON response: missing 'esearchresult'")
				return data
			except retryable as e:
				if attempt == max_retries - 1:
					raise
				wait_time = (self.sleep_time * (attempt + 1)) + (10 if self._is_429(e) else 0)
				print(f"Error in {context} ({type(e).__name__}). Retrying in {wait_time}s... (Attempt {attempt+1}/{max_retries})")
				sleep(wait_time)

	def _throttle(self):
		# Pace requests ourselves instead of discovering NCBI's limit via 429s,
		# which cost a full backoff each time they are hit. Locked so concurrent
		# efetch workers share one global request budget.
		if self.min_request_interval <= 0:
			return
		with self._throttle_lock:
			elapsed = time.time() - self._last_request
			if elapsed < self.min_request_interval:
				sleep(self.min_request_interval - elapsed)
			self._last_request = time.time()

	def _search_term_encoded(self):
		base_term = f"txid{self.taxid}[Organism:exp]"
		if not self.entrez_filter:
			return base_term
		return quote_plus(f"{base_term} AND ({self.entrez_filter})")

	def _parse_json(self, resp, context):
		"""Decode an E-utilities JSON body, tolerating stray control characters.

		NCBI intermittently returns bodies with unescaped control characters that
		json.loads rejects. Strip them and retry before giving up, and log the raw
		body so a failure is diagnosable rather than an opaque JSONDecodeError.
		"""
		try:
			return resp.json()
		except ValueError:
			pass

		sanitised = "".join(ch for ch in resp.text if ord(ch) >= 32)
		try:
			data = json.loads(sanitised)
			print(f"[ncbi] recovered malformed JSON for {context} by stripping control characters")
			return data
		except ValueError as e:
			body = resp.text[:500].replace("\n", "\\n")
			print(f"[ncbi] undecodable response for {context} "
				f"(status={resp.status_code}, {len(resp.text)} bytes): {body!r}")
			raise ValueError(f"Undecodable JSON response for {context}: {e}") from e

	def get_record_count(self):
		search_url = f"{self.base_url}esearch.fcgi?db=nucleotide&term={self._search_term_encoded()}&retmode=json&email={self.email}"
		response = self._get(search_url)
		response.raise_for_status()
		data = self._parse_json(response, "record count")
		return int(data['esearchresult']['count'])

	def fetch_ids(self):
		start_all = time.time()

		# 1) Initialize search history and retrieve WebEnv, QueryKey and count
		t0 = time.time()
		hist_url = (
			f"{self.base_url}esearch.fcgi?db=nucleotide"
			f"&term={self._search_term_encoded()}"
			f"&retmax=0&idtype=acc"
			f"&usehistory=y&email={self.email}&retmode=json"
		)
		hist = self._esearch_json_with_retry(hist_url, "fetch_ids ESearch history")["esearchresult"]
		
		count = int(hist["count"])
		print(f"[fetch_ids] ESearch history → count={count:,} took {time.time() - t0:.1f}s")
		
		if count == 0:
			return []

		webenv   = hist["webenv"]
		querykey = hist["querykey"]

		# 2) Paginate in chunks, advancing by actual returned page length
		all_accs = []
		start = 0
		while start < count:
			t1 = time.time()
			page_url = (
				f"{self.base_url}esearch.fcgi?db=nucleotide"
				f"&WebEnv={webenv}&query_key={querykey}"
				f"&retstart={start}&retmax={self.batch_size}"
				f"&idtype=acc&retmode=json"
			)
			data = self._esearch_json_with_retry(page_url, f"fetch_ids chunk at retstart={start}")
			page = data["esearchresult"]["idlist"]
			if not page:
				break
			page_bare = [acc.split('.', 1)[0] for acc in page]
			all_accs.extend(page_bare)
			print(f"[fetch_ids] chunk {start:,}-{min(start+len(page), count):,} "
				f"({len(page):,} records) took {time.time() - t1:.1f}s")
			start += len(page)

		print(f"[fetch_ids] total accs fetched: {len(all_accs):,} in {time.time() - start_all:.1f}s")
		return all_accs

	def fetch_accs(self):
		start_all = time.time()

		t0 = time.time()
		hist_url = (
			f"{self.base_url}esearch.fcgi?db=nucleotide"
			f"&term={self._search_term_encoded()}"
			f"&retmax=0&idtype=acc"
			f"&usehistory=y&email={self.email}&retmode=json"
		)
		hist = self._esearch_json_with_retry(hist_url, "fetch_accs ESearch history")["esearchresult"]

		count = int(hist["count"])
		print(f"[fetch_accs] ESearch history → count={count:,} took {time.time() - t0:.1f}s")

		if count == 0:
			return []

		webenv = hist["webenv"]
		querykey = hist["querykey"]

		all_accs = []
		start = 0
		while start < count:
			t1 = time.time()
			page_url = (
				f"{self.base_url}esearch.fcgi?db=nucleotide"
				f"&WebEnv={webenv}&query_key={querykey}"
				f"&retstart={start}&retmax={self.batch_size}"
				f"&idtype=acc&retmode=json"
			)
			data = self._esearch_json_with_retry(page_url, f"fetch_accs chunk at retstart={start}")
			page = data["esearchresult"]["idlist"]
			if not page:
				break
			all_accs.extend(page)
			print(f"[fetch_accs] chunk {start:,}-{min(start+len(page), count):,} ({len(page):,} records) took {time.time() - t1:.1f}s")
			start += len(page)

		print(f"[fetch_accs] total accs fetched: {len(all_accs):,} in {time.time() - start_all:.1f}s")
		return all_accs
	def iter_accs(self):
		t0 = time.time()
		hist_url = (
			f"{self.base_url}esearch.fcgi?db=nucleotide"
			f"&term={self._search_term_encoded()}"
			f"&retmax=0&idtype=acc"
			f"&usehistory=y&email={self.email}&retmode=json"
		)
		hist = self._esearch_json_with_retry(hist_url, "iter_accs ESearch history")["esearchresult"]

		count = int(hist["count"])
		print(f"[fetch_accs] ESearch history → count={count:,} took {time.time() - t0:.1f}s")

		if count == 0:
			return

		webenv = hist["webenv"]
		querykey = hist["querykey"]

		start = 0
		while start < count:
			t1 = time.time()
			page_url = (
				f"{self.base_url}esearch.fcgi?db=nucleotide"
				f"&WebEnv={webenv}&query_key={querykey}"
				f"&retstart={start}&retmax={self.batch_size}"
				f"&idtype=acc&retmode=json"
			)
			data = self._esearch_json_with_retry(page_url, f"iter_accs chunk at retstart={start}")
			page = data["esearchresult"]["idlist"]
			if not page:
				break
			print(f"[fetch_accs] chunk {start:,}-{min(start+len(page), count):,} ({len(page):,} records) took {time.time() - t1:.1f}s")
			yield page
			start += len(page)


	def _download_one_efetch_batch(self, i, chunk, batch_n, total_ids):
		ids_str = ",".join(chunk)
		url = (
			f"{self.base_url}efetch.fcgi?db=nucleotide"
			f"&id={ids_str}"
			f"&retmode=xml&email={self.email}"
		)
		max_retries = 8
		timeout_exc = getattr(requests.exceptions, "Timeout", TimeoutError)
		req_exc = getattr(requests.exceptions, "RequestException", Exception)
		for attempt in range(max_retries):
			try:
				resp = self._get(url)
				resp.raise_for_status()
				body = resp.text or ""
				if body.strip() and "</GBSet>" not in body:
					raise ValueError("Truncated GenBank XML response: missing </GBSet>")
				self.save_data(resp.text, i + batch_n)
				break
			except (
				requests.exceptions.ChunkedEncodingError,
				requests.exceptions.ConnectionError,
				requests.exceptions.HTTPError,
				timeout_exc,
				req_exc,
				ValueError,
			) as e:
				if attempt == max_retries - 1:
					print(f"Failed to download batch starting with index {i} after {max_retries} attempts.")
					raise e
				status_code = getattr(getattr(e, "response", None), "status_code", None)
				is_backoff_status = status_code in (429, 500, 502, 503, 504)
				wait_time = (max(1, self.sleep_time) * (attempt + 1)) + (10 if is_backoff_status else 0)
				print(f"Network error ({type(e).__name__}, status={status_code}) for batch starting with index {i}. Retrying in {wait_time}s... (Attempt {attempt+1}/{max_retries})")
				sleep(wait_time)
		print(f"Downloaded XML {i+1:,}–{min(i+batch_n, total_ids):,}")
		if self.sleep_time > 0:
			sleep(self.sleep_time)

	def fetch_genbank_data(self, ids):
		# Use a separate internal batch size just for efetch
		batch_n = self.efetch_batch_size
		if self.test_run and not self.update_file:
			ids = random.sample(ids, k=min(self.test_fetch_count, len(ids)))
			if self.efetch_batch_size == 100 and self.test_fetch_count <= 100:
				batch_n = 10

		if self.ref_list and not self.update_file:
			try:
				ref_list = load_reference_file_accessions(self.ref_list)
				insdc_refs = [acc for acc in ref_list if not accession_utils.is_gisaid_accession(acc)]
				skipped_refs = [acc for acc in ref_list if accession_utils.is_gisaid_accession(acc)]
				ids.extend(insdc_refs)
				print(f"Added {len(insdc_refs)} IDs from reference list.")
				if skipped_refs:
					print(f"[info] Skipped {len(skipped_refs)} GISAID reference ID(s) during GenBank fetch (e.g. {', '.join(skipped_refs[:5])}).")
			except Exception as e:
				print(f"Warning: Could not read reference list: {e}")
		
		# Remove duplicates deterministically
		ids = list(dict.fromkeys(ids))
	
		max_ids = len(ids)
		batches = [(i, ids[i:i+batch_n]) for i in range(0, max_ids, batch_n)]
		if self.efetch_workers <= 1 or len(batches) <= 1:
			for i, chunk in batches:
				self._download_one_efetch_batch(i, chunk, batch_n, max_ids)
		else:
			print(f"[efetch] Downloading {max_ids:,} IDs across {len(batches):,} batches with {self.efetch_workers} workers (batch_size={batch_n})")
			with ThreadPoolExecutor(max_workers=self.efetch_workers) as pool:
				futures = [pool.submit(self._download_one_efetch_batch, i, chunk, batch_n, max_ids) for i, chunk in batches]
				for fut in as_completed(futures):
					fut.result()

	def save_data(self, data, batch_size):
		with self._save_lock:
			os.makedirs(self.output_dir, exist_ok=True)
			os.makedirs(join(self.output_dir, self.base_dir), exist_ok=True)

			base_filename = join(self.output_dir, self.base_dir, f"batch-{batch_size}")
			filename = f"{base_filename}.xml"
			counter = 1
			while os.path.exists(filename):
				filename = f"{base_filename}_{counter}.xml"
				counter += 1

			tmp_filename = f"{filename}.tmp"
			with open(tmp_filename, "w", encoding="utf-8") as file:
				file.write(data)
			os.replace(tmp_filename, filename)
			print(f"Data written to: {filename}")

	def download(self):
		if self.test_run:
			try:
				search_url = (
					f"{self.base_url}esearch.fcgi?db=nucleotide"
					f"&term={self._search_term_encoded()}"
					f"&retmax={self.test_fetch_count}&idtype=acc"
					f"&usehistory=y&email={self.email}&retmode=json"
				)
				response = self._get(search_url)
				response.raise_for_status()
				data = self._parse_json(response, "test-run ESearch")
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
					excluded_primary.update(self._canonical_exclusion_keys(row[0] for row in cursor.fetchall()))
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
	parser.add_argument('--entrez_filter', default=None, help='Optional Entrez query filter ANDed with txid<TAXID>[Organism:exp] (e.g. "29500:30200[SLEN] NOT patent[Division] NOT synthetic[Division]")')
	parser.add_argument('--efetch_batch_size', default=100, type=int, help='Number of accessions per efetch XML request (default: 100)')
	parser.add_argument('--efetch_workers', default=1, type=int, help='Concurrent efetch worker threads sharing the NCBI rate limiter (default: 1)')
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
		test_fetch_count = args.test_fetch_count,
		entrez_filter = args.entrez_filter,
		efetch_batch_size = args.efetch_batch_size,
		efetch_workers = args.efetch_workers,
	)
	fetcher.update_log = args.update_log

	if args.update:
		fetcher.update(args.update, meta_acc_col=args.meta_acc_col)
	else:
		fetcher.download()

