import os
import re
import csv
import shlex
import subprocess
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
from FastaHandler import RemoveRedundantSequence
from ExportRefListFromUpdateDb import load_master_accessions, load_master_accessions_from_file
import accession_utils
import segment_utils


class NextalignNotAvailable(RuntimeError):
	"""The ``nextalign`` executable is not on PATH.

	Kept distinct from an alignment failure so the operator is pointed at the
	conda environment rather than at the sequence data. Running through the
	shell used to turn this into exit status 127, which the relaxation loop
	retried five times and then reported as "extreme sequence divergence".
	"""


class NextalignAlignment:
	#: Bases that can be compared position by position. Everything else - the
	#: gap character and every IUPAC ambiguity code - means "not known here",
	#: which is not the same as "does not match". Scoring ``N`` as a mismatch
	#: made a sequence whose called bases all agree with the reference score 0%
	#: on a poly-N 3' end, and via :meth:`nextalign_master` that aborted the run.
	COMPARABLE_BASES = frozenset('ACGT')

	#: A query is judged against the final 30% of its aligned bases, where the
	#: conserved downstream genes sit.
	CONSERVED_WINDOW_FRACTION = 0.70

	#: Below this identity in that window the alignment is treated as
	#: frame-shifted.
	MIN_DOWNSTREAM_IDENTITY = 0.50

	def __init__(self, gb_matrix, query_dir, ref_dir, ref_fa_file, master_seq_dir, tmp_dir, master_ref, nextalign_dir, reference_alignment, update_db=None, max_threads=10):
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
		self.max_threads = max_threads

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
		"""Delegates to :mod:`accession_utils` - the single authority on
		accession identity.

		This used to be ``os.path.basename(path).split('.')[0]``, which
		truncates at the *first* dot rather than at a recognised file
		extension. On the data this pipeline actually ships that happens to
		give the right answer, because every reference list and every update DB
		holds the **bare** accession - so this was a latent hazard, not a live
		defect. It produced the wrong answer for four shapes that do occur:
		``NC_001542.1.fasta`` (a versioned download), ``X.aligned.fasta``, any
		dot-leading filename (which became ``""``, an empty accession used as a
		directory name and as a command-line value), and a group name
		containing a dot.

		The canonical form is bare, matching ``meta_data.primary_accession``.
		The versioned spelling lives only in ``meta_data.accession_version``,
		where GenBankFetcher uses it to spot a revised record.
		"""
		return accession_utils.accession_from_filename(file_path)

	@staticmethod
	def _normalize_segment(segment_value):
		"""Delegates to :mod:`segment_utils` - the single normalisation authority.

		This used to scrape every digit out of the string, which turned ``4.0`` into
		segment 40 and, worse, inverted the polymerase segments: ``PB2`` (segment 1)
		became ``2`` and ``PB1`` (segment 2) became ``1``. Missing still maps to
		``""`` here, because callers of this particular helper depend on it.
		"""
		normalised = segment_utils.normalise_segment(segment_value)
		if normalised is None or normalised.casefold() in segment_utils.PANDAS_NULL_TOKENS:
			return ""
		return normalised

	def get_master_list(self):
		"""Master accessions in canonical (bare) form.

		Normalised because the other side of the lookup in :meth:`process` is
		:meth:`path_to_basename`, which is also canonical. An update DB or a
		hand-written reference list carrying ``NC_001542.1`` would otherwise be
		compared against a directory named ``NC_001542`` and silently miss.
		"""
		return [
			accession_utils.normalise_accession(master)
			for master in self._raw_master_list()
			if accession_utils.normalise_accession(master)
		]

	def _raw_master_list(self):
		if self.update_db:
			masters = load_master_accessions(self.update_db)
			if masters:
				return masters
		if os.path.isfile(self.master_ref):
			try:
				return load_master_accessions_from_file(self.master_ref)
			except Exception as exc:
				# Was a bare `except:`, which also swallowed KeyboardInterrupt
				# and returned an empty master list without a word.
				print(f"WARNING: could not read master reference list {self.master_ref}: {exc}")
				return []
		else:
			return [x.strip() for x in self.master_ref.split(',') if x.strip()]

	@staticmethod
	def _listdir_files(directory):
		"""Real, non-hidden files in ``directory``.

		``os.listdir`` also returns editor swap files, ``.DS_Store`` and
		nextflow/snakemake bookkeeping entries, each of which used to be
		launched as its own alignment job.
		"""
		if not os.path.isdir(directory):
			print(f"WARNING: directory does not exist: {directory}")
			return []
		return sorted(
			name for name in os.listdir(directory)
			if not name.startswith('.') and os.path.isfile(join(directory, name))
		)

	@staticmethod
	def _accessions_in_fasta(fasta_path):
		"""Sequence ids in a FASTA, used to attribute a whole-group failure."""
		try:
			return [record.id for record in SeqIO.parse(fasta_path, "fasta")]
		except (OSError, ValueError) as exc:
			print(f"WARNING: could not read sequence ids from {fasta_path}: {exc}")
			return []

	def _validate_alignment_integrity(self, aligned_fasta_path):
		"""
		Evaluates the structural validity of the alignment by checking nucleotide
		identity in the final 30% of the alignment window (conserved downstream genes)
		for all query sequences in the alignment against the reference sequence.

		Returns ``(bool_passed, min_calculated_identity)``, where the identity is
		the worst score among the queries that could actually be scored. A query
		with nothing comparable in the window - because the reference is all gap
		there, or the query's 3' end is all ``N`` - is reported as unvalidated
		rather than as a total mismatch; those two cases used to return 0.0 and
		abort the pipeline.
		"""
		if not os.path.exists(aligned_fasta_path):
			return False, 0.0

		records = list(SeqIO.parse(aligned_fasta_path, "fasta"))
		if len(records) < 2:
			return False, 0.0

		ref_seq = str(records[0].seq).upper()
		total_columns = len(ref_seq)

		min_identity = 1.0
		all_passed = True

		for r_idx in range(1, len(records)):
			query_seq = str(records[r_idx].seq).upper()

			# An alignment is rectangular by definition. A record of a different
			# length means nextalign produced something malformed - the loop used
			# to index the query with the reference's length and raise IndexError.
			if len(query_seq) != total_columns:
				print(f"[validation] Sequence {records[r_idx].id} is {len(query_seq)} columns but reference {records[0].id} is {total_columns}. Alignment is not rectangular.")
				all_passed = False
				min_identity = 0.0
				continue

			non_gap_indices = [i for i in range(total_columns) if query_seq[i] != '-']
			if not non_gap_indices:
				print(f"[validation] Sequence {records[r_idx].id} is entirely gap against reference {records[0].id}.")
				all_passed = False
				min_identity = 0.0
				continue

			start_idx = int(len(non_gap_indices) * self.CONSERVED_WINDOW_FRACTION)
			target_indices = non_gap_indices[start_idx:]

			matches = 0
			valid_positions = 0

			for i in target_indices:
				n1 = ref_seq[i]
				n2 = query_seq[i]
				if n1 not in self.COMPARABLE_BASES or n2 not in self.COMPARABLE_BASES:
					continue
				valid_positions += 1
				if n1 == n2:
					matches += 1

			if valid_positions == 0:
				# Nothing to compare: an insertion past the reference's 3' end,
				# or an ambiguous tail. Not evidence of a frame shift.
				print(f"[validation] Sequence {records[r_idx].id} has no comparable bases in the conserved window against reference {records[0].id}; not validated.")
				continue

			identity = matches / valid_positions

			if identity < min_identity:
				min_identity = identity

			if identity < self.MIN_DOWNSTREAM_IDENTITY:
				all_passed = False
				print(f"[validation] Sequence {records[r_idx].id} failed frame validation against reference {records[0].id}. Downstream Identity: {identity:.2%}")

		return all_passed, min_identity

	def _run_nextalign(self, command, description):
		"""Execute nextalign without a shell and return its exit code.

		The argv list is handed to subprocess unchanged, matching every other
		external-tool call in this repo. It used to be flattened with
		``os.system(" ".join(command))``, which word-split any path containing a
		space, executed anything after a ``;``, dropped empty arguments
		entirely, and returned a wait status that was printed as if it were an
		exit code.
		"""
		print(f"{description}: {shlex.join(command)}")
		try:
			completed = subprocess.run(command, check=False)
		except FileNotFoundError as exc:
			raise NextalignNotAvailable(
				f"nextalign is not on PATH; cannot run: {shlex.join(command)}"
			) from exc
		return completed.returncode

	def _nextalign_command(self, profile, query_acc_path, ref_acc_path, output_subdir, accession):
		return [
			'nextalign', 'run',
			'--min-seeds', str(profile["min_seeds"]),
			'--seed-spacing', str(profile["seed_spacing"]),
			'--min-match-rate', str(profile["min_match_rate"]),
			'--input-ref', ref_acc_path,
			'--output-all', output_subdir,
			'--output-basename', f'{accession}',
			'--include-reference',
			'--jobs', str(self.max_threads),
			query_acc_path
		]

	def nextalign_master(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(ref_acc_path)
		output_subdir = join(query_aln_op, f'{accession}')
		aligned_fasta = join(output_subdir, f'{accession}.aligned.fasta')

		# Loop dynamically through relaxation parameters
		for attempt, profile in enumerate(self.relaxation_profiles, 1):
			command = self._nextalign_command(profile, query_acc_path, ref_acc_path, output_subdir, accession)
			return_code = self._run_nextalign(command, f"Executing Master Alignment (Attempt {attempt})")

			if return_code == 0:
				passed, identity = self._validate_alignment_integrity(aligned_fasta)
				if passed:
					print(f"SUCCESS: {accession} aligned cleanly (Downstream Identity: {identity:.2%}) using seeds={profile['min_seeds']}.")
					return
				else:
					print(f"WARNING: Nextalign completed for {accession} but failed frame validation. Downstream Identity: {identity:.2%}. Frame-shift suspected near hypervariable locus.")
			else:
				print(f"WARNING: Nextalign exited with code {return_code} for {accession} on attempt {attempt}.")

		# If all relaxation attempts fail, crash hard to preserve database parity
		print(f"\nCRITICAL ERROR: Fundamental alignment failure for reference sequence: {accession}")
		print(f"Reason: Out-of-frame insertions/deletions or extreme sequence divergence could not be resolved across all heuristic profiles.")
		print(f"Aborting execution to protect SQLite coordinate lookups.")
		raise RuntimeError(f"Alignment Frame Failure: {accession} aborted pipeline execution.")

	def nextalign_query(self, query_acc_path, ref_acc_path, query_aln_op):
		"""Align one query group. Returns True when nextalign succeeded.

		The return value matters: when nextalign dies before writing an
		errors.csv there is nothing for :meth:`update_gb_matrix` to read, so
		without a signal here the whole group passed into the tree build as
		though it had aligned cleanly.
		"""
		accession = self.path_to_basename(query_acc_path)
		# Use baseline profile for intra-group queries as they are tightly clustered
		base_profile = self.relaxation_profiles[0]
		command = self._nextalign_command(
			base_profile, query_acc_path, ref_acc_path,
			join(query_aln_op, f'{accession}'), accession,
		)

		return_code = self._run_nextalign(command, f"Executing Query Alignment for {accession}")
		if return_code != 0:
			print(f"Query alignment error: {accession} exited with code {return_code}")
			return False
		return True

	@staticmethod
	def _read_nextalign_errors(error_file_path, reference_accession):
		"""``{accession: [error, ...]}`` from a nextalign errors.csv.

		Falls back to reading the first two columns positionally when the file
		does not carry the ``seqName``/``errors`` header. That fallback used to
		be unreachable: the seek that rewinds to the top sat *outside* the
		``with`` block, so it ran against a closed file and every unexpected or
		absent header raised ``ValueError`` and took the run down.
		"""
		failed = {}
		with open(error_file_path, newline='', encoding='utf-8') as handle:
			reader = csv.DictReader(handle)
			fieldnames = reader.fieldnames or []
			if {'seqName', 'errors'}.issubset(set(fieldnames)):
				for row in reader:
					acc = str(row.get('seqName') or '').strip()
					error = str(row.get('errors') or '').strip()
					if not acc or not error or acc == reference_accession:
						continue
					failed.setdefault(acc, []).append(error)
				return failed

			handle.seek(0)
			for row in csv.reader(handle):
				if len(row) < 2:
					continue
				acc = str(row[0]).strip()
				error = str(row[1]).strip()
				if not acc or not error or acc == reference_accession:
					continue
				if acc == 'seqName':
					# A header we did not recognise, not a sequence.
					continue
				failed.setdefault(acc, []).append(error)
		return failed

	@staticmethod
	def _sanitise_criteria(text):
		"""Collapse whitespace so an exclusion criterion cannot break the TSV.

		The matrix is read downstream both with ``csv.DictReader`` and with a
		bare ``line.strip().split('\\t')`` (ValidateMatrix.py), so quoting an
		embedded tab satisfies the first reader while silently giving the second
		the wrong number of fields.
		"""
		return re.sub(r'\s+', ' ', str(text)).strip()

	@staticmethod
	def _split_criteria(value):
		return [part.strip() for part in str(value or '').split(';') if part.strip()]

	def _apply_exclusions(self, gB_matrix_file, failed_accessions):
		"""Read the matrix and return ``(fieldnames, rows)`` with exclusions applied.

		Nothing is written here, deliberately. The matrix used to be truncated
		by ``open(..., 'w')`` before the rows were serialised, so a row with more
		fields than the header raised part-way through ``writerows`` and left the
		pipeline's central table holding only the rows written before the raise.
		"""
		with open(gB_matrix_file, newline='', encoding='utf-8') as csvfile:
			reader = csv.DictReader(csvfile, delimiter='\t', restval='')
			if not reader.fieldnames:
				raise ValueError(f"gB matrix has no header row: {gB_matrix_file}")

			# A copy. Appending to reader.fieldnames while the reader is still
			# parsing materialises the new keys on every row with restval, so
			# row.get('exclusion_status', '0') returned None - the key exists -
			# and every row was written with a blank status instead of '0'.
			fieldnames = list(reader.fieldnames)
			for column in ('exclusion_status', 'exclusion_criteria'):
				if column not in fieldnames:
					fieldnames.append(column)

			updated_rows = []
			for line_number, row in enumerate(reader, start=2):
				overflow = row.pop(None, None)
				if overflow:
					raise ValueError(
						f"{gB_matrix_file} line {line_number} has more fields than the "
						f"{len(reader.fieldnames)}-column header; trailing values: {overflow!r}. "
						f"The matrix has not been modified."
					)

				gi = row.get('gi_number')
				acc_type = str(row.get('accession_type') or '').strip().lower()
				status = row.get('exclusion_status') or '0'
				criteria = self._split_criteria(row.get('exclusion_criteria'))

				if gi in failed_accessions and acc_type not in {'reference', 'master', 'exclusion_list'}:
					status = '1'
					for error in failed_accessions[gi]:
						criterion = self._sanitise_criteria(error)
						# Appended only if new, so a resumed run does not grow
						# the field by one copy of every criterion per pass.
						if criterion and criterion not in criteria:
							criteria.append(criterion)

				row['exclusion_status'] = status
				row['exclusion_criteria'] = '; '.join(criteria)
				updated_rows.append(row)

		return fieldnames, updated_rows

	@staticmethod
	def _write_matrix_atomically(gB_matrix_file, fieldnames, rows):
		"""Serialise to a sibling temp file, then ``os.replace`` it into place.

		``os.replace`` is atomic within a filesystem, so an interrupted or
		failing write leaves the previous matrix intact rather than truncated.
		"""
		tmp_path = gB_matrix_file + '.tmp'
		try:
			with open(tmp_path, 'w', newline='', encoding='utf-8') as csvfile:
				writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
				writer.writeheader()
				writer.writerows(rows)
			os.replace(tmp_path, gB_matrix_file)
		finally:
			if os.path.exists(tmp_path):
				os.remove(tmp_path)

	def update_gb_matrix(self, alignment_dir: list, gB_matrix_file, extra_failures=None):
		"""Mark sequences nextalign could not align as excluded.

		``extra_failures`` carries ``{accession: [reason]}`` for groups that
		never produced an errors.csv because nextalign itself failed to run.
		"""
		failed_accessions = {}

		for each_aln_type in alignment_dir:
			if os.path.basename(os.path.normpath(each_aln_type)) != "query_aln":
				continue
			if not os.path.isdir(each_aln_type):
				continue
			for each_aln in sorted(os.listdir(each_aln_type)):
				error_file_path = join(each_aln_type, each_aln, each_aln + ".errors.csv")
				if os.path.exists(error_file_path):
					for acc, errors in self._read_nextalign_errors(error_file_path, each_aln).items():
						failed_accessions.setdefault(acc, []).extend(errors)

		for acc, errors in (extra_failures or {}).items():
			failed_accessions.setdefault(acc, []).extend(errors)

		fieldnames, updated_rows = self._apply_exclusions(gB_matrix_file, failed_accessions)
		self._write_matrix_atomically(gB_matrix_file, fieldnames, updated_rows)

	def _align_query_groups(self, query_aln_output_dir):
		"""Align every query group against its like-named reference file.

		Returns ``{accession: [reason]}`` for groups that could not be aligned
		at all, so those sequences can be excluded rather than passing silently
		into the tree build.
		"""
		unaligned = {}
		for each_query_file in self._listdir_files(self.query_dir):
			query_path = join(self.query_dir, each_query_file)
			ref_path = join(self.ref_dir, each_query_file)

			if not os.path.isfile(ref_path):
				reason = f"No reference sequence file for query group {each_query_file}"
				print(f"WARNING: {reason}; group not aligned.")
				self._record_group_failure(unaligned, query_path, reason)
				continue

			if not self.nextalign_query(query_path, ref_path, query_aln_output_dir):
				self._record_group_failure(
					unaligned, query_path,
					f"nextalign failed to run for query group {each_query_file}",
				)

		return unaligned

	@classmethod
	def _record_group_failure(cls, unaligned, query_path, reason):
		for accession in cls._accessions_in_fasta(query_path):
			unaligned.setdefault(accession, []).append(reason)

	def process(self):
		query_aln_output_dir = join(self.tmp_dir, self.nextalign_dir, "query_aln")
		ref_aln_output_dir = join(self.tmp_dir, self.nextalign_dir, "reference_aln")

		os.makedirs(query_aln_output_dir, exist_ok=True)
		os.makedirs(ref_aln_output_dir, exist_ok=True)

		unaligned = self._align_query_groups(query_aln_output_dir)

		if self.reference_alignment:
			self.update_gb_matrix([query_aln_output_dir], self.gb_matrix, unaligned)
			return

		for each_ref in self._listdir_files(self.master_seq_dir):
			self.nextalign_master(self.ref_fa_file, join(self.master_seq_dir, each_ref), ref_aln_output_dir)

		masters = self.get_master_list()
		for master in masters:
			input_seq = join(ref_aln_output_dir, master, master + ".aligned.fasta")
			if os.path.exists(input_seq):
				output_seq = input_seq
				unique_seqs = RemoveRedundantSequence(input_seq, output_seq)
				unique_seqs.remove_redundant_fasta()
			else:
				# Used to be an unconditional `if os.path.exists`, so a master
				# whose alignment was written under a different spelling of the
				# accession skipped redundancy removal without a word.
				print(f"WARNING: no master alignment at {input_seq}; redundancy removal skipped for {master}.")

		self.update_gb_matrix([query_aln_output_dir, ref_aln_output_dir], self.gb_matrix, unaligned)

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
	parser.add_argument('--max_threads', type=int, help='Maximum number of threads to use', default=10)

	args = parser.parse_args()

	processor = NextalignAlignment(args.gB_matrix, args.query_dir, args.ref_dir, args.ref_fa_file, args.master_seq_dir, args.tmp_dir, args.master_ref, args.nextalign_dir, args.ref_alignment_file, args.update_db, args.max_threads)
	processor.process()
