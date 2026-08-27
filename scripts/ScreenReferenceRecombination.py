"""Screen the curated reference set for recombination with UShER/RIPPLES.

Why only the reference set
--------------------------
RIPPLES detects recombination on a mutation-annotated tree (MAT) by trying to
explain a branch's mutations with two parents instead of one. Its cost per
branch scales with the *number of breakpoint pairs* it must try, which grows
quadratically with the mutations on that branch. For SARS-CoV-2, where branches
carry a handful of mutations, this is cheap. For a diverse virus like HCV, where
inter-genotype branches carry thousands, it is not: on a 33k-sample HCV MAT the
first branch alone required 19,341,090 breakpoint pairs, and there were 1,790
branches to test.

The reference set is the tractable and the *important* case. It is small (238
sequences for HCV), it changes rarely, and it is the foundation every genotype
call rests on - a recombinant hiding among the references mis-genotypes every
query that matches it. Screening queries at scale is a different problem needing
a different method; see `info_help/reference_recombination_screening.md`.

Nothing here is virus-specific. The reference set, the master accession and the
alignments all come from the database.
"""

import argparse
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import time


# ---------------------------------------------------------------------------
# Output formats, taken from the RIPPLES binary itself rather than guessed
# ---------------------------------------------------------------------------

#: Header RIPPLES writes to recombination.tsv.
RECOMBINATION_COLUMNS = (
	"recomb_node_id",
	"breakpoint-1_interval",
	"breakpoint-2_interval",
	"donor_node_id",
	"donor_is_sibling",
	"donor_parsimony",
	"acceptor_node_id",
	"acceptor_is_sibling",
	"acceptor_parsimony",
	"original_parsimony",
	"min_starting_parsimony",
	"recomb_parsimony",
)

#: Header RIPPLES writes to descendants.tsv.
DESCENDANTS_COLUMNS = ("node_id", "descendants")

#: Column names as stored in the database. Hyphens are not usable unquoted in
#: SQL, so the two interval columns are renamed; everything else is verbatim.
DB_COLUMN_FOR = {
	"breakpoint-1_interval": "breakpoint_1_interval",
	"breakpoint-2_interval": "breakpoint_2_interval",
}

#: Hard ceiling on threads for every external tool this script launches.
#: This runs on a shared machine - RIPPLES and UShER both default to taking
#: every core they can see (224 here), which starves other users. Requests above
#: this are clamped, not honoured, and the clamp is announced.
MAX_THREADS = 16

RECOMBINATION_TABLE = "reference_recombination"
RECOMBINATION_STATUS_COLUMN = "recombination_status"

#: Values written to meta_data.recombination_status.
STATUS_RECOMBINANT = "recombinant_reference"
STATUS_SCREENED = "screened_no_evidence"


# ---------------------------------------------------------------------------
# Parsing - pure functions, no subprocesses, so they are directly testable
# ---------------------------------------------------------------------------

def parse_interval(text):
	"""Parse a RIPPLES breakpoint interval into (start, end) integers.

	RIPPLES writes intervals as ``1234-5678``. A single coordinate, an empty
	field, or anything unparseable yields (None, None) rather than raising:
	a malformed interval must not lose the rest of the event.
	"""
	if text is None:
		return (None, None)
	text = str(text).strip()
	if not text:
		return (None, None)
	match = re.match(r'^\s*(\d+)\s*-\s*(\d+)\s*$', text)
	if match:
		return (int(match.group(1)), int(match.group(2)))
	if re.match(r'^\d+$', text):
		value = int(text)
		return (value, value)
	return (None, None)


def _read_tsv(path, expected_columns, label):
	"""Read a RIPPLES TSV, tolerating the '#'-prefixed header or none at all."""
	if not path or not os.path.exists(path):
		return []
	rows = []
	with open(path, "r", encoding="utf-8", errors="replace") as handle:
		header = None
		for line in handle:
			line = line.rstrip("\n\r")
			if not line.strip():
				continue
			fields = line.split("\t")
			if header is None:
				if fields[0].lstrip("#").strip() == expected_columns[0]:
					header = [f.lstrip("#").strip() for f in fields]
					continue
				# No header written (RIPPLES omits it when there are no events)
				header = list(expected_columns)
			if len(fields) < len(expected_columns):
				fields = fields + [""] * (len(expected_columns) - len(fields))
			rows.append(dict(zip(header, fields)))
	return rows


def parse_recombination_tsv(path):
	"""Parse recombination.tsv into a list of dicts with parsed intervals."""
	events = []
	for row in _read_tsv(path, RECOMBINATION_COLUMNS, "recombination"):
		event = {DB_COLUMN_FOR.get(k, k): v for k, v in row.items()}
		start1, end1 = parse_interval(row.get("breakpoint-1_interval"))
		start2, end2 = parse_interval(row.get("breakpoint-2_interval"))
		event["breakpoint_1_start"] = start1
		event["breakpoint_1_end"] = end1
		event["breakpoint_2_start"] = start2
		event["breakpoint_2_end"] = end2
		events.append(event)
	return events


def parse_descendants_tsv(path):
	"""Parse descendants.tsv into {node_id: [sample, ...]}.

	RIPPLES reports events against internal nodes. The descendants file is the
	only thing that maps a node back to the accessions it covers, so without it
	an event cannot be attributed to any sequence.
	"""
	mapping = {}
	for row in _read_tsv(path, DESCENDANTS_COLUMNS, "descendants"):
		node = (row.get("node_id") or "").strip()
		if not node:
			continue
		raw = (row.get("descendants") or "").strip()
		samples = [s.strip() for s in re.split(r'[,\s]+', raw) if s.strip()]
		mapping.setdefault(node, [])
		mapping[node].extend(samples)
	return mapping


def attribute_events_to_accessions(events, descendants):
	"""Expand node-level events into one row per (accession, event).

	An event on an internal node applies to every accession beneath it. Events
	whose node has no descendants listed are kept with a NULL accession rather
	than dropped - losing a detected event because the descendants file was
	incomplete would be worse than reporting it unattributed.
	"""
	rows = []
	for event in events:
		node = (event.get("recomb_node_id") or "").strip()
		samples = descendants.get(node) or []
		if not samples:
			row = dict(event)
			row["primary_accession"] = None
			rows.append(row)
			continue
		for sample in samples:
			row = dict(event)
			row["primary_accession"] = sample
			rows.append(row)
	return rows


def chunk_bounds(n_branches, n_chunks):
	"""Split branch indices into contiguous [start, end) chunks for -S/-E.

	RIPPLES parallelises poorly *within* a branch - measured on HCV, going from
	8 to 48 threads did not shorten the first branch - so the way to use a big
	machine is to run several RIPPLES processes over disjoint branch ranges.
	"""
	if n_branches <= 0 or n_chunks <= 0:
		return []
	n_chunks = min(n_chunks, n_branches)
	base, extra = divmod(n_branches, n_chunks)
	bounds = []
	start = 0
	for i in range(n_chunks):
		size = base + (1 if i < extra else 0)
		bounds.append((start, start + size))
		start += size
	return bounds


# ---------------------------------------------------------------------------
# Command construction - pure, so the CLI contract is testable without binaries
# ---------------------------------------------------------------------------

def clamp_threads(requested, maximum=MAX_THREADS):
	"""Clamp a thread request into [1, maximum].

	Returns the clamped value; callers announce it when it differs from the
	request so a silently-reduced run is never mistaken for the one asked for.
	"""
	try:
		value = int(requested)
	except (TypeError, ValueError):
		return 1
	return max(1, min(value, int(maximum)))


def build_fatovcf_command(fasta_path, vcf_path, binary="faToVcf"):
	return [binary, fasta_path, vcf_path]


def build_usher_command(tree_path, vcf_path, output_mat, outdir, threads, binary="usher"):
	return [
		binary,
		"-t", tree_path,
		"-v", vcf_path,
		"-o", output_mat,
		"-d", outdir,
		"-T", str(clamp_threads(threads)),
	]


def build_ripples_command(mat_path, outdir, threads, branch_length=3,
                          min_range=1000, max_range=10000000,
                          parsimony_improvement=3, num_descendants=3,
                          start_index=None, end_index=None,
                          samples_file=None, binary="ripples"):
	"""Assemble a RIPPLES invocation.

	Defaults differ from RIPPLES' own in one place: ``num_descendants`` is 3
	rather than 10. A curated reference set has only a handful of sequences per
	subtype, so a threshold of 10 descendants would skip most of the tree.
	"""
	command = [
		binary,
		"-i", mat_path,
		"-d", outdir,
		"-T", str(clamp_threads(threads)),
		"-l", str(branch_length),
		"-r", str(min_range),
		"-R", str(max_range),
		"-p", str(parsimony_improvement),
		"-n", str(num_descendants),
	]
	if start_index is not None:
		command += ["-S", str(start_index)]
	if end_index is not None:
		command += ["-E", str(end_index)]
	if samples_file:
		command += ["-s", samples_file]
	return command


# ---------------------------------------------------------------------------
# Reading the reference set out of the database
# ---------------------------------------------------------------------------

REFERENCE_TYPES = ("reference", "master")


def reference_records_from_db(conn):
	"""Return (master_accession, {accession: alignment}) for the reference set.

	Raises ValueError if there is no master or no reference alignments - both
	mean the database cannot support this screen, and failing loudly beats
	writing an empty result that looks like 'no recombination found'.
	"""
	placeholders = ",".join("?" * len(REFERENCE_TYPES))
	master_row = conn.execute(
		"SELECT primary_accession FROM meta_data "
		"WHERE lower(COALESCE(accession_type,'')) = 'master' LIMIT 1"
	).fetchone()
	if not master_row:
		raise ValueError("no accession_type='master' row in meta_data")
	master = master_row[0]

	rows = conn.execute(
		"SELECT sa.primary_accession, sa.alignment FROM sequence_alignment AS sa "
		"JOIN meta_data AS md ON md.primary_accession = sa.primary_accession "
		f"WHERE lower(COALESCE(md.accession_type,'')) IN ({placeholders})",
		REFERENCE_TYPES,
	).fetchall()
	sequences = {a: s for a, s in rows if s}
	if not sequences:
		raise ValueError("no alignments found for reference/master accessions")
	if master not in sequences:
		raise ValueError(f"master {master} has no row in sequence_alignment")

	widths = {len(s) for s in sequences.values()}
	if len(widths) != 1:
		raise ValueError(
			f"reference alignments are ragged: widths {sorted(widths)}. "
			"faToVcf requires a rectangular alignment."
		)
	return master, sequences


def write_reference_fasta(master, sequences, path):
	"""Write the reference alignment with the master FIRST.

	faToVcf treats the first record as the VCF reference, so the order is not
	cosmetic - putting anything else first would express every variant against
	the wrong coordinate frame.
	"""
	with open(path, "w", encoding="utf-8") as handle:
		handle.write(f">{master}\n{sequences[master]}\n")
		for accession in sorted(sequences):
			if accession == master:
				continue
			handle.write(f">{accession}\n{sequences[accession]}\n")
	return path


# ---------------------------------------------------------------------------
# Newick handling, for deriving a starting tree from the database
# ---------------------------------------------------------------------------

_NEWICK_TOKEN = re.compile(r'[(),;]|[^(),;]+')


class _Node:
	__slots__ = ("name", "children")

	def __init__(self, name=None):
		self.name = name
		self.children = []


def _parse_newick(text):
	root = _Node()
	current = root
	stack = []
	for token in _NEWICK_TOKEN.findall(text.strip()):
		if token == '(':
			node = _Node()
			current.children.append(node)
			stack.append(current)
			current = node
		elif token == ',':
			parent = stack[-1]
			node = _Node()
			parent.children.append(node)
			current = node
		elif token == ')':
			current = stack.pop()
		elif token == ';':
			break
		else:
			label = token.split(':')[0].strip().strip("'\"")
			if label and not current.children:
				current.name = label
	return root


def _restrict(node, keep):
	if not node.children:
		return _Node(node.name) if node.name in keep else None
	kept = [k for k in (_restrict(c, keep) for c in node.children) if k is not None]
	if not kept:
		return None
	if len(kept) == 1:
		return kept[0]
	parent = _Node()
	parent.children = kept
	return parent


def _write_newick(node):
	if not node.children:
		return node.name or ""
	return "(" + ",".join(_write_newick(c) for c in node.children) + ")"


def prune_newick_to(newick, keep):
	"""Restrict a newick string to `keep`, suppressing unary nodes.

	Branch lengths and support values are discarded: UShER re-derives branch
	lengths from the VCF, and only the topology is used as a starting point.
	"""
	pruned = _restrict(_parse_newick(newick), set(keep))
	if pruned is None:
		raise ValueError("none of the requested tips are present in the tree")
	return _write_newick(pruned) + ";"


def starting_tree_from_db(conn, keep):
	"""Pick the best available tree in the DB and prune it to the reference set.

	Selection is by **coverage first**, then by source. A tree that contains
	more of the reference set is always preferred, because a missing reference
	would simply be absent from the screen. Among trees with equal coverage an
	IQ-TREE wins over an UShER tree, since UShER branch lengths are parsimony
	placements - though only topology is used here either way.

	In practice the UShER tree usually wins on HCV: it holds all 238 references
	where the IQ-TREE, built from cluster representatives, holds 219. Pass an
	explicit --tree to override this entirely.
	"""
	rows = conn.execute("SELECT name, source, newick FROM trees").fetchall()
	if not rows:
		raise ValueError("no rows in the trees table")

	def rank(row):
		source = (row[1] or "").lower()
		return 0 if "iqtree" in source else 1

	best = None
	for row in sorted(rows, key=rank):
		try:
			pruned = prune_newick_to(row[2], keep)
		except ValueError:
			continue
		covered = set(re.findall(r'[(,]([^(),;]+)', pruned)) & set(keep)
		if best is None or len(covered) > best[0]:
			best = (len(covered), row[0], pruned)
		if len(covered) == len(keep):
			break
	if best is None:
		raise ValueError("no tree in the database covers the reference set")
	return best[1], best[2]


# ---------------------------------------------------------------------------
# Database output
# ---------------------------------------------------------------------------

def _quote_identifier(name):
	return '"' + str(name).replace('"', '""') + '"'


def ensure_recombination_tables(conn):
	"""Create the results table and the meta_data status column if absent."""
	stored = [DB_COLUMN_FOR.get(c, c) for c in RECOMBINATION_COLUMNS]
	extra = [
		"primary_accession", "breakpoint_1_start", "breakpoint_1_end",
		"breakpoint_2_start", "breakpoint_2_end", "detected_by", "detected_at",
	]
	columns = ", ".join(_quote_identifier(c) + " TEXT" for c in extra + stored)
	conn.execute(f"CREATE TABLE IF NOT EXISTS {RECOMBINATION_TABLE} ({columns})")
	conn.execute(
		f"CREATE UNIQUE INDEX IF NOT EXISTS idx_{RECOMBINATION_TABLE}_event ON "
		f"{RECOMBINATION_TABLE} (primary_accession, recomb_node_id, "
		f"breakpoint_1_interval, breakpoint_2_interval)"
	)
	existing = {row[1].lower() for row in conn.execute("PRAGMA table_info(meta_data)")}
	if RECOMBINATION_STATUS_COLUMN.lower() not in existing:
		conn.execute(
			f"ALTER TABLE meta_data ADD COLUMN "
			f"{_quote_identifier(RECOMBINATION_STATUS_COLUMN)} TEXT"
		)


def store_results(conn, rows, screened_accessions, detected_by="ripples", detected_at=None):
	"""Write events and set meta_data.recombination_status.

	Idempotent: re-running the same screen replaces the same rows rather than
	accumulating duplicates, which is what the unique index enforces.
	"""
	ensure_recombination_tables(conn)
	detected_at = detected_at or time.strftime("%Y-%m-%d %H:%M:%S")
	stored = [DB_COLUMN_FOR.get(c, c) for c in RECOMBINATION_COLUMNS]
	extra = [
		"primary_accession", "breakpoint_1_start", "breakpoint_1_end",
		"breakpoint_2_start", "breakpoint_2_end", "detected_by", "detected_at",
	]
	columns = extra + stored
	placeholders = ",".join("?" * len(columns))
	payload = []
	for row in rows:
		record = dict(row)
		record["detected_by"] = detected_by
		record["detected_at"] = detected_at
		payload.append([
			None if record.get(c) is None else str(record.get(c)) for c in columns
		])
	conn.executemany(
		f"INSERT OR REPLACE INTO {RECOMBINATION_TABLE} "
		f"({','.join(_quote_identifier(c) for c in columns)}) VALUES ({placeholders})",
		payload,
	)

	flagged = sorted({r["primary_accession"] for r in rows if r.get("primary_accession")})
	status_column = _quote_identifier(RECOMBINATION_STATUS_COLUMN)
	if screened_accessions:
		marks = ",".join("?" * len(screened_accessions))
		conn.execute(
			f"UPDATE meta_data SET {status_column} = ? "
			f"WHERE primary_accession IN ({marks})",
			[STATUS_SCREENED] + list(screened_accessions),
		)
	if flagged:
		marks = ",".join("?" * len(flagged))
		conn.execute(
			f"UPDATE meta_data SET {status_column} = ? "
			f"WHERE primary_accession IN ({marks})",
			[STATUS_RECOMBINANT] + flagged,
		)
	conn.commit()
	return len(payload), len(flagged)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def _require_binary(name):
	path = shutil.which(name)
	if not path:
		raise FileNotFoundError(
			f"{name} not found on PATH. It ships with the UShER suite "
			f"(conda install -c bioconda usher)."
		)
	return path


def count_long_branches(log_text):
	"""Pull 'Found N long branches' out of a RIPPLES log.

	This is the only number that tells you how much work the run implies, and
	it is worth surfacing before committing hours to it.
	"""
	match = re.search(r'Found\s+(\d+)\s+long branches', log_text or "")
	return int(match.group(1)) if match else None


class ReferenceRecombinationScreen:
	def __init__(self, db=None, alignment=None, tree=None, outdir=".", threads=8,
	             branch_length=3, min_range=1000, max_range=10000000,
	             parsimony_improvement=3, num_descendants=3,
	             start_index=None, end_index=None, chunks=1,
	             write_db=False, skip_build=False, timeout=None):
		self.db = db
		self.alignment = alignment
		self.tree = tree
		self.outdir = outdir
		requested = threads
		self.threads = clamp_threads(threads)
		if self.threads != requested:
			print(f"[threads] clamped {requested} -> {self.threads} "
			      f"(MAX_THREADS={MAX_THREADS}; this is a shared machine)")
		self.branch_length = branch_length
		self.min_range = min_range
		self.max_range = max_range
		self.parsimony_improvement = parsimony_improvement
		self.num_descendants = num_descendants
		self.start_index = start_index
		self.end_index = end_index
		self.chunks = max(1, int(chunks or 1))
		self.write_db = write_db
		self.skip_build = skip_build
		self.timeout = timeout
		self.accessions = []

	# -- inputs ------------------------------------------------------------
	def prepare_inputs(self):
		os.makedirs(self.outdir, exist_ok=True)
		fasta = os.path.join(self.outdir, "reference_set.fasta")
		tree_path = os.path.join(self.outdir, "reference_start.nwk")

		if self.db:
			conn = sqlite3.connect(f"file:{self.db}?mode=ro", uri=True)
			try:
				master, sequences = reference_records_from_db(conn)
				write_reference_fasta(master, sequences, fasta)
				self.accessions = sorted(sequences)
				if self.tree:
					shutil.copyfile(self.tree, tree_path)
				else:
					name, newick = starting_tree_from_db(conn, set(sequences))
					print(f"[tree] derived starting tree from '{name}' in the database")
					with open(tree_path, "w", encoding="utf-8") as handle:
						handle.write(newick + "\n")
			finally:
				conn.close()
		else:
			if not (self.alignment and self.tree):
				raise ValueError("without --db you must supply both --alignment and --tree")
			shutil.copyfile(self.alignment, fasta)
			shutil.copyfile(self.tree, tree_path)
			with open(fasta, "r", encoding="utf-8") as handle:
				self.accessions = [
					line[1:].strip().split()[0] for line in handle if line.startswith(">")
				]
		print(f"[input] {len(self.accessions)} reference sequences -> {fasta}")
		return fasta, tree_path

	# -- MAT ---------------------------------------------------------------
	def build_mat(self, fasta, tree_path):
		mat = os.path.join(self.outdir, "reference.pb")
		if self.skip_build and os.path.exists(mat):
			print(f"[mat] reusing {mat}")
			return mat
		vcf = os.path.join(self.outdir, "reference_set.vcf")
		_require_binary("faToVcf")
		_require_binary("usher")
		subprocess.run(build_fatovcf_command(fasta, vcf), check=True)
		subprocess.run(
			build_usher_command(tree_path, vcf, "reference.pb", self.outdir, self.threads),
			check=True,
		)
		print(f"[mat] built {mat}")
		return mat

	# -- RIPPLES -----------------------------------------------------------
	def run_ripples(self, mat):
		_require_binary("ripples")
		results = os.path.join(self.outdir, "ripples")
		os.makedirs(results, exist_ok=True)
		command = build_ripples_command(
			mat, results, self.threads,
			branch_length=self.branch_length, min_range=self.min_range,
			max_range=self.max_range, parsimony_improvement=self.parsimony_improvement,
			num_descendants=self.num_descendants,
			start_index=self.start_index, end_index=self.end_index,
		)
		log_path = os.path.join(self.outdir, "ripples.log")
		print("[ripples] " + " ".join(command))
		started = time.time()
		with open(log_path, "w", encoding="utf-8") as log:
			try:
				subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT,
				               timeout=self.timeout)
			except subprocess.TimeoutExpired:
				print(f"[ripples] TIMED OUT after {self.timeout}s - partial results kept")
		elapsed = time.time() - started
		with open(log_path, "r", encoding="utf-8", errors="replace") as log:
			text = log.read()
		branches = count_long_branches(text)
		if branches is not None:
			print(f"[ripples] {branches} long branches to test; ran {elapsed:.0f}s")
			if self.chunks == 1 and branches > 20:
				per_chunk = chunk_bounds(branches, min(16, branches))
				print(f"[ripples] this is slow serially. To split it across "
				      f"{len(per_chunk)} concurrent runs use --start_index/--end_index, e.g.")
				for lo, hi in per_chunk[:3]:
					print(f"             --start_index {lo} --end_index {hi}")
				print("             ... (see the docs for a ready-made loop)")
		return results

	# -- results -----------------------------------------------------------
	def collect(self, results_dir):
		events = parse_recombination_tsv(os.path.join(results_dir, "recombination.tsv"))
		descendants = parse_descendants_tsv(os.path.join(results_dir, "descendants.tsv"))
		rows = attribute_events_to_accessions(events, descendants)
		print(f"[result] {len(events)} recombination events covering "
		      f"{len({r['primary_accession'] for r in rows if r.get('primary_accession')})} accessions")
		for row in rows[:20]:
			print(f"   {row.get('primary_accession') or '(unattributed)':12s} "
			      f"node={row.get('recomb_node_id')} "
			      f"bp1={row.get('breakpoint_1_interval')} bp2={row.get('breakpoint_2_interval')} "
			      f"donor={row.get('donor_node_id')} acceptor={row.get('acceptor_node_id')}")
		if self.write_db:
			if not self.db:
				raise ValueError("--write_db needs --db")
			conn = sqlite3.connect(self.db)
			try:
				written, flagged = store_results(conn, rows, self.accessions)
				print(f"[db] wrote {written} rows to {RECOMBINATION_TABLE}; "
				      f"{flagged} accessions marked {STATUS_RECOMBINANT}")
			finally:
				conn.close()
		return rows

	def run(self):
		fasta, tree_path = self.prepare_inputs()
		mat = self.build_mat(fasta, tree_path)
		results = self.run_ripples(mat)
		return self.collect(results)


def parse_args(argv=None):
	parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
	parser.add_argument("--db", help="SQLite database holding the reference set")
	parser.add_argument("--alignment", help="reference alignment FASTA (instead of --db)")
	parser.add_argument("--tree", help="starting tree; derived from --db when omitted")
	parser.add_argument("--outdir", default="reference_recombination")
	parser.add_argument("--threads", type=int, default=8,
	                    help=f"threads for usher/ripples; clamped to {MAX_THREADS}")
	parser.add_argument("--branch_length", type=int, default=3,
	                    help="RIPPLES -l: minimum mutations on a branch to test it")
	parser.add_argument("--min_range", type=int, default=1000, help="RIPPLES -r")
	parser.add_argument("--max_range", type=int, default=10000000, help="RIPPLES -R")
	parser.add_argument("--parsimony_improvement", type=int, default=3, help="RIPPLES -p")
	parser.add_argument("--num_descendants", type=int, default=3, help="RIPPLES -n")
	parser.add_argument("--start_index", type=int, default=None, help="RIPPLES -S")
	parser.add_argument("--end_index", type=int, default=None, help="RIPPLES -E")
	parser.add_argument("--chunks", type=int, default=1,
	                    help="advisory only; suppresses the split-it-up hint when >1")
	parser.add_argument("--write_db", action="store_true")
	parser.add_argument("--skip_build", action="store_true",
	                    help="reuse an existing reference.pb in --outdir")
	parser.add_argument("--timeout", type=float, default=None,
	                    help="seconds before RIPPLES is stopped; partial results are kept")
	return parser.parse_args(argv)


def main(argv=None):
	args = parse_args(argv)
	screen = ReferenceRecombinationScreen(
		db=args.db, alignment=args.alignment, tree=args.tree, outdir=args.outdir,
		threads=args.threads, branch_length=args.branch_length,
		min_range=args.min_range, max_range=args.max_range,
		parsimony_improvement=args.parsimony_improvement,
		num_descendants=args.num_descendants, start_index=args.start_index,
		end_index=args.end_index, chunks=args.chunks, write_db=args.write_db,
		skip_build=args.skip_build, timeout=args.timeout,
	)
	screen.run()
	return 0


if __name__ == "__main__":
	sys.exit(main())
