#!/usr/bin/env python3
"""Hunt for recombination across a V-gTK database, then report on it.

This is the driver around :mod:`ScreenReferenceRecombination`, which does the
actual UShER/RIPPLES work. Two things are added here, both of which the screen
deliberately leaves to its caller:

1. **The chunked run is executed, not merely suggested.** RIPPLES does not
   parallelise *within* a branch - measured on HCV, 8 threads and 48 threads
   placed the first branch in the same time - so the only way to use a big
   machine is to run several processes over disjoint branch ranges. The screen
   prints the ranges and stops. This module builds the MAT once, discovers how
   many long branches there are, splits them, runs the chunks concurrently
   under a fixed thread budget, and merges the results back into one answer.

2. **A report.** The screen leaves its findings in a database table. What a
   reviewer actually needs is which references are implicated, where the
   breakpoints fall, how far the two-parent explanation beats the one-parent
   one, and - just as important - how much of the reference set was screened
   and came back clean. A screened-clean reference and a never-screened one are
   different states and the report keeps them apart.

Discovering the branch count needs care. RIPPLES only announces
``Found N long branches`` once it has loaded the MAT and started work, so the
count is probed by launching RIPPLES, waiting for that line, and stopping it
again. The probe costs seconds against a run that costs hours.

The thread ceiling of the screen (16) is inherited and applies to the *whole*
fan-out, not per chunk: this runs on a shared machine, and 8 chunks that each
grabbed 16 threads would take 128 cores.

Usage
-----
    # full screen, 8 chunks x 2 threads, results written back to the database
    python scripts/recombination_hunter_db.py --db path/to.db --chunks 8 --write_db

    # regenerate the report from a database that has already been screened
    python scripts/recombination_hunter_db.py --db path/to.db --report_only
"""

import argparse
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

# The screen lives beside this file. Importing it by path rather than relying on
# the caller's cwd means the script works from anywhere, which matters because
# Nextflow invokes scripts by absolute path from a task work directory.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ScreenReferenceRecombination as screen  # noqa: E402


MAX_THREADS = screen.MAX_THREADS
RECOMBINATION_TABLE = screen.RECOMBINATION_TABLE
STATUS_RECOMBINANT = screen.STATUS_RECOMBINANT
STATUS_SCREENED = screen.STATUS_SCREENED

#: 8 x 2 saturates the 16-thread ceiling. More chunks beat more threads per
#: chunk for the reason in the module docstring, so the default leans that way.
DEFAULT_CHUNKS = 8
DEFAULT_THREADS_PER_CHUNK = 2

#: How long to wait for RIPPLES to announce its branch count before giving up
#: on chunking and falling back to a single serial run.
PROBE_TIMEOUT = 1800

#: Which RIPPLES binary to drive. Both ship with the UShER suite and take an
#: identical flag surface, so this is a drop-in choice.
#:
#: The default is the fast one, because the exhaustive one cannot finish on the
#: data this tool exists to screen. RIPPLES was built for SARS-CoV-2, where a
#: branch carries a handful of mutations. On a reference set spanning all eight
#: HCV genotypes, 76% of the 9,644 bp genome is a variable site: the parsimony
#: tree carries 223,252 mutations over 475 branches, the median branch carries
#: 422 and only 7 branches carry fewer than 30. A `-l 3` "long branch" therefore
#: selects essentially the whole tree, and plain `ripples` responds by
#: enumerating ~2.6 million breakpoint pairs per branch at ~2.3 pairs/second -
#: about 13 days per branch, ~250 days for 156 branches. Measured here: six days
#: of wall clock finished 0 branches, reaching pair 1,582,696 of 2,591,226 on
#: the first. `ripples-fast` completed all 156 in 65 seconds.
#:
#: The trade is real and is recorded in the report and in `detected_by`: the
#: fast binary prunes the search rather than enumerating it, so it is not
#: guaranteed to agree with the exhaustive scan - and on data like this there is
#: no head-to-head to check it against, because the exhaustive scan never
#: finishes a single branch.
DEFAULT_RIPPLES_BINARY = "ripples-fast"

#: Metadata worth seeing beside a flagged reference. Every one is optional -
#: the columns are probed at runtime because schemas differ between builds.
CONTEXT_COLUMNS = (
	"nearest_reference_genotype",
	"nearest_reference_subtype",
	"country",
	"collection_date",
	"real_length",
)


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------

def plan_thread_budget(chunks, threads_per_chunk, ceiling=MAX_THREADS):
	"""Fit ``chunks x threads_per_chunk`` inside the ceiling.

	Threads per chunk are given up before chunks are, because concurrency across
	branches is what actually buys wall-clock here. Returns the adjusted pair;
	the caller announces any change, so a quietly-shrunk run is never mistaken
	for the one that was asked for.
	"""
	chunks = max(1, int(chunks or 1))
	threads_per_chunk = max(1, int(threads_per_chunk or 1))
	ceiling = max(1, int(ceiling))

	if chunks > ceiling:
		# More chunks than cores allowed: one thread each, and drop the excess.
		return ceiling, 1
	while chunks * threads_per_chunk > ceiling and threads_per_chunk > 1:
		threads_per_chunk -= 1
	return chunks, threads_per_chunk


def probe_long_branches(mat, outdir, threads, ripples_kwargs,
                        binary=DEFAULT_RIPPLES_BINARY,
                        timeout=PROBE_TIMEOUT, poll=2.0):
	"""Launch RIPPLES only to read ``Found N long branches``, then stop it.

	There is no way to ask RIPPLES for this number without starting it, and the
	number is what decides how the work should be split. The probe is stopped as
	soon as the line appears, so it costs seconds rather than the hours the real
	run will take.

	Returns the branch count, or ``None`` if RIPPLES never announced one - in
	which case the caller should fall back to a single unchunked run rather than
	guess a split.
	"""
	probe_dir = os.path.join(outdir, "probe")
	results_dir = os.path.join(probe_dir, "ripples")
	os.makedirs(results_dir, exist_ok=True)
	log_path = os.path.join(probe_dir, "ripples.probe.log")

	command = screen.build_ripples_command(mat, results_dir, threads,
	                                      binary=binary, **ripples_kwargs)
	print("[probe] " + " ".join(command))

	deadline = time.time() + timeout
	count = None
	with open(log_path, "w", encoding="utf-8") as log:
		process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
		try:
			while True:
				if os.path.exists(log_path):
					with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
						count = screen.count_long_branches(handle.read())
					if count is not None:
						break
				if process.poll() is not None:
					# Exited before we saw the line; the whole log is all we get.
					with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
						count = screen.count_long_branches(handle.read())
					break
				if time.time() > deadline:
					print(f"[probe] no branch count after {timeout:.0f}s - giving up on chunking")
					break
				time.sleep(poll)
		finally:
			if process.poll() is None:
				process.terminate()
				try:
					process.wait(timeout=30)
				except subprocess.TimeoutExpired:
					process.kill()
					process.wait()

	# The probe's partial output is not a result - it covers an arbitrary prefix
	# of the branches and would look like a complete screen if left behind.
	shutil.rmtree(results_dir, ignore_errors=True)
	return count


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def run_ripples(binary, mat, outdir, threads, ripples_kwargs,
                start_index=None, end_index=None, timeout=None):
	"""Run one RIPPLES process and return its results directory.

	The screen's own ``run_ripples`` hardcodes the ``ripples`` binary, so this
	reimplements the few lines around it rather than monkey-patching a module
	that other callers share.
	"""
	screen._require_binary(binary)
	results = os.path.join(outdir, "ripples")
	os.makedirs(results, exist_ok=True)
	command = screen.build_ripples_command(
		mat, results, threads,
		start_index=start_index, end_index=end_index, binary=binary,
		**ripples_kwargs,
	)
	log_path = os.path.join(outdir, "ripples.log")
	print("[ripples] " + " ".join(command))
	with open(log_path, "w", encoding="utf-8") as log:
		try:
			subprocess.run(command, check=True, stdout=log,
			               stderr=subprocess.STDOUT, timeout=timeout)
		except subprocess.TimeoutExpired:
			# Partial output is still output; discarding it would turn a slow
			# run into a silently clean one.
			print(f"[ripples] TIMED OUT after {timeout}s - partial results kept")
	return results


def run_chunk(binary, mat, chunk_dir, lo, hi, threads, ripples_kwargs, timeout=None):
	"""Run one RIPPLES process over the half-open branch range [lo, hi)."""
	os.makedirs(chunk_dir, exist_ok=True)
	started = time.time()
	run_ripples(binary, mat, chunk_dir, threads, ripples_kwargs,
	            start_index=lo, end_index=hi, timeout=timeout)
	elapsed = time.time() - started
	print(f"[chunk] branches [{lo},{hi}) finished in {elapsed / 60:.1f} min")
	return chunk_dir


def merge_chunks(chunk_dirs):
	"""Combine per-chunk RIPPLES output into one event list and node mapping.

	Branch ranges are disjoint so events should not repeat, but the merge
	de-duplicates anyway: a re-run of a single chunk into an existing directory
	would otherwise double-count, and a duplicated event would inflate the
	report's headline number.
	"""
	events, descendants, seen = [], {}, set()
	for chunk_dir in chunk_dirs:
		results = os.path.join(chunk_dir, "ripples")
		for event in screen.parse_recombination_tsv(os.path.join(results, "recombination.tsv")):
			key = (
				event.get("recomb_node_id"),
				event.get("breakpoint_1_interval"),
				event.get("breakpoint_2_interval"),
				event.get("donor_node_id"),
				event.get("acceptor_node_id"),
			)
			if key in seen:
				continue
			seen.add(key)
			events.append(event)
		for node, samples in screen.parse_descendants_tsv(
			os.path.join(results, "descendants.tsv")
		).items():
			bucket = descendants.setdefault(node, [])
			for sample in samples:
				if sample not in bucket:
					bucket.append(sample)
	return events, descendants


# ---------------------------------------------------------------------------
# Database reads for the report
# ---------------------------------------------------------------------------

def _existing_columns(conn, table, wanted):
	"""Intersect wanted columns with what the table actually has."""
	try:
		present = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
	except sqlite3.Error:
		return []
	return [c for c in wanted if c in present]


def gather_context(conn, accessions):
	"""Fetch the metadata shown beside each flagged accession.

	Missing columns are skipped rather than fatal: this runs against databases
	built by several profiles and the report is worth having either way.
	"""
	accessions = [a for a in accessions if a]
	if not accessions:
		return {}
	columns = _existing_columns(conn, "meta_data", CONTEXT_COLUMNS)
	if not columns:
		return {}
	marks = ",".join("?" * len(accessions))
	select = ",".join(["primary_accession"] + columns)
	rows = conn.execute(
		f"SELECT {select} FROM meta_data WHERE primary_accession IN ({marks})",
		accessions,
	)
	return {row[0]: dict(zip(columns, row[1:])) for row in rows}


def gather_status_counts(conn):
	"""Count reference-set rows by recombination_status.

	NULL is reported as 'not screened' and kept distinct from
	'screened_no_evidence'. Treating the two as one would turn an absence of
	evidence into evidence of absence, which is the whole point of the column.
	"""
	column = screen.RECOMBINATION_STATUS_COLUMN
	if column not in _existing_columns(conn, "meta_data", [column]):
		return None
	rows = conn.execute(
		f"SELECT COALESCE({column}, '(not screened)'), COUNT(*) "
		f"FROM meta_data "
		f"WHERE lower(COALESCE(accession_type,'')) IN ('reference','master') "
		f"GROUP BY 1 ORDER BY 2 DESC"
	)
	return list(rows)


def count_screened(conn):
	"""How many reference rows carry a recombination_status at all.

	Zero means the screen has never run against this database. That is a
	different statement from "the screen ran and found nothing", and the report
	must not blur them - reporting an unscreened set as clean is the one failure
	mode that would actively mislead.
	"""
	column = screen.RECOMBINATION_STATUS_COLUMN
	if column not in _existing_columns(conn, "meta_data", [column]):
		return 0
	row = conn.execute(
		f"SELECT COUNT(*) FROM meta_data "
		f"WHERE {column} IS NOT NULL "
		f"AND lower(COALESCE(accession_type,'')) IN ('reference','master')"
	).fetchone()
	return int(row[0]) if row else 0


SCREEN_RUNS_TABLE = "recombination_screen_runs"


def record_screen_run(conn, meta, rows, accessions, branches_tested, binary):
	"""Append a provenance row describing this screen.

	``meta_data.recombination_status`` records *that* a sequence was screened,
	and ``reference_recombination`` records what was found. Neither records
	*how* - and when nothing is found the events table is empty, so the method
	disappears entirely. That matters here: `ripples-fast` prunes the breakpoint
	search instead of enumerating it, so "no evidence" from a pruned search is a
	weaker statement than "no evidence" from an exhaustive one, and a reader six
	months from now cannot tell the two apart from a status column alone.

	Rows accumulate rather than replace, so re-screening leaves a history.
	"""
	conn.execute(
		f"CREATE TABLE IF NOT EXISTS {SCREEN_RUNS_TABLE} ("
		f"screened_at TEXT, ripples_binary TEXT, search_mode TEXT, "
		f"references_screened INTEGER, branches_tested INTEGER, "
		f"events_found INTEGER, accessions_flagged INTEGER, "
		f"branch_length INTEGER, min_range INTEGER, max_range INTEGER, "
		f"parsimony_improvement INTEGER, num_descendants INTEGER, "
		f"chunks TEXT, note TEXT)"
	)
	flagged = len({r["primary_accession"] for r in rows if r.get("primary_accession")})
	events = len({r.get("recomb_node_id") for r in rows}) if rows else 0
	pruned = str(binary).endswith("-fast")
	note = (
		"pruned search: not guaranteed identical to the exhaustive scan, and no "
		"head-to-head is available on trees this divergent because exhaustive "
		"ripples does not finish a single branch"
		if pruned else
		"exhaustive breakpoint-pair enumeration"
	)
	conn.execute(
		f"INSERT INTO {SCREEN_RUNS_TABLE} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
		(
			time.strftime("%Y-%m-%d %H:%M:%S"), binary,
			"pruned" if pruned else "exhaustive",
			len(accessions), branches_tested, events, flagged,
			meta.get("branch_length (-l)"), meta.get("min_range (-r)"),
			meta.get("max_range (-R)"), meta.get("parsimony_improvement (-p)"),
			meta.get("num_descendants (-n)"), str(meta.get("chunks x threads")),
			note,
		),
	)
	conn.commit()
	return note


def read_events_from_db(conn):
	"""Read stored events back out, for --report_only."""
	try:
		columns = [row[1] for row in conn.execute(f"PRAGMA table_info({RECOMBINATION_TABLE})")]
	except sqlite3.Error:
		return []
	if not columns:
		return []
	rows = conn.execute(f"SELECT {','.join(columns)} FROM {RECOMBINATION_TABLE}")
	return [dict(zip(columns, row)) for row in rows]


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def _as_int(value):
	try:
		return int(value)
	except (TypeError, ValueError):
		return None


def parsimony_gain(row):
	"""How much better the two-parent explanation is. None when unscored.

	This is the number an event should be believed in proportion to; it is what
	``--parsimony_improvement`` puts a floor under.
	"""
	original = _as_int(row.get("original_parsimony"))
	recombinant = _as_int(row.get("recomb_parsimony"))
	if original is None or recombinant is None:
		return None
	return original - recombinant


def breakpoint_histogram(rows, bins=40, width=48):
	"""Text histogram of breakpoint midpoints along the genome.

	Clustering of breakpoints in one region is the thing worth eyeballing - a
	real recombination hotspot looks different from calls scattered uniformly,
	which is more often a sign the parsimony floor is too low.
	"""
	points = []
	for row in rows:
		for prefix in ("breakpoint_1", "breakpoint_2"):
			start = _as_int(row.get(f"{prefix}_start"))
			end = _as_int(row.get(f"{prefix}_end"))
			if start is not None and end is not None:
				points.append((start + end) / 2.0)
	if not points:
		return []
	lo, hi = min(points), max(points)
	if hi <= lo:
		hi = lo + 1
	counts = [0] * bins
	for point in points:
		index = int((point - lo) / (hi - lo) * bins)
		counts[min(index, bins - 1)] += 1
	peak = max(counts) or 1
	lines = []
	for i, count in enumerate(counts):
		if not count:
			continue
		start = lo + (hi - lo) * i / bins
		end = lo + (hi - lo) * (i + 1) / bins
		bar = "#" * max(1, int(round(count / peak * width)))
		lines.append(f"{int(start):>7}-{int(end):<7} {bar} {count}")
	return lines


def build_report(db, rows, accessions, context, status_counts, meta, screened=True,
                 branches_tested=None):
	"""Render the whole report as Markdown."""
	flagged = sorted({r["primary_accession"] for r in rows if r.get("primary_accession")})
	unattributed = [r for r in rows if not r.get("primary_accession")]
	events = {r.get("recomb_node_id") for r in rows}

	out = []
	add = out.append

	add(f"# Recombination screen - {os.path.basename(db) if db else 'reference set'}")
	add("")
	if not screened:
		add(f"**Not screened.** This database has no recombination results stored: "
		    f"none of its {len(accessions)} reference sequences carry a "
		    f"`{screen.RECOMBINATION_STATUS_COLUMN}`. This is *not* a clean result - "
		    f"nothing has been tested yet. Run the screen without `--report_only`.")
	elif branches_tested == 0:
		add(f"**Nothing was tested.** No branch in the {len(accessions)}-sequence "
		    f"reference tree carried enough mutations to reach the "
		    f"`--branch_length` threshold, so RIPPLES had no candidate to examine. "
		    f"This is *not* a clean result - lower the threshold, or check that "
		    f"the mutation-annotated tree really holds the reference set.")
	elif not rows:
		add(f"**No recombination detected.** {len(accessions)} reference sequences were "
		    f"screened"
		    + (f" across {branches_tested} long branch(es)" if branches_tested else "")
		    + f" and RIPPLES found no branch better explained by two parents "
		      f"than by one, at the thresholds below.")
	else:
		add(f"**{len(events)} candidate event(s) across {len(flagged)} reference "
		    f"sequence(s)**, out of {len(accessions)} screened.")
	add("")

	# -- how it was run -----------------------------------------------------
	add("## How this was run")
	add("")
	add("| setting | value |")
	add("|---|---|")
	for key, value in meta.items():
		add(f"| {key} | {value} |")
	add("")

	# -- flagged references -------------------------------------------------
	if flagged:
		add("## Flagged references")
		add("")
		add("Ordered by parsimony gain - how much better the two-parent explanation "
		    "fits than the single-parent one. Treat the top of this table as the "
		    "strongest candidates, not as a pass/fail line.")
		add("")
		header = ["accession", "genotype", "subtype", "node", "breakpoint 1",
		          "breakpoint 2", "donor", "acceptor", "orig", "recomb", "gain"]
		add("| " + " | ".join(header) + " |")
		add("|" + "---|" * len(header))
		ranked = sorted(
			[r for r in rows if r.get("primary_accession")],
			key=lambda r: (parsimony_gain(r) is None, -(parsimony_gain(r) or 0)),
		)
		for row in ranked:
			accession = row["primary_accession"]
			info = context.get(accession, {})
			gain = parsimony_gain(row)
			add("| " + " | ".join(str(x) for x in [
				accession,
				info.get("nearest_reference_genotype") or "-",
				info.get("nearest_reference_subtype") or "-",
				row.get("recomb_node_id") or "-",
				row.get("breakpoint_1_interval") or "-",
				row.get("breakpoint_2_interval") or "-",
				row.get("donor_node_id") or "-",
				row.get("acceptor_node_id") or "-",
				row.get("original_parsimony") or "-",
				row.get("recomb_parsimony") or "-",
				"-" if gain is None else gain,
			]) + " |")
		add("")

	if unattributed:
		add("## Events with no listed descendants")
		add("")
		add(f"{len(unattributed)} event(s) were called on nodes that RIPPLES did not "
		    f"list descendants for. They are kept rather than dropped - losing a "
		    f"detection because the descendants file was incomplete would be worse "
		    f"than reporting it unattributed.")
		add("")
		for row in unattributed:
			add(f"- node `{row.get('recomb_node_id')}` "
			    f"bp1 `{row.get('breakpoint_1_interval')}` "
			    f"bp2 `{row.get('breakpoint_2_interval')}`")
		add("")

	# -- breakpoints --------------------------------------------------------
	histogram = breakpoint_histogram(rows)
	if histogram:
		add("## Breakpoint distribution")
		add("")
		add("Midpoint of every reported breakpoint interval, along the alignment.")
		add("")
		add("```")
		out.extend(histogram)
		add("```")
		add("")

	# -- coverage -----------------------------------------------------------
	if status_counts:
		add("## Screen coverage")
		add("")
		add("| status | references |")
		add("|---|---|")
		for status, count in status_counts:
			add(f"| {status} | {count} |")
		add("")
		add("`(not screened)` and `screened_no_evidence` are deliberately distinct. "
		    "A query that treated a missing screen as a clean one would be wrong.")
		add("")

	# -- caveats ------------------------------------------------------------
	binary = str(meta.get("ripples binary") or "")
	if binary.endswith("-fast"):
		add("## How the search was done")
		add("")
		add(f"This screen used `{binary}`, which **prunes** the breakpoint search "
		    f"rather than enumerating it. That is what makes it finish: on a "
		    f"reference set spanning this much divergence, exhaustive `ripples` "
		    f"evaluates ~2.6 million breakpoint pairs per branch at ~2.3/second, "
		    f"about 13 days per branch.")
		add("")
		add("The consequence is that a pruned search is not guaranteed to agree "
		    "with the exhaustive one, and here there is no head-to-head to check "
		    "it against - exhaustive `ripples` does not finish a single branch on "
		    "this data. **An empty result is therefore weaker evidence than it "
		    "looks.**")
		add("")

	if not rows:
		return "\n".join(out)

	add("## Before acting on this")
	add("")
	add("These are **candidates**. RIPPLES proposes a two-parent explanation that "
	    "fits the mutations better than one parent; that is evidence, not proof.")
	add("")
	add("Before removing anything from the reference set:")
	add("")
	add("1. Confirm with an independent method - 3SEQ, RDP, or a windowed identity "
	    "scan against genotype prototypes.")
	add("2. Check the breakpoint falls somewhere biologically plausible. For HCV "
	    "that is usually the NS2/NS3 junction.")
	add("3. Weigh the parsimony gain. A gain barely above "
	    "`--parsimony_improvement` is a weak call by construction.")
	add("")
	add("A recombinant left in the reference set mis-genotypes every query that "
	    "matches it, so a confirmed hit matters - but so does not removing a good "
	    "reference on a marginal score.")
	add("")
	return "\n".join(out)


def write_events_tsv(path, rows):
	"""Write the flat per-(accession, event) table beside the Markdown."""
	columns = [
		"primary_accession", "recomb_node_id",
		"breakpoint_1_interval", "breakpoint_2_interval",
		"breakpoint_1_start", "breakpoint_1_end",
		"breakpoint_2_start", "breakpoint_2_end",
		"donor_node_id", "donor_is_sibling", "donor_parsimony",
		"acceptor_node_id", "acceptor_is_sibling", "acceptor_parsimony",
		"original_parsimony", "min_starting_parsimony", "recomb_parsimony",
		"parsimony_gain",
	]
	with open(path, "w", encoding="utf-8") as handle:
		handle.write("\t".join(columns) + "\n")
		for row in rows:
			record = dict(row)
			record["parsimony_gain"] = parsimony_gain(row)
			handle.write("\t".join(
				"" if record.get(c) is None else str(record.get(c)) for c in columns
			) + "\n")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

class RecombinationHunter:
	def __init__(self, db, outdir="recombination_hunt", tree=None, alignment=None,
	             chunks=DEFAULT_CHUNKS, threads_per_chunk=DEFAULT_THREADS_PER_CHUNK,
	             branch_length=3, min_range=1000, max_range=10000000,
	             parsimony_improvement=3, num_descendants=3,
	             write_db=False, skip_build=False, timeout=None,
	             probe_timeout=PROBE_TIMEOUT,
	             ripples_binary=DEFAULT_RIPPLES_BINARY):
		self.db = db
		self.outdir = outdir
		self.tree = tree
		self.alignment = alignment
		requested = (chunks, threads_per_chunk)
		self.chunks, self.threads_per_chunk = plan_thread_budget(chunks, threads_per_chunk)
		if (self.chunks, self.threads_per_chunk) != requested:
			print(f"[budget] {requested[0]}x{requested[1]} threads exceeds the "
			      f"{MAX_THREADS}-thread ceiling; using "
			      f"{self.chunks}x{self.threads_per_chunk}")
		self.write_db = write_db
		self.skip_build = skip_build
		self.timeout = timeout
		self.probe_timeout = probe_timeout
		self.ripples_binary = ripples_binary
		self.accessions = []
		# Everything RIPPLES itself is parameterised by, passed explicitly rather
		# than left to the binary's defaults, which change between versions.
		self.ripples_kwargs = dict(
			branch_length=branch_length,
			min_range=min_range,
			max_range=max_range,
			parsimony_improvement=parsimony_improvement,
			num_descendants=num_descendants,
		)

	def build_mat(self, base, fasta, tree_path):
		"""Build the MAT, insisting on knowing where it actually landed.

		UShER treats ``-o`` as a path relative to the *process* working
		directory, while ``-d`` only redirects its auxiliary output. Asking for
		``reference.pb`` inside an ``--outdir`` therefore writes the MAT next to
		wherever the command happened to be launched, and the screen then hands
		RIPPLES a path to a file that does not exist.

		RIPPLES does not treat that as an error. It warns ``Tree found empty``,
		finds zero long branches and exits 0 - so the run reports a clean
		reference set without having tested a single branch. That is the one
		outcome this whole tool exists to prevent, so the MAT is built with an
		absolute ``-o`` and then checked before anything is allowed to use it.
		"""
		mat = os.path.abspath(os.path.join(self.outdir, "reference.pb"))
		if self.skip_build and os.path.exists(mat) and os.path.getsize(mat) > 0:
			print(f"[mat] reusing {mat}")
			return mat

		vcf = os.path.join(self.outdir, "reference_set.vcf")
		screen._require_binary("faToVcf")
		screen._require_binary("usher")
		subprocess.run(screen.build_fatovcf_command(fasta, vcf), check=True)
		subprocess.run(
			screen.build_usher_command(tree_path, vcf, mat, self.outdir, base.threads),
			check=True,
		)

		if not os.path.exists(mat) or os.path.getsize(mat) == 0:
			# Do not fall back to hunting for it: if UShER put the MAT somewhere
			# unexpected we want to be told, not to screen a file we cannot name.
			raise SystemExit(
				f"UShER reported success but no usable MAT is at {mat}. "
				f"Refusing to run RIPPLES against a missing tree - it would "
				f"report 'no recombination' without testing anything."
			)
		print(f"[mat] built {mat} ({os.path.getsize(mat) / 1e6:.1f} MB)")
		return mat

	def _base_screen(self):
		return screen.ReferenceRecombinationScreen(
			db=self.db,
			alignment=self.alignment,
			tree=self.tree,
			outdir=self.outdir,
			threads=min(MAX_THREADS, self.chunks * self.threads_per_chunk),
			skip_build=self.skip_build,
			**self.ripples_kwargs,
		)

	def run(self):
		os.makedirs(self.outdir, exist_ok=True)
		started = time.time()

		base = self._base_screen()
		fasta, tree_path = base.prepare_inputs()
		mat = self.build_mat(base, fasta, tree_path)
		self.accessions = base.accessions

		# One probe decides the split for the whole run.
		branches = probe_long_branches(
			mat, self.outdir, self.threads_per_chunk, self.ripples_kwargs,
			binary=self.ripples_binary, timeout=self.probe_timeout,
		)
		if branches == 0:
			# With a MAT that has been checked, this is a real answer rather than
			# a broken run - but it is "nothing met the threshold", which is not
			# the same claim as "screened and clean", and the report says so.
			print("[plan] no branch carries >= "
			      f"{self.ripples_kwargs['branch_length']} mutations; nothing to test")
			bounds = []
		elif branches is None:
			print("[plan] branch count unknown - running a single unchunked RIPPLES")
			bounds = [(None, None)]
		elif self.chunks == 1:
			# One chunk is the whole run, so say so by omitting -S/-E entirely.
			# Both binaries mark those flags EXPERIMENTAL; there is no reason to
			# exercise them when nothing is being split.
			bounds = [(None, None)]
			print(f"[plan] {branches} long branches in a single unchunked run "
			      f"x {self.threads_per_chunk} thread(s)")
		else:
			bounds = screen.chunk_bounds(branches, self.chunks) or [(None, None)]
			print(f"[plan] {branches} long branches over {len(bounds)} chunk(s) "
			      f"x {self.threads_per_chunk} thread(s)")

		chunk_dirs = []
		if bounds:
			chunk_dirs = self._run_bounds(mat, bounds)

		events, descendants = merge_chunks(chunk_dirs) if chunk_dirs else ([], {})
		rows = screen.attribute_events_to_accessions(events, descendants)
		elapsed = time.time() - started
		print(f"[result] {len(events)} event(s), "
		      f"{len({r['primary_accession'] for r in rows if r.get('primary_accession')})} "
		      f"accession(s) flagged, in {elapsed / 3600:.2f} h")

		meta = {
			"database": self.db or "(alignment input)",
			"reference sequences": len(self.accessions),
			"long branches tested": branches if branches is not None else "unknown",
			"chunks x threads": f"{len(bounds)} x {self.threads_per_chunk}",
			"ripples binary": self.ripples_binary,
			"wall clock": f"{elapsed / 3600:.2f} h",
			"branch_length (-l)": self.ripples_kwargs["branch_length"],
			"num_descendants (-n)": self.ripples_kwargs["num_descendants"],
			"min_range (-r)": self.ripples_kwargs["min_range"],
			"max_range (-R)": self.ripples_kwargs["max_range"],
			"parsimony_improvement (-p)": self.ripples_kwargs["parsimony_improvement"],
			"finished": time.strftime("%Y-%m-%d %H:%M:%S"),
		}

		if self.write_db and self.db:
			conn = sqlite3.connect(self.db)
			try:
				written, flagged = screen.store_results(
					conn, rows, self.accessions, detected_by=self.ripples_binary
				)
				record_screen_run(conn, meta, rows, self.accessions,
				                  branches, self.ripples_binary)
				print(f"[db] wrote {written} row(s) to {RECOMBINATION_TABLE}; "
				      f"{flagged} marked {STATUS_RECOMBINANT}; "
				      f"provenance in {SCREEN_RUNS_TABLE}")
			finally:
				conn.close()

		self.branches_tested = branches
		return rows, meta

	def _run_bounds(self, mat, bounds):
		"""Run one RIPPLES process per branch range, concurrently."""
		chunk_dirs = []
		with ThreadPoolExecutor(max_workers=len(bounds)) as pool:
			futures = []
			for index, (lo, hi) in enumerate(bounds):
				chunk_dir = os.path.join(self.outdir, f"chunk_{index:03d}")
				futures.append(pool.submit(
					run_chunk, self.ripples_binary, mat, chunk_dir, lo, hi,
					self.threads_per_chunk, self.ripples_kwargs, self.timeout,
				))
			for future in futures:
				chunk_dirs.append(future.result())
		return chunk_dirs

	def collect(self):
		"""Merge chunk output that already exists, without re-running RIPPLES.

		A chunked screen is hours long and chunks fail independently. Being able
		to merge, store and report from whatever is on disk means a lost chunk
		costs one chunk, not the whole run - and it lets the results be written
		to the database later, when nothing else holds it open.
		"""
		chunk_dirs = sorted(
			os.path.join(self.outdir, name)
			for name in os.listdir(self.outdir)
			if name.startswith("chunk_")
			and os.path.isdir(os.path.join(self.outdir, name))
		)
		if not chunk_dirs:
			raise SystemExit(f"no chunk_* directories under {self.outdir}")
		events, descendants = merge_chunks(chunk_dirs)
		rows = screen.attribute_events_to_accessions(events, descendants)
		print(f"[collect] {len(chunk_dirs)} chunk(s) -> {len(events)} event(s)")

		if self.db:
			conn = sqlite3.connect(f"file:{self.db}?mode=ro", uri=True)
			try:
				self.accessions = [
					r[0] for r in conn.execute(
						"SELECT primary_accession FROM meta_data "
						"WHERE lower(COALESCE(accession_type,'')) IN ('reference','master')"
					)
				]
			finally:
				conn.close()

		if self.write_db and self.db:
			conn = sqlite3.connect(self.db)
			try:
				written, flagged = screen.store_results(
					conn, rows, self.accessions, detected_by=self.ripples_binary
				)
				print(f"[db] wrote {written} row(s) to {RECOMBINATION_TABLE}; "
				      f"{flagged} marked {STATUS_RECOMBINANT}")
			finally:
				conn.close()

		return rows, {
			"database": self.db or "(alignment input)",
			"reference sequences": len(self.accessions),
			"chunks merged": len(chunk_dirs),
			"ripples binary": self.ripples_binary,
			"source": f"existing chunk output under {self.outdir} (--collect_only)",
			"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
		}

	def report(self, rows, meta, screened=True):
		context, status_counts = {}, None
		if self.db:
			conn = sqlite3.connect(f"file:{self.db}?mode=ro", uri=True)
			try:
				context = gather_context(
					conn, [r.get("primary_accession") for r in rows]
				)
				status_counts = gather_status_counts(conn)
			finally:
				conn.close()

		text = build_report(self.db, rows, self.accessions, context, status_counts,
		                    meta, screened=screened,
		                    branches_tested=getattr(self, "branches_tested", None))
		markdown = os.path.join(self.outdir, "recombination_report.md")
		tsv = os.path.join(self.outdir, "recombination_events.tsv")
		with open(markdown, "w", encoding="utf-8") as handle:
			handle.write(text)
		write_events_tsv(tsv, rows)
		print(f"[report] {markdown}")
		print(f"[report] {tsv}")
		return markdown, tsv


def parse_args(argv=None):
	parser = argparse.ArgumentParser(
		description="Run the UShER/RIPPLES recombination screen over a V-gTK "
		            "database in parallel chunks, and report the result."
	)
	parser.add_argument("--db", help="SQLite database holding the reference set")
	parser.add_argument("--alignment", help="reference alignment FASTA (instead of --db)")
	parser.add_argument("--tree", help="starting tree; derived from --db when omitted")
	parser.add_argument("--outdir", default="recombination_hunt")
	parser.add_argument("--chunks", type=int, default=DEFAULT_CHUNKS,
	                    help=f"concurrent RIPPLES processes (default {DEFAULT_CHUNKS})")
	parser.add_argument("--threads_per_chunk", type=int, default=DEFAULT_THREADS_PER_CHUNK,
	                    help=f"threads each chunk gets; chunks x threads is capped "
	                         f"at {MAX_THREADS} (default {DEFAULT_THREADS_PER_CHUNK})")
	parser.add_argument("--branch_length", type=int, default=3, help="RIPPLES -l")
	parser.add_argument("--min_range", type=int, default=1000, help="RIPPLES -r")
	parser.add_argument("--max_range", type=int, default=10000000, help="RIPPLES -R")
	parser.add_argument("--parsimony_improvement", type=int, default=3, help="RIPPLES -p")
	parser.add_argument("--num_descendants", type=int, default=3, help="RIPPLES -n")
	parser.add_argument("--write_db", action="store_true",
	                    help="store results in the database as well as reporting")
	parser.add_argument("--skip_build", action="store_true",
	                    help="reuse an existing reference.pb in --outdir")
	parser.add_argument("--timeout", type=float, default=None,
	                    help="seconds before each chunk is stopped; partial results kept")
	parser.add_argument("--probe_timeout", type=float, default=PROBE_TIMEOUT,
	                    help="seconds to wait for the RIPPLES branch count")
	parser.add_argument("--ripples_binary", default=DEFAULT_RIPPLES_BINARY,
	                    help=f"which RIPPLES binary to drive (default "
	                         f"{DEFAULT_RIPPLES_BINARY}). Use 'ripples' for the "
	                         f"exhaustive search - only viable on trees whose "
	                         f"branches carry few mutations.")
	parser.add_argument("--collect_only", action="store_true",
	                    help="merge/store/report existing chunk_* output without "
	                         "re-running RIPPLES")
	parser.add_argument("--report_only", action="store_true",
	                    help="rebuild the report from an already-screened database")
	return parser.parse_args(argv)


def main(argv=None):
	args = parse_args(argv)
	if not args.db and not args.alignment:
		raise SystemExit("one of --db or --alignment is required")
	if args.report_only and not args.db:
		raise SystemExit("--report_only needs --db")

	hunter = RecombinationHunter(
		db=args.db,
		alignment=args.alignment,
		tree=args.tree,
		outdir=args.outdir,
		chunks=args.chunks,
		threads_per_chunk=args.threads_per_chunk,
		branch_length=args.branch_length,
		min_range=args.min_range,
		max_range=args.max_range,
		parsimony_improvement=args.parsimony_improvement,
		num_descendants=args.num_descendants,
		write_db=args.write_db,
		skip_build=args.skip_build,
		timeout=args.timeout,
		probe_timeout=args.probe_timeout,
		ripples_binary=args.ripples_binary,
	)

	if args.report_only:
		conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
		try:
			rows = read_events_from_db(conn)
			screened = count_screened(conn) > 0
			hunter.accessions = [
				r[0] for r in conn.execute(
					"SELECT primary_accession FROM meta_data "
					"WHERE lower(COALESCE(accession_type,'')) IN ('reference','master')"
				)
			]
		finally:
			conn.close()
		os.makedirs(args.outdir, exist_ok=True)
		meta = {
			"database": args.db,
			"reference sequences": len(hunter.accessions),
			"source": f"existing {RECOMBINATION_TABLE} table (--report_only)",
			"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
		}
		meta["screened"] = "yes" if screened else "no - never run"
		hunter.report(rows, meta, screened=screened)
		return 0

	if args.collect_only:
		rows, meta = hunter.collect()
	else:
		rows, meta = hunter.run()
	hunter.report(rows, meta)
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
