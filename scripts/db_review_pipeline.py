#!/usr/bin/env python3
"""Screen a database for recombination, midpoint-root its trees, validate it.

Each step is separate and each writes into ONE new database, so the source is
never modified and the result is a single file carrying every output:

    1. copy           source DB -> output DB (everything below is in-place on
                      the copy, so a failed run never damages the input)
    2. recombination  UShER/RIPPLES screen via recombination_hunter_db, storing
                      reference_recombination, meta_data.recombination_status
                      and recombination_screen_runs
    3. reroot         midpoint-root every tree row via TreeReRoot, writing each
                      back with ReplaceUsherTreeInDb
    4. validate       ValidateDbTree over the finished database

Steps are selectable with --steps, so a failed or slow one can be re-run alone
against the database the previous steps already produced.

Usage
-----
    python scripts/db_review_pipeline.py \\
        --db  test_out/HCV_XML_full_plus_update/HCV_full_new_Usher.db \\
        --out test_out/HCV_XML_full_plus_update/HCV_full_reviewed.db \\
        --workdir dev/db_review --threads 5
"""

import argparse
import os
import shutil
import sqlite3
import subprocess
import sys
import time

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS_DIR)

ALL_STEPS = ("copy", "recombination", "reroot", "validate")


def _say(message):
	print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def _run(command, label):
	"""Run a step, failing loudly rather than continuing on a broken database."""
	_say(f"$ {' '.join(str(c) for c in command)}")
	started = time.time()
	result = subprocess.run(command)
	elapsed = time.time() - started
	if result.returncode != 0:
		raise SystemExit(
			f"[{label}] exited {result.returncode} after {elapsed / 60:.1f} min. "
			f"The output database is left as it is so the failure can be "
			f"inspected; re-run with --steps to resume from here."
		)
	_say(f"{label} finished in {elapsed / 60:.1f} min")


# ---------------------------------------------------------------------------
# 1. copy
# ---------------------------------------------------------------------------

def step_copy(source, output, force=False):
	if os.path.exists(output):
		if not force:
			raise SystemExit(
				f"{output} already exists. Pass --force to overwrite, or point "
				f"--out somewhere else. Refusing to silently replace a database."
			)
		_say(f"overwriting existing {output}")
	if not os.path.exists(source):
		raise SystemExit(f"source database not found: {source}")
	size = os.path.getsize(source) / 1e9
	_say(f"copying {source} -> {output} ({size:.1f} GB)")
	started = time.time()
	shutil.copyfile(source, output)
	_say(f"copied in {time.time() - started:.0f}s")


# ---------------------------------------------------------------------------
# 2. recombination
# ---------------------------------------------------------------------------

def step_recombination(output, workdir, threads, binary, extra=()):
	_run(
		[
			sys.executable, os.path.join(SCRIPTS_DIR, "recombination_hunter_db.py"),
			"--db", output,
			"--outdir", os.path.join(workdir, "recombination"),
			# One process holding the whole thread budget: chunking exists to
			# parallelise the exhaustive binary, and ripples-fast finishes in
			# about a minute, so splitting it only adds moving parts.
			"--chunks", "1",
			"--threads_per_chunk", str(threads),
			"--ripples_binary", binary,
			"--write_db",
			*extra,
		],
		"recombination",
	)


# ---------------------------------------------------------------------------
# 3. midpoint rooting
# ---------------------------------------------------------------------------

def tree_rows(db):
	"""Every tree row, with the selectors ReplaceUsherTreeInDb matches on."""
	conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
	try:
		return [
			{
				"name": row[0] or "",
				"source": row[1] or "",
				"segment_key": row[2] or "",
				"segment": row[3] or "",
				"tips": (row[4] or "").count(",") + 1 if row[4] else 0,
			}
			for row in conn.execute(
				"SELECT name, source, segment_key, segment, newick FROM trees"
			)
		]
	finally:
		conn.close()


def step_reroot(output, workdir, order_node="increase", only_source=None):
	rows = tree_rows(output)
	if not rows:
		_say("no tree rows to reroot")
		return
	tree_dir = os.path.join(workdir, "trees")
	os.makedirs(tree_dir, exist_ok=True)

	conn = sqlite3.connect(f"file:{output}?mode=ro", uri=True)
	try:
		for row in rows:
			if only_source and row["source"] != only_source:
				_say(f"skipping {row['name']} (source {row['source']})")
				continue
			_say(f"rerooting {row['name']} ({row['tips']:,} tips, source={row['source']})")
			newick = conn.execute(
				"SELECT newick FROM trees WHERE COALESCE(name,'')=? "
				"AND COALESCE(source,'')=? AND COALESCE(segment_key,'')=? "
				"AND COALESCE(segment,'')=?",
				(row["name"], row["source"], row["segment_key"], row["segment"]),
			).fetchone()[0]

			raw = os.path.join(tree_dir, f"{row['name']}.original.nwk")
			rooted = os.path.join(tree_dir, f"{row['name']}.midpoint.nwk")
			with open(raw, "w", encoding="utf-8") as handle:
				handle.write(newick.strip() + "\n")

			_run(
				[
					sys.executable, os.path.join(SCRIPTS_DIR, "TreeReRoot.py"),
					"--input_tree", raw,
					"--output_tree", rooted,
					"--order_node", order_node,
				],
				f"reroot:{row['name']}",
			)
			_run(
				[
					sys.executable, os.path.join(SCRIPTS_DIR, "ReplaceUsherTreeInDb.py"),
					"--db", output,
					"--tree", rooted,
					"--source", row["source"],
					"--name", row["name"],
					"--segment-key", row["segment_key"],
					"--segment", row["segment"],
				],
				f"store:{row['name']}",
			)
	finally:
		conn.close()


# ---------------------------------------------------------------------------
# 4. validate
# ---------------------------------------------------------------------------

def step_validate(output, workdir, extra=()):
	report_dir = os.path.join(workdir, "validation")
	os.makedirs(report_dir, exist_ok=True)
	_run(
		[
			sys.executable, os.path.join(SCRIPTS_DIR, "ValidateDbTree.py"),
			"--db", output,
			"--outdir", report_dir,
			*extra,
		],
		"validate",
	)
	return report_dir


# ---------------------------------------------------------------------------

def summarise(output):
	"""What the finished database now carries, read back from the file itself."""
	conn = sqlite3.connect(f"file:{output}?mode=ro", uri=True)
	lines = []
	try:
		def table_exists(name):
			return bool(list(conn.execute(
				"SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (name,)
			)))

		if table_exists("recombination_screen_runs"):
			for row in conn.execute(
				"SELECT screened_at, ripples_binary, search_mode, branches_tested, "
				"events_found, accessions_flagged FROM recombination_screen_runs "
				"ORDER BY screened_at DESC LIMIT 1"
			):
				lines.append(
					f"recombination: {row[4]} event(s) over {row[3]} branch(es), "
					f"{row[5]} reference(s) flagged, via {row[1]} ({row[2]} search) at {row[0]}"
				)
		columns = {r[1] for r in conn.execute("PRAGMA table_info(meta_data)")}
		if "recombination_status" in columns:
			for status, count in conn.execute(
				"SELECT COALESCE(recombination_status,'(not screened)'), COUNT(*) "
				"FROM meta_data WHERE lower(COALESCE(accession_type,'')) "
				"IN ('reference','master') GROUP BY 1 ORDER BY 2 DESC"
			):
				lines.append(f"  {status}: {count}")
		for name, source, created in conn.execute(
			"SELECT name, source, created_at FROM trees"
		):
			lines.append(f"tree {name} (source={source}) last written {created}")
	finally:
		conn.close()
	return lines


def parse_args(argv=None):
	parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
	parser.add_argument("--db", required=True, help="source database (never modified)")
	parser.add_argument("--out", required=True, help="new database carrying every output")
	parser.add_argument("--workdir", default="db_review",
	                    help="intermediates, reports and logs")
	parser.add_argument("--threads", type=int, default=5,
	                    help="thread budget for the recombination screen (default 5)")
	parser.add_argument("--ripples_binary", default="ripples-fast")
	parser.add_argument("--order_node", default="increase",
	                    choices=["increase", "decrease", "none"])
	parser.add_argument("--reroot_source", default=None,
	                    help="only reroot trees with this source (default: all)")
	parser.add_argument("--steps", default=",".join(ALL_STEPS),
	                    help=f"comma-separated subset of {','.join(ALL_STEPS)}")
	parser.add_argument("--force", action="store_true",
	                    help="overwrite --out if it already exists")
	parser.add_argument("--validator_arg", action="append", default=[],
	                    help="extra flag passed through to ValidateDbTree "
	                         "(repeatable)")
	return parser.parse_args(argv)


def main(argv=None):
	args = parse_args(argv)
	steps = [s.strip() for s in args.steps.split(",") if s.strip()]
	unknown = [s for s in steps if s not in ALL_STEPS]
	if unknown:
		raise SystemExit(f"unknown step(s) {unknown}; choose from {list(ALL_STEPS)}")

	os.makedirs(args.workdir, exist_ok=True)
	overall = time.time()
	_say(f"steps: {', '.join(steps)}")

	if "copy" in steps:
		step_copy(args.db, args.out, force=args.force)
	elif not os.path.exists(args.out):
		raise SystemExit(
			f"{args.out} does not exist and 'copy' is not in --steps. Run the "
			f"copy step first, or point --out at the database to work on."
		)

	if "recombination" in steps:
		step_recombination(args.out, args.workdir, args.threads, args.ripples_binary)
	if "reroot" in steps:
		step_reroot(args.out, args.workdir, args.order_node, args.reroot_source)
	if "validate" in steps:
		step_validate(args.out, args.workdir, args.validator_arg)

	_say(f"all steps done in {(time.time() - overall) / 60:.1f} min -> {args.out}")
	for line in summarise(args.out):
		print("   " + line)
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
