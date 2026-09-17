#!/usr/bin/env python3
"""Surgically add new references to a COPY of an existing database.

********************************************************************************
*  WORK IN PROGRESS - NOT EXPECTED TO WORK YET.                                *
*  Only `plan`, `export-backbone` and `verify-unchanged` are implemented and    *
*  tested. `fetch`, `align` and `insert` are partial sketches. Do not use the  *
*  output database for anything that matters.                                 *
********************************************************************************

WHY THIS EXISTS
===============
New, hyper-divergent lineages are discovered or evolve. Each needs a reference
so its sequences are aligned and typed against something close. A full build
with the new reference list is the complete answer - it rebuilds trees and
reassigns every clade - but it changes everything.

This tool is the opposite: the smallest change that puts a new reference into a
database, and a proof that nothing that already existed was altered.

    - the original database is never modified; a copy is written
    - existing references, masters, segments and genotypes must not change
    - the existing backbone's column layout is kept: new references are added to
      it with MAFFT --add --keeplength, so no stored alignment moves
    - trees, clade assignments and mutation calls are NOT rebuilt; run a full
      build for those
    - every row of every table in the original must be byte-identical in the copy

It is deliberately separate from the main pipeline: update runs refuse a
reference list that adds references (ValidateRefListAgainstDb.py).

KNOWN OPEN QUESTION
===================
A reference added here is not a tip in the stored UShER tree, so a later normal
update run on the new database fails ValidateRefListAgainstDb's tree check. The
added references are recorded in a ``reference_update_log`` table so that check
can decide how to treat them; that decision has not been made.

Subcommands
===========
plan              diff the database's references against a new list -> new_references.tsv
export-backbone   write the database's per-segment backbone (refset_<segment>_aln.fasta)
fetch             download the new references' records            [SKETCH]
align             add them to the backbone, keeping its columns   [SKETCH]
insert            write a copy of the database with the new rows  [SKETCH]
verify-unchanged  prove every original row is unchanged in the copy
"""

import argparse
import csv
import hashlib
import os
import shutil
import sqlite3
import subprocess
import sys

import pandas as pd

import accession_utils
import segment_utils
from ExportRefListFromUpdateDb import load_reference_file_table

REFERENCE_TYPES = ("master", "reference")
LOG_TABLE = "reference_update_log"


def _key(value):
	value = str(value or "").strip()
	return accession_utils.normalise_accession(value) or value


def _segment(value):
	label = segment_utils.normalise_segment(value)
	if label is None or label.casefold() in segment_utils.PANDAS_NULL_TOKENS:
		return ""
	return label


# ---------------------------------------------------------------------------
# plan
# ---------------------------------------------------------------------------

class PlanError(ValueError):
	"""The new list asks for more than an addition."""


def _list_references(ref_list):
	table = load_reference_file_table(ref_list)
	rows = {}
	for _, row in table.iterrows():
		kind = str(row.get("accession_type", "")).strip().lower()
		if kind not in REFERENCE_TYPES:
			continue
		rows[_key(row["primary_accession"])] = {
			"accession_type": kind,
			"segment": _segment(row.get("segment", "")),
			"genotype": str(row.get("genotype", "") or "").strip(),
			"subtype": str(row.get("subtype", "") or "").strip(),
		}
	return rows


def _db_references(db):
	conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
	try:
		meta = pd.read_sql_query("SELECT * FROM meta_data", conn).fillna("")
	finally:
		conn.close()
	meta = meta[meta["accession_type"].astype(str).str.strip().str.lower().isin(REFERENCE_TYPES)]
	return {
		_key(row["primary_accession"]): {
			"accession_type": str(row["accession_type"]).strip().lower(),
			"segment": _segment(row.get("segment", "")),
		}
		for _, row in meta.iterrows()
	}


def plan(db, ref_list):
	"""New references to add, as dicts. Raises PlanError for anything but a pure addition."""
	in_db, in_list = _db_references(db), _list_references(ref_list)
	problems = []
	removed = sorted(set(in_db) - set(in_list))
	if removed:
		problems.append(f"{len(removed)} reference(s) in the database are not in the new list: {removed[:10]}")
	for acc in sorted(set(in_db) & set(in_list)):
		old, new = in_db[acc], in_list[acc]
		if old["accession_type"] != new["accession_type"]:
			problems.append(f"{acc}: {old['accession_type']} in the database, {new['accession_type']} in the list")
		if old["segment"] and new["segment"] and old["segment"] != new["segment"]:
			problems.append(f"{acc}: segment {old['segment']} in the database, {new['segment']} in the list")
	added = {acc: row for acc, row in in_list.items() if acc not in in_db}
	if any(row["accession_type"] == "master" for row in added.values()):
		problems.append("a new master changes every coordinate in its segment; that needs a full build")
	if problems:
		raise PlanError("Not a surgical addition:\n  " + "\n  ".join(problems))
	return [dict(primary_accession=acc, **row) for acc, row in sorted(added.items())]


def write_plan(rows, path):
	columns = ["primary_accession", "accession_type", "segment", "genotype", "subtype"]
	with open(path, "w", newline="") as handle:
		writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
		writer.writeheader()
		writer.writerows(rows)


# ---------------------------------------------------------------------------
# export-backbone
# ---------------------------------------------------------------------------

def export_backbone(db, output_dir):
	"""Write refset_<segment>_aln.fasta from the database, exactly as update mode derives it."""
	from PadAlignment import PadAlignment
	pad = PadAlignment(None, ".", ".", ".", False, update_db=db)
	return pad.export_update_backbones(output_dir)


# ---------------------------------------------------------------------------
# fetch / align / insert  [SKETCHES]
# ---------------------------------------------------------------------------

def fetch(accessions, output_fasta, email):
	"""[SKETCH] FASTA for the new references. Metadata (host, date, country) is not fetched yet."""
	from Bio import Entrez
	Entrez.email = email
	with Entrez.efetch(db="nuccore", id=",".join(accessions), rettype="fasta", retmode="text") as handle:
		text = handle.read()
	with open(output_fasta, "w") as out:
		for line in text.splitlines():
			out.write((">" + _key(line[1:].split()[0]) if line.startswith(">") else line.lower()) + "\n")
	# TODO: GenBank XML for these accessions through GenBankParser.py, so the
	# meta_data row carries the same fields as every other reference.


def align(backbone_fasta, new_fasta, output_fasta, threads=2):
	"""[SKETCH] Add new references to a backbone without changing its columns.

	--keeplength drops any column the new sequence would insert. Those bases must
	be written to the insertions table (TODO: read --mapout), exactly as update
	mode does for queries.
	"""
	with open(output_fasta, "w") as out:
		subprocess.run(["mafft", "--add", new_fasta, "--keeplength", "--mapout", "--thread", str(threads),
						backbone_fasta], stdout=out, check=True)


def insert(db, new_db, plan_rows, aligned_fasta):
	"""[SKETCH] Copy the database and add the new references' rows. Never touches ``db``."""
	if os.path.abspath(db) == os.path.abspath(new_db):
		raise ValueError("refusing to write over the original database")
	shutil.copyfile(db, new_db)
	raise NotImplementedError(
		"insert is not written yet: meta_data, sequences, sequence_alignment, insertions, features "
		f"and the {LOG_TABLE} table for the new references")


# ---------------------------------------------------------------------------
# verify-unchanged
# ---------------------------------------------------------------------------

def _row_digests(conn, table):
	columns = [row[1] for row in conn.execute(f'PRAGMA table_info("{table}")')]
	counts = {}
	for row in conn.execute(f'SELECT * FROM "{table}"'):
		digest = hashlib.sha256(repr(tuple(row)).encode()).hexdigest()
		counts[digest] = counts.get(digest, 0) + 1
	return columns, counts


def verify_unchanged(original_db, new_db):
	"""``[problems]``: every table, column and row of the original must be present unchanged in the copy.

	Additions are allowed - new rows, new tables - but nothing may be removed or
	edited, and existing tables may not gain columns (that would change every row).
	"""
	problems = []
	old = sqlite3.connect(f"file:{original_db}?mode=ro", uri=True)
	new = sqlite3.connect(f"file:{new_db}?mode=ro", uri=True)
	try:
		old_tables = [r[0] for r in old.execute("SELECT name FROM sqlite_master WHERE type='table'")]
		new_tables = {r[0] for r in new.execute("SELECT name FROM sqlite_master WHERE type='table'")}
		for table in old_tables:
			if table not in new_tables:
				problems.append(f"{table}: table missing")
				continue
			old_cols, old_rows = _row_digests(old, table)
			new_cols, new_rows = _row_digests(new, table)
			if old_cols != new_cols:
				problems.append(f"{table}: columns changed {old_cols} -> {new_cols}")
				continue
			lost = sum(max(0, n - new_rows.get(d, 0)) for d, n in old_rows.items())
			if lost:
				problems.append(f"{table}: {lost} original row(s) removed or changed")
	finally:
		old.close()
		new.close()
	return problems


# ---------------------------------------------------------------------------

def main(argv=None):
	parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
	sub = parser.add_subparsers(dest="command", required=True)

	p = sub.add_parser("plan")
	p.add_argument("--db", required=True)
	p.add_argument("--ref_list", required=True)
	p.add_argument("--output", default="new_references.tsv")

	p = sub.add_parser("export-backbone")
	p.add_argument("--db", required=True)
	p.add_argument("--output_dir", default="db_backbone")

	p = sub.add_parser("fetch")
	p.add_argument("--plan", required=True)
	p.add_argument("--output", default="new_references.fasta")
	p.add_argument("--email", required=True)

	p = sub.add_parser("align")
	p.add_argument("--backbone", required=True)
	p.add_argument("--new_fasta", required=True)
	p.add_argument("--output", required=True)
	p.add_argument("--threads", type=int, default=2)

	p = sub.add_parser("insert")
	p.add_argument("--db", required=True)
	p.add_argument("--new_db", required=True)
	p.add_argument("--plan", required=True)
	p.add_argument("--aligned", required=True)

	p = sub.add_parser("verify-unchanged")
	p.add_argument("--original_db", required=True)
	p.add_argument("--new_db", required=True)

	args = parser.parse_args(argv)
	print("[AddReferencesToDb] WORK IN PROGRESS - not expected to work yet", file=sys.stderr)

	if args.command == "plan":
		try:
			rows = plan(args.db, args.ref_list)
		except PlanError as exc:
			print(f"ERROR: {exc}", file=sys.stderr)
			return 1
		write_plan(rows, args.output)
		print(f"{len(rows)} new reference(s) -> {args.output}")
		return 0 if rows else 1
	if args.command == "export-backbone":
		written = export_backbone(args.db, args.output_dir)
		print(f"backbone -> {written}")
		return 0 if written else 1
	if args.command == "fetch":
		accessions = pd.read_csv(args.plan, sep="\t", dtype=str)["primary_accession"].tolist()
		fetch(accessions, args.output, args.email)
		return 0
	if args.command == "align":
		align(args.backbone, args.new_fasta, args.output, args.threads)
		return 0
	if args.command == "insert":
		insert(args.db, args.new_db, pd.read_csv(args.plan, sep="\t", dtype=str).to_dict("records"), args.aligned)
		return 0
	if args.command == "verify-unchanged":
		problems = verify_unchanged(args.original_db, args.new_db)
		for problem in problems:
			print(f"CHANGED: {problem}", file=sys.stderr)
		print("verify-unchanged:", "FAIL" if problems else "every original row is unchanged")
		return 1 if problems else 0
	return 2


if __name__ == "__main__":
	sys.exit(main())
