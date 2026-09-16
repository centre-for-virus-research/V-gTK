#!/usr/bin/env python3
"""Fail an update run at the start if its reference list has outgrown the update DB.

Two checks, both against the ref list's master/reference rows:

1. **meta_data** - every reference must already be stored as a master or
   reference.

2. **the stored UShER tree** - every reference must be a tip. Clade assignment
   gives each query the clade of its nearest labelled *reference tip*, and
   update mode places new sequences onto the tree already in the DB. A reference
   added to the list after that tree was built is absent from it, so clades would
   be assigned as if it did not exist - silently. UShER trees hold every placed
   sample, so a complete build has every reference in them.

   IQ-TREE trees are deliberately not checked: they hold only cluster
   representatives, and references are routinely missing from them (13 of 28 in
   the RABV update test DB, 220 of 238 in HCV_test).

A DB with no stored UShER tree (a tree-free build) has nothing to check against,
and passes check 2 with a note.
"""

import re
import sqlite3
from argparse import ArgumentParser

import pandas as pd

import accession_utils
from ExportRefListFromUpdateDb import load_reference_file_table

#: A tip label is whatever follows '(' or ',' up to the next ':', ',' or ')'.
_TIP = re.compile(r"[(,]\s*('(?:[^']|'')*'|[^():,;\s]+)")


def _key(value):
	value = str(value).strip().strip("'")
	return accession_utils.normalise_accession(value) or value


def newick_tips(newick):
	"""Tip labels of a Newick string, normalised to bare accessions."""
	return {_key(label) for label in _TIP.findall(newick or "")}


def usher_trees(conn):
	"""``[(name, tips)]`` for every stored UShER tree; ``[]`` when there is no trees table."""
	tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
	if "trees" not in tables:
		return []
	rows = conn.execute("SELECT name, source, newick FROM trees").fetchall()
	return [(name, newick_tips(newick)) for name, source, newick in rows
			if "usher" in f"{source} {name}".lower() and newick]


def main(args):
	conn = sqlite3.connect(args.db)
	try:
		df = pd.read_sql_query("SELECT primary_accession, accession_type FROM meta_data", conn)
		trees = usher_trees(conn)
	finally:
		conn.close()

	db_refs = set(
		df[df["accession_type"].fillna("").str.lower().isin(["master", "reference"])]["primary_accession"].astype(str).str.strip().tolist()
	)

	ref_df = load_reference_file_table(args.ref_list)
	if ref_df.shape[1] < 2:
		raise ValueError("ref_list must have at least 2 columns: accession and type")
	ref_df["primary_accession"] = ref_df["primary_accession"].fillna("").astype(str).str.strip()
	ref_df["accession_type"] = ref_df["accession_type"].fillna("").astype(str).str.lower().str.strip()
	required = set(ref_df[ref_df["accession_type"].isin(["master", "reference"])] ["primary_accession"].tolist())
	missing = sorted([x for x in required if x and x not in db_refs])
	if missing:
		raise ValueError("Reference mismatch: ref_list entries not present in DB as master/reference: " + ", ".join(missing))
	print("Reference list matches DB master/reference set")

	if not trees:
		print("No UShER tree stored in the update DB: reference tips not checked")
		return
	# A segmented DB stores one UShER tree per segment, so a reference only has
	# to be a tip in one of them.
	tips = set().union(*(t for _, t in trees))
	not_in_tree = sorted(x for x in required if x and _key(x) not in tips)
	if not_in_tree:
		shown = ", ".join(not_in_tree[:20]) + (f" ... (+{len(not_in_tree) - 20} more)" if len(not_in_tree) > 20 else "")
		raise ValueError(
			f"Reference mismatch: {len(not_in_tree)} ref_list reference(s) are not tips in the update DB's UShER "
			f"tree ({', '.join(name for name, _ in trees)}): {shown}. The stored tree predates these references, "
			"so clade assignment would ignore them. Rebuild the database with this reference list, "
			"or use the reference list it was built with.")
	print(f"All {len(required)} references are tips in the stored UShER tree(s)")


if __name__ == "__main__":
	parser = ArgumentParser(description="Validate that ref_list master/reference entries exist in update DB and its UShER tree")
	parser.add_argument("--ref_list", required=True)
	parser.add_argument("--db", required=True)
	args = parser.parse_args()
	main(args)
