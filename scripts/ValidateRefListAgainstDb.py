#!/usr/bin/env python3
import sqlite3
import pandas as pd
from argparse import ArgumentParser

from ExportRefListFromUpdateDb import load_reference_file_table


def main(args):
	conn = sqlite3.connect(args.db)
	try:
		df = pd.read_sql_query("SELECT primary_accession, accession_type FROM meta_data", conn)
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


if __name__ == "__main__":
	parser = ArgumentParser(description="Validate that ref_list master/reference entries exist in update DB")
	parser.add_argument("--ref_list", required=True)
	parser.add_argument("--db", required=True)
	args = parser.parse_args()
	main(args)
