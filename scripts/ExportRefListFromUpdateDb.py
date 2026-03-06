#!/usr/bin/env python3
import os
import sqlite3
from argparse import ArgumentParser

import pandas as pd


def _segment_norm(value):
	if value is None:
		return "0"
	s = str(value).strip()
	if not s or s.lower() == "nan":
		return "0"
	digits = ''.join(ch for ch in s if ch.isdigit())
	return digits if digits else s


def load_reference_rows(update_db):
	if not os.path.isfile(update_db):
		raise FileNotFoundError(f"Update DB not found: {update_db}")

	conn = sqlite3.connect(update_db)
	try:
		df = pd.read_sql_query("SELECT * FROM meta_data", conn)
	finally:
		conn.close()

	if "primary_accession" not in df.columns:
		raise ValueError("Update DB table meta_data must contain column: primary_accession")
	if "accession_type" not in df.columns:
		df["accession_type"] = "query"
	if "segment" not in df.columns:
		df["segment"] = "0"

	if "accession_type" not in df.columns:
		df["accession_type"] = "query"

	df["segment"] = df["segment"].apply(_segment_norm)
	df["primary_accession"] = df["primary_accession"].fillna("").astype(str).str.strip()
	df = df[df["primary_accession"] != ""]

	return df[
		df["accession_type"].fillna("").astype(str).str.lower().isin(["master", "reference"])
	][["primary_accession", "accession_type", "segment"]].drop_duplicates()


def load_reference_dict(update_db):
	refs = load_reference_rows(update_db)
	return {
		str(row["primary_accession"]).strip(): str(row["accession_type"]).strip()
		for _, row in refs.iterrows()
	}


def load_master_accessions(update_db):
	refs = load_reference_rows(update_db)
	masters = refs[refs["accession_type"].fillna("").astype(str).str.lower() == "master"]
	if masters.empty:
		return refs["primary_accession"].astype(str).str.strip().tolist()
	return masters["primary_accession"].astype(str).str.strip().tolist()


def export_ref_list(update_db, output_path):
	refs = load_reference_rows(update_db)
	refs.to_csv(output_path, sep="\t", index=False, header=False)


if __name__ == "__main__":
	parser = ArgumentParser(description="Export reference list from an existing update DB")
	parser.add_argument("--db", required=True, help="Existing SQLite DB path")
	parser.add_argument("-o", "--output", required=True, help="Output TSV path")
	args = parser.parse_args()
	export_ref_list(args.db, args.output)