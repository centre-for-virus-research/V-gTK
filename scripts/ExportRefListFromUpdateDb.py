#!/usr/bin/env python3
import os
import sqlite3
from argparse import ArgumentParser

import pandas as pd


REFLIST_COLUMN_ALIASES = {
	"primary_accession": "primary_accession",
	"accession": "primary_accession",
	"accession_id": "primary_accession",
	"acc": "primary_accession",
	"id": "primary_accession",
	"accession_type": "accession_type",
	"type": "accession_type",
	"status": "accession_type",
	"segment": "segment",
	"genotype": "genotype",
	"subtype": "subtype",
}


def _normalize_reflist_text(value):
	if value is None:
		return ""
	return str(value).strip()


def _default_reflist_column_names(count):
	base = ["primary_accession", "accession_type", "segment", "genotype", "subtype"]
	if count <= len(base):
		return base[:count]
	return base + [f"col{i}" for i in range(len(base), count)]


def _deduplicate_column_names(columns):
	seen = {}
	result = []
	for col in columns:
		base = col if col else "column"
		idx = seen.get(base, 0)
		result.append(base if idx == 0 else f"{base}_{idx}")
		seen[base] = idx + 1
	return result


def _looks_like_header_row(values):
	normalized = [_normalize_reflist_text(v).lower() for v in values]
	if not normalized:
		return False
	if normalized[0] not in {"primary_accession", "accession", "accession_id", "acc", "id"}:
		return False
	return any(value in {"accession_type", "type", "status", "segment", "genotype", "subtype"} for value in normalized[1:])


def load_reference_file_table(ref_list_path):
	if not os.path.isfile(ref_list_path):
		raise FileNotFoundError(f"Reference list not found: {ref_list_path}")

	df = pd.read_csv(ref_list_path, sep="\t", header=None, dtype=str, keep_default_na=False).fillna("")
	if df.empty:
		return pd.DataFrame(columns=["primary_accession", "accession_type", "segment", "genotype", "subtype"])

	if _looks_like_header_row(df.iloc[0].tolist()):
		header = [_normalize_reflist_text(v).lower() for v in df.iloc[0].tolist()]
		columns = _deduplicate_column_names([REFLIST_COLUMN_ALIASES.get(col, col or f"col{i}") for i, col in enumerate(header)])
		df = df.iloc[1:].reset_index(drop=True)
		df.columns = columns
	else:
		df.columns = _default_reflist_column_names(df.shape[1])

	for col in df.columns:
		df[col] = df[col].map(_normalize_reflist_text)

	if "primary_accession" not in df.columns and len(df.columns) > 0:
		df = df.rename(columns={df.columns[0]: "primary_accession"})

	for required in ["primary_accession", "accession_type", "segment", "genotype", "subtype"]:
		if required not in df.columns:
			df[required] = ""

	df = df[df["primary_accession"] != ""].reset_index(drop=True)
	return df


def load_reference_file_dict(ref_list_path):
	refs = load_reference_file_table(ref_list_path)
	refs = refs[refs["accession_type"].fillna("").astype(str).str.strip() != ""]
	return {
		str(row["primary_accession"]).strip(): str(row["accession_type"]).strip()
		for _, row in refs.iterrows()
		if str(row["primary_accession"]).strip()
	}


def load_reference_file_accessions(ref_list_path, allowed_types=None):
	refs = load_reference_file_table(ref_list_path)
	if allowed_types is not None:
		allowed = {str(v).strip().lower() for v in allowed_types}
		refs = refs[refs["accession_type"].fillna("").astype(str).str.lower().isin(allowed)]
	return refs["primary_accession"].astype(str).str.strip().drop_duplicates().tolist()


def load_master_accessions_from_file(ref_list_path):
	refs = load_reference_file_table(ref_list_path)
	if refs.empty:
		return []
	masters = refs[refs["accession_type"].fillna("").astype(str).str.lower() == "master"]
	if not masters.empty:
		return masters["primary_accession"].astype(str).str.strip().drop_duplicates().tolist()
	non_excluded = refs[refs["accession_type"].fillna("").astype(str).str.lower() != "exclusion_list"]
	if not non_excluded.empty:
		return non_excluded["primary_accession"].astype(str).str.strip().drop_duplicates().tolist()
	return refs["primary_accession"].astype(str).str.strip().drop_duplicates().tolist()


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