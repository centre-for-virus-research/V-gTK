#!/usr/bin/env python3
import os
import sqlite3
import pandas as pd
from argparse import ArgumentParser


def _pick_col(df, candidates):
	for c in candidates:
		if c in df.columns:
			return c
	return None


def _segment_norm(value):
	if value is None:
		return "0"
	s = str(value).strip()
	if not s:
		return "0"
	digits = ''.join([x for x in s if x.isdigit()])
	return digits if digits else s


def main(args):
	if not os.path.isfile(args.db):
		raise FileNotFoundError(f"DB not found: {args.db}")

	os.makedirs(args.output_dir, exist_ok=True)
	ref_dir = os.path.join(args.output_dir, "ref_backbones")
	tree_dir = os.path.join(args.output_dir, "trees")
	existing_ids_dir = os.path.join(args.output_dir, "existing_ids")
	os.makedirs(ref_dir, exist_ok=True)
	os.makedirs(tree_dir, exist_ok=True)
	os.makedirs(existing_ids_dir, exist_ok=True)

	conn = sqlite3.connect(args.db)
	try:
		df_meta = pd.read_sql_query("SELECT * FROM meta_data", conn)
		df_aln = pd.read_sql_query("SELECT * FROM sequence_alignment", conn)
		df_trees = pd.read_sql_query("SELECT * FROM trees", conn) if _table_exists(conn, "trees") else pd.DataFrame()
	finally:
		conn.close()

	if "primary_accession" not in df_meta.columns:
		raise ValueError("meta_data must contain primary_accession")
	if "accession_type" not in df_meta.columns:
		df_meta["accession_type"] = "query"
	if "segment" not in df_meta.columns:
		df_meta["segment"] = "0"

	df_meta["segment"] = df_meta["segment"].apply(_segment_norm)
	df_meta["accession_type"] = df_meta["accession_type"].fillna("query").astype(str)

	# exported reference list from DB
	df_refs = df_meta[df_meta["accession_type"].str.lower().isin(["master", "reference"])][["primary_accession", "accession_type", "segment"]].drop_duplicates()
	ref_list_path = os.path.join(args.output_dir, "ref_list_from_db.tsv")
	df_refs.to_csv(ref_list_path, sep="\t", index=False, header=False)

	# existing IDs by segment for update placement exclusion (place new tips only)
	excluded_ids = set()
	conn_excl = sqlite3.connect(args.db)
	try:
		meta_cols = {row[1] for row in conn_excl.execute("PRAGMA table_info(meta_data)").fetchall()}
		if "exclusion_status" in meta_cols:
			df_excl_meta = pd.read_sql_query(
				"SELECT primary_accession FROM meta_data WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) <> '' "
				"AND LOWER(TRIM(COALESCE(CAST(exclusion_status AS TEXT), ''))) NOT IN ('', '0', 'false', 'no', 'na', 'none', 'nan')",
				conn_excl,
			)
			excluded_ids.update(df_excl_meta["primary_accession"].fillna("").astype(str).str.strip().tolist())
		if _table_exists(conn_excl, "excluded_accessions"):
			df_excl = pd.read_sql_query("SELECT primary_accession FROM excluded_accessions", conn_excl)
			excluded_ids.update(df_excl["primary_accession"].fillna("").astype(str).str.strip().tolist())
	finally:
		conn_excl.close()

	for segment, part in df_meta.groupby("segment"):
		accessions = (
			part["primary_accession"]
			.fillna("")
			.astype(str)
			.str.strip()
		)
		accessions = sorted([x for x in accessions if x and x not in excluded_ids])
		ids_path = os.path.join(existing_ids_dir, f"segment_{segment}_ids.txt")
		with open(ids_path, "w", encoding="utf-8") as handle:
			for acc in accessions:
				handle.write(acc + "\n")

	# join with sequence alignment to export per-segment backbones
	if "primary_accession" not in df_aln.columns and "sequence_id" in df_aln.columns:
		df_aln["primary_accession"] = df_aln["sequence_id"]
	if "alignment" not in df_aln.columns:
		align_col = _pick_col(df_aln, ["aligned_seq", "sequence", "aln", "alignment_seq"])
		if not align_col:
			raise ValueError("sequence_alignment must contain alignment-like column")
		df_aln["alignment"] = df_aln[align_col]
	if "alignment_name" not in df_aln.columns:
		df_aln["alignment_name"] = "0"

	df_join = df_aln.merge(
		df_meta[["primary_accession", "accession_type", "segment"]],
		on="primary_accession",
		how="left",
		suffixes=("_aln", "_meta"),
	)
	if "segment" not in df_join.columns:
		if "segment_meta" in df_join.columns:
			df_join["segment"] = df_join["segment_meta"]
		elif "segment_aln" in df_join.columns:
			df_join["segment"] = df_join["segment_aln"]
		elif "segment_y" in df_join.columns:
			df_join["segment"] = df_join["segment_y"]
		elif "segment_x" in df_join.columns:
			df_join["segment"] = df_join["segment_x"]
		else:
			df_join["segment"] = "0"
	df_join["segment"] = df_join["segment"].apply(_segment_norm)
	df_join = df_join[df_join["accession_type"].fillna("").str.lower().isin(["master", "reference"])]

	manifest_rows = []
	for segment, part in df_join.groupby("segment"):
		rows = []
		for _, r in part.drop_duplicates(subset=["primary_accession", "alignment_name"]).iterrows():
			acc = str(r["primary_accession"]).strip()
			seq = str(r["alignment"]).strip()
			if not acc or not seq:
				continue
			rows.append((acc, seq))
		if not rows:
			continue
		out_fa = os.path.join(ref_dir, f"refset_{segment}_aln.fasta")
		with open(out_fa, "w", encoding="utf-8") as handle:
			for acc, seq in rows:
				handle.write(f">{acc}\n{seq}\n")
		manifest_rows.append({"segment": segment, "path": out_fa})

	pd.DataFrame(manifest_rows).to_csv(os.path.join(args.output_dir, "ref_backbone_manifest.tsv"), sep="\t", index=False)

	# export preferred trees by segment (usher > iqtree > veryfasttree)
	tree_manifest = []
	if not df_trees.empty and "newick" in df_trees.columns:
		df_trees["source"] = df_trees.get("source", "unknown")
		if "segment" not in df_trees.columns:
			df_trees["segment"] = "0"
		df_trees["segment"] = df_trees["segment"].apply(_segment_norm)
		priority = {"usher": 0, "iqtree": 1, "veryfasttree": 2}
		df_trees["p"] = df_trees["source"].map(priority).fillna(9)
		for segment, part in df_trees.sort_values(["p"]).groupby("segment", sort=False):
			row = part.iloc[0]
			nwk = str(row.get("newick", "")).strip()
			if not nwk:
				continue
			out_nwk = os.path.join(tree_dir, f"segment_{segment}.nwk")
			with open(out_nwk, "w", encoding="utf-8") as handle:
				handle.write(nwk + "\n")
			tree_manifest.append({"segment": segment, "source": row.get("source", "unknown"), "path": out_nwk})
	tree_manifest_path = os.path.join(args.output_dir, "tree_manifest.tsv")
	pd.DataFrame(tree_manifest, columns=["segment", "source", "path"]).to_csv(tree_manifest_path, sep="\t", index=False)

	print(f"Exported update assets under: {args.output_dir}")


def _table_exists(conn, table):
	row = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)).fetchone()
	return row is not None


if __name__ == "__main__":
	parser = ArgumentParser(description="Export DB-backed update assets (ref list, backbones, trees)")
	parser.add_argument("--db", required=True, help="Existing SQLite DB path")
	parser.add_argument("--output_dir", required=True, help="Output directory")
	args = parser.parse_args()
	main(args)
