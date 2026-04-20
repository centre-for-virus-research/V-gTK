#!/usr/bin/env python3
import os
import re
import csv
import shutil
import sqlite3
import sys
import pandas as pd
from Bio import SeqIO
from os.path import join, normpath
from datetime import datetime
from argparse import ArgumentParser


class CreateSqliteDB:
	def __init__(
		self,
		meta_data,
		features,
		pad_aln,
		gene_info,
		m49_countries,
		m49_interm_region,
		m49_regions,
		m49_sub_regions,
		proj_settings,
		fasta_sequence_file,
		insertions,
		host_taxa_file,
		base_dir,
		output_dir,
		db_name,
		db_status,
		tree_file=None,
		iqtree_file=None,
		usher_tree=None,
		cluster_tsv=None,
		cluster_min_seq_id=None,
		filtered_ids_file=None,
		filtered_details_file=None,
		tree_manifest=None,
		update=False,
		update_db=None,
		batch_id=None,
	):
		self.meta_data = meta_data
		self.features = features
		self.pad_aln = pad_aln
		self.gene_info = gene_info
		self.m49_countries = m49_countries
		self.m49_interm_region = m49_interm_region
		self.m49_regions = m49_regions
		self.m49_sub_regions = m49_sub_regions
		self.proj_settings = proj_settings
		self.fasta_sequence_file = fasta_sequence_file
		self.insertions = insertions
		self.host_taxa_file = host_taxa_file
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.db_name = db_name
		self.db_status = db_status
		self.tree_file = tree_file
		self.iqtree_file = iqtree_file
		self.usher_tree = usher_tree
		self.cluster_tsv = cluster_tsv
		self.cluster_min_seq_id = cluster_min_seq_id
		self.filtered_ids_file = filtered_ids_file
		self.filtered_details_file = filtered_details_file
		self.tree_manifest = tree_manifest
		self.update = bool(update)
		self.update_db = update_db
		self.batch_id = batch_id or datetime.now().strftime("batch_%Y%m%d_%H%M%S")

	@staticmethod
	def _read_tree_file(tree_path):
		if not tree_path:
			return None
		try:
			with open(tree_path, "r", encoding="utf-8") as handle:
				return handle.read().strip()
		except FileNotFoundError:
			return None

	@staticmethod
	def _load_tree_manifest(manifest_path):
		if not manifest_path or not os.path.isfile(manifest_path):
			return []
		rows = []
		with open(manifest_path, "r", encoding="utf-8") as f:
			reader = csv.DictReader(f, delimiter='\t')
			for row in reader:
				path_val = (row.get("path") or "").strip()
				if not path_val:
					continue
				rows.append({
					"source": (row.get("source") or "").strip() or "unknown",
					"name": (row.get("name") or "").strip() or None,
					"segment_key": (row.get("segment_key") or "").strip() or None,
					"path": path_val,
				})
		return rows

	@staticmethod
	def _segment_from_key(segment_key):
		if not segment_key:
			return None
		key = str(segment_key).strip()
		if not key:
			return None
		if key.isdigit():
			return key

		patterns = [
			r"(?:^|[_-])segment[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])seg[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])refset[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])(\d+)(?:$|[_-])",
		]
		for pat in patterns:
			m = re.search(pat, key, flags=re.IGNORECASE)
			if m:
				return m.group(1)
		return None

	def _load_filtered_ids(self):
		if not self.filtered_ids_file:
			return set()
		try:
			with open(self.filtered_ids_file, "r", encoding="utf-8") as f:
				return {line.strip() for line in f if line.strip()}
		except FileNotFoundError:
			return set()

	@staticmethod
	def _require_file(path, label):
		if not path or not os.path.isfile(path):
			raise FileNotFoundError(f"{label} file not found: {path}")

	@staticmethod
	def _read_tsv_required(path, required_columns, label, dtype=None):
		df = pd.read_csv(path, sep="\t", dtype=dtype)
		missing = [c for c in required_columns if c not in df.columns]
		if missing:
			raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
		return df

	@staticmethod
	def _read_csv_required(path, required_columns, label, dtype=None):
		df = pd.read_csv(path, sep=",", dtype=dtype)
		missing = [c for c in required_columns if c not in df.columns]
		if missing:
			raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
		return df

	@staticmethod
	def _ensure_primary_accession(df, label, aliases=None):
		aliases = aliases or []
		if "primary_accession" in df.columns:
			return df
		for alias in aliases:
			if alias in df.columns:
				df["primary_accession"] = df[alias]
				return df
		raise ValueError(f"{label} is missing required columns: primary_accession")

	@staticmethod
	def _normalize_alignment_columns(df, label):
		if "primary_accession" not in df.columns and "sequence_id" in df.columns:
			df["primary_accession"] = df["sequence_id"]
		elif "primary_accession" in df.columns and "sequence_id" not in df.columns:
			df["sequence_id"] = df["primary_accession"]
		else:
			missing = [c for c in ["primary_accession", "sequence_id"] if c not in df.columns]
			if missing:
				raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
		return df

	def _load_filtered_details(self):
		reasons = {}
		if not self.filtered_details_file or not os.path.isfile(self.filtered_details_file):
			return reasons
		try:
			df = pd.read_csv(self.filtered_details_file, sep="\t", dtype=str).fillna("")
		except Exception:
			return reasons
		if "seq_name" not in df.columns:
			return reasons
		for _, row in df.iterrows():
			seq = str(row.get("seq_name", "")).strip()
			if not seq:
				continue
			err = str(row.get("error", "")).strip()
			ref = str(row.get("reference", "")).strip()
			if err:
				reason = f"alignment_filtering: {err}"
			elif ref:
				reason = f"alignment_filtering: filtered against reference {ref}"
			else:
				reason = "alignment_filtering"
			reasons[seq] = reason
		return reasons

	def _add_cluster_column(self, df_meta_data):
		if not self.cluster_tsv:
			return df_meta_data
		try:
			cluster_df = pd.read_csv(self.cluster_tsv, sep="\t", header=None, dtype=str)
		except FileNotFoundError:
			return df_meta_data
		if cluster_df.shape[1] < 2:
			return df_meta_data
		cluster_df = cluster_df.iloc[:, :2]
		cluster_df.columns = ["cluster_rep", "member"]
		cluster_map = dict(zip(cluster_df["member"], cluster_df["cluster_rep"]))
		try:
			min_id = float(self.cluster_min_seq_id) if self.cluster_min_seq_id is not None else None
		except (TypeError, ValueError):
			min_id = None
		col_name = f"cluster_{int(round(min_id * 100))}pct" if min_id is not None else "cluster"
		if "primary_accession" in df_meta_data.columns:
			df_meta_data[col_name] = df_meta_data["primary_accession"].map(cluster_map)
		return df_meta_data

	def load_fasta(self):
		fasta_data = []
		for record in SeqIO.parse(self.fasta_sequence_file, "fasta"):
			fasta_data.append({"header": record.id, "sequence": str(record.seq)})
		return pd.DataFrame(fasta_data)

	@staticmethod
	def _normalize_db_status(db_status):
		s = (db_status or "").strip().lower()
		if s in {"new", "new db", "create", "created", "fresh"}:
			return "new db"
		if s in {"modified", "update", "updated", "changed"}:
			return "last updated"
		return db_status

	def _db_path(self):
		if self.update and not self.update_db:
			raise ValueError("--update requires --update_db path")
		return join(self.base_dir, self.output_dir, self.db_name + ".db")

	@staticmethod
	def _paths_equivalent(path_a, path_b):
		if not path_a or not path_b:
			return False
		return os.path.realpath(os.path.abspath(path_a)) == os.path.realpath(os.path.abspath(path_b))

	def _prepare_update_target_db(self, db_path):
		if not self.update:
			return
		if not self.update_db:
			raise ValueError("--update requires --update_db path")
		if not os.path.isfile(self.update_db):
			raise FileNotFoundError(f"update_db file not found: {self.update_db}")
		if self._paths_equivalent(db_path, self.update_db):
			return
		os.makedirs(os.path.dirname(db_path), exist_ok=True)
		shutil.copyfile(self.update_db, db_path)

	@staticmethod
	def _table_exists(conn, table):
		row = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)).fetchone()
		return row is not None

	@staticmethod
	def _table_columns(conn, table):
		if not CreateSqliteDB._table_exists(conn, table):
			return []
		rows = conn.execute(f"PRAGMA table_info({table})").fetchall()
		return [r[1] for r in rows]

	@staticmethod
	def _normalize_key_series(s):
		return s.fillna("").astype(str).str.strip()

	@staticmethod
	def _cluster_placeholder_columns(columns):
		return [c for c in columns if re.fullmatch(r"cluster_\d+pct", str(c or "").strip())]

	def _fetch_existing_keys(self, conn, table, key_cols):
		if not self._table_exists(conn, table):
			return set()
		df = pd.read_sql_query(f"SELECT {', '.join(key_cols)} FROM {table}", conn)
		for c in key_cols:
			if c not in df.columns:
				raise ValueError(f"DB table '{table}' is missing expected key column '{c}'")
			df[c] = self._normalize_key_series(df[c])
		if len(key_cols) == 1:
			return set(df[key_cols[0]].tolist())
		return set(map(tuple, df[key_cols].itertuples(index=False, name=None)))

	@staticmethod
	def _dedupe_incoming_df(df, key_cols):
		before = len(df)
		df2 = df.drop_duplicates(subset=key_cols, keep="first")
		return df2, (before - len(df2))

	def _align_df_to_existing_schema(self, conn, table, df):
		"""Align incoming dataframe to an existing SQLite table schema.

		- Drops incoming columns absent from existing table (warn)
		- Fails if any existing table columns are missing in incoming dataframe
		- Reorders to existing table column order
		"""
		existing_cols = self._table_columns(conn, table)
		if not existing_cols:
			return df

		incoming_cols = list(df.columns)
		extra_cols = [c for c in incoming_cols if c not in existing_cols]
		if extra_cols:
			print(
				f"[CreateSqliteDB][warn] Incoming '{table}' has columns not present in existing DB schema; "
				f"dropping: {extra_cols}"
			)
			df = df.drop(columns=extra_cols)

		missing_cols = [c for c in existing_cols if c not in df.columns]
		if self.update and table == "meta_data":
			for cluster_col in self._cluster_placeholder_columns(missing_cols):
				df[cluster_col] = "NA- see tree"
			missing_cols = [c for c in existing_cols if c not in df.columns]
		if missing_cols:
			raise ValueError(
				f"Incoming '{table}' dataframe is missing columns required by existing DB schema: {missing_cols}"
			)

		return df[existing_cols]

	def _resolve_key_cols_for_existing_schema(self, conn, table, key_cols, incoming_cols):
		existing_cols = set(self._table_columns(conn, table))
		if not existing_cols:
			return key_cols

		resolved = [c for c in key_cols if c in existing_cols and c in incoming_cols]
		if resolved != key_cols:
			print(
				f"[CreateSqliteDB][warn] Adjusted key columns for '{table}' to match existing DB schema: "
				f"requested={key_cols}, resolved={resolved}"
			)

		if resolved:
			return resolved

		fallback_candidates = [
			"primary_accession", "accession", "header", "name", "m49_code", "id"
		]
		for c in fallback_candidates:
			if c in existing_cols and c in incoming_cols:
				print(
					f"[CreateSqliteDB][warn] Falling back to key column '{c}' for table '{table}' "
					"due to schema mismatch."
				)
				return [c]

		raise ValueError(
			f"Cannot resolve compatible key columns for table '{table}'. "
			f"Requested keys={key_cols}; existing columns={sorted(existing_cols)}; "
			f"incoming columns={sorted(list(incoming_cols))}"
		)

	def _infer_key_cols(self, table, df):
		cols = list(df.columns)
		if table == "meta_data":
			return ["primary_accession", "segment"] if "segment" in cols else ["primary_accession"]
		if table == "sequences":
			return ["header"]
		if table == "sequence_alignment":
			key = ["primary_accession"]
			if "alignment_name" in cols:
				key.append("alignment_name")
			if "segment" in cols:
				key.append("segment")
			return key
		if table == "features":
			key = ["accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end", "cds_start", "cds_end", "product"]
			if "segment" in cols:
				key.append("segment")
			missing = [c for c in key if c not in cols]
			if not missing:
				return key
			if "primary_accession" in cols and "segment" in cols:
				return ["primary_accession", "segment"]
			if "primary_accession" in cols:
				return ["primary_accession"]
			if "accession" in cols and "segment" in cols:
				return ["accession", "segment"]
			if "accession" in cols:
				return ["accession"]
			return cols[:1]
		if table == "insertions":
			key = [c for c in ["primary_accession", "reference", "insertion"] if c in cols]
			if "segment" in cols:
				key.append("segment")
			return key if key else cols[:1]
		if table == "host_taxa":
			for c in ["taxonomy_id", "host_taxa_id", "tax_id", "id"]:
				if c in cols:
					return [c]
			return [cols[0]] if cols else []
		if table == "genes":
			for c in ["name", "gene_name", "description"]:
				if c in cols:
					return [c]
			return cols[:1]
		if table in ["m49_country", "m49_intermediate", "m49_regions", "m49_sub_regions", "project_settings"]:
			return cols[:1] if cols else []
		return cols[:1] if cols else []

	def merge_table_append_nonredundant(self, conn, df, table, key_cols=None, update_exclusions=None):
		if df is None:
			return 0
		df = df.copy()
		key_cols = key_cols if key_cols is not None else self._infer_key_cols(table, df)
		if not key_cols:
			if not self.update:
				df.to_sql(table, conn, if_exists="replace", index=False)
			else:
				df.to_sql(table, conn, if_exists="append" if self._table_exists(conn, table) else "replace", index=False)
			return len(df)

		if self.update and self._table_exists(conn, table):
			df = self._align_df_to_existing_schema(conn, table, df)
			key_cols = self._resolve_key_cols_for_existing_schema(conn, table, key_cols, set(df.columns))

		for c in key_cols:
			if c not in df.columns:
				raise ValueError(f"Incoming '{table}' dataframe missing key column '{c}'")
			df[c] = self._normalize_key_series(df[c])
		df, dropped_internal = self._dedupe_incoming_df(df, key_cols)
		if dropped_internal:
			print(f"[CreateSqliteDB] Dropped {dropped_internal} duplicate incoming rows in '{table}' by key {key_cols}")
		if not self.update:
			df.to_sql(table, conn, if_exists="replace", index=False)
			return len(df)
		if not self._table_exists(conn, table):
			df.to_sql(table, conn, if_exists="replace", index=False)
			print(f"[CreateSqliteDB] Created table '{table}' with {len(df)} rows (table did not exist)")
			return len(df)
		existing_keys = self._fetch_existing_keys(conn, table, key_cols)
		upsert_tables = {"meta_data", "sequence_alignment", "features", "insertions", "sequences"}
		if table in upsert_tables:
			if len(key_cols) == 1:
				key = key_cols[0]
				keys_to_replace = sorted(set(df[key].tolist()))
				if keys_to_replace:
					placeholders = ",".join(["?"] * len(keys_to_replace))
					conn.execute(f"DELETE FROM {table} WHERE {key} IN ({placeholders})", keys_to_replace)
			else:
				where_clause = " AND ".join([f"{c}=?" for c in key_cols])
				for key_tuple in set(map(tuple, df[key_cols].itertuples(index=False, name=None))):
					conn.execute(f"DELETE FROM {table} WHERE {where_clause}", tuple(key_tuple))
			df.to_sql(table, conn, if_exists="append", index=False)
			print(f"[CreateSqliteDB] Upserted {len(df)} rows into '{table}' by key {key_cols}")
			return len(df)

		if len(key_cols) == 1:
			key = key_cols[0]
			new_mask = ~df[key].isin(existing_keys)
			dup_keys = df.loc[~new_mask, key].tolist()
		else:
			incoming_keys = list(map(tuple, df[key_cols].itertuples(index=False, name=None)))
			keep_flags = [k not in existing_keys for k in incoming_keys]
			new_mask = pd.Series(keep_flags, index=df.index)
			dup_keys = [k for k, keep in zip(incoming_keys, keep_flags) if not keep]
		df_new = df.loc[new_mask].copy()
		if not df_new.empty:
			df_new.to_sql(table, conn, if_exists="append", index=False)
			print(f"[CreateSqliteDB] Appended {len(df_new)} new rows into '{table}' (non-redundant)")
		else:
			print(f"[CreateSqliteDB] No new rows to append into '{table}' (all duplicates)")
		if update_exclusions is not None and dup_keys:
			now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
			for k in dup_keys:
				update_exclusions.append({"batch_id": self.batch_id, "table_name": table, "key": str(k), "reason": "duplicate_key_in_db", "date": now_str})
		return len(df_new)

	@staticmethod
	def _table_row_count(conn, table):
		if not CreateSqliteDB._table_exists(conn, table):
			return 0
		row = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()
		return int(row[0]) if row else 0

	def create_db(self):
		self._require_file(self.meta_data, "meta_data")
		self._require_file(self.features, "features")
		self._require_file(self.pad_aln, "pad_aln")
		self._require_file(self.gene_info, "gene_info")
		self._require_file(self.m49_countries, "m49_countries")
		self._require_file(self.m49_interm_region, "m49_interm_region")
		self._require_file(self.m49_regions, "m49_regions")
		self._require_file(self.m49_sub_regions, "m49_sub_regions")
		self._require_file(self.proj_settings, "proj_settings")
		self._require_file(self.fasta_sequence_file, "fasta_sequences")
		self._require_file(self.insertions, "insertions")
		self._require_file(self.host_taxa_file, "host_taxa_file")

		excluded_records = []
		update_exclusions = []
		filtered_ids = self._load_filtered_ids()
		filtered_details = self._load_filtered_details()
		if filtered_ids:
			print(f"[CreateSqliteDB] Excluding {len(filtered_ids)} filtered sequences from DB")
			for fid in filtered_ids:
				excluded_records.append({"primary_accession": fid, "reason": filtered_details.get(fid, "alignment_filtering")})

		df_meta_data = self._read_tsv_required(self.meta_data, ["primary_accession"], "meta_data", dtype=str)
		acc_type_col = "accession_type" if "accession_type" in df_meta_data.columns else None
		if acc_type_col:
			acc_type_norm = df_meta_data[acc_type_col].fillna("").str.strip().str.lower()
			is_ref_or_master = acc_type_norm.isin(["reference", "master"])
		else:
			is_ref_or_master = pd.Series(False, index=df_meta_data.index)
		if filtered_ids and "primary_accession" in df_meta_data.columns:
			before_count = len(df_meta_data)
			remove_mask = df_meta_data["primary_accession"].isin(filtered_ids) & (~is_ref_or_master)
			df_meta_data = df_meta_data[~remove_mask]
			after_count = len(df_meta_data)
			if before_count != after_count:
				print(f"[CreateSqliteDB] Removed {before_count - after_count} filtered non-reference sequences from meta_data")
			if acc_type_col:
				ref_or_master_ids = set(
					df_meta_data[
						df_meta_data[acc_type_col].fillna("").str.strip().str.lower().isin(["master", "reference"])
					]["primary_accession"].astype(str).str.strip().tolist()
				)
				excluded_records = [
					r for r in excluded_records
					if str(r.get("primary_accession", "")).strip() not in ref_or_master_ids
				]
		is_ref_or_master = is_ref_or_master.reindex(df_meta_data.index, fill_value=False)
		exclusion_reason = pd.Series("", index=df_meta_data.index, dtype=str)
		for col in ["exclusion", "exclusion_criteria"]:
			if col in df_meta_data.columns:
				col_values = df_meta_data[col].fillna("").astype(str).str.strip()
				exclusion_reason = exclusion_reason.where(exclusion_reason != "", col_values)

		status_mask = pd.Series(False, index=df_meta_data.index)
		if "exclusion_status" in df_meta_data.columns:
			status_values = df_meta_data["exclusion_status"].fillna("").astype(str).str.strip()
			status_mask = ~status_values.str.lower().isin(["", "0", "false", "no", "na", "none", "nan"])
			fallback_reason = "metadata_exclusion"
			exclusion_reason = exclusion_reason.mask((exclusion_reason == "") & status_mask, fallback_reason)

		exclusion_mask = exclusion_reason.astype(str).str.strip() != ""
		excluded_rows = df_meta_data[exclusion_mask]
		if not excluded_rows.empty:
			print(f"[CreateSqliteDB] Found {len(excluded_rows)} rows with exclusions in meta_data")
			for idx, row in excluded_rows.iterrows():
				acc = str(row.get("primary_accession", "")).strip()
				reason = str(exclusion_reason.loc[idx]).strip()
				if acc:
					excluded_records.append({"primary_accession": acc, "reason": reason})
			remove_mask = exclusion_mask & (~is_ref_or_master)
			df_meta_data = df_meta_data[~remove_mask]
			retained_refs = (exclusion_mask & is_ref_or_master).sum()
			if retained_refs:
				print(f"[CreateSqliteDB] Retained {retained_refs} reference/master rows despite exclusion flags")

		df_meta_data = self._add_cluster_column(df_meta_data)
		valid_accessions = set(df_meta_data["primary_accession"].dropna().astype(str).str.strip()) if "primary_accession" in df_meta_data.columns else set()

		df_features = self._read_tsv_required(self.features, [], "features")
		if valid_accessions and "accession" in df_features.columns:
			df_features = df_features[df_features["accession"].astype(str).str.strip().isin(valid_accessions)]
		elif valid_accessions and "primary_accession" in df_features.columns:
			df_features = df_features[df_features["primary_accession"].astype(str).str.strip().isin(valid_accessions)]

		df_aln = self._read_tsv_required(self.pad_aln, [], "pad_aln")
		df_aln = self._normalize_alignment_columns(df_aln, "pad_aln")
		if valid_accessions and "sequence_id" in df_aln.columns:
			df_aln = df_aln[df_aln["sequence_id"].astype(str).str.strip().isin(valid_accessions)]
		elif valid_accessions and "primary_accession" in df_aln.columns:
			df_aln = df_aln[df_aln["primary_accession"].astype(str).str.strip().isin(valid_accessions)]

		df_gene = self._read_tsv_required(self.gene_info, [], "gene_info")
		df_m49_country = self._read_csv_required(self.m49_countries, ["m49_code"], "m49_countries", dtype={"m49_code": str})
		df_m49_interm = self._read_csv_required(self.m49_interm_region, [], "m49_interm_region")
		df_m49_region = self._read_csv_required(self.m49_regions, [], "m49_regions")
		df_m49_sub_region = self._read_csv_required(self.m49_sub_regions, [], "m49_sub_regions")
		df_proj_setting = self._read_tsv_required(self.proj_settings, [], "proj_settings")
		
		df_insertions = self._read_tsv_required(self.insertions, [], "insertions")
		df_insertions = self._ensure_primary_accession(df_insertions, "insertions", aliases=["accession", "sequence_id"])
		if valid_accessions and "primary_accession" in df_insertions.columns:
			df_insertions = df_insertions[df_insertions["primary_accession"].astype(str).str.strip().isin(valid_accessions)]
		
		df_host_taxa = self._read_tsv_required(self.host_taxa_file, [], "host_taxa_file", dtype=str)
		host_meta_col = next((c for c in ["host_taxa_id", "taxonomy_id", "host_tax_id"] if c in df_meta_data.columns), None)
		if host_meta_col and not df_host_taxa.empty:
			valid_taxa = set(df_meta_data[host_meta_col].dropna().astype(str).str.strip())
			taxa_col = next((c for c in ["taxa_id", "taxonomy_id", "host_taxa_id", "tax_id", "id"] if c in df_host_taxa.columns), None)
			if taxa_col:
				df_host_taxa = df_host_taxa[df_host_taxa[taxa_col].astype(str).str.strip().isin(valid_taxa)]
		
		df_fasta_sequences = self.load_fasta()
		if valid_accessions and "header" in df_fasta_sequences.columns:
			df_fasta_sequences = df_fasta_sequences[df_fasta_sequences["header"].astype(str).str.strip().isin(valid_accessions)]

		db_path = self._db_path()
		os.makedirs(os.path.dirname(db_path), exist_ok=True)
		self._prepare_update_target_db(db_path)
		conn = sqlite3.connect(db_path)
		cursor = conn.cursor()
		cursor.execute("PRAGMA foreign_keys = ON;")
		cursor.execute("CREATE TABLE IF NOT EXISTS trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS info (creation_type TEXT, date TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS excluded_accessions (primary_accession TEXT, reason TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS update_exclusions (batch_id TEXT, table_name TEXT, key TEXT, reason TEXT, date TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS update_batches (batch_id TEXT PRIMARY KEY, started_at TEXT, finished_at TEXT, update_db TEXT, mode TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS update_table_deltas (batch_id TEXT, table_name TEXT, before_count INTEGER, after_count INTEGER, delta INTEGER);")

		now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		cursor.execute("INSERT OR REPLACE INTO update_batches (batch_id, started_at, finished_at, update_db, mode) VALUES (?, ?, ?, ?, ?)", (self.batch_id, now_str, None, self.update_db if self.update else None, "update" if self.update else "create"))

		tables_for_delta = ["meta_data", "features", "sequence_alignment", "genes", "sequences", "insertions", "host_taxa"]
		before_counts = {t: self._table_row_count(conn, t) for t in tables_for_delta}

		self.merge_table_append_nonredundant(conn, df_meta_data, "meta_data", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_features, "features", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_aln, "sequence_alignment", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_gene, "genes", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_m49_country, "m49_country", ["m49_code"], update_exclusions)
		self.merge_table_append_nonredundant(conn, df_m49_interm, "m49_intermediate", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_m49_region, "m49_regions", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_m49_sub_region, "m49_sub_regions", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_proj_setting, "project_settings", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_fasta_sequences, "sequences", ["header"], update_exclusions)
		self.merge_table_append_nonredundant(conn, df_insertions, "insertions", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_host_taxa, "host_taxa", None, update_exclusions)

		if excluded_records:
			df_excluded = pd.DataFrame(excluded_records).drop_duplicates(subset=["primary_accession"])
			
			df_excluded['base_reason'] = df_excluded['reason'].apply(lambda x: str(x).split(':')[0].strip() if ':' in str(x) else str(x).strip())
			reason_counts = df_excluded['base_reason'].value_counts()
			
			summary_lines = []
			summary_lines.append("\n" + "="*80)
			summary_lines.append(f"[CreateSqliteDB] DB Construction Summary:")
			summary_lines.append("="*80)
			summary_lines.append(f"Successfully included in DB : {len(df_meta_data)}")
			summary_lines.append(f"Sequences filtered out      :{len(df_excluded)}")
			summary_lines.append("-" * 80)
			summary_lines.append(f"{'Category':<45} | {'Count':<10}")
			summary_lines.append("-" * 80)
			for reason, count in reason_counts.items():
				summary_lines.append(f"{reason:<45} | {count:<10}")
			
			summary_lines.append("\n[CreateSqliteDB] Detailed Exclusion Reasons (Top 10):")
			detailed_counts = df_excluded['reason'].value_counts().head(10)
			for reason, count in detailed_counts.items():
				summary_lines.append(f"  - {count:<4} : {reason}")
			summary_lines.append("="*80 + "\n")
			
			summary_str = "\n".join(summary_lines)
			print(summary_str)
			
			summary_path = os.path.join(self.base_dir, self.output_dir, "db_summary.txt")
			with open(summary_path, "w", encoding="utf-8") as sf:
				sf.write(summary_str)

			if self.update and self._table_exists(conn, "excluded_accessions"):
				df_excluded[['primary_accession', 'reason']].to_sql("excluded_accessions", conn, if_exists="append", index=False)
			else:
				df_excluded[['primary_accession', 'reason']].to_sql("excluded_accessions", conn, if_exists="replace", index=False)
			print(f"[CreateSqliteDB] excluded_accessions now has >= {len(df_excluded)} newly added records")

		if update_exclusions:
			pd.DataFrame(update_exclusions).to_sql("update_exclusions", conn, if_exists="append", index=False)
			print(f"[CreateSqliteDB] Logged {len(update_exclusions)} update-mode duplicate keys into update_exclusions")

		accession_to_segment = {}
		if "primary_accession" in df_meta_data.columns and "segment" in df_meta_data.columns:
			for _, row in df_meta_data[["primary_accession", "segment"]].dropna().iterrows():
				acc = str(row["primary_accession"]).strip()
				seg = str(row["segment"]).strip()
				if acc and seg and acc not in accession_to_segment:
					accession_to_segment[acc] = seg

		tree_records = []
		for name, source, tree_path in [("veryfasttree", "veryfasttree", self.tree_file), ("iqtree", "iqtree", self.iqtree_file), ("usher", "usher", self.usher_tree)]:
			newick = self._read_tree_file(tree_path)
			if newick:
				tree_records.append({"name": name, "source": source, "segment_key": None, "segment": None, "newick": newick})

		for entry in self._load_tree_manifest(self.tree_manifest):
			newick = self._read_tree_file(entry["path"])
			if not newick:
				continue
			seg_key = entry.get("segment_key")
			seg_num = accession_to_segment.get(seg_key) if seg_key else None
			if seg_num is None:
				seg_num = self._segment_from_key(seg_key)
			name = entry.get("name") or f"{entry['source']}_{seg_key or 'tree'}"
			tree_records.append({"name": name, "source": entry["source"], "segment_key": seg_key, "segment": seg_num, "newick": newick})

		if tree_records:
			df_tree = pd.DataFrame(tree_records)
			df_tree["created_at"] = now_str
			if not self.update:
				df_tree.to_sql("trees", conn, if_exists="replace", index=False)
			else:
				if self._table_exists(conn, "trees"):
					for _, tr in df_tree.iterrows():
						conn.execute(
							"DELETE FROM trees WHERE COALESCE(source, '')=? AND COALESCE(name, '')=? AND COALESCE(segment_key, '')=? AND COALESCE(segment, '')=?",
							(
								str(tr.get("source") or "").strip(),
								str(tr.get("name") or "").strip(),
								str(tr.get("segment_key") or "").strip(),
								str(tr.get("segment") or "").strip(),
							),
						)
				df_tree.to_sql("trees", conn, if_exists="append", index=False)

		creation_type = self._normalize_db_status(self.db_status if self.db_status else ("last updated" if self.update else "new db"))
		pd.DataFrame([{"creation_type": creation_type, "date": now_str}]).to_sql("info", conn, if_exists="append", index=False)

		after_counts = {t: self._table_row_count(conn, t) for t in tables_for_delta}
		deltas = []
		for table_name in tables_for_delta:
			before = before_counts.get(table_name, 0)
			after = after_counts.get(table_name, 0)
			deltas.append({"batch_id": self.batch_id, "table_name": table_name, "before_count": int(before), "after_count": int(after), "delta": int(after - before)})
		pd.DataFrame(deltas).to_sql("update_table_deltas", conn, if_exists="append", index=False)
		cursor.execute("UPDATE update_batches SET finished_at=? WHERE batch_id=?", (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.batch_id))

		conn.commit()
		conn.close()


def process(args):
	db_creator = CreateSqliteDB(
		args.meta_data,
		args.features,
		args.pad_aln,
		args.gene_info,
		args.m49_countries,
		args.m49_interm_region,
		args.m49_regions,
		args.m49_sub_regions,
		args.proj_settings,
		args.fasta_sequences,
		args.insertion_file,
		args.host_taxa_file,
		args.base_dir,
		args.output_dir,
		args.db_name,
		args.db_status,
		args.tree_file,
		args.iqtree_file,
		args.usher_tree,
		args.cluster_tsv,
		args.cluster_min_seq_id,
		args.filtered_ids,
		args.filtered_details,
		args.tree_manifest,
		update=args.update,
		update_db=args.update_db,
		batch_id=args.batch_id,
	)
	db_creator.create_db()


if __name__ == "__main__":
	parser = ArgumentParser(description='Creating sqlite DB')
	parser.add_argument('-m', '--meta_data', help='Meta data table', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-b', '--base_dir', help='Base directory', default="tmp")
	parser.add_argument('-o', '--output_dir', help='tmp directory where the database is stored', default="SqliteDB")
	parser.add_argument('-rf', '--features', help='Features table', default="tmp/Tables/features.tsv")
	parser.add_argument('-p', '--pad_aln', help='Padded alignment file', default="tmp/Tables/sequence_alignment.tsv")
	parser.add_argument('-g', '--gene_info', help='Gene table', default="generic/rabv/Tables/gene_info.csv")
	parser.add_argument('-mc', '--m49_countries', help='M49 countries', default="assets/m49_country.csv")
	parser.add_argument('-mir', '--m49_interm_region', help='M49 intermediate regions', default="assets/m49_intermediate_region.csv")
	parser.add_argument('-mr', '--m49_regions', help='M49 regions', default="assets/m49_region.csv")
	parser.add_argument('-msr', '--m49_sub_regions', help='M49 sub-regions', default="assets/m49_sub_region.csv")
	parser.add_argument('-s', '--proj_settings', help='Project settings', default="tmp/Software_info/software_info.tsv")
	parser.add_argument('-fa', '--fasta_sequences', help='Fasta sequences', default="tmp/GenBank-matrix/sequences.fa")
	parser.add_argument('-i', '--insertion_file', help='Nextalign insertion file', default="tmp/Tables/insertions.tsv")
	parser.add_argument('-ht', '--host_taxa_file', help='Host Taxanomy file', default="tmp/HostTaxa/Host_taxa.tsv")
	parser.add_argument('-d', '--db_name', help='Name of the Sqlite database', default="gdb")
	parser.add_argument('-ds', '--db_status', help='Database status: "new db" (default) or "last modified"/"last updated". Determines info.creation_type.', default="new db")
	parser.add_argument('-t', '--tree_file', help='VeryFastTree Newick file', default=None)
	parser.add_argument('-it', '--iqtree_file', help='IQ-TREE Newick file', default=None)
	parser.add_argument('-ut', '--usher_tree', help='UShER output Newick file', default=None)
	parser.add_argument('--tree_manifest', help='TSV manifest with columns: source, name, segment_key, path', default=None)
	parser.add_argument('-ct', '--cluster_tsv', help='MMseqs clustering TSV (rep\tmember)', default=None)
	parser.add_argument('-ci', '--cluster_min_seq_id', help='MMseqs min sequence identity used for clustering', default=None)
	parser.add_argument('-fi', '--filtered_ids', help='File with filtered sequence IDs (one per line) to exclude from DB', default=None)
	parser.add_argument('-fd', '--filtered_details', help='TSV with filtered sequence details (seq_name, reference, error, warnings)', default=None)
	parser.add_argument('--update', action='store_true', help='Append/update into an existing DB instead of replacing tables')
	parser.add_argument('--update_db', default=None, help='Path to existing DB file for update mode')
	parser.add_argument('--batch_id', default=None, help='Optional update batch identifier for audit logging')
	args = parser.parse_args()
	if args.update_db:
		args.update_db = normpath(args.update_db)
	try:
		process(args)
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		sys.exit(2)
