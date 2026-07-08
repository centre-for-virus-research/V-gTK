#!/usr/bin/env python3
import os
import re
import csv
import sys
import time
import shutil
import sqlite3
import read_file
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from datetime import datetime
from os.path import join, normpath
from argparse import ArgumentParser
from collections import defaultdict


"""
Update mode:
  python scripts/CreateSqliteDB.py --meta_data tmp/Update/GenBank-matrix/gB_matrix_raw.tsv \
    --features tmp/Update/Tables/features.tsv --pad_aln tmp/Update/Tables/sequence_alignment.tsv \
    --fasta_sequences tmp/Update/GenBank-matrix/sequences.fa --host_taxa_file tmp/Update/HostTaxa/Host_taxa.tsv --update
"""

 
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
        base_dir,
        output_dir,
        db_name,
        db_status,
        host_taxa_file,
        host_lineage_file,
        host_children_file,
        host_lineage_lookup_file,
        db_file=None,
        tree_file=None,
        iqtree_file=None,
        usher_tree=None,
        tree_dir=None,
        cluster_tsv=None,
        cluster_min_seq_id=None,
        filtered_ids_file=None,
        filtered_details_file=None,
        tree_manifest=None,
        update=False,
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
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.db_file = db_file
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
        self.tree_dir = tree_dir
        self.host_taxa_file = host_taxa_file
        self.host_lineage_file = host_lineage_file
        self.host_children_file = host_children_file
        self.host_lineage_lookup_file = host_lineage_lookup_file
        self.update = bool(update)


    @staticmethod
    def _timestamp_for_backup() -> str:
        return datetime.now().strftime("%Y%m%d_%H%M%S")

    def _backup_existing_db_if_update(self, db_path: str) -> str | None:

        if not self.update:
            return None

        if not os.path.isfile(db_path):
            print(f"[CreateSqliteDB] Update mode enabled but DB does not exist yet, so no backup created: {db_path}")
            return None

        db_dir = os.path.dirname(db_path)
        db_file = os.path.basename(db_path)
        db_stem, db_ext = os.path.splitext(db_file)
        timestamp = self._timestamp_for_backup()

        backup_name = f"{db_stem}_backup_{timestamp}{db_ext or '.db'}"
        backup_path = os.path.join(db_dir, backup_name)

        # extra safety: avoid overwrite if somehow same timestamp/path exists
        counter = 1
        while os.path.exists(backup_path):
            backup_name = f"{db_stem}_backup_{timestamp}_{counter}{db_ext or '.db'}"
            backup_path = os.path.join(db_dir, backup_name)
            counter += 1

        #subprocess.run(["cp", db_path, backup_path], check=True)
        shutil.copy2(db_path, backup_path)
        print(f"[CreateSqliteDB] Backup created: {backup_path}")
        return backup_path

    @staticmethod
    def _coerce_nullable_int_columns(df: pd.DataFrame, cols) -> pd.DataFrame:
        """
        Safely coerce columns to pandas nullable Int64.
        Any blank/whitespace/non-numeric becomes <NA>.
        """
        out = df.copy()
        blankish = {"", "nan", "none", "na", "n/a", "-", "null"}

        for c in cols:
            if c not in out.columns:
                continue

            # normalize to string dtype, strip whitespace
            s = out[c].astype("string").str.strip()

            # blankish -> NA
            s = s.mask(s.isna() | s.str.lower().isin(blankish), pd.NA)

            # numeric conversion (gives float + NaN; never leaves '' behind)
            num = pd.to_numeric(s, errors="coerce")

            # now safe to cast to nullable Int64
            out[c] = num.astype("Int64")

        return out
    

    # ----------------------
    # IO helpers (kept)
    # ----------------------
    @staticmethod
    def _read_delimited(path: str, dtype=None) -> pd.DataFrame:
        return pd.read_csv(path, sep=None, engine="python", dtype=dtype)

    def load_trees_from_dir(self, tree_dir: str) -> pd.DataFrame:
        manifest_path = join(tree_dir, "meta_data.tsv")
        if not os.path.exists(manifest_path):
            raise FileNotFoundError(f"Tree manifest not found: {manifest_path}")

        df_manifest = self._read_delimited(manifest_path, dtype=str)

        required = {"chromosome", "segment_number", "tree_type", "tree_name", "tree_model", "description"}
        missing = required - set(df_manifest.columns)
        if missing:
            raise ValueError(
                f"Tree manifest missing columns: {sorted(missing)}. Found: {list(df_manifest.columns)}"
            )

        rows = []
        for _, r in df_manifest.iterrows():
            # NOTE: your old code references r["tree_path"] but required columns list does not include it.
            # We support either "tree_path" OR "tree_name" as the file name column.
            tree_name = str(r.get("tree_path", "")).strip() or str(r.get("tree_name", "")).strip()
            if not tree_name:
                raise ValueError("Tree manifest must contain 'tree_path' or 'tree_name' with file name")

            tree_path = join(tree_dir, tree_name)
            if not os.path.exists(tree_path):
                raise FileNotFoundError(f"Tree file listed in manifest not found: {tree_path}")

            tree_newick = Path(tree_path).read_text(encoding="utf-8").strip()

            out = dict(r)
            out["tree_path"] = tree_path
            out["newick"] = tree_newick
            rows.append(out)

        return pd.DataFrame(rows)

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
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                path_val = (row.get("path") or "").strip()
                if not path_val:
                    continue
                rows.append(
                    {
                        "source": (row.get("source") or "").strip() or "unknown",
                        "name": (row.get("name") or "").strip() or None,
                        "segment_key": (row.get("segment_key") or "").strip() or None,
                        "path": path_val,
                    }
                )
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

    def _load_filtered_ids(self) -> set:
        """Load the set of sequence IDs that were filtered during alignment."""
        if not self.filtered_ids_file:
            return set()
        try:
            with open(self.filtered_ids_file, "r", encoding="utf-8") as f:
                return {line.strip() for line in f if line.strip()}
        except FileNotFoundError:
            return set()

    @staticmethod
    def _require_file(path: str, label: str):
        if not path or not os.path.isfile(path):
            raise FileNotFoundError(f"{label} file not found: {path}")

    @staticmethod
    def _read_tsv_required(path: str, required_columns, label: str, dtype=None):
        df = pd.read_csv(path, sep="\t", dtype=dtype)
        missing = [c for c in required_columns if c not in df.columns]
        if missing:
            raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
        return df

    @staticmethod
    def _read_csv_required(path: str, required_columns, label: str, dtype=None):
        df = pd.read_csv(path, sep=",", dtype=dtype)
        missing = [c for c in required_columns if c not in df.columns]
        if missing:
            raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
        return df

    @staticmethod
    def _ensure_primary_accession(df: pd.DataFrame, label: str, aliases=None) -> pd.DataFrame:
        if aliases is None:
            aliases = []
        if "primary_accession" in df.columns:
            return df
        for alias in aliases:
            if alias in df.columns:
                df["primary_accession"] = df[alias]
                return df
        raise ValueError(f"{label} is missing required columns: primary_accession")

    @staticmethod
    def _normalize_alignment_columns(df: pd.DataFrame, label: str) -> pd.DataFrame:
        if "primary_accession" not in df.columns and "sequence_id" in df.columns:
            df["primary_accession"] = df["sequence_id"]
        elif "primary_accession" in df.columns and "sequence_id" not in df.columns:
            df["sequence_id"] = df["primary_accession"]
        else:
            missing = [c for c in ["primary_accession", "sequence_id"] if c not in df.columns]
            if missing:
                raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
        return df

    def _load_filtered_details(self) -> dict:
        """Load filtered sequence reasons from filtered_sequences.tsv if available."""
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

        if min_id is not None:
            pct = int(round(min_id * 100))
            col_name = f"cluster_{pct}pct"
        else:
            col_name = "cluster"

        if "primary_accession" in df_meta_data.columns:
            df_meta_data[col_name] = df_meta_data["primary_accession"].map(cluster_map)
        return df_meta_data

    def load_fasta(self):
        fasta_data = []
        for record in SeqIO.parse(self.fasta_sequence_file, "fasta"):
            fasta_data.append({"header": record.id, "sequence": str(record.seq)})
        return pd.DataFrame(fasta_data)

    @staticmethod
    def _normalize_db_status(db_status: str):
        s = (db_status or "").strip().lower()
        if s in {"new", "new db", "create", "created", "fresh"}:
            return "new db"
        if s in {"modified", "update", "updated", "changed"}:
            return "last updated"
        return db_status

    # ----------------------
    # Update-mode merge (non-redundant) helpers
    # ----------------------
    def _db_path(self) -> str:
        
        if self.db_file:
            return normpath(self.db_file)

        if self.update:
            bd = normpath(self.base_dir)
            if bd.endswith(normpath("Update")):
                parent = os.path.dirname(bd)
            else:
                parent = bd
            return join(parent, self.output_dir, self.db_name + ".db")

        return join(self.base_dir, self.output_dir, self.db_name + ".db")

    @staticmethod
    def _table_exists(conn, table: str) -> bool:
        row = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (table,),
        ).fetchone()
        return row is not None

    @staticmethod
    def _table_columns(conn, table: str) -> list[str]:
        try:
            rows = conn.execute(f"PRAGMA table_info({table})").fetchall()
        except Exception:
            return []
        return [r[1] for r in rows]

    def _migrate_legacy_tree_table(self, conn, created_at: str) -> None:
        """
        Older runs accidentally created a table named 'tree'.
        Move compatible rows into the single canonical 'trees' table, then drop 'tree'.
        """
        if not self._table_exists(conn, "tree"):
            return

        legacy_cols = set(self._table_columns(conn, "tree"))
        required = {"newick"}
        if not required.issubset(legacy_cols):
            print("[CreateSqliteDB][WARN] Dropping legacy 'tree' table because it has no newick column")
            conn.execute("DROP TABLE IF EXISTS tree;")
            return

        df_legacy = pd.read_sql_query("SELECT * FROM tree", conn).fillna("")
        migrated = 0

        for _, row in df_legacy.iterrows():
            newick = str(row.get("newick", "")).strip()
            if not newick:
                continue

            name = str(row.get("name", "")).strip()
            if not name:
                name = str(row.get("tree_name", "")).strip()
            if not name:
                name = str(row.get("tree_path", "")).strip()
            if not name:
                name = "tree_dir_tree"

            source = str(row.get("source", "")).strip()
            if not source:
                source = "tree_dir"
            tree_type = str(row.get("tree_type", "")).strip()
            tree_model = str(row.get("tree_model", "")).strip()

            segment_key = str(row.get("segment_key", "")).strip()
            if not segment_key:
                segment_key = str(row.get("chromosome", "")).strip()
            segment = str(row.get("segment", "")).strip()
            if not segment:
                segment = str(row.get("segment_number", "")).strip()

            existing = conn.execute(
                """
                SELECT 1 FROM trees
                WHERE COALESCE(source, '') = ?
                  AND COALESCE(name, '') = ?
                  AND COALESCE(segment_key, '') = ?
                  AND COALESCE(segment, '') = ?
                LIMIT 1;
                """,
                (source, name, segment_key, segment),
            ).fetchone()

            if existing:
                continue

            conn.execute(
                """
                INSERT INTO trees
                    (name, source, tree_type, tree_model, segment_key, segment, newick, created_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    name,
                    source,
                    tree_type or None,
                    tree_model or None,
                    segment_key or None,
                    segment or None,
                    newick,
                    created_at,
                ),
            )
            migrated += 1

        conn.execute("DROP TABLE IF EXISTS tree;")
        if migrated:
            print(f"[CreateSqliteDB] Migrated {migrated} legacy rows from 'tree' into 'trees'")
        print("[CreateSqliteDB] Removed legacy table 'tree'; using only 'trees'")

    def _ensure_trees_table(self, conn) -> None:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS trees (
                name TEXT,
                source TEXT,
                tree_type TEXT,
                tree_model TEXT,
                segment_key TEXT,
                segment TEXT,
                newick TEXT,
                created_at TEXT,
                description TEXT
            );
            """
        )

        existing_cols = set(self._table_columns(conn, "trees"))
        for col, col_type in [
            ("tree_type", "TEXT"),
            ("tree_model", "TEXT"),
            ("description", "TEXT"),
        ]:
            if col not in existing_cols:
                conn.execute(f"ALTER TABLE trees ADD COLUMN {col} {col_type};")

    @staticmethod
    def _normalize_key_series(s: pd.Series) -> pd.Series:
        if pd.api.types.is_string_dtype(s) or s.dtype == object:
            return s.astype("string").fillna("").str.strip()
        # leave numeric/bool/datetime columns unchanged
        return s

    def _fetch_existing_keys(self, conn, table: str, key_cols: list[str]) -> set:
        if not self._table_exists(conn, table):
            return set()

        cols_sql = ", ".join(key_cols)
        df = pd.read_sql_query(f"SELECT {cols_sql} FROM {table}", conn)

        for c in key_cols:
            if c not in df.columns:
                raise ValueError(f"DB table '{table}' is missing expected key column '{c}'")

        for c in key_cols:
            df[c] = self._normalize_key_series(df[c])

        if len(key_cols) == 1:
            return set(df[key_cols[0]].tolist())

        return set(map(tuple, df[key_cols].itertuples(index=False, name=None)))

    @staticmethod
    def _dedupe_incoming_df(df: pd.DataFrame, key_cols: list[str]) -> tuple[pd.DataFrame, int]:
        before = len(df)
        df2 = df.drop_duplicates(subset=key_cols, keep="first")
        return df2, (before - len(df2))

    def _infer_key_cols(self, table: str, df: pd.DataFrame) -> list[str]:
        """
        Conservative key inference.
        If a required key column is missing, we fall back to "first columns" but WARN loudly.
        """
        cols = list(df.columns)

        if table == "meta_data":
            return ["primary_accession"]

        if table == "sequences":
            return ["header"]

        if table == "sequence_alignment":
            if "sequence_id" in cols and "alignment_name" in cols:
                return ["sequence_id", "alignment_name"]
            if "primary_accession" in cols and "alignment_name" in cols:
                return ["primary_accession", "alignment_name"]
            if "sequence_id" in cols:
                return ["sequence_id"]
            # fallback
            print("[CreateSqliteDB][WARN] sequence_alignment key inference fallback to first 2 columns", file=sys.stderr)
            return cols[:2] if len(cols) >= 2 else cols[:1]
        
        if table == "features":
            key = ["accession", "reference_accession", "aln_start", "aln_end", "cds_start", "cds_end", "product"]
            missing = [c for c in key if c not in cols]
            if missing:
                raise ValueError(f"features.tsv is missing required key columns: {missing}")
            return key

        if table == "host_taxa":
            for c in ["host_taxa_id", "tax_id", "id"]:
                if c in cols:
                    return [c]
            print("[CreateSqliteDB][WARN] host_taxa key inference fallback to first column", file=sys.stderr)
            return [cols[0]]

        if table == "host_children":
            for a, b in [("parent_taxa_id", "child_taxa_id"), ("parent", "child")]:
                if a in cols and b in cols:
                    return [a, b]
            print("[CreateSqliteDB][WARN] host_children key inference fallback to first 2 columns", file=sys.stderr)
            return cols[:2] if len(cols) >= 2 else cols[:1]

        if table == "host_lineage":
            key = []
            for c in ["host_taxa_id", "tax_id"]:
                if c in cols:
                    key.append(c)
                    break
            for c in ["rank", "level", "lineage_rank"]:
                if c in cols:
                    key.append(c)
                    break
            if len(key) >= 1:
                return key
            print("[CreateSqliteDB][WARN] host_lineage key inference fallback to first 2 columns", file=sys.stderr)
            return cols[:2] if len(cols) >= 2 else cols[:1]

        if table == "host_lineage_lookup":
            for c in ["host_taxa_id", "tax_id", "query", "name"]:
                if c in cols:
                    return [c]
            print("[CreateSqliteDB][WARN] host_lineage_lookup key inference fallback to first column", file=sys.stderr)
            return [cols[0]]

        if table == "trees":
            # try stable composite if present
            for key in [["tree_path"], ["source", "name"], ["tree_name", "tree_model", "segment_number"]]:
                if all(k in cols for k in key):
                    return key
            return ["tree_path"] if "tree_path" in cols else cols[:2]

        # default fallback
        print(f"[CreateSqliteDB][WARN] {table} key inference fallback to first column", file=sys.stderr)
        return [cols[0]] if cols else []

    def merge_table_append_nonredundant(
        self,
        conn,
        df: pd.DataFrame,
        table: str,
        key_cols: list[str] | None = None,
        update_exclusions: list[dict] | None = None,
    ) -> int:
        """
        Normal mode: replace table with df.
        Update mode: append only rows whose key(s) are not already present in DB.
        Also dedupe within incoming df itself.
        """
        if df is None:
            return 0
        df = df.copy()

        if key_cols is None:
            key_cols = self._infer_key_cols(table, df)

        for c in key_cols:
            if c not in df.columns:
                raise ValueError(f"Incoming '{table}' dataframe missing key column '{c}'")

        # normalize keys
        for c in key_cols:
            df[c] = self._normalize_key_series(df[c])

        # dedupe incoming itself
        df, dropped_internal = self._dedupe_incoming_df(df, key_cols)
        if dropped_internal:
            print(f"[CreateSqliteDB] Dropped {dropped_internal} duplicate incoming rows in '{table}' by key {key_cols}")

        # normal mode: replace
        if not self.update:
            df.to_sql(table, conn, if_exists="replace", index=False)
            return len(df)

        # update mode: create if missing
        if not self._table_exists(conn, table):
            df.to_sql(table, conn, if_exists="replace", index=False)
            print(f"[CreateSqliteDB] Created table '{table}' with {len(df)} rows (table did not exist)")
            return len(df)

        existing_keys = self._fetch_existing_keys(conn, table, key_cols)

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
                update_exclusions.append(
                    {
                        "table_name": table,
                        "key": str(k),
                        "reason": "duplicate_key_in_db",
                        "date": now_str,
                    }
                )

        return len(df_new)

    def write_table_as_is(self, conn, df: pd.DataFrame, table: str) -> int:
        if df is None:
            return 0

        df = df.copy()

        if not self.update:
            df.to_sql(table, conn, if_exists="replace", index=False)
            print(f"[CreateSqliteDB] Wrote {len(df)} rows into '{table}' as-is")
        else:
            df.to_sql(table, conn, if_exists="append", index=False)
            print(f"[CreateSqliteDB] Appended {len(df)} rows into '{table}' as-is")

        return len(df)
    # ----------------------
    # DB creation
    # ----------------------
    def create_db(self):
        # Validate input files (always from args paths; in update mode these are under tmp/Update)
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
        self._require_file(self.host_taxa_file, "host_taxa_file")
        self._require_file(self.host_lineage_file, "host_lineage_file")
        self._require_file(self.host_children_file, "host_children_file")
        self._require_file(self.host_lineage_lookup_file, "host_lineage_lookup_file")

        excluded_records = []          # your existing excluded_accessions table (primarily meta_data exclusions)
        update_exclusions = []         # new: duplicates for any table in update mode

        # filtered sequence IDs handling (kept)
        filtered_ids = self._load_filtered_ids()
        filtered_details = self._load_filtered_details()
        if filtered_ids:
            print(f"[CreateSqliteDB] Excluding {len(filtered_ids)} filtered sequences from DB")
            for fid in filtered_ids:
                excluded_records.append(
                    {"primary_accession": fid, "reason": filtered_details.get(fid, "alignment_filtering")}
                )

        # Load incoming dataframes
        df_meta_data = self._read_tsv_required(self.meta_data, ["primary_accession"], "meta_data", dtype=str)
        df_meta_data_updates = df_meta_data.copy()

        df_meta_data = self._coerce_nullable_int_columns(df_meta_data,["length", "exclusion_status", "a", "t", "g", "c", "n", "real_length", "country_validated", "host_taxa_id"])

        # Track reference/master rows so they are always retained in meta_data
        acc_type_col = "accession_type" if "accession_type" in df_meta_data.columns else None
        if acc_type_col:
            acc_type_norm = df_meta_data[acc_type_col].fillna("").str.strip().str.lower()
            is_ref_or_master = acc_type_norm.isin(["reference", "master"])
        else:
            is_ref_or_master = pd.Series(False, index=df_meta_data.index)

        # Exclude filtered sequences from incoming meta_data, but never remove reference/master rows
        if filtered_ids and "primary_accession" in df_meta_data.columns:
            before_count = len(df_meta_data)
            remove_mask = df_meta_data["primary_accession"].isin(filtered_ids) & (~is_ref_or_master)
            df_meta_data = df_meta_data[~remove_mask]
            after_count = len(df_meta_data)
            if before_count != after_count:
                print(f"[CreateSqliteDB] Removed {before_count - after_count} filtered non-reference sequences from incoming meta_data")

        is_ref_or_master = is_ref_or_master.reindex(df_meta_data.index, fill_value=False)

        # Collect exclusions from meta_data (e.g. invalid division)
        if "exclusion" in df_meta_data.columns:
            exclusion_mask = df_meta_data["exclusion"].notna() & (df_meta_data["exclusion"] != "")
            excluded_rows = df_meta_data[exclusion_mask]
            if not excluded_rows.empty:
                print(f"[CreateSqliteDB] Found {len(excluded_rows)} rows with exclusions in incoming meta_data")
                for _, row in excluded_rows.iterrows():
                    acc = row.get("primary_accession", "")
                    if acc:
                        excluded_records.append({"primary_accession": acc, "reason": row["exclusion"]})

                # Remove excluded rows from incoming meta_data, but never remove reference/master rows
                remove_mask = exclusion_mask & (~is_ref_or_master)
                df_meta_data = df_meta_data[~remove_mask]
                retained_refs = (exclusion_mask & is_ref_or_master).sum()
                if retained_refs:
                    print(f"[CreateSqliteDB] Retained {retained_refs} reference/master rows despite exclusion flags")

        df_meta_data = self._add_cluster_column(df_meta_data)

        df_features = self._read_tsv_required(self.features, [], "features", dtype=str)
        df_features = self._coerce_nullable_int_columns(df_features,cols=["aln_start", "aln_end", "cds_start", "cds_end"])
        

        df_aln = self._read_tsv_required(self.pad_aln, [], "pad_aln")
        df_aln = self._normalize_alignment_columns(df_aln, "pad_aln")
        df_gene = self._read_tsv_required(self.gene_info, [], "gene_info")

        df_gene = self._coerce_nullable_int_columns(df_gene,cols=["start", "end"])

        df_m49_country = self._read_csv_required(
            self.m49_countries, ["m49_code"], "m49_countries", dtype={"m49_code": str}
        )
        df_m49_country = self._coerce_nullable_int_columns(df_m49_country,cols=["m49_code", "is_ldc", "is_lldc","is_sids"])
        
        df_m49_interm = self._read_csv_required(self.m49_interm_region, [], "m49_interm_region")
        df_m49_interm = self._coerce_nullable_int_columns(df_m49_interm,cols=["m49_code"])
        
        df_m49_region = self._read_csv_required(self.m49_regions, [], "m49_regions")
        df_m49_region = self._coerce_nullable_int_columns(
            df_m49_region,
            cols=["m49_code"])

        df_m49_sub_region = self._read_csv_required(self.m49_sub_regions, [], "m49_sub_regions")
        df_m49_sub_region = self._coerce_nullable_int_columns(df_m49_sub_region,cols=["m49_code"])
        
        df_proj_setting = self._read_tsv_required(self.proj_settings, [], "proj_settings")

        df_host_taxa = self._read_tsv_required(self.host_taxa_file, [], "host_taxa_file", dtype=str)
        df_host_taxa = self._coerce_nullable_int_columns(df_host_taxa,cols=["taxa_id"])
        
        df_host_lineage = self._read_tsv_required(self.host_lineage_file, [], "host_lineage_file", dtype=str)
        df_host_lineage = self._coerce_nullable_int_columns(df_host_lineage,cols=["host_taxa_id"])
        

        
        df_host_children = self._read_tsv_required(self.host_children_file, [], "host_children_file", dtype=str)
        df_host_children = self._coerce_nullable_int_columns(df_host_children,cols=["parent_taxa_id", "child_taxa_id"])
        
        df_host_lineage_lookup = self._read_tsv_required(
            self.host_lineage_lookup_file, [], "host_lineage_lookup_file", dtype=str
        )
        df_host_lineage_lookup = self._coerce_nullable_int_columns(df_host_lineage_lookup,cols=["lineage_taxa_id", "desc_taxa_id"])
        

        df_fasta_sequences = self.load_fasta()

        # Trees (optional)
        df_trees = None
        if self.tree_dir:
            df_trees = self.load_trees_from_dir(self.tree_dir)
        
        db_path = self._db_path()
        db_dir = os.path.dirname(db_path)

        if db_dir:
            os.makedirs(db_dir, exist_ok=True)

        self._backup_existing_db_if_update(db_path)

        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Ensure tables for info/trees exist when we need them
        cursor.execute("PRAGMA foreign_keys = ON;")
        self._ensure_trees_table(conn)
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS info (creation_type TEXT, date TEXT);"
        )
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS excluded_accessions (primary_accession TEXT, reason TEXT);"
        )
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS update_exclusions (table_name TEXT, key TEXT, reason TEXT, date TEXT);"
        )

        if self.update:
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS updates (
                    primary_accession TEXT,
                    updated_at TEXT
                );
                """
            )
        if self.update:
            self._log_updates(conn, df_meta_data_updates)
        # ----------------------
        # Write/merge tables
        # ----------------------
        # NOTE: normal mode = replace, update mode = append non-redundant.
        # meta_data


        self.merge_table_append_nonredundant(
            conn, df_meta_data, "meta_data", ["primary_accession"], None
        )

        # features
        features_key = [
            "accession",
            "reference_accession",
            "aln_start",
            "aln_end",
            "cds_start",
            "cds_end",
            "product",
        ]
        self.merge_table_append_nonredundant(conn, df_features, "features", features_key, update_exclusions)
      
        #self.merge_table_append_nonredundant(conn, df_features, "features", None, update_exclusions)

        # sequence_alignment
        # prefer composite (sequence_id, alignment_name) if available
        aln_key = ["sequence_id", "alignment_name"] if "alignment_name" in df_aln.columns else ["sequence_id"]
        self.merge_table_append_nonredundant(conn, df_aln, "sequence_alignment", aln_key, update_exclusions)

        # genes and M49 and project_settings
        #self.merge_table_append_nonredundant(conn, df_gene, "gene_info", None, update_exclusions)
        gene_key = ["accession", "category", "start", "end", "gene_name"]
        self.merge_table_append_nonredundant(conn, df_gene, "genes", gene_key, update_exclusions)

        self.merge_table_append_nonredundant(conn, df_m49_country, "m49_country", ["m49_code"], update_exclusions)
        self.merge_table_append_nonredundant(conn, df_m49_interm, "m49_intermediate", None, update_exclusions)
        self.merge_table_append_nonredundant(conn, df_m49_region, "m49_regions", None, update_exclusions)
        self.merge_table_append_nonredundant(conn, df_m49_sub_region, "m49_sub_regions", None, update_exclusions)
        self.merge_table_append_nonredundant(conn, df_proj_setting, "project_settings", None, update_exclusions)

        # sequences
        self.merge_table_append_nonredundant(conn, df_fasta_sequences, "sequences", ["header"], update_exclusions)

        # host tables
        self.write_table_as_is(conn, df_host_taxa, "host_taxa")
        self.write_table_as_is(conn, df_host_lineage, "host_lineage")
        self.write_table_as_is(conn, df_host_children, "host_children")
        self.write_table_as_is(conn, df_host_lineage_lookup, "host_lineage_lookup")


        # ----------------------
        # Exclusion tables (your existing behavior)
        # ----------------------
        if excluded_records:
            df_excluded = pd.DataFrame(excluded_records).drop_duplicates(subset=["primary_accession"])
            # append in update mode so you don't wipe previous exclusions
            if self.update and self._table_exists(conn, "excluded_accessions"):
                df_excluded.to_sql("excluded_accessions", conn, if_exists="append", index=False)
            else:
                df_excluded.to_sql("excluded_accessions", conn, if_exists="replace", index=False)
            print(f"[CreateSqliteDB] excluded_accessions now has >= {len(df_excluded)} newly added records")

        if update_exclusions:
            df_upd_excl = pd.DataFrame(update_exclusions)
            # always append
            df_upd_excl.to_sql("update_exclusions", conn, if_exists="append", index=False)
            print(f"[CreateSqliteDB] Logged {len(df_upd_excl)} update-mode duplicate keys into update_exclusions")

        # ----------------------
        # Trees (single Newick files + manifest) -> insert into trees table (kept behavior; fixed creation)
        # ----------------------
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # If an older run created the accidental singular "tree" table,
        # move compatible rows into "trees" and remove the duplicate table.
        self._migrate_legacy_tree_table(conn, now_str)

        # build accession->segment mapping from current DB meta_data (best effort)
        accession_to_segment = {}
        try:
            md_cols = pd.read_sql_query("PRAGMA table_info(meta_data);", conn)
            md_has_segment = "segment" in md_cols["name"].tolist()
        except Exception:
            md_has_segment = False

        if md_has_segment:
            try:
                df_md_seg = pd.read_sql_query("SELECT primary_accession, segment FROM meta_data", conn)
                for _, row in df_md_seg.dropna().iterrows():
                    acc = str(row["primary_accession"]).strip()
                    seg = str(row["segment"]).strip()
                    if acc and seg and acc not in accession_to_segment:
                        accession_to_segment[acc] = seg
            except Exception:
                pass

        tree_records = []

        # Add --tree_dir rows to the same canonical trees table.
        # Manifest columns are preserved as separate fields where possible:
        #   tree_name -> name
        #   chromosome -> segment_key
        #   segment_number -> segment
        if df_trees is not None and not df_trees.empty:
            for _, row in df_trees.iterrows():
                newick = str(row.get("newick", "")).strip()
                if not newick:
                    continue

                tree_name = str(row.get("tree_name", "")).strip()
                tree_type = str(row.get("tree_type", "")).strip()
                chromosome = str(row.get("chromosome", "")).strip()
                segment_number = str(row.get("segment_number", "")).strip()
                tree_model = str(row.get("tree_model", "")).strip()
                description = str(row.get("description", "")).strip()

                tree_records.append(
                    {
                        "name": tree_name or "tree_dir_tree",
                        "source": "tree_dir",
                        "tree_type": tree_type or None,
                        "tree_model": tree_model or None,
                        "segment_key": chromosome or None,
                        "segment": segment_number or self._segment_from_key(chromosome),
                        "newick": newick,
                        "description": description or None,
                    }
                )

        for name, source, tree_path in [
            ("veryfasttree", "veryfasttree", self.tree_file),
            ("iqtree", "iqtree", self.iqtree_file),
            ("usher", "usher", self.usher_tree),
        ]:
            newick = self._read_tree_file(tree_path)
            if newick:
                tree_records.append(
                    {
                        "name": name,
                        "source": source,
                        "tree_type": None,
                        "tree_model": None,
                        "segment_key": None,
                        "segment": None,
                        "newick": newick,
                    }
                )

        for entry in self._load_tree_manifest(self.tree_manifest):
            newick = self._read_tree_file(entry["path"])
            if not newick:
                continue
            seg_key = entry.get("segment_key")
            seg_num = accession_to_segment.get(seg_key) if seg_key else None
            if seg_num is None:
                seg_num = self._segment_from_key(seg_key)
            nm = entry.get("name") or f"{entry['source']}_{seg_key or 'tree'}"
            tree_records.append(
                {
                    "name": nm,
                    "source": entry["source"],
                    "tree_type": None,
                    "tree_model": None,
                    "segment_key": seg_key,
                    "segment": seg_num,
                    "newick": newick,
                }
            )

        # insert trees non-redundantly in update mode based on (source,name,segment_key,segment)
        if tree_records:
            df_tree_ins = pd.DataFrame(tree_records)
            df_tree_ins["created_at"] = now_str
            # create a temporary key column for dedupe and merge
            df_tree_ins["__k"] = (
                df_tree_ins["source"].fillna("").astype(str).str.strip()
                + "|"
                + df_tree_ins["name"].fillna("").astype(str).str.strip()
                + "|"
                + df_tree_ins["tree_type"].fillna("").astype(str).str.strip()
                + "|"
                + df_tree_ins["tree_model"].fillna("").astype(str).str.strip()
                + "|"
                + df_tree_ins["segment_key"].fillna("").astype(str).str.strip()
                + "|"
                + df_tree_ins["segment"].fillna("").astype(str).str.strip()
            )

            if not self.update:
                # replace by clearing trees table then insert
                cursor.execute("DELETE FROM trees;")
                for _, tr in df_tree_ins.iterrows():
                    cursor.execute(
                        """
                        INSERT INTO trees
                            (name, source, tree_type, tree_model, segment_key, segment, newick, created_at, description)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            tr["name"],
                            tr["source"],
                            tr["tree_type"],
                            tr["tree_model"],
                            tr["segment_key"],
                            tr["segment"],
                            tr["newick"],
                            tr["created_at"],
                            tr.get("description"),
                        ),
                    )
            else:
                # fetch existing keys
                cursor.execute("SELECT source, name, tree_type, tree_model, segment_key, segment FROM trees;")
                existing = set()
                for s, n, tt, tm, sk, sg in cursor.fetchall():
                    existing.add(
                        (
                            str(s or "").strip(),
                            str(n or "").strip(),
                            str(tt or "").strip(),
                            str(tm or "").strip(),
                            str(sk or "").strip(),
                            str(sg or "").strip(),
                        )
                    )
                appended = 0
                for _, tr in df_tree_ins.iterrows():
                    key = (
                        str(tr["source"] or "").strip(),
                        str(tr["name"] or "").strip(),
                        str(tr["tree_type"] or "").strip(),
                        str(tr["tree_model"] or "").strip(),
                        str(tr["segment_key"] or "").strip(),
                        str(tr["segment"] or "").strip(),
                    )
                    if key in existing:
                        continue
                    cursor.execute(
                        """
                        INSERT INTO trees
                            (name, source, tree_type, tree_model, segment_key, segment, newick, created_at, description)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            tr["name"],
                            tr["source"],
                            tr["tree_type"],
                            tr["tree_model"],
                            tr["segment_key"],
                            tr["segment"],
                            tr["newick"],
                            tr["created_at"],
                            tr.get("description"),
                        ),
                    )
                    appended += 1
                if appended:
                    print(f"[CreateSqliteDB] Appended {appended} new trees (non-redundant)")

        # ----------------------
        # info table
        # ----------------------
        creation_type = self._normalize_db_status(self.db_status)
        cursor.execute(
            "INSERT INTO info (creation_type, date) VALUES (?, ?)",
            (creation_type, now_str),
        )

        # ----------------------
        # Indexes (Dana’s list, safe)
        # ----------------------
        def _has_column(conn_, table: str, col: str) -> bool:
            try:
                rows = conn_.execute(f"PRAGMA table_info({table})").fetchall()
            except Exception:
                return False
            return any(r[1] == col for r in rows)

        def _create_index(table: str, col: str, index_name: str):
            if self._table_exists(conn, table) and _has_column(conn, table, col):
                conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON {table}({col})")

        # meta_data
        _create_index("meta_data", "isolate", "idx_meta_data_isolate")
        _create_index("meta_data", "primary_accession", "idx_meta_data_primary_accession")
        _create_index("meta_data", "host", "idx_meta_data_host")
        _create_index("meta_data", "length", "idx_meta_data_length")
        _create_index("meta_data", "accession_type", "idx_meta_data_accession_type")
        _create_index("meta_data", "strain", "idx_meta_data_strain")
        # _create_index("meta_data", "segment", "idx_meta_data_segment")

        # sequences
        _create_index("sequences", "header", "idx_sequences_header")

        # sequence_alignment
        _create_index("sequence_alignment", "sequence_id", "idx_sequence_alignment_sequence_id")
        _create_index("sequence_alignment", "alignment_name", "idx_sequence_alignment_alignment_name")

        # features
        _create_index("features", "accession", "idx_features_accession")
        _create_index("features", "reference_accession", "idx_features_reference_accession")

        # cluster_members (ONLY if that table exists in your DB build)
        _create_index("cluster_members", "primary_accession", "idx_cluster_members_primary_accession")

        conn.commit()
        conn.close()

    # adding new sequences to the update table
    def _log_updates(self, conn, df_meta_data: pd.DataFrame) -> int:

        if not self.update or df_meta_data is None or df_meta_data.empty:
            return 0

        if "primary_accession" not in df_meta_data.columns:
            raise ValueError("meta_data is missing required column: primary_accession")

        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS updates (
                primary_accession TEXT,
                updated_at TEXT
            );
            """
        )

        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        df_updates = (
            df_meta_data[["primary_accession"]]
            .copy()
            .dropna(subset=["primary_accession"])
        )
        df_updates["primary_accession"] = df_updates["primary_accession"].astype(str).str.strip()
        df_updates = df_updates[df_updates["primary_accession"] != ""]
        df_updates = df_updates.drop_duplicates(subset=["primary_accession"], keep="first")
        df_updates["updated_at"] = now_str

        if df_updates.empty:
            return 0

        df_updates.to_sql("updates", conn, if_exists="append", index=False)
        print(f"[CreateSqliteDB] Logged {len(df_updates)} accessions into updates table")

        return len(df_updates)
    
def process(args):
    db_creator = CreateSqliteDB(
        meta_data=args.meta_data,
        features=args.features,
        pad_aln=args.pad_aln,
        gene_info=args.gene_info,
        m49_countries=args.m49_countries,
        m49_interm_region=args.m49_interm_region,
        m49_regions=args.m49_regions,
        m49_sub_regions=args.m49_sub_regions,
        proj_settings=args.proj_settings,
        fasta_sequence_file=args.fasta_sequences,
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        db_name=args.db_name,
        db_status=args.db_status,
        host_taxa_file=args.host_taxa_file,
        host_lineage_file=args.host_lineage_file,
        host_children_file=args.host_children_file,
        host_lineage_lookup_file=args.host_lineage_lookup_file,
        db_file=args.db_file,
        tree_file=args.tree_file,
        iqtree_file=args.iqtree_file,
        usher_tree=args.usher_tree,
        tree_dir=args.tree_dir,
        cluster_tsv=args.cluster_tsv,
        cluster_min_seq_id=args.cluster_min_seq_id,
        filtered_ids_file=args.filtered_ids,
        filtered_details_file=args.filtered_details,
        tree_manifest=args.tree_manifest,
        update=args.update,
    )
    db_creator.create_db()

if __name__ == "__main__":
    parser = ArgumentParser(description="Creating sqlite DB")
    parser.add_argument("-m", "--meta_data", help="Meta data table", default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
    parser.add_argument("-b", "--base_dir", help="Base directory", default="tmp")
    parser.add_argument("-o", "--output_dir", help="tmp directory where the database is stored", default="SqliteDB")
    parser.add_argument("-rf", "--features", help="Features table", default="tmp/Tables/features.tsv")
    parser.add_argument("-p", "--pad_aln", help="Padded alignment file", default="tmp/Tables/sequence_alignment.tsv")
    parser.add_argument("-g", "--gene_info", help="Gene table", default="tmp/Tables/genes.tsv")
    parser.add_argument("-mc", "--m49_countries", help="M49 countries", default="assets/m49_country.csv")
    parser.add_argument(
        "-mir", "--m49_interm_region", help="M49 intermediate regions", default="assets/m49_intermediate_region.csv"
    )
    parser.add_argument("-mr", "--m49_regions", help="M49 regions", default="assets/m49_region.csv")
    parser.add_argument("-msr", "--m49_sub_regions", help="M49 sub-regions", default="assets/m49_sub_region.csv")
    parser.add_argument("-s", "--proj_settings", help="Project settings", default="tmp/Software_info/software_info.tsv")
    parser.add_argument("-fa", "--fasta_sequences", help="Fasta sequences", default="tmp/GenBank-matrix/sequences.fa")
    parser.add_argument(
        "--db_file",
        help="Full path to an existing or new SQLite DB file. If provided, this overrides --base_dir/--output_dir/--db_name.",
        default=None,
    )
    parser.add_argument("-d", "--db_name", help="Name of the Sqlite database", default="gdb")
    parser.add_argument(
        "-ds",
        "--db_status",
        help='Database status: "new db" (default) or "last modified"/"last updated". Determines info.creation_type.',
        default="new db",
    )
    parser.add_argument("-t", "--tree_file", help="VeryFastTree Newick file", default=None)
    parser.add_argument("-it", "--iqtree_file", help="IQ-TREE Newick file", default=None)
    parser.add_argument("-ut", "--usher_tree", help="UShER output Newick file", default=None)
    parser.add_argument("--tree_manifest", help="TSV manifest with columns: source, name, segment_key, path", default=None)
    parser.add_argument("-ct", "--cluster_tsv", help="MMseqs clustering TSV (rep\\tmember)", default=None)
    parser.add_argument("-ci", "--cluster_min_seq_id", help="MMseqs min sequence identity used for clustering", default=None)
    parser.add_argument("-fi", "--filtered_ids", help="File with filtered sequence IDs (one per line) to exclude from DB", default=None)
    parser.add_argument(
        "-fd",
        "--filtered_details",
        help="TSV with filtered sequence details (seq_name, reference, error, warnings)",
        default=None,
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="If enabled, merge tables into the existing DB (append-only, non-redundant).",
    )
    parser.add_argument(
        "--tree_dir",
        help="Directory containing tree files and a manifest meta_data.tsv (chromosome, segment_number, tree_type, tree_name, tree_model, description).",
        default=None,
    )
    parser.add_argument("-ht", "--host_taxa_file", help="Host Taxanomy file", default="tmp/HostTaxa/Host_taxa.tsv")
    parser.add_argument("-hl", "--host_lineage_file", help="Host Lineage file", default="tmp/HostTaxa/Host_taxa_lineage.tsv")
    parser.add_argument("-hc", "--host_children_file", help="Host Children file", default="tmp/HostTaxa/Host_taxa_children.tsv")
    parser.add_argument(
        "-hll",
        "--host_lineage_lookup_file",
        help="Host Lineage lookup file",
        default="tmp/HostTaxa/Host_taxa_lineage_lookup.tsv",
    )

    args = parser.parse_args()

    # Keep your existing behavior: in update mode, base_dir -> base_dir/Update (avoid Update/Update)
    if args.update:
        if not normpath(args.base_dir).endswith(normpath("Update")):
            args.base_dir = join(args.base_dir, "Update")

    try:
        process(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
