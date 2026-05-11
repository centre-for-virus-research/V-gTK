#!/usr/bin/env python3
import os
import csv
import re
import shutil
import sqlite3
from datetime import datetime
from argparse import ArgumentParser
from typing import List, Dict, Tuple


class BlastAlignment:
    COMMENT_COL = "comment"   # the free-text comment column name (both files)
    _BLANKISH = {"", "-", "na", "n/a", "NA", "Na", "N/A"}

    def __init__(
        self,
        gb_matrix: str,
        curated_file: str,
        base_dir: str = "tmp",
        output_dir: str = "Curated",
        output_file: str = "gB_matrix.tsv",
        update: bool = False,
        db_file: str = None,
        primary_accession: str = "primary_accession",
    ) -> None:
        self.gb_matrix = gb_matrix
        self.curated_file = curated_file
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.output_file = output_file
        self.update = bool(update)
        self.db_file = db_file
        self.primary_accession = primary_accession

    # ---------------------- Public API ----------------------

    def process(self) -> None:
        # Load source from TSV or DB
        if self.gb_matrix:
            gb_header, gb_rows = self._read_tsv(self.gb_matrix, "GenBank matrix")
            source_label = self.gb_matrix
        else:
            if not self.update:
                raise ValueError("--gb_matrix is required unless --update is enabled with --db_file.")
            gb_header, gb_rows = self._read_meta_data_from_db("meta_data")
            source_label = f"{self.db_file}:meta_data"

        cur_header, cur_rows = self._read_tsv(self.curated_file, "curated file")

        # Required columns
        self._require_column(gb_header, self.primary_accession, "GenBank matrix / DB meta_data")
        self._require_column(cur_header, self.primary_accession, "curated file")

        # Optional columns
        comment_col = self.COMMENT_COL if self.COMMENT_COL in gb_header else None
        curator_col = "curator" if "curator" in gb_header else None

        # Update mode requirements
        if self.update and not self.db_file:
            raise ValueError("--db_file is required when --update is enabled.")

        # Index source matrix by primary accession
        gb_index: Dict[str, Dict[str, str]] = {}
        for row in gb_rows:
            k = (row.get(self.primary_accession, "") or "").strip()
            if k and k not in gb_index:
                gb_index[k] = row

        # Columns to update (excluding comment and curator)
        candidate_cols = [
            c for c in cur_header
            if c not in (self.COMMENT_COL, "curator")
        ]

        updated_rows = 0
        updated_cells = 0
        missing = 0
        ignored_cols = set()

        # Track rows that actually changed, for DB syncing
        rows_for_db_sync: List[Dict[str, str]] = []

        for crow in cur_rows:
            key = (crow.get(self.primary_accession, "") or "").strip()
            if not key:
                continue

            target = gb_index.get(key)
            if target is None:
                missing += 1
                continue

            row_changes = 0
            change_notes: List[str] = []

            # --- 1) Update metadata fields ---
            for col in candidate_cols:
                if col == self.primary_accession:
                    continue

                if col not in gb_header:
                    ignored_cols.add(col)
                    continue

                new_val_raw = crow.get(col, "")
                new_val = (new_val_raw if new_val_raw is not None else "").strip()
                if self._is_blankish(new_val):
                    continue

                old_val_raw = target.get(col, "")
                old_val = (old_val_raw if old_val_raw is not None else "").strip()

                if new_val != old_val:
                    target[col] = new_val
                    updated_cells += 1
                    row_changes += 1
                    change_notes.append(f"{col}: '{old_val}' -> '{new_val}'")

            # --- 2) Only append comment/curator if real metadata changed ---
            if row_changes > 0:
                updated_rows += 1

                # Handle comment
                if comment_col:
                    existing_comment = (target.get(comment_col, "") or "").strip()
                    curated_comment = (crow.get(comment_col, "") or "").strip()
                    existing_comment = self._strip_unwanted_fragments(existing_comment)

                    if not self._is_blankish(curated_comment):
                        existing_comment = self._append_unique(existing_comment, curated_comment)

                    if change_notes:
                        change_blob = "; ".join(change_notes)
                        existing_comment = self._append_unique(existing_comment, change_blob)

                    target[comment_col] = self._clean_separators(existing_comment)

                # Handle curator
                if curator_col:
                    existing_curator = (target.get(curator_col, "") or "").strip()
                    curated_curator = (crow.get(curator_col, "") or "").strip()
                    if not self._is_blankish(curated_curator):
                        target[curator_col] = self._append_unique(existing_curator, curated_curator)

                # Store updated row for DB sync
                rows_for_db_sync.append(dict(target))

        # --- Write output TSV only if source TSV path was provided ---
        out_path = None
        if self.gb_matrix:
            out_dir = os.path.join(self.base_dir, self.output_dir)
            os.makedirs(out_dir, exist_ok=True)

            out_path = os.path.join(out_dir, self.output_file)
            self._write_tsv(out_path, gb_header, gb_rows)

            # overwrite original input matrix
            self._write_tsv(self.gb_matrix, gb_header, gb_rows)

        # --- Sync curated rows into DB only in update mode ---
        db_inserted = db_updated = db_skipped = 0
        if self.update and rows_for_db_sync:
            self._backup_existing_db_if_update(self.db_file)
            db_inserted, db_updated, db_skipped = self._sync_curated_rows_to_db(
                gb_header=gb_header,
                rows_to_sync=rows_for_db_sync,
                table_name="meta_data",
            )

        # --- Summary ---
        print("\n=== Curation Summary ===")
        print(f"Source               : {source_label}")
        print(f"Original rows        : {len(gb_rows)}")
        print(f"Curated rows         : {len(cur_rows)}")
        print(f"Rows updated         : {updated_rows}")
        print(f"Cells updated        : {updated_cells}")
        print(f"Curated keys missing : {missing}")
        if ignored_cols:
            print(f"Ignored curated cols : {sorted(ignored_cols)} (not in source header)")
        if out_path:
            print(f"Output               : {out_path}")
        else:
            print("Output               : skipped TSV write (DB-only update mode)")
        if self.update:
            print(f"DB file              : {self.db_file}")
            print(f"DB inserted          : {db_inserted}")
            print(f"DB updated           : {db_updated}")
            print(f"DB skipped           : {db_skipped}")
        print("========================\n")

    # ---------------------- Utilities ----------------------

    @staticmethod
    def _timestamp_for_backup() -> str:
        return datetime.now().strftime("%Y%m%d_%H%M%S")

    def _backup_existing_db_if_update(self, db_path: str) -> str | None:
        if not self.update:
            return None

        if not db_path:
            raise ValueError("--db_file is required when --update is enabled.")

        if not os.path.isfile(db_path):
            raise FileNotFoundError(f"Update mode enabled but DB file does not exist: {db_path}")

        db_dir = os.path.dirname(db_path)
        db_file = os.path.basename(db_path)
        db_stem, db_ext = os.path.splitext(db_file)
        timestamp = self._timestamp_for_backup()

        backup_name = f"{db_stem}_backup_{timestamp}{db_ext or '.db'}"
        backup_path = os.path.join(db_dir, backup_name)

        counter = 1
        while os.path.exists(backup_path):
            backup_name = f"{db_stem}_backup_{timestamp}_{counter}{db_ext or '.db'}"
            backup_path = os.path.join(db_dir, backup_name)
            counter += 1

        shutil.copy2(db_path, backup_path)
        print(f"[Curator] Backup created: {backup_path}")
        return backup_path

    @staticmethod
    def _table_exists(conn, table: str) -> bool:
        row = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (table,),
        ).fetchone()
        return row is not None

    def _read_meta_data_from_db(self, table_name: str = "meta_data") -> Tuple[List[str], List[Dict[str, str]]]:
        if not self.db_file:
            raise ValueError("--db_file is required to read meta_data from DB.")

        if not os.path.isfile(self.db_file):
            raise FileNotFoundError(f"DB file not found: {self.db_file}")

        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            if not self._table_exists(conn, table_name):
                raise ValueError(f"Table '{table_name}' does not exist in DB: {self.db_file}")

            cur = conn.execute(f"SELECT * FROM {table_name}")
            rows_raw = cur.fetchall()
            header = [desc[0] for desc in cur.description] if cur.description else []

            rows: List[Dict[str, str]] = []
            for rec in rows_raw:
                row_dict = {}
                for col in header:
                    val = rec[col]
                    row_dict[col] = "" if val is None else str(val)
                rows.append(row_dict)

            return header, rows
        finally:
            conn.close()

    def _sync_curated_rows_to_db(
        self,
        gb_header: List[str],
        rows_to_sync: List[Dict[str, str]],
        table_name: str = "meta_data",
    ) -> Tuple[int, int, int]:
        """
        Upsert curated rows into SQLite meta_data table using self.primary_accession as key.
        - If key exists: UPDATE matching row
        - If key does not exist: INSERT new row

        Returns:
            inserted_count, updated_count, skipped_count
        """
        if not self.update:
            return 0, 0, 0

        if not self.db_file:
            raise ValueError("--db_file is required when --update is enabled.")

        if self.primary_accession not in gb_header:
            raise ValueError(
                f"Primary accession column '{self.primary_accession}' is missing in source header."
            )

        conn = sqlite3.connect(self.db_file)
        try:
            if not self._table_exists(conn, table_name):
                raise ValueError(f"Table '{table_name}' does not exist in DB: {self.db_file}")

            db_cols_info = conn.execute(f"PRAGMA table_info({table_name})").fetchall()
            db_columns = [row[1] for row in db_cols_info]

            if self.primary_accession not in db_columns:
                raise ValueError(
                    f"Primary accession column '{self.primary_accession}' is missing in DB table '{table_name}'."
                )

            usable_cols = [c for c in gb_header if c in db_columns]

            inserted_count = 0
            updated_count = 0
            skipped_count = 0

            for row in rows_to_sync:
                key_val = (row.get(self.primary_accession, "") or "").strip()
                if not key_val:
                    skipped_count += 1
                    continue

                exists = conn.execute(
                    f"SELECT 1 FROM {table_name} WHERE {self.primary_accession} = ? LIMIT 1",
                    (key_val,),
                ).fetchone()

                if exists:
                    update_cols = [c for c in usable_cols if c != self.primary_accession]
                    if update_cols:
                        set_clause = ", ".join([f"{c} = ?" for c in update_cols])
                        values = [row.get(c, "") for c in update_cols]
                        values.append(key_val)

                        conn.execute(
                            f"UPDATE {table_name} SET {set_clause} WHERE {self.primary_accession} = ?",
                            values,
                        )
                        updated_count += 1
                    else:
                        skipped_count += 1
                else:
                    insert_cols = usable_cols
                    placeholders = ", ".join(["?"] * len(insert_cols))
                    col_sql = ", ".join(insert_cols)
                    values = [row.get(c, "") for c in insert_cols]

                    conn.execute(
                        f"INSERT INTO {table_name} ({col_sql}) VALUES ({placeholders})",
                        values,
                    )
                    inserted_count += 1

            conn.commit()
            print(
                f"[Curator] DB sync completed for table '{table_name}': "
                f"{inserted_count} inserted, {updated_count} updated, {skipped_count} skipped"
            )
            return inserted_count, updated_count, skipped_count

        finally:
            conn.close()

    @staticmethod
    def _is_blankish(val: str) -> bool:
        return (
            (val is None)
            or (val.strip() == "")
            or (val.strip() in BlastAlignment._BLANKISH)
            or (val.strip().lower() in {"na", "n/a"})
        )

    @staticmethod
    def _split_comment_items(s: str) -> List[str]:
        parts = [p.strip() for p in s.split(";") if p.strip()]
        return parts

    def _append_unique(self, existing: str, new_item: str) -> str:
        if self._is_blankish(new_item):
            return existing
        existing_items = self._split_comment_items(existing)
        existing_norm = {re.sub(r"\s+", " ", it.lower()) for it in existing_items}
        new_norm = re.sub(r"\s+", " ", new_item.strip().lower())
        if new_norm not in existing_norm:
            existing_items.append(new_item.strip())
        return "; ".join(existing_items)

    @staticmethod
    def _clean_separators(s: str) -> str:
        s = re.sub(r"(^\s*NA\s*;?\s*|\s*;?\s*NA\s*$)", "", s, flags=re.IGNORECASE)
        s = re.sub(r"\s*;\s*;\s*", "; ", s)
        return s.strip(" ;")

    @staticmethod
    def _strip_unwanted_fragments(s: str) -> str:
        if not s:
            return s
        items = [
            p for p in BlastAlignment._split_comment_items(s)
            if not p.lower().strip().startswith("curator:")
        ]
        return "; ".join(items)

    # ---------- I/O helpers ----------

    def _read_tsv(self, path: str, label: str) -> Tuple[List[str], List[Dict[str, str]]]:
        if not path:
            raise ValueError(f"{label} path is missing.")
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found: {path}")
        with open(path, "r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            header = reader.fieldnames or []
            rows: List[Dict[str, str]] = []
            for rec in reader:
                for k in list(rec.keys()):
                    rec[k] = "" if rec[k] is None else str(rec[k])
                rows.append(rec)
        return header, rows

    def _write_tsv(self, path: str, header: List[str], rows: List[Dict[str, str]]) -> None:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=header,
                delimiter="\t",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator="\n",
            )
            writer.writeheader()
            for rec in rows:
                writer.writerow({c: rec.get(c, "") for c in header})

    @staticmethod
    def _require_column(header: List[str], col: str, label: str) -> None:
        if col not in header:
            raise ValueError(f"Required column '{col}' missing in {label}.")


# ---------------------- CLI Wrapper ----------------------

def _build_argparser() -> ArgumentParser:
    p = ArgumentParser(description=(
        "Replace fields in gB_matrix or DB meta_data from curated TSV by primary accession. "
        "Curator/comment fields are appended only if other metadata fields differ."
    ))
    p.add_argument(
        '-g', '--gb_matrix',
        help='GenBank matrix file (TSV). Optional in update mode.',
        default=None
    )
    p.add_argument(
        '-c', '--curated_file',
        help='Curated file (TSV)',
        default='generic/curation.tsv'
    )
    p.add_argument('-b', '--base_dir', help='Base directory for outputs', default='tmp')
    p.add_argument('-t', '--output_dir', help='Output subdirectory', default='Curated')
    p.add_argument('-o', '--output_file', help='Output filename (TSV)', default='gB_matrix_raw.tsv')

    p.add_argument(
        '--update',
        action='store_true',
        help='If enabled, also sync curated rows into the SQLite meta_data table.'
    )
    p.add_argument(
        '--db_file',
        default=None,
        help='Full path to the SQLite DB file. Required when --update is enabled.'
    )
    p.add_argument(
        '--primary_accession',
        default='primary_accession',
        help='Column name used as key for matching and DB update.'
    )

    return p


def main() -> None:
    args = _build_argparser().parse_args()

    if args.update and not args.db_file:
        raise ValueError("--db_file is required when --update is enabled.")

    if (not args.update) and (not args.gb_matrix):
        raise ValueError("--gb_matrix is required in normal mode.")

    BlastAlignment(
        gb_matrix=args.gb_matrix,
        curated_file=args.curated_file,
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        output_file=args.output_file,
        update=args.update,
        db_file=args.db_file,
        primary_accession=args.primary_accession,
    ).process()


if __name__ == "__main__":
    main()