#!/usr/bin/env python3
import os
import sqlite3
import argparse
from os.path import join
from typing import Dict, List, Optional, Set, Tuple, Any


class ValidateDB:
    """
    Validate consistency across V-gTK SQLite tables.

    Rules implemented:
      - meta_data (after exclusions) defines expected accessions.
      - sequence_alignment.primary_accession must match expected accessions.
      - features.accession must match expected accessions.
      - host_taxa.taxa_id must contain all host_taxa_id values used in meta_data (after exclusions).

    Notes:
      - Exclusion is controlled by --exclude_value (default "1").
      - NA-like values are ignored for host_taxa_id checks.
    """

    NA_SET = {"", "NA", "Na", "N/A", "na", "n/a", "-", None}

    def __init__(
        self,
        accession_column: str,
        base_dir: str,
        output_dir: str,
        output_file: str,
        db: str,
        exclusion_column: str = "exclusion_status",
        exclude_value: str = "1",
    ) -> None:
        self.accession_column = accession_column
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.output_file = output_file
        self.db = db
        self.exclusion_column = exclusion_column
        self.exclude_value = str(exclude_value)

        os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)

    # ---------------- SQLite helpers ----------------

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db)
        conn.row_factory = sqlite3.Row
        return conn

    def table_exists(self, table_name: str) -> bool:
        conn = self._connect()
        try:
            row = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
                (table_name,),
            ).fetchone()
            return row is not None
        finally:
            conn.close()

    def get_table_columns(self, table_name: str) -> List[str]:
        conn = self._connect()
        try:
            rows = conn.execute(f"PRAGMA table_info({table_name})").fetchall()
            return [r["name"] for r in rows]
        finally:
            conn.close()

    def fetch_distinct_values(
        self,
        table: str,
        column: str,
        where_sql: Optional[str] = None,
        params: Tuple[Any, ...] = (),
    ) -> Set[str]:
        conn = self._connect()
        try:
            sql = f"SELECT DISTINCT {column} AS v FROM {table}"
            if where_sql:
                sql += f" WHERE {where_sql}"
            rows = conn.execute(sql, params).fetchall()
            out: Set[str] = set()
            for r in rows:
                v = r["v"]
                if v is None:
                    continue
                s = str(v).strip()
                if s:
                    out.add(s)
            return out
        finally:
            conn.close()

    def fetch_count(
        self,
        table: str,
        where_sql: Optional[str] = None,
        params: Tuple[Any, ...] = (),
    ) -> int:
        conn = self._connect()
        try:
            sql = f"SELECT COUNT(*) AS n FROM {table}"
            if where_sql:
                sql += f" WHERE {where_sql}"
            row = conn.execute(sql, params).fetchone()
            return int(row["n"]) if row else 0
        finally:
            conn.close()

    # ---------------- Core logic ----------------

    def _meta_where_nonexcluded(self, meta_cols: List[str]) -> Tuple[Optional[str], Tuple[Any, ...]]:
        """
        Build a WHERE clause to select non-excluded meta_data rows.

        Default: include rows where exclusion_column != exclude_value OR exclusion_column IS NULL
        If exclusion_column doesn't exist, no filter.
        """
        if self.exclusion_column not in meta_cols:
            return None, ()

        # treat NULL as not excluded
        where_sql = f"({self.exclusion_column} IS NULL OR CAST({self.exclusion_column} AS TEXT) != ?)"
        return where_sql, (self.exclude_value,)

    def get_expected_accessions(self) -> Dict[str, Any]:
        """
        Return expected accession set from meta_data after exclusions.
        """
        if not self.table_exists("meta_data"):
            raise RuntimeError("Table 'meta_data' not found in DB.")

        meta_cols = self.get_table_columns("meta_data")
        if self.accession_column not in meta_cols:
            raise RuntimeError(
                f"meta_data does not contain accession column '{self.accession_column}'. "
                f"Available columns: {meta_cols}"
            )

        where_sql, params = self._meta_where_nonexcluded(meta_cols)

        expected = self.fetch_distinct_values(
            table="meta_data",
            column=self.accession_column,
            where_sql=where_sql,
            params=params,
        )

        total_rows = self.fetch_count("meta_data")
        nonexcluded_rows = self.fetch_count("meta_data", where_sql=where_sql, params=params) if where_sql else total_rows

        return {
            "expected_accessions": expected,
            "meta_total_rows": total_rows,
            "meta_nonexcluded_rows": nonexcluded_rows,
            "meta_has_exclusion": self.exclusion_column in meta_cols,
        }

    def compare_sets(self, expected: Set[str], observed: Set[str]) -> Dict[str, Any]:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        return {
            "ok": (len(missing) == 0 and len(extra) == 0),
            "expected_count": len(expected),
            "observed_count": len(observed),
            "missing": missing,
            "extra": extra,
        }

    def validate_sequence_alignment(self, expected_accessions: Set[str]) -> Dict[str, Any]:
        table = "sequence_alignment"
        if not self.table_exists(table):
            return {"ok": False, "error": f"Table '{table}' not found."}

        cols = self.get_table_columns(table)
        # try common column names if your schema differs
        candidate_cols = ["primary_accession", "accession", "sequence_id"]
        col = next((c for c in candidate_cols if c in cols), None)
        if col is None:
            return {
                "ok": False,
                "error": f"Table '{table}' missing any of columns {candidate_cols}. Columns present: {cols}",
            }

        observed = self.fetch_distinct_values(table=table, column=col)
        res = self.compare_sets(expected_accessions, observed)
        res.update({"table": table, "column": col})
        return res

    def validate_features(self, expected_accessions: Set[str]) -> Dict[str, Any]:
        table = "features"
        if not self.table_exists(table):
            return {"ok": False, "error": f"Table '{table}' not found."}

        cols = self.get_table_columns(table)
        col = "accession" if "accession" in cols else None
        if col is None:
            return {"ok": False, "error": f"Table '{table}' missing column 'accession'. Columns present: {cols}"}

        observed = self.fetch_distinct_values(table=table, column=col)
        res = self.compare_sets(expected_accessions, observed)
        res.update({"table": table, "column": col})
        return res

    def validate_host_taxa(self) -> Dict[str, Any]:
        """
        Compare host_taxa_id used in meta_data (after exclusions) vs host_taxa.taxa_id.
        """
        if not self.table_exists("meta_data"):
            return {"ok": False, "error": "Table 'meta_data' not found."}
        if not self.table_exists("host_taxa"):
            return {"ok": False, "error": "Table 'host_taxa' not found."}

        meta_cols = self.get_table_columns("meta_data")
        if "host_taxa_id" not in meta_cols:
            return {"ok": False, "error": "meta_data missing column 'host_taxa_id'."}

        where_sql, params = self._meta_where_nonexcluded(meta_cols)
        meta_hosts = self.fetch_distinct_values(
            table="meta_data",
            column="host_taxa_id",
            where_sql=where_sql,
            params=params,
        )

        # drop NA-like
        meta_hosts = {h for h in meta_hosts if (h not in self.NA_SET and str(h).strip() not in self.NA_SET)}

        host_cols = self.get_table_columns("host_taxa")
        col = "taxa_id" if "taxa_id" in host_cols else None
        if col is None:
            return {"ok": False, "error": f"host_taxa missing column 'taxa_id'. Columns present: {host_cols}"}

        observed = self.fetch_distinct_values(table="host_taxa", column=col)

        res = self.compare_sets(meta_hosts, observed)
        res.update({"table": "host_taxa", "column": col, "expected_count_source": "meta_data.host_taxa_id"})
        return res

    # ---------------- Reporting ----------------

    def _format_result(self, title: str, res: Dict[str, Any], show_n: int = 25) -> str:
        lines = [f"[{title}]"]
        if "error" in res:
            lines.append(f"  OK: False")
            lines.append(f"  ERROR: {res['error']}")
            return "\n".join(lines)

        lines.append(f"  OK: {res.get('ok')}")
        lines.append(f"  Table: {res.get('table')}  Column: {res.get('column')}")
        lines.append(f"  Expected: {res.get('expected_count')}  Observed: {res.get('observed_count')}")
        miss = res.get("missing", [])
        extra = res.get("extra", [])
        lines.append(f"  Missing: {len(miss)}  Extra: {len(extra)}")

        if miss:
            lines.append(f"  Missing examples (up to {show_n}):")
            for x in miss[:show_n]:
                lines.append(f"    - {x}")
        if extra:
            lines.append(f"  Extra examples (up to {show_n}):")
            for x in extra[:show_n]:
                lines.append(f"    - {x}")
        return "\n".join(lines)

    def validate(self) -> None:
        meta_info = self.get_expected_accessions()
        expected_accessions = meta_info["expected_accessions"]

        # Run validations
        seq_res = self.validate_sequence_alignment(expected_accessions)
        feat_res = self.validate_features(expected_accessions)
        host_res = self.validate_host_taxa()

        # Build summary
        lines: List[str] = []
        lines.append("=== ValidateDB Summary ===")
        lines.append(f"DB: {self.db}")
        lines.append(f"meta_data total rows: {meta_info['meta_total_rows']}")
        if meta_info["meta_has_exclusion"]:
            lines.append(
                f"meta_data non-excluded rows: {meta_info['meta_nonexcluded_rows']} "
                f"(excluding {self.exclusion_column} == {self.exclude_value})"
            )
        else:
            lines.append("meta_data exclusion column not found; using all rows.")
        lines.append(f"Expected distinct accessions (from meta_data.{self.accession_column}): {len(expected_accessions)}")
        lines.append("")

        lines.append(self._format_result("sequence_alignment vs meta_data", seq_res))
        lines.append("")
        lines.append(self._format_result("features vs meta_data", feat_res))
        lines.append("")
        lines.append(self._format_result("host_taxa vs meta_data.host_taxa_id", host_res))
        lines.append("")
        lines.append("===========================")

        report = "\n".join(lines)

        # Print to console
        print(report)

        # Write to file
        out_path = join(self.base_dir, self.output_dir, self.output_file)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(report + "\n")

        print(f"\nReport written to: {out_path}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate a V-gTK style SQLite DB for basic consistency checks.")
    parser.add_argument("-id", "--accession_column", default="primary_accession",
                        help="Accession column name in meta_data used as the expected set")
    parser.add_argument("-b", "--base_dir", default="tmp", help="Base directory")
    parser.add_argument("-o", "--output_dir", default="ValidateDB", help="Output directory")
    parser.add_argument("-f", "--output_file", default="Validation_summary.txt", help="Output report filename")
    parser.add_argument("--db", required=True, help="Path to SQLite DB")
    parser.add_argument("--exclusion_column", default="exclusion_status",
                        help="Exclusion column in meta_data (default: exclusion_status)")
    parser.add_argument("--exclude_value", default="1",
                        help="Value in exclusion_column that means excluded (default: 1)")
    args = parser.parse_args()

    ValidateDB(
        accession_column=args.accession_column,
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        output_file=args.output_file,
        db=args.db,
        exclusion_column=args.exclusion_column,
        exclude_value=str(args.exclude_value),
    ).validate()
