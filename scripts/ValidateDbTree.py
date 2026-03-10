import argparse
import os
import sqlite3
from io import StringIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

ACCESSION_COLUMNS = {
    "primary_accession",
    "accession",
    "accession_id",
    "accession_version",
    "sequence_id",
    "sequence_name",
    "seq_name",
    "sample",
    "query",
    "reference",
    "header",
}

CLUSTER_PLACEHOLDER_VALUES = {
    "na-see-tree",
    "na- see tree",
    "na - see tree",
}


def get_table_columns(conn, table_name):
    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info({table_name})")
    return [row[1] for row in cursor.fetchall()]


def collect_table_values(conn, table_name, columns):
    values = set()
    cursor = conn.cursor()
    for col in columns:
        cursor.execute(f"SELECT DISTINCT {col} FROM {table_name}")
        values.update([r[0] for r in cursor.fetchall() if r[0]])
    return values


def find_cluster_column(columns):
    for col in columns:
        if col == "cluster" or col.startswith("cluster_"):
            return col
    return None


def is_cluster_placeholder(value):
    if value is None:
        return True
    norm = str(value).strip().lower()
    if not norm:
        return True
    return norm in CLUSTER_PLACEHOLDER_VALUES


def fetch_tree_newick(conn):
    cursor = conn.cursor()
    cursor.execute(
        "SELECT name, newick FROM trees WHERE newick IS NOT NULL AND length(newick) > 0 ORDER BY CASE WHEN source='usher' THEN 0 WHEN source='iqtree' THEN 1 ELSE 2 END, name LIMIT 1"
    )
    row = cursor.fetchone()
    if not row:
        return None, None
    return row[0], row[1]


def fetch_trees(conn, source=None):
    cols = get_table_columns(conn, "trees")
    select_cols = ["name", "source", "newick"]
    if "segment_key" in cols:
        select_cols.append("segment_key")
    if "segment" in cols:
        select_cols.append("segment")

    cursor = conn.cursor()
    query = f"SELECT {', '.join(select_cols)} FROM trees WHERE newick IS NOT NULL AND length(newick) > 0"
    params = []
    if source:
        query += " AND source = ?"
        params.append(source)
    cursor.execute(query, params)

    rows = []
    for record in cursor.fetchall():
        row = {k: v for k, v in zip(select_cols, record)}
        row.setdefault("segment_key", None)
        row.setdefault("segment", None)
        rows.append(row)
    return rows


def fetch_table_column(conn, table_name, column_name):
    cursor = conn.cursor()
    cursor.execute(f"SELECT {column_name} FROM {table_name}")
    return [r[0] for r in cursor.fetchall()]


def table_exists(conn, table_name):
    cursor = conn.cursor()
    cursor.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1",
        (table_name,),
    )
    return cursor.fetchone() is not None


def main():
    parser = argparse.ArgumentParser(description="Validate SQLite DB contents against tree and plot the tree")
    parser.add_argument("--db", required=True, help="Path to SQLite DB")
    parser.add_argument("--outdir", required=True, help="Output directory for report and plot")
    parser.add_argument(
        "--test-mode",
        action="store_true",
        help="Relax strict failures for known test-only edge cases (e.g., no placeable queries)",
    )
    parser.add_argument(
        "--expect-segment-trees",
        action="store_true",
        help="In segmented test mode, require at least one USHER tree per segment in meta_data",
    )
    parser.add_argument(
        "--check-update-integrity",
        action="store_true",
        help="Validate update audit tables and update-mode integrity constraints",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    conn = sqlite3.connect(args.db)
    try:
        tree_name, newick = fetch_tree_newick(conn)
        if not newick:
            raise SystemExit("No tree found in DB (trees table is empty)")

        tree = Phylo.read(StringIO(newick), "newick")
        tree_terminals = {t.name for t in tree.get_terminals() if t.name}

        meta_accessions = [x for x in fetch_table_column(conn, "meta_data", "primary_accession") if x]
        seq_headers = [x for x in fetch_table_column(conn, "sequences", "header") if x]

        meta_set = set(meta_accessions)
        seq_set = set(seq_headers)
        excluded_accessions = set()
        if table_exists(conn, "excluded_accessions"):
            excluded_accessions = set(
                x for x in fetch_table_column(conn, "excluded_accessions", "primary_accession") if x
            )
        expected_meta_set = meta_set - excluded_accessions
        missing_in_sequences = sorted(expected_meta_set - seq_set)
        missing_in_meta = sorted(tree_terminals - meta_set)

        meta_columns = get_table_columns(conn, "meta_data")
        cluster_col = find_cluster_column(meta_columns)
        centroid_set = set()
        missing_centroids_in_tree = []
        extra_in_tree = []
        meta_expected_in_tree = set(expected_meta_set)
        if cluster_col:
            cursor = conn.cursor()
            cursor.execute(f"SELECT primary_accession, {cluster_col} FROM meta_data")
            cluster_rows = cursor.fetchall()
            centroid_set = {
                str(cluster_val).strip()
                for _, cluster_val in cluster_rows
                if not is_cluster_placeholder(cluster_val)
            }
            missing_centroids_in_tree = sorted(centroid_set - tree_terminals)
            extra_in_tree = sorted(tree_terminals - centroid_set)

        missing_in_tree = sorted(meta_expected_in_tree - tree_terminals)

        # Basic completeness checks
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM meta_data")
        meta_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM sequences")
        seq_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM trees WHERE source='usher'")
        usher_tree_count = cursor.fetchone()[0]

        segment_count = None
        meta_columns_for_seg = get_table_columns(conn, "meta_data")
        if "segment" in meta_columns_for_seg:
            segment_values = [
                x for x in fetch_table_column(conn, "meta_data", "segment")
                if x is not None and str(x).strip() != ""
            ]
            segment_count = len(set(segment_values))

        segmented_validation_ok = False
        if args.expect_segment_trees and segment_count is not None and segment_count > 1:
            usher_trees = fetch_trees(conn, source="usher")

            excluded_accessions = set()
            if table_exists(conn, "excluded_accessions"):
                excluded_accessions = set(
                    x for x in fetch_table_column(conn, "excluded_accessions", "primary_accession") if x
                )

            cursor.execute("SELECT primary_accession, segment FROM meta_data")
            meta_rows = cursor.fetchall()

            accession_to_segment = {}
            accessions_by_segment = {}
            segments_in_meta = set()
            for acc, seg in meta_rows:
                if acc:
                    accession_to_segment[str(acc)] = str(seg).strip() if seg is not None else None
                if seg is None:
                    continue
                seg_str = str(seg).strip()
                if not seg_str:
                    continue
                segments_in_meta.add(seg_str)
                if acc and str(acc) not in excluded_accessions:
                    accessions_by_segment.setdefault(seg_str, set()).add(str(acc))

            trees_by_segment = {}
            for tr in usher_trees:
                seg = tr.get("segment")
                seg_key = tr.get("segment_key")
                if (seg is None or str(seg).strip() == "") and seg_key:
                    seg = accession_to_segment.get(str(seg_key))
                if seg is None or str(seg).strip() == "":
                    continue
                seg_str = str(seg).strip()
                t = Phylo.read(StringIO(tr["newick"]), "newick")
                terms = {x.name for x in t.get_terminals() if x.name}
                trees_by_segment.setdefault(seg_str, set()).update(terms)

            missing_segment_trees = sorted(s for s in segments_in_meta if s not in trees_by_segment)
            per_segment_missing = {}
            for seg in sorted(segments_in_meta):
                expected_accessions = accessions_by_segment.get(seg, set())
                tree_terms = trees_by_segment.get(seg, set())
                m = sorted(expected_accessions - tree_terms)
                if m:
                    per_segment_missing[seg] = m

            print(f"[info] Segmented validation: meta segments={len(segments_in_meta)} tree segments={len(trees_by_segment)}")
            if missing_segment_trees:
                print(f"[info] Missing tree segments: {', '.join(missing_segment_trees[:10])}")
            if per_segment_missing:
                first_seg = sorted(per_segment_missing.keys())[0]
                print(f"[info] Segment {first_seg} missing accessions in tree (first 10): {', '.join(per_segment_missing[first_seg][:10])}")

            if missing_segment_trees:
                raise SystemExit(
                    "Validation failed: segmented run expects at least one UShER tree per segment "
                    f"(missing segments={', '.join(missing_segment_trees[:10])})"
                )
            if per_segment_missing:
                raise SystemExit("Validation failed: per-segment UShER tree is missing unfiltered accessions")

            segmented_validation_ok = True

        report_path = os.path.join(args.outdir, "db_tree_validation.txt")
        missing_tree_path = os.path.join(args.outdir, "missing_in_tree.txt")
        missing_seq_path = os.path.join(args.outdir, "missing_in_sequences.txt")
        missing_meta_path = os.path.join(args.outdir, "missing_in_meta.txt")
        missing_centroids_path = os.path.join(args.outdir, "missing_centroids_in_tree.txt")
        extra_in_tree_path = os.path.join(args.outdir, "extra_in_tree.txt")

        with open(report_path, "w", encoding="utf-8") as report:
            report.write(f"Tree source: {tree_name}\n")
            report.write(f"Meta rows: {meta_count}\n")
            report.write(f"Sequence rows: {seq_count}\n")
            report.write(f"Tree terminals: {len(tree_terminals)}\n")
            report.write(f"UShER tree rows: {usher_tree_count}\n")
            if segment_count is not None:
                report.write(f"Segment count in meta_data: {segment_count}\n")
            if cluster_col:
                report.write(f"Cluster column: {cluster_col}\n")
                report.write(f"Centroid count: {len(centroid_set)}\n")
                report.write(f"Missing centroids in tree: {len(missing_centroids_in_tree)}\n")
                report.write(f"Extra nodes in tree: {len(extra_in_tree)}\n")
            report.write(f"Missing in tree (meta_data -> tree): {len(missing_in_tree)}\n")
            report.write(f"Missing in sequences (meta_data -> sequences): {len(missing_in_sequences)}\n")
            report.write(f"Missing in meta_data (tree -> meta_data): {len(missing_in_meta)}\n")

            if missing_in_tree:
                report.write("\nMissing in tree (first 50):\n")
                report.write("\n".join(missing_in_tree[:50]) + "\n")
            if missing_in_sequences:
                report.write("\nMissing in sequences (first 50):\n")
                report.write("\n".join(missing_in_sequences[:50]) + "\n")
            if missing_in_meta:
                report.write("\nMissing in meta_data (first 50):\n")
                report.write("\n".join(missing_in_meta[:50]) + "\n")

            if cluster_col and missing_centroids_in_tree:
                report.write("\nMissing centroids in tree (first 50):\n")
                report.write("\n".join(missing_centroids_in_tree[:50]) + "\n")
            if cluster_col and extra_in_tree:
                report.write("\nExtra nodes in tree (first 50):\n")
                report.write("\n".join(extra_in_tree[:50]) + "\n")

            # Table-level coverage summary
            report.write("\nTable-level accession coverage:\n")
            for table_name in ["meta_data", "sequences", "sequence_alignment", "features", "insertions", "host_taxa"]:
                cols = get_table_columns(conn, table_name)
                acc_cols = [c for c in cols if c in ACCESSION_COLUMNS]
                if not acc_cols:
                    report.write(f"- {table_name}: no accession-like columns found\n")
                    continue
                values = collect_table_values(conn, table_name, acc_cols)
                missing_from_table = sorted(tree_terminals - values)
                report.write(
                    f"- {table_name}: columns={','.join(acc_cols)} values={len(values)} missing_from_table={len(missing_from_table)}\n"
                )

            if args.check_update_integrity:
                report.write("\nUpdate integrity checks:\n")
                if table_exists(conn, "update_batches") and table_exists(conn, "update_table_deltas"):
                    cursor.execute("SELECT batch_id, mode, started_at, finished_at FROM update_batches ORDER BY rowid DESC LIMIT 1")
                    last_batch = cursor.fetchone()
                    if last_batch:
                        batch_id, mode, started_at, finished_at = last_batch
                        report.write(f"- last_batch_id: {batch_id}\n")
                        report.write(f"- last_batch_mode: {mode}\n")
                        report.write(f"- last_batch_started: {started_at}\n")
                        report.write(f"- last_batch_finished: {finished_at}\n")
                        cursor.execute(
                            "SELECT table_name, before_count, after_count, delta FROM update_table_deltas WHERE batch_id=? ORDER BY table_name",
                            (batch_id,),
                        )
                        rows = cursor.fetchall()
                        for table_name, before_count, after_count, delta in rows:
                            report.write(f"  - {table_name}: before={before_count} after={after_count} delta={delta}\n")
                else:
                    report.write("- update audit tables missing (update_batches/update_table_deltas)\n")

                if table_exists(conn, "features") and table_exists(conn, "meta_data"):
                    feature_cols = get_table_columns(conn, "features")
                    meta_cols = get_table_columns(conn, "meta_data")
                    if "segment" in feature_cols and "segment" in meta_cols and "accession" in feature_cols:
                        cursor.execute(
                            """
                            SELECT COUNT(*)
                            FROM features f
                            JOIN meta_data m ON m.primary_accession = f.accession
                            WHERE COALESCE(TRIM(f.segment), '') != COALESCE(TRIM(m.segment), '')
                            """
                        )
                        seg_mismatch_count = cursor.fetchone()[0]
                        report.write(f"- feature_segment_mismatch: {seg_mismatch_count}\n")

        # Write full lists for debugging
        if missing_in_tree:
            with open(missing_tree_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_tree) + "\n")
        if missing_in_sequences:
            with open(missing_seq_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_sequences) + "\n")
        if missing_in_meta:
            with open(missing_meta_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_meta) + "\n")

        if cluster_col and missing_centroids_in_tree:
            with open(missing_centroids_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_centroids_in_tree) + "\n")
        if cluster_col and extra_in_tree:
            with open(extra_in_tree_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(extra_in_tree) + "\n")

        # Print summary to stdout for quick visibility in logs
        print(f"[info] Tree source: {tree_name}")
        print(f"[info] Meta rows: {meta_count}, Sequence rows: {seq_count}, Tree terminals: {len(tree_terminals)}")
        print(f"[info] UShER tree rows: {usher_tree_count}")
        if segment_count is not None:
            print(f"[info] Segment count in meta_data: {segment_count}")
        if cluster_col:
            print(f"[info] Cluster column: {cluster_col}, Centroid count: {len(centroid_set)}")
            print(f"[info] Missing centroids in tree: {len(missing_centroids_in_tree)}")
            print(f"[info] Extra nodes in tree: {len(extra_in_tree)}")
        print(f"[info] Missing in tree: {len(missing_in_tree)}")
        print(f"[info] Missing in sequences: {len(missing_in_sequences)}")
        print(f"[info] Missing in meta_data: {len(missing_in_meta)}")
        if missing_in_tree:
            print(f"[info] Missing in tree (first 10): {', '.join(missing_in_tree[:10])}")
        if missing_in_sequences:
            print(f"[info] Missing in sequences (first 10): {', '.join(missing_in_sequences[:10])}")
        if missing_in_meta:
            print(f"[info] Missing in meta_data (first 10): {', '.join(missing_in_meta[:10])}")
        if cluster_col and missing_centroids_in_tree:
            print(f"[info] Missing centroids in tree (first 10): {', '.join(missing_centroids_in_tree[:10])}")
        if cluster_col and extra_in_tree:
            print(f"[info] Extra nodes in tree (first 10): {', '.join(extra_in_tree[:10])}")

        # Plot tree
        fig = plt.figure(figsize=(12, 18))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, do_show=False, axes=ax)
        plot_path = os.path.join(args.outdir, "db_tree.png")
        plt.tight_layout()
        fig.savefig(plot_path, dpi=150)

        overlap_count = len(meta_set & tree_terminals)
        disjoint_sets = overlap_count == 0 and len(meta_set) > 0 and len(tree_terminals) > 0

        # Test-mode exception: disjoint USHER tree/meta can happen when no query sequences were placeable
        # and the resulting tree contains only reference nodes.
        if args.test_mode and tree_name == "usher" and disjoint_sets:
            print(
                "[warn] Test mode: tree terminals and meta_data are disjoint; "
                "treating as no-placeable-queries scenario and not failing validation."
            )
            return

        if args.expect_segment_trees and segment_count is not None and segment_count > 1 and not segmented_validation_ok:
            if usher_tree_count < segment_count:
                raise SystemExit(
                    "Validation failed: segmented run expects at least one UShER tree per segment "
                    f"(segments={segment_count}, usher_trees={usher_tree_count})"
                )

        if args.check_update_integrity:
            if not table_exists(conn, "update_batches") or not table_exists(conn, "update_table_deltas"):
                raise SystemExit("Validation failed: update audit tables are missing")
            cursor.execute("SELECT batch_id FROM update_batches ORDER BY rowid DESC LIMIT 1")
            row = cursor.fetchone()
            if row:
                batch_id = row[0]
                cursor.execute("SELECT COUNT(*) FROM update_table_deltas WHERE batch_id=?", (batch_id,))
                if cursor.fetchone()[0] == 0:
                    raise SystemExit("Validation failed: no update_table_deltas rows found for latest batch")
            if table_exists(conn, "features") and table_exists(conn, "meta_data"):
                feature_cols = get_table_columns(conn, "features")
                meta_cols = get_table_columns(conn, "meta_data")
                if "segment" in feature_cols and "segment" in meta_cols and "accession" in feature_cols:
                    cursor.execute(
                        """
                        SELECT COUNT(*)
                        FROM features f
                        JOIN meta_data m ON m.primary_accession = f.accession
                        WHERE COALESCE(TRIM(f.segment), '') != COALESCE(TRIM(m.segment), '')
                        """
                    )
                    if cursor.fetchone()[0] > 0:
                        raise SystemExit("Validation failed: segment contamination detected in features table")

        # Fail if any validation issues detected
        if tree_name == "usher":
            if segmented_validation_ok:
                return
            if missing_in_tree or missing_in_sequences or missing_in_meta:
                raise SystemExit("Validation failed: UShER tree does not match accessions")
            if cluster_col and extra_in_tree:
                print(
                    "[info] UShER tree contains non-centroid terminals "
                    f"({len(extra_in_tree)} extra vs centroid context); full accession coverage passed."
                )
        elif cluster_col:
            if missing_centroids_in_tree or extra_in_tree:
                raise SystemExit("Validation failed: IQ-TREE does not match centroid set")
        else:
            if missing_in_tree or missing_in_sequences or missing_in_meta:
                raise SystemExit("Validation failed: missing accessions detected")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
