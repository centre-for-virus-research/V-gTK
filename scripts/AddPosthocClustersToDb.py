#!/usr/bin/env python3
"""Run MMseqs2 clustering post-hoc on a finished RABV-gTK alignment and attach cluster_<N>pct to SQLite meta_data.

Intended for large treeless builds (e.g. SARS-CoV-2 with `tree_free = true` and
`skip_clustering = true`) where clustering is deferred until after the primary
database build completes.
"""

import argparse
import os
import sqlite3
import sys
import tempfile

from accession_utils import normalise_accession
from ExportRefListFromUpdateDb import load_master_accessions_from_file, load_reference_dict
from MMseqsClustering import run_mmseqs_clustering


def _resolve_master_ids(db_path, ref_list_path=None):
    if ref_list_path and os.path.isfile(ref_list_path):
        return load_master_accessions_from_file(ref_list_path)
    ref_dict = load_reference_dict(db_path)
    return [acc for acc, rtype in ref_dict.items() if str(rtype).strip().lower() == "master"]


def _load_cluster_tsv(tsv_path):
    cluster_map = {}
    with open(tsv_path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0].strip(), parts[1].strip()
            if not rep or not member:
                continue
            cluster_map[member] = rep
            bare = normalise_accession(member)
            if bare and bare not in cluster_map:
                cluster_map[bare] = rep
    return cluster_map


def attach_clusters_to_db(db_path, cluster_tsv, min_seq_id):
    col_name = f"cluster_{int(round(float(min_seq_id) * 100))}pct"
    cluster_map = _load_cluster_tsv(cluster_tsv)
    if not cluster_map:
        raise ValueError(f"No cluster assignments found in {cluster_tsv}")

    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()
        cur.execute("PRAGMA temp_store = MEMORY;")
        cur.execute("PRAGMA cache_size = -67108864;")
        cur.execute("PRAGMA mmap_size = 274877906944;")
        cols = [row[1] for row in cur.execute("PRAGMA table_info(meta_data)").fetchall()]
        if not cols:
            raise ValueError(f"Database {db_path} has no meta_data table")
        if col_name not in cols:
            cur.execute(f'ALTER TABLE meta_data ADD COLUMN "{col_name}" TEXT')

        cur.execute("CREATE TEMP TABLE _tmp_posthoc_clusters (acc TEXT PRIMARY KEY, rep TEXT)")
        cur.executemany(
            "INSERT OR REPLACE INTO _tmp_posthoc_clusters (acc, rep) VALUES (?, ?)",
            list(cluster_map.items()),
        )
        cur.execute(
            f'''
            UPDATE meta_data
            SET "{col_name}" = (
                SELECT rep FROM _tmp_posthoc_clusters
                WHERE _tmp_posthoc_clusters.acc = meta_data.primary_accession
            )
            WHERE EXISTS (
                SELECT 1 FROM _tmp_posthoc_clusters
                WHERE _tmp_posthoc_clusters.acc = meta_data.primary_accession
            )
            '''
        )
        updated = cur.rowcount
        conn.commit()
        return col_name, len(cluster_map), updated
    finally:
        conn.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--db", required=True, help="Path to target RABV-gTK SQLite .db file")
    parser.add_argument("--msa", required=True, help="Path to deduplicated merged MSA FASTA (*_merged_MSA_dedup.fasta)")
    parser.add_argument("--ref-list", default=None, help="Optional reference list TSV (otherwise master IDs are read from --db)")
    parser.add_argument("--work-dir", default=None, help="Directory for MMseqs intermediate/output files (default: temporary dir)")
    parser.add_argument("--min-seq-id", type=float, default=0.999, help="MMseqs minimum sequence identity (default: 0.999)")
    parser.add_argument("--min-completeness", type=float, default=0.98, help="Minimum completeness for step-1 representatives (default: 0.98)")
    parser.add_argument("--threads", type=int, default=48, help="Threads for MMseqs2 (default: 48)")
    parser.add_argument("--max-memory", default="1000G",
                        help="MMseqs2 --split-memory-limit cap (default: 1000G = ~1 TB)")
    parser.add_argument("--exact", action="store_true",
                        help="Use full all-vs-all 'mmseqs cluster' instead of default linear-time 'mmseqs linclust'")
    args = parser.parse_args(argv)

    if not os.path.isfile(args.db):
        raise FileNotFoundError(f"SQLite DB not found: {args.db}")
    if not os.path.isfile(args.msa):
        raise FileNotFoundError(f"MSA FASTA not found: {args.msa}")

    master_ids = _resolve_master_ids(args.db, args.ref_list)
    base_name = os.path.splitext(os.path.basename(args.msa))[0]

    def _run_in_dir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
        run_mmseqs_clustering(
            args.msa,
            out_dir,
            args.min_seq_id,
            args.threads,
            strip_gaps=True,
            sort_by_quality=True,
            two_step=True,
            min_completeness=args.min_completeness,
            fast=not args.exact,
            reference_ids=master_ids,
            split_memory_limit=args.max_memory,
        )
        cluster_tsv = os.path.join(out_dir, base_name, f"{base_name}_clusters.tsv")
        col_name, n_members, n_updated = attach_clusters_to_db(args.db, cluster_tsv, args.min_seq_id)
        print(f"[posthoc-mmseqs] Populated meta_data.{col_name} in {args.db}: "
              f"{n_members:,} clustered sequences -> {n_updated:,} meta_data rows updated.")

    if args.work_dir:
        _run_in_dir(args.work_dir)
    else:
        with tempfile.TemporaryDirectory(prefix="posthoc_mmseqs_") as tmp:
            _run_in_dir(tmp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
