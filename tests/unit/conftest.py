import sqlite3
from pathlib import Path

import pytest


@pytest.fixture()
def basic_update_db(tmp_path: Path) -> Path:
    db_path = tmp_path / "basic_update_seed.db"
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
            [
                ("REF1", "master", "1"),
                ("REF2", "reference", "2"),
                ("Q_OLD", "query", "1"),
            ],
        )

        cur.execute("CREATE TABLE excluded_accessions (primary_accession TEXT, reason TEXT)")
        cur.execute("INSERT INTO excluded_accessions(primary_accession, reason) VALUES ('Q_EXCL', 'manual')")

        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, alignment_name, alignment, segment) VALUES (?, ?, ?, ?)",
            [
                ("REF1", "REF1", "ATGC", "1"),
                ("REF2", "REF2", "A-GC", "2"),
                ("Q_OLD", "REF1", "AT-T", "1"),
            ],
        )

        cur.execute(
            "CREATE TABLE features (accession TEXT, master_ref_accession TEXT, reference_accession TEXT, aln_start TEXT, aln_end TEXT, cds_start TEXT, cds_end TEXT, product TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO features(accession, master_ref_accession, reference_accession, aln_start, aln_end, cds_start, cds_end, product, segment) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("Q_OLD", "REF1", "REF1", "1", "4", "1", "4", "P", "1"),
            ],
        )

        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.executemany(
            "INSERT INTO sequences(header, sequence) VALUES (?, ?)",
            [("REF1", "ATGC"), ("REF2", "AAGC"), ("Q_OLD", "ATTT")],
        )

        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, segment_key TEXT, newick TEXT)")
        cur.executemany(
            "INSERT INTO trees(name, source, segment, segment_key, newick) VALUES (?, ?, ?, ?, ?)",
            [
                ("usher_seg1", "usher", "1", "seg_1", "(REF1:0.1,REF2:0.1,Q_OLD:0.2);"),
                ("usher_seg2", "usher", "2", "seg_2", "(REF2:0.1);"),
                ("iqtree_seg2", "iqtree", "2", "seg_2", "(REF2:0.1);"),
            ],
        )

        cur.execute("CREATE TABLE update_batches (batch_id TEXT, mode TEXT, started_at TEXT, finished_at TEXT)")
        cur.execute(
            "INSERT INTO update_batches(batch_id, mode, started_at, finished_at) VALUES ('batch_seed', 'update', '2026-01-01', '2026-01-01')"
        )
        cur.execute(
            "CREATE TABLE update_table_deltas (batch_id TEXT, table_name TEXT, before_count INTEGER, after_count INTEGER, delta INTEGER)"
        )
        cur.execute(
            "INSERT INTO update_table_deltas(batch_id, table_name, before_count, after_count, delta) VALUES ('batch_seed', 'meta_data', 2, 3, 1)"
        )

        conn.commit()
    finally:
        conn.close()

    return db_path