import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "ValidateDbTree.py"


def _add_required_alignment_tables(db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS insertions (primary_accession TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS host_taxa (taxonomy_id TEXT)")
        cur.execute("INSERT INTO insertions(primary_accession) VALUES ('REF1')")
        cur.execute("INSERT INTO host_taxa(taxonomy_id) VALUES ('111')")
        cur.execute("INSERT INTO host_taxa(taxonomy_id) VALUES ('222')")
        conn.commit()
    finally:
        conn.close()


def _add_mutation_tables(db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE sequence_mutations (
                primary_accession TEXT,
                mutation_id TEXT,
                protein_name TEXT,
                segment TEXT,
                aa_position INTEGER,
                alt_residue TEXT,
                combination_id TEXT
            )
            """
        )
        cur.executemany(
            "INSERT INTO sequence_mutations VALUES (?, ?, ?, ?, ?, ?, ?)",
            [
                ("REF1", "NS3:107I", "NS3", "1", 107, "I", "comboA"),
                ("REF1", "NS3:155K", "NS3", "1", 155, "K", "comboA"),
                ("Q_OLD", "NS5A:30Y", "NS5A", "1", 30, "Y", ""),
            ],
        )
        cur.execute(
            """
            CREATE TABLE sequence_drug_resistance (
                primary_accession TEXT,
                combination_id TEXT,
                combination_status TEXT,
                mutations_detected INTEGER,
                mutations_required INTEGER,
                resistance_category TEXT,
                drug TEXT
            )
            """
        )
        cur.execute(
            "INSERT INTO sequence_drug_resistance VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("REF1", "comboA", "complete", 2, 2, "III", "drugA"),
        )
        conn.commit()
    finally:
        conn.close()


def _add_compact_mutation_tables(db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE mutation_catalog (
                mutation_id TEXT,
                protein_name TEXT,
                segment TEXT,
                aa_position TEXT,
                alt_residue TEXT,
                reference_accession TEXT,
                mutation_type TEXT,
                signature_id TEXT,
                signature_kind TEXT,
                combination_id TEXT,
                combination_size TEXT,
                resistance_category TEXT,
                drug TEXT
            )
            """
        )
        cur.executemany(
            "INSERT INTO mutation_catalog VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("NS3:107I", "NS3", "1", "107", "I", "REF1", "snp", "sig_ns3_107i", "single", "", "", "", ""),
                ("comboA_part1", "NS3", "1", "107", "I", "REF1", "snp", "comboA_sig", "combination", "comboA", "2", "III", "drugA"),
                ("comboA_part2", "NS3", "1", "155", "K", "REF1", "snp", "comboA_sig", "combination", "comboA", "2", "III", "drugA"),
            ],
        )
        cur.execute(
            """
            CREATE TABLE sequence_relevant_mutation_summary (
                primary_accession TEXT,
                relevant_mutations_present TEXT,
                total_relevant_mutation_count INTEGER
            )
            """
        )
        cur.executemany(
            "INSERT INTO sequence_relevant_mutation_summary VALUES (?, ?, ?)",
            [
                ("REF1", "NS3:107I;comboA_part1;comboA_part2", 3),
                ("Q_OLD", "NS3:107I", 1),
            ],
        )
        cur.execute(
            """
            CREATE TABLE completed_signatures_only (
                primary_accession TEXT,
                signature_id TEXT,
                signature_kind TEXT
            )
            """
        )
        cur.executemany(
            "INSERT INTO completed_signatures_only VALUES (?, ?, ?)",
            [
                ("REF1", "sig_ns3_107i", "single"),
                ("REF1", "comboA_sig", "combination"),
                ("Q_OLD", "sig_ns3_107i", "single"),
            ],
        )
        conn.commit()
    finally:
        conn.close()


def _create_segmented_subsample_db(db_path: Path, tree_source: str = "usher"):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, segment TEXT, exclusion_status TEXT, exclusion_criteria TEXT)"
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, segment, exclusion_status, exclusion_criteria) VALUES (?, ?, ?, ?)",
            [
                ("SEG1_A", "1", "", ""),
                ("SEG1_B", "1", "", ""),
                ("SEG2_A", "2", "", ""),
                ("SEG2_B", "2", "", ""),
            ],
        )
        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.executemany(
            "INSERT INTO sequences(header, sequence) VALUES (?, ?)",
            [("SEG1_A", "ATGC"), ("SEG1_B", "ATGT"), ("SEG2_A", "ATGA"), ("SEG2_B", "ATGG")],
        )
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, segment) VALUES (?, ?, ?, ?)",
            [
                ("SEG1_A", "SEG1_A", "SEG1_A", "1"),
                ("SEG1_B", "SEG1_B", "SEG1_B", "1"),
                ("SEG2_A", "SEG2_A", "SEG2_A", "2"),
                ("SEG2_B", "SEG2_B", "SEG2_B", "2"),
            ],
        )
        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
        cur.executemany(
            "INSERT INTO trees(name, source, newick, segment, segment_key) VALUES (?, ?, ?, ?, ?)",
            [
                (f"{tree_source}_seg1", tree_source, "(SEG1_A:0.1);", "1", "refset_1"),
                (f"{tree_source}_seg2", tree_source, "(SEG2_A:0.1);", "2", "refset_2"),
            ],
        )
        cur.execute("CREATE TABLE update_batches (batch_id TEXT, mode TEXT, started_at TEXT, finished_at TEXT)")
        cur.execute("INSERT INTO update_batches(batch_id, mode, started_at, finished_at) VALUES ('b1', 'update', '2026-01-01', '2026-01-01')")
        cur.execute("CREATE TABLE update_table_deltas (batch_id TEXT, table_name TEXT, before_count INTEGER, after_count INTEGER, delta INTEGER)")
        cur.execute("INSERT INTO update_table_deltas(batch_id, table_name, before_count, after_count, delta) VALUES ('b1', 'meta_data', 0, 4, 4)")
        cur.execute("CREATE TABLE features (accession TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO features(accession, segment) VALUES (?, ?)",
            [("SEG1_A", "1"), ("SEG1_B", "1"), ("SEG2_A", "2"), ("SEG2_B", "2")],
        )
        cur.execute("CREATE TABLE insertions (primary_accession TEXT)")
        cur.executemany(
            "INSERT INTO insertions(primary_accession) VALUES (?)",
            [("SEG1_A",), ("SEG1_B",), ("SEG2_A",), ("SEG2_B",)],
        )
        cur.execute("CREATE TABLE host_taxa (primary_accession TEXT)")
        cur.executemany(
            "INSERT INTO host_taxa(primary_accession) VALUES (?)",
            [("SEG1_A",), ("SEG1_B",), ("SEG2_A",), ("SEG2_B",)],
        )
        conn.commit()
    finally:
        conn.close()


def _create_segmented_iqtree_centroid_db(db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT, cluster_98pct TEXT, exclusion_status TEXT, exclusion_criteria TEXT)"
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment, cluster_98pct, exclusion_status, exclusion_criteria) VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("SEG1_A", "query", "1", "SEG1_A", "", ""),
                ("SEG1_B", "query", "1", "SEG1_A", "", ""),
                ("SEG1_REF", "reference", "1", "", "1", "reference_context"),
                ("SEG2_A", "query", "2", "SEG2_A", "", ""),
                ("SEG2_B", "query", "2", "SEG2_A", "", ""),
                ("SEG2_MASTER", "master", "2", "", "1", "reference_context"),
            ],
        )
        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.executemany(
            "INSERT INTO sequences(header, sequence) VALUES (?, ?)",
            [("SEG1_A", "ATGC"), ("SEG1_B", "ATGT"), ("SEG1_REF", "ATGA"), ("SEG2_A", "ATGA"), ("SEG2_B", "ATGG"), ("SEG2_MASTER", "ATGC")],
        )
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, segment) VALUES (?, ?, ?, ?)",
            [
                ("SEG1_A", "SEG1_A", "SEG1_A", "1"),
                ("SEG1_B", "SEG1_B", "SEG1_B", "1"),
                ("SEG2_A", "SEG2_A", "SEG2_A", "2"),
                ("SEG2_B", "SEG2_B", "SEG2_B", "2"),
            ],
        )
        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
        cur.executemany(
            "INSERT INTO trees(name, source, newick, segment, segment_key) VALUES (?, ?, ?, ?, ?)",
            [
                ("iqtree_seg1", "iqtree", "(SEG1_A:0.1,SEG1_REF:0.2);", "1", "refset_1"),
                ("iqtree_seg2", "iqtree", "(SEG2_A:0.1,SEG2_MASTER:0.2);", "2", "refset_2"),
            ],
        )
        cur.execute("CREATE TABLE update_batches (batch_id TEXT, mode TEXT, started_at TEXT, finished_at TEXT)")
        cur.execute("INSERT INTO update_batches(batch_id, mode, started_at, finished_at) VALUES ('b1', 'update', '2026-01-01', '2026-01-01')")
        cur.execute("CREATE TABLE update_table_deltas (batch_id TEXT, table_name TEXT, before_count INTEGER, after_count INTEGER, delta INTEGER)")
        cur.execute("INSERT INTO update_table_deltas(batch_id, table_name, before_count, after_count, delta) VALUES ('b1', 'meta_data', 0, 4, 4)")
        cur.execute("CREATE TABLE features (accession TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO features(accession, segment) VALUES (?, ?)",
            [("SEG1_A", "1"), ("SEG1_B", "1"), ("SEG2_A", "2"), ("SEG2_B", "2")],
        )
        cur.execute("CREATE TABLE insertions (primary_accession TEXT)")
        cur.executemany(
            "INSERT INTO insertions(primary_accession) VALUES (?)",
            [("SEG1_A",), ("SEG1_B",), ("SEG2_A",), ("SEG2_B",)],
        )
        cur.execute("CREATE TABLE host_taxa (primary_accession TEXT)")
        cur.executemany(
            "INSERT INTO host_taxa(primary_accession) VALUES (?)",
            [("SEG1_A",), ("SEG1_B",), ("SEG2_A",), ("SEG2_B",)],
        )
        conn.commit()
    finally:
        conn.close()


def test_validate_db_tree_update_integrity_passes_on_seed_db(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    _add_mutation_tables(basic_update_db)
    outdir = tmp_path / "out"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(outdir),
            "--check-update-integrity",
                "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "[sequence_alignment vs meta_data]" in report
    assert "[features vs meta_data]" in report
    assert "[host_taxa vs meta_data]" in report
    assert "Accepted distinct accessions:" in report
    assert "Filtered distinct accessions:" in report
    assert "Update integrity checks:" in report
    assert "last_batch_id:" in report
    assert "[feature projection integrity]" in report
    assert "[sequence_mutations integrity]" in report
    assert "[sequence_drug_resistance integrity]" in report


def test_validate_db_tree_update_integrity_allows_host_taxa_used_only_by_excluded_rows(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)

    conn = sqlite3.connect(str(basic_update_db))
    try:
        cur = conn.cursor()
        cur.execute("UPDATE meta_data SET host_taxa_id='9796' WHERE primary_accession='Q_EXCL'")
        cur.execute("INSERT INTO host_taxa(taxonomy_id) VALUES ('9796')")
        conn.commit()
    finally:
        conn.close()

    outdir = tmp_path / "out_host_taxa_excluded"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(outdir),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "[host_taxa vs meta_data]" in report
    assert "Table: host_taxa  Column: taxonomy_id" in report
    assert "Expected: 3  Observed: 3" in report
    assert "Extra: 0" in report


def test_validate_db_tree_update_integrity_fails_for_out_of_span_feature_projection(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)

    conn = sqlite3.connect(str(basic_update_db))
    try:
        cur = conn.cursor()
        cur.execute("DELETE FROM features WHERE accession='REF1'")
        cur.execute("DELETE FROM features WHERE accession='Q_OLD'")
        cur.executemany(
            "INSERT INTO features(accession, master_ref_accession, reference_accession, aln_start, aln_end, cds_start, cds_end, cds_start_OG_seq, cds_end_OG_seq, product, segment) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("REF1", "REF1", "REF1", "1", "4", "1", "2", "1", "2", "P1", "1"),
                ("REF1", "REF1", "REF1", "1", "4", "3", "4", "3", "4", "P2", "1"),
                ("Q_OLD", "REF1", "REF1", "1", "2", "1", "2", "1", "2", "P1", "1"),
                ("Q_OLD", "REF1", "REF1", "1", "2", "2", "2", "2", "2", "P2", "1"),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(outdir),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "feature projection integrity check failed" in result.stderr

    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "[feature projection integrity]" in report
    assert "offending_rows: 1" in report
    assert "mismatched_cds_rows: 1" in report
    assert "accession=Q_OLD product=P2" in report
    assert "reason=cds_span_mismatch" in report


def test_validate_db_tree_update_integrity_fails_for_clipped_cds_projection_with_overlapping_alignment(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)

    conn = sqlite3.connect(str(basic_update_db))
    try:
        cur = conn.cursor()
        cur.execute("DELETE FROM features WHERE accession='REF1'")
        cur.execute("DELETE FROM features WHERE accession='Q_OLD'")
        cur.executemany(
            "INSERT INTO features(accession, master_ref_accession, reference_accession, aln_start, aln_end, cds_start, cds_end, cds_start_OG_seq, cds_end_OG_seq, product, segment) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("REF1", "REF1", "REF1", "1", "4", "1", "2", "1", "2", "P1", "1"),
                ("REF1", "REF1", "REF1", "1", "4", "3", "4", "3", "4", "P2", "1"),
                # Alignment span still overlaps P1 on REF1, but cds span was incorrectly clipped to the covered region.
                ("Q_OLD", "REF1", "REF1", "1", "2", "1", "1", "1", "1", "P1", "1"),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    outdir = tmp_path / "out_clipped"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(outdir),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "feature projection integrity check failed" in result.stderr

    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "[feature projection integrity]" in report
    assert "offending_rows: 1" in report
    assert "mismatched_cds_rows: 1" in report
    assert "accession=Q_OLD product=P1" in report
    assert "cds_start=1 cds_end=1" in report
    assert "reason=cds_span_mismatch" in report


def test_validate_db_tree_update_integrity_passes_with_compact_mutation_tables(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    _add_compact_mutation_tables(basic_update_db)
    outdir = tmp_path / "out"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(outdir),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "[sequence_relevant_mutation_summary integrity]" in report
    assert "[completed_signatures_only integrity]" in report


def test_validate_db_tree_fails_when_features_do_not_match_expected_accessions(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DELETE FROM features WHERE accession='REF2'")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "DB consistency checks failed for features vs meta_data" in result.stderr


def test_validate_db_tree_fails_without_trees_by_default(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DELETE FROM trees")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "No tree found in DB (trees table is empty)" in result.stderr


def test_validate_db_tree_allows_tree_free_db_with_flag(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DELETE FROM trees")
        conn.commit()
    finally:
        conn.close()

    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(outdir),
            "--allow-no-trees",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "Tree status: not present (allowed by --allow-no-trees)" in report
    assert (outdir / "db_tree.png").exists()


def test_validate_db_tree_tree_free_mode_still_checks_consistency(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DELETE FROM trees")
        conn.execute("DELETE FROM features WHERE accession='REF2'")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--allow-no-trees",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "DB consistency checks failed for features vs meta_data" in result.stderr


def test_validate_db_tree_update_integrity_fails_without_audit_tables(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DROP TABLE update_batches")
        conn.execute("DROP TABLE update_table_deltas")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "update audit tables are missing" in result.stderr


def test_validate_db_tree_update_integrity_detects_segment_contamination(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("UPDATE features SET segment='9' WHERE accession='Q_OLD'")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "segment contamination detected" in result.stderr


def test_validate_db_tree_detects_mutation_drug_resistance_count_mismatch(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    _add_mutation_tables(basic_update_db)

    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute(
            "UPDATE sequence_drug_resistance SET mutations_detected=1, combination_status='partial' WHERE primary_accession='REF1' AND combination_id='comboA'"
        )
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "sequence_drug_resistance integrity" in result.stderr


def test_validate_db_tree_detects_mutation_orphan_accession(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    _add_mutation_tables(basic_update_db)

    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute(
            "INSERT INTO sequence_mutations VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("MISSING_ACC", "NS3:999X", "NS3", "1", 999, "X", ""),
        )
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "sequence_mutations integrity" in result.stderr


def test_validate_db_tree_update_integrity_ignores_placeholder_cluster_assignments_for_excluded_rows(tmp_path: Path):
    db_path = tmp_path / "placeholder_cluster.db"
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, cluster_98pct TEXT, segment TEXT, exclusion_status TEXT, exclusion_criteria TEXT)"
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, cluster_98pct, segment, exclusion_status, exclusion_criteria) VALUES (?, ?, ?, ?, ?)",
            [
                ("REF1", "REF1", "", "", ""),
                ("Q_PLACEABLE", "REF1", "", "", ""),
                ("Q_UNPLACED", "NA- see tree", "", "", ""),
            ],
        )
        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.executemany(
            "INSERT INTO sequences(header, sequence) VALUES (?, ?)",
            [("REF1", "ATGC"), ("Q_PLACEABLE", "ATGT"), ("Q_UNPLACED", "ATGA")],
        )
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, segment) VALUES (?, ?, ?, ?)",
            [("REF1", "REF1", "REF1", ""), ("Q_PLACEABLE", "Q_PLACEABLE", "REF1", "")],
        )
        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
        cur.execute(
            "INSERT INTO trees(name, source, newick, segment, segment_key) VALUES (?, ?, ?, ?, ?)",
            ("usher", "usher", "(REF1:0.1,Q_PLACEABLE:0.2);", "", None),
        )
        cur.execute("UPDATE meta_data SET exclusion_status='1', exclusion_criteria='alignment_filtering' WHERE primary_accession='Q_UNPLACED'")
        cur.execute("CREATE TABLE update_batches (batch_id TEXT, mode TEXT, started_at TEXT, finished_at TEXT)")
        cur.execute("INSERT INTO update_batches(batch_id, mode, started_at, finished_at) VALUES ('b1', 'update', '2026-01-01', '2026-01-01')")
        cur.execute("CREATE TABLE update_table_deltas (batch_id TEXT, table_name TEXT, before_count INTEGER, after_count INTEGER, delta INTEGER)")
        cur.execute("INSERT INTO update_table_deltas(batch_id, table_name, before_count, after_count, delta) VALUES ('b1', 'meta_data', 2, 3, 1)")
        cur.execute("CREATE TABLE features (accession TEXT, segment TEXT)")
        cur.executemany("INSERT INTO features(accession, segment) VALUES (?, ?)", [("REF1", ""), ("Q_PLACEABLE", "")])
        cur.execute("CREATE TABLE insertions (primary_accession TEXT)")
        cur.execute("CREATE TABLE host_taxa (primary_accession TEXT)")
        conn.commit()
    finally:
        conn.close()

    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(outdir),
            "--check-update-integrity",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "Missing centroids in tree: 0" in report
    assert "Non-centroid UShER leaves: 1" in report
    assert "Missing in tree (meta_data -> tree): 0" in report


def test_validate_db_tree_update_integrity_fails_when_latest_batch_is_noop(tmp_path: Path, basic_update_db: Path):
    _add_required_alignment_tables(basic_update_db)
    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DELETE FROM update_table_deltas WHERE batch_id='batch_seed'")
        conn.executemany(
            "INSERT INTO update_table_deltas(batch_id, table_name, before_count, after_count, delta) VALUES (?, ?, ?, ?, ?)",
            [
                ("batch_seed", "meta_data", 3, 3, 0),
                ("batch_seed", "sequences", 3, 3, 0),
                ("batch_seed", "sequence_alignment", 3, 3, 0),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(basic_update_db),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "latest update batch made no DB changes" in result.stderr


def test_validate_db_tree_segmented_test_mode_allows_subsampled_trees(tmp_path: Path):
    db_path = tmp_path / "segmented_subsample.db"
    _create_segmented_subsample_db(db_path)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "allowing per-segment UShER trees to be subsampled" in result.stdout


def test_validate_db_tree_segmented_test_mode_allows_subsampled_iqtree_trees(tmp_path: Path):
    db_path = tmp_path / "segmented_iqtree_subsample.db"
    _create_segmented_subsample_db(db_path, tree_source="iqtree")

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--segment-tree-source",
            "iqtree",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "allowing per-segment IQ-TREE trees to be subsampled" in result.stdout


def test_validate_db_tree_segmented_iqtree_centroid_only_comparison(tmp_path: Path):
    db_path = tmp_path / "segmented_iqtree_centroid.db"
    _create_segmented_iqtree_centroid_db(db_path)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
            "--segment-tree-source",
            "iqtree",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "IQ-TREE segmented validation: comparing per-segment centroid nodes only" in result.stdout


def test_validate_db_tree_segmented_non_test_mode_rejects_subsampled_trees(tmp_path: Path):
    db_path = tmp_path / "segmented_subsample_fail.db"
    _create_segmented_subsample_db(db_path)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(tmp_path / "out"),
            "--check-update-integrity",
            "--expect-segment-trees",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "per-segment UShER tree is missing unfiltered accessions" in result.stderr


def test_validate_db_tree_segmented_test_mode_skips_strict_consistency_without_update_checks(tmp_path: Path):
    db_path = tmp_path / "segmented_relaxed_consistency.db"
    _create_segmented_subsample_db(db_path)

    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute("DELETE FROM features WHERE accession='SEG2_B'")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(tmp_path / "out"),
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "skipping strict DB consistency failures for segmented validation run" in result.stdout


def test_validate_db_tree_ignores_reference_and_master_alignment_rows_in_consistency_check(tmp_path: Path):
    db_path = tmp_path / "alignment_refs_ok.db"
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, exclusion_status TEXT, exclusion_criteria TEXT)"
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, exclusion_status, exclusion_criteria) VALUES (?, ?, ?, ?)",
            [
                ("QUERY1", "query", "", ""),
                ("REF1", "reference", "1", "reference_context"),
                ("MASTER1", "master", "1", "master_context"),
            ],
        )
        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.execute("INSERT INTO sequences(header, sequence) VALUES ('QUERY1', 'ATGC')")
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, segment) VALUES (?, ?, ?, ?)",
            [
                ("QUERY1", "QUERY1", "REF1", "1"),
                ("REF1", "REF1", "REF1", "1"),
                ("MASTER1", "MASTER1", "MASTER1", "1"),
            ],
        )
        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
        cur.execute(
            "INSERT INTO trees(name, source, newick, segment, segment_key) VALUES ('usher_tree', 'usher', '(QUERY1:0.1);', '1', 'refset_1')"
        )
        cur.execute("CREATE TABLE features (accession TEXT)")
        cur.execute("INSERT INTO features(accession) VALUES ('QUERY1')")
        cur.execute("CREATE TABLE insertions (primary_accession TEXT)")
        cur.execute("INSERT INTO insertions(primary_accession) VALUES ('QUERY1')")
        cur.execute("CREATE TABLE host_taxa (other_col TEXT)")
        cur.execute("INSERT INTO host_taxa(other_col) VALUES ('placeholder')")
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(tmp_path / "out"),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "sequence_alignment vs meta_data: ok=True" in result.stdout


def test_validate_db_tree_segmented_allows_excluded_reference_context_rows(tmp_path: Path):
    db_path = tmp_path / "segmented_excluded_refs.db"
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT, exclusion_status TEXT, exclusion_criteria TEXT)"
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment, exclusion_status, exclusion_criteria) VALUES (?, ?, ?, ?, ?)",
            [
                ("SEG1_MASTER", "master", "1", "1", "reference_context"),
                ("SEG2_REF", "reference", "2", "1", "reference_context"),
                ("SEG1_Q", "query", "1", "", ""),
                ("SEG2_Q", "query", "2", "", ""),
            ],
        )
        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.executemany(
            "INSERT INTO sequences(header, sequence) VALUES (?, ?)",
            [("SEG1_Q", "ATGC"), ("SEG2_Q", "ATGA")],
        )
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, segment) VALUES (?, ?, ?, ?)",
            [
                ("SEG1_MASTER", "SEG1_MASTER", "SEG1_MASTER", "1"),
                ("SEG2_REF", "SEG2_REF", "SEG2_REF", "2"),
                ("SEG1_Q", "SEG1_Q", "SEG1_MASTER", "1"),
                ("SEG2_Q", "SEG2_Q", "SEG2_REF", "2"),
            ],
        )
        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
        cur.executemany(
            "INSERT INTO trees(name, source, newick, segment, segment_key) VALUES (?, ?, ?, ?, ?)",
            [
                ("usher_seg1", "usher", "(SEG1_Q:0.1);", "1", "refset_1"),
                ("usher_seg2", "usher", "(SEG2_Q:0.1);", "2", "refset_2"),
            ],
        )
        cur.execute("CREATE TABLE features (accession TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO features(accession, segment) VALUES (?, ?)",
            [("SEG1_Q", "1"), ("SEG2_Q", "2")],
        )
        cur.execute("CREATE TABLE insertions (primary_accession TEXT)")
        cur.executemany(
            "INSERT INTO insertions(primary_accession) VALUES (?)",
            [("SEG1_Q",), ("SEG2_Q",)],
        )
        cur.execute("CREATE TABLE host_taxa (primary_accession TEXT)")
        cur.executemany(
            "INSERT INTO host_taxa(primary_accession) VALUES (?)",
            [("SEG1_Q",), ("SEG2_Q",)],
        )
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(tmp_path / "out"),
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "Segmented validation: meta segments=2 tree segments=2" in result.stdout
    assert "sequence_alignment vs meta_data: ok=True" in result.stdout


def test_validate_db_tree_segmented_ignores_excluded_reference_alignment_rows(tmp_path: Path):
    db_path = tmp_path / "segmented_excluded_reference_alignment.db"
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT, exclusion_status TEXT, exclusion_criteria TEXT)"
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment, exclusion_status, exclusion_criteria) VALUES (?, ?, ?, ?, ?)",
            [
                ("SEG1_A", "query", "1", "", ""),
                ("SEG2_A", "query", "2", "", ""),
                ("REF_EXCL", "reference", "1", "1", "trimmed_from_update_assets"),
            ],
        )
        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.executemany(
            "INSERT INTO sequences(header, sequence) VALUES (?, ?)",
            [("SEG1_A", "ATGC"), ("SEG2_A", "ATGA")],
        )
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, segment) VALUES (?, ?, ?, ?)",
            [
                ("SEG1_A", "SEG1_A", "SEG1_A", "1"),
                ("SEG2_A", "SEG2_A", "SEG2_A", "2"),
                ("REF_EXCL", "REF_EXCL", "REF_EXCL", "1"),
            ],
        )
        cur.execute("CREATE TABLE features (accession TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO features(accession, segment) VALUES (?, ?)",
            [("SEG1_A", "1"), ("SEG2_A", "2")],
        )
        cur.execute("CREATE TABLE insertions (primary_accession TEXT)")
        cur.executemany(
            "INSERT INTO insertions(primary_accession) VALUES (?)",
            [("SEG1_A",), ("SEG2_A",)],
        )
        cur.execute("CREATE TABLE host_taxa (primary_accession TEXT)")
        cur.executemany(
            "INSERT INTO host_taxa(primary_accession) VALUES (?)",
            [("SEG1_A",), ("SEG2_A",)],
        )
        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
        cur.executemany(
            "INSERT INTO trees(name, source, newick, segment, segment_key) VALUES (?, ?, ?, ?, ?)",
            [
                ("usher_seg1", "usher", "(SEG1_A:0.1);", "1", "refset_1"),
                ("usher_seg2", "usher", "(SEG2_A:0.1);", "2", "refset_2"),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(tmp_path / "out"),
            "--expect-segment-trees",
            "--test-mode",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "Segmented validation: meta segments=2 tree segments=2" in result.stdout
    assert "sequence_alignment vs meta_data: ok=True" in result.stdout


def test_validate_db_tree_prefers_largest_usher_tree_row(tmp_path: Path):
    db_path = tmp_path / "multiple_usher_rows.db"
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, cluster_95pct TEXT, exclusion_status TEXT, exclusion_criteria TEXT)"
        )
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, cluster_95pct, exclusion_status, exclusion_criteria) VALUES (?, ?, ?, ?)",
            [("A", "A", "", ""), ("B", "B", "", ""), ("C", "C", "", "")],
        )
        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        cur.executemany(
            "INSERT INTO sequences(header, sequence) VALUES (?, ?)",
            [("A", "ATGC"), ("B", "ATGT"), ("C", "ATGA")],
        )
        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, segment) VALUES (?, ?, ?, ?)",
            [("A", "A", "A", ""), ("B", "B", "B", ""), ("C", "C", "C", "")],
        )
        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
        cur.executemany(
            "INSERT INTO trees(name, source, newick, segment, segment_key) VALUES (?, ?, ?, ?, ?)",
            [
                ("usher_small", "usher", "(A:0.1);", "", None),
                ("usher_large", "usher", "(A:0.1,B:0.2,C:0.3);", "", None),
                ("iqtree_large", "iqtree", "(A:0.1,B:0.2,C:0.3,D:0.4);", "", None),
            ],
        )
        cur.execute("CREATE TABLE insertions (primary_accession TEXT)")
        cur.executemany("INSERT INTO insertions(primary_accession) VALUES (?)", [("A",), ("B",), ("C",)])
        cur.execute("CREATE TABLE host_taxa (primary_accession TEXT)")
        cur.executemany("INSERT INTO host_taxa(primary_accession) VALUES (?)", [("A",), ("B",), ("C",)])
        cur.execute("CREATE TABLE features (accession TEXT)")
        cur.executemany("INSERT INTO features(accession) VALUES (?)", [("A",), ("B",), ("C",)])
        conn.commit()
    finally:
        conn.close()

    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--db",
            str(db_path),
            "--outdir",
            str(outdir),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "Tree terminals: 3" in report
    assert "Tree source: usher" in report
    assert "Tree name: usher_large" in report
