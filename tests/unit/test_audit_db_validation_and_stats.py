"""Carpet-bomb audit tests for SQLite DB layer, validators, mutations, and stats math.

Includes direct coverage of TestPipelineOutput, SoftwareVersion, growth_core,
growth_stats, and tree_growth_metrics alongside strict xfail bug reproducers
(AUDIT-DB-*, AUDIT-VAL-*, AUDIT-MUT-*, AUDIT-STAT-*).
"""

import argparse
import os
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import CreateSqliteDB as cdb_mod
import SoftwareVersion as sv_mod
import TestPipelineOutput as tpo_mod
import ValidateDbTree as vdt_mod
import VerifyMutations as vm_mod
import growth_core as gc_mod
import growth_stats as gs_mod
import protein_alleles as pal_mod
import tree_growth_metrics as tgm_mod
import within_lineage_growth as wlg_mod


def _make_minimal_db(db_path: Path, *, segmented: bool = False, with_trees: bool = True):
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute(
        "CREATE TABLE meta_data (primary_accession TEXT, segment TEXT, accession_type TEXT, "
        "exclusion_status TEXT, exclusion_criteria TEXT, host_taxa_id TEXT, collection_date TEXT, collection_year TEXT)"
    )
    cur.execute(
        "CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT, segment TEXT)"
    )
    cur.execute(
        "CREATE TABLE sequences (header TEXT, sequence TEXT, segment TEXT)"
    )
    cur.execute(
        "CREATE TABLE features (accession TEXT, cds_start_OG_seq TEXT, cds_end_OG_seq TEXT, "
        "aln_start TEXT, aln_end TEXT, product TEXT, segment TEXT)"
    )
    cur.execute("CREATE TABLE insertions (primary_accession TEXT, reference TEXT, insertion TEXT, segment TEXT)")
    cur.execute("CREATE TABLE genes (description TEXT, display_name TEXT, name TEXT, parent_name TEXT)")
    cur.execute("CREATE TABLE host_taxa (taxa_id TEXT, name TEXT, name_type TEXT, taxonomy_level TEXT, common_name TEXT)")
    cur.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, newick TEXT)")
    cur.execute("CREATE TABLE info (creation_type TEXT, creation_date TEXT)")
    cur.execute("INSERT INTO info VALUES ('new db', '2026-11-24')")

    segs = ["1", "2"] if segmented else ["1"]
    for seg in segs:
        m = f"M{seg}"
        q = f"Q{seg}"
        cur.execute("INSERT INTO meta_data VALUES (?, ?, 'master', '0', '', '9606', '2024-01-01', '2024')", (m, seg))
        cur.execute("INSERT INTO meta_data VALUES (?, ?, 'query', '0', '', '9606', '2024-02-01', '2024')", (q, seg))
        cur.execute("INSERT INTO sequence_alignment VALUES (?, ?, ?, 'ATGCGT', ?)", (m, m, m, seg))
        cur.execute("INSERT INTO sequence_alignment VALUES (?, ?, ?, 'ATGCGT', ?)", (q, q, m, seg))
        cur.execute("INSERT INTO sequences VALUES (?, 'ATGCGT', ?)", (m, seg))
        cur.execute("INSERT INTO sequences VALUES (?, 'ATGCGT', ?)", (q, seg))
        cur.execute("INSERT INTO features VALUES (?, '1', '6', '1', '6', 'P', ?)", (m, seg))
        cur.execute("INSERT INTO features VALUES (?, '1', '6', '1', '6', 'P', ?)", (q, seg))
        if with_trees:
            cur.execute("INSERT INTO trees VALUES (?, 'usher', ?, ?)", (f"usher_{seg}", seg, f"({m}:0.1,{q}:0.1);"))
    cur.execute("INSERT INTO host_taxa VALUES ('9606', 'Homo sapiens', 'scientific', 'species', 'human')")
    cur.execute("INSERT INTO genes VALUES ('P', 'P', 'P', 'whole_genome')")
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# Passing tests raising coverage on TestPipelineOutput, SoftwareVersion, and stats
# ---------------------------------------------------------------------------

def test_test_pipeline_output_non_segmented_passes_on_valid_inputs(tmp_path, monkeypatch):
    blast = tmp_path / "query_uniq_tophits.tsv"
    blast.write_text("Q1\tM1\t99.0\t+\n", encoding="utf-8")
    mat = tmp_path / "gB_matrix_validated.tsv"
    mat.write_text("primary_accession\texclusion_status\texclusion_criteria\nM1\t0\t\nQ1\t0\t\n", encoding="utf-8")
    db = tmp_path / "valid.db"
    _make_minimal_db(db, segmented=False)
    out = tmp_path / "nonseg_out.txt"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "TestPipelineOutput.py",
            "--mode", "non_segmented",
            "--blast_hits", str(blast),
            "--gb_matrix", str(mat),
            "--sqlite_db", str(db),
            "--output", str(out),
        ],
    )
    rc = tpo_mod.main()
    assert rc == 0
    assert out.exists()
    assert "PASS" in out.read_text(encoding="utf-8")


def test_test_pipeline_output_segmented_passes_on_valid_inputs(tmp_path, monkeypatch):
    ablast = tmp_path / "query_uniq_tophit_annotated.tsv"
    ablast.write_text("Q1\tM1\t99.0\t+\t1\nQ2\tM2\t99.0\t+\t2\n", encoding="utf-8")
    vmat = tmp_path / "gB_matrix_validated_segment.tsv"
    vmat.write_text(
        "primary_accession\tsegment\tsegment_validated\texclusion_status\n"
        "M1\t1\t1\t0\nQ1\t1\t1\t0\nM2\t2\t2\t0\nQ2\t2\t2\t0\n",
        encoding="utf-8",
    )
    pmat = tmp_path / "gB_matrix_pivoted_segments.tsv"
    pmat.write_text(
        "Parsed_strain\t1\t2\tComplete_status\n"
        "strainA\tQ1\tQ2\tComplete\n",
        encoding="utf-8",
    )
    db = tmp_path / "valid_seg.db"
    _make_minimal_db(db, segmented=True)
    out = tmp_path / "seg_out.txt"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "TestPipelineOutput.py",
            "--mode", "segmented",
            "--annotated_blast", str(ablast),
            "--validated_matrix", str(vmat),
            "--pivoted_matrix", str(pmat),
            "--sqlite_db", str(db),
            "--output", str(out),
        ],
    )
    rc = tpo_mod.main()
    assert rc == 0
    assert out.exists()
    assert "PASS" in out.read_text(encoding="utf-8")


def test_test_pipeline_output_fails_on_empty_blast_with_active_queries(tmp_path, monkeypatch):
    blast = tmp_path / "query_uniq_tophits.tsv"
    blast.write_text("", encoding="utf-8")
    mat = tmp_path / "gB_matrix.tsv"
    mat.write_text("primary_accession\texclusion_status\nQ1\t0\n", encoding="utf-8")
    db = tmp_path / "valid.db"
    _make_minimal_db(db, segmented=False)
    out = tmp_path / "out.txt"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "TestPipelineOutput.py",
            "--mode", "non_segmented",
            "--blast_hits", str(blast),
            "--gb_matrix", str(mat),
            "--sqlite_db", str(db),
            "--output", str(out),
        ],
    )
    rc = tpo_mod.main()
    assert rc == 1
    assert "FAIL" in out.read_text(encoding="utf-8")


def test_software_version_checker_writes_tsv(tmp_path, monkeypatch):
    monkeypatch.setattr(
        sv_mod.subprocess,
        "run",
        lambda *a, **kw: argparse.Namespace(returncode=0, stdout="tool 1.2.3\n", stderr=""),
    )
    checker = sv_mod.SoftwareVersionChecker(str(tmp_path), "Software_info", "software_info.tsv")
    checker.save_versions_to_tsv()
    out_file = tmp_path / "Software_info" / "software_info.tsv"
    assert out_file.exists()
    txt = out_file.read_text(encoding="utf-8")
    assert "Python" in txt and "vgtk" in txt


def test_growth_core_and_stats_basic_behaviours():
    dt, prec = gc_mod.parse_collection_date("2024-03-15")
    assert prec == "day" and dt.year == 2024 and dt.month == 3
    p_poisson = gs_mod.poisson_growth(
        np.array([0.0, 1.0, 2.0, 3.0]),
        np.array([10.0, 12.0, 14.0, 16.0]),
        np.array([100.0, 100.0, 100.0, 100.0]),
    )
    assert p_poisson["converged"] is True
    assert gs_mod.chi2_sf(0.0, 1) == 1.0
    tree = tgm_mod.parse_newick("((A:1,B:1):1,(C:1,D:1):1);")
    assert len(tree.tip_counts()) > 0


# ---------------------------------------------------------------------------
# Promoted regression tests for fixed DB & Validator bugs (WP-02)
# ---------------------------------------------------------------------------

def test_create_sqlite_db_strain_retirement_cleans_current_mutation_tables(tmp_path):
    src_db = tmp_path / "src.db"
    _make_minimal_db(src_db, segmented=False)
    conn = sqlite3.connect(str(src_db))
    for tbl in (
        "sequence_relevant_mutation_summary",
        "completed_signatures_only",
        "sequence_mutation_calls",
        "sequence_reading_frame_scan",
    ):
        conn.execute(f"CREATE TABLE {tbl} (primary_accession TEXT, val TEXT)")
        conn.execute(f"INSERT INTO {tbl} VALUES ('OLD_RETIRED', 'hit')")
    conn.execute("INSERT INTO meta_data VALUES ('OLD_RETIRED', '1', 'query', '0', '', '9606', '2024-01-01', '2024')")
    conn.commit()
    conn.close()

    out_db = tmp_path / "out.db"
    import shutil
    shutil.copy2(src_db, out_db)
    helper = cdb_mod.CreateSqliteDB.__new__(cdb_mod.CreateSqliteDB)
    helper.update = True
    helper.update_db = str(src_db)
    helper.is_segmented = "N"
    helper.batch_id = "test_batch"

    conn = sqlite3.connect(str(out_db))
    helper._apply_update_db_strain_retirements(
        conn,
        superseded_update_list=[
            {
                "retired_accession": "OLD_RETIRED",
                "segment": "1",
                "winner_accession": "Q1",
                "reason": "duplicate",
            }
        ],
        update_winner_patches={},
        update_exclusions=[],
        now_str="2026-11-24 00:00:00",
    )
    conn.commit()
    orphans = {
        tbl: conn.execute(f"SELECT COUNT(*) FROM {tbl} WHERE primary_accession='OLD_RETIRED'").fetchone()[0]
        for tbl in (
            "sequence_relevant_mutation_summary",
            "completed_signatures_only",
            "sequence_mutation_calls",
            "sequence_reading_frame_scan",
        )
    }
    conn.close()
    assert all(v == 0 for v in orphans.values()), f"Orphan rows left behind: {orphans}"


def test_create_sqlite_db_preserves_influenza_na_in_features_and_insertions(tmp_path):
    feat = tmp_path / "features.tsv"
    feat.write_text(
        "accession\tcds_start_OG_seq\tcds_end_OG_seq\taln_start\taln_end\tproduct\tsegment\n"
        "CY001\t1\t1410\t1\t1410\tNA\t6\n",
        encoding="utf-8",
    )
    src = Path(cdb_mod.__file__).read_text(encoding="utf-8")
    assert 'self._read_tsv_required(self.pad_aln, [], "pad_aln", dtype=str, keep_na_strings=True)' in src
    assert 'self._read_tsv_required(self.features, [], "features", dtype=str, keep_na_strings=True)' in src
    assert 'self._read_tsv_required(self.insertions, [], "insertions", dtype=str, keep_na_strings=True)' in src

    helper = cdb_mod.CreateSqliteDB.__new__(cdb_mod.CreateSqliteDB)
    df = helper._read_tsv_required(str(feat), ["accession", "product"], "features", keep_na_strings=True)
    assert df.loc[0, "product"] == "NA"


def test_create_sqlite_db_load_fasta_handles_empty_fasta(tmp_path):
    empty_fa = tmp_path / "empty.fa"
    empty_fa.write_text("", encoding="utf-8")
    helper = cdb_mod.CreateSqliteDB.__new__(cdb_mod.CreateSqliteDB)
    helper.fasta_sequence_file = str(empty_fa)
    df = helper.load_fasta()
    assert "header" in df.columns and "sequence" in df.columns


def test_create_sqlite_db_atomic_staging_preserves_original_on_crash(tmp_path, monkeypatch):
    """AUDIT-DB-04: Mid-build crash during in-place update leaves target DB intact with no .tmp.* left."""
    out_dir = tmp_path / "SqliteDB"
    out_dir.mkdir(parents=True)
    target_db = out_dir / "upd.db"
    _make_minimal_db(target_db, segmented=False)

    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "aln.tsv"
    m49_c = tmp_path / "m49_c.csv"
    m49_i = tmp_path / "m49_i.csv"
    m49_r = tmp_path / "m49_r.csv"
    m49_s = tmp_path / "m49_s.csv"
    proj = tmp_path / "proj.tsv"
    ins = tmp_path / "ins.tsv"
    host = tmp_path / "host.tsv"
    fa = tmp_path / "seqs.fa"
    meta.write_text("primary_accession\texclusion\tsegment\nQ_NEW\t0\t1\n", encoding="utf-8")
    features.write_text("accession\tmaster_ref_accession\treference_accession\taln_start\taln_end\tcds_start\tcds_end\tcds_start_OG_seq\tcds_end_OG_seq\tproduct\tsegment\nQ_NEW\tM1\tM1\t1\t6\t1\t6\t1\t6\tP\t1\n", encoding="utf-8")
    aln.write_text("primary_accession\talignment_name\talignment\tsegment\nQ_NEW\tM1\tATGCGT\t1\n", encoding="utf-8")
    m49_c.write_text("m49_code,name\n001,World\n", encoding="utf-8")
    m49_i.write_text("code,name\nX,Inter\n", encoding="utf-8")
    m49_r.write_text("code,name\nY,Region\n", encoding="utf-8")
    m49_s.write_text("code,name\nZ,Sub\n", encoding="utf-8")
    proj.write_text("Software\tVersion\nPython\t3.10\n", encoding="utf-8")
    ins.write_text("primary_accession\treference\tinsertion\tsegment\n", encoding="utf-8")
    host.write_text("taxa_id\tname\n9606\tHomo sapiens\n", encoding="utf-8")
    fa.write_text(">Q_NEW\nATGCGT\n", encoding="utf-8")

    creator = cdb_mod.CreateSqliteDB(
        meta_data=str(meta), features=str(features), pad_aln=str(aln), gene_info=None,
        m49_countries=str(m49_c), m49_interm_region=str(m49_i), m49_regions=str(m49_r),
        m49_sub_regions=str(m49_s), proj_settings=str(proj), fasta_sequence_file=str(fa),
        insertions=str(ins), host_taxa_file=str(host),
        base_dir=str(tmp_path), output_dir="SqliteDB", db_name="upd",
        db_status="last updated", update=True, update_db=str(target_db), batch_id="crash_test",
    )

    def _boom(*_args, **_kwargs):
        raise RuntimeError("simulated mid-update crash")

    monkeypatch.setattr(creator, "merge_table_append_nonredundant", _boom)
    with pytest.raises(RuntimeError, match="simulated mid-update crash"):
        creator.create_db()

    # Original DB untouched & uncorrupted, no leftover .tmp.* files
    assert list(out_dir.glob("*.tmp.*")) == []
    conn = sqlite3.connect(str(target_db))
    rows = {r[0] for r in conn.execute("SELECT primary_accession FROM meta_data").fetchall()}
    assert "Q_NEW" not in rows and "M1" in rows
    if cdb_mod.CreateSqliteDB._table_exists(conn, "update_batches"):
        unfinished = conn.execute("SELECT COUNT(*) FROM update_batches WHERE finished_at IS NULL OR TRIM(finished_at)=''").fetchone()[0]
        assert unfinished == 0
    conn.close()


def test_validate_db_tree_writes_fail_status_when_validation_fails(tmp_path):
    db = tmp_path / "broken.db"
    _make_minimal_db(db, segmented=False)
    conn = sqlite3.connect(str(db))
    # Delete sequence_alignment row for Q1 so validation MUST fail
    conn.execute("DELETE FROM sequence_alignment WHERE primary_accession='Q1'")
    conn.commit()
    conn.close()

    outdir = tmp_path / "val_out"
    outdir.mkdir()
    with pytest.raises(SystemExit):
        vdt_mod.main(["--db", str(db), "--outdir", str(outdir)])
    report = (outdir / "db_tree_validation.txt").read_text(encoding="utf-8")
    assert "Validation status: PASS" not in report
    assert "Validation status: FAIL" in report


def test_validate_feature_projection_integrity_catches_out_of_window_cds(tmp_path):
    db = tmp_path / "clamped.db"
    conn = sqlite3.connect(str(db))
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, exclusion_status TEXT, segment TEXT)")
    conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT, segment TEXT)")
    conn.execute(
        "CREATE TABLE features (accession TEXT, master_ref_accession TEXT, product TEXT, "
        "aln_start TEXT, aln_end TEXT, cds_start TEXT, cds_end TEXT, segment TEXT)"
    )
    conn.execute("INSERT INTO meta_data VALUES ('M1', 'master', '0', '1')")
    conn.execute("INSERT INTO meta_data VALUES ('Q1', 'query', '0', '1')")
    conn.execute("INSERT INTO sequence_alignment VALUES ('M1', 'M1', 'ATGAAACCCGGG', '1')")
    conn.execute("INSERT INTO sequence_alignment VALUES ('Q1', 'M1', '---AAACCC---', '1')")
    conn.execute("INSERT INTO features VALUES ('M1', 'M1', 'P', '1', '12', '1', '12', '1')")
    # Deliberately store out-of-window cds_start/cds_end (1..12 while aln is 4..9):
    conn.execute("INSERT INTO features VALUES ('Q1', 'M1', 'P', '4', '9', '1', '12', '1')")
    conn.commit()
    res = vdt_mod.validate_feature_projection_integrity(conn)
    conn.close()
    assert vdt_mod.result_failed(res)


def test_validate_sequence_alignment_vs_meta_handles_null_exclusion_status(tmp_path):
    db = tmp_path / "null_excl.db"
    conn = sqlite3.connect(str(db))
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, exclusion_status TEXT)")
    conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    conn.execute("INSERT INTO meta_data VALUES ('M1', 'master', NULL)")
    conn.execute("INSERT INTO sequence_alignment VALUES ('M1', 'M1', 'ATGC')")
    conn.commit()
    res = vdt_mod.validate_sequence_alignment_vs_meta(conn, {"M1"}, "primary_accession", "exclusion_status", "1")
    conn.close()
    assert not vdt_mod.result_failed(res)


def test_validate_collection_dates_records_finding_on_invalid_calendar_date(tmp_path):
    db = tmp_path / "bad_date.db"
    conn = sqlite3.connect(str(db))
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, exclusion_status TEXT, collection_date TEXT, collection_year TEXT)")
    conn.execute("INSERT INTO meta_data VALUES ('Q1', '0', '2020-02-30', '2020')")
    conn.commit()
    res = vdt_mod.validate_collection_dates(conn, ["primary_accession", "exclusion_status", "collection_date", "collection_year"])
    conn.close()
    assert vdt_mod.result_failed(res)


def test_verify_mutations_regex_and_translator_accept_stop_and_deletion():
    vm_src = Path(vm_mod.__file__).read_text(encoding="utf-8")
    assert r"([A-Z_X])" not in vm_src
    assert vm_mod.translate_codon("TAA") == "*"
    assert vm_mod.translate_codon("---") == "-"


def test_protein_alleles_reads_sequence_mutation_calls(tmp_path):
    db = tmp_path / "mut.db"
    conn = sqlite3.connect(str(db))
    conn.execute(
        "CREATE TABLE sequence_mutation_calls "
        "(primary_accession TEXT, mutation_id TEXT, protein_name TEXT, aa_position INTEGER, alt_residue TEXT, call_status TEXT)"
    )
    conn.execute("INSERT INTO sequence_mutation_calls VALUES ('Q1', 'NS3:80K', 'NS3', 80, 'K', 'emitted')")
    conn.commit()
    prots = pal_mod.list_proteins(conn)
    alleles = pal_mod.catalog_alleles(conn, "NS3")
    conn.close()
    assert not prots.empty and len(alleles) == 1


def test_firth_likelihood_ratio_p_is_time_unit_invariant_on_flat_data():
    t_years = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    t_days = t_years * 365.25
    trials = np.full(6, 20.0)
    y = np.full(6, 0.5)
    X_yr = np.column_stack([np.ones(6), t_years])
    X_day = np.column_stack([np.ones(6), t_days])
    p_yr, _ = gs_mod.likelihood_ratio_p(X_yr, y, column=1, family="binomial", trials=trials, firth=True)
    p_day, _ = gs_mod.likelihood_ratio_p(X_day, y, column=1, family="binomial", trials=trials, firth=True)
    assert abs(p_yr - p_day) < 1e-6
    assert p_yr > 0.99


def test_growth_core_chi2_sf_matches_growth_stats_precision():
    # Exact chi2_sf(3.841458820694124, df=1) == 0.0500000
    x = 3.841458820694124
    got = gc_mod.chi2_sf(x, 1)
    assert abs(got - 0.05) < 1e-6


def test_skyline_accumulates_waiting_time_across_sampling_intervals():
    # Three intervals leading into one coalescence at tau=3.0:
    # (0..1, k=2 => w=1), (1..2, k=3 => w=3), (2..3, k=4 => w=6) => total Ne = 1 + 3 + 6 = 10.0
    intervals = [
        (0.0, 1.0, 2),
        (1.0, 2.0, 3),
        (2.0, 3.0, 4),
    ]
    pts = tgm_mod.skyline(intervals, [3.0])
    assert len(pts) == 1
    assert abs(float(pts[0]["ne"]) - 10.0) < 1e-9


def test_within_lineage_growth_excludes_qc_failed_rows(tmp_path):
    db = tmp_path / "wlg.db"
    conn = sqlite3.connect(str(db))
    conn.execute(
        "CREATE TABLE meta_data (primary_accession TEXT, collection_date TEXT, collection_year TEXT, "
        "country TEXT, host TEXT, exclusion_status TEXT, nearest_reference_genotype TEXT, nearest_reference_subtype TEXT)"
    )
    conn.execute("INSERT INTO meta_data VALUES ('Q_KEEP', '2024-01-01', '2024', 'UK', 'Human', '0', '1', '1a')")
    conn.execute("INSERT INTO meta_data VALUES ('Q_DROP', '2024-01-02', '2024', 'UK', 'Human', '1', '1', '1a')")
    conn.commit()
    conn.close()

    with wlg_mod.WithinLineageGrowthAnalyzer(db_path=str(db), protein="NS3") as analyzer:
        cohort = analyzer.load_cohort()
    assert set(cohort["primary_accession"]) == {"Q_KEEP"}
