import sqlite3
import shutil
from pathlib import Path

import pandas as pd
import pytest

from CreateSqliteDB import CreateSqliteDB


REPO_ROOT = Path(__file__).resolve().parents[2]
REAL_UPDATE_DB = REPO_ROOT / "test_data" / "RABV_test" / "rabv-jul0425.db"


@pytest.fixture()
def real_update_db_copy(tmp_path: Path) -> Path:
    if not REAL_UPDATE_DB.exists():
        pytest.skip(f"Real update DB not found at {REAL_UPDATE_DB}")
    dst = tmp_path / "rabv-jul0425.copy.db"
    shutil.copyfile(REAL_UPDATE_DB, dst)
    return dst


def _write_tsv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _write_csv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def _inputs(tmp_path: Path, suffix: str, aln_a: str = "ATGC"):
    meta = tmp_path / f"meta_{suffix}.tsv"
    features = tmp_path / f"features_{suffix}.tsv"
    aln = tmp_path / f"sequence_alignment_{suffix}.tsv"
    gene = tmp_path / f"gene_{suffix}.tsv"
    m49_country = tmp_path / f"m49_country_{suffix}.csv"
    m49_inter = tmp_path / f"m49_inter_{suffix}.csv"
    m49_region = tmp_path / f"m49_region_{suffix}.csv"
    m49_sub = tmp_path / f"m49_sub_{suffix}.csv"
    proj = tmp_path / f"software_{suffix}.tsv"
    insertions = tmp_path / f"insertions_{suffix}.tsv"
    host_taxa = tmp_path / f"host_{suffix}.tsv"
    fasta = tmp_path / f"seqs_{suffix}.fa"

    _write_tsv(meta, [["A", "", "1"]], ["primary_accession", "exclusion", "segment"])
    _write_tsv(features, [["A", "M", "R", "1", "10", "1", "10", "P", "1"]], ["accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end", "cds_start", "cds_end", "product", "segment"])
    _write_tsv(aln, [["A", "R", aln_a, "1"]], ["primary_accession", "alignment_name", "alignment", "segment"])
    _write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    _write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    _write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    _write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    _write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    _write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    _write_tsv(insertions, [["A", "R", "ins:5:A", "1"]], ["primary_accession", "reference", "insertion", "segment"])
    _write_tsv(host_taxa, [["A", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">A\nATGC\n", encoding="utf-8")

    return {
        "meta": meta,
        "features": features,
        "aln": aln,
        "gene": gene,
        "m49_country": m49_country,
        "m49_inter": m49_inter,
        "m49_region": m49_region,
        "m49_sub": m49_sub,
        "proj": proj,
        "insertions": insertions,
        "host_taxa": host_taxa,
        "fasta": fasta,
    }


def _build_db(tmp_path: Path, inp: dict, update=False, update_db=None, filtered_ids_file=None, iqtree_file=None, usher_tree=None):
    db = CreateSqliteDB(
        meta_data=str(inp["meta"]),
        features=str(inp["features"]),
        pad_aln=str(inp["aln"]),
        gene_info=str(inp["gene"]),
        m49_countries=str(inp["m49_country"]),
        m49_interm_region=str(inp["m49_inter"]),
        m49_regions=str(inp["m49_region"]),
        m49_sub_regions=str(inp["m49_sub"]),
        proj_settings=str(inp["proj"]),
        fasta_sequence_file=str(inp["fasta"]),
        insertions=str(inp["insertions"]),
        host_taxa_file=str(inp["host_taxa"]),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="upd_db",
        db_status="last updated" if update else "new db",
        iqtree_file=str(iqtree_file) if iqtree_file else None,
        usher_tree=str(usher_tree) if usher_tree else None,
        update=update,
        update_db=str(update_db) if update_db else None,
        batch_id="batch_test",
        filtered_ids_file=str(filtered_ids_file) if filtered_ids_file else None,
    )
    db.create_db()
    return tmp_path / "SqliteDB" / "upd_db.db"


def _write_realdb_compatible_lookup_tables(inp: dict):
    pd.DataFrame(
        [["id1", "001", "World", "World Full", "0", "0", "0", "dev", "r1", "sr1", "ir1"]],
        columns=[
            "id",
            "m49_code",
            "display_name",
            "full_name",
            "is_ldc",
            "is_lldc",
            "is_sids",
            "development_status",
            "m49_region_id",
            "m49_sub_region_id",
            "m49_intermediate_region_id",
        ],
    ).to_csv(inp["m49_country"], index=False)

    pd.DataFrame(
        [["ir1", "100", "Inter", "sr1"]],
        columns=["id", "m49_code", "display_name", "m49_sub_region_id"],
    ).to_csv(inp["m49_inter"], index=False)

    pd.DataFrame(
        [["r1", "10", "Region"]],
        columns=["id", "m49_code", "display_name"],
    ).to_csv(inp["m49_region"], index=False)

    pd.DataFrame(
        [["sr1", "11", "SubRegion", "r1"]],
        columns=["id", "m49_code", "display_name", "m49_region_id"],
    ).to_csv(inp["m49_sub"], index=False)


def _write_realdb_compatible_meta(inp: dict, db_path: Path, rows: list[dict]):
    conn = sqlite3.connect(str(db_path))
    try:
        base = pd.read_sql_query("SELECT * FROM meta_data LIMIT 1", conn)
    finally:
        conn.close()

    if base.empty:
        pytest.skip("Real update DB has empty meta_data table")

    base_row = base.iloc[0].to_dict()
    out_rows = []
    for override in rows:
        row = dict(base_row)
        row.update(override)
        out_rows.append(row)

    pd.DataFrame(out_rows, columns=base.columns).to_csv(inp["meta"], sep="\t", index=False)


def _write_realdb_compatible_genes(inp: dict, db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        base = pd.read_sql_query("SELECT * FROM genes LIMIT 1", conn)
    finally:
        conn.close()

    if base.empty:
        # fallback to current minimal schema if real DB has no rows
        pd.DataFrame([["geneA", "Gene A", "Gene A", "Genome"]], columns=["name", "description", "display_name", "parent_name"]).to_csv(
            inp["gene"], sep="\t", index=False
        )
        return

    row = base.iloc[0].to_dict()
    row.update({"name": "geneA", "description": "Gene A", "display_name": "Gene A", "parent_name": "Genome"})
    pd.DataFrame([row], columns=base.columns).to_csv(inp["gene"], sep="\t", index=False)


def test_create_sqlite_db_update_mode_upserts_without_growth(tmp_path: Path):
    initial = _inputs(tmp_path, "initial", aln_a="ATGC")
    _build_db(tmp_path, initial, update=False)

    db_path = tmp_path / "SqliteDB" / "upd_db.db"
    assert db_path.exists()

    update_inputs = _inputs(tmp_path, "update", aln_a="AT--")
    _build_db(tmp_path, update_inputs, update=True, update_db=db_path)

    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute("SELECT alignment FROM sequence_alignment WHERE primary_accession='A' AND segment='1'")
        assert cur.fetchone()[0] == "AT--"
        cur.execute("SELECT COUNT(*) FROM sequence_alignment WHERE primary_accession='A' AND segment='1'")
        assert cur.fetchone()[0] == 1

        # Re-run same update; row count should remain stable (idempotent net effect)
        _build_db(tmp_path, update_inputs, update=True, update_db=db_path)
        cur.execute("SELECT COUNT(*) FROM sequence_alignment WHERE primary_accession='A' AND segment='1'")
        assert cur.fetchone()[0] == 1
    finally:
        conn.close()


def test_create_sqlite_db_filtered_ids_do_not_exclude_reference_rows(tmp_path: Path):
    inp = _inputs(tmp_path, "filtered_refs", aln_a="ATGC")

    # Override meta_data to include explicit master/query typing
    pd.DataFrame(
        [
            ["REF1", "", "1", "master"],
            ["Q1", "", "1", "query"],
        ],
        columns=["primary_accession", "exclusion", "segment", "accession_type"],
    ).to_csv(inp["meta"], sep="\t", index=False)

    # Keep alignment/sequence rows for both accessions to make filtering behavior explicit
    pd.DataFrame(
        [["REF1", "REF1", "ATGC", "1"], ["Q1", "REF1", "AT--", "1"]],
        columns=["primary_accession", "alignment_name", "alignment", "segment"],
    ).to_csv(inp["aln"], sep="\t", index=False)
    inp["fasta"].write_text(">REF1\nATGC\n>Q1\nATTT\n", encoding="utf-8")

    filtered_ids = tmp_path / "filtered_ids.txt"
    filtered_ids.write_text("REF1\nQ1\n", encoding="utf-8")

    _build_db(tmp_path, inp, update=False, filtered_ids_file=filtered_ids)

    conn = sqlite3.connect(str(tmp_path / "SqliteDB" / "upd_db.db"))
    try:
        cur = conn.cursor()
        cur.execute("SELECT primary_accession FROM meta_data ORDER BY primary_accession")
        kept = [row[0] for row in cur.fetchall()]
        assert kept == ["REF1"]

        cur.execute("SELECT primary_accession, reason FROM excluded_accessions ORDER BY primary_accession")
        excluded = cur.fetchall()
        assert excluded == [("Q1", "alignment_filtering")]
    finally:
        conn.close()


def test_update_mode_logs_duplicate_non_upsert_keys_against_real_db(tmp_path: Path, real_update_db_copy: Path):
    inp = _inputs(tmp_path, "real_dup_keys", aln_a="ATGC")

    conn = sqlite3.connect(str(real_update_db_copy))
    try:
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM update_exclusions")
        seed_update_exclusions = cur.fetchone()[0]
        cur.execute("SELECT m49_code FROM m49_country WHERE m49_code IS NOT NULL AND TRIM(m49_code) != '' LIMIT 1")
        row = cur.fetchone()
    finally:
        conn.close()

    if not row:
        pytest.skip("No m49_code values found in real update DB")
    existing_code = str(row[0]).strip()

    _write_realdb_compatible_meta(
        inp,
        real_update_db_copy,
        [
            {"primary_accession": "A", "segment": "1", "accession_type": "query", "exclusion_status": ""},
            {"primary_accession": "B", "segment": "1", "accession_type": "query", "exclusion_status": ""},
        ],
    )

    # Intentionally collide with existing m49_country code and add one synthetic new code.
    pd.DataFrame(
        [
            ["id_dup", existing_code, "World_duplicate", "World_duplicate_full", "0", "0", "0", "dev", "r1", "sr1", "ir1"],
            ["id_new", "ZZZ", "Synthetic Test Region", "Synthetic Full", "0", "0", "0", "dev", "r1", "sr1", "ir1"],
        ],
        columns=[
            "id",
            "m49_code",
            "display_name",
            "full_name",
            "is_ldc",
            "is_lldc",
            "is_sids",
            "development_status",
            "m49_region_id",
            "m49_sub_region_id",
            "m49_intermediate_region_id",
        ],
    ).to_csv(inp["m49_country"], index=False)

    _write_realdb_compatible_lookup_tables(inp)
    # overwrite country again with our duplicate/new rows after compatibility helper
    pd.DataFrame(
        [
            ["id_dup", existing_code, "World_duplicate", "World_duplicate_full", "0", "0", "0", "dev", "r1", "sr1", "ir1"],
            ["id_new", "ZZZ", "Synthetic Test Region", "Synthetic Full", "0", "0", "0", "dev", "r1", "sr1", "ir1"],
        ],
        columns=[
            "id",
            "m49_code",
            "display_name",
            "full_name",
            "is_ldc",
            "is_lldc",
            "is_sids",
            "development_status",
            "m49_region_id",
            "m49_sub_region_id",
            "m49_intermediate_region_id",
        ],
    ).to_csv(inp["m49_country"], index=False)

    pd.DataFrame(
        [["1001", "Host A", "species", "root;A"], ["1002", "Host B", "species", "root;B"]],
        columns=["taxonomy_id", "scientific_name", "rank", "lineage"],
    ).to_csv(inp["host_taxa"], sep="\t", index=False)
    _write_realdb_compatible_genes(inp, real_update_db_copy)

    out_db = _build_db(tmp_path, inp, update=True, update_db=real_update_db_copy)

    conn = sqlite3.connect(str(out_db))
    try:
        cur = conn.cursor()
        cur.execute(
            "SELECT reason FROM update_exclusions WHERE table_name='m49_country' ORDER BY rowid DESC LIMIT 5"
        )
        reasons = [r[0] for r in cur.fetchall()]
        assert "duplicate_key_in_db" in reasons

        cur.execute("SELECT COUNT(*) FROM m49_country WHERE m49_code='ZZZ'")
        assert cur.fetchone()[0] == 1
    finally:
        conn.close()

    conn = sqlite3.connect(str(real_update_db_copy))
    try:
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM update_exclusions")
        assert cur.fetchone()[0] == seed_update_exclusions
    finally:
        conn.close()


def test_update_mode_filtered_ids_keeps_real_master_even_when_listed(tmp_path: Path, real_update_db_copy: Path):
    conn = sqlite3.connect(str(real_update_db_copy))
    try:
        cur = conn.cursor()
        cur.execute(
            "SELECT primary_accession, segment FROM meta_data "
            "WHERE lower(coalesce(accession_type,''))='master' LIMIT 1"
        )
        row = cur.fetchone()
    finally:
        conn.close()

    if not row:
        pytest.skip("No master accession found in real update DB")

    master_acc, master_segment = row[0], str(row[1]) if row[1] is not None else "1"

    inp = _inputs(tmp_path, "real_master_filtered", aln_a="ATGC")

    _write_realdb_compatible_meta(
        inp,
        real_update_db_copy,
        [
            {"primary_accession": master_acc, "segment": master_segment, "accession_type": "master", "exclusion_status": ""},
            {"primary_accession": "Q_REAL_FILTER", "segment": master_segment, "accession_type": "query", "exclusion_status": ""},
        ],
    )

    pd.DataFrame(
        [
            [master_acc, master_acc, "ATGC", master_segment],
            ["Q_REAL_FILTER", master_acc, "AT--", master_segment],
        ],
        columns=["primary_accession", "alignment_name", "alignment", "segment"],
    ).to_csv(inp["aln"], sep="\t", index=False)

    pd.DataFrame(
        [
            [master_acc, master_acc, master_acc, "1", "4", "1", "4", "P", master_segment],
            ["Q_REAL_FILTER", master_acc, master_acc, "1", "4", "1", "4", "P", master_segment],
        ],
        columns=[
            "accession",
            "master_ref_accession",
            "reference_accession",
            "aln_start",
            "aln_end",
            "cds_start",
            "cds_end",
            "product",
            "segment",
        ],
    ).to_csv(inp["features"], sep="\t", index=False)

    pd.DataFrame(
        [[master_acc, master_acc, "ins:2:A", master_segment], ["Q_REAL_FILTER", master_acc, "ins:3:T", master_segment]],
        columns=["primary_accession", "reference", "insertion", "segment"],
    ).to_csv(inp["insertions"], sep="\t", index=False)

    pd.DataFrame(
        [["2001", "Host Master", "species", "root;master"], ["2002", "Host Query", "species", "root;query"]],
        columns=["taxonomy_id", "scientific_name", "rank", "lineage"],
    ).to_csv(inp["host_taxa"], sep="\t", index=False)
    _write_realdb_compatible_genes(inp, real_update_db_copy)

    _write_realdb_compatible_lookup_tables(inp)

    inp["fasta"].write_text(f">{master_acc}\nATGC\n>Q_REAL_FILTER\nATTT\n", encoding="utf-8")

    filtered_ids = tmp_path / "filtered_real_ids.txt"
    filtered_ids.write_text(f"{master_acc}\nQ_REAL_FILTER\n", encoding="utf-8")

    out_db = _build_db(tmp_path, inp, update=True, update_db=real_update_db_copy, filtered_ids_file=filtered_ids)

    conn = sqlite3.connect(str(out_db))
    try:
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM meta_data WHERE primary_accession=?", (master_acc,))
        assert cur.fetchone()[0] >= 1

        cur.execute("SELECT reason FROM excluded_accessions WHERE primary_accession='Q_REAL_FILTER' ORDER BY rowid DESC LIMIT 1")
        q_row = cur.fetchone()
        assert q_row is not None
        assert q_row[0] == "alignment_filtering"

        cur.execute("SELECT COUNT(*) FROM excluded_accessions WHERE primary_accession=?", (master_acc,))
        assert cur.fetchone()[0] == 0
    finally:
        conn.close()


def test_update_mode_fails_on_missing_existing_table_columns(tmp_path: Path, real_update_db_copy: Path):
    inp = _inputs(tmp_path, "missing_cols_realdb", aln_a="ATGC")

    # Deliberately incomplete for existing real DB m49_country schema.
    pd.DataFrame(
        [["001", "World"]],
        columns=["m49_code", "display_name"],
    ).to_csv(inp["m49_country"], index=False)

    with pytest.raises(ValueError, match="missing columns required by existing DB schema"):
        _build_db(tmp_path, inp, update=True, update_db=real_update_db_copy)


def test_update_mode_autofills_missing_cluster_98pct_with_placeholder(tmp_path: Path, real_update_db_copy: Path):
    conn = sqlite3.connect(str(real_update_db_copy))
    try:
        cols = [row[1] for row in conn.execute("PRAGMA table_info(meta_data)").fetchall()]
    finally:
        conn.close()

    if "cluster_98pct" not in cols:
        pytest.skip("Real update DB meta_data does not include cluster_98pct")

    inp = _inputs(tmp_path, "cluster_placeholder_realdb", aln_a="ATGC")
    _write_realdb_compatible_lookup_tables(inp)
    _write_realdb_compatible_genes(inp, real_update_db_copy)

    acc = "Q_CLUSTER_NA"
    _write_realdb_compatible_meta(
        inp,
        real_update_db_copy,
        [
            {"primary_accession": acc, "segment": "1", "accession_type": "query", "exclusion_status": ""},
        ],
    )

    meta_df = pd.read_csv(inp["meta"], sep="\t", dtype=str)
    if "cluster_98pct" in meta_df.columns:
        meta_df = meta_df.drop(columns=["cluster_98pct"])
    meta_df.to_csv(inp["meta"], sep="\t", index=False)

    pd.DataFrame(
        [["3001", "Host Cluster", "species", "root;cluster"]],
        columns=["taxonomy_id", "scientific_name", "rank", "lineage"],
    ).to_csv(inp["host_taxa"], sep="\t", index=False)

    out_db = _build_db(tmp_path, inp, update=True, update_db=real_update_db_copy)

    conn = sqlite3.connect(str(out_db))
    try:
        row = conn.execute(
            "SELECT cluster_98pct FROM meta_data WHERE primary_accession=? ORDER BY rowid DESC LIMIT 1",
            (acc,),
        ).fetchone()
        assert row is not None
        assert row[0] == "NA- see tree"
    finally:
        conn.close()


def test_update_mode_copies_seed_db_to_named_output_without_mutating_seed(tmp_path: Path):
    initial = _inputs(tmp_path, "seed_initial", aln_a="ATGC")
    seed_output = _build_db(tmp_path, initial, update=False)
    seed_db = tmp_path / "seed.db"
    shutil.copyfile(seed_output, seed_db)
    seed_output.unlink()

    update_inputs = _inputs(tmp_path, "seed_update", aln_a="AT--")
    out_db = _build_db(tmp_path, update_inputs, update=True, update_db=seed_db)

    assert out_db.exists()
    assert out_db != seed_db

    seed_conn = sqlite3.connect(str(seed_db))
    try:
        seed_alignment = seed_conn.execute(
            "SELECT alignment FROM sequence_alignment WHERE primary_accession='A' AND segment='1'"
        ).fetchone()[0]
    finally:
        seed_conn.close()

    out_conn = sqlite3.connect(str(out_db))
    try:
        out_alignment = out_conn.execute(
            "SELECT alignment FROM sequence_alignment WHERE primary_accession='A' AND segment='1'"
        ).fetchone()[0]
    finally:
        out_conn.close()

    assert seed_alignment == "ATGC"
    assert out_alignment == "AT--"


def test_update_mode_replaces_existing_usher_tree_with_same_key(tmp_path: Path):
    initial = _inputs(tmp_path, "tree_initial", aln_a="ATGC")
    seed_tree = tmp_path / "seed_tree.nwk"
    seed_tree.write_text("(A:0.1);\n", encoding="utf-8")
    seed_db = _build_db(tmp_path, initial, update=False, usher_tree=seed_tree)

    update_inputs = _inputs(tmp_path, "tree_update", aln_a="AT--")
    updated_tree = tmp_path / "updated_tree.nwk"
    updated_tree.write_text("(A:0.1,B:0.2);\n", encoding="utf-8")
    out_db = _build_db(tmp_path, update_inputs, update=True, update_db=seed_db, usher_tree=updated_tree)

    conn = sqlite3.connect(str(out_db))
    try:
        rows = conn.execute("SELECT name, source, newick FROM trees WHERE source='usher'").fetchall()
        assert len(rows) == 1
        assert rows[0][0] == "usher"
        assert rows[0][2].strip() == "(A:0.1,B:0.2);"
    finally:
        conn.close()


def test_update_mode_uses_meta_exclusion_criteria_when_filtered_files_are_empty(tmp_path: Path):
    initial = _inputs(tmp_path, "meta_excl_initial", aln_a="ATGC")
    seed_db = _build_db(tmp_path, initial, update=False)

    update_inp = _inputs(tmp_path, "meta_excl_update", aln_a="ATGC")
    pd.DataFrame(
        [
            ["REF1", "", "1", "master", "", ""],
            ["Q_MISSING", "", "1", "query", "Unable to align during update", "1"],
        ],
        columns=[
            "primary_accession",
            "exclusion",
            "segment",
            "accession_type",
            "exclusion_criteria",
            "exclusion_status",
        ],
    ).to_csv(update_inp["meta"], sep="\t", index=False)

    pd.DataFrame(
        [["REF1", "REF1", "ATGC", "1"]],
        columns=["primary_accession", "alignment_name", "alignment", "segment"],
    ).to_csv(update_inp["aln"], sep="\t", index=False)

    pd.DataFrame(
        [["REF1", "REF1", "REF1", "1", "4", "1", "4", "P", "1"]],
        columns=[
            "accession",
            "master_ref_accession",
            "reference_accession",
            "aln_start",
            "aln_end",
            "cds_start",
            "cds_end",
            "product",
            "segment",
        ],
    ).to_csv(update_inp["features"], sep="\t", index=False)

    pd.DataFrame(
        [["REF1", "REF1", "ins:2:A", "1"]],
        columns=["primary_accession", "reference", "insertion", "segment"],
    ).to_csv(update_inp["insertions"], sep="\t", index=False)

    update_inp["fasta"].write_text(">REF1\nATGC\n>Q_MISSING\nATTT\n", encoding="utf-8")

    out_db = _build_db(tmp_path, update_inp, update=True, update_db=seed_db)

    conn = sqlite3.connect(str(out_db))
    try:
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM meta_data WHERE primary_accession='Q_MISSING'")
        assert cur.fetchone()[0] == 0

        cur.execute("SELECT reason FROM excluded_accessions WHERE primary_accession='Q_MISSING' ORDER BY rowid DESC LIMIT 1")
        row = cur.fetchone()
        assert row is not None
        assert row[0] == "Unable to align during update"
    finally:
        conn.close()
