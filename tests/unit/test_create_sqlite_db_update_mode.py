import sqlite3
from pathlib import Path

import pandas as pd

from CreateSqliteDB import CreateSqliteDB


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


def _build_db(tmp_path: Path, inp: dict, update=False, update_db=None):
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
        update=update,
        update_db=str(update_db) if update_db else None,
        batch_id="batch_test",
    )
    db.create_db()


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
