import sqlite3
import re
from pathlib import Path

import pandas as pd

from CreateSqliteDB import CreateSqliteDB


def write_tsv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def write_csv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def test_create_sqlite_db_exclusions_clusters_and_trees(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    cluster_tsv = tmp_path / "clusters.tsv"
    filtered_ids = tmp_path / "filtered_ids.txt"
    iqtree = tmp_path / "iqtree.treefile"
    usher = tmp_path / "final-tree.nh"

    write_tsv(
        meta,
        [
            ["A", ""],
            ["B", ""],
            ["C", "manual exclusion reason"],
        ],
        ["primary_accession", "exclusion"],
    )
    write_tsv(features, [["A", "P"]], ["primary_accession", "feature"])
    write_tsv(aln, [["A", "ATGC"]], ["primary_accession", "aligned_seq"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["A", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["A", "host1"]], ["primary_accession", "host"])

    fasta.write_text(">A\nATGC\n>B\nATGA\n", encoding="utf-8")
    cluster_tsv.write_text("REP_A\tA\nREP_B\tB\n", encoding="utf-8")
    filtered_ids.write_text("B\n", encoding="utf-8")
    iqtree.write_text("(A:0.1,B:0.2);\n", encoding="utf-8")
    usher.write_text("(A:0.2,C:0.3);\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="testdb",
        db_status="updated",
        tree_file=None,
        iqtree_file=str(iqtree),
        usher_tree=str(usher),
        cluster_tsv=str(cluster_tsv),
        cluster_min_seq_id="0.95",
        filtered_ids_file=str(filtered_ids),
    )

    db.create_db()

    db_path = tmp_path / "SqliteDB" / "testdb.db"
    assert db_path.exists()

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute(
        "SELECT primary_accession, cluster_95pct, exclusion_status, exclusion_criteria FROM meta_data ORDER BY primary_accession"
    )
    meta_rows = cur.fetchall()
    assert meta_rows == [
        ("A", "REP_A", "", ""),
        ("B", "REP_B", "1", "alignment_filtering"),
        ("C", None, "1", "manual exclusion reason"),
    ]

    cur.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='excluded_accessions'")
    assert cur.fetchone()[0] == 0

    cur.execute("SELECT name, source, newick FROM trees ORDER BY source")
    trees = cur.fetchall()
    assert len(trees) == 2
    assert trees[0][1] == "iqtree"
    assert trees[1][1] == "usher"

    cur.execute("SELECT creation_type FROM info ORDER BY rowid DESC LIMIT 1")
    creation_type = cur.fetchone()[0]
    assert creation_type == "last updated"

    conn.close()


def test_create_sqlite_db_stores_iqtree_when_usher_missing(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    iqtree = tmp_path / "iqtree.treefile"

    write_tsv(meta, [["A", ""]], ["primary_accession", "exclusion"])
    write_tsv(features, [["A", "P"]], ["primary_accession", "feature"])
    write_tsv(aln, [["A", "ATGC"]], ["primary_accession", "aligned_seq"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["A", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["A", "host1"]], ["primary_accession", "host"])

    fasta.write_text(">A\nATGC\n", encoding="utf-8")
    iqtree.write_text("(A:0.1);\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="iqtree_only",
        db_status="new db",
        iqtree_file=str(iqtree),
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "iqtree_only.db")
    cur = conn.cursor()
    cur.execute("SELECT source, newick FROM trees")
    trees = cur.fetchall()
    conn.close()

    assert trees == [("iqtree", "(A:0.1);")]


def test_create_sqlite_db_uses_filtered_details_reason(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    filtered_ids = tmp_path / "filtered_ids.txt"
    filtered_details = tmp_path / "filtered_sequences.tsv"

    write_tsv(meta, [["A", ""], ["B", ""]], ["primary_accession", "exclusion"])
    write_tsv(features, [["A", "P"]], ["primary_accession", "feature"])
    write_tsv(aln, [["A", "ATGC"]], ["primary_accession", "aligned_seq"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["A", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["A", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">A\nATGC\n>B\nATGA\n", encoding="utf-8")

    filtered_ids.write_text("B\n", encoding="utf-8")
    write_tsv(
        filtered_details,
        [["B", "EU747327", "reference not present in master-projected reference_aln; query cannot be projected into merged segment alignment", ""]],
        ["seq_name", "reference", "error", "warnings"],
    )

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="testdb2",
        db_status="new db",
        filtered_ids_file=str(filtered_ids),
        filtered_details_file=str(filtered_details),
    )
    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "testdb2.db")
    cur = conn.cursor()
    cur.execute(
        "SELECT primary_accession, exclusion_status, exclusion_criteria FROM meta_data WHERE primary_accession='B'"
    )
    rows = cur.fetchall()
    cur.execute("SELECT COUNT(*) FROM trees")
    tree_count = cur.fetchone()[0]
    conn.close()

    assert rows == [
        (
            "B",
            "1",
            "alignment_filtering: reference not present in master-projected reference_aln; query cannot be projected into merged segment alignment",
        )
    ]
    assert tree_count == 0


def test_create_sqlite_db_summary_collapses_per_sequence_alignment_errors(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    filtered_ids = tmp_path / "filtered_ids.txt"
    filtered_details = tmp_path / "filtered_sequences.tsv"

    write_tsv(meta, [["Q1", ""], ["Q2", ""]], ["primary_accession", "exclusion"])
    write_tsv(features, [["Q1", "P"]], ["primary_accession", "feature"])
    write_tsv(aln, [["Q1", "ATGC"]], ["primary_accession", "aligned_seq"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["Q1", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["Q1", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">Q1\nATGC\n>Q2\nATGA\n", encoding="utf-8")

    filtered_ids.write_text("Q1\nQ2\n", encoding="utf-8")
    write_tsv(
        filtered_details,
        [
            ["Q1", "REF1", "In sequence #30 'Q1': Unable to align: not enough matches. Details: number of seed matches: 0.", ""],
            ["Q2", "REF1", "In sequence #42 'Q2': Unable to align: not enough matches. Details: number of seed matches: 0.", ""],
        ],
        ["seq_name", "reference", "error", "warnings"],
    )

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="testdb_summary",
        db_status="new db",
        filtered_ids_file=str(filtered_ids),
        filtered_details_file=str(filtered_details),
    )
    db.create_db()

    summary = (tmp_path / "SqliteDB" / "db_summary.txt").read_text(encoding="utf-8")
    assert "alignment_filtering: Unable to align: not enough matches. Details: number of seed matches: 0." in summary
    assert "alignment_filtering                           |" not in summary
    assert "  - 2    : alignment_filtering: Unable to align: not enough matches. Details: number of seed matches: 0." in summary
    assert "In sequence #30 'Q1'" not in summary
    assert "In sequence #42 'Q2'" not in summary


def test_create_sqlite_db_preserves_host_taxa_for_excluded_metadata(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"

    write_tsv(
        meta,
        [
            ["A", "", "9606"],
            ["B", "manual exclusion reason", "9796"],
        ],
        ["primary_accession", "exclusion", "host_taxa_id"],
    )
    write_tsv(features, [["A", "P"]], ["primary_accession", "feature"])
    write_tsv(aln, [["A", "ATGC"]], ["primary_accession", "aligned_seq"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["A", "none"]], ["primary_accession", "insertions"])
    write_tsv(
        host_taxa,
        [["9606", "human", "generic", "species"], ["9796", "horse", "generic", "species"]],
        ["taxa_id", "name", "name_type", "taxonomy_level"],
    )
    fasta.write_text(">A\nATGC\n>B\nATGA\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="testdb_host_taxa",
        db_status="updated",
        tree_file=None,
        iqtree_file=None,
        usher_tree=None,
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "testdb_host_taxa.db")
    cur = conn.cursor()
    cur.execute("SELECT taxa_id FROM host_taxa ORDER BY taxa_id")
    rows = cur.fetchall()
    conn.close()

    assert rows == [("9606",), ("9796",)]


def test_create_sqlite_db_maps_tree_manifest_segment_from_refset_key(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    tree_manifest = tmp_path / "tree_manifest.tsv"
    seg_tree = tmp_path / "seg1.treefile"

    write_tsv(meta, [["A", "", "1"]], ["primary_accession", "exclusion", "segment"])
    write_tsv(features, [["A", "P"]], ["primary_accession", "feature"])
    write_tsv(aln, [["A", "ATGC"]], ["primary_accession", "aligned_seq"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["A", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["A", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">A\nATGC\n", encoding="utf-8")

    seg_tree.write_text("(A:0.1);\n", encoding="utf-8")
    write_tsv(
        tree_manifest,
        [["usher", "usher_refset_1", "refset_1_aln_merged_MSA", str(seg_tree)]],
        ["source", "name", "segment_key", "path"],
    )

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="testdb_manifest_segment",
        db_status="new db",
        tree_manifest=str(tree_manifest),
    )
    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "testdb_manifest_segment.db")
    cur = conn.cursor()
    cur.execute("SELECT name, source, segment_key, segment FROM trees WHERE source='usher' LIMIT 1")
    row = cur.fetchone()
    conn.close()

    assert row == ("Usher_tree_full_segment_1", "usher", None, "1")


def test_create_sqlite_db_preserves_multi_segment_rows_and_tree_segments(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    tree_manifest = tmp_path / "tree_manifest.tsv"
    seg1_tree = tmp_path / "seg1.treefile"
    seg2_tree = tmp_path / "seg2.treefile"

    write_tsv(
        meta,
        [
            ["MASTER1", "master", "1"],
            ["MASTER2", "master", "2"],
            ["Q_SEG1", "query", "1"],
            ["Q_SEG2", "query", "2"],
        ],
        ["primary_accession", "accession_type", "segment"],
    )
    write_tsv(
        features,
        [
            ["MASTER1", "MASTER1", "MASTER1", "1", "4", "1", "4", "1", "4", "P1", "1"],
            ["MASTER2", "MASTER2", "MASTER2", "1", "4", "1", "4", "1", "4", "P2", "2"],
            ["Q_SEG1", "MASTER1", "MASTER1", "1", "4", "1", "4", "1", "4", "P1", "1"],
            ["Q_SEG2", "MASTER2", "MASTER2", "1", "4", "1", "4", "1", "4", "P2", "2"],
        ],
        [
            "accession",
            "master_ref_accession",
            "reference_accession",
            "aln_start",
            "aln_end",
            "cds_start",
            "cds_end",
            "cds_start_OG_seq",
            "cds_end_OG_seq",
            "product",
            "segment",
        ],
    )
    write_tsv(
        aln,
        [
            ["MASTER1", "MASTER1", "ATGC", "1"],
            ["MASTER2", "MASTER2", "ATGC", "2"],
            ["Q_SEG1", "MASTER1", "AT-T", "1"],
            ["Q_SEG2", "MASTER2", "A-GC", "2"],
        ],
        ["primary_accession", "alignment_name", "alignment", "segment"],
    )
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(
        insertions,
        [["Q_SEG1", "MASTER1", "ins:2:A", "1"], ["Q_SEG2", "MASTER2", "ins:3:T", "2"]],
        ["primary_accession", "reference", "insertion", "segment"],
    )
    write_tsv(host_taxa, [["MASTER1", "host1"], ["MASTER2", "host2"]], ["primary_accession", "host"])
    fasta.write_text(
        ">MASTER1\nATGC\n>MASTER2\nATGC\n>Q_SEG1\nATTT\n>Q_SEG2\nAGGC\n",
        encoding="utf-8",
    )

    seg1_tree.write_text("(MASTER1:0.1,Q_SEG1:0.2);\n", encoding="utf-8")
    seg2_tree.write_text("(MASTER2:0.1,Q_SEG2:0.2);\n", encoding="utf-8")
    write_tsv(
        tree_manifest,
        [
            ["usher", "usher_refset_1", "refset_1_aln_merged_MSA", str(seg1_tree)],
            ["usher", "usher_refset_2", "refset_2_aln_merged_MSA", str(seg2_tree)],
        ],
        ["source", "name", "segment_key", "path"],
    )

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="testdb_segmented_rows",
        db_status="new db",
        tree_manifest=str(tree_manifest),
    )
    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "testdb_segmented_rows.db")
    cur = conn.cursor()
    cur.execute("SELECT primary_accession, segment FROM meta_data WHERE accession_type='query' ORDER BY primary_accession")
    meta_rows = cur.fetchall()
    cur.execute("SELECT primary_accession, segment FROM sequence_alignment WHERE primary_accession LIKE 'Q_%' ORDER BY primary_accession")
    aln_rows = cur.fetchall()
    cur.execute("SELECT accession, segment FROM features WHERE accession LIKE 'Q_%' ORDER BY accession")
    feature_rows = cur.fetchall()
    cur.execute("SELECT segment_key, segment FROM trees ORDER BY segment_key")
    tree_rows = cur.fetchall()
    conn.close()

    assert meta_rows == [("Q_SEG1", "1"), ("Q_SEG2", "2")]
    assert aln_rows == [("Q_SEG1", "1"), ("Q_SEG2", "2")]
    assert feature_rows == [("Q_SEG1", "1"), ("Q_SEG2", "2")]
    assert tree_rows == [("refset_1_aln_merged_MSA", "1"), ("refset_2_aln_merged_MSA", "2")]


def test_create_sqlite_db_preserves_segmented_rows_across_tables(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    tree_manifest = tmp_path / "tree_manifest.tsv"
    usher_seg1 = tmp_path / "usher_seg1.tree"
    usher_seg2 = tmp_path / "usher_seg2.tree"

    write_tsv(
        meta,
        [
            ["MASTER1", "master", "1", ""],
            ["MASTER2", "master", "2", ""],
            ["Q1", "query", "1", ""],
            ["Q2", "query", "2", ""],
        ],
        ["primary_accession", "accession_type", "segment", "exclusion"],
    )
    write_tsv(
        features,
        [
            ["MASTER1", "MASTER1", "MASTER1", "1", "4", "1", "4", "1", "4", "P1", "1"],
            ["MASTER2", "MASTER2", "MASTER2", "1", "4", "1", "4", "1", "4", "P2", "2"],
            ["Q1", "MASTER1", "MASTER1", "1", "4", "1", "4", "1", "4", "P1", "1"],
            ["Q2", "MASTER2", "MASTER2", "1", "4", "1", "4", "1", "4", "P2", "2"],
        ],
        [
            "accession",
            "master_ref_accession",
            "reference_accession",
            "aln_start",
            "aln_end",
            "cds_start",
            "cds_end",
            "cds_start_OG_seq",
            "cds_end_OG_seq",
            "product",
            "segment",
        ],
    )
    write_tsv(
        aln,
        [
            ["MASTER1", "MASTER1", "ATGC", "1"],
            ["MASTER2", "MASTER2", "ATGA", "2"],
            ["Q1", "MASTER1", "AT-T", "1"],
            ["Q2", "MASTER2", "A-GA", "2"],
        ],
        ["primary_accession", "alignment_name", "alignment", "segment"],
    )
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(
        insertions,
        [["Q1", "MASTER1", "ins:2:T", "1"], ["Q2", "MASTER2", "ins:3:G", "2"]],
        ["primary_accession", "reference", "insertion", "segment"],
    )
    write_tsv(host_taxa, [["MASTER1", "host1"], ["MASTER2", "host2"]], ["primary_accession", "host"])
    fasta.write_text(
        ">MASTER1\nATGC\n"
        ">MASTER2\nATGA\n"
        ">Q1\nATTT\n"
        ">Q2\nAGGA\n",
        encoding="utf-8",
    )

    usher_seg1.write_text("(MASTER1:0.1,Q1:0.2);\n", encoding="utf-8")
    usher_seg2.write_text("(MASTER2:0.1,Q2:0.2);\n", encoding="utf-8")
    write_tsv(
        tree_manifest,
        [
            ["usher", "usher_refset_1", "refset_1_aln_merged_MSA", str(usher_seg1)],
            ["usher", "usher_refset_2", "refset_2_aln_merged_MSA", str(usher_seg2)],
        ],
        ["source", "name", "segment_key", "path"],
    )

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="testdb_segmented_tables",
        db_status="new db",
        tree_manifest=str(tree_manifest),
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "testdb_segmented_tables.db")
    cur = conn.cursor()
    cur.execute("SELECT primary_accession, segment FROM meta_data ORDER BY primary_accession")
    meta_rows = cur.fetchall()
    cur.execute("SELECT primary_accession, segment FROM sequence_alignment ORDER BY primary_accession")
    aln_rows = cur.fetchall()
    cur.execute("SELECT accession, segment FROM features ORDER BY accession")
    feature_rows = cur.fetchall()
    cur.execute("SELECT primary_accession, segment FROM insertions ORDER BY primary_accession")
    insertion_rows = cur.fetchall()
    cur.execute("SELECT name, segment_key, segment FROM trees ORDER BY segment")
    tree_rows = cur.fetchall()
    conn.close()

    assert meta_rows == [("MASTER1", "1"), ("MASTER2", "2"), ("Q1", "1"), ("Q2", "2")]
    assert aln_rows == [("MASTER1", "1"), ("MASTER2", "2"), ("Q1", "1"), ("Q2", "2")]
    assert feature_rows == [("MASTER1", "1"), ("MASTER2", "2"), ("Q1", "1"), ("Q2", "2")]
    assert insertion_rows == [("Q1", "1"), ("Q2", "2")]
    assert tree_rows == [
        ("usher_refset_1", "refset_1_aln_merged_MSA", "1"),
        ("usher_refset_2", "refset_2_aln_merged_MSA", "2"),
    ]


def test_create_sqlite_db_summary_reports_segment_inclusion_exclusion_counts(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"

    write_tsv(
        meta,
        [
            ["REF_SEG1", "", "1", "reference"],
            ["REF_SEG2", "", "2", "reference"],
            ["Q_SEG1_OK", "", "1", "query"],
            ["Q_SEG1_EX", "manual exclusion reason", "1", "query"],
            ["Q_SEG2_OK", "", "2", "query"],
            ["Q_SEG2_EX", "manual exclusion reason", "2", "query"],
        ],
        ["primary_accession", "exclusion", "segment", "accession_type"],
    )
    write_tsv(features, [["Q_SEG1_OK", "P1", "1"], ["Q_SEG2_OK", "P2", "2"]], ["primary_accession", "feature", "segment"])
    write_tsv(aln, [["Q_SEG1_OK", "ATGC", "1"], ["Q_SEG2_OK", "ATGC", "2"]], ["primary_accession", "aligned_seq", "segment"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["Q_SEG1_OK", "none", "1"], ["Q_SEG2_OK", "none", "2"]], ["primary_accession", "insertions", "segment"])
    write_tsv(host_taxa, [["Q_SEG1_OK", "host1"]], ["primary_accession", "host"])
    fasta.write_text(
        ">REF_SEG1\nATGC\n"
        ">REF_SEG2\nATGC\n"
        ">Q_SEG1_OK\nATGC\n"
        ">Q_SEG1_EX\nATGC\n"
        ">Q_SEG2_OK\nATGC\n"
        ">Q_SEG2_EX\nATGC\n",
        encoding="utf-8",
    )

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="testdb_segment_summary",
        db_status="new db",
    )
    db.create_db()

    summary = (tmp_path / "SqliteDB" / "db_summary.txt").read_text(encoding="utf-8")
    assert "[CreateSqliteDB] Query Sequences by Segment (passed/failed QC):" in summary
    assert "Query sequences passing QC" in summary

    segment_block = summary.split("[CreateSqliteDB] Query Sequences by Segment (passed/failed QC):", 1)[1]
    assert re.search(r"(?m)^1\s*\|\s*1\s*\|\s*1\s*$", segment_block)
    assert re.search(r"(?m)^2\s*\|\s*1\s*\|\s*1\s*$", segment_block)


def test_create_sqlite_db_preserves_projected_gff_feature_columns(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"

    write_tsv(
        meta,
        [["NC_004102", "master"], ["EU781827", "reference"], ["seq1", ""]],
        ["primary_accession", "accession_type"],
    )
    write_tsv(
        features,
        [
            ["NC_004102", "NC_004102", "NC_004102", "1", "9", "4", "9", "4", "9", "NS3"],
            ["EU781827", "NC_004102", "EU781827", "1", "9", "4", "9", "4", "9", "NS3"],
            ["seq1", "NC_004102", "EU781827", "1", "9", "4", "9", "4", "9", "NS3"],
        ],
        ["accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end", "cds_start", "cds_end", "cds_start_OG_seq", "cds_end_OG_seq", "product"],
    )
    write_tsv(
        aln,
        [["NC_004102", "NC_004102", "ATGC"], ["EU781827", "EU781827", "ATGC"], ["seq1", "EU781827", "ATGC"]],
        ["primary_accession", "alignment_name", "aligned_seq"],
    )
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["NC_004102", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["NC_004102", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">NC_004102\nATGC\n>EU781827\nATGC\n>seq1\nATGC\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="projected_features",
        db_status="new db",
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "projected_features.db")
    cur = conn.cursor()
    cur.execute("PRAGMA table_info(features)")
    columns = [row[1] for row in cur.fetchall()]
    cur.execute(
        "SELECT accession, master_ref_accession, reference_accession, cds_start, cds_end, product FROM features ORDER BY accession"
    )
    rows = cur.fetchall()
    conn.close()

    assert columns == [
        'accession',
        'master_ref_accession',
        'reference_accession',
        'aln_start',
        'aln_end',
        'cds_start',
        'cds_end',
        'cds_start_OG_seq',
        'cds_end_OG_seq',
        'product',
        'segment',
    ]
    assert rows == [
        ('EU781827', 'NC_004102', 'EU781827', '4', '9', 'NS3'),
        ('NC_004102', 'NC_004102', 'NC_004102', '4', '9', 'NS3'),
        ('seq1', 'NC_004102', 'EU781827', '4', '9', 'NS3'),
    ]


def test_create_sqlite_db_raises_when_meta_file_missing(tmp_path: Path):
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"

    write_tsv(features, [["A", "P"]], ["primary_accession", "feature"])
    write_tsv(aln, [["A", "ATGC"]], ["primary_accession", "aligned_seq"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["A", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["A", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">A\nATGC\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(tmp_path / "missing_meta.tsv"),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="bad",
        db_status="new db",
    )

    try:
        db.create_db()
        assert False, "Expected FileNotFoundError"
    except FileNotFoundError as exc:
        assert "meta_data file not found" in str(exc)


def test_create_sqlite_db_raises_when_alignment_missing_primary_accession(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"

    write_tsv(meta, [["A", ""]], ["primary_accession", "exclusion"])
    write_tsv(features, [["A", "P"]], ["primary_accession", "feature"])
    write_tsv(aln, [["ATGC"]], ["aligned_seq"])
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["A", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["A", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">A\nATGC\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="bad2",
        db_status="new db",
    )

    try:
        db.create_db()
        assert False, "Expected ValueError"
    except ValueError as exc:
        assert "pad_aln is missing required columns" in str(exc)


def test_create_sqlite_db_adds_nearest_reference_metadata_columns(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    reference_tsv = tmp_path / "reference.tsv"

    write_tsv(
        meta,
        [["NC_004102", "master"], ["EU781827", "reference"], ["seq1", ""], ["seq2", ""]],
        ["primary_accession", "accession_type"],
    )
    write_tsv(features, [["NC_004102", "P"]], ["primary_accession", "feature"])
    write_tsv(
        aln,
        [["NC_004102", "NC_004102", "ATGC"], ["EU781827", "EU781827", "ATGC"], ["seq1", "NC_004102", "ATGC"], ["seq2", "EU781827", "ATGC"]],
        ["primary_accession", "alignment_name", "aligned_seq"],
    )
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["NC_004102", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["NC_004102", "host1"]], ["primary_accession", "host"])
    fasta.write_text(
        ">NC_004102\nATGC\n>EU781827\nATGC\n>seq1\nATGC\n>seq2\nATGC\n",
        encoding="utf-8",
    )
    write_tsv(
        reference_tsv,
        [["NC_004102", "master", "1", "1", "a"], ["EU781827", "reference", "1", "1", "NA"]],
        ["primary_accession", "status", "segment", "genotype", "subtype"],
    )

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="with_reference_metadata",
        db_status="new db",
        reference_tsv=str(reference_tsv),
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "with_reference_metadata.db")
    cur = conn.cursor()
    cur.execute(
        "SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype FROM meta_data ORDER BY primary_accession"
    )
    rows = cur.fetchall()
    conn.close()

    assert rows == [
        ("EU781827", "1", ""),
        ("NC_004102", "1", "a"),
        ("seq1", "1", "a"),
        ("seq2", "1", ""),
    ]


def test_create_sqlite_db_update_adds_reference_metadata_columns_to_existing_db(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    reference_tsv = tmp_path / "reference.tsv"

    write_tsv(meta, [["NC_004102", "master"], ["seq1", ""]], ["primary_accession", "accession_type"])
    write_tsv(features, [["NC_004102", "P"]], ["primary_accession", "feature"])
    write_tsv(
        aln,
        [["NC_004102", "NC_004102", "ATGC"], ["seq1", "NC_004102", "ATGC"]],
        ["primary_accession", "alignment_name", "aligned_seq"],
    )
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["NC_004102", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["NC_004102", "host1"]], ["primary_accession", "host"])
    fasta.write_text(
        ">NC_004102\nATGC\n>seq1\nATGC\n",
        encoding="utf-8",
    )
    write_tsv(
        reference_tsv,
        [["NC_004102", "master", "1", "1", "a"]],
        ["primary_accession", "status", "segment", "genotype", "subtype"],
    )

    base_db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="update_reference_metadata",
        db_status="new db",
    )
    base_db.create_db()

    updated_db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="update_reference_metadata",
        db_status="last updated",
        reference_tsv=str(reference_tsv),
        update=True,
        update_db=str(tmp_path / "SqliteDB" / "update_reference_metadata.db"),
        batch_id="test_batch",
    )
    updated_db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "update_reference_metadata.db")
    cur = conn.cursor()
    cur.execute("PRAGMA table_info(meta_data)")
    columns = [row[1] for row in cur.fetchall()]
    cur.execute(
        "SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype FROM meta_data ORDER BY primary_accession"
    )
    rows = cur.fetchall()
    conn.close()

    assert "nearest_reference_genotype" in columns
    assert "nearest_reference_subtype" in columns
    assert rows == [("NC_004102", "1", "a"), ("seq1", "1", "a")]


def test_create_sqlite_db_allows_reference_tsv_without_subtype_column(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    reference_tsv = tmp_path / "reference.tsv"

    write_tsv(meta, [["NC_004102", "master"], ["seq1", ""]], ["primary_accession", "accession_type"])
    write_tsv(features, [["NC_004102", "P"]], ["primary_accession", "feature"])
    write_tsv(
        aln,
        [["NC_004102", "NC_004102", "ATGC"], ["seq1", "NC_004102", "ATGC"]],
        ["primary_accession", "alignment_name", "aligned_seq"],
    )
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["NC_004102", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["NC_004102", "host1"]], ["primary_accession", "host"])
    fasta.write_text(
        ">NC_004102\nATGC\n>seq1\nATGC\n",
        encoding="utf-8",
    )
    write_tsv(
        reference_tsv,
        [["NC_004102", "master", "1", "1"]],
        ["primary_accession", "status", "segment", "genotype"],
    )

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="genotype_only_reference_metadata",
        db_status="new db",
        reference_tsv=str(reference_tsv),
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "genotype_only_reference_metadata.db")
    cur = conn.cursor()
    cur.execute(
        "SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype FROM meta_data ORDER BY primary_accession"
    )
    rows = cur.fetchall()
    conn.close()

    assert rows == [("NC_004102", "1", ""), ("seq1", "1", "")]


def test_create_sqlite_db_accepts_headerless_reference_tsv(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    reference_tsv = tmp_path / "reference.tsv"

    write_tsv(meta, [["NC_004102", "master"], ["seq1", ""]], ["primary_accession", "accession_type"])
    write_tsv(features, [["NC_004102", "P"]], ["primary_accession", "feature"])
    write_tsv(
        aln,
        [["NC_004102", "NC_004102", "ATGC"], ["seq1", "NC_004102", "ATGC"]],
        ["primary_accession", "alignment_name", "aligned_seq"],
    )
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["NC_004102", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["NC_004102", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">NC_004102\nATGC\n>seq1\nATGC\n", encoding="utf-8")
    reference_tsv.write_text("NC_004102\tmaster\t1\t1\ta\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="headerless_reference_metadata",
        db_status="new db",
        reference_tsv=str(reference_tsv),
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "headerless_reference_metadata.db")
    cur = conn.cursor()
    cur.execute(
        "SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype FROM meta_data ORDER BY primary_accession"
    )
    rows = cur.fetchall()
    conn.close()

    assert rows == [("NC_004102", "1", "a"), ("seq1", "1", "a")]


def test_create_sqlite_db_with_none_gene_info(tmp_path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "aln.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"

    write_tsv(meta, [["NC_004102", "master"], ["seq1", ""]], ["primary_accession", "accession_type"])
    write_tsv(
        features,
        [
            ["NC_004102", "NC_004102", "NC_004102", "1", "9", "1", "4", "1", "4", "G"],
            ["NC_004102", "NC_004102", "NC_004102", "1", "9", "5", "9", "5", "9", "L"]
        ],
        ["accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end", "cds_start", "cds_end", "cds_start_OG_seq", "cds_end_OG_seq", "product"]
    )
    write_tsv(
        aln,
        [["NC_004102", "NC_004102", "ATGC"], ["seq1", "NC_004102", "ATGC"]],
        ["primary_accession", "alignment_name", "aligned_seq"],
    )
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["NC_004102", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["NC_004102", "host1"]], ["primary_accession", "host"])
    fasta.write_text(">NC_004102\nATGC\n>seq1\nATGC\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=None,
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="test_none_gene_info",
        db_status="new db",
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "test_none_gene_info.db")
    cur = conn.cursor()
    cur.execute("SELECT description, display_name, name, parent_name FROM genes ORDER BY name")
    rows = cur.fetchall()
    conn.close()

    assert set(rows) == {
        ("G", "G", "G", "whole_genome"),
        ("L", "L", "L", "whole_genome"),
        ("Whole genome", "Whole genome", "whole_genome", "NULL")
    }


def test_nearest_reference_prefers_tree_over_blast_hit(tmp_path: Path):
    meta = tmp_path / "meta.tsv"
    features = tmp_path / "features.tsv"
    aln = tmp_path / "sequence_alignment.tsv"
    gene = tmp_path / "gene.tsv"
    m49_country = tmp_path / "m49_country.csv"
    m49_inter = tmp_path / "m49_inter.csv"
    m49_region = tmp_path / "m49_region.csv"
    m49_sub = tmp_path / "m49_sub.csv"
    proj = tmp_path / "software.tsv"
    insertions = tmp_path / "insertions.tsv"
    host_taxa = tmp_path / "host.tsv"
    fasta = tmp_path / "seqs.fa"
    reference_tsv = tmp_path / "reference.tsv"
    usher_tree = tmp_path / "final-tree.nh"

    write_tsv(
        meta,
        [
            ["REF1", "reference"],
            ["REF2", "reference"],
            ["seqTree", "query"],
            ["seqBlast", "query"],
        ],
        ["primary_accession", "accession_type"],
    )
    write_tsv(features, [["REF1", "P"]], ["primary_accession", "feature"])
    # BLAST top hit (alignment_name) points seqTree at REF2 (genotype 2), but the
    # tree places it next to REF1 (genotype 1). seqBlast is absent from the tree.
    write_tsv(
        aln,
        [
            ["REF1", "REF1", "ATGC"],
            ["REF2", "REF2", "ATGC"],
            ["seqTree", "REF2", "ATGC"],
            ["seqBlast", "REF1", "ATGC"],
        ],
        ["primary_accession", "alignment_name", "aligned_seq"],
    )
    write_tsv(gene, [["geneA", "Gene A"]], ["name", "description"])
    write_csv(m49_country, [["001", "World"]], ["m49_code", "name"])
    write_csv(m49_inter, [["X", "Inter"]], ["code", "name"])
    write_csv(m49_region, [["Y", "Region"]], ["code", "name"])
    write_csv(m49_sub, [["Z", "SubRegion"]], ["code", "name"])
    write_tsv(proj, [["Python", "3.11"]], ["Software", "Version"])
    write_tsv(insertions, [["REF1", "none"]], ["primary_accession", "insertions"])
    write_tsv(host_taxa, [["REF1", "host1"]], ["primary_accession", "host"])
    fasta.write_text(
        ">REF1\nATGC\n>REF2\nATGC\n>seqTree\nATGC\n>seqBlast\nATGC\n",
        encoding="utf-8",
    )
    write_tsv(
        reference_tsv,
        [["REF1", "reference", "1", "1", "a"], ["REF2", "reference", "2", "2", "b"]],
        ["primary_accession", "status", "segment", "genotype", "subtype"],
    )
    # seqTree clusters with REF1 (genotype 1) in the tree; seqBlast is not a tip.
    usher_tree.write_text("((REF1:0.1,seqTree:0.1):0.2,REF2:0.2);\n", encoding="utf-8")

    db = CreateSqliteDB(
        meta_data=str(meta),
        features=str(features),
        pad_aln=str(aln),
        gene_info=str(gene),
        m49_countries=str(m49_country),
        m49_interm_region=str(m49_inter),
        m49_regions=str(m49_region),
        m49_sub_regions=str(m49_sub),
        proj_settings=str(proj),
        fasta_sequence_file=str(fasta),
        insertions=str(insertions),
        host_taxa_file=str(host_taxa),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="tree_pref",
        db_status="new db",
        reference_tsv=str(reference_tsv),
        usher_tree=str(usher_tree),
    )

    db.create_db()

    conn = sqlite3.connect(tmp_path / "SqliteDB" / "tree_pref.db")
    cur = conn.cursor()
    cur.execute(
        "SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype FROM meta_data ORDER BY primary_accession"
    )
    rows = cur.fetchall()
    conn.close()

    assert rows == [
        ("REF1", "1", "a"),
        ("REF2", "2", "b"),
        # Absent from the tree -> falls back to the BLAST hit REF1.
        ("seqBlast", "1", "a"),
        # Tree neighbour REF1 wins over the BLAST hit REF2.
        ("seqTree", "1", "a"),
    ]

