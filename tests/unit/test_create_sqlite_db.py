import sqlite3
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
    cur.execute("SELECT source, segment_key, segment FROM trees WHERE source='usher' LIMIT 1")
    row = cur.fetchone()
    conn.close()

    assert row == ("usher", "refset_1_aln_merged_MSA", "1")


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
            ["NC_004102", "NC_004102", "NC_004102", "1", "9", "4", "9", "NS3"],
            ["EU781827", "NC_004102", "EU781827", "1", "9", "4", "9", "NS3"],
            ["seq1", "NC_004102", "EU781827", "1", "9", "4", "9", "NS3"],
        ],
        ["accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end", "cds_start", "cds_end", "product"],
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
        'product',
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
