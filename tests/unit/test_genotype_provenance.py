"""Provenance for the genotype/subtype labels stored in meta_data.

``nearest_reference_genotype`` / ``nearest_reference_subtype`` are resolved
through a chain of fallbacks - the sequence's own curated reference-list entry,
then the phylogenetic tree neighbourhood, then an EPA-ng placement, then the
BLAST top hit - and until now the database recorded only the answer.  A stored
``6a`` on a reference the curators deliberately left as "genotype 6, subtype not
assigned" was indistinguishable from a curated ``6a``: the letter came from
whatever the sequence happened to align against.  69 of the 238 references in
``generic/hcv/ref_list_subtype_genotype.txt`` are in that state, and 437 of the
497 rows in ``generic/rabv/ref_list_clades.txt`` carry "NA" for both fields.

``genotype_origin`` / ``subtype_origin`` name the source that produced each
value, from the fixed vocabulary::

    curated_reflist | tree_usher | tree_iqtree | epa_placement | blast_tophit
    | gisaid_declared | ncbi_declared | unresolved

Filling a curated "NA" subtype by inference is still allowed - that behaviour is
deliberate - and the origin column is what makes it visible.
"""

import csv
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from CreateSqliteDB import CreateSqliteDB


REPO_ROOT = Path(__file__).resolve().parents[2]
HCV_DB = REPO_ROOT / "test_out" / "HCV_OM_test" / "HCV_OM_test.db"

VOCABULARY = {
    "curated_reflist", "tree_usher", "tree_iqtree", "epa_placement",
    "blast_tophit", "gisaid_declared", "ncbi_declared", "unresolved",
}


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _ref_list(path: Path, rows) -> Path:
    """Reference list in the shipped 5-column shape:
    primary_accession / accession_type / segment / genotype / subtype."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for row in rows:
            writer.writerow(row)
    return path


def _db(tmp_path: Path, **kwargs) -> CreateSqliteDB:
    """CreateSqliteDB carrying only the clade-relevant inputs.

    ``_add_reference_columns`` touches none of the other input files, so
    placeholders are enough and no build has to run.  Keyword arguments only, so
    a new constructor parameter cannot silently shift a positional here.
    """
    base = dict(
        meta_data="unused", features="unused", pad_aln="unused", gene_info=None,
        m49_countries="unused", m49_interm_region="unused", m49_regions="unused",
        m49_sub_regions="unused", proj_settings="unused", fasta_sequence_file="unused",
        insertions="unused", host_taxa_file="unused", base_dir=str(tmp_path),
        output_dir="SqliteDB", db_name="provenance", db_status="new db",
    )
    base.update(kwargs)
    return CreateSqliteDB(**base)


def _meta(accessions, **columns):
    frame = pd.DataFrame({"primary_accession": list(accessions)})
    for name, values in columns.items():
        frame[name] = values
    return frame


def _aln(pairs):
    if not pairs:
        return pd.DataFrame(columns=["primary_accession", "alignment_name"])
    return pd.DataFrame(pairs, columns=["primary_accession", "alignment_name"])


def _labels(df, accession):
    hit = df[df["primary_accession"] == accession].iloc[0]
    return (
        hit["nearest_reference_genotype"],
        hit["genotype_origin"],
        hit["nearest_reference_subtype"],
        hit["subtype_origin"],
    )


# --------------------------------------------------------------------------
# one test per source in the vocabulary
# --------------------------------------------------------------------------

def test_a_curated_reference_records_curated_reflist_for_both_fields(tmp_path: Path):
    """A reference's own ref_list row is the top of the chain and must say so.

    Real trigger: a curated dataset must not drift on rebuild.  With the origin
    recorded, "this label was curated" is a query rather than an assumption.
    """
    ref_list = _ref_list(tmp_path / "ref.txt", [["R1a", "reference", "1", "1", "1a"]])
    out = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(
        _meta(["R1a"]), _aln([])
    )
    assert _labels(out, "R1a") == ("1", "curated_reflist", "1a", "curated_reflist")


def test_a_query_placed_by_the_usher_tree_records_tree_usher(tmp_path: Path):
    """UShER holds every placed sample, so it is the preferred inference tier."""
    ref_list = _ref_list(tmp_path / "ref.txt", [
        ["R1a", "reference", "1", "1", "1a"],
        ["R2b", "reference", "1", "2", "2b"],
    ])
    usher = tmp_path / "usher.nh"
    usher.write_text("((Q1:0.1,R1a:0.1):0.1,R2b:0.5);\n", encoding="utf-8")
    out = _db(tmp_path, reference_tsv=str(ref_list), usher_tree=str(usher))._add_reference_columns(
        _meta(["Q1"]), _aln([("Q1", "R2b")])
    )
    assert _labels(out, "Q1") == ("1", "tree_usher", "1a", "tree_usher")


def test_a_query_placed_only_by_the_iqtree_backbone_records_tree_iqtree(tmp_path: Path):
    """The IQ-TREE backbone holds cluster representatives only, so it is a
    distinct - and distinguishable - tier from UShER."""
    ref_list = _ref_list(tmp_path / "ref.txt", [
        ["R1a", "reference", "1", "1", "1a"],
        ["R2b", "reference", "1", "2", "2b"],
    ])
    iqtree = tmp_path / "iqtree.treefile"
    iqtree.write_text("((Q1:0.1,R1a:0.1):0.1,R2b:0.5);\n", encoding="utf-8")
    out = _db(tmp_path, reference_tsv=str(ref_list), iqtree_file=str(iqtree))._add_reference_columns(
        _meta(["Q1"]), _aln([("Q1", "R2b")])
    )
    assert _labels(out, "Q1") == ("1", "tree_iqtree", "1a", "tree_iqtree")


def test_an_epa_ng_placement_records_epa_placement(tmp_path: Path):
    """EPA-ng only runs when the run produced no usable tree; a database built
    that way must not look like a tree-derived one."""
    ref_list = _ref_list(tmp_path / "ref.txt", [["R2b", "reference", "1", "2", "2b"]])
    clade_tsv = tmp_path / "clade_assignments.tsv"
    clade_tsv.write_text(
        "primary_accession\tgenotype\tsubtype\nQ1\t3\t3c\n", encoding="utf-8"
    )
    out = _db(
        tmp_path, reference_tsv=str(ref_list), clade_assignments=str(clade_tsv)
    )._add_reference_columns(_meta(["Q1"]), _aln([("Q1", "R2b")]))
    assert _labels(out, "Q1") == ("3", "epa_placement", "3c", "epa_placement")


def test_the_blast_top_hit_records_blast_tophit(tmp_path: Path):
    """--tree_free skips phylogenetics and the EPA-ng fallback, so alignment_name
    is the only source left and the weakest evidence in the chain."""
    ref_list = _ref_list(tmp_path / "ref.txt", [["R2b", "reference", "1", "2", "2b"]])
    out = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(
        _meta(["Q1"]), _aln([("Q1", "R2b")])
    )
    assert _labels(out, "Q1") == ("2", "blast_tophit", "2b", "blast_tophit")


def test_a_declared_gisaid_genotype_records_gisaid_declared(tmp_path: Path):
    """Vendor metadata is a statement, not an inference, and is ranked below
    everything this pipeline derives itself."""
    ref_list = _ref_list(tmp_path / "ref.txt", [["R2b", "reference", "1", "2", "2b"]])
    meta = _meta(["Q1"], gisaid_genotype=["7"], gisaid_subtype=["7x"])
    out = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(meta, _aln([]))
    assert _labels(out, "Q1") == ("7", "gisaid_declared", "7x", "gisaid_declared")


def test_a_declared_genbank_genotype_records_ncbi_declared(tmp_path: Path):
    """A bare `genotype`/`subtype` column comes off the GenBank matrix."""
    ref_list = _ref_list(tmp_path / "ref.txt", [["R2b", "reference", "1", "2", "2b"]])
    meta = _meta(["Q1"], genotype=["4"], subtype=["4r"])
    out = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(meta, _aln([]))
    assert _labels(out, "Q1") == ("4", "ncbi_declared", "4r", "ncbi_declared")


def test_gisaid_is_preferred_over_genbank_and_both_lose_to_inference(tmp_path: Path):
    """Precedence, end to end, on one accession that every source can speak to."""
    ref_list = _ref_list(tmp_path / "ref.txt", [
        ["R1a", "reference", "1", "1", "1a"],
        ["R2b", "reference", "1", "2", "2b"],
    ])
    usher = tmp_path / "usher.nh"
    usher.write_text("((Q1:0.1,R1a:0.1):0.1,R2b:0.5);\n", encoding="utf-8")
    meta = _meta(["Q1"], gisaid_genotype=["7"], genotype=["4"])

    with_tree = _db(
        tmp_path, reference_tsv=str(ref_list), usher_tree=str(usher)
    )._add_reference_columns(meta.copy(), _aln([("Q1", "R2b")]))
    assert _labels(with_tree, "Q1")[:2] == ("1", "tree_usher")

    without_tree = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(
        meta.copy(), _aln([])
    )
    assert _labels(without_tree, "Q1")[:2] == ("7", "gisaid_declared")


def test_a_sequence_no_source_can_label_is_unresolved_not_blank(tmp_path: Path):
    """"unresolved" is a value, not an empty cell.

    A blank origin would be indistinguishable from "this database predates the
    column", and a column that is 100% empty is exactly what
    validate_meta_column_population is meant to catch.
    """
    ref_list = _ref_list(tmp_path / "ref.txt", [["R2b", "reference", "1", "2", "2b"]])
    out = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(
        _meta(["Q_ORPHAN"]), _aln([])
    )
    assert _labels(out, "Q_ORPHAN") == ("", "unresolved", "", "unresolved")


def test_every_origin_written_is_in_the_agreed_vocabulary(tmp_path: Path):
    """The vocabulary is a contract with downstream consumers; nothing else may
    be written into either column."""
    ref_list = _ref_list(tmp_path / "ref.txt", [
        ["R1a", "reference", "1", "1", "1a"],
        ["R_NA", "reference", "1", "6", "NA"],
    ])
    usher = tmp_path / "usher.nh"
    usher.write_text("((Q1:0.1,R1a:0.1):0.1,R_NA:0.5);\n", encoding="utf-8")
    meta = _meta(["R1a", "R_NA", "Q1", "Q_ORPHAN"], gisaid_subtype=["", "", "", "9z"])
    out = _db(tmp_path, reference_tsv=str(ref_list), usher_tree=str(usher))._add_reference_columns(
        meta, _aln([("Q1", "R1a"), ("R_NA", "R1a")])
    )
    written = set(out["genotype_origin"]) | set(out["subtype_origin"])
    assert written <= VOCABULARY, f"origins outside the agreed vocabulary: {written - VOCABULARY}"


# --------------------------------------------------------------------------
# the point of the columns: making an inference visible
# --------------------------------------------------------------------------

def test_a_curated_na_subtype_filled_by_inference_says_where_the_letter_came_from(tmp_path: Path):
    """The deliberate case: curated subtype "NA" is still treated as missing.

    Real trigger: 69 HCV references are curated genotype-known / subtype-NA and
    come out of a build carrying a subtype letter taken from the reference they
    aligned against.  That behaviour is kept on purpose - but the pair of origins
    now says "genotype curated, subtype inferred from the BLAST hit", which is
    what makes the fabricated letter recognisable.
    """
    ref_list = _ref_list(tmp_path / "ref.txt", [
        ["R1a", "reference", "1", "1", "1a"],
        ["R_NA", "reference", "1", "6", "NA"],
    ])
    out = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(
        _meta(["R_NA"]), _aln([("R_NA", "R1a")])
    )
    genotype, genotype_origin, subtype, subtype_origin = _labels(out, "R_NA")
    assert (genotype, genotype_origin) == ("6", "curated_reflist")
    assert (subtype, subtype_origin) == ("1a", "blast_tophit"), (
        "the inferred subtype must be attributed to the source that invented it"
    )


def test_genotype_and_subtype_are_taken_from_the_same_tree(tmp_path: Path):
    """A label assembled out of two neighbourhoods is supported by neither.

    Real trigger: segmented builds stage one UShER tree per segment plus an
    IQ-TREE backbone, and the same accession is routinely in several of them.
    Here UShER knows only a genotype-2 reference (subtype "NA") and the IQ-TREE
    backbone only a genotype-1 one; merging the fields independently produced
    (2, 1b), a pair no reference carries.
    """
    ref_list = _ref_list(tmp_path / "ref.txt", [
        ["R_GENO_ONLY", "reference", "1", "2", "NA"],
        ["R_SUBTYPED", "reference", "1", "1", "1b"],
    ])
    usher = tmp_path / "usher.nh"
    usher.write_text("(Q1:0.1,R_GENO_ONLY:0.1);\n", encoding="utf-8")
    iqtree = tmp_path / "iqtree.treefile"
    iqtree.write_text("(Q1:0.1,R_SUBTYPED:0.1);\n", encoding="utf-8")
    out = _db(
        tmp_path, reference_tsv=str(ref_list), usher_tree=str(usher), iqtree_file=str(iqtree)
    )._add_reference_columns(_meta(["Q1"]), _aln([]))

    genotype, genotype_origin, subtype, subtype_origin = _labels(out, "Q1")
    assert (genotype, genotype_origin) == ("2", "tree_usher")
    assert subtype == "", f"subtype {subtype!r} came from a different tree than the genotype"
    assert subtype_origin == "unresolved"


def test_serotype_is_not_read_as_a_declared_subtype(tmp_path: Path):
    """Influenza's serotype is H1N1, which is not a subtype in this sense.

    Every one of the 518 rows in the shipped IAV database carries a serotype, so
    a loose match on "*type" columns would relabel the whole database from vendor
    metadata and call it a nearest-reference result.
    """
    ref_list = _ref_list(tmp_path / "ref.txt", [["R2b", "reference", "1", "2", "2b"]])
    meta = _meta(["Q1"], serotype=["H1N1"], serotype_validated=["H1N1"])
    out = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(meta, _aln([]))
    assert _labels(out, "Q1") == ("", "unresolved", "", "unresolved")


def test_the_origin_columns_are_not_fed_back_in_on_a_rebuild(tmp_path: Path):
    """meta_data re-exported from a database already carries the label columns.

    They must not be mistaken for vendor-declared genotypes, or the second build
    would attribute everything to ncbi_declared and the chain would freeze.
    """
    ref_list = _ref_list(tmp_path / "ref.txt", [["R2b", "reference", "1", "2", "2b"]])
    meta = _meta(
        ["Q1"],
        nearest_reference_genotype=["9"],
        nearest_reference_subtype=["9z"],
        genotype_origin=["blast_tophit"],
        subtype_origin=["blast_tophit"],
    )
    out = _db(tmp_path, reference_tsv=str(ref_list))._add_reference_columns(meta, _aln([]))
    assert _labels(out, "Q1") == ("", "unresolved", "", "unresolved")


# --------------------------------------------------------------------------
# end to end, through create_db
# --------------------------------------------------------------------------

def _pipeline_inputs(directory: Path, accessions, ref_rows, aln_pairs):
    directory.mkdir(parents=True, exist_ok=True)
    files = {name: directory / filename for name, filename in {
        "meta": "meta_data.tsv", "features": "features.tsv", "aln": "sequence_alignment.tsv",
        "gene": "gene_info.tsv", "country": "m49_country.csv", "interm": "m49_intermediate.csv",
        "region": "m49_region.csv", "sub_region": "m49_sub_region.csv",
        "software": "software_info.tsv", "insertions": "insertions.tsv",
        "host_taxa": "host_taxa.tsv", "fasta": "sequences.fa", "ref_list": "ref_list.txt",
    }.items()}

    def _tsv(path, rows, columns):
        pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)

    def _csv(path, rows, columns):
        pd.DataFrame(rows, columns=columns).to_csv(path, index=False)

    _tsv(files["meta"], [[acc, "1", kind] for acc, kind in accessions],
         ["primary_accession", "segment", "accession_type"])
    _tsv(files["features"],
         [[acc, "M", "R", "1", "10", "1", "10", "1", "10", "P", "1"] for acc, _ in accessions],
         ["accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end",
          "cds_start", "cds_end", "cds_start_OG_seq", "cds_end_OG_seq", "product", "segment"])
    _tsv(files["aln"], [[q, r, "ATGC", "1"] for q, r in aln_pairs],
         ["primary_accession", "alignment_name", "alignment", "segment"])
    _tsv(files["gene"], [["geneA", "Gene A"]], ["name", "description"])
    _csv(files["country"], [["001", "World"]], ["m49_code", "name"])
    _csv(files["interm"], [["X", "Inter"]], ["code", "name"])
    _csv(files["region"], [["Y", "Region"]], ["code", "name"])
    _csv(files["sub_region"], [["Z", "SubRegion"]], ["code", "name"])
    _tsv(files["software"], [["MAFFT", "v7.4"]], ["Software", "Version"])
    _tsv(files["insertions"], [[accessions[0][0], "R", "ins:5:A", "1"]],
         ["primary_accession", "reference", "insertion", "segment"])
    _tsv(files["host_taxa"], [[accessions[0][0], "host1"]], ["primary_accession", "host"])
    files["fasta"].write_text(
        "".join(f">{acc}\nATGC\n" for acc, _ in accessions), encoding="utf-8"
    )
    _ref_list(files["ref_list"], ref_rows)
    return files


def _build(base_dir: Path, files: dict, *, db_name: str, **kwargs) -> Path:
    creator = CreateSqliteDB(
        meta_data=str(files["meta"]), features=str(files["features"]),
        pad_aln=str(files["aln"]), gene_info=str(files["gene"]),
        m49_countries=str(files["country"]), m49_interm_region=str(files["interm"]),
        m49_regions=str(files["region"]), m49_sub_regions=str(files["sub_region"]),
        proj_settings=str(files["software"]), fasta_sequence_file=str(files["fasta"]),
        insertions=str(files["insertions"]), host_taxa_file=str(files["host_taxa"]),
        base_dir=str(base_dir), output_dir="SqliteDB", db_name=db_name,
        db_status=None, reference_tsv=str(files["ref_list"]), **kwargs
    )
    creator.create_db()
    return base_dir / "SqliteDB" / f"{db_name}.db"


def test_a_built_database_carries_both_origin_columns_populated(tmp_path: Path):
    """The columns have to survive the whole build, not just the resolver."""
    files = _pipeline_inputs(
        tmp_path / "in1",
        [("R1a", "reference"), ("Q1", "query")],
        [["R1a", "reference", "1", "1", "1a"]],
        [("Q1", "R1a"), ("R1a", "R1a")],
    )
    db_path = _build(tmp_path, files, db_name="provenance")

    conn = sqlite3.connect(str(db_path))
    try:
        columns = {row[1] for row in conn.execute("PRAGMA table_info(meta_data)")}
        assert {"genotype_origin", "subtype_origin"} <= columns
        stored = {
            row[0]: row[1:]
            for row in conn.execute(
                "SELECT primary_accession, nearest_reference_genotype, genotype_origin, "
                "nearest_reference_subtype, subtype_origin FROM meta_data"
            )
        }
    finally:
        conn.close()

    assert stored["R1a"] == ("1", "curated_reflist", "1a", "curated_reflist")
    assert stored["Q1"] == ("1", "blast_tophit", "1a", "blast_tophit")


def test_an_update_adds_the_origin_columns_to_a_database_that_lacks_them(tmp_path: Path):
    """Every published database predates these columns.

    ``_ensure_update_columns`` has to ALTER them in on the first update, or the
    incoming frame is rejected as carrying columns the schema does not have.
    """
    seed_files = _pipeline_inputs(
        tmp_path / "in1", [("R1a", "reference")], [["R1a", "reference", "1", "1", "1a"]],
        [("R1a", "R1a")],
    )
    seed = _build(tmp_path, seed_files, db_name="seed")
    conn = sqlite3.connect(str(seed))
    try:
        conn.execute("ALTER TABLE meta_data DROP COLUMN genotype_origin")
        conn.execute("ALTER TABLE meta_data DROP COLUMN subtype_origin")
        conn.commit()
        assert "genotype_origin" not in {r[1] for r in conn.execute("PRAGMA table_info(meta_data)")}
    finally:
        conn.close()

    update_files = _pipeline_inputs(
        tmp_path / "in2", [("R1a", "reference"), ("Q1", "query")],
        [["R1a", "reference", "1", "1", "1a"]], [("Q1", "R1a"), ("R1a", "R1a")],
    )
    updated = _build(tmp_path, update_files, db_name="updated", update=True, update_db=str(seed))

    conn = sqlite3.connect(str(updated))
    try:
        stored = dict(
            conn.execute("SELECT primary_accession, genotype_origin FROM meta_data").fetchall()
        )
    finally:
        conn.close()
    assert stored == {"R1a": "curated_reflist", "Q1": "blast_tophit"}


# --------------------------------------------------------------------------
# the shipped HCV database
# --------------------------------------------------------------------------

@pytest.mark.skipif(not HCV_DB.exists(), reason=f"HCV test DB not present: {HCV_DB}")
def test_shipped_hcv_db_shows_why_the_origin_columns_are_needed():
    """Evidence, on a real build, that the answer alone is not enough.

    HCV_OM_test.db carries a genotype for references and queries alike and there
    is no way to tell, from the database, which of them was curated and which was
    inferred from a BLAST hit.  This is read-only: it pins the state the new
    columns exist to describe, and it will keep passing once the database is
    rebuilt with them (the assertion is about what the OLD columns cannot say).
    """
    conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
    try:
        columns = {row[1] for row in conn.execute("PRAGMA table_info(meta_data)")}
        labelled = conn.execute(
            "SELECT COUNT(*) FROM meta_data "
            "WHERE TRIM(COALESCE(nearest_reference_genotype, '')) <> ''"
        ).fetchone()[0]
    finally:
        conn.close()

    assert labelled > 0, "no genotype labels to reason about"
    if not {"genotype_origin", "subtype_origin"} <= columns:
        pytest.skip("shipped HCV DB predates the origin columns; rebuild to populate them")
    conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
    try:
        origins = {
            row[0] for row in conn.execute(
                "SELECT DISTINCT genotype_origin FROM meta_data"
            )
        } | {
            row[0] for row in conn.execute("SELECT DISTINCT subtype_origin FROM meta_data")
        }
    finally:
        conn.close()
    assert origins <= VOCABULARY
