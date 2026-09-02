"""CreateSqliteDB's own defences against a versioned accession reaching an identity column.

The pipeline's identity key is the BARE accession ("NC_001542"). The versioned
spelling ("NC_001542.1") is legitimate in exactly one place: meta_data.accession_version,
where GenBankFetcher uses it to notice that GenBank has revised a record.

CreateSqliteDB used to perform no accession normalisation whatsoever - whatever
arrived was stored verbatim after .str.strip(). That was correct only because
every upstream producer happens to agree on bare, and it broke in three silent
ways the moment one of them did not:

(a) a versioned accession in features.tsv or the padded alignment, against a bare
    meta_data, is absent from valid_pairs and the row is filtered away with no
    message - for features that can empty the table;
(b) sequences.header comes straight from Bio.SeqIO's record.id, so a versioned
    FASTA left a record with metadata and no sequence;
(c) on an --update run the UNIQUE indexes are over the raw column, so bare and
    versioned are two distinct keys: the upsert found nothing to conflict with
    and inserted a DUPLICATE row instead of overwriting. The run exited 0.

Every test below fails on the pre-fix code in one of those three ways.
"""

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from CreateSqliteDB import CreateSqliteDB


VERSION_WARNING_MARKER = "arrived carrying a GenBank version suffix"


def _write_tsv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _write_csv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def _inputs(
    tmp_path: Path,
    suffix: str,
    meta_accession="KX148218",
    accession_version="KX148218.1",
    locus="KX148218",
    feature_accession="KX148218",
    reference_accession="NC_001542",
    aln_accession="KX148218",
    alignment_name="NC_001542",
    alignment="ATGC",
    fasta_header="KX148218",
    sequence="ATGC",
):
    """One query record, one segment, spelled independently in each table.

    Every accession is a separate argument on purpose: the failures this file
    pins all come from two tables disagreeing about how to spell the SAME
    record, so a test has to be able to make exactly one of them versioned.
    """
    paths = {
        name: tmp_path / f"{name}_{suffix}{ext}"
        for name, ext in [
            ("meta", ".tsv"),
            ("features", ".tsv"),
            ("aln", ".tsv"),
            ("gene", ".tsv"),
            ("proj", ".tsv"),
            ("insertions", ".tsv"),
            ("host_taxa", ".tsv"),
            ("m49_country", ".csv"),
            ("m49_inter", ".csv"),
            ("m49_region", ".csv"),
            ("m49_sub", ".csv"),
            ("fasta", ".fa"),
        ]
    }

    _write_tsv(
        paths["meta"],
        [[meta_accession, accession_version, locus, "", "1", "query"]],
        ["primary_accession", "accession_version", "locus", "exclusion", "segment", "accession_type"],
    )
    _write_tsv(
        paths["features"],
        [[feature_accession, reference_accession, reference_accession, "1", "4", "1", "4", "1", "4", "P", "1"]],
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
    _write_tsv(
        paths["aln"],
        [[aln_accession, alignment_name, alignment, "1"]],
        ["primary_accession", "alignment_name", "alignment", "segment"],
    )
    _write_tsv(paths["gene"], [["geneA", "Gene A"]], ["name", "description"])
    _write_csv(paths["m49_country"], [["001", "World"]], ["m49_code", "name"])
    _write_csv(paths["m49_inter"], [["X", "Inter"]], ["code", "name"])
    _write_csv(paths["m49_region"], [["Y", "Region"]], ["code", "name"])
    _write_csv(paths["m49_sub"], [["Z", "SubRegion"]], ["code", "name"])
    _write_tsv(paths["proj"], [["Python", "3.11"]], ["Software", "Version"])
    _write_tsv(paths["insertions"], [["A", "R", "ins:5:A", "1"]], ["primary_accession", "reference", "insertion", "segment"])
    _write_tsv(paths["host_taxa"], [["A", "host1"]], ["primary_accession", "host"])
    paths["fasta"].write_text(f">{fasta_header}\n{sequence}\n", encoding="utf-8")

    return paths


def _build_db(tmp_path: Path, inp: dict, db_name: str, update=False, update_db=None) -> Path:
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
        db_name=db_name,
        db_status="last updated" if update else "new db",
        update=update,
        update_db=str(update_db) if update_db else None,
        batch_id="batch_identity",
    )
    db.create_db()
    return tmp_path / "SqliteDB" / f"{db_name}.db"


def _rows(db_path: Path, sql):
    conn = sqlite3.connect(str(db_path))
    try:
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# (a) the valid_pairs membership filter
# ---------------------------------------------------------------------------


def test_a_versioned_features_row_joins_instead_of_vanishing(tmp_path: Path):
    """features.tsv spelled the record 'KX148218.1' while meta_data spelled it
    'KX148218'. valid_pairs is built from meta_data, so the tuple missed and the
    ONLY feature row was dropped - leaving an empty features table and no error."""
    inp = _inputs(tmp_path, "featver", feature_accession="KX148218.1")

    db_path = _build_db(tmp_path, inp, "featver")

    assert _rows(db_path, "SELECT accession FROM features") == [("KX148218",)]


def test_a_versioned_fasta_header_keeps_its_sequence_row(tmp_path: Path):
    """record.id is whatever the FASTA carried, and sequences.header is filtered
    against meta_data's bare accessions - so the sequence disappeared while its
    meta_data row stayed behind, describing a record with no sequence."""
    inp = _inputs(tmp_path, "fastaver", fasta_header="KX148218.1", sequence="ATGG")

    db_path = _build_db(tmp_path, inp, "fastaver")

    assert _rows(db_path, "SELECT header, sequence FROM sequences") == [("KX148218", "ATGG")]


def test_a_versioned_alignment_row_joins_instead_of_vanishing(tmp_path: Path):
    """Same filter, applied to the padded alignment. A dropped alignment row is
    the worst of the three: the record still has metadata and a sequence, so
    nothing downstream looks obviously broken."""
    inp = _inputs(tmp_path, "alnver", aln_accession="KX148218.1")

    db_path = _build_db(tmp_path, inp, "alnver")

    assert _rows(db_path, "SELECT primary_accession, alignment FROM sequence_alignment") == [
        ("KX148218", "ATGC")
    ]


def test_sequence_id_is_normalised_with_the_primary_accession_it_mirrors(tmp_path: Path):
    """_normalize_alignment_columns makes primary_accession a verbatim copy of
    sequence_id, so normalising one and not the other would leave a single row
    disagreeing with itself about its own name."""
    inp = _inputs(tmp_path, "seqid", aln_accession="KX148218.1")

    db_path = _build_db(tmp_path, inp, "seqid")

    assert _rows(db_path, "SELECT primary_accession, sequence_id FROM sequence_alignment") == [
        ("KX148218", "KX148218")
    ]


def test_the_features_reference_columns_are_normalised_too(tmp_path: Path):
    """master_ref_accession and reference_accession are joined against meta_data
    by VerifyMutations, not by this file, so a version leak in either survives
    the valid_pairs filter here and only surfaces much later as a feature map
    with no reference."""
    inp = _inputs(tmp_path, "featref", reference_accession="NC_001542.1")

    db_path = _build_db(tmp_path, inp, "featref")

    assert _rows(db_path, "SELECT master_ref_accession, reference_accession FROM features") == [
        ("NC_001542", "NC_001542")
    ]


def test_the_meta_data_identity_columns_are_normalised(tmp_path: Path):
    """primary_accession and locus are the same value by definition; a producer
    that versions one usually versions both."""
    inp = _inputs(
        tmp_path,
        "metaver",
        meta_accession="KX148218.1",
        locus="KX148218.1",
        feature_accession="KX148218.1",
        aln_accession="KX148218.1",
        fasta_header="KX148218.1",
    )

    db_path = _build_db(tmp_path, inp, "metaver")

    assert _rows(db_path, "SELECT primary_accession, locus FROM meta_data") == [
        ("KX148218", "KX148218")
    ]


# ---------------------------------------------------------------------------
# accession_version - the one column that MUST keep its version
# ---------------------------------------------------------------------------


def test_accession_version_survives_verbatim(tmp_path: Path):
    """The revision marker. Normalise it and every record looks unrevised for
    ever, so an update run never re-fetches a corrected sequence."""
    inp = _inputs(tmp_path, "accver", meta_accession="KX148218.1", accession_version="KX148218.1")

    db_path = _build_db(tmp_path, inp, "accver")

    assert _rows(db_path, "SELECT primary_accession, accession_version FROM meta_data") == [
        ("KX148218", "KX148218.1")
    ]


def test_a_null_identity_column_stays_null(tmp_path: Path):
    """locus and gi_number were never .fillna('')'d before this change, so a
    blanket fillna in the normaliser would rewrite every NULL in every existing
    database into an empty string."""
    inp = _inputs(tmp_path, "nulls", locus="")

    db_path = _build_db(tmp_path, inp, "nulls")

    assert _rows(db_path, "SELECT locus FROM meta_data") == [(None,)]


# ---------------------------------------------------------------------------
# (c) update mode: overwrite, not duplicate
# ---------------------------------------------------------------------------


def test_a_versioned_incoming_accession_updates_the_stored_row(tmp_path: Path):
    """The UNIQUE indexes are over the raw column, so 'KX148218' and 'KX148218.2'
    were two distinct keys: INSERT ... ON CONFLICT found no conflict and simply
    inserted. Every upsert table grew a second, half-populated row for a record
    it already held, nothing logged it, and the run exited 0."""
    seed = _inputs(tmp_path, "seed", alignment="ATGC", sequence="ATGC")
    seed_db = _build_db(tmp_path, seed, "seed_db")
    assert _rows(seed_db, "SELECT COUNT(*) FROM meta_data") == [(1,)]

    revised = _inputs(
        tmp_path,
        "revised",
        meta_accession="KX148218.2",
        accession_version="KX148218.2",
        locus="KX148218.2",
        feature_accession="KX148218.2",
        aln_accession="KX148218.2",
        fasta_header="KX148218.2",
        alignment="ATGA",
        sequence="ATGA",
    )
    merged_db = _build_db(tmp_path, revised, "merged_db", update=True, update_db=seed_db)

    assert _rows(merged_db, "SELECT primary_accession, accession_version FROM meta_data") == [
        ("KX148218", "KX148218.2")
    ]
    assert _rows(merged_db, "SELECT primary_accession, alignment FROM sequence_alignment") == [
        ("KX148218", "ATGA")
    ]
    assert _rows(merged_db, "SELECT header, sequence FROM sequences") == [("KX148218", "ATGA")]
    assert _rows(merged_db, "SELECT COUNT(*) FROM features") == [(1,)]


# ---------------------------------------------------------------------------
# report, do not silently repair
# ---------------------------------------------------------------------------


def test_the_repair_is_reported_with_table_column_and_sample(tmp_path: Path, capsys):
    """A version leak is a bug in whatever produced the file. Absorbing it here
    without a word is how it would survive to the next release, so the warning
    has to name the table, the column and the offending value."""
    inp = _inputs(tmp_path, "warn", feature_accession="KX148218.1")

    _build_db(tmp_path, inp, "warn")

    out = capsys.readouterr().out
    assert VERSION_WARNING_MARKER in out
    assert "features.accession" in out
    assert "KX148218.1 -> KX148218" in out
    assert "accession_version" in out


def test_bare_input_prints_no_warning_and_is_stored_untouched(tmp_path: Path, capsys):
    """Every shipped reference list and every SQLite fixture is already bare, so
    on today's data this whole mechanism must be a no-op that says nothing."""
    inp = _inputs(tmp_path, "quiet")

    db_path = _build_db(tmp_path, inp, "quiet")

    assert VERSION_WARNING_MARKER not in capsys.readouterr().out
    assert _rows(db_path, "SELECT primary_accession, locus, accession_version FROM meta_data") == [
        ("KX148218", "KX148218", "KX148218.1")
    ]
    assert _rows(db_path, "SELECT accession FROM features") == [("KX148218",)]
    assert _rows(db_path, "SELECT header FROM sequences") == [("KX148218",)]


# ---------------------------------------------------------------------------
# conservatism - a label is not an accession
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "identifier",
    [
        "EPI_ISL_402124",       # GISAID, not GenBank
        "NC_002016_4_plus",     # a segmented-run group name
        "KX148218.2021",        # four digits is a year, not a version
    ],
)
def test_labels_pass_through_verbatim_and_quietly(tmp_path: Path, capsys, identifier):
    """This is why normalisation goes through accession_utils rather than
    split('.')[0]: a token that is not accession-shaped is a label, and
    truncating one would break exactly the join it was meant to protect."""
    suffix = identifier.replace("/", "_").replace(".", "_")
    inp = _inputs(
        tmp_path,
        suffix,
        meta_accession=identifier,
        locus=identifier,
        feature_accession=identifier,
        aln_accession=identifier,
        fasta_header=identifier,
    )

    db_path = _build_db(tmp_path, inp, f"label_{suffix}")

    assert VERSION_WARNING_MARKER not in capsys.readouterr().out
    assert _rows(db_path, "SELECT primary_accession FROM meta_data") == [(identifier,)]
    assert _rows(db_path, "SELECT accession FROM features") == [(identifier,)]
    assert _rows(db_path, "SELECT header FROM sequences") == [(identifier,)]


# ---------------------------------------------------------------------------
# Half-normalisation: the comparisons that read a raw file and match it against
# an already-normalised primary_accession.
#
# Normalising one side of a comparison and leaving the other raw is worse than
# normalising neither: the QC exclusion list, the filtered-reason lookup and the
# cluster assignment all match by equality, so a versioned entry that used to
# match would silently stop matching. Each set below is kept a SUPERSET of what
# it held before, so these can only ever match more, never less.
# ---------------------------------------------------------------------------


def _bare_builder(**attrs):
    """A CreateSqliteDB with just the attributes one loader needs."""
    builder = CreateSqliteDB.__new__(CreateSqliteDB)
    for key, value in attrs.items():
        setattr(builder, key, value)
    return builder


def test_filtered_ids_carry_both_spellings(tmp_path: Path):
    ids_file = tmp_path / "filtered_sequences_ids.txt"
    ids_file.write_text("PV547761.1\nNC_001542\n", encoding="utf-8")

    loaded = _bare_builder(filtered_ids_file=str(ids_file))._load_filtered_ids()

    assert "PV547761" in loaded, "the QC exclusion would no longer match the normalised accession"
    assert "PV547761.1" in loaded, "the original spelling must stay - the set may only grow"
    assert "NC_001542" in loaded


def test_filtered_ids_are_unchanged_for_bare_input(tmp_path: Path):
    """The no-op guarantee: every file this pipeline writes today is bare."""
    ids_file = tmp_path / "filtered_sequences_ids.txt"
    ids_file.write_text("NC_001542\nKX148218\nNC_001542_1\n", encoding="utf-8")

    loaded = _bare_builder(filtered_ids_file=str(ids_file))._load_filtered_ids()

    assert loaded == {"NC_001542", "KX148218", "NC_001542_1"}


def test_composite_accession_segment_ids_are_left_alone(tmp_path: Path):
    """An accession may itself contain an underscore, so splitting a composite
    apart is guesswork; it stays comparable because the other side rebuilds the
    composite from an already-normalised accession."""
    ids_file = tmp_path / "filtered_sequences_ids.txt"
    ids_file.write_text("NC_001542_4\n", encoding="utf-8")

    assert _bare_builder(filtered_ids_file=str(ids_file))._load_filtered_ids() == {"NC_001542_4"}


def test_filtered_details_are_reachable_by_the_bare_accession(tmp_path: Path):
    details = tmp_path / "filtered_details.tsv"
    details.write_text(
        "seq_name\terror\treference\nPV547761.1\tnot enough matches\t\n", encoding="utf-8"
    )

    reasons = _bare_builder(filtered_details_file=str(details))._load_filtered_details()

    assert "PV547761" in reasons
    assert reasons["PV547761"] == reasons["PV547761.1"]


def test_filtered_details_bare_input_is_untouched(tmp_path: Path):
    details = tmp_path / "filtered_details.tsv"
    details.write_text("seq_name\terror\treference\nNC_001542\tboom\t\n", encoding="utf-8")

    reasons = _bare_builder(filtered_details_file=str(details))._load_filtered_details()

    assert list(reasons) == ["NC_001542"]


def test_cluster_assignment_survives_a_versioned_cluster_tsv(tmp_path: Path):
    """primary_accession is normalised before the map is applied, so a versioned
    member column would map to nothing and drop every cluster assignment."""
    cluster_tsv = tmp_path / "clusters.tsv"
    cluster_tsv.write_text("REP1\tPV547761.1\nREP1\tNC_001542\n", encoding="utf-8")

    builder = _bare_builder(cluster_tsv=str(cluster_tsv), cluster_min_seq_id=0.95)
    df = pd.DataFrame({"primary_accession": ["PV547761", "NC_001542"]})

    out = builder._add_cluster_column(df)

    assert list(out["cluster_95pct"]) == ["REP1", "REP1"]


def test_cluster_assignment_is_unchanged_for_a_bare_cluster_tsv(tmp_path: Path):
    cluster_tsv = tmp_path / "clusters.tsv"
    cluster_tsv.write_text("REP1\tNC_001542\nREP2\tKX148218\n", encoding="utf-8")

    builder = _bare_builder(cluster_tsv=str(cluster_tsv), cluster_min_seq_id=0.95)
    df = pd.DataFrame({"primary_accession": ["NC_001542", "KX148218"]})

    out = builder._add_cluster_column(df)

    assert list(out["cluster_95pct"]) == ["REP1", "REP2"]


def test_insertions_is_an_identity_table():
    """insertions.primary_accession is filtered against valid_pairs with the
    same membership idiom as features, so it needs the same canonical spelling."""
    assert "insertions" in CreateSqliteDB.IDENTITY_COLUMNS
    assert "primary_accession" in CreateSqliteDB.IDENTITY_COLUMNS["insertions"]


def test_accession_version_is_never_an_identity_column():
    """The revision marker. Normalising it would make every record look
    unrevised for ever, which is the opposite failure to the one being fixed."""
    for columns in CreateSqliteDB.IDENTITY_COLUMNS.values():
        assert "accession_version" not in columns
