"""Edge cases in ``--update`` (incremental) mode, across all viruses.

Update mode is the path that runs every week against a database that already
holds years of curated data, and it is the one path where a mistake is not a
crash but a quiet corruption of rows nobody re-derives.  Two mechanisms carry
almost all the risk:

* ``merge_table_append_nonredundant`` writes the *upsert* tables with
  ``INSERT OR REPLACE``, which only replaces when a UNIQUE index makes the row
  conflict, and replaces the **whole** row when it does; and
* ``_create_table_unique_indexes`` creates those indexes only in update mode,
  only for three of the five upsert tables, and only when every key column
  happens to be present.

Every test below names the real-world trigger in its docstring.
"""

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from CreateSqliteDB import CreateSqliteDB
from ExportRefListFromUpdateDb import export_ref_list, load_reference_file_table


REPO_ROOT = Path(__file__).resolve().parents[2]
# The RABV database produced by the shipped --update integration run.  Read-only:
# it is evidence about what update mode actually wrote, not a fixture to mutate.
REAL_UPDATE_DB = REPO_ROOT / "test_out" / "update_test" / "rabv-jul0425-update-test.db"


def _ro(db_path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)


# --------------------------------------------------------------------------
# Minimal but realistic pipeline inputs.  Every helper writes under tmp_path.
# --------------------------------------------------------------------------

FEATURE_COLUMNS = [
    "accession", "master_ref_accession", "reference_accession", "aln_start",
    "aln_end", "cds_start", "cds_end", "cds_start_OG_seq", "cds_end_OG_seq",
    "product", "segment",
]
LEGACY_FEATURE_COLUMNS = [
    "accession", "master_ref_accession", "reference_accession", "aln_start",
    "aln_end", "cds_start", "cds_end", "product", "segment",
]


def _tsv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _csv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def _inputs(
    directory: Path,
    *,
    meta_rows=None,
    meta_columns=None,
    feature_rows=None,
    feature_columns=None,
    aln_rows=None,
    insertion_rows=None,
    sequence_rows=None,
    software_rows=None,
    segment="1",
):
    """Write one complete set of CreateSqliteDB inputs into ``directory``."""
    directory.mkdir(parents=True, exist_ok=True)
    files = {
        key: directory / name
        for key, name in {
            "meta": "meta_data.tsv",
            "features": "features.tsv",
            "aln": "sequence_alignment.tsv",
            "gene": "gene_info.tsv",
            "country": "m49_country.csv",
            "interm": "m49_intermediate.csv",
            "region": "m49_region.csv",
            "sub_region": "m49_sub_region.csv",
            "software": "software_info.tsv",
            "insertions": "insertions.tsv",
            "host_taxa": "host_taxa.tsv",
            "fasta": "sequences.fa",
        }.items()
    }

    _tsv(
        files["meta"],
        meta_rows if meta_rows is not None else [["A", segment, "Bos taurus"]],
        meta_columns or ["primary_accession", "segment", "host"],
    )
    _tsv(
        files["features"],
        feature_rows
        if feature_rows is not None
        else [["A", "M", "R", "1", "10", "1", "10", "1", "10", "P", segment]],
        feature_columns or FEATURE_COLUMNS,
    )
    _tsv(
        files["aln"],
        aln_rows if aln_rows is not None else [["A", "R", "ATGC", segment]],
        ["primary_accession", "alignment_name", "alignment", "segment"],
    )
    _tsv(files["gene"], [["geneA", "Gene A"]], ["name", "description"])
    _csv(files["country"], [["001", "World"]], ["m49_code", "name"])
    _csv(files["interm"], [["X", "Inter"]], ["code", "name"])
    _csv(files["region"], [["Y", "Region"]], ["code", "name"])
    _csv(files["sub_region"], [["Z", "SubRegion"]], ["code", "name"])
    _tsv(
        files["software"],
        software_rows if software_rows is not None else [["MAFFT", "v7.4"]],
        ["Software", "Version"],
    )
    _tsv(
        files["insertions"],
        insertion_rows if insertion_rows is not None else [["A", "R", "ins:5:A", segment]],
        ["primary_accession", "reference", "insertion", "segment"],
    )
    _tsv(files["host_taxa"], [["A", "host1"]], ["primary_accession", "host"])
    rows = sequence_rows if sequence_rows is not None else [("A", "ATGC")]
    files["fasta"].write_text("".join(f">{h}\n{s}\n" for h, s in rows), encoding="utf-8")
    return files


def _build(
    base_dir: Path,
    files: dict,
    *,
    db_name: str,
    update: bool = False,
    update_db: Path | None = None,
    batch_id: str = "batch_test",
    db_status: str = "new db",
    cluster_tsv: Path | None = None,
    cluster_min_seq_id: str | None = None,
) -> Path:
    creator = CreateSqliteDB(
        meta_data=str(files["meta"]),
        features=str(files["features"]),
        pad_aln=str(files["aln"]),
        gene_info=str(files["gene"]),
        m49_countries=str(files["country"]),
        m49_interm_region=str(files["interm"]),
        m49_regions=str(files["region"]),
        m49_sub_regions=str(files["sub_region"]),
        proj_settings=str(files["software"]),
        fasta_sequence_file=str(files["fasta"]),
        insertions=str(files["insertions"]),
        host_taxa_file=str(files["host_taxa"]),
        base_dir=str(base_dir),
        output_dir="SqliteDB",
        db_name=db_name,
        db_status=db_status,
        cluster_tsv=str(cluster_tsv) if cluster_tsv else None,
        cluster_min_seq_id=cluster_min_seq_id,
        update=update,
        update_db=str(update_db) if update_db else None,
        batch_id=batch_id,
    )
    creator.create_db()
    return base_dir / "SqliteDB" / f"{db_name}.db"


def _rows(db_path: Path, sql: str):
    conn = sqlite3.connect(str(db_path))
    try:
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


def _count(db_path: Path, table: str) -> int:
    return _rows(db_path, f"SELECT COUNT(*) FROM {table}")[0][0]


# --------------------------------------------------------------------------
# 1. Upsert tables with no unique index
# --------------------------------------------------------------------------


@pytest.mark.xfail(
    reason="sequences and insertions are in merge_table_append_nonredundant's "
           "upsert_tables and are written with INSERT OR REPLACE, but "
           "_create_table_unique_indexes only builds indexes for features, "
           "sequence_alignment and meta_data. With no UNIQUE index nothing can "
           "conflict, so INSERT OR REPLACE degrades to a plain INSERT and the "
           "row is appended alongside the old one.",
    strict=False,
)
def test_update_replaces_sequence_text_instead_of_appending_a_second_row(tmp_path: Path):
    """A resubmitted GenBank record must not leave two sequences under one header.

    Real trigger: GenBank revises a record between weekly builds (an accession
    gains 40 bases at the 5' end).  The update batch re-ships that accession, so
    ``sequences`` ends up holding both the old and the new string for the same
    (header, segment).  Every downstream ``SELECT sequence FROM sequences WHERE
    header = ?`` then returns whichever row SQLite reaches first - silently the
    pre-revision sequence half the time.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1", sequence_rows=[("A", "ATGC")]), db_name="seed")

    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", sequence_rows=[("A", "ATGCGGGGGGGG")]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    stored = _rows(updated, "SELECT header, segment, sequence FROM sequences ORDER BY sequence")
    assert stored == [("A", "1", "ATGCGGGGGGGG")], (
        f"expected one row carrying the revised sequence, got {stored}"
    )


@pytest.mark.xfail(
    reason="insertions has the same defect as sequences: it is in upsert_tables "
           "but _create_table_unique_indexes never creates an index over its key "
           "(primary_accession, reference, insertion, segment), so INSERT OR "
           "REPLACE appends.",
    strict=False,
)
def test_update_does_not_accumulate_insertions_for_the_same_accession(tmp_path: Path):
    """Nextalign insertion calls change when the reference alignment changes.

    Real trigger: the reference backbone is re-padded between builds, so
    nextalign reports the same accession's insertion at a different offset.  The
    old call is never retired, so ``insertions`` grows a contradictory second
    row and any per-accession insertion report double-counts.
    """
    seed = _build(
        tmp_path,
        _inputs(tmp_path / "in1", insertion_rows=[["A", "R", "ins:5:A", "1"]]),
        db_name="seed",
    )

    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", insertion_rows=[["A", "R", "ins:7:A", "1"]]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    stored = _rows(updated, "SELECT primary_accession, insertion FROM insertions")
    assert stored == [("A", "ins:7:A")], f"expected only the re-called insertion, got {stored}"


@pytest.mark.xfail(
    reason="Re-running an identical update is not idempotent: sequences and "
           "insertions gain one duplicate row per run because they are upserted "
           "without a unique index.",
    strict=False,
)
def test_rerunning_an_identical_update_leaves_the_db_unchanged(tmp_path: Path):
    """Weekly updates get re-run after a failure downstream of CreateSqliteDB.

    Real trigger: the pipeline dies in ANNOTATE_MUTATIONS or the tree step and
    the operator resubmits the same batch.  An incremental merge must be safe to
    repeat; here each repeat silently inflates ``sequences`` and ``insertions``,
    so row counts drift away from the number of real sequences with no error and
    no audit entry.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    batch = _inputs(tmp_path / "in2")

    first = _build(tmp_path, batch, db_name="u1", update=True, update_db=seed, batch_id="b1")
    second = _build(tmp_path, batch, db_name="u2", update=True, update_db=first, batch_id="b2")

    tables = ["meta_data", "features", "sequence_alignment", "sequences", "insertions", "host_taxa"]
    after_first = {t: _count(first, t) for t in tables}
    after_second = {t: _count(second, t) for t in tables}
    assert after_second == after_first, (
        f"replaying the same update changed row counts: {after_first} -> {after_second}"
    )


@pytest.mark.xfail(
    reason="_create_table_unique_indexes builds idx_features_upsert only when the "
           "features table already has cds_start_OG_seq and cds_end_OG_seq. On a "
           "database whose features table predates those columns it silently "
           "creates no index at all, while merge_table_append_nonredundant still "
           "routes features through INSERT OR REPLACE - so every update appends.",
    strict=False,
)
def test_update_of_legacy_features_schema_does_not_duplicate_rows(tmp_path: Path):
    """An older DB whose features table has only cds_start/cds_end must still upsert.

    Real trigger: a database built before the OG-coordinate columns existed is
    handed to ``--update``.  ``_infer_key_cols`` correctly falls back to
    ``cds_start``/``cds_end``, so the merge believes it has a key, but the guard
    in ``_create_table_unique_indexes`` requires the OG columns and quietly does
    nothing.  Re-supplying a byte-identical feature row then doubles it, and the
    features table doubles again on every subsequent update.
    """
    legacy_row = [["A", "M", "R", "1", "10", "1", "10", "P", "1"]]
    seed = _build(
        tmp_path,
        _inputs(tmp_path / "in1", feature_rows=legacy_row, feature_columns=LEGACY_FEATURE_COLUMNS),
        db_name="seed",
    )
    assert _count(seed, "features") == 1

    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", feature_rows=legacy_row, feature_columns=LEGACY_FEATURE_COLUMNS),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    assert _count(updated, "features") == 1, (
        "re-supplying an identical feature row duplicated it: "
        f"{_rows(updated, 'SELECT accession, cds_start, cds_end, product, segment FROM features')}"
    )


@pytest.mark.xfail(
    reason="SQLite treats NULLs as distinct in a UNIQUE index, so a pre-existing "
           "row whose key column is NULL never conflicts with the incoming row. "
           "INSERT OR REPLACE then appends and the stale row survives forever.",
    strict=False,
)
def test_update_supersedes_a_row_whose_stored_segment_is_null(tmp_path: Path):
    """A NULL in a key column makes the upsert index unable to see the old row.

    Real trigger: a database whose meta_data predates the ``segment`` column.
    ``_ensure_update_columns`` ALTERs the column in, so every existing row gets
    ``segment IS NULL``; ``idx_metadata_upsert`` on (primary_accession, segment)
    is then created happily because NULLs never collide.  The update inserts a
    second row for the same accession and the old, now-unreachable metadata is
    still returned by ``SELECT * FROM meta_data WHERE primary_accession = ?``.
    """
    seed = _build(
        tmp_path,
        _inputs(tmp_path / "in1", segment="4", meta_rows=[["A", "4", "Bos taurus"]]),
        db_name="seed",
    )
    conn = sqlite3.connect(str(seed))
    try:
        conn.execute("UPDATE meta_data SET segment = NULL")
        conn.commit()
    finally:
        conn.close()

    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", segment="4", meta_rows=[["A", "4", "Homo sapiens"]]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    stored = _rows(updated, "SELECT primary_accession, segment, host FROM meta_data")
    assert stored == [("A", "4", "Homo sapiens")], (
        f"the NULL-segment row was not superseded: {stored}"
    )


# --------------------------------------------------------------------------
# 2. INSERT OR REPLACE replaces the WHOLE row
# --------------------------------------------------------------------------


@pytest.mark.xfail(
    reason="INSERT OR REPLACE deletes the conflicting row and inserts the incoming "
           "one wholesale, so a column that is blank in this batch overwrites a "
           "populated value with NULL. There is no per-column 'keep existing when "
           "incoming is empty' step anywhere in merge_table_append_nonredundant.",
    strict=False,
)
def test_update_with_a_blank_cell_does_not_erase_curated_metadata(tmp_path: Path):
    """A thinner batch must not blank fields an earlier batch had filled.

    Real trigger: ``_read_meta_data_tsv`` maps 'NA', 'None' and '' to NaN for
    every column outside NA_IS_A_VALUE_COLUMNS - and a real HCV matrix carries
    274k ``host_validated='NA'`` cells.  When an accession is re-shipped in a
    later batch whose host validation came back 'NA' (or whose GISAID merge was
    not supplied), the previously validated host is silently replaced with NULL.
    """
    columns = ["primary_accession", "segment", "host", "host_validated"]
    seed = _build(
        tmp_path,
        _inputs(tmp_path / "in1", meta_columns=columns, meta_rows=[["A", "1", "Bos taurus", "Bos taurus"]]),
        db_name="seed",
    )
    assert _rows(seed, "SELECT host, host_validated FROM meta_data") == [("Bos taurus", "Bos taurus")]

    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", meta_columns=columns, meta_rows=[["A", "1", "", "NA"]]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    stored = _rows(updated, "SELECT host, host_validated FROM meta_data")
    assert stored == [("Bos taurus", "Bos taurus")], (
        f"a blank incoming cell erased curated metadata: {stored}"
    )


@pytest.mark.xfail(
    reason="_align_df_to_existing_schema fabricates the literal string "
           "'NA- see tree' for any cluster_<N>pct column the incoming frame lacks, "
           "and the meta_data upsert then writes that placeholder over the real "
           "cluster representative stored by the previous build.",
    strict=False,
)
def test_update_without_cluster_tsv_preserves_existing_cluster_assignments(tmp_path: Path):
    """Skipping MMseqs on an update must not destroy the stored clustering.

    Real trigger: the seed build ran MMseqs and stored ``cluster_98pct`` for
    every accession.  The next update runs without ``--cluster_tsv`` (clustering
    disabled, or ``cluster_min_seq_id`` unset so the column would be named
    ``cluster`` instead).  ``_align_df_to_existing_schema`` fills the missing
    column with a placeholder rather than leaving it out, so every re-supplied
    accession loses its cluster representative.
    """
    cluster_tsv = tmp_path / "clusters.tsv"
    cluster_tsv.write_text("A\tA\n", encoding="utf-8")

    seed = _build(
        tmp_path,
        _inputs(tmp_path / "in1"),
        db_name="seed",
        cluster_tsv=cluster_tsv,
        cluster_min_seq_id="0.98",
    )
    assert _rows(seed, "SELECT cluster_98pct FROM meta_data") == [("A",)]

    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2"),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    stored = _rows(updated, "SELECT primary_accession, cluster_98pct FROM meta_data")
    assert stored == [("A", "A")], (
        f"the stored cluster representative was overwritten with a placeholder: {stored}"
    )


# --------------------------------------------------------------------------
# 3. _dedupe_incoming_df: 'first' is file order, and nothing records the loss
# --------------------------------------------------------------------------


@pytest.mark.xfail(
    reason="_dedupe_incoming_df uses drop_duplicates(keep='first'), so which of "
           "two same-key incoming rows survives is decided by their order in the "
           "TSV rather than by content, and the discarded row is only printed to "
           "stdout - nothing is written to update_exclusions.",
    strict=False,
)
def test_incoming_duplicate_keys_do_not_depend_on_row_order(tmp_path: Path):
    """Two rows for one accession must not resolve by whichever was written first.

    Real trigger: the gB matrix carries the same accession twice - one row from
    the GenBank XML and one from a GISAID/external merge, or the same accession
    present in two XML batches.  Whichever line pandas read first wins, so a row
    with an empty host can silently beat a fully populated one, and the losing
    row leaves no trace in ``update_exclusions``.
    """
    results = {}
    for label, meta_rows in (
        ("blank_first", [["A", "1", ""], ["A", "1", "Homo sapiens"]]),
        ("blank_last", [["A", "1", "Homo sapiens"], ["A", "1", ""]]),
    ):
        work = tmp_path / label
        seed = _build(work, _inputs(work / "in1"), db_name="seed")
        updated = _build(
            work,
            _inputs(work / "in2", meta_rows=meta_rows),
            db_name="updated",
            update=True,
            update_db=seed,
            batch_id="batch_update",
        )
        results[label] = _rows(updated, "SELECT host FROM meta_data")

    assert results["blank_first"] == results["blank_last"], (
        f"row order changed which duplicate survived: {results}"
    )


@pytest.mark.xfail(
    reason="update_exclusions is only populated for the non-upsert append path. "
           "Rows dropped by _dedupe_incoming_df, and rows replaced in the upsert "
           "tables, are never recorded, so the audit table cannot account for the "
           "difference between rows submitted and rows stored.",
    strict=False,
)
def test_dropped_incoming_duplicates_are_recorded_in_update_exclusions(tmp_path: Path):
    """The audit table exists to explain missing rows; it must see this drop.

    Real trigger: an operator reconciling 'we submitted 4,000 accessions but
    meta_data only grew by 3,997' queries ``update_exclusions``.  For meta_data
    it is always empty, because the only writer is the non-upsert branch, so the
    three lost rows are invisible outside the process stdout.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", meta_rows=[["A", "1", "Bos taurus"], ["A", "1", "Homo sapiens"]]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    logged = _rows(updated, "SELECT table_name, key, reason FROM update_exclusions WHERE table_name='meta_data'")
    assert logged, "the duplicate incoming meta_data row was dropped without an audit entry"


# --------------------------------------------------------------------------
# 4. Bookkeeping tables: deltas and provenance
# --------------------------------------------------------------------------


@pytest.mark.xfail(
    reason="update_table_deltas records after_count - before_count. A replaced row "
           "yields delta 0 (indistinguishable from 'nothing happened') and a row "
           "duplicated by the index-less upsert yields delta 1 (indistinguishable "
           "from 'one new sequence'), so the delta cannot be read as a count of "
           "new sequences - which is exactly how the run summary prints it.",
    strict=False,
)
def test_update_delta_counts_new_sequences_not_row_growth(tmp_path: Path):
    """'Sequences added to DB' is read as a count of new accessions.

    Real trigger: the operator checks ``update_table_deltas`` after a weekly run.
    Re-supplying an accession whose sequence was revised adds nothing new, yet
    the ``sequences`` delta reads 1 (a duplicate row was appended) while the
    ``meta_data`` delta reads 0 (the row was genuinely replaced).  The two
    tables disagree about the same event, and the misleading one is the one the
    summary line quotes.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1", sequence_rows=[("A", "ATGC")]), db_name="seed")
    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", sequence_rows=[("A", "ATGCGGGG")]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    deltas = dict(
        _rows(updated, "SELECT table_name, delta FROM update_table_deltas WHERE batch_id='batch_update'")
    )
    assert deltas["sequences"] == deltas["meta_data"] == 0, (
        "no new accession was submitted, yet the deltas claim growth: "
        f"meta_data={deltas['meta_data']}, sequences={deltas['sequences']}"
    )


@pytest.mark.xfail(
    reason="project_settings, genes, host_taxa and the m49_* tables take the "
           "non-upsert append-non-redundant path keyed on their first column, so "
           "an existing key is skipped entirely. A changed tool version - or the "
           "'Time of creation' row - can never be refreshed by an update.",
    strict=False,
)
def test_update_refreshes_tool_versions_in_project_settings(tmp_path: Path):
    """The provenance table must describe the run that produced the current data.

    Real trigger: MAFFT is upgraded between builds.  The update runs with the new
    binary, but ``project_settings`` is keyed on Software alone and 'MAFFT' is
    already present, so the row is skipped and logged as ``duplicate_key_in_db``.
    The database then attributes its newest alignments to the old version - and
    the same mechanism freezes 'Time of creation' at the seed build's timestamp.
    """
    seed = _build(
        tmp_path,
        _inputs(tmp_path / "in1", software_rows=[["MAFFT", "v7.4"]]),
        db_name="seed",
    )
    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", software_rows=[["MAFFT", "v7.526"]]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    stored = dict(_rows(updated, "SELECT Software, Version FROM project_settings"))
    assert stored["MAFFT"] == "v7.526", f"tool version was not refreshed by the update: {stored}"


@pytest.mark.xfail(
    reason="create_db only falls back to 'last updated' when db_status is falsy, "
           "and vgtk-init.nf never passes -ds, so args.db_status keeps its "
           "'new db' default. Every --update run therefore stamps "
           "info.creation_type='new db'.",
    strict=False,
)
def test_update_run_records_itself_as_an_update_in_info(tmp_path: Path):
    """``info.creation_type`` is the only human-facing 'is this a fresh build?' flag.

    Real trigger: the shipped nextflow DB process builds its command line without
    ``-ds``, so an incremental run is indistinguishable from a fresh one in the
    published database.  ``update_batches.mode`` says 'update' at the same moment
    ``info`` says 'new db'.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2"),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
        db_status="new db",  # the value vgtk-init.nf leaves in place
    )

    modes = dict(_rows(updated, "SELECT batch_id, mode FROM update_batches"))
    assert modes["batch_update"] == "update"
    creation_types = [row[0] for row in _rows(updated, "SELECT creation_type FROM info")]
    assert creation_types[-1] == "last updated", (
        f"an --update run recorded itself in info as {creation_types[-1]!r}"
    )


# --------------------------------------------------------------------------
# 5. Update assets exported back out of the DB
# --------------------------------------------------------------------------


@pytest.mark.xfail(
    reason="load_reference_rows hard-selects only "
           "[primary_accession, accession_type, segment], so export_ref_list "
           "writes a 3-column headerless TSV. Reloading it through "
           "load_reference_file_table backfills genotype and subtype as empty "
           "strings, which is what CladeAssignment and _load_reference_lookup "
           "read clade labels from.",
    strict=False,
)
def test_exported_ref_list_round_trips_genotype_and_subtype(tmp_path: Path):
    """Regenerating a ref_list from an update DB must not blank the clade labels.

    Real trigger: the operator loses the curated ref_list and regenerates it with
    ``ExportRefListFromUpdateDb.py --db ... -o ref_list.tsv`` before the next
    update.  The exported file has three columns; every reference's genotype and
    subtype come back empty, so the next run assigns blank clades to everything
    while reporting no error.
    """
    db_path = tmp_path / "refs.db"
    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, "
            "segment TEXT, genotype TEXT, subtype TEXT)"
        )
        conn.executemany(
            "INSERT INTO meta_data VALUES (?, ?, ?, ?, ?)",
            [
                ("R1", "master", "1", "1a", "1a.1"),
                ("R2", "reference", "1", "2b", "2b.3"),
                ("Q1", "query", "1", "", ""),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    out = tmp_path / "ref_list_from_db.tsv"
    export_ref_list(str(db_path), str(out))
    reloaded = load_reference_file_table(str(out)).set_index("primary_accession")

    assert reloaded.loc["R1", "genotype"] == "1a", (
        f"genotype was lost by the export/reload round trip: {out.read_text()!r}"
    )
    assert reloaded.loc["R2", "subtype"] == "2b.3"


# --------------------------------------------------------------------------
# 6. Behaviour that is correct - regression tests
# --------------------------------------------------------------------------


def test_update_replaces_meta_data_and_alignment_rows_in_place(tmp_path: Path):
    """The three indexed upsert tables must not grow when a key is re-supplied.

    This is the invariant the whole upsert design rests on, and it is what makes
    the missing indexes on ``sequences``/``insertions`` a defect rather than a
    style choice: with ``idx_metadata_upsert`` and ``idx_seq_alignment_upsert``
    present, re-supplying an accession updates it in place.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1", aln_rows=[["A", "R", "ATGC", "1"]]), db_name="seed")
    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", meta_rows=[["A", "1", "Homo sapiens"]], aln_rows=[["A", "R", "AT--", "1"]]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    assert _rows(updated, "SELECT primary_accession, segment, host FROM meta_data") == [
        ("A", "1", "Homo sapiens")
    ]
    assert _rows(updated, "SELECT primary_accession, alignment_name, alignment FROM sequence_alignment") == [
        ("A", "R", "AT--")
    ]
    assert _count(updated, "features") == 1


def test_update_creates_unique_indexes_only_for_three_of_five_upsert_tables(tmp_path: Path):
    """Pin the exact index set, because it is what decides replace-vs-append.

    ``merge_table_append_nonredundant`` routes five tables through
    ``INSERT OR REPLACE`` but ``_create_table_unique_indexes`` knows about three.
    This test documents that asymmetry so that adding an index for ``sequences``
    or ``insertions`` (the fix) has to come here and flip the assertion.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    updated = _build(
        tmp_path, _inputs(tmp_path / "in2"), db_name="updated", update=True, update_db=seed, batch_id="b"
    )

    indexed = {
        row[0]
        for row in _rows(
            updated,
            "SELECT tbl_name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'",
        )
    }
    assert indexed == {"meta_data", "features", "sequence_alignment"}
    upsert_tables = {"meta_data", "sequence_alignment", "features", "insertions", "sequences"}
    assert upsert_tables - indexed == {"sequences", "insertions"}


def test_update_seed_db_is_copied_not_mutated(tmp_path: Path):
    """``--update_db`` names an input; the run must write to ``-d``/``-o`` instead.

    Real trigger: the published database is passed as the update seed while it is
    still being served.  ``_prepare_update_target_db`` copies it to the output
    path first, so a failed or duplicated run cannot corrupt the live file.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    seed_bytes_before = seed.read_bytes()

    updated = _build(
        tmp_path,
        _inputs(tmp_path / "in2", meta_rows=[["A", "1", "Homo sapiens"]]),
        db_name="updated",
        update=True,
        update_db=seed,
        batch_id="batch_update",
    )

    assert updated != seed
    assert seed.read_bytes() == seed_bytes_before
    assert _rows(seed, "SELECT host FROM meta_data") == [("Bos taurus",)]
    assert _rows(updated, "SELECT host FROM meta_data") == [("Homo sapiens",)]


def test_case_only_column_mismatch_fails_loudly_instead_of_dropping_the_column(tmp_path: Path):
    """A DB column that differs only in case must never silently lose its data.

    Real trigger: the segment/Segment collision that killed an influenza build.
    ``_ensure_update_columns`` matches case-insensitively so it adds nothing, and
    ``_align_df_to_existing_schema`` then drops the incoming ``segment``/``host``
    as 'extra'.  It must raise before writing rather than store a frame whose
    values went nowhere.
    """
    db_path = tmp_path / "mixed_case.db"
    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute('CREATE TABLE meta_data ("primary_accession" TEXT, "Segment" TEXT, "Host" TEXT)')
        conn.execute("INSERT INTO meta_data VALUES ('A', '1', 'Bos taurus')")
        conn.commit()

        merger = CreateSqliteDB.__new__(CreateSqliteDB)
        merger.update = True
        merger.update_db = str(db_path)
        merger.batch_id = "b"

        incoming = pd.DataFrame([{"primary_accession": "A", "segment": "1", "host": "Homo sapiens"}])
        with pytest.raises(ValueError) as excinfo:
            merger.merge_table_append_nonredundant(conn, incoming, "meta_data")

        assert "missing columns required by existing DB schema" in str(excinfo.value)
        # nothing was written
        assert conn.execute("SELECT host FROM meta_data").fetchall() == [("Bos taurus",)]
    finally:
        conn.close()


@pytest.mark.xfail(
    reason="The ValueError raised for a case-only column mismatch blames namespaced "
           "external columns and recommends rebuilding the database, never "
           "mentioning that 'segment' and 'Segment' are the same SQLite identifier. "
           "assert_no_case_insensitive_duplicates already produces the right "
           "wording for the frame-internal version of this collision.",
    strict=False,
)
def test_case_only_column_mismatch_names_the_case_collision(tmp_path: Path):
    """The operator acts on this message; pointing at the wrong cause costs a rebuild.

    Real trigger: the same segment/Segment collision.  The message says
    "A database built before that change cannot be updated in place - rebuild
    it", so the operator discards a good database instead of renaming one column.
    """
    db_path = tmp_path / "mixed_case.db"
    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute('CREATE TABLE meta_data ("primary_accession" TEXT, "Segment" TEXT)')
        conn.execute("INSERT INTO meta_data VALUES ('A', '1')")
        conn.commit()

        merger = CreateSqliteDB.__new__(CreateSqliteDB)
        merger.update = True
        merger.update_db = str(db_path)
        merger.batch_id = "b"

        incoming = pd.DataFrame([{"primary_accession": "A", "segment": "1"}])
        with pytest.raises(ValueError) as excinfo:
            merger.merge_table_append_nonredundant(conn, incoming, "meta_data")

        message = str(excinfo.value).lower()
        assert "case" in message, f"diagnostic does not mention the case collision: {excinfo.value}"
    finally:
        conn.close()


# --------------------------------------------------------------------------
# 7. Evidence from the shipped --update integration database
# --------------------------------------------------------------------------


@pytest.mark.skipif(
    not REAL_UPDATE_DB.exists(), reason=f"RABV update-mode DB not found at {REAL_UPDATE_DB}"
)
def test_real_update_db_has_no_unique_index_on_sequences_or_insertions():
    """Confirm on real output that the two upsert tables are unprotected.

    The synthetic tests above show INSERT OR REPLACE degrading to INSERT; this
    one shows the shipped RABV update database is in exactly that state, so the
    next update against it will start duplicating rows rather than replacing.
    """
    conn = _ro(REAL_UPDATE_DB)
    try:
        indexed = {
            row[0]
            for row in conn.execute(
                "SELECT tbl_name FROM sqlite_master WHERE type='index' AND sql IS NOT NULL"
            )
        }
    finally:
        conn.close()

    assert "meta_data" in indexed and "sequence_alignment" in indexed
    assert "sequences" not in indexed
    assert "insertions" not in indexed


@pytest.mark.xfail(
    reason="The shipped update-mode database records info.creation_type='new db' "
           "for both its create batch and its update batch, because vgtk-init.nf "
           "never passes -ds and create_db only falls back when db_status is falsy.",
    strict=False,
)
@pytest.mark.skipif(
    not REAL_UPDATE_DB.exists(), reason=f"RABV update-mode DB not found at {REAL_UPDATE_DB}"
)
def test_real_update_db_records_its_update_batch_as_an_update():
    """A published database should say how it was last built.

    Real trigger: this database was produced by the repository's own --update
    integration run.  ``update_batches`` holds one 'create' and one 'update'
    batch, while ``info`` holds two rows that both say 'new db'.
    """
    conn = _ro(REAL_UPDATE_DB)
    try:
        modes = [row[0] for row in conn.execute("SELECT mode FROM update_batches ORDER BY started_at")]
        creation_types = [row[0] for row in conn.execute("SELECT creation_type FROM info ORDER BY date")]
    finally:
        conn.close()

    assert "update" in modes, "this database was not built by an update run"
    assert creation_types[-1] == "last updated", (
        f"an update run stamped info.creation_type={creation_types[-1]!r}"
    )


@pytest.mark.xfail(
    reason="project_settings is keyed on its first column via the non-upsert "
           "append path, so the update batch's rows were all rejected as "
           "duplicate_key_in_db and 'Time of creation' still reports the seed "
           "build's timestamp.",
    strict=False,
)
@pytest.mark.skipif(
    not REAL_UPDATE_DB.exists(), reason=f"RABV update-mode DB not found at {REAL_UPDATE_DB}"
)
def test_real_update_db_project_settings_reflect_the_latest_batch():
    """Provenance frozen at the seed build misattributes every later sequence.

    Real trigger: the update batch shipped its own software_info.tsv, including a
    fresh 'Time of creation'.  Every row was skipped and logged into
    ``update_exclusions``; the database still advertises the seed run's
    timestamp and tool versions for data added afterwards.
    """
    conn = _ro(REAL_UPDATE_DB)
    try:
        skipped = conn.execute(
            "SELECT COUNT(*) FROM update_exclusions WHERE table_name='project_settings'"
        ).fetchone()[0]
        creation = conn.execute(
            "SELECT Version FROM project_settings WHERE Software LIKE 'Time of creation%'"
        ).fetchone()
        last_batch = conn.execute(
            "SELECT started_at FROM update_batches WHERE mode='update' ORDER BY started_at DESC LIMIT 1"
        ).fetchone()
    finally:
        conn.close()

    if creation is None or last_batch is None:
        pytest.skip("database does not carry a 'Time of creation' row or an update batch")

    from datetime import datetime, timezone

    stored = datetime.fromtimestamp(float(creation[0]), tz=timezone.utc)
    batch_started = datetime.strptime(last_batch[0], "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)

    assert skipped == 0, f"{skipped} project_settings rows were dropped as duplicate keys"
    assert stored >= batch_started, (
        f"'Time of creation' ({stored}) predates the update batch ({batch_started})"
    )
