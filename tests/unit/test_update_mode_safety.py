"""Safety properties of ``--update``: idempotence, non-destruction, honesty.

Update mode runs weekly against a database holding years of curated data, so its
failure mode is not a crash but a quiet corruption of rows nobody re-derives.
Three mechanisms carry the risk, and each has a section below:

1. ``INSERT OR REPLACE`` only replaces when a UNIQUE index makes the row
   conflict.  ``sequences`` and ``insertions`` are written that way and had no
   index at all, so replaying an identical update appended one row per record.
2. ``INSERT OR REPLACE`` replaces the WHOLE row, so any column that is absent or
   blank in a batch erases what an earlier run stored.  A run without
   ``--cluster_tsv`` wrote the literal string ``NA- see tree`` over the cluster
   representative MMseqs had computed - 100 such rows in the shipped
   ``test_out/update_test/rabv-jul0425-update-test.db``.
3. ``info.creation_type`` is the only human-facing record of how a database was
   built, and every ``--update`` run stamped it "new db".
"""

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from CreateSqliteDB import CreateSqliteDB


REPO_ROOT = Path(__file__).resolve().parents[2]
#: The RABV database produced by the shipped --update integration run. Read-only.
REAL_UPDATE_DB = REPO_ROOT / "test_out" / "update_test" / "rabv-jul0425-update-test.db"

UPSERT_TABLES = {"meta_data", "sequence_alignment", "features", "insertions", "sequences"}

FEATURE_COLUMNS = [
    "accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end",
    "cds_start", "cds_end", "cds_start_OG_seq", "cds_end_OG_seq", "product", "segment",
]


# --------------------------------------------------------------------------
# helpers - one complete set of CreateSqliteDB inputs under tmp_path
# --------------------------------------------------------------------------

def _tsv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _csv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def _inputs(directory: Path, *, meta_rows=None, meta_columns=None, feature_rows=None,
            aln_rows=None, insertion_rows=None, sequence_rows=None, segment="1"):
    directory.mkdir(parents=True, exist_ok=True)
    files = {key: directory / name for key, name in {
        "meta": "meta_data.tsv", "features": "features.tsv", "aln": "sequence_alignment.tsv",
        "gene": "gene_info.tsv", "country": "m49_country.csv", "interm": "m49_intermediate.csv",
        "region": "m49_region.csv", "sub_region": "m49_sub_region.csv",
        "software": "software_info.tsv", "insertions": "insertions.tsv",
        "host_taxa": "host_taxa.tsv", "fasta": "sequences.fa",
    }.items()}

    _tsv(files["meta"],
         meta_rows if meta_rows is not None else [["A", segment, "Bos taurus"]],
         meta_columns or ["primary_accession", "segment", "host"])
    _tsv(files["features"],
         feature_rows if feature_rows is not None
         else [["A", "M", "R", "1", "10", "1", "10", "1", "10", "P", segment]],
         FEATURE_COLUMNS)
    _tsv(files["aln"], aln_rows if aln_rows is not None else [["A", "R", "ATGC", segment]],
         ["primary_accession", "alignment_name", "alignment", "segment"])
    _tsv(files["gene"], [["geneA", "Gene A"]], ["name", "description"])
    _csv(files["country"], [["001", "World"]], ["m49_code", "name"])
    _csv(files["interm"], [["X", "Inter"]], ["code", "name"])
    _csv(files["region"], [["Y", "Region"]], ["code", "name"])
    _csv(files["sub_region"], [["Z", "SubRegion"]], ["code", "name"])
    _tsv(files["software"], [["MAFFT", "v7.4"]], ["Software", "Version"])
    _tsv(files["insertions"],
         insertion_rows if insertion_rows is not None else [["A", "R", "ins:5:A", segment]],
         ["primary_accession", "reference", "insertion", "segment"])
    _tsv(files["host_taxa"], [["A", "host1"]], ["primary_accession", "host"])
    rows = sequence_rows if sequence_rows is not None else [("A", "ATGC")]
    files["fasta"].write_text("".join(f">{h}\n{s}\n" for h, s in rows), encoding="utf-8")
    return files


def _build(base_dir: Path, files: dict, *, db_name: str, update: bool = False,
           update_db: Path | None = None, batch_id: str = "batch_test",
           db_status=None, cluster_tsv: Path | None = None,
           cluster_min_seq_id: str | None = None) -> Path:
    creator = CreateSqliteDB(
        meta_data=str(files["meta"]), features=str(files["features"]),
        pad_aln=str(files["aln"]), gene_info=str(files["gene"]),
        m49_countries=str(files["country"]), m49_interm_region=str(files["interm"]),
        m49_regions=str(files["region"]), m49_sub_regions=str(files["sub_region"]),
        proj_settings=str(files["software"]), fasta_sequence_file=str(files["fasta"]),
        insertions=str(files["insertions"]), host_taxa_file=str(files["host_taxa"]),
        base_dir=str(base_dir), output_dir="SqliteDB", db_name=db_name, db_status=db_status,
        cluster_tsv=str(cluster_tsv) if cluster_tsv else None,
        cluster_min_seq_id=cluster_min_seq_id, update=update,
        update_db=str(update_db) if update_db else None, batch_id=batch_id,
    )
    creator.create_db()
    return base_dir / "SqliteDB" / f"{db_name}.db"


def _rows(db_path: Path, sql: str, params=()):
    conn = sqlite3.connect(str(db_path))
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


def _counts(db_path: Path, tables=UPSERT_TABLES):
    return {t: _rows(db_path, f"SELECT COUNT(*) FROM {t}")[0][0] for t in sorted(tables)}


def _indexed_tables(db_path: Path):
    return {
        row[0] for row in _rows(
            db_path, "SELECT tbl_name FROM sqlite_master WHERE type='index' AND sql IS NOT NULL"
        )
    }


def _index_columns(db_path: Path, index_name: str):
    conn = sqlite3.connect(str(db_path))
    try:
        return [row[2] for row in conn.execute(f'PRAGMA index_info("{index_name}")')]
    finally:
        conn.close()


# ==========================================================================
# 1. the two missing unique indexes
# ==========================================================================

def test_sequences_and_insertions_get_a_unique_index_over_their_upsert_key(tmp_path: Path):
    """The key has to be the one ``_infer_key_cols`` resolved, or ON CONFLICT
    has no target and the write silently degrades to an append."""
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")

    assert UPSERT_TABLES <= _indexed_tables(seed)
    assert _index_columns(seed, "idx_sequences_upsert") == ["header", "segment"]
    assert _index_columns(seed, "idx_insertions_upsert") == [
        "primary_accession", "reference", "insertion", "segment"
    ]


def test_replaying_an_identical_update_leaves_every_row_count_unchanged(tmp_path: Path):
    """Idempotence is the property update mode is supposed to have.

    Real trigger: a weekly run is re-submitted after a cluster failure, or the
    same batch is applied to two mirrors.  Without an index over their keys,
    ``sequences`` and ``insertions`` grew by one row per record every time.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    first = _build(tmp_path, _inputs(tmp_path / "in2"), db_name="first",
                   update=True, update_db=seed, batch_id="b1")
    before = _counts(first)

    second = _build(tmp_path, _inputs(tmp_path / "in3"), db_name="second",
                    update=True, update_db=first, batch_id="b2")

    assert _counts(second) == before, "replaying an identical update changed the row counts"


def test_a_resupplied_sequence_is_replaced_not_appended(tmp_path: Path):
    """The other half of the same defect: the new sequence text must win.

    Appending left the old row in place, and every ``SELECT sequence FROM
    sequences WHERE header = ?`` after that returned whichever row SQLite reached
    first - i.e. the superseded one.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1", sequence_rows=[("A", "ATGC")]),
                  db_name="seed")
    updated = _build(tmp_path, _inputs(tmp_path / "in2", sequence_rows=[("A", "TTTT")]),
                     db_name="updated", update=True, update_db=seed, batch_id="b")

    assert _rows(updated, "SELECT header, sequence FROM sequences") == [("A", "TTTT")]


def test_insertions_do_not_accumulate_for_the_same_accession(tmp_path: Path):
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    updated = _build(tmp_path, _inputs(tmp_path / "in2"), db_name="updated",
                     update=True, update_db=seed, batch_id="b")

    assert _rows(updated, "SELECT primary_accession, reference, insertion FROM insertions") == [
        ("A", "R", "ins:5:A")
    ]


def test_a_legacy_features_table_still_gets_an_index_over_the_key_it_has(tmp_path: Path):
    """``cds_start_OG_seq``/``cds_end_OG_seq`` are recent columns.

    ``_infer_key_cols`` already falls back to ``cds_start``/``cds_end`` for a
    database that predates them; the index has to follow the same fallback or
    that database gets no index at all while still being routed through the
    upsert.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    conn = sqlite3.connect(str(seed))
    try:
        conn.execute("DROP INDEX idx_features_upsert")
        conn.execute("ALTER TABLE features DROP COLUMN cds_start_OG_seq")
        conn.execute("ALTER TABLE features DROP COLUMN cds_end_OG_seq")
        conn.commit()
    finally:
        conn.close()

    legacy_features = [["A", "M", "R", "1", "10", "1", "10", "P", "1"]]
    files = _inputs(tmp_path / "in2")
    _tsv(files["features"], legacy_features,
         ["accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end",
          "cds_start", "cds_end", "product", "segment"])

    updated = _build(tmp_path, files, db_name="updated", update=True, update_db=seed, batch_id="b")

    assert _index_columns(updated, "idx_features_upsert") == [
        "accession", "cds_start", "cds_end", "product", "segment"
    ]
    assert _rows(updated, "SELECT COUNT(*) FROM features") == [(1,)]


def test_a_stored_null_key_column_is_normalised_so_the_index_can_see_it(tmp_path: Path, capsys):
    """SQLite treats NULLs as DISTINCT in a UNIQUE index.

    Real trigger: a database whose ``sequences`` table predates the ``segment``
    column.  ``_ensure_update_columns`` ALTERs it in, so every existing row has
    ``segment IS NULL``; the index is then created happily because NULLs never
    collide, and the update appends beside rows it can never reach again.
    Incoming keys are normalised to ``''``, so stored NULLs are too.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    conn = sqlite3.connect(str(seed))
    try:
        conn.execute("DROP INDEX idx_sequences_upsert")
        conn.execute("UPDATE sequences SET segment = NULL")
        conn.commit()
    finally:
        conn.close()

    updated = _build(tmp_path, _inputs(tmp_path / "in2", sequence_rows=[("A", "TTTT")]),
                     db_name="updated", update=True, update_db=seed, batch_id="b")

    message = capsys.readouterr().out
    assert "NULL value(s) in the upsert key" in message
    assert _rows(updated, "SELECT segment FROM sequences WHERE segment IS NULL") == []
    assert _rows(updated, "SELECT header, segment, sequence FROM sequences") == [("A", "1", "TTTT")]


def test_pre_existing_duplicate_rows_fail_with_a_message_naming_them(tmp_path: Path):
    """The first update against an old database is the first time the constraint
    exists, so the CREATE is where a pre-existing violation surfaces.

    "UNIQUE constraint failed: sequences.header, sequences.segment" on its own
    says nothing about which rows, or what to do; the operator needs both.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    conn = sqlite3.connect(str(seed))
    try:
        conn.execute("DROP INDEX idx_sequences_upsert")
        conn.execute("INSERT INTO sequences (header, sequence, segment) VALUES ('A', 'GGGG', '1')")
        conn.commit()
    finally:
        conn.close()

    with pytest.raises(ValueError) as excinfo:
        _build(tmp_path, _inputs(tmp_path / "in2"), db_name="updated",
               update=True, update_db=seed, batch_id="b")

    message = str(excinfo.value)
    assert "idx_sequences_upsert" in message
    assert "'A'" in message, f"the offending key is not named: {message}"
    assert "de-duplicate" in message.lower()


@pytest.mark.xfail(
    reason="Normalising a stored NULL key to '' lets the index see the row, but "
           "('A', '') is still a different key from ('A', '4'), so a segmented "
           "database whose rows predate the segment column keeps one unreachable "
           "row per accession. Superseding it would mean matching on a subset of "
           "the key, which cannot be told apart from a genuinely different segment.",
    strict=False,
)
def test_a_null_segment_row_is_superseded_by_the_batch_that_names_its_segment(tmp_path: Path):
    seed = _build(tmp_path, _inputs(tmp_path / "in1", segment="4"), db_name="seed")
    conn = sqlite3.connect(str(seed))
    try:
        conn.execute("DROP INDEX idx_metadata_upsert")
        conn.execute("UPDATE meta_data SET segment = NULL")
        conn.commit()
    finally:
        conn.close()

    updated = _build(tmp_path, _inputs(tmp_path / "in2", segment="4",
                                       meta_rows=[["A", "4", "Homo sapiens"]]),
                     db_name="updated", update=True, update_db=seed, batch_id="b")

    assert _rows(updated, "SELECT primary_accession, segment, host FROM meta_data") == [
        ("A", "4", "Homo sapiens")
    ]


@pytest.mark.skipif(not REAL_UPDATE_DB.exists(), reason=f"not found: {REAL_UPDATE_DB}")
def test_the_shipped_update_db_can_take_the_new_indexes(tmp_path: Path):
    """Applying a new constraint to published data is only safe if the data
    already satisfies it - so check the real database rather than assuming.

    Read-only: the shipped file is copied to tmp_path first.
    """
    import shutil

    copy = tmp_path / "rabv-update.copy.db"
    shutil.copyfile(REAL_UPDATE_DB, copy)

    creator = CreateSqliteDB.__new__(CreateSqliteDB)
    creator.update = True
    conn = sqlite3.connect(str(copy))
    try:
        created = {}
        for table in sorted(UPSERT_TABLES):
            columns = CreateSqliteDB._table_columns(conn, table)
            key_cols = creator._infer_key_cols(table, pd.DataFrame(columns=columns))
            created[table] = creator._create_table_unique_indexes(conn, table, key_cols)
        conn.commit()
    finally:
        conn.close()

    assert created["sequences"] == "idx_sequences_upsert"
    assert created["insertions"] == "idx_insertions_upsert"
    assert all(created.values()), f"no index could be built for: {sorted(k for k, v in created.items() if not v)}"


# ==========================================================================
# 2. a column the batch does not carry must keep its stored value
# ==========================================================================

def test_an_update_without_cluster_tsv_preserves_the_stored_cluster_rep(tmp_path: Path):
    """Skipping MMseqs on an update must not destroy the stored clustering.

    Real trigger: the seed build ran MMseqs and stored ``cluster_98pct`` for
    every accession; the next update runs without ``--cluster_tsv``.  The absent
    column was fabricated as the literal string "NA- see tree" and written over
    the real representative - 100 rows in the shipped update-mode database carry
    that placeholder today.
    """
    cluster_tsv = tmp_path / "clusters.tsv"
    cluster_tsv.write_text("REP_A\tA\n", encoding="utf-8")
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed",
                  cluster_tsv=cluster_tsv, cluster_min_seq_id="0.98")
    assert _rows(seed, "SELECT cluster_98pct FROM meta_data") == [("REP_A",)]

    updated = _build(tmp_path, _inputs(tmp_path / "in2"), db_name="updated",
                     update=True, update_db=seed, batch_id="b")

    assert _rows(updated, "SELECT cluster_98pct FROM meta_data") == [("REP_A",)]


def test_a_new_accession_still_gets_the_cluster_placeholder(tmp_path: Path):
    """The placeholder is what ValidateDbTree recognises as "not clustered", so
    it still has to be written for rows that have nothing to preserve."""
    cluster_tsv = tmp_path / "clusters.tsv"
    cluster_tsv.write_text("REP_A\tA\n", encoding="utf-8")
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed",
                  cluster_tsv=cluster_tsv, cluster_min_seq_id="0.98")

    files = _inputs(tmp_path / "in2",
                    meta_rows=[["A", "1", "Bos taurus"], ["NEW1", "1", "Bos taurus"]],
                    aln_rows=[["A", "R", "ATGC", "1"], ["NEW1", "R", "ATGC", "1"]],
                    sequence_rows=[("A", "ATGC"), ("NEW1", "ATGC")])
    updated = _build(tmp_path, files, db_name="updated", update=True, update_db=seed, batch_id="b")

    stored = dict(_rows(updated, "SELECT primary_accession, cluster_98pct FROM meta_data"))
    assert stored == {"A": "REP_A", "NEW1": "NA- see tree"}


def test_a_column_absent_from_the_batch_keeps_its_stored_value(tmp_path: Path):
    """The general form of the cluster hazard.

    Real trigger: a column added by an optional merge step (GISAID metadata, host
    validation) that a later, thinner batch does not run.  INSERT OR REPLACE
    rewrites the whole row, so the curated value went to NULL.
    """
    seed = _build(tmp_path,
                  _inputs(tmp_path / "in1",
                          meta_columns=["primary_accession", "segment", "host", "host_validated"],
                          meta_rows=[["A", "1", "Bos taurus", "Bos taurus"]]),
                  db_name="seed")

    updated = _build(tmp_path,
                     _inputs(tmp_path / "in2",
                             meta_columns=["primary_accession", "segment", "host"],
                             meta_rows=[["A", "1", "Homo sapiens"]]),
                     db_name="updated", update=True, update_db=seed, batch_id="b")

    assert _rows(updated, "SELECT host, host_validated FROM meta_data") == [
        ("Homo sapiens", "Bos taurus")
    ]


def test_a_blank_incoming_cell_does_not_erase_curated_metadata(tmp_path: Path):
    """A blank cell means "this batch says nothing", not "delete that".

    Real trigger: ``_read_meta_data_tsv`` maps 'NA', 'None' and '' to NaN for
    every column outside NA_IS_A_VALUE_COLUMNS, and a real HCV matrix carries
    274k ``host_validated='NA'`` cells.  Re-shipping such an accession replaced
    the previously validated host with NULL.
    """
    columns = ["primary_accession", "segment", "host", "host_validated"]
    seed = _build(tmp_path,
                  _inputs(tmp_path / "in1", meta_columns=columns,
                          meta_rows=[["A", "1", "Bos taurus", "Bos taurus"]]),
                  db_name="seed")

    updated = _build(tmp_path,
                     _inputs(tmp_path / "in2", meta_columns=columns,
                             meta_rows=[["A", "1", "", "NA"]]),
                     db_name="updated", update=True, update_db=seed, batch_id="b")

    assert _rows(updated, "SELECT host, host_validated FROM meta_data") == [
        ("Bos taurus", "Bos taurus")
    ]


def test_a_blank_alignment_does_not_erase_the_stored_one(tmp_path: Path):
    """Same rule in the table where the value is the actual science."""
    seed = _build(tmp_path, _inputs(tmp_path / "in1", aln_rows=[["A", "R", "ATGC", "1"]]),
                  db_name="seed")
    updated = _build(tmp_path, _inputs(tmp_path / "in2", aln_rows=[["A", "R", "", "1"]]),
                     db_name="updated", update=True, update_db=seed, batch_id="b")

    assert _rows(updated, "SELECT alignment FROM sequence_alignment") == [("ATGC",)]


def test_exclusion_status_can_still_be_cleared_by_a_later_batch(tmp_path: Path):
    """The deliberate exemption from "a blank never overwrites".

    Real trigger: a sequence excluded by one batch (a bad alignment against a
    reference that has since been replaced) passes QC in the next.  If blanks
    could never overwrite, nothing could ever un-exclude it and the row would
    stay hidden from every downstream query for good.
    """
    columns = ["primary_accession", "segment", "host", "exclusion_status", "exclusion_criteria"]
    seed = _build(tmp_path,
                  _inputs(tmp_path / "in1", meta_columns=columns,
                          meta_rows=[["A", "1", "Bos taurus", "1", "alignment_filtering"]]),
                  db_name="seed")
    assert _rows(seed, "SELECT exclusion_status FROM meta_data") == [("1",)]

    updated = _build(tmp_path,
                     _inputs(tmp_path / "in2", meta_columns=columns,
                             meta_rows=[["A", "1", "Bos taurus", "", ""]]),
                     db_name="updated", update=True, update_db=seed, batch_id="b")

    status, criteria = _rows(updated, "SELECT exclusion_status, exclusion_criteria FROM meta_data")[0]
    assert (status or "") == "", f"an excluded row could not be un-excluded: {status!r}"
    assert (criteria or "") == ""


def test_a_missing_column_in_a_non_upsert_table_still_fails_loudly(tmp_path: Path):
    """The tolerant path is for the upsert tables only.

    ``m49_country`` and friends are appended wholesale, so a batch that has lost
    a column is a broken input rather than a partial update, and must not be
    written at all.
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    files = _inputs(tmp_path / "in2")
    _csv(files["country"], [["001", "World"]], ["m49_code", "display_name"])

    with pytest.raises(ValueError, match="missing columns required by existing DB schema"):
        _build(tmp_path, files, db_name="updated", update=True, update_db=seed, batch_id="b")


# ==========================================================================
# 3. db_status: recording, and using, how the database was built
# ==========================================================================

def test_a_fresh_build_records_new_db(tmp_path: Path):
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    assert [row[0] for row in _rows(seed, "SELECT creation_type FROM info")] == ["new db"]


def test_an_update_records_last_updated_even_when_told_new_db(tmp_path: Path, capsys):
    """The run mode is the authority, and the contradiction is reported.

    Real trigger: vgtk-init.nf never passed ``-ds``, so ``args.db_status`` kept
    its "new db" default and ``if db_status`` was always true - every update run
    stamped itself a fresh build.  The shipped update-mode database has two info
    rows and both say "new db".
    """
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed", db_status="new db")
    updated = _build(tmp_path, _inputs(tmp_path / "in2"), db_name="updated", update=True,
                     update_db=seed, batch_id="b", db_status="new db")

    creation_types = [row[0] for row in _rows(updated, "SELECT creation_type FROM info")]
    assert creation_types == ["new db", "last updated"], creation_types
    assert "--db_status says 'new db'" in capsys.readouterr().out


def test_the_info_table_keeps_its_two_column_shape(tmp_path: Path):
    """Backwards compatibility: consumers read (creation_type, date)."""
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    conn = sqlite3.connect(str(seed))
    try:
        assert [row[1] for row in conn.execute("PRAGMA table_info(info)")] == ["creation_type", "date"]
    finally:
        conn.close()


def test_history_accumulates_so_a_database_can_say_how_it_was_built(tmp_path: Path):
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    first = _build(tmp_path, _inputs(tmp_path / "in2"), db_name="first", update=True,
                   update_db=seed, batch_id="b1")
    second = _build(tmp_path, _inputs(tmp_path / "in3"), db_name="second", update=True,
                    update_db=first, batch_id="b2")

    assert [row[0] for row in _rows(second, "SELECT creation_type FROM info")] == [
        "new db", "last updated", "last updated"
    ]


def test_updating_something_that_is_not_a_pipeline_database_is_refused(tmp_path: Path):
    """``--update_db`` naming the wrong file must not silently produce a new one.

    Real trigger: a mistyped path that happens to exist, or a database from a
    different project.  Without meta_data there is nothing to update
    incrementally, and every "duplicate key" guard downstream is inert.
    """
    stray = tmp_path / "not_a_pipeline.db"
    conn = sqlite3.connect(str(stray))
    try:
        conn.execute("CREATE TABLE something_else (x TEXT)")
        conn.commit()
    finally:
        conn.close()

    with pytest.raises(ValueError, match="no meta_data table"):
        _build(tmp_path, _inputs(tmp_path / "in1"), db_name="updated", update=True,
               update_db=stray, batch_id="b")


def test_updating_a_database_with_no_recorded_history_warns_but_proceeds(tmp_path: Path, capsys):
    """A legacy database is a legitimate state, not an error - but say so."""
    seed = _build(tmp_path, _inputs(tmp_path / "in1"), db_name="seed")
    conn = sqlite3.connect(str(seed))
    try:
        conn.execute("DELETE FROM info")
        conn.commit()
    finally:
        conn.close()

    updated = _build(tmp_path, _inputs(tmp_path / "in2"), db_name="updated", update=True,
                     update_db=seed, batch_id="b")

    assert "no recorded build history" in capsys.readouterr().out
    assert [row[0] for row in _rows(updated, "SELECT creation_type FROM info")] == ["last updated"]


@pytest.mark.skipif(not REAL_UPDATE_DB.exists(), reason=f"not found: {REAL_UPDATE_DB}")
def test_the_shipped_update_db_shows_the_history_that_was_lost():
    """Read-only evidence for the defect this section fixes: the database was
    built by one create batch and one update batch, and both info rows say
    "new db".  A rebuild with the fix in place produces "new db", "last updated".
    """
    conn = sqlite3.connect(f"file:{REAL_UPDATE_DB}?mode=ro", uri=True)
    try:
        modes = [row[0] for row in conn.execute("SELECT mode FROM update_batches ORDER BY started_at")]
        creation_types = [row[0] for row in conn.execute("SELECT creation_type FROM info ORDER BY date")]
    finally:
        conn.close()

    assert "update" in modes, "this database was not built by an update run"
    assert len(creation_types) == len(modes), (
        "info and update_batches disagree about how many runs touched this database"
    )
    if creation_types[-1] == "last updated":
        pytest.skip("shipped DB has been rebuilt with the fix")
    assert creation_types == ["new db", "new db"], creation_types
