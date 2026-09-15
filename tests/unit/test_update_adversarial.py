"""Update mode, under the conditions that actually corrupt a database.

An update is the dangerous operation: the rebuild path starts from nothing and
is self-consistent by construction, whereas an update writes new rows beside
years of old ones and every disagreement between the two becomes a silent,
permanent inconsistency in a database nobody rebuilds.

The scenarios here are the ones that a real change to the pipeline produces:

  * a coordinate correction lands, so the same (accession, product) now has
    DIFFERENT coordinates from the rows already stored. Those columns are part
    of the features upsert key, so the corrected row does not replace the old
    one - it joins it, and the gene now exists twice with two answers;
  * the stored alignment gets wider between runs, so historic rows are in a
    narrower column space than the incoming ones;
  * the master for a segment changes between runs;
  * a segment appears that the database has never seen;
  * the same batch is applied twice.

Where a scenario is currently unsafe the test records precisely what happens,
so the cost of the migration decision is visible rather than folklore.
"""

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from CreateSqliteDB import CreateSqliteDB


FEATURE_COLUMNS = [
    "accession", "master_ref_accession", "reference_accession", "aln_start",
    "aln_end", "cds_start", "cds_end", "cds_start_OG_seq", "cds_end_OG_seq",
    "product", "segment",
]


def _tsv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _csv(path: Path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def make_inputs(tmp_path: Path, suffix: str, feature_rows=None, meta_rows=None,
                alignment="ATGC", accession="A", segment="1"):
    """The full input set CreateSqliteDB requires, with the parts we vary."""
    paths = {name: tmp_path / f"{name}_{suffix}" for name in (
        "meta.tsv", "features.tsv", "aln.tsv", "gene.tsv", "m49c.csv", "m49i.csv",
        "m49r.csv", "m49s.csv", "proj.tsv", "insertions.tsv", "host.tsv", "seqs.fa",
    )}
    if feature_rows is None:
        feature_rows = [[accession, "M", "R", "1", "10", "1", "10", "1", "10", "P", segment]]
    if meta_rows is None:
        meta_rows = [[accession, "", segment]]

    _tsv(paths["meta.tsv"], meta_rows, ["primary_accession", "exclusion", "segment"])
    _tsv(paths["features.tsv"], feature_rows, FEATURE_COLUMNS)
    _tsv(paths["aln.tsv"], [[accession, "R", alignment, segment]],
         ["primary_accession", "alignment_name", "alignment", "segment"])
    _tsv(paths["gene.tsv"], [["geneA", "Gene A"]], ["name", "description"])
    _csv(paths["m49c.csv"], [["001", "World"]], ["m49_code", "name"])
    _csv(paths["m49i.csv"], [["X", "Inter"]], ["code", "name"])
    _csv(paths["m49r.csv"], [["Y", "Region"]], ["code", "name"])
    _csv(paths["m49s.csv"], [["Z", "SubRegion"]], ["code", "name"])
    _tsv(paths["proj.tsv"], [["Python", "3.11"]], ["Software", "Version"])
    _tsv(paths["insertions.tsv"], [[accession, "R", "ins:5:A", segment]],
         ["primary_accession", "reference", "insertion", "segment"])
    _tsv(paths["host.tsv"], [[accession, "host1"]], ["primary_accession", "host"])
    paths["seqs.fa"].write_text(f">{accession}\n{alignment.replace('-', '')}\n", encoding="utf-8")
    return paths


def build_db(tmp_path: Path, paths, db_name, update=False, update_db=None,
             batch_id="batch_test"):
    db = CreateSqliteDB(
        meta_data=str(paths["meta.tsv"]),
        features=str(paths["features.tsv"]),
        pad_aln=str(paths["aln.tsv"]),
        gene_info=str(paths["gene.tsv"]),
        m49_countries=str(paths["m49c.csv"]),
        m49_interm_region=str(paths["m49i.csv"]),
        m49_regions=str(paths["m49r.csv"]),
        m49_sub_regions=str(paths["m49s.csv"]),
        proj_settings=str(paths["proj.tsv"]),
        fasta_sequence_file=str(paths["seqs.fa"]),
        insertions=str(paths["insertions.tsv"]),
        host_taxa_file=str(paths["host.tsv"]),
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name=db_name,
        db_status="last updated" if update else "new db",
        update=update,
        update_db=str(update_db) if update_db else None,
        batch_id=batch_id,
    )
    db.create_db()
    return tmp_path / "SqliteDB" / f"{db_name}.db"


def features_of(db_path, accession="A", product="P"):
    """Coordinates as TEXT, which is how the shipped schema stores them.

    Normalised to str because the columns are TEXT in a fresh build but arrive
    as INTEGER through some update paths, and a test that cared about the
    difference would be testing sqlite's type affinity rather than the pipeline.
    """
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        rows = conn.execute(
            "SELECT cds_start, cds_end, cds_start_OG_seq, cds_end_OG_seq "
            "FROM features WHERE accession=? AND product=?", (accession, product)
        ).fetchall()
    finally:
        conn.close()
    return [tuple("" if value is None else str(value) for value in row) for row in rows]


# ---------------------------------------------------------------------------
# 1. a coordinate correction arriving as an update
# ---------------------------------------------------------------------------


def test_corrected_coordinates_do_not_replace_the_rows_they_correct(tmp_path):
    """The migration hazard, demonstrated rather than described.

    features is upserted on a key that INCLUDES cds_start_OG_seq and
    cds_end_OG_seq (CreateSqliteDB._infer_key_cols). So when a coordinate fix
    changes those values, the corrected row does not match the stored one and is
    inserted alongside it. The gene then has two rows with two different
    answers, and nothing downstream knows which to believe: AnnotateMutations
    folds features into a per-product map and choose_feature_entry breaks the
    tie by preferring the SHORTEST span.
    """
    first = make_inputs(tmp_path, "v1")
    db = build_db(tmp_path, first, "adv_db")
    assert features_of(db) == [("1", "10", "1", "10")]

    # The same accession and product, now with corrected coordinates.
    corrected = make_inputs(
        tmp_path, "v2",
        feature_rows=[["A", "M", "R", "1", "10", "4", "13", "4", "13", "P", "1"]],
    )
    updated = build_db(tmp_path, corrected, "adv_db_2", update=True, update_db=db,
                       batch_id="batch_correction")

    rows = features_of(updated)
    assert len(rows) == 2, (
        "a coordinate correction currently ADDS a row rather than replacing it; "
        "if this becomes 1 the upsert key has been changed and the migration "
        "question is settled"
    )
    assert {(r[2], r[3]) for r in rows} == {("1", "10"), ("4", "13")}


def test_identical_coordinates_do_replace_cleanly(tmp_path):
    """The control: when nothing about the key changes, the row is replaced."""
    first = make_inputs(tmp_path, "s1")
    db = build_db(tmp_path, first, "same_db")

    again = make_inputs(tmp_path, "s2")
    updated = build_db(tmp_path, again, "same_db_2", update=True, update_db=db,
                       batch_id="batch_same")

    assert len(features_of(updated)) == 1


def test_applying_the_same_batch_twice_is_idempotent(tmp_path):
    """Re-running an update must not multiply rows."""
    first = make_inputs(tmp_path, "i1")
    db = build_db(tmp_path, first, "idem_db")

    once = build_db(tmp_path, make_inputs(tmp_path, "i2"), "idem_2",
                    update=True, update_db=db, batch_id="batch_x")
    twice = build_db(tmp_path, make_inputs(tmp_path, "i3"), "idem_3",
                     update=True, update_db=once, batch_id="batch_x")

    assert len(features_of(twice)) == 1


# ---------------------------------------------------------------------------
# 2. the shape of the database changes underneath the update
# ---------------------------------------------------------------------------


def test_update_into_a_db_without_the_og_columns(tmp_path):
    """An older database predates cds_*_OG_seq; the columns must be added."""
    first = make_inputs(tmp_path, "o1")
    db = build_db(tmp_path, first, "old_db")

    conn = sqlite3.connect(str(db))
    try:
        # The upsert index names the OG columns, so it has to go first - which is
        # itself worth knowing: a migration that drops or rewrites these columns
        # has to rebuild idx_features_upsert, not just ALTER the table.
        conn.execute("DROP INDEX IF EXISTS idx_features_upsert")
        conn.execute("ALTER TABLE features DROP COLUMN cds_start_OG_seq")
        conn.execute("ALTER TABLE features DROP COLUMN cds_end_OG_seq")
        conn.commit()
        cols = [r[1] for r in conn.execute("PRAGMA table_info(features)")]
        assert "cds_start_OG_seq" not in cols
    finally:
        conn.close()

    updated = build_db(tmp_path, make_inputs(tmp_path, "o2"), "old_db_2",
                       update=True, update_db=db, batch_id="batch_schema")

    conn = sqlite3.connect(f"file:{updated}?mode=ro", uri=True)
    try:
        cols = [r[1] for r in conn.execute("PRAGMA table_info(features)")]
    finally:
        conn.close()
    assert "cds_start_OG_seq" in cols and "cds_end_OG_seq" in cols


def test_backbone_widening_between_runs(tmp_path):
    """The stored alignment gets wider, so old rows are in a narrower space.

    Historic coordinates were computed against the old width. Nothing in the
    features table records which width a row belongs to, so both generations sit
    in one table indistinguishable from each other.
    """
    first = make_inputs(tmp_path, "w1", alignment="ATGC")
    db = build_db(tmp_path, first, "wide_db")

    wider = make_inputs(tmp_path, "w2", alignment="AT--GC",
                        feature_rows=[["A", "M", "R", "1", "6", "1", "6", "1", "4", "P", "1"]])
    updated = build_db(tmp_path, wider, "wide_db_2", update=True, update_db=db,
                       batch_id="batch_wide")

    conn = sqlite3.connect(f"file:{updated}?mode=ro", uri=True)
    try:
        widths = {len(r[0]) for r in conn.execute(
            "SELECT alignment FROM sequence_alignment WHERE primary_accession='A'")}
    finally:
        conn.close()
    # One row per accession per alignment_name/segment: the widened one wins.
    assert widths == {6}


def test_a_segment_the_database_has_never_seen(tmp_path):
    """A new segment arriving in an update must not disturb the existing one."""
    first = make_inputs(tmp_path, "g1", accession="A", segment="1")
    db = build_db(tmp_path, first, "seg_db")

    new_segment = make_inputs(tmp_path, "g2", accession="B", segment="4",
                              feature_rows=[["B", "M4", "R4", "1", "10", "1", "10", "1", "10", "P", "4"]],
                              meta_rows=[["B", "", "4"]])
    updated = build_db(tmp_path, new_segment, "seg_db_2", update=True, update_db=db,
                       batch_id="batch_newseg")

    conn = sqlite3.connect(f"file:{updated}?mode=ro", uri=True)
    try:
        rows = conn.execute(
            "SELECT accession, segment FROM features ORDER BY accession").fetchall()
    finally:
        conn.close()
    assert ("A", "1") in rows and ("B", "4") in rows


def test_master_for_a_segment_changes_between_runs(tmp_path):
    """A different master means a different coordinate frame for the same gene.

    Nothing rejects this, and the features table keeps both rows keyed by their
    own accession - so a segment ends up annotated against two frames at once.
    """
    first = make_inputs(tmp_path, "m1",
                        feature_rows=[["A", "OLD_MASTER", "R", "1", "10", "1", "10", "1", "10", "P", "1"]])
    db = build_db(tmp_path, first, "master_db")

    swapped = make_inputs(tmp_path, "m2",
                          feature_rows=[["A", "NEW_MASTER", "R", "1", "10", "1", "10", "1", "10", "P", "1"]])
    updated = build_db(tmp_path, swapped, "master_db_2", update=True, update_db=db,
                       batch_id="batch_master")

    conn = sqlite3.connect(f"file:{updated}?mode=ro", uri=True)
    try:
        masters = {r[0] for r in conn.execute(
            "SELECT master_ref_accession FROM features WHERE accession='A'")}
    finally:
        conn.close()
    # master_ref_accession is not part of the upsert key, so the row is replaced
    # and the old frame silently disappears rather than being reconciled.
    assert masters == {"NEW_MASTER"}


@pytest.mark.xfail(
    strict=True,
    reason="Nothing records which alignment width, master, or pipeline version a "
           "features row was computed under. After any coordinate change the table "
           "holds two generations of rows that are indistinguishable, so neither a "
           "consumer nor a migration can tell which rows need recomputing.",
)
def test_features_rows_record_their_provenance(tmp_path):
    first = make_inputs(tmp_path, "p1")
    db = build_db(tmp_path, first, "prov_db")

    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        cols = [r[1] for r in conn.execute("PRAGMA table_info(features)")]
    finally:
        conn.close()
    assert any(c in cols for c in ("batch_id", "pipeline_version", "alignment_width"))
