"""The meta_data column contract, and how it legitimately varies by run mode.

Three separate incidents in this repo came from schema drift going unnoticed:

* ``segment_name`` for neuraminidase was silently NULLed because a per-column
  rule was applied file-wide;
* the ``genes`` table shipped with a nameless Neuraminidase row for months;
* four update-mode tests silently skipped against a database that never existed,
  and when finally pointed at a real one they failed on fixtures that had drifted
  from every tracked schema.

Each was invisible because nothing asserted what a finished database is supposed
to contain. These tests pin the contract: which columns every run mode must
produce, and which are legitimately mode-specific. They skip cleanly when the
reference databases are absent (CI does not build them), so they are a guard for
local work rather than a CI gate.
"""

import re
import sqlite3
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT = REPO_ROOT / "test_out"

#: One reference database per run mode.
RUN_MODE_DBS = {
    "rabv_tree_free": OUT / "basic_test_treefree" / "rabv-jul0425.db",
    "rabv_update": OUT / "update_test" / "rabv-jul0425-update-test.db",
    "hcv_mutations": OUT / "HCV_OM_test" / "HCV_OM_test.db",
    "iav_segmented": OUT / "IAV_DB" / "iav-db.db",
}

#: Columns every run mode must produce. Measured across all four reference DBs.
#: A column dropping out of this set means some run mode stopped writing it.
#:
#: The MMseqs cluster column is deliberately NOT listed: its name encodes the
#: profile's clustering identity (mmseqs_min_seq_id 0.95 -> cluster_95pct,
#: 0.98 -> cluster_98pct), so it is not a fixed name. That every run mode has
#: exactly one such column is asserted separately, by pattern.
UNIVERSAL_META_COLUMNS = frozenset({
    "a", "c", "g", "n", "t",  # ATGCN composition counts
    "accession_type", "accession_version", "authors", "cds_info",
    "collection_date", "collection_date_validated", "collection_day",
    "collection_mon", "collection_year", "comment", "country", "country_validated",
    "create_date", "data_source", "db_xref", "definition", "division",
    "exclusion_criteria", "exclusion_status", "genes", "geo_loc", "gi_number",
    "host", "host_scientific_name", "host_taxa_id", "host_validated", "isolate",
    "isolation_source", "journal", "length", "locus", "mol_type", "molecule_type",
    "nearest_reference_genotype", "nearest_reference_subtype", "organism",
    "position", "primary_accession", "pubmed_id", "real_length",
    "reference_number", "segment", "serotype", "source", "strain", "strandedness",
    "taxonomy", "title", "topology", "update_date",
})

#: Columns produced by newer code that the older reference builds predate.
#: They are permitted anywhere but required nowhere, so a stale reference DB
#: does not fail the contract while a fresh one still declares them. Once every
#: reference DB has been rebuilt these belong in UNIVERSAL_META_COLUMNS.
#: Measured: present and fully populated (338/338) in the HCV build; absent from
#: the rabv and IAV builds, which predate the genotype-provenance work.
RECENTLY_ADDED_META_COLUMNS = frozenset({
    "genotype_origin", "subtype_origin",
})

#: Columns only a segmented (influenza) build produces, because the processes
#: that write them - VALIDATE_SEGMENT and VALIDATE_STRAIN - are gated on
#: is_segmented / is_flu in vgtk-init.nf.
SEGMENTED_ONLY_META_COLUMNS = frozenset({
    "Parsed_strain", "closest_reference", "exclusion", "segment_validated",
    "serotype_validated",
    # The segment_utils namespacing: `segment` is the numeric segment, while
    # `segment_name` is the gene label ('NA' for neuraminidase - a real value,
    # not a null) and `segment_source` records which stream supplied it. Only a
    # segmented build writes them; measured 498/518 populated in the IAV build
    # and absent from rabv and HCV.
    "segment_name", "segment_source",
})

#: Tables only a mutation_catalog run produces (HCV).
MUTATION_ONLY_TABLES = frozenset({
    "mutation_catalog", "sequence_relevant_mutation_summary",
    "completed_signatures_only",
})

#: Tables every finished database must have, whatever the run mode.
CORE_TABLES = frozenset({
    "meta_data", "features", "sequence_alignment", "sequences", "insertions",
    "genes", "host_taxa", "trees", "info",
})


def _connect(path):
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True)


def _columns(conn, table):
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}


def _tables(conn):
    return {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}


#: Artifacts built before a fix landed, so they cannot satisfy the current
#: contract. Rebuilding clears them; until then they would report a bug that no
#: longer exists. Keyed by run mode -> what is stale about it.
KNOWN_STALE = {
    "rabv_tree_free": (
        "built 2026-05-06, before CreateSqliteDB learned to fill blank segments "
        "with '1' for a non-segmented virus. Its meta_data.segment is 100% empty. "
        "Current code handles this (verified by "
        "test_non_segmented_build_edge_cases.py); rebuild the DB to clear this skip."
    ),
}


def _require(path):
    if not path.exists():
        pytest.skip(f"reference database not built here: {path}")
    return path


def _require_fresh(mode):
    if mode in KNOWN_STALE:
        pytest.skip(f"{mode}: {KNOWN_STALE[mode]}")
    return _require(RUN_MODE_DBS[mode])


@pytest.mark.parametrize("mode", sorted(RUN_MODE_DBS))
class TestEveryRunModeHonoursTheContract:
    def test_core_tables_all_present(self, mode):
        conn = _connect(_require(RUN_MODE_DBS[mode]))
        missing = CORE_TABLES - _tables(conn)
        conn.close()
        assert not missing, f"{mode} is missing core tables: {sorted(missing)}"

    def test_universal_meta_columns_all_present(self, mode):
        """A column vanishing from one run mode is exactly the drift that made
        the update-mode fixtures rot unnoticed."""
        conn = _connect(_require(RUN_MODE_DBS[mode]))
        missing = UNIVERSAL_META_COLUMNS - _columns(conn, "meta_data")
        conn.close()
        assert not missing, f"{mode} lost universal meta_data columns: {sorted(missing)}"

    def test_exactly_one_mmseqs_cluster_column(self, mode):
        """The name varies with the profile's clustering identity, but every run
        mode must produce exactly one - two would mean a stale column survived a
        threshold change, and none would mean clustering did not run."""
        conn = _connect(_require(RUN_MODE_DBS[mode]))
        cluster_cols = sorted(
            c for c in _columns(conn, "meta_data") if re.fullmatch(r"cluster_\d+pct", c))
        conn.close()
        assert len(cluster_cols) == 1, (
            f"{mode} has cluster columns {cluster_cols}; expected exactly one")

    def test_no_unexpected_extra_meta_columns(self, mode):
        """New columns are fine, but they must be declared here deliberately.

        This is the check that would have caught `Segment` arriving alongside
        `segment` from a vendor export.
        """
        conn = _connect(_require(RUN_MODE_DBS[mode]))
        cols = _columns(conn, "meta_data")
        conn.close()
        known = UNIVERSAL_META_COLUMNS | SEGMENTED_ONLY_META_COLUMNS
        # cluster_<N>pct is named from the profile's clustering identity, so it
        # is matched by pattern rather than declared by name.
        unexpected = {
            c for c in cols - known - RECENTLY_ADDED_META_COLUMNS
            if not c.startswith("gisaid_") and not re.fullmatch(r"cluster_\d+pct", c)
        }
        assert not unexpected, (
            f"{mode} has undeclared meta_data columns {sorted(unexpected)}. "
            f"If deliberate, add them to the contract in this module."
        )

    def test_no_two_columns_differ_only_by_case(self, mode):
        """SQLite cannot represent it, so a DB containing it is impossible -
        but the guard is cheap and states the invariant."""
        conn = _connect(_require(RUN_MODE_DBS[mode]))
        for table in sorted(CORE_TABLES & _tables(conn)):
            cols = [row[1] for row in conn.execute(f"PRAGMA table_info({table})")]
            folded = [c.casefold() for c in cols]
            assert len(folded) == len(set(folded)), f"{mode}.{table} has a case collision"
        conn.close()

    def test_meta_data_key_is_unique(self, mode):
        """(primary_accession, segment) is the upsert key. A duplicate means an
        update would replace an arbitrary one of them."""
        conn = _connect(_require(RUN_MODE_DBS[mode]))
        dups = conn.execute(
            "SELECT COUNT(*) FROM (SELECT primary_accession, segment FROM meta_data "
            "GROUP BY 1, 2 HAVING COUNT(*) > 1)"
        ).fetchone()[0]
        conn.close()
        assert dups == 0, f"{mode} has {dups} duplicate (primary_accession, segment) keys"


class TestModeSpecificColumnsStayModeSpecific:
    def test_segmented_columns_absent_from_non_segmented_builds(self):
        for mode in ("rabv_tree_free", "rabv_update", "hcv_mutations"):
            path = RUN_MODE_DBS[mode]
            if not path.exists():
                continue
            conn = _connect(path)
            present = SEGMENTED_ONLY_META_COLUMNS & _columns(conn, "meta_data")
            conn.close()
            assert not present, (
                f"{mode} is non-segmented but carries segmented-only columns {sorted(present)}"
            )

    def test_segmented_build_has_them(self):
        conn = _connect(_require(RUN_MODE_DBS["iav_segmented"]))
        cols = _columns(conn, "meta_data")
        conn.close()
        missing = SEGMENTED_ONLY_META_COLUMNS - cols
        assert not missing, f"segmented build is missing {sorted(missing)}"

    def test_mutation_tables_only_where_a_catalog_was_supplied(self):
        for mode, path in RUN_MODE_DBS.items():
            if not path.exists():
                continue
            conn = _connect(path)
            present = MUTATION_ONLY_TABLES & _tables(conn)
            conn.close()
            if mode == "hcv_mutations":
                assert present == MUTATION_ONLY_TABLES, f"{mode} missing {sorted(MUTATION_ONLY_TABLES - present)}"
            else:
                assert not present, f"{mode} has mutation tables without a catalog: {sorted(present)}"

    def test_non_segmented_builds_use_segment_one_throughout(self):
        """`_force_segment_one_df` should leave no other value behind."""
        for mode in ("rabv_tree_free", "rabv_update", "hcv_mutations"):
            path = RUN_MODE_DBS[mode]
            if not path.exists() or mode in KNOWN_STALE:
                continue
            conn = _connect(path)
            values = {r[0] for r in conn.execute("SELECT DISTINCT segment FROM meta_data")}
            conn.close()
            assert values <= {"1"}, f"{mode} is non-segmented but has segments {sorted(values)}"


class TestFeatureKeyDistinguishesSplicedGenes:
    """The features upsert key must keep the CDS coordinates.

    Influenza M2 and NEP are spliced and PA-X is produced by a ribosomal
    frameshift, so each has TWO CDS intervals under one product name. In a real
    IAV database that is 136 (accession, product, segment) groups holding two
    rows each. The upsert index is on
    (accession, cds_start_OG_seq, cds_end_OG_seq, product, segment), which tells
    the exons apart - dropping the coordinates from that key would silently
    collapse every spliced gene to one interval.
    """

    def test_product_alone_is_not_unique_for_spliced_genes(self):
        conn = _connect(_require(RUN_MODE_DBS["iav_segmented"]))
        groups = conn.execute(
            "SELECT COUNT(*) FROM (SELECT accession, product, segment FROM features "
            "GROUP BY 1, 2, 3 HAVING COUNT(*) > 1)"
        ).fetchone()[0]
        products = {
            r[0] for r in conn.execute(
                "SELECT DISTINCT product FROM (SELECT accession, product, segment FROM features "
                "GROUP BY 1, 2, 3 HAVING COUNT(*) > 1)"
            )
        }
        conn.close()
        assert groups > 0, "expected spliced influenza genes to share a product name"
        assert products <= {"PA-X protein", "matrix protein 2", "nuclear export protein"}, (
            f"unexpected products with multiple CDS intervals: {sorted(products)}"
        )

    def test_the_real_key_including_coordinates_is_unique(self):
        conn = _connect(_require(RUN_MODE_DBS["iav_segmented"]))
        dups = conn.execute(
            "SELECT COUNT(*) FROM (SELECT accession, cds_start_OG_seq, cds_end_OG_seq, "
            "product, segment FROM features GROUP BY 1, 2, 3, 4, 5 HAVING COUNT(*) > 1)"
        ).fetchone()[0]
        conn.close()
        assert dups == 0, f"{dups} duplicate rows on the features upsert key"


class TestNoColumnIsWhollyEmpty:
    """A column that is 100% NULL means something stopped writing it.

    The shipped IAV database stored the neuraminidase gene as
    ('Neuraminidase', None, '', 'whole_genome') for months without anyone
    noticing, because nothing asserted that populated columns stay populated.
    """

    #: Legitimately sparse or optional; excluded from the all-empty check.
    ALLOWED_EMPTY = frozenset({
        "comment", "db_xref", "cds_info", "isolation_source", "geo_loc",
        "serotype", "serotype_validated", "Parsed_strain", "segment_validated",
        "closest_reference", "exclusion", "collection_date_validated",
        "nearest_reference_genotype", "nearest_reference_subtype",
        "collection_day", "collection_mon", "collection_year", "collection_date",
        "host_taxa_id", "host_scientific_name", "host_validated",
        "country_validated", "isolate", "strain", "journal", "title",
        "pubmed_id", "authors", "reference_number", "position",
        # Only populated when something was actually excluded: 0/338 in the HCV
        # build and 0/78 in rabv_tree_free, 1/228 in rabv_update, 22/518 in IAV.
        # An empty column here means "nothing was excluded", not a broken write.
        "exclusion_criteria",
    })

    @pytest.mark.parametrize("mode", sorted(RUN_MODE_DBS))
    def test_expected_columns_are_populated(self, mode):
        conn = _connect(_require_fresh(mode))
        cols = [r[1] for r in conn.execute("PRAGMA table_info(meta_data)")]
        checked = [c for c in cols if c not in self.ALLOWED_EMPTY]
        selects = ", ".join(
            f"SUM(CASE WHEN TRIM(COALESCE(CAST(\"{c}\" AS TEXT), '')) != '' THEN 1 ELSE 0 END)"
            for c in checked
        )
        counts = conn.execute(f"SELECT {selects} FROM meta_data").fetchone()
        total = conn.execute("SELECT COUNT(*) FROM meta_data").fetchone()[0]
        conn.close()
        if total == 0:
            pytest.skip("empty meta_data")
        empty = [c for c, n in zip(checked, counts) if n == 0]
        assert not empty, f"{mode}: columns are 100% empty: {empty}"
