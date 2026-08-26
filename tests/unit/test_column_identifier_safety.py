"""SQLite identifier hazards at the DB-write boundary.

SQLite column identifiers are case-insensitive; pandas column labels are not. A
frame carrying both ``segment`` and ``Segment`` is legal in pandas and illegal in
SQLite, and a real GenBank+GISAID influenza build produced exactly that: GISAID
ships a declared ``Segment`` column against the pipeline's canonical ``segment``.
The build died with ``duplicate column name: Segment`` after the trees were
already computed, and the error named neither the two columns nor the input.

The spaced-identifier case is worse than an error: GISAID also ships ten headers
of the form ``HA INSDC_Upload``, and an unquoted
``ALTER TABLE t ADD COLUMN HA INSDC_Upload TEXT`` does not fail - SQLite reads
``HA`` as the column name and ``INSDC_Upload TEXT`` as its type, silently
creating the wrong column.
"""

import sqlite3

import pandas as pd
import pytest

from CreateSqliteDB import CreateSqliteDB


class TestCaseInsensitiveCollisionGuard:
    def test_collision_is_rejected_and_names_both_columns(self):
        frame = pd.DataFrame({"primary_accession": ["A1"]})
        frame["segment"] = "4"
        frame.insert(2, "Segment", "HA")

        with pytest.raises(ValueError) as excinfo:
            CreateSqliteDB.assert_no_case_insensitive_duplicates(frame.columns, "incoming meta_data")

        message = str(excinfo.value)
        assert "incoming meta_data" in message
        assert "'segment'" in message and "'Segment'" in message
        # The positions matter: a 101-column matrix is not greppable by eye.
        assert "col 2" in message and "col 3" in message

    def test_clean_columns_are_accepted(self):
        CreateSqliteDB.assert_no_case_insensitive_duplicates(
            ["primary_accession", "segment", "segment_name", "gisaid_segment"], "meta_data"
        )

    def test_the_namespaced_form_of_the_real_collision_is_accepted(self):
        """What the merge now produces must pass the guard."""
        CreateSqliteDB.assert_no_case_insensitive_duplicates(
            ["primary_accession", "segment", "segment_declared", "gisaid_genotype"], "meta_data"
        )

    def test_guard_catches_what_sqlite_would_have_caught_cryptically(self):
        """Cross-check: the guard rejects exactly the frames sqlite rejects."""
        frame = pd.DataFrame({"primary_accession": ["A1"]})
        frame["segment"] = "4"
        frame.insert(2, "Segment", "HA")

        with sqlite3.connect(":memory:") as conn:
            with pytest.raises(sqlite3.OperationalError, match="duplicate column name"):
                frame.to_sql("meta_data", conn, if_exists="replace", index=False)

        with pytest.raises(ValueError):
            CreateSqliteDB.assert_no_case_insensitive_duplicates(frame.columns, "meta_data")


class TestIdentifierQuoting:
    @pytest.mark.parametrize("name", [
        "HA INSDC_Upload", "NA INSDC_Upload", "PB1 INSDC_Upload",
        'weird"quote', "select", "group by",
    ])
    def test_quoted_identifier_creates_exactly_that_column(self, name):
        conn = sqlite3.connect(":memory:")
        conn.execute('CREATE TABLE t ("primary_accession" TEXT)')
        conn.execute(f"ALTER TABLE t ADD COLUMN {CreateSqliteDB._quote_identifier(name)} TEXT")

        columns = [row[1] for row in conn.execute("PRAGMA table_info(t)")]
        assert name in columns
        conn.close()

    def test_unquoted_spaced_identifier_would_have_created_the_wrong_column(self):
        """Documents the hazard the quoting exists to prevent."""
        conn = sqlite3.connect(":memory:")
        conn.execute('CREATE TABLE t ("primary_accession" TEXT)')
        conn.execute("ALTER TABLE t ADD COLUMN HA INSDC_Upload TEXT")  # deliberately unquoted

        schema = {row[1]: row[2] for row in conn.execute("PRAGMA table_info(t)")}
        assert "HA INSDC_Upload" not in schema
        assert schema["HA"] == "INSDC_Upload TEXT"  # silently wrong, no error raised
        conn.close()


class TestEnsureUpdateColumns:
    """`_ensure_update_columns` re-read the schema once, outside its loop."""

    @staticmethod
    def _updater():
        instance = CreateSqliteDB.__new__(CreateSqliteDB)
        instance.update = True
        return instance

    def test_repeated_label_does_not_issue_the_same_alter_twice(self):
        conn = sqlite3.connect(":memory:")
        conn.execute('CREATE TABLE meta_data ("primary_accession" TEXT)')
        frame = pd.DataFrame([["A1", "x", "y"]], columns=["primary_accession", "genotype", "genotype"])

        # The duplicate label is caught up front rather than reaching ALTER TABLE.
        with pytest.raises(ValueError, match="differ only by case"):
            self._updater()._ensure_update_columns(conn, "meta_data", frame)
        conn.close()

    def test_column_differing_only_by_case_from_the_existing_schema_is_not_re_added(self):
        conn = sqlite3.connect(":memory:")
        conn.execute('CREATE TABLE meta_data ("primary_accession" TEXT, "Segment" TEXT)')
        frame = pd.DataFrame({"primary_accession": ["A1"], "segment": ["4"]})

        self._updater()._ensure_update_columns(conn, "meta_data", frame)

        columns = [row[1] for row in conn.execute("PRAGMA table_info(meta_data)")]
        assert columns == ["primary_accession", "Segment"]
        conn.close()

    def test_genuinely_new_columns_are_added(self):
        conn = sqlite3.connect(":memory:")
        conn.execute('CREATE TABLE meta_data ("primary_accession" TEXT)')
        frame = pd.DataFrame({
            "primary_accession": ["A1"], "segment": ["4"],
            "segment_name": ["HA"], "segment_source": ["blast_inferred"],
        })

        self._updater()._ensure_update_columns(conn, "meta_data", frame)

        columns = [row[1] for row in conn.execute("PRAGMA table_info(meta_data)")]
        assert columns == ["primary_accession", "segment", "segment_name", "segment_source"]
        conn.close()

    def test_spaced_vendor_column_is_added_under_its_real_name(self):
        conn = sqlite3.connect(":memory:")
        conn.execute('CREATE TABLE meta_data ("primary_accession" TEXT)')
        frame = pd.DataFrame({"primary_accession": ["A1"], "HA INSDC_Upload": ["EPI1"]})

        self._updater()._ensure_update_columns(conn, "meta_data", frame)

        schema = {row[1]: row[2] for row in conn.execute("PRAGMA table_info(meta_data)")}
        assert "HA INSDC_Upload" in schema
        assert "HA" not in schema
        conn.close()
