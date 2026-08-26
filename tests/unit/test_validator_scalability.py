"""The DB validator has to survive a production-sized database.

``validate_meta_column_population`` originally issued one ``SELECT COUNT(*)``
per meta_data column. That is one full table scan per column, and meta_data is
wide. On a real influenza build - 106 columns, 3.9 M rows, 28 GB - it performed
roughly 3 TB of logical reads and had not finished after half an hour, at which
point it had to be killed. The invariant it checks is worth having, but a check
that cannot run on a real database provides nothing.

These tests pin the cost model rather than the wall-clock time, so they stay
meaningful on any machine: the number of queries must not grow with the number
of columns, and the number of table scans must not grow with the number of rows.
"""

import sqlite3

import pytest

import ValidateDbTree as V


class CountingConnection:
    """A sqlite3 connection that records every statement executed through it."""

    def __init__(self, conn):
        self._conn = conn
        self.statements = []

    def execute(self, sql, *args, **kwargs):
        self.statements.append(sql)
        return self._conn.execute(sql, *args, **kwargs)

    def cursor(self):
        return CountingCursor(self._conn.cursor(), self.statements)

    def __getattr__(self, name):
        return getattr(self._conn, name)


class CountingCursor:
    def __init__(self, cursor, sink):
        self._cursor = cursor
        self._sink = sink

    def execute(self, sql, *args, **kwargs):
        self._sink.append(sql)
        return self._cursor.execute(sql, *args, **kwargs)

    def __getattr__(self, name):
        return getattr(self._cursor, name)


def _make_db(n_columns, n_rows, populated=True):
    """A meta_data table of the requested shape."""
    conn = sqlite3.connect(":memory:")
    cols = ["primary_accession"] + [f"col_{i}" for i in range(n_columns - 1)]
    conn.execute(
        "CREATE TABLE meta_data (" + ", ".join(f'"{c}" TEXT' for c in cols) + ")"
    )
    row = tuple(f"v{i}" if populated else "" for i in range(n_columns))
    placeholders = ", ".join("?" for _ in cols)
    conn.executemany(
        f"INSERT INTO meta_data VALUES ({placeholders})", [row] * n_rows
    )
    conn.commit()
    return conn, cols


def _flagged(result):
    """Every column the check reported as entirely empty.

    Core and reference-label columns are surfaced as findings (with `examples`);
    everything else is listed under the top-level `empty_other_columns`.
    """
    flagged = set(result.get("empty_other_columns") or [])
    for finding in result.get("findings", []):
        flagged.update(finding.get("examples") or [])
    return flagged


def _scan_count(statements):
    """How many statements read the whole meta_data table."""
    return sum(1 for s in statements if "FROM meta_data" in s)


class TestColumnPopulationCostDoesNotGrowWithColumns:
    """The regression that killed a 28 GB validation run."""

    @pytest.mark.parametrize("n_columns", [8, 40, 120])
    def test_scan_count_is_bounded_regardless_of_width(self, n_columns):
        conn, cols = _make_db(n_columns, n_rows=5)
        counting = CountingConnection(conn)
        V.validate_meta_column_population(counting, cols)
        scans = _scan_count(counting.statements)
        # Chunked at 64 columns per statement, so a 120-column table takes two.
        assert scans <= 4, (
            f"{n_columns} columns caused {scans} scans of meta_data; the check "
            f"must not scan once per column"
        )

    def test_a_wide_table_costs_no_more_scans_than_a_narrow_one_per_column(self):
        narrow_conn, narrow_cols = _make_db(8, n_rows=5)
        wide_conn, wide_cols = _make_db(120, n_rows=5)

        narrow = CountingConnection(narrow_conn)
        wide = CountingConnection(wide_conn)
        V.validate_meta_column_population(narrow, narrow_cols)
        V.validate_meta_column_population(wide, wide_cols)

        narrow_scans = _scan_count(narrow.statements)
        wide_scans = _scan_count(wide.statements)
        # 15x the columns must not mean anything like 15x the scans.
        assert wide_scans <= narrow_scans + 2, (
            f"narrow={narrow_scans} scans, wide={wide_scans} scans - cost is "
            f"growing with column count"
        )

    def test_scan_count_does_not_grow_with_rows(self):
        small_conn, cols = _make_db(20, n_rows=5)
        big_conn, _ = _make_db(20, n_rows=2000)
        small = CountingConnection(small_conn)
        big = CountingConnection(big_conn)
        V.validate_meta_column_population(small, cols)
        V.validate_meta_column_population(big, cols)
        assert _scan_count(small.statements) == _scan_count(big.statements)


class TestColumnPopulationStillDetectsWhatItIsFor:
    """Speed is worthless if the check stopped working."""

    def test_an_entirely_empty_column_is_flagged(self):
        conn = sqlite3.connect(":memory:")
        conn.execute('CREATE TABLE meta_data ("primary_accession" TEXT, "host" TEXT, "organism" TEXT)')
        conn.executemany("INSERT INTO meta_data VALUES (?, ?, ?)",
                         [("A", "", "Influenza A"), ("B", "", "Influenza A")])
        conn.commit()
        cols = ["primary_accession", "host", "organism"]
        result = V.validate_meta_column_population(conn, cols)
        assert "host" in _flagged(result)
        assert "organism" not in _flagged(result)
        assert "primary_accession" not in _flagged(result)

    def test_a_column_with_one_value_is_not_flagged(self):
        """The check is 'entirely empty', not 'mostly empty'."""
        conn = sqlite3.connect(":memory:")
        conn.execute('CREATE TABLE meta_data ("primary_accession" TEXT, "serotype" TEXT)')
        conn.executemany("INSERT INTO meta_data VALUES (?, ?)",
                         [("A", ""), ("B", ""), ("C", "H1N1")])
        conn.commit()
        result = V.validate_meta_column_population(conn, ["primary_accession", "serotype"])
        assert "serotype" not in _flagged(result)

    def test_chunk_boundary_columns_are_not_dropped(self):
        """Chunking must cover every column, including across a chunk edge."""
        n = 130
        conn = sqlite3.connect(":memory:")
        cols = [f"col_{i}" for i in range(n)]
        conn.execute("CREATE TABLE meta_data (" + ", ".join(f'"{c}" TEXT' for c in cols) + ")")
        # Every column blank except the two either side of the 64-column boundary.
        row = ["" for _ in cols]
        row[63] = "x"
        row[64] = "y"
        conn.execute(f"INSERT INTO meta_data VALUES ({', '.join('?' for _ in cols)})", row)
        conn.commit()
        result = V.validate_meta_column_population(conn, cols)
        flagged = _flagged(result)
        assert "col_63" not in flagged and "col_64" not in flagged
        assert "col_0" in flagged and f"col_{n - 1}" in flagged
        assert len(flagged) == n - 2
