"""Tests for the reference-set recombination screen.

None of these run UShER or RIPPLES. The screen's correctness lives in three
places that are all pure functions - how the reference alignment is extracted,
how RIPPLES output is parsed, and how results reach the database - and those are
what is tested here. The subprocess calls are thin wrappers whose only real
contract is the command line they build, which is asserted directly.
"""

import sqlite3
import sys
from pathlib import Path

import pytest

import ScreenReferenceRecombination as SRR


REPO_ROOT = Path(__file__).resolve().parents[2]


# --------------------------------------------------------------------------
# Thread ceiling - a shared machine, and both tools grab every core by default
# --------------------------------------------------------------------------

class TestThreadCeiling:
    def test_ceiling_is_sixteen(self):
        assert SRR.MAX_THREADS == 16

    @pytest.mark.parametrize("requested,expected", [
        (1, 1), (8, 8), (16, 16), (17, 16), (48, 16), (224, 16),
        (0, 1), (-4, 1), (None, 1), ("nonsense", 1),
    ])
    def test_clamped(self, requested, expected):
        assert SRR.clamp_threads(requested) == expected

    def test_ripples_command_never_exceeds_the_ceiling(self):
        command = SRR.build_ripples_command("m.pb", "out", 224)
        assert command[command.index("-T") + 1] == "16"

    def test_usher_command_never_exceeds_the_ceiling(self):
        command = SRR.build_usher_command("t.nwk", "v.vcf", "o.pb", "d", 224)
        assert command[command.index("-T") + 1] == "16"

    def test_the_class_clamps_and_says_so(self, capsys):
        screen = SRR.ReferenceRecombinationScreen(threads=64)
        assert screen.threads == 16
        assert "clamped" in capsys.readouterr().out

    def test_a_reasonable_request_is_left_alone_and_silent(self, capsys):
        screen = SRR.ReferenceRecombinationScreen(threads=8)
        assert screen.threads == 8
        assert "clamped" not in capsys.readouterr().out


# --------------------------------------------------------------------------
# Breakpoint intervals
# --------------------------------------------------------------------------

class TestParseInterval:
    @pytest.mark.parametrize("text,expected", [
        ("1234-5678", (1234, 5678)),
        (" 1234 - 5678 ", (1234, 5678)),
        ("42", (42, 42)),
        ("0-0", (0, 0)),
    ])
    def test_valid(self, text, expected):
        assert SRR.parse_interval(text) == expected

    @pytest.mark.parametrize("text", ["", "   ", None, "NA", "abc-def", "1234-", "-5678", "1-2-3"])
    def test_unparseable_is_none_not_an_exception(self, text):
        """A malformed interval must not cost us the rest of the event."""
        assert SRR.parse_interval(text) == (None, None)


# --------------------------------------------------------------------------
# RIPPLES output
# --------------------------------------------------------------------------

HEADER = "#" + "\t".join(SRR.RECOMBINATION_COLUMNS)
EVENT = "\t".join([
    "node_42", "3200-3260", "7100-7180", "node_7", "y", "120",
    "node_19", "n", "133", "410", "260", "253",
])


@pytest.fixture
def ripples_dir(tmp_path):
    (tmp_path / "recombination.tsv").write_text(HEADER + "\n" + EVENT + "\n")
    (tmp_path / "descendants.tsv").write_text(
        "#node_id\tdescendants\nnode_42\tMK548365,KU746825\n")
    return tmp_path


class TestParseRecombination:
    def test_header_and_one_event(self, ripples_dir):
        events = SRR.parse_recombination_tsv(ripples_dir / "recombination.tsv")
        assert len(events) == 1
        e = events[0]
        assert e["recomb_node_id"] == "node_42"
        assert e["donor_node_id"] == "node_7"
        assert e["acceptor_node_id"] == "node_19"

    def test_intervals_are_split_into_numbers(self, ripples_dir):
        e = SRR.parse_recombination_tsv(ripples_dir / "recombination.tsv")[0]
        assert (e["breakpoint_1_start"], e["breakpoint_1_end"]) == (3200, 3260)
        assert (e["breakpoint_2_start"], e["breakpoint_2_end"]) == (7100, 7180)

    def test_hyphenated_columns_are_renamed_for_sql(self, ripples_dir):
        """`breakpoint-1_interval` cannot be an unquoted SQL identifier."""
        e = SRR.parse_recombination_tsv(ripples_dir / "recombination.tsv")[0]
        assert "breakpoint_1_interval" in e
        assert "breakpoint-1_interval" not in e

    def test_missing_file_is_empty_not_an_error(self, tmp_path):
        assert SRR.parse_recombination_tsv(tmp_path / "nope.tsv") == []

    def test_empty_file_is_empty(self, tmp_path):
        """RIPPLES creates the file up front and leaves it empty when it finds
        nothing - that is 'no evidence', not a failure."""
        p = tmp_path / "recombination.tsv"; p.write_text("")
        assert SRR.parse_recombination_tsv(p) == []

    def test_header_only_is_empty(self, tmp_path):
        p = tmp_path / "recombination.tsv"; p.write_text(HEADER + "\n")
        assert SRR.parse_recombination_tsv(p) == []

    def test_a_headerless_file_is_still_read(self, tmp_path):
        p = tmp_path / "recombination.tsv"; p.write_text(EVENT + "\n")
        events = SRR.parse_recombination_tsv(p)
        assert len(events) == 1 and events[0]["recomb_node_id"] == "node_42"

    def test_short_rows_are_padded_rather_than_dropped(self, tmp_path):
        p = tmp_path / "recombination.tsv"
        p.write_text(HEADER + "\nnode_9\t100-200\t300-400\n")
        events = SRR.parse_recombination_tsv(p)
        assert len(events) == 1
        assert events[0]["recomb_node_id"] == "node_9"
        assert events[0]["donor_node_id"] == ""

    def test_column_set_matches_the_binary(self):
        """These names were read out of the RIPPLES binary. If a version bump
        changes them this test is the tripwire."""
        assert SRR.RECOMBINATION_COLUMNS[0] == "recomb_node_id"
        assert SRR.RECOMBINATION_COLUMNS[-1] == "recomb_parsimony"
        assert len(SRR.RECOMBINATION_COLUMNS) == 12


class TestParseDescendants:
    def test_comma_separated(self, ripples_dir):
        d = SRR.parse_descendants_tsv(ripples_dir / "descendants.tsv")
        assert d == {"node_42": ["MK548365", "KU746825"]}

    def test_whitespace_separated(self, tmp_path):
        p = tmp_path / "descendants.tsv"
        p.write_text("#node_id\tdescendants\nnode_1\tA B  C\n")
        assert SRR.parse_descendants_tsv(p) == {"node_1": ["A", "B", "C"]}

    def test_missing_file_is_empty(self, tmp_path):
        assert SRR.parse_descendants_tsv(tmp_path / "nope.tsv") == {}

    def test_repeated_node_accumulates(self, tmp_path):
        p = tmp_path / "descendants.tsv"
        p.write_text("#node_id\tdescendants\nnode_1\tA\nnode_1\tB\n")
        assert SRR.parse_descendants_tsv(p) == {"node_1": ["A", "B"]}


class TestAttribution:
    def test_one_row_per_descendant(self, ripples_dir):
        events = SRR.parse_recombination_tsv(ripples_dir / "recombination.tsv")
        desc = SRR.parse_descendants_tsv(ripples_dir / "descendants.tsv")
        rows = SRR.attribute_events_to_accessions(events, desc)
        assert {r["primary_accession"] for r in rows} == {"MK548365", "KU746825"}
        assert all(r["recomb_node_id"] == "node_42" for r in rows)

    def test_an_event_with_no_descendants_is_kept_unattributed(self, ripples_dir):
        """Dropping a detected event because descendants.tsv was incomplete
        would be worse than reporting it without an accession."""
        events = SRR.parse_recombination_tsv(ripples_dir / "recombination.tsv")
        rows = SRR.attribute_events_to_accessions(events, {})
        assert len(rows) == 1 and rows[0]["primary_accession"] is None

    def test_no_events_gives_no_rows(self):
        assert SRR.attribute_events_to_accessions([], {"n": ["A"]}) == []


# --------------------------------------------------------------------------
# Chunking - the only practical way to run this on a diverse virus
# --------------------------------------------------------------------------

class TestChunkBounds:
    def test_covers_every_branch_exactly_once(self):
        bounds = SRR.chunk_bounds(156, 16)
        assert sum(b - a for a, b in bounds) == 156
        assert bounds[0][0] == 0 and bounds[-1][1] == 156
        for (_, end), (start, _) in zip(bounds, bounds[1:]):
            assert end == start

    def test_sizes_differ_by_at_most_one(self):
        sizes = [b - a for a, b in SRR.chunk_bounds(156, 16)]
        assert max(sizes) - min(sizes) <= 1

    def test_more_chunks_than_branches_is_capped(self):
        bounds = SRR.chunk_bounds(3, 16)
        assert len(bounds) == 3
        assert all(b - a == 1 for a, b in bounds)

    @pytest.mark.parametrize("n,k", [(0, 4), (10, 0), (-1, 4)])
    def test_degenerate_inputs(self, n, k):
        assert SRR.chunk_bounds(n, k) == []

    def test_indices_reach_ripples_as_S_and_E(self):
        cmd = SRR.build_ripples_command("m.pb", "o", 8, start_index=10, end_index=20)
        assert cmd[cmd.index("-S") + 1] == "10"
        assert cmd[cmd.index("-E") + 1] == "20"

    def test_indices_are_omitted_when_not_set(self):
        cmd = SRR.build_ripples_command("m.pb", "o", 8)
        assert "-S" not in cmd and "-E" not in cmd


class TestCountLongBranches:
    def test_reads_the_number(self):
        assert SRR.count_long_branches("Found 156 long branches\nCompleted") == 156

    def test_absent_is_none(self):
        assert SRR.count_long_branches("nothing here") is None
        assert SRR.count_long_branches("") is None
        assert SRR.count_long_branches(None) is None


# --------------------------------------------------------------------------
# Defaults chosen deliberately
# --------------------------------------------------------------------------

class TestDefaults:
    def test_num_descendants_is_lowered_from_the_ripples_default(self):
        """RIPPLES defaults to -n 10, tuned for SARS-CoV-2. A curated reference
        set has a handful of sequences per subtype, so 10 would skip most of
        the tree."""
        cmd = SRR.build_ripples_command("m.pb", "o", 8)
        assert cmd[cmd.index("-n") + 1] == "3"

    def test_every_ripples_parameter_is_passed_explicitly(self):
        """Never rely on a tool's defaults - they change between versions."""
        cmd = SRR.build_ripples_command("m.pb", "o", 8)
        for flag in ("-i", "-d", "-T", "-l", "-r", "-R", "-p", "-n"):
            assert flag in cmd

    def test_fatovcf_command(self):
        assert SRR.build_fatovcf_command("a.fa", "b.vcf") == ["faToVcf", "a.fa", "b.vcf"]


# --------------------------------------------------------------------------
# Reading the reference set out of a database
# --------------------------------------------------------------------------

def _make_db(path, rows, alignments):
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment TEXT)")
    conn.execute("CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, "
                 "segment TEXT, newick TEXT, created_at TEXT)")
    conn.executemany("INSERT INTO meta_data VALUES (?,?)", rows)
    conn.executemany("INSERT INTO sequence_alignment VALUES (?,?)", alignments)
    conn.commit()
    return conn


@pytest.fixture
def small_db(tmp_path):
    path = tmp_path / "t.db"
    conn = _make_db(
        str(path),
        [("REF0", "master"), ("REF1", "reference"), ("REF2", "reference"),
         ("Q1", "query"), ("Q2", "query")],
        [("REF0", "ACGTACGT"), ("REF1", "ACGTACGA"), ("REF2", "ACGTTCGT"),
         ("Q1", "ACGTACGT"), ("Q2", "ACGTACGT")],
    )
    conn.execute("INSERT INTO trees VALUES ('iq','iqtree',NULL,'1',"
                 "'((REF0:0.1,REF1:0.1):0.1,(REF2:0.1,Q1:0.1):0.1);','now')")
    conn.commit()
    return path, conn


class TestReferenceRecordsFromDb:
    def test_returns_master_and_reference_alignments_only(self, small_db):
        _, conn = small_db
        master, seqs = SRR.reference_records_from_db(conn)
        assert master == "REF0"
        assert set(seqs) == {"REF0", "REF1", "REF2"}, "queries must not be included"

    def test_no_master_is_an_error(self, tmp_path):
        conn = _make_db(str(tmp_path / "a.db"),
                        [("R1", "reference")], [("R1", "ACGT")])
        with pytest.raises(ValueError, match="master"):
            SRR.reference_records_from_db(conn)

    def test_no_reference_alignments_is_an_error(self, tmp_path):
        conn = _make_db(str(tmp_path / "b.db"), [("M", "master")], [])
        with pytest.raises(ValueError):
            SRR.reference_records_from_db(conn)

    def test_ragged_alignments_are_refused(self, tmp_path):
        """faToVcf needs a rectangular alignment; a ragged one would silently
        produce nonsense coordinates."""
        conn = _make_db(str(tmp_path / "c.db"),
                        [("M", "master"), ("R", "reference")],
                        [("M", "ACGT"), ("R", "ACGTAC")])
        with pytest.raises(ValueError, match="ragged"):
            SRR.reference_records_from_db(conn)


class TestWriteReferenceFasta:
    def test_master_is_written_first(self, small_db, tmp_path):
        """faToVcf takes the FIRST record as the VCF reference, so the order
        determines the coordinate frame - it is not cosmetic."""
        _, conn = small_db
        master, seqs = SRR.reference_records_from_db(conn)
        out = tmp_path / "r.fasta"
        SRR.write_reference_fasta(master, seqs, out)
        names = [l[1:].strip() for l in out.read_text().splitlines() if l.startswith(">")]
        assert names[0] == "REF0"
        assert sorted(names) == ["REF0", "REF1", "REF2"]

    def test_every_sequence_appears_once(self, small_db, tmp_path):
        _, conn = small_db
        master, seqs = SRR.reference_records_from_db(conn)
        out = tmp_path / "r.fasta"
        SRR.write_reference_fasta(master, seqs, out)
        text = out.read_text()
        assert text.count(">REF0") == 1


# --------------------------------------------------------------------------
# Starting tree
# --------------------------------------------------------------------------

class TestPruneNewick:
    def test_restricts_to_requested_tips(self):
        out = SRR.prune_newick_to("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);", {"A", "C"})
        assert set("".join(c for c in out if c.isalpha())) == {"A", "C"}

    def test_unary_nodes_are_suppressed(self):
        out = SRR.prune_newick_to("((A:0.1,B:0.1):0.2,C:0.3);", {"A", "C"})
        assert out == "(A,C);"

    def test_branch_lengths_and_support_are_dropped(self):
        """UShER re-derives branch lengths from the VCF; only topology is used."""
        out = SRR.prune_newick_to("((A:0.1,B:0.2)95:0.3,C:0.4);", {"A", "B", "C"})
        assert ":" not in out and "95" not in out

    def test_no_overlap_is_an_error(self):
        with pytest.raises(ValueError):
            SRR.prune_newick_to("((A,B),C);", {"X", "Y"})

    def test_quoted_labels_are_handled(self):
        out = SRR.prune_newick_to("(('A':0.1,'B':0.1),C);", {"A", "B"})
        assert out == "(A,B);"


class TestStartingTreeFromDb:
    def test_prefers_iqtree_over_usher(self, small_db):
        _, conn = small_db
        conn.execute("INSERT INTO trees VALUES ('ush','usher',NULL,'1',"
                     "'((REF0,REF1),(REF2,Q1));','now')")
        conn.commit()
        name, _ = SRR.starting_tree_from_db(conn, {"REF0", "REF1", "REF2"})
        assert name == "iq"

    def test_result_is_pruned_to_the_reference_set(self, small_db):
        _, conn = small_db
        _, newick = SRR.starting_tree_from_db(conn, {"REF0", "REF1", "REF2"})
        assert "Q1" not in newick

    def test_no_trees_is_an_error(self, tmp_path):
        conn = _make_db(str(tmp_path / "d.db"), [("M", "master")], [("M", "ACGT")])
        with pytest.raises(ValueError, match="trees"):
            SRR.starting_tree_from_db(conn, {"M"})


# --------------------------------------------------------------------------
# Writing results back
# --------------------------------------------------------------------------

@pytest.fixture
def db_for_writing(tmp_path):
    path = str(tmp_path / "w.db")
    conn = _make_db(path,
                    [("REF0", "master"), ("REF1", "reference"), ("REF2", "reference")],
                    [("REF0", "ACGT"), ("REF1", "ACGT"), ("REF2", "ACGT")])
    return conn


ROWS = [
    {"primary_accession": "REF1", "recomb_node_id": "node_42",
     "breakpoint_1_interval": "3200-3260", "breakpoint_2_interval": "7100-7180",
     "breakpoint_1_start": 3200, "breakpoint_1_end": 3260,
     "breakpoint_2_start": 7100, "breakpoint_2_end": 7180,
     "donor_node_id": "node_7", "acceptor_node_id": "node_19"},
]


class TestStoreResults:
    def test_creates_table_and_status_column(self, db_for_writing):
        SRR.store_results(db_for_writing, ROWS, ["REF0", "REF1", "REF2"])
        names = {r[0] for r in db_for_writing.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        assert SRR.RECOMBINATION_TABLE in names
        cols = {r[1] for r in db_for_writing.execute("PRAGMA table_info(meta_data)")}
        assert SRR.RECOMBINATION_STATUS_COLUMN in cols

    def test_flagged_and_screened_are_distinguished(self, db_for_writing):
        """A sequence that was screened and came back clean is a different
        thing from one that was never screened - the column must say which."""
        SRR.store_results(db_for_writing, ROWS, ["REF0", "REF1", "REF2"])
        status = dict(db_for_writing.execute(
            f"SELECT primary_accession, {SRR.RECOMBINATION_STATUS_COLUMN} FROM meta_data"))
        assert status["REF1"] == SRR.STATUS_RECOMBINANT
        assert status["REF0"] == SRR.STATUS_SCREENED
        assert status["REF2"] == SRR.STATUS_SCREENED

    def test_rerunning_does_not_duplicate(self, db_for_writing):
        SRR.store_results(db_for_writing, ROWS, ["REF0", "REF1", "REF2"])
        SRR.store_results(db_for_writing, ROWS, ["REF0", "REF1", "REF2"])
        n = db_for_writing.execute(
            f"SELECT COUNT(*) FROM {SRR.RECOMBINATION_TABLE}").fetchone()[0]
        assert n == 1, "the unique index should make a re-run idempotent"

    def test_parsed_coordinates_are_stored_alongside_the_raw_interval(self, db_for_writing):
        SRR.store_results(db_for_writing, ROWS, ["REF1"])
        row = db_for_writing.execute(
            f"SELECT breakpoint_1_interval, breakpoint_1_start, breakpoint_1_end "
            f"FROM {SRR.RECOMBINATION_TABLE}").fetchone()
        assert row == ("3200-3260", "3200", "3260")

    def test_no_events_still_marks_everything_screened(self, db_for_writing):
        """An empty result is a real finding and must be recorded as such,
        otherwise a clean screen is indistinguishable from no screen."""
        written, flagged = SRR.store_results(db_for_writing, [], ["REF0", "REF1"])
        assert (written, flagged) == (0, 0)
        status = dict(db_for_writing.execute(
            f"SELECT primary_accession, {SRR.RECOMBINATION_STATUS_COLUMN} FROM meta_data"))
        assert status["REF0"] == SRR.STATUS_SCREENED

    def test_unattributed_events_are_still_stored(self, db_for_writing):
        rows = [dict(ROWS[0], primary_accession=None)]
        written, flagged = SRR.store_results(db_for_writing, rows, ["REF0"])
        assert written == 1 and flagged == 0

    def test_running_twice_is_safe_on_the_status_column(self, db_for_writing):
        SRR.store_results(db_for_writing, ROWS, ["REF0", "REF1"])
        SRR.store_results(db_for_writing, [], ["REF0", "REF1"])
        status = dict(db_for_writing.execute(
            f"SELECT primary_accession, {SRR.RECOMBINATION_STATUS_COLUMN} FROM meta_data"))
        assert status["REF1"] == SRR.STATUS_SCREENED, (
            "a clean re-run must clear a stale recombinant flag")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

class TestCli:
    def test_defaults(self):
        args = SRR.parse_args([])
        assert args.threads == 8
        assert args.num_descendants == 3
        assert args.write_db is False

    def test_write_db_without_db_is_refused(self, tmp_path):
        screen = SRR.ReferenceRecombinationScreen(write_db=True, db=None,
                                                  outdir=str(tmp_path))
        with pytest.raises(ValueError, match="--db"):
            screen.collect(str(tmp_path))

    def test_alignment_without_tree_is_refused(self, tmp_path):
        screen = SRR.ReferenceRecombinationScreen(alignment="a.fa", tree=None,
                                                  outdir=str(tmp_path))
        with pytest.raises(ValueError, match="--alignment and --tree"):
            screen.prepare_inputs()

    def test_missing_binary_names_the_install(self):
        with pytest.raises(FileNotFoundError, match="usher"):
            SRR._require_binary("definitely-not-a-real-binary-xyz")


# --------------------------------------------------------------------------
# The script is wired into the repo
# --------------------------------------------------------------------------

def test_script_is_executable_as_a_module():
    assert (REPO_ROOT / "scripts" / "ScreenReferenceRecombination.py").exists()


def test_documentation_exists():
    """The screen is slow and parameter-sensitive; it is not usable without
    the written guidance on chunking and on why queries need a different
    method entirely."""
    assert (REPO_ROOT / "info_help" / "reference_recombination_screening.md").exists()


# --------------------------------------------------------------------------
# Tripwire against real RIPPLES output
# --------------------------------------------------------------------------

REAL_OUTPUT = REPO_ROOT / "test_data" / "unit" / "reference_recombination"


class TestAgainstRealRipplesOutput:
    """These fixtures were produced by the installed RIPPLES, not hand-written.

    The column names in RECOMBINATION_COLUMNS were read out of the binary. If a
    RIPPLES upgrade renames or reorders a column, the parser would keep running
    and quietly mis-attribute every field. These tests are the tripwire.
    """

    def test_header_matches_the_constant_exactly(self):
        path = REAL_OUTPUT / "real_empty_recombination.tsv"
        header = path.read_text().splitlines()[0].lstrip("#").split("\t")
        assert tuple(h.strip() for h in header) == SRR.RECOMBINATION_COLUMNS

    def test_descendants_header_matches(self):
        path = REAL_OUTPUT / "real_empty_descendants.tsv"
        header = path.read_text().splitlines()[0].lstrip("#").split("\t")
        assert tuple(h.strip() for h in header) == SRR.DESCENDANTS_COLUMNS

    def test_a_real_empty_result_parses_to_nothing(self):
        """RIPPLES writes the header and no rows when it finds nothing. That is
        'screened, no evidence' - it must not raise, and must not invent a row."""
        assert SRR.parse_recombination_tsv(REAL_OUTPUT / "real_empty_recombination.tsv") == []
        assert SRR.parse_descendants_tsv(REAL_OUTPUT / "real_empty_descendants.tsv") == {}

    def test_every_stored_column_is_a_legal_sql_identifier(self):
        """`breakpoint-1_interval` is not usable unquoted; the rename must cover
        every hyphenated name the binary emits."""
        import re as _re
        for column in SRR.RECOMBINATION_COLUMNS:
            stored = SRR.DB_COLUMN_FOR.get(column, column)
            assert _re.match(r'^[A-Za-z_][A-Za-z0-9_]*$', stored), stored
