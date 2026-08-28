"""Tests for the chunked recombination driver and its report.

Nothing here launches UShER or RIPPLES. What this module adds over
``ScreenReferenceRecombination`` is planning, merging and reporting, and all
three are pure enough to test directly. The one behaviour worth calling out is
the report's refusal to describe an unscreened database as clean - that is the
single failure mode that would actively mislead a reader, so it is asserted
rather than assumed.
"""

import sqlite3
from pathlib import Path

import pytest

import recombination_hunter_db as HUNT


# --------------------------------------------------------------------------
# Thread budget - the ceiling covers the whole fan-out, not each chunk
# --------------------------------------------------------------------------

class TestThreadBudget:
    def test_ceiling_is_inherited_from_the_screen(self):
        import ScreenReferenceRecombination as SRR
        assert HUNT.MAX_THREADS == SRR.MAX_THREADS == 16

    @pytest.mark.parametrize("chunks,threads", [
        (8, 2), (4, 4), (16, 1), (1, 16), (2, 8), (3, 5),
    ])
    def test_never_exceeds_the_ceiling(self, chunks, threads):
        got_chunks, got_threads = HUNT.plan_thread_budget(chunks, threads)
        assert got_chunks * got_threads <= HUNT.MAX_THREADS

    def test_default_saturates_the_ceiling(self):
        chunks, threads = HUNT.plan_thread_budget(
            HUNT.DEFAULT_CHUNKS, HUNT.DEFAULT_THREADS_PER_CHUNK
        )
        assert (chunks, threads) == (8, 2)
        assert chunks * threads == HUNT.MAX_THREADS

    def test_threads_are_surrendered_before_chunks(self):
        # Concurrency across branches is what buys wall-clock, so an
        # over-subscribed request should keep its chunks and lose its threads.
        assert HUNT.plan_thread_budget(8, 8) == (8, 2)
        assert HUNT.plan_thread_budget(4, 16) == (4, 4)

    def test_more_chunks_than_the_ceiling_are_dropped(self):
        assert HUNT.plan_thread_budget(32, 4) == (16, 1)

    @pytest.mark.parametrize("chunks,threads", [
        (0, 0), (-1, -1), (None, None),
    ])
    def test_degenerate_requests_become_one_by_one(self, chunks, threads):
        assert HUNT.plan_thread_budget(chunks, threads) == (1, 1)


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------

class TestParsimonyGain:
    def test_gain_is_original_minus_recombinant(self):
        assert HUNT.parsimony_gain(
            {"original_parsimony": "40", "recomb_parsimony": "22"}
        ) == 18

    @pytest.mark.parametrize("row", [
        {"original_parsimony": None, "recomb_parsimony": "22"},
        {"original_parsimony": "40", "recomb_parsimony": None},
        {"original_parsimony": "", "recomb_parsimony": ""},
        {"original_parsimony": "n/a", "recomb_parsimony": "22"},
        {},
    ])
    def test_unscored_events_are_none_not_zero(self, row):
        # Zero would sort as a real-but-weak call; None keeps it out of the
        # ranking instead of quietly claiming no improvement.
        assert HUNT.parsimony_gain(row) is None


class TestBreakpointHistogram:
    def test_no_intervals_gives_no_histogram(self):
        assert HUNT.breakpoint_histogram([]) == []

    def test_unparseable_intervals_are_skipped(self):
        rows = [{"breakpoint_1_start": None, "breakpoint_1_end": None,
                 "breakpoint_2_start": None, "breakpoint_2_end": None}]
        assert HUNT.breakpoint_histogram(rows) == []

    def test_every_breakpoint_is_counted_once(self):
        rows = [{
            "breakpoint_1_start": 100, "breakpoint_1_end": 200,
            "breakpoint_2_start": 5000, "breakpoint_2_end": 5100,
        }]
        lines = HUNT.breakpoint_histogram(rows)
        assert sum(int(line.split()[-1]) for line in lines) == 2

    def test_identical_coordinates_do_not_divide_by_zero(self):
        rows = [{
            "breakpoint_1_start": 500, "breakpoint_1_end": 500,
            "breakpoint_2_start": 500, "breakpoint_2_end": 500,
        }]
        assert len(HUNT.breakpoint_histogram(rows)) == 1


# --------------------------------------------------------------------------
# Merging chunk output
# --------------------------------------------------------------------------

RECOMB_HEADER = "\t".join([
    "recomb_node_id", "breakpoint-1_interval", "breakpoint-2_interval",
    "donor_node_id", "donor_is_sibling", "donor_parsimony",
    "acceptor_node_id", "acceptor_is_sibling", "acceptor_parsimony",
    "original_parsimony", "min_starting_parsimony", "recomb_parsimony",
])


def _write_chunk(root: Path, name: str, events, descendants):
    results = root / name / "ripples"
    results.mkdir(parents=True)
    lines = [RECOMB_HEADER]
    for node, bp1, bp2 in events:
        lines.append("\t".join([node, bp1, bp2, "dn", "0", "5", "an", "0", "6",
                                "40", "38", "22"]))
    (results / "recombination.tsv").write_text("\n".join(lines) + "\n")
    (results / "descendants.tsv").write_text(
        "node_id\tdescendants\n"
        + "".join(f"{n}\t{','.join(s)}\n" for n, s in descendants.items())
    )
    return root / name


class TestMergeChunks:
    def test_events_from_every_chunk_are_kept(self, tmp_path):
        a = _write_chunk(tmp_path, "chunk_000", [("n1", "10-20", "30-40")], {"n1": ["A"]})
        b = _write_chunk(tmp_path, "chunk_001", [("n2", "50-60", "70-80")], {"n2": ["B"]})
        events, descendants = HUNT.merge_chunks([a, b])
        assert {e["recomb_node_id"] for e in events} == {"n1", "n2"}
        assert descendants == {"n1": ["A"], "n2": ["B"]}

    def test_the_same_event_seen_twice_is_counted_once(self, tmp_path):
        # Ranges are disjoint so this should not happen, but re-running a single
        # chunk into an existing directory would otherwise inflate the headline.
        a = _write_chunk(tmp_path, "chunk_000", [("n1", "10-20", "30-40")], {"n1": ["A"]})
        b = _write_chunk(tmp_path, "chunk_001", [("n1", "10-20", "30-40")], {"n1": ["A"]})
        events, descendants = HUNT.merge_chunks([a, b])
        assert len(events) == 1
        assert descendants["n1"] == ["A"]

    def test_a_chunk_that_produced_nothing_is_not_an_error(self, tmp_path):
        a = _write_chunk(tmp_path, "chunk_000", [("n1", "10-20", "30-40")], {"n1": ["A"]})
        empty = tmp_path / "chunk_001" / "ripples"
        empty.mkdir(parents=True)
        events, _ = HUNT.merge_chunks([a, tmp_path / "chunk_001"])
        assert len(events) == 1

    def test_intervals_are_parsed_into_integers(self, tmp_path):
        a = _write_chunk(tmp_path, "chunk_000", [("n1", "3200-3260", "5100-5180")], {})
        events, _ = HUNT.merge_chunks([a])
        assert events[0]["breakpoint_1_start"] == 3200
        assert events[0]["breakpoint_2_end"] == 5180


# --------------------------------------------------------------------------
# Database reads
# --------------------------------------------------------------------------

def _make_db(path: Path, status_column=True, statuses=()):
    conn = sqlite3.connect(str(path))
    columns = "primary_accession TEXT, accession_type TEXT"
    if status_column:
        columns += ", recombination_status TEXT"
    conn.execute(f"CREATE TABLE meta_data ({columns})")
    rows = [("REF1", "master"), ("REF2", "reference"), ("Q1", "query")]
    if status_column:
        lookup = dict(statuses)
        conn.executemany(
            "INSERT INTO meta_data VALUES (?,?,?)",
            [(a, t, lookup.get(a)) for a, t in rows],
        )
    else:
        conn.executemany("INSERT INTO meta_data VALUES (?,?)", rows)
    conn.commit()
    conn.close()
    return path


class TestDatabaseReads:
    def test_no_status_column_means_never_screened(self, tmp_path):
        db = _make_db(tmp_path / "a.db", status_column=False)
        conn = sqlite3.connect(str(db))
        assert HUNT.count_screened(conn) == 0
        conn.close()

    def test_a_column_of_nulls_still_means_never_screened(self, tmp_path):
        db = _make_db(tmp_path / "b.db")
        conn = sqlite3.connect(str(db))
        assert HUNT.count_screened(conn) == 0
        conn.close()

    def test_queries_do_not_count_towards_screen_coverage(self, tmp_path):
        db = _make_db(tmp_path / "c.db",
                      statuses=[("REF1", "screened_no_evidence"),
                                ("Q1", "screened_no_evidence")])
        conn = sqlite3.connect(str(db))
        assert HUNT.count_screened(conn) == 1
        conn.close()

    def test_missing_results_table_reads_as_no_events(self, tmp_path):
        db = _make_db(tmp_path / "d.db")
        conn = sqlite3.connect(str(db))
        assert HUNT.read_events_from_db(conn) == []
        conn.close()

    def test_status_counts_separate_null_from_clean(self, tmp_path):
        db = _make_db(tmp_path / "e.db",
                      statuses=[("REF1", "screened_no_evidence")])
        conn = sqlite3.connect(str(db))
        counts = dict(HUNT.gather_status_counts(conn))
        conn.close()
        assert counts["screened_no_evidence"] == 1
        assert counts["(not screened)"] == 1


# --------------------------------------------------------------------------
# The report
# --------------------------------------------------------------------------

FINDING = {
    "primary_accession": "EF407452", "recomb_node_id": "node_12",
    "breakpoint_1_interval": "3200-3260", "breakpoint_2_interval": "5100-5180",
    "breakpoint_1_start": 3200, "breakpoint_1_end": 3260,
    "breakpoint_2_start": 5100, "breakpoint_2_end": 5180,
    "donor_node_id": "d", "acceptor_node_id": "a",
    "original_parsimony": "40", "recomb_parsimony": "22",
}


class TestReport:
    def test_unscreened_is_never_described_as_clean(self):
        text = HUNT.build_report("x.db", [], ["A"] * 238, {}, None, {}, screened=False)
        assert "Not screened" in text
        assert "No recombination detected" not in text
        assert "not* a clean result" in text or "not a clean result" in text.replace("*", "")

    def test_screened_and_empty_reports_a_clean_result(self):
        text = HUNT.build_report("x.db", [], ["A"] * 238, {}, None, {}, screened=True)
        assert "No recombination detected" in text
        assert "238 reference sequences" in text

    def test_a_clean_report_does_not_lecture_about_acting_on_findings(self):
        text = HUNT.build_report("x.db", [], ["A"] * 2, {}, None, {}, screened=True)
        assert "Before acting on this" not in text

    def test_findings_appear_with_their_context_and_gain(self):
        context = {"EF407452": {"nearest_reference_genotype": "2",
                                "nearest_reference_subtype": "2k"}}
        text = HUNT.build_report("x.db", [FINDING], ["A"] * 238, context, None, {})
        assert "EF407452" in text
        assert "2k" in text
        assert "| 18 |" in text          # 40 - 22
        assert "Before acting on this" in text

    def test_events_are_ranked_by_gain(self):
        weak = dict(FINDING, primary_accession="WEAK",
                    original_parsimony="30", recomb_parsimony="28")
        text = HUNT.build_report("x.db", [weak, FINDING], ["A"] * 238, {}, None, {})
        assert text.index("EF407452") < text.index("WEAK")

    def test_unattributed_events_are_reported_not_dropped(self):
        orphan = dict(FINDING, primary_accession=None)
        text = HUNT.build_report("x.db", [orphan], ["A"] * 238, {}, None, {})
        assert "no listed descendants" in text
        assert "node_12" in text

    def test_run_settings_are_recorded(self):
        meta = {"branch_length (-l)": 3, "chunks x threads": "8 x 2"}
        text = HUNT.build_report("x.db", [FINDING], ["A"], {}, None, meta)
        assert "| branch_length (-l) | 3 |" in text
        assert "| chunks x threads | 8 x 2 |" in text


class TestEventsTsv:
    def test_gain_is_written_as_its_own_column(self, tmp_path):
        out = tmp_path / "events.tsv"
        HUNT.write_events_tsv(out, [FINDING])
        header, row = out.read_text().splitlines()
        assert "parsimony_gain" in header.split("\t")
        assert row.split("\t")[header.split("\t").index("parsimony_gain")] == "18"

    def test_missing_fields_are_empty_not_the_string_none(self, tmp_path):
        out = tmp_path / "events.tsv"
        HUNT.write_events_tsv(out, [{"primary_accession": "A"}])
        row = out.read_text().splitlines()[1]
        assert "None" not in row


# --------------------------------------------------------------------------
# The MAT must exist before RIPPLES is allowed near it
#
# UShER resolves -o against the process working directory while -d only moves
# its auxiliary output, so a MAT asked for inside --outdir is written next to
# wherever the command was launched. RIPPLES then reads a file that is not
# there, warns "Tree found empty", finds zero long branches and exits 0 - a
# clean bill of health for a screen that never ran. Both halves are tested:
# the path is absolute, and a missing result is fatal.
# --------------------------------------------------------------------------

class TestMatBuild:
    def _hunter(self, tmp_path, **kwargs):
        return HUNT.RecombinationHunter(db=None, outdir=str(tmp_path), **kwargs)

    def test_output_path_is_absolute(self, tmp_path, monkeypatch):
        seen = {}

        def fake_run(command, **_):
            if command[0] == "usher":
                seen["mat"] = command[command.index("-o") + 1]
                Path(seen["mat"]).write_bytes(b"protobuf")
            return None

        monkeypatch.setattr(HUNT.screen, "_require_binary", lambda name: name)
        monkeypatch.setattr(HUNT.subprocess, "run", fake_run)

        hunter = self._hunter(tmp_path)
        base = hunter._base_screen()
        mat = hunter.build_mat(base, "in.fa", "in.nwk")

        assert Path(seen["mat"]).is_absolute()
        assert Path(mat).is_absolute()
        assert Path(mat).parent == tmp_path

    def test_a_missing_mat_is_fatal_not_an_empty_screen(self, tmp_path, monkeypatch):
        # UShER "succeeds" but writes nothing where we asked.
        monkeypatch.setattr(HUNT.screen, "_require_binary", lambda name: name)
        monkeypatch.setattr(HUNT.subprocess, "run", lambda *a, **k: None)

        hunter = self._hunter(tmp_path)
        base = hunter._base_screen()
        with pytest.raises(SystemExit) as excinfo:
            hunter.build_mat(base, "in.fa", "in.nwk")
        assert "missing tree" in str(excinfo.value)

    def test_a_zero_byte_mat_is_fatal(self, tmp_path, monkeypatch):
        def fake_run(command, **_):
            if command[0] == "usher":
                Path(command[command.index("-o") + 1]).write_bytes(b"")

        monkeypatch.setattr(HUNT.screen, "_require_binary", lambda name: name)
        monkeypatch.setattr(HUNT.subprocess, "run", fake_run)

        hunter = self._hunter(tmp_path)
        with pytest.raises(SystemExit):
            hunter.build_mat(hunter._base_screen(), "in.fa", "in.nwk")


class TestNothingTestedIsNotCleanEither:
    def test_zero_branches_is_reported_as_untested(self):
        text = HUNT.build_report("x.db", [], ["A"] * 238, {}, None, {},
                                 screened=True, branches_tested=0)
        assert "Nothing was tested" in text
        assert "No recombination detected" not in text

    def test_a_real_clean_screen_says_how_much_was_tested(self):
        text = HUNT.build_report("x.db", [], ["A"] * 238, {}, None, {},
                                 screened=True, branches_tested=156)
        assert "No recombination detected" in text
        assert "156 long branch" in text
