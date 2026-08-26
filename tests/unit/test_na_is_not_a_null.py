"""Influenza's neuraminidase is literally named "NA", which pandas nulls by default.

pandas' default NA sentinel list contains the string ``"NA"``. Influenza segment 6
and its gene are both named NA, so every default-configured read erases them. This
was not theoretical: the shipped IAV database stored the neuraminidase gene as

    ('Neuraminidase', None, '', 'whole_genome')

a gene with no name, and the checked-in GISAID fixture has eight rows with a blank
``Segment`` sitting exactly where neuraminidase belongs.

The suite previously codified the bug as intended behaviour - the only three
occurrences of the string "NA" in 14k lines of tests all asserted that it *becomes*
empty or is dropped. Those assertions are correct for the columns they cover
(``host_validated``/``collection_date_validated`` use 'NA' as a deliberate "not
applicable" encoding, 274k cells of it in a real HCV matrix), which is exactly why
the fix has to be per-site rather than a blanket ``keep_default_na=False``.
"""

import csv

import pandas as pd
import pytest

import gisaid_tidy
import segment_utils
from CreateSqliteDB import CreateSqliteDB

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
GENE_INFO = REPO_ROOT / "generic" / "influenza" / "Tables" / "gene_info.tsv"


class TestGeneInfoKeepsTheNaGene:
    def test_the_shipped_flu_gene_table_still_contains_an_NA_gene(self):
        """If this ever changes, the rest of this class is testing nothing."""
        rows = list(csv.DictReader(open(GENE_INFO, encoding="utf-8"), delimiter="\t"))
        names = {r["name"] for r in rows}
        assert "NA" in names, "gene_info.tsv no longer has an NA gene"

    def test_neuraminidase_keeps_its_name_when_read(self):
        frame = CreateSqliteDB._read_tsv_required(
            str(GENE_INFO), [], "gene_info", keep_na_strings=True
        )
        row = frame[frame["description"] == "Neuraminidase"].iloc[0]
        assert row["name"] == "NA"
        assert row["display_name"] == "NA"

    def test_the_default_read_still_erases_it(self):
        """Documents precisely why the opt-in exists."""
        frame = CreateSqliteDB._read_tsv_required(str(GENE_INFO), [], "gene_info")
        row = frame[frame["description"] == "Neuraminidase"].iloc[0]
        assert pd.isna(row["name"])

    def test_no_other_gene_is_affected_by_the_opt_in(self):
        default = CreateSqliteDB._read_tsv_required(str(GENE_INFO), [], "gene_info")
        kept = CreateSqliteDB._read_tsv_required(
            str(GENE_INFO), [], "gene_info", keep_na_strings=True
        )
        assert len(default) == len(kept)
        for column in ("description", "display_name", "name"):
            differing = [
                (a, b) for a, b in zip(default[column], kept[column])
                if not (a == b or (pd.isna(a) and b in ("NA", "NULL", "None", "nan")))
            ]
            assert not differing, f"{column} changed beyond the NA sentinels: {differing}"


class TestGisaidTidyKeepsNaSegments:
    def test_na_segment_survives_the_read(self, tmp_path):
        source = tmp_path / "metadata.tsv"
        source.write_text(
            "Segment_Id\tSegment\tNote\n"
            "EPI1\tHA\tx\n"
            "EPI2\tNA\t\n"      # neuraminidase, plus a genuinely blank cell
            "EPI3\t\tz\n",      # genuinely blank segment
            encoding="utf-8",
        )
        frame = gisaid_tidy._restore_blank_as_missing(
            pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
        )
        assert frame.loc[1, "Segment"] == "NA"

    def test_a_genuinely_blank_cell_is_still_missing(self, tmp_path):
        """`dropna(how='all')` and `pd.notnull(...)` gates must keep working."""
        source = tmp_path / "metadata.tsv"
        source.write_text("Segment_Id\tSegment\nEPI3\t\n", encoding="utf-8")
        frame = gisaid_tidy._restore_blank_as_missing(
            pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
        )
        assert pd.isna(frame.loc[0, "Segment"])
        assert pd.notnull(frame.loc[0, "Segment"]) is False

    def test_an_all_blank_row_is_still_dropped(self, tmp_path):
        source = tmp_path / "metadata.tsv"
        source.write_text("Segment_Id\tSegment\nEPI1\tHA\n\t\n", encoding="utf-8")
        frame = gisaid_tidy._restore_blank_as_missing(
            pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
        )
        assert len(frame.dropna(how="all")) == 1


class TestNaIsNotUniversallyMissing:
    """The judgement is per column, and the shared helper reflects that."""

    def test_na_is_a_segment_name_for_a_name_keyed_virus(self):
        assert segment_utils.is_unavailable("NA", {"ha", "na", "pb1"}) is False

    def test_na_is_missing_for_a_number_keyed_virus(self):
        assert segment_utils.is_unavailable("NA", {"1", "2", "6"}) is True

    @pytest.mark.parametrize("token", ["nan", "none", "<na>", "null"])
    def test_stringified_nulls_are_never_real_labels(self, token):
        assert token in segment_utils.PANDAS_NULL_TOKENS

    def test_na_is_not_in_the_stringified_null_set(self):
        """`NA` is deliberately absent: it is a real influenza name."""
        assert "na" not in segment_utils.PANDAS_NULL_TOKENS
