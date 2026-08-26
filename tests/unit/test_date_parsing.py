"""GenBank collection_date is free text, and every form of it must round-trip.

A single HCV database carries six different date formats:

    empty            49.6%
    DD-Mon-YYYY      28.8%
    year only        17.3%
    ISO yyyy-mm-dd    2.0%
    Mon-YYYY          1.7%
    other             0.5%   (includes ranges such as '2016/2019')
    ISO yyyy-mm       0.1%

The parser inferred which components were present with three regexes, and its
month test matched only alphabetic month names. So ``15-Mar-2020`` kept its day
and month while ``2020-03-15`` was reduced to a year - the ISO form that INSDC is
standardising on was the one losing the most precision. 5,598 rows in that one
database had a full date in the source and a blank ``collection_day`` and
``collection_mon`` in the built database.

Two further inputs produced partial dates that were misparses rather than genuine
partial information: ``8/25/2022`` gave a day with no month, and ``15-Mar-20``
gave a day and month with no year at all.
"""

import pytest

from date_utils import split_date_components


def _components(value):
    result = split_date_components(value)
    return result["day"], result["month"], result["year"]


class TestIsoDatesKeepFullPrecision:
    """The regression that mattered: ISO dates were reduced to a bare year."""

    @pytest.mark.parametrize("value,expected", [
        ("2020-03-15", (15, 3, 2020)),
        ("2018-02-21", (21, 2, 2018)),
        ("2019-01-18", (18, 1, 2019)),
        ("2020-3-5", (5, 3, 2020)),      # sloppy but unambiguous
        ("2020-12-31", (31, 12, 2020)),
        ("2020-01-01", (1, 1, 2020)),
    ])
    def test_full_iso_date(self, value, expected):
        assert _components(value) == expected

    @pytest.mark.parametrize("value,expected", [
        ("2020-03", ("", 3, 2020)),
        ("2020-12", ("", 12, 2020)),
    ])
    def test_iso_year_month(self, value, expected):
        assert _components(value) == expected

    def test_iso_and_legacy_forms_of_the_same_date_agree(self):
        """The whole bug in one assertion."""
        assert _components("2013-04-08") == _components("08-Apr-2013")


class TestLegacyFormsStillWork:
    """These already worked and must keep working."""

    @pytest.mark.parametrize("value,expected", [
        ("15-Mar-2020", (15, 3, 2020)),
        ("08-Apr-2013", (8, 4, 2013)),
        ("Mar-2020", ("", 3, 2020)),
        ("Mar 2020", ("", 3, 2020)),
        ("2020", ("", "", 2020)),
    ])
    def test_form(self, value, expected):
        assert _components(value) == expected


class TestNoComponentIsInvented:
    @pytest.mark.parametrize("value", ["2020", "2013", "1999"])
    def test_a_year_only_date_never_becomes_january_first(self, value):
        """dateutil is given a 1900-01-01 default; if that leaked through, every
        year-only record would claim to be sampled on 1 January."""
        day, month, year = _components(value)
        assert day == "" and month == ""
        assert year == int(value)

    @pytest.mark.parametrize("value", ["", "   ", "unknown", "NA", "circa 2015", "n/a"])
    def test_unparseable_input_yields_nothing(self, value):
        assert _components(value) == ("", "", "")

    @pytest.mark.parametrize("value", ["2016/2019", "1999/2000"])
    def test_a_year_range_is_not_silently_collapsed_to_one_year(self, value):
        """A range spans years; picking one would be a fabrication. Better to
        record nothing than to record a specific year that is possibly wrong."""
        assert _components(value) == ("", "", "")


class TestInvalidCalendarDates:
    """An impossible date should lose the impossible parts, not the whole record."""

    @pytest.mark.parametrize("value,year", [
        ("2020-02-30", 2020),   # February never has 30 days
        ("2020-13-01", 2020),   # no month 13
        ("2020-00-00", 2020),
    ])
    def test_year_survives_but_day_and_month_do_not(self, value, year):
        day, month, got_year = _components(value)
        assert (day, month) == ("", "")
        assert got_year == year

    def test_a_leap_day_in_a_leap_year_is_valid(self):
        assert _components("2020-02-29") == (29, 2, 2020)

    def test_a_leap_day_in_a_non_leap_year_is_not(self):
        day, month, year = _components("2021-02-29")
        assert (day, month) == ("", "")
        assert year == 2021


class TestPartialResultsAreNeverMisleading:
    """A day with no month, or a day and month with no year, is a misparse.

    Downstream code cannot distinguish a genuine partial date from a mangled one,
    so these must not be emitted at all.
    """

    def test_a_day_is_never_returned_without_a_month(self):
        for value in ["8/25/2022", "1/2/2022", "25/12/2022"]:
            day, month, _ = _components(value)
            assert not (day != "" and month == ""), f"{value} gave day={day} with no month"

    def test_a_date_is_never_returned_without_a_year(self):
        for value in ["15-Mar-20", "Mar-20", "15-Mar"]:
            day, month, year = _components(value)
            assert not (year == "" and (day != "" or month != "")), (
                f"{value} gave day={day} month={month} with no year"
            )

    def test_every_returned_component_set_is_a_valid_prefix(self):
        """Only year, year+month, or year+month+day are meaningful."""
        samples = [
            "2020-03-15", "2020-03", "2020", "15-Mar-2020", "Mar-2020", "8/25/2022",
            "15-Mar-20", "2016/2019", "2020-02-30", "unknown", "", "NA", "Mar 2020",
        ]
        for value in samples:
            day, month, year = _components(value)
            present = (year != "", month != "", day != "")
            assert present in {
                (False, False, False),
                (True, False, False),
                (True, True, False),
                (True, True, True),
            }, f"{value} gave an incoherent combination day={day} month={month} year={year}"


class TestParserNeverRaises:
    """GenBankParser calls this on every record; an exception would abort a build."""

    @pytest.mark.parametrize("value", [
        "", "   ", "?", "not a date", "2020-", "-", "//", "9" * 400,
        "2020-03-15T00:00:00Z", "between 2010 and 2012", "0000-00-00",
        "Feb-30-2020", "32-Jan-2020", "2020/03/15",
    ])
    def test_hostile_input(self, value):
        result = split_date_components(value)
        assert set(result) == {"day", "month", "year"}
