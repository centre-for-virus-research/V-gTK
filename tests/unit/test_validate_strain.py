import csv
from pathlib import Path

import pytest

from ValidateStrain import (
    clean_strain_name,
    compare_strain_year_to_collection,
    curate_serotype,
    extract_serotype,
    parse_strain_from_definition,
    process_matrix,
    strain_year_from_name,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "validate_strain_edge"


def write_tsv(path: Path, rows, header):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv_dicts(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_parse_and_clean_strain():
    definition = "Influenza A virus (A/chicken/Chile/11(H5N1)) segment"
    parsed = parse_strain_from_definition(definition)
    assert parsed == "A/chicken/Chile/11"

    assert clean_strain_name("A/chicken Chile///11___") == "A/chicken_Chile///11"


def test_extract_serotype_and_curation_edge_cases():
    assert extract_serotype({"serotype": " H3N2 ", "definition": "x"}) == "H3N2"
    assert extract_serotype({"serotype": "", "definition": "virus (H10N7)"}) == "H10N7"
    assert extract_serotype({"definition": "virus without subtype"}) == ""

    replacements = {"HXNY": "unknown"}
    assert curate_serotype("A/H1N1", replacements) == "H1N1"
    assert curate_serotype("h5n1", replacements) == "H5N1"
    assert curate_serotype("H7", replacements) == "H7"
    assert curate_serotype("N2", replacements) == "N2"
    assert curate_serotype("H?N3", replacements) == "N3"
    assert curate_serotype("H9N?", replacements) == "H9"
    assert curate_serotype("HXNY", replacements) == "unknown"
    assert curate_serotype("???", replacements) == "unknown"


def test_process_matrix_end_to_end(tmp_path: Path):
    input_tsv = tmp_path / "gB_matrix.tsv"
    output_tsv = tmp_path / "gB_matrix_validated.tsv"

    write_tsv(
        input_tsv,
        [
            ["ACC1", "A/duck/CA/1", "", "Influenza A virus (A/duck/CA/1(H3N2))"],
            ["ACC2", "", "", "Influenza A virus (A/chicken/PE/99(H5N1))"],
            ["ACC3", "", "A/H1N1", "Influenza A virus (A/human/US/20(H1N1))"],
            ["ACC4", "", "", "no subtype info"],
        ],
        header=["primary_accession", "strain", "serotype", "definition"],
    )

    process_matrix(
        str(input_tsv),
        str(output_tsv),
        overwrite=False,
        serotype_replacements={"UNKNOWNVALUE": "unknown"},
    )

    rows = read_tsv_dicts(output_tsv)

    assert rows[0]["Parsed_strain"] == "A/duck/CA/1"
    assert rows[0]["serotype_validated"] == "H3N2"

    assert rows[1]["Parsed_strain"] == "A/chicken/PE/99"
    assert rows[1]["serotype_validated"] == "H5N1"

    assert rows[2]["Parsed_strain"] == "A/human/US/20"
    assert rows[2]["serotype_validated"] == "H1N1"

    assert rows[3]["Parsed_strain"] == ""
    assert rows[3]["serotype_validated"] == "unknown"


def test_process_matrix_overwrite_mode(tmp_path: Path):
    input_tsv = tmp_path / "gB_matrix.tsv"
    write_tsv(
        input_tsv,
        [["ACCX", "", "", "Influenza A virus (A/goose/CN/7(H9N2))"]],
        header=["primary_accession", "strain", "serotype", "definition"],
    )

    process_matrix(str(input_tsv), str(input_tsv), overwrite=True, serotype_replacements={})
    rows = read_tsv_dicts(input_tsv)

    assert rows[0]["Parsed_strain"] == "A/goose/CN/7"
    assert rows[0]["serotype_validated"] == "H9N2"


def test_edge_dataset_matches_expected_output(tmp_path: Path):
    input_tsv = tmp_path / "input.tsv"
    output_tsv = tmp_path / "output.tsv"

    input_tsv.write_text((DATA_DIR / "input.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    replacements = {}
    with (DATA_DIR / "replacements.tsv").open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            replacements[row["raw_value"]] = row["curated_value"]

    process_matrix(
        str(input_tsv),
        str(output_tsv),
        overwrite=False,
        serotype_replacements=replacements,
    )

    assert read_tsv_dicts(output_tsv) == read_tsv_dicts(DATA_DIR / "expected.tsv")


# ---------------------------------------------------------------------------
# Strain name vs collection date
#
# A flu strain name ends in the collection year, which makes it an independent
# witness to the collection date -- and the date is the field that gets corrupted.
# A rebuilt influenza database carried 277 rows dated 1905: the collection year
# typed into a cell Excel had formatted as a date, so the year came back as a
# 1900-epoch day-serial (serial 2021 is 1905-07-13, serial 2025 is 1905-07-17).
# The strain name was correct in every one of them. Two further records cleared a
# 2010 plausibility floor unnoticed -- a G.1.1 sequence dated 2011 and a J.2
# sequence dated 2012 -- and surfaced as the first sampling date of their lineage
# in a manuscript figure. That is the failure this check exists to prevent.
# ---------------------------------------------------------------------------


THIS_YEAR = 2026


def verdict_for(strain, collection_year):
    row = {"Parsed_strain": strain, "collection_year": collection_year}
    return compare_strain_year_to_collection(row, THIS_YEAR)[0]


@pytest.mark.parametrize("strain,collection_year", [
    ("A/Lisboa/12/2021", 1905),             # Excel serial 2021
    ("A/Baltimore/JH-1150/2025", 1905),     # Excel serial 2025
    ("A/Lisboa/322/2022", 2011),            # the G.1.1 record
    ("A/England/225020837/2022", 2002),
    ("A/_South_Africa/K056330/2023", 2005),
])
def test_a_date_years_from_its_own_strain_name_is_flagged(strain, collection_year):
    assert verdict_for(strain, collection_year) == "mismatch"


@pytest.mark.parametrize("strain,collection_year", [
    ("A/Brisbane/80/2015", 2015),
    ("A/New_Zealand/14CF0730/2014", 2014),
    ("A/Darwin/6/2021", 2022),   # collected December, named for the reporting season
    ("A/Darwin/6/2022", 2021),   # and the same offset the other way
])
def test_an_agreeing_or_season_offset_date_is_not_flagged(strain, collection_year):
    assert verdict_for(strain, collection_year) == "ok"


@pytest.mark.parametrize("strain", [
    "A/duck_Chile//01",                        # '01' is an isolate number
    "A/gull/Peru/88",                          # could be 1988, could be isolate 88
    "A/Anas_platyrhynchos/Camargue/091863/09",  # a real two-digit year, still ambiguous
])
def test_a_short_final_field_is_never_read_as_a_year(strain):
    """Two-digit tails are left unresolved rather than guessed.

    Flu names do carry two-digit years, but nothing in the string separates those
    from an isolate number, and guessing puts healthy records into the conflicts
    report -- the fastest way to get the report ignored. Real corrupt records all
    carried an unambiguous four-digit year, so nothing is lost by declining here.
    """
    assert verdict_for(strain, 2015) == "ambiguous_strain_year"


@pytest.mark.parametrize("strain", [
    "A/chicken/China/15830",   # four digits, but 15830 is not a year
    "A/swine/Iowa/101165",
    "",
])
def test_a_name_with_no_year_is_recorded_as_uncheckable_not_as_a_pass(strain):
    assert verdict_for(strain, 2014) == "no_strain_year"


def test_a_four_digit_isolate_number_before_1918_is_not_a_year():
    assert strain_year_from_name("A/swine/Iowa/1890", THIS_YEAR) is None


def test_a_missing_collection_year_is_reported_separately_from_a_mismatch():
    assert verdict_for("A/Brisbane/1/2010", "") == "no_collection_year"


def test_the_check_never_rewrites_the_collection_date():
    """It flags; it does not correct.

    A mismatch means one of the two fields is wrong and nothing here can tell
    which -- in the 53 rows dated 1905-07-14 the strain name was corrupt too.
    """
    row = {"Parsed_strain": "A/Lisboa/12/2021", "collection_year": 1905,
           "collection_date": "1905-07-13"}
    compare_strain_year_to_collection(row, THIS_YEAR)
    assert row["collection_date"] == "1905-07-13"
    assert row["collection_year"] == 1905
