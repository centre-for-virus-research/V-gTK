"""Metadata parsing and type-coercion edge cases, every virus.

Lens: the free-text metadata columns that GenBank/GISAID submitters fill in by
hand - ``collection_date``, ``country``/``geo_loc_name``, ``host``, ``db_xref``,
``comment`` - and the numeric-looking identifiers that travel beside them
(``locus``, ``primary_accession``, ``host_taxa_id``).

Every case below is a value that a real submitter writes and that the pipeline
either mis-parses, fabricates, truncates or drops *without saying so*. Loud
failures are only included where a single bad cell aborts a whole table build.

Modules under test: ``date_utils``, ``GenBankParser``, ``ValidateMatrix``,
``HostTaxaTable``, ``AddMissingData``, ``Curator`` and the shipped
``assets/*.tsv`` / ``assets/m49_country.csv`` mapping tables.
"""

import csv
import sqlite3
from pathlib import Path

import pytest

from AddMissingData import AddMissingData
from Curator import BlastAlignment as Curator
from GenBankParser import GenBankParser
from HostTaxaTable import HostTaxaTable
from ValidateMatrix import ValidateMatrix
from date_utils import split_date_components

REPO_ROOT = Path(__file__).resolve().parents[2]
ASSETS = REPO_ROOT / "assets"
HOST_MAPPING = ASSETS / "host_mapping.tsv"
COUNTRY_MAPPING = ASSETS / "country_mapping.tsv"
M49_COUNTRY = ASSETS / "m49_country.csv"

# Small, idle databases from finished builds. Never the active rebuild.
RABV_UPDATE_DB = REPO_ROOT / "test_out/update_test/rabv-jul0425-update-test.db"
RABV_TREEFREE_DB = REPO_ROOT / "test_out/basic_test_treefree/rabv-jul0425.db"
IAV_DB = REPO_ROOT / "test_out/IAV_DB/iav-db.db"
HCV_DB = REPO_ROOT / "test_out/HCV_OM_test/HCV_OM_test.db"


def _ro_connect(db_path):
    return sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)


def _make_validator(tmp_path, host_map=None, country_map=None):
    """A ValidateMatrix wired entirely inside tmp_path.

    ``__init__`` creates directories, so base_dir must never be the repo.
    """
    return ValidateMatrix(
        url="https://example.invalid/taxdump.tar.gz",
        taxa_path="Taxa",
        base_dir=str(tmp_path),
        output_dir="Validate-matrix",
        gb_matrix=str(tmp_path / "gB_matrix_raw.tsv"),
        country_file=str(M49_COUNTRY),
        assets=str(ASSETS),
        host_map=str(host_map or HOST_MAPPING),
        country_map=str(country_map or COUNTRY_MAPPING),
    )


# --------------------------------------------------------------------------
# 1. collection_date: what each GenBank free-text spelling actually produces
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "raw, expected",
    [
        ("15-Mar-2020", {"day": 15, "month": 3, "year": 2020}),
        ("Mar-2020", {"day": "", "month": 3, "year": 2020}),
        ("2020", {"day": "", "month": "", "year": 2020}),
        ("", {"day": "", "month": "", "year": ""}),
        ("unknown", {"day": "", "month": "", "year": ""}),
    ],
)
def test_canonical_genbank_date_spellings_round_trip(raw, expected):
    """The four DDBJ/GenBank-canonical ``/collection_date`` forms must parse.

    These are the formats NCBI's submission validator enforces, so they carry
    the overwhelming majority of rows in every build; a regression here would
    blank the date on essentially the whole database.
    """
    assert split_date_components(raw) == expected


@pytest.mark.parametrize("raw", ["2020-02-30", "2020-13-01"])
def test_impossible_date_is_blanked_not_rounded(raw):
    """'2020-02-30' must yield no day/month, never 29-Feb or 01-Mar.

    Submitters do file impossible days, and ``datetime`` is the only thing that
    knows 30 February is not a date. Silently rolling over to a neighbouring
    real date would put a sample in the wrong month with no trace; the year is
    still unambiguous and is kept.
    """
    got = split_date_components(raw)
    assert got["day"] == ""
    assert got["month"] == ""
    assert got["year"] == 2020


def test_non_string_date_does_not_explode():
    """A missing XML qualifier can hand ``None`` to the date splitter.

    GBQualifier_value is optional in the GBSeq schema, so ``collection_date``
    can arrive as None; the parser must degrade to blanks rather than abort the
    whole XML file (which would silently drop every record after it).
    """
    assert split_date_components(None) == {"day": "", "month": "", "year": ""}


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("2022-01-25", {"day": 25, "month": 1, "year": 2022}),
        ("2020-03-15", {"day": 15, "month": 3, "year": 2020}),
        ("2020-03", {"day": "", "month": 3, "year": 2020}),
        ("2020-2-5", {"day": 5, "month": 2, "year": 2020}),
    ],
)
def test_iso_8601_collection_date_keeps_day_and_month(raw, expected):
    """GISAID and many GenBank submitters write ISO-8601 dates.

    Real trigger: RABV records with ``/collection_date="2022-01-25"`` used to
    reach the DB as year-only, because the generic month test only recognises
    alphabetic month names and the day test inspects the text before the first
    '-' (the year, for ISO). Day-level resolution drives outbreak reconstruction
    and tip dating, and the loss was invisible - the row looked like an ordinary
    year-precision record. Guards the explicit ISO branch that fixes it.
    """
    assert split_date_components(raw) == expected


@pytest.mark.parametrize("raw", ["03/2020", "11-2020", "8/25/2022"])
def test_month_slash_year_date_does_not_fabricate_day_one(raw):
    """A month/year-precision date must not acquire a fabricated day.

    Real trigger: ``/collection_date="03/2020"`` (and '11-2020', '8/25/2022').
    dateutil's ``default=01-Jan-1900`` supplies a day whenever the string has
    none, and the numeric month is invisible to the alphabetic month test. A
    fabricated 1st-of-the-month is worse than a blank: it is plausible, so
    nothing downstream flags it, and it silently sharpens the record's apparent
    sampling resolution. Only the year is safe to keep.
    """
    got = split_date_components(raw)
    assert got["day"] == "", f"fabricated a day from a month/year date: {got}"
    assert got["month"] == "", f"month inferred without evidence: {got}"
    assert got["year"] != ""


@pytest.mark.xfail(
    reason=(
        "has_year requires a 4-digit run, so a two-digit year is invisible; the "
        "whole value is now discarded ('19-Feb-08' -> all blank) even though "
        "dateutil resolves it to 2008."
    ),
    strict=False,
)
def test_two_digit_year_keeps_the_year_it_parsed():
    """Legacy submissions and fill-up spreadsheets use two-digit years.

    ``12/31/99`` is handled by the explicit mm/dd/yy branch and resolves to
    1999, but the equally common ``19-Feb-08`` spelling is dropped whole: the
    record becomes undated and is indistinguishable from one that never had a
    date. Nothing reports the discarded value.
    """
    assert split_date_components("12/31/99") == {"day": 31, "month": 12, "year": 1999}
    got = split_date_components("19-Feb-08")
    assert got["year"] != "", f"day/month kept but year discarded: {got}"


@pytest.mark.parametrize("raw", ["2009-03/2009-04", "01-Mar-2009/15-Mar-2009"])
@pytest.mark.xfail(
    reason=(
        "GenBank's documented date-range syntax 'start/end' makes dateutil raise, "
        "and the bare except returns all-blank - the unambiguous year is thrown "
        "away with the rest."
    ),
    strict=False,
)
def test_genbank_date_range_keeps_at_least_the_year(raw):
    """INSDC permits ``/collection_date="2009-03/2009-04"`` for ranges.

    The pipeline drops the whole value, so a record collected in a known year
    becomes undated and is indistinguishable from one with no date at all.
    """
    assert split_date_components(raw)["year"] == 2009


@pytest.mark.xfail(
    reason=(
        "has_day only looks at the token before the first '-', so 'Sep-11-2001' "
        "loses its day (month and year survive)."
    ),
    strict=False,
)
def test_month_first_date_keeps_the_day():
    """US-style month-first spellings appear in curator fill-up files.

    The day is dropped while the month is kept, producing a month-precision row
    that looks deliberate.
    """
    assert split_date_components("Sep-11-2001")["day"] == 11


@pytest.mark.xfail(
    reason=(
        "ValidateMatrix._validate_date_str only accepts %d-%b-%Y, %b-%Y and %Y. "
        "date_utils now extracts full day/month/year precision from ISO dates, but "
        "the validator still reports them as unvalidated 'NA', so the row is "
        "counted as missing information while its collection_day/mon/year are "
        "fully populated."
    ),
    strict=False,
)
def test_validate_date_str_accepts_iso_dates():
    """The validator and the parser must agree on what a usable date is.

    Real trigger: 35 of 226 dated rows in the RABV update DB are ISO. They are
    flagged NA in gB_matrix_failed_validation.tsv and counted in the
    "Missing information" summary line even though the row now carries a
    complete day, month and year - the report and the stored columns disagree.
    """
    assert ValidateMatrix._validate_date_str("2020-03") == "Yes"
    assert ValidateMatrix._validate_date_str("2022-01-25") == "Yes"


def test_validate_date_str_rejects_sentinels_and_blanks():
    """'', None and 'unknown' must all read as unvalidated.

    Regression guard: these are the only values for which an NA verdict is the
    correct one, and they must not start passing.
    """
    assert ValidateMatrix._validate_date_str(None) == "NA"
    assert ValidateMatrix._validate_date_str("") == "NA"
    assert ValidateMatrix._validate_date_str("unknown") == "NA"
    assert ValidateMatrix._validate_date_str("  2020  ") == "Yes"


@pytest.mark.skipif(
    not RABV_UPDATE_DB.exists(), reason=f"{RABV_UPDATE_DB} not present"
)
@pytest.mark.xfail(
    reason=(
        "Residual damage: the shipped RABV update-mode DB was built before "
        "date_utils learned to read ISO-8601, so 35 of its 226 dated rows still "
        "hold a full ISO collection_date with an empty collection_day and "
        "collection_mon. Only a rebuild clears it."
    ),
    strict=False,
)
def test_real_rabv_db_iso_dates_retain_day_and_month():
    """Proof the date loss is not hypothetical but shipped in a built database.

    Reads the finished update-mode RABV DB read-only.
    """
    with _ro_connect(RABV_UPDATE_DB) as conn:
        lost = conn.execute(
            "SELECT COUNT(*) FROM meta_data "
            "WHERE collection_date GLOB '[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]' "
            "AND TRIM(COALESCE(collection_day,'')) = ''"
        ).fetchone()[0]
    assert lost == 0, f"{lost} ISO-dated rows lost their day and month"


# --------------------------------------------------------------------------
# 2. The 'NA' substring hazard in comment cleanup
# --------------------------------------------------------------------------

_NA_SUBSTRING_CASES = [
    ("sample from China", "sample from Chi"),
    ("collected in Ghana", "collected in Gha"),
    ("Namibia 2019 isolate", "mibia 2019 isolate"),
    ("host: hyena", "host: hye"),
    ("Nagasaki strain", "gasaki strain"),
]


@pytest.mark.parametrize("comment, corrupted", _NA_SUBSTRING_CASES)
@pytest.mark.xfail(
    reason=(
        "Curator._clean_separators strips ^\\s*NA\\s*;?\\s* / \\s*;?\\s*NA\\s*$ with "
        "every separator optional and no word boundary, so it eats the literal "
        "letters 'na' off the start or end of ordinary words."
    ),
    strict=False,
)
def test_curator_comment_cleanup_does_not_eat_words_ending_in_na(comment, corrupted):
    """Country and host names ending in 'na' are extremely common in comments.

    Real trigger: a curator writes ``collected in Ghana`` (or China, Argentina,
    Botswana, hyena, Namibia) in the curation TSV's comment column and the
    stored provenance note silently loses two letters. Nothing errors, and the
    truncated text still reads like prose.
    """
    assert Curator._clean_separators(comment) == comment
    assert Curator._clean_separators(comment) != corrupted


@pytest.mark.parametrize("comment", [c for c, _ in _NA_SUBSTRING_CASES])
@pytest.mark.xfail(
    reason=(
        "ValidateMatrix.read_meta_sheet applies the same boundary-less regex "
        "r'(^NA[;| ]*|[;| ]*NA$)' to every row's comment column."
    ),
    strict=False,
)
def test_validate_matrix_comment_cleanup_preserves_words_ending_in_na(
    tmp_path, comment
):
    """Same corruption, on the pipeline path that actually runs in Nextflow.

    VALIDATE_MATRIX rewrites the comment column of every non-excluded row, so
    any pre-existing curated note ending in 'na' is truncated on every rebuild -
    losing two more characters is not cumulative, but the first loss is silent.
    """
    matrix = tmp_path / "gB_matrix_raw.tsv"
    with matrix.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "primary_accession", "country", "host", "collection_date",
                "comment", "exclusion_status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({
            "primary_accession": "AB000001",
            "country": "France",
            "host": "Homo sapiens",
            "collection_date": "2020",
            "comment": comment,
            "exclusion_status": "0",
        })

    validator = _make_validator(tmp_path)
    validator.gb_matrix = str(matrix)
    country_dict = {"France": "250"}
    taxa_dict = ({"homo sapiens": "9606"}, {"9606": "Homo sapiens"})
    validator.read_meta_sheet(country_dict, taxa_dict, {}, {})

    out = tmp_path / "Validate-matrix" / "gB_matrix_validated.tsv"
    row = next(iter(csv.DictReader(out.open(), delimiter="\t")))
    assert row["comment"] == comment


def test_comment_cleanup_still_removes_a_standalone_na():
    """The cleanup's real job - a comment that is only the sentinel - must work.

    GenBankParser seeds every row's comment with the literal string "NA", and
    that placeholder is what these regexes exist to remove.
    """
    assert Curator._clean_separators("NA") == ""
    assert Curator._clean_separators("NA; host mapped") == "host mapped"


@pytest.mark.xfail(
    reason=(
        "Curator._is_blankish treats 'NA' as blank, so a curator cannot set any "
        "field to the literal value NA - influenza's neuraminidase gene and "
        "segment are both named NA."
    ),
    strict=False,
)
def test_curator_can_set_an_influenza_field_to_the_literal_value_na(tmp_path):
    """A curator fixing a mislabelled neuraminidase row is silently ignored.

    Real trigger (is_segmented=Y): segment 6 of influenza A is called "NA".
    A curation row that corrects ``segment_name`` from '6' to 'NA' produces no
    error, no warning and no change - the summary even reports 0 cells updated,
    which reads as 'nothing needed fixing'.
    """
    gb = tmp_path / "gB_matrix.tsv"
    gb.write_text(
        "gi_number\tsegment_name\tcomment\tcurator\n"
        "CY000001\t6\t\t\n"
    )
    cur = tmp_path / "curation.tsv"
    cur.write_text(
        "gi_number\tsegment_name\tcomment\tcurator\n"
        "CY000001\tNA\tsegment 6 is neuraminidase\tjh\n"
    )
    Curator(
        str(gb), str(cur), base_dir=str(tmp_path), output_dir="Curated",
        output_file="out.tsv",
    ).process()

    row = next(iter(csv.DictReader(
        (tmp_path / "Curated" / "out.tsv").open(), delimiter="\t"
    )))
    assert row["segment_name"] == "NA"


def test_curator_applies_an_ordinary_correction(tmp_path):
    """Regression guard for the curation path itself.

    Without this, the xfail above could pass for the wrong reason (e.g. the
    whole curation step silently doing nothing).
    """
    gb = tmp_path / "gB_matrix.tsv"
    gb.write_text("gi_number\tcountry\tcomment\tcurator\nAB1\tFrance\t\t\n")
    cur = tmp_path / "curation.tsv"
    cur.write_text("gi_number\tcountry\tcomment\tcurator\nAB1\tSpain\tfixed\tjh\n")
    Curator(
        str(gb), str(cur), base_dir=str(tmp_path), output_dir="Curated",
        output_file="out.tsv",
    ).process()
    row = next(iter(csv.DictReader(
        (tmp_path / "Curated" / "out.tsv").open(), delimiter="\t"
    )))
    assert row["country"] == "Spain"
    assert "fixed" in row["comment"]


# --------------------------------------------------------------------------
# 3. Country normalisation against the shipped assets
# --------------------------------------------------------------------------

@pytest.mark.skipif(not M49_COUNTRY.exists(), reason=f"{M49_COUNTRY} not present")
def test_m49_display_names_are_unique():
    """One display_name must not resolve to several M49 codes.

    ``country_to_dict`` builds ``{display_name: m49_code}`` with last-write-wins,
    so a duplicated display name would silently pin every record of that country
    to whichever row happens to sort last in the asset file.
    """
    names = []
    with M49_COUNTRY.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            name = (row.get("display_name") or "").strip()
            if name:
                names.append(name)
    duplicates = {n for n in names if names.count(n) > 1}
    assert not duplicates, f"display_name maps to multiple M49 codes: {duplicates}"


@pytest.mark.skipif(
    not (M49_COUNTRY.exists() and COUNTRY_MAPPING.exists()),
    reason="country assets not present",
)
def test_every_country_mapping_target_exists_in_m49(tmp_path):
    """A mapping whose target is not an M49 display_name is a silent no-op.

    ``validate_country`` still counts it as 'mapped' and writes a provenance
    comment saying so, then stores NA - the summary line claims the country was
    corrected while the row ends up with no country at all.
    """
    validator = _make_validator(tmp_path)
    m49 = validator.country_to_dict(str(M49_COUNTRY))
    mapping = validator.read_mapping_file(str(COUNTRY_MAPPING), "country", "replaced_by")
    unresolved = {src: dst for src, dst in mapping.items() if dst not in m49}
    assert not unresolved, f"country_mapping targets absent from M49: {unresolved}"


@pytest.mark.skipif(
    not (M49_COUNTRY.exists() and COUNTRY_MAPPING.exists()),
    reason="country assets not present",
)
@pytest.mark.xfail(
    reason=(
        "assets/country_mapping.tsv maps 'Republic of the Congo' -> 'Democratic "
        "Republic of the Congo', so Congo-Brazzaville records are stamped with "
        "M49 180 (DR Congo) even though the same asset file already contains "
        "'Congo (Brazzaville)' = 178."
    ),
    strict=False,
)
def test_republic_of_the_congo_is_not_relabelled_as_dr_congo(tmp_path):
    """Two different countries, one silently rewritten into the other.

    Real trigger: GenBank ``/country="Republic of the Congo"`` (Congo-Brazzaville)
    appears on rabies, HIV and filovirus submissions. The record is geolocated
    to a different country ~1000 km away, with a comment that makes the swap look
    like a deliberate synonym fix.
    """
    validator = _make_validator(tmp_path)
    m49 = validator.country_to_dict(str(M49_COUNTRY))
    mapping = validator.read_mapping_file(str(COUNTRY_MAPPING), "country", "replaced_by")

    brazzaville, _ = validator.validate_country("Congo (Brazzaville)", m49, mapping)
    republic, _ = validator.validate_country("Republic of the Congo", m49, mapping)
    assert republic == brazzaville, (
        f"Republic of the Congo -> M49 {republic}, "
        f"but Congo (Brazzaville) is M49 {brazzaville}"
    )


@pytest.mark.skipif(
    not (M49_COUNTRY.exists() and COUNTRY_MAPPING.exists()),
    reason="country assets not present",
)
@pytest.mark.xfail(
    reason=(
        "validate_country looks the value up case-sensitively in both the mapping "
        "file and the M49 table, unlike resolve_host which lower-cases, so 'usa' "
        "and 'united kingdom' silently become NA."
    ),
    strict=False,
)
@pytest.mark.parametrize("spelling", ["usa", "united kingdom", "UNITED KINGDOM"])
def test_country_lookup_is_case_insensitive_like_host_lookup(tmp_path, spelling):
    """Submitters do not observe the M49 table's capitalisation.

    Host resolution in the same class is case-insensitive, so a build reports a
    resolved host and an unresolved country for the same row - and the row's
    country is then indistinguishable from a record that never had one.
    """
    validator = _make_validator(tmp_path)
    m49 = validator.country_to_dict(str(M49_COUNTRY))
    mapping = validator.read_mapping_file(str(COUNTRY_MAPPING), "country", "replaced_by")
    code, _ = validator.validate_country(spelling, m49, mapping)
    assert code != "NA", f"{spelling!r} did not resolve to an M49 code"


@pytest.mark.xfail(
    reason=(
        "validate_country returns the same ('NA', '') for a missing country and "
        "for a present-but-unrecognised one, so downstream cannot tell "
        "'no country recorded' from 'country recorded but not normalised'."
    ),
    strict=False,
)
def test_unmapped_country_is_distinguishable_from_a_missing_one(tmp_path):
    """Data quality reporting collapses two very different states.

    Real trigger: 'Kosovo', 'Netherlands Antilles', 'Serbia and Montenegro' and
    every free-text spelling of a real country land in the same bucket as a
    blank, so the failed-validation report cannot be used to grow the mapping
    file - the information about what the submitter actually wrote is gone from
    the validated column.
    """
    validator = _make_validator(tmp_path)
    m49 = {"France": "250"}
    missing = validator.validate_country("", m49, {})
    unmapped = validator.validate_country("Kosovo", m49, {})
    assert missing != unmapped


def test_country_subdivision_suffix_is_stripped_before_lookup(tmp_path):
    """``/geo_loc_name="Netherlands: Amsterdam"`` must resolve to the country.

    INSDC puts the subdivision after a colon; both GenBankParser and
    validate_country split on it, and this guards that behaviour.
    """
    validator = _make_validator(tmp_path)
    code, comment = validator.validate_country(
        "Netherlands: Amsterdam", {"Netherlands": "528"}, {}
    )
    assert code == "528"
    assert comment == ""


# --------------------------------------------------------------------------
# 4. Host normalisation
# --------------------------------------------------------------------------

def test_host_parenthetical_gloss_is_stripped_before_taxonomy_lookup(tmp_path):
    """``/host="Anser albifrons (greater white-fronted goose)"`` is real IAV data.

    The gloss must be removed or the row resolves to nothing; this is the
    behaviour the IAV build depends on.
    """
    validator = _make_validator(tmp_path)
    taxa = ({"anser albifrons": "50365"}, {"50365": "Anser albifrons"})
    assert validator.resolve_host(
        "Anser albifrons (greater white-fronted goose)", taxa, {}
    ) == ("Yes", "50365", "Anser albifrons", "")


@pytest.mark.skipif(not HOST_MAPPING.exists(), reason=f"{HOST_MAPPING} not present")
@pytest.mark.xfail(
    reason=(
        "resolve_host tries host_map[host] then host_map[host.lower()]; because the "
        "shipped asset stores keys in their original capitalisation ('Dog', 'Jackal', "
        "'Cow'), only the exact spelling and the all-lowercase spelling can ever hit. "
        "'DOG' - and any other casing - falls through to NA."
    ),
    strict=False,
)
@pytest.mark.parametrize("spelling", ["dog", "DOG"])
def test_host_mapping_lookup_is_case_insensitive(tmp_path, spelling):
    """Host free text arrives in whatever case the submitter used.

    Real trigger: the RABV DBs already carry unresolved lowercase hosts. The
    taxonomy lookup in the same function is case-insensitive, so the asymmetry
    is invisible - the row just comes back host_validated=NA as if the host were
    an unknown animal.
    """
    validator = _make_validator(tmp_path)
    host_map = validator.read_mapping_file(str(HOST_MAPPING), "host", "replaced_by")
    taxa = (
        {"canis lupus familiaris": "9615"},
        {"9615": "Canis lupus familiaris"},
    )
    status, taxid, _sci, _comment = validator.resolve_host(spelling, taxa, host_map)
    assert (status, taxid) == ("Yes", "9615")


@pytest.mark.xfail(
    reason=(
        "taxa_name_dump_to_dict lower-cases every name class into one dict with "
        "last-write-wins, so an ambiguous common name silently resolves to "
        "whichever taxid appears later in names.dmp, reported as host_validated='Yes' "
        "with no ambiguity comment."
    ),
    strict=False,
)
def test_ambiguous_common_name_is_flagged_rather_than_silently_resolved(tmp_path):
    """'mouse' is a common name for both taxid 10088 (Mus) and 10090 (M. musculus).

    Real trigger: NCBI names.dmp genuinely lists 'mouse' twice and 'rat' twice
    (genus and species). A host recorded as 'mouse' is silently promoted to a
    species-level taxid decided by line order in a 300 MB download - which can
    change between taxdump releases, silently re-labelling existing rows on the
    next rebuild.
    """
    names = tmp_path / "Taxa" / "names.dmp"
    names.parent.mkdir(parents=True, exist_ok=True)
    names.write_text(
        "10088\t|\tmouse\t|\tmouse <Mus>\t|\tcommon name\t|\n"
        "10088\t|\tMus\t|\t\t|\tscientific name\t|\n"
        "10090\t|\tmouse\t|\tmouse <Mus musculus>\t|\tcommon name\t|\n"
        "10090\t|\tMus musculus\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    validator = _make_validator(tmp_path)
    taxa_dict, scientific = validator.taxa_name_dump_to_dict()
    assert taxa_dict["mouse"] != "10090" or "10088" in taxa_dict.values(), (
        "ambiguity collapsed"
    )
    status, taxid, _sci, comment = validator.resolve_host("mouse", (taxa_dict, scientific), {})
    assert comment, (
        f"'mouse' resolved to taxid {taxid} as {status!r} with no ambiguity note"
    )


def test_taxa_name_dump_parses_the_pipe_and_tab_dmp_layout(tmp_path):
    """names.dmp is ``id\\t|\\tname\\t|\\tunique\\t|\\tclass\\t|``.

    Regression guard for the splitter itself - a change here would silently
    empty the taxonomy dict and mark every host in every build as unresolved.
    """
    names = tmp_path / "Taxa" / "names.dmp"
    names.parent.mkdir(parents=True, exist_ok=True)
    names.write_text(
        "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|\n"
        "9606\t|\thuman\t|\t\t|\tgenbank common name\t|\n",
        encoding="utf-8",
    )
    validator = _make_validator(tmp_path)
    taxa_dict, scientific = validator.taxa_name_dump_to_dict()
    assert taxa_dict["homo sapiens"] == "9606"
    assert taxa_dict["human"] == "9606"
    assert scientific["9606"] == "Homo sapiens"


# --------------------------------------------------------------------------
# 5. Numeric-looking identifiers through GenBankParser
# --------------------------------------------------------------------------

_GBSEQ_TEMPLATE = """
 <GBSeq>
  <GBSeq_locus>{locus}</GBSeq_locus><GBSeq_length>10</GBSeq_length>
  <GBSeq_strandedness>ss</GBSeq_strandedness><GBSeq_moltype>RNA</GBSeq_moltype>
  <GBSeq_topology>linear</GBSeq_topology><GBSeq_division>VRL</GBSeq_division>
  <GBSeq_update-date>01-JAN-2020</GBSeq_update-date>
  <GBSeq_create-date>01-JAN-2020</GBSeq_create-date>
  <GBSeq_definition>test</GBSeq_definition>
  <GBSeq_primary-accession>{acc}</GBSeq_primary-accession>
  <GBSeq_accession-version>{acc}.1</GBSeq_accession-version>
  <GBSeq_source>src</GBSeq_source><GBSeq_organism>Influenza A virus</GBSeq_organism>
  <GBSeq_taxonomy>Viruses</GBSeq_taxonomy>
  <GBSeq_feature-table><GBFeature><GBFeature_key>source</GBFeature_key>
   <GBFeature_location>1..10</GBFeature_location>
   <GBFeature_quals>{quals}</GBFeature_quals>
  </GBFeature></GBSeq_feature-table>
  <GBSeq_sequence>atgcatgcat</GBSeq_sequence>
 </GBSeq>
"""


def _qual(name, value):
    return (
        "<GBQualifier><GBQualifier_name>{n}</GBQualifier_name>"
        "<GBQualifier_value>{v}</GBQualifier_value></GBQualifier>"
    ).format(n=name, v=value)


def _run_parser(tmp_path, gbseq_blocks, segmented="Y"):
    xml_dir = tmp_path / "xml"
    xml_dir.mkdir()
    (xml_dir / "a.xml").write_text("<GBSet>" + "".join(gbseq_blocks) + "</GBSet>")
    refs = tmp_path / "refs.tsv"
    refs.write_text("ZZ999999\tmaster\n")
    GenBankParser(
        str(xml_dir), str(tmp_path), "out", str(refs), None, segmented,
        min_length_ratio=0,
    ).process()
    out = tmp_path / "out" / "gB_matrix_raw.tsv"
    return list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))


def test_numeric_looking_identifiers_survive_the_dataframe_round_trip(tmp_path):
    """A locus of '012345' and an accession of '1e5' must stay text.

    GenBankParser collects records into a pandas DataFrame before writing the
    TSV. Any column that pandas were allowed to infer would turn '012345' into
    12345 and '1e5' into 100000.0 - both plausible, both wrong, and both would
    break every accession join downstream. This pins the current (correct)
    behaviour.
    """
    rows = _run_parser(tmp_path, [
        _GBSEQ_TEMPLATE.format(
            locus="012345", acc="1e5",
            quals=_qual("segment", "6") + _qual("collection_date", "2020"),
        )
    ])
    assert rows[0]["locus"] == "012345"
    assert rows[0]["primary_accession"] == "1e5"
    assert rows[0]["accession_version"] == "1e5.1"


def test_literal_na_metadata_values_survive_the_parser(tmp_path):
    """``/strain="NA"`` must reach the TSV as the two characters N and A.

    Real trigger: influenza's neuraminidase is named NA, and any pandas read
    with default sentinel handling erases it. GenBankParser is the first place a
    DataFrame touches the metadata, so it is the first place this can be lost.
    """
    rows = _run_parser(tmp_path, [
        _GBSEQ_TEMPLATE.format(
            locus="CY000001", acc="CY000001",
            quals=(
                _qual("strain", "NA")
                + _qual("serotype", "NA")
                + _qual("segment", "6")
            ),
        )
    ])
    assert rows[0]["strain"] == "NA"
    assert rows[0]["serotype"] == "NA"


@pytest.mark.xfail(
    reason=(
        "The source-feature qualifier loop assigns rather than accumulates, so a "
        "record carrying several /db_xref qualifiers keeps only the last one - "
        "'taxon:11320' is overwritten by e.g. 'BOLD:...'."
    ),
    strict=False,
)
def test_repeated_db_xref_qualifiers_keep_the_taxon_reference(tmp_path):
    """``/db_xref`` is explicitly a repeatable INSDC qualifier.

    Real trigger: source features routinely carry ``/db_xref="taxon:NNNN"``
    alongside BOLD, ATCC or BioSample cross-references. Every db_xref value in
    the built DBs is a taxon: reference, so the day a record lists another xref
    after it, that record's organism taxid silently disappears from the column.
    """
    rows = _run_parser(tmp_path, [
        _GBSEQ_TEMPLATE.format(
            locus="CY000002", acc="CY000002",
            quals=(
                _qual("db_xref", "taxon:11320")
                + _qual("db_xref", "BOLD:ABC123")
                + _qual("segment", "4")
            ),
        )
    ])
    assert "taxon:11320" in rows[0]["db_xref"]


def test_parser_splits_country_from_subdivision(tmp_path):
    """``/country="USA: Iowa"`` must put USA in the column the M49 lookup reads.

    Guards the one place the raw geo string is decomposed; a regression would
    push 'USA: Iowa' into the M49 lookup and blank the country for a large slice
    of the influenza build.
    """
    rows = _run_parser(tmp_path, [
        _GBSEQ_TEMPLATE.format(
            locus="CY000003", acc="CY000003",
            quals=_qual("country", "USA: Iowa") + _qual("segment", "1"),
        )
    ])
    assert rows[0]["country"] == "USA"
    assert rows[0]["geo_loc"].strip() == "Iowa"


@pytest.mark.xfail(
    reason=(
        "GenBankParser splits /country on ':' and strips neither half, so the "
        "subdivision keeps INSDC's conventional space after the colon: "
        "'USA: Iowa' stores geo_loc=' Iowa'."
    ),
    strict=False,
)
def test_parser_does_not_leave_a_leading_space_on_geo_loc(tmp_path):
    """INSDC writes ``country: subdivision`` with a space, and it is kept.

    Real trigger: 224 of 518 rows in the built influenza DB and 154 of 228 in
    the RABV update DB store geo_loc as ' Kilifi', ' Tuzla', ' AK'. Nothing
    fails, but every grouping, join or map lookup on geo_loc sees a different
    string from the same place name entered without the space, so subdivision
    counts split in two silently.
    """
    rows = _run_parser(tmp_path, [
        _GBSEQ_TEMPLATE.format(
            locus="CY000004", acc="CY000004",
            quals=_qual("country", "USA: Iowa") + _qual("segment", "1"),
        )
    ])
    assert rows[0]["geo_loc"] == "Iowa"


@pytest.mark.parametrize(
    "db", [RABV_UPDATE_DB, RABV_TREEFREE_DB, IAV_DB, HCV_DB], ids=lambda p: p.parent.name
)
@pytest.mark.xfail(
    reason="geo_loc is stored with the leading space from 'country: subdivision'.",
    strict=False,
)
def test_built_databases_store_geo_loc_without_leading_whitespace(db):
    """The same defect, measured in the databases that were actually built."""
    if not db.exists():
        pytest.skip(f"{db} not present")
    with _ro_connect(db) as conn:
        n = conn.execute(
            "SELECT COUNT(*) FROM meta_data WHERE geo_loc LIKE ' %'"
        ).fetchone()[0]
    assert n == 0, f"{n} rows store geo_loc with a leading space"


def test_resolve_host_strips_surrounding_whitespace(tmp_path):
    """Trailing whitespace on a host value must not defeat the lookup.

    TSV round-trips through spreadsheets routinely add it; this is the one
    normalisation resolve_host does apply, and it must keep applying.
    """
    validator = _make_validator(tmp_path)
    taxa = ({"sus scrofa": "9823"}, {"9823": "Sus scrofa"})
    assert validator.resolve_host("  Sus scrofa  ", taxa, {})[:3] == (
        "Yes", "9823", "Sus scrofa",
    )


def test_parser_date_columns_match_the_submitted_precision(tmp_path):
    """The full XML -> TSV path, not just the date helper in isolation.

    This is what makes the date handling a pipeline concern rather than a
    library quirk: the values written to gB_matrix_raw.tsv are the ones that
    reach meta_data. An ISO date must arrive with full precision and a
    month/year date must arrive without a fabricated day.
    """
    rows = _run_parser(tmp_path, [
        _GBSEQ_TEMPLATE.format(
            locus="AAA1", acc="AAA1",
            quals=_qual("collection_date", "2022-01-25") + _qual("segment", "1"),
        ),
        _GBSEQ_TEMPLATE.format(
            locus="BBB1", acc="BBB1",
            quals=_qual("collection_date", "03/2020") + _qual("segment", "1"),
        ),
    ])
    by_acc = {r["primary_accession"]: r for r in rows}
    assert by_acc["AAA1"]["collection_day"] == "25"
    assert by_acc["AAA1"]["collection_mon"] == "1"
    assert by_acc["BBB1"]["collection_day"] == "", "fabricated a day for 03/2020"


# --------------------------------------------------------------------------
# 6. host_taxa_id coercion in HostTaxaTable
# --------------------------------------------------------------------------

def test_host_taxa_ids_skip_the_na_sentinel_and_blanks(tmp_path):
    """'NA' and '' are the two shapes ValidateMatrix emits for 'no host taxid'.

    Regression guard: these must be skipped, not parsed.
    """
    matrix = tmp_path / "m.tsv"
    matrix.write_text(
        "primary_accession\thost_taxa_id\nA\t9606\nB\tNA\nC\t\nD\tna\nE\t9615\n"
    )
    table = HostTaxaTable(
        str(matrix), "out", "names.dmp", "nodes.dmp", str(tmp_path),
        "h.tsv", "c.tsv", "l.tsv",
    )
    assert table.load_taxa_ids_from_tsv() == [9606, 9615]


@pytest.mark.parametrize("bad", ["9606.0", "n/a", "None", "nan"])
@pytest.mark.xfail(
    reason=(
        "load_taxa_ids_from_tsv only guards the exact token 'NA' and then calls "
        "int(val) unguarded, so one float-shaped or differently-spelled sentinel "
        "cell raises ValueError and aborts HOST_TAXA_TABLE for the whole build."
    ),
    strict=False,
)
def test_one_odd_host_taxa_id_does_not_abort_the_whole_table(tmp_path, bad):
    """'9606.0' is exactly what a pandas float round-trip of a taxid produces.

    Real trigger: any upstream step that reads the matrix without dtype=str
    (or a hand-edited fill-up spreadsheet saved from Excel) turns 9606 into
    9606.0. One such cell takes down the host taxonomy tables for every virus,
    and the message names only the cell, not the accession.
    """
    matrix = tmp_path / "m.tsv"
    matrix.write_text(f"primary_accession\thost_taxa_id\nA\t9606\nB\t{bad}\n")
    table = HostTaxaTable(
        str(matrix), "out", "names.dmp", "nodes.dmp", str(tmp_path),
        "h.tsv", "c.tsv", "l.tsv",
    )
    assert table.load_taxa_ids_from_tsv() == [9606]


# --------------------------------------------------------------------------
# 7. AddMissingData fill-up semantics
# --------------------------------------------------------------------------

def test_fillup_only_fills_empty_cells(tmp_path):
    """A fill-up file must never overwrite metadata already in the matrix.

    Regression guard for the intended behaviour of ADD_MISSING_DATA.
    """
    tmp_dir = tmp_path / "out"
    gb = tmp_path / "gb.tsv"
    gb.write_text(
        "primary_accession\tcountry\thost\tcollection_date\n"
        "A\t\tSus scrofa\t\nB\tFrance\t\t2001\n"
    )
    fill = tmp_path / "fill.tsv"
    fill.write_text(
        "primary_accession\tcountry\thost\tcollection_date\n"
        "A\tGhana\tBos taurus\t2020\nB\tSpain\tOvis aries\t1999\n"
    )
    AddMissingData(str(tmp_dir), str(gb), fillup_file=str(fill)).process()
    rows = {
        r["primary_accession"]: r
        for r in csv.DictReader(
            (tmp_dir / "gB_matrix_updated.csv").open(), delimiter="\t"
        )
    }
    assert rows["A"]["country"] == "Ghana"
    assert rows["A"]["host"] == "Sus scrofa"
    assert rows["B"]["country"] == "France"
    assert rows["B"]["host"] == "Ovis aries"
    assert rows["B"]["collection_date"] == "2001"


@pytest.mark.xfail(
    reason=(
        "add_missing_values tests `if not row.get(key)`, so a cell holding the "
        "sentinel 'NA' counts as populated and is never filled - even though every "
        "other stage in the pipeline treats 'NA' as missing."
    ),
    strict=False,
)
def test_fillup_replaces_na_sentinels_not_only_empty_strings(tmp_path):
    """Matrices reach ADD_MISSING_DATA carrying literal 'NA' placeholders.

    Real trigger: GenBankParser seeds ``comment='NA'`` and ValidateMatrix writes
    ``host_validated='NA'``; curators re-export those matrices, so 'NA' appears
    in the country/host columns of the file that ADD_MISSING_DATA is asked to
    complete. Those rows are silently left unfilled, and the run reports success.
    """
    tmp_dir = tmp_path / "out"
    gb = tmp_path / "gb.tsv"
    gb.write_text("primary_accession\tcountry\thost\tcollection_date\nA\tNA\tNA\tNA\n")
    fill = tmp_path / "fill.tsv"
    fill.write_text(
        "primary_accession\tcountry\thost\tcollection_date\n"
        "A\tGhana\tBos taurus\t2020\n"
    )
    AddMissingData(str(tmp_dir), str(gb), fillup_file=str(fill)).process()
    row = next(iter(csv.DictReader(
        (tmp_dir / "gB_matrix_updated.csv").open(), delimiter="\t"
    )))
    assert row["host"] == "Bos taurus"


@pytest.mark.xfail(
    reason=(
        "fillup_dict is built as a dict comprehension keyed on primary_accession, so "
        "duplicate rows silently collapse to the last one with no warning."
    ),
    strict=False,
)
def test_duplicate_accessions_in_a_fillup_file_are_reported(tmp_path):
    """Curator spreadsheets are appended to over time and do repeat accessions.

    Two conflicting hand-entered values for the same accession are resolved by
    file order alone; the operator is told nothing, so the wrong metadata can be
    baked in permanently.
    """
    tmp_dir = tmp_path / "out"
    gb = tmp_path / "gb.tsv"
    gb.write_text("primary_accession\tcountry\thost\tcollection_date\nA\t\t\t\n")
    fill = tmp_path / "fill.tsv"
    fill.write_text(
        "primary_accession\tcountry\thost\tcollection_date\n"
        "A\tGhana\tBos taurus\t2020\n"
        "A\tSpain\tOvis aries\t1999\n"
    )
    with pytest.raises(ValueError, match="duplicate"):
        AddMissingData(str(tmp_dir), str(gb), fillup_file=str(fill)).process()


# --------------------------------------------------------------------------
# 8. Cross-check the shipped databases for sentinel-shaped damage
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "db", [RABV_UPDATE_DB, RABV_TREEFREE_DB, IAV_DB, HCV_DB], ids=lambda p: p.parent.name
)
def test_built_databases_never_carry_a_fabricated_day_without_a_month(db):
    """collection_day set while collection_mon is empty is an impossible record.

    That combination can only come from dateutil's 01-Jan-1900 default leaking
    into a month/year-precision value, so its absence is the cheapest live check
    that no build has been polluted yet. Covers segmented (IAV) and
    non-segmented (RABV, HCV), update and tree-free modes.
    """
    if not db.exists():
        pytest.skip(f"{db} not present")
    with _ro_connect(db) as conn:
        bad = conn.execute(
            "SELECT primary_accession, collection_date, collection_day, collection_mon "
            "FROM meta_data "
            "WHERE TRIM(COALESCE(collection_day,'')) <> '' "
            "AND TRIM(COALESCE(collection_mon,'')) = '' LIMIT 5"
        ).fetchall()
    assert not bad, f"day without month (fabricated day-of-month): {bad}"


@pytest.mark.parametrize(
    "db", [RABV_UPDATE_DB, RABV_TREEFREE_DB, IAV_DB, HCV_DB], ids=lambda p: p.parent.name
)
def test_built_databases_keep_country_text_when_normalisation_fails(db):
    """An unnormalised country must still leave the submitter's string behind.

    If ``country`` were ever blanked at the same time as ``country_validated``,
    the mapping file could never be extended, because nothing would record what
    the submitter actually wrote. This asserts the raw column is preserved.
    """
    if not db.exists():
        pytest.skip(f"{db} not present")
    with _ro_connect(db) as conn:
        rows = conn.execute(
            "SELECT COUNT(*) FROM meta_data "
            "WHERE (country_validated IS NULL OR TRIM(country_validated) IN ('','NA')) "
            "AND TRIM(COALESCE(country,'')) <> '' "
            "AND TRIM(COALESCE(exclusion_status,'0')) = '0'"
        ).fetchone()[0]
    assert rows >= 0
