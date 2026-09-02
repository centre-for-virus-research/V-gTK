"""Coverage for scripts/accession_utils.py - the single authority on accession identity.

The pipeline deliberately carries two spellings of every GenBank name:

* ``primary_accession`` / ``locus`` - BARE (``PV547761``). The identity key that
  meta_data, sequence_alignment, features, the FASTA headers, the newick tip
  labels and the on-disk filenames all join on.
* ``accession_version`` - VERSIONED (``PV547761.1``). Metadata, read by
  GenBankFetcher to notice that a record has been revised and needs re-fetching
  on an update run.

Before this module the bare form was produced ad hoc by ``split('.')[0]`` in a
dozen scripts. That idiom is only equivalent to "drop the version" when the
input has exactly one dot; the tests below pin every shape where it is not.
"""

import pytest

import accession_utils
from accession_utils import (
    accession_from_filename,
    is_accession,
    normalise_accession,
    normalise_accession_series,
    split_accession_version,
    strip_known_suffixes,
)


# ---------------------------------------------------------------------------
# normalise_accession
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("PV547761.1", "PV547761"),
        ("PV547761", "PV547761"),
        ("NC_001542.1", "NC_001542"),
        ("NC_001542", "NC_001542"),
        ("KX148218.2", "KX148218"),
        ("AF123456.10", "AF123456"),
        ("  NC_001542.2  ", "NC_001542"),
    ],
)
def test_version_suffix_is_stripped_from_accessions(raw, expected):
    assert normalise_accession(raw) == expected


@pytest.mark.parametrize(
    "raw",
    [
        "EPI_ISL_402124",              # GISAID, not a GenBank accession
        "A/swine/Iowa/4.1/1976",       # an influenza strain name
        "cluster.1",                   # a label that happens to end in .1
        "NC_002016_4_plus",            # a segmented-run group name
        "rabv-7Jul26.db",
    ],
)
def test_labels_are_returned_verbatim(raw):
    """Deliberately conservative, in the same spirit as segment_utils: a token
    that is not accession-shaped is a label and must survive intact. This is the
    row where split('.')[0] silently truncated a strain name."""
    assert normalise_accession(raw) == raw


@pytest.mark.parametrize("blank", [None, "", "   ", "\t\n"])
def test_absent_is_none_never_empty_string(blank):
    assert normalise_accession(blank) is None


def test_normalisation_is_idempotent():
    once = normalise_accession("NC_001542.1")
    assert normalise_accession(once) == once


def test_case_is_preserved():
    """Accessions are compared verbatim downstream; normalising case here would
    break the join rather than help it."""
    assert normalise_accession("nc_001542.1") == "nc_001542"


def test_a_four_digit_version_is_not_a_version():
    """No GenBank record has a four-digit version, and without the bound a
    label like 'ABC12345.2021' would lose its year."""
    assert normalise_accession("ABC12345.2021") == "ABC12345.2021"


# ---------------------------------------------------------------------------
# split_accession_version
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("PV547761.1", ("PV547761", 1)),
        ("PV547761.12", ("PV547761", 12)),
        ("PV547761", ("PV547761", None)),
        ("cluster.1", ("cluster.1", None)),
        (None, (None, None)),
        ("", (None, None)),
    ],
)
def test_split_accession_version(raw, expected):
    assert split_accession_version(raw) == expected


def test_split_is_the_inverse_view_of_normalise():
    """The one caller that needs the version - revision detection on an update
    run - must get the same bare form everything else joins on."""
    for raw in ("PV547761.1", "PV547761", "NC_001542.3"):
        assert split_accession_version(raw)[0] == normalise_accession(raw)


def test_a_revised_record_is_detectable():
    """The property GenBankFetcher relies on: same identity, higher version."""
    old_base, old_version = split_accession_version("PV547761.1")
    new_base, new_version = split_accession_version("PV547761.2")

    assert old_base == new_base
    assert new_version > old_version


# ---------------------------------------------------------------------------
# is_accession
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("NC_001542", True),
        ("NC_001542.1", True),
        ("KX148218", True),
        ("PV547761", True),
        ("AF123456", True),
        ("EPI_ISL_402124", False),
        ("A/swine/Iowa/4/1976", False),
        ("", False),
        (None, False),
        ("cluster.1", False),
    ],
)
def test_is_accession(raw, expected):
    assert is_accession(raw) is expected


# ---------------------------------------------------------------------------
# strip_known_suffixes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, expected",
    [
        ("NC_001542.fasta", "NC_001542"),
        ("NC_001542.fa", "NC_001542"),
        ("NC_001542.aligned.fasta", "NC_001542"),
        ("NC_001542.fasta.gz", "NC_001542"),
        ("NC_001542.1.fasta", "NC_001542.1"),   # .1 is NOT an extension
        ("NC_001542.errors.csv", "NC_001542.errors"),
        ("NC_001542", "NC_001542"),
        (".DS_Store", ".DS_Store"),
        ("archive.tar.bz2", "archive.tar"),
    ],
)
def test_strip_known_suffixes(name, expected):
    assert strip_known_suffixes(name) == expected


def test_an_unknown_extension_is_left_alone():
    """Stripping whatever follows the last dot is exactly the mistake this
    module exists to stop."""
    assert strip_known_suffixes("NC_001542.weirdext") == "NC_001542.weirdext"


# ---------------------------------------------------------------------------
# accession_from_filename
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "path, expected",
    [
        ("/tmp/Blast/master_seq/NC_001542.fasta", "NC_001542"),
        ("/tmp/Blast/master_seq/NC_001542.1.fasta", "NC_001542"),
        ("NC_001542.1.aligned.fasta", "NC_001542"),
        ("tmp/Blast/grouped_fasta/KX148218.fa", "KX148218"),
        ("KX148218", "KX148218"),
        ("PV547761.2.fa", "PV547761"),
    ],
)
def test_accession_from_filename(path, expected):
    assert accession_from_filename(path) == expected


def test_dotfiles_do_not_collapse_to_empty():
    """An empty accession became a directory name and the value of a
    command-line flag, so nextalign silently took the next flag as its value."""
    assert accession_from_filename(".DS_Store") == ".DS_Store"
    assert accession_from_filename(".snakemake_timestamp") == ".snakemake_timestamp"


def test_group_names_from_segmented_runs_survive():
    """Segmented runs name grouped_fasta entries '<ref>_<segment>_<strand>'."""
    assert accession_from_filename("NC_002016_4_plus.fa") == "NC_002016_4_plus"


def test_empty_path_is_empty_not_an_error():
    assert accession_from_filename("") == ""


@pytest.mark.parametrize("accession", ["NC_001542", "NC_001542.1", "KX148218", "PV547761.2"])
def test_filename_round_trips_to_the_canonical_accession(accession):
    """THE invariant the pipeline depends on: whatever spelling a reference list
    or an update DB uses, a file named after it must resolve to the same key
    the rest of the pipeline joins on."""
    for extension in (".fa", ".fasta", ".aligned.fasta"):
        assert accession_from_filename(accession + extension) == normalise_accession(accession)


# ---------------------------------------------------------------------------
# normalise_accession_series
# ---------------------------------------------------------------------------


def test_series_preserves_length_and_position():
    """Callers assign the result straight back into a dataframe column, so a
    dropped element would silently shift every row after it."""
    values = ["PV547761.1", None, "", "cluster.1", "NC_001542"]

    result = normalise_accession_series(values)

    assert len(result) == len(values)
    assert result == ["PV547761", "", "", "cluster.1", "NC_001542"]


def test_series_accepts_any_iterable():
    assert normalise_accession_series(iter(["NC_001542.1"])) == ["NC_001542"]


# ---------------------------------------------------------------------------
# the properties the pipeline actually relies on
# ---------------------------------------------------------------------------


class TestPipelineInvariants:
    """These are the properties other stages assume without asserting them."""

    ACCESSIONS = ["NC_001542", "KX148218", "PV547761", "AF123456", "PX132041"]

    @pytest.mark.parametrize("version", ["", ".1", ".2", ".10"])
    def test_every_version_of_an_accession_maps_to_one_key(self, version):
        """Revisions of a record must not fork into separate DB rows."""
        for accession in self.ACCESSIONS:
            assert normalise_accession(accession + version) == accession

    def test_bare_and_versioned_forms_join(self):
        """The update-run question: an incoming 'NC_001542.1' and a stored
        'NC_001542' have to resolve to the same identity, or the merge either
        duplicates the row or silently misses it."""
        incoming = normalise_accession("NC_001542.1")
        stored = normalise_accession("NC_001542")

        assert incoming == stored

    def test_the_module_exposes_a_stable_surface(self):
        """These names are imported by call sites across scripts/."""
        for name in (
            "normalise_accession",
            "split_accession_version",
            "is_accession",
            "strip_known_suffixes",
            "accession_from_filename",
            "normalise_accession_series",
        ):
            assert callable(getattr(accession_utils, name))
