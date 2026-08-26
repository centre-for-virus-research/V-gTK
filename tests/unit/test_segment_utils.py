"""Regression tests for segment identity.

Every case here is drawn from a real failure in a GenBank+GISAID influenza
build, not invented. Before ``segment_utils`` existed, six scripts each carried
their own ``_normalize_segment`` and three of the disagreements were silent data
corruption rather than mere inconsistency.
"""

import pytest

import segment_utils
from BlastAlignment import BlastAlignment
from CalcAlignmentCord import CalculateAlignmentCoordinates
from GenBankParser import GenBankParser
from NextalignAlignment import NextalignAlignment
from PadAlignment import PadAlignment
from CreateSqliteDB import CreateSqliteDB
from ExportRefListFromUpdateDb import _segment_norm as reflist_segment_norm
from ExportUpdateAssets import _segment_norm as update_assets_segment_norm
from SegmentPivotTable import normalise_segment_label
from UsherPlacement import UsherPlacement

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

FLU_SEGMENT_NAMES = REPO_ROOT / "generic" / "influenza" / "segment_names_iav.tsv"

#: Influenza A segment number for each gene/segment name. The polymerase pair is
#: the trap: PB2 is segment 1 and PB1 is segment 2, so any rule that reads the
#: digit out of the name gets both of them exactly backwards.
TRUE_SEGMENT_FOR_NAME = {
    "PB2": "1", "PB1": "2", "PA": "3", "HA": "4",
    "NP": "5", "NA": "6", "MP": "7", "NS": "8",
    "M1": "7", "M2": "7", "NS1": "8", "NEP": "8",
}


class TestNormaliseSegment:
    @pytest.mark.parametrize("raw,expected", [
        ("4", "4"),
        ("04", "4"),        # zero-padded
        ("4.0", "4"),       # float-ish; the old digit-scrape made this '40'
        (" 4 ", "4"),
        ("Segment 4", "4"),  # GenBank /segment= freehand
        ("segment 7", "7"),
        ("RNA 4", "4"),
        ("seg-6", "6"),
        ("sgt_3", "3"),
    ])
    def test_numeric_forms_collapse_to_the_integer(self, raw, expected):
        assert segment_utils.normalise_segment(raw) == expected

    @pytest.mark.parametrize("raw", ["L", "M", "S", "HA", "NA", "PB1", "PB2", "NS1"])
    def test_labels_survive_verbatim(self, raw):
        assert segment_utils.normalise_segment(raw) == raw

    @pytest.mark.parametrize("raw", [None, "", "   "])
    def test_absent_is_none_never_zero_or_empty_string(self, raw):
        # BlastAlignment used to return "0" for a missing segment, which reads
        # downstream as a real segment rather than as an absent one.
        assert segment_utils.normalise_segment(raw) is None

    def test_a_gene_name_is_never_mistaken_for_its_embedded_digit(self):
        """The PB1/PB2 inversion, guarded directly."""
        assert segment_utils.normalise_segment("PB1") == "PB1"
        assert segment_utils.normalise_segment("PB2") == "PB2"
        assert segment_utils.normalise_segment("PB1") != "1"
        assert segment_utils.normalise_segment("PB2") != "2"

    def test_a_non_integral_decimal_is_not_silently_truncated(self):
        # '2.5' is not segment 2; it is malformed and must stay visible as such.
        assert segment_utils.normalise_segment("2.5") == "2.5"


class TestSegmentNameAsset:
    @pytest.fixture(scope="class")
    def flu(self):
        return segment_utils.load_segment_names(str(FLU_SEGMENT_NAMES))

    def test_asset_exists(self):
        assert FLU_SEGMENT_NAMES.is_file(), f"missing asset: {FLU_SEGMENT_NAMES}"

    @pytest.mark.parametrize("name,segment", sorted(TRUE_SEGMENT_FOR_NAME.items()))
    def test_every_gene_name_maps_to_its_true_segment(self, flu, name, segment):
        _, by_name = flu
        assert segment_utils.name_to_segment(name, by_name) == segment

    def test_name_lookup_is_case_insensitive(self, flu):
        _, by_name = flu
        assert segment_utils.name_to_segment("ha", by_name) == "4"
        assert segment_utils.name_to_segment("Ha", by_name) == "4"

    def test_neuraminidase_survives_as_a_real_name(self, flu):
        """`NA` is a segment name, not a null. It must round-trip."""
        by_segment, by_name = flu
        assert segment_utils.name_to_segment("NA", by_name) == "6"
        assert segment_utils.segment_to_name("6", by_segment) == "NA"

    @pytest.mark.parametrize("segment,name", [
        ("1", "PB2"), ("2", "PB1"), ("3", "PA"), ("4", "HA"),
        ("5", "NP"), ("6", "NA"), ("7", "MP"), ("8", "NS"),
    ])
    def test_segment_to_name_round_trips(self, flu, segment, name):
        by_segment, by_name = flu
        assert segment_utils.segment_to_name(segment, by_segment) == name
        assert segment_utils.name_to_segment(name, by_name) == segment

    def test_already_named_segments_pass_through(self, flu):
        """A virus keyed on L/S already has names; do not blank them."""
        by_segment, _ = flu
        assert segment_utils.segment_to_name("L", by_segment) == "L"
        assert segment_utils.segment_to_name("S", by_segment) == "S"

    def test_unknown_numeric_segment_has_no_name(self, flu):
        by_segment, _ = flu
        assert segment_utils.segment_to_name("99", by_segment) == ""

    def test_missing_asset_degrades_quietly(self):
        by_segment, by_name = segment_utils.load_segment_names("/nonexistent/path.tsv")
        assert by_segment == {} and by_name == {}


class TestAvailabilityIsContextual:
    """`NA` means neuraminidase to influenza and 'missing' to everyone else.

    Which one it is cannot be decided from the string, so it is decided from the
    reference list's declared label set.
    """

    NUMERIC = {"1", "2", "3", "4", "5", "6", "7", "8"}
    NAMED = {"pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"}

    def test_na_is_missing_when_the_virus_keys_on_numbers(self):
        assert segment_utils.is_unavailable("NA", self.NUMERIC) is True

    def test_na_is_a_real_segment_when_the_virus_keys_on_names(self):
        assert segment_utils.is_unavailable("NA", self.NAMED) is False

    def test_real_numeric_segment_is_available(self):
        assert segment_utils.is_unavailable("6", self.NUMERIC) is False

    def test_without_a_label_set_the_conservative_default_applies(self):
        assert segment_utils.is_unavailable("NA") is True
        assert segment_utils.is_unavailable("4") is False


class TestAllNormalisersNowAgree:
    """The ten former copies must not drift apart again.

    Each call site keeps its own *missing* sentinel because its callers depend on
    it, so missing values are excluded here; what must agree is the
    normalisation of a present value.
    """

    PRESENT_CASES = ["4", "04", "4.0", " 4 ", "segment 7", "RNA 4",
                     "PB1", "PB2", "M1", "NS1", "HA", "L", "S", "8"]

    @staticmethod
    def _genbank(value):
        parser = GenBankParser.__new__(GenBankParser)
        parser.is_segmented_virus = "Y"
        return parser._normalize_segment_value(value)

    #: Every former copy of the normaliser, with the missing-value sentinel its
    #: own callers depend on. Ten call sites, one rule.
    ALL_NORMALISERS = None  # populated below the class body

    @pytest.mark.parametrize("raw", PRESENT_CASES)
    def test_every_normaliser_produces_the_same_canonical_value(self, raw):
        expected = segment_utils.normalise_segment(raw)
        for label, fn, _sentinel in self.ALL_NORMALISERS:
            assert fn(raw) == expected, f"{label} disagreed on {raw!r}"

    @pytest.mark.parametrize("gene,segment", sorted(TRUE_SEGMENT_FOR_NAME.items()))
    def test_no_normaliser_turns_a_gene_name_into_a_wrong_segment_number(self, gene, segment):
        """The regression that mattered most: PB2->2 and PB1->1, inverted.

        A normaliser may keep the name, but it must never emit a *different*
        segment number than the one that gene actually belongs to.
        """
        for label, fn, _sentinel in self.ALL_NORMALISERS:
            produced = fn(gene)
            if produced is not None and str(produced).isdigit():
                assert str(produced) == segment, (
                    f"{label}: {gene} belongs to segment {segment} but it produced {produced}"
                )

    @pytest.mark.parametrize("missing", [None, "", "nan"])
    def test_each_call_site_keeps_its_own_missing_sentinel(self, missing):
        for label, fn, sentinel in self.ALL_NORMALISERS:
            assert fn(missing) == sentinel, f"{label} changed its missing sentinel"


def _genbank_normaliser(value):
    parser = GenBankParser.__new__(GenBankParser)
    parser.is_segmented_virus = "Y"
    return parser._normalize_segment_value(value)


TestAllNormalisersNowAgree.ALL_NORMALISERS = (
    ("GenBankParser", _genbank_normaliser, ""),
    ("BlastAlignment", BlastAlignment._normalize_segment, "0"),
    ("PadAlignment", PadAlignment._normalize_segment, None),
    ("NextalignAlignment", NextalignAlignment._normalize_segment, ""),
    ("CalcAlignmentCord", CalculateAlignmentCoordinates._normalize_segment_value, ""),
    ("SegmentPivotTable", normalise_segment_label, None),
    ("CreateSqliteDB", CreateSqliteDB._normalize_segment_value, ""),
    ("ExportUpdateAssets", update_assets_segment_norm, "0"),
    ("ExportRefListFromUpdateDb", reflist_segment_norm, "0"),
    ("UsherPlacement", UsherPlacement._normalize_segment, "0"),
)
