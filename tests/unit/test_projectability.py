"""Tests for the shared projectability predicate.

The defect these guard against: PadAlignment and CollectFilteredSequences each
carried their own answer to "can this reference be projected?". PadAlignment
answered it from the per-segment backbones (refset_<N>_aln.fasta, all references
present); CollectFilteredSequences answered it from Nextalign/reference_aln/
(only references nextalign could align to a *master*). On influenza those
disagree by 87 of 117 HA references, and the stricter, wrong answer became
--skip_ids and deleted 633,988 sequences. See MISSING_H3_report.md.
"""

from pathlib import Path

import pytest

import projectability


def write_fasta(path: Path, records, crlf: bool = False):
    """Write a FASTA. ``crlf=True`` reproduces the real refset files, which are
    CRLF-terminated - a naive header parse yields ids with a trailing \\r."""
    path.parent.mkdir(parents=True, exist_ok=True)
    eol = "\r\n" if crlf else "\n"
    body = eol.join(f">{name}{eol}{seq}" for name, seq in records) + eol
    path.write_bytes(body.encode("utf-8"))


def make_ref_list(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join("\t".join(r) + "\n" for r in rows), encoding="utf-8")


# --------------------------------------------------------------------------
# read_fasta_ids
# --------------------------------------------------------------------------

def test_read_fasta_ids_strips_version_and_description(tmp_path):
    fasta = tmp_path / "refs.fasta"
    write_fasta(fasta, [("AB573800.1 Influenza A virus segment 4", "ACGT"), ("MZ198323", "ACGT")])
    assert projectability.read_fasta_ids(str(fasta)) == {"AB573800", "MZ198323"}


def test_read_fasta_ids_is_crlf_safe(tmp_path):
    """generic/influenza/ref_set_aligned/*.fasta are CRLF; ids must not carry \\r."""
    fasta = tmp_path / "crlf.fasta"
    write_fasta(fasta, [("MZ198323", "ACGT"), ("DQ864721", "ACGT")], crlf=True)
    ids = projectability.read_fasta_ids(str(fasta))
    assert ids == {"MZ198323", "DQ864721"}
    assert not any(i.endswith("\r") for i in ids)


def test_read_fasta_ids_missing_file_is_empty_not_error(tmp_path):
    assert projectability.read_fasta_ids(str(tmp_path / "nope.fasta")) == set()


# --------------------------------------------------------------------------
# normalise_optional_path
# --------------------------------------------------------------------------

@pytest.mark.parametrize("value", [None, "", "null", "NULL", "none", "UNSET", "  "])
def test_unset_sentinels_normalise_to_none(value):
    """Nextflow interpolates absent params as the literal 'null'/'UNSET'."""
    assert projectability.normalise_optional_path(value) is None


def test_real_path_survives_normalisation():
    assert projectability.normalise_optional_path(" /tmp/x ") == "/tmp/x"


# --------------------------------------------------------------------------
# find_precomputed_reference_alignment
# --------------------------------------------------------------------------

def test_canonical_refset_name_is_preferred(tmp_path):
    write_fasta(tmp_path / "refset_4_aln.fasta", [("A", "ACGT")])
    write_fasta(tmp_path / "refset_6_aln.fasta", [("B", "ACGT")])
    resolved = projectability.find_precomputed_reference_alignment(str(tmp_path), "4")
    assert Path(resolved).name == "refset_4_aln.fasta"


def test_segment_normalisation_applies(tmp_path):
    """'04' and '4.0' and 'segment 4' are all segment 4 (segment_utils contract)."""
    write_fasta(tmp_path / "refset_4_aln.fasta", [("A", "ACGT")])
    for spelling in ("4", "04", "4.0", "segment 4"):
        resolved = projectability.find_precomputed_reference_alignment(str(tmp_path), spelling)
        assert resolved is not None, spelling
        assert Path(resolved).name == "refset_4_aln.fasta"


def test_sole_fasta_used_when_segment_unknown(tmp_path):
    write_fasta(tmp_path / "only_one.fasta", [("A", "ACGT")])
    assert projectability.find_precomputed_reference_alignment(str(tmp_path), None) is not None


def test_ambiguous_dir_with_unknown_segment_resolves_nothing(tmp_path):
    write_fasta(tmp_path / "a.fasta", [("A", "ACGT")])
    write_fasta(tmp_path / "b.fasta", [("B", "ACGT")])
    assert projectability.find_precomputed_reference_alignment(str(tmp_path), None) is None


def test_missing_dir_resolves_nothing(tmp_path):
    assert projectability.find_precomputed_reference_alignment(str(tmp_path / "absent"), "4") is None


# --------------------------------------------------------------------------
# projectable_reference_ids - the predicate itself
# --------------------------------------------------------------------------

def test_precomputed_backbones_beat_reference_aln(tmp_path):
    """THE regression. A reference nextalign could not align to the master is
    still projectable when the per-segment backbone contains it.

    This is influenza segment 4 in miniature: master AB573800 is H1, the H3 and
    H7 references are too divergent for nextalign to align against it, but
    refset_4_aln.fasta holds all three and PadAlignment projects all three.
    """
    ref_dir = tmp_path / "ref_set_aligned"
    write_fasta(ref_dir / "refset_4_aln.fasta",
                [("MASTER_H1", "ACGT"), ("REF_H3", "ACGT"), ("REF_H7", "ACGT")])

    nextalign = tmp_path / "Nextalign"
    # nextalign only managed the master itself against the master
    write_fasta(nextalign / "reference_aln" / "MASTER_H1" / "MASTER_H1.aligned.fasta",
                [("MASTER_H1", "ACGT")])

    ids, source = projectability.projectable_reference_ids(
        nextalign_dir=str(nextalign),
        precomputed_ref_dir=str(ref_dir),
        master_segment_map={"MASTER_H1": "4"},
    )
    assert ids == {"MASTER_H1", "REF_H3", "REF_H7"}
    assert source == "precomputed"


def test_reference_aln_used_when_no_precomputed_dir(tmp_path):
    """Non-segmented RABV/HCV behaviour is unchanged: reference_aln is the answer."""
    nextalign = tmp_path / "Nextalign"
    write_fasta(nextalign / "reference_aln" / "MASTER" / "MASTER.aligned.fasta",
                [("MASTER", "ACGT"), ("REF_OK", "ACGT")])

    ids, source = projectability.projectable_reference_ids(
        nextalign_dir=str(nextalign), precomputed_ref_dir=None, master_segment_map={})
    assert ids == {"MASTER", "REF_OK"}
    assert source == "reference_aln"


def test_falls_back_per_master_when_a_segment_backbone_is_missing(tmp_path):
    """Segment 4 has a backbone, segment 6 does not: use each master's best source."""
    ref_dir = tmp_path / "ref_set_aligned"
    write_fasta(ref_dir / "refset_4_aln.fasta", [("M4", "ACGT"), ("REF_H3", "ACGT")])

    nextalign = tmp_path / "Nextalign"
    write_fasta(nextalign / "reference_aln" / "M6" / "M6.aligned.fasta",
                [("M6", "ACGT"), ("REF_N2", "ACGT")])

    ids, source = projectability.projectable_reference_ids(
        nextalign_dir=str(nextalign),
        precomputed_ref_dir=str(ref_dir),
        master_segment_map={"M4": "4", "M6": "6"},
    )
    assert ids == {"M4", "REF_H3", "M6", "REF_N2"}
    assert source == "mixed"


def test_union_spans_all_segments(tmp_path):
    """A reference is projectable if ANY master's backbone holds it - PadAlignment
    iterates masters, so the union is the correct set, not a per-segment slice."""
    ref_dir = tmp_path / "ref_set_aligned"
    write_fasta(ref_dir / "refset_4_aln.fasta", [("M4", "ACGT"), ("HA_REF", "ACGT")])
    write_fasta(ref_dir / "refset_6_aln.fasta", [("M6", "ACGT"), ("NA_REF", "ACGT")])

    ids, _ = projectability.projectable_reference_ids(
        precomputed_ref_dir=str(ref_dir), master_segment_map={"M4": "4", "M6": "6"})
    assert ids == {"M4", "HA_REF", "M6", "NA_REF"}


# --------------------------------------------------------------------------
# load_master_segment_map
# --------------------------------------------------------------------------

def test_master_segment_map_reads_only_master_rows(tmp_path):
    ref_list = tmp_path / "ref_list_refmast.txt"
    make_ref_list(ref_list, [
        ("AB573800", "master", "4"),
        ("MZ198323", "reference", "4"),
        ("AB472016", "master", "6"),
        ("IBV_X", "exclusion_list", "B"),
    ])
    assert projectability.load_master_segment_map(str(ref_list)) == {
        "AB573800": "4", "AB472016": "6"}


def test_master_segment_map_empty_for_unsegmented_ref_list(tmp_path):
    """generic/rabv and generic/hcv ref lists are 2-column with no segment; an
    empty map must read as 'not segmented', which routes to reference_aln."""
    ref_list = tmp_path / "ref_list.txt"
    make_ref_list(ref_list, [("NC_001542", "master"), ("KJ1", "reference")])
    assert projectability.load_master_segment_map(str(ref_list)) == {}


def test_master_segment_map_missing_file_is_empty(tmp_path):
    assert projectability.load_master_segment_map(str(tmp_path / "nope.txt")) == {}
