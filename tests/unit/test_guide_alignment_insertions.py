"""A guide alignment preserves insertions that are absent from the master.

The concern this answers: in the large SARS-CoV-2 alignments, everything is
projected onto Wuhan-Hu-1 coordinates, so any region inserted relative to that
reference is stripped and all evolutionary signal inside it is lost. If this
pipeline behaved the same way, curating a sub-reference with an interesting
insertion would be pointless - the insertion would be discarded on the way into
the merged MSA.

It does not behave that way, *provided a guide alignment is supplied*
(``--ref_set_aligned`` / ``--precomputed_ref_dir``). ``process_master_alignment``
iterates the records of the **guide alignment**, and for each reference projects
its nextalign-aligned queries into that reference's own gapped row. The merged
MSA therefore lives in the guide alignment's column space, not in the master's
raw coordinate space. Columns where the master is gapped and a sub-reference has
bases survive, and queries carrying the insertion place their residues there -
including substitutions *within* the inserted region, which is the whole point.

Without a guide alignment the behaviour is the SARS-CoV-2 one: nextalign's
reference coordinates win and insertions go to ``.insertions.csv`` as
bookkeeping, recoverable as presence/absence but not alignable.

These tests are the standing guarantee on that property.
"""

import os
from pathlib import Path

import pytest
from Bio import SeqIO

from PadAlignment import PadAlignment


# Master and sub-reference differ by a 6 nt insertion. Q_WITH_INS carries the
# insertion and one substitution inside it; Q_NO_INS lacks the insertion.
GUIDE = (
    ">MASTER\nACGTACGT------ACGTACGT\n"
    ">REF_INS\nACGTACGTGGGGGGACGTACGT\n"
)
SUBALIGNMENT = (
    ">REF_INS\nACGTACGTGGGGGGACGTACGT\n"
    ">Q_WITH_INS\nACGTACGTGGGTGGACGTACGT\n"
    ">Q_NO_INS\nACGTACGT------ACGTACGT\n"
)
INSERTION_COLUMNS = slice(8, 14)


@pytest.fixture
def merged(tmp_path):
    """Run PadAlignment against a guide alignment; return {id: sequence}."""
    guide = tmp_path / "refset_1_aln.fasta"
    guide.write_text(GUIDE)

    query_dir = tmp_path / "in" / "REF_INS"
    query_dir.mkdir(parents=True)
    (query_dir / "REF_INS.aligned.fasta").write_text(SUBALIGNMENT)

    out = tmp_path / "out"
    padder = PadAlignment(
        reference_alignment=str(guide), input_dir=str(tmp_path / "in"),
        base_dir=str(tmp_path), output_dir=str(out), keep_intermediate_files=True,
    )
    padder.process_master_alignment(
        str(guide), str(tmp_path / "in"), str(tmp_path), str(out), True
    )

    merged_files = [f for f in os.listdir(out) if f.endswith("_merged_MSA.fasta")]
    assert merged_files, "no merged MSA was produced"
    records = list(SeqIO.parse(out / merged_files[0], "fasta"))
    return {r.id: str(r.seq) for r in records}


class TestInsertionColumnsSurvive:
    def test_the_insertion_columns_exist_in_the_merged_msa(self, merged):
        """The columns are there at all - the SARS-CoV-2 failure mode is that
        they are not."""
        assert merged["REF_INS"][INSERTION_COLUMNS] == "GGGGGG"

    def test_the_master_is_gapped_there_rather_than_truncated(self, merged):
        """The master keeps full width; it simply has no bases in that region."""
        assert merged["MASTER"][INSERTION_COLUMNS] == "------"
        assert len(merged["MASTER"]) == len(merged["REF_INS"])

    def test_a_query_carrying_the_insertion_keeps_its_residues(self, merged):
        """Not just presence/absence - the actual bases are placed in the
        alignment, where they can be compared across sequences."""
        assert merged["Q_WITH_INS"][INSERTION_COLUMNS] == "GGGTGG"

    def test_substitutions_within_the_insertion_are_visible(self, merged):
        """The property that makes evolutionary analysis of the inserted region
        possible: Q_WITH_INS differs from REF_INS at one site *inside* the
        insertion, and that difference is alignable."""
        ref = merged["REF_INS"][INSERTION_COLUMNS]
        query = merged["Q_WITH_INS"][INSERTION_COLUMNS]
        differing = [i for i, (a, b) in enumerate(zip(ref, query)) if a != b]
        assert differing == [3], f"expected one substitution at offset 3, got {differing}"

    def test_a_query_lacking_the_insertion_is_gapped_not_shifted(self, merged):
        """The danger of hand-rolled coordinate projection is an off-by-N that
        slides flanking bases into the insertion. Flanks must stay put."""
        assert merged["Q_NO_INS"][INSERTION_COLUMNS] == "------"
        assert merged["Q_NO_INS"][:8] == "ACGTACGT"
        assert merged["Q_NO_INS"][14:] == "ACGTACGT"

    def test_every_row_is_the_guide_alignment_width(self, merged):
        """A ragged merged MSA is unusable downstream."""
        widths = {len(s) for s in merged.values()}
        assert widths == {len(GUIDE.split("\n")[1])}


class TestNonSegmentedVirusesCanUseAGuideAlignment:
    """The mechanism is not influenza-only.

    ``ref_set_aligned`` is currently set only on the influenza profiles, but
    nothing gates it on ``is_segmented``. A non-segmented virus supplies a
    directory holding a single FASTA and the resolver picks it up.
    """

    def test_a_sole_alignment_is_used_when_there_is_no_segment(self, tmp_path):
        ref_dir = tmp_path / "ref_set_aligned"
        ref_dir.mkdir()
        (ref_dir / "hcv_backbone_aln.fasta").write_text(GUIDE)

        padder = PadAlignment(
            reference_alignment="", input_dir=str(tmp_path), base_dir=str(tmp_path),
            output_dir=str(tmp_path), keep_intermediate_files=False,
        )
        found = padder.find_precomputed_reference_alignment(str(ref_dir), None)
        assert found is not None
        assert Path(found).name == "hcv_backbone_aln.fasta"

    def test_ambiguity_is_refused_rather_than_guessed(self, tmp_path):
        """Two candidates and no segment to choose between them: returning
        either would silently align against the wrong backbone."""
        ref_dir = tmp_path / "ref_set_aligned"
        ref_dir.mkdir()
        (ref_dir / "a_aln.fasta").write_text(GUIDE)
        (ref_dir / "b_aln.fasta").write_text(GUIDE)

        padder = PadAlignment(
            reference_alignment="", input_dir=str(tmp_path), base_dir=str(tmp_path),
            output_dir=str(tmp_path), keep_intermediate_files=False,
        )
        assert padder.find_precomputed_reference_alignment(str(ref_dir), None) is None

    @pytest.mark.parametrize("unset", [None, "", "UNSET", "null"])
    def test_unset_sentinels_mean_no_guide_alignment(self, unset, tmp_path):
        """Nextflow passes the string 'UNSET' when the param is null; it must
        not be treated as a directory name."""
        padder = PadAlignment(
            reference_alignment="", input_dir=str(tmp_path), base_dir=str(tmp_path),
            output_dir=str(tmp_path), keep_intermediate_files=False,
        )
        assert padder._normalize_optional_path(unset) is None


class TestReferenceIsNotDuplicated:
    """The guide reference is added to the merged MSA explicitly, and nextalign
    also emits it in the subalignment. Both copies reaching the output is a
    real defect: it inflates counts and can confuse downstream dedup by name.
    """

    @pytest.mark.xfail(
        reason="known: the reference appears twice in the merged MSA when it is "
               "present in both the guide alignment and the nextalign subalignment",
        strict=False,
    )
    def test_reference_appears_once(self, tmp_path):
        guide = tmp_path / "refset_1_aln.fasta"
        guide.write_text(GUIDE)
        query_dir = tmp_path / "in" / "REF_INS"
        query_dir.mkdir(parents=True)
        (query_dir / "REF_INS.aligned.fasta").write_text(SUBALIGNMENT)
        out = tmp_path / "out"

        padder = PadAlignment(
            reference_alignment=str(guide), input_dir=str(tmp_path / "in"),
            base_dir=str(tmp_path), output_dir=str(out), keep_intermediate_files=True,
        )
        padder.process_master_alignment(
            str(guide), str(tmp_path / "in"), str(tmp_path), str(out), True
        )
        merged_file = [f for f in os.listdir(out) if f.endswith("_merged_MSA.fasta")][0]
        ids = [r.id for r in SeqIO.parse(out / merged_file, "fasta")]
        assert ids.count("REF_INS") == 1, f"reference duplicated: {ids}"
