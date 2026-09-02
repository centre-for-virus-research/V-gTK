"""Guards on the accession CALL SITES that still used ``split('.')[0]``.

``scripts/accession_utils.py`` is the single authority on accession identity,
and ``tests/unit/test_accession_utils.py`` pins its rules. This file pins the
four *call sites* that were still deciding identity for themselves, because a
helper with perfect coverage is worth nothing while the production path does not
call it - and each of these sites was wrong in its own way:

``GenBankSequenceSubmitter``
    ``path_to_basename`` cut a filename at its first dot, and one line cut a
    whole *path* at its first dot, which loses the directory.
``CollectFilteredSequences``
    wrote ``filtered_sequences_ids.txt`` in a mixture of spellings. That file is
    read verbatim by ``PadAlignment --skip_ids``, which compares it against
    version-stripped record ids, so anything non-bare in it skipped nothing.
``CheckSegmentRetention``
    version-stripped three of the four id streams it subtracts from each other.
    The fourth was the reference list, so on versioned input the guard's
    subtractions removed nothing and the guard passed segments it exists to fail.
``projectability``
    truncated any backbone header id that is not a plain accession.

Every accession this repo ships is bare, so all of these are no-ops on current
data - which is exactly why they survived. The cases below are the shapes where
the old idiom and the canonical rule disagree.
"""

import subprocess
from pathlib import Path

import pytest

import accession_utils
import projectability
from CheckSegmentRetention import check, load_excused
from CollectFilteredSequences import (
    _normalize_accession,
    collect_filtered_sequences,
    write_filtered_ids_only,
)
from GenBankSequenceSubmitter import GenBankSequenceSubmitter


# ---------------------------------------------------------------------------
# GenBankSequenceSubmitter - both sites parse a FILENAME, so both delegate to
# accession_from_filename rather than to normalise_accession.
# ---------------------------------------------------------------------------


class TestGenBankSequenceSubmitter:

    @pytest.mark.parametrize(
        "path",
        [
            "/tmp/x/NC_001542.fasta",
            "/tmp/x/NC_001542.1.fasta",
            "NC_001542.1.aligned.fasta",
            "NC_001542_aln.fasta",
            "NC_002016_4_plus.fa",
            "KX148218",
        ],
    )
    def test_path_to_basename_is_the_shared_authority(self, path):
        """One rule for "what accession is this file named after", not two.

        NextalignAlignment.path_to_basename delegates to the same function; the
        two scripts name files for each other's stages, so a private rule here
        is a join waiting to break.
        """
        assert GenBankSequenceSubmitter.path_to_basename(path) == \
            accession_utils.accession_from_filename(path)

    def test_a_dotfile_no_longer_collapses_to_an_empty_accession(self):
        """``.split('.')[0]`` returned ``''`` for any dot-leading name.

        The call sites append ``"_aln.fasta"`` to this value, so every stray
        dotfile in a scanned directory produced the same output name,
        ``_aln.fasta``, silently overwriting the previous one.
        """
        assert GenBankSequenceSubmitter.path_to_basename(".DS_Store") == ".DS_Store"
        assert GenBankSequenceSubmitter.path_to_basename(".snakemake_timestamp") != ""

    def test_mafft_query_sequences_keeps_the_directory_of_a_dotted_path(self, tmp_path, monkeypatch):
        """The worse of the two sites: it split a PATH, not a filename.

        ``process()`` passes ``<reference_alignments>/<accession>`` with no
        extension and this line rebuilds the ``_aln.fasta`` sibling that
        ``extract_matching_sequences()`` wrote. Cutting the path at its first dot
        threw the directory away as soon as any parent component contained one -
        a ``--tmp_dir`` of ``./tmp`` or ``../run.v2`` is enough - and mafft was
        handed a path that does not exist.
        """
        run_dir = tmp_path / "run.v2" / "reference_alignments"
        run_dir.mkdir(parents=True)
        out_dir = tmp_path / "query_ref_alignment"
        out_dir.mkdir()
        captured = {}

        def fake_run(command, **kwargs):
            captured["command"] = command
            return subprocess.CompletedProcess(command, 0)

        monkeypatch.setattr(subprocess, "run", fake_run)

        submitter = GenBankSequenceSubmitter(
            "seqs", str(tmp_path), "Table2asn", "meta.tsv", "template.sbt",
            "ref.gff3", "vgtk.db", 30,
        )
        submitter.mafft_query_sequences(
            str(tmp_path / "grouped" / "NC_001542.fasta"),
            str(run_dir / "NC_001542"),
            str(out_dir),
        )

        assert captured["command"][-1] == str(run_dir / "NC_001542_aln.fasta")
        assert (out_dir / "NC_001542_aln.fasta").exists()

    def test_a_dot_free_path_is_rebuilt_exactly_as_before(self, tmp_path, monkeypatch):
        """The no-op half of the same site: on the paths this repo actually
        produces, the new construction is byte-identical to the old idiom."""
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        legacy_input = "/data/tmp/analysis/abc-123/reference_alignments/NC_001542"
        captured = {}

        def fake_run(command, **kwargs):
            captured["command"] = command
            return subprocess.CompletedProcess(command, 0)

        monkeypatch.setattr(subprocess, "run", fake_run)

        submitter = GenBankSequenceSubmitter(
            "seqs", str(tmp_path), "Table2asn", "meta.tsv", "template.sbt",
            "ref.gff3", "vgtk.db", 30,
        )
        submitter.mafft_query_sequences("query.fa", legacy_input, str(out_dir))

        assert captured["command"][-1] == legacy_input.split(".")[0] + "_aln.fasta"


# ---------------------------------------------------------------------------
# CollectFilteredSequences - a FASTA header / seqName cell is an ACCESSION with
# a description after it, not a filename: first token, then normalise_accession.
# ---------------------------------------------------------------------------


class TestCollectFilteredSequences:

    @pytest.mark.parametrize(
        "raw, expected",
        [
            ("NC_001542.1 Rabies lyssavirus, complete genome", "NC_001542"),
            ("NC_001542 Rabies lyssavirus", "NC_001542"),
            ("  KX148218.2  ", "KX148218"),
            ("KX148218", "KX148218"),
        ],
    )
    def test_first_token_still_wins_and_the_version_still_goes(self, raw, expected):
        """The ``.split()[0]`` half of the old chain was doing real work - these
        inputs are FASTA header lines - and must survive the fix."""
        assert _normalize_accession(raw) == expected

    def test_a_strain_name_is_no_longer_truncated_mid_name(self):
        """Influenza strain names routinely carry a dot. ``split('.')[0]`` turned
        this one into ``A/swine/Iowa/4``, a different sequence's id."""
        assert _normalize_accession("A/swine/Iowa/4.1/1976 H1N1") == "A/swine/Iowa/4.1/1976"

    @pytest.mark.parametrize("blank", ["", "   ", "\t\n"])
    def test_an_empty_header_is_empty_not_an_exception(self, blank):
        """``"".split()[0]`` raised IndexError inside the per-file try/except,
        which abandoned the parse of every remaining record in that file."""
        assert _normalize_accession(blank) == ""

    def _errors_csv(self, nextalign_dir, ref_dir_name, rows):
        path = nextalign_dir / "query_aln" / ref_dir_name / "chunk.errors.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        body = "seqName,errors,warnings\n" + "".join(f"{n},{e},\n" for n, e in rows)
        path.write_text(body, encoding="utf-8")

    def test_ids_file_carries_the_spelling_padalignment_compares_against(self, tmp_path):
        """The consequence, end to end.

        ``PadAlignment._load_skip_ids`` keeps each line's first token VERBATIM,
        then drops the version from each record id before the lookup::

            ids.add(line.strip().split()[0])
            ...
            [r for r in seqs if r.id.split('.')[0] not in self.skip_ids]

        So the writer has to emit the bare form. It used to emit the raw
        ``seqName`` cell - description, version and all - and every such id
        skipped nothing at all.
        """
        nextalign_dir = tmp_path / "Nextalign"
        self._errors_csv(nextalign_dir, "NC_001542", [("PV547761.1 some description", "frame shift")])

        filtered = collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "out.tsv"))
        ids_file = write_filtered_ids_only(filtered, str(tmp_path / "out.tsv"))
        lines = Path(ids_file).read_text(encoding="utf-8").split()

        assert lines == ["PV547761"]

        # PadAlignment's side of the join, reproduced exactly.
        skip_ids = {line.strip().split()[0] for line in Path(ids_file).read_text(
            encoding="utf-8").splitlines() if line.strip()}
        assert "PV547761.1".split(".")[0] in skip_ids

    def test_the_reference_column_is_canonical_too(self, tmp_path):
        """The other two collectors normalise the ``query_aln/<ref>`` directory
        name; this branch did not, so one run could report the same reference
        under two spellings."""
        nextalign_dir = tmp_path / "Nextalign"
        self._errors_csv(nextalign_dir, "NC_001542.1", [("KX148218", "too many Ns")])

        filtered = collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "out.tsv"))

        assert filtered["KX148218"]["reference"] == "NC_001542"

    def test_bare_input_is_untouched(self, tmp_path):
        """The no-op proof: on the spelling every dataset in this repo uses, the
        writer emits exactly what it emitted before."""
        nextalign_dir = tmp_path / "Nextalign"
        self._errors_csv(nextalign_dir, "REF_A", [("SEQ_1", "frame shift")])

        filtered = collect_filtered_sequences(str(nextalign_dir), str(tmp_path / "out.tsv"))

        assert filtered == {
            "SEQ_1": {"reference": "REF_A", "error": "frame shift", "warnings": ""}
        }


# ---------------------------------------------------------------------------
# CheckSegmentRetention - four id streams, one rule.
# ---------------------------------------------------------------------------


class TestCheckSegmentRetention:

    def _write(self, tmp_path, ref_seg_rows, merged_headers, tophit_rows, excused=()):
        ref_seg = tmp_path / "ref_seq_seg.tsv"
        ref_seg.write_text("".join(f"{a}\t{s}\n" for a, s in ref_seg_rows), encoding="utf-8")

        merged = tmp_path / "refset_4_aln_merged_MSA.fasta"
        merged.write_text("".join(f">{h}\nACGT\n" for h in merged_headers), encoding="utf-8")

        tophit = tmp_path / "query_uniq_tophit_annotated.tsv"
        tophit.write_text("".join("\t".join(r) + "\n" for r in tophit_rows), encoding="utf-8")

        filtered_ids = tmp_path / "filtered_sequences_ids.txt"
        filtered_ids.write_text("".join(f"{a}\n" for a in excused), encoding="utf-8")

        return [str(merged)], str(tophit), str(ref_seg), str(filtered_ids)

    def test_a_versioned_reference_list_no_longer_disarms_the_guard(self, tmp_path):
        """The reference row PadAlignment inserts into every block must be
        subtracted from ``retained``.

        ``read_fasta_ids`` gave ``NC_001542`` while ``load_reference_segments``
        gave ``NC_001542.1``, so the subtraction removed nothing and the backbone
        row was counted as a retained query. Here that turns a real 3-of-4 loss
        into a reported 4 of 4 and the guard returns no failure at all.
        """
        queries = ["KX148218", "KX148219", "KX148220", "KX148221"]
        merged, tophit, ref_seg, filtered = self._write(
            tmp_path,
            ref_seg_rows=[("NC_001542.1", "4")],
            merged_headers=["NC_001542"] + queries[:3],
            tophit_rows=[(q, "NC_001542", "99.0", "plus", "4") for q in queries],
        )

        failures = check(merged, tophit, ref_seg, filtered, min_retention=0.98)

        assert len(failures) == 1
        segment, retained, expected, ratio, _examples = failures[0]
        assert (segment, retained, expected) == ("4", 3, 4)
        assert ratio == pytest.approx(0.75)

    def test_a_versioned_reference_column_still_maps_to_its_segment(self, tmp_path):
        """The fallback lookup, for a tophit file written before column 5 existed.

        ``ref_segments.get(reference)`` was a raw-vs-raw lookup. A versioned
        reference column missed the map, every row was dropped as segment-less,
        and ``expected`` came out 0 - which this guard scores as ratio 1.0. A
        segment that lost every single query reported perfect health.
        """
        queries = ["KX148218", "KX148219", "KX148220", "KX148221"]
        merged, tophit, ref_seg, filtered = self._write(
            tmp_path,
            ref_seg_rows=[("NC_001542", "4")],
            merged_headers=["NC_001542"],
            tophit_rows=[(q, "NC_001542.1", "99.0", "plus") for q in queries],
        )

        failures = check(merged, tophit, ref_seg, filtered, min_retention=0.98)

        assert len(failures) == 1
        segment, retained, expected, ratio, _examples = failures[0]
        assert (segment, retained, expected, ratio) == ("4", 0, 4, 0.0)

    def test_excused_ids_are_read_the_way_they_are_written(self, tmp_path):
        """``load_excused`` reads what CollectFilteredSequences writes, so both
        must apply the same rule. A label id the writer now preserves would have
        been truncated here, failing to excuse the sequence it names and tripping
        the guard on a loss that was already accounted for."""
        path = tmp_path / "filtered_sequences_ids.txt"
        path.write_text("PV547761.1\nKX148218\ncluster.1\n", encoding="utf-8")

        assert load_excused(str(path)) == {"PV547761", "KX148218", "cluster.1"}

    def test_a_healthy_bare_segment_still_passes(self, tmp_path):
        """The no-op proof for this script."""
        queries = [f"KX1482{i:02d}" for i in range(20)]
        merged, tophit, ref_seg, filtered = self._write(
            tmp_path,
            ref_seg_rows=[("NC_001542", "4")],
            merged_headers=["NC_001542"] + queries,
            tophit_rows=[(q, "NC_001542", "99.0", "plus", "4") for q in queries],
        )

        assert check(merged, tophit, ref_seg, filtered, min_retention=0.98) == []


# ---------------------------------------------------------------------------
# projectability - backbone headers decide which queries survive.
# ---------------------------------------------------------------------------


class TestProjectability:

    def test_backbone_ids_lose_a_version_but_not_a_label(self, tmp_path):
        """These ids are membership-tested against ``query_aln/<ref>`` directory
        names. Truncating one at a dot it does not own means its whole query
        group is declared unprojectable and deleted - the exact failure mode this
        module was written to prevent."""
        fasta = tmp_path / "refset_4_aln.fasta"
        fasta.write_text(
            ">NC_001542.1 Rabies lyssavirus\nACGT\n"
            ">KX148218\nACGT\n"
            ">cluster.1\nACGT\n",
            encoding="utf-8",
        )

        assert projectability.read_fasta_ids(str(fasta)) == {
            "NC_001542", "KX148218", "cluster.1",
        }
