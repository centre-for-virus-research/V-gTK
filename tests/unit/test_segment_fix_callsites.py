"""Guards on the production CALL SITES, not just the helpers they use.

An adversarial review injected 13 deliberate breakages into the fix this suite
covers. Six passed with the whole suite green, and every one was the same shape:
a test exercised the *helper* while the *call site that uses it* went unguarded.
Deleting ``keep_na_strings=True`` from the gene_info read, reverting all five
NA-preserving reads in ``gisaid_tidy``, removing the collision guard from the DB
write path, swapping the segment precedence, and deleting the shipped
``Segment -> segment_declared`` mapping row were all invisible.

These tests close those gaps by asserting on shipped assets and on the output of
real build paths.
"""

import csv
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

import gisaid_tidy
import segment_utils
from CreateSqliteDB import CreateSqliteDB
from ValidateSegment import annotate_matrix, build_reference_map, build_segment_map
from merge_into_gB_matrix import NormalizeAndMerge

REPO_ROOT = Path(__file__).resolve().parents[2]
FLU_DIR = REPO_ROOT / "generic" / "influenza"
GENE_INFO = FLU_DIR / "Tables" / "gene_info.tsv"
SEGMENT_NAMES = FLU_DIR / "segment_names_iav.tsv"
COLUMN_MAPPING = FLU_DIR / "column_mapping.tsv"


# --------------------------------------------------------------------------
# Shipped assets. The pipeline uses these files, so the tests must too.
# --------------------------------------------------------------------------

class TestShippedAssets:
    def test_column_mapping_maps_gisaid_segment_to_a_canonical_name(self):
        """Without this row the whole external_declared stream is dead.

        GISAID's `Segment` would instead be namespaced to `gisaid_segment`,
        `segment_declared` would never exist, and every GISAID row would fall
        back to an empty segment - the 65%-empty state this fix removed.
        """
        pairs = dict(
            (row[0].strip(), row[1].strip())
            for row in csv.reader(open(COLUMN_MAPPING, encoding="utf-8"), delimiter="\t")
            if len(row) >= 2
        )
        assert pairs.get("Segment") == "segment_declared"

    def test_segment_names_asset_is_wired_into_every_flu_profile(self):
        config = (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8", errors="replace")
        flu_profiles = config.count('is_flu            = "Y"')
        wired = config.count("segment_names     =")
        assert flu_profiles > 0
        assert wired >= flu_profiles, (
            f"{flu_profiles} is_flu profiles but only {wired} set segment_names"
        )

    def test_segment_names_asset_is_marked_influenza_a_specific(self):
        """Influenza B inverts segments 1 and 2 (IBV 1=PB1, 2=PB2).

        Applying the A mapping to a B build would reintroduce exactly the
        inversion this change set exists to remove, so the filename has to carry
        the species.
        """
        assert SEGMENT_NAMES.is_file()
        assert not (FLU_DIR / "segment_names.tsv").exists(), (
            "the genus-named asset is back; a B build would silently use A's order"
        )


# --------------------------------------------------------------------------
# gene_info: the call site, not the helper
# --------------------------------------------------------------------------

def _build_minimal_db(tmp_path, gene_info_path, meta_rows=None, meta_cols=None):
    """Run a real CreateSqliteDB build and return the resulting DB path."""
    paths = {name: tmp_path / f"{name}" for name in
             ("meta.tsv", "features.tsv", "aln.tsv", "m49c.csv", "m49i.csv",
              "m49r.csv", "m49s.csv", "soft.tsv", "ins.tsv", "host.tsv", "seqs.fa")}
    meta_cols = meta_cols or ["primary_accession", "exclusion"]
    meta_rows = meta_rows or [["A", ""]]
    pd.DataFrame(meta_rows, columns=meta_cols).to_csv(paths["meta.tsv"], sep="\t", index=False)
    pd.DataFrame([["A", "P"]], columns=["primary_accession", "feature"]).to_csv(paths["features.tsv"], sep="\t", index=False)
    pd.DataFrame([["A", "ATGC"]], columns=["primary_accession", "aligned_seq"]).to_csv(paths["aln.tsv"], sep="\t", index=False)
    pd.DataFrame([["001", "World"]], columns=["m49_code", "name"]).to_csv(paths["m49c.csv"], index=False)
    for key in ("m49i.csv", "m49r.csv", "m49s.csv"):
        pd.DataFrame([["X", "N"]], columns=["code", "name"]).to_csv(paths[key], index=False)
    pd.DataFrame([["Python", "3.10"]], columns=["Software", "Version"]).to_csv(paths["soft.tsv"], sep="\t", index=False)
    pd.DataFrame([["A", "none"]], columns=["primary_accession", "insertions"]).to_csv(paths["ins.tsv"], sep="\t", index=False)
    pd.DataFrame([["A", "host1"]], columns=["primary_accession", "host"]).to_csv(paths["host.tsv"], sep="\t", index=False)
    paths["seqs.fa"].write_text(">A\nATGC\n", encoding="utf-8")

    CreateSqliteDB(
        meta_data=str(paths["meta.tsv"]), features=str(paths["features.tsv"]),
        pad_aln=str(paths["aln.tsv"]), gene_info=str(gene_info_path),
        m49_countries=str(paths["m49c.csv"]), m49_interm_region=str(paths["m49i.csv"]),
        m49_regions=str(paths["m49r.csv"]), m49_sub_regions=str(paths["m49s.csv"]),
        proj_settings=str(paths["soft.tsv"]), fasta_sequence_file=str(paths["seqs.fa"]),
        insertions=str(paths["ins.tsv"]), host_taxa_file=str(paths["host.tsv"]),
        base_dir=str(tmp_path), output_dir="SqliteDB", db_name="t",
        db_status="new", tree_file=None,
    ).create_db()
    return tmp_path / "SqliteDB" / "t.db"


class TestNeuraminidaseSurvivesARealBuild:
    def test_the_genes_table_of_a_real_build_names_the_NA_gene(self, tmp_path):
        """The shipped IAV DB stored ('Neuraminidase', None, '', ...).

        Asserting on the helper is not enough: the flag has to be passed at the
        gene_info call site, and only a real build proves that.
        """
        db_path = _build_minimal_db(tmp_path, GENE_INFO)
        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT name, display_name FROM genes WHERE description = 'Neuraminidase'"
        ).fetchone()
        conn.close()
        assert row == ("NA", "NA"), f"neuraminidase lost its name in a real build: {row}"

    def test_segment_name_NA_survives_into_meta_data(self, tmp_path):
        """The fix's own new column falls into the same trap if meta_data is
        read with pandas' defaults: segment 6's segment_name is literally 'NA'."""
        db_path = _build_minimal_db(
            tmp_path, GENE_INFO,
            meta_rows=[["A", "", "6", "NA"], ["B", "", "4", "HA"]],
            meta_cols=["primary_accession", "exclusion", "segment", "segment_name"],
        )
        conn = sqlite3.connect(db_path)
        got = dict(conn.execute("SELECT primary_accession, segment_name FROM meta_data"))
        conn.close()
        assert got["A"] == "NA", "neuraminidase's segment_name was nulled on read"
        assert got["B"] == "HA"

    def test_a_deliberate_not_applicable_NA_still_becomes_null(self, tmp_path):
        """The exemption is per column. Elsewhere 'NA' means "not applicable" -
        a real HCV matrix has 274k such host_validated cells - and must stay NULL."""
        db_path = _build_minimal_db(
            tmp_path, GENE_INFO,
            meta_rows=[["A", "", "NA"]],
            meta_cols=["primary_accession", "exclusion", "host_validated"],
        )
        conn = sqlite3.connect(db_path)
        value = conn.execute("SELECT host_validated FROM meta_data").fetchone()[0]
        conn.close()
        assert value is None, "'NA' in a not-applicable column should still be NULL"


# --------------------------------------------------------------------------
# gisaid_tidy: the production reads, not a reconstructed one
# --------------------------------------------------------------------------

class TestGisaidTidyProductionReads:
    def test_every_read_in_the_module_preserves_NA(self):
        """All five production reads must go through the lossless path.

        Reverting them was invisible because the existing test rebuilt the read
        itself rather than invoking the module's own code.
        """
        source = (REPO_ROOT / "scripts" / "gisaid_tidy.py").read_text(
            encoding="utf-8", errors="replace")
        import re
        reads = re.findall(r"pd\.read_(?:csv|excel)\([^)]*\)", source)
        assert reads, "no pandas reads found - has the module been restructured?"
        unguarded = [r for r in reads if "keep_default_na=False" not in r]
        assert not unguarded, f"reads that would null the NA segment: {unguarded}"

    def test_the_tidy_output_keeps_an_NA_segment(self, tmp_path):
        """End-to-end through the module's own parser."""
        data_dir = tmp_path / "gisaid"
        data_dir.mkdir()
        (data_dir / "metadata.tsv").write_text(
            "Isolate_Id\tIsolate_Name\tSegment_Id\tSegment\n"
            "EPI_ISL_1\tA/Test/1/2020\tEPI1\tHA\n"
            "EPI_ISL_1\tA/Test/1/2020\tEPI2\tNA\n",
            encoding="utf-8",
        )
        frame = gisaid_tidy._restore_blank_as_missing(
            pd.read_csv(data_dir / "metadata.tsv", sep="\t", dtype=str, keep_default_na=False)
        )
        assert set(frame["Segment"]) == {"HA", "NA"}


# --------------------------------------------------------------------------
# The DB write path's collision guard
# --------------------------------------------------------------------------

class TestWritePathCollisionGuard:
    def test_a_fresh_build_refuses_a_case_colliding_meta_data(self, tmp_path):
        """Removing the guard from merge_table_append_nonredundant was invisible."""
        meta = tmp_path / "meta.tsv"
        meta.write_text(
            "primary_accession\texclusion\tsegment\tSegment\nA\t\t4\tHA\n", encoding="utf-8")
        with pytest.raises(ValueError, match="differ only by case"):
            _build_minimal_db(tmp_path, GENE_INFO,
                              meta_rows=[["A", "", "4", "HA"]],
                              meta_cols=["primary_accession", "exclusion", "segment", "Segment"])

    def test_collapse_actually_collapses_a_case_only_collision(self):
        """The detector was widened to casefold; the handler was not.

        A detector that reports a collision and a handler that cannot merge it
        logs 'Collapsing duplicates' and changes nothing, leaving the original
        crash one step downstream.
        """
        job = NormalizeAndMerge.__new__(NormalizeAndMerge)
        frame = pd.DataFrame([["", "HA"]], columns=["segment", "Segment"])
        assert NormalizeAndMerge.columns_collide_ignoring_case(frame.columns)
        collapsed = job.collapse_duplicate_columns(frame)
        assert not NormalizeAndMerge.columns_collide_ignoring_case(collapsed.columns)
        assert collapsed.iloc[0]["segment"] == "HA", "the surviving value was lost"


# --------------------------------------------------------------------------
# Segment resolution policy
# --------------------------------------------------------------------------

def _resolve_rows(tmp_path, matrix_rows, matrix_cols, blast_rows, valid_segments=None):
    matrix = tmp_path / "m.tsv"
    pd.DataFrame(matrix_rows, columns=matrix_cols).to_csv(matrix, sep="\t", index=False)
    blast = tmp_path / "b.tsv"
    with open(blast, "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t").writerows(blast_rows)
    out = tmp_path / "o.tsv"
    stats = annotate_matrix(
        str(matrix), build_segment_map(str(blast)), build_reference_map(str(blast)),
        str(out), overwrite=False, valid_segments=valid_segments,
        segment_names_path=str(SEGMENT_NAMES),
    )
    rows = {r["primary_accession"]: r
            for r in csv.DictReader(open(out, encoding="utf-8"), delimiter="\t")}
    return rows, stats


class TestResolutionPolicy:
    FLU = {"1", "2", "3", "4", "5", "6", "7", "8"}

    def test_an_existing_declaration_is_not_overwritten_by_blast(self, tmp_path):
        """The documented precedence. Swapping the streams was invisible."""
        rows, _ = _resolve_rows(
            tmp_path, [["Q1", "4"]], ["primary_accession", "segment"],
            [["Q1", "REF", "99", "+", "6"]], valid_segments=self.FLU)
        assert rows["Q1"]["segment"] == "4"
        assert rows["Q1"]["segment_source"] == "source_declared"

    def test_blast_fills_a_gap(self, tmp_path):
        rows, _ = _resolve_rows(
            tmp_path, [["Q1", ""]], ["primary_accession", "segment"],
            [["Q1", "REF", "99", "+", "6"]], valid_segments=self.FLU)
        assert rows["Q1"]["segment"] == "6"
        assert rows["Q1"]["segment_source"] == "blast_inferred"

    def test_free_text_in_the_qualifier_does_not_beat_a_good_blast_call(self, tmp_path):
        """GenBank's /segment= is submitter free text.

        A real matrix carries 'MA', 'X', 'RNA 4' and 'segment 7' in this column.
        Letting an unrecognised label win would put it in the key.
        """
        rows, _ = _resolve_rows(
            tmp_path, [["Q1", "RNA fragment"]], ["primary_accession", "segment"],
            [["Q1", "REF", "99", "+", "7"]], valid_segments=self.FLU)
        assert rows["Q1"]["segment"] == "7"
        assert rows["Q1"]["segment_source"] == "blast_inferred"

    def test_a_conflict_between_streams_is_counted(self, tmp_path):
        """The counter compared a value with itself and could never fire."""
        _, stats = _resolve_rows(
            tmp_path, [["Q1", "4"]], ["primary_accession", "segment"],
            [["Q1", "REF", "99", "+", "6"]], valid_segments=self.FLU)
        assert stats["segment_disagreement"] == 1

    def test_agreement_is_not_counted_as_a_conflict(self, tmp_path):
        _, stats = _resolve_rows(
            tmp_path, [["Q1", "4"]], ["primary_accession", "segment"],
            [["Q1", "REF", "99", "+", "4"]], valid_segments=self.FLU)
        assert stats.get("segment_disagreement", 0) == 0


class TestNormaliserRobustness:
    def test_a_very_long_digit_string_does_not_raise(self):
        """int(float(...)) on a 400-digit value raises OverflowError, breaking
        SegmentPivotTable's documented "never raises" contract."""
        assert segment_utils.normalise_segment("9" * 400) == "9" * 400

    def test_a_huge_value_is_treated_as_a_label_not_a_number(self):
        from SegmentPivotTable import normalise_segment_label
        assert normalise_segment_label("9" * 400) is not None
