"""End-to-end guard for the GenBank+GISAID influenza path.

This file exists because the unit suite could not have caught the bug that broke
a real build. Every other module here is a hermetic single-script test on 2-6
column synthetic fixtures; nothing ran ``merge_into_gB_matrix.py``'s output into
``ValidateSegment.py`` and on into ``CreateSqliteDB.py``. The only code path that
produces both a canonical ``segment`` and a vendor ``Segment`` is exactly that
chain, and it had no test at either end - nor in CI, where ``segmented_test`` was
commented out of the Nextflow matrix.

The production failure being guarded:

    ERROR: duplicate column name: Segment

raised after the trees were already built, because GISAID ships a declared
``Segment`` column of gene names against the pipeline's canonical numeric
``segment``, and SQLite identifiers are case-insensitive.
"""

import csv
import subprocess
import sys

import pandas as pd
import pytest

import segment_utils
from CreateSqliteDB import CreateSqliteDB
from ValidateSegment import annotate_matrix, build_reference_map, build_segment_map
from merge_into_gB_matrix import NormalizeAndMerge

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "scripts"
FLU_SEGMENT_NAMES = REPO_ROOT / "generic" / "influenza" / "segment_names_iav.tsv"

#: The four canonical columns GISAID's headers collide with, case-insensitively.
#: Four are neutralised by the mapping file; `Segment` was the one with no entry.
BASE_COLUMNS = [
    "primary_accession", "accession_version", "locus", "gi_number",
    "segment", "host", "collection_date", "authors", "update_date",
    "sequence", "real_length", "dataset_source", "exclusion", "accession_type",
]


def _write_tsv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    return path


@pytest.fixture
def gisaid_bundle(tmp_path):
    """A GenBank base matrix plus a GISAID export shaped like the real one."""
    base = _write_tsv(
        tmp_path / "gB_matrix.tsv", BASE_COLUMNS,
        [["NC_007373", "NC_007373.1", "NC_007373", "NC_007373", "1", "Human",
          "2005-01-01", "Someone", "2005-01-02", "ATGC", "4", "genbank", "", "master"]],
    )

    # Field names taken verbatim from test_data/gisaid_data/metadata.tsv, including
    # the space-bearing "<SEG> INSDC_Upload" headers and the NA-segment row.
    gisaid = _write_tsv(
        tmp_path / "metadata.tsv",
        ["Segment_Id", "Segment", "Host", "Collection_Date", "Genotype",
         "HA INSDC_Upload", "Zip_Code"],
        [
            ["EPI_HA", "HA", "Human", "2020-01-01", "clade1", "CY000001", "G1"],
            ["EPI_NA", "NA", "Human", "2020-01-02", "clade1", "", "G1"],
            ["EPI_PB1", "PB1", "Human", "2020-01-03", "clade1", "", "G1"],
            ["EPI_PB2", "PB2", "Human", "2020-01-04", "clade1", "", "G1"],
        ],
    )

    fasta = tmp_path / "all_nuc.fas"
    fasta.write_text(
        "".join(f">{acc}\nATGCATGC\n" for acc in ("EPI_HA", "EPI_NA", "EPI_PB1", "EPI_PB2")),
        encoding="utf-8",
    )

    mapping = _write_tsv(
        tmp_path / "column_mapping.tsv", ["source", "target"],
        [["Host", "host"], ["Collection_Date", "collection_date"],
         ["Segment", "segment_declared"]],
    )
    # read_mapping() has no header concept, so drop the header row we just wrote.
    lines = mapping.read_text(encoding="utf-8").splitlines()[1:]
    mapping.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {"base": base, "gisaid": gisaid, "fasta": fasta, "mapping": mapping,
            "out": tmp_path / "gB_matrix_merged.tsv", "tmp": tmp_path}


def _merge(bundle):
    NormalizeAndMerge(
        gb_matrix=str(bundle["base"]),
        input_tsv=str(bundle["gisaid"]),
        input_fasta=str(bundle["fasta"]),
        mapping_file=str(bundle["mapping"]),
        output=str(bundle["out"]),
        key="Segment_Id",
        dataset_source="gisaid",
        log_file=str(bundle["tmp"] / "merge.log"),
    ).process()
    return bundle["out"]


class TestMergeProducesADbSafeMatrix:
    def test_no_two_columns_differ_only_by_case(self, gisaid_bundle):
        """The exact production failure, guarded at the merge boundary."""
        header = list(pd.read_csv(_merge(gisaid_bundle), sep="\t", dtype=str, nrows=0).columns)
        lowered = [c.casefold() for c in header]
        assert len(lowered) == len(set(lowered)), (
            f"case-collision survived the merge: "
            f"{[c for c in header if lowered.count(c.casefold()) > 1]}"
        )

    def test_the_db_writer_guard_accepts_the_merged_matrix(self, gisaid_bundle):
        header = list(pd.read_csv(_merge(gisaid_bundle), sep="\t", dtype=str, nrows=0).columns)
        CreateSqliteDB.assert_no_case_insensitive_duplicates(header, "merged meta_data")

    def test_the_merged_matrix_actually_reaches_sqlite(self, gisaid_bundle):
        """Before the fix this raised OperationalError: duplicate column name."""
        import sqlite3

        frame = pd.read_csv(_merge(gisaid_bundle), sep="\t", dtype=str, keep_default_na=False)
        conn = sqlite3.connect(":memory:")
        frame.to_sql("meta_data", conn, if_exists="replace", index=False)
        assert conn.execute("SELECT COUNT(*) FROM meta_data").fetchone()[0] == 5
        conn.close()

    def test_declared_segment_is_kept_under_its_canonical_name(self, gisaid_bundle):
        rows = pd.read_csv(_merge(gisaid_bundle), sep="\t", dtype=str, keep_default_na=False)
        assert "segment_declared" in rows.columns
        assert "Segment" not in rows.columns
        declared = rows.set_index("primary_accession")["segment_declared"].to_dict()
        assert declared["EPI_HA"] == "HA"
        assert declared["EPI_PB1"] == "PB1"

    def test_neuraminidase_is_not_read_as_a_null(self, gisaid_bundle):
        """`NA` is a segment name. pandas' default sentinels erase it."""
        rows = pd.read_csv(_merge(gisaid_bundle), sep="\t", dtype=str, keep_default_na=False)
        declared = rows.set_index("primary_accession")["segment_declared"].to_dict()
        assert declared["EPI_NA"] == "NA", "neuraminidase was nulled on read"

    def test_spaced_vendor_headers_become_legal_identifiers(self, gisaid_bundle):
        header = list(pd.read_csv(_merge(gisaid_bundle), sep="\t", dtype=str, nrows=0).columns)
        assert "HA INSDC_Upload" not in header
        assert "gisaid_ha_insdc_upload" in header
        assert not [c for c in header if " " in c], f"spaces survive in {header}"

    def test_unmapped_vendor_columns_are_namespaced_not_dropped(self, gisaid_bundle):
        header = list(pd.read_csv(_merge(gisaid_bundle), sep="\t", dtype=str, nrows=0).columns)
        assert "gisaid_genotype" in header and "gisaid_zip_code" in header

    def test_mapped_columns_are_not_namespaced(self, gisaid_bundle):
        header = list(pd.read_csv(_merge(gisaid_bundle), sep="\t", dtype=str, nrows=0).columns)
        assert "host" in header and "collection_date" in header
        assert "gisaid_host" not in header


class TestSegmentIsResolvedForGisaidRows:
    """Before the fix, 65% of a real merged matrix left VALIDATE_SEGMENT with
    `segment` empty and had to be backfilled by a hand-written rescue script."""

    def _run(self, tmp_path, merged, blast_rows):
        blast = tmp_path / "annotated_hits.tsv"
        with open(blast, "w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, delimiter="\t").writerows(blast_rows)
        out = tmp_path / "validated.tsv"
        annotate_matrix(
            str(merged), build_segment_map(str(blast)), build_reference_map(str(blast)),
            str(out), overwrite=False,
            segment_names_path=str(FLU_SEGMENT_NAMES),
        )
        rows = list(csv.DictReader(open(out, encoding="utf-8"), delimiter="\t"))
        return {r["primary_accession"]: r for r in rows}

    def test_declared_names_resolve_to_the_right_segment_numbers(self, gisaid_bundle):
        """With no BLAST hit at all, the submitter's declaration carries the row."""
        merged = _merge(gisaid_bundle)
        rows = self._run(gisaid_bundle["tmp"], merged, [])

        assert rows["EPI_HA"]["segment"] == "4"
        assert rows["EPI_NA"]["segment"] == "6"
        # The trap: digit-scraping made PB1 segment 1 and PB2 segment 2, inverted.
        assert rows["EPI_PB1"]["segment"] == "2"
        assert rows["EPI_PB2"]["segment"] == "1"

    def test_resolution_records_where_the_value_came_from(self, gisaid_bundle):
        merged = _merge(gisaid_bundle)
        rows = self._run(gisaid_bundle["tmp"], merged, [])
        assert rows["EPI_HA"]["segment_source"] == "external_declared"
        assert rows["NC_007373"]["segment_source"] == "source_declared"

    def test_blast_inference_takes_precedence_over_a_declaration(self, gisaid_bundle):
        merged = _merge(gisaid_bundle)
        rows = self._run(gisaid_bundle["tmp"], merged,
                         [["EPI_HA", "REF_HA", "99", "+", "4"]])
        assert rows["EPI_HA"]["segment"] == "4"
        assert rows["EPI_HA"]["segment_source"] == "blast_inferred"

    def test_segment_name_is_derived_for_every_resolved_row(self, gisaid_bundle):
        merged = _merge(gisaid_bundle)
        rows = self._run(gisaid_bundle["tmp"], merged, [])
        assert rows["EPI_HA"]["segment_name"] == "HA"
        assert rows["EPI_NA"]["segment_name"] == "NA"
        assert rows["EPI_PB1"]["segment_name"] == "PB1"
        assert rows["NC_007373"]["segment_name"] == "PB2"  # segment 1

    def test_the_key_column_never_holds_a_gene_name(self, gisaid_bundle):
        """SegmentPivotTable treats a segment value of 'na' as "no segment".

        If gene names were used as identifiers, every neuraminidase row would be
        silently dropped from the completeness table and every flu isolate would
        report Incomplete.
        """
        merged = _merge(gisaid_bundle)
        rows = self._run(gisaid_bundle["tmp"], merged, [])
        for accession, row in rows.items():
            segment = row["segment"]
            if segment:
                assert segment.isdigit(), f"{accession} has a non-numeric segment {segment!r}"

    def test_an_untranslatable_declaration_does_not_leak_into_the_key(self, tmp_path, gisaid_bundle):
        """A vendor vocabulary we have no mapping for must not become the key."""
        merged = _merge(gisaid_bundle)
        frame = pd.read_csv(merged, sep="\t", dtype=str, keep_default_na=False)
        frame.loc[frame["primary_accession"] == "EPI_HA", "segment_declared"] = "WEIRD_SEG"
        patched = tmp_path / "patched.tsv"
        frame.to_csv(patched, sep="\t", index=False)

        rows = self._run(tmp_path, patched, [])
        assert rows["EPI_HA"]["segment"] == ""
        assert rows["EPI_HA"]["segment_source"] == ""


class TestValidateSegmentCli:
    def test_cli_accepts_the_segment_names_asset(self, tmp_path, gisaid_bundle):
        """Guards the wiring vgtk-init.nf depends on."""
        merged = _merge(gisaid_bundle)
        blast = tmp_path / "hits.tsv"
        blast.write_text("EPI_HA\tREF\t99\t+\t4\n", encoding="utf-8")
        out = tmp_path / "cli_out.tsv"

        result = subprocess.run(
            [sys.executable, str(SCRIPTS / "ValidateSegment.py"),
             "-g", str(merged), "-s", str(blast), "-o", str(out),
             "--segment_names", str(FLU_SEGMENT_NAMES)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert out.is_file()
        assert "Segment resolution:" in result.stdout
