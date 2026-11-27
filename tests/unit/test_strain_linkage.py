"""Extensive TDD unit and integration tests for Workstream 2:
Segmented virus strain-name linkage table (`strain_segment_linkage`) and
AIV/IAV + general segmented strain canonicalisation.
"""

from __future__ import annotations

import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from strain_linkage import (
    build_strain_segment_linkage,
    load_location_aliases,
    normalise_segmented_strain_key,
    write_strain_linkage_tables,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATE_DB_TREE = REPO_ROOT / "scripts" / "ValidateDbTree.py"
STRAIN_LOC_ALIASES = REPO_ROOT / "generic" / "influenza" / "strain_location_aliases.tsv"


class TestSegmentedStrainNormalisation:
    def test_aiv_real_world_spellings_normalise_identically(self):
        aliases = load_location_aliases(str(STRAIN_LOC_ALIASES)) if STRAIN_LOC_ALIASES.exists() else {"ny": "new_york"}
        # Real EPI_ISL_1089 variants observed in live iav-db.db
        variants = [
            ("A/Chicken/New_York/6777-3/97", "1997"),
            ("A/Chicken/New_York/6777-3/97_H7N2NSB", "1997"),
            ("A/chicken/NY/6777-3/1997", "1997"),
            ("A/Chicken/NY/6777-3/1997", "1997"),
            ("A/Chicken/New York/6777-3/1997 (H7N2)", "1997"),
        ]
        keys = {normalise_segmented_strain_key(s, collection_year=y, location_aliases=aliases) for s, y in variants}
        assert len(keys) == 1
        assert next(iter(keys)) == "a/chicken/new_york/6777-3/1997"

    def test_prefixes_suffixes_double_slashes_and_two_digit_years(self):
        assert (
            normalise_segmented_strain_key("A/swine//Ohio/02026/2008", collection_year="2008")
            == "a/swine/ohio/02026/2008"
        )
        assert (
            normalise_segmented_strain_key("_A/duck/Vietnam/12/05-original", collection_year="2005")
            == "a/duck/vietnam/12/2005"
        )
        assert (
            normalise_segmented_strain_key("partial/A/mallard/Alaska/44/2019-HA", collection_year="2019")
            == "a/mallard/alaska/44/2019"
        )
        assert (
            normalise_segmented_strain_key("Recombinant_A/PR/8/34(H1N1)", collection_year="1934")
            == "a/pr/8/1934"
        )
        assert (
            normalise_segmented_strain_key("A/barnacle_goose/Sweden/SVA231107SZ0383/FB004550/H-2023", collection_year="2023")
            == "a/barnacle_goose/sweden/sva231107sz0383/fb004550/h-2023"
        )

    def test_non_flu_segmented_virus_strain_preserved(self):
        assert normalise_segmented_strain_key("LASV/Josiah/SL/1976") == "lasv/josiah/sl/1976"
        assert normalise_segmented_strain_key("Unknown") == ""
        assert normalise_segmented_strain_key("N/A") == ""


class TestStrainSegmentLinkageGrouping:
    def _make_row(self, acc, seg, strain="", iso_id="", alt="", src="genbank", passage="", year="2024", rlen=1700, excl=""):
        return {
            "primary_accession": acc,
            "segment": str(seg),
            "Parsed_strain": strain,
            "strain": strain,
            "isolate": "",
            "Isolate_Id": iso_id,
            "alt_accessions": alt,
            "dataset_source": src,
            "Passage_History": passage,
            "collection_year": year,
            "collection_date": f"{year}-01-01" if year else "",
            "real_length": str(rlen),
            "exclusion_status": excl,
            "exclusion_criteria": "excluded" if excl == "1" else "",
            "accession_type": "query",
        }

    def test_reunites_intra_db_split_isolate_ids_like_abudhabi(self):
        # Mirrors real A/ABUDHABI/2373/2025: segs 1,2,4,5,6,7,8 have EPI_ISL_20133490,
        # seg 3 won by EPI_ISL_20178480 with EPI_ISL_20133490 in alt_accessions.
        rows = []
        for s in (1, 2, 4, 5, 6, 7, 8):
            rows.append(
                self._make_row(
                    f"EPI_A_{s}", s, "A/ABUDHABI/2373/2025",
                    iso_id="EPI_ISL_20133490",
                    alt=f"gisaid:EPI_B_{s};EPI_ISL_20178480",
                    src="gisaid", passage="cs", year="2025",
                )
            )
        rows.append(
            self._make_row(
                "EPI_B_3", 3, "A/ABUDHABI/2373/2025",
                iso_id="EPI_ISL_20178480",
                alt="gisaid:EPI_A_3;EPI_ISL_20133490",
                src="gisaid", passage="SIAT1", year="2025",
            )
        )
        linkage, alternates, seg_order = build_strain_segment_linkage(rows, expected_segments=["1","2","3","4","5","6","7","8"])
        assert len(linkage) == 1
        rec = linkage[0]
        assert rec["complete_status"] == "Complete"
        assert rec["segments_present"] == 8
        assert rec["segment_3"] == "EPI_B_3"
        assert rec["segment_4"] == "EPI_A_4"
        assert "EPI_ISL_20178480" in rec["alt_isolate_ids"] or "EPI_ISL_20133490" in rec["alt_isolate_ids"]
        assert alternates == []

    def test_keeps_distinct_passaged_genomes_separate_like_aichi(self):
        # Mirrors real A/AICHI/2/2018: MDCK 1 (EPI_ISL_312775) vs Original (EPI_ISL_312776),
        # non-identical sequences so NOT in each other's alt_accessions.
        rows = []
        for s in range(1, 9):
            rows.append(self._make_row(f"EPI_M_{s}", s, "A/AICHI/2/2018", iso_id="EPI_ISL_312775", src="gisaid", passage="MDCK 1", year="2018"))
            rows.append(self._make_row(f"EPI_O_{s}", s, "A/AICHI/2/2018", iso_id="EPI_ISL_312776", src="gisaid", passage="Original", year="2018"))
        linkage, alternates, _ = build_strain_segment_linkage(rows, expected_segments=[str(i) for i in range(1, 9)])
        assert len(linkage) == 2
        by_iso = {r["isolate_id"]: r for r in linkage}
        assert by_iso["EPI_ISL_312775"]["complete_status"] == "Complete"
        assert by_iso["EPI_ISL_312775"]["segment_4"] == "EPI_M_4"
        assert by_iso["EPI_ISL_312776"]["complete_status"] == "Complete"
        assert by_iso["EPI_ISL_312776"]["segment_4"] == "EPI_O_4"
        assert alternates == []

    def test_recovers_missing_a_prefix_and_merges_isolate_id_less_genbank_segments(self):
        # Mirrors EPI_ISL_1080 (missing A/) + partial GenBank without Isolate_Id filling remaining segments
        rows = [
            self._make_row("EPI_1", 1, "A/guinea_fowl/New_York/13820-3/1995", iso_id="EPI_ISL_1080", src="gisaid", year="1995"),
            self._make_row("PQ_2", 2, "Guinea_Fowl/New_York/13820-3/95_H7N2NSB", iso_id="", src="genbank", year="1995"),
            self._make_row("PQ_3", 3, "A/Guinea_Fowl/NY/13820-3/1995", iso_id="", src="genbank", year="1995"),
        ]
        linkage, _, _ = build_strain_segment_linkage(
            rows,
            expected_segments=["1", "2", "3"],
            location_aliases={"ny": "new_york"},
        )
        assert len(linkage) == 1
        rec = linkage[0]
        assert rec["complete_status"] == "Complete"
        assert rec["segment_1"] == "EPI_1"
        assert rec["segment_2"] == "PQ_2"
        assert rec["segment_3"] == "PQ_3"
        assert rec["dataset_sources"] == "genbank;gisaid"

    def test_ambiguous_multi_sequence_cell_picks_deterministic_winner_and_records_alternate(self):
        # Mirrors EPI_ISL_10186 having two different kept GenBank sequences on segment 1
        rows = [
            self._make_row("DQ492837", 1, "A/duck/Test/1/2006", iso_id="EPI_ISL_10186", src="genbank", year="2006", rlen=2280),
            self._make_row("DQ492833", 1, "A/duck/Test/1/2006", iso_id="EPI_ISL_10186", src="genbank", year="2006", rlen=2200),
            self._make_row("DQ493362", 2, "A/duck/Test/1/2006", iso_id="EPI_ISL_10186", src="genbank", year="2006", rlen=2274),
        ]
        linkage, alternates, _ = build_strain_segment_linkage(rows, expected_segments=["1", "2"])
        assert len(linkage) == 1
        rec = linkage[0]
        assert rec["segment_1"] == "DQ492837"  # longer sequence wins deterministically
        assert rec["segment_2"] == "DQ493362"
        assert rec["ambiguous_segments"] == "1"
        assert len(alternates) == 1
        assert alternates[0]["alternate_accession"] == "DQ492833"
        assert alternates[0]["selected_accession"] == "DQ492837"
        assert alternates[0]["segment"] == "1"

    def test_excludes_excluded_rows_and_handles_letter_segments_and_blank_singletons(self):
        rows = [
            self._make_row("LASV_L", "L", "LASV/Josiah/1976", src="genbank", year="1976"),
            self._make_row("LASV_S", "S", "LASV/Josiah/1976", src="genbank", year="1976"),
            self._make_row("EXCL_S", "S", "LASV/Josiah/1976", src="genbank", year="1976", excl="1"),
            self._make_row("ORPHAN_L", "L", "", src="genbank", year="1976"),
        ]
        linkage, alternates, seg_order = build_strain_segment_linkage(rows, expected_segments=["L", "S"])
        assert seg_order == ["L", "S"]
        assert len(linkage) == 2
        by_acc = {r["strain_accession"]: r for r in linkage}
        josiah = [r for r in linkage if r["canonical_strain_key"] == "lasv/josiah/1976"][0]
        assert josiah["segment_L"] == "LASV_L"
        assert josiah["segment_S"] == "LASV_S"
        assert josiah["complete_status"] == "Complete"
        assert "acc:ORPHAN_L" in by_acc
        assert by_acc["acc:ORPHAN_L"]["segment_L"] == "ORPHAN_L"
        assert by_acc["acc:ORPHAN_L"]["complete_status"] == "Incomplete"
        assert alternates == []

    def test_write_and_validate_db_tree_linkage_invariant(self, tmp_path: Path):
        db = tmp_path / "seg_linkage.db"
        conn = sqlite3.connect(db)
        conn.executescript(
            """
            CREATE TABLE meta_data (
                primary_accession TEXT, accession_type TEXT, dataset_source TEXT,
                segment TEXT, Parsed_strain TEXT, strain TEXT, isolate TEXT,
                Isolate_Id TEXT, alt_accessions TEXT, Passage_History TEXT,
                collection_year TEXT, collection_date TEXT, real_length TEXT,
                exclusion_status TEXT, exclusion_criteria TEXT, host_taxa_id TEXT
            );
            CREATE TABLE sequence_alignment (primary_accession TEXT, segment TEXT, alignment TEXT);
            CREATE TABLE sequences (header TEXT, segment TEXT, sequence TEXT);
            CREATE TABLE features (accession TEXT, segment TEXT, product TEXT, cds_start INT, cds_end INT);
            CREATE TABLE host_taxa (taxa_id TEXT, name TEXT, name_type TEXT, taxonomy_level TEXT, common_name TEXT);
            CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT);
            """
        )
        conn.execute(
            "INSERT INTO meta_data VALUES ('A1','query','genbank','1','A/Mallard/NY/1/2024','A/Mallard/NY/1/2024','','','','Original','2024','2024-01-01','2280','','','8839')"
        )
        conn.execute(
            "INSERT INTO meta_data VALUES ('A2','query','genbank','2','A/Mallard/New_York/1/24','A/Mallard/New_York/1/24','','','','Original','2024','2024-01-01','2274','','','8839')"
        )
        conn.execute("INSERT INTO host_taxa VALUES ('8839','Anas platyrhynchos','scientific','species','mallard')")
        conn.execute("INSERT INTO sequence_alignment VALUES ('A1','1','ATGCGT')")
        conn.execute("INSERT INTO sequence_alignment VALUES ('A2','2','ATGCGT')")
        conn.execute("INSERT INTO sequences VALUES ('A1','1','ATGCGT')")
        conn.execute("INSERT INTO sequences VALUES ('A2','2','ATGCGT')")
        conn.execute("INSERT INTO features VALUES ('A1','1','PB2',1,6)")
        conn.execute("INSERT INTO features VALUES ('A2','2','PB1',1,6)")
        conn.execute("INSERT INTO trees VALUES ('usher_seg1','usher','segment_1','1','(A1:0.01);','2026-11-25')")
        conn.execute("INSERT INTO trees VALUES ('usher_seg2','usher','segment_2','2','(A2:0.01);','2026-11-25')")
        conn.commit()

        summary = write_strain_linkage_tables(conn, expected_segments=["1", "2"], location_aliases={"ny": "new_york"})
        conn.commit()
        assert summary["strain_rows"] == 1
        assert summary["complete_strains"] == 1
        row = conn.execute("SELECT segment_1, segment_2, complete_status FROM strain_segment_linkage").fetchone()
        assert row == ("A1", "A2", "Complete")
        conn.close()

        # ValidateDbTree must pass on clean linkage table
        res_ok = subprocess.run(
            [sys.executable, str(VALIDATE_DB_TREE), "--db", str(db), "--outdir", str(tmp_path / "out_ok"), "--strict-invariants"],
            capture_output=True, text=True,
        )
        assert res_ok.returncode == 0, res_ok.stdout + res_ok.stderr

        # Corrupt segment_2 to point to wrong-segment accession A1 -> ValidateDbTree must FAIL
        conn = sqlite3.connect(db)
        conn.execute("UPDATE strain_segment_linkage SET segment_2 = 'A1'")
        conn.commit()
        conn.close()
        res_bad = subprocess.run(
            [sys.executable, str(VALIDATE_DB_TREE), "--db", str(db), "--outdir", str(tmp_path / "out_bad"), "--strict-invariants"],
            capture_output=True, text=True,
        )
        assert res_bad.returncode != 0
        assert "strain_segment_linkage" in (res_bad.stdout + res_bad.stderr).lower()
