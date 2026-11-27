"""Unit and end-to-end invariant tests for cross-database dual-upload deduplication.

Covers:
  1. ``scripts/cross_db_dedup.py`` helpers (strain/isolate normalization across
     flu & unsegmented viruses, passage compatibility, UTR overhang sequence
     matching, coverage-parity stub protection, date-precision coalescing, and
     iterative ``O(N)`` Newick tip pruning/renaming).
  2. ``scripts/merge_into_gB_matrix.py`` build-time & ``--update_db`` deduplication
     (including ``--base_fasta`` when ``GenBankParser.py`` omits ``sequence`` from
     ``gB_matrix_raw.tsv``, Pathoplexus ``insdcAccessionBase``, uploader INSDC
     typos, and ``FilterAndExtractSequences.py`` tombstone passthrough).
  3. ``scripts/DedupStrainAlignment.py``, ``scripts/UsherPlacement.py`` starter-tree
     pruning, ``scripts/CreateSqliteDB.py`` fresh & reverse-arrival ``--update``, and
     ``scripts/ValidateDbTree.py`` duplicate-strain invariant checks.
"""

from __future__ import annotations

import csv
import sqlite3
import subprocess
import sys
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from Bio import Phylo

from cross_db_dedup import (
    aligned_sequences_agree,
    classify_passage,
    coalesce_winner_metadata,
    extract_canonical_isolate,
    extract_insdc_links,
    normalize_isolate_key,
    passages_compatible,
    prune_and_rename_newick,
    unaligned_sequences_agree,
)
from merge_into_gB_matrix import NormalizeAndMerge

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATE_DB_TREE = REPO_ROOT / "scripts" / "ValidateDbTree.py"
DEDUP_STRAIN_SCRIPT = REPO_ROOT / "scripts" / "DedupStrainAlignment.py"


def _write_tsv(path: Path, header, rows) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    with path.open("w", encoding="utf-8") as fh:
        for hdr, seq in records.items():
            fh.write(f">{hdr}\n{seq}\n")
    return path


# ---------------------------------------------------------------------------
# 1. Domain helper tests
# ---------------------------------------------------------------------------

class TestCrossDbDedupHelpers:
    def test_normalize_isolate_key_handles_flu_and_unsegmented_labels(self):
        assert (
            normalize_isolate_key("A/Wisconsin/67/2022 (H1N1)", collection_year="2022")
            == "a/wisconsin/67/2022"
        )
        assert (
            normalize_isolate_key("A / Michigan / 45 / 15(H1N1)", collection_year="2015")
            == "a/michigan/45/2015"
        )
        # Non-segmented virus isolate (e.g. RABV / Mpox / CCHF) survives intact
        assert normalize_isolate_key("RABV/Dog/Kenya/Nairobi-042/2023") == "rabv/dog/kenya/nairobi-042/2023"
        # Placeholders must never match each other
        assert normalize_isolate_key("Unknown") == ""
        assert normalize_isolate_key("N/A") == ""

    def test_extract_canonical_isolate_falls_back_to_genbank_definition(self):
        row = {
            "strain": "",
            "definition": "Influenza A virus (A/Texas/50/2012(H3N2)) segment 4 hemagglutinin (HA) gene, complete cds",
            "collection_date": "2012-04-10",
        }
        assert extract_canonical_isolate(row) == "a/texas/50/2012"

    def test_passage_compatibility_separates_original_from_egg_and_cell(self):
        assert classify_passage("Original") == "orig"
        assert classify_passage("Clinical specimen") == "orig"
        assert classify_passage("E2") == "egg"
        assert classify_passage("SPF Egg 1") == "egg"
        assert classify_passage("MDCK2/SIAT1") == "cell"
        assert classify_passage("") == "unknown"

        assert not passages_compatible("orig", "egg")
        assert not passages_compatible("orig", "cell")
        assert not passages_compatible("egg", "cell")
        assert passages_compatible("orig", "unknown")
        assert passages_compatible("orig", "orig")

    def test_extract_insdc_links_is_segment_scoped_and_supports_pathoplexus(self):
        flu_row = {
            "segment_declared": "HA",
            "gisaid_pb2_insdc_upload": "PX111111.1",
            "gisaid_ha_insdc_upload": "PX222222.1; PX222223",
        }
        assert extract_insdc_links(flu_row) == ["PX222222", "PX222223"]

        ppx_row = {
            "pathoplexus_insdcaccessionbase": "OQ987654.1",
        }
        assert extract_insdc_links(ppx_row) == ["OQ987654"]

    def test_unaligned_sequences_agree_handles_utr_overhang_and_blocks_stubs(self):
        core = ("ATGCGTAC" * 75)  # 600 bp (>=95% of 624 bp with 24 bp UTR)
        with_utr = "AGCAAAAGCAGG" + core + "CTTGTTTCTACT"  # 624 bp (12 bp 5' + 12 bp 3')
        ok, mm, _ = unaligned_sequences_agree(core, with_utr, max_mm=0)
        assert ok and mm == 0

        # Opposite-end overhangs: GenBank has 5' UTR, GISAID has 3' UTR
        gb_5p = "AGCAAAAGCAGG" + core
        gi_3p = core + "CTTGTTTCTACT"
        ok2, mm2, _ = unaligned_sequences_agree(gb_5p, gi_3p, max_mm=0)
        assert ok2 and mm2 == 0

        # Truncated GenBank stub (300 bp vs 600 bp) must be rejected
        stub = core[:300]
        ok_stub, _, _ = unaligned_sequences_agree(stub, core, max_mm=0)
        assert not ok_stub

    def test_coalesce_winner_metadata_updates_split_date_columns_in_lockstep(self):
        winner = {
            "primary_accession": "PX393448",
            "collection_date": "2024-07",
            "collection_year": "2024",
            "collection_mon": "7",
            "collection_day": "",
            "Originating_Lab": "",
            "alt_accessions": "",
        }
        retired = {
            "primary_accession": "EPI4717132",
            "Isolate_Id": "EPI_ISL_19283741",
            "collection_date": "2024-07-19",
            "Originating_Lab": "WHO Collaborating Centre",
        }
        upd = coalesce_winner_metadata(winner, retired, retired_source="gisaid")
        assert upd["collection_date"] == "2024-07-19"
        assert upd["collection_year"] == "2024"
        assert upd["collection_mon"] == "7"
        assert upd["collection_day"] == "19"
        assert upd["Originating_Lab"] == "WHO Collaborating Centre"
        assert upd["alt_accessions"] == "gisaid:EPI4717132;EPI_ISL_19283741"

        # Contradictory month must never overwrite winner's date
        bad_retired = dict(retired, collection_date="2024-08-19")
        upd_bad = coalesce_winner_metadata(winner, bad_retired, retired_source="gisaid")
        assert "collection_date" not in upd_bad


# ---------------------------------------------------------------------------
# 2. Build-time & --update_db merge_into_gB_matrix tests
# ---------------------------------------------------------------------------

class TestMergeCrossSourceDeduplication:
    def test_merge_deduplicates_with_base_fasta_and_guards_typos_stubs_and_passages(self, tmp_path: Path):
        """Simulates real GenBankParser output (NO `sequence` column in gB_matrix_raw.tsv,
        sequences supplied via --base_fasta)."""
        seq_ha_1 = "ATGC" * 100          # 400 bp exact
        seq_ha_2 = "GCTA" * 100          # 400 bp exact
        seq_ha_2_utr = "AAATTT" + seq_ha_2 + "CCCGGG"  # 412 bp (12 bp UTR overhang)
        seq_typo = "TTTT" * 100          # completely different sequence (uploader typo)
        seq_stub = seq_ha_1[:180]        # 180 bp truncated GenBank stub (<95% of 400 bp)
        seq_egg_mut = ("ATGC" * 99) + "ATGG"  # 1 mismatch + Egg passage

        # Real GenBankParser drops `sequence` from gB_matrix_raw.tsv!
        base_cols = [
            "primary_accession", "accession_version", "locus", "gi_number",
            "segment", "strain", "Passage_History", "collection_date",
            "collection_year", "collection_mon", "collection_day",
            "Originating_Lab", "real_length", "exclusion", "accession_type",
        ]
        base_tsv = _write_tsv(
            tmp_path / "gB_matrix_raw.tsv",
            base_cols,
            [
                # 1. Explicit INSDC target with year-only date (should be enriched to 2024-03-15)
                ["PX100001", "PX100001.1", "PX100001", "PX100001", "4", "A/Sydney/10/2024", "Original", "2024", "2024", "", "", "", "400", "", ""],
                # 2. Unlinked same-strain target
                ["PX100002", "PX100002.1", "PX100002", "PX100002", "4", "A/Sydney/11/2024", "", "2024-04-01", "2024", "4", "1", "", "400", "", ""],
                # 3. Truncated GenBank stub (180 bp) - GISAID full-length 400 bp must NOT be retired!
                ["PX100003", "PX100003.1", "PX100003", "PX100003", "4", "A/Sydney/12/2024", "Original", "2024-04-02", "2024", "4", "2", "", "180", "", ""],
            ],
        )
        base_fa = _write_fasta(
            tmp_path / "sequences.fa",
            {
                "PX100001": seq_ha_1,
                "PX100002": seq_ha_2,
                "PX100003": seq_stub,
            },
        )

        gisaid_tsv = _write_tsv(
            tmp_path / "gisaid_meta.tsv",
            [
                "Segment_Id", "Isolate_Id", "Isolate_Name", "Segment",
                "Passage_History", "Collection_Date", "Originating_Lab", "HA INSDC_Upload",
            ],
            [
                # Should retire (Tier 1 explicit INSDC link) & enrich PX100001 date + lab
                ["EPI_DUP_LINK", "EPI_ISL_1", "A/Sydney/10/2024", "HA", "Original", "2024-03-15", "VIDRL Melbourne", "PX100001.1"],
                # Should retire (Tier 2 unlinked same strain + 0mm with UTR overhang)
                ["EPI_DUP_STRAIN", "EPI_ISL_2", "A/Sydney/11/2024 (H3N2)", "HA", "Original", "2024-04-01", "VIDRL Melbourne", ""],
                # Uploader pasted PX100001 by mistake on an unrelated sequence -> MUST SURVIVE
                ["EPI_BAD_LINK", "EPI_ISL_3", "A/Perth/99/2024", "HA", "Original", "2024-05-01", "PathWest", "PX100001"],
                # Points to truncated stub PX100003 -> MUST SURVIVE
                ["EPI_FULL_OVER_STUB", "EPI_ISL_4", "A/Sydney/12/2024", "HA", "Original", "2024-04-02", "VIDRL", "PX100003"],
                # Egg-passaged isolate vs Original PX100001 -> MUST SURVIVE
                ["EPI_EGG_PASSAGE", "EPI_ISL_5", "A/Sydney/10/2024", "HA", "E2", "2024-03-15", "VIDRL", ""],
            ],
        )
        gisaid_fa = _write_fasta(
            tmp_path / "gisaid.fa",
            {
                "EPI_DUP_LINK": seq_ha_1,
                "EPI_DUP_STRAIN": seq_ha_2_utr,
                "EPI_BAD_LINK": seq_typo,
                "EPI_FULL_OVER_STUB": seq_ha_1,
                "EPI_EGG_PASSAGE": seq_egg_mut,
            },
        )
        mapping = tmp_path / "col_map.tsv"
        mapping.write_text(
            "Isolate_Name\tstrain\n"
            "Collection_Date\tcollection_date\n"
            "Segment\tsegment_declared\n"
            "Passage_History\tPassage_History\n"
            "Originating_Lab\tOriginating_Lab\n"
            "Isolate_Id\tIsolate_Id\n",
            encoding="utf-8",
        )

        out_tsv = tmp_path / "out" / "gB_matrix_merged.tsv"
        dedup_fa = tmp_path / "out" / "gisaid_kept.fa"

        NormalizeAndMerge(
            gb_matrix=str(base_tsv),
            input_tsv=str(gisaid_tsv),
            input_fasta=str(gisaid_fa),
            base_fasta=str(base_fa),
            dedup_fasta_out=str(dedup_fa),
            mapping_file=str(mapping),
            output=str(out_tsv),
            key="Segment_Id",
            dataset_source="gisaid",
            log_file=str(tmp_path / "merge.log"),
        ).process()

        df = pd.read_csv(out_tsv, sep="\t", dtype=str, keep_default_na=False)
        by_acc = df.set_index("primary_accession").to_dict(orient="index")

        # Retired tombstones remain in meta_data with exclusion_status='1'
        assert by_acc["EPI_DUP_LINK"]["exclusion_status"] == "1"
        assert by_acc["EPI_DUP_LINK"]["exclusion_criteria"] == "cross_db_duplicate_of:PX100001"
        assert by_acc["EPI_DUP_LINK"]["canonical_accession"] == "PX100001"

        assert by_acc["EPI_DUP_STRAIN"]["exclusion_status"] == "1"
        assert by_acc["EPI_DUP_STRAIN"]["canonical_accession"] == "PX100002"

        # Option B: Full-length GISAID (400 bp) beats truncated GenBank stub PX100003 (180 bp),
        # keeping EPI_FULL_OVER_STUB active with dual-source tag & retiring the redundant stub.
        assert by_acc["PX100003"]["exclusion_status"] == "1"
        assert by_acc["PX100003"]["canonical_accession"] == "EPI_FULL_OVER_STUB"
        assert by_acc["EPI_FULL_OVER_STUB"]["dataset_source"] == "genbank;gisaid"
        assert by_acc["EPI_FULL_OVER_STUB"]["sequence_source"] == "gisaid"

        # Protected / non-duplicate / winning rows stay active
        for kept_id in ("PX100001", "PX100002", "EPI_BAD_LINK", "EPI_FULL_OVER_STUB", "EPI_EGG_PASSAGE"):
            assert by_acc[kept_id].get("exclusion_status", "") != "1", f"{kept_id} should not be excluded"
        assert by_acc["PX100001"]["dataset_source"] == "genbank;gisaid"
        assert by_acc["PX100002"]["dataset_source"] == "genbank;gisaid"

        # Metadata coalescing + split-date lockstep on PX100001
        assert by_acc["PX100001"]["collection_date"] == "2024-03-15"
        assert by_acc["PX100001"]["collection_year"] == "2024"
        assert by_acc["PX100001"]["collection_mon"] == "3"
        assert by_acc["PX100001"]["collection_day"] == "15"
        assert by_acc["PX100001"]["Originating_Lab"] == "VIDRL Melbourne"
        assert "gisaid:EPI_DUP_LINK" in by_acc["PX100001"]["alt_accessions"]
        assert "EPI_ISL_1" in by_acc["PX100001"]["alt_accessions"]

        # Deduplicated FASTA omits only the two retired GISAID accessions
        kept_fa_text = dedup_fa.read_text(encoding="utf-8")
        assert ">EPI_DUP_LINK" not in kept_fa_text
        assert ">EPI_DUP_STRAIN" not in kept_fa_text
        assert ">EPI_BAD_LINK" in kept_fa_text
        assert ">EPI_FULL_OVER_STUB" in kept_fa_text
        assert ">EPI_EGG_PASSAGE" in kept_fa_text

    def test_merge_deduplicates_against_existing_update_db_and_pathoplexus(self, tmp_path: Path):
        """An incremental --update_db run catches external records whose INSDC
        counterpart was already committed in an earlier DB build (unsegmented virus)."""
        seq_rabv = "ACGTACGT" * 50  # 400 bp
        prev_db = tmp_path / "prev_rabv.db"
        conn = sqlite3.connect(str(prev_db))
        conn.execute(
            """
            CREATE TABLE meta_data (
                primary_accession TEXT, strain TEXT, collection_date TEXT,
                real_length TEXT, dataset_source TEXT, exclusion_status TEXT, accession_type TEXT
            )
            """
        )
        conn.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
        conn.execute(
            "INSERT INTO meta_data VALUES ('OQ555001', 'RABV/Fox/Poland/12/2023', '2023-09-01', '400', 'genbank', '0', '')"
        )
        conn.execute("INSERT INTO sequences VALUES ('OQ555001', ?)", (seq_rabv,))
        conn.commit()
        conn.close()

        # Empty new GenBank delta in this update run
        base_tsv = _write_tsv(
            tmp_path / "gB_delta.tsv",
            ["primary_accession", "accession_version", "locus", "gi_number", "strain", "real_length", "exclusion"],
            [],
        )
        ppx_tsv = _write_tsv(
            tmp_path / "ppx.tsv",
            ["accession", "isolate_name", "insdcAccessionBase", "collection_date"],
            [
                ["PP_0001.1", "RABV/Fox/Poland/12/2023", "OQ555001", "2023-09-01"],
                ["PP_0002.1", "RABV/Fox/Poland/99/2023", "", "2023-10-01"],
            ],
        )
        ppx_fa = _write_fasta(
            tmp_path / "ppx.fa",
            {"PP_0001.1": seq_rabv, "PP_0002.1": "TTTTGGGG" * 50},
        )
        mapping = tmp_path / "map.tsv"
        mapping.write_text("isolate_name\tstrain\ncollection_date\tcollection_date\n", encoding="utf-8")

        out_tsv = tmp_path / "out" / "merged.tsv"
        NormalizeAndMerge(
            gb_matrix=str(base_tsv),
            input_tsv=str(ppx_tsv),
            input_fasta=str(ppx_fa),
            update_db=str(prev_db),
            mapping_file=str(mapping),
            output=str(out_tsv),
            key="accession",
            dataset_source="pathoplexus",
            log_file=str(tmp_path / "ppx.log"),
        ).process()

        df = pd.read_csv(out_tsv, sep="\t", dtype=str, keep_default_na=False).set_index("primary_accession")
        assert df.loc["PP_0001.1", "exclusion_status"] == "1"
        assert df.loc["PP_0001.1", "canonical_accession"] == "OQ555001"
        assert df.loc["PP_0002.1"].get("exclusion_status", "") != "1"

    def test_merge_cli_matches_nextflow_merge_gisaid_contract(self, tmp_path: Path):
        """Guards the exact CLI flags MERGE_GISAID passes in vgtk-init.nf, including
        always writing --dedup_report even when 0 duplicates are detected."""
        merge_script = REPO_ROOT / "scripts" / "merge_into_gB_matrix.py"
        seq_a = "ATGCGTAC" * 50
        seq_b = "TTAACCGG" * 50
        base_tsv = _write_tsv(
            tmp_path / "gB.tsv",
            ["primary_accession", "accession_version", "locus", "gi_number", "segment", "strain", "real_length", "exclusion"],
            [["PX900001", "PX900001.1", "PX900001", "PX900001", "4", "A/Cork/1/2024", "400", ""]],
        )
        base_fa = _write_fasta(tmp_path / "sequences.fa", {"PX900001": seq_a})
        gisaid_tsv = _write_tsv(
            tmp_path / "metadata.tsv",
            ["Segment_Id", "Isolate_Name", "Segment"],
            [["EPI_UNIQ_1", "A/Dublin/99/2024", "HA"]],
        )
        gisaid_fa = _write_fasta(tmp_path / "all_nuc.fas", {"EPI_UNIQ_1": seq_b})
        mapping = tmp_path / "col_map.tsv"
        mapping.write_text("Isolate_Name\tstrain\nSegment\tsegment_declared\n", encoding="utf-8")

        out_merged = tmp_path / "gB_matrix_merged.tsv"
        out_fa = tmp_path / "gisaid_dedup.fas"
        out_rep = tmp_path / "gisaid_cross_source_dedup_report.tsv"

        proc = subprocess.run(
            [
                sys.executable,
                str(merge_script),
                "-g", str(base_tsv),
                "--base_fasta", str(base_fa),
                "-t", str(gisaid_tsv),
                "-f", str(gisaid_fa),
                "-o", str(out_merged),
                "--dedup_fasta_out", str(out_fa),
                "--dedup_report", str(out_rep),
                "-k", "Segment_Id",
                "--dataset_source", "gisaid",
                "-m", str(mapping),
                "--log_file", str(tmp_path / "merge_cli.log"),
            ],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"
        assert out_merged.is_file()
        assert out_fa.is_file()
        assert out_rep.is_file(), "MERGE_GISAID requires --dedup_report to exist even when 0 duplicates occur"
        rep_df = pd.read_csv(out_rep, sep="\t", dtype=str)
        assert "winner_accession" in rep_df.columns and len(rep_df) == 0


# ---------------------------------------------------------------------------
# 3. Iterative O(N) Newick tip pruning, renaming & cherry collapse
# ---------------------------------------------------------------------------

class TestNewickPruneAndRename:
    def test_newick_pruner_collapses_polytomy_cherry_and_preserves_branch_length(self):
        nwk = "((A:0.01,(PX1:0,EPI1:0)node_99:0.004)node_10:0.005,B:0.03);"
        pruned, stats = prune_and_rename_newick(nwk, {"EPI1"})
        assert stats == {"orig_tips": 4, "kept_tips": 3, "pruned_tips": 1, "renamed_tips": 0}
        tree = Phylo.read(StringIO(pruned), "newick")
        tips = {t.name: t.branch_length for t in tree.get_terminals()}
        assert set(tips.keys()) == {"A", "PX1", "B"}
        assert pytest.approx(tips["PX1"], abs=1e-9) == 0.004

    def test_newick_renames_retired_centroid_and_prunes_second_tip_mapping_to_same_winner(self):
        iq4 = (
            "(((NC_SEG4:0.01,REF_SEG4:0.02):0.01,"
            "(EPI_S4_CENT:0.005,EPI_S4_SING:0.001)node_c:0.01):0.01,"
            "(EPI_S4_TYPO:0.05,EPI_S4_EGG:0.02):0.01);"
        )
        pruned, stats = prune_and_rename_newick(
            iq4,
            prune_leaves=set(),
            rename_leaves={"EPI_S4_CENT": "GB_S4_WIN", "EPI_S4_SING": "GB_S4_WIN"},
        )
        assert stats == {"orig_tips": 6, "kept_tips": 5, "pruned_tips": 1, "renamed_tips": 1}
        tree = Phylo.read(StringIO(pruned), "newick")
        tip_names = [t.name for t in tree.get_terminals()]
        assert len(tip_names) == len(set(tip_names)), f"Duplicate tips after rename: {tip_names}"
        assert set(tip_names) == {"NC_SEG4", "REF_SEG4", "GB_S4_WIN", "EPI_S4_TYPO", "EPI_S4_EGG"}
        tips = {t.name: t.branch_length for t in tree.get_terminals()}
        # node_c collapsed when EPI_S4_SING was pruned, adding 0.01 + 0.005 = 0.015 onto GB_S4_WIN
        assert pytest.approx(tips["GB_S4_WIN"], abs=1e-9) == 0.015


class TestDedupStrainAlignmentCliAndTiers:
    def test_cli_preserves_protected_refs_and_handles_tier1_and_intra_pools(self, tmp_path: Path):
        seq_ref = "ATG" * 80
        seq_gb = "ATG" * 79 + "ATC"
        seq_gb_1mm = "ATG" * 78 + "ATAATC"
        seq_gis = "ATG" * 78 + "ACCATC"

        ref_list = tmp_path / "ref_list.txt"
        ref_list.write_text("NC_REF4\tmaster\t4\nREF_DUP_SAME_STRAIN\treference\t4\n", encoding="utf-8")

        meta_in = tmp_path / "meta.tsv"
        pd.DataFrame([
            # Protected master & reference sharing the same strain: neither may ever be retired!
            {"primary_accession": "NC_REF4", "accession_type": "master", "dataset_source": "genbank", "segment": "4", "Parsed_strain": "A/Ref/1/2020", "Passage_History": "Original", "collection_date": "2020-01-01", "gisaid_ha_insdc_upload": ""},
            {"primary_accession": "REF_DUP_SAME_STRAIN", "accession_type": "query", "dataset_source": "genbank", "segment": "4", "Parsed_strain": "A/Ref/1/2020", "Passage_History": "Original", "collection_date": "2020-01-01", "gisaid_ha_insdc_upload": ""},
            # Intra-GenBank duplicate pair (PX200002 has full date)
            {"primary_accession": "PX200001", "accession_type": "query", "dataset_source": "genbank", "segment": "4", "Parsed_strain": "A/Madrid/10/2024", "Passage_History": "Original", "collection_date": "2024-03", "gisaid_ha_insdc_upload": ""},
            {"primary_accession": "PX200002", "accession_type": "query", "dataset_source": "genbank", "segment": "4", "Parsed_strain": "A/Madrid/10/2024", "Passage_History": "Original", "collection_date": "2024-03-11", "gisaid_ha_insdc_upload": ""},
            # Explicit Tier 1 INSDC link (1 mismatch <= 2mm, redirected from PX200001 -> PX200002!)
            {"primary_accession": "EPI_LINKED", "accession_type": "query", "dataset_source": "gisaid", "segment": "4", "Parsed_strain": "A/Madrid/10_variant_label/2024", "Passage_History": "Original", "collection_date": "2024-03-11", "gisaid_ha_insdc_upload": "PX200001.1"},
            # Intra-GISAID duplicate pair
            {"primary_accession": "EPI_ONLY_1", "accession_type": "query", "dataset_source": "gisaid", "segment": "4", "Parsed_strain": "A/Lima/5/2024", "Passage_History": "MDCK", "collection_date": "2024-04-02", "gisaid_ha_insdc_upload": ""},
            {"primary_accession": "EPI_ONLY_2", "accession_type": "query", "dataset_source": "gisaid", "segment": "4", "Parsed_strain": "A/Lima/5/2024", "Passage_History": "C1", "collection_date": "2024-04", "gisaid_ha_insdc_upload": ""},
        ]).to_csv(meta_in, sep="\t", index=False)

        msa_in = tmp_path / "refset_4_merged_MSA.fasta"
        _write_fasta(
            msa_in,
            {
                "NC_REF4": seq_ref,
                "REF_DUP_SAME_STRAIN": seq_ref,
                "PX200001": seq_gb,
                "PX200002": seq_gb,
                "EPI_LINKED": seq_gb_1mm,
                "EPI_ONLY_1": seq_gis,
                "EPI_ONLY_2": seq_gis,
            },
        )

        out_meta = tmp_path / "out" / "dedup_meta.tsv"
        out_msa_dir = tmp_path / "out" / "clean_msa"
        out_ret_upd = tmp_path / "out" / "retired_upd.txt"
        report_tsv = tmp_path / "out" / "report.tsv"

        proc = subprocess.run(
            [
                sys.executable,
                str(DEDUP_STRAIN_SCRIPT),
                "--meta_data", str(meta_in),
                "--msa", str(msa_in),
                "--ref_list", str(ref_list),
                "--out_meta", str(out_meta),
                "--out_msa_dir", str(out_msa_dir),
                "--out_retired_update_ids", str(out_ret_upd),
                "--report_tsv", str(report_tsv),
            ],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"

        df = pd.read_csv(out_meta, sep="\t", dtype=str, keep_default_na=False).set_index("primary_accession")
        # Protected ref_list entries survive even though REF_DUP_SAME_STRAIN was marked 'query' in meta_data
        assert df.loc["NC_REF4", "exclusion_status"] in ("", "0")
        assert df.loc["REF_DUP_SAME_STRAIN", "exclusion_status"] in ("", "0")
        assert df.loc["PX200002", "exclusion_status"] in ("", "0")
        assert df.loc["PX200001", "exclusion_status"] == "1"
        assert df.loc["EPI_LINKED", "exclusion_status"] == "1"
        assert df.loc["EPI_LINKED", "canonical_accession"] == "PX200002"
        assert df.loc["EPI_ONLY_1", "exclusion_status"] in ("", "0")
        assert df.loc["EPI_ONLY_2", "exclusion_status"] == "1"

        clean_fa = (out_msa_dir / msa_in.name).read_text(encoding="utf-8")
        for kept in ("NC_REF4", "REF_DUP_SAME_STRAIN", "PX200002", "EPI_ONLY_1"):
            assert f">{kept}" in clean_fa
        for dropped in ("PX200001", "EPI_LINKED", "EPI_ONLY_2"):
            assert dropped not in clean_fa
        assert report_tsv.is_file()


class TestValidateDbTreeCheckStrainDuplicatesCli:
    def test_cli_check_strain_duplicates_fails_on_single_source_unretired_duplicate(self, tmp_path: Path):
        seq_m = "ATG" * 40
        seq_q = "ATG" * 39 + "ATC"
        db_path = tmp_path / "single_src_dup.db"
        with sqlite3.connect(str(db_path)) as conn:
            conn.execute(
                "CREATE TABLE meta_data (primary_accession TEXT, dataset_source TEXT, accession_type TEXT, "
                "exclusion_status TEXT, strain TEXT, Passage_History TEXT, collection_date TEXT, real_length TEXT)"
            )
            conn.executemany(
                "INSERT INTO meta_data VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    ("NC_M", "genbank", "master", "0", "RABV/Master/1990", "Original", "1990-01-01", "120"),
                    ("OQ_DUP1", "genbank", "query", "0", "RABV/Dog/Kenya/12/2024", "Original", "2024-05-01", "120"),
                    ("OQ_DUP2", "genbank", "query", "0", "RABV/Dog/Kenya/12/2024", "Original", "2024-05-01", "120"),
                ],
            )
            conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
            conn.executemany(
                "INSERT INTO sequence_alignment VALUES (?, 'NC_M', ?)",
                [("NC_M", seq_m), ("OQ_DUP1", seq_q), ("OQ_DUP2", seq_q)],
            )
            conn.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
            conn.executemany(
                "INSERT INTO sequences VALUES (?, ?)",
                [("NC_M", seq_m), ("OQ_DUP1", seq_q), ("OQ_DUP2", seq_q)],
            )
            conn.execute(
                "CREATE TABLE features (accession TEXT, master_ref_accession TEXT, reference_accession TEXT, "
                "aln_start TEXT, aln_end TEXT, cds_start TEXT, cds_end TEXT, cds_start_OG_seq TEXT, cds_end_OG_seq TEXT, product TEXT)"
            )
            conn.executemany(
                "INSERT INTO features VALUES (?, 'NC_M', 'NC_M', '1', '120', '1', '120', '1', '120', 'N')",
                [("NC_M",), ("OQ_DUP1",), ("OQ_DUP2",)],
            )
            conn.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT)")
            conn.execute("INSERT INTO trees VALUES ('usher', 'usher', '(NC_M:0.01,(OQ_DUP1:0,OQ_DUP2:0):0.002);')")
            conn.commit()

        proc = subprocess.run(
            [
                sys.executable,
                str(VALIDATE_DB_TREE),
                "--db", str(db_path),
                "--outdir", str(tmp_path / "val_dup"),
                "--check-strain-duplicates",
            ],
            capture_output=True,
            text=True,
        )
        assert proc.returncode != 0
        assert "same-strain duplicate sequences detected" in (proc.stdout + proc.stderr)


# ---------------------------------------------------------------------------
# 4. Update-mode both arrival orders, UsherPlacement, & CreateSqliteDB
# ---------------------------------------------------------------------------


class TestUpdateModeBothArrivalOrdersAndAlignmentGate:
    def test_merge_into_gb_matrix_handles_intra_genbank_and_intra_gisaid(self, tmp_path: Path):
        from merge_into_gB_matrix import NormalizeAndMerge

        seq = "ATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCC" * 5
        gb_tsv = tmp_path / "gb.tsv"
        pd.DataFrame([
            {"primary_accession": "OQ100001", "segment": "4", "strain": "A/Texas/50/2024", "Passage_History": "Original", "collection_date": "2024-02"},
            {"primary_accession": "OQ100002", "segment": "4", "strain": "A/Texas/50/2024(H3N2)", "Passage_History": "Clinical", "collection_date": "2024-02-14"},
        ]).to_csv(gb_tsv, sep="\t", index=False)
        gb_fa = tmp_path / "gb.fa"
        gb_fa.write_text(f">OQ100001\n{seq}\n>OQ100002\n{seq}\n", encoding="utf-8")

        gis_tsv = tmp_path / "gis.tsv"
        pd.DataFrame([
            {"Segment_Id": "EPI900001", "Segment": "4", "Isolate_Name": "A/Sydney/77/2024", "Passage_History": "MDCK1", "Collection_Date": "2024-03"},
            {"Segment_Id": "EPI900002", "Segment": "4", "Isolate_Name": "A/Sydney/77/2024", "Passage_History": "C1", "Collection_Date": "2024-03-19"},
        ]).to_csv(gis_tsv, sep="\t", index=False)
        gis_fa = tmp_path / "gis.fa"
        gis_fa.write_text(f">EPI900001\n{seq}\n>EPI900002\n{seq}\n", encoding="utf-8")

        col_map = tmp_path / "col_map.tsv"
        col_map.write_text("Segment_Id\tprimary_accession\nSegment\tsegment\nIsolate_Name\tstrain\nCollection_Date\tcollection_date\nPassage_History\tPassage_History\n", encoding="utf-8")

        out_tsv = tmp_path / "merged.tsv"
        dedup_fa = tmp_path / "gis_dedup.fa"
        nm = NormalizeAndMerge(
            gb_matrix=str(gb_tsv),
            input_tsv=str(gis_tsv),
            input_fasta=str(gis_fa),
            mapping_file=str(col_map),
            output=str(out_tsv),
            key="Segment_Id",
            dataset_source="gisaid",
            base_fasta=str(gb_fa),
            dedup_fasta_out=str(dedup_fa),
            log_file=str(tmp_path / "merge.log"),
        )
        nm.process()

        merged = pd.read_csv(out_tsv, sep="\t", dtype=str, keep_default_na=False)
        # OQ100002 has full date 2024-02-14 vs OQ100001 2024-02 -> OQ100002 wins intra-GenBank!
        by_acc = merged.set_index("primary_accession")
        assert by_acc.loc["OQ100002", "exclusion_status"] in ("", "0")
        assert by_acc.loc["OQ100001", "exclusion_status"] == "1"
        assert "intra_db_duplicate_of:OQ100002" in by_acc.loc["OQ100001", "exclusion_criteria"]

        # EPI900002 has full date 2024-03-19 vs EPI900001 2024-03 -> EPI900002 wins intra-GISAID,
        # EPI900001 is kept as an exclusion_status='1' tombstone in merged.tsv and omitted from dedup_fa!
        assert by_acc.loc["EPI900002", "exclusion_status"] in ("", "0")
        assert by_acc.loc["EPI900001", "exclusion_status"] == "1"
        assert "intra_db_duplicate_of:EPI900002" in by_acc.loc["EPI900001", "exclusion_criteria"]
        fa_text = dedup_fa.read_text(encoding="utf-8")
        assert ">EPI900002" in fa_text
        assert "EPI900001" not in fa_text

    def test_dedup_strain_alignment_and_usher_reverse_arrival_order(self, tmp_path: Path):
        from DedupStrainAlignment import deduplicate_strain_alignments
        from UsherPlacement import UsherPlacement

        seq_ref = "ATG" * 80
        seq_dup = "ATG" * 79 + "ATC"
        seq_other = "ATG" * 78 + "ACCATC"

        # Seed update_db where GISAID arrived FIRST (EPI_FIRST is active on the starter tree)
        upd_db = tmp_path / "seed_update.db"
        with sqlite3.connect(str(upd_db)) as conn:
            conn.execute(
                "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, dataset_source TEXT, "
                "segment TEXT, Parsed_strain TEXT, Passage_History TEXT, collection_date TEXT, "
                "Lineage TEXT, exclusion_status TEXT, real_length TEXT)"
            )
            conn.executemany(
                "INSERT INTO meta_data VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    ("REF_M4", "master", "genbank", "4", "A/Ref/1/2020", "", "2020-01-01", "", "0", "240"),
                    ("EPI_FIRST", "query", "gisaid", "4", "A/Berlin/99/2025", "Original", "2025-01-20", "J.2", "0", "240"),
                    ("GB_KEEP", "query", "genbank", "4", "A/Oslo/11/2025", "Original", "2025-01-10", "", "0", "240"),
                ],
            )
            conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, segment TEXT, alignment TEXT)")
            conn.executemany(
                "INSERT INTO sequence_alignment VALUES (?, ?, ?)",
                [
                    ("REF_M4", "4", seq_ref),
                    ("EPI_FIRST", "4", seq_dup),
                    ("GB_KEEP", "4", seq_other),
                ],
            )
            conn.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, segment_key TEXT, newick TEXT)")
            conn.execute(
                "INSERT INTO trees VALUES ('usher_seg4', 'usher', '4', 'seg_4', '(REF_M4:0.001,(EPI_FIRST:0.0,GB_KEEP:0.003):0.002);')"
            )
            conn.commit()

        # Later update batch brings in GenBank upload PQ_LATER for A/Berlin/99/2025 (with coarser date 2025-01)
        # PLUS an incoming duplicate EPI_DUP_OF_KEEP for A/Oslo/11/2025
        batch_meta = pd.DataFrame([
            {"primary_accession": "REF_M4", "accession_type": "master", "dataset_source": "genbank", "segment": "4", "Parsed_strain": "A/Ref/1/2020", "Passage_History": "", "collection_date": "2020-01-01", "exclusion_status": "0"},
            {"primary_accession": "PQ_LATER", "accession_type": "query", "dataset_source": "genbank", "segment": "4", "Parsed_strain": "A/Berlin/99/2025", "Passage_History": "Clinical", "collection_date": "2025-01", "exclusion_status": "0"},
            {"primary_accession": "EPI_DUP_OF_KEEP", "accession_type": "query", "dataset_source": "gisaid", "segment": "4", "Parsed_strain": "A/Oslo/11/2025", "Passage_History": "Original", "collection_date": "2025-01-10", "exclusion_status": "0"},
        ])
        msa_in = tmp_path / "segment_4_merged_MSA.fasta"
        msa_in.write_text(
            f">REF_M4\n{seq_ref}\n>PQ_LATER\n{seq_dup}\n>EPI_DUP_OF_KEEP\n{seq_other}\n",
            encoding="utf-8",
        )
        clean_dir = tmp_path / "clean_msa"
        dedup_df, retired_upd_ids, pairs = deduplicate_strain_alignments(
            df_meta=batch_meta,
            msa_files=[msa_in],
            out_msa_dir=clean_dir,
            update_db_path=str(upd_db),
        )

        # Reverse arrival order: EPI_FIRST in update_db is superseded by PQ_LATER!
        assert retired_upd_ids == {"EPI_FIRST"}
        by_acc = dedup_df.set_index("primary_accession")
        assert by_acc.loc["PQ_LATER", "exclusion_status"] == "0"
        # PQ_LATER inherited full collection_date and Lineage from superseded EPI_FIRST!
        assert by_acc.loc["PQ_LATER", "collection_date"] == "2025-01-20"
        assert by_acc.loc["PQ_LATER", "Lineage"] == "J.2"
        # Forward arrival order: EPI_DUP_OF_KEEP is retired in batch and omitted from clean MSA!
        assert by_acc.loc["EPI_DUP_OF_KEEP", "exclusion_status"] == "1"
        clean_msa_text = (clean_dir / msa_in.name).read_text(encoding="utf-8")
        assert ">PQ_LATER" in clean_msa_text
        assert "EPI_DUP_OF_KEEP" not in clean_msa_text

        # UsherPlacement prunes EPI_FIRST from the starter tree via --retire_ids_file
        retire_file = tmp_path / "retired_upd.txt"
        retire_file.write_text("EPI_FIRST\n", encoding="utf-8")
        up = UsherPlacement(
            padded_aln=str(clean_dir / msa_in.name),
            output_dir=str(tmp_path / "usher_out"),
            update_db=str(upd_db),
            retire_ids_file=str(retire_file),
        )
        tree_out, ids_out = up.prepare_update_assets()
        assert "EPI_FIRST" not in Path(tree_out).read_text(encoding="utf-8")
        assert "EPI_FIRST" not in Path(ids_out).read_text(encoding="utf-8")
        assert "GB_KEEP" in Path(tree_out).read_text(encoding="utf-8")

    def test_validate_db_tree_flags_unretired_strain_duplicates(self, tmp_path: Path):
        import ValidateDbTree as vdt

        seq = "ATG" * 90
        bad_db = tmp_path / "unretired.db"
        with sqlite3.connect(str(bad_db)) as conn:
            conn.execute(
                "CREATE TABLE meta_data (primary_accession TEXT, dataset_source TEXT, accession_type TEXT, "
                "exclusion_status TEXT, segment TEXT, strain TEXT, Passage_History TEXT, real_length TEXT)"
            )
            conn.executemany(
                "INSERT INTO meta_data VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    ("OQ_1", "genbank", "query", "0", "4", "A/Paris/10/2024", "Original", "270"),
                    ("EPI_1", "gisaid", "query", "0", "4", "A/Paris/10/2024", "Original", "270"),
                    ("OQ_2A", "genbank", "query", "0", "4", "A/Rome/20/2024", "MDCK", "270"),
                    ("OQ_2B", "genbank", "query", "0", "4", "A/Rome/20/2024", "C1", "270"),
                ],
            )
            conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, segment TEXT, alignment TEXT)")
            conn.executemany(
                "INSERT INTO sequence_alignment VALUES (?, ?, ?)",
                [("OQ_1", "4", seq), ("EPI_1", "4", seq), ("OQ_2A", "4", seq), ("OQ_2B", "4", seq)],
            )
            conn.commit()

        with sqlite3.connect(f"file:{bad_db}?mode=ro", uri=True) as ro_conn:
            meta_cols = vdt.get_table_columns(ro_conn, "meta_data")
            res = vdt.validate_strain_sequence_duplicates(ro_conn, meta_cols)
            assert res["ok"] is False
            codes = {f["code"]: f["count"] for f in res["findings"]}
            assert codes["cross_db_duplicate_strain_sequences"] == 1
            assert codes["intra_db_duplicate_strain_sequences"] == 1

    def test_create_sqlite_db_update_reverse_arrival_retires_and_prunes_gisaid(self, tmp_path: Path):
        from CreateSqliteDB import CreateSqliteDB

        seq_m = "ATG" * 60
        seq_q = "ATG" * 59 + "ATC"

        def _make_bundle(sub: Path, meta_rows, aln_rows, feat_rows, fa_records, nwk_text):
            sub.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(meta_rows).to_csv(sub / "meta.tsv", sep="\t", index=False)
            pd.DataFrame(aln_rows).to_csv(sub / "aln.tsv", sep="\t", index=False)
            pd.DataFrame(feat_rows).to_csv(sub / "feat.tsv", sep="\t", index=False)
            pd.DataFrame([{"name": "HA", "description": "Hemagglutinin"}]).to_csv(sub / "gene.tsv", sep="\t", index=False)
            pd.DataFrame([{"m49_code": "001", "name": "World"}]).to_csv(sub / "c.csv", index=False)
            pd.DataFrame([{"code": "X", "name": "Inter"}]).to_csv(sub / "i.csv", index=False)
            pd.DataFrame([{"code": "Y", "name": "Region"}]).to_csv(sub / "r.csv", index=False)
            pd.DataFrame([{"code": "Z", "name": "SubRegion"}]).to_csv(sub / "sr.csv", index=False)
            pd.DataFrame([{"Software": "Python", "Version": "3.12"}]).to_csv(sub / "sw.tsv", sep="\t", index=False)
            pd.DataFrame(columns=["primary_accession", "reference", "insertion", "segment"]).to_csv(sub / "ins.tsv", sep="\t", index=False)
            pd.DataFrame([{"taxonomy_id": "9606", "scientific_name": "Homo sapiens"}]).to_csv(sub / "ht.tsv", sep="\t", index=False)
            (sub / "seqs.fa").write_text(
                "".join(f">{k}\n{v}\n" for k, v in fa_records),
                encoding="utf-8",
            )
            (sub / "usher.nwk").write_text(nwk_text, encoding="utf-8")
            (sub / "manifest.tsv").write_text(
                f"source\tname\tsegment_key\tpath\nusher\tusher_seg4\tseg_4\t{sub / 'usher.nwk'}\n",
                encoding="utf-8",
            )
            return sub

        # Batch 1 (create): GISAID upload EPI_EARLY arrives first
        b1 = _make_bundle(
            tmp_path / "b1",
            meta_rows=[
                {"primary_accession": "REF4", "accession_type": "master", "dataset_source": "genbank", "segment": "4", "strain": "A/Ref/1/2020", "Passage_History": "", "collection_date": "2020-01-01", "exclusion_status": "0", "host_taxa_id": "9606"},
                {"primary_accession": "EPI_EARLY", "accession_type": "query", "dataset_source": "gisaid", "segment": "4", "strain": "A/Zurich/42/2025", "Passage_History": "Original", "collection_date": "2025-02-18", "exclusion_status": "0", "host_taxa_id": "9606"},
            ],
            aln_rows=[
                {"primary_accession": "REF4", "alignment_name": "REF4", "alignment": seq_m, "segment": "4"},
                {"primary_accession": "EPI_EARLY", "alignment_name": "REF4", "alignment": seq_q, "segment": "4"},
            ],
            feat_rows=[
                {"accession": "REF4", "master_ref_accession": "REF4", "reference_accession": "REF4", "aln_start": "1", "aln_end": "180", "cds_start": "1", "cds_end": "180", "cds_start_OG_seq": "1", "cds_end_OG_seq": "180", "product": "HA", "segment": "4"},
                {"accession": "EPI_EARLY", "master_ref_accession": "REF4", "reference_accession": "REF4", "aln_start": "1", "aln_end": "180", "cds_start": "1", "cds_end": "180", "cds_start_OG_seq": "1", "cds_end_OG_seq": "180", "product": "HA", "segment": "4"},
            ],
            fa_records=[("REF4", seq_m), ("EPI_EARLY", seq_q)],
            nwk_text="(REF4:0.001,EPI_EARLY:0.002);",
        )
        CreateSqliteDB(
            meta_data=str(b1 / "meta.tsv"), features=str(b1 / "feat.tsv"), pad_aln=str(b1 / "aln.tsv"),
            gene_info=str(b1 / "gene.tsv"), m49_countries=str(b1 / "c.csv"), m49_interm_region=str(b1 / "i.csv"),
            m49_regions=str(b1 / "r.csv"), m49_sub_regions=str(b1 / "sr.csv"), proj_settings=str(b1 / "sw.tsv"),
            fasta_sequence_file=str(b1 / "seqs.fa"), insertions=str(b1 / "ins.tsv"), host_taxa_file=str(b1 / "ht.tsv"),
            base_dir=str(tmp_path), output_dir="db1", db_name="flu", db_status="new db",
            tree_manifest=str(b1 / "manifest.tsv"), is_segmented=True,
        ).create_db()
        seed_db = tmp_path / "db1" / "flu.db"

        # Batch 2 (--update): GenBank upload PQ_ARRIVED_LATER arrives for the same strain!
        b2 = _make_bundle(
            tmp_path / "b2",
            meta_rows=[
                {"primary_accession": "REF4", "accession_type": "master", "dataset_source": "genbank", "segment": "4", "strain": "A/Ref/1/2020", "Passage_History": "", "collection_date": "2020-01-01", "exclusion_status": "0", "host_taxa_id": "9606"},
                {"primary_accession": "PQ_ARRIVED_LATER", "accession_type": "query", "dataset_source": "genbank", "segment": "4", "strain": "A/Zurich/42/2025", "Passage_History": "Original", "collection_date": "2025-02", "exclusion_status": "0", "host_taxa_id": "9606"},
            ],
            aln_rows=[
                {"primary_accession": "REF4", "alignment_name": "REF4", "alignment": seq_m, "segment": "4"},
                {"primary_accession": "PQ_ARRIVED_LATER", "alignment_name": "REF4", "alignment": seq_q, "segment": "4"},
            ],
            feat_rows=[
                {"accession": "REF4", "master_ref_accession": "REF4", "reference_accession": "REF4", "aln_start": "1", "aln_end": "180", "cds_start": "1", "cds_end": "180", "cds_start_OG_seq": "1", "cds_end_OG_seq": "180", "product": "HA", "segment": "4"},
                {"accession": "PQ_ARRIVED_LATER", "master_ref_accession": "REF4", "reference_accession": "REF4", "aln_start": "1", "aln_end": "180", "cds_start": "1", "cds_end": "180", "cds_start_OG_seq": "1", "cds_end_OG_seq": "180", "product": "HA", "segment": "4"},
            ],
            fa_records=[("REF4", seq_m), ("PQ_ARRIVED_LATER", seq_q)],
            # Suppose even an unpruned tree arrived with both EPI_EARLY and PQ_ARRIVED_LATER: CreateSqliteDB prunes EPI_EARLY!
            nwk_text="(REF4:0.001,(EPI_EARLY:0.0,PQ_ARRIVED_LATER:0.0):0.002);",
        )
        CreateSqliteDB(
            meta_data=str(b2 / "meta.tsv"), features=str(b2 / "feat.tsv"), pad_aln=str(b2 / "aln.tsv"),
            gene_info=str(b2 / "gene.tsv"), m49_countries=str(b2 / "c.csv"), m49_interm_region=str(b2 / "i.csv"),
            m49_regions=str(b2 / "r.csv"), m49_sub_regions=str(b2 / "sr.csv"), proj_settings=str(b2 / "sw.tsv"),
            fasta_sequence_file=str(b2 / "seqs.fa"), insertions=str(b2 / "ins.tsv"), host_taxa_file=str(b2 / "ht.tsv"),
            base_dir=str(tmp_path), output_dir="db2", db_name="flu", db_status="last updated",
            tree_manifest=str(b2 / "manifest.tsv"), update=True, update_db=str(seed_db), is_segmented=True,
        ).create_db()
        updated_db = tmp_path / "db2" / "flu.db"

        with sqlite3.connect(str(updated_db)) as conn:
            meta = pd.read_sql_query("SELECT * FROM meta_data", conn).set_index("primary_accession")
            assert meta.loc["EPI_EARLY", "exclusion_status"] == "1"
            assert meta.loc["EPI_EARLY", "canonical_accession"] == "PQ_ARRIVED_LATER"
            assert meta.loc["PQ_ARRIVED_LATER", "exclusion_status"] == "0"
            # Full day-precision collection_date coalesced from superseded EPI_EARLY!
            assert meta.loc["PQ_ARRIVED_LATER", "collection_date"] == "2025-02-18"
            aln_accs = {r[0] for r in conn.execute("SELECT primary_accession FROM sequence_alignment")}
            assert aln_accs == {"REF4", "PQ_ARRIVED_LATER"}
            nwk = conn.execute("SELECT newick FROM trees WHERE source='usher'").fetchone()[0]
            assert "EPI_EARLY" not in nwk
            assert "PQ_ARRIVED_LATER" in nwk
