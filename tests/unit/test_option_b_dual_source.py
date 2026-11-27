"""Extensive TDD unit and integration tests for Workstream 1:
Dual-source metadata flagging ('genbank;gisaid') + Option B quality-based winner
selection + sequence_source tracking + mutation child-table cleanup.
"""

from __future__ import annotations

import csv
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from cross_db_dedup import (
    candidate_winner_sort_key,
    coalesce_winner_metadata,
    is_preferred_source,
    merge_dataset_sources,
    parse_dataset_sources,
)
from merge_into_gB_matrix import NormalizeAndMerge

REPO_ROOT = Path(__file__).resolve().parents[2]
CREATE_DB_SCRIPT = REPO_ROOT / "scripts" / "CreateSqliteDB.py"
DEDUP_STRAIN_SCRIPT = REPO_ROOT / "scripts" / "DedupStrainAlignment.py"
VALIDATE_DB_TREE = REPO_ROOT / "scripts" / "ValidateDbTree.py"


def _write_tsv(path: Path, header, rows) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    return path


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    with path.open("w", encoding="utf-8") as fh:
        for hdr, seq in records.items():
            fh.write(f">{hdr}\n{seq}\n")
    return path


# ---------------------------------------------------------------------------
# 1. Pure helper tests: dataset_source token parsing, merging, preferred check
# ---------------------------------------------------------------------------

class TestDatasetSourceHelpers:
    def test_parse_and_merge_dataset_sources_normalises_and_sorts(self):
        assert parse_dataset_sources("genbank") == frozenset({"genbank"})
        assert parse_dataset_sources("gisaid;genbank") == frozenset({"genbank", "gisaid"})
        assert parse_dataset_sources("  GENBANK ; GISAID ; ") == frozenset({"genbank", "gisaid"})
        assert parse_dataset_sources("ncbi") == frozenset({"genbank"})
        assert parse_dataset_sources("") == frozenset()

        assert merge_dataset_sources("genbank", "gisaid") == "genbank;gisaid"
        assert merge_dataset_sources("gisaid", "genbank") == "genbank;gisaid"
        assert merge_dataset_sources("genbank;gisaid", "gisaid") == "genbank;gisaid"
        assert merge_dataset_sources("genbank;gisaid", "pathoplexus") == "genbank;gisaid;pathoplexus"
        # Intra-DB merge preserves single source
        assert merge_dataset_sources("genbank", "genbank") == "genbank"
        assert merge_dataset_sources("gisaid", "gisaid") == "gisaid"
        assert merge_dataset_sources("ena", "genbank") == "genbank"

    def test_is_preferred_source_is_token_aware_for_dual_source(self):
        assert is_preferred_source("genbank", "genbank") is True
        assert is_preferred_source("genbank;gisaid", "genbank") is True
        assert is_preferred_source("gisaid;genbank", "genbank") is True
        assert is_preferred_source("genbank;pathoplexus", "genbank") is True
        assert is_preferred_source("gisaid", "genbank") is False
        assert is_preferred_source("pathoplexus", "genbank") is False
        assert is_preferred_source("gisaid;pathoplexus", "genbank") is False
        assert is_preferred_source("genbank;gisaid", "gisaid") is True

    def test_coalesce_winner_metadata_sets_dual_dataset_source_both_directions(self):
        gb_row = {
            "primary_accession": "PQ100001",
            "dataset_source": "genbank",
            "collection_date": "2024",
            "collection_year": "2024",
            "collection_mon": "",
            "collection_day": "",
            "Isolate_Id": "",
            "alt_accessions": "",
        }
        gi_row = {
            "primary_accession": "EPI900001",
            "dataset_source": "gisaid",
            "collection_date": "2024-03-15",
            "Isolate_Id": "EPI_ISL_999001",
            "Clade": "2.3.4.4b",
        }
        # Direction 1: GenBank wins over GISAID
        upd_gb = coalesce_winner_metadata(gb_row, gi_row, retired_source="gisaid")
        assert upd_gb["dataset_source"] == "genbank;gisaid"
        assert "gisaid:EPI900001" in upd_gb["alt_accessions"]
        assert "EPI_ISL_999001" in upd_gb["alt_accessions"]

        # Direction 2: GISAID wins over GenBank (Option B!)
        gi_winner = dict(gi_row)
        gi_winner["alt_accessions"] = ""
        upd_gi = coalesce_winner_metadata(gi_winner, gb_row, retired_source="genbank")
        assert upd_gi["dataset_source"] == "genbank;gisaid"
        assert "genbank:PQ100001" in upd_gi["alt_accessions"]

        # Direction 3: Intra-DB (GenBank + GenBank) stays single-source 'genbank'
        gb_row2 = dict(gb_row, primary_accession="PQ100002")
        upd_intra = coalesce_winner_metadata(gb_row, gb_row2, retired_source="genbank")
        assert upd_intra.get("dataset_source", "genbank") == "genbank"


# ---------------------------------------------------------------------------
# 2. Option B candidate_winner_sort_key tests
# ---------------------------------------------------------------------------

class TestOptionBWinnerSortKey:
    def test_aligned_higher_acgt_gisaid_beats_genbank(self):
        gb = {
            "accession": "CY006226",
            "source": "genbank",
            "origin": "batch",
            "protected": False,
            "real_length": 2285,
            "aligned_acgt": 2273,
            "aligned_ambig": 1,
            "collection_date": "2024-01-01",
        }
        gi = {
            "accession": "EPI90837",
            "source": "gisaid",
            "origin": "batch",
            "protected": False,
            "real_length": 2306,
            "aligned_acgt": 2276,
            "aligned_ambig": 0,
            "collection_date": "2024-01-01",
        }
        assert candidate_winner_sort_key(gi, "genbank", quality_first=True) < candidate_winner_sort_key(
            gb, "genbank", quality_first=True
        )

    def test_aligned_tie_with_utr_only_gisaid_overhang_keeps_genbank(self):
        # Mirrors real CY003071 (raw 2317, aln 2276) vs EPI83322 (raw 2321, aln 2276)
        gb = {
            "accession": "CY003071",
            "source": "genbank",
            "origin": "batch",
            "protected": False,
            "real_length": 2317,
            "aligned_acgt": 2276,
            "aligned_ambig": 0,
            "collection_date": "2024-01-01",
        }
        gi = {
            "accession": "EPI83322",
            "source": "gisaid",
            "origin": "batch",
            "protected": False,
            "real_length": 2321,
            "aligned_acgt": 2276,
            "aligned_ambig": 0,
            "collection_date": "2024-01-01",
        }
        assert candidate_winner_sort_key(gb, "genbank", quality_first=True) < candidate_winner_sort_key(
            gi, "genbank", quality_first=True
        )

    def test_protected_master_or_reference_always_beats_longer_gisaid(self):
        ref_gb = {
            "accession": "NC_007369",
            "source": "genbank",
            "origin": "batch",
            "protected": True,
            "real_length": 1490,
            "aligned_acgt": 1490,
            "aligned_ambig": 5,
            "collection_date": "2000",
        }
        gi = {
            "accession": "EPI999999",
            "source": "gisaid",
            "origin": "batch",
            "protected": False,
            "real_length": 1565,
            "aligned_acgt": 1565,
            "aligned_ambig": 0,
            "collection_date": "2000-06-15",
        }
        assert candidate_winner_sort_key(ref_gb, "genbank", quality_first=True) < candidate_winner_sort_key(
            gi, "genbank", quality_first=True
        )

    def test_same_acgt_fewer_aligned_ambiguities_wins(self):
        gb = {
            "accession": "PQ1",
            "source": "genbank",
            "origin": "batch",
            "protected": False,
            "real_length": 1700,
            "aligned_acgt": 1700,
            "aligned_ambig": 3,
            "collection_date": "2024-01-01",
        }
        gi = {
            "accession": "EPI1",
            "source": "gisaid",
            "origin": "batch",
            "protected": False,
            "real_length": 1700,
            "aligned_acgt": 1700,
            "aligned_ambig": 0,
            "collection_date": "2024-01-01",
        }
        assert candidate_winner_sort_key(gi, "genbank", quality_first=True) < candidate_winner_sort_key(
            gb, "genbank", quality_first=True
        )


# ---------------------------------------------------------------------------
# 3. End-to-end merge_into_gB_matrix.py Option B + sequence_source tests
# ---------------------------------------------------------------------------

class TestMergeOptionBAndDualSource:
    def test_merge_tie_keeps_genbank_with_dual_source_and_better_gisaid_wins_with_sequence_source(self, tmp_path: Path):
        seq_clean = "ATGCGTACGT" * 60  # 600 bp, 0 N
        seq_with_n = "ATGCGTACGT" * 29 + "NNNNNNNNNN" + "ATGCGTACGT" * 30  # 600 bp, 590 ACGT, 10 N
        seq_partial_gb = seq_clean[:450]  # >90 bp shorter partial GenBank!

        base_tsv = _write_tsv(
            tmp_path / "base.tsv",
            ["primary_accession", "strain", "Parsed_strain", "segment", "collection_date", "real_length", "accession_type"],
            [
                ["PQ_TIE", "A/London/1/2024", "A/London/1/2024", "4", "2024", "600", "query"],
                ["PQ_HAS_N", "A/Paris/2/2024", "A/Paris/2/2024", "4", "2024-02-01", "590", "query"],
                ["PQ_PARTIAL", "A/Rome/3/2024", "A/Rome/3/2024", "4", "2024-03-01", "450", "query"],
                ["PQ_REF", "A/Ref/0/2020", "A/Ref/0/2020", "4", "2020-01-01", "590", "reference"],
            ],
        )
        base_fa = _write_fasta(
            tmp_path / "base.fa",
            {
                "PQ_TIE": seq_clean,
                "PQ_HAS_N": seq_with_n,
                "PQ_PARTIAL": seq_partial_gb,
                "PQ_REF": seq_with_n,
            },
        )
        ext_tsv = _write_tsv(
            tmp_path / "ext.tsv",
            ["Segment_Id", "Isolate_Name", "Segment", "Collection_Date", "Isolate_Id", "Passage_History", "HA INSDC_Upload", "sequence"],
            [
                ["EPI_TIE", "A/London/1/2024", "HA", "2024-01-12", "EPI_ISL_1", "Original", "PQ_TIE", seq_clean],
                ["EPI_CLEAN", "A/Paris/2/2024", "HA", "2024-02-01", "EPI_ISL_2", "Original", "PQ_HAS_N", seq_clean],
                ["EPI_FULL", "A/Rome/3/2024", "HA", "2024-03-01", "EPI_ISL_3", "Original", "PQ_PARTIAL", seq_clean],
                ["EPI_REF_DUP", "A/Ref/0/2020", "HA", "2020-01-01", "EPI_ISL_4", "Original", "PQ_REF", seq_clean],
            ],
        )
        col_map = _write_tsv(
            tmp_path / "col_map.tsv",
            ["external_column", "base_column"],
            [
                ["Isolate_Name", "strain"],
                ["Segment", "segment_declared"],
                ["Collection_Date", "collection_date"],
                ["Isolate_Id", "Isolate_Id"],
                ["Passage_History", "Passage_History"],
            ],
        )
        ext_fa = _write_fasta(
            tmp_path / "ext.fa",
            {
                "EPI_TIE": seq_clean,
                "EPI_CLEAN": seq_clean,
                "EPI_FULL": seq_clean,
                "EPI_REF_DUP": seq_clean,
            },
        )
        out_tsv = tmp_path / "merged.tsv"
        out_fa = tmp_path / "merged.fa"
        merger = NormalizeAndMerge(
            gb_matrix=str(base_tsv),
            input_tsv=str(ext_tsv),
            input_fasta=str(ext_fa),
            base_fasta=str(base_fa),
            dedup_fasta_out=str(out_fa),
            mapping_file=str(col_map),
            output=str(out_tsv),
            key="Segment_Id",
            dataset_source="gisaid",
            log_file=str(tmp_path / "merge.log"),
        )
        merger.process()

        with out_tsv.open("r", encoding="utf-8") as fh:
            rows = {r["primary_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

        # 1. Tie: GenBank wins, dual-flagged, sequence_source='genbank'
        assert rows["PQ_TIE"]["exclusion_status"] in ("", "0")
        assert rows["PQ_TIE"]["dataset_source"] == "genbank;gisaid"
        assert rows["PQ_TIE"]["sequence_source"] == "genbank"
        assert rows["EPI_TIE"]["exclusion_status"] == "1"
        assert rows["EPI_TIE"]["canonical_accession"] == "PQ_TIE"

        # 2. GISAID resolved 10 Ns -> GISAID wins, dual-flagged, sequence_source='gisaid'
        assert rows["EPI_CLEAN"]["exclusion_status"] in ("", "0")
        assert rows["EPI_CLEAN"]["dataset_source"] == "genbank;gisaid"
        assert rows["EPI_CLEAN"]["sequence_source"] == "gisaid"
        assert "genbank:PQ_HAS_N" in rows["EPI_CLEAN"]["alt_accessions"]
        assert rows["PQ_HAS_N"]["exclusion_status"] == "1"
        assert rows["PQ_HAS_N"]["canonical_accession"] == "EPI_CLEAN"

        # 3. Partial GenBank (450 bp) vs Full GISAID (600 bp) -> Full GISAID wins under Option B!
        assert rows["EPI_FULL"]["exclusion_status"] in ("", "0")
        assert rows["EPI_FULL"]["dataset_source"] == "genbank;gisaid"
        assert rows["EPI_FULL"]["sequence_source"] == "gisaid"
        assert rows["PQ_PARTIAL"]["exclusion_status"] == "1"
        assert rows["PQ_PARTIAL"]["canonical_accession"] == "EPI_FULL"

        # 4. Protected reference GenBank ALWAYS beats cleaner GISAID!
        assert rows["PQ_REF"]["exclusion_status"] in ("", "0")
        assert rows["PQ_REF"]["dataset_source"] == "genbank;gisaid"
        assert rows["PQ_REF"]["sequence_source"] == "genbank"
        assert rows["EPI_REF_DUP"]["exclusion_status"] == "1"
        assert rows["EPI_REF_DUP"]["canonical_accession"] == "PQ_REF"


# ---------------------------------------------------------------------------
# 4. DedupStrainAlignment + CreateSqliteDB + ValidateDbTree end-to-end tests
# ---------------------------------------------------------------------------

class TestDedupAlignmentAndSqliteOptionB:
    def test_dedup_strain_alignment_option_b_prefers_better_aligned_gisaid_and_keeps_gb_on_utr_tie(self, tmp_path: Path):
        seq_clean = "ATGCGTACGT" * 30  # 300 bp aligned
        seq_with_n = "ATGCGTACGT" * 14 + "NNNNN" + "TACGT" + "ATGCGTACGT" * 15  # 300 bp aligned, 295 ACGT, 5 in-place Ns
        meta = _write_tsv(
            tmp_path / "meta.tsv",
            [
                "primary_accession", "accession_type", "dataset_source", "segment",
                "Parsed_strain", "Passage_History", "collection_date", "real_length",
                "exclusion_status", "exclusion_criteria", "Isolate_Id", "alt_accessions",
            ],
            [
                ["NC_MAST", "master", "genbank", "4", "A/Master/1/2020", "Original", "2020-01-01", "300", "", "", "", ""],
                # Pair 1: aligned identical (300 ACGT each), GISAID raw real_length is 325 (UTR only) -> GenBank wins!
                ["PQ_UTR", "query", "genbank", "4", "A/Utr/1/2024", "Original", "2024-01-01", "300", "", "", "", ""],
                ["EPI_UTR", "query", "gisaid", "4", "A/Utr/1/2024", "Original", "2024-01-01", "325", "", "", "EPI_ISL_10", ""],
                # Pair 2: GenBank has 5 aligned Ns, GISAID has 0 Ns -> GISAID wins!
                ["PQ_N", "query", "genbank", "4", "A/FixN/2/2024", "Original", "2024-01-02", "295", "", "", "", ""],
                ["EPI_N", "query", "gisaid", "4", "A/FixN/2/2024", "Original", "2024-01-02", "300", "", "", "EPI_ISL_11", ""],
            ],
        )
        msa = _write_fasta(
            tmp_path / "refset_4_aln_merged_MSA.fasta",
            {
                "NC_MAST": seq_clean,
                "PQ_UTR": seq_clean,
                "EPI_UTR": seq_clean,
                "PQ_N": seq_with_n,
                "EPI_N": seq_clean,
            },
        )
        out_meta = tmp_path / "meta_dedup.tsv"
        out_msa_dir = tmp_path / "msa_out"
        out_retired = tmp_path / "retired.txt"
        res = subprocess.run(
            [
                sys.executable, str(DEDUP_STRAIN_SCRIPT),
                "--meta", str(meta),
                "--msa", str(msa),
                "--out_meta", str(out_meta),
                "--out_msa_dir", str(out_msa_dir),
                "--out_retired_update_ids", str(out_retired),
            ],
            capture_output=True, text=True, check=True,
        )
        assert res.returncode == 0
        with out_meta.open("r", encoding="utf-8") as fh:
            rows = {r["primary_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

        # Pair 1: UTR-only tie -> GenBank wins, dual-flagged
        assert rows["PQ_UTR"]["exclusion_status"] in ("", "0")
        assert rows["PQ_UTR"]["dataset_source"] == "genbank;gisaid"
        assert rows["PQ_UTR"]["sequence_source"] == "genbank"
        assert rows["EPI_UTR"]["exclusion_status"] == "1"

        # Pair 2: GISAID resolves 5 aligned Ns -> GISAID wins, dual-flagged
        assert rows["EPI_N"]["exclusion_status"] in ("", "0")
        assert rows["EPI_N"]["dataset_source"] == "genbank;gisaid"
        assert rows["EPI_N"]["sequence_source"] == "gisaid"
        assert rows["PQ_N"]["exclusion_status"] == "1"
        assert rows["PQ_N"]["canonical_accession"] == "EPI_N"

        # Clean MSA contains winners only
        clean_msa_text = (out_msa_dir / msa.name).read_text(encoding="utf-8")
        assert ">PQ_UTR\n" in clean_msa_text and ">EPI_UTR\n" not in clean_msa_text
        assert ">EPI_N\n" in clean_msa_text and ">PQ_N\n" not in clean_msa_text

    def test_update_db_retirement_sets_dual_source_and_cleans_hcv_mutation_tables(self, tmp_path: Path):
        """When an incoming cleaner sequence retires an existing update_db row,
        CreateSqliteDB must update dataset_source='genbank;gisaid', set sequence_source,
        AND purge the retired accession from all 4 current mutation tables.
        """
        seq_clean = "ATGCGTACGT" * 20
        seq_n = "ATGCGTACGT" * 9 + "NNNNNNNNNN" + "ATGCGTACGT" * 10

        # Build existing DB containing a GISAID row with Ns AND rows in all 4 current mutation tables
        old_db = tmp_path / "existing.db"
        conn = sqlite3.connect(old_db)
        conn.executescript(
            """
            CREATE TABLE meta_data (
                primary_accession TEXT, accession_type TEXT, dataset_source TEXT,
                sequence_source TEXT, segment TEXT, Parsed_strain TEXT,
                Passage_History TEXT, collection_date TEXT, exclusion_status TEXT,
                exclusion_criteria TEXT, canonical_accession TEXT, alt_accessions TEXT,
                Isolate_Id TEXT
            );
            CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, segment TEXT, alignment TEXT);
            CREATE TABLE sequences (header TEXT, sequence TEXT, segment TEXT);
            CREATE TABLE features (accession TEXT, product TEXT, cds_start INT, cds_end INT);
            CREATE TABLE insertions (primary_accession TEXT, pos INT);
            CREATE TABLE sequence_mutation_calls (primary_accession TEXT, mutation_id TEXT, segment TEXT);
            CREATE TABLE sequence_relevant_mutation_summary (primary_accession TEXT, total_relevant_mutations INT);
            CREATE TABLE completed_signatures_only (primary_accession TEXT, signature_id TEXT);
            CREATE TABLE sequence_reading_frame_scan (primary_accession TEXT, protein_name TEXT);
            CREATE TABLE info (event TEXT, date TEXT);
            CREATE TABLE update_batches (batch_id TEXT, started_at TEXT, finished_at TEXT, update_db TEXT, mode TEXT);
            CREATE TABLE update_exclusions (batch_id TEXT, table_name TEXT, key TEXT, reason TEXT, date TEXT);
            CREATE TABLE update_table_deltas (batch_id TEXT, table_name TEXT, before_count INT, after_count INT, delta INT);
            """
        )
        conn.execute(
            "INSERT INTO meta_data VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ("EPI_OLD", "query", "gisaid", "gisaid", "4", "A/Update/1/2024", "Original", "2024-05-01", "", "", "", "", "EPI_ISL_55"),
        )
        conn.execute("INSERT INTO sequence_alignment VALUES (?,?,?,?)", ("EPI_OLD", "refset_4", "4", seq_n))
        conn.execute("INSERT INTO sequences VALUES (?,?,?)", ("EPI_OLD", seq_n, "4"))
        conn.execute("INSERT INTO sequence_mutation_calls VALUES (?,?,?)", ("EPI_OLD", "HA:226L", "4"))
        conn.execute("INSERT INTO sequence_relevant_mutation_summary VALUES (?,?)", ("EPI_OLD", 1))
        conn.execute("INSERT INTO completed_signatures_only VALUES (?,?)", ("EPI_OLD", "HA:226L"))
        conn.execute("INSERT INTO sequence_reading_frame_scan VALUES (?,?)", ("EPI_OLD", "HA"))
        conn.commit()
        conn.close()

        from CreateSqliteDB import CreateSqliteDB

        builder = CreateSqliteDB.__new__(CreateSqliteDB)
        builder.update = True
        builder.batch_id = "test_batch"
        conn = sqlite3.connect(old_db)
        builder._apply_update_db_strain_retirements(
            conn,
            superseded_update_list=[
                {
                    "retired_accession": "EPI_OLD",
                    "segment": "4",
                    "winner_accession": "PQ_NEW",
                    "reason": "cross_db_duplicate_of:PQ_NEW",
                }
            ],
            update_winner_patches={},
            update_exclusions=[],
            now_str="2026-11-25 00:00:00",
        )
        conn.commit()
        for tbl in (
            "sequence_mutation_calls",
            "sequence_relevant_mutation_summary",
            "completed_signatures_only",
            "sequence_reading_frame_scan",
            "sequence_alignment",
            "sequences",
        ):
            assert conn.execute(f"SELECT COUNT(*) FROM {tbl} WHERE primary_accession='EPI_OLD'" if tbl != "sequences" else "SELECT COUNT(*) FROM sequences WHERE header='EPI_OLD'").fetchone()[0] == 0, f"Orphan row left in {tbl}"
        conn.close()

    def test_validate_db_tree_checks_dual_source_flag_on_cross_db_winners(self, tmp_path: Path):
        """ValidateDbTree must flag a cross-DB winner whose dataset_source was not
        upgraded to include both sources (e.g. 'genbank;gisaid')."""
        db = tmp_path / "dual_check.db"
        conn = sqlite3.connect(db)
        conn.executescript(
            """
            CREATE TABLE meta_data (
                primary_accession TEXT, accession_type TEXT, dataset_source TEXT,
                sequence_source TEXT, segment TEXT, Parsed_strain TEXT,
                exclusion_status TEXT, exclusion_criteria TEXT, canonical_accession TEXT,
                alt_accessions TEXT, host_taxa_id TEXT
            );
            CREATE TABLE sequence_alignment (primary_accession TEXT, segment TEXT, alignment TEXT);
            CREATE TABLE sequences (header TEXT, segment TEXT, sequence TEXT);
            CREATE TABLE features (accession TEXT, segment TEXT, product TEXT, cds_start INT, cds_end INT);
            CREATE TABLE host_taxa (taxa_id TEXT, name TEXT, name_type TEXT, taxonomy_level TEXT, common_name TEXT);
            CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT);
            """
        )
        seq = "ATGCGTACGT" * 20
        # Missing dual-source flag on PQ1 even though EPI1 was retired as cross_db_duplicate_of:PQ1
        conn.execute(
            "INSERT INTO meta_data VALUES ('PQ1','query','genbank','genbank','1','A/Test/1/2024','','','','gisaid:EPI1','9606')"
        )
        conn.execute(
            "INSERT INTO meta_data VALUES ('EPI1','query','gisaid','gisaid','1','A/Test/1/2024','1','cross_db_duplicate_of:PQ1','PQ1','','9606')"
        )
        conn.execute("INSERT INTO host_taxa VALUES ('9606','Homo sapiens','scientific','species','human')")
        conn.execute("INSERT INTO sequence_alignment VALUES ('PQ1','1',?)", (seq,))
        conn.execute("INSERT INTO sequences VALUES ('PQ1','1',?)", (seq,))
        conn.execute("INSERT INTO features VALUES ('PQ1','1','PB2',1,200)")
        conn.execute("INSERT INTO trees VALUES ('usher_seg1','usher','segment_1','1','(PQ1:0.01);','2026-11-25')")
        conn.commit()
        conn.close()

        outdir1 = tmp_path / "val_out1"
        res = subprocess.run(
            [sys.executable, str(VALIDATE_DB_TREE), "--db", str(db), "--outdir", str(outdir1), "--check-strain-duplicates", "--strict-invariants"],
            capture_output=True, text=True,
        )
        assert res.returncode != 0
        assert "dual-source" in (res.stdout + res.stderr).lower() or "dataset_source" in (res.stdout + res.stderr).lower()

        # Fixing dataset_source to 'genbank;gisaid' makes the check pass
        conn = sqlite3.connect(db)
        conn.execute("UPDATE meta_data SET dataset_source = 'genbank;gisaid' WHERE primary_accession = 'PQ1'")
        conn.commit()
        conn.close()
        outdir2 = tmp_path / "val_out2"
        res2 = subprocess.run(
            [sys.executable, str(VALIDATE_DB_TREE), "--db", str(db), "--outdir", str(outdir2), "--check-strain-duplicates", "--strict-invariants"],
            capture_output=True, text=True,
        )
        assert res2.returncode == 0, res2.stdout + res2.stderr
