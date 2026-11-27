"""Extensive TDD unit and integration tests for Workstream 4:
IAV mutation catalog conversion from dev/flu_raw_tables/, spliced multi-exon CDS
support (M2, NEP, PA-X), NP product canonicalisation, influenza profile columns,
serotype_validated H/N gating, pdm09 NS1 truncation handling, chunked per-segment
streaming parity, and segment-aware VerifyMutations.py.
"""

from __future__ import annotations

import csv
import sqlite3
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import AnnotateMutations as AM
import BuildIavMutationCatalog as BIC
import reading_frame as rf
import VerifyMutations as VM

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_FLU_DIR = REPO_ROOT / "dev" / "flu_raw_tables"


# ---------------------------------------------------------------------------
# 1. Raw flu CSV catalog conversion tests
# ---------------------------------------------------------------------------

class TestBuildIavMutationCatalog:
    def test_parses_all_eight_raw_tables_without_dropping_na_segment(self, tmp_path: Path):
        if not RAW_FLU_DIR.is_dir():
            pytest.skip("Local dev/flu_raw_tables not present")
        out_cat = tmp_path / "iav_catalog.tsv"
        out_pub = tmp_path / "iav_publications.csv"
        stats = BIC.build_iav_catalog(
            raw_dir=str(RAW_FLU_DIR),
            out_catalog_tsv=str(out_cat),
            out_publications_csv=str(out_pub),
        )
        assert out_cat.is_file()
        assert out_pub.is_file()

        df = pd.read_csv(out_cat, sep="\t", dtype=str, keep_default_na=False)
        pubs = pd.read_csv(out_pub, dtype=str, keep_default_na=False)

        # All 10 curated proteins MUST be present (especially NA — NOT dropped as null!)
        proteins = set(df["protein_name"])
        assert proteins == {"HA", "NA", "PB2", "PB1", "PA", "NP", "M1", "M2", "NS1", "NEP"}
        assert len(df[df["protein_name"] == "NA"]) >= 44
        assert len(df[df["protein_name"] == "M2"]) >= 18
        assert len(df[df["protein_name"] == "NEP"]) >= 11
        assert len(df[df["protein_name"] == "HA"]) >= 540

        # HCV contract required columns must all exist
        for col in (
            "mutation_id", "protein_name", "segment", "aa_position", "alt_residue",
            "reference_accession", "mutation_type", "signature_id", "signature_kind",
            "combination_id", "combination_size", "phenotype", "genotype",
            "resistance_category", "drug", "evidence_id", "data_source",
            "evidence_url", "evidence_label", "relevant_genotypes", "wild_type_residues",
            "serotypes_tested", "h3_numbering", "h5_numbering",
        ):
            assert col in df.columns, f"Missing HCV-compatible column {col}"

        # H3-to-AB573800 coordinate mapping landmarks
        q226l = df[(df["protein_name"] == "HA") & (df["h3_numbering"] == "Q226L")]
        assert not q226l.empty
        assert set(q226l["aa_position"]) == {"223"}  # H3 226 -> AB573800 master pos 223!
        assert set(q226l["alt_residue"]) == {"L"}

        g228s = df[(df["protein_name"] == "HA") & (df["h3_numbering"] == "G228S")]
        assert not g228s.empty
        assert set(g228s["aa_position"]) == {"225"}

        # HA deletion L133∆ maps cleanly to alt_residue='-' (alongside L133V/R/K substitutions)
        l133del = df[(df["protein_name"] == "HA") & (df["h3_numbering"].str.contains("L133∆|L133Δ", regex=True))]
        assert not l133del.empty
        assert set(l133del["alt_residue"]) == {"-"}

        # NS1 stop mutations survive cleanly
        k216stop = df[(df["protein_name"] == "NS1") & (df["aa_position"] == "216") & (df["alt_residue"] == "*")]
        assert not k216stop.empty

        # M2 PMID-only row ("No DOI - PMID - 11252668") resolves to pmid:11252668
        m2_93 = df[(df["protein_name"] == "M2") & (df["aa_position"] == "93")]
        assert not m2_93.empty
        assert any("11252668" in eid for eid in m2_93["evidence_id"])

        # Publications table has 1 row per unique evidence_id matching catalog
        assert len(pubs) == pubs["evidence_id"].nunique()
        cat_eids = {e for e in df["evidence_id"] if e}
        assert cat_eids.issubset(set(pubs["evidence_id"]))

        # Wild-type scoping protects seasonal human H3N2 from PB2:E627K flood while firing on H5N1
        e627k = df[(df["protein_name"] == "PB2") & (df["aa_position"] == "627") & (df["alt_residue"] == "K")]
        assert not e627k.empty
        assert any("H3N2:K" in wt for wt in e627k["wild_type_residues"])
        # PB2:M631L is NOT restricted only to H10N7 (must fire on H5N1 dairy cattle!)
        m631l = df[(df["protein_name"] == "PB2") & (df["aa_position"] == "631") & (df["alt_residue"] == "L")]
        assert not m631l.empty
        assert set(m631l["relevant_genotypes"]) == {""}


# ---------------------------------------------------------------------------
# 2. Multi-exon spliced/frameshifted CDS & reading_frame tests
# ---------------------------------------------------------------------------

class TestSplicedCdsAndReadingFrame:
    def test_m2_spliced_two_exon_codon_grid_and_split_codon_9(self):
        # M2 master layout: exon 1 = 1..26 (8 codons + 2 bases of codon 9),
        # intron = 27..714, exon 2 = 715..982 (1st base finishes codon 9, then codons 10..97)
        coord_map = {pos: pos - 1 for pos in range(1, 983)}
        exons = [(1, 26), (715, 982)]
        grid = rf.codon_grid(coord_map, 1, 982, exons=exons)
        assert grid[1] == (0, 1, 2)
        assert grid[8] == (21, 22, 23)
        # Codon 9 splits across splice junction: nt 25 (col 24), nt 26 (col 25), nt 715 (col 714)!
        assert grid[9] == (24, 25, 714)
        # Codon 10 starts at nt 716 (col 715)
        assert grid[10] == (715, 716, 717)
        # Codon 31 (S31N!) is 21 codons after codon 10 -> nt 716 + 21*3 = 779 (cols 778,779,780)
        assert grid[31] == (778, 779, 780)

    def test_reading_frame_does_not_falsely_displace_m2_exon2(self):
        # Construct realistic 982-nt segment 7 sequence with valid M2 exons (1..26 and 715..982)
        # Exon 1: ATG + 7x GCT (A) + AG (first 2 nt of codon 9 -> AGC = S) = 26 nt
        # Intron: 688 nt of 'TAA' stops (in M1 reading frame, irrelevant to M2!)
        # Exon 2: C (3rd nt of codon 9) + 21x GCT (codons 10..30) + AAT (codon 31 = N!) + 66x GCT + TAA = 268 nt
        exon1 = "ATG" + "GCT" * 7 + "AG"
        intron = "TAA" * 229 + "T"  # 688 nt -> exon2 starts at 1-based nt 715 (0-based 714)
        exon2 = "C" + "GCT" * 21 + "AAT" + "GCT" * 66 + "TAA"
        seq = exon1 + intron + exon2
        assert len(exon1) == 26
        assert len(seq) == 982

        coord_map = {pos: pos - 1 for pos in range(1, len(seq) + 1)}
        grid = rf.codon_grid(coord_map, 1, len(seq), exons=[(1, 26), (715, 982)])
        rm = rf.RecordMapping(seq, seq)
        pf = rf.ProteinFrame(rm, grid, AM.translate_codon, "X", master_alignment=seq)
        # Codon 31 in exon 2 MUST read 'N' in-frame, NOT guarded/displaced!
        res31, verdict31 = pf.read(grid[31])
        assert res31 == "N"
        assert verdict31 in (rf.IN_FRAME, rf.DEFER)
        # Codon 2 in exon 1 reads 'A' cleanly
        res2, _ = pf.read(grid[2])
        assert res2 == "A"

    def test_pdm09_ns1_natural_terminal_truncation_stop_does_not_guard_upstream(self):
        # 230-aa NS1 master, pdm09 query has natural stop TAA at codon 220 + a 1-codon aligner split at codon 150
        master_codons = ["ATG"] + ["GCT"] * 228 + ["TAA"]
        master_seq = "".join(master_codons)
        # Query has P42S (codon 42 = TCT), natural stop at codon 220 (TAA)
        query_codons = list(master_codons)
        query_codons[41] = "TCT"  # pos 42 = S
        query_codons[219] = "TAA"  # pos 220 = * (natural pdm09 truncation)
        query_seq = "".join(query_codons)
        # Simulate minor aligner artefact at codon 150 (1-base gap + shifted base)
        aln_q = query_seq[:447] + "-" + query_seq[447:450] + query_seq[451:]
        # Master padded with 1 gap at col 447 to keep equal length
        aln_m = master_seq[:447] + "-" + master_seq[447:]
        # Wait: make simple alignment with 1 deleted base at codon 150 (split codon) so split is non-empty
        aln_q2 = query_seq[:448] + "-" + query_seq[449:]
        coord_map = {pos: pos - 1 for pos in range(1, len(master_seq) + 1)}
        grid = rf.codon_grid(coord_map, 1, len(master_seq))
        rm = rf.RecordMapping(aln_q2, query_seq)
        pf = rf.ProteinFrame(rm, grid, AM.translate_codon, "X", master_alignment=master_seq)
        # Codon 42 upstream of natural terminal truncation MUST still read 'S', NOT be guarded as frameshift!
        res42, verdict42 = pf.read(grid[42])
        assert res42 == "S", f"pdm09 NS1 pos 42 wrongly guarded: {(res42, verdict42)}"


# ---------------------------------------------------------------------------
# 3. End-to-end multi-segment IAV AnnotateMutations + VerifyMutations tests
# ---------------------------------------------------------------------------

class TestIavEndToEndAnnotationAndStreaming:
    def _create_synthetic_iav_db(self, db_path: Path) -> None:
        conn = sqlite3.connect(db_path)
        conn.executescript(
            """
            CREATE TABLE meta_data (
                primary_accession TEXT, accession_type TEXT, segment TEXT,
                serotype_validated TEXT, nearest_reference_genotype TEXT,
                nearest_reference_subtype TEXT, exclusion_status TEXT, exclusion_criteria TEXT
            );
            CREATE TABLE genes (name TEXT, display_name TEXT, description TEXT);
            CREATE TABLE features (
                accession TEXT, segment TEXT, product TEXT, feature_type TEXT,
                cds_start INT, cds_end INT, cds_start_OG_seq INT, cds_end_OG_seq INT
            );
            CREATE TABLE sequence_alignment (primary_accession TEXT, segment TEXT, alignment TEXT);
            CREATE TABLE sequences (header TEXT, segment TEXT, sequence TEXT);
            """
        )
        # Genes table matching generic/influenza/Tables/gene_info.tsv (plus NP alias)
        conn.executemany(
            "INSERT INTO genes VALUES (?,?,?)",
            [
                ("PB2", "PB2", "Polymerase PB2"),
                ("NP", "NP", "Nucleoprotein"),
                ("M2", "M2", "Matrix protein 2"),
            ],
        )

        # Segment 1 (PB2): 15 codons (45 nt), pos 10 = E (GAA) in master -> K (AAA) in H5N1 & H3N2
        pb2_m = "ATG" + "GCT" * 8 + "GAA" + "GCT" * 4 + "TAA"
        pb2_k = "ATG" + "GCT" * 8 + "AAA" + "GCT" * 4 + "TAA"
        # Segment 5 (NP): 15 codons (45 nt), product literally 'nucleocapsid protein' (like NC_007369!)
        np_m = "ATG" + "GCT" * 4 + "TAT" + "GCT" * 8 + "TAA"  # pos 6 = Y
        np_h = "ATG" + "GCT" * 4 + "CAT" + "GCT" * 8 + "TAA"  # pos 6 = H
        # Segment 7 (M2 spliced): exon 1 = 1..26 (8 codons + AG), intron = 27..44 (18 nt), exon 2 = 45..72 (C + 8 codons + TAA)
        # Codon 9 = AGC (S), codon 10 = GCT (A), codon 11 = AGT (S) -> AAT (N) in query!
        m2_m = "ATG" + "GCT" * 7 + "AG" + ("TAA" * 6) + "C" + "GCT" + "AGT" + "GCT" * 6 + "TAA"
        m2_n = "ATG" + "GCT" * 7 + "AG" + ("TAA" * 6) + "C" + "GCT" + "AAT" + "GCT" * 6 + "TAA"

        meta_rows = [
            ("MAST_PB2", "master", "1", "H5N1", "", "", "", ""),
            ("Q_H5_PB2", "query", "1", "H5N1", "", "", "", ""),
            ("Q_H3_PB2", "query", "1", "H3N2", "", "", "", ""),
            ("MAST_NP", "master", "5", "H1N1", "", "", "", ""),
            ("Q_H5_NP", "query", "5", "H5N1", "", "", "", ""),
            ("MAST_M2", "master", "7", "H1N1", "", "", "", ""),
            ("Q_H5_M2", "query", "7", "H5N1", "", "", "", ""),
        ]
        conn.executemany("INSERT INTO meta_data VALUES (?,?,?,?,?,?,?,?)", meta_rows)

        feat_rows = [
            ("MAST_PB2", "1", "polymerase PB2", "CDS", 1, 45, 1, 45),
            ("MAST_NP", "5", "nucleocapsid protein", "CDS", 1, 45, 1, 45),
            # TWO spliced exons for M2 (just like NC_002016.gff3!)
            ("MAST_M2", "7", "matrix protein 2", "CDS", 1, 26, 1, 26),
            ("MAST_M2", "7", "matrix protein 2", "CDS", 45, 72, 45, 72),
        ]
        conn.executemany("INSERT INTO features VALUES (?,?,?,?,?,?,?,?)", feat_rows)

        seqs = [
            ("MAST_PB2", "1", pb2_m),
            ("Q_H5_PB2", "1", pb2_k),
            ("Q_H3_PB2", "1", pb2_k),
            ("MAST_NP", "5", np_m),
            ("Q_H5_NP", "5", np_h),
            ("MAST_M2", "7", m2_m),
            ("Q_H5_M2", "7", m2_n),
        ]
        conn.executemany("INSERT INTO sequence_alignment VALUES (?,?,?)", seqs)
        conn.executemany("INSERT INTO sequences VALUES (?,?,?)", seqs)
        conn.commit()
        conn.close()

    def _write_synthetic_catalog(self, path: Path, pub_path: Path) -> None:
        rows = [
            {
                "mutation_id": "PB2:10K", "protein_name": "PB2", "segment": "1",
                "aa_position": "10", "alt_residue": "K", "reference_accession": "MAST_PB2",
                "mutation_type": "aminoAcidSimplePolymorphism", "signature_id": "PB2:10K",
                "signature_kind": "single", "combination_id": "", "combination_size": "",
                "phenotype": "mammalian_adaptation", "genotype": "all",
                "resistance_category": "", "drug": "", "drug_category": "", "drug_producer": "",
                "any_in_vitro_evidence": "", "in_vitro_max_ec50_midpoint": "",
                "any_in_vivo_evidence": "1", "in_vivo_baseline": "", "in_vivo_treatment_emergent": "",
                "evidence_id": "doi:10.1000/pb2", "data_source": "publication",
                "evidence_url": "https://doi.org/10.1000/pb2", "evidence_label": "Test PB2 2024",
                "evidence_type": "literature", "regimens": "", "linked_evidence_ids": "", "finding_ids": "",
                "relevant_genotypes": "", "wild_type_residues": "H3N2:K:99.0;H5N1:E:99.0",
                "serotypes_tested": "H5N1; H3N2",
            },
            {
                "mutation_id": "NP:6H", "protein_name": "NP", "segment": "5",
                "aa_position": "6", "alt_residue": "H", "reference_accession": "MAST_NP",
                "mutation_type": "aminoAcidSimplePolymorphism", "signature_id": "NP:6H",
                "signature_kind": "single", "combination_id": "", "combination_size": "",
                "phenotype": "btn3a3_evasion", "genotype": "H5N1",
                "resistance_category": "", "drug": "", "drug_category": "", "drug_producer": "",
                "any_in_vitro_evidence": "1", "in_vitro_max_ec50_midpoint": "",
                "any_in_vivo_evidence": "", "in_vivo_baseline": "", "in_vivo_treatment_emergent": "",
                "evidence_id": "doi:10.1000/np", "data_source": "publication",
                "evidence_url": "https://doi.org/10.1000/np", "evidence_label": "Test NP 2024",
                "evidence_type": "literature", "regimens": "", "linked_evidence_ids": "", "finding_ids": "",
                "relevant_genotypes": "H5", "wild_type_residues": "H5N1:Y:99.0",
                "serotypes_tested": "H5N1",
            },
            {
                "mutation_id": "M2:11N", "protein_name": "M2", "segment": "7",
                "aa_position": "11", "alt_residue": "N", "reference_accession": "MAST_M2",
                "mutation_type": "aminoAcidSimplePolymorphism", "signature_id": "M2:11N",
                "signature_kind": "single", "combination_id": "", "combination_size": "",
                "phenotype": "mammalian_adaptation", "genotype": "all",
                "resistance_category": "", "drug": "", "drug_category": "", "drug_producer": "",
                "any_in_vitro_evidence": "", "in_vitro_max_ec50_midpoint": "",
                "any_in_vivo_evidence": "", "in_vivo_baseline": "", "in_vivo_treatment_emergent": "",
                "evidence_id": "doi:10.1000/m2", "data_source": "publication",
                "evidence_url": "https://doi.org/10.1000/m2", "evidence_label": "Test M2 2024",
                "evidence_type": "literature", "regimens": "", "linked_evidence_ids": "", "finding_ids": "",
                "relevant_genotypes": "", "wild_type_residues": "",
                "serotypes_tested": "H5N1",
            },
        ]
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        pd.DataFrame(
            [
                {"evidence_id": "doi:10.1000/pb2", "citation": "Test PB2 2024", "doi": "10.1000/pb2", "url": "https://doi.org/10.1000/pb2"},
                {"evidence_id": "doi:10.1000/np", "citation": "Test NP 2024", "doi": "10.1000/np", "url": "https://doi.org/10.1000/np"},
                {"evidence_id": "doi:10.1000/m2", "citation": "Test M2 2024", "doi": "10.1000/m2", "url": "https://doi.org/10.1000/m2"},
            ]
        ).to_csv(pub_path, index=False)

    def test_iav_annotation_streaming_spliced_np_serotype_gating_and_verify_mutations(self, tmp_path: Path):
        db_full = tmp_path / "iav_full.db"
        db_chunked = tmp_path / "iav_chunked.db"
        self._create_synthetic_iav_db(db_full)
        self._create_synthetic_iav_db(db_chunked)
        cat = tmp_path / "cat.tsv"
        pub = tmp_path / "pub.csv"
        self._write_synthetic_catalog(cat, pub)

        for db, chunk_sz in ((db_full, 0), (db_chunked, 1)):
            cmd = [
                sys.executable, str(REPO_ROOT / "scripts" / "AnnotateMutations.py"),
                "--db", str(db),
                "--catalog", str(cat),
                "--publications", str(pub),
                "--virus", "influenza",
                "--chunk_size", str(chunk_sz),
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            assert res.returncode == 0, f"AnnotateMutations failed (chunk={chunk_sz}):\n{res.stdout}\n{res.stderr}"

        # Verify all 6 target tables (clinical_trials optional) + exact calls!
        conn = sqlite3.connect(db_chunked)
        tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        for t in (
            "mutation_catalog", "publications",
            "sequence_relevant_mutation_summary", "completed_signatures_only",
            "sequence_mutation_calls", "sequence_reading_frame_scan",
        ):
            assert t in tables, f"Missing table {t}"

        # Influenza profile MUST preserve evidence/phenotype/wild_type columns in mutation_catalog!
        cat_cols = {r[1] for r in conn.execute("PRAGMA table_info(mutation_catalog)")}
        for col in ("evidence_id", "evidence_url", "phenotype", "wild_type_residues", "relevant_genotypes", "serotypes_tested"):
            assert col in cat_cols, f"mutation_catalog missing {col}"

        summary = {
            r[0]: r[1]
            for r in conn.execute("SELECT primary_accession, relevant_mutations_present FROM sequence_relevant_mutation_summary")
        }
        # PB2:10K fires on H5N1, suppressed as wild-type on H3N2 via serotype_validated fallback!
        assert "PB2:10K" in summary.get("Q_H5_PB2", "")
        assert "PB2:10K" not in summary.get("Q_H3_PB2", "")
        # NP:6H resolves 'nucleocapsid protein' -> 'NP' AND matches relevant_genotypes='H5' on H5N1!
        assert "NP:6H" in summary.get("Q_H5_NP", "")
        # M2:11N resolves across 2-exon splice junction!
        assert "M2:11N" in summary.get("Q_H5_M2", "")

        # Chunked streaming parity against unchunked run
        conn_full = sqlite3.connect(db_full)
        for tbl in ("sequence_relevant_mutation_summary", "completed_signatures_only", "sequence_mutation_calls", "sequence_reading_frame_scan"):
            rows_c = sorted(conn.execute(f"SELECT * FROM {tbl}").fetchall())
            rows_f = sorted(conn_full.execute(f"SELECT * FROM {tbl}").fetchall())
            assert rows_c == rows_f, f"Chunked vs unchunked mismatch in {tbl}"
        conn.close()
        conn_full.close()

        # Segment-aware VerifyMutations.py MUST verify multi-segment IAV DB cleanly!
        vres = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "VerifyMutations.py"), "--db", str(db_chunked), "--virus", "influenza"],
            capture_output=True, text=True,
        )
        assert vres.returncode == 0, f"VerifyMutations failed on IAV DB:\n{vres.stdout}\n{vres.stderr}"
