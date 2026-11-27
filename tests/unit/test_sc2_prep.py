"""Unit tests for SARS-CoV-2 DB prep: GFF/frameshift handling, EpiCoV GISAID tidy/merge, QC filters, and Pango reference builder."""

import csv
from pathlib import Path

import pandas as pd
import pytest

import accession_utils
import AddPosthocClustersToDb as posthoc_mmseqs
import BuildReferenceAlignment as bra
import BuildSc2ReferenceSet as sc2_refs
from FilterAndExtractSequences import FilterAndExtractSequences
from GenBankFetcher import GenBankFetcher
from GenBankParser import GenBankParser
from GffToDictionary import GffDictionary
from gisaid_tidy import GISAIDTidy, parse_gisaid_location
from merge_into_gB_matrix import NormalizeAndMerge
from ValidateMatrix import ValidateMatrix

REPO_ROOT = Path(__file__).resolve().parents[2]
SC2_GFF = REPO_ROOT / "generic" / "sc2" / "Gff" / "NC_045512.gff3"
HCV_GFF = REPO_ROOT / "generic" / "hcv" / "Gff" / "NC_004102.gff3"
SC2_GENE_INFO = REPO_ROOT / "generic" / "sc2" / "Tables" / "gene_info.tsv"
SC2_COL_MAP = REPO_ROOT / "generic" / "sc2" / "column_mapping.tsv"


def test_sc2_gff_retains_mature_peptides_and_downstream_orfs_without_pp1a_duplicates():
    gff_dict = GffDictionary(str(SC2_GFF)).gff_dict
    cds = gff_dict["CDS"]
    products = [entry["product"] for entry in cds]

    # Polyproteins themselves must be replaced by mature peptides
    assert "ORF1ab polyprotein" not in products
    assert "ORF1a polyprotein" not in products

    # Shared pp1a/pp1ab mature peptides (nsp1..nsp10) must appear once, not twice
    assert products.count("leader protein") == 1
    assert products.count("nsp3") == 1
    assert products.count("3C-like proteinase") == 1
    assert products.count("nsp10") == 1

    # Both pp1a-specific nsp11 and pp1ab-specific nsp12..nsp16 must survive
    assert "nsp11" in products
    assert "RNA-dependent RNA polymerase" in products
    assert "helicase" in products
    assert "2'-O-ribose methyltransferase" in products

    # Standalone structural and accessory ORFs downstream of ORF1ab must all survive
    for downstream in (
        "surface glycoprotein",
        "ORF3a protein",
        "envelope protein",
        "membrane glycoprotein",
        "ORF6 protein",
        "ORF7a protein",
        "ORF7b",
        "ORF8 protein",
        "nucleocapsid phosphoprotein",
        "ORF10 protein",
    ):
        assert downstream in products, f"Missing downstream SARS-CoV-2 CDS: {downstream}"

    # Every product emitted by NC_045512.gff3 must have an entry in generic/sc2/Tables/gene_info.tsv
    gene_df = pd.read_csv(SC2_GENE_INFO, sep="\t", dtype=str, keep_default_na=False)
    known_names = set(gene_df["name"])
    assert set(products) <= known_names


def test_hcv_gff_behaviour_is_unchanged():
    products = [entry["product"] for entry in GffDictionary(str(HCV_GFF)).gff_dict["CDS"]]
    assert "polyprotein" not in products
    assert "protein F" not in products
    assert products[0] == "core protein"
    assert products[-1] == "RNA-dependent RNA polymerase NS5B"
    assert len(products) == 10


def test_choose_codon_regions_trims_frameshift_overlap_to_keep_orf1b():
    # Synthetic master mimicking SARS-CoV-2 ORF1a (1..180, 60 codons, longer) + frameshifted ORF1b (174..308, 45 codons)
    orf1a = "ATG" + "GCT" * 58 + "TAA"          # 180 nt (positions 1..180)
    tail_1b = "GGT" * 42 + "TGA"                # 129 nt (positions 181..309)
    master = orf1a + tail_1b                    # 309 nt; 174..308 is 135 nt = 45 codons in -1 frame
    cds = [
        {"start": 1, "end": 180, "strand": "+", "phase": "0", "product": "ORF1a"},
        {"start": 174, "end": 308, "strand": "+", "phase": "0", "product": "ORF1ab"},
    ]
    chosen, _ = bra.choose_codon_regions(cds, master)
    assert len(chosen) == 2
    assert (chosen[0]["start"], chosen[0]["end"]) == (1, 180)
    # Trimmed in whole-codon steps past 180 -> 174 + 3*3 = 183..308 (42 codons)
    assert (chosen[1]["start"], chosen[1]["end"]) == (183, 308)
    assert (chosen[1]["end"] - chosen[1]["start"] + 1) % 3 == 0


def test_pango_major_parent_derivation():
    alias_key = {
        "A": "",
        "B": "",
        "BA": "B.1.1.529",
        "JN": "B.1.1.529.2.86.1",
        "KP": "B.1.1.529.2.86.1.1.11.1",
        "XBB": ["BJ.1", "BM.1.1.1"],
        "XFG": ["LF.7", "LP.8.1.2"],
    }
    assert sc2_refs.derive_pango_major_parent("A", alias_key) == "A"
    assert sc2_refs.derive_pango_major_parent("B.1.1.7", alias_key) == "B.1"
    assert sc2_refs.derive_pango_major_parent("BA.5.2.1", alias_key) == "BA.5"
    assert sc2_refs.derive_pango_major_parent("JN.1.4", alias_key) == "JN.1"
    assert sc2_refs.derive_pango_major_parent("KP.3.1.1", alias_key) == "KP.3"
    assert sc2_refs.derive_pango_major_parent("XBB.1.5", alias_key) == "XBB"
    assert sc2_refs.derive_pango_major_parent("XFG.3.1", alias_key) == "XFG"


def test_epicov_gisaid_tidy_and_merge_with_hybrid_refs(tmp_path: Path):
    raw_dir = tmp_path / "gisaid_raw"
    tidy_dir = tmp_path / "gisaid_tidy"
    raw_dir.mkdir()

    # Mix of bulk tar.xz FASTA header (no EPI_ISL in header) and explicit EPI_ISL header
    (raw_dir / "metadata.tsv").write_text(
        "Virus name\tAccession ID\tCollection date\tSubmission date\tLocation\tHost\n"
        "hCoV-19/Belgium/UZA-UA-48355442/2023\tEPI_ISL_17000001\t2023-01-10\t2023-01-20\tEurope / Belgium / Antwerp\tHuman\n"
        "hCoV-19/England/QRY-1/2023\tEPI_ISL_17000002\t2023-02-05\t2023-02-12\tEurope / United Kingdom / England\tHuman\n",
        encoding="utf-8",
    )
    (raw_dir / "sequences.fasta").write_text(
        ">hCoV-19/Belgium/UZA-UA-48355442/2023|2023-01-10|2023-01-20\n"
        "ATGCATGCATGC\n"
        ">hCoV-19/England/QRY-1/2023|EPI_ISL_17000002|2023-02-05\n"
        "ATGCNNATGCAT\n",
        encoding="utf-8",
    )

    ref_list = tmp_path / "refs.tsv"
    ref_list.write_text(
        "primary_accession\taccession_type\tsegment\tgenotype\tsubtype\n"
        "NC_045512\tmaster\t1\tB\tB\n"
        "EPI_ISL_17000001\treference\t1\tXBZ\tXBZ\n",
        encoding="utf-8",
    )

    tidy = GISAIDTidy(
        data_dir=str(raw_dir),
        output_dir=str(tidy_dir),
        filetype="tsv",
        log_file=str(tmp_path / "tidy.log"),
        ref_list=str(ref_list),
    )
    tidy.process()

    tidy_meta = pd.read_csv(tidy_dir / "metadata.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert set(tidy_meta["Accession_Id"]) == {"EPI_ISL_17000001", "EPI_ISL_17000002"}
    assert set(tidy_meta["Country"]) == {"Belgium", "United Kingdom"}

    base_gb = tmp_path / "gB_matrix_validated.tsv"
    base_gb.write_text(
        "primary_accession\taccession_version\tlocus\tgi_number\taccession_type\tsegment\tdivision\t"
        "length\treal_length\ta\tt\tg\tc\tn\tcountry\tgeo_loc\thost\tstrain\t"
        "collection_date\tcollection_day\tcollection_mon\tcollection_year\texclusion_status\texclusion_criteria\n"
        "NC_045512\tNC_045512.2\tNC_045512\tNC_045512\tmaster\t1\tVRL\t"
        "12\t12\t3\t3\t3\t3\t0\tChina\tWuhan\tHomo sapiens\tWuhan-Hu-1\t"
        "2019-12-01\t1\t12\t2019\t0\t\n",
        encoding="utf-8",
    )

    merged_out = tmp_path / "gB_matrix_merged.tsv"
    NormalizeAndMerge(
        gb_matrix=str(base_gb),
        input_tsv=str(tidy_dir / "metadata.tsv"),
        input_fasta=str(tidy_dir / "all_nuc.fas"),
        mapping_file=str(SC2_COL_MAP),
        output=str(merged_out),
        key="Accession_Id",
        dataset_source="gisaid",
        log_file=str(tmp_path / "merge.log"),
        ref_list=str(ref_list),
        require_refs=True,
    ).process()

    merged_df = pd.read_csv(merged_out, sep="\t", dtype=str, keep_default_na=False).set_index("primary_accession")
    assert merged_df.loc["EPI_ISL_17000001", "accession_type"] == "reference"
    assert merged_df.loc["EPI_ISL_17000002", "accession_type"] == "query"
    assert merged_df.loc["EPI_ISL_17000001", "exclusion_status"] == "0"
    assert merged_df.loc["EPI_ISL_17000001", "segment"] == "1"
    assert merged_df.loc["EPI_ISL_17000002", "n"] == "2"
    assert merged_df.loc["EPI_ISL_17000001", "collection_year"] == "2023"


def test_filter_and_extract_enforces_max_n_and_ambig_on_queries_only(tmp_path: Path):
    matrix = tmp_path / "gB_matrix.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    matrix.write_text(
        "division\tgi_number\tlength\treal_length\tn\texclusion_status\texclusion_criteria\taccession_type\n"
        "VRL\tNC_045512\t10\t8\t2\t0\t\tmaster\n"
        "VRL\tQ_CLEAN\t10\t10\t0\t0\t\tquery\n"
        "VRL\tQ_HIGH_N\t10\t8\t2\t0\t\tquery\n"
        "VRL\tQ_WOBBLE\t10\t8\t0\t0\t\tquery\n",
        encoding="utf-8",
    )
    seqs.write_text(
        ">NC_045512\nATGCNNATGC\n"
        ">Q_CLEAN\nATGCATGCAT\n"
        ">Q_HIGH_N\nATGCNNATGC\n"
        ">Q_WOBBLE\nATGCRWATGC\n",
        encoding="utf-8",
    )
    refs.write_text("NC_045512\tmaster\t1\n", encoding="utf-8")

    FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=10,
        real_length=8,
        prop_ambigious=None,
        segmented_virus="N",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
        max_n_count=1,
        max_ambig_fraction=0.1,
    ).process()

    ref_fa = (tmp_path / "Sequences" / "ref_seq.fa").read_text(encoding="utf-8")
    query_fa = (tmp_path / "Sequences" / "query_seq.fa").read_text(encoding="utf-8")
    assert ">NC_045512" in ref_fa
    assert ">Q_CLEAN" in query_fa
    assert ">Q_HIGH_N" not in query_fa
    assert ">Q_WOBBLE" not in query_fa


def test_genbank_fetcher_skips_gisaid_refs_and_applies_entrez_filter(tmp_path: Path, monkeypatch):
    assert accession_utils.is_gisaid_accession("EPI_ISL_402124")
    assert accession_utils.is_gisaid_accession("EPI317")
    assert not accession_utils.is_gisaid_accession("NC_045512")

    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text(
        "primary_accession\taccession_type\tsegment\n"
        "NC_045512\tmaster\t1\n"
        "EPI_ISL_402124\treference\t1\n",
        encoding="utf-8",
    )
    fetcher = GenBankFetcher(
        taxid="2697049",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
        ref_list=str(ref_file),
        entrez_filter="29500:30200[SLEN]",
    )
    assert "29500%3A30200%5BSLEN%5D" in fetcher._search_term_encoded()

    called_urls = []

    class _Resp:
        status_code = 200
        text = "<GBSet></GBSet>"
        def raise_for_status(self): pass

    monkeypatch.setattr("GenBankFetcher.requests.get", lambda url: (called_urls.append(url), _Resp())[1])
    monkeypatch.setattr(fetcher, "save_data", lambda *_: None)
    fetcher.fetch_genbank_data(["OZ494601"])
    assert any("NC_045512" in u for u in called_urls)
    assert not any("EPI_ISL_402124" in u for u in called_urls)


def test_gisaid_tidy_supports_tar_xz_and_ignores_lineages_and_csv_when_tsv_exists(tmp_path: Path):
    import io
    import tarfile

    raw_dir = tmp_path / "SC2_prep"
    tidy_dir = tmp_path / "tidy"
    raw_dir.mkdir()

    # Decoy files matching dev/SC2_prep/ layout (lineages.csv + incomplete metadata.csv)
    (raw_dir / "lineages.csv").write_text("taxon,lineage\nhCoV-19/Belgium/A/2023,JN.1\n", encoding="utf-8")
    (raw_dir / "metadata.csv").write_text(
        "accession_id,virus_name,collection_date\nEPI_ISL_99999999,hCoV-19/Nowhere/DECOY/2023,2023-01-01\n",
        encoding="utf-8",
    )

    # Colleague's full GISAID metadata_tsv_*.tar.xz archive
    tsv_bytes = (
        "Virus name\tAccession ID\tCollection date\tSubmission date\tLocation\tHost\n"
        "hCoV-19/Scotland/CVR-1/2026\tEPI_ISL_18000001\t2026-11-10\t2026-11-12\tEurope / United Kingdom / Scotland\tHuman\n"
    ).encode("utf-8")
    tar_path = raw_dir / "metadata_tsv_2026_11_15.tar.xz"
    with tarfile.open(tar_path, mode="w:xz") as tf:
        info = tarfile.TarInfo(name="metadata.tsv")
        info.size = len(tsv_bytes)
        tf.addfile(info, io.BytesIO(tsv_bytes))

    (raw_dir / "sequences.fasta").write_text(
        ">hCoV-19/Scotland/CVR-1/2026|2026-11-10|2026-11-12\nATGCATGCATGC\n",
        encoding="utf-8",
    )

    GISAIDTidy(
        data_dir=str(raw_dir),
        output_dir=str(tidy_dir),
        filetype="tsv",
        log_file=str(tmp_path / "tidy.log"),
    ).process()

    tidy_meta = pd.read_csv(tidy_dir / "metadata.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(tidy_meta["Accession_Id"]) == ["EPI_ISL_18000001"]
    assert list(tidy_meta["Country"]) == ["United Kingdom"]
    assert list(tidy_meta["Host"]) == ["Human"]


def test_build_sc2_refs_strict_qc_rejects_high_n_and_verifies_gisaid_fasta(tmp_path: Path):
    meta_tsv = tmp_path / "metadata.tsv"
    meta_tsv.write_text(
        "Virus name\tAccession ID\tCollection date\tSequence length\tN-Content\tIs complete?\n"
        "hCoV-19/USA/BAD-N/2024\tEPI_ISL_1001\t2024-01-01\t29850\t0.02\tTrue\n"
        "hCoV-19/USA/GOOD-1/2024\tEPI_ISL_1002\t2024-01-02\t29860\t0.001\tTrue\n"
        "hCoV-19/USA/SNEAKY-FASTA-N/2024\tEPI_ISL_1003\t2024-01-03\t29860\t0.0001\tTrue\n",
        encoding="utf-8",
    )
    taxon_to_lin = {
        "USA/BAD-N/2024": "JN.1",
        "USA/GOOD-1/2024": "JN.1",
        "USA/SNEAKY-FASTA-N/2024": "KP.2",
    }
    cands = sc2_refs.resolve_from_gisaid_metadata(
        str(meta_tsv), taxon_to_lin, min_length=29500, max_n=150, max_ambig_fraction=0.005
    )
    # EPI_ISL_1001 rejected immediately on metadata N-Content=0.02 (>0.005)
    assert [c["accession"] for c in cands["JN.1"]] == ["EPI_ISL_1002"]

    # Sequence-level verification catches EPI_ISL_1003 having 300 Ns in sequences.fasta
    fasta = tmp_path / "sequences.fasta"
    clean_seq = "ATGC" * 7465          # 29,860 ATGC, 0 N
    dirty_seq = "ATGC" * 7390 + "N" * 300  # 29,860 bp, 300 Ns (>150)
    fasta.write_text(
        f">hCoV-19/USA/GOOD-1/2024|2024-01-02|2024-01-05\n{clean_seq}\n"
        f">hCoV-19/USA/SNEAKY-FASTA-N/2024|2024-01-03|2024-01-05\n{dirty_seq}\n",
        encoding="utf-8",
    )
    verified = sc2_refs.verify_gisaid_candidates_from_fasta(
        str(fasta), cands, min_length=29500, min_real_length=29350, max_n=150, max_ambig_fraction=0.005
    )
    assert "JN.1" in verified and verified["JN.1"]["accession"] == "EPI_ISL_1002"
    assert "KP.2" not in verified


def test_posthoc_clusters_attach_to_sqlite_db_and_iso_dates_validate(tmp_path: Path):
    import sqlite3

    assert ValidateMatrix._validate_date_str("2026-11-15") == "Yes"
    assert ValidateMatrix._validate_date_str("2026-11") == "Yes"
    assert ValidateMatrix._validate_date_str("15-Nov-2026") == "Yes"

    db_path = tmp_path / "sc2.db"
    conn = sqlite3.connect(db_path)
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    conn.executemany(
        "INSERT INTO meta_data VALUES (?, ?)",
        [("NC_045512", "master"), ("EPI_ISL_18000001", "query"), ("OZ123456", "query")],
    )
    conn.commit()
    conn.close()

    cluster_tsv = tmp_path / "clusters.tsv"
    cluster_tsv.write_text(
        "NC_045512\tNC_045512\n"
        "NC_045512\tEPI_ISL_18000001\n"
        "OZ123456\tOZ123456.1\n",
        encoding="utf-8",
    )
    col_name, n_members, n_updated = posthoc_mmseqs.attach_clusters_to_db(str(db_path), str(cluster_tsv), 0.999)
    assert col_name == "cluster_100pct"
    assert n_updated == 3

    conn = sqlite3.connect(db_path)
    rows = dict(conn.execute("SELECT primary_accession, cluster_100pct FROM meta_data").fetchall())
    conn.close()
    assert rows == {
        "NC_045512": "NC_045512",
        "EPI_ISL_18000001": "NC_045512",
        "OZ123456": "OZ123456",
    }
