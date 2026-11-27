"""Carpet-bomb audit tests covering all 10 scripts in scripts/divergent_mode/.

Lifts all 7 previously unmeasured (0% imported) divergent_mode scripts into
coverage via real synthetic data, and pins verified bugs (AUDIT-DIV-*) with
strict xfail reproducers.
"""

import csv
import importlib.util
from pathlib import Path
import sys

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
DM_DIR = REPO_ROOT / "scripts" / "divergent_mode"


def _load_dm(name: str):
    path = DM_DIR / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"dm_{name}", str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Passing tests exercising all 10 divergent_mode scripts
# ---------------------------------------------------------------------------

def test_split_segments_and_build_protein_profiles_helpers(tmp_path):
    ss = _load_dm("split_segments")
    bpp = _load_dm("BuildProteinProfiles")

    rec_l = SeqRecord(Seq("A" * 6000), id="ACC1.1", description="Virus segment L complete")
    rec_l.annotations["organism"] = "Test bunyavirus"
    rec_l.features.append(
        SeqFeature(FeatureLocation(0, 6000), type="source", qualifiers={"segment": ["L"], "strain": ["S1"]})
    )
    seg_re = ss.segment_regex({"L", "M", "S"})
    assert ss.segment_of(rec_l, {"L", "M", "S"}, seg_re) == "L"
    assert ss.source_qual(rec_l, "strain") == "S1"
    assert bpp.segment_of(rec_l) == "L"

    feat = SeqFeature(FeatureLocation(0, 300), type="CDS", qualifiers={"product": ["polymerase P2"]})
    assert bpp.canonical_name(feat) == "P2"


def test_extract_cds_and_codon_aware_align_happy_path():
    edc = _load_dm("ExtractDivergentCds")
    caa = _load_dm("CodonAwareAlign")

    # 300-nt complete ORF starting ATG ... TAA (0-based lo=4, hi=304)
    nt_str = "ATG" + ("GCT" * 98) + "TAA"
    seq = Seq("TTTT" + nt_str + "AAAA")
    start, end, strand, flag = edc.resolve_orf(seq, 4, 304, strand=1, max_extend=400, allow_partial=False)
    assert flag == "complete" and strand == 1 and str(seq[start:end]) == nt_str

    nt_recs = [SeqRecord(Seq("ATGGCT"), id="TAX1")]
    aln_prots = [SeqRecord(Seq("M-A"), id="TAX1")]
    codon_recs = caa.back_project(nt_recs, aln_prots)
    assert str(codon_recs[0].seq) == "ATG---GCT"
    rep = caa.check_codon_integrity(codon_recs)
    assert rep == []


def test_select_genes_concatenate_and_summarise_end_to_end(tmp_path):
    sbg = _load_dm("SelectBackboneGenes")
    cca = _load_dm("ConcatenateCodonAlignments")
    ssb = _load_dm("SummariseSegmentBackbones")

    codon_dir = tmp_path / "codons"
    codon_dir.mkdir()
    for gene, seq in (("G1", "ATGAAATAA"), ("G2", "ATGCCCGGG")):
        (codon_dir / f"{gene}.codon.fasta").write_text(
            f">T1\n{seq}\n>T2\n{seq}\n>T3\n{seq}\n>T4\n{seq}\n",
            encoding="utf-8",
        )
    report = tmp_path / "cds_extraction_report.tsv"
    with open(report, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["group", "segment", "accession", "status"], delimiter="\t")
        w.writeheader()
        for g in ("G1", "G2"):
            for t in ("T1", "T2", "T3", "T4"):
                w.writerow({"group": g, "segment": "L", "accession": t, "status": "ok"})

    chosen = tmp_path / "chosen.txt"
    status = tmp_path / "status.tsv"
    summary = tmp_path / "summary.tsv"
    orig_argv = sys.argv[:]
    try:
        sys.argv = [
            "SelectBackboneGenes.py",
            "--codon-dir", str(codon_dir),
            "--report", str(report),
            "--segments", "L",
            "--min-coverage", "0.9",
            "-o", str(chosen),
            "--status", str(status),
            "--summary", str(summary),
        ]
        sbg.main()
    finally:
        sys.argv = orig_argv
    assert chosen.read_text(encoding="utf-8").strip() == "G1,G2"

    concat_out = tmp_path / "concat"
    try:
        sys.argv = [
            "ConcatenateCodonAlignments.py",
            "--codon-dir", str(codon_dir),
            "--genes", "G1,G2",
            "-o", str(concat_out),
        ]
        cca.main()
    finally:
        sys.argv = orig_argv
    assert (concat_out / "backbone.codon.fasta").exists()
    assert (concat_out / "backbone.gene_coordinates.tsv").exists()

    bb_sum = tmp_path / "backbone_summary.tsv"
    try:
        sys.argv = [
            "SummariseSegmentBackbones.py",
            "--status", str(status),
            "--segments", "L,S",
            "--built", "L",
            "-o", str(bb_sum),
        ]
        ssb.main()
    finally:
        sys.argv = orig_argv
    txt = bb_sum.read_text(encoding="utf-8")
    assert "L\tY" in txt and "S\tN" in txt


def test_partition_backbone_and_graft_helpers(tmp_path):
    pbd = _load_dm("PartitionByDivergence")
    bdb = _load_dm("BuildDivergentBackbone")
    gdc = _load_dm("GraftDenseClusters")

    d, cmp_n = pbd.p_distance("ACGTACGT", "ACGTAAGT")
    assert cmp_n == 8 and abs(d - 0.125) < 1e-9

    assert bdb.translate_gapped("ATGGCT---TAA") == "MA-*"
    coords_tsv = tmp_path / "coords.tsv"
    coords_tsv.write_text("gene\tnt_start\tnt_end\nG1\t1\t9\n", encoding="utf-8")
    coords = bdb.load_coords(str(coords_tsv))
    assert bdb.locate(coords, 4) == ("G1", 2, 1)

    with pytest.raises(ValueError):
        gdc.reject_option_like(rep="-BAD_ID")


# ---------------------------------------------------------------------------
# Strict xfail reproducers for verified Divergent Mode bugs
# ---------------------------------------------------------------------------

@pytest.mark.xfail(strict=True, reason="AUDIT-DIV-01: BuildProteinProfiles.segment_of hardcodes L/M/S and returns '' for RNA1/1..8")
def test_build_protein_profiles_supports_non_lms_segments():
    bpp = _load_dm("BuildProteinProfiles")
    rec = SeqRecord(Seq("A" * 1000), id="ACC1", description="segment RNA1")
    rec.features.append(SeqFeature(FeatureLocation(0, 1000), type="source", qualifiers={"segment": ["RNA1"]}))
    assert bpp.segment_of(rec) == "RNA1"


@pytest.mark.xfail(strict=True, reason="AUDIT-DIV-02: ExtractDivergentCds rejects 5'-partial (partial5) CDS ending in stop codon as internal_stop")
def test_extract_divergent_cds_accepts_partial5_orf_with_terminal_stop():
    edc = _load_dm("ExtractDivergentCds")
    # Sequence missing 5' ATG (runs straight into TAA at end): partial5, NOT internal_stop
    seq = Seq(("GCT" * 20) + "TAA")
    nt, flag = edc.resolve_orf(seq, 1, len(seq), strand=1, max_extend=30, allow_partial=True)
    assert "partial5" in flag and "partial3" not in flag
    aa = str(nt.translate(table=11))
    # Exact production check at ExtractDivergentCds.py:179:
    body = aa[:-1] if flag == "complete" else aa
    assert "*" not in body


@pytest.mark.xfail(strict=True, reason="AUDIT-DIV-03: PartitionByDivergence auto cutoff drops between-species distances > auto_ceiling before scanning gaps")
def test_partition_by_divergence_detects_band_straddling_auto_ceiling(tmp_path):
    pbd = _load_dm("PartitionByDivergence")
    # Build 4 sequences of 300 nt:
    # A, B identical (d=0.0); C, D identical (d=0.0); A/B vs C/D differ at 150/300 = 0.50 (> auto_ceiling 0.35!)
    aln = tmp_path / "aln.fasta"
    s1 = "A" * 300
    s2 = ("C" * 150) + ("A" * 150)
    aln.write_text(f">A\n{s1}\n>B\n{s1}\n>C\n{s2}\n>D\n{s2}\n", encoding="utf-8")
    outdir = tmp_path / "part"
    orig_argv = sys.argv[:]
    try:
        sys.argv = [
            "PartitionByDivergence.py",
            "-a", str(aln),
            "-o", str(outdir),
            "--cutoff", "auto",
            "--auto-ceiling", "0.35",
        ]
        pbd.main()
    finally:
        sys.argv = orig_argv
    import json
    meta = json.loads((outdir / "partition_cutoff.json").read_text(encoding="utf-8"))
    assert meta["mode"] == "auto"


@pytest.mark.xfail(strict=True, reason="AUDIT-DIV-04: BuildDivergentBackbone.translate_gapped emits 'X' on '?'/IUPAC codons instead of '-'")
def test_build_divergent_backbone_treats_ambiguous_question_mark_as_gap():
    bdb = _load_dm("BuildDivergentBackbone")
    assert bdb.translate_gapped("???") == "-"


@pytest.mark.xfail(strict=True, reason="AUDIT-DIV-05: SummariseSegmentBackbones writes status='ok' when downstream step failed (backbone_built=N)")
def test_summarise_segment_backbones_flags_downstream_failure(tmp_path):
    ssb = _load_dm("SummariseSegmentBackbones")
    status = tmp_path / "status_M.tsv"
    with open(status, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["segment", "stage", "status", "detail", "n_genes_scored", "best_gene", "best_coverage", "chosen_genes"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerow({
            "segment": "M", "stage": "selection", "status": "ok", "detail": "",
            "n_genes_scored": "1", "best_gene": "G1", "best_coverage": "1.0", "chosen_genes": "G1",
        })
    out = tmp_path / "summary.tsv"
    orig_argv = sys.argv[:]
    try:
        sys.argv = [
            "SummariseSegmentBackbones.py",
            "--status", str(status),
            "--segments", "M",
            "--built", "",
            "--allow-no-backbone",
            "-o", str(out),
        ]
        ssb.main()
    finally:
        sys.argv = orig_argv
    rows = list(csv.DictReader(open(out, encoding="utf-8"), delimiter="\t"))
    assert rows[0]["backbone_built"] == "N" and rows[0]["status"] != "ok"


@pytest.mark.xfail(strict=True, reason="AUDIT-DIV-06: GraftDenseClusters divides UShER substitution counts by padded supermatrix width instead of rep ungapped sites")
def test_graft_dense_clusters_rescales_by_ungapped_rep_sites(tmp_path, monkeypatch):
    gdc = _load_dm("GraftDenseClusters")
    # Supermatrix width = 1000, but rep R1 only has 100 ungapped bases (900 gaps)
    aln = tmp_path / "aln.fasta"
    seq_r = ("A" * 100) + ("-" * 900)
    seq_q = "G" + ("A" * 99) + ("-" * 900)
    aln.write_text(f">R1\n{seq_r}\n>Q1\n{seq_q}\n>R2\n{seq_r}\n", encoding="utf-8")
    bb = tmp_path / "bb.treefile"
    bb.write_text("(R1:0.1,R2:0.1);", encoding="utf-8")
    assign = tmp_path / "assign.tsv"
    assign.write_text("cluster\taccession\trole\nC1\tR1\trepresentative\nC1\tQ1\tmember\nC2\tR2\trepresentative\n", encoding="utf-8")

    sub_tree = tmp_path / "sub.nh"
    sub_tree.write_text("(R1:0.0,Q1:1.0);", encoding="utf-8")  # 1 substitution over 100 ungapped sites -> 0.01 subs/site
    monkeypatch.setattr(gdc, "place_cluster", lambda *a, **kw: (str(sub_tree), 1))

    outdir = tmp_path / "graft_out"
    orig_argv = sys.argv[:]
    try:
        sys.argv = [
            "GraftDenseClusters.py",
            "--backbone-tree", str(bb),
            "-a", str(aln),
            "--assignment", str(assign),
            "-o", str(outdir),
        ]
        gdc.main()
    finally:
        sys.argv = orig_argv
    from Bio import Phylo
    tree = Phylo.read(outdir / "grafted.treefile", "newick")
    q1 = next(t for t in tree.get_terminals() if t.name == "Q1")
    # Current buggy code gives 1.0 / 1000 = 0.001 instead of 1.0 / 100 = 0.01
    assert abs(q1.branch_length - 0.01) < 1e-6
