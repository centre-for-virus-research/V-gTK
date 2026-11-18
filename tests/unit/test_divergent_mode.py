"""Unit tests for the opt-in divergent run mode.

Two things are being protected here:

* the invariants the mode exists to provide - a codon alignment that cannot
  contain a frameshift, and an ORF resolver that never emits a partial codon;
* the guarantee that turning the mode off leaves the legacy path alone.
"""
import importlib.util
import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DM_DIR = os.path.join(REPO_ROOT, "scripts", "divergent_mode")


def _load(name):
    path = os.path.join(DM_DIR, f"{name}.py")
    spec = importlib.util.spec_from_file_location(f"dm_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


Bio = pytest.importorskip("Bio")
extract = _load("ExtractDivergentCds")
codon = _load("CodonAwareAlign")
partition = _load("PartitionByDivergence")

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402


# --------------------------------------------------------------------------
# ORF resolution
# --------------------------------------------------------------------------

# resolve_orf rejects anything under 90 nt as too short to be a gene, so the
# fixtures below use a 30-codon body.
BODY = "AAACCCGGGTTT"[:12] * 10          # 120 nt, no in-frame stop
ORF = "ATG" + BODY + "TAA"               # 126 nt


def test_resolve_orf_finds_complete_orf_on_plus_strand():
    seq = Seq("CCCC" + ORF + "CCCC")
    got = extract.resolve_orf(seq, 4, 4 + len(ORF), 1)
    assert got is not None
    start, end, strand, flag = got
    assert (strand, flag) == (1, "complete")
    assert str(seq[start:end]).startswith("ATG")
    assert str(seq[start:end]).endswith("TAA")
    assert (end - start) % 3 == 0


def test_resolve_orf_handles_reverse_strand():
    seq = Seq("CCCC") + Seq(ORF).reverse_complement() + Seq("CCCC")
    got = extract.resolve_orf(seq, 4, 4 + len(ORF), -1)
    assert got is not None
    start, end, strand, flag = got
    assert strand == -1
    recovered = str(seq[start:end].reverse_complement())
    assert recovered.startswith("ATG") and recovered.endswith("TAA")


def test_resolve_orf_flags_three_prime_truncation_and_keeps_frame():
    """A record ending inside the gene has no stop; keep it, flagged, in frame."""
    seq = Seq("ATG" + "AAACCCGGGTTT" * 8 + "AC")  # deliberately not a multiple of 3
    got = extract.resolve_orf(seq, 0, 60, 1, allow_partial=True)
    assert got is not None
    start, end, _, flag = got
    assert "partial3" in flag
    assert (end - start) % 3 == 0, "a partial ORF must still be whole codons"


def test_resolve_orf_rejects_truncation_when_complete_required():
    seq = Seq("ATG" + "AAACCCGGGTTT" * 8 + "AC")
    assert extract.resolve_orf(seq, 0, 60, 1, allow_partial=False) is None


def test_resolve_orf_rejects_a_hit_too_short_to_be_a_gene():
    seq = Seq("ATG" + "AAACCC" + "TAA")
    assert extract.resolve_orf(seq, 0, 9, 1) is None


# --------------------------------------------------------------------------
# Codon back-projection - the frameshift guarantee
# --------------------------------------------------------------------------

def test_back_project_emits_three_gap_characters_per_protein_gap():
    nt = [SeqRecord(Seq("ATGAAACCC"), id="a"), SeqRecord(Seq("ATGCCC"), id="b")]
    prot = [SeqRecord(Seq("MKP"), id="a"), SeqRecord(Seq("M-P"), id="b")]
    out = codon.back_project(nt, prot)
    by_id = {r.id: str(r.seq) for r in out}
    assert by_id["a"] == "ATGAAACCC"
    assert by_id["b"] == "ATG---CCC"
    assert len({len(s) for s in by_id.values()}) == 1


def test_back_project_output_never_contains_a_frameshift_gap():
    nt = [
        SeqRecord(Seq("ATGAAACCCGGG"), id="a"),
        SeqRecord(Seq("ATGGGG"), id="b"),
        SeqRecord(Seq("ATGAAAGGG"), id="c"),
    ]
    prot = [
        SeqRecord(Seq("MKPG"), id="a"),
        SeqRecord(Seq("M--G"), id="b"),
        SeqRecord(Seq("MK-G"), id="c"),
    ]
    out = codon.back_project(nt, prot)
    assert codon.check_codon_integrity(out) == []


def test_check_codon_integrity_catches_a_planted_frameshift():
    bad = [
        SeqRecord(Seq("ATGAAACCC"), id="a"),
        SeqRecord(Seq("AT-AAACCC"), id="b"),  # single-base gap, out of phase
    ]
    problems = codon.check_codon_integrity(bad)
    assert any("frameshift" in p[1] for p in problems), problems


def test_back_project_refuses_a_protein_with_no_nucleotide_record():
    with pytest.raises(KeyError):
        codon.back_project(
            [SeqRecord(Seq("ATG"), id="a")],
            [SeqRecord(Seq("M"), id="missing")],
        )


# --------------------------------------------------------------------------
# Distance partitioning
# --------------------------------------------------------------------------

def test_p_distance_ignores_columns_gapped_in_either_sequence():
    p, compared = partition.p_distance("ACGT", "ACG-")
    assert compared == 3
    assert p == 0.0


def test_p_distance_counts_real_differences():
    p, compared = partition.p_distance("AAAA", "AAAT")
    assert compared == 4
    assert p == pytest.approx(0.25)


def test_p_distance_on_no_overlap_reports_maximum():
    p, compared = partition.p_distance("----", "ACGT")
    assert compared == 0
    assert p == 1.0


# --------------------------------------------------------------------------
# Gene selection is per segment, and a skipped segment is never silent
#
# The defect these cover: selection used to be scoped to one seed segment's
# cluster representatives while scoring genes from every segment. Genes on the
# other segments share no accession with that partition, so they scored zero
# coverage and dropped out - a family where all three segments carry a core
# gene got one backbone and no word about the other two.
# --------------------------------------------------------------------------

SELECT = os.path.join(DM_DIR, "SelectBackboneGenes.py")
SUMMARISE = os.path.join(DM_DIR, "SummariseSegmentBackbones.py")


def _selection_fixture(tmp_path):
    """Two segments, each with its own accessions, genes and partition."""
    codon = tmp_path / "codon"
    codon.mkdir()
    accs = {"L": [f"L{i}" for i in range(4)], "S": [f"S{i}" for i in range(4)]}
    genes = {"L_RdRp": ("L", accs["L"], 600), "S_N": ("S", accs["S"], 300)}
    for gene, (seg, members, width) in genes.items():
        with open(codon / f"{gene}.codon.fasta", "w") as fh:
            for a in members:
                fh.write(f">{a}\n{'ACG' * (width // 3)}\n")

    report = tmp_path / "report.tsv"
    with open(report, "w") as fh:
        fh.write("segment\tgroup\taccession\tstatus\n")
        for gene, (seg, members, _) in genes.items():
            for a in members:
                fh.write(f"{seg}\t{gene}\t{a}\tok\n")

    for seg, members in accs.items():
        with open(tmp_path / f"assign_{seg}.tsv", "w") as fh:
            fh.write("cluster\taccession\trole\trepresentative\n")
            for i, a in enumerate(members):
                role = "representative" if i % 2 == 0 else "dense_member"
                fh.write(f"C001\t{a}\t{role}\t{members[0]}\n")
    return codon, report


def _run_select(tmp_path, segment, assignment, *extra):
    codon, report = _selection_fixture(tmp_path)
    out = tmp_path / f"{segment}_genes.txt"
    status = tmp_path / f"{segment}_status.tsv"
    proc = subprocess.run(
        [sys.executable, "-I", SELECT, "--codon-dir", str(codon),
         "--report", str(report), "-o", str(out), "--status", str(status),
         "--segments", segment, "--assignment", str(tmp_path / assignment),
         *extra],
        capture_output=True, text=True, timeout=300)
    return proc, out, status


def test_selection_scoped_to_a_segment_picks_that_segments_gene(tmp_path):
    proc, out, _ = _run_select(tmp_path, "S", "assign_S.tsv")
    assert proc.returncode == 0, proc.stderr[-800:]
    assert out.read_text().strip() == "S_N"


def test_selection_refuses_a_partition_built_on_another_segment(tmp_path):
    """The original defect, made unrepresentable.

    Scoring segment S against segment L's clusters silently returns "no
    usable gene on S" when nothing about S was ever measured.
    """
    proc, out, _ = _run_select(tmp_path, "S", "assign_L.tsv", "--allow-empty")
    assert proc.returncode != 0, "mis-scoped selection must abort, not skip"
    assert "no accession from segment S" in (proc.stderr + proc.stdout)
    assert not out.exists() or not out.read_text().strip()


def test_selection_refuses_several_segments_against_one_partition(tmp_path):
    proc, _, _ = _run_select(tmp_path, "L,S", "assign_L.tsv")
    assert proc.returncode != 0
    assert "exactly one" in (proc.stderr + proc.stdout)


def test_segment_with_no_qualifying_gene_skips_instead_of_killing_the_run(tmp_path):
    """One unusable segment must not take down the segments that work."""
    proc, out, status = _run_select(
        tmp_path, "S", "assign_S.tsv", "--allow-empty", "--min-coverage", "1.1")
    assert proc.returncode == 0, proc.stderr[-800:]
    assert out.read_text().strip() == "", "a skipped segment emits no genes"
    body = status.read_text()
    assert "below_min_coverage" in body, body


def test_skipping_is_fatal_without_allow_empty(tmp_path):
    proc, _, _ = _run_select(tmp_path, "S", "assign_S.tsv",
                             "--min-coverage", "1.1")
    assert proc.returncode != 0


def _coverage_fixture(tmp_path):
    """One segment: a dense cluster of five, plus four lone taxa.

    `M_CORE` is in every cluster, `M_SPIKE` only in the dense cluster and one
    lone taxon - the shape of Cystoviridae's host-attachment gene, which is
    universal inside phi6 and absent from almost everything else.
    """
    dense = [f"M{i}" for i in range(5)]
    lone = [f"M{i}" for i in range(5, 9)]
    codon = tmp_path / "codon"
    codon.mkdir()
    members = {"M_CORE": dense + lone, "M_SPIKE": dense + [lone[0]]}
    for gene, accs in members.items():
        with open(codon / f"{gene}.codon.fasta", "w") as fh:
            for a in accs:
                fh.write(f">{a}\n{'ACG' * 100}\n")

    report = tmp_path / "report.tsv"
    with open(report, "w") as fh:
        fh.write("segment\tgroup\taccession\tstatus\n")
        for gene, accs in members.items():
            for a in accs:
                fh.write(f"M\t{gene}\t{a}\tok\n")

    manifest = tmp_path / "segment_manifest.tsv"
    with open(manifest, "w") as fh:
        fh.write("accession\tsegment\tlength\tkept\n")
        for a in dense + lone:
            fh.write(f"{a}\tM\t4000\tY\n")

    # The partition is built on M_SPIKE's alignment, so it only knows the
    # taxa that carry M_SPIKE - exactly the circularity being guarded against.
    assignment = tmp_path / "assign_M.tsv"
    with open(assignment, "w") as fh:
        fh.write("cluster\taccession\trole\trepresentative\n")
        for a in dense:
            role = "representative" if a == dense[0] else "dense_member"
            fh.write(f"C001\t{a}\t{role}\t{dense[0]}\n")
        fh.write(f"C002\t{lone[0]}\trepresentative\t{lone[0]}\n")
    return codon, report, manifest, assignment


def _coverage(tmp_path, use_manifest=True):
    codon, report, manifest, assignment = _coverage_fixture(tmp_path)
    out = tmp_path / "genes.txt"
    summary = tmp_path / "selection.tsv"
    cmd = [sys.executable, "-I", SELECT, "--codon-dir", str(codon),
           "--report", str(report), "-o", str(out), "--summary", str(summary),
           "--segments", "M", "--assignment", str(assignment),
           "--min-coverage", "0.9", "--allow-empty"]
    if use_manifest:
        cmd += ["--manifest", str(manifest)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-800:]
    rows = [dict(zip(summary.read_text().splitlines()[0].split("\t"),
                     line.split("\t")))
            for line in summary.read_text().splitlines()[1:]]
    return out, {r["gene"]: r for r in rows}


def test_coverage_counts_a_dense_cluster_once(tmp_path):
    """Five near-identical isolates are one unit of evidence, not five."""
    _, rows = _coverage(tmp_path)
    assert rows["M_CORE"]["n_taxa"] == "5", "9 genomes collapse to 5 clusters"
    assert rows["M_CORE"]["n_genomes"] == "9", "raw genome count is kept too"
    assert float(rows["M_CORE"]["coverage"]) == pytest.approx(1.0)


def test_coverage_denominator_is_not_the_partitions_own_membership(tmp_path):
    """The defect that survived the first per-segment fix.

    A partition is built from one gene's alignment, so taxa missing that gene
    are missing from the partition. Scoring against partition membership
    measures the seed gene over the taxa that have it and returns 1.00 every
    time - which turned Cystoviridae's M_P3 (10% of divergent taxa) into a
    backbone gene at full coverage.
    """
    genes, rows = _coverage(tmp_path)
    assert float(rows["M_SPIKE"]["coverage"]) == pytest.approx(0.4), rows
    assert genes.read_text().strip() == "M_CORE", \
        "only the gene present across the divergent taxa may carry a backbone"


def test_coverage_without_a_manifest_falls_back_to_the_report(tmp_path):
    _, rows = _coverage(tmp_path, use_manifest=False)
    assert rows["M_CORE"]["n_genomes"] == "9"


def test_a_skipped_segment_still_publishes_its_gene_ranking(tmp_path):
    """The ranking is the evidence for the skip; it must not vanish with it."""
    codon, report, manifest, assignment = _coverage_fixture(tmp_path)
    summary = tmp_path / "selection.tsv"
    proc = subprocess.run(
        [sys.executable, "-I", SELECT, "--codon-dir", str(codon),
         "--report", str(report), "-o", str(tmp_path / "genes.txt"),
         "--summary", str(summary), "--manifest", str(manifest),
         "--segments", "M", "--assignment", str(assignment),
         "--min-coverage", "1.1", "--allow-empty"],
        capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-800:]
    body = summary.read_text().splitlines()
    assert len(body) > 1, "a skipped segment must still rank its genes"
    assert all(line.endswith("\tN") for line in body[1:]), \
        "nothing may be marked chosen on a skipped segment"


def _write_status(path, segment, status, stage="selection", genes=""):
    with open(path, "w") as fh:
        fh.write("segment\tstage\tstatus\tdetail\tn_genes_scored\tbest_gene\t"
                 "best_coverage\tchosen_genes\n")
        fh.write(f"{segment}\t{stage}\t{status}\t\t1\tX\t0.5\t{genes}\n")


def test_summary_names_every_requested_segment_including_the_skipped(tmp_path):
    ok = tmp_path / "L.tsv"
    skipped = tmp_path / "M.tsv"
    _write_status(ok, "L", "ok", genes="L_RdRp")
    _write_status(skipped, "M", "below_min_coverage")
    out = tmp_path / "backbone_summary.tsv"
    proc = subprocess.run(
        [sys.executable, "-I", SUMMARISE, "--status", str(ok), str(skipped),
         "--segments", "L,M,S", "--built", "L", "-o", str(out)],
        capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr[-800:]
    rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
    assert [r[0] for r in rows] == ["L", "M", "S"], \
        "every requested segment needs a row, built or not"
    built = {r[0]: r[1] for r in rows}
    assert built == {"L": "Y", "M": "N", "S": "N"}
    # A segment that never reached selection is still reported.
    assert any("not_run" in r for r in rows[2])


def test_summary_fails_when_no_segment_produced_a_backbone(tmp_path):
    skipped = tmp_path / "M.tsv"
    _write_status(skipped, "M", "below_min_coverage")
    out = tmp_path / "backbone_summary.tsv"
    proc = subprocess.run(
        [sys.executable, "-I", SUMMARISE, "--status", str(skipped),
         "--segments", "M", "--built", "", "-o", str(out)],
        capture_output=True, text=True, timeout=120)
    assert proc.returncode != 0, "a run with no backbone at all has no result"


def test_workflow_fans_out_over_segments():
    """The subworkflow must drive selection from a channel of segments."""
    sub = _read("subworkflows/local/divergent_mode/divergent_mode.nf")
    assert "Channel.fromList" in sub and "divergent_segments" in sub, \
        "gene selection must run once per segment"
    assert "DM_SEGMENT_SUMMARY(" in sub, \
        "a run must report what each segment did"
    # The shared directories are consumed once per segment, so they have to be
    # value channels; a queue channel is drained by the first segment.
    assert ".first()" in sub


def test_per_segment_outputs_do_not_overwrite_each_other():
    mod = _read("modules/local/divergent_mode/main.nf")
    for proc in ("DM_CONCATENATE", "DM_BACKBONE_ML", "DM_GRAFT",
                 "DM_SELECT_BACKBONE_GENES"):
        body = mod[mod.index(f"process {proc} "):]
        body = body[:body.index("shell:")]
        assert "by_segment/${segment}" in body, \
            f"{proc} must publish under its own segment directory"


def test_split_segments_is_not_hard_coded_to_lms():
    """Segment names are a parameter; the mode has to generalise past L/M/S."""
    src = open(os.path.join(DM_DIR, "split_segments.py"), encoding="utf-8").read()
    assert '"--segments"' in src
    assert "segment_regex" in src


# --------------------------------------------------------------------------
# Non-regression: the flag is off, and the legacy path is untouched
# --------------------------------------------------------------------------

def _read(path):
    with open(os.path.join(REPO_ROOT, path), encoding="utf-8") as fh:
        return fh.read()


def test_divergent_mode_defaults_to_off():
    cfg = _read("nextflow.config")
    assert "run_divergent_mode        = false" in cfg, \
        "divergent mode must default to false"


def test_every_divergent_param_is_allowed_by_the_param_guard():
    """Any param in the config that the guard does not know aborts the run."""
    import re
    cfg = _read("nextflow.config")
    nf = _read("vgtk-init.nf")
    declared = set(re.findall(r"^\s*(divergent_\w+|run_divergent_mode)\s*=", cfg, re.M))
    assert declared, "no divergent params found in nextflow.config"
    for name in declared:
        assert f'"{name}"' in nf, f"{name} missing from scriptDefinedParams"


def test_divergent_branch_returns_before_any_legacy_process():
    """The opt-in branch must short-circuit, not run alongside the legacy path."""
    nf = _read("vgtk-init.nf")
    start = nf.index("workflow {")
    branch = nf.index("if( params.run_divergent_mode", start)
    legacy = nf.index("TEST_DEPENDENCIES()", start)
    assert branch < legacy, "divergent branch must come before the legacy path"
    body = nf[branch:legacy]
    assert "DIVERGENT_MODE(" in body
    assert "return" in body, "the divergent branch must return"


def test_divergent_modules_live_in_their_own_namespace():
    for path in ("modules/local/divergent_mode/main.nf",
                 "subworkflows/local/divergent_mode/divergent_mode.nf"):
        assert os.path.exists(os.path.join(REPO_ROOT, path)), f"missing {path}"


def test_divergent_processes_do_not_collide_with_legacy_process_names():
    import re
    legacy = set(re.findall(r"^process\s+(\w+)", _read("vgtk-init.nf"), re.M))
    new = set(re.findall(r"^process\s+(\w+)",
                         _read("modules/local/divergent_mode/main.nf"), re.M))
    assert new, "no processes found in the divergent module"
    assert not (legacy & new), f"name collision: {sorted(legacy & new)}"


def test_extraction_defaults_to_a_single_tblastn_pass():
    """Query enrichment lost genes on Cystoviridae; it must stay off by default."""
    cfg = _read("nextflow.config")
    assert "divergent_extract_iterations = 1" in cfg


@pytest.mark.parametrize("script", [
    "BuildProteinProfiles", "ExtractDivergentCds", "CodonAwareAlign",
    "ConcatenateCodonAlignments", "PartitionByDivergence",
    "BuildDivergentBackbone", "GraftDenseClusters", "SelectBackboneGenes",
    "SummariseSegmentBackbones", "split_segments",
])
def test_every_divergent_script_compiles_and_offers_help(script):
    path = os.path.join(DM_DIR, f"{script}.py")
    assert os.path.exists(path), f"missing {path}"
    proc = subprocess.run([sys.executable, "-I", path, "--help"],
                          capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr[-800:]


# --------------------------------------------------------------------------
# Grafting: the unified tree must carry one branch-length unit
# --------------------------------------------------------------------------

def test_graft_rescales_parsimony_branch_lengths_by_default():
    """UShER counts substitutions; the ML backbone measures substitutions/site.

    Grafting without converting yields one tree in two units, where a
    within-cluster branch of 119 substitutions dwarfs a between-genus branch of
    1.2 subs/site and every downstream distance silently reads the mixture as
    comparable.
    """
    src = open(os.path.join(DM_DIR, "GraftDenseClusters.py"), encoding="utf-8").read()
    assert "rescale_subtree_branches" in src
    assert "--no-rescale-subtree-branches" in src, "the conversion must be opt-out"
    assert "aln_width" in src


def test_grafted_tree_branch_lengths_are_on_one_scale():
    """No branch may look like a raw substitution count."""
    pytest.importorskip("Bio.Phylo")
    from Bio import Phylo
    import glob as _glob
    trees = sorted(_glob.glob(os.path.join(
        REPO_ROOT, "test_out", "divergent_test", "divergent_mode",
        "by_segment", "*", "grafted", "grafted.treefile")))
    if not trees:
        pytest.skip("run the divergent_test profile first")
    lengths = []
    for tree_path in trees:
        tree = Phylo.read(tree_path, "newick")
        lengths += [c.branch_length for c in tree.find_clades()
                    if c.branch_length is not None]
    assert lengths, "tree has no branch lengths"
    # A substitutions/site branch on a few-kb alignment cannot plausibly
    # exceed ~3; a raw parsimony count on a dense cluster runs to the hundreds.
    assert max(lengths) < 3.0, f"branch length {max(lengths)} looks like a raw count"


# --------------------------------------------------------------------------
# The auto cut-off judges a band against the spread, not its own lower edge
#
# The defect these cover: the band test was `gap >= max(0.01, 2*lo)`, i.e.
# `hi >= 3*lo`. Because candidate bands are capped at --auto-ceiling, no band
# whose lower edge exceeded ceiling/3 could ever qualify, however empty it
# was. A species diverse enough to spread its own isolates past that bound was
# forced to the 1.5% fallback and shattered into singletons - which is what
# Puumala orthohantavirus does on segment L (isolates 8-16% apart, then a
# clean empty band to 23.4% before any other species).
# --------------------------------------------------------------------------

PARTITION = os.path.join(DM_DIR, "PartitionByDivergence.py")
_ALN_N = 2000


def _aln_with_distances(tmp_path, muts, name="aln.fasta"):
    """One taxon per entry in `muts`, mutated over disjoint blocks.

    Disjoint blocks make the pairwise p-distance exact: taxa i and j differ at
    every site either one of them mutated, so d(i,j) = (m_i + m_j) / N.
    """
    path = tmp_path / name
    start = 0
    with open(path, "w") as fh:
        for i, m in enumerate(muts):
            seq = ["A"] * _ALN_N
            for p in range(start, start + m):
                seq[p] = "C"
            start += m
            fh.write(f">t{i}\n{''.join(seq)}\n")
    assert start <= _ALN_N, "mutation blocks must stay disjoint"
    return path


def _run_partition(tmp_path, muts, *extra):
    aln = _aln_with_distances(tmp_path, muts)
    out = tmp_path / "part"
    proc = subprocess.run(
        [sys.executable, "-I", PARTITION, "-a", str(aln), "-o", str(out),
         *extra],
        capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-800:]
    import json
    with open(out / "partition_cutoff.json") as fh:
        return json.load(fh), proc


def test_auto_cutoff_accepts_a_clean_band_at_high_cluster_scale(tmp_path):
    """The Hantaviridae segment L shape: a real cluster that is not tight.

    Four taxa 12% apart, then an empty band to 24%. The old rule needed a gap
    of 2*0.12 = 0.24 and rejected it.
    """
    report, _ = _run_partition(tmp_path, [120, 120, 120, 120, 360, 360, 360])
    assert report["mode"] == "auto", report
    assert report["band_lo"] == pytest.approx(0.12)
    assert report["band_hi"] == pytest.approx(0.24)
    assert report["cutoff"] == pytest.approx(0.18)


def test_auto_cutoff_at_high_scale_keeps_the_cluster_together(tmp_path):
    """The consequence that matters: the dense cluster survives as one."""
    _run_partition(tmp_path, [120, 120, 120, 120, 360, 360, 360])
    rows = [l.split("\t") for l in
            (tmp_path / "part" / "partition_assignment.tsv")
            .read_text().splitlines()[1:]]
    sizes = sorted({r[0]: int(r[5]) for r in rows}.values(), reverse=True)
    assert sizes[0] == 4, f"the 12%-apart cluster was split: {sizes}"


def test_auto_cutoff_still_falls_back_on_a_continuum(tmp_path):
    """The guard the band test exists for must survive the fix."""
    report, _ = _run_partition(tmp_path, [40 * i for i in range(1, 10)])
    assert report["mode"] == "auto_fallback", report
    assert report["cutoff"] == 0.015
    assert "of the range" in report["reason"]


def test_auto_cutoff_unchanged_on_a_tight_cluster(tmp_path):
    """The Cystoviridae shape: a very tight cluster, far from everything.

    Both the old and the new rule accept this one, at the same cut-off.
    """
    report, _ = _run_partition(tmp_path, [4, 4, 4, 300, 300, 300])
    assert report["mode"] == "auto"
    assert report["cutoff"] == pytest.approx((0.004 + 0.152) / 2)


def test_auto_cutoff_fallback_records_why(tmp_path):
    """A fallback that does not say why cannot be checked by the reader."""
    report, proc = _run_partition(tmp_path, [40 * i for i in range(1, 10)])
    assert report["reason"]
    assert report["band_frac"] < report["min_band_frac"]
    assert "falling back" in (proc.stdout + proc.stderr)
