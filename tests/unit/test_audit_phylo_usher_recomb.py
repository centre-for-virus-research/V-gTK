"""Carpet-bomb audit tests for phylogenetics, UShER, rooting, clades, recombination, and orphans.

Adds passing coverage for ComputeRobinsonFoulds and MMseqsClustering, plus strict
xfail reproducers (AUDIT-PHY-*).
"""

import argparse
from io import StringIO
import os
import sqlite3
from pathlib import Path

from Bio import Phylo
import pytest

import AddReferencesToDb as artd_mod
import BuildTreeManifest as btm_mod
import CladeAssignment as ca_mod
import ComputeRobinsonFoulds as crf_mod
import Curator as cur_mod
import MMseqsClustering as mm_mod
import PrepareUsherResume as pur_mod
import ScreenReferenceRecombination as srr_mod
import TextFileHandler as tfh_mod
import TreeReRoot as trr_mod
import UsherPlacement as up_mod
import clade_from_tree as cft_mod


# ---------------------------------------------------------------------------
# Passing coverage lifts: ComputeRobinsonFoulds & MMseqsClustering
# ---------------------------------------------------------------------------

def test_compute_robinson_foulds_basic_cases(tmp_path):
    t1 = tmp_path / "t1.nwk"
    t2 = tmp_path / "t2.nwk"
    t3 = tmp_path / "t3.nwk"
    t1.write_text("((A,B),(C,(D,E)));", encoding="utf-8")
    t2.write_text("((A,B),(C,(D,E)));", encoding="utf-8")
    t3.write_text("((A,C),(B,(D,E)));", encoding="utf-8")

    res_same = crf_mod.compute_rf(crf_mod.load_tree(str(t1)), crf_mod.load_tree(str(t2)))
    assert res_same["rf"] == 0 and res_same["normalized_rf"] == 0.0

    res_diff = crf_mod.compute_rf(crf_mod.load_tree(str(t1)), crf_mod.load_tree(str(t3)))
    assert res_diff["rf"] > 0 and res_diff["normalized_rf"] > 0.0

    empty = tmp_path / "empty.nwk"
    empty.write_text("   ", encoding="utf-8")
    with pytest.raises(ValueError):
        crf_mod.load_tree(str(empty))


def test_mmseqs_clustering_completeness_and_sorting_helpers(tmp_path):
    aln = tmp_path / "aln.fasta"
    aln.write_text(
        ">M1\n" + ("A" * 100) + "\n"
        ">Q_SHORT\n" + ("A" * 20) + ("-" * 80) + "\n",
        encoding="utf-8",
    )
    complete = tmp_path / "complete.fa"
    rem = tmp_path / "rem.fa"
    n_c, n_r = mm_mod.split_by_completeness(str(aln), 0.9, str(complete), str(rem), reference_length=100)
    assert n_c == 1 and n_r == 1

    sorted_fa = tmp_path / "sorted.fa"
    mm_mod.sort_fasta_by_informative_length(str(aln), str(sorted_fa))
    headers = [l[1:].strip() for l in sorted_fa.read_text(encoding="utf-8").splitlines() if l.startswith(">")]
    assert headers == ["M1", "Q_SHORT"]


# ---------------------------------------------------------------------------
# Strict xfail reproducers for Phylo, UShER, Rooting, Clades & Recombination
# ---------------------------------------------------------------------------

@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-01: ScreenReferenceRecombination.build_mat passes relative 'reference.pb' to usher instead of outdir/reference.pb")
def test_screen_reference_recombination_writes_mat_inside_outdir(tmp_path, monkeypatch):
    outdir = tmp_path / "recomb_out"
    outdir.mkdir()
    fa = tmp_path / "refs.fa"
    fa.write_text(">R1\nATGC\n", encoding="utf-8")
    nwk = tmp_path / "tree.nwk"
    nwk.write_text("(R1:0.1);", encoding="utf-8")

    captured_cmds = []

    def fake_run(cmd, **kw):
        captured_cmds.append(cmd)
        return argparse.Namespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(srr_mod, "_require_binary", lambda *_: None)
    monkeypatch.setattr(srr_mod.subprocess, "run", fake_run)

    scr = srr_mod.ReferenceRecombinationScreen(outdir=str(outdir), skip_build=False)
    scr.build_mat(str(fa), str(nwk))
    usher_cmd = captured_cmds[1]
    out_idx = usher_cmd.index("-o") + 1
    assert os.path.isabs(usher_cmd[out_idx]) and str(outdir) in usher_cmd[out_idx]


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-02: CladeAssignment EPA-ng fallback crashes with TypeError when gb_matrix=None (how vgtk-init.nf calls it)")
def test_clade_assignment_alignment_works_without_gb_matrix(tmp_path):
    padded = tmp_path / "padded.fa"
    padded.write_text(">R1\nATGC\n>Q1\nATGC\n", encoding="utf-8")
    ref_tsv = tmp_path / "refs.txt"
    ref_tsv.write_text("R1\treference\t1\tG1\t1a\n", encoding="utf-8")
    ca = ca_mod.CladeAssignment(
        major_clade=None,
        minor_clade=None,
        padded_alignment=str(padded),
        base_dir=str(tmp_path),
        output_dir="ca_out",
        gb_matrix=None,
        reference_tsv=str(ref_tsv),
    )
    ref_fa, q_fa = ca.alignment()
    assert os.path.isfile(ref_fa) and os.path.isfile(q_fa)


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-03: ComputeRobinsonFoulds gives normalized_rf=1.0 on identical unrooted 4-taxon topologies")
def test_compute_robinson_foulds_even_taxa_rooting_invariant():
    t_a = Phylo.read(StringIO("(A,(B,(C,D)));"), "newick")
    t_b = Phylo.read(StringIO("(C,(D,(A,B)));"), "newick")
    res = crf_mod.compute_rf(t_a, t_b)
    assert res["rf"] == 0 and res["normalized_rf"] == 0.0


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-04: clade_from_tree swallows RecursionError and returns {} on deep caterpillar tree (depth 1200)")
def test_clade_from_tree_handles_deep_caterpillar_tree():
    # Build depth-1200 caterpillar tree
    nwk = "R1:0.01"
    for i in range(1, 1200):
        nwk = f"({nwk},Q{i}:0.01)"
    nwk += ";"
    calls = cft_mod.assign_labels_from_tree(nwk, {"R1": "G1"}, {"R1": "1a"})
    assert calls.get("Q1") == ("G1", "1a")


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-05: TreeReRoot.find_outgroup rejects monophyletic {O1, O2} when input root lies inside the outgroup")
def test_tree_reroot_handles_root_inside_monophyletic_outgroup():
    tree = Phylo.read(StringIO("(O1:0.1,(O2:0.1,(A:0.2,B:0.2,C:0.2):0.5):0.1);"), "newick")
    rerooter = trr_mod.TreeReRooter(tree=tree, outgroup=["O1", "O2"])
    clade = rerooter.find_outgroup()
    assert {t.name for t in clade.get_terminals()} == {"O1", "O2"}


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-06: PrepareUsherResume omits collapsed identical sequences recorded in identical_sequence_groups.tsv")
def test_prepare_usher_resume_retains_collapsed_identical_sequences(tmp_path):
    run_dir = tmp_path / "usher_run"
    c1 = run_dir / "chunk_0001"
    c2 = run_dir / "chunk_0002"
    c1.mkdir(parents=True)
    c2.mkdir(parents=True)
    (c1 / "final-tree.nh").write_text("(M1:0.1,Q1:0.1);", encoding="utf-8")
    (c2 / "chunk_0002.fasta").write_text(">Q2\nATGC\n", encoding="utf-8")
    (run_dir / "identical_sequence_groups.tsv").write_text(
        "representative\tmembers\nQ1\tQ1,Q1_DUP\n", encoding="utf-8"
    )
    padded = tmp_path / "padded.fa"
    padded.write_text(">M1\nATGC\n>Q1\nATGC\n>Q1_DUP\nATGC\n>Q2\nATGC\n", encoding="utf-8")
    resume_dir = tmp_path / "resume"
    pur_mod.prepare_resume(str(run_dir), str(padded), str(resume_dir))
    resume_ids = pur_mod.read_ids_from_fasta(str(resume_dir / "resume_alignment.fasta"))
    assert "Q1_DUP" in resume_ids


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-07: UsherPlacement._write_newick truncates branch lengths < 5e-6 to 0.00000")
def test_usher_placement_preserves_small_branch_lengths(tmp_path):
    tree_path = tmp_path / "in.nh"
    tree_path.write_text("(M1:1e-6,Q1:2e-6);", encoding="utf-8")
    out_path = tmp_path / "out.nh"
    tree = up_mod.UsherPlacement._read_newick(str(tree_path))
    up_mod.UsherPlacement._write_newick(tree, str(out_path))
    txt = out_path.read_text(encoding="utf-8")
    assert "0.00000" not in txt


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-08: MMseqsClustering.split_by_completeness can drop a master with >10% Ns into remainder")
def test_mmseqs_split_by_completeness_never_drops_master_or_reference(tmp_path):
    aln = tmp_path / "aln.fa"
    aln.write_text(
        ">MASTER\n" + ("A" * 80) + ("N" * 20) + "\n"
        ">Q1\n" + ("A" * 100) + "\n",
        encoding="utf-8",
    )
    comp = tmp_path / "comp.fa"
    rem = tmp_path / "rem.fa"
    mm_mod.split_by_completeness(
        str(aln), 0.9, str(comp), str(rem), reference_length=100, reference_ids={"MASTER"}
    )
    assert "MASTER" in comp.read_text(encoding="utf-8")


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-09: BuildTreeManifest falls back to chunk_0001/final-tree.nh on an aborted chunked run")
def test_build_tree_manifest_ignores_intermediate_chunk_trees(tmp_path):
    usher_root = tmp_path / "usher_inputs"
    seg1 = usher_root / "Usher_refset_1"
    (seg1 / "chunk_0001").mkdir(parents=True)
    (seg1 / "chunk_0001" / "final-tree.nh").write_text("(M1:1,Q1:1);", encoding="utf-8")
    rows = btm_mod.collect_tree_rows(str(usher_root), "usher")
    assert rows == []


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-10: AddReferencesToDb.verify_unchanged compares SQLite internal sqlite_sequence table")
def test_add_references_to_db_verify_unchanged_ignores_sqlite_sequence(tmp_path):
    db1 = tmp_path / "db1.db"
    db2 = tmp_path / "db2.db"
    for idx, p in enumerate((db1, db2), 1):
        conn = sqlite3.connect(str(p))
        conn.execute("CREATE TABLE items (id INTEGER PRIMARY KEY AUTOINCREMENT, val TEXT)")
        conn.execute("INSERT INTO items (val) VALUES ('a')")
        if idx == 2:
            # Simulate an insert+delete that advances sqlite_sequence while items rows stay identical
            conn.execute("INSERT INTO items (val) VALUES ('temp')")
            conn.execute("DELETE FROM items WHERE val='temp'")
        conn.commit()
        conn.close()
    problems = artd_mod.verify_unchanged(str(db1), str(db2))
    assert problems == []


@pytest.mark.xfail(strict=True, reason="AUDIT-PHY-11: TextFileHandler raises NameError and Curator treats Influenza 'NA' as blankish")
def test_orphan_scripts_text_file_handler_and_curator_na(tmp_path):
    f = tmp_path / "t.tsv"
    f.write_text("a\tb\n1\t2\n", encoding="utf-8")
    loader = tfh_mod.TextFileLoader(str(f), "\t")
    assert not loader.load_gb_matrix().empty
    assert not cur_mod.BlastAlignment._is_blankish("NA")
