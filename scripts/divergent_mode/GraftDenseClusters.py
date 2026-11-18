#!/usr/bin/env python3
"""Place dense clusters by parsimony and graft them onto the ML backbone.

Each cluster's members are placed with UShER against their own representative,
which is the only coordinate system they are close enough to share. The
resulting subtree then replaces that representative's tip in the ML backbone,
so deep splits keep their maximum-likelihood branch lengths while dense
clusters keep the resolution parsimony gives them cheaply.

This is the part that a single-reference parsimony run gets wrong. UShER
assumes the samples it places are a short distance from the tree it is given;
across cystovirus genera that assumption fails, and every placement collapses
towards the root. Restricting each UShER run to one low-divergence cluster
keeps it inside its working range.

faToVcf is called with `-includeNoAltN` for the reason documented at
scripts/UsherPlacement.py: without it, a site with no ALT allele in the batch
is dropped from the VCF, and UShER then reads every sample as carrying the
reference base at that site instead of as missing data.
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys

from Bio import Phylo, SeqIO


def run(cmd, **kw):
    proc = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {proc.stdout[-800:]}\n{proc.stderr[-800:]}")
    return proc


def reject_option_like(**values):
    """faToVcf parses a leading '-' as an option; it has no '--' separator."""
    for name, value in values.items():
        if value is not None and str(value).startswith("-"):
            raise ValueError(f"{name} {value!r} would be read as a faToVcf option")


def place_cluster(rep, members, recs, workdir, threads):
    """UShER-place one cluster's members against its representative."""
    cdir = os.path.join(workdir, rep.replace(".", "_"))
    os.makedirs(cdir, exist_ok=True)

    ordered = [rep] + [m for m in members if m != rep]
    fa = os.path.join(cdir, "cluster.fa")
    SeqIO.write([recs[a] for a in ordered], fa, "fasta")

    seed = os.path.join(cdir, "seed.nwk")
    with open(seed, "w") as fh:
        fh.write(f"({rep});\n")

    vcf = os.path.join(cdir, "cluster.vcf")
    reject_option_like(rep=rep, fa=fa, vcf=vcf)
    run(["faToVcf", "-includeNoAltN", f"-ref={rep}", fa, vcf])

    outdir = os.path.join(cdir, "usher_out")
    run(["usher", "-t", seed, "-v", vcf, "-d", outdir, "-T", str(threads)])

    for name in ("final-tree.nh", "uncondensed-final-tree.nh"):
        path = os.path.join(outdir, name)
        if os.path.isfile(path):
            variable_sites = sum(
                1 for line in open(vcf) if line and not line.startswith("#")
            )
            return path, variable_sites
    raise FileNotFoundError(f"usher wrote no tree in {outdir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone-tree", required=True)
    ap.add_argument("-a", "--alignment", required=True)
    ap.add_argument("--assignment", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--prefix", default="grafted")
    ap.add_argument("--no-rescale-subtree-branches", dest="rescale_subtree_branches",
                    action="store_false",
                    help="leave UShER parsimony branch lengths as substitution "
                         "counts instead of converting to substitutions/site")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    workdir = os.path.join(args.outdir, "work")
    os.makedirs(workdir, exist_ok=True)

    recs = {r.id: r for r in SeqIO.parse(args.alignment, "fasta")}
    aln_width = len(next(iter(recs.values())).seq) if recs else 0
    rows = list(csv.DictReader(open(args.assignment), delimiter="\t"))

    clusters = {}
    for r in rows:
        clusters.setdefault(r["cluster"], {"rep": None, "members": []})
        clusters[r["cluster"]]["members"].append(r["accession"])
        if r["role"] == "representative":
            clusters[r["cluster"]]["rep"] = r["accession"]

    tree = Phylo.read(args.backbone_tree, "newick")
    # Strip IQ-TREE's "/support" suffix so tip lookups match accessions.
    for clade in tree.find_clades():
        if clade.name and "/" in clade.name:
            clade.name = clade.name.split("/")[0]

    report = []
    grafted = 0
    for cname, info in sorted(clusters.items()):
        rep, members = info["rep"], info["members"]
        if rep is None:
            continue
        if len(members) == 1:
            report.append({"cluster": cname, "representative": rep,
                           "n_members": 1, "action": "backbone_tip_only",
                           "variable_sites": "", "subtree_tips": 1})
            continue

        sub_path, variable_sites = place_cluster(rep, members, recs, workdir, args.threads)
        subtree = Phylo.read(sub_path, "newick")
        tips = subtree.count_terminals()

        # UShER reports branch lengths as a COUNT of substitutions; the ML
        # backbone reports them as substitutions per site. Grafting without
        # converting produces one tree carrying two different units, where a
        # within-cluster branch of 119 substitutions dwarfs a between-genus
        # branch of 1.2 subs/site - and any downstream distance, plot or
        # rooting silently reads the mixture as comparable. Divide by the
        # alignment width to put the subtree on the backbone's scale.
        if args.rescale_subtree_branches and aln_width:
            for clade in subtree.find_clades():
                if clade.branch_length:
                    clade.branch_length = clade.branch_length / aln_width

        target = next((c for c in tree.get_terminals() if c.name == rep), None)
        if target is None:
            report.append({"cluster": cname, "representative": rep,
                           "n_members": len(members), "action": "rep_not_in_backbone",
                           "variable_sites": variable_sites, "subtree_tips": tips})
            continue

        # Replace the representative tip with the placed subtree, keeping the
        # backbone branch that leads to it.
        target.clades = subtree.root.clades
        target.name = None
        grafted += 1
        report.append({"cluster": cname, "representative": rep,
                       "n_members": len(members), "action": "grafted",
                       "variable_sites": variable_sites, "subtree_tips": tips})
        print(f"{cname}: grafted {tips} tips at {rep} "
              f"({variable_sites} variable sites)", flush=True)

    out_tree = os.path.join(args.outdir, f"{args.prefix}.treefile")
    Phylo.write(tree, out_tree, "newick")

    rep_path = os.path.join(args.outdir, f"{args.prefix}_graft_report.tsv")
    with open(rep_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(report[0]), delimiter="\t")
        w.writeheader()
        w.writerows(report)

    total_tips = tree.count_terminals()
    print(f"\ngrafted {grafted} cluster(s); unified tree has {total_tips} tips")
    print(f"tree   -> {out_tree}")
    print(f"report -> {rep_path}")
    if total_tips != len(recs):
        print(f"WARNING: alignment has {len(recs)} taxa but tree has {total_tips} tips")


if __name__ == "__main__":
    main()
