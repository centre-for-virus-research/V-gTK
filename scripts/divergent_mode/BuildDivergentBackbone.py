#!/usr/bin/env python3
"""ML backbone with joint ancestral state reconstruction, and PCP extraction.

Runs IQ-TREE over the backbone representatives only. Keeping the dense cluster
members out is the point: they carry almost no extra signal about the deep
splits, and including 48 near-identical phi6 sequences biases branch length
estimation towards that one cluster.

Ancestral sequences come from IQ-TREE's `-asr`, which writes marginal
posterior probabilities per site in a `.state` file. Each branch of the tree
then yields a parent-child pair: the ancestral sequence at the parent node and
the sequence at the child, with the amino-acid substitutions between them
resolved back to a gene and a codon index via the supermatrix coordinate map.
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys

from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_iqtree(aln, prefix, model, threads, mem, extra=None, seed=12345):
    exe = shutil.which("iqtree3") or shutil.which("iqtree2") or shutil.which("iqtree")
    if exe is None:
        sys.exit("no iqtree binary on PATH")
    cmd = [exe, "-s", aln, "--prefix", prefix, "-m", model,
           "-nt", str(threads), "-seed", str(seed), "-redo"]
    if mem:
        cmd += ["-mem", mem]
    if extra:
        cmd += list(extra)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"iqtree failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    return exe


def read_state_file(path):
    """Ancestral sequences keyed by internal node name, from IQ-TREE .state."""
    seqs = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts or parts[0] == "Node":
                continue
            node, site, state = parts[0], int(parts[1]), parts[2]
            seqs.setdefault(node, {})[site] = state
    out = {}
    for node, sites in seqs.items():
        n = max(sites)
        out[node] = "".join(sites.get(i, "-") for i in range(1, n + 1))
    return out


def load_coords(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rows:
        r["nt_start"] = int(r["nt_start"])
        r["nt_end"] = int(r["nt_end"])
    return rows


def locate(coords, nt_pos):
    """Map a 1-based supermatrix nt position to (gene, codon index, codon pos)."""
    for r in coords:
        if r["nt_start"] <= nt_pos <= r["nt_end"]:
            off = nt_pos - r["nt_start"]
            return r["gene"], off // 3 + 1, off % 3 + 1
    return "", 0, 0


def translate_gapped(seq):
    """Codon-wise translation that keeps alignment columns; '-' for gap codons."""
    out = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        if "-" in codon or "N" in codon:
            out.append("-")
        else:
            try:
                out.append(str(Seq(codon).translate(table=11)))
            except Exception:
                out.append("X")
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--alignment", required=True)
    ap.add_argument("--assignment", required=True)
    ap.add_argument("--coords", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--model", default="GTR+I+G4")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--mem", default="")
    ap.add_argument("--prefix", default="backbone")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--reuse", action="store_true",
                    help="reuse an existing IQ-TREE run if its outputs are present")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    reps = [r["accession"] for r in csv.DictReader(open(args.assignment), delimiter="\t")
            if r["role"] == "representative"]
    if len(reps) < 4:
        sys.exit(f"only {len(reps)} backbone representatives; need at least 4")

    recs = {r.id: r for r in SeqIO.parse(args.alignment, "fasta")}
    sub = [recs[a] for a in reps if a in recs]
    sub_path = os.path.join(args.outdir, f"{args.prefix}_representatives.fasta")
    SeqIO.write(sub, sub_path, "fasta")
    print(f"backbone: {len(sub)} representatives x {len(sub[0].seq)} nt")

    prefix = os.path.join(args.outdir, args.prefix)
    extra = ["-asr"]
    if args.bootstrap and len(sub) >= 4:
        extra += ["-B", str(args.bootstrap)]
    if args.reuse and os.path.exists(f"{prefix}.state") and os.path.exists(f"{prefix}.treefile"):
        print("reusing existing IQ-TREE outputs")
    else:
        run_iqtree(sub_path, prefix, args.model, args.threads, args.mem, extra)

    tree_path = f"{prefix}.treefile"
    state_path = f"{prefix}.state"
    if not os.path.exists(state_path):
        sys.exit(f"iqtree produced no .state file at {state_path}")

    ancestors = read_state_file(state_path)
    tree = Phylo.read(tree_path, "newick")
    coords = load_coords(args.coords)

    anc_path = os.path.join(args.outdir, f"{args.prefix}_ancestral.fasta")
    SeqIO.write([SeqRecord(Seq(s), id=n, description="") for n, s in sorted(ancestors.items())],
                anc_path, "fasta")

    tip_seq = {r.id: str(r.seq).upper() for r in sub}
    all_seq = dict(tip_seq)
    all_seq.update({k: v.upper() for k, v in ancestors.items()})

    # IQ-TREE writes internal labels as "Node7/97" - the node name it used in
    # the .state file, then the bootstrap support. Biopython keeps the whole
    # string as the clade name, so matching it against the .state keys fails on
    # every supported node and only the root survives. Split the two apart.
    def node_name(clade):
        if not clade.name:
            return None, None
        name, _, support = clade.name.partition("/")
        return name, (support or None)

    pcp_rows = []
    for clade in tree.find_clades(order="level"):
        parent_name, _ = node_name(clade)
        if not parent_name or parent_name not in all_seq:
            continue
        for child in clade.clades:
            child_name, child_support = node_name(child)
            if not child_name or child_name not in all_seq:
                continue
            p_nt, c_nt = all_seq[parent_name], all_seq[child_name]
            p_aa, c_aa = translate_gapped(p_nt), translate_gapped(c_nt)
            subs = []
            for idx, (x, y) in enumerate(zip(p_aa, c_aa), start=1):
                if x == y or x == "-" or y == "-":
                    continue
                gene, codon, _ = locate(coords, (idx - 1) * 3 + 1)
                subs.append(f"{gene}:{x}{codon}{y}")
            nt_subs = sum(1 for x, y in zip(p_nt, c_nt)
                          if x != y and x not in "-N" and y not in "-N")
            pcp_rows.append({
                "parent": parent_name,
                "child": child_name,
                "child_is_tip": "Y" if child_name in tip_seq else "N",
                "branch_length": round(child.branch_length or 0.0, 6),
                "support": child_support if child_support is not None
                           else (child.confidence if child.confidence is not None else ""),
                "nt_substitutions": nt_subs,
                "aa_substitutions": len(subs),
                "substitution_list": ";".join(subs[:200]),
            })

    pcp_path = os.path.join(args.outdir, f"{args.prefix}_parent_child_pairs.tsv")
    with open(pcp_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pcp_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(pcp_rows)

    print(f"tree       -> {tree_path}")
    print(f"ancestral  -> {anc_path} ({len(ancestors)} internal nodes)")
    print(f"PCPs       -> {pcp_path} ({len(pcp_rows)} branches)")
    print(f"total aa substitutions across backbone: "
          f"{sum(r['aa_substitutions'] for r in pcp_rows)}")


if __name__ == "__main__":
    main()
