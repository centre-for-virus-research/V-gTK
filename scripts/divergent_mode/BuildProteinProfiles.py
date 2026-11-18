#!/usr/bin/env python3
"""Build labelled protein orthologue groups from annotated reference genomes.

Product strings in cystovirus GenBank records are not usable as gene keys:
the same phrase ("prohead", "hypothetical protein") covers different genes and
`/gene` is mostly absent. So the groups are derived from the sequences —
all-vs-all blastp, single-linkage clustering on reciprocal coverage — and only
the *naming* uses annotation, taken from a seed genome whose notes carry the
canonical P-numbers.

Output is one protein FASTA per orthologue group, which is the tblastn query
set that ExtractDivergentCds.py uses to pull CDS out of unannotated genomes.

Generalises to any divergent family: supply annotated references and one seed
accession whose annotation names the genes.
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PNUM_RE = re.compile(r"\bP(\d+[a-z]?)\b")


def canonical_name(feat):
    """Pull a canonical P-number out of a CDS feature's annotation."""
    for key in ("note", "product", "gene"):
        for val in feat.qualifiers.get(key, []):
            m = PNUM_RE.search(val)
            if m:
                return f"P{m.group(1)}"
    return ""


def segment_of(rec):
    for feat in rec.features:
        if feat.type == "source":
            for val in feat.qualifiers.get("segment", []):
                val = val.strip().upper()
                if val in {"L", "M", "S"}:
                    return val
    m = re.search(r"segment\s+([LMS])\b", rec.description or "", re.IGNORECASE)
    return m.group(1).upper() if m else ""


def extract_proteins(gb_path, accessions):
    """Translated CDS for the given accessions, tagged with segment and name."""
    out = []
    wanted = set(accessions)
    for rec in SeqIO.parse(gb_path, "genbank"):
        if rec.id not in wanted:
            continue
        seg = segment_of(rec)
        for idx, feat in enumerate(f for f in rec.features if f.type == "CDS"):
            trans = feat.qualifiers.get("translation", [None])[0]
            if not trans:
                try:
                    trans = str(feat.extract(rec.seq).translate(table=11, to_stop=True))
                except Exception:
                    continue
            if not trans or len(trans) < 30:
                continue
            pid = f"{rec.id}|{seg}|cds{idx}"
            out.append(
                SeqRecord(
                    Seq(trans),
                    id=pid,
                    description=canonical_name(feat),
                )
            )
    return out


def all_vs_all(proteins, workdir, threads, evalue):
    faa = os.path.join(workdir, "all_proteins.faa")
    SeqIO.write(proteins, faa, "fasta")
    subprocess.run(
        ["makeblastdb", "-in", faa, "-dbtype", "prot", "-out", os.path.join(workdir, "db")],
        check=True, capture_output=True, text=True,
    )
    tsv = os.path.join(workdir, "all_vs_all.tsv")
    subprocess.run(
        [
            "blastp", "-query", faa, "-db", os.path.join(workdir, "db"),
            "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
            "-evalue", str(evalue), "-num_threads", str(threads),
            "-max_target_seqs", "5000", "-out", tsv,
        ],
        check=True, capture_output=True, text=True,
    )
    return tsv


def cluster(proteins, tsv, min_cov, min_ident):
    """Single-linkage clustering on hits clearing coverage + identity."""
    parent = {p.id: p.id for p in proteins}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    seg_of = {p.id: p.id.split("|")[1] for p in proteins}
    with open(tsv) as fh:
        for line in fh:
            q, s, pident, length, qlen, slen, _, _ = line.rstrip("\n").split("\t")
            if q == s or q not in parent or s not in parent:
                continue
            # Orthologues of a segmented virus live on the same segment; a cross
            # segment hit is a shared domain, not the same gene.
            if seg_of[q] != seg_of[s]:
                continue
            length, qlen, slen = int(length), int(qlen), int(slen)
            cov = length / max(qlen, slen)
            if cov >= min_cov and float(pident) >= min_ident:
                union(q, s)

    groups = {}
    for p in proteins:
        groups.setdefault(find(p.id), []).append(p)
    return groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-g", "--genbank", required=True, help="multi-record GenBank flatfile")
    ap.add_argument("-a", "--accessions", required=True,
                    help="comma-separated reference accessions to mine")
    ap.add_argument("-s", "--seed", default="",
                    help="comma-separated accessions whose annotation names groups")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-coverage", type=float, default=0.5)
    ap.add_argument("--min-identity", type=float, default=20.0)
    args = ap.parse_args()

    accessions = [a.strip() for a in args.accessions.split(",") if a.strip()]
    seeds = {a.strip() for a in args.seed.split(",") if a.strip()}

    proteins = extract_proteins(args.genbank, accessions)
    if not proteins:
        sys.exit("no CDS translations found for the requested accessions")
    print(f"extracted {len(proteins)} proteins from {len(accessions)} references")

    os.makedirs(args.outdir, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix="orthogroups_")
    try:
        tsv = all_vs_all(proteins, workdir, args.threads, args.evalue)
        groups = cluster(proteins, tsv, args.min_coverage, args.min_identity)
    finally:
        pass  # workdir cleaned below after reading nothing further

    rows = []
    used = set()
    for members in sorted(groups.values(), key=lambda m: -len(m)):
        seg = members[0].id.split("|")[1]
        # Name from a seed genome's annotation where one is in the group.
        names = [m.description for m in members
                 if m.id.split("|")[0] in seeds and m.description]
        if names:
            base = sorted(set(names))[0]
        else:
            base = "OG"
        label = f"{seg}_{base}"
        n = 2
        while label in used:
            label = f"{seg}_{base}_{n}"
            n += 1
        used.add(label)

        path = os.path.join(args.outdir, f"{label}.faa")
        SeqIO.write(members, path, "fasta")
        rows.append(
            {
                "group": label,
                "segment": seg,
                "n_members": len(members),
                "min_len": min(len(m.seq) for m in members),
                "max_len": max(len(m.seq) for m in members),
                "members": ";".join(m.id for m in members),
            }
        )

    man = os.path.join(args.outdir, "orthogroup_manifest.tsv")
    with open(man, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    import shutil
    shutil.rmtree(workdir, ignore_errors=True)

    for r in rows:
        print(f"  {r['group']:<18} n={r['n_members']:<3} len {r['min_len']}-{r['max_len']}")
    print(f"\n{len(rows)} orthologue groups -> {args.outdir}")


if __name__ == "__main__":
    main()
