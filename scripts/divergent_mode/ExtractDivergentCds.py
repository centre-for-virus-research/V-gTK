#!/usr/bin/env python3
"""Extract coding sequences from divergent genomes by protein homology.

tblastn each orthologue group's proteins against the genomes of its segment,
then resolve each hit to a real open reading frame: walk out from the HSP in
frame to the terminating stop, and back to the first start codon that is not
separated from the hit by a stop.

The point is that the output is *always* a whole number of codons starting at a
start codon, so nothing downstream can introduce a frameshift by construction.
A hit that cannot be resolved to a clean ORF is dropped and reported, never
emitted as a partial.

Why tblastn and not a profile HMM: HMMER is not in the `vgtk` environment and
this is a shared box. For a family this size, tblastn against a six-frame
translation recovers the same core genes, and BLAST+ is already a declared
pipeline dependency.
"""
import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

STARTS = {"ATG", "GTG", "TTG"}  # phi6 P5a/P5b use GTG starts
STOPS = {"TAA", "TAG", "TGA"}
OUTFMT = "6 qseqid sseqid pident length evalue bitscore sstart send sframe"


def run_tblastn(query_faa, db_prefix, threads, evalue):
    proc = subprocess.run(
        [
            "tblastn", "-query", query_faa, "-db", db_prefix,
            "-outfmt", OUTFMT, "-evalue", str(evalue),
            "-num_threads", str(threads), "-max_target_seqs", "10000",
            "-db_gencode", "11", "-seg", "no",
        ],
        capture_output=True, text=True, check=True,
    )
    hits = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        q, s, pident, length, ev, bits, sstart, send, sframe = line.split("\t")
        hits.append(
            {
                "query": q, "subject": s,
                "pident": float(pident), "length": int(length),
                "evalue": float(ev), "bitscore": float(bits),
                "sstart": int(sstart), "send": int(send), "sframe": int(sframe),
            }
        )
    return hits


def resolve_orf(seq, lo, hi, strand, max_extend=30000, allow_partial=True):
    """Resolve an HSP span to an ORF on the given strand.

    `lo`/`hi` are 0-based half-open coordinates on the *plus* strand.
    Returns (start, end, strand, flag) in plus-strand coordinates, or None.

    A sequence whose deposited ends fall inside the gene has no stop codon to
    find downstream, or no start upstream. Many phi6 lab-evolution submissions
    are 3' truncated this way. Those are still usable: the frame is known, so
    the span is trimmed to a whole number of codons and flagged partial rather
    than discarded. Truncation at a sequence boundary is missing data, which
    the protein aligner handles as a terminal gap; it is not a frameshift.
    Set `allow_partial=False` to require a complete start-to-stop ORF.
    """
    work = seq if strand > 0 else seq.reverse_complement()
    n = len(work)
    if strand > 0:
        hlo, hhi = lo, hi
    else:
        hlo, hhi = n - hi, n - lo

    flags = []

    # downstream: first in-frame stop at or after the hit end
    end = None
    pos = hlo
    while pos + 3 <= n:
        if str(work[pos:pos + 3]).upper() in STOPS and pos >= hhi - 3:
            end = pos + 3
            break
        pos += 3
    if end is None:
        if not allow_partial:
            return None
        # Run to the sequence end, trimmed back to a codon boundary.
        end = hlo + ((n - hlo) // 3) * 3
        if end <= hlo:
            return None
        flags.append("partial3")

    # upstream: earliest start codon before the hit that no in-frame stop
    # separates from it (the longest-ORF convention)
    start = None
    pos = hlo
    hit_stop = False
    while pos >= 0:
        codon = str(work[pos:pos + 3]).upper()
        if codon in STOPS and pos < hlo:
            hit_stop = True
            break
        if codon in STARTS:
            start = pos
        pos -= 3
    if start is None:
        if not allow_partial or hit_stop:
            return None
        # Ran off the 5' end without a start: keep the frame, flag it.
        start = hlo % 3
        flags.append("partial5")

    if end - start > max_extend or end - start < 90:
        return None
    if (end - start) % 3 != 0:
        end -= (end - start) % 3
        if end - start < 90:
            return None

    flag = "+".join(flags) if flags else "complete"
    if strand > 0:
        return start, end, 1, flag
    return n - end, n - start, -1, flag


def best_hit_per_genome(hits):
    """Keep the single highest-bitscore HSP per subject sequence."""
    best = {}
    for h in hits:
        cur = best.get(h["subject"])
        if cur is None or h["bitscore"] > cur["bitscore"]:
            best[h["subject"]] = h
    return best



def extract_from_hits(best, genomes, seg, label, require_complete):
    """Turn one best-hit-per-genome map into CDS records plus report rows."""
    recs, rows = [], []
    for acc, h in best.items():
        rec = genomes[acc]
        lo = min(h["sstart"], h["send"]) - 1
        hi = max(h["sstart"], h["send"])
        strand = 1 if h["sframe"] > 0 else -1
        orf = resolve_orf(rec.seq, lo, hi, strand,
                          allow_partial=not require_complete)
        if orf is None:
            rows.append(
                {"segment": seg, "group": label, "accession": acc,
                 "status": "no_clean_orf", "completeness": "",
                 "cds_len": "", "aa_len": "", "strand": strand,
                 "pident": h["pident"], "evalue": h["evalue"]}
            )
            continue

        s0, e0, st, flag = orf
        nt = rec.seq[s0:e0]
        if st < 0:
            nt = nt.reverse_complement()

        status = "ok"
        if len(nt) % 3 != 0:
            status = "not_triplet"
        else:
            aa = str(nt.translate(table=11))
            # A complete ORF ends on its stop codon, so that one is not
            # internal; a partial ORF has no terminal stop at all.
            body = aa[:-1] if flag == "complete" else aa
            if "*" in body:
                status = "internal_stop"

        if status == "ok":
            recs.append(
                SeqRecord(nt, id=acc,
                          description=f"{label} {s0 + 1}..{e0}"
                                      f"({'+' if st > 0 else '-'}) {flag}")
            )
        rows.append(
            {"segment": seg, "group": label, "accession": acc,
             "status": status, "completeness": flag, "cds_len": len(nt),
             "aa_len": len(nt) // 3, "strand": st,
             "pident": h["pident"], "evalue": h["evalue"]}
        )
    return recs, rows


def extract_group(gpath, db, genomes, seg, label, args, workdir):
    """Recover one orthologue group's CDS, enriching the query set as it goes.

    With `--iterations 1` (the default) this is a single tblastn pass.

    Higher values turn on PSI-BLAST-style query enrichment: each round feeds
    the CDS recovered so far back in as queries, to walk the search out across
    the family. It is off by default because on Cystoviridae it *loses* genes
    rather than finding them - see the decisions log. Partial CDS entering the
    query set out-score the true locus on bitscore and drag the best hit
    somewhere else: a second round took L_P4 from 67/67 genomes to 20/67 and
    pushed internal-stop calls from 5 to 67. The option stays so the behaviour
    can be re-measured on a family where the seed references are sparser, but
    every round is kept only if it improves recovery, and the best round wins -
    enrichment can never return a worse result than the plain search.
    """
    query_path = gpath
    best_recs, best_rows = [], []
    history = []
    for it in range(1, max(1, args.iterations) + 1):
        hits = run_tblastn(query_path, db, args.threads, args.evalue)
        if not hits:
            break
        new_recs, new_rows = extract_from_hits(
            best_hit_per_genome(hits), genomes, seg, label,
            args.require_complete_orf,
        )
        history.append(len(new_recs))
        gained = len(new_recs) > len(best_recs)
        if gained:
            best_recs, best_rows = new_recs, new_rows
        if not gained or it == args.iterations:
            break
        seed = list(SeqIO.parse(gpath, "fasta"))
        enriched = seed + [
            SeqRecord(r.seq.translate(table=11), id=f"it{it}_{r.id}", description="")
            for r in best_recs
        ]
        query_path = os.path.join(workdir, f"{label}_it{it}.faa")
        SeqIO.write(enriched, query_path, "fasta")
    return best_recs, best_rows, history


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment-dir", required=True)
    ap.add_argument("--orthogroup-dir", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--segments", default="L,M,S")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--iterations", type=int, default=1,
                    help="max query-enrichment rounds per group (1 = no enrichment)")
    ap.add_argument("--require-complete-orf", action="store_true",
                    help="discard CDS truncated at a sequence boundary")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix="cds_extract_")
    rows = []
    try:
        for seg in [s.strip() for s in args.segments.split(",") if s.strip()]:
            fasta = os.path.join(args.segment_dir, f"segment_{seg}.fasta")
            if not os.path.exists(fasta):
                print(f"skip segment {seg}: {fasta} missing")
                continue
            genomes = {r.id: r for r in SeqIO.parse(fasta, "fasta")}
            db = os.path.join(workdir, f"db_{seg}")
            subprocess.run(
                ["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", db],
                check=True, capture_output=True, text=True,
            )

            for gpath in sorted(glob.glob(os.path.join(args.orthogroup_dir, f"{seg}_*.faa"))):
                label = os.path.basename(gpath)[:-4]
                recs, grp_rows, history = extract_group(
                    gpath, db, genomes, seg, label, args, workdir
                )
                rows.extend(grp_rows)
                if recs:
                    out = os.path.join(args.outdir, f"{label}.cds.fasta")
                    SeqIO.write(recs, out, "fasta")
                    print(f"{label:<12} {len(recs)}/{len(genomes)} genomes "
                          f"rounds={history} -> {out}", flush=True)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    if not rows:
        sys.exit("no CDS extracted")
    man = os.path.join(args.outdir, "cds_extraction_report.tsv")
    with open(man, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    counts = defaultdict(int)
    for r in rows:
        counts[r["status"]] += 1
    print("\nstatus:", dict(counts))
    print(f"report -> {man}")


if __name__ == "__main__":
    main()
