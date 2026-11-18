#!/usr/bin/env python3
"""Split taxa into a divergent ML backbone and dense parsimony clusters.

Pairwise divergence is computed over the codon supermatrix, ignoring columns
where either sequence has a gap, so a taxon missing a gene is not scored as
divergent on that account.

Taxa within `--cutoff` of each other form a cluster; each cluster contributes
one representative to the ML backbone and the rest are placed by parsimony
against that representative later. A cluster of one is simply a backbone taxon.

The representative is the cluster member with the lowest mean distance to the
rest of its cluster - its medoid - and, among ties, the one with the fewest
ambiguous sites, so the backbone is anchored on the most complete sequence
available rather than on whichever accession sorted first.
"""
import argparse
import csv
import itertools
import os
import sys

from Bio import SeqIO

GAPS = set("-nN?")


def p_distance(a, b):
    """Uncorrected p-distance over columns where both have a called base."""
    same = compared = 0
    for x, y in zip(a, b):
        if x in GAPS or y in GAPS:
            continue
        compared += 1
        if x == y:
            same += 1
    if compared == 0:
        return 1.0, 0
    return 1.0 - (same / compared), compared


def jc69(p):
    """Jukes-Cantor correction. Saturates, so report both."""
    import math
    if p >= 0.75:
        return float("inf")
    return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--alignment", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--cutoff", default="auto",
                    help="max p-distance within a dense cluster, or 'auto' to "
                         "place it in the widest empty band of the distance "
                         "distribution")
    ap.add_argument("--auto-ceiling", type=float, default=0.35,
                    help="auto mode ignores gaps above this distance")
    ap.add_argument("--auto-min-band-frac", type=float, default=0.15,
                    help="auto mode accepts the widest empty band only if it "
                         "spans at least this fraction of the observed "
                         "distance range; guards against a continuum")
    ap.add_argument("--min-compared", type=int, default=200,
                    help="minimum overlapping columns to trust a distance")
    ap.add_argument("--prefix", default="partition")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(args.alignment, "fasta")}
    ids = sorted(recs)
    if len(ids) < 2:
        sys.exit("need at least two sequences")

    dist = {}
    for a, b in itertools.combinations(ids, 2):
        p, compared = p_distance(recs[a], recs[b])
        dist[(a, b)] = dist[(b, a)] = (p, compared)

    def d(a, b):
        return 0.0 if a == b else dist[(a, b)][0]

    # Pick the cutoff from the data when asked. Within-population and
    # between-species divergence are usually separated by a band that contains
    # no pairs at all; the midpoint of the widest such band is a cutoff that
    # does not depend on knowing the family in advance. Falls back to 1.5% if
    # the distribution has no clear gap, which is the signal that the dataset
    # is a divergence continuum rather than clusters.
    #
    # The band is judged against the spread of the data, not against its own
    # lower edge. Measuring it against the lower edge assumes dense clusters
    # are tight, and that assumption is unsatisfiable above a certain scale:
    # requiring gap >= 2*lo means hi >= 3*lo, so no band whose lower edge
    # exceeds auto_ceiling/3 could ever qualify however empty it was. A
    # species diverse enough to spread its own members past that bound — e.g.
    # Puumala orthohantavirus on segment L, whose isolates sit 8-16% apart —
    # was therefore forced to the fallback and shattered into singletons.
    auto_report = {}
    if str(args.cutoff).strip().lower() == "auto":
        observed = sorted({dist[(a, b)][0]
                           for a, b in itertools.combinations(ids, 2)
                           if dist[(a, b)][1] >= args.min_compared})
        below = [v for v in observed if v <= args.auto_ceiling]
        best_gap, best_lo, best_hi = 0.0, None, None
        for lo, hi in zip(below, below[1:]):
            if hi - lo > best_gap:
                best_gap, best_lo, best_hi = hi - lo, lo, hi
        spread = (observed[-1] - observed[0]) if len(observed) > 1 else 0.0
        gap_frac = (best_gap / spread) if spread > 0 else 0.0
        # A band counts if it is wide in absolute terms and holds a real share
        # of the observed range. A continuum has no such gap at any scale.
        if (best_lo is not None and best_gap >= 0.01
                and gap_frac >= args.auto_min_band_frac):
            cutoff = (best_lo + best_hi) / 2.0
            auto_report = {"mode": "auto", "band_lo": best_lo,
                           "band_hi": best_hi, "band_gap": best_gap,
                           "band_frac": gap_frac,
                           "min_band_frac": args.auto_min_band_frac,
                           "cutoff": cutoff}
            print(f"auto cutoff: widest empty band {best_lo:.4f}..{best_hi:.4f} "
                  f"({gap_frac:.1%} of range) -> cutoff {cutoff:.4f}")
        else:
            cutoff = 0.015
            if best_lo is None:
                reason = "no gap in the distance distribution"
            elif best_gap < 0.01:
                reason = f"widest gap {best_gap:.4f} below the 0.01 floor"
            else:
                reason = (f"widest band {best_lo:.4f}..{best_hi:.4f} spans "
                          f"{gap_frac:.1%} of the range, under the "
                          f"{args.auto_min_band_frac:.0%} bar")
            auto_report = {"mode": "auto_fallback", "cutoff": cutoff,
                           "reason": reason, "band_lo": best_lo,
                           "band_hi": best_hi, "band_gap": best_gap,
                           "band_frac": gap_frac,
                           "min_band_frac": args.auto_min_band_frac}
            print(f"auto cutoff: {reason}; falling back to 1.5%")
    else:
        cutoff = float(args.cutoff)
        auto_report = {"mode": "fixed", "cutoff": cutoff}
    args.cutoff = cutoff

    # single-linkage clustering at the cutoff
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in itertools.combinations(ids, 2):
        p, compared = dist[(a, b)]
        if compared < args.min_compared:
            continue
        if p <= args.cutoff:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

    clusters = {}
    for i in ids:
        clusters.setdefault(find(i), []).append(i)

    rows = []
    reps = []
    for idx, members in enumerate(
        sorted(clusters.values(), key=lambda m: (-len(m), m[0])), 1
    ):
        if len(members) == 1:
            rep = members[0]
        else:
            def score(m):
                mean_d = sum(d(m, o) for o in members if o != m) / (len(members) - 1)
                ambiguous = sum(1 for ch in recs[m] if ch in GAPS)
                return (mean_d, ambiguous, m)
            rep = min(members, key=score)
        reps.append(rep)
        spread = max((d(x, y) for x, y in itertools.combinations(members, 2)), default=0.0)
        for m in members:
            rows.append({
                "cluster": f"C{idx:03d}",
                "accession": m,
                "role": "representative" if m == rep else "dense_member",
                "representative": rep,
                "dist_to_rep": round(d(m, rep), 6),
                "cluster_size": len(members),
                "cluster_max_spread": round(spread, 6),
            })

    out = os.path.join(args.outdir, f"{args.prefix}_assignment.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    mat = os.path.join(args.outdir, f"{args.prefix}_distances.tsv")
    with open(mat, "w") as fh:
        fh.write("accession\t" + "\t".join(ids) + "\n")
        for a in ids:
            fh.write(a + "\t" + "\t".join(f"{d(a, b):.6f}" for b in ids) + "\n")

    sizes = sorted((len(m) for m in clusters.values()), reverse=True)
    print(f"cutoff {args.cutoff:.3%}: {len(ids)} taxa -> {len(clusters)} clusters")
    print(f"  backbone representatives : {len(reps)}")
    print(f"  dense members to graft   : {len(ids) - len(reps)}")
    print(f"  cluster sizes            : {sizes[:12]}{'...' if len(sizes) > 12 else ''}")
    import json
    with open(os.path.join(args.outdir, f"{args.prefix}_cutoff.json"), "w") as fh:
        json.dump({**auto_report, "n_taxa": len(ids), "n_clusters": len(clusters),
                   "n_representatives": len(reps),
                   "n_dense_members": len(ids) - len(reps)}, fh, indent=2)
    print(f"assignment -> {out}")
    print(f"distances  -> {mat}")


if __name__ == "__main__":
    main()
