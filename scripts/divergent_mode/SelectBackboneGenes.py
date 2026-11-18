#!/usr/bin/env python3
"""Choose which genes carry the ML backbone, for one segment at a time.

A gene only helps the deep topology if it is recoverable across the divergent
taxa, not merely across whatever dense cluster dominates the dataset. Scoring
raw taxon counts picks the opposite genes: on Cystoviridae every gene is
present in all 48 phi6 isolates, so raw counts rank them all near-equal while
coverage outside phi6 ranges from 100% down to 10%.

So coverage is measured over *dense clusters* where a partition is supplied -
each cluster counts once, however many near-identical isolates it holds - and
a gene is kept when it clears `--min-coverage` on those. Ties are broken by
alignment length, since a longer gene carries more signal.

The denominator is every genome of this segment that entered extraction, not
only those the partition happens to contain. A partition is built from one
seed gene's alignment, so a genome missing that gene is missing from the
partition too; scoring against the partition's membership would measure each
gene over the taxa that already have the seed gene and return 1.00 for the
seed gene every time. Genomes outside the partition are counted as clusters of
one, which is what a divergent taxon with no near neighbour is.

One segment per call. Segments of a segmented virus have their own accessions,
their own dense clusters and - where the family reassorts - their own
topology, so a partition built on one segment's alignment says nothing about
another's. Calling this with several segments and one partition scores every
gene outside that partition's segment at zero coverage, which reads as "those
segments have no usable gene" when nothing of the sort was measured. The
`--assignment` scope is therefore checked against the requested segment and
the run aborts if they do not match.

Writes the chosen gene names as one comma-separated line, which is what the
Nextflow concatenation step consumes. With `--allow-empty` a segment carrying
no qualifying gene writes an empty file and exits 0, so one unusable segment
does not take down a run over a family where the others are fine; the caller
is expected to read `--status` and say so out loud.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

from Bio import SeqIO


def write_summary(path, scored, chosen_names):
    """The full ranking, including for a segment that is about to be skipped.

    A skipped segment is where this table matters most: it is the evidence for
    the skip, and the number you would lower --min-coverage to.
    """
    if not path or not scored:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["segment", "gene", "coverage", "n_present", "n_taxa",
                        "n_genomes_present", "n_genomes", "nt_columns",
                        "chosen"],
            delimiter="\t",
        )
        w.writeheader()
        for r in scored:
            w.writerow({**r, "chosen": "Y" if r["gene"] in chosen_names else "N"})


def write_status(path, row):
    if not path:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["segment", "stage", "status", "detail", "n_genes_scored",
                        "best_gene", "best_coverage", "chosen_genes"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerow(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--codon-dir", required=True)
    ap.add_argument("--report", required=True, help="cds_extraction_report.tsv")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--segments", default="L,M,S",
                    help="segment to score; one per call (a comma list is "
                         "accepted only without --assignment)")
    ap.add_argument("--min-coverage", type=float, default=0.9)
    ap.add_argument("--max-genes", type=int, default=6)
    ap.add_argument("--summary", default="")
    ap.add_argument("--status", default="",
                    help="one-row TSV describing the outcome for this segment")
    ap.add_argument("--pick-seed", action="store_true",
                    help="write only the single highest-coverage gene, used to "
                         "seed the partition that scopes the real selection")
    ap.add_argument("--allow-empty", action="store_true",
                    help="write an empty gene list and exit 0 when no gene "
                         "qualifies, instead of aborting the run")
    ap.add_argument("--assignment", default="",
                    help="partition_assignment.tsv; collapses each dense "
                         "cluster to one unit of coverage (recommended)")
    ap.add_argument("--manifest", default="",
                    help="segment_manifest.tsv; the coverage denominator is "
                         "every genome kept for this segment. Without it the "
                         "denominator is whatever the extraction report saw, "
                         "which omits genomes where no group hit at all and "
                         "so overstates coverage.")
    args = ap.parse_args()

    segments = [s.strip() for s in args.segments.split(",") if s.strip()]
    stage = "seed" if args.pick_seed else "selection"
    seg_label = ",".join(segments)

    def give_up(status, detail, scored=None):
        """Either abort, or write an empty list and let the caller carry on."""
        scored = scored or []
        write_summary(args.summary, scored, set())
        best = scored[0] if scored else None
        write_status(args.status, {
            "segment": seg_label, "stage": stage, "status": status,
            "detail": detail, "n_genes_scored": len(scored),
            "best_gene": best["gene"] if best else "",
            "best_coverage": f"{best['coverage']:.4f}" if best else "",
            "chosen_genes": "",
        })
        if args.allow_empty:
            open(args.output, "w").close()
            print(f"[skip] segment {seg_label}: {detail}")
            sys.exit(0)
        sys.exit(detail)

    rows = list(csv.DictReader(open(args.report), delimiter="\t"))

    # Scope the coverage denominator to cluster representatives when a
    # partition is available. Without it, a dataset dominated by one dense
    # cluster scores every gene as near-universal: on Cystoviridae all 48 phi6
    # isolates carry every gene, which hides that coverage outside phi6 ranges
    # from 100% (L_P4) down to 10% (M_P3).
    cluster_of = {}
    if args.assignment:
        if len(segments) != 1:
            sys.exit(
                "--assignment collapses one segment's dense clusters, so "
                f"--segments must name exactly one; got '{args.segments}'. "
                "Run this once per segment."
            )
        for r in csv.DictReader(open(args.assignment), delimiter="\t"):
            cluster_of[r["accession"]] = r["cluster"]
        if not cluster_of:
            sys.exit(f"no cluster assignments in {args.assignment}")

        # The partition must be this segment's own. A partition built on
        # another segment shares no accessions with this one, so every gene
        # here would score zero coverage and the segment would look unusable
        # when it was simply never measured.
        seg_accs = {r["accession"] for r in rows if r["segment"] == segments[0]}
        if not (set(cluster_of) & seg_accs):
            sys.exit(
                f"partition {args.assignment} has no accession from segment "
                f"{segments[0]} ({len(cluster_of)} clustered accessions, none "
                "of this segment's). Gene selection must be scoped to a "
                "partition built on the same segment, one segment per call."
            )
    else:
        print("WARNING: no --assignment given; coverage counts every isolate "
              "separately and a dominant dense cluster will skew the choice")

    # Taxa per segment: every genome kept for the segment where a manifest is
    # available, otherwise everything the extractor saw. The manifest is the
    # honest universe - a genome where no group hit at all never reaches the
    # extraction report, and leaving it out inflates every gene's coverage.
    seg_taxa = defaultdict(set)
    if args.manifest:
        for r in csv.DictReader(open(args.manifest), delimiter="\t"):
            if r.get("kept") == "Y" and r.get("segment") in segments:
                seg_taxa[r["segment"]].add(r["accession"])
    for r in rows:
        if r["segment"] in segments and r["accession"] not in \
                seg_taxa.get(r["segment"], set()):
            # Only widen the universe when no manifest pinned it down.
            if not args.manifest:
                seg_taxa[r["segment"]].add(r["accession"])

    def unit(acc):
        """One unit of coverage: a dense cluster, or a lone taxon."""
        return cluster_of.get(acc, f"singleton:{acc}")

    seg_units = {seg: {unit(a) for a in accs} for seg, accs in seg_taxa.items()}
    seg_raw = {seg: len(accs) for seg, accs in seg_taxa.items()}

    present = defaultdict(set)
    present_raw = defaultdict(set)
    for r in rows:
        if r["status"] == "ok" and r["accession"] in seg_taxa.get(r["segment"], ()):
            present[(r["segment"], r["group"])].add(unit(r["accession"]))
            present_raw[(r["segment"], r["group"])].add(r["accession"])

    if args.assignment:
        seg = segments[0]
        print(f"segment {seg}: scoring coverage over "
              f"{len(seg_units.get(seg, ()))} dense clusters "
              f"({seg_raw.get(seg, 0)} genomes)")

    widths = {}
    for path in glob.glob(os.path.join(args.codon_dir, "*.codon.fasta")):
        label = os.path.basename(path)[: -len(".codon.fasta")]
        rec = next(SeqIO.parse(path, "fasta"), None)
        if rec is not None:
            widths[label] = len(rec.seq)

    scored = []
    for (seg, group), units in present.items():
        if seg not in segments or group not in widths:
            continue
        total = len(seg_units.get(seg, ()))
        cov = len(units) / total if total else 0.0
        scored.append({"segment": seg, "gene": group, "coverage": round(cov, 4),
                       "n_present": len(units), "n_taxa": total,
                       "n_genomes_present": len(present_raw[(seg, group)]),
                       "n_genomes": seg_raw.get(seg, 0),
                       "nt_columns": widths[group]})

    scored.sort(key=lambda r: (-r["coverage"], -r["nt_columns"]))

    if not scored:
        give_up("no_gene_recovered",
                f"no aligned gene on segment {seg_label}: nothing was "
                "recovered, or every group fell below --min-taxa in the codon "
                "alignment step")

    if args.pick_seed:
        # The seed only has to be the gene most likely to contain every taxon,
        # so the clustering it feeds sees the whole dataset. Coverage first,
        # length as the tie-break.
        chosen = scored[:1]
    else:
        chosen = [r for r in scored
                  if r["coverage"] >= args.min_coverage][: args.max_genes]

    if not chosen:
        best = scored[0]
        give_up(
            "below_min_coverage",
            f"no gene on segment {seg_label} reaches --min-coverage "
            f"{args.min_coverage}; best is {best['gene']} at "
            f"{best['coverage']:.2f}. Lower --divergent_min_gene_coverage or "
            "widen the reference set.",
            scored,
        )

    with open(args.output, "w") as fh:
        fh.write(",".join(r["gene"] for r in chosen) + "\n")

    write_summary(args.summary, scored, {r["gene"] for r in chosen})

    write_status(args.status, {
        "segment": seg_label, "stage": stage, "status": "ok", "detail": "",
        "n_genes_scored": len(scored),
        "best_gene": scored[0]["gene"],
        "best_coverage": f"{scored[0]['coverage']:.4f}",
        "chosen_genes": ",".join(r["gene"] for r in chosen),
    })

    print(f"backbone genes: {', '.join(r['gene'] for r in chosen)}")
    for r in chosen:
        print(f"  {r['gene']:<10} coverage={r['coverage']:.2f} "
              f"({r['n_present']}/{r['n_taxa']} clusters, "
              f"{r['n_genomes_present']}/{r['n_genomes']} genomes) "
              f"{r['nt_columns']} nt")


if __name__ == "__main__":
    main()
