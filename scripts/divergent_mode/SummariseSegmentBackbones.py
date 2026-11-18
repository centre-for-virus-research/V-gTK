#!/usr/bin/env python3
"""One row per segment: did it get a backbone, and if not, why not.

Gene selection runs once per segment, and a segment can legitimately drop out
of a run - no gene recovered at all, or none recoverable across the divergent
taxa. That is the right answer on a segmented dsRNA family, where only the
RdRp segment is conserved family-wide. It is a red flag on a negative-sense
segmented family, where every segment should carry a core gene.

Either way it must be said out loud. The failure this file exists to prevent
is a run that builds one sound backbone, ignores the other segments and
reports nothing about them, so the output looks complete and is not.

Reads the per-segment status files written by SelectBackboneGenes.py and the
list of segments that reached a grafted tree, and writes one table. Exits
non-zero if no segment produced a backbone at all, because then there is no
result.
"""
import argparse
import csv
import os
import sys

STAGE_ORDER = {"seed": 0, "selection": 1}

ADVICE = {
    "no_gene_recovered":
        "no orthologue group survived extraction and alignment on this "
        "segment; check cds_extraction_report.tsv and --divergent_min_taxa",
    "below_min_coverage":
        "genes exist but none is recoverable across the divergent taxa; "
        "expected on a dsRNA family outside the RdRp segment, otherwise "
        "widen the reference set or lower --divergent_min_gene_coverage",
    "not_run":
        "no status file; the segment never reached gene selection",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", nargs="*", default=[],
                    help="per-segment status TSVs from SelectBackboneGenes.py")
    ap.add_argument("--segments", required=True,
                    help="comma-separated segments the run was asked for")
    ap.add_argument("--built", default="",
                    help="comma-separated segments that produced a backbone")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--allow-no-backbone", action="store_true",
                    help="exit 0 even when no segment produced a backbone")
    args = ap.parse_args()

    segments = [s.strip() for s in args.segments.split(",") if s.strip()]
    built = {s.strip() for s in args.built.split(",") if s.strip()}

    # Keep the furthest stage each segment reached: a selection-stage status
    # supersedes the seed-stage one for the same segment.
    latest = {}
    for path in args.status:
        if not os.path.exists(path):
            continue
        for row in csv.DictReader(open(path), delimiter="\t"):
            seg = row.get("segment", "")
            if not seg:
                continue
            prev = latest.get(seg)
            if prev is None or STAGE_ORDER.get(row.get("stage", ""), -1) >= \
                    STAGE_ORDER.get(prev.get("stage", ""), -1):
                latest[seg] = row

    out_rows = []
    for seg in segments:
        row = latest.get(seg, {})
        status = row.get("status", "not_run")
        if seg in built:
            status = "ok"
        out_rows.append({
            "segment": seg,
            "backbone_built": "Y" if seg in built else "N",
            "stage_reached": row.get("stage", "none"),
            "status": status,
            "genes": row.get("chosen_genes", ""),
            "n_genes_scored": row.get("n_genes_scored", "0"),
            "best_gene": row.get("best_gene", ""),
            "best_coverage": row.get("best_coverage", ""),
            "note": "" if seg in built else ADVICE.get(status, row.get("detail", "")),
        })

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"{len(built)}/{len(segments)} segments produced a backbone")
    for r in out_rows:
        mark = "ok  " if r["backbone_built"] == "Y" else "SKIP"
        print(f"  [{mark}] {r['segment']:<4} {r['status']:<20} "
              f"{r['genes'] or '-'}")
        if r["note"]:
            print(f"           {r['note']}")
    print(f"summary -> {args.output}")

    missing = [r["segment"] for r in out_rows if r["backbone_built"] == "N"]
    if missing and built:
        print(f"[warn] no backbone for segment(s): {', '.join(missing)}. "
              "This is the correct answer only if those segments carry no gene "
              "recoverable across the divergent taxa - see the note column.",
              file=sys.stderr)

    if not built and not args.allow_no_backbone:
        sys.exit("no segment produced a backbone; there is no result to "
                 "report. See the note column in " + args.output)


if __name__ == "__main__":
    main()
