#!/usr/bin/env python3
"""Split a multi-segment GenBank flatfile into per-segment FASTA sets.

Keeps only records that carry a recognised segment assignment and whose length
clears a per-segment completeness floor, so partial gene submissions do not
enter the backbone. Writes one FASTA per segment plus a manifest TSV.

Segment names come from `--segments` rather than being hard-coded, so a family
labelled 1/2/3 or RNA1..RNA4 works the same way as an L/M/S bunyavirus.
"""
import argparse
import csv
import os
import re

from Bio import SeqIO

DEFAULT_SEGMENTS = ["L", "M", "S"]

# Expected complete lengths for the Cystoviridae testbed: L ~6.4 kb, M ~4.0 kb,
# S ~2.9 kb. Floors sit in the empty gap between the partial and complete modes
# of the observed length distribution (see dev/plans_and_progress notes).
# Override per family with --floors; Hantaviridae needs L:6000,M:3000,S:1500.
DEFAULT_FLOORS = {"L": 5500, "M": 3000, "S": 2000}


def segment_regex(names):
    """Match 'segment X' for any of the segment names this run expects."""
    alt = "|".join(re.escape(n) for n in sorted(names, key=len, reverse=True))
    return re.compile(rf"segment\s+({alt})\b", re.IGNORECASE)


def segment_of(rec, names, seg_re):
    for feat in rec.features:
        if feat.type != "source":
            continue
        for val in feat.qualifiers.get("segment", []):
            val = val.strip().upper()
            if val in names:
                return val
            m = seg_re.search(val)
            if m:
                return m.group(1).upper()
    m = seg_re.search(rec.description or "")
    return m.group(1).upper() if m else ""


def source_qual(rec, key):
    for feat in rec.features:
        if feat.type == "source":
            vals = feat.qualifiers.get(key, [])
            if vals:
                return vals[0]
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument(
        "--segments",
        default=",".join(DEFAULT_SEGMENTS),
        help="segment names to keep, e.g. L,M,S",
    )
    ap.add_argument(
        "--floors",
        default=",".join(f"{k}:{v}" for k, v in DEFAULT_FLOORS.items()),
        help="per-segment minimum length, e.g. L:5500,M:3000,S:2000",
    )
    args = ap.parse_args()

    names = {s.strip().upper() for s in args.segments.split(",") if s.strip()}
    if not names:
        raise SystemExit("--segments is empty")

    # A floor for a segment this run does not want is ignored; a segment with
    # no floor keeps everything, which is the honest default for a family
    # whose complete length is not known yet.
    floors = {seg: DEFAULT_FLOORS.get(seg, 0) for seg in names}
    for item in args.floors.split(","):
        if not item.strip():
            continue
        seg, val = item.split(":")
        seg = seg.strip().upper()
        if seg in floors:
            floors[seg] = int(val)

    seg_re = segment_regex(names)
    os.makedirs(args.outdir, exist_ok=True)
    kept = {seg: [] for seg in names}
    manifest = []

    with open(args.input) as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            seg = segment_of(rec, names, seg_re)
            length = len(rec.seq)
            reason = ""
            if not seg:
                reason = "no_segment"
            elif seg not in names:
                reason = "segment_not_requested"
            elif length < floors.get(seg, 0):
                reason = "below_floor"

            manifest.append(
                {
                    "accession": rec.id,
                    "segment": seg,
                    "length": length,
                    "organism": rec.annotations.get("organism", ""),
                    "strain": source_qual(rec, "strain") or source_qual(rec, "isolate"),
                    "n_cds": sum(1 for f in rec.features if f.type == "CDS"),
                    "kept": "N" if reason else "Y",
                    "drop_reason": reason,
                }
            )
            if not reason:
                rec.description = rec.annotations.get("organism", "")
                kept[seg].append(rec)

    for seg, recs in kept.items():
        path = os.path.join(args.outdir, f"segment_{seg}.fasta")
        SeqIO.write(recs, path, "fasta")
        print(f"segment {seg}: {len(recs)} records -> {path}")

    man_path = os.path.join(args.outdir, "segment_manifest.tsv")
    with open(man_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(manifest[0]), delimiter="\t")
        w.writeheader()
        w.writerows(manifest)
    print(f"manifest -> {man_path}")


if __name__ == "__main__":
    main()
