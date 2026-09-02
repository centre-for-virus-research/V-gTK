#!/usr/bin/env python3
"""Guard: assert each segment kept the sequences BLAST assigned to it.

Nothing in this pipeline used to compare "how many sequences did BLAST assign to
segment N" against "how many reached segment N's merged alignment". Because of
that, a defect in the projectability test deleted 366,623 of segment 4's 572,876
HA sequences and 267,350 of segment 6's 516,485 NA sequences - 99.98% of all H3
and 99.97% of all N1 - and the run still exited 0. The numbers were sitting in
db_summary.txt the whole time:

    Segment | Passed QC | Failed QC
       1    |  462975   |    925
       4    |  205397   | 367894      <- 35.8% retained
       6    |  248412   | 268067      <- 48.1% retained
       7    |  492067   |   1044

Every healthy segment retains >99.8%. This script turns that into a hard stop.

It runs between PAD_ALIGNMENT and DEDUP_ALIGNMENT, so a bad merge is caught
before hours of clustering and tree building are spent on it.

Accounting
----------
retained  = unique query accessions in refset_<N>_aln_merged_MSA.fasta
assigned  = rows in query_uniq_tophit_annotated.tsv whose reference maps to segment N
excused   = assigned accessions present in filtered_sequences_ids.txt

    ratio = retained / (assigned - excused)

"Excused" covers sequences nextalign genuinely could not align (too short,
ambiguous large indels). Those are legitimate losses and must not trip the guard.
Sequences dropped by a *rule* rather than by alignment failure are not excused -
which is the whole point, since that is what went wrong.

Reference rows are excluded from ``retained``: PadAlignment inserts each backbone
row into its own block, so every reference that also has a sub-alignment appears
twice in the merged FASTA (113 duplicate ids in segment 4). Counting unique
accessions and subtracting references avoids both double counting and an
inflated ratio on a segment whose queries all vanished.
"""

import argparse
import os
import re
import sys

import accession_utils


def read_fasta_ids(path):
    """Unique canonical accessions in a FASTA. CRLF-safe.

    All four id streams this guard subtracts from each other - retained,
    assigned, excused and reference - now go through :mod:`accession_utils`.
    Three of them used to canonicalise with ``.split(".")[0]`` and the fourth,
    :func:`load_reference_segments`, not at all, so a reference list spelled
    ``NC_001542.1`` left ``retained`` and ``assigned`` bare while
    ``reference_ids`` stayed versioned. Both subtractions below then removed
    nothing: the backbone row PadAlignment inserts into every block was counted
    as a retained query, the ratio came out too high, and the guard passed a
    segment it exists to fail. Half-normalising is worse than not normalising,
    because it looks like it works.
    """
    ids = set()
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            token = line[1:].strip().split()
            if token:
                ids.add(accession_utils.normalise_accession(token[0]))
    return ids


def load_reference_segments(ref_seg_path):
    """{accession: segment} from a reference->segment TSV.

    Takes the LAST column as the segment, matching BlastAlignment.py's convention,
    so both shapes in this repo work: ref_seq_seg.tsv (accession, segment) and
    ref_list_refmast.txt (accession, accession_type, segment). Rows whose segment
    is not a number (the influenza B/C/D exclusion_list entries) are kept as-is and
    simply never match a refset_<N> file.

    Keys are canonicalised because they are the other side of two subtractions
    and one lookup, all of which are canonical. This is the side that was left
    verbatim, which is what made the guard a partial no-op on a versioned
    reference list. Segment values are not touched - they are not accessions.
    """
    mapping = {}
    with open(ref_seg_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = [c.strip() for c in line.rstrip("\n").rstrip("\r").split("\t")]
            if len(parts) >= 2 and parts[0]:
                mapping[accession_utils.normalise_accession(parts[0])] = parts[-1]
    return mapping


def load_assigned_counts(tophit_path, ref_segments):
    """{segment: set(query accessions)} from query_uniq_tophit_annotated.tsv.

    Column 5 carries the segment when present; fall back to the reference->segment
    map so the check still works on a file written before that column existed.

    The *reference* is canonicalised as well as the query. Only the query used
    to be, so on a file whose reference column carries a version the fallback
    lookup missed the map, the row was dropped as segment-less, and a segment
    whose sequences had all vanished reported ``expected = 0`` - which this
    guard scores as ratio 1.0 and passes.
    """
    assigned = {}
    with open(tophit_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if len(parts) < 2:
                continue
            query = accession_utils.normalise_accession(parts[0]) or ""
            reference = accession_utils.normalise_accession(parts[1]) or ""
            segment = parts[4].strip() if len(parts) >= 5 and parts[4].strip() else ref_segments.get(reference, "")
            if not segment:
                continue
            assigned.setdefault(segment, set()).add(query)
    return assigned


def load_excused(filtered_ids_path):
    """Canonical accessions from filtered_sequences_ids.txt.

    The first whitespace token is kept deliberately: CollectFilteredSequences
    writes one id per line, but a line that ever carried a description must not
    excuse a sequence called ``"ACC description"``. The version is dropped by
    :mod:`accession_utils` rather than at the first dot, so this reader and that
    writer now agree on one spelling instead of two that happened to coincide.
    """
    if not filtered_ids_path or not os.path.isfile(filtered_ids_path):
        return set()
    excused = set()
    with open(filtered_ids_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            token = line.strip().split()
            if token:
                excused.add(accession_utils.normalise_accession(token[0]))
    return excused


def segment_of_merged_file(path):
    """refset_4_aln_merged_MSA.fasta -> '4'. Returns None when unrecognised."""
    match = re.search(r"refset_(\d+)_", os.path.basename(path))
    return match.group(1) if match else None


def check(merged_files, tophit_path, ref_seg_path, filtered_ids_path=None,
          min_retention=0.98):
    ref_segments = load_reference_segments(ref_seg_path)
    reference_ids = set(ref_segments)
    assigned = load_assigned_counts(tophit_path, ref_segments)
    excused = load_excused(filtered_ids_path)

    rows, failures = [], []
    for merged in sorted(merged_files):
        segment = segment_of_merged_file(merged)
        if segment is None:
            print(f"[warn] Cannot infer segment from {merged}; skipping retention check for it")
            continue
        retained_ids = read_fasta_ids(merged) - reference_ids
        expected_ids = assigned.get(segment, set()) - reference_ids - excused
        expected = len(expected_ids)
        retained = len(retained_ids)
        ratio = (retained / expected) if expected else 1.0
        rows.append((segment, retained, expected, ratio))
        if expected and ratio < min_retention:
            failures.append((segment, retained, expected, ratio, sorted(expected_ids - retained_ids)[:5]))

    width = max((len(s) for s, *_ in rows), default=7)
    print(f"{'segment':<{width}}  {'retained':>10}  {'expected':>10}  {'ratio':>7}")
    for segment, retained, expected, ratio in rows:
        print(f"{segment:<{width}}  {retained:>10}  {expected:>10}  {ratio:>7.4f}")

    for segment, retained, expected, ratio, examples in failures:
        print(
            f"ERROR: segment {segment} retained {retained}/{expected} sequences "
            f"({ratio:.4f} < {min_retention}). Missing examples: {', '.join(examples) or 'n/a'}",
            file=sys.stderr,
        )
    return failures


def main():
    parser = argparse.ArgumentParser(
        description="Fail the run if a segment's merged alignment lost sequences BLAST assigned to it"
    )
    parser.add_argument("-m", "--merged", nargs="+", required=True,
                        help="refset_<N>_aln_merged_MSA.fasta files to check")
    parser.add_argument("-t", "--tophit", required=True,
                        help="query_uniq_tophit_annotated.tsv from BLAST_ALIGNMENT")
    parser.add_argument("-s", "--ref_seg", required=True,
                        help="TSV mapping reference accession to segment (ref_seq_seg.tsv)")
    parser.add_argument("-f", "--filtered_ids", default=None,
                        help="filtered_sequences_ids.txt; these are excused from the ratio")
    parser.add_argument("--min_retention", type=float, default=0.98,
                        help="Minimum retained/expected per segment (default: 0.98)")
    args = parser.parse_args()

    failures = check(
        merged_files=args.merged,
        tophit_path=args.tophit,
        ref_seg_path=args.ref_seg,
        filtered_ids_path=args.filtered_ids,
        min_retention=args.min_retention,
    )
    if failures:
        sys.exit(4)
    print("[guard] All segments passed the retention check.")


if __name__ == "__main__":
    main()
