#!/usr/bin/env python3
"""
Collect all sequences that were filtered/failed during nextalign alignment.
Reads .errors.csv files from nextalign output and writes a summary of filtered sequences.
"""

import os
import csv
import sys
import argparse
from pathlib import Path

# Sibling import. This script is always invoked as `python scripts/CollectFilteredSequences.py`
# (vgtk-init.nf:626) so scripts/ is sys.path[0], and tests/conftest.py puts it on
# sys.path too. The explicit insert covers a caller that imports this module from
# elsewhere without doing either.
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

import projectability


def _normalize_accession(acc: str) -> str:
    if not acc:
        return ""
    return acc.strip().split()[0].split(".")[0]


def _validate_nextalign_dir(nextalign_dir: str) -> Path:
    path = Path(nextalign_dir)
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f"Nextalign directory not found: {nextalign_dir}")
    return path


def collect_unprojectable_queries(nextalign_dir: str, projectable_ids=None) -> dict:
    """
    Collect queries whose BLAST reference cannot be projected into the merged
    segment alignment. These queries would be dropped by PadAlignment and must be
    excluded downstream.

    ``projectable_ids`` is the authoritative set of reference accessions that
    PadAlignment will actually project, obtained from :mod:`projectability` by
    resolving the same per-segment backbones PadAlignment opens. Pass it whenever
    the caller knows the backbone layout - i.e. whenever ``--precomputed_ref_dir``
    or ``--ref_list`` is available.

    When it is ``None`` the historical behaviour is kept exactly: the set is the
    union of ``reference_aln/*/*.aligned.fasta`` headers. That is correct for the
    single-master, non-segmented builds (RABV, HCV), which supply no per-segment
    backbones and whose references are close enough to their master for nextalign
    to align them all.

    It was *not* correct for segmented influenza, and this is the defect this
    argument exists to close: influenza segment 4's master is AB573800 (H1N1), so
    nextalign could not align 87 of the 117 HA references to it and every query in
    those 87 groups - 366,623 sequences, including 242,686 of 242,732 H3 - was
    marked unprojectable and deleted, even though refset_4_aln.fasta contains all
    117 references and projects each one correctly. See MISSING_H3_report.md.

    Returns dict:
        {seq_name: {"reference": ref_id, "error": reason, "warnings": ""}}
    """
    filtered = {}
    nextalign_path = Path(nextalign_dir)

    if projectable_ids is None:
        aligned_reference_ids = set()
        for ref_aln in nextalign_path.glob("reference_aln/*/*.aligned.fasta"):
            try:
                with open(ref_aln, "r", encoding="utf-8") as handle:
                    for line in handle:
                        if line.startswith(">"):
                            aligned_reference_ids.add(_normalize_accession(line[1:]))
            except Exception as e:
                print(f"Warning: Could not parse {ref_aln}: {e}")
    else:
        aligned_reference_ids = set(projectable_ids)

    for query_ref_dir in nextalign_path.glob("query_aln/*"):
        if not query_ref_dir.is_dir():
            continue
        ref_id = _normalize_accession(query_ref_dir.name)
        if ref_id in aligned_reference_ids:
            continue

        aln_file = query_ref_dir / f"{query_ref_dir.name}.aligned.fasta"
        if not aln_file.exists():
            continue

        try:
            with open(aln_file, "r", encoding="utf-8") as handle:
                for line in handle:
                    if not line.startswith(">"):
                        continue
                    seq_id = _normalize_accession(line[1:])
                    if not seq_id or seq_id == ref_id:
                        continue
                    filtered[seq_id] = {
                        "reference": ref_id,
                        "error": (
                            "reference not present in master-projected reference_aln; "
                            "query cannot be projected into merged segment alignment"
                        ),
                        "warnings": "",
                    }
        except Exception as e:
            print(f"Warning: Could not parse {aln_file}: {e}")

    return filtered


def collect_high_gap_sequences(nextalign_dir: str, max_gap_proportion: float = 0.5) -> dict:
    """
    Scan nextalign query_aln aligned FASTA files and flag query sequences that have
    a gap proportion exceeding max_gap_proportion. Reference sequences (matching the
    directory name) are always skipped.

    Returns dict: {seq_name: {"reference": ref_id, "error": reason, "warnings": ""}}
    """
    filtered = {}
    nextalign_path = Path(nextalign_dir)

    for query_ref_dir in sorted(nextalign_path.glob("query_aln/*")):
        if not query_ref_dir.is_dir():
            continue
        ref_id = _normalize_accession(query_ref_dir.name)
        aln_file = query_ref_dir / f"{query_ref_dir.name}.aligned.fasta"
        if not aln_file.exists():
            continue

        try:
            with open(aln_file, "r", encoding="utf-8") as handle:
                current_id = None
                seq_parts: list = []
                for line in handle:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if current_id and current_id != ref_id and seq_parts:
                            seq = "".join(seq_parts)
                            if seq:
                                gap_prop = seq.count("-") / len(seq)
                                if gap_prop > max_gap_proportion:
                                    filtered[current_id] = {
                                        "reference": ref_id,
                                        "error": (
                                            f"high_gap_proportion:{gap_prop:.3f} > {max_gap_proportion} "
                                            f"against reference {ref_id}"
                                        ),
                                        "warnings": "",
                                    }
                        current_id = _normalize_accession(line[1:])
                        seq_parts = []
                    else:
                        seq_parts.append(line)
                # process last record
                if current_id and current_id != ref_id and seq_parts:
                    seq = "".join(seq_parts)
                    if seq:
                        gap_prop = seq.count("-") / len(seq)
                        if gap_prop > max_gap_proportion:
                            filtered[current_id] = {
                                "reference": ref_id,
                                "error": (
                                    f"high_gap_proportion:{gap_prop:.3f} > {max_gap_proportion} "
                                    f"against reference {ref_id}"
                                ),
                                "warnings": "",
                            }
        except Exception as e:
            print(f"Warning: Could not parse {aln_file}: {e}")

    return filtered


def count_query_sequences(nextalign_dir: str) -> int:
    """Total query sequences presented to PadAlignment, for the G1 ratio.

    Counts headers in every ``query_aln/<ref>/<ref>.aligned.fasta`` excluding the
    reference's own row, which PadAlignment re-inserts from the backbone rather
    than carrying through. Header-only pass, so it is IO-bound and cheap relative
    to the gap filter, which already streams the same files in full.
    """
    total = 0
    for query_ref_dir in Path(nextalign_dir).glob("query_aln/*"):
        if not query_ref_dir.is_dir():
            continue
        ref_id = _normalize_accession(query_ref_dir.name)
        aln_file = query_ref_dir / f"{query_ref_dir.name}.aligned.fasta"
        if not aln_file.exists():
            continue
        try:
            with open(aln_file, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(">") and _normalize_accession(line[1:]) != ref_id:
                        total += 1
        except OSError as exc:
            print(f"Warning: Could not count {aln_file}: {exc}")
    return total


def resolve_projectable_ids(nextalign_dir: str, precomputed_ref_dir=None, ref_list=None):
    """Return the reference accessions PadAlignment can project, or ``None``.

    ``None`` means "caller supplied nothing to resolve backbones with", and the
    unprojectable filter then falls back to its historical ``reference_aln``
    behaviour. That keeps RABV/HCV and every direct caller byte-identical.

    A resolved-but-empty set is a hard error, not a quiet zero: an empty
    projectable set marks *every* query unprojectable, which is precisely the
    failure that deleted 16% of the influenza database.
    """
    precomputed_ref_dir = projectability.normalise_optional_path(precomputed_ref_dir)
    ref_list = projectability.normalise_optional_path(ref_list)
    if not precomputed_ref_dir and not ref_list:
        return None

    master_segment_map = projectability.load_master_segment_map(ref_list)
    ids, source = projectability.projectable_reference_ids(
        nextalign_dir=nextalign_dir,
        precomputed_ref_dir=precomputed_ref_dir,
        master_segment_map=master_segment_map,
    )
    if not ids:
        raise ValueError(
            "Resolved zero projectable reference accessions from "
            f"precomputed_ref_dir={precomputed_ref_dir!r} ref_list={ref_list!r}. "
            "Refusing to mark every query unprojectable."
        )
    print(
        f"[projectability] {len(ids)} projectable reference(s) resolved from {source} "
        f"({len(master_segment_map)} master->segment mapping(s))"
    )
    return ids


def collect_filtered_sequences(nextalign_dir: str, output_file: str, max_gap_proportion: float = 0.5,
                               precomputed_ref_dir=None, ref_list=None) -> dict:
    """
    Scan nextalign directory for .errors.csv files and collect filtered sequence IDs.

    ``precomputed_ref_dir`` / ``ref_list`` are optional and, when given, make the
    unprojectable filter ask PadAlignment's own question instead of the
    ``reference_aln`` proxy. See :func:`collect_unprojectable_queries`.

    Returns dict with structure:
        {seq_name: {"reference": ref_id, "error": error_message, "warnings": warnings}}
    """
    filtered = {}
    nextalign_path = _validate_nextalign_dir(nextalign_dir)
    
    # Only query_aln failures correspond to candidate query sequences. Reference
    # alignment failures are diagnostics for subject/reference rows and must not
    # be turned into downstream query exclusions.
    for errors_file in nextalign_path.glob("query_aln/*/*.errors.csv"):
        # Get reference ID from parent directory name
        ref_id = errors_file.parent.name
        
        try:
            with open(errors_file, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                required_cols = {"seqName", "errors", "warnings"}
                if not reader.fieldnames or not required_cols.issubset(set(reader.fieldnames)):
                    raise ValueError(
                        f"Malformed nextalign errors file {errors_file}: required columns are seqName, errors, warnings"
                    )
                for row in reader:
                    seq_name = row.get("seqName", "").strip()
                    errors = row.get("errors", "").strip()
                    warnings = row.get("warnings", "").strip()
                    
                    # Only include if there's an actual error (not just warnings)
                    if seq_name and seq_name != ref_id and errors:
                        filtered[seq_name] = {
                            "reference": ref_id,
                            "error": errors,
                            "warnings": warnings,
                        }
        except Exception as e:
            if isinstance(e, ValueError):
                raise
            print(f"Warning: Could not read {errors_file}: {e}")
    
    # Add sequences that would be silently dropped later during PadAlignment
    # because their query reference is absent from every backbone that run opens.
    projectable_ids = resolve_projectable_ids(
        nextalign_dir, precomputed_ref_dir=precomputed_ref_dir, ref_list=ref_list
    )
    unprojectable = collect_unprojectable_queries(nextalign_dir, projectable_ids=projectable_ids)
    for seq_name, info in unprojectable.items():
        if seq_name in filtered:
            continue
        filtered[seq_name] = info

    # Add sequences with too many gaps in their reference-aligned FASTA.
    high_gap = collect_high_gap_sequences(nextalign_dir, max_gap_proportion=max_gap_proportion)
    gap_added = 0
    for seq_name, info in high_gap.items():
        if seq_name not in filtered:
            filtered[seq_name] = info
            gap_added += 1
    if gap_added:
        print(f"[gap_filter] Added {gap_added} sequence(s) with gap proportion > {max_gap_proportion}")

    return filtered


def write_summary(filtered: dict, output_file: str):
    """Write a compact summary for visibility in published outputs."""
    summary_file = output_file.replace(".tsv", "_summary.txt")
    if summary_file == output_file:
        summary_file = output_file + ".summary.txt"

    unprojectable = {
        seq: info for seq, info in filtered.items()
        if "cannot be projected into merged segment alignment" in info.get("error", "")
    }

    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("filtered_sequences_total\t{}\n".format(len(filtered)))
        f.write("unprojectable_queries_total\t{}\n".format(len(unprojectable)))
        if unprojectable:
            refs = sorted({info.get("reference", "") for info in unprojectable.values() if info.get("reference", "")})
            f.write("unprojectable_reference_count\t{}\n".format(len(refs)))
            f.write("unprojectable_reference_examples\t{}\n".format(",".join(refs[:20])))
            seq_examples = sorted(unprojectable.keys())[:20]
            f.write("unprojectable_sequence_examples\t{}\n".format(",".join(seq_examples)))

    return summary_file


def write_filtered_list(filtered: dict, output_file: str):
    """Write filtered sequences to TSV file."""
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else ".", exist_ok=True)
    
    with open(output_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["seq_name", "reference", "error", "warnings"])
        for seq_name, info in sorted(filtered.items()):
            writer.writerow([
                seq_name,
                info["reference"],
                info["error"],
                info["warnings"],
            ])


def write_filtered_ids_only(filtered: dict, output_file: str):
    """Write just the filtered sequence IDs, one per line (for easy exclusion)."""
    ids_file = output_file.replace(".tsv", "_ids.txt")
    if ids_file == output_file:
        ids_file = output_file + ".ids.txt"
    
    with open(ids_file, "w", encoding="utf-8") as f:
        for seq_name in sorted(filtered.keys()):
            f.write(seq_name + "\n")
    
    return ids_file


def main():
    parser = argparse.ArgumentParser(
        description="Collect sequences filtered by nextalign due to errors"
    )
    parser.add_argument(
        "-n", "--nextalign_dir",
        required=True,
        help="Path to Nextalign output directory containing query_aln/ and reference_aln/"
    )
    parser.add_argument(
        "-o", "--output",
        default="filtered_sequences.tsv",
        help="Output TSV file with filtered sequences and reasons"
    )
    parser.add_argument(
        "-b", "--base_dir",
        default=".",
        help="Base directory for output"
    )
    parser.add_argument(
        "--max_gap_proportion",
        type=float,
        default=0.5,
        help="Maximum allowed gap proportion in query_aln aligned FASTA; sequences exceeding this are filtered (default: 0.5)"
    )
    parser.add_argument(
        "--precomputed_ref_dir",
        default=None,
        help=(
            "Directory of per-segment backbone alignments (refset_<segment>_aln.fasta), the same "
            "value PadAlignment receives. When given, projectability is decided against these "
            "backbones instead of Nextalign/reference_aln/, which is the only correct question "
            "for segmented viruses. Omit for non-segmented builds."
        )
    )
    parser.add_argument(
        "--ref_list",
        default=None,
        help=(
            "Ref list TSV (accession, type, segment), the same value PadAlignment receives as -m. "
            "Used only to learn which masters exist and which segment each covers, so the right "
            "backbone file is opened per master. Never used as the reference set itself."
        )
    )
    parser.add_argument(
        "--max_filtered_fraction",
        type=float,
        default=1.0,
        help=(
            "Abort if the filtered fraction of all query sequences exceeds this. The pipeline "
            "passes 0.02; the script default of 1.0 keeps direct/legacy callers unguarded. "
            "Set to 0 to disable the check entirely."
        )
    )

    args = parser.parse_args()

    output_path = os.path.join(args.base_dir, args.output)

    try:
        print(f"Scanning {args.nextalign_dir} for filtered sequences...")
        filtered = collect_filtered_sequences(
            args.nextalign_dir,
            output_path,
            max_gap_proportion=args.max_gap_proportion,
            precomputed_ref_dir=args.precomputed_ref_dir,
            ref_list=args.ref_list,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    print(f"Found {len(filtered)} filtered sequences")

    # Guard G1. Losing a sixth of the database used to be a stderr line and exit 0.
    total_queries = count_query_sequences(args.nextalign_dir)
    if args.max_filtered_fraction and 0 < args.max_filtered_fraction < 1.0 and total_queries:
        fraction = len(filtered) / total_queries
        print(f"[guard] filtered {len(filtered)}/{total_queries} query sequences ({fraction:.4f})")
        if fraction > args.max_filtered_fraction:
            print(
                f"ERROR: filtered fraction {fraction:.4f} exceeds --max_filtered_fraction "
                f"{args.max_filtered_fraction}. {len(filtered)} of {total_queries} query "
                f"sequences would be dropped. Refusing to continue; inspect "
                f"{output_path.replace('.tsv', '_summary.txt')}.",
                file=sys.stderr,
            )
            sys.exit(3)


    if filtered:
        write_filtered_list(filtered, output_path)
        ids_file = write_filtered_ids_only(filtered, output_path)
        summary_file = write_summary(filtered, output_path)
        print(f"Wrote filtered sequences to: {output_path}")
        print(f"Wrote filtered IDs to: {ids_file}")
        print(f"Wrote filtered summary to: {summary_file}")

        unprojectable_count = sum(
            1
            for info in filtered.values()
            if "cannot be projected into merged segment alignment" in info.get("error", "")
        )
        if unprojectable_count:
            print(
                "\n[ALERT] {} query sequence(s) were dropped because their BLAST-hit references "
                "were not present in master-projected reference alignments. "
                "See filtered_sequences_summary.txt and filtered_sequences.tsv for details.\n".format(unprojectable_count),
                file=sys.stderr,
            )
        
        # Print summary
        print("\nFiltered sequences:")
        for seq_name, info in sorted(filtered.items()):
            # Truncate long error messages for display
            error_short = info["error"][:80] + "..." if len(info["error"]) > 80 else info["error"]
            print(f"  {seq_name}: {error_short}")
    else:
        # Write empty files so downstream processes have something to read
        write_filtered_list({}, output_path)
        write_filtered_ids_only({}, output_path)
        write_summary({}, output_path)
        print("No filtered sequences found")


if __name__ == "__main__":
    main()
