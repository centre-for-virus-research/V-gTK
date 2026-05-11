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


def _normalize_accession(acc: str) -> str:
    if not acc:
        return ""
    return acc.strip().split()[0].split(".")[0]


def _validate_nextalign_dir(nextalign_dir: str) -> Path:
    path = Path(nextalign_dir)
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f"Nextalign directory not found: {nextalign_dir}")
    return path


def collect_unprojectable_queries(nextalign_dir: str) -> dict:
    """
    Collect queries aligned to references that are not present in any reference_aln
    master-projected alignment. These queries are dropped by PadAlignment and must
    be excluded downstream.

    Returns dict:
        {seq_name: {"reference": ref_id, "error": reason, "warnings": ""}}
    """
    filtered = {}
    nextalign_path = Path(nextalign_dir)

    aligned_reference_ids = set()
    for ref_aln in nextalign_path.glob("reference_aln/*/*.aligned.fasta"):
        try:
            with open(ref_aln, "r", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith(">"):
                        aligned_reference_ids.add(_normalize_accession(line[1:]))
        except Exception as e:
            print(f"Warning: Could not parse {ref_aln}: {e}")

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


def collect_filtered_sequences(nextalign_dir: str, output_file: str) -> dict:
    """
    Scan nextalign directory for .errors.csv files and collect filtered sequence IDs.
    
    Returns dict with structure:
        {seq_name: {"reference": ref_id, "error": error_message, "warnings": warnings}}
    """
    filtered = {}
    nextalign_path = _validate_nextalign_dir(nextalign_dir)
    
    # Search in both query_aln and reference_aln subdirectories
    for errors_file in nextalign_path.rglob("*.errors.csv"):
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
                    if seq_name and errors:
                        filtered[seq_name] = {
                            "reference": ref_id,
                            "error": errors,
                            "warnings": warnings,
                        }
        except Exception as e:
            if isinstance(e, ValueError):
                raise
            print(f"Warning: Could not read {errors_file}: {e}")
    
    # Add sequences that are silently dropped later during PadAlignment because
    # their query reference was not alignable in reference_aln.
    unprojectable = collect_unprojectable_queries(nextalign_dir)
    for seq_name, info in unprojectable.items():
        if seq_name in filtered:
            continue
        filtered[seq_name] = info

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
    
    args = parser.parse_args()
    
    output_path = os.path.join(args.base_dir, args.output)
    
    try:
        print(f"Scanning {args.nextalign_dir} for filtered sequences...")
        filtered = collect_filtered_sequences(args.nextalign_dir, output_path)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    
    print(f"Found {len(filtered)} filtered sequences")
    
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
