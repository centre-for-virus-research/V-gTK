import os
import subprocess
import shutil
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_ranges(ranges, alignment_length):
    trimmed_ranges = []
    for r in ranges:
        start_end = r.split(":")
        start = int(start_end[0]) if start_end[0] else 1
        end = int(start_end[1]) if len(start_end) > 1 and start_end[1] else alignment_length
        trimmed_ranges.append((start - 1, end))
    return trimmed_ranges


def trim_alignment(input_fasta, output_fasta, ranges):
    trimmed_seqs = []
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        trimmed_seq_parts = []
        alignment_length = len(seq_record.seq)
        trimmed_ranges = parse_ranges(ranges, alignment_length)
        for start, end in trimmed_ranges:
            trimmed_seq_parts.append(str(seq_record.seq[start:end]))
        trimmed_seq = Seq("".join(trimmed_seq_parts))
        trimmed_record = SeqRecord(trimmed_seq, id=seq_record.id, description=seq_record.description)
        trimmed_seqs.append(trimmed_record)
    SeqIO.write(trimmed_seqs, output_fasta, "fasta")


def apply_trimming(cds_file, input_dir, trim_output_dir):
    trimmed_files = {}
    os.makedirs(trim_output_dir, exist_ok=True)
    with open(cds_file, mode='r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        for row in csv_reader:
            input_fasta_name = row['input_fasta']
            input_fasta = os.path.join(input_dir, input_fasta_name)
            basename = os.path.splitext(input_fasta_name)[0]
            output_fasta_name = f"{basename}_trimmed.fas"
            output_fasta = os.path.join(trim_output_dir, output_fasta_name)

            ranges = row['ranges'].split(",")
            trim_alignment(input_fasta, output_fasta, ranges)
            trimmed_files[input_fasta_name] = output_fasta
    return trimmed_files


GAP_CHARACTERS = "-."


def strip_alignment_gaps(input_fasta, output_fasta):
    """Write an unaligned copy of input_fasta, preserving sequence IDs.

    MMseqs is not alignment-aware: it treats '-' as ordinary sequence content.
    In a padded MSA that both breaks the k-mer prefilter (gaps interrupt the
    k-mers used to find candidate neighbours) and drags computed identities
    down, so the same sequences fragment into far more clusters aligned than
    unaligned. Measured on refset_3 (influenza PA) at --min-seq-id 0.95:
    1,137 clusters from a 100k padded MSA vs 442 from the same sequences
    unaligned.
    """
    translation = str.maketrans("", "", GAP_CHARACTERS)
    written = 0
    empty = 0
    with open(output_fasta, "w") as handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequence = str(record.seq).translate(translation)
            if not sequence:
                empty += 1
                continue
            handle.write(f">{record.id}\n{sequence}\n")
            written += 1
    if empty:
        print(f"[warn] Skipped {empty} sequence(s) that were entirely gaps")
    return written


def write_aligned_representatives(alignment_fasta, cluster_tsv, output_fasta):
    """Emit cluster representatives using their ALIGNED sequences.

    Clustering runs on unaligned sequences, but the representatives feed
    VeryFastTree and IQ-TREE, which need equal-length input - IQ_TREE in
    vgtk-init.nf explicitly aborts on a ragged alignment. So the representative
    IDs come from MMseqs while the sequences come from the original padded MSA.
    """
    representative_ids = set()
    with open(cluster_tsv) as handle:
        for line in handle:
            if line.strip():
                representative_ids.add(line.split("\t")[0].strip())

    written = set()
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(alignment_fasta, "fasta"):
            seq_id = str(record.id).strip()
            if seq_id in representative_ids and seq_id not in written:
                out_handle.write(f">{record.id}\n{str(record.seq)}\n")
                written.add(seq_id)

    missing = representative_ids - written
    if missing:
        raise ValueError(
            f"{len(missing)} cluster representative(s) missing from {alignment_fasta}, "
            f"e.g. {sorted(missing)[:5]}"
        )
    return len(written)


def run_mmseqs_clustering(input_fasta, output_dir, min_seq_id, threads=8, strip_gaps=True, max_seqs=None):
    base_name = os.path.splitext(os.path.basename(input_fasta))[0]
    mmseqs_dir = os.path.join(output_dir, base_name)
    segments_db_dir = os.path.join(mmseqs_dir, "segments_DB")
    tmp_dir = os.path.join(mmseqs_dir, "tmp_mmseq2")

    os.makedirs(segments_db_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    os.environ["TMPDIR"] = tmp_dir

    db_path = os.path.join(segments_db_dir, f"{base_name}_db")
    cluster_path = os.path.join(segments_db_dir, f"{base_name}_cluster")
    cluster_seq_path = os.path.join(segments_db_dir, f"{base_name}_cluster_seq")
    cluster_rep_path = os.path.join(segments_db_dir, f"{base_name}_cluster_rep")

    # Cluster the unaligned sequences; keep the padded MSA for the output below.
    #
    # The scratch copy deliberately lives under segments_DB with a non-.fasta
    # extension and is deleted as soon as the DB is built. Several consumers
    # (IQ_TREE and VERY_FAST_TREE in vgtk-init.nf, vgtk-rabv.sh, and
    # UsherPlacement.resolve_non_update_assets) fall back to a generic `*.fasta`
    # glob over this directory, so any stray unaligned FASTA left here can be
    # picked up in place of the representatives - and VeryFastTree has no
    # equal-length guard to catch it.
    cluster_source = input_fasta
    if strip_gaps:
        cluster_source = os.path.join(segments_db_dir, f"{base_name}_ungapped.seq")
        count = strip_alignment_gaps(input_fasta, cluster_source)
        print(f"[info] Clustering {count} gap-stripped sequences (representatives stay aligned)")

    cluster_cmd = [
        "mmseqs", "cluster",
        "--min-seq-id", str(min_seq_id),
        db_path, cluster_path, tmp_dir,
        "--threads", str(threads), '--cov-mode', '2', '--cluster-mode', '2'
    ]
    if max_seqs is not None:
        # Prefilter candidates kept per query. The default can strand a sequence
        # as its own cluster when a large near-identical mass crowds its true
        # neighbours out of the candidate list.
        cluster_cmd += ["--max-seqs", str(max_seqs)]

    try:
        subprocess.run(["mmseqs", "createdb", cluster_source, db_path, "--threads", str(threads)], check=True)
    finally:
        if strip_gaps and os.path.isfile(cluster_source):
            # The DB now holds the sequences, so drop the scratch copy - in a
            # finally block so a failed run cannot leave an unaligned FASTA
            # behind for the generic `*.fasta` fallbacks to pick up later.
            os.remove(cluster_source)

    subprocess.run(cluster_cmd, check=True)
    tsv_output = os.path.join(mmseqs_dir, f"{base_name}_clusters.tsv")
    subprocess.run(["mmseqs", "createtsv", db_path, db_path, cluster_path, tsv_output, "--threads", str(threads)], check=True)

    if not strip_gaps:
        # segments_DB/<base>_cluster_seq.fasta is read by nothing in this repo -
        # only reachable via the generic `*.fasta` fallbacks - and when clustering
        # gap-stripped input it would be an unaligned copy of every sequence, the
        # exact thing those fallbacks must not find. Skip it (also saves several
        # GB per segment); --keep-gaps still produces it for backwards compat.
        subprocess.run(["mmseqs", "createseqfiledb", db_path, cluster_path, cluster_seq_path, "--threads", str(threads)], check=True)
        subprocess.run(["mmseqs", "result2flat", db_path, db_path, cluster_seq_path, f"{cluster_seq_path}.fasta"], check=True)

    rep_fasta = os.path.join(mmseqs_dir, f"{base_name}_cluster_rep.fasta")
    if strip_gaps:
        # convert2fasta would emit the ungapped sequences MMseqs clustered on,
        # which IQ-TREE rejects as ragged, so rebuild from the padded MSA.
        count = write_aligned_representatives(input_fasta, tsv_output, rep_fasta)
        print(f"[info] Wrote {count} aligned cluster representatives to {rep_fasta}")
    else:
        subprocess.run(["mmseqs", "createsubdb", cluster_path, db_path, cluster_rep_path], check=True)
        subprocess.run(["mmseqs", "convert2fasta", cluster_rep_path, rep_fasta], check=True)

    shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MMseqs2 clustering on FASTA files, with optional trimming.")
    parser.add_argument("-i", "--input_dir", help="Directory containing padded nucleotide alignments.", default="tmp/Pad-Alignment")
    parser.add_argument("-o", "--output_dir", help="Directory where outputs will be saved", default="tmp/MMseqClusters")
    parser.add_argument("--min-seq-id", type=float, default=0.95, help="Minimum sequence identity for clustering (default: 0.95)")
    parser.add_argument("--trim_cds_file", help="Optional: TSV file with ranges for trimming sequences before clustering. If omitted, no trimming is applied.", default=None)
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for MMseqs2 clustering (default: 8)")
    parser.add_argument("--keep-gaps", action="store_true",
                        help="Cluster the padded alignment as-is instead of stripping gaps first. MMseqs treats "
                             "'-' as sequence content, which inflates the cluster count; only use this to "
                             "reproduce pre-existing results.")
    parser.add_argument("--max-seqs", type=int, default=None,
                        help="Prefilter candidates kept per query (MMseqs default). Raise it if a large "
                             "near-identical group appears to be crowding true neighbours out of the candidate list.")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.trim_cds_file and os.path.isfile(args.trim_cds_file):
        print(f"[info] Trimming enabled. Using CDS file: {args.trim_cds_file}")
        trim_output_dir = os.path.join(args.output_dir, "trimmed_fastas")
        trimmed_files = apply_trimming(args.trim_cds_file, args.input_dir, trim_output_dir)
        if not trimmed_files:
            print("[warn] No trimmed outputs were produced; falling back to original FASTA files.")
        else:
            for original_filename, trimmed_fasta in trimmed_files.items():
                print(f"Clustering trimmed file: {trimmed_fasta}")
                run_mmseqs_clustering(trimmed_fasta, args.output_dir, args.min_seq_id, args.threads,
                                      strip_gaps=not args.keep_gaps, max_seqs=args.max_seqs)
            print("All processing completed.")
            exit(0)

    if args.trim_cds_file and not os.path.isfile(args.trim_cds_file):
        print(f"[warn] Provided --trim_cds_file does not exist: {args.trim_cds_file}. Proceeding without trimming.")

    # Default path: no trimming
    print("[info] Trimming not requested. Clustering original FASTA files.")
    fasta_files = [f for f in os.listdir(args.input_dir) if f.endswith(".fas") or f.endswith(".fasta")]
    if not fasta_files:
        print(f"No FASTA files found in {args.input_dir}")
    else:
        for fasta_file in fasta_files:
            input_fasta_path = os.path.join(args.input_dir, fasta_file)
            print(f"Clustering original file: {input_fasta_path}")
            run_mmseqs_clustering(input_fasta_path, args.output_dir, args.min_seq_id, args.threads,
                                  strip_gaps=not args.keep_gaps, max_seqs=args.max_seqs)

    print("All processing completed.")



