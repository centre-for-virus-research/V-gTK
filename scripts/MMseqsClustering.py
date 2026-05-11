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


def run_mmseqs_clustering(input_fasta, output_dir, min_seq_id, threads=8):
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

    subprocess.run(["mmseqs", "createdb", input_fasta, db_path], check=True)
    subprocess.run([
        "mmseqs", "cluster",
        "--min-seq-id", str(min_seq_id),
        db_path, cluster_path, tmp_dir,
        "--threads", str(threads)
    ], check=True)
    tsv_output = os.path.join(mmseqs_dir, f"{base_name}_clusters.tsv")
    subprocess.run(["mmseqs", "createtsv", db_path, db_path, cluster_path, tsv_output], check=True)

    subprocess.run(["mmseqs", "createseqfiledb", db_path, cluster_path, cluster_seq_path], check=True)
    subprocess.run(["mmseqs", "result2flat", db_path, db_path, cluster_seq_path, f"{cluster_seq_path}.fasta"], check=True)

    subprocess.run(["mmseqs", "createsubdb", cluster_path, db_path, cluster_rep_path], check=True)
    subprocess.run(["mmseqs", "convert2fasta", cluster_rep_path, os.path.join(mmseqs_dir, f"{base_name}_cluster_rep.fasta")], check=True)

    shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MMseqs2 clustering on FASTA files, with optional trimming.")
    parser.add_argument("-i", "--input_dir", help="Directory containing padded nucleotide alignments.", default="tmp/Pad-Alignment")
    parser.add_argument("-o", "--output_dir", help="Directory where outputs will be saved", default="tmp/MMseqClusters")
    parser.add_argument("--min-seq-id", type=float, default=0.95, help="Minimum sequence identity for clustering (default: 0.95)")
    parser.add_argument("--trim_cds_file", help="Optional: TSV file with ranges for trimming sequences before clustering. If omitted, no trimming is applied.", default=None)
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for MMseqs2 clustering (default: 8)")

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
                run_mmseqs_clustering(trimmed_fasta, args.output_dir, args.min_seq_id, args.threads)
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
            run_mmseqs_clustering(input_fasta_path, args.output_dir, args.min_seq_id, args.threads)

    print("All processing completed.")



