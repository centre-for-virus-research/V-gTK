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


INFORMATIVE_BASES = ("A", "C", "G", "T", "U")


def informative_length(sequence):
    """Unambiguous bases: everything that is not a gap, an N, or an IUPAC
    ambiguity code. str.count runs in C, so this stays cheap over a full
    segment alignment."""
    seq = sequence.upper()
    return sum(seq.count(base) for base in INFORMATIVE_BASES)


def sort_fasta_by_informative_length(input_fasta, output_fasta):
    """Rewrite a FASTA in descending order of informative (non-N, non-gap) length.

    --cluster-mode 2 is greedy incremental: it sorts by length and makes the
    longest sequence of each group the representative. But MMseqs counts N
    toward length (verified: ACGTNNNNNNACGT is length 14), so a half-N sequence
    outranks a shorter clean one and becomes the representative that ends up as
    a tree tip. Ordering the input by informative length, with createdb
    --shuffle 0, makes the tie-break prefer the highest-quality sequence
    instead of an arbitrary one.

    Sorting is done through an external sort on a temporary TSV so memory stays
    flat regardless of segment size.
    """
    tmp_tsv = output_fasta + ".sort.tsv"
    sorted_tsv = output_fasta + ".sorted.tsv"
    try:
        with open(tmp_tsv, "w") as handle:
            for record in SeqIO.parse(input_fasta, "fasta"):
                sequence = str(record.seq)
                handle.write(f"{informative_length(sequence)}\t{record.id}\t{sequence}\n")

        env = dict(os.environ, LC_ALL="C")
        with open(sorted_tsv, "w") as handle:
            subprocess.run(["sort", "-t", "\t", "-k1,1nr", "-s", tmp_tsv],
                           check=True, stdout=handle, env=env)

        # Write to a scratch path and swap it in only on success. This is
        # normally called with output_fasta == input_fasta, so writing in place
        # would truncate the input if the sort failed.
        staged = output_fasta + ".ordered"
        written = 0
        with open(sorted_tsv) as src, open(staged, "w") as out:
            for line in src:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                out.write(f">{parts[1]}\n{parts[2]}\n")
                written += 1
        if written == 0:
            os.remove(staged)
            raise ValueError(f"informative-length sort produced no sequences from {input_fasta}")
        os.replace(staged, output_fasta)
        return written
    finally:
        for path in (tmp_tsv, sorted_tsv):
            if os.path.isfile(path):
                os.remove(path)


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


def split_by_completeness(alignment_fasta, min_completeness, complete_out, remainder_out):
    """Split the alignment into near-complete sequences and everything else.

    Completeness is informative bases / alignment width, so it measures how much
    of the reference-coordinate alignment a sequence actually covers. Both
    outputs are written unaligned, ready for MMseqs.

    The point is to keep fragments out of the backbone. Measured on influenza HA:
    a one-pass clustering put 35% of its representatives at >25% gaps, because
    every fragment without a same-threshold partner becomes its own cluster and
    therefore its own tree tip. Clustering only the >=90% complete sequences
    gives 0% such tips while discarding ~9% of the input - and step 2 hands
    those discarded sequences back a cluster assignment.
    """
    translation = str.maketrans("", "", GAP_CHARACTERS)
    n_complete = n_remainder = 0
    with open(complete_out, "w") as complete_handle, open(remainder_out, "w") as remainder_handle:
        for record in SeqIO.parse(alignment_fasta, "fasta"):
            aligned = str(record.seq)
            width = len(aligned)
            if width == 0:
                continue
            sequence = aligned.translate(translation)
            if not sequence:
                continue
            if informative_length(aligned) / width >= min_completeness:
                complete_handle.write(f">{record.id}\n{sequence}\n")
                n_complete += 1
            else:
                remainder_handle.write(f">{record.id}\n{sequence}\n")
                n_remainder += 1
    return n_complete, n_remainder


def assign_remainder_to_representatives(remainder_fasta, representatives_fasta, tmp_dir,
                                        threads, min_seq_id):
    """Search the held-back sequences against the step-1 representatives and
    return {member_id: representative_id} for those meeting min_seq_id.

    Assigning members to EXISTING representatives (rather than letting them form
    new clusters) keeps the set of cluster centroids identical to the set of
    tree tips, which ValidateDbTree relies on.
    """
    os.makedirs(tmp_dir, exist_ok=True)
    hits_path = os.path.join(tmp_dir, "remainder_hits.m8")

    # The representative FASTA is aligned (it feeds the tree), but searching
    # against gapped sequences reintroduces the gaps-as-X problem this whole
    # design exists to avoid - measured cost on HA: 482/707 assignments instead
    # of 696/707. Search against an unaligned copy.
    ungapped_reps = os.path.join(tmp_dir, "representatives_ungapped.seq")
    strip_alignment_gaps(representatives_fasta, ungapped_reps)

    subprocess.run([
        "mmseqs", "easy-search", remainder_fasta, ungapped_reps, hits_path,
        os.path.join(tmp_dir, "search_tmp"),
        "--threads", str(threads),
        # REQUIRED for nucleotide input: without it mmseqs silently produces no
        # output at all rather than failing loudly.
        "--search-type", "3",
        "--cov-mode", "2", "-c", "0.8",
        "--max-seqs", "10",
        "--format-output", "query,target,pident",
    ], check=True)

    threshold = float(min_seq_id) * 100.0
    best = {}
    if os.path.isfile(hits_path):
        with open(hits_path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                query, target = parts[0], parts[1]
                try:
                    pident = float(parts[2])
                except ValueError:
                    continue
                if pident < threshold:
                    continue
                if query not in best or pident > best[query][1]:
                    best[query] = (target, pident)
    return {query: target for query, (target, _) in best.items()}


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


def run_mmseqs_clustering(input_fasta, output_dir, min_seq_id, threads=8, strip_gaps=True, max_seqs=None,
                          sort_by_quality=True, two_step=False, min_completeness=0.9):
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
    remainder_source = None
    if two_step:
        # Step 1 clusters only the near-complete sequences, so fragments cannot
        # become representatives (and therefore tree tips). Step 2 gives them a
        # cluster assignment afterwards.
        cluster_source = os.path.join(segments_db_dir, f"{base_name}_complete.seq")
        remainder_source = os.path.join(segments_db_dir, f"{base_name}_remainder.seq")
        n_complete, n_remainder = split_by_completeness(
            input_fasta, min_completeness, cluster_source, remainder_source)
        print(f"[info] Two-step clustering: {n_complete} sequences >= {min_completeness:.0%} complete "
              f"form the backbone, {n_remainder} held back for assignment")
        if n_complete == 0:
            raise ValueError(
                f"No sequences in {input_fasta} reach {min_completeness:.0%} completeness; "
                f"lower --min-completeness or disable --two-step")
        if n_remainder == 0:
            remainder_source = None
    elif strip_gaps:
        cluster_source = os.path.join(segments_db_dir, f"{base_name}_ungapped.seq")
        count = strip_alignment_gaps(input_fasta, cluster_source)
        print(f"[info] Clustering {count} gap-stripped sequences (representatives stay aligned)")
    if (strip_gaps or two_step) and sort_by_quality:
        count = sort_fasta_by_informative_length(cluster_source, cluster_source)
        print(f"[info] Ordered {count} sequences by descending informative length so the "
              f"greedy clustering prefers complete sequences as representatives")

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
        # No --threads: createdb is IO-bound and mmseqs <= v14 rejects the flag
        # outright ("Unrecognized parameter"), so passing it breaks the run on any
        # host where an older mmseqs is first on PATH.
        createdb_cmd = ["mmseqs", "createdb", cluster_source, db_path]
        if (strip_gaps or two_step) and sort_by_quality:
            # createdb shuffles by default (--shuffle 1), which would discard the
            # informative-length ordering established above.
            createdb_cmd += ["--shuffle", "0"]
        subprocess.run(createdb_cmd, check=True)
    finally:
        if (strip_gaps or two_step) and os.path.isfile(cluster_source):
            # The DB now holds the sequences, so drop the scratch copy - in a
            # finally block so a failed run cannot leave an unaligned FASTA
            # behind for the generic `*.fasta` fallbacks to pick up later.
            os.remove(cluster_source)

    subprocess.run(cluster_cmd, check=True)
    tsv_output = os.path.join(mmseqs_dir, f"{base_name}_clusters.tsv")
    subprocess.run(["mmseqs", "createtsv", db_path, db_path, cluster_path, tsv_output, "--threads", str(threads)], check=True)

    if not strip_gaps and not two_step:
        # segments_DB/<base>_cluster_seq.fasta is read by nothing in this repo -
        # only reachable via the generic `*.fasta` fallbacks - and when clustering
        # gap-stripped input it would be an unaligned copy of every sequence, the
        # exact thing those fallbacks must not find. Skip it (also saves several
        # GB per segment); --keep-gaps still produces it for backwards compat.
        subprocess.run(["mmseqs", "createseqfiledb", db_path, cluster_path, cluster_seq_path, "--threads", str(threads)], check=True)
        subprocess.run(["mmseqs", "result2flat", db_path, db_path, cluster_seq_path, f"{cluster_seq_path}.fasta"], check=True)

    rep_fasta = os.path.join(mmseqs_dir, f"{base_name}_cluster_rep.fasta")
    if two_step:
        # Representatives come from step 1 only, so the tree gets complete
        # sequences. Write them before step 2 so the search has a target DB.
        count = write_aligned_representatives(input_fasta, tsv_output, rep_fasta)
        print(f"[info] Backbone: {count} aligned representatives from complete sequences")

        if remainder_source:
            # Deliberately the same identity as step 1. Step 2 assigns fragments to
            # clusters step 1 already defined; a looser threshold here would label a
            # sequence as belonging to a min_seq_id cluster on a weaker match, so the
            # cluster column would mean two different things. Sequences that match
            # nothing get a NULL cluster, which is the honest answer - they stay in
            # the database and are still placed on the tree by UShER.
            assignments = assign_remainder_to_representatives(
                remainder_source, rep_fasta, os.path.join(mmseqs_dir, "step2_tmp"),
                threads, min_seq_id)
            with open(tsv_output, "a") as handle:
                for member, representative in sorted(assignments.items()):
                    handle.write(f"{representative}\t{member}\n")
            held = sum(1 for _ in SeqIO.parse(remainder_source, "fasta"))
            print(f"[info] Step 2: assigned {len(assignments)}/{held} held-back sequences to a "
                  f"backbone cluster at >= {float(min_seq_id):.0%} identity; "
                  f"{held - len(assignments)} left unassigned (NULL cluster)")
            shutil.rmtree(os.path.join(mmseqs_dir, "step2_tmp"), ignore_errors=True)
        for path in (remainder_source,):
            if path and os.path.isfile(path):
                os.remove(path)
    elif strip_gaps:
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
    parser.add_argument("--two-step", action="store_true",
                        help="Two-step clustering. Step 1 clusters only sequences at least "
                             "--min-completeness covered, so fragments never become representatives "
                             "(and never become tree tips). Step 2 searches the held-back sequences "
                             "against those representatives and assigns each to its best match, so "
                             "they keep a cluster label without polluting the backbone.")
    parser.add_argument("--min-completeness", type=float, default=0.9,
                        help="Fraction of the alignment width a sequence must cover with unambiguous "
                             "bases to join step 1 (default: 0.9)")
    parser.add_argument("--no-quality-sort", action="store_true",
                        help="Do not order clustering input by informative (non-N, non-gap) length. "
                             "By default the input is sorted longest-informative-first and createdb "
                             "--shuffle is disabled, so --cluster-mode 2 picks the most complete "
                             "sequence of each cluster as its representative rather than an N-padded one.")
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
                                      strip_gaps=not args.keep_gaps, max_seqs=args.max_seqs,
                                      sort_by_quality=not args.no_quality_sort,
                                      two_step=args.two_step, min_completeness=args.min_completeness)
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
                                  strip_gaps=not args.keep_gaps, max_seqs=args.max_seqs,
                                      sort_by_quality=not args.no_quality_sort,
                                      two_step=args.two_step, min_completeness=args.min_completeness)

    print("All processing completed.")



