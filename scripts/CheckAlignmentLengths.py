import argparse
from collections import Counter

from Bio import SeqIO


def collect_lengths(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in: {fasta_file}")
    lengths = {record.id: len(str(record.seq)) for record in records}
    return lengths


def check_uniform_lengths(fasta_file):
    lengths = collect_lengths(fasta_file)
    counts = Counter(lengths.values())
    unique_lengths = sorted(counts.keys())
    min_len = min(unique_lengths)
    max_len = max(unique_lengths)
    is_uniform = len(unique_lengths) == 1

    non_majority = []
    if not is_uniform:
        majority_len = counts.most_common(1)[0][0]
        non_majority = [seq_id for seq_id, length in lengths.items() if length != majority_len]

    return {
        "is_uniform": is_uniform,
        "count": len(lengths),
        "min_len": min_len,
        "max_len": max_len,
        "unique_lengths": unique_lengths,
        "non_majority_ids": non_majority,
    }


def main():
    parser = argparse.ArgumentParser(description="Check that all sequences in a FASTA alignment have the same length.")
    parser.add_argument("-i", "--input_fasta", required=True, help="Path to input FASTA alignment file")
    parser.add_argument("--show_examples", type=int, default=10, help="How many non-majority IDs to print when mismatch exists")
    args = parser.parse_args()

    result = check_uniform_lengths(args.input_fasta)

    print(f"records={result['count']}")
    print(f"min_len={result['min_len']}")
    print(f"max_len={result['max_len']}")
    print(f"unique_lengths={','.join(str(x) for x in result['unique_lengths'])}")

    if result["is_uniform"]:
        print("status=PASS all sequences are the same length")
        return

    print("status=FAIL sequence lengths are not uniform")
    if result["non_majority_ids"]:
        example_ids = result["non_majority_ids"][: args.show_examples]
        print("example_non_majority_ids=" + ",".join(example_ids))
    raise SystemExit(1)


if __name__ == "__main__":
    main()
