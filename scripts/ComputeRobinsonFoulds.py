import argparse
from io import StringIO

from Bio import Phylo


def load_tree(path):
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read().strip()
    if not text:
        raise ValueError(f"Empty tree file: {path}")
    return Phylo.read(StringIO(text), "newick")


def terminal_names(tree):
    return {t.name for t in tree.get_terminals() if t.name}


def bipartitions(tree, shared_taxa):
    all_leaves = set(shared_taxa)
    splits = set()
    if len(all_leaves) < 4:
        return splits

    for clade in tree.find_clades(order="postorder"):
        if clade is tree.root:
            continue
        leaves = {t.name for t in clade.get_terminals() if t.name in all_leaves}
        if len(leaves) < 2:
            continue
        complement = all_leaves - leaves
        if len(complement) < 2:
            continue
        canonical = frozenset(leaves) if len(leaves) <= len(complement) else frozenset(complement)
        splits.add(canonical)
    return splits


def compute_rf(tree_a, tree_b):
    taxa_a = terminal_names(tree_a)
    taxa_b = terminal_names(tree_b)
    shared = taxa_a & taxa_b
    if len(shared) < 4:
        return {
            "shared_taxa": len(shared),
            "splits_a": 0,
            "splits_b": 0,
            "rf": 0,
            "max_rf": 0,
            "normalized_rf": 0.0,
            "unique_to_a": len(taxa_a - shared),
            "unique_to_b": len(taxa_b - shared),
        }

    s_a = bipartitions(tree_a, shared)
    s_b = bipartitions(tree_b, shared)
    rf = len(s_a - s_b) + len(s_b - s_a)
    max_rf = max(0, 2 * (len(shared) - 3))
    norm = (rf / max_rf) if max_rf else 0.0

    return {
        "shared_taxa": len(shared),
        "splits_a": len(s_a),
        "splits_b": len(s_b),
        "rf": rf,
        "max_rf": max_rf,
        "normalized_rf": norm,
        "unique_to_a": len(taxa_a - shared),
        "unique_to_b": len(taxa_b - shared),
    }


def main():
    parser = argparse.ArgumentParser(description="Compute Robinson-Foulds distance between two Newick trees")
    parser.add_argument("--tree-a", required=True, help="Path to first Newick tree")
    parser.add_argument("--tree-b", required=True, help="Path to second Newick tree")
    parser.add_argument("--max-rf", type=int, default=None, help="Fail if RF distance is greater than this value")
    parser.add_argument(
        "--max-normalized-rf",
        type=float,
        default=None,
        help="Fail if normalized RF (RF/maxRF) is greater than this value",
    )
    args = parser.parse_args()

    tree_a = load_tree(args.tree_a)
    tree_b = load_tree(args.tree_b)
    metrics = compute_rf(tree_a, tree_b)

    print(f"shared_taxa={metrics['shared_taxa']}")
    print(f"splits_tree_a={metrics['splits_a']}")
    print(f"splits_tree_b={metrics['splits_b']}")
    print(f"rf_distance={metrics['rf']}")
    print(f"max_rf={metrics['max_rf']}")
    print(f"normalized_rf={metrics['normalized_rf']:.6f}")
    print(f"unique_taxa_tree_a={metrics['unique_to_a']}")
    print(f"unique_taxa_tree_b={metrics['unique_to_b']}")

    if args.max_rf is not None and metrics["rf"] > args.max_rf:
        raise SystemExit(f"RF distance {metrics['rf']} exceeds threshold {args.max_rf}")

    if args.max_normalized_rf is not None and metrics["normalized_rf"] > args.max_normalized_rf:
        raise SystemExit(
            f"Normalized RF {metrics['normalized_rf']:.6f} exceeds threshold {args.max_normalized_rf}"
        )


if __name__ == "__main__":
    main()