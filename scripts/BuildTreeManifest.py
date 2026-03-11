#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Union


SOURCE_CONFIG = {
    "iqtree": {
        "patterns": ("*.treefile",),
        "prefixes": (r"^IQTree_MMseqClusters_",),
    },
    "usher": {
        "patterns": ("final-tree.nh", "uncondensed-final-tree.nh"),
        "prefixes": (r"^Usher_MMseqClusters_", r"^UsherUpdate_"),
    },
}


def _strip_source_prefixes(name: str, source: str) -> str:
    value = name
    for pattern in SOURCE_CONFIG[source]["prefixes"]:
        value = re.sub(pattern, "", value)
    return value


def _find_first_tree_file(directory: Path, source: str) -> Optional[Path]:
    for pattern in SOURCE_CONFIG[source]["patterns"]:
        matches = sorted(directory.rglob(pattern))
        if matches:
            return matches[0]
    return None


def collect_tree_rows(root_dir: Optional[Union[str, Path]], source: str) -> List[Dict[str, str]]:
    if not root_dir:
        return []

    root = Path(root_dir)
    if not root.is_dir():
        return []

    candidate_dirs = sorted(child for child in root.iterdir() if child.is_dir())
    if not candidate_dirs:
        candidate_dirs = [root]

    rows = []
    for candidate in candidate_dirs:
        tree_path = _find_first_tree_file(candidate, source)
        if tree_path is None:
            continue

        segment_key = ""
        if candidate != root:
            segment_key = _strip_source_prefixes(candidate.name, source)

        rows.append(
            {
                "source": source,
                "name": f"{source}_{segment_key or 'tree'}",
                "segment_key": segment_key,
                "path": str(tree_path),
            }
        )

    return rows


def build_tree_manifest(
    output_path: Union[str, Path],
    iqtree_dir: Optional[Union[str, Path]] = None,
    usher_dir: Optional[Union[str, Path]] = None,
) -> List[Dict[str, str]]:
    rows = []
    rows.extend(collect_tree_rows(iqtree_dir, "iqtree"))
    rows.extend(collect_tree_rows(usher_dir, "usher"))

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "name", "segment_key", "path"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a tree manifest from staged IQ-TREE and UShER output directories.")
    parser.add_argument("--output", required=True, help="Output TSV manifest path")
    parser.add_argument("--iqtree-dir", default=None, help="Root directory containing staged IQ-TREE outputs")
    parser.add_argument("--usher-dir", default=None, help="Root directory containing staged UShER outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_tree_manifest(args.output, iqtree_dir=args.iqtree_dir, usher_dir=args.usher_dir)


if __name__ == "__main__":
    main()
