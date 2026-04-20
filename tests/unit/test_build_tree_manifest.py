# pyright: reportMissingImports=false
from pathlib import Path

from BuildTreeManifest import build_tree_manifest


def test_build_tree_manifest_collects_segmented_iqtree_and_usher_outputs(tmp_path: Path):
    iqtree_root = tmp_path / "iqtree_inputs"
    usher_root = tmp_path / "usher_inputs"
    manifest_path = tmp_path / "tree_manifest.tsv"

    iqtree_dir = iqtree_root / "IQTree_MMseqClusters_refset_1_aln_merged_MSA_dedup_cluster_input"
    usher_dir = usher_root / "Usher_MMseqClusters_refset_1_aln_merged_MSA_dedup_cluster_input"
    update_dir = usher_root / "UsherUpdate_refset_2_aln_merged_MSA_dedup_cluster_input"

    iqtree_dir.mkdir(parents=True)
    usher_dir.mkdir(parents=True)
    update_dir.mkdir(parents=True)

    (iqtree_dir / "iqtree.treefile").write_text("(A:0.1);\n", encoding="utf-8")
    (usher_dir / "uncondensed-final-tree.nh").write_text("(A:0.2);\n", encoding="utf-8")
    (update_dir / "final-tree.nh").write_text("(B:0.3);\n", encoding="utf-8")

    rows = build_tree_manifest(manifest_path, iqtree_dir=iqtree_root, usher_dir=usher_root)

    expected_rows = [
        {
            "source": "iqtree",
            "name": "iqtree_refset_1_aln_merged_MSA_dedup_cluster_input",
            "segment_key": "refset_1_aln_merged_MSA_dedup_cluster_input",
            "path": str(iqtree_dir / "iqtree.treefile"),
        },
        {
            "source": "usher",
            "name": "usher_refset_1_aln_merged_MSA_dedup_cluster_input",
            "segment_key": "refset_1_aln_merged_MSA_dedup_cluster_input",
            "path": str(usher_dir / "uncondensed-final-tree.nh"),
        },
        {
            "source": "usher",
            "name": "usher_refset_2_aln_merged_MSA_dedup_cluster_input",
            "segment_key": "refset_2_aln_merged_MSA_dedup_cluster_input",
            "path": str(update_dir / "final-tree.nh"),
        },
    ]

    assert sorted(rows, key=lambda row: (row["source"], row["segment_key"], row["path"])) == sorted(
        expected_rows,
        key=lambda row: (row["source"], row["segment_key"], row["path"]),
    )

    manifest_lines = manifest_path.read_text(encoding="utf-8").strip().splitlines()
    assert manifest_lines[0] == "source\tname\tsegment_key\tpath"
    assert len(manifest_lines) == 4


def test_build_tree_manifest_prefers_top_level_usher_tree_over_chunk_outputs(tmp_path: Path):
    usher_root = tmp_path / "usher_inputs"
    manifest_path = tmp_path / "tree_manifest.tsv"

    update_dir = usher_root / "UsherUpdate_refset_0_aln_merged_MSA_dedup"
    chunk_dir = update_dir / "chunk_0001"
    chunk_dir.mkdir(parents=True)

    (chunk_dir / "uncondensed-final-tree.nh").write_text("(CHUNK_A:0.1,CHUNK_B:0.2);\n", encoding="utf-8")
    (update_dir / "uncondensed-final-tree.nh").write_text("(FINAL_A:0.1,FINAL_B:0.2,FINAL_C:0.3);\n", encoding="utf-8")

    rows = build_tree_manifest(manifest_path, usher_dir=usher_root)

    assert rows == [
        {
            "source": "usher",
            "name": "usher_refset_0_aln_merged_MSA_dedup",
            "segment_key": "refset_0_aln_merged_MSA_dedup",
            "path": str(update_dir / "uncondensed-final-tree.nh"),
        }
    ]
