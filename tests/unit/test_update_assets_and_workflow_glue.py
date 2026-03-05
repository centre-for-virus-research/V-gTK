from pathlib import Path
from types import SimpleNamespace

from ExportUpdateAssets import main as export_update_assets_main

def test_export_update_assets_emits_manifest_and_existing_ids(tmp_path: Path, basic_update_db: Path):
    out_dir = tmp_path / "UpdateAssets"

    export_update_assets_main(SimpleNamespace(db=str(basic_update_db), output_dir=str(out_dir)))

    assert (out_dir / "ref_backbones" / "refset_1_aln.fasta").exists()
    assert (out_dir / "ref_backbones" / "refset_2_aln.fasta").exists()
    assert (out_dir / "tree_manifest.tsv").exists()
    assert (out_dir / "existing_ids" / "segment_1_ids.txt").exists()

    ids_seg1 = (out_dir / "existing_ids" / "segment_1_ids.txt").read_text(encoding="utf-8").splitlines()
    assert "REF1" in ids_seg1
    assert "Q_EXCL" not in ids_seg1


def test_export_update_assets_prefers_usher_tree_source(tmp_path: Path, basic_update_db: Path):
    out_dir = tmp_path / "UpdateAssets"
    export_update_assets_main(SimpleNamespace(db=str(basic_update_db), output_dir=str(out_dir)))

    tree_manifest = (out_dir / "tree_manifest.tsv").read_text(encoding="utf-8")
    assert "usher" in tree_manifest
    assert "segment_1.nwk" in tree_manifest


def test_export_update_assets_writes_tree_manifest_header_when_no_trees(tmp_path: Path, basic_update_db: Path):
    import sqlite3

    conn = sqlite3.connect(str(basic_update_db))
    try:
        conn.execute("DELETE FROM trees")
        conn.commit()
    finally:
        conn.close()

    out_dir = tmp_path / "UpdateAssets"
    export_update_assets_main(SimpleNamespace(db=str(basic_update_db), output_dir=str(out_dir)))
    content = (out_dir / "tree_manifest.tsv").read_text(encoding="utf-8").strip()
    assert content == "segment\tsource\tpath"


def test_workflow_update_glue_contains_skip_and_guard():
    workflow_path = Path(__file__).resolve().parents[2] / "vgtk-init.nf"
    text = workflow_path.read_text(encoding="utf-8")

    assert "def UPDATE_MODE = params.update_db" in text
    assert "update DB path matches output DB path" in text
    assert "if( UPDATE_MODE )" in text and "USHER_UPDATE_PLACEMENT" in text
    assert "MMSEQS_CLUSTERING(cluster_input_ch)" in text
