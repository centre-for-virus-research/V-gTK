import sqlite3
from pathlib import Path
from types import SimpleNamespace

from ExportUpdateAssets import main as export_update_assets_main


def _build_update_db(db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
        cur.executemany(
            "INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
            [
                ("REF1", "master", "1"),
                ("REF2", "reference", "2"),
                ("Q1", "query", "1"),
            ],
        )
        cur.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
        cur.executemany(
            "INSERT INTO sequence_alignment(primary_accession, alignment_name, alignment) VALUES (?, ?, ?)",
            [
                ("REF1", "REF1", "ATGC"),
                ("REF2", "REF2", "A-GC"),
            ],
        )
        cur.execute("CREATE TABLE trees (source TEXT, segment TEXT, newick TEXT)")
        cur.executemany(
            "INSERT INTO trees(source, segment, newick) VALUES (?, ?, ?)",
            [
                ("usher", "1", "(REF1:0.1,Q1:0.2);"),
                ("iqtree", "2", "(REF2:0.1);"),
            ],
        )
        cur.execute("CREATE TABLE excluded_accessions (primary_accession TEXT, reason TEXT)")
        cur.execute("INSERT INTO excluded_accessions(primary_accession, reason) VALUES ('Q1', 'x')")
        conn.commit()
    finally:
        conn.close()


def test_export_update_assets_emits_manifest_and_existing_ids(tmp_path: Path):
    db_path = tmp_path / "prev.db"
    out_dir = tmp_path / "UpdateAssets"
    _build_update_db(db_path)

    export_update_assets_main(SimpleNamespace(db=str(db_path), output_dir=str(out_dir)))

    assert (out_dir / "ref_backbones" / "refset_1_aln.fasta").exists()
    assert (out_dir / "ref_backbones" / "refset_2_aln.fasta").exists()
    assert (out_dir / "tree_manifest.tsv").exists()
    assert (out_dir / "existing_ids" / "segment_1_ids.txt").exists()

    ids_seg1 = (out_dir / "existing_ids" / "segment_1_ids.txt").read_text(encoding="utf-8").splitlines()
    assert "REF1" in ids_seg1
    assert "Q1" not in ids_seg1


def test_workflow_update_glue_contains_skip_and_guard():
    workflow_path = Path(__file__).resolve().parents[2] / "vgtk-init.nf"
    text = workflow_path.read_text(encoding="utf-8")

    assert "def UPDATE_MODE = params.update_db" in text
    assert "update DB path matches output DB path" in text
    assert "if( UPDATE_MODE )" in text and "USHER_UPDATE_PLACEMENT" in text
    assert "MMSEQS_CLUSTERING(cluster_input_ch)" in text
