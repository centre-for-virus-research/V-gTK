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
    assert "def TREE_FREE_MODE = params.tree_free.toString().toBoolean()" in text
    assert "def BASE_TREE_ONLY_MODE = params.base_tree_only.toString().toBoolean()" in text
    assert "params.tree_free=true with params.update_db: skipping tree placement and retaining any existing trees already stored in the update DB" in text
    assert "params.base_tree_only=true: skipping USHER placement and storing IQ-TREE outputs in the DB" in text
    assert "update DB path matches output DB path" in text
    assert "if( UPDATE_MODE )" in text and "USHER_PLACEMENT(usher_input_ch)" in text
    assert "if( TREE_FREE_MODE )" in text
    assert "if( BASE_TREE_ONLY_MODE )" in text
    assert 'EXTRA_ARGS="${EXTRA_ARGS} --allow-no-trees"' in text
    assert 'EXTRA_ARGS="${EXTRA_ARGS} --segment-tree-source iqtree"' in text
    assert "MMSEQS_CLUSTERING(cluster_input_ch)" in text
    assert "REF_FASTA_FROM_UPDATE_DB" not in text
    assert "REF_LIST_FROM_UPDATE_DB" not in text
    assert "blast_ref_fasta = REF_FASTA_FROM_UPDATE_DB.out.ref_fasta_from_db" not in text
    assert "UPDATE_BLAST_ARGS=\"--update_db !{params.update_db}\"" in text
    assert "UsherPlacement.py" in text
    assert "process USHER_UPDATE_PLACEMENT" not in text
    assert "db_ref_backbones" not in text
    assert "process HOST_TAXA_TABLE" in text
    assert "HOST_TAXA_TABLE.out.host_taxa" in text
    assert "GENERATE_TABLES.out.host_taxa" not in text
    assert 'for USHER_DIR in usher_inputs/*; do' in text
    assert 'USHER_FILE=$(find -L usher_inputs -name "final-tree.nh" -print -quit || true)' not in text
    assert 'find -L mmseq_inputs -type f -name "*_clusters.tsv" -print0 | sort -z | xargs -0 -r cat > "$MERGED_CLUSTER_TSV"' in text
    assert 'path "Taxa/names.dmp", emit: taxa_names' in text
    assert 'path "Taxa/nodes.dmp", emit: taxa_nodes' in text
    assert 'if [ "!{params.mutation_catalog}" != "null" ] && [ -n "!{params.mutation_catalog}" ]; then' in text
    assert "AnnotateMutations.py" in text
    assert "--mutation_catalog !{params.mutation_catalog}" in text
    assert "--catalog_column_profile ${CATALOG_PROFILE}" in text
    assert "--db !{params.db_name}.db" in text
    # TEST_SUBSAMPLE_CLUSTER_INPUT is invoked in test mode to speed up CI
    assert "TEST_SUBSAMPLE_CLUSTER_INPUT(" in text


def test_workflow_genbank_parser_passes_segmented_flag():
    workflow_path = Path(__file__).resolve().parents[2] / "vgtk-init.nf"
    text = workflow_path.read_text(encoding="utf-8")

    assert "GenBankParser.py -r !{ref_list_path} -d !{gen_bank_XML} -o . -b . -s !{params.is_segmented} ${extra}" in text


def test_nextflow_config_defines_mutation_defaults():
    config_path = Path(__file__).resolve().parents[2] / "nextflow.config"
    text = config_path.read_text(encoding="utf-8")

    assert "mutation_catalog  = null" in text
    assert "mutation_virus    = null" in text
    assert "tree_free       = false" in text
    assert "base_tree_only  = false" in text
    assert "segmented_base_tree" in text
