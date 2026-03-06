from pathlib import Path

from UsherPlacement import UsherPlacement


def test_prepare_update_assets_exports_tree_and_existing_ids(tmp_path: Path, basic_update_db: Path):
	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(">REF1\nATGC\n>Q_NEW\nATGT\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		update_db=str(basic_update_db),
	)

	tree_path, ids_path = processor.prepare_update_assets()

	assert Path(tree_path).exists()
	assert Path(ids_path).exists()
	assert "REF1" in Path(ids_path).read_text(encoding="utf-8")
	assert "Q_EXCL" not in Path(ids_path).read_text(encoding="utf-8")
	assert "Q_OLD" in Path(tree_path).read_text(encoding="utf-8")