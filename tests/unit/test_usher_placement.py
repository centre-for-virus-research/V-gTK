import sqlite3
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


def test_update_mode_uses_tree_terminals_as_existing_ids_and_rehydrates_missing_db_alignments(tmp_path: Path):
	db_path = tmp_path / "update.db"
	conn = sqlite3.connect(str(db_path))
	try:
		cur = conn.cursor()
		cur.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)")
		cur.executemany(
			"INSERT INTO meta_data(primary_accession, accession_type, segment) VALUES (?, ?, ?)",
			[("REF1", "master", "1"), ("Q_OLD", "query", "1"), ("Q_NEW", "query", "1")],
		)
		cur.execute(
			"CREATE TABLE sequence_alignment (primary_accession TEXT, sequence_id TEXT, alignment_name TEXT, alignment TEXT, segment TEXT)"
		)
		cur.executemany(
			"INSERT INTO sequence_alignment(primary_accession, sequence_id, alignment_name, alignment, segment) VALUES (?, ?, ?, ?, ?)",
			[("REF1", "REF1", "REF1", "ATGC", "1"), ("Q_OLD", "Q_OLD", "REF1", "ATGT", "1")],
		)
		cur.execute("CREATE TABLE trees (name TEXT, source TEXT, segment TEXT, segment_key TEXT, newick TEXT)")
		cur.execute(
			"INSERT INTO trees(name, source, segment, segment_key, newick) VALUES (?, ?, ?, ?, ?)",
			("usher_seg1", "usher", "1", "seg_1", "(REF1:0.1);"),
		)
		conn.commit()
	finally:
		conn.close()

	msa = tmp_path / "segment_1_dedup.fasta"
	msa.write_text(">REF1\nATGC\n>Q_NEW\nATGA\n", encoding="utf-8")

	processor = UsherPlacement(
		padded_aln=str(msa),
		output_dir=str(tmp_path / "usher_out"),
		update_db=str(db_path),
	)

	tree_path, ids_path = processor.prepare_update_assets()
	assert Path(ids_path).read_text(encoding="utf-8").strip().splitlines() == ["REF1"]

	merged_path = processor.build_update_alignment_input(tree_path)
	merged_text = Path(merged_path).read_text(encoding="utf-8")
	assert ">Q_NEW" in merged_text
	assert ">Q_OLD" in merged_text