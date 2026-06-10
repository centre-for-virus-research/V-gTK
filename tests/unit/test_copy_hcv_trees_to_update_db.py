import sqlite3
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEV_DIR = REPO_ROOT / "dev"

if str(DEV_DIR) not in sys.path:
	sys.path.insert(0, str(DEV_DIR))

from copy_hcv_trees_to_update_db import copy_trees_between_dbs  # type: ignore[reportMissingImports]


def _create_tree_db(db_path: Path, rows):
	conn = sqlite3.connect(str(db_path))
	try:
		conn.execute(
			"CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT)"
		)
		if rows:
			conn.executemany(
				"INSERT INTO trees(name, source, segment_key, segment, newick, created_at) VALUES (?, ?, ?, ?, ?, ?)",
				rows,
			)
		conn.commit()
	finally:
		conn.close()


def test_copy_trees_between_dbs_copies_source_rows_into_new_output(tmp_path: Path):
	source_db = tmp_path / "source.db"
	template_db = tmp_path / "template.db"
	output_db = tmp_path / "output.db"

	_create_tree_db(
		source_db,
		[
			("usher", "usher", "", "", "(A:0.1,B:0.2);", "2026-05-01 00:00:00"),
			("segment_tree", "usher", "segA", "A", "(C:0.3,D:0.4);", "2026-05-02 00:00:00"),
		],
	)
	_create_tree_db(
		template_db,
		[
			("placeholder", "other", "", "", "(OLD:1.0);", "2025-01-01 00:00:00"),
		],
	)

	returned_path, tree_count = copy_trees_between_dbs(source_db, template_db, output_db)

	assert returned_path == output_db
	assert tree_count == 2
	assert template_db.exists()

	conn = sqlite3.connect(str(output_db))
	try:
		rows = conn.execute(
			"SELECT name, source, COALESCE(segment_key, ''), COALESCE(segment, ''), newick, created_at FROM trees ORDER BY name"
		).fetchall()
	finally:
		conn.close()

	assert rows == [
		("segment_tree", "usher", "segA", "A", "(C:0.3,D:0.4);", "2026-05-02 00:00:00"),
		("usher", "usher", "", "", "(A:0.1,B:0.2);", "2026-05-01 00:00:00"),
	]


def test_copy_trees_between_dbs_raises_for_missing_source(tmp_path: Path):
	template_db = tmp_path / "template.db"
	_create_tree_db(template_db, [])

	try:
		copy_trees_between_dbs(tmp_path / "missing.db", template_db, tmp_path / "output.db")
		raise AssertionError("Expected FileNotFoundError for missing source DB")
	except FileNotFoundError as exc:
		assert "Source DB does not exist" in str(exc)