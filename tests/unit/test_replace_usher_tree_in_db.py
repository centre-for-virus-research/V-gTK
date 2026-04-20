import sqlite3
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"

if str(SCRIPTS_DIR) not in sys.path:
	sys.path.insert(0, str(SCRIPTS_DIR))

from ReplaceUsherTreeInDb import replace_tree  # type: ignore[reportMissingImports]


def test_replace_tree_updates_only_selected_top_level_usher_row(tmp_path: Path):
	db_path = tmp_path / "test.db"
	replacement_tree = tmp_path / "replacement.nwk"
	replacement_tree.write_text("(A:0.1,B:0.2,C:0.3);\n", encoding="utf-8")

	conn = sqlite3.connect(str(db_path))
	try:
		conn.execute(
			"CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT)"
		)
		conn.executemany(
			"INSERT INTO trees(name, source, segment_key, segment, newick, created_at) VALUES (?, ?, ?, ?, ?, ?)",
			[
				("usher", "usher", None, None, "(OLD:0.1);", "2026-01-01 00:00:00"),
				(
					"usher_AF009606.aligned_merged_MSA_dedup",
					"usher",
					"AF009606.aligned_merged_MSA_dedup",
					None,
					"(SEGMENT_OLD:0.1);",
					"2026-01-01 00:00:00",
				),
			],
		)
		conn.commit()
	finally:
		conn.close()

	replace_tree(db_path=db_path, tree_path=replacement_tree)

	conn = sqlite3.connect(str(db_path))
	try:
		rows = conn.execute(
			"SELECT name, source, COALESCE(segment_key, ''), COALESCE(segment, ''), newick FROM trees ORDER BY name"
		).fetchall()
	finally:
		conn.close()

	assert rows == [
		(
			"usher",
			"usher",
			"",
			"",
			"(A:0.1,B:0.2,C:0.3);",
		),
		(
			"usher_AF009606.aligned_merged_MSA_dedup",
			"usher",
			"AF009606.aligned_merged_MSA_dedup",
			"",
			"(SEGMENT_OLD:0.1);",
		),
	]


def test_replace_tree_raises_when_selected_row_is_missing(tmp_path: Path):
	db_path = tmp_path / "test.db"
	replacement_tree = tmp_path / "replacement.nwk"
	replacement_tree.write_text("(A:0.1,B:0.2);\n", encoding="utf-8")

	conn = sqlite3.connect(str(db_path))
	try:
		conn.execute(
			"CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT)"
		)
		conn.execute(
			"INSERT INTO trees(name, source, segment_key, segment, newick, created_at) VALUES (?, ?, ?, ?, ?, ?)",
			("usher_other", "usher", "seg1", "1", "(OLD:0.1);", "2026-01-01 00:00:00"),
		)
		conn.commit()
	finally:
		conn.close()

	try:
		replace_tree(db_path=db_path, tree_path=replacement_tree)
		raise AssertionError("Expected replace_tree to raise when no matching row exists")
	except ValueError as exc:
		assert "Expected exactly one matching tree row" in str(exc)