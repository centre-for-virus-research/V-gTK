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

def test_replace_tree_allows_insert_when_row_missing(tmp_path: Path):
	db_path = tmp_path / "test.db"
	tree_file = tmp_path / "new_tree.nwk"
	tree_file.write_text("(A:0.1,B:0.2);\n")

	conn = sqlite3.connect(str(db_path))
	try:
		conn.execute(
			"CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT)"
		)
		conn.commit()
	finally:
		conn.close()

	action = replace_tree(db_path=db_path, tree_path=tree_file, source="usher", name="usher_seg4", segment_key="seg4_key", segment="4", allow_insert=True)
	assert action == "inserted"

	conn = sqlite3.connect(str(db_path))
	try:
		rows = conn.execute("SELECT name, source, segment_key, segment, newick FROM trees").fetchall()
	finally:
		conn.close()

	assert rows == [("usher_seg4", "usher", "seg4_key", "4", "(A:0.1,B:0.2);")]


def _trees_db(path: Path, rows=()):
	conn = sqlite3.connect(str(path))
	try:
		conn.execute(
			"CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT)"
		)
		conn.executemany(
			"INSERT INTO trees(name, source, segment_key, segment, newick, created_at) VALUES (?, ?, ?, ?, ?, ?)",
			rows,
		)
		conn.commit()
	finally:
		conn.close()


def test_missing_row_error_names_the_selector_and_the_way_out(tmp_path: Path):
	"""The message was truncated mid-sentence - "Expected exactly one matching
	tree row, " and nothing more. Finalising the influenza DB means inserting
	eight UShER rows into a table that holds none, so the one error you are
	actually going to hit has to say which selector missed and that
	--allow-insert is the answer."""
	db_path = tmp_path / "test.db"
	tree_file = tmp_path / "t.nwk"
	tree_file.write_text("(A:0.1,B:0.2);\n", encoding="utf-8")
	_trees_db(db_path)

	try:
		replace_tree(db_path=db_path, tree_path=tree_file, source="usher",
					 name="usher_refset_4", segment_key="refset_4", segment="4")
		raise AssertionError("expected a ValueError")
	except ValueError as exc:
		message = str(exc)
		assert "found 0" in message
		assert "usher_refset_4" in message
		assert "--allow-insert" in message


def test_ambiguous_selector_error_does_not_suggest_inserting(tmp_path: Path):
	"""Two matches is the opposite problem - inserting a third would make it
	worse, so the hint has to differ."""
	db_path = tmp_path / "test.db"
	tree_file = tmp_path / "t.nwk"
	tree_file.write_text("(A:0.1,B:0.2);\n", encoding="utf-8")
	_trees_db(db_path, [
		("usher_seg4", "usher", "refset_4", "4", "(OLD:0.1);", "2026-01-01 00:00:00"),
		("usher_seg4", "usher", "refset_4", "4", "(OLD:0.2);", "2026-01-01 00:00:00"),
	])

	try:
		replace_tree(db_path=db_path, tree_path=tree_file, source="usher",
					 name="usher_seg4", segment_key="refset_4", segment="4")
		raise AssertionError("expected a ValueError")
	except ValueError as exc:
		message = str(exc)
		assert "found 2" in message
		assert "--allow-insert" not in message


def test_allow_insert_still_updates_a_row_that_does_exist(tmp_path: Path):
	"""--allow-insert must not mean "always insert" - running the finalisation
	twice would otherwise leave two rows per segment."""
	db_path = tmp_path / "test.db"
	tree_file = tmp_path / "t.nwk"
	tree_file.write_text("(A:0.1,B:0.9);\n", encoding="utf-8")
	_trees_db(db_path, [
		("usher_seg4", "usher", "refset_4", "4", "(OLD:0.1);", "2026-01-01 00:00:00"),
	])

	action = replace_tree(db_path=db_path, tree_path=tree_file, source="usher",
						  name="usher_seg4", segment_key="refset_4", segment="4",
						  allow_insert=True)
	assert action == "updated"

	conn = sqlite3.connect(str(db_path))
	try:
		rows = conn.execute("SELECT newick FROM trees").fetchall()
	finally:
		conn.close()
	assert rows == [("(A:0.1,B:0.9);",)]
