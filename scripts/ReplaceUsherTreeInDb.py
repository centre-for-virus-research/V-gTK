#!/usr/bin/env python3
import argparse
import sqlite3
from datetime import datetime
from pathlib import Path


def parse_args():
	parser = argparse.ArgumentParser(description="Replace or insert a stored tree row in a SQLite DB with a supplied Newick file.")
	parser.add_argument("--db", required=True, help="Path to the SQLite DB to update")
	parser.add_argument("--tree", required=True, help="Path to the replacement Newick file")
	parser.add_argument("--source", default="usher", help="Tree source (default: usher)")
	parser.add_argument("--name", default="usher", help="Tree name (default: usher)")
	parser.add_argument("--segment-key", default="", help="Optional segment_key selector for the tree row")
	parser.add_argument("--segment", default="", help="Optional segment selector for the tree row")
	parser.add_argument("--allow-insert", action="store_true", help="Insert a new tree row if none exists")
	return parser.parse_args()


def read_newick(tree_path):
	text = Path(tree_path).read_text(encoding="utf-8").strip()
	if not text:
		raise ValueError(f"Replacement tree file is empty: {tree_path}")
	return text


def replace_tree(db_path, tree_path, source="usher", name="usher", segment_key="", segment="", allow_insert=False):
	newick = read_newick(tree_path)
	conn = sqlite3.connect(str(db_path))
	try:
		cursor = conn.cursor()
		cursor.execute(
			"""
			SELECT COUNT(*)
			FROM trees
			WHERE COALESCE(source, '')=?
			  AND COALESCE(name, '')=?
			  AND COALESCE(segment_key, '')=?
			  AND COALESCE(segment, '')=?
			""",
			(source.strip(), name.strip(), segment_key.strip(), segment.strip()),
		)
		match_count = cursor.fetchone()[0]
		if match_count == 0 and allow_insert:
			cursor.execute(
				"""
				INSERT INTO trees (name, source, segment_key, segment, newick, created_at)
				VALUES (?, ?, ?, ?, ?, ?)
				""",
				(
					name.strip(),
					source.strip(),
					segment_key.strip() or None,
					segment.strip() or None,
					newick,
					datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
				),
			)
			conn.commit()
			return "inserted"
		elif match_count != 1:
			raise ValueError(
				"Expected exactly one matching tree row, "
			)

		cursor.execute(
			"""
			UPDATE trees
			SET newick=?, created_at=?
			WHERE COALESCE(source, '')=?
			  AND COALESCE(name, '')=?
			  AND COALESCE(segment_key, '')=?
			  AND COALESCE(segment, '')=?
			""",
			(
				newick,
				datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
				source.strip(),
				name.strip(),
				segment_key.strip(),
				segment.strip(),
			),
		)
		conn.commit()
		return "updated"
	finally:
		conn.close()


def main():
	args = parse_args()
	action = replace_tree(
		db_path=args.db,
		tree_path=args.tree,
		source=args.source,
		name=args.name,
		segment_key=args.segment_key,
		segment=args.segment,
		allow_insert=args.allow_insert,
	)
	print(
	)


if __name__ == "__main__":
	main()
