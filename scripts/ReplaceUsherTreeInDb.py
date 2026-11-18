#!/usr/bin/env python3
import argparse
import sqlite3
from datetime import datetime
from pathlib import Path

import usher_assets


def parse_args():
	parser = argparse.ArgumentParser(description="Replace or insert a stored tree row in a SQLite DB with a supplied Newick file.")
	parser.add_argument("--db", required=True, help="Path to the SQLite DB to update")
	parser.add_argument("--tree", required=True, help="Path to the replacement Newick file")
	parser.add_argument("--source", default="usher", help="Tree source (default: usher)")
	parser.add_argument("--name", default="usher", help="Tree name (default: usher)")
	parser.add_argument("--segment-key", default="", help="Optional segment_key selector for the tree row")
	parser.add_argument("--segment", default="", help="Optional segment selector for the tree row")
	parser.add_argument("--allow-insert", action="store_true", help="Insert a new tree row if none exists")
	parser.add_argument("--usher-dir", default=None, help="UsherPlacement output dir whose usher.pb and sidecars belong to --tree; stored with the row")
	parser.add_argument("--keep-usher-assets", action="store_true", help="The new tree is the old one rewritten (e.g. rerooted): keep its stored usher.pb if the tip set is unchanged")
	return parser.parse_args()


def read_newick(tree_path):
	text = Path(tree_path).read_text(encoding="utf-8").strip()
	if not text:
		raise ValueError(f"Replacement tree file is empty: {tree_path}")
	return text


def _sync_usher_bundle(conn, key, tree_path, old_newick, newick, created_at, usher_dir, keep_usher_assets):
	"""Keep usher_assets pointing at the tree now in the row.

	A supplied --usher-dir replaces the bundle (it must verify against --tree);
	otherwise the old bundle is relinked or marked stale, never left claiming
	to describe the new tree.
	"""
	if usher_dir:
		bundle_tree = Path(usher_dir) / Path(tree_path).name
		bundle = usher_assets.collect_bundle(str(bundle_tree)) if bundle_tree.is_file() else None
		if bundle is None or bundle_tree.read_bytes().strip() != Path(tree_path).read_bytes().strip():
			raise ValueError(f"{usher_dir} holds no usher.pb verified against {tree_path}")
		usher_assets.store_bundle(conn, key, bundle, newick, created_at)
		print(f"[ReplaceUsherTreeInDb] Stored {bundle['status']} UShER bundle from {usher_dir}")
		return
	status = usher_assets.retarget_bundle(conn, key, old_newick, newick, created_at, keep_usher_assets)
	if status == "stale":
		print("[ReplaceUsherTreeInDb] Warning: stored usher.pb no longer matches this tree; marked stale (pass --usher-dir to replace it)")
	elif status:
		print(f"[ReplaceUsherTreeInDb] Stored usher.pb {status} to the new tree")


def replace_tree(db_path, tree_path, source="usher", name="usher", segment_key="", segment="", allow_insert=False, usher_dir=None, keep_usher_assets=False):
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
		key = {"name": name, "source": source, "segment_key": segment_key, "segment": segment}
		created_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
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
					created_at,
				),
			)
			_sync_usher_bundle(conn, key, tree_path, None, newick, created_at, usher_dir, keep_usher_assets)
			conn.commit()
			return "inserted"
		elif match_count != 1:
			selector = (
				f"source={source.strip()!r}, name={name.strip()!r}, "
				f"segment_key={segment_key.strip()!r}, segment={segment.strip()!r}"
			)
			hint = (
				" Pass --allow-insert to create it."
				if match_count == 0
				else " Narrow the selector so it matches one row."
			)
			raise ValueError(
				f"Expected exactly one matching tree row, found {match_count} for {selector}.{hint}"
			)

		old_newick = cursor.execute(
			"""
			SELECT newick FROM trees
			WHERE COALESCE(source, '')=?
			  AND COALESCE(name, '')=?
			  AND COALESCE(segment_key, '')=?
			  AND COALESCE(segment, '')=?
			""",
			(source.strip(), name.strip(), segment_key.strip(), segment.strip()),
		).fetchone()[0]
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
				created_at,
				source.strip(),
				name.strip(),
				segment_key.strip(),
				segment.strip(),
			),
		)
		_sync_usher_bundle(conn, key, tree_path, old_newick, newick, created_at, usher_dir, keep_usher_assets)
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
		usher_dir=args.usher_dir,
		keep_usher_assets=args.keep_usher_assets,
	)
	print(
		f"{action} tree row: source={args.source!r}, name={args.name!r}, "
		f"segment={args.segment!r}, segment_key={args.segment_key!r} "
		f"<- {args.tree}"
	)


if __name__ == "__main__":
	main()
