#!/usr/bin/env python3
"""The UShER mutation-annotated tree (usher.pb) and its sidecars, kept in the
DB beside the trees row they describe.

A trees row holds only a Newick string, so the mutations behind an UShER tree
used to live in a loose usher.pb that nothing tied to the row: rewriting the
row (ReplaceUsherTreeInDb, an update run) left whatever .pb was lying around
silently describing a different tree. Here every usher_assets row carries the
same (name, source, segment_key, segment) key as its trees row, plus the
sha256 of that row's Newick and its created_at, so a reader can always tell
whether the bundle still matches the tree.

UsherPlacement writes usher_provenance.json next to its outputs, recording the
sha256 of the .pb and of every tree file it produced. collect_bundle only
accepts a .pb whose provenance names the exact tree file being stored.

status values:
	verified    provenance confirms the .pb was produced with this tree file
	unverified  no provenance (output from before it existed); the .pb was
	            found next to the tree but cannot be proven to match it
	relinked    the trees row was rewritten without changing its tip set
	            (e.g. midpoint rerooting); mutations still hold, the .pb root
	            may differ from the stored Newick
	stale       the trees row was replaced by a tree the .pb does not describe

Extract a bundle back to files with:
	python usher_assets.py --db x.db --name usher_seg4 --out_dir dir
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import sqlite3
from datetime import datetime

PROVENANCE_FILE = "usher_provenance.json"
PB_FILE = "usher.pb"
TREE_FILES = ("final-tree.nh", "uncondensed-final-tree.nh")
# column -> file. Stored gzipped: sample paths for a large tree run to GBs.
SIDECARS = {
	"sample_paths_gz": "sample-paths.tsv",
	"placement_stats_gz": "placement_stats.tsv",
	"identical_sequence_groups_gz": "identical_sequence_groups.tsv",
	"excluded_sequences_gz": "excluded_sequences.tsv",
}
KEY_COLUMNS = ("name", "source", "segment_key", "segment")

SCHEMA = (
	"CREATE TABLE IF NOT EXISTS usher_assets ("
	"name TEXT, source TEXT, segment_key TEXT, segment TEXT, "
	"tree_created_at TEXT, newick_sha256 TEXT, "
	"pb BLOB, pb_sha256 TEXT, "
	+ "".join(f"{col} BLOB, " for col in SIDECARS)
	+ "provenance TEXT, status TEXT, note TEXT, created_at TEXT);"
)

_KEY_WHERE = " AND ".join(f"COALESCE({col}, '')=?" for col in KEY_COLUMNS)
_TIP_RE = re.compile(r"[(,]\s*([^(),:;\[\]\s]+)")


def sha256_bytes(data):
	return hashlib.sha256(data).hexdigest()


def sha256_text(text):
	return sha256_bytes((text or "").strip().encode("utf-8"))


def sha256_file(path):
	digest = hashlib.sha256()
	with open(path, "rb") as handle:
		for block in iter(lambda: handle.read(1 << 20), b""):
			digest.update(block)
	return digest.hexdigest()


def _read_bytes(path):
	with open(path, "rb") as handle:
		return handle.read()


def tip_labels(newick):
	"""Tip labels of a Newick string. Internal labels follow ')' and are skipped."""
	return set(_TIP_RE.findall(newick or ""))


def _key_values(key):
	return tuple(str(key.get(col) or "").strip() for col in KEY_COLUMNS)


def write_provenance(output_dir, pb_produced=True, **fields):
	"""Pin the .pb in output_dir to the tree files it was produced with.

	Call once every output is final (after identical-sequence expansion), so
	the recorded tree hashes are those of the files a DB build will read.
	A run that produced no .pb records pb_sha256 = None, which stops a .pb
	left over from an earlier run in the same directory from being bundled.
	"""
	pb_path = os.path.join(output_dir, PB_FILE)
	record = {
		"created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
		"pb_sha256": sha256_file(pb_path) if pb_produced and os.path.isfile(pb_path) else None,
		"tree_sha256": {
			name: sha256_file(os.path.join(output_dir, name))
			for name in TREE_FILES
			if os.path.isfile(os.path.join(output_dir, name))
		},
		"sidecar_sha256": {
			name: sha256_file(os.path.join(output_dir, name))
			for name in SIDECARS.values()
			if os.path.isfile(os.path.join(output_dir, name))
		},
	}
	record.update(fields)
	path = os.path.join(output_dir, PROVENANCE_FILE)
	with open(path, "w", encoding="utf-8") as handle:
		json.dump(record, handle, indent=2, sort_keys=True)
	return path


def collect_bundle(tree_path):
	"""The .pb and sidecars that belong to tree_path, or None.

	None when there is no .pb, or when provenance shows the .pb (or the tree
	file) is not the one that run produced: storing it would recreate the
	mismatch this module exists to prevent.
	"""
	if not tree_path:
		return None
	out_dir = os.path.dirname(os.path.abspath(tree_path))
	pb_path = os.path.join(out_dir, PB_FILE)
	if not os.path.isfile(pb_path):
		return None

	pb_sha = sha256_file(pb_path)
	provenance = None
	prov_path = os.path.join(out_dir, PROVENANCE_FILE)
	if os.path.isfile(prov_path):
		with open(prov_path, encoding="utf-8") as handle:
			provenance = json.load(handle)

	if provenance is None:
		status, note = "unverified", f"no {PROVENANCE_FILE}; {PB_FILE} found beside the tree file"
	else:
		recorded_tree = (provenance.get("tree_sha256") or {}).get(os.path.basename(tree_path))
		if provenance.get("pb_sha256") != pb_sha:
			print(f"[usher_assets] Not bundling {pb_path}: it is not the .pb recorded in {prov_path}")
			return None
		if recorded_tree != sha256_file(tree_path):
			print(f"[usher_assets] Not bundling {pb_path}: {tree_path} changed after the .pb was produced")
			return None
		status, note = "verified", None

	bundle = {"pb": _read_bytes(pb_path), "pb_sha256": pb_sha, "status": status, "note": note,
	          "provenance": json.dumps(provenance, sort_keys=True) if provenance else None}
	for col, name in SIDECARS.items():
		path = os.path.join(out_dir, name)
		bundle[col] = gzip.compress(_read_bytes(path)) if os.path.isfile(path) else None
	return bundle


def ensure_schema(conn):
	conn.execute(SCHEMA)


def delete_bundle(conn, key):
	ensure_schema(conn)
	conn.execute(f"DELETE FROM usher_assets WHERE {_KEY_WHERE}", _key_values(key))


def store_bundle(conn, key, bundle, newick, tree_created_at):
	"""Replace the bundle for key with one linked to this Newick."""
	delete_bundle(conn, key)
	columns = list(KEY_COLUMNS) + ["tree_created_at", "newick_sha256", "pb", "pb_sha256"] + list(SIDECARS) + [
		"provenance", "status", "note", "created_at"]
	values = [key.get(col) for col in KEY_COLUMNS] + [tree_created_at, sha256_text(newick), bundle["pb"], bundle["pb_sha256"]]
	values += [bundle.get(col) for col in SIDECARS]
	values += [bundle.get("provenance"), bundle["status"], bundle.get("note"), datetime.now().strftime("%Y-%m-%d %H:%M:%S")]
	conn.execute(
		f"INSERT INTO usher_assets ({', '.join(columns)}) VALUES ({', '.join('?' for _ in columns)})", values)


def retarget_bundle(conn, key, old_newick, new_newick, tree_created_at, keep):
	"""Follow a trees row whose Newick was rewritten in place.

	keep=True claims the rewrite preserved the tree (e.g. rerooting). The claim
	is checked on tip sets; when it holds the bundle is relinked to the new
	Newick, otherwise it is marked stale. Returns the new status, or None when
	the row has no bundle.
	"""
	ensure_schema(conn)
	if conn.execute(f"SELECT COUNT(*) FROM usher_assets WHERE {_KEY_WHERE}", _key_values(key)).fetchone()[0] == 0:
		return None
	if keep and tip_labels(old_newick) == tip_labels(new_newick):
		status, note = "relinked", "trees row rewritten with the same tip set; .pb root may differ from the stored Newick"
	else:
		status = "stale"
		note = ("trees row rewritten with a different tip set" if keep
		        else "trees row replaced; the .pb describes the previous tree")
	conn.execute(
		f"UPDATE usher_assets SET status=?, note=?, newick_sha256=?, tree_created_at=? WHERE {_KEY_WHERE}",
		(status, note, sha256_text(new_newick), tree_created_at) + _key_values(key))
	return status


def reconcile(conn):
	"""Drop bundles whose trees row is gone; mark stale any whose row changed.

	A safety net for writers that do not know about usher_assets.
	"""
	ensure_schema(conn)
	trees = {}
	for row in conn.execute("SELECT name, source, segment_key, segment, newick FROM trees"):
		trees.setdefault(tuple(str(v or "").strip() for v in row[:4]), []).append(row[4])
	removed = marked = 0
	for rowid, *key_vals, newick_sha, status in conn.execute(
		"SELECT rowid, name, source, segment_key, segment, newick_sha256, status FROM usher_assets").fetchall():
		key = tuple(str(v or "").strip() for v in key_vals)
		if key not in trees:
			conn.execute("DELETE FROM usher_assets WHERE rowid=?", (rowid,))
			removed += 1
		elif status != "stale" and newick_sha not in {sha256_text(n) for n in trees[key]}:
			conn.execute(
				"UPDATE usher_assets SET status='stale', note=? WHERE rowid=?",
				("trees row changed without updating its UShER bundle", rowid))
			marked += 1
	return removed, marked


def extract_bundle(conn, key, out_dir):
	"""Write a stored bundle back out as usher.pb plus sidecars. Returns the row's status."""
	columns = ["pb", "provenance", "status"] + list(SIDECARS)
	row = conn.execute(f"SELECT {', '.join(columns)} FROM usher_assets WHERE {_KEY_WHERE}", _key_values(key)).fetchone()
	if row is None:
		raise LookupError(f"No usher_assets row for {dict(zip(KEY_COLUMNS, _key_values(key)))}")
	record = dict(zip(columns, row))
	os.makedirs(out_dir, exist_ok=True)
	with open(os.path.join(out_dir, PB_FILE), "wb") as handle:
		handle.write(record["pb"])
	if record["provenance"]:
		with open(os.path.join(out_dir, PROVENANCE_FILE), "w", encoding="utf-8") as handle:
			handle.write(record["provenance"])
	for col, name in SIDECARS.items():
		if record[col] is not None:
			with open(os.path.join(out_dir, name), "wb") as handle:
				handle.write(gzip.decompress(record[col]))
	return record["status"]


def main():
	parser = argparse.ArgumentParser(description="Extract a stored UShER bundle (usher.pb and sidecars) from a DB.")
	parser.add_argument("--db", required=True)
	parser.add_argument("--name", default="")
	parser.add_argument("--source", default="usher")
	parser.add_argument("--segment_key", default="")
	parser.add_argument("--segment", default="")
	parser.add_argument("--out_dir", required=True)
	args = parser.parse_args()
	conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
	try:
		status = extract_bundle(conn, vars(args), args.out_dir)
	finally:
		conn.close()
	print(f"[usher_assets] Wrote bundle to {args.out_dir} (status: {status})")
	if status != "verified":
		print(f"[usher_assets] Warning: bundle status is '{status}'")


if __name__ == "__main__":
	main()
