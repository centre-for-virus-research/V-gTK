"""The UShER mutation-annotated tree (usher.pb) travels in the DB with its tree.

Before usher_assets the .pb stayed in the UsherPlacement output directory and
nothing tied it to a trees row, so a replaced row left a .pb describing a
different tree. These tests pin the link: bundles are keyed like trees rows,
refused when provenance says they belong to another tree, and relinked or
marked stale whenever the row they describe is rewritten.
"""

import gzip
import json
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

import usher_assets
from CreateSqliteDB import CreateSqliteDB
from ReplaceUsherTreeInDb import replace_tree
from UsherPlacement import UsherPlacement


TREE = "((A:0.1,B:0.2):0.1,C:0.3);\n"


def make_usher_dir(path: Path, tree=TREE, pb=b"PB-v1", provenance=True):
	path.mkdir(parents=True, exist_ok=True)
	(path / "final-tree.nh").write_text(tree, encoding="utf-8")
	(path / "usher.pb").write_bytes(pb)
	(path / "sample-paths.tsv").write_text("A\tnode_1:A1G\n", encoding="utf-8")
	(path / "placement_stats.tsv").write_text("A\t1\t1\t\n", encoding="utf-8")
	(path / "identical_sequence_groups.tsv").write_text("anchor_id\tmember_id\tanchor_requires_placement\n", encoding="utf-8")
	if provenance:
		usher_assets.write_provenance(str(path), placed=1)
	return path / "final-tree.nh"


def trees_db(tmp_path: Path, rows):
	db = tmp_path / "t.db"
	conn = sqlite3.connect(db)
	conn.execute("CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT)")
	conn.executemany("INSERT INTO trees VALUES (?, ?, ?, ?, ?, ?)", rows)
	conn.commit()
	conn.close()
	return db


KEY = {"name": "usher", "source": "usher", "segment_key": None, "segment": "1"}


def asset_rows(db):
	conn = sqlite3.connect(db)
	try:
		return conn.execute("SELECT name, status, newick_sha256, pb FROM usher_assets").fetchall()
	finally:
		conn.close()


def test_collect_bundle_verified_and_round_trips(tmp_path: Path):
	tree = make_usher_dir(tmp_path / "u")
	bundle = usher_assets.collect_bundle(str(tree))
	assert bundle["status"] == "verified"

	db = trees_db(tmp_path, [("usher", "usher", None, "1", TREE.strip(), "t0")])
	conn = sqlite3.connect(db)
	usher_assets.store_bundle(conn, KEY, bundle, TREE, "t0")
	status = usher_assets.extract_bundle(conn, {"name": "usher", "source": "usher", "segment": "1"}, str(tmp_path / "x"))
	conn.close()

	assert status == "verified"
	assert (tmp_path / "x" / "usher.pb").read_bytes() == b"PB-v1"
	assert (tmp_path / "x" / "sample-paths.tsv").read_text() == "A\tnode_1:A1G\n"
	assert json.loads((tmp_path / "x" / "usher_provenance.json").read_text())["placed"] == 1


def test_collect_bundle_refuses_pb_not_from_this_tree(tmp_path: Path):
	"""The original failure: a tree rewritten (or a .pb swapped) after placement."""
	tree = make_usher_dir(tmp_path / "u")
	tree.write_text("((A:0.1,C:0.2):0.1,B:0.3);\n", encoding="utf-8")
	assert usher_assets.collect_bundle(str(tree)) is None

	tree = make_usher_dir(tmp_path / "v")
	(tmp_path / "v" / "usher.pb").write_bytes(b"PB-from-another-run")
	assert usher_assets.collect_bundle(str(tree)) is None


def test_run_without_placement_disowns_leftover_pb(tmp_path: Path):
	tree = make_usher_dir(tmp_path / "u")
	usher_assets.write_provenance(str(tmp_path / "u"), pb_produced=False, placed=0)
	assert usher_assets.collect_bundle(str(tree)) is None


def test_collect_bundle_without_provenance_is_unverified(tmp_path: Path):
	tree = make_usher_dir(tmp_path / "u", provenance=False)
	assert usher_assets.collect_bundle(str(tree))["status"] == "unverified"
	(tmp_path / "u" / "usher.pb").unlink()
	assert usher_assets.collect_bundle(str(tree)) is None


def test_reconcile_drops_orphans_and_flags_changed_rows(tmp_path: Path):
	db = trees_db(tmp_path, [("usher", "usher", None, "1", TREE.strip(), "t0")])
	bundle = usher_assets.collect_bundle(str(make_usher_dir(tmp_path / "u")))
	conn = sqlite3.connect(db)
	usher_assets.store_bundle(conn, KEY, bundle, TREE, "t0")
	usher_assets.store_bundle(conn, {**KEY, "name": "gone"}, bundle, TREE, "t0")
	conn.execute("UPDATE trees SET newick='(Z:1,Y:1);'")
	assert usher_assets.reconcile(conn) == (1, 1)
	conn.commit()
	conn.close()
	assert [(r[0], r[1]) for r in asset_rows(db)] == [("usher", "stale")]


@pytest.mark.parametrize(
	"new_tree, keep, expected",
	[
		("(C:0.3,(A:0.1,B:0.2):0.1);\n", True, "relinked"),
		("(C:0.3,(A:0.1,D:0.2):0.1);\n", True, "stale"),
		("(C:0.3,(A:0.1,B:0.2):0.1);\n", False, "stale"),
	],
)
def test_replace_tree_never_leaves_a_pb_claiming_the_wrong_tree(tmp_path: Path, new_tree, keep, expected):
	db = trees_db(tmp_path, [("usher", "usher", None, "1", TREE.strip(), "t0")])
	conn = sqlite3.connect(db)
	usher_assets.store_bundle(conn, KEY, usher_assets.collect_bundle(str(make_usher_dir(tmp_path / "u"))), TREE, "t0")
	conn.commit()
	conn.close()

	replacement = tmp_path / "new.nwk"
	replacement.write_text(new_tree, encoding="utf-8")
	replace_tree(db, replacement, name="usher", segment="1", keep_usher_assets=keep)

	(name, status, newick_sha, _), = asset_rows(db)
	assert status == expected
	assert newick_sha == usher_assets.sha256_text(new_tree)


def test_replace_tree_with_usher_dir_stores_the_new_bundle(tmp_path: Path):
	db = trees_db(tmp_path, [("usher", "usher", None, "1", "(OLD:1);", "t0")])
	new_tree = make_usher_dir(tmp_path / "u2", pb=b"PB-v2")
	replace_tree(db, new_tree, name="usher", segment="1", usher_dir=str(tmp_path / "u2"))
	(_, status, _, pb), = asset_rows(db)
	assert (status, pb) == ("verified", b"PB-v2")

	with pytest.raises(ValueError):
		replace_tree(db, tmp_path / "u2" / "final-tree.nh", name="usher", segment="1", usher_dir=str(tmp_path / "nope"))


def test_promote_merges_per_batch_stats_and_paths(tmp_path: Path):
	"""Each batch writes placement stats for its own samples; the last batch alone
	used to be promoted, losing every earlier placement."""
	processor = UsherPlacement(padded_aln=str(tmp_path / "a.fa"), output_dir=str(tmp_path / "out"))
	(tmp_path / "out").mkdir()
	chunks = []
	for idx, sample in enumerate(["S1", "S2"], start=1):
		chunk = tmp_path / "out" / f"chunk_{idx:04d}"
		chunk.mkdir()
		(chunk / "placement_stats.tsv").write_text(f"{sample}\t1\t1\t\n")
		(chunk / "mutation-paths.txt").write_text(f"{sample}\tnode_1:A1G\n")
		(chunk / "usher.pb").write_bytes(f"pb{idx}".encode())
		chunks.append(str(chunk))
	processor.promote_final_usher_outputs(chunks[-1], chunks)
	out = tmp_path / "out"
	assert (out / "placement_stats.tsv").read_text() == "S1\t1\t1\t\nS2\t1\t1\t\n"
	assert (out / "mutation-paths.txt").read_text().count("\n") == 2
	assert (out / "usher.pb").read_bytes() == b"pb2"


def _write_inputs(tmp_path: Path):
	def tsv(name, rows, cols):
		pd.DataFrame(rows, columns=cols).to_csv(tmp_path / name, sep="\t", index=False)
		return str(tmp_path / name)

	def csv(name, rows, cols):
		pd.DataFrame(rows, columns=cols).to_csv(tmp_path / name, index=False)
		return str(tmp_path / name)

	(tmp_path / "seqs.fa").write_text(">A\nATGC\n>B\nATGA\n>C\nATGG\n", encoding="utf-8")
	return dict(
		meta_data=tsv("meta.tsv", [["A", ""], ["B", ""], ["C", ""]], ["primary_accession", "exclusion"]),
		features=tsv("features.tsv", [["A", "P"]], ["primary_accession", "feature"]),
		pad_aln=tsv("aln.tsv", [["A", "ATGC"]], ["primary_accession", "aligned_seq"]),
		gene_info=tsv("gene.tsv", [["g", "G"]], ["name", "description"]),
		m49_countries=csv("c.csv", [["001", "World"]], ["m49_code", "name"]),
		m49_interm_region=csv("i.csv", [["X", "I"]], ["code", "name"]),
		m49_regions=csv("r.csv", [["Y", "R"]], ["code", "name"]),
		m49_sub_regions=csv("s.csv", [["Z", "S"]], ["code", "name"]),
		proj_settings=tsv("sw.tsv", [["Python", "3"]], ["Software", "Version"]),
		fasta_sequence_file=str(tmp_path / "seqs.fa"),
		insertions=tsv("ins.tsv", [["A", "none"]], ["primary_accession", "insertions"]),
		host_taxa_file=tsv("host.tsv", [["A", "h"]], ["primary_accession", "host"]),
		base_dir=str(tmp_path),
		output_dir="SqliteDB",
		db_name="testdb",
		db_status="new",
		tree_file=None,
		iqtree_file=None,
	)


def test_create_db_bundles_the_pb_with_its_usher_row(tmp_path: Path):
	tree = make_usher_dir(tmp_path / "Usher_seg1")
	CreateSqliteDB(usher_tree=str(tree), **_write_inputs(tmp_path)).create_db()

	conn = sqlite3.connect(tmp_path / "SqliteDB" / "testdb.db")
	trees = conn.execute("SELECT name, segment, newick, created_at FROM trees WHERE source='usher'").fetchall()
	assets = conn.execute(
		"SELECT name, segment, newick_sha256, tree_created_at, pb, status, sample_paths_gz FROM usher_assets").fetchall()
	conn.close()

	assert len(trees) == 1 and len(assets) == 1
	(t_name, t_seg, newick, created_at), (a_name, a_seg, sha, a_created, pb, status, paths) = trees[0], assets[0]
	assert (a_name, a_seg, a_created) == (t_name, t_seg, created_at)
	assert sha == usher_assets.sha256_text(newick)
	assert (pb, status) == (b"PB-v1", "verified")
	assert gzip.decompress(paths) == b"A\tnode_1:A1G\n"


def test_create_db_does_not_bundle_a_pb_for_another_tree(tmp_path: Path):
	tree = make_usher_dir(tmp_path / "Usher_seg1")
	tree.write_text("((A:0.1,C:0.2):0.1,B:0.3);\n", encoding="utf-8")
	CreateSqliteDB(usher_tree=str(tree), **_write_inputs(tmp_path)).create_db()

	conn = sqlite3.connect(tmp_path / "SqliteDB" / "testdb.db")
	assert conn.execute("SELECT COUNT(*) FROM trees WHERE source='usher'").fetchone()[0] == 1
	assert conn.execute("SELECT COUNT(*) FROM usher_assets").fetchone()[0] == 0
	conn.close()
