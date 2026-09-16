#!/usr/bin/env python3
"""Check a reference backbone holds every reference, and a cached one is still current.

WHY THIS EXISTS
===============
BUILD_REFERENCE_ALIGNMENT aligns every reference with MAFFT. For HCV's 238
genomes that is ~14 minutes on 8 threads, repeated identically on every fresh
build of the same reference list. So a profile can point ``params.ref_set_aligned``
at a backbone built once and kept under ``generic/<virus>/``.

A cached backbone goes stale without any error: add a reference to the list and
the cache silently lacks its row; change ``ref_aln_min_insertion_support`` and
the cache keeps columns the build would have dropped. This check turns both
into a failure that says how to rebuild.

WHAT IS CHECKED
===============
Every backbone - cached, hand-curated or built in this run:

``coverage``               every master/reference in the list must be a row of the
                           backbone file the pipeline will open for its segment
                           (projectability.find_precomputed_reference_alignment,
                           the resolver PadAlignment uses). A missing reference is
                           a FAIL: its queries would be aligned against a backbone
                           that does not hold it. Rows not in the list are a WARN.

A cached backbone, additionally against ``backbone_manifest.tsv``, written
beside it by ``--write_manifest``:

``reference_accessions``   the set of accessions in the reference list (after the
                           builder's own normalisation) must equal the set in the
                           backbone FASTAs, and must match the manifest. FAIL.
``min_insertion_support``  must equal the run's setting. FAIL.
``builder_sha256``         BuildReferenceAlignment.py's hash when the cache was
                           built. A different hash means the builder has changed
                           since: WARN, because the cache may no longer be what
                           the builder produces.

A directory with no manifest (hand-curated, or built in this run) gets the
coverage check only.

References are never added to a backbone or a database after the fact: a
changed reference list means a rebuild.

Usage:
    CheckReferenceBackbone.py --backbone_dir DIR --ref_list LIST --min_insertion_support 2
    CheckReferenceBackbone.py --backbone_dir DIR --ref_list LIST --min_insertion_support 2 \\
        --write_manifest --source "dev/nf_runs/HCV_test/out/ref_set_aligned (HCV_test, 15 Sep 2026)"
Exit code 0 = usable, 1 = stale or broken.
"""

import argparse
import datetime
import glob
import hashlib
import os
import sys

import accession_utils
import projectability
import segment_utils
from ExportRefListFromUpdateDb import load_reference_file_table

MANIFEST = "backbone_manifest.tsv"
BUILDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BuildReferenceAlignment.py")


def _key(value):
	value = str(value).strip()
	return accession_utils.normalise_accession(value) or value


def reference_rows(ref_list):
	"""``[(accession, segment)]`` for master/reference rows, accessions normalised as the builder does.

	``exclusion_list`` rows (influenza B/C/D in the flu list) name sequences to keep
	out of the build, not references, so they are not expected in the backbone.
	"""
	table = load_reference_file_table(ref_list)
	types = table["accession_type"].fillna("").astype(str).str.strip().str.lower()
	segments = table["segment"] if "segment" in table.columns else [""] * len(table)
	return [(_key(acc), str(seg or "").strip())
			for acc, seg, kind in zip(table["primary_accession"], segments, types)
			if str(acc).strip() and kind != "exclusion_list"]


def reference_accessions(ref_list):
	"""Unique master/reference accessions in the reference list."""
	return {acc for acc, _ in reference_rows(ref_list)}


def read_backbone_file(path):
	with open(path) as handle:
		return [_key(line[1:].split()[0]) for line in handle if line.startswith(">") and line[1:].split()]


def backbone_accessions(backbone_dir):
	"""``{fasta name: [accessions in file order]}`` for every refset_*_aln.fasta."""
	return {os.path.basename(path): read_backbone_file(path)
			for path in sorted(glob.glob(os.path.join(backbone_dir, "refset_*_aln.fasta")))}


def coverage(backbone_dir, ref_list, is_segmented="N"):
	"""``(missing, extra)``: references absent from the file their segment resolves to, and unlisted rows.

	Unsegmented builds have one backbone, so a reference only has to be in it.
	Segmented builds must have each reference in its own segment's file; a
	reference with no segment in the list only has to be in some file.
	"""
	files = backbone_accessions(backbone_dir)
	everywhere = {acc for accs in files.values() for acc in accs}
	rows = reference_rows(ref_list)
	cache = {}
	missing = []
	for acc, segment in rows:
		label = segment_utils.normalise_segment(segment) if segment else None
		if str(is_segmented).upper() != "Y" or not label:
			if acc not in everywhere:
				missing.append(acc)
			continue
		if label not in cache:
			path = projectability.find_precomputed_reference_alignment(backbone_dir, label)
			cache[label] = (os.path.basename(path), set(read_backbone_file(path))) if path else (None, set())
		name, members = cache[label]
		if acc not in members:
			missing.append(f"{acc} (segment {label}: {name or 'no backbone file'})")
	listed = {acc for acc, _ in rows}
	return sorted(set(missing)), sorted(everywhere - listed)


def digest_accessions(accessions):
	return hashlib.sha256("\n".join(sorted(accessions)).encode()).hexdigest()


def file_sha256(path):
	with open(path, "rb") as handle:
		return hashlib.sha256(handle.read()).hexdigest()


def read_manifest(backbone_dir):
	path = os.path.join(backbone_dir, MANIFEST)
	if not os.path.isfile(path):
		return None
	with open(path) as handle:
		return dict(line.rstrip("\n").split("\t", 1) for line in handle if "\t" in line)


def check(backbone_dir, ref_list, min_insertion_support, is_segmented="N"):
	"""``(failures, warnings, notes)``."""
	failures, warnings, notes = [], [], []
	files = backbone_accessions(backbone_dir)
	if not files:
		return [f"no refset_*_aln.fasta in {backbone_dir}"], warnings, notes
	in_backbone = {acc for accs in files.values() for acc in accs}

	missing, extra = coverage(backbone_dir, ref_list, is_segmented)
	if missing:
		failures.append(f"{len(missing)} reference(s) from the list are not in the backbone: "
						+ ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))
	manifest = read_manifest(backbone_dir)
	if manifest is None:
		if extra:
			warnings.append(f"{len(extra)} backbone row(s) are not in the reference list: {extra[:10]}")
		notes.append(f"{len(in_backbone)} rows across {len(files)} backbone file(s), no manifest "
					 "(hand-curated or built in this run): coverage checked")
		return failures, warnings, notes

	# A cache was built from exactly its list, so an unlisted row means the list shrank.
	if extra:
		failures.append(f"{len(extra)} backbone row(s) are not in the reference list: {extra[:10]}")
	if not missing and not extra and manifest.get("reference_accessions_sha256") != digest_accessions(reference_accessions(ref_list)):
		failures.append("reference accessions differ from those the cache was built for")
	if str(manifest.get("min_insertion_support")) != str(min_insertion_support):
		failures.append(
			f"cache built with min_insertion_support={manifest.get('min_insertion_support')}, "
			f"this run uses {min_insertion_support}")
	if manifest.get("builder_sha256") != file_sha256(BUILDER):
		warnings.append(
			f"BuildReferenceAlignment.py has changed since the cache was made on "
			f"{manifest.get('cached_on', '?')}; the cache may not match what it now produces")
	notes.append(f"{len(in_backbone)} references across {len(files)} backbone file(s); "
				 f"cached {manifest.get('cached_on', '?')} from {manifest.get('source', '?')}")
	return failures, warnings, notes


def write_manifest(backbone_dir, ref_list, min_insertion_support, source):
	files = backbone_accessions(backbone_dir)
	in_backbone = {acc for accs in files.values() for acc in accs}
	in_list = reference_accessions(ref_list)
	if in_backbone != in_list:
		raise ValueError(f"backbone and reference list disagree ({len(in_list ^ in_backbone)} accessions); "
						 "refusing to write a manifest for a mismatched cache")
	rows = {
		"cached_on": datetime.date.today().isoformat(),
		"source": source,
		"reference_list": os.path.basename(ref_list),
		"reference_count": str(len(in_list)),
		"reference_accessions_sha256": digest_accessions(in_list),
		"min_insertion_support": str(min_insertion_support),
		"builder_sha256": file_sha256(BUILDER),
		"files": ";".join(files),
	}
	with open(os.path.join(backbone_dir, MANIFEST), "w") as handle:
		for key, value in rows.items():
			handle.write(f"{key}\t{value}\n")


def main(argv=None):
	parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
	parser.add_argument("--backbone_dir", required=True)
	parser.add_argument("--ref_list", required=True)
	parser.add_argument("--min_insertion_support", required=True)
	parser.add_argument("--is_segmented", default="N", help="Y: each reference must be in its own segment's backbone file")
	parser.add_argument("--write_manifest", action="store_true", help="record this backbone's provenance")
	parser.add_argument("--source", default="", help="with --write_manifest: where the backbone came from")
	args = parser.parse_args(argv)

	if args.write_manifest:
		write_manifest(args.backbone_dir, args.ref_list, args.min_insertion_support, args.source)
		print(f"Wrote {os.path.join(args.backbone_dir, MANIFEST)}")
		return 0

	failures, warnings, notes = check(args.backbone_dir, args.ref_list, args.min_insertion_support, args.is_segmented)
	for note in notes:
		print(f"[backbone] {note}")
	for warning in warnings:
		print(f"[backbone][warn] {warning}")
	if failures:
		for failure in failures:
			print(f"[backbone][FAIL] {failure}", file=sys.stderr)
		print("[backbone] References are never added to a backbone after the fact. For a cache, rebuild: run the "
			  "profile with --ref_set_aligned null (the backbone is published to <publish_dir>/ref_set_aligned), "
			  "copy it over the cache, then run CheckReferenceBackbone.py --write_manifest. For a hand-curated "
			  "backbone, add the missing references to it or use a reference list that matches.", file=sys.stderr)
		return 1
	return 0


if __name__ == "__main__":
	sys.exit(main())
