#!/usr/bin/env python3
'''
Clade assignment for query sequences.

Two strategies, chosen automatically:

  1. TREE-BASED (preferred). If an IQ-TREE or UShER tree is supplied, every query
     inherits the clade of its nearest labelled reference tip in that tree. This
     is consistent with the phylogeny the rest of the pipeline produces: a query
     cannot be assigned a clade that disagrees with all of its own tree
     neighbours. UShER is preferred over IQ-TREE because it holds every placed
     sample, whereas the IQ-TREE backbone only holds cluster representatives.

  2. EPA-ng FALLBACK. Only when no tree is available:
       - build a reference-only tree with iqtree2,
       - place queries onto it with epa-ng,
       - assign clades with `gappa examine assign --per-query-results`,
     taking gappa's BEST hit per query (highest accumulated likelihood weight),
     not whichever row happens to appear last in the file.

Either way the query's major clade (genotype) and minor clade (subtype) are
written back into the GenBank matrix.
'''

import os
import csv
import subprocess
from os.path import join
from argparse import ArgumentParser

import read_file
from clade_from_tree import assign_labels_from_tree


def _clean(value):
	text = "" if value is None else str(value).strip()
	if text.lower() in {"", "na", "n/a", "nan", "none", "null"}:
		return ""
	return text


class CladeAssignment:
	def __init__(self, major_clade, minor_clade, padded_alignment, base_dir, output_dir, gb_matrix,
				 threads, iqtree_model, iqtree_tree=None, usher_tree=None):
		self.major_clade = major_clade
		self.minor_clade = minor_clade
		self.padded_alignment = padded_alignment
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.gb_matrix = gb_matrix
		self.threads = str(threads)
		self.iqtree_model = iqtree_model
		self.iqtree_tree = self._normalize_optional_path(iqtree_tree)
		self.usher_tree = self._normalize_optional_path(usher_tree)

	@staticmethod
	def _normalize_optional_path(path_value):
		if path_value is None:
			return None
		text = str(path_value).strip()
		if not text or text.lower() in {"null", "unset", "none"}:
			return None
		return text

	# ------------------------------------------------------------------ inputs

	@staticmethod
	def _read_taxon_file(path):
		"""Reference taxon file: <accession>\\t<clade>. Tolerates a header row and
		extra columns; keeps the first non-empty clade seen per accession."""
		mapping = {}
		if not path or not os.path.isfile(path):
			return mapping
		with open(path, newline='', encoding='utf-8') as handle:
			for row in csv.reader(handle, delimiter='\t'):
				if len(row) < 2:
					continue
				acc = _clean(row[0])
				clade = _clean(row[1])
				if acc.lower() in {"name", "accession", "primary_accession"}:
					continue
				if acc and clade and acc not in mapping:
					mapping[acc] = clade
		return mapping

	@staticmethod
	def _matrix_key(row):
		"""Accession key that matches tree tips / FASTA headers. Prefer
		primary_accession; fall back to gi_number then locus."""
		for field in ("primary_accession", "gi_number", "locus"):
			value = _clean(row.get(field))
			if value:
				return value
		return ""

	def _query_and_reference_accessions(self):
		"""Split matrix accessions into query vs reference, skipping excluded rows."""
		queries, references = set(), set()
		with open(self.gb_matrix, newline='', encoding='utf-8') as handle:
			for row in csv.DictReader(handle, delimiter='\t'):
				if _clean(row.get("exclusion_status")) == "1":
					continue
				acc = self._matrix_key(row)
				if not acc:
					continue
				if _clean(row.get("accession_type")).lower() == "query":
					queries.add(acc)
				else:
					references.add(acc)
		return queries, references

	# -------------------------------------------------------------- tree path

	def _select_tree(self):
		"""Return (path, source) for the most trustworthy available tree."""
		if self.usher_tree and os.path.isfile(self.usher_tree):
			return self.usher_tree, "usher"
		if self.iqtree_tree and os.path.isfile(self.iqtree_tree):
			return self.iqtree_tree, "iqtree"
		return None, None

	def assign_from_tree(self, tree_path, reference_major, reference_minor):
		"""Assign each query tip the clade of its nearest labelled reference."""
		reference_labels = {}
		for acc in set(reference_major) | set(reference_minor):
			reference_labels[acc] = {
				"genotype": reference_major.get(acc, ""),
				"subtype": reference_minor.get(acc, ""),
			}
		assignments = assign_labels_from_tree(tree_path, reference_labels)
		major = {acc: labels.get("genotype", "") for acc, labels in assignments.items()}
		minor = {acc: labels.get("subtype", "") for acc, labels in assignments.items()}
		return major, minor

	# ------------------------------------------------------------- EPA-ng path

	def alignment(self):
		"""Split the padded alignment into query-only and reference-only FASTAs
		for the EPA-ng fallback."""
		output_path = join(self.base_dir, self.output_dir)
		os.makedirs(output_path, exist_ok=True)

		queries, _references = self._query_and_reference_accessions()
		query_path = join(output_path, "query_aln.fa")
		reference_path = join(output_path, "reference_aln.fa")
		with open(query_path, "w", encoding='utf-8') as query_out, \
				open(reference_path, "w", encoding='utf-8') as reference_out:
			for header, seq in read_file.fasta(self.padded_alignment):
				acc = _clean(header).split()[0] if header else ""
				target = query_out if acc in queries else reference_out
				target.write(">" + acc + "\n" + seq.strip() + "\n")
		return query_path, reference_path

	@staticmethod
	def _resolve_iqtree_binary():
		"""The pipeline environment ships iqtree3; older installs have iqtree2."""
		for candidate in ("iqtree3", "iqtree2", "iqtree"):
			if shutil.which(candidate):
				return candidate
		raise FileNotFoundError("No iqtree binary found on PATH (looked for iqtree3, iqtree2, iqtree)")

	def iqtree(self, ref_fa, prefix):
		command = [self._resolve_iqtree_binary(), '-s', ref_fa, '-m', self.iqtree_model,
				   '-nt', self.threads, '-pre', prefix]
		subprocess.run(command, check=True)

	def epa_ng(self, ref_fa, query_fa, prefix):
		command = [
			'epa-ng', '--redo', '-m', self.iqtree_model,
			'-t', prefix + '.treefile', '-s', ref_fa, '-q', query_fa,
			'--outdir', join(self.base_dir, self.output_dir),
		]
		subprocess.run(command, check=True)

	def _gappa_assign(self, taxon_file, out_subdir):
		out_dir = join(self.base_dir, self.output_dir, out_subdir)
		command = [
			'gappa', 'examine', 'assign',
			'--jplace-path', join(self.base_dir, self.output_dir, 'epa_result.jplace'),
			'--taxon-file', taxon_file,
			'--out-dir', out_dir,
			'--per-query-results',
		]
		subprocess.run(command, check=True)
		return join(out_dir, 'per_query.tsv')

	@staticmethod
	def parse_gappa_per_query(per_query_file):
		"""Return {query: best_taxopath}. gappa lists multiple candidate taxopaths
		per query, sorted best-first by accumulated LWR; keep the highest-weight
		row, not whichever appears last."""
		best = {}
		if not per_query_file or not os.path.isfile(per_query_file):
			return best
		with open(per_query_file, 'r', newline='', encoding='utf-8') as handle:
			reader = csv.DictReader(handle, delimiter='\t')
			for row in reader:
				name = _clean(row.get('name'))
				taxopath = _clean(row.get('taxopath'))
				if not name or not taxopath:
					continue
				# aLWR (accumulated likelihood weight ratio) is gappa's confidence;
				# fall back to LWR, then to file order, when the column is absent.
				try:
					score = float(row.get('aLWR') or row.get('LWR') or 0.0)
				except (TypeError, ValueError):
					score = 0.0
				current = best.get(name)
				if current is None or score > current[0]:
					best[name] = (score, taxopath)
		return {name: taxopath for name, (_score, taxopath) in best.items()}

	@staticmethod
	def _leaf_taxon(taxopath):
		"""gappa taxopaths are ';'-separated; the clade label is the deepest rank."""
		parts = [p for p in str(taxopath).split(';') if _clean(p)]
		return _clean(parts[-1]) if parts else ""

	def assign_with_epa_ng(self):
		query_fa, ref_fa = self.alignment()
		prefix = join(self.base_dir, self.output_dir, 'ref_tree')
		self.iqtree(ref_fa, prefix)
		self.epa_ng(ref_fa, query_fa, prefix)
		major_pq = self._gappa_assign(self.major_clade, 'gappa_major_clades')
		minor_pq = self._gappa_assign(self.minor_clade, 'gappa_minor_clades')
		major = {name: self._leaf_taxon(tp) for name, tp in self.parse_gappa_per_query(major_pq).items()}
		minor = {name: self._leaf_taxon(tp) for name, tp in self.parse_gappa_per_query(minor_pq).items()}
		return major, minor

	# ----------------------------------------------------------------- output

	def write_clades_to_gb_matrix(self, major_clades, minor_clades):
		temp_file = self.gb_matrix + '.tmp'
		with open(self.gb_matrix, 'r', newline='', encoding='utf-8') as infile, \
				open(temp_file, 'w', newline='', encoding='utf-8') as outfile:
			reader = csv.DictReader(infile, delimiter='\t')
			base_fields = [f for f in reader.fieldnames if f not in ('major_clade', 'minor_clade')]
			new_fields = base_fields + ['major_clade', 'minor_clade']
			writer = csv.DictWriter(outfile, fieldnames=new_fields, delimiter='\t')
			writer.writeheader()
			for row in reader:
				acc = self._matrix_key(row)
				row = {k: v for k, v in row.items() if k in base_fields}
				row['major_clade'] = major_clades.get(acc, '')
				row['minor_clade'] = minor_clades.get(acc, '')
				writer.writerow(row)
		os.replace(temp_file, self.gb_matrix)
		print(f"Updated file in place: {self.gb_matrix}")

	def process(self):
		os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)
		reference_major = self._read_taxon_file(self.major_clade)
		reference_minor = self._read_taxon_file(self.minor_clade)

		tree_path, source = self._select_tree()
		if tree_path:
			print(f"[info] Assigning clades from the {source} tree: {tree_path}")
			query_major, query_minor = self.assign_from_tree(tree_path, reference_major, reference_minor)
		else:
			print("[info] No IQ-TREE/UShER tree available; falling back to EPA-ng placement.")
			query_major, query_minor = self.assign_with_epa_ng()

		# References keep their curated clade; queries take the assigned one.
		major = dict(reference_major)
		major.update({acc: clade for acc, clade in query_major.items() if clade})
		minor = dict(reference_minor)
		minor.update({acc: clade for acc, clade in query_minor.items() if clade})
		self.write_clades_to_gb_matrix(major, minor)


if __name__ == "__main__":
	parser = ArgumentParser(description='Assign clades (genotype/subtype) to query sequences')
	parser.add_argument('-c', '--major_clade', help='Major clade (genotype) reference TSV', default='generic/rabv/major_clades.tsv')
	parser.add_argument('-s', '--minor_clade', help='Minor clade (subtype) reference TSV', default='generic/rabv/minor_clades.tsv')
	parser.add_argument('-p', '--padded_alignment', help='Padded alignment FASTA', default='tmp/Pad-alignment/NC_001542.aligned_merged_MSA.fasta')
	parser.add_argument('-b', '--base_dir', help='Base directory', default='tmp')
	parser.add_argument('-o', '--output_dir', help='Output directory', default='CladeAssignment')
	parser.add_argument('-g', '--gb_matrix', help='GenBank matrix TSV', default='tmp/GenBank-matrix/gB_matrix_raw.tsv')
	parser.add_argument('-t', '--threads', help='Threads for iqtree (EPA-ng fallback)', default='7')
	parser.add_argument('-m', '--iqtree_model', help='iqtree model (EPA-ng fallback)', default='GTR+G')
	parser.add_argument('-it', '--iqtree_tree', help='IQ-TREE newick with query + reference tips', default=None)
	parser.add_argument('-ut', '--usher_tree', help='UShER newick with query + reference tips', default=None)
	args = parser.parse_args()

	CladeAssignment(
		args.major_clade,
		args.minor_clade,
		args.padded_alignment,
		args.base_dir,
		args.output_dir,
		args.gb_matrix,
		args.threads,
		args.iqtree_model,
		iqtree_tree=args.iqtree_tree,
		usher_tree=args.usher_tree,
	).process()
