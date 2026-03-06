#!/usr/bin/env python3
import argparse
import os
import re
import sqlite3
import subprocess


class UsherPlacement:
	def __init__(self, padded_aln, output_dir, mmseq_cluster_dir=None, iqtree_dir=None, update_db=None, threads=1, test_mode=False):
		self.padded_aln = padded_aln
		self.output_dir = output_dir
		self.mmseq_cluster_dir = self._normalize_optional_path(mmseq_cluster_dir)
		self.iqtree_dir = self._normalize_optional_path(iqtree_dir)
		self.update_db = self._normalize_optional_path(update_db)
		self.threads = max(1, int(threads))
		self.test_mode = str(test_mode).strip() == "1" if not isinstance(test_mode, bool) else test_mode

	@staticmethod
	def _normalize_optional_path(path_value):
		if path_value is None:
			return None
		path_str = str(path_value).strip()
		if not path_str or path_str.lower() == "null" or path_str == "UNSET":
			return None
		return path_str

	@staticmethod
	def _normalize_segment(value):
		if value is None:
			return "0"
		s = str(value).strip()
		if not s:
			return "0"
		digits = ''.join(ch for ch in s if ch.isdigit())
		return digits if digits else s

	@staticmethod
	def _find_first(path_root, pattern_suffix=None, exact_name=None):
		for root, _, files in os.walk(path_root):
			for name in files:
				if exact_name and name == exact_name:
					return os.path.join(root, name)
				if pattern_suffix and name.endswith(pattern_suffix):
					return os.path.join(root, name)
		return None

	@staticmethod
	def _read_ids_from_fasta(fasta_path):
		ids = []
		with open(fasta_path, "r", encoding="utf-8") as handle:
			for line in handle:
				if line.startswith(">"):
					ids.append(line[1:].strip().split()[0])
		return ids

	def _segment_from_padded_alignment(self):
		name = os.path.basename(self.padded_aln)
		name = re.sub(r"_dedup\.fasta$", "", name)
		match = re.search(r"(\d+)", name)
		return match.group(1) if match else "0"

	def prepare_update_assets(self):
		if not self.update_db:
			raise ValueError("update_db is required for update-mode USHER placement")
		if not os.path.isfile(self.update_db):
			raise FileNotFoundError(f"Update DB not found: {self.update_db}")
		os.makedirs(self.output_dir, exist_ok=True)

		segment = self._segment_from_padded_alignment()
		tree_out = os.path.join(self.output_dir, f"tree_segment_{segment}.nwk")
		ids_out = os.path.join(self.output_dir, f"existing_ids_segment_{segment}.txt")

		conn = sqlite3.connect(self.update_db)
		try:
			cur = conn.cursor()

			tree_written = False
			if cur.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='trees'").fetchone():
				rows = cur.execute(
					"SELECT COALESCE(segment, '0') AS segment, COALESCE(source, '') AS source, newick "
					"FROM trees WHERE newick IS NOT NULL AND TRIM(newick) != ''"
				).fetchall()
				priority = {"usher": 0, "iqtree": 1, "veryfasttree": 2}
				best = None
				for row in rows:
					row_seg = self._normalize_segment(row[0])
					if row_seg != segment:
						continue
					source = str(row[1]).strip().lower()
					score = priority.get(source, 9)
					if best is None or score < best[0]:
						best = (score, row[2])
				if best and best[1]:
					with open(tree_out, "w", encoding="utf-8") as handle:
						handle.write(str(best[1]).strip() + "\n")
					tree_written = True

			if not tree_written:
				raise ValueError(f"Missing tree for segment {segment} in update DB")

			excluded = set()
			if cur.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='excluded_accessions'").fetchone():
				for (acc,) in cur.execute("SELECT primary_accession FROM excluded_accessions"):
					if acc:
						excluded.add(str(acc).strip())

			ids = []
			for acc, segv in cur.execute("SELECT primary_accession, COALESCE(segment, '0') FROM meta_data"):
				if not acc:
					continue
				acc = str(acc).strip()
				if not acc or acc in excluded:
					continue
				if self._normalize_segment(segv) == segment:
					ids.append(acc)
		finally:
			conn.close()

		with open(ids_out, "w", encoding="utf-8") as handle:
			for acc in sorted(set(ids)):
				handle.write(acc + "\n")

		return tree_out, ids_out

	def resolve_non_update_assets(self):
		if not self.mmseq_cluster_dir or not os.path.isdir(self.mmseq_cluster_dir):
			raise FileNotFoundError(f"MMseqs cluster directory not found: {self.mmseq_cluster_dir}")
		if not self.iqtree_dir or not os.path.isdir(self.iqtree_dir):
			raise FileNotFoundError(f"IQ-TREE directory not found: {self.iqtree_dir}")

		cluster_rep = self._find_first(self.mmseq_cluster_dir, pattern_suffix="_cluster_rep.fasta")
		if not cluster_rep:
			cluster_rep = self._find_first(self.mmseq_cluster_dir, pattern_suffix=".fasta")
			if cluster_rep:
				print(f"[warn] *_cluster_rep.fasta not found; falling back to {cluster_rep}")
		if not cluster_rep:
			raise FileNotFoundError(f"No FASTA found in {self.mmseq_cluster_dir}")

		tree_file = self._find_first(self.iqtree_dir, pattern_suffix=".treefile")
		if not tree_file:
			raise FileNotFoundError(f"No treefile found in {self.iqtree_dir}")

		return cluster_rep, tree_file

	def resolve_reference_id(self, cluster_rep=None):
		if cluster_rep:
			cluster_ids = self._read_ids_from_fasta(cluster_rep)
			if cluster_ids:
				ref_id = cluster_ids[0]
				if ref_id in set(self._read_ids_from_fasta(self.padded_aln)):
					return ref_id
				print(f"[warn] Reference ID '{ref_id}' is not present in {self.padded_aln}; using first MSA ID instead.")

		aln_ids = self._read_ids_from_fasta(self.padded_aln)
		if not aln_ids:
			raise ValueError(f"Could not resolve reference ID from {self.padded_aln}")
		return aln_ids[0]

	def write_ids_file(self, rel_name, ids):
		path = os.path.join(self.output_dir, rel_name)
		with open(path, "w", encoding="utf-8") as handle:
			for item in ids:
				handle.write(str(item).strip() + "\n")
		return path

	def run_command(self, command, retry_message=None):
		try:
			subprocess.run(command, check=True)
			return True
		except subprocess.CalledProcessError:
			if retry_message:
				print(retry_message)
			return False

	def build_vcf(self, ref_id, exclude_ids_file=None):
		vcf_path = os.path.join(self.output_dir, "all_samples.vcf")
		base_cmd = ["faToVcf", f"-ref={ref_id}"]
		if exclude_ids_file:
			cmd = base_cmd + [f"-excludeFile={exclude_ids_file}", self.padded_aln, vcf_path]
			if self.run_command(cmd, "[warn] faToVcf failed with -excludeFile; retrying without exclude filter."):
				return vcf_path
		cmd = base_cmd + [self.padded_aln, vcf_path]
		subprocess.run(cmd, check=True)
		return vcf_path

	def run(self):
		os.makedirs(self.output_dir, exist_ok=True)

		cluster_rep = None
		if self.update_db:
			tree_file, existing_ids_file = self.prepare_update_assets()
		else:
			cluster_rep, tree_file = self.resolve_non_update_assets()
			centroid_ids = self._read_ids_from_fasta(cluster_rep)
			self.write_ids_file("centroid_ids.txt", centroid_ids)
			self.write_ids_file("aln_ids.txt", self._read_ids_from_fasta(self.padded_aln))
			existing_ids_file = self.write_ids_file("exclude_ids.txt", centroid_ids[1:])

		ref_id = self.resolve_reference_id(cluster_rep=cluster_rep)
		aln_ids = self._read_ids_from_fasta(self.padded_aln)
		if len(aln_ids) <= 1:
			raise ValueError(f"Need at least 2 sequences in {self.padded_aln}, found {len(aln_ids)}")

		exclude_ids_file = None
		if self.update_db:
			if os.path.isfile(existing_ids_file):
				existing_ids = [x for x in self._read_text_lines(existing_ids_file) if x != ref_id]
				exclude_ids_file = self.write_ids_file("exclude_ids.txt", existing_ids)
		else:
			if os.path.isfile(existing_ids_file):
				exclude_ids = self._read_text_lines(existing_ids_file)
				exclude_count = len(exclude_ids)
				if exclude_count < (len(aln_ids) - 1):
					exclude_ids_file = existing_ids_file
				else:
					mode = " in test mode" if self.test_mode else ""
					print(f"[warn] Exclude list would remove all non-reference sequences ({exclude_count}/{len(aln_ids)}){mode}; running faToVcf without -excludeFile.")

		vcf_path = self.build_vcf(ref_id, exclude_ids_file=exclude_ids_file)

		usher_help = subprocess.run(["usher", "--help"], capture_output=True, text=True, check=False)
		usher_cmd = [
			"usher",
			"-v", vcf_path,
			"-t", tree_file,
			"-d", self.output_dir,
			"-o", os.path.join(self.output_dir, "usher.pb"),
			"-C", "-u",
		]
		if " -T " in (usher_help.stdout + usher_help.stderr):
			usher_cmd.extend(["-T", str(self.threads)])

		env = os.environ.copy()
		for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
			env[var] = str(self.threads)
		subprocess.run(usher_cmd, check=True, env=env)

	@staticmethod
	def _read_text_lines(path):
		with open(path, "r", encoding="utf-8") as handle:
			return [line.strip() for line in handle if line.strip()]


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Run USHER placement for either standard or update mode")
	parser.add_argument("--padded_aln", required=True, help="Padded alignment FASTA")
	parser.add_argument("--output_dir", required=True, help="Output directory")
	parser.add_argument("--mmseq_cluster_dir", default=None, help="MMseqs cluster directory for standard mode")
	parser.add_argument("--iqtree_dir", default=None, help="IQ-TREE directory for standard mode")
	parser.add_argument("--update_db", default=None, help="Existing update DB for update-mode tree and ID extraction")
	parser.add_argument("--threads", default=1, type=int, help="Thread count")
	parser.add_argument("--test_mode", default="0", help="Whether test mode is enabled (1/0)")
	args = parser.parse_args()

	UsherPlacement(
		padded_aln=args.padded_aln,
		output_dir=args.output_dir,
		mmseq_cluster_dir=args.mmseq_cluster_dir,
		iqtree_dir=args.iqtree_dir,
		update_db=args.update_db,
		threads=args.threads,
		test_mode=args.test_mode,
	).run()