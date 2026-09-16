#!/usr/bin/env python3
"""Build the reference backbone alignment: every reference, every shared indel.

WHY THIS EXISTS
===============

Without a backbone, each reference is aligned to the master by Nextalign and
anything the master lacks is stripped to the insertions table. Differences
between major clades then have no columns, so no tree can see them. Measured on
HCV: 35,537 reference bases (1.6%) were thrown away that way.

This step builds one alignment per segment holding all references, and
PadAlignment slots each query into its reference's row of it. Only new,
non-reference sequences still lose their private insertions to the insertions
table. See info_help/guide_alignment_and_insertions.md.

WHAT "HIGH QUALITY" MEANS HERE
==============================

1. **Whole records.** The skeleton alignment is built from full reference
   records, so the master's row is its entire record and the 5' trim that
   CalcAlignmentCord has to correct for is zero.
2. **Codon-aware genes.** A plain nucleotide aligner puts gaps anywhere. On HCV,
   MAFFT left 8 of 34 master gap runs inside genes out of frame, which shifts
   every codon read downstream of them. So each master CDS is re-aligned as
   protein and back-translated: gaps inside a gene are whole codons.
3. **The most accurate MAFFT strategy the input size allows** (L-INS-i when it
   fits the budget, FFT-NS-i otherwise), for both the skeleton and the proteins.
4. **Rare insertions filtered, never lost.** An insertion column carried by
   fewer than ``--min_insertion_support`` references widens every row in the
   database for one sequence. Those columns are dropped and their bases written
   to ``dropped_insertions.tsv`` in the insertions-table format. On the curated
   HCV alignment 326 of 523 insertion columns were single-reference.

It uses only what every NCBI build already has - the reference FASTA, the
reference list and each master's GFF - so it works for any simple genome.

OUTPUTS (in --output_dir)
=========================

``refset_<segment>_aln.fasta``
    One per master, master row first, canonical accessions. ``<segment>`` is the
    master's segment from the reference list, or ``0`` when it has none, which is
    exactly what ``projectability.find_precomputed_reference_alignment`` resolves.
``dropped_insertions.tsv``
    ``primary_accession, reference, insertion, segment`` - the bases of each
    reference that are not in its row, as ``<master position>:<bases>``.
``build_report.tsv``
    One row per segment: sizes, strategies, frame checks, what was dropped.

The step fails rather than writing a backbone that breaks one of its guarantees:
the master row must degap to the master record, no master gap run inside a
codon region may be out of frame, and every reference base must be either in
its row or recorded as dropped.
"""

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter

from Bio.Data import CodonTable

import accession_utils
import segment_utils
from ExportRefListFromUpdateDb import load_reference_file_table

#: The shared server caps any single job at 14 threads.
MAX_THREADS = 14

#: Pairwise DP cells above which L-INS-i is too slow and FFT-NS-i is used.
#: ``pairs * length^2``. 28 RABV references of 12 kb (5e10) take L-INS-i; 238
#: HCV genomes (2.6e12) do not. Proteins are three times shorter, so an HCV
#: polyprotein alignment (2.5e11) still fits.
LINSI_BUDGET_NUCLEOTIDE = 1e11
LINSI_BUDGET_PROTEIN = 3e11

#: Above this many sequences even the iterative refinement is dropped.
FFTNS2_ABOVE = 2000

#: A reading frame is accepted with at most this many internal stops per codon
#: (with a floor of one, so a single sequencing error does not demote a gene).
MAX_INTERNAL_STOP_RATE = 0.005

#: k-mer length used to assign references to masters and to orient them.
KMER = 9

_STANDARD = CodonTable.unambiguous_dna_by_id[1]
_CODON_TO_AA = dict(_STANDARD.forward_table)
for _stop in _STANDARD.stop_codons:
	_CODON_TO_AA[_stop] = "*"

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN-", "TGCAYRSWMKVHDBN-")


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def clean_sequence(text):
	"""Uppercase DNA with whitespace and gap characters removed, U read as T."""
	return "".join(text.split()).upper().replace("U", "T").replace("-", "").replace(".", "")


def reverse_complement(seq):
	return seq.translate(_COMPLEMENT)[::-1]


def translate(nucleotides):
	"""Standard-code translation. Anything not a clean codon is ``X``."""
	return "".join(
		_CODON_TO_AA.get(nucleotides[i:i + 3], "X")
		for i in range(0, len(nucleotides) - len(nucleotides) % 3, 3)
	)


def read_fasta(path):
	"""``[(canonical accession, sequence)]`` in file order, first copy of each id kept."""
	records = []
	seen = set()
	header = None
	chunks = []

	def flush():
		if header is None:
			return
		token = header.split()[0] if header.split() else ""
		acc = accession_utils.normalise_accession(token) or token
		if not acc:
			return
		if acc in seen:
			print(f"[warn] Duplicate reference {acc} in the reference FASTA; keeping the first copy.")
			return
		seen.add(acc)
		records.append((acc, clean_sequence("".join(chunks))))

	with open(path, "r", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			line = line.rstrip("\r\n")
			if line.startswith(">"):
				flush()
				header = line[1:].strip()
				chunks = []
			else:
				chunks.append(line)
	flush()
	return records


def write_fasta(path, records):
	with open(path, "w", encoding="utf-8") as handle:
		for acc, seq in records:
			handle.write(f">{acc}\n{seq}\n")


def read_aligned_fasta(path):
	"""``{id: aligned sequence}`` from a MAFFT output, uppercased."""
	out = {}
	header = None
	chunks = []
	with open(path, "r", encoding="utf-8") as handle:
		for line in handle:
			line = line.rstrip("\r\n")
			if line.startswith(">"):
				if header is not None:
					out[header] = "".join(chunks).upper()
				header = line[1:].strip().split()[0]
				chunks = []
			else:
				chunks.append(line.strip())
	if header is not None:
		out[header] = "".join(chunks).upper()
	return out


def kmers(seq, k=KMER):
	return {seq[i:i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i:i + k]}


def containment(a, b):
	if not a or not b:
		return 0.0
	return len(a & b) / min(len(a), len(b))


# ---------------------------------------------------------------------------
# Grouping references by master
# ---------------------------------------------------------------------------

def segment_label(value):
	normalised = segment_utils.normalise_segment(value)
	if normalised is None or normalised.casefold() in segment_utils.PANDAS_NULL_TOKENS:
		return ""
	return normalised


def group_references(records, ref_list, is_segmented):
	"""``[(label, master, [member accessions])]`` with the master first in each group.

	Unsegmented builds, and any build with a single master, put every reference
	in one group. Segmented builds group by the reference list's segment column;
	a reference with no usable segment goes to the master it shares most k-mers
	with, because leaving it out would orphan every query BLAST assigned to it.
	"""
	table = load_reference_file_table(ref_list)
	table["acc"] = table["primary_accession"].map(
		lambda v: accession_utils.normalise_accession(v) or str(v).strip())
	types = table["accession_type"].astype(str).str.strip().str.lower()
	sequences = dict(records)

	masters = []
	for _, row in table[types == "master"].iterrows():
		if row["acc"] in sequences and row["acc"] not in [m for m, _ in masters]:
			masters.append((row["acc"], segment_label(row.get("segment", ""))))
		elif row["acc"] not in sequences:
			print(f"[warn] Master {row['acc']} has no sequence in the reference FASTA; skipped.")
	if not masters:
		raise ValueError("No master accession from the reference list has a sequence in the reference FASTA")

	if str(is_segmented).upper() != "Y" or len(masters) == 1:
		master, label = masters[0]
		if len(masters) > 1:
			print(f"[warn] {len(masters)} masters in an unsegmented build; using {master} as the backbone master.")
		members = [master] + [acc for acc, _ in records if acc != master]
		return [(label or "0", master, members)]

	by_label = {}
	for master, label in masters:
		if not label:
			raise ValueError(f"Segmented build: master {master} has no segment in the reference list")
		if label in by_label:
			print(f"[warn] Masters {by_label[label]} and {master} share segment {label}; keeping {by_label[label]}.")
			continue
		by_label[label] = master

	ref_segment = {}
	for _, row in table.iterrows():
		label = segment_label(row.get("segment", ""))
		if label in by_label:
			ref_segment.setdefault(row["acc"], label)

	groups = {label: [master] for label, master in by_label.items()}
	master_kmers = None
	unassigned = []
	for acc, seq in records:
		if acc in by_label.values():
			continue
		label = ref_segment.get(acc)
		if label is None:
			if master_kmers is None:
				master_kmers = {lab: kmers(sequences[m]) for lab, m in by_label.items()}
			query = kmers(seq)
			rc = kmers(reverse_complement(seq))
			label = max(by_label, key=lambda lab: max(containment(query, master_kmers[lab]),
												   containment(rc, master_kmers[lab])))
			unassigned.append((acc, label))
		groups[label].append(acc)
	if unassigned:
		preview = ", ".join(f"{a}->{l}" for a, l in unassigned[:10])
		print(f"[info] {len(unassigned)} reference(s) had no segment in the reference list and were "
			  f"assigned to the most similar master: {preview}")
	return [(label, by_label[label], groups[label]) for label in sorted(by_label, key=str)]


def orient_to_master(records, master):
	"""Reverse-complement any reference that matches its master better that way."""
	master_kmers = kmers(records[master])
	flipped = []
	for acc, seq in list(records.items()):
		if acc == master:
			continue
		forward = containment(kmers(seq), master_kmers)
		reverse = containment(kmers(reverse_complement(seq)), master_kmers)
		if reverse > forward and reverse > 0.05:
			records[acc] = reverse_complement(seq)
			flipped.append(acc)
	if flipped:
		print(f"[warn] Reverse-complemented {len(flipped)} reference(s) to match master {master}: "
			  f"{', '.join(flipped[:10])}. PadAlignment cannot project their queries, which are "
			  f"aligned in the reference's original orientation.")
	return flipped


# ---------------------------------------------------------------------------
# Master annotation
# ---------------------------------------------------------------------------

def find_gff(gff_paths, master):
	for path in gff_paths or []:
		stem = os.path.basename(path)
		if accession_utils.accession_from_filename(stem) == master or stem.startswith(master):
			return path
	return None


def read_cds(gff_path):
	"""CDS lines of a GFF3: ``[{start, end, strand, phase, product}]`` and the region length."""
	cds = []
	region_length = None
	with open(gff_path, "r", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			if line.startswith("#") or not line.strip():
				continue
			parts = line.rstrip("\r\n").split("\t")
			if len(parts) < 9:
				continue
			feature = parts[2]
			try:
				start, end = int(parts[3]), int(parts[4])
			except ValueError:
				continue
			if feature == "region" and region_length is None:
				region_length = end
			if feature != "CDS":
				continue
			product = ""
			for attribute in parts[8].split(";"):
				if attribute.startswith("product="):
					product = attribute.split("=", 1)[1]
			cds.append({"start": start, "end": end, "strand": parts[6], "phase": parts[7], "product": product})
	return cds, region_length


def choose_codon_regions(cds, master_seq):
	"""Non-overlapping, plus-strand CDS that translate cleanly from the master.

	Longest first, so a polyprotein beats the short frameshift product nested in
	it and influenza M1 beats the spliced M2 exons. Everything else - minus strand,
	wrapped, not a whole number of codons, internal stops in the master - stays
	nucleotide-aligned rather than being forced into a frame it does not have.
	"""
	usable = []
	skipped = []
	for entry in cds:
		start, end = entry["start"], entry["end"]
		reason = None
		if entry["strand"] != "+":
			reason = "not plus strand"
		elif start > end or start < 1 or end > len(master_seq):
			reason = "outside the master record"
		elif entry["phase"] not in ("0", "."):
			reason = f"phase {entry['phase']}"
		elif (end - start + 1) % 3:
			reason = "not a whole number of codons"
		else:
			protein = translate(master_seq[start - 1:end])
			if "*" in protein[:-1]:
				reason = "internal stop in the master"
		if reason:
			skipped.append((entry, reason))
		else:
			usable.append(entry)

	chosen = []
	for entry in sorted(usable, key=lambda e: (-(e["end"] - e["start"]), e["start"])):
		if all(entry["end"] < c["start"] or entry["start"] > c["end"] for c in chosen):
			chosen.append(entry)
	return sorted(chosen, key=lambda e: e["start"]), skipped


# ---------------------------------------------------------------------------
# MAFFT
# ---------------------------------------------------------------------------

def mafft_strategy(n_sequences, max_length, protein, override="auto"):
	"""``(label, mafft arguments)`` for the most accurate strategy that fits the budget."""
	strategies = {
		"linsi": ("L-INS-i", ["--localpair", "--maxiterate", "1000"]),
		"fftnsi": ("FFT-NS-i", ["--retree", "2", "--maxiterate", "2"]),
		"fftns2": ("FFT-NS-2", ["--retree", "2", "--maxiterate", "0"]),
	}
	if override != "auto":
		return strategies[override]
	pairs = n_sequences * (n_sequences - 1) / 2
	budget = LINSI_BUDGET_PROTEIN if protein else LINSI_BUDGET_NUCLEOTIDE
	if pairs * max_length ** 2 <= budget:
		return strategies["linsi"]
	if n_sequences <= FFTNS2_ABOVE:
		return strategies["fftnsi"]
	return strategies["fftns2"]


def run_mafft(arguments, input_path, output_path, threads, log_path):
	command = ["mafft", "--thread", str(threads)] + arguments + [input_path]
	with open(output_path, "w", encoding="utf-8") as out, open(log_path, "a", encoding="utf-8") as log:
		log.write("$ " + " ".join(command) + "\n")
		log.flush()
		try:
			completed = subprocess.run(command, stdout=out, stderr=log, check=False)
		except FileNotFoundError as exc:
			raise RuntimeError("mafft is not on PATH; activate the vgtk environment") from exc
	if completed.returncode != 0:
		raise RuntimeError(f"mafft failed (exit {completed.returncode}); see {log_path}")


# ---------------------------------------------------------------------------
# Codon regions
# ---------------------------------------------------------------------------

def best_frame(nucleotides):
	"""``(offset, codon count, protein, internal stops)`` for the least-broken frame."""
	best = None
	for offset in (0, 1, 2):
		n_codons = (len(nucleotides) - offset) // 3
		if n_codons <= 0:
			continue
		protein = translate(nucleotides[offset:offset + 3 * n_codons])
		stops = protein[:-1].count("*")
		key = (stops, offset)
		if best is None or key < best[0]:
			best = (key, (offset, n_codons, protein, stops))
	if best is None:
		return 0, 0, "", 0
	return best[1]


def frame_is_usable(n_codons, stops, master_codons):
	if n_codons < max(10, master_codons // 10):
		return False
	return stops <= max(1, int(n_codons * MAX_INTERNAL_STOP_RATE))


def back_translate(aligned_protein, coding_nucleotides):
	pieces = []
	index = 0
	for residue in aligned_protein:
		if residue == "-":
			pieces.append("---")
		else:
			pieces.append(coding_nucleotides[index:index + 3])
			index += 3
	if index != len(coding_nucleotides):
		raise ValueError("Back-translation consumed a different number of codons than the protein holds")
	return "".join(pieces)


def tuck_partial_codons(left, block, right):
	"""Move the 0-2 bases outside a row's frame into the gap codon beside its own codons.

	Returns ``(left, block, right)``. A partial genome that stops mid-codon leaves
	1-2 bases after its last whole codon. Parked in the columns at the block's
	edge, they sit hundreds of columns from the rest of the row, usually in
	master-gap columns: the row's aligned span then ends far past the last
	master base it really covers, which is exactly what
	ValidateDbTree's feature projection check measured on HCV_test (13 partial
	references). Placed in the empty codon slot next to the row's last (or
	first) codon they stay adjacent, and the column filter still treats that
	slot as one whole codon. Only when the row's codons already reach the edge of
	the block do the bases stay in the edge columns, where they are adjacent anyway.
	"""
	residues = [i for i, char in enumerate(block) if char != "-"]
	if not residues:
		return left, block, right
	if right:
		end = (residues[-1] // 3 + 1) * 3
		if end + 3 <= len(block) and block[end:end + 3] == "---":
			block = block[:end] + right.ljust(3, "-") + block[end + 3:]
			right = ""
	if left:
		start = (residues[0] // 3) * 3
		if start >= 3 and block[start - 3:start] == "---":
			block = block[:start - 3] + left.rjust(3, "-") + block[start:]
			left = ""
	return left, block, right


def align_codon_region(slices, master, work_dir, threads, log_path, protein_override="auto"):
	"""Codon alignment of one gene from each row's nucleotide slice.

	Returns ``(pieces, broken)`` where ``pieces`` is ``{acc: (left, block, right)}``:
	``block`` is the codon alignment, ``left``/``right`` the 0-2 bases outside the
	chosen frame, which the caller places in their own columns.
	"""
	master_codons = len(slices[master]) // 3
	good = {}
	broken = []
	for acc, nucleotides in slices.items():
		if not nucleotides:
			continue
		offset, n_codons, protein, stops = best_frame(nucleotides)
		if acc == master:
			offset, n_codons, protein = 0, master_codons, translate(nucleotides)
		elif not frame_is_usable(n_codons, stops, master_codons):
			broken.append(acc)
			continue
		good[acc] = (offset, n_codons, protein)

	protein_in = os.path.join(work_dir, "protein.fa")
	protein_out = os.path.join(work_dir, "protein.aln.fa")
	write_fasta(protein_in, [(acc, protein.replace("*", "X")) for acc, (_, _, protein) in good.items()])
	if len(good) > 1:
		label, arguments = mafft_strategy(len(good), max(len(p) for _, _, p in good.values()), True, protein_override)
		run_mafft(["--amino", "--anysymbol", "--quiet"] + arguments, protein_in, protein_out, threads, log_path)
		aligned = read_aligned_fasta(protein_out)
	else:
		label = "single sequence"
		aligned = {acc: protein.replace("*", "X") for acc, (_, _, protein) in good.items()}

	pieces = {}
	for acc, (offset, n_codons, _) in good.items():
		nucleotides = slices[acc]
		coding = nucleotides[offset:offset + 3 * n_codons]
		pieces[acc] = tuck_partial_codons(nucleotides[:offset], back_translate(aligned[acc], coding),
										  nucleotides[offset + 3 * n_codons:])

	if broken:
		block_in = os.path.join(work_dir, "codon_block.fa")
		add_in = os.path.join(work_dir, "broken.fa")
		add_out = os.path.join(work_dir, "codon_block_added.fa")
		write_fasta(block_in, [(acc, block) for acc, (_, block, _) in pieces.items()])
		write_fasta(add_in, [(acc, slices[acc]) for acc in broken])
		run_mafft(["--nuc", "--quiet", "--keeplength", "--add", add_in], block_in, add_out, threads, log_path)
		added = read_aligned_fasta(add_out)
		for acc in broken:
			pieces[acc] = ("", added[acc], "")
	return pieces, broken, label


# ---------------------------------------------------------------------------
# Assembly, filtering and accounting
# ---------------------------------------------------------------------------

def master_columns(master_row):
	"""0-based column of each 1-based master position."""
	positions = {}
	count = 0
	for column, char in enumerate(master_row):
		if char != "-":
			count += 1
			positions[count] = column
	return positions


def assemble(skeleton, master, regions, codon_pieces):
	"""Stitch skeleton noncoding columns and codon regions into one alignment.

	Returns ``(rows, units)`` where ``units`` lists ``(first column, width, kind)``
	for the filter: ``kind`` is ``"codon"`` for a whole-codon unit, else
	``"column"``. Keeping codon units whole is what stops the filter from ever
	putting a master gap run out of frame.
	"""
	accessions = list(skeleton)
	columns = master_columns(skeleton[master])
	rows = {acc: [] for acc in accessions}
	units = []
	width = 0

	def add_columns(texts, kind):
		nonlocal width
		span = len(next(iter(texts.values()))) if texts else 0
		if span == 0:
			return
		for acc in accessions:
			rows[acc].append(texts[acc])
		step = 3 if kind == "codon" else 1
		for offset in range(0, span, step):
			units.append((width + offset, step, kind))
		width += span

	cursor = 0
	for region, pieces in zip(regions, codon_pieces):
		first = columns[region["start"]]
		last = columns[region["end"]]
		add_columns({acc: skeleton[acc][cursor:first] for acc in accessions}, "column")
		left_width = max(len(pieces.get(acc, ("", "", ""))[0]) for acc in accessions)
		right_width = max(len(pieces.get(acc, ("", "", ""))[2]) for acc in accessions)
		block_width = len(pieces[master][1])
		add_columns({acc: pieces.get(acc, ("", "", ""))[0].rjust(left_width, "-") for acc in accessions}, "column")
		add_columns({acc: pieces.get(acc, ("", "-" * block_width, ""))[1] for acc in accessions}, "codon")
		add_columns({acc: pieces.get(acc, ("", "", ""))[2].ljust(right_width, "-") for acc in accessions}, "column")
		cursor = last + 1
	add_columns({acc: skeleton[acc][cursor:] for acc in accessions}, "column")
	return {acc: "".join(parts) for acc, parts in rows.items()}, units


def filter_columns(rows, master, units, min_support):
	"""Drop all-gap units, and insertion units carried by too few references.

	An insertion unit is one where the master is entirely gap. Its support is the
	number of other rows with at least one base in it.
	"""
	accessions = [acc for acc in rows if acc != master]
	keep = []
	dropped_all_gap = 0
	dropped_low_support = 0
	master_row = rows[master]
	for first, width, _ in units:
		if any(master_row[c] != "-" for c in range(first, first + width)):
			keep.append((first, width))
			continue
		support = sum(1 for acc in accessions
					  if any(rows[acc][c] != "-" for c in range(first, first + width)))
		if support == 0:
			dropped_all_gap += 1
		elif support < min_support:
			dropped_low_support += 1
		else:
			keep.append((first, width))
	kept_columns = [c for first, width in keep for c in range(first, first + width)]
	filtered = {acc: "".join(row[c] for c in kept_columns) for acc, row in rows.items()}
	return filtered, dropped_all_gap, dropped_low_support


def removed_bases(record, row):
	"""Runs of ``record`` missing from ``row``: ``[(row column after which, bases)]``.

	The row's bases are matched to the record greedily from the left; a row that
	is not a subsequence of its record raises, because that means the backbone
	holds bases the reference does not have.
	"""
	runs = []
	record_index = 0
	last_column = -1
	pending = []
	for column, char in enumerate(row):
		if char == "-":
			continue
		while record_index < len(record) and record[record_index] != char:
			pending.append(record[record_index])
			record_index += 1
		if record_index >= len(record):
			raise ValueError("aligned row is not a subsequence of its record")
		if pending:
			runs.append((last_column, "".join(pending)))
			pending = []
		record_index += 1
		last_column = column
	if record_index < len(record):
		runs.append((last_column, record[record_index:]))
	return runs


def master_gap_runs_out_of_frame(master_row, regions):
	"""Gap runs inside a codon region whose length is not a whole number of codons."""
	columns = master_columns(master_row)
	bad = []
	for region in regions:
		segment = master_row[columns[region["start"]]:columns[region["end"]] + 1]
		run = 0
		for char in segment + "X":
			if char == "-":
				run += 1
			else:
				if run % 3:
					bad.append((region.get("product", ""), run))
				run = 0
	return bad


# ---------------------------------------------------------------------------
# One segment
# ---------------------------------------------------------------------------

def build_segment(label, master, members, sequences, gff_path, work_dir, threads,
				  min_support, skeleton_override="auto", protein_override="auto"):
	started = time.time()
	log_path = os.path.join(work_dir, "mafft.log")
	records = {acc: sequences[acc] for acc in members}
	flipped = orient_to_master(records, master)
	report = {
		"segment": label, "master": master, "references": len(records) - 1,
		"reverse_complemented": ",".join(flipped),
	}

	skeleton_in = os.path.join(work_dir, "skeleton.fa")
	skeleton_out = os.path.join(work_dir, "skeleton.aln.fa")
	write_fasta(skeleton_in, list(records.items()))
	if len(records) > 1:
		strategy, arguments = mafft_strategy(len(records), max(map(len, records.values())), False, skeleton_override)
		run_mafft(["--nuc", "--quiet"] + arguments, skeleton_in, skeleton_out, threads, log_path)
		skeleton = read_aligned_fasta(skeleton_out)
	else:
		strategy = "single sequence"
		skeleton = dict(records)
	report["skeleton_strategy"] = strategy
	skeleton = {acc: skeleton[acc] for acc in records}  # master first, input order

	regions = []
	if gff_path:
		cds, region_length = read_cds(gff_path)
		if region_length and region_length != len(records[master]):
			print(f"[warn] GFF {os.path.basename(gff_path)} describes {region_length} bases but master "
				  f"{master} has {len(records[master])}; genes aligned as nucleotides only.")
		else:
			regions, skipped = choose_codon_regions(cds, records[master])
			for entry, reason in skipped:
				print(f"[info] {master} CDS {entry['start']}-{entry['end']} ({entry['product'] or 'unnamed'}) "
					  f"not codon-aligned: {reason}.")
	else:
		print(f"[warn] No GFF for master {master}; segment {label} aligned as nucleotides only.")

	codon_pieces = []
	broken_total = Counter()
	protein_labels = set()
	columns = master_columns(skeleton[master])
	for index, region in enumerate(regions):
		first, last = columns[region["start"]], columns[region["end"]]
		slices = {acc: row[first:last + 1].replace("-", "") for acc, row in skeleton.items()}
		region_dir = os.path.join(work_dir, f"cds_{index + 1}")
		os.makedirs(region_dir, exist_ok=True)
		pieces, broken, protein_label = align_codon_region(slices, master, region_dir, threads, log_path, protein_override)
		protein_labels.add(protein_label)
		broken_total.update(broken)
		codon_pieces.append(pieces)
		if broken:
			print(f"[info] {master} {region['product'] or 'CDS'}: {len(broken)} reference(s) with a broken "
				  f"reading frame added by nucleotide ({', '.join(broken[:5])}).")

	rows, units = assemble(skeleton, master, regions, codon_pieces)
	rows, dropped_all_gap, dropped_low_support = filter_columns(rows, master, units, min_support)

	if rows[master].replace("-", "") != records[master]:
		raise AssertionError(f"Master {master} row does not degap to its record")
	out_of_frame = master_gap_runs_out_of_frame(rows[master], regions)
	if out_of_frame:
		raise AssertionError(f"Master {master} has gap runs out of frame inside genes: {out_of_frame[:5]}")
	width = len(rows[master])
	if any(len(row) != width for row in rows.values()):
		raise AssertionError("Backbone rows are not all the same width")

	positions = []
	count = 0
	for char in rows[master]:
		if char != "-":
			count += 1
		positions.append(count)
	insertions = []
	recorded_bases = 0
	for acc, row in rows.items():
		if acc == master:
			continue
		runs = removed_bases(records[acc], row)
		if not runs:
			continue
		recorded_bases += sum(len(bases) for _, bases in runs)
		text = ";".join(f"{positions[col] if col >= 0 else 0}:{bases}" for col, bases in runs)
		insertions.append((acc, master, text))

	report.update({
		"protein_strategy": ",".join(sorted(protein_labels)) or "none",
		"width": width,
		"master_length": len(records[master]),
		"master_gap_columns": rows[master].count("-"),
		"codon_regions": len(regions),
		"broken_frame_references": len(broken_total),
		"insertion_units_dropped_low_support": dropped_low_support,
		"all_gap_units_dropped": dropped_all_gap,
		"reference_bases_recorded_as_insertions": recorded_bases,
		"seconds": round(time.time() - started, 1),
	})
	ordered = [(master, rows[master])] + [(acc, row) for acc, row in rows.items() if acc != master]
	return ordered, insertions, report


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

REPORT_COLUMNS = [
	"segment", "master", "references", "width", "master_length", "master_gap_columns",
	"codon_regions", "skeleton_strategy", "protein_strategy", "broken_frame_references",
	"insertion_units_dropped_low_support", "all_gap_units_dropped",
	"reference_bases_recorded_as_insertions", "reverse_complemented", "seconds",
]


def build(ref_fasta, ref_list, gff_paths, output_dir, threads=4, min_insertion_support=2,
		  is_segmented="N", keep_work=False, skeleton_strategy="auto", protein_strategy="auto"):
	threads = max(1, min(int(threads), MAX_THREADS))
	if min_insertion_support < 1:
		raise ValueError("--min_insertion_support must be at least 1")
	if shutil.which("mafft") is None:
		raise RuntimeError("mafft is not on PATH; activate the vgtk environment")

	records = read_fasta(ref_fasta)
	if not records:
		raise ValueError(f"No sequences in {ref_fasta}")
	sequences = dict(records)
	groups = group_references(records, ref_list, is_segmented)

	os.makedirs(output_dir, exist_ok=True)
	work_root = tempfile.mkdtemp(prefix="_work_", dir=output_dir)
	reports = []
	all_insertions = []
	try:
		for label, master, members in groups:
			print(f"[info] Segment {label}: master {master}, {len(members) - 1} reference(s).")
			segment_dir = os.path.join(work_root, f"segment_{label}")
			os.makedirs(segment_dir, exist_ok=True)
			rows, insertions, report = build_segment(
				label, master, members, sequences, find_gff(gff_paths, master), segment_dir,
				threads, min_insertion_support, skeleton_strategy, protein_strategy)
			write_fasta(os.path.join(output_dir, f"refset_{label}_aln.fasta"), rows)
			segment_value = label if str(is_segmented).upper() == "Y" else ""
			all_insertions.extend((acc, ref, text, segment_value) for acc, ref, text in insertions)
			reports.append(report)
			print("[info] " + ", ".join(f"{k}={report.get(k, '')}" for k in REPORT_COLUMNS))
	finally:
		if keep_work:
			print(f"[info] Intermediate files kept in {work_root}")
		else:
			shutil.rmtree(work_root, ignore_errors=True)

	with open(os.path.join(output_dir, "dropped_insertions.tsv"), "w", newline="", encoding="utf-8") as handle:
		writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
		writer.writerow(["primary_accession", "reference", "insertion", "segment"])
		writer.writerows(all_insertions)
	with open(os.path.join(output_dir, "build_report.tsv"), "w", newline="", encoding="utf-8") as handle:
		writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t", lineterminator="\n")
		writer.writeheader()
		writer.writerows(reports)
	return reports


def main(argv=None):
	parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
	parser.add_argument("--ref_fasta", required=True, help="Reference and master sequences (BLAST_ALIGNMENT's ref_seq_filtered.fa)")
	parser.add_argument("--ref_list", required=True, help="Reference list TSV (accession, type, segment, ...)")
	parser.add_argument("--gff", nargs="*", default=[], help="Master GFF3 file(s), named by accession")
	parser.add_argument("--output_dir", default="ref_set_aligned", help="Where refset_<segment>_aln.fasta files are written")
	parser.add_argument("--threads", type=int, default=4, help=f"MAFFT threads (capped at {MAX_THREADS})")
	parser.add_argument("--min_insertion_support", type=int, default=2,
						help="Keep an insertion column only if at least this many references have bases in it. "
							 "1 keeps every insertion. Dropped bases go to dropped_insertions.tsv.")
	parser.add_argument("--is_segmented", default="N", help="Y for a segmented virus")
	parser.add_argument("--skeleton_strategy", default="auto", choices=["auto", "linsi", "fftnsi", "fftns2"],
						help="Force the nucleotide MAFFT strategy (default: most accurate that fits)")
	parser.add_argument("--protein_strategy", default="auto", choices=["auto", "linsi", "fftnsi", "fftns2"],
						help="Force the protein MAFFT strategy (default: most accurate that fits)")
	parser.add_argument("--keep_work", action="store_true", help="Keep intermediate alignments for debugging")
	args = parser.parse_args(argv)
	try:
		build(args.ref_fasta, args.ref_list, args.gff, args.output_dir, args.threads,
			  args.min_insertion_support, args.is_segmented, args.keep_work,
			  args.skeleton_strategy, args.protein_strategy)
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		return 2
	return 0


if __name__ == "__main__":
	sys.exit(main())
