"""Project each master's annotated CDS features onto every row of its alignment.

THREE COORDINATE SPACES
=======================

Almost every defect this module has shipped has been one space being handed to
something expecting another, so they are named here and referred to by number
throughout:

  1. **Master genome coordinates.** 1-based positions along the master's own
     ungapped genome. A GFF is written in this space, so ``cds['start']`` and
     ``cds['end']`` are space 1, and so is the covered span returned by
     ``CalculateGenomeCoordinates`` (it reports *master residue numbers*, not
     columns - see CalcGenomeCords.py).

  2. **Alignment columns.** 1-based column numbers in the padded MSA. Every row
     of one alignment file shares this space, and that is the only reason a
     master's annotation can be transferred to a query at all.
     ``master_coord_to_aln_pos`` is the sole bridge from space 1 to space 2.
     ``get_gap_ranges`` and ``count_gaps_before_position`` speak space 2 and
     nothing else.

  3. **Per-sequence ungapped coordinates.** 1-based positions along one row's
     own residues with its gaps removed. Every row has its own space 3, and it
     coincides with space 2 only for a row carrying no gaps before the position
     being converted.

Space 1 and space 2 are the same number only while the master's own row is
ungapped. The moment a guide alignment introduces a column the master lacks
(``ref_set_aligned``; see info_help/guide_alignment_and_insertions.md) they
diverge by the number of master gaps, and any arithmetic that confuses them is
silently wrong by exactly that much.

WHAT EACH OUTPUT COLUMN MEANS
=============================

The column names are not a reliable guide to the space - ``aln_start`` is a
master coordinate while ``cds_start`` is an alignment column - so:

  ``aln_start`` / ``aln_end``
      Space 1. The master-coordinate span this row actually covers, used to
      clamp features to the part of the genome a partial sequence reached.

  ``cds_start`` / ``cds_end``
      Space 2. The feature's alignment columns. Identical for every row of one
      alignment, because they describe the master's feature, not the row's.

  ``cds_start_OG_seq`` / ``cds_end_OG_seq``
      Space 3. The same feature expressed in this row's own numbering, so that
      slicing the row's raw (unaligned) sequence between them yields the CDS.

THE 5' TRIM
===========

There is a fourth thing that is emphatically NOT a coordinate space: the master's
alignment row need not begin at the first base of the master's record. The row
comes from the guide alignment, whose curated copy of a reference is often
trimmed - measured on the influenza build, by 27, 10, 24, 69, 45, 25 and 26 bases
on segments 1, 2, 3, 4, 5, 7 and 8, with segment 6 alone flush.

The GFF is written in RECORD coordinates. Counting residues down the alignment
row numbers them from the start of the row. Those two agree only when the trim is
zero, so ``compute_trim_offset`` locates the row inside the record and the map is
keyed by ``residue_index + trim``. The covered span from
``CalculateGenomeCoordinates`` is shifted by the same amount, because it counts
master residues down the row as well and the span clamp has to compare like with
like.

Left uncorrected this produced a row mixing two spaces: ``og_start`` was the
record coordinate (it survived as a pass-through, the conversion being a no-op
when the master's gaps all sit downstream) while ``og_end`` was clamped to the
row's length - truncating the CDS by the trim and dropping the true stop codon.

THE INVARIANT THAT MATTERS DOWNSTREAM
=====================================

``cds_*_OG_seq`` is in the row's own numbering (space 3), which is the numbering
of the alignment row as stored in ``sequence_alignment`` - the TRIMMED row. It is
NOT the record coordinate whenever the trim is non-zero. For a flush master
(trim 0) the two coincide, and space 3 equals space 1: the master's ungapped
residue numbering *is* the master genome, so its ``cds_*_OG_seq`` comes back out
as its plain GFF coordinate whatever gaps its row carries.

``AnnotateMutations.resolve_master_coordinate_space`` depends on precisely that.
It builds ``build_alignment_coordinate_map(master_alignment)``, which is keyed by
the master's ungapped position (space 3 == space 1), and looks up
``cds_start + (aa_pos - 1) * 3`` in it to find a catalogue residue's columns. Feed
it a master ``cds_start`` that is short by the master's gap count and every
catalogue position resolves to the wrong codon - with no error, because the
lookup still succeeds.

THE TWO DOWNSTREAM READERS DISAGREE ON PURPOSE
==============================================

They read different columns under different conventions, and both are right:

  ``AnnotateMutations`` prefers ``cds_start_OG_seq``/``cds_end_OG_seq`` (space 3)
  and resolves them through an ungapped-keyed coordinate map. For the master's
  own row that is space 1, per the invariant above.

  ``protein_alleles`` reads ``cds_start``/``cds_end`` (space 2) and treats them
  as columns outright (``codon_columns(..., coord_space='column')``). Its own
  note records how that was discovered: "the HA master's `cds_end` is 1782,
  which is exactly the alignment width, while its ungapped length is 1650", and
  ``detect_coord_space`` exists to tell the two readings apart by translating
  the master both ways and preferring the one with no internal stops.

So a change to the ``OG`` columns must leave ``cds_start``/``cds_end`` alone,
and vice versa. They are not two spellings of one number.

KNOWN RESIDUAL LIMITATION
=========================

``count_gaps_before_position`` clamps backwards at BOTH ends. For ``og_end``
that is right - the last real base at or before the feature end. For
``og_start`` it is not: a CDS start falling inside one of this row's gaps
resolves to the last residue *before* the CDS rather than the first one inside
it, and a row whose leading padding covers the start column can still produce
``og_start = 0``. In the pipeline the span clamp in
``recalculate_cds_coordinates_with_span`` hides the leading-padding case, because
the span begins at the first master coordinate the row actually covers; an
*internal* deletion spanning a CDS start is not covered by that clamp. Tracked
as xfails in tests/unit/test_alignment_coordinate_integrity.py.
"""

import os
import re
import sys
import sqlite3
import pandas as pd
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
from GffToDictionary import GffDictionary
from CalcGenomeCords import CalculateGenomeCoordinates 
from ExportRefListFromUpdateDb import load_master_accessions_from_file, load_reference_file_table
import segment_utils

def alignment_covered_columns(sequence):
	"""``(first, last)`` 1-based ALIGNMENT COLUMNS carrying a residue, or None.

	Everything outside this is padding added to square the alignment up: columns
	the submitter never reported, not sequence missing from the virus.
	"""
	first = last = None
	for index, base in enumerate(sequence, start=1):
		if base != "-":
			if first is None:
				first = index
			last = index
	if first is None:
		return None
	return (first, last)


class CalculateAlignmentCoordinates:

	def __init__(self, paded_alignment, master_gff, tmp_dir, output_dir, output_file, master_accession, blast_uniq_hits, update_db=None, update_scope_tsv=None, segment_map_tsv=None, master_seq_dir=None):
		self.paded_alignment = paded_alignment
		self.master_gff = master_gff
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.output_file = output_file
		self.master_accession = master_accession
		self.blast_uniq_hits = blast_uniq_hits
		self.update_db = update_db
		self.update_scope_tsv = update_scope_tsv
		self.segment_map_tsv = segment_map_tsv
		self.master_seq_dir = master_seq_dir
		self._master_records = None
		self._trim_warned = set()

	def load_master_records(self):
		"""``{accession: full record sequence}`` for every master, if supplied.

		Indexed by FASTA record id rather than by filename, because the only
		contract on this directory is that NextalignAlignment iterates whatever
		files are in it.
		"""
		if self._master_records is not None:
			return self._master_records
		records = {}
		if self.master_seq_dir and os.path.isdir(self.master_seq_dir):
			for name in sorted(os.listdir(self.master_seq_dir)):
				path = join(self.master_seq_dir, name)
				if not os.path.isfile(path):
					continue
				try:
					for record in SeqIO.parse(path, "fasta"):
						records.setdefault(str(record.id).strip(), str(record.seq).replace("-", "").upper())
				except Exception as exc:
					print(f"[warn] Could not read master sequence file {path}: {exc}", file=sys.stderr)
		self._master_records = records
		return records

	def compute_trim_offset(self, master, master_alignment, gff_dict=None):
		"""How far into the master's RECORD its alignment row begins, in bases.

		The master's GFF is written in full-record coordinates, but the row the
		master occupies in the merged MSA comes from the guide alignment, whose
		copy of that reference is frequently trimmed at the 5' end - measured on
		the influenza build, by 27, 10, 24, 69, 45, 25 and 26 bases on segments
		1, 2, 3, 4, 5, 7 and 8 (segment 6 alone is flush). Counting residues in
		the trimmed row therefore numbers them differently from the GFF, and
		every annotated coordinate lands that many bases out.

		Returns the offset to ADD to a residue index of the alignment row to get
		the record coordinate. 0 when the row is flush with the record, when no
		master sequence was supplied, or when the row cannot be located - the
		last of which is reported rather than assumed.
		"""
		degapped = master_alignment.replace("-", "").upper()
		if not degapped:
			return 0

		record = self.load_master_records().get(master)
		if record:
			index = record.find(degapped)
			if index >= 0:
				return index
			# Not contiguous: nextalign strips insertions relative to the
			# reference, so the row can be a subsequence rather than a substring.
			# Anchor on the 5' end, which is all the offset depends on.
			for probe in (60, 40, 30, 20):
				if len(degapped) >= probe:
					index = record.find(degapped[:probe])
					if index >= 0:
						return index
			if master not in self._trim_warned:
				self._trim_warned.add(master)
				print(f"[warn] Could not locate the alignment row of master {master} within its "
					  f"record ({len(degapped)} aligned bases vs {len(record)} in the record); "
					  f"assuming no 5' trim, so CDS coordinates may be offset.", file=sys.stderr)
			return 0

		# No master sequence to compare against. The GFF's region line still
		# says how long the record is, so a mismatch can at least be reported.
		if gff_dict and gff_dict.get('region'):
			try:
				record_length = int(gff_dict['region'][0]['end'])
			except (KeyError, IndexError, TypeError, ValueError):
				record_length = None
			if record_length and record_length != len(degapped) and master not in self._trim_warned:
				self._trim_warned.add(master)
				print(f"[warn] Master {master} has {len(degapped)} bases in the alignment but its GFF "
					  f"describes a record of {record_length}. Without --master_seq_dir the 5' trim "
					  f"cannot be located and CDS coordinates may be offset.", file=sys.stderr)
		return 0

	@staticmethod
	def _normalize_segment_value(value):
		"""Delegates to :mod:`segment_utils` - the single normalisation authority.

		This used to scrape every digit out of the string, which turned ``4.0`` into
		segment 40 and, worse, inverted the polymerase segments: ``PB2`` (segment 1)
		became ``2`` and ``PB1`` (segment 2) became ``1``. Missing still maps to
		``""`` here, because callers of this particular helper depend on it.
		"""
		normalised = segment_utils.normalise_segment(value)
		if normalised is None or normalised.casefold() in segment_utils.PANDAS_NULL_TOKENS:
			return ""
		return normalised

	@staticmethod
	def _infer_segment_from_alignment_name(fasta_file):
		name = os.path.basename(fasta_file)
		patterns = [
			r"(?:^|[_-])refset[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])segment[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])seg[_-]?(\d+)(?:$|[_-])",
		]
		for pattern in patterns:
			match = re.search(pattern, name, flags=re.IGNORECASE)
			if match:
				return match.group(1)
		return ""

	def load_existing_feature_accessions(self):
		if not self.update_db:
			return set()
		if not os.path.isfile(self.update_db):
			raise FileNotFoundError(f"Update DB not found: {self.update_db}")
		conn = sqlite3.connect(self.update_db)
		try:
			row = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='features'").fetchone()
			if row is None:
				return set()
			cols = [r[1] for r in conn.execute("PRAGMA table_info(features)").fetchall()]
			if "accession" not in cols:
				raise ValueError("Update DB features table is missing required column: accession")
			df = pd.read_sql_query("SELECT accession FROM features WHERE accession IS NOT NULL", conn)
			return set(df["accession"].astype(str).str.strip().tolist())
		finally:
			conn.close()

	def load_update_scope_accessions(self):
		if not self.update_scope_tsv or not os.path.isfile(self.update_scope_tsv):
			return set()
		df = pd.read_csv(self.update_scope_tsv, sep='\t', dtype=str)
		if 'primary_accession' not in df.columns:
			raise ValueError(f"Update scope TSV is missing required column: primary_accession ({self.update_scope_tsv})")
		return set(df['primary_accession'].fillna('').astype(str).str.strip().tolist())

	def load_segment_map(self):
		if not self.segment_map_tsv or not os.path.isfile(self.segment_map_tsv):
			return {}
		df = pd.read_csv(self.segment_map_tsv, sep='\t', dtype=str)
		if 'primary_accession' not in df.columns or 'segment' not in df.columns:
			raise ValueError(f"Segment map TSV is missing required columns: primary_accession, segment ({self.segment_map_tsv})")
		df['primary_accession'] = df['primary_accession'].fillna('').astype(str).str.strip()
		df['segment'] = df['segment'].fillna('').astype(str).str.strip().map(self._normalize_segment_value)
		return dict(zip(df['primary_accession'], df['segment']))

	def get_master_list(self):
		if os.path.isfile(self.master_accession):
			try:
				return load_master_accessions_from_file(self.master_accession)
			except:
				return []
		else:
			return [x.strip() for x in self.master_accession.split(',') if x.strip()]

	def load_master_segment_map(self):
		if not os.path.isfile(self.master_accession):
			return {}
		try:
			ref_df = load_reference_file_table(self.master_accession)
		except Exception:
			return {}
		if ref_df.empty or 'accession_type' not in ref_df.columns:
			return {}
		masters = ref_df[
			ref_df['accession_type'].fillna('').astype(str).str.strip().str.lower() == 'master'
		].copy()
		if masters.empty or 'segment' not in masters.columns:
			return {}
		masters['segment'] = masters['segment'].map(self._normalize_segment_value)
		masters['primary_accession'] = masters['primary_accession'].fillna('').astype(str).str.strip()
		masters = masters[(masters['segment'] != '') & (masters['primary_accession'] != '')]
		masters = masters.drop_duplicates(subset=['segment'], keep='first')
		return dict(zip(masters['segment'], masters['primary_accession']))

	def resolve_master_for_alignment(self, fasta_file, masters, master_segment_map=None):
		for master in masters:
			if fasta_file.startswith(master):
				return master

		segment = self._infer_segment_from_alignment_name(fasta_file)
		if segment and master_segment_map:
			master = master_segment_map.get(segment)
			if master:
				return master

		if len(masters) == 1:
			return masters[0]

		return None

	def get_gff_for_master(self, master):
		if isinstance(self.master_gff, list):
			for gff in self.master_gff:
				if master in os.path.basename(gff):
					return gff
		elif isinstance(self.master_gff, str):
			if master in os.path.basename(self.master_gff):
				return self.master_gff
		return None

	def get_gap_ranges(self, sequence):
		"""Inclusive 1-based ALIGNMENT COLUMN ranges (space 2) of each gap run."""
		gap_ranges = []
		start = None

		for i, char in enumerate(sequence):
			if char == '-':
				if start is None:
					start = i + 1
			else:
				if start is not None:
					gap_ranges.append([start, i])
					start = None

		if start is not None:
			gap_ranges.append([start, len(sequence)])

		return gap_ranges

	def count_gaps_before_position(self, gap_ranges, position):
		"""Gap characters at or before ``position`` in one row.

		``position`` is an ALIGNMENT COLUMN (space 2). It cannot be a master
		genome coordinate: ``gap_ranges`` are column ranges, so comparing a
		space-1 number against them is not an off-by-one, it is a comparison
		between two different number lines that happen to agree only while the
		master's row is ungapped.

		Subtracting the result from the same alignment column converts it to
		that row's own ungapped numbering (space 2 -> space 3), which is the
		only thing this is used for.

		A position landing *inside* a gap counts that gap up to and including
		itself, so the answer is the last real residue at or before it - what a
		feature boundary falling in a deletion should clamp to.
		"""
		count = 0
		for start, end in gap_ranges:
			if end < position:
				count += (end - start + 1)
			elif start <= position <= end:
				count += (position - start + 1)
		return count

	def count_gaps_strictly_before(self, gap_ranges, position):
		"""Gap characters strictly before ``position``, which is an ALIGNMENT COLUMN.

		The companion to :meth:`count_gaps_before_position`, and the difference
		between them is the difference between clamping a feature boundary
		backwards and clamping it forwards.

		``position - count_gaps_strictly_before(...)`` is the row's own index for
		the first residue at or AFTER the column - what a feature *start* wants.
		``position - count_gaps_before_position(...)`` is the index of the last
		residue at or BEFORE it - what a feature *end* wants.

		Using the inclusive count for a start was a defect in its own right: a CDS
		beginning inside one of this row's gaps resolved to the last residue
		before the feature, so the slice began early and the frame shifted; and a
		row whose leading padding covered the start column produced ``0``, which
		is not a valid 1-based coordinate at all. This form cannot go below 1,
		because it counts only whole gap columns that precede the position.
		"""
		count = 0
		for start, end in gap_ranges:
			if end < position:
				count += (end - start + 1)
			elif start < position <= end:
				count += (position - start)
		return count

	def format_genome_coverage(self, query_alignment, master_coord_to_aln_pos, feature_start, feature_end):
		try:
			feature_start = int(feature_start)
			feature_end = int(feature_end)
		except (TypeError, ValueError):
			return "NA"

		if feature_start > feature_end:
			feature_start, feature_end = feature_end, feature_start

		aln_indexes = []
		for master_pos in range(feature_start, feature_end + 1):
			if master_pos in master_coord_to_aln_pos:
				# Convert 1-based alignment position to 0-based string index
				aln_indexes.append(master_coord_to_aln_pos[master_pos] - 1)

		if not aln_indexes:
			return "NA"

		covered = 0
		for aln_index in aln_indexes:
			if aln_index >= len(query_alignment):
				continue

			base = query_alignment[aln_index]
			if base != '-' and base.upper() != 'N':
				covered += 1

		coverage = (covered / len(aln_indexes)) * 100
		return f"{coverage:.2f}"

	def recalculate_cds_coordinates_with_span(self, sequence_id, gap_ranges, cds_list, start_offset, genome_cord_start=None, genome_cord_end=None, master_coord_to_aln_pos=None):
		"""One row's view of every master CDS, in all three coordinate spaces.

		The order of operations is the whole point, so it is spelled out:

		  a. ``cds['start']`` / ``cds['end']`` arrive as master genome
		     coordinates (space 1).
		  b. They are clamped to ``genome_cord_start``/``genome_cord_end``, the
		     master-coordinate span this row covers - also space 1, so the
		     comparison is legitimate. A feature falling entirely outside the
		     span is dropped rather than recorded as a coordinate the row never
		     reached.
		  c. The clamped coordinates are converted to alignment columns
		     (space 1 -> space 2) through ``master_coord_to_aln_pos``. **This
		     has to happen before any gap arithmetic.**
		  d. Only then is the row's own gap count subtracted, converting the
		     alignment column to the row's ungapped numbering (space 2 ->
		     space 3).

		Doing (d) before (c) - subtracting a gap count measured in columns from
		a master genome coordinate - is the defect this ordering exists to
		prevent. It shifted every ``og`` coordinate upstream by the number of
		gaps in the master's own row, which is zero for an unsegmented build
		with no guide alignment (RABV, HCV) and non-zero for influenza segments
		2, 4, 6, 7 and 8. Those gap counts are multiples of three, so the CDS
		window stayed in frame and the damage never announced itself: the
		recorded start simply sat one or more codons early, translation began on
		a UTR codon, and the real stop codon fell outside the recorded end.

		``master_coord_to_aln_pos`` is optional. Without it, spaces 1 and 2 are
		assumed to coincide - correct only for an ungapped master, and the
		behaviour every caller that omits it already relies on.
		"""
		adjusted_coords = []
		clamp_to_span = genome_cord_start not in (None, "NA") and genome_cord_end not in (None, "NA")
		span_start = None
		span_end = None
		if clamp_to_span:
			span_start = int(str(genome_cord_start))
			span_end = int(str(genome_cord_end))

		for cds in cds_list:
			cds_start = int(cds['start'])
			cds_end = int(cds['end'])
			if clamp_to_span:
				assert span_start is not None and span_end is not None
				overlap_start = max(cds_start, span_start)
				overlap_end = min(cds_end, span_end)
				if overlap_start > overlap_end:
					continue
			else:
				overlap_start = cds_start
				overlap_end = cds_end

			# Convert to alignment positions BEFORE counting gaps.
			# count_gaps_before_position() measures the query's gaps in alignment
			# columns, so it can only be handed an alignment position.
			# overlap_start/overlap_end are master *genome* coordinates, and the two
			# part company the moment the master's own row carries gaps - which is
			# exactly what a guide alignment does to every segment with an insertion
			# (influenza segments 2, 4, 6, 7 and 8; segments 3 and 5 have none).
			# Subtracting the query's gap count from a master coordinate shifted every
			# OG coordinate upstream by the master's gap count, so the recorded CDS
			# start landed before the real ATG and the true stop codon fell outside
			# the recorded end. Those counts are multiples of three, so the window
			# stayed in frame and the damage showed up only as a translation that
			# began on a UTR codon - an internal stop whenever that codon was one.
			if master_coord_to_aln_pos is not None:
				aln_start = master_coord_to_aln_pos.get(overlap_start, overlap_start)
				aln_end = master_coord_to_aln_pos.get(overlap_end, overlap_end)
			else:
				aln_start = overlap_start
				aln_end = overlap_end

			# Start clamps FORWARDS (first residue at or after the column), end
			# clamps BACKWARDS (last residue at or before it). Using the inclusive
			# count for both put a start that landed inside a gap one or more
			# residues before the feature, and could return 0.
			adj_start = aln_start - self.count_gaps_strictly_before(gap_ranges, aln_start)
			adj_end = aln_end - self.count_gaps_before_position(gap_ranges, aln_end)

			adjusted_entry = {
				'start': aln_start,
				'end': aln_end,
				'og_start': adj_start,
				'og_end': adj_end,
				'feature_start': cds_start,
				'feature_end': cds_end,
				'product': cds['product'],
			}
			if adjusted_entry not in adjusted_coords:
				adjusted_coords.append(adjusted_entry)

		return adjusted_coords

	def load_blast_hits(self):
		if not self.blast_uniq_hits or not os.path.isfile(self.blast_uniq_hits):
			raise FileNotFoundError(f"BLAST unique hits file not found: {self.blast_uniq_hits}")
		acc_dict = {}
		for i in open(self.blast_uniq_hits):
			parts = i.strip().split('\t')
			if len(parts) != 4:
				raise ValueError(f"Malformed BLAST hits row in {self.blast_uniq_hits}: {i.strip()}")
			query, ref, score, strand = parts
			acc_dict[query] = ref
		return acc_dict

	def load_historical_alignment_lengths(self):
		if not self.update_db or not os.path.isfile(self.update_db):
			return {}
		conn = sqlite3.connect(self.update_db)
		try:
			row = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='sequence_alignment'").fetchone()
			if row is None:
				return {}
			
			cols = [r[1] for r in conn.execute("PRAGMA table_info(sequence_alignment)").fetchall()]
			aln_col = None
			for candidate in ["alignment", "aligned_seq", "sequence", "aln", "alignment_seq"]:
				if candidate in cols:
					aln_col = candidate
					break
			
			if not aln_col:
				return {}
				
			id_col = "primary_accession" if "primary_accession" in cols else "sequence_id" if "sequence_id" in cols else None
			if not id_col:
				return {}

			df = pd.read_sql_query(f"SELECT {id_col} as acc, LENGTH({aln_col}) as aln_len FROM sequence_alignment", conn)
			return dict(zip(df["acc"].astype(str), df["aln_len"]))
		except Exception as e:
			print(f"[warn] Could not load historic alignment lengths: {e}", file=sys.stderr)
			return {}
		finally:
			conn.close()

	def find_gaps_in_fasta(self):
		os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)
		update_scope_accessions = self.load_update_scope_accessions()
		existing_features = self.load_existing_feature_accessions()
		segment_map = self.load_segment_map()
		historical_lengths = self.load_historical_alignment_lengths()

		fasta_file_dir = self.paded_alignment
		if not fasta_file_dir or not os.path.isdir(fasta_file_dir):
			raise FileNotFoundError(f"Padded alignment directory not found: {fasta_file_dir}")
		fasta_files = [f for f in os.listdir(fasta_file_dir) if os.path.isfile(join(fasta_file_dir, f))]
		if not fasta_files:
			raise ValueError(f"No alignment files found in directory: {fasta_file_dir}")

		blast_dict = self.load_blast_hits()
		masters = self.get_master_list()
		master_segment_map = self.load_master_segment_map()
		if not masters:
			raise ValueError("No master accession could be resolved from --master_accession")

		header = [
			"accession",
			"master_ref_accession",
			"reference_accession",
			"aln_start",
			"aln_end",
			"cds_start",
			"cds_end",
			"cds_start_OG_seq",
			"cds_end_OG_seq",
			"product",
			"genome_coverage"
		]
		if segment_map:
			header.append("segment")

		with open(join(self.tmp_dir, self.output_dir, self.output_file), "w") as out_f:
			out_f.write("\t".join(header) + "\n")

			for fasta_file in fasta_files:
				current_master = self.resolve_master_for_alignment(fasta_file, masters, master_segment_map)
				if not current_master:
					print(f"Could not determine master for {fasta_file}. Skipping.")
					continue

				gff_file = self.get_gff_for_master(current_master)
				if not gff_file:
					print(f"No GFF found for master {current_master}. Skipping.")
					continue

				gff_dict = GffDictionary(gff_file).gff_dict
				cds_list = gff_dict['CDS']

				calc = CalculateGenomeCoordinates(join(fasta_file_dir, fasta_file), current_master)
				genome_coords = calc.extract_alignment_coordinates()
				
				fasta_records = list(SeqIO.parse(join(fasta_file_dir, fasta_file), "fasta"))
				master_record = next((r for r in fasta_records if r.id == current_master), None)
				if not master_record:
					master_record = fasta_records[0]
				master_alignment = str(master_record.seq)

				# The GFF is in record coordinates; this row may start partway
				# into the record. Everything derived from the row is therefore
				# shifted into record space by the same offset, so that the CDS
				# coordinates, the coordinate map and the covered span are all
				# in one space and can legitimately be compared.
				trim_offset = self.compute_trim_offset(current_master, master_alignment, gff_dict)

				master_coord_to_aln_pos = {}
				master_res_count = 0
				for align_pos, base in enumerate(master_alignment, start=1):
					if base != "-":
						master_res_count += 1
						master_coord_to_aln_pos[master_res_count + trim_offset] = align_pos

				for record in fasta_records:
					record_id = str(record.id).strip()
					sequence = str(record.seq)
					current_len = len(sequence)
					
					old_len = historical_lengths.get(record_id)
					backbone_expanded = old_len is not None and current_len > old_len

					if update_scope_accessions and record_id not in update_scope_accessions:
						if not backbone_expanded:
							continue
						else:
							print(f"[info] Forcing coordinate recalculation for historic ID {record_id} due to backbone expansion ({old_len} -> {current_len})")

					if update_scope_accessions and record_id in existing_features and not backbone_expanded:
						continue

					if record.id in genome_coords:
						master_acc, genome_cord_start, genome_cord_end = genome_coords[record.id]
						# CalculateGenomeCoordinates counts the master's residues
						# down the alignment row, so its span is in trimmed-row
						# numbering. Lift it into record space, because it is
						# about to be compared against GFF coordinates, which are
						# in record space, and the clamp has to compare like with
						# like. This shifted span stays internal.
						if trim_offset:
							if genome_cord_start not in (None, "NA"):
								genome_cord_start = int(genome_cord_start) + trim_offset
							if genome_cord_end not in (None, "NA"):
								genome_cord_end = int(genome_cord_end) + trim_offset
					else:
						genome_cord_start, genome_cord_end = "NA", "NA"

					# What gets REPORTED as aln_start/aln_end is the row's covered
					# span in alignment columns - the same space as cds_start and
					# cds_end. The two are compared directly by
					# ValidateDbTree.validate_feature_projection_integrity, which
					# clamps the master's cds columns by this span, so reporting
					# one in record coordinates and the other in columns makes
					# that check compare two different number lines. They agree
					# only for a master that is both flush and ungapped, which is
					# why this held for RABV and HCV and broke on influenza.
					covered_columns = alignment_covered_columns(sequence)
					if covered_columns is None:
						aln_span_start, aln_span_end = "NA", "NA"
					else:
						aln_span_start, aln_span_end = covered_columns

					gaps = self.get_gap_ranges(sequence)

					if gaps and gaps[0][0] == 1:
						start_offset = gaps[0][1] + 1
					else:
						start_offset = 1

					adjusted = self.recalculate_cds_coordinates_with_span(
						record.id,
						gaps,
						cds_list,
						start_offset,
						genome_cord_start=genome_cord_start,
						genome_cord_end=genome_cord_end,
						master_coord_to_aln_pos=master_coord_to_aln_pos,
					)

					for each_cords in adjusted:
						reference_acc = blast_dict[record.id] if record.id in blast_dict else current_master
						if segment_map:
							record_segment = segment_map.get(record.id, "")
							master_segment = segment_map.get(current_master, record_segment)
							ref_segment = segment_map.get(reference_acc, record_segment)
							if record_segment and master_segment and record_segment != master_segment:
								raise ValueError(f"Segment mismatch for {record.id}: record={record_segment}, master={master_segment}")
							if record_segment and ref_segment and record_segment != ref_segment:
								raise ValueError(f"Segment mismatch for {record.id}: record={record_segment}, ref={ref_segment}")

						# Execute coverage evaluation using tracking maps
						genome_coverage = self.format_genome_coverage(
							sequence,
							master_coord_to_aln_pos,
							each_cords.get('feature_start'),
							each_cords.get('feature_end')
						)

						data = [
							record.id,
							current_master,
							reference_acc,
							str(aln_span_start),
							str(aln_span_end),
							str(each_cords['start']),
							str(each_cords['end']),
							str(each_cords['og_start']),
							str(each_cords['og_end']),
							each_cords['product'],
							str(genome_coverage)
						]
						if segment_map:
							data.append(segment_map.get(record.id, ""))
						out_f.write('\t'.join(data) + "\n")

if __name__ == "__main__":
	parser = ArgumentParser(description='Calculates the genome and cds coordinates for a given sequences')
	parser.add_argument('-i', '--paded_alignment', help='Sequence file directory, it can be single or multiple fasta sequence files.', required=True)
	parser.add_argument('-b', '--tmp_dir', help='Base directory', default="tmp")
	parser.add_argument('-d', '--output_dir', help='Output directory where processed data and results are stored', default='Tables')
	parser.add_argument('-o', '--output_file', help='Output file name', default='features.tsv')
	parser.add_argument('-m', '--master_accession', help='Master accession', required=True)
	parser.add_argument('-bh', '--blast_uniq_hits', help='Blast unique hits file', default='tmp/Blast/query_uniq_tophits.tsv')
	parser.add_argument('-g', '--master_gff', help='Master GFF3 file(s)', required=True, nargs='+')
	parser.add_argument('--update_db', help='Existing DB path; when set, only emit feature rows for accessions not already in DB features table', default=None)
	parser.add_argument('--update_scope_tsv', help='TSV with primary_accession column; when set, only recalculate coordinates for these accessions', default=None)
	parser.add_argument('--segment_map_tsv', help='TSV with primary_accession and segment columns for segment-consistency checks', default=None)
	parser.add_argument('--master_seq_dir', help="Directory of master reference FASTA records (BLAST_ALIGNMENT's master_seq). Lets the 5' trim between a master's GFF, which is in record coordinates, and its row in the merged MSA be located. Without it a trim is reported but cannot be corrected.", default=None)
	args = parser.parse_args()

	processor = CalculateAlignmentCoordinates(args.paded_alignment, args.master_gff, args.tmp_dir, args.output_dir, args.output_file, args.master_accession, args.blast_uniq_hits, args.update_db, args.update_scope_tsv, args.segment_map_tsv, args.master_seq_dir)
	try:
		processor.find_gaps_in_fasta()
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		sys.exit(2)