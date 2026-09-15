#!/usr/bin/env python3
"""Reading a cohort and its protein alleles out of a V-gTK database.

The growth-rate estimators need two things from the database and nothing else:
a table of sequences with a numeric sampling time, and, for each of them,
which variant of a protein residue they carry. Getting those two things right
is most of the work, so it is kept here, away from the statistics.

Three routes to the alleles are supported, in decreasing order of how much the
database has to already know:

``features``
    The protein is an annotated product in `features`, so every sequence has
    its own CDS start and the residue is read straight out of the alignment.
    Exact, and the route that works for a virus with one protein per feature -
    RABV N or G, an influenza segment.

``catalog``
    The protein only appears in `sequence_mutations`, i.e. it was annotated
    against a mutation catalogue. Carriers are the recorded calls. This is
    exact for the sequences that were called and says *nothing* about the rest,
    so the denominator is restricted to sequences that were evaluated for that
    protein. Reporting a frequency against the whole database instead would
    divide a real numerator by an unrelated denominator.

``alignment``
    The protein is a sub-peptide of a polyprotein that the database does not
    annotate separately (HCV NS5A inside the polyprotein, say). The protein's
    start in the master's coordinates has to be supplied, or calibrated against
    the catalogue calls - and because that calibration is an inference, its
    agreement is always reported and it is never selected automatically.
"""

import sqlite3

import numpy as np
import pandas as pd

from AnnotateMutations import (DELETION_RESIDUE, UNKNOWN_RESIDUE,
							   build_alignment_coordinate_map, translate_codon)

#: meta_data columns that may hold a genotype, best first. The provenance-aware
#: columns come first because a build that has them also has genotype_origin
#: beside them; a bare `genotype` column may be a vendor declaration.
GENOTYPE_COLUMNS = ('nearest_reference_genotype', 'genotype', 'gisaid_genotype',
					'ncbi_genotype', 'genbank_genotype')
SUBTYPE_COLUMNS = ('nearest_reference_subtype', 'subtype', 'gisaid_subtype',
				   'ncbi_subtype', 'genbank_subtype')

#: Symbols the vectorised translator distinguishes. Everything that is not a
#: base or a gap collapses onto the last slot and translates to X, which is the
#: same collapse AnnotateMutations.translate_codon makes.
_SYMBOLS = 'ACGT-N'
_SYMBOL_COUNT = len(_SYMBOLS)


def _build_symbol_table():
	lookup = np.full(256, _SYMBOL_COUNT - 1, dtype=np.uint8)
	for index, symbol in enumerate(_SYMBOLS):
		lookup[ord(symbol)] = index
		lookup[ord(symbol.lower())] = index
	lookup[ord('U')] = _SYMBOLS.index('T')
	lookup[ord('u')] = _SYMBOLS.index('T')
	return lookup


def _build_codon_table():
	"""Residue for every symbol triple, filled from translate_codon itself.

	Deriving the table rather than writing one out keeps this translator and the
	annotator's in step by construction - a change to how the annotator reads an
	ambiguous or gapped codon cannot silently disagree with the frequencies
	reported here.
	"""
	table = np.empty(_SYMBOL_COUNT ** 3, dtype='<U1')
	for i, a in enumerate(_SYMBOLS):
		for j, b in enumerate(_SYMBOLS):
			for k, c in enumerate(_SYMBOLS):
				codon = a + b + c
				if codon == '---':
					residue = DELETION_RESIDUE
				else:
					residue = translate_codon(codon) or UNKNOWN_RESIDUE
				table[i * _SYMBOL_COUNT ** 2 + j * _SYMBOL_COUNT + k] = residue
	return table


_SYMBOL_LOOKUP = _build_symbol_table()
_CODON_TABLE = _build_codon_table()


_GAP_SYMBOL = _SYMBOLS.index('-')


def covered_spans(block):
	"""First and last non-gap column of every row in a stacked alignment.

	Everything outside that span is padding added to square the alignment up -
	sequence the submitter never reported, not sequence the virus is missing.
	"""
	non_gap = block != ord('-')
	has_any = non_gap.any(axis=1)
	width = block.shape[1]
	first = np.argmax(non_gap, axis=1)
	last = width - 1 - np.argmax(non_gap[:, ::-1], axis=1)
	# A row with no bases at all gets an empty span, so nothing is ever inside it.
	return np.where(has_any, first, width), np.where(has_any, last, -1)


def translate_columns(block, offsets, covered=None):
	"""Translate one codon out of a stack of aligned sequences.

	`block` is a ``(sequences x columns)`` uint8 array of ASCII, `offsets` the
	three column indices of the codon. Vectorised because the whole point of
	this route is genotyping every sequence in the database at a site, and doing
	that a codon at a time in Python is minutes per position.

	`covered` is the ``(first, last)`` pair from :func:`covered_spans`. With it,
	an all-gap codon counts as a deletion only when it sits strictly *inside*
	the sequence's covered span; outside, it is padding and translates to X.
	This is the same distinction AnnotateMutations draws, and it is not a
	nicety: without it every sequence that stops short of the CDS end acquires
	"deletion" alleles across the tail. Sequencing became more complete over
	time, so those deletions correlate with date and come out as the fastest
	growing alleles in the run - an artefact of coverage, not of the virus.
	"""
	a, b, c = (_SYMBOL_LOOKUP[block[:, offset]].astype(np.int32) for offset in offsets)
	residues = _CODON_TABLE[a * _SYMBOL_COUNT ** 2 + b * _SYMBOL_COUNT + c]
	if covered is None:
		return residues
	first, last = covered
	all_gap = (a == _GAP_SYMBOL) & (b == _GAP_SYMBOL) & (c == _GAP_SYMBOL)
	inside = (min(offsets) > first) & (max(offsets) < last)
	return np.where(all_gap & ~inside, UNKNOWN_RESIDUE, residues)


# ---------------------------------------------------------------------------
# Database introspection
# ---------------------------------------------------------------------------

def table_columns(conn, table):
	try:
		return [row[1] for row in conn.execute('PRAGMA table_info("%s")' % table)]
	except sqlite3.DatabaseError:
		return []


def table_exists(conn, table):
	row = conn.execute("SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?",
					   (table,)).fetchone()
	return row is not None


def resolve_label_column(conn, requested, candidates):
	"""Pick the meta_data column holding a genotype (or subtype) label."""
	columns = table_columns(conn, 'meta_data')
	if requested:
		if requested not in columns:
			raise ValueError('meta_data has no column %r (available: %s)'
							 % (requested, ', '.join(columns)))
		return requested
	lowered = {column.lower(): column for column in columns}
	for candidate in candidates:
		if candidate in lowered:
			return lowered[candidate]
	return None


# ---------------------------------------------------------------------------
# Cohort
# ---------------------------------------------------------------------------

#: Days per month used to place a month-only date in the middle of its month.
_MONTH_STARTS = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)


def _as_int(value):
	try:
		text = str(value).strip()
		if not text:
			return None
		return int(float(text))
	except (TypeError, ValueError):
		return None


def decimal_date(year, month=None, day=None):
	"""Decimal year for a possibly partial collection date.

	A year-only date becomes mid-year and a month-only date mid-month, never
	1 January. Defaulting to the start of the period puts every imprecise date
	systematically earlier than the precise ones around it, and since imprecise
	dates are commoner in older records that bias is *correlated with time* -
	exactly the axis every estimator here regresses on.
	"""
	year = _as_int(year)
	if year is None:
		return float('nan'), 'none'
	month = _as_int(month)
	day = _as_int(day)
	leap = (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
	length = 366.0 if leap else 365.0
	if month is None or not 1 <= month <= 12:
		return year + 0.5, 'year'
	start = _MONTH_STARTS[month - 1] + (1 if (leap and month > 2) else 0)
	if day is None or not 1 <= day <= 31:
		next_start = (_MONTH_STARTS[month] + (1 if (leap and month >= 2) else 0)
					  if month < 12 else length)
		return year + ((start + next_start) / 2.0) / length, 'month'
	return year + (start + day - 0.5) / length, 'day'


PRECISION_RANK = {'none': 0, 'year': 1, 'month': 2, 'day': 3}


def load_cohort(conn, genotype_column=None, subtype_column=None, require_date=True,
				drop_excluded=True, min_precision='year', accession_types=None,
				min_year=None, max_year=None, segment=None):
	"""One row per sequence, with a decimal sampling date and its labels."""
	columns = table_columns(conn, 'meta_data')
	if not columns:
		raise ValueError('this database has no meta_data table')
	genotype_column = resolve_label_column(conn, genotype_column, GENOTYPE_COLUMNS)
	subtype_column = resolve_label_column(conn, subtype_column, SUBTYPE_COLUMNS)

	wanted = ['primary_accession', 'collection_year', 'collection_mon', 'collection_day',
			  'accession_type', 'country_validated', 'country', 'host_scientific_name',
			  'host', 'segment', 'exclusion_status', 'strain']
	select = [column for column in dict.fromkeys(wanted) if column in columns]
	for extra in (genotype_column, subtype_column):
		if extra and extra not in select:
			select.append(extra)
	frame = pd.read_sql_query('SELECT %s FROM meta_data' % ', '.join('"%s"' % c for c in select), conn)

	frame['primary_accession'] = frame['primary_accession'].astype(str).str.strip()
	frame = frame[frame['primary_accession'] != '']
	frame = frame.drop_duplicates(subset=['primary_accession'], keep='first')

	dates = [decimal_date(row.get('collection_year'), row.get('collection_mon'),
						  row.get('collection_day'))
			 for row in frame.to_dict('records')]
	frame['decimal_date'] = [value for value, _ in dates]
	frame['date_precision'] = [precision for _, precision in dates]

	frame['genotype'] = (frame[genotype_column].astype(str).str.strip()
						 if genotype_column else '')
	frame['subtype'] = (frame[subtype_column].astype(str).str.strip()
						if subtype_column else '')
	for column in ('genotype', 'subtype'):
		frame[column] = frame[column].replace({'nan': '', 'None': '', 'NA': '', 'na': ''})

	if drop_excluded and table_exists(conn, 'excluded_accessions'):
		excluded = {str(row[0]).strip() for row
					in conn.execute('SELECT primary_accession FROM excluded_accessions')}
		if excluded:
			frame = frame[~frame['primary_accession'].isin(excluded)]
	if drop_excluded and 'exclusion_status' in frame.columns:
		status = frame['exclusion_status'].astype(str).str.strip().str.lower()
		frame = frame[~status.isin({'excluded', 'exclude', 'true', 'yes', '1'})]

	if accession_types:
		frame = frame[frame.get('accession_type', '').astype(str).str.strip().isin(set(accession_types))]
	if segment is not None and 'segment' in frame.columns:
		frame = frame[frame['segment'].astype(str).str.strip() == str(segment).strip()]
	if require_date:
		frame = frame[np.isfinite(frame['decimal_date'])]
		floor = PRECISION_RANK.get(min_precision, 1)
		frame = frame[frame['date_precision'].map(PRECISION_RANK).fillna(0) >= floor]
	if min_year is not None:
		frame = frame[frame['decimal_date'] >= float(min_year)]
	if max_year is not None:
		frame = frame[frame['decimal_date'] < float(max_year) + 1.0]

	frame = frame.reset_index(drop=True)
	frame.attrs['genotype_column'] = genotype_column
	frame.attrs['subtype_column'] = subtype_column
	return frame


def genotype_labels_from_reference_list(conn, reference_tsv, tree_newick):
	"""Genotype/subtype per accession, inherited from tree neighbours.

	A database built before genotype columns existed - or one whose genotypes
	were never resolved - has the information anyway: the curated reference list
	names the genotype of each reference, and the tree says which references
	each query sits with. This reuses the pipeline's own clade assignment so a
	genotype derived here is the same genotype the pipeline would have stored.
	"""
	from ExportRefListFromUpdateDb import load_reference_file_table
	from clade_from_tree import assign_labels_from_tree

	references = load_reference_file_table(reference_tsv)
	labels = {}
	for _, row in references.iterrows():
		accession = str(row.get('primary_accession', '')).strip()
		if not accession:
			continue
		labels[accession] = {'genotype': str(row.get('genotype', '')).strip(),
							 'subtype': str(row.get('subtype', '')).strip()}
	if not labels:
		raise ValueError('reference list %r carried no genotypes' % reference_tsv)
	assigned = assign_labels_from_tree(tree_newick, labels)
	merged = dict(labels)
	for accession, entry in assigned.items():
		merged.setdefault(accession, {'genotype': '', 'subtype': ''})
		merged[accession] = {'genotype': str(entry.get('genotype', '')).strip(),
							 'subtype': str(entry.get('subtype', '')).strip()}
	return merged


# ---------------------------------------------------------------------------
# Alleles
# ---------------------------------------------------------------------------

def list_proteins(conn):
	"""Proteins this database can be asked about, and where they come from."""
	rows = []
	if table_exists(conn, 'sequence_mutations'):
		for name, sites, calls in conn.execute(
				'SELECT protein_name, COUNT(DISTINCT aa_position), '
				'COUNT(DISTINCT primary_accession) FROM sequence_mutations '
				'WHERE protein_name IS NOT NULL AND TRIM(protein_name) != "" '
				'GROUP BY protein_name'):
			rows.append({'protein': str(name), 'source': 'catalog',
						 'n_sites': int(sites), 'n_sequences': int(calls)})
	if table_exists(conn, 'features'):
		for product, count in conn.execute(
				'SELECT product, COUNT(*) FROM features '
				'WHERE product IS NOT NULL AND TRIM(product) != "" GROUP BY product'):
			rows.append({'protein': str(product), 'source': 'features',
						 'n_sites': 0, 'n_sequences': int(count)})
	return pd.DataFrame(rows, columns=['protein', 'source', 'n_sites', 'n_sequences'])


def catalog_alleles(conn, protein):
	"""Recorded residue calls for one protein, de-duplicated.

	`sequence_mutations` carries one row per (sequence, mutation, combination),
	so a sequence appearing in sixty drug combinations has sixty identical rows
	for the same residue. Counting those as sixty carriers would weight a
	sequence by how many combinations its catalogue happens to list.
	"""
	frame = pd.read_sql_query(
		'SELECT DISTINCT primary_accession, mutation_id, protein_name, aa_position, '
		'alt_residue FROM sequence_mutations WHERE protein_name = ?', conn, params=(protein,))
	if frame.empty:
		return frame
	frame['primary_accession'] = frame['primary_accession'].astype(str).str.strip()
	frame['alt_residue'] = frame['alt_residue'].astype(str).str.strip()
	frame['aa_position'] = pd.to_numeric(frame['aa_position'], errors='coerce')
	return frame.dropna(subset=['aa_position'])


def master_accession(conn, segment=None):
	"""The master reference, for one segment if the build has several.

	`segment` is not optional in practice for a segmented virus. An influenza
	database carries one master per segment - eight of them - and the first row
	SQLite happens to return is whichever segment it stored first. Resolving HA
	coordinates against segment 2's master produces residues for every sequence
	and every one of them is wrong, with nothing in the output to say so.
	"""
	if segment is not None:
		row = conn.execute(
			"SELECT primary_accession FROM meta_data "
			"WHERE TRIM(LOWER(COALESCE(accession_type,''))) = 'master' "
			"AND TRIM(COALESCE(segment,'')) = ? LIMIT 1",
			(str(segment).strip(),)).fetchone()
		if row:
			return str(row[0]).strip()
		return None
	row = conn.execute("SELECT primary_accession FROM meta_data "
					   "WHERE TRIM(LOWER(COALESCE(accession_type,''))) = 'master' "
					   "LIMIT 1").fetchone()
	return str(row[0]).strip() if row else None


def count_masters(conn):
	"""How many masters the build has, so a caller can insist on a segment."""
	row = conn.execute("SELECT COUNT(*) FROM meta_data "
					   "WHERE TRIM(LOWER(COALESCE(accession_type,''))) = 'master'").fetchone()
	return int(row[0]) if row else 0


def _fetch_alignment(conn, accession, segment=None):
	if segment is not None:
		row = conn.execute(
			"SELECT alignment FROM sequence_alignment WHERE primary_accession = ? "
			"AND TRIM(COALESCE(segment,'')) = ? LIMIT 1",
			(accession, str(segment).strip())).fetchone()
		if row and row[0]:
			return str(row[0])
	row = conn.execute('SELECT alignment FROM sequence_alignment WHERE primary_accession = ? '
					   'LIMIT 1', (accession,)).fetchone()
	return str(row[0]) if row and row[0] else None


def feature_cds_start(conn, protein, accession):
	"""CDS start of an annotated product for one accession, if it has one."""
	row = conn.execute('SELECT cds_start FROM features WHERE accession = ? AND '
					   'LOWER(TRIM(product)) = LOWER(TRIM(?)) LIMIT 1',
					   (accession, protein)).fetchone()
	if not row or row[0] in (None, ''):
		return None
	try:
		return int(float(row[0]))
	except (TypeError, ValueError):
		return None


def calibrate_protein_start(conn, protein, master=None, max_anchors=600, seed=0):
	"""Infer where a protein starts, by testing candidate frames against the calls.

	For a protein the database annotates only through a mutation catalogue -
	a sub-peptide of a polyprotein - nothing records where it begins in the
	master's coordinates, but the recorded calls do constrain it: at the true
	start every call should agree with the residue the alignment carries.

	The agreement at the winning start is returned with it and is the whole
	point. A clean frame agrees with nearly every call; a protein whose
	numbering is genotype-relative, or whose region carries indels, does not,
	and the caller is expected to refuse a weak calibration rather than quietly
	genotype a database against a guess. Nothing selects this route
	automatically.
	"""
	master = master or master_accession(conn)
	if not master:
		raise ValueError('no master accession in meta_data; supply one explicitly')
	master_alignment = _fetch_alignment(conn, master)
	if not master_alignment:
		raise ValueError('master %s has no stored alignment' % master)

	calls = catalog_alleles(conn, protein)
	calls = calls[calls['alt_residue'].str.len() == 1]
	if calls.empty:
		raise ValueError('no catalogue calls for protein %r to calibrate against' % protein)
	if len(calls) > max_anchors:
		calls = calls.sample(n=max_anchors, random_state=seed)

	accessions = sorted(set(calls['primary_accession']))
	placeholders = ','.join('?' * len(accessions))
	stored = dict(conn.execute(
		'SELECT primary_accession, alignment FROM sequence_alignment '
		'WHERE primary_accession IN (%s)' % placeholders, accessions))
	usable = [accession for accession in accessions if stored.get(accession)]
	if not usable:
		raise ValueError('none of the calling sequences for %r has an alignment' % protein)

	width = len(master_alignment)
	block = np.frombuffer(''.join(str(stored[a]) for a in usable).encode('ascii'),
						  dtype=np.uint8).reshape(len(usable), width)
	row_of = {accession: index for index, accession in enumerate(usable)}
	calls = calls[calls['primary_accession'].isin(row_of)]
	rows = np.array([row_of[a] for a in calls['primary_accession']])
	positions = calls['aa_position'].to_numpy(dtype=int)
	residues = calls['alt_residue'].to_numpy(dtype='<U1')

	coord_map = build_alignment_coordinate_map(master_alignment)
	highest = int(positions.max())
	best = {'protein': protein, 'master': master, 'start': None,
			'agreement': float('nan'), 'n_anchors': int(len(calls)), 'n_scored': 0,
			'runner_up_agreement': float('nan')}
	scores = []
	for start in range(1, len(coord_map) - highest * 3 + 1):
		columns = np.array([[coord_map.get(start + (p - 1) * 3 + k, -1) for k in range(3)]
							for p in np.unique(positions)])
		if np.any(columns < 0):
			continue
		column_of = {position: columns[index] for index, position
					 in enumerate(np.unique(positions))}
		hit = scored = 0
		for position in np.unique(positions):
			mask = positions == position
			offsets = column_of[position]
			codons = block[rows[mask]][:, offsets]
			a, b, c = (_SYMBOL_LOOKUP[codons[:, i]].astype(np.int32) for i in range(3))
			called = _CODON_TABLE[a * _SYMBOL_COUNT ** 2 + b * _SYMBOL_COUNT + c]
			informative = (called != UNKNOWN_RESIDUE)
			scored += int(np.sum(informative))
			hit += int(np.sum(called[informative] == residues[mask][informative]))
		if scored:
			scores.append((hit / scored, scored, start))
	if not scores:
		return best
	scores.sort(reverse=True)
	best.update({'start': int(scores[0][2]), 'agreement': float(scores[0][0]),
				 'n_scored': int(scores[0][1]),
				 'runner_up_agreement': float(scores[1][0]) if len(scores) > 1 else float('nan')})
	return best


#: How `features.cds_start` / `cds_end` are to be read against a stored
#: alignment. The two are the same thing only when the master carries no gaps,
#: which is why this went unnoticed on HCV.
COORD_SPACES = ('auto', 'column', 'ungapped')


def codon_columns(master_alignment, start, positions, coord_space='column'):
	"""Alignment columns of each requested codon, under one coordinate reading.

	``column``   - `start` is an alignment column, so codon n sits at columns
				   ``start + 3(n-1)`` onwards. This is what the influenza builds
				   store: the HA master's `cds_end` is 1782, which is exactly the
				   alignment width, while its ungapped length is 1650.
	``ungapped`` - `start` counts bases in the master's own sequence, ignoring
				   gaps, and each is then mapped to the column holding it.

	On a gapless master the two agree exactly. On the influenza HA master, which
	carries 132 gaps, reading the ungapped way puts the frame 18 codons out and
	runs 44 codons off the end of the sequence.
	"""
	if coord_space not in ('column', 'ungapped'):
		raise ValueError('coord_space must be column or ungapped')
	if coord_space == 'ungapped':
		coord_map = build_alignment_coordinate_map(master_alignment)

		def resolve(nucleotide):
			return coord_map.get(nucleotide, -1)
	else:
		width = len(master_alignment)

		def resolve(nucleotide):
			return nucleotide - 1 if 1 <= nucleotide <= width else -1

	offsets = {}
	for position in sorted({int(p) for p in positions}):
		columns = [resolve(int(start) + (position - 1) * 3 + k) for k in range(3)]
		if min(columns) >= 0:
			offsets[position] = columns
	return offsets


def _as_block(text):
	# 'replace' rather than strict: a stray byte in one stored alignment should
	# translate to X, not abort the whole run.
	return np.frombuffer(text.encode('ascii', 'replace'), dtype=np.uint8).reshape(1, -1)


def detect_coord_space(master_alignment, start, end):
	"""Work out which reading of the CDS coordinates actually yields a protein.

	Translates the master's own CDS both ways and prefers the reading with no
	internal stop codons, then the one that leaves fewest codons unresolved. A
	coding sequence read in the right frame has no internal stops; one read in
	the wrong frame has them every few dozen codons, so the test is decisive and
	needs nothing but the master itself.

	Returns ``(coord_space, scores)``.
	"""
	n_codons = max(0, (int(end) - int(start) + 1) // 3)
	scores = {}
	if n_codons < 10:
		return 'column', scores
	block = _as_block(master_alignment)
	for space in ('column', 'ungapped'):
		offsets = codon_columns(master_alignment, start, range(1, n_codons + 1), space)
		text = ''.join(str(translate_columns(block, columns)[0])
					   for columns in offsets.values())
		scores[space] = {
			'resolved': len(offsets),
			'missing': n_codons - len(offsets),
			'internal_stops': (text[:-1].count('*') if text else n_codons),
			'unknown': text.count(UNKNOWN_RESIDUE),
			'protein_head': text[:40],
		}

	def rank(space):
		entry = scores[space]
		return (entry['internal_stops'], entry['missing'] + entry['unknown'])

	return min(('column', 'ungapped'), key=rank), scores


def resolve_coord_space(conn, master, start, end, segment=None, requested='auto'):
	"""``(coord_space, scores)`` for one master's CDS, honouring an explicit choice.

	Kept beside the detection so a caller never has to reach for the master's
	alignment itself just to ask the question.
	"""
	if requested in ('column', 'ungapped'):
		return requested, {}
	alignment = _fetch_alignment(conn, master, segment=segment)
	if not alignment:
		return 'column', {}
	return detect_coord_space(alignment, start, end)


def alignment_residues(conn, positions, start, accessions=None, master=None,
					   segment=None, coord_space='column', batch_size=20000):
	"""Residue at each requested protein position, for every sequence.

	Returns ``(frame, info)``: a DataFrame with one row per accession and one
	column per position, and a dict recording which master and coordinate
	reading were used and how many rows were skipped.

	`segment` is pushed into the SQL rather than filtered in Python. On a
	segmented build that is not a tidiness point: the influenza database stores
	all eight segments in one table, so without it every HA query drags the
	other seven segments' alignments off disk to throw them away - eight times
	the I/O for the same answer, against a 30 GB file.

	Rows whose alignment is a different width from the master's are counted and
	skipped. They are in a different coordinate space, and reading fixed columns
	out of them would invent residues rather than fail.
	"""
	master = master or master_accession(conn, segment=segment)
	master_alignment = _fetch_alignment(conn, master, segment=segment) if master else None
	if not master_alignment:
		raise ValueError('no master alignment to resolve protein coordinates against '
						 '(segment %r, master %r)' % (segment, master))
	offsets = codon_columns(master_alignment, start, positions, coord_space)

	width = len(master_alignment)
	# Columns where the master itself is a gap are insertion columns relative to
	# it: nearly every sequence is a gap there, and that is a property of the
	# alignment rather than evidence about the reading frame.
	master_gaps = sorted(position for position, columns in offsets.items()
						 if set(master_alignment[column] for column in columns) == {'-'})
	info = {'master': master, 'segment': segment, 'coord_space': coord_space,
			'master_width': width, 'positions_requested': len(set(int(p) for p in positions)),
			'positions_resolved': len(offsets), 'rows_read': 0, 'rows_used': 0,
			'rows_wrong_width': 0, 'master_gap_positions': master_gaps,
			'last_position': max(offsets) if offsets else None}
	if not offsets:
		return pd.DataFrame(columns=['primary_accession']), info

	query = 'SELECT primary_accession, alignment FROM sequence_alignment'
	params = ()
	if segment is not None:
		query += " WHERE TRIM(COALESCE(segment,'')) = ?"
		params = (str(segment).strip(),)
	wanted = None
	if accessions is not None:
		wanted = {str(a).strip() for a in accessions}

	names, chunks, frames = [], [], []

	def flush():
		if not names:
			return
		block = np.frombuffer(''.join(chunks).encode('ascii', 'replace'),
							  dtype=np.uint8).reshape(len(names), width)
		covered = covered_spans(block)
		data = {'primary_accession': list(names)}
		for position, columns in offsets.items():
			data[position] = translate_columns(block, columns, covered=covered)
		frames.append(pd.DataFrame(data))
		names.clear()
		chunks.clear()

	for accession, alignment in conn.execute(query, params):
		info['rows_read'] += 1
		accession = str(accession).strip()
		if wanted is not None and accession not in wanted:
			continue
		text = str(alignment or '')
		if len(text) != width:
			info['rows_wrong_width'] += 1
			continue
		names.append(accession)
		chunks.append(text)
		info['rows_used'] += 1
		if len(names) >= batch_size:
			flush()
	flush()
	if not frames:
		return pd.DataFrame(columns=['primary_accession']), info
	return pd.concat(frames, ignore_index=True), info
