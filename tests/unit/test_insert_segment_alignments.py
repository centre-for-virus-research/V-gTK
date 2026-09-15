"""InsertSegmentAlignments.py - the targeted per-segment write to sequence_alignment.

The behaviours pinned here are the ones that fail *silently* if they regress:
an upsert that appends instead of replacing, a reference losing its
self-referential alignment_name, or a record in the wrong coordinate frame being
stored next to correct ones.
"""

import sqlite3

import pytest

import InsertSegmentAlignments as isa


WIDTH = 12


def _make_db(path, rows=()):
	conn = sqlite3.connect(str(path))
	conn.execute('CREATE TABLE IF NOT EXISTS "sequence_alignment" ('
				 '"primary_accession" TEXT, "alignment_name" TEXT, "alignment" TEXT, '
				 '"segment" TEXT, "sequence_id" TEXT)')
	conn.executemany('INSERT INTO sequence_alignment VALUES (?,?,?,?,?)', rows)
	conn.commit()
	conn.close()


def _make_fasta(path, records):
	path.write_text(''.join('>%s\n%s\n' % (name, seq) for name, seq in records))


def _make_tophit(path, rows):
	# query, reference, percent identity, strand, segment - no header
	path.write_text(''.join('%s\t%s\t99.0\tplus\t%s\n' % row for row in rows))


@pytest.fixture
def build(tmp_path):
	"""A database, a two-query-plus-one-reference MSA, and its tophit table."""
	db = tmp_path / 'test.db'
	aln = tmp_path / 'aln.fasta'
	tophit = tmp_path / 'tophit.tsv'
	_make_db(db)
	_make_fasta(aln, [('Q00001', 'ACGT' * 3), ('Q00002', 'ACG-' * 3), ('REF001', 'AAAA' * 3)])
	# REF001 is deliberately absent: a reference was never a query.
	_make_tophit(tophit, [('Q00001', 'REF001', '4'), ('Q00002', 'REF001', '4'),
						  ('Q09999', 'REF777', '6')])
	return db, aln, tophit


def _run(db, aln, tophit, *extra):
	return isa.main(['--db', str(db), '--padded_aln', str(aln), '--tophit', str(tophit),
					 '--segment', '4', '--assert_width', str(WIDTH)] + list(extra))


def _rows(db, segment='4'):
	conn = sqlite3.connect(str(db))
	try:
		return conn.execute('SELECT primary_accession, alignment_name, alignment, segment, '
							'sequence_id FROM sequence_alignment WHERE segment = ? '
							'ORDER BY primary_accession', (segment,)).fetchall()
	finally:
		conn.close()


def test_inserts_one_row_per_record(build):
	db, aln, tophit = build
	_run(db, aln, tophit)
	rows = _rows(db)
	assert [r[0] for r in rows] == ['Q00001', 'Q00002', 'REF001']
	assert all(len(r[2]) == WIDTH for r in rows)
	assert all(r[3] == '4' for r in rows)


def test_alignment_name_is_the_reference_and_sequence_id_mirrors_accession(build):
	db, aln, tophit = build
	_run(db, aln, tophit)
	by_accession = {r[0]: r for r in _rows(db)}
	assert by_accession['Q00001'][1] == 'REF001'
	# A reference has no tophit row and is aligned against itself - the live DB
	# holds exactly this for the segment-4 master (AB573800|AB573800|4|1782).
	assert by_accession['REF001'][1] == 'REF001'
	assert all(row[4] == row[0] for row in _rows(db))


def test_tophit_lookup_is_segment_scoped(build):
	"""A query listed only under segment 6 must not pick up a segment-4 reference."""
	mapping = isa.load_tophit(str(build[2]), '4')
	assert 'Q09999' not in mapping
	assert isa.load_tophit(str(build[2]), '6') == {'Q09999': 'REF777'}


def test_upsert_replaces_rather_than_appending(build):
	"""The regression that silently doubles the table.

	INSERT OR REPLACE needs a UNIQUE index to conflict on; without one SQLite
	degrades it to a plain INSERT. Running twice must not grow the table.
	"""
	db, aln, tophit = build
	_run(db, aln, tophit)
	first = _rows(db)
	_run(db, aln, tophit)
	assert _rows(db) == first


def test_upsert_overwrites_a_stale_alignment(build):
	db, aln, tophit = build
	_make_db(db, [('Q00001', 'REF001', 'TTTT' * 3, '4', 'Q00001')])
	_run(db, aln, tophit)
	stored = {r[0]: r[2] for r in _rows(db)}
	assert stored['Q00001'] == 'ACGT' * 3


def test_replace_segment_clears_only_its_own_segment(build):
	db, aln, tophit = build
	_make_db(db, [('OLD001', 'REFOLD', 'GGGG' * 3, '4', 'OLD001'),
				  ('KEEP01', 'REFSIX', 'CCCC' * 3, '6', 'KEEP01')])
	_run(db, aln, tophit, '--mode', 'replace-segment')
	assert [r[0] for r in _rows(db)] == ['Q00001', 'Q00002', 'REF001']
	assert [r[0] for r in _rows(db, segment='6')] == ['KEEP01']


def test_record_of_the_wrong_width_is_refused(build):
	"""A differently-framed record must never be stored beside correct ones."""
	db, aln, tophit = build
	_make_fasta(aln, [('Q00001', 'ACGT' * 3), ('Q00002', 'ACGTACGT')])
	_run(db, aln, tophit)
	assert [r[0] for r in _rows(db)] == ['Q00001']


def test_dry_run_writes_nothing(build):
	db, aln, tophit = build
	_run(db, aln, tophit, '--dry_run')
	assert _rows(db) == []


def test_dry_run_cannot_modify_the_database(build, tmp_path):
	"""Opened read-only, so pointing a dry run at the live DB is safe.

	Setting journal_mode alone would write to the file, which is why this is
	asserted rather than assumed.
	"""
	db, aln, tophit = build
	before = (tmp_path / 'test.db').stat().st_mtime_ns
	_run(db, aln, tophit, '--dry_run')
	assert (tmp_path / 'test.db').stat().st_mtime_ns == before


def test_refuses_a_segment_with_no_tophit_rows(build):
	db, aln, tophit = build
	with pytest.raises(SystemExit):
		isa.main(['--db', str(db), '--padded_aln', str(aln), '--tophit', str(tophit),
				  '--segment', '7'])
