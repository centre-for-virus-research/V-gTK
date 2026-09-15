#!/usr/bin/env python3
"""Insert one segment's padded alignment rows into `sequence_alignment`.

Nothing else in `scripts/` does a *targeted, per-segment* write to this table.
`CreateSqliteDB.py` builds it wholesale as part of a full database build, and
`PadAlignment --update_db` / `UsherPlacement --update_db` only ever read it. When
a rebuild produces a corrected alignment for one segment - as the treepatch
rebuild does for HA and NA - there is no way to put it back without rebuilding
the whole 29 GB database. This script is that way.

Written against `treepatch/FINALISE_DB.md` section 1.5, which specifies it.

What a row is
-------------
One record of the padded merged MSA becomes one row:

    primary_accession  the bare accession, from the FASTA header
    alignment_name     the REFERENCE this query was aligned against, looked up
                       in query_uniq_tophit_annotated.tsv. Not the refset name.
                       A reference or master is aligned against itself, so its
                       alignment_name is its own accession - which is exactly
                       what the live database holds for AB573800.
    alignment          the gapped sequence, asserted to --assert_width
    segment            --segment, as TEXT, matching meta_data.segment
    sequence_id        a verbatim copy of primary_accession, which is what
                       CreateSqliteDB._normalize_alignment_columns produces

Why the UNIQUE index is not optional
------------------------------------
`CreateSqliteDB.UPSERT_TABLES` includes `sequence_alignment` with natural key
(primary_accession, alignment_name, segment), and the established write is
INSERT OR REPLACE on that key. But INSERT OR REPLACE needs something to
conflict *on*: with no unique index it degrades to a plain INSERT and appends a
duplicate instead of replacing. That is a bug CreateSqliteDB has already been
bitten by and documents at UPSERT_INDEX_NAMES. The live influenza database
carries no indexes at all, so this script creates the same index
(`idx_seq_alignment_upsert`) under the same name before writing.

Two modes
---------
--mode upsert (default)
    INSERT OR REPLACE on the natural key. Faithful to CreateSqliteDB. Note the
    sharp edge: if a stored row for the same accession has a *different*
    alignment_name - because the tophit assignment changed between runs - the
    key differs, so the old row is not replaced and the accession ends up with
    two alignment rows in different coordinate frames. The script counts that
    afterwards and says so.

--mode replace-segment
    DELETE every row for the segment, then INSERT. Guarantees exactly one row
    per accession for that segment, so the final count is just the number of
    records in the MSA. Correct when the padded MSA is the new truth for the
    whole segment and is a superset of what is stored, which is the treepatch
    case. Slower by one table scan.

Usage
-----
    python scripts/InsertSegmentAlignments.py \
        --db "$RDB" \
        --padded_aln treepatch/work/pad/refset_4_aln_merged_MSA.fasta \
        --tophit .../inputs/query_uniq_tophit_annotated.tsv \
        --segment 4 --assert_width 1782 --mode replace-segment
"""

import argparse
import os
import sqlite3
import sys
import time


def log(message):
	print('[%s] %s' % (time.strftime('%H:%M:%S'), message), flush=True)


def bare_accession(name):
	"""'PV123012.1 some description' -> 'PV123012'.

	The same normalisation CreateSqliteDB.IDENTITY_COLUMNS enforces, so the
	accession joins against meta_data.primary_accession.
	"""
	return str(name).strip().split()[0].split('.')[0] if str(name).strip() else ''


def load_tophit(path, segment):
	"""{query accession: reference accession} for one segment.

	Mirrors `load_closest_reference()` in treepatch/bin/tp_common.py, which is
	the existing parser for this file: no header, tab separated, and the columns
	are query, reference, percent identity, strand, segment.
	"""
	segment = str(segment).strip()
	mapping = {}
	with open(path, 'r', encoding='utf-8', errors='replace') as handle:
		for line in handle:
			parts = line.rstrip('\n').rstrip('\r').split('\t')
			if len(parts) >= 5 and parts[4].strip() == segment:
				mapping[bare_accession(parts[0])] = parts[1].strip()
	return mapping


def iter_fasta(path):
	"""Stream (header, sequence) pairs.

	Hand rolled rather than Biopython because these files are ~1 GB and this
	only needs the header line and the joined sequence.
	"""
	name, chunks = None, []
	with open(path, 'r', encoding='utf-8', errors='replace') as handle:
		for line in handle:
			if line.startswith('>'):
				if name is not None:
					yield name, ''.join(chunks)
				name, chunks = line[1:].strip(), []
			elif name is not None:
				chunks.append(line.strip())
	if name is not None:
		yield name, ''.join(chunks)


def ensure_upsert_index(conn):
	"""Create the UNIQUE index INSERT OR REPLACE needs, or explain why it cannot.

	Same name and same key columns as CreateSqliteDB.UPSERT_INDEX_NAMES, so a
	database that already carries it is left alone rather than given a second
	identical index.
	"""
	try:
		conn.execute(
			'CREATE UNIQUE INDEX IF NOT EXISTS idx_seq_alignment_upsert '
			'ON sequence_alignment (primary_accession, alignment_name, segment)')
		return True
	except sqlite3.IntegrityError:
		duplicates = conn.execute(
			'SELECT primary_accession, alignment_name, segment, COUNT(*) c '
			'FROM sequence_alignment GROUP BY 1,2,3 HAVING c > 1 LIMIT 5').fetchall()
		rendered = '; '.join('%s/%s/%s x%d' % row for row in duplicates) or '(none listed)'
		raise SystemExit(
			'Cannot create the UNIQUE index on sequence_alignment'
			'(primary_accession, alignment_name, segment): the database already holds\n'
			'rows sharing a key, e.g. %s\n'
			'Without the index INSERT OR REPLACE silently appends instead of replacing,\n'
			'so fix the duplicates or run with --mode replace-segment.' % rendered)


def main(argv=None):
	parser = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--db', required=True, help='sqlite database to write to')
	parser.add_argument('--padded_aln', required=True,
						help='padded merged MSA for this segment, e.g. refset_4_aln_merged_MSA.fasta')
	parser.add_argument('--tophit', required=True,
						help='query_uniq_tophit_annotated.tsv, source of alignment_name')
	parser.add_argument('--segment', required=True, help='segment this alignment belongs to')
	parser.add_argument('--assert_width', type=int,
						help='required alignment width; a record of any other width is refused')
	parser.add_argument('--mode', choices=('upsert', 'replace-segment'), default='upsert',
						help='upsert on the natural key (default), or clear the segment first')
	parser.add_argument('--batch_size', type=int, default=20000,
						help='rows per transaction (default 20000)')
	parser.add_argument('--dry_run', action='store_true',
						help='read and validate the alignment, write nothing')
	args = parser.parse_args(argv)

	for path in (args.db, args.padded_aln, args.tophit):
		if not os.path.exists(path):
			raise SystemExit('no such file: %s' % path)

	segment = str(args.segment).strip()

	log('loading tophit references for segment %s' % segment)
	tophit = load_tophit(args.tophit, segment)
	log('  %d query->reference assignments' % len(tophit))
	if not tophit:
		raise SystemExit('no tophit rows for segment %s - wrong file or wrong segment?' % segment)

	# A dry run must not be able to modify the database it is pointed at, and
	# setting journal_mode alone would: it writes to the file. Opening read-only
	# makes "validate this against the live DB" a safe thing to ask for.
	if args.dry_run:
		conn = sqlite3.connect('file:%s?mode=ro' % os.path.abspath(args.db), uri=True)
	else:
		conn = sqlite3.connect(args.db)
		conn.execute('PRAGMA journal_mode=WAL')
		conn.execute('PRAGMA synchronous=NORMAL')

	before = conn.execute('SELECT COUNT(*) FROM sequence_alignment').fetchone()[0]
	before_segment = conn.execute(
		"SELECT COUNT(*) FROM sequence_alignment WHERE TRIM(COALESCE(segment,'')) = ?",
		(segment,)).fetchone()[0]
	log('sequence_alignment holds %d rows, %d for segment %s'
		% (before, before_segment, segment))

	if args.dry_run:
		log('dry run: validating %s' % args.padded_aln)
	elif args.mode == 'replace-segment':
		log('replace-segment: deleting the %d existing segment-%s rows' % (before_segment, segment))
		conn.execute("DELETE FROM sequence_alignment WHERE TRIM(COALESCE(segment,'')) = ?", (segment,))
		conn.commit()
		ensure_upsert_index(conn)
	else:
		ensure_upsert_index(conn)

	statement = ('INSERT OR REPLACE INTO sequence_alignment '
				 '(primary_accession, alignment_name, alignment, segment, sequence_id) '
				 'VALUES (?, ?, ?, ?, ?)')

	read = written = self_named = bad_width = blank = 0
	widths = {}
	batch = []
	log('streaming %s' % args.padded_aln)
	for header, sequence in iter_fasta(args.padded_aln):
		read += 1
		accession = bare_accession(header)
		if not accession or not sequence:
			blank += 1
			continue
		widths[len(sequence)] = widths.get(len(sequence), 0) + 1
		if args.assert_width and len(sequence) != args.assert_width:
			bad_width += 1
			continue
		# A reference or master has no tophit row because it was never a query;
		# it is aligned against itself, which is how the live DB stores AB573800.
		reference = tophit.get(accession)
		if not reference:
			reference = accession
			self_named += 1
		batch.append((accession, reference, sequence, segment, accession))
		if len(batch) >= args.batch_size:
			if not args.dry_run:
				conn.executemany(statement, batch)
				conn.commit()
			written += len(batch)
			batch = []
			if written % (args.batch_size * 10) == 0:
				log('  %d records read, %d written' % (read, written))
	if batch:
		if not args.dry_run:
			conn.executemany(statement, batch)
			conn.commit()
		written += len(batch)

	log('read %d records, wrote %d' % (read, written))
	log('  %d had no tophit row and were named against themselves (references/masters)' % self_named)
	if blank:
		log('  %d records skipped: empty header or sequence' % blank)
	if bad_width:
		log('  %d records REFUSED for a width other than %d' % (bad_width, args.assert_width))
	log('  widths seen: %s' % sorted(widths.items(), key=lambda kv: -kv[1])[:5])

	if args.assert_width and bad_width:
		log('WARNING: refused records mean this MSA is not a single coordinate frame.')

	if args.dry_run:
		log('dry run: nothing written')
		conn.close()
		return 0

	after = conn.execute('SELECT COUNT(*) FROM sequence_alignment').fetchone()[0]
	after_segment = conn.execute(
		"SELECT COUNT(*) FROM sequence_alignment WHERE TRIM(COALESCE(segment,'')) = ?",
		(segment,)).fetchone()[0]
	log('sequence_alignment now holds %d rows (%+d), %d for segment %s (%+d)'
		% (after, after - before, after_segment, segment, after_segment - before_segment))

	# An accession with two alignment rows for one segment is the failure mode
	# described in the module docstring: alignment_residues() would read whichever
	# row SQLite returned first, and the two can be in different frames.
	doubled = conn.execute(
		"SELECT COUNT(*) FROM (SELECT primary_accession FROM sequence_alignment "
		"WHERE TRIM(COALESCE(segment,'')) = ? GROUP BY 1 HAVING COUNT(*) > 1)",
		(segment,)).fetchone()[0]
	if doubled:
		log('WARNING: %d accessions now have more than one segment-%s alignment row.' % (doubled, segment))
		log('         Re-run with --mode replace-segment to collapse them.')
	else:
		log('  one alignment row per accession for segment %s' % segment)

	conn.close()
	return 0


if __name__ == '__main__':
	sys.exit(main())
