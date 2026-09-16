#!/usr/bin/env python3
"""Read-only audit of stored alignments, projected features and mutation evidence.

Uses residue-column lists and bisect, independently of CalcAlignmentCord and
AnnotateMutations. OG coordinates index the *degapped alignment*, not the raw
record. Raw records may contain insertions omitted from the alignment.

Coding diagnostics apply to complete, contiguous, forward-strand features with
the standard genetic code. They are warnings: a nucleotide projection need not
preserve the query's reading frame. This is not a test of biological homology,
nor an implementation of spliced/minus-strand/circular feature translation.

Reading-frame checks map each alignment residue back to the raw record. A
reference codon should come from three adjacent raw bases (not split by a
stripped insertion) in the feature's majority raw frame. Bases with more than
one possible raw position (e.g. a homopolymer beside an insertion) are left
unassessed, never flagged. The frame comes from the feature itself, not the
record's own CDS annotation, so a feature off frame along most of its length is
not detected, and an evenly split feature is reported as undetermined.

A call the annotator read from the record rather than the columns
(``codon_read = 'record_corrected'``) cannot be checked against the columns.
It is re-read from the raw record with reading_frame.ProteinFrame and must
match. That confirms the stored call is what the reader says; whether the
reader is right is what tests/unit/test_reading_frame.py's chaos tests
establish against simulated ground truth.
"""

import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import sqlite3

from Bio.Seq import Seq

import reading_frame


def positive_integer(value):
    """Reject fractional coordinates instead of silently truncating them."""
    text = str(value).strip()
    if not re.fullmatch(r"[0-9]+(?:\.0+)?", text):
        raise ValueError(f"Invalid coordinate: {value!r}")
    number = int(text.split('.')[0])
    if number < 1:
        raise ValueError(f"Non-positive coordinate: {value!r}")
    return number


def unique_embedding(short, long):
    """Index in long of each base of short, or None where several are possible.

    Returns None outright when short is not a subsequence of long (removing
    bases, e.g. nextalign insertions, cannot yield it). The leftmost and
    rightmost greedy embeddings bound every valid placement of a base; only
    where they agree is its raw position certain.
    """
    early, offset = [], 0
    for base in short:
        offset = long.find(base, offset)
        if offset == -1:
            return None
        early.append(offset)
        offset += 1
    late, offset = [], len(long)
    for base in reversed(short):
        offset = long.rfind(base, 0, offset)
        late.append(offset)
    late.reverse()
    return [left if left == right else None for left, right in zip(early, late)]


def runs(numbers):
    """[3, 4, 5, 9] -> [[3, 5], [9, 9]]"""
    grouped = []
    for number in numbers:
        if grouped and number == grouped[-1][1] + 1:
            grouped[-1][1] = number
        else:
            grouped.append([number, number])
    return grouped


def observed_residue(alignment, columns, covered):
    codon = ''.join(alignment[i] for i in columns)
    if codon == '---':
        # Terminal padding is missing coverage, not evidence of a deletion.
        if covered and covered[0] < columns[0] and columns[-1] < covered[-1]:
            return '-'
        return 'X'
    if len(codon) != 3 or set(codon) - set('ACGT'):
        return 'X'
    return str(Seq(codon).translate())


def translate_triplet(codon):
    if len(codon) != 3 or set(codon) - set('ACGT'):
        return 'X'
    return str(Seq(codon).translate())


def product_matches(product, protein):
    """Exact name or a complete token (NS3 in 'protease/helicase protein NS3').

    Callers require a unique matching feature; ambiguous names are reported.
    No virus-specific coordinate/feature lookup from the annotator is reused.
    """
    return product.casefold() == protein.casefold() or bool(
        re.search(r'(?<!\w)' + re.escape(protein) + r'(?!\w)', product, re.I)
    )


def audit_connection(conn):
    counts = Counter()
    findings = {}
    warnings = {}
    skipped = []

    def flag(code, example, warning=False):
        target = warnings if warning else findings
        item = target.setdefault(code, {'count': 0, 'examples': []})
        item['count'] += 1
        if len(item['examples']) < 10:
            item['examples'].append(example)

    def finish():
        return {'ok': not findings, 'counts': dict(counts), 'errors': findings,
                'warnings': warnings, 'skipped': skipped}

    def rows(table):
        cursor = conn.execute(f'SELECT * FROM "{table}"')
        columns = [column[0] for column in cursor.description]
        return [dict(zip(columns, row)) for row in cursor]

    def key(row, accession='primary_accession'):
        segment = str(row.get('segment') or '').strip()
        if re.fullmatch(r'[0-9]+\.0+', segment):
            segment = segment.split('.')[0]
        return str(row.get(accession) or '').strip(), segment

    frame_anchors = {}

    def raw_codon_start(identity, columns):
        """First raw position of a reference codon in this row, or why it has none.

        'gap': a column is gapped here (deletion or coverage, not a frame
        question); 'ambiguous': a base has more than one possible raw source;
        'split': the three bases are not adjacent in the raw record.
        """
        alignment = alignments[identity]
        if columns[-1] >= len(alignment) or any(alignment[i] == '-' for i in columns):
            return 'gap'
        residues = residue_columns[identity]
        positions = [embeddings[identity][bisect_left(residues, i)] for i in columns]
        if None in positions:
            return 'ambiguous'
        if positions != list(range(positions[0], positions[0] + 3)):
            return 'split'
        return positions[0]

    def check_feature_frame(identity, product, example, master_row, master_columns):
        """Walk the master's codon grid and find codons off the feature's majority raw frame.

        The majority, not the first readable codon, anchors the frame: a local
        shift at a feature boundary would otherwise flag every correct codon
        after it (measured on HCV mature peptides). With a tie, the codon runs
        are reported against the earlier-seen frame and calls are undetermined.
        """
        try:
            first, last = (positive_integer(master_row[k])
                           for k in ('cds_start_OG_seq', 'cds_end_OG_seq'))
        except ValueError:
            return
        grid = master_columns[first - 1:last]
        split, phases = [], []
        for number, offset in enumerate(range(0, len(grid) - 2, 3), 1):
            start = raw_codon_start(identity, grid[offset:offset + 3])
            if start == 'split':
                split.append(number)
            elif isinstance(start, str):
                counts[f'frame_codons_{start}'] += 1
            else:
                phases.append((number, start % 3))
        top = Counter(phase for _, phase in phases).most_common(2)
        frame = top[0][0] if top else None
        tie = len(top) == 2 and top[0][1] == top[1][1]
        shifted = [number for number, phase in phases if phase != frame]
        frame_anchors[(identity, product)] = (frame, tie)
        counts['frame_features_checked'] += 1
        counts['frame_features_tied'] += tie
        counts['frame_codons_in_frame'] += len(phases) - len(shifted)
        counts['frame_codons_split'] += len(split)
        counts['frame_codons_out_of_frame'] += len(shifted)
        if split:
            flag('codon_split_by_raw_insertion',
                 {**example, 'count': len(split), 'codons': split[:20]}, warning=True)
        if shifted:
            flag('reading_frame_shift_in_raw',
                 {**example, 'count': len(shifted), 'codon_runs': runs(shifted)[:20],
                  'frame_tie': tie}, warning=True)

    tables ={r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    integrity = [r[0] for r in conn.execute('PRAGMA integrity_check')]
    if integrity != ['ok']:
        flag('sqlite_integrity', integrity)
    required = {
        'sequence_alignment': {'primary_accession', 'alignment', 'alignment_name'},
        'features': {'accession', 'master_ref_accession', 'product', 'aln_start',
                     'aln_end', 'cds_start', 'cds_end', 'cds_start_OG_seq', 'cds_end_OG_seq'},
    }
    for table, columns in required.items():
        present = {r[1] for r in conn.execute(f'PRAGMA table_info("{table}")')}
        if not columns <= present:
            flag('missing_schema', {'table': table, 'columns': sorted(columns - present)})
    if findings:
        return finish()

    alignments = {}
    residue_columns = {}
    widths = defaultdict(set)
    for row in rows('sequence_alignment'):
        counts['alignment_rows'] += 1
        identity = key(row)
        if identity in alignments:
            flag('duplicate_alignment', identity)
        text = str(row['alignment'] or '').upper()
        alignments[identity] = text
        residue_columns[identity] = [i for i, base in enumerate(text) if base != '-']
        widths[identity[1]].add(len(text))
        if not residue_columns[identity]:
            flag('empty_alignment', identity)
        if set(text) - set('ACGTRYSWKMBDHVNU-'):
            flag('invalid_alignment_alphabet', identity)
    if not alignments:
        flag('no_alignments', 'sequence_alignment is empty')
    for segment, lengths in widths.items():
        if len(lengths) > 1:
            flag('ragged_alignment', {'segment': segment, 'widths': sorted(lengths)})
    for row in rows('sequence_alignment'):
        if (str(row['alignment_name']), key(row)[1]) not in alignments:
            flag('missing_alignment_reference', key(row))

    # Raw record position of each alignment residue (None where ambiguous).
    embeddings = {}
    oriented = {}
    if 'sequences' in tables:
        raw = {key(r, 'header'): re.sub(r'\s+', '', str(r['sequence'] or '')).upper().replace('U', 'T')
               for r in rows('sequences')}
        for identity, alignment in alignments.items():
            sequence = raw.get(identity)
            if sequence is None:
                flag('missing_raw_sequence', identity)
                continue
            ungapped = alignment.replace('-', '').replace('U', 'T')
            embedding = unique_embedding(ungapped, sequence)
            if ungapped == sequence:
                counts['raw_exact_matches'] += 1
            elif embedding is not None:
                counts['raw_forward_subsequences'] += 1
            else:
                sequence = str(Seq(sequence).reverse_complement())
                embedding = unique_embedding(ungapped, sequence)
                if embedding is None:
                    flag('alignment_not_raw_subsequence', identity)
                    continue
                counts['raw_reverse_subsequences'] += 1
            embeddings[identity] = embedding
            oriented[identity] = sequence
    else:
        skipped.append('raw sequence comparison and reading-frame checks: sequences table absent')

    features = rows('features')
    counts['feature_rows'] = len(features)
    if not features:
        flag('no_features', 'features is empty')
    master_features = defaultdict(list)
    row_masters = defaultdict(set)
    for row in features:
        identity = key(row, 'accession')
        master = str(row['master_ref_accession'])
        row_masters[identity].add(master)
        if row['accession'] == master:
            master_features[(master, identity[1])].append(row)

    for row in features:
        identity = key(row, 'accession')
        example = {'accession': identity[0], 'segment': identity[1], 'product': row['product']}
        alignment = alignments.get(identity)
        columns = residue_columns.get(identity)
        if alignment is None or not columns:
            flag('feature_without_alignment', example)
            continue
        try:
            start, end, og_start, og_end, span_start, span_end = (
                positive_integer(row[k]) for k in ('cds_start', 'cds_end',
                    'cds_start_OG_seq', 'cds_end_OG_seq', 'aln_start', 'aln_end')
            )
        except ValueError as exc:
            flag('invalid_feature_coordinate', {**example, 'reason': str(exc)})
            continue
        if not 1 <= start <= end <= len(alignment):
            flag('feature_alignment_bounds', example)
            continue
        if (span_start, span_end) != (columns[0] + 1, columns[-1] + 1):
            flag('covered_span_mismatch', example)
        expected = (bisect_left(columns, start - 1) + 1, bisect_right(columns, end - 1))
        if not 1 <= og_start <= og_end <= len(columns):
            flag('feature_og_bounds', example)
        if (og_start, og_end) != expected:
            flag('feature_og_mismatch', {**example, 'stored': [og_start, og_end],
                                       'expected': list(expected)})
        fragment = alignment[start - 1:end].replace('-', '')
        if not fragment:
            flag('feature_without_residues', example)
        if fragment != alignment.replace('-', '')[og_start - 1:og_end]:
            flag('feature_slice_mismatch', example)
        master_key = (str(row['master_ref_accession']), identity[1])
        candidates = [r for r in master_features[master_key] if r['product'] == row['product']]
        master_columns = residue_columns.get(master_key, [])
        if not candidates or not master_columns:
            flag('missing_master_feature', example)
            continue
        if identity in embeddings and len(candidates) == 1:
            check_feature_frame(identity, row['product'], example, candidates[0], master_columns)
        # Coverage is measured where BOTH rows have residues; query-only
        # insertion columns do not define a master-record coordinate.
        shared = [i for i in master_columns if i < len(alignment) and alignment[i] != '-']
        expected_spans = []
        full_spans = []
        for candidate in candidates:
            try:
                left, right = (positive_integer(candidate[k]) for k in ('cds_start', 'cds_end'))
            except ValueError:
                continue
            full_spans.append((left, right))
            if shared:
                expected_spans.append((max(left, shared[0] + 1), min(right, shared[-1] + 1)))
        if (start, end) not in expected_spans:
            flag('master_projection_mismatch', example)
        if (start, end) in full_spans and fragment and not set(fragment) - set('ACGTRYSWKMBDHVNU'):
            counts['complete_feature_slices'] += 1
            if len(fragment) % 3:
                flag('coding_length_not_multiple_of_three', example, warning=True)
            protein = str(Seq(fragment[:len(fragment) // 3 * 3]).translate())
            if '*' in protein[:-1]:
                flag('internal_stop_in_degapped_feature', example, warning=True)
            reference_columns = [i for i in master_columns if start - 1 <= i < end]
            stop_positions = []
            # Read the reference's triplets without joining across query gaps.
            # This diagnoses problems a degapped-ORF check alone cannot separate.
            for offset in range(0, len(reference_columns) - 3, 3):
                triplet = reference_columns[offset:offset + 3]
                codon = ''.join(alignment[i] for i in triplet)
                if not set(codon) - set('ACGT'):
                    counts['concrete_internal_reference_codons'] += 1
                    if codon in {'TAA', 'TAG', 'TGA'}:
                        stop_positions.append(offset // 3 + 1)
                elif '-' in codon and codon != '---':
                    counts['partial_gap_internal_reference_codons'] += 1
                else:
                    counts['unknown_or_deleted_internal_reference_codons'] += 1
            if stop_positions:
                flag('internal_stop_in_reference_codons',
                     {**example, 'positions': stop_positions}, warning=True)

    if 'sequence_mutation_calls' not in tables:
        skipped.append('mutation evidence: sequence_mutation_calls table absent')
        return finish()

    call_columns = {r[1] for r in conn.execute('PRAGMA table_info(sequence_mutation_calls)')}
    if not {'primary_accession', 'protein_name', 'aa_position', 'alt_residue',
            'observed_residue'} <= call_columns:
        flag('mutation_schema', sorted(call_columns))
        return finish()
    for row in rows('sequence_mutation_calls'):
        counts['mutation_call_rows'] += 1
        identity = key(row)
        example = {'accession': identity[0], 'segment': identity[1],
                   'mutation_id': row.get('mutation_id')}
        masters = row_masters.get(identity, set())
        if identity not in alignments or len(masters) != 1:
            flag('unresolved_mutation_master', example)
            continue
        master_key = (next(iter(masters)), identity[1])
        matches = [r for r in master_features[master_key]
                   if product_matches(r['product'], row['protein_name'])]
        if len(matches) != 1 or master_key not in residue_columns:
            flag('unresolved_mutation_feature', example)
            continue
        try:
            position = positive_integer(row['aa_position'])
            start = positive_integer(matches[0]['cds_start_OG_seq'])
            end = positive_integer(matches[0]['cds_end_OG_seq'])
        except ValueError:
            flag('invalid_mutation_coordinate', example)
            continue
        offset = start - 1 + (position - 1) * 3
        columns = residue_columns[master_key][offset:offset + 3]
        if offset + 3 > end or len(columns) != 3 or columns[-1] >= len(alignments[identity]):
            flag('mutation_outside_feature', example)
            continue
        frame, tie = frame_anchors.get((identity, matches[0]['product']), (None, False))
        if row.get('codon_read') == 'record_corrected':
            # Read from the record by construction, so the columns cannot
            # confirm it; a fresh read of the raw record has to.
            reread = None
            if identity in oriented and master_key in alignments:
                master_columns = residue_columns[master_key]
                coord_map = {i + 1: column for i, column in enumerate(master_columns)}
                grid = reading_frame.codon_grid(coord_map, start, end)
                frame_reader = reading_frame.ProteinFrame(
                    reading_frame.RecordMapping(alignments[identity], oriented[identity]), grid,
                    translate_triplet, 'X', alignments[master_key])
                reread = frame_reader.read(tuple(columns))
            stored = row['observed_residue']
            if reread != (stored, reading_frame.CORRECTED) or stored != row['alt_residue']:
                flag('record_corrected_residue_unverified',
                     {**example, 'aa_position': position, 'stored': stored,
                      'reread': list(reread) if reread else None})
            else:
                counts['record_corrected_residues_verified'] += 1
            continue
        actual = observed_residue(alignments[identity], columns, residue_columns[identity])
        # Unknown/missing sequence cannot support any stored call evidence.
        if actual == 'X' or actual != row['observed_residue'] or actual != row['alt_residue']:
            flag('mutation_residue_mismatch', {**example, 'actual': actual,
                                             'stored': row['observed_residue']})
        else:
            counts['mutation_residues_verified'] += 1
        if identity in embeddings:
            start = raw_codon_start(identity, columns)
            if isinstance(start, str):
                status = start
            elif frame is None:
                status = 'unanchored'
            elif tie:
                status = 'frame_undetermined'
            else:
                status = 'out_of_frame' if start % 3 != frame else 'in_frame'
            counts[f'call_codons_{status}'] += 1
            if status in {'split', 'out_of_frame', 'frame_undetermined'}:
                # Without a call_status column every stored call counts as evidence.
                emitted = row.get('call_status') in (None, 'emitted')
                flag(f"{'emitted' if emitted else 'suppressed'}_call_not_one_raw_codon",
                     {**example, 'aa_position': position, 'call_status': row.get('call_status'),
                      'problem': status}, warning=True)
    return finish()


def audit_database(path):
    path = Path(path).resolve()
    conn = sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)
    try:
        conn.execute('PRAGMA query_only=ON')
        conn.execute('BEGIN')  # All checks see the same snapshot.
        return audit_connection(conn)
    finally:
        conn.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, required=True)
    parser.add_argument('--json', type=Path, help='Write the complete audit report')
    parser.add_argument('--strict-coding', action='store_true',
                        help='Also fail on coding warnings (use for synthetic intact ORFs)')
    args = parser.parse_args(argv)
    report = audit_database(args.db)
    if args.json:
        if args.json.resolve() == args.db.resolve():
            parser.error('--json must differ from --db')
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    return 0 if report['ok'] and not (args.strict_coding and report['warnings']) else 1


if __name__ == '__main__':
    raise SystemExit(main())
