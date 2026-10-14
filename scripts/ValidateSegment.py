import csv
import argparse
import os
import shutil
from collections import Counter

import segment_utils

'''Annotate gB_matrix with BLAST results. Includes overwrite input file (default) and overwrite exclusion options'''


#: A submitter-declared segment that BLAST contradicts is only overridden when
#: BLAST is clearly right. Below either floor the call is not trusted, so the row
#: is excluded rather than labelled one way and aligned another. The floors come
#: from the 301 influenza rows that disagreed in the 2026 DB: every one of the 234
#: kept rows had >= 83% identity, and the 3 that failed these floors were 164-340
#: bp fragments or 83% hits.
MIN_CONFLICT_PIDENT = 85.0
MIN_CONFLICT_LENGTH = 300


def _row_length(row):
    for key in ('real_length', 'length'):
        try:
            return int(float(row.get(key)))
        except (TypeError, ValueError):
            continue
    return None


def _finalise_segment(row, blast_call, by_segment, by_name, stats, valid_segments=None,
                      blast_identity=None, min_pident=MIN_CONFLICT_PIDENT,
                      min_length=MIN_CONFLICT_LENGTH, conflicts=None):
    """Write the canonical ``segment``, its display name, and its provenance.

    ``segment`` stays the single definitive key: normalised, and numeric wherever
    the virus keys on numbers. ``segment_name`` is a derived label only - never a
    key - so that the gene vocabulary (HA, NA, PB1) is available for display and
    for cross-checking a submitter's declaration without ever being mistaken for
    an identifier. Keeping them apart matters: SegmentPivotTable treats the
    literal 'na' as "no segment assigned", so a gene name used as an identifier
    would silently drop every neuraminidase row.

    Returns the exclusion reason when a submitter/BLAST conflict could not be
    resolved, else ''.
    """
    resolved, source, candidates = resolve_segment(row, blast_call, by_name, valid_segments)
    row.setdefault('segment_submitted', '')
    unresolved = ''

    # Compare the STREAMS against each other. The previous version re-derived the
    # prior value from the same field stream 1 had already used, so it compared a
    # value with itself: when stream 1 won it could never fire, and when a later
    # stream won it fired on a value that had merely been rejected as missing -
    # reporting a "disagreement" that was nothing of the kind.
    distinct = {value for value in candidates.values() if value}
    if len(distinct) > 1:
        stats['segment_disagreement'] += 1
        stats[f'conflict:{"|".join(f"{k}={v}" for k, v in sorted(candidates.items()))}'] += 1

    # The alignment is built from the BLAST call, never from the submitter's
    # label, so when the two disagree the column the database is keyed on has to
    # follow BLAST or the row sits in one segment's table and another's tree.
    # GenBank genome assemblies (OX/OY/OZ) are the usual source: their
    # "segment: N" is the contig order, not the gene, and BLAST is a clean 1-8
    # permutation across each assembly. The old rule let the declaration win and
    # counted the conflict without acting on it.
    declared, blasted = candidates.get('source_declared'), candidates.get('blast_inferred')
    if source == 'source_declared' and blasted and blasted != declared \
            and (valid_segments is None or blasted.casefold() in valid_segments):
        length = _row_length(row)
        confident = ((blast_identity is None or blast_identity >= min_pident)
                     and (length is None or length >= min_length))
        outcome = 'blast_overrode_declared' if confident else 'conflict_unresolved'
        stats[f'conflict_outcome:{outcome}'] += 1
        if conflicts is not None:
            conflicts.append((row.get('primary_accession', ''), declared, blasted,
                              '' if blast_identity is None else blast_identity,
                              '' if length is None else length, outcome))
        if confident:
            row['segment_submitted'] = declared
            resolved, source = blasted, 'blast_overrode_declared'
        else:
            source = 'conflict_unresolved'
            unresolved = (f"segment conflict: submitter declares {declared}, BLAST says {blasted} "
                          f"but not confidently enough to relabel "
                          f"(identity {blast_identity}, length {length})")

    row['segment'] = resolved
    row['segment_name'] = segment_utils.segment_to_name(resolved, by_segment)
    row['segment_source'] = source
    stats[f'source:{source or "unresolved"}'] += 1
    return unresolved

def build_segment_map(file_path):
    segment_map = {}
    with open(file_path, 'r') as ann_file:
        reader = csv.reader(ann_file, delimiter='\t')
        for row in reader:
            segment_map[row[0]] = row[4]
    return segment_map

def build_identity_map(file_path):
    """{accession: BLAST percent identity} from the annotated top-hit file."""
    identity_map = {}
    with open(file_path, 'r') as ann_file:
        for row in csv.reader(ann_file, delimiter='\t'):
            try:
                identity_map[row[0]] = float(row[2])
            except (IndexError, ValueError):
                continue
    return identity_map


def build_reference_map(file_path):
    reference_map = {}
    with open(file_path, 'r') as ann_file:
        reader = csv.reader(ann_file, delimiter='\t')
        for row in reader:
            reference_map[row[0]] = row[1]
    return reference_map


def _is_reference_like(row):
    acc_type = str(row.get('accession_type', '')).strip().lower()
    return acc_type in {'master', 'reference', 'excluded', 'exclusion_list'}

def load_valid_segments(ref_list_path):
    """The segment labels this build targets, from all non-exclusion_list rows.

    A superset of the master rows, so it can never exclude a segment the pipeline
    is actually building. For influenza it is exactly {1..8} - the B/C/D decoys in
    the segment column are all exclusion_list entries - so passing a flu reference
    list reproduces the hardcoded 1-8 range check exactly.
    """
    from ExportRefListFromUpdateDb import load_reference_file_table
    from SegmentPivotTable import normalise_segment_label

    refs = load_reference_file_table(ref_list_path)
    acc_type = refs['accession_type'].astype(str).str.strip().str.lower()
    labels = set()
    for value in refs[acc_type != 'exclusion_list']['segment']:
        label = normalise_segment_label(value)
        if label is not None:
            labels.add(label.casefold())
    return labels or None


def resolve_segment(row, blast_call, by_name, valid_segments=None):
    """Decide the canonical ``segment`` for one row, and say where it came from.

    Three independent streams can carry the segment, and before this function
    only the first ever reached the ``segment`` column - which is why 65% of a
    GenBank+GISAID influenza matrix came out of this script with ``segment``
    empty and had to be backfilled by hand:

    1. the source's own declaration (GenBank ``/segment=``), already in ``segment``
    2. this pipeline's BLAST top-hit inference, i.e. ``blast_call``
    3. an external submitter's declared segment *name* in ``segment_declared``
       (GISAID ships ``HA``/``NA``/``PB1``... where we key on ``1``..``8``)

    Precedence follows the established policy of the manual backfill: an existing
    value is never overwritten, gaps are filled from the inference, and anything
    still empty falls back to the external declaration. Names are resolved to
    canonical numbers via the segment-name asset, so a row declaring ``PB1``
    becomes segment ``2`` rather than - as the old digit-scraping normalisers did -
    segment ``1``.

    Returns ``(segment, source)``; ``source`` is ``''`` when nothing resolved.
    """
    def _resolve(value, translated_only):
        """Normalise, then translate a name to its canonical number if possible.

        ``translated_only`` refuses a label we could not translate. It guards the
        external stream: a vendor vocabulary must never leak into the key column
        just because we had no mapping for it. The pipeline's own streams are
        allowed to keep an untranslatable label, because that is how genuinely
        name-keyed viruses (Lassa's L/S) are represented.
        """
        normalised = segment_utils.normalise_segment(value)
        if normalised is None:
            return None
        if normalised.isdigit():
            return normalised
        mapped = segment_utils.name_to_segment(normalised, by_name)
        if mapped is not None:
            return mapped
        return None if translated_only else normalised

    # Every stream is evaluated before one is chosen. Returning on the first hit
    # would mean the later streams are never computed, and then a source-vs-BLAST
    # conflict - the thing most worth reporting - is invisible by construction.
    candidates, acceptable = {}, {}
    for value, source, translated_only in (
        (row.get('segment'), 'source_declared', False),
        (blast_call, 'blast_inferred', False),
        (row.get('segment_declared'), 'external_declared', True),
    ):
        if value in (None, '', 'not found'):
            continue
        resolved = _resolve(value, translated_only)
        if resolved is None or resolved.casefold() in segment_utils.FALLBACK_MISSING_TOKENS:
            continue
        candidates[source] = resolved
        # The reference list is authoritative about which labels this build has.
        # Without this check any free text in GenBank's /segment= qualifier -
        # 'X', 'RNA fragment', 'unknown' - is accepted as the canonical key and
        # beats a good BLAST call, because /segment= is submitter free text.
        if valid_segments is not None and resolved.casefold() not in valid_segments:
            continue
        acceptable[source] = resolved

    for source in ('source_declared', 'blast_inferred', 'external_declared'):
        if source in acceptable:
            return acceptable[source], source, candidates
    return '', '', candidates


def annotate_matrix(matrix_file, segment_map, reference_map, output_file, overwrite=False, overwrite_exclusions=False, valid_segments=None, segment_names_path=None, identity_map=None, min_pident=MIN_CONFLICT_PIDENT, min_length=MIN_CONFLICT_LENGTH, conflicts_out=None):
    tmp_output = output_file if not overwrite else output_file + ".tmp"
    by_segment, by_name = segment_utils.load_segment_names(segment_names_path)
    stats = Counter()
    conflicts = []
    identity_map = identity_map or {}

    with open(matrix_file, 'r') as matrix, open(tmp_output, 'w', newline='') as output:
        reader = csv.DictReader(matrix, delimiter='\t')

        fieldnames = reader.fieldnames.copy()
        # 'segment' is normally already present, but a matrix can reach here without
        # it (non-segmented builds, and the minimal fixtures in the test suite).
        # _finalise_segment always writes it, so it must always be declared.
        for added in ('segment', 'segment_validated', 'closest_reference', 'exclusion',
                      'segment_name', 'segment_source', 'segment_submitted'):
            if added not in fieldnames:
                fieldnames.append(added)

        writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            primary_accession = row['primary_accession']

            if _is_reference_like(row):
                row['segment_validated'] = row.get('segment_validated', '') or row.get('segment', '')
                row['closest_reference'] = row.get('closest_reference', '') or primary_accession
                _finalise_segment(row, segment_map.get(primary_accession), by_segment, by_name, stats, valid_segments)
                writer.writerow(row)
                continue

            row['segment_validated'] = segment_map.get(primary_accession, 'not found')
            unresolved = _finalise_segment(
                row, segment_map.get(primary_accession), by_segment, by_name, stats, valid_segments,
                blast_identity=identity_map.get(primary_accession), min_pident=min_pident,
                min_length=min_length, conflicts=conflicts)
            row['closest_reference'] = reference_map.get(primary_accession, 'not found')

            if overwrite_exclusions or not row.get('exclusion'):
                if row['segment_validated'] == 'not found':
                    row['exclusion'] = 'not significant BLAST hit'
                elif valid_segments is not None:
                    # Segmented virus of any kind: the reference list defines the
                    # valid labels, so L/M/S survive instead of being stamped
                    # 'non IAV genomic sequence' and dropped from the database.
                    from SegmentPivotTable import normalise_segment_label

                    label = normalise_segment_label(row['segment_validated'])
                    if label is None or label.casefold() not in valid_segments:
                        row['exclusion'] = 'segment not in reference segment set'
                    else:
                        row['exclusion'] = ''
                else:
                    try:
                        segment_value = int(row['segment_validated'])
                        if not 1 <= segment_value <= 8:
                            row['exclusion'] = 'non IAV genomic sequence'
                        else:
                            row['exclusion'] = ''
                    except ValueError:
                        row['exclusion'] = 'non IAV genomic sequence'

            # A conflict that could not be resolved is excluded, but never over
            # an exclusion the row already carries.
            if unresolved and not row.get('exclusion'):
                row['exclusion'] = unresolved

            writer.writerow(row)

    if overwrite:
        shutil.move(tmp_output, output_file)

    if conflicts_out:
        with open(conflicts_out, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['primary_accession', 'declared', 'blast', 'blast_pident', 'length', 'outcome'])
            writer.writerows(conflicts)

    print(f"Annotated file written to {output_file}")
    resolved_total = sum(v for k, v in stats.items() if k.startswith('source:') and k != 'source:unresolved')
    print(f"Segment resolution: {resolved_total} resolved, "
          f"{stats.get('source:unresolved', 0)} unresolved")
    for key in sorted(k for k in stats if k.startswith('source:')):
        print(f"  - {key.split(':', 1)[1]}: {stats[key]}")
    if stats.get('segment_disagreement'):
        print(f"  [warn] {stats['segment_disagreement']} row(s) where the submitter's segment "
              f"disagreed with BLAST: "
              f"{stats.get('conflict_outcome:blast_overrode_declared', 0)} relabelled from BLAST, "
              f"{stats.get('conflict_outcome:conflict_unresolved', 0)} excluded as unresolved")
    return stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate gB_matrix.tsv with segment information from query_tophit_unique_annotated.tsv')
    parser.add_argument('-g', '--gb_matrix', required=True, help='Path to gB_matrix.tsv file')
    parser.add_argument('-s', '--blast_segment', required=True, help='Path to query_tophit_unique_annotated.tsv file (BLAST output)')
    parser.add_argument('-o', '--output_file', help='Path to output file (optional, defaults to overwriting input file)')
    parser.add_argument('-r', '--ref_list', default=None, help='Reference list TSV; its non-exclusion rows define the valid segment labels (without it, the legacy influenza 1-8 range check is used)')
    parser.add_argument('--overwrite_exclusions', action='store_true', help='Force overwrite existing exclusion values')
    parser.add_argument('--segment_names', default=None,
                        help='Segment-name asset TSV (segment/segment_name/aliases). Populates the derived '
                             'segment_name column and lets a submitter-declared segment name (GISAID ships HA/NA/PB1) '
                             'be translated to the canonical segment number.')
    parser.add_argument('--min_conflict_pident', type=float, default=MIN_CONFLICT_PIDENT,
                        help='When the submitter\'s segment disagrees with BLAST, BLAST wins only at or above this '
                             'percent identity; below it the row is excluded (default %(default)s)')
    parser.add_argument('--min_conflict_length', type=int, default=MIN_CONFLICT_LENGTH,
                        help='...and only for sequences at least this long (default %(default)s)')
    parser.add_argument('--conflicts_out', default=None,
                        help='Write every submitter-vs-BLAST disagreement and its outcome to this TSV')
    args = parser.parse_args()

    matrix_file = args.gb_matrix
    annotated_file = args.blast_segment
    output_file = args.output_file if args.output_file else matrix_file
    overwrite = args.output_file is None

    segment_map = build_segment_map(annotated_file)
    reference_map = build_reference_map(annotated_file)
    valid_segments = load_valid_segments(args.ref_list) if args.ref_list else None
    if valid_segments:
        print(f"Valid segment labels from {args.ref_list}: {sorted(valid_segments)}")

    annotate_matrix(matrix_file, segment_map, reference_map, output_file, overwrite=overwrite,
                    overwrite_exclusions=args.overwrite_exclusions, valid_segments=valid_segments,
                    segment_names_path=args.segment_names,
                    identity_map=build_identity_map(annotated_file),
                    min_pident=args.min_conflict_pident, min_length=args.min_conflict_length,
                    conflicts_out=args.conflicts_out)
