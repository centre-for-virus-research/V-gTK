"""Tabulate isolates/strains against the segments they carry, for any segmented virus.

Generalisation of the former influenza-only FluPivotTable. Two things used to be
hardwired to influenza A:

  * the expected segment set was the module constant ['1'..'8'], so a 2-segment
    arenavirus or an 11-segment rotavirus could never be "Complete";
  * the grouping key was ``Parsed_strain``, a column produced only by
    ValidateStrain.py, which the workflow runs only when ``params.is_flu == "Y"``.

Here the expected segment set is derived from the reference list (the master rows
define one segment each), and the grouping key is elected from the columns the
matrix actually carries: Parsed_strain for flu, otherwise the GenBank ``isolate``
or ``strain`` qualifier.

With no reference list and no explicit segment set, the influenza 1-8 default is
used, so the legacy CLI reproduces the old output exactly.
"""

import csv
import re
from argparse import ArgumentParser
from collections import OrderedDict
from os import makedirs
from os.path import dirname, exists
import segment_utils

# The influenza A/B genome. Retained as the fallback when no reference list is
# supplied, which is what keeps the legacy CLI byte-compatible.
DEFAULT_REQUIRED_SEGMENTS = ['1', '2', '3', '4', '5', '6', '7', '8']
required_segments = DEFAULT_REQUIRED_SEGMENTS  # backwards-compatible alias

# Ordered candidates for the isolate grouping key. Parsed_strain only exists on
# flu runs; isolate and strain come from GenBankParser and exist for every virus.
DEFAULT_ISOLATE_KEY_COLUMNS = ('Parsed_strain', 'isolate', 'strain')

# Values in an `exclusion` column that mean "not excluded".
DEFAULT_EXCLUSION_FALSE_VALUES = ('', '0', 'false', 'no', 'none', 'na', 'nan', '-')

# Segment values that mean "no segment could be assigned".
UNAVAILABLE_SEGMENT_VALUES = frozenset(
    {'', 'not found', 'na', 'n/a', 'none', 'unknown', 'nan', '-'}
)

RESERVED_OUTPUT_COLUMNS = (
    'Complete_status',
    'isolate_key_source',
    'Segments_present',
    'Missing_segments',
    'Duplicate_segments',
)

EXTENDED_COLUMNS = (
    'isolate_key_source',
    'Segments_present',
    'Missing_segments',
    'Duplicate_segments',
)

# Length-bounded: without it, int(float(text)) on a long digit string raises
# OverflowError, breaking this module's documented "never raises" contract.
_NUMERIC_SEGMENT_RE = re.compile(r'^-?\d{1,9}(\.0+)?$')


def normalise_segment_label(value):
    """Canonicalise a segment label, or return None if no segment is available.

    Numeric labels collapse to their integer form ('4', '4.0', ' 4 ', '08' -> '4')
    so the influenza matrix keeps pivoting as it always has. Non-numeric labels
    (L, M, S, HA, ...) survive verbatim. Never raises - the old
    ``str(int(float(v)))`` killed the process on the first letter segment.

    Normalisation itself is delegated to :mod:`segment_utils` so that this and
    the five other former copies cannot drift apart again. What stays local is
    the *availability* judgement: this function is the only one that treats
    'na', 'unknown', 'not found' and friends as "no segment", which is right for
    a completeness table keyed on numeric segments but would silently discard
    every neuraminidase row of a virus keyed on segment names.
    """
    if value is None:
        return None
    text = str(value).strip()
    if text.casefold() in UNAVAILABLE_SEGMENT_VALUES:
        return None
    # This local numeric branch is kept because it accepts forms the shared
    # normaliser deliberately rejects as labels - notably a bare '-1', and digit
    # strings longer than the shared normaliser's length bound. Both then agree.
    if _NUMERIC_SEGMENT_RE.match(text):
        return str(int(float(text)))
    return segment_utils.normalise_segment(text)


def natural_segment_sort(labels):
    """Numeric segments ascending, then non-numeric ones case-insensitively.

    Keeps an 11-segment virus as 1..11 rather than the lexicographic 1,10,11,2...
    """
    def sort_key(label):
        try:
            return (0, float(label), '')
        except (TypeError, ValueError):
            return (1, 0.0, str(label).casefold())

    return sorted(labels, key=sort_key)


def _unique_segments(values):
    seen = OrderedDict()
    for value in values:
        label = normalise_segment_label(value)
        if label is not None:
            seen.setdefault(label, None)
    return list(seen)


def required_segments_from_ref_list(ref_list_path, source='master'):
    """Derive the expected segment set from the reference list.

    Returns ``(segments, source_label)``; ``segments`` is empty when the reference
    list carries no segment column (e.g. a 2-column list), leaving the caller to
    fall back.

    Master rows are authoritative: ValidateRefList guarantees exactly one master
    per real segment, so the set is exactly the segments this build targets. It
    also drops influenza's B/C/D decoy rows, which are ``exclusion_list`` entries
    parked in the segment column - counting those would make every flu isolate
    permanently Incomplete.
    """
    # Imported lazily so this module stays stdlib-only unless a ref_list is used.
    from ExportRefListFromUpdateDb import load_reference_file_table

    refs = load_reference_file_table(ref_list_path)
    acc_type = refs['accession_type'].astype(str).str.strip().str.lower()

    if source == 'all':
        segments = _unique_segments(refs['segment'])
        return natural_segment_sort(segments), 'ref_list:all'

    if source == 'all_non_exclusion':
        segments = _unique_segments(refs[acc_type != 'exclusion_list']['segment'])
        return natural_segment_sort(segments), 'ref_list:non_exclusion'

    segments = _unique_segments(refs[acc_type == 'master']['segment'])
    if segments:
        return natural_segment_sort(segments), 'ref_list:master'

    # Legacy reference lists without master rows.
    segments = _unique_segments(refs[acc_type != 'exclusion_list']['segment'])
    return natural_segment_sort(segments), 'ref_list:non_exclusion'


def elect_isolate_key_column(rows, candidates, min_coverage=0.5):
    """Pick one column to group every row by, and report each candidate's coverage.

    One column is elected for the whole file rather than falling back per row:
    cross-segment grouping only works if every segment of a specimen is keyed the
    same way. Keying segment L by `isolate` and segment S by `strain` would never
    join them, and could collide two different namespaces on an equal string.
    """
    coverages = OrderedDict()
    total = len(rows)
    for column in candidates:
        if rows and column not in rows[0]:
            continue
        if not rows:
            coverages[column] = 0.0
            continue
        populated = sum(1 for row in rows if str(row.get(column) or '').strip())
        coverages[column] = populated / total if total else 0.0

    for column, coverage in coverages.items():
        if coverage >= min_coverage:
            return column, coverage, coverages

    best = None
    for column, coverage in coverages.items():
        if coverage > 0 and (best is None or coverage > coverages[best]):
            best = column
    if best is not None:
        print(
            f"[WARNING] elected isolate key column '{best}' with only "
            f"{coverages[best] * 100:.0f}% coverage"
        )
        return best, coverages[best], coverages

    return None, 0.0, coverages


def _normalise_key(value, mode):
    text = str(value or '').strip()
    if mode in ('collapse', 'casefold'):
        text = ' '.join(text.split())
    return text


def _is_excluded(value, false_values):
    return str(value or '').strip().casefold() not in false_values


def _resolve_required_segments(required, ref_list, reflist_segment_source):
    if required:
        return [str(segment) for segment in required], 'argument'
    if ref_list:
        segments, source = required_segments_from_ref_list(ref_list, reflist_segment_source)
        if segments:
            return segments, source
    print(
        '[WARNING] no reference list segments available; falling back to the '
        'influenza segment set 1-8'
    )
    return list(DEFAULT_REQUIRED_SEGMENTS), 'default'


def _write_summary(summary_file, summary):
    if not summary_file:
        return
    output_dir = dirname(summary_file)
    if output_dir and not exists(output_dir):
        makedirs(output_dir)
    with open(summary_file, mode='w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['key', 'value'])
        for key, value in summary.items():
            writer.writerow([key, value])


def pivot_data(
    input_file,
    output_file,
    required_segments=None,
    *,
    ref_list=None,
    reflist_segment_source='master',
    isolate_key_columns=None,
    isolate_column_name=None,
    min_key_coverage=0.5,
    blank_isolate='accession',
    key_normalise='strip',
    segment_column='segment_validated',
    exclusion_column='exclusion',
    exclusion_false_values=DEFAULT_EXCLUSION_FALSE_VALUES,
    exclusion_status_column=None,
    accession_types=None,
    include_unexpected_segments=False,
    extended_columns=False,
    summary_file=None,
    require_columns=None,
):
    """Pivot the long-form matrix to one row per isolate, one column per segment.

    Returns the summary dict that is also written to ``summary_file``.
    """
    candidates = tuple(isolate_key_columns or DEFAULT_ISOLATE_KEY_COLUMNS)
    false_values = {str(v).strip().casefold() for v in exclusion_false_values}
    allowed_types = (
        {str(t).strip().lower() for t in accession_types} if accession_types else None
    )

    with open(input_file, mode='r', newline='') as csv_file:
        reader = csv.DictReader(csv_file, delimiter='\t')
        header = list(reader.fieldnames or [])
        rows = list(reader)

    needed = set(require_columns) if require_columns else {'primary_accession', segment_column}
    missing = sorted(needed - set(header))
    if missing:
        raise ValueError(
            'Input file is missing one or more required columns: '
            + ', '.join(repr(column) for column in missing)
        )
    if not require_columns and not any(column in header for column in candidates):
        raise ValueError(
            'Input file is missing one or more required columns: no isolate key '
            'column found among ' + ', '.join(repr(column) for column in candidates)
        )

    n_rows_read = len(rows)
    n_rows_excluded = 0
    groupable = []
    for row in rows:
        if allowed_types is not None:
            if str(row.get('accession_type') or '').strip().lower() not in allowed_types:
                continue
        if exclusion_column and exclusion_column in header:
            if _is_excluded(row.get(exclusion_column), false_values):
                n_rows_excluded += 1
                continue
        if exclusion_status_column and exclusion_status_column in header:
            if _is_excluded(row.get(exclusion_status_column), false_values):
                n_rows_excluded += 1
                continue
        groupable.append(row)

    segments, segments_source = _resolve_required_segments(
        required_segments, ref_list, reflist_segment_source
    )
    print(f"Required segments ({segments_source}): {segments}")

    key_column, key_coverage, key_coverages = elect_isolate_key_column(
        groupable, candidates, min_key_coverage
    )
    key_header = isolate_column_name or key_column or 'isolate_key'

    segment_lookup = {label.casefold(): label for label in segments}
    reserved = set(RESERVED_OUTPUT_COLUMNS) | {key_header}
    clashes = sorted(set(segments) & reserved)
    if clashes:
        raise ValueError(
            'Segment label(s) collide with reserved output column names: '
            + ', '.join(clashes)
        )

    data_wide = OrderedDict()
    key_sources = {}
    n_rows_blank_key = 0
    n_rows_no_segment = 0
    n_rows_unexpected = 0
    unexpected_labels = OrderedDict()
    n_duplicate_cells = 0

    for row in groupable:
        primary_accession = str(row.get('primary_accession') or '').strip()

        raw_key = _normalise_key(row.get(key_column) if key_column else '', key_normalise)
        key_source = key_column or 'none'
        if not raw_key:
            n_rows_blank_key += 1
            if blank_isolate == 'drop':
                continue
            if blank_isolate == 'accession':
                raw_key = primary_accession
                key_source = 'accession_fallback'
            else:  # 'group' - the legacy behaviour, everything blank in one bucket
                key_source = 'blank'
        group_key = raw_key.casefold() if key_normalise == 'casefold' else raw_key

        label = normalise_segment_label(row.get(segment_column))
        if label is None:
            n_rows_no_segment += 1
            print(f"[WARNING] Empty {segment_column} for: {primary_accession}, {raw_key}")
            continue

        canonical = segment_lookup.get(label.casefold())
        if canonical is None:
            n_rows_unexpected += 1
            unexpected_labels.setdefault(label, None)
            if not include_unexpected_segments:
                continue
            canonical = label
            segment_lookup[label.casefold()] = label

        if group_key not in data_wide:
            data_wide[group_key] = (raw_key, OrderedDict())
            key_sources[group_key] = key_source
        accessions = data_wide[group_key][1].setdefault(canonical, [])
        if primary_accession and primary_accession not in accessions:
            if accessions:
                n_duplicate_cells += 1
                print(
                    f"[WARNING] duplicate accessions for {raw_key} segment "
                    f"{canonical}: {','.join(accessions + [primary_accession])}"
                )
            accessions.append(primary_accession)

    extra_segments = [label for label in segment_lookup.values() if label not in segments]
    output_segments = list(segments) + natural_segment_sort(extra_segments)

    complete_count = 0
    incomplete_count = 0

    output_dir = dirname(output_file)
    if output_dir and not exists(output_dir):
        makedirs(output_dir)

    with open(output_file, mode='w', newline='') as csv_output_file:
        fieldnames = [key_header] + output_segments
        if extended_columns:
            fieldnames += list(EXTENDED_COLUMNS)
        fieldnames.append('Complete_status')
        writer = csv.DictWriter(csv_output_file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for group_key, (raw_key, loci_dict) in data_wide.items():
            row = {key_header: raw_key}
            for segment in output_segments:
                row[segment] = ','.join(loci_dict.get(segment, []))

            present = [s for s in segments if loci_dict.get(s)]
            complete = len(present) == len(segments)

            if extended_columns:
                duplicates = [
                    f"{s}:{len(loci_dict[s])}"
                    for s in output_segments
                    if len(loci_dict.get(s, [])) > 1
                ]
                row['isolate_key_source'] = key_sources.get(group_key, '')
                row['Segments_present'] = ','.join(present)
                row['Missing_segments'] = ','.join(s for s in segments if not loci_dict.get(s))
                row['Duplicate_segments'] = ';'.join(duplicates)

            row['Complete_status'] = 'Complete' if complete else 'Incomplete'

            if complete:
                complete_count += 1
            else:
                incomplete_count += 1

            writer.writerow(row)

    if n_rows_unexpected and groupable and n_rows_unexpected > 0.2 * len(groupable):
        print(
            f"[WARNING] {n_rows_unexpected}/{len(groupable)} rows carry a segment "
            f"outside the expected set {segments}: {list(unexpected_labels)}. If the "
            'reference list has more than three columns, check BlastAlignment.py:466 '
            "(`segment = parts[-1]`) - it reads the LAST column, not the segment one."
        )

    print(f"Number of Complete genomes: {complete_count}")
    print(f"Number of Incomplete genomes: {incomplete_count}")

    summary = OrderedDict([
        ('required_segments', ','.join(segments)),
        ('required_segments_source', segments_source),
        ('ref_list', ref_list or ''),
        ('isolate_key_column', key_column or ''),
        ('isolate_key_coverage', f"{key_coverage:.4f}"),
        ('isolate_key_candidates',
         ','.join(f"{c}={v:.4f}" for c, v in key_coverages.items())),
        ('n_rows_read', n_rows_read),
        ('n_rows_excluded', n_rows_excluded),
        ('n_rows_blank_key', n_rows_blank_key),
        ('n_rows_no_segment', n_rows_no_segment),
        ('n_rows_unexpected_segment', n_rows_unexpected),
        ('unexpected_segment_labels', ','.join(unexpected_labels)),
        ('n_duplicate_isolate_segment_cells', n_duplicate_cells),
        ('n_isolates', complete_count + incomplete_count),
        ('n_complete', complete_count),
        ('n_incomplete', incomplete_count),
    ])
    _write_summary(summary_file, summary)
    return summary


if __name__ == "__main__":
    parser = ArgumentParser(
        description='Compile the per-isolate accession numbers per segment and flag '
                    'which isolates carry a complete genome. Requires a deduplicated '
                    'GenBank matrix.'
    )
    parser.add_argument('-g', '--gb_matrix', help='Genbank matrix (deduplicated) to be processed', default='./tmp/Validate-matrix/gB_matrix_validated.tsv')
    parser.add_argument('-o', '--output_file', help='Output file path for the pivot strains matrix', default='tmp/Validate-matrix/gB_matrix_cast_strains.tsv')
    parser.add_argument('-r', '--ref_list', default=None,
                        help='Reference list TSV; its master rows define the expected segment set')
    parser.add_argument('--required_segments', default=None,
                        help='Comma-separated expected segment set; overrides --ref_list')
    parser.add_argument('--reflist_segment_source', default='master',
                        choices=['master', 'all_non_exclusion', 'all'],
                        help='Which reference list rows define the expected segment set')
    parser.add_argument('-k', '--isolate_key', default=','.join(DEFAULT_ISOLATE_KEY_COLUMNS),
                        help='Ordered candidate columns for the isolate grouping key')
    parser.add_argument('--isolate_column_name', default=None,
                        help='Name for the first output column (defaults to the elected key column)')
    parser.add_argument('--min_key_coverage', type=float, default=0.5,
                        help='Minimum non-blank fraction for a key column to be elected')
    parser.add_argument('--blank_isolate', default='accession', choices=['accession', 'group', 'drop'],
                        help='What to do with a row whose isolate key is blank')
    parser.add_argument('--key_normalise', default='strip', choices=['strip', 'collapse', 'casefold'],
                        help='Whitespace/case normalisation applied to the isolate key')
    parser.add_argument('--segment_column', default='segment_validated',
                        help='Matrix column carrying the segment label')
    parser.add_argument('--exclusion_column', default='exclusion',
                        help="Matrix column carrying the exclusion reason ('' disables the filter)")
    parser.add_argument('--exclusion_false_values', default=','.join(DEFAULT_EXCLUSION_FALSE_VALUES),
                        help='Comma-separated values in the exclusion column that mean "not excluded"')
    parser.add_argument('--exclusion_status_column', default=None,
                        help='Additionally drop rows excluded by this QC column (e.g. exclusion_status)')
    parser.add_argument('--accession_types', default=None,
                        help='Comma-separated accession_type values to keep (default: all)')
    parser.add_argument('--include_unexpected_segments', action='store_true',
                        help='Emit segments outside the expected set as trailing columns')
    parser.add_argument('--extended_columns', action='store_true',
                        help='Add isolate_key_source/Segments_present/Missing_segments/Duplicate_segments')
    parser.add_argument('--summary_file', default=None,
                        help='Path for the run summary sidecar (default: <output_file>.summary.tsv)')
    args = parser.parse_args()

    output_path = args.output_file
    summary_path = args.summary_file or (output_path + '.summary.tsv')

    pivot_data(
        args.gb_matrix,
        output_path,
        [s.strip() for s in args.required_segments.split(',') if s.strip()] if args.required_segments else None,
        ref_list=args.ref_list,
        reflist_segment_source=args.reflist_segment_source,
        isolate_key_columns=[c.strip() for c in args.isolate_key.split(',') if c.strip()],
        isolate_column_name=args.isolate_column_name,
        min_key_coverage=args.min_key_coverage,
        blank_isolate=args.blank_isolate,
        key_normalise=args.key_normalise,
        segment_column=args.segment_column,
        exclusion_column=args.exclusion_column,
        exclusion_false_values=args.exclusion_false_values.split(','),
        exclusion_status_column=args.exclusion_status_column,
        accession_types=[t.strip() for t in args.accession_types.split(',')] if args.accession_types else None,
        include_unexpected_segments=args.include_unexpected_segments,
        extended_columns=args.extended_columns,
        summary_file=summary_path,
    )

    print(f"Pivoted data has been written to {output_path}")
