import csv
import argparse
import os
import shutil

'''Annotate gB_matrix with BLAST results. Includes overwrite input file (default) and overwrite exclusion options'''

def build_segment_map(file_path):
    segment_map = {}
    with open(file_path, 'r') as ann_file:
        reader = csv.reader(ann_file, delimiter='\t')
        for row in reader:
            segment_map[row[0]] = row[4]
    return segment_map

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


def annotate_matrix(matrix_file, segment_map, reference_map, output_file, overwrite=False, overwrite_exclusions=False, valid_segments=None):
    tmp_output = output_file if not overwrite else output_file + ".tmp"

    with open(matrix_file, 'r') as matrix, open(tmp_output, 'w', newline='') as output:
        reader = csv.DictReader(matrix, delimiter='\t')

        fieldnames = reader.fieldnames.copy()
        if 'segment_validated' not in fieldnames:
            fieldnames.append('segment_validated')
        if 'closest_reference' not in fieldnames:
            fieldnames.append('closest_reference')
        if 'exclusion' not in fieldnames:
            fieldnames.append('exclusion')

        writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            primary_accession = row['primary_accession']

            if _is_reference_like(row):
                row['segment_validated'] = row.get('segment_validated', '') or row.get('segment', '')
                row['closest_reference'] = row.get('closest_reference', '') or primary_accession
                writer.writerow(row)
                continue

            row['segment_validated'] = segment_map.get(primary_accession, 'not found')
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

            writer.writerow(row)

    if overwrite:
        shutil.move(tmp_output, output_file)

    print(f"Annotated file written to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate gB_matrix.tsv with segment information from query_tophit_unique_annotated.tsv')
    parser.add_argument('-g', '--gb_matrix', required=True, help='Path to gB_matrix.tsv file')
    parser.add_argument('-s', '--blast_segment', required=True, help='Path to query_tophit_unique_annotated.tsv file (BLAST output)')
    parser.add_argument('-o', '--output_file', help='Path to output file (optional, defaults to overwriting input file)')
    parser.add_argument('-r', '--ref_list', default=None, help='Reference list TSV; its non-exclusion rows define the valid segment labels (without it, the legacy influenza 1-8 range check is used)')
    parser.add_argument('--overwrite_exclusions', action='store_true', help='Force overwrite existing exclusion values')
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

    annotate_matrix(matrix_file, segment_map, reference_map, output_file, overwrite=overwrite, overwrite_exclusions=args.overwrite_exclusions, valid_segments=valid_segments)
