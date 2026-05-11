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

def annotate_matrix(matrix_file, segment_map, reference_map, output_file, overwrite=False, overwrite_exclusions=False):
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

            row['segment_validated'] = segment_map.get(primary_accession, 'not found')
            row['closest_reference'] = reference_map.get(primary_accession, 'not found')

            if overwrite_exclusions or not row.get('exclusion'):
                if row['segment_validated'] == 'not found':
                    row['exclusion'] = 'not significant BLAST hit'
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
    parser.add_argument('--overwrite_exclusions', action='store_true', help='Force overwrite existing exclusion values')
    args = parser.parse_args()

    matrix_file = args.gb_matrix
    annotated_file = args.blast_segment
    output_file = args.output_file if args.output_file else matrix_file
    overwrite = args.output_file is None

    segment_map = build_segment_map(annotated_file)
    reference_map = build_reference_map(annotated_file)

    annotate_matrix(matrix_file, segment_map, reference_map, output_file, overwrite=overwrite, overwrite_exclusions=args.overwrite_exclusions)
