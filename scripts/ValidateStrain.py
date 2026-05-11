import os
import csv
import re
import shutil
from argparse import ArgumentParser

def parse_strain_from_definition(definition: str) -> str:
    """Extract strain name from definition field using regex patterns."""
    pattern = r"(?<=\().*?(?=\))"
    second_pattern = r"^(.*?)(?=\(H|\(mixed)"

    match = re.search(pattern, definition)
    parsed_text = match.group(0) if match else ''
    if parsed_text:
        second_match = re.match(second_pattern, parsed_text)
        if second_match:
            parsed_text = second_match.group(1)
    return parsed_text

def clean_strain_name(strain: str) -> str:
    """Clean strain name by removing virus subtype patterns and normalizing characters."""
    strain = re.sub(r'\(H\d{1}N\d{1}\)', '', strain)
    strain = re.sub(r'[^a-zA-Z0-9/_-]', '_', strain)
    strain = re.sub(r'_+', '_', strain)
    strain = re.sub(r'[_/]+$', '', strain)
    return strain

def extract_serotype(row):
    """Extract serotype from column or parse from definition."""
    if 'serotype' in row and row['serotype'].strip():
        return row['serotype'].strip()
    definition = row.get('definition', '')
    match = re.search(r'\((H\d{1,2}N\d{1,2})\)', definition)
    if match:
        return match.group(1)
    return ''

def load_serotype_replacements(replacement_file):
    """Load serotype replacements from a TSV file (two columns)."""
    replacements = {}
    with open(replacement_file, encoding="utf-8") as fin:
        reader = csv.reader(fin, delimiter='\t')
        for row in reader:
            # Salta header si lo hay
            if len(row) == 0 or row[0].startswith('#') or row[0].lower() == 'raw_value':
                continue
            if len(row) < 2:
                continue
            raw, curated = row[0].strip(), row[1].strip()
            if raw:
                replacements[raw] = curated
    return replacements

def curate_serotype(raw_serotype, replacements=None):
    """Normalize serotype string"""
    if not raw_serotype:
        return ''
    # Normalize case and remove odd characters
    serotype = raw_serotype.upper().replace(" ", "").replace("-", "").replace("/", "")

    # Clean GISAID serotype "A / H1N1" or "A/H3N2" -> return HxNy (or Hx / Ny)
    m = re.search(r'^[ABCD]\s*/\s*(H\d{1,2}(?:N\d{1,2})?|N\d{1,2})$', raw_serotype.strip(), re.IGNORECASE)
    if m:
        return m.group(1).upper()
        
    # Step 1: Regex extraction (priority order)
    match = re.search(r'H\d{1,2}N\d{1,2}', serotype)
    if match:
        return match.group(0)
    match = re.search(r'H\d{1,2}', serotype)
    if match:
        return match.group(0)
    match = re.search(r'N\d{1,2}', serotype)
    if match:
        return match.group(0)
    match = re.search(r'H\?N(\d{1,2})', serotype)
    if match:
        return f'N{match.group(1)}'
    match = re.search(r'H(\d{1,2})N\?', serotype)
    if match:
        return f'H{match.group(1)}'
    match = re.search(r'H(\d{1,2})NX', serotype)
    if match:
        return f'H{match.group(1)}'

    # Step 2: replacements using diccionary (provided by user)
    if replacements and serotype in replacements:
        return replacements[serotype]

    # Step 3: Mark as unknown
    return 'unknown'

def process_matrix(input_path: str, output_path: str, overwrite: bool = False, serotype_replacements=None):
    tmp_output = output_path if not overwrite else output_path + ".tmp"

    from_strain = 0
    from_definition = 0
    empty_count = 0

    from_serotype_col = 0
    from_serotype_def = 0
    empty_serotype = 0

    row_count = 0

    with open(input_path, 'r', encoding='utf-8', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        
        fieldnames = list(reader.fieldnames)
        if 'Parsed_strain' not in fieldnames:
            fieldnames.append('Parsed_strain')
        if 'serotype_validated' not in fieldnames:
            fieldnames.append('serotype_validated')

        with open(tmp_output, 'w', encoding='utf-8', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            for idx, row in enumerate(reader, start=1):
                try:
                    if not row:
                        print(f"[WARNING] Empty row at line {idx + 1}. Skipping.")
                        continue

                    # --- Strain ---
                    parsed_strain = ''
                    if 'strain' in row and row['strain'].strip():
                        parsed_strain = row['strain'].strip()
                        from_strain += 1
                    else:
                        definition = row.get('definition', '').strip()
                        parsed_strain = parse_strain_from_definition(definition)
                        if parsed_strain:
                            from_definition += 1
                        else:
                            empty_count += 1
                    row['Parsed_strain'] = clean_strain_name(parsed_strain)

                    # --- Serotype ---
                    raw_serotype = extract_serotype(row)
                    if raw_serotype:
                        if 'serotype' in row and row['serotype'].strip():
                            from_serotype_col += 1
                        else:
                            from_serotype_def += 1
                        row['serotype_validated'] = curate_serotype(raw_serotype, replacements=serotype_replacements)
                    else:
                        empty_serotype += 1
                        row['serotype_validated'] = 'unknown'

                    writer.writerow(row)
                    row_count += 1

                except Exception as e:
                    print(f"[ERROR] Exception at row {idx + 1}: {e}")
                    print(f"Row content: {row}")
                    continue

    if overwrite:
        shutil.move(tmp_output, output_path)

    print("Strain and serotype parsing finished.")
    print(f"- Total rows written: {row_count}")
    print(f"- Parsed_strain from 'strain': {from_strain}")
    print(f"- Parsed_strain from 'definition': {from_definition}")
    print(f"- Empty Parsed_strain values: {empty_count}")
    print(f"- serotype_validated from 'serotype' column: {from_serotype_col}")
    print(f"- serotype_validated parsed from 'definition': {from_serotype_def}")
    print(f"- Empty serotype_validated values: {empty_serotype} - including exclusion TRUE rows")

if __name__ == "__main__":
    parser = ArgumentParser(description='Extract or validate Parsed_strain and serotype_validated fields from gB_matrix based on strain/serotype or definition columns.')
    parser.add_argument('-d', '--tmp_dir', help='Path to save validated gB_matrix', default='tmp/Validate-matrix')
    parser.add_argument('-g', '--gb_matrix', help='Genbank matrix file', default='tmp/Validate-matrix/gB_matrix_validated.tsv')
    parser.add_argument('-o', '--output_file', help='Output file path (optional, defaults to overwriting input file)', default=None)
    parser.add_argument('-m','--mapping_file', help='TSV file with raw and curated serotype values (default: ./generic-influenza/serotype_mapping.tsv)', default='./generic-influenza/serotype_mapping.tsv')
    args = parser.parse_args()

    os.makedirs(args.tmp_dir, exist_ok=True)

    input_path = args.gb_matrix
    output_path = args.output_file if args.output_file else input_path
    overwrite = args.output_file is None

    serotype_replacements = load_serotype_replacements(args.mapping_file)
    process_matrix(input_path, output_path, overwrite=overwrite, serotype_replacements=serotype_replacements)
