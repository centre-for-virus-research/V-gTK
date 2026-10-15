import os
import csv
import datetime
import re
import shutil
from argparse import ArgumentParser
from collections import Counter

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

#: Influenza strain names carry the collection year as their last '/'-delimited field:
#: A/Lisboa/12/2021. That makes the strain name an independent witness to the collection
#: date, which is worth having because the date is the field that gets corrupted. A GISAID
#: export where the year was typed into a cell Excel had formatted as a date comes back as
#: a 1905 day-serial -- serial 2021 is 1905-07-13, serial 2025 is 1905-07-17 -- and 277
#: rows of one rebuilt influenza database carried exactly that. The strain name was right
#: in every one of them. Nothing here rewrites a date: a mismatch means one of the two
#: fields is wrong and the script cannot tell which, so it flags and reports.
STRAIN_YEAR_RE = re.compile(r'^\d{4}$')

#: A final field of one or two digits is ambiguous and is never read as a year. Flu names
#: really do carry two-digit years (``A/.../09``), but the last field of
#: ``A/duck/Chile/01`` is an isolate number and nothing in the string separates the two
#: cases. Guessing would put a steady trickle of healthy records into the conflicts
#: report, which is the fastest way to get the report ignored.
SHORT_TAIL_RE = re.compile(r'^\d{1,2}$')

#: Influenza sequencing starts here; a four-digit final field earlier than this is an
#: isolate number that happens to have four digits, not a year.
EARLIEST_PLAUSIBLE_YEAR = 1918


def strain_year_from_name(parsed_strain: str, this_year: int):
    """Return the four-digit collection year a flu strain name claims, or None.

    None is returned for anything that is not unambiguously a year, rather than a
    guess -- a strain whose last field is an isolate number must not be read as a date.
    """
    if not parsed_strain:
        return None
    last = parsed_strain.rstrip('/_').split('/')[-1].strip()
    if not STRAIN_YEAR_RE.match(last):
        return None
    value = int(last)
    if not (EARLIEST_PLAUSIBLE_YEAR <= value <= this_year):
        return None
    return value


def compare_strain_year_to_collection(row, this_year: int):
    """Classify one row's strain name against its collection year.

    Returns (verdict, strain_year, collection_year). The verdict is one of ``ok``,
    ``mismatch``, ``ambiguous_strain_year``, ``no_strain_year`` or
    ``no_collection_year``; only ``mismatch`` is a defect, the rest record that the
    cross-check could not be run on that row.
    """
    parsed_strain = row.get('Parsed_strain', '') or ''
    strain_year = strain_year_from_name(parsed_strain, this_year)
    raw_collection = str(row.get('collection_year', '')).strip()
    collection_year = int(raw_collection) if raw_collection.isdigit() else None

    if strain_year is None:
        last = parsed_strain.rstrip('/_').split('/')[-1].strip()
        if SHORT_TAIL_RE.match(last):
            return 'ambiguous_strain_year', None, collection_year
        return 'no_strain_year', None, collection_year
    if collection_year is None:
        return 'no_collection_year', strain_year, collection_year
    # A sample collected in December and named for the season it was reported in is a
    # routine one-year offset, not a corrupt date. Only a wider gap is flagged.
    if abs(strain_year - collection_year) <= 1:
        return 'ok', strain_year, collection_year
    return 'mismatch', strain_year, collection_year


def process_matrix(input_path: str, output_path: str, overwrite: bool = False, serotype_replacements=None,
                   conflicts_out: str = None):
    tmp_output = output_path if not overwrite else output_path + ".tmp"
    this_year = datetime.date.today().year
    agreement_counts = Counter()
    conflicts = []

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
        if 'strain_year' not in fieldnames:
            fieldnames.append('strain_year')
        if 'strain_date_agreement' not in fieldnames:
            fieldnames.append('strain_date_agreement')

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

                    # --- Strain name vs collection date ---
                    verdict, strain_year, collection_year = compare_strain_year_to_collection(
                        row, this_year)
                    row['strain_year'] = strain_year if strain_year is not None else ''
                    row['strain_date_agreement'] = verdict
                    agreement_counts[verdict] += 1
                    if verdict == 'mismatch':
                        conflicts.append({
                            'accession': row.get('primary_accession', '') or row.get('locus', ''),
                            'Parsed_strain': row['Parsed_strain'],
                            'collection_date': row.get('collection_date', ''),
                            'collection_year': collection_year,
                            'strain_year': strain_year,
                            'years_apart': strain_year - collection_year,
                        })

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

    if conflicts_out:
        conflict_fields = ['accession', 'Parsed_strain', 'collection_date',
                           'collection_year', 'strain_year', 'years_apart']
        with open(conflicts_out, 'w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=conflict_fields, delimiter='\t')
            writer.writeheader()
            # Worst first: the bigger the gap the more likely it is a corrupt date
            # rather than a season-boundary sample.
            for conflict in sorted(conflicts, key=lambda c: -abs(c['years_apart'])):
                writer.writerow(conflict)

    print("Strain name vs collection date:")
    print(f"- Agreeing (within one year): {agreement_counts['ok']}")
    print(f"- MISMATCHED: {agreement_counts['mismatch']}")
    print(f"- No usable year in the strain name: {agreement_counts['no_strain_year']}")
    print(f"- Strain name ends in an ambiguous 1-2 digit field: "
          f"{agreement_counts['ambiguous_strain_year']}")
    print(f"- No collection year to compare against: {agreement_counts['no_collection_year']}")
    if conflicts:
        excel = [c for c in conflicts if str(c['collection_date']).startswith('190')]
        if excel:
            print(f"  NOTE: {len(excel)} mismatched rows carry a date in the 1900s, which is "
                  "what a collection year becomes when read as an Excel 1900 day-serial.")
        if conflicts_out:
            print(f"  Written to {conflicts_out}")

if __name__ == "__main__":
    parser = ArgumentParser(description='Extract or validate Parsed_strain and serotype_validated fields from gB_matrix based on strain/serotype or definition columns.')
    parser.add_argument('-d', '--tmp_dir', help='Path to save validated gB_matrix', default='tmp/Validate-matrix')
    parser.add_argument('-g', '--gb_matrix', help='Genbank matrix file', default='tmp/Validate-matrix/gB_matrix_validated.tsv')
    parser.add_argument('-o', '--output_file', help='Output file path (optional, defaults to overwriting input file)', default=None)
    parser.add_argument('-m','--mapping_file', help='TSV file with raw and curated serotype values (default: ./generic-influenza/serotype_mapping.tsv)', default='./generic-influenza/serotype_mapping.tsv')
    parser.add_argument('--conflicts_out', help='TSV listing rows whose strain-name year disagrees with their collection year. Flagged, never corrected.', default=None)
    args = parser.parse_args()

    os.makedirs(args.tmp_dir, exist_ok=True)

    input_path = args.gb_matrix
    output_path = args.output_file if args.output_file else input_path
    overwrite = args.output_file is None

    serotype_replacements = load_serotype_replacements(args.mapping_file)
    process_matrix(input_path, output_path, overwrite=overwrite, serotype_replacements=serotype_replacements,
                   conflicts_out=args.conflicts_out)
