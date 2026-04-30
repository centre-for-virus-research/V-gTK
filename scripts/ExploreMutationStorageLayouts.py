#!/usr/bin/env python3

import argparse
import sqlite3
from pathlib import Path

import pandas as pd


REQUIRED_SEQUENCE_MUTATION_COLUMNS = {
    'primary_accession',
    'mutation_id',
    'protein_name',
    'segment',
    'aa_position',
    'alt_residue',
    'combination_id',
}

REQUIRED_MUTATION_CATALOG_COLUMNS = {
    'mutation_id',
    'protein_name',
    'segment',
    'aa_position',
    'alt_residue',
    'reference_accession',
    'signature_id',
    'signature_kind',
    'combination_id',
}


def _require_columns(df, required_columns, table_name):
    missing = sorted(required_columns - set(df.columns))
    if missing:
        raise ValueError(f"Table '{table_name}' is missing required columns: {', '.join(missing)}")


def _normalize_text(value):
    if pd.isna(value):
        return ''
    return str(value).strip()


def _normalize_df_strings(df):
    normalized = df.copy()
    for column in normalized.columns:
        normalized[column] = normalized[column].apply(_normalize_text)
    return normalized


def _format_human_bytes(byte_count):
    units = ['B', 'KB', 'MB', 'GB', 'TB']
    size = float(byte_count)
    unit_index = 0
    while size >= 1024 and unit_index < len(units) - 1:
        size /= 1024
        unit_index += 1
    if unit_index == 0:
        return f'{int(size):,} {units[unit_index]}'
    return f'{size:,.2f} {units[unit_index]}'


def _table_exists(conn, table_name):
    cursor = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name = ?",
        (table_name,),
    )
    return cursor.fetchone() is not None


def infer_mutation_catalog_from_sequence_mutations(sequence_mutations):
    singles = sequence_mutations[sequence_mutations['combination_id'] == ''].copy()
    single_records = []
    if not singles.empty:
        singles = singles.sort_values(['mutation_id'])
        for _, row in singles.drop_duplicates(subset=['mutation_id']).iterrows():
            single_records.append(
                {
                    'mutation_id': row['mutation_id'],
                    'protein_name': row['protein_name'],
                    'segment': row['segment'],
                    'aa_position': row['aa_position'],
                    'alt_residue': row['alt_residue'],
                    'reference_accession': '',
                    'mutation_type': '',
                    'signature_id': row['mutation_id'],
                    'signature_kind': 'single',
                    'combination_id': '',
                    'combination_size': '',
                    'resistance_category': '',
                    'drug': '',
                }
            )

    combination_records = []
    combos = sequence_mutations[sequence_mutations['combination_id'] != ''].copy()
    if not combos.empty:
        combos = combos.sort_values(['combination_id', 'mutation_id'])
        for combination_id, group in combos.groupby('combination_id', sort=False):
            component_mutations = []
            for mutation_id in group['mutation_id'].tolist():
                if mutation_id and mutation_id not in component_mutations:
                    component_mutations.append(mutation_id)
            first_row = group.iloc[0]
            for mutation_id in component_mutations:
                component_row = group[group['mutation_id'] == mutation_id].iloc[0]
                combination_records.append(
                    {
                        'mutation_id': mutation_id,
                        'protein_name': component_row['protein_name'],
                        'segment': component_row['segment'],
                        'aa_position': str(component_row['aa_position']),
                        'alt_residue': component_row['alt_residue'],
                        'reference_accession': '',
                        'mutation_type': '',
                        'signature_id': combination_id,
                        'signature_kind': 'combination',
                        'combination_id': combination_id,
                        'combination_size': str(len(component_mutations)),
                        'resistance_category': '',
                        'drug': '',
                    }
                )

    inferred = pd.DataFrame(single_records + combination_records)
    if inferred.empty:
        inferred = pd.DataFrame(columns=sorted(REQUIRED_MUTATION_CATALOG_COLUMNS | {'combination_size', 'mutation_type', 'resistance_category', 'drug'}))
    return inferred


def load_mutation_tables(db_path, mutation_catalog_tsv=None):
    conn = sqlite3.connect(str(db_path))
    try:
        sequence_mutations = pd.read_sql_query('SELECT * FROM sequence_mutations', conn)
        has_mutation_catalog_table = _table_exists(conn, 'mutation_catalog')
        if has_mutation_catalog_table:
            mutation_catalog = pd.read_sql_query('SELECT * FROM mutation_catalog', conn)
            catalog_source = 'db_table'
        else:
            mutation_catalog = None
            catalog_source = 'missing'
    finally:
        conn.close()

    _require_columns(sequence_mutations, REQUIRED_SEQUENCE_MUTATION_COLUMNS, 'sequence_mutations')

    if mutation_catalog is None and mutation_catalog_tsv:
        mutation_catalog = pd.read_csv(mutation_catalog_tsv, sep='\t', dtype=str)
        catalog_source = 'tsv'

    if mutation_catalog is None:
        mutation_catalog = infer_mutation_catalog_from_sequence_mutations(_normalize_df_strings(sequence_mutations))
        catalog_source = 'inferred_from_sequence_mutations'

    _require_columns(mutation_catalog, REQUIRED_MUTATION_CATALOG_COLUMNS, 'mutation_catalog')

    return _normalize_df_strings(sequence_mutations), _normalize_df_strings(mutation_catalog), catalog_source


def build_combination_catalog(mutation_catalog):
    combo_rows = mutation_catalog[
        (mutation_catalog['signature_kind'] == 'combination')
        & (mutation_catalog['combination_id'] != '')
    ].copy()
    if combo_rows.empty:
        return pd.DataFrame(
            columns=[
                'signature_id',
                'combination_id',
                'protein_name',
                'segment',
                'reference_accession',
                'drug',
                'resistance_category',
                'required_mutation_ids',
                'required_count',
                'combination_size',
            ]
        )

    if 'component_order' in combo_rows.columns:
        sort_columns = ['combination_id', 'component_order', 'mutation_id']
    else:
        sort_columns = ['combination_id', 'mutation_id']
    combo_rows = combo_rows.sort_values(sort_columns)

    records = []
    for combination_id, group in combo_rows.groupby('combination_id', sort=False):
        required_mutations = []
        for mutation_id in group['mutation_id'].tolist():
            if mutation_id and mutation_id not in required_mutations:
                required_mutations.append(mutation_id)

        first_row = group.iloc[0]
        signature_id = first_row['signature_id'] or combination_id
        combination_size = first_row.get('combination_size', '')
        if not combination_size:
            combination_size = str(len(required_mutations))
        records.append(
            {
                'signature_id': signature_id,
                'combination_id': combination_id,
                'protein_name': first_row['protein_name'],
                'segment': first_row['segment'],
                'reference_accession': first_row['reference_accession'],
                'drug': first_row.get('drug', ''),
                'resistance_category': first_row.get('resistance_category', ''),
                'required_mutation_ids': ';'.join(required_mutations),
                'required_count': len(required_mutations),
                'combination_size': combination_size,
            }
        )

    return pd.DataFrame(records)


def build_single_catalog(mutation_catalog):
    single_rows = mutation_catalog[mutation_catalog['signature_kind'] == 'single'].copy()
    if single_rows.empty:
        return pd.DataFrame(
            columns=[
                'signature_id',
                'mutation_id',
                'protein_name',
                'segment',
                'aa_position',
                'alt_residue',
                'reference_accession',
                'drug',
                'resistance_category',
            ]
        )

    single_rows = single_rows.sort_values(['mutation_id'])
    deduped = single_rows.drop_duplicates(subset=['mutation_id']).copy()
    deduped['signature_id'] = deduped['signature_id'].where(deduped['signature_id'] != '', deduped['mutation_id'])
    return deduped[
        [
            'signature_id',
            'mutation_id',
            'protein_name',
            'segment',
            'aa_position',
            'alt_residue',
            'reference_accession',
            'drug',
            'resistance_category',
        ]
    ].reset_index(drop=True)


def build_long_without_combination_components(sequence_mutations):
    rows = sequence_mutations[sequence_mutations['combination_id'] == ''].copy()
    return rows.drop_duplicates().reset_index(drop=True)


def build_combination_signature_summary(sequence_mutations, mutation_catalog):
    combination_catalog = build_combination_catalog(mutation_catalog)
    combo_hits = sequence_mutations[sequence_mutations['combination_id'] != ''].copy()
    if combo_hits.empty or combination_catalog.empty:
        return pd.DataFrame(
            columns=[
                'primary_accession',
                'signature_id',
                'signature_kind',
                'protein_name',
                'segment',
                'reference_accession',
                'drug',
                'resistance_category',
                'required_mutation_ids',
                'mutations_present',
                'mutations_missing',
                'present_count',
                'required_count',
                'combination_status',
            ]
        )

    combo_hits = combo_hits.sort_values(['primary_accession', 'combination_id', 'mutation_id'])
    records = []
    combo_lookup = {
        row['combination_id']: row
        for _, row in combination_catalog.iterrows()
    }

    grouped = combo_hits.groupby(['primary_accession', 'combination_id'], sort=False)
    for (primary_accession, combination_id), group in grouped:
        combo_metadata = combo_lookup.get(combination_id)
        if combo_metadata is None:
            continue

        present_mutations = []
        for mutation_id in group['mutation_id'].tolist():
            if mutation_id and mutation_id not in present_mutations:
                present_mutations.append(mutation_id)

        required_mutations = [
            mutation_id for mutation_id in combo_metadata['required_mutation_ids'].split(';') if mutation_id
        ]
        missing_mutations = [mutation_id for mutation_id in required_mutations if mutation_id not in present_mutations]
        records.append(
            {
                'primary_accession': primary_accession,
                'signature_id': combo_metadata['signature_id'],
                'signature_kind': 'combination',
                'protein_name': combo_metadata['protein_name'],
                'segment': combo_metadata['segment'],
                'reference_accession': combo_metadata['reference_accession'],
                'drug': combo_metadata['drug'],
                'resistance_category': combo_metadata['resistance_category'],
                'required_mutation_ids': combo_metadata['required_mutation_ids'],
                'mutations_present': ';'.join(present_mutations),
                'mutations_missing': ';'.join(missing_mutations),
                'present_count': len(present_mutations),
                'required_count': int(combo_metadata['required_count']),
                'combination_status': 'complete' if not missing_mutations else 'partial',
            }
        )

    return pd.DataFrame(records)


def build_signature_summary(sequence_mutations, mutation_catalog):
    single_catalog = build_single_catalog(mutation_catalog)
    combination_summary = build_combination_signature_summary(sequence_mutations, mutation_catalog)

    single_hits = sequence_mutations[sequence_mutations['combination_id'] == ''].copy()
    if not single_hits.empty:
        single_hits = single_hits.drop_duplicates(subset=['primary_accession', 'mutation_id']).copy()
        single_hits = single_hits.merge(single_catalog, on='mutation_id', how='left', suffixes=('', '_catalog'))
        single_hits['signature_id'] = single_hits['signature_id'].where(single_hits['signature_id'] != '', single_hits['mutation_id'])
        single_hits['signature_kind'] = 'single'
        single_hits['required_mutation_ids'] = single_hits['mutation_id']
        single_hits['mutations_present'] = single_hits['mutation_id']
        single_hits['mutations_missing'] = ''
        single_hits['present_count'] = 1
        single_hits['required_count'] = 1
        single_hits['combination_status'] = 'complete'
        single_hits = single_hits[
            [
                'primary_accession',
                'signature_id',
                'signature_kind',
                'protein_name_catalog',
                'segment_catalog',
                'reference_accession',
                'drug',
                'resistance_category',
                'required_mutation_ids',
                'mutations_present',
                'mutations_missing',
                'present_count',
                'required_count',
                'combination_status',
            ]
        ].rename(
            columns={
                'protein_name_catalog': 'protein_name',
                'segment_catalog': 'segment',
            }
        )
    else:
        single_hits = pd.DataFrame(columns=combination_summary.columns)

    if combination_summary.empty:
        return single_hits.reset_index(drop=True)
    if single_hits.empty:
        return combination_summary.reset_index(drop=True)
    return pd.concat([single_hits, combination_summary], ignore_index=True).sort_values(
        ['primary_accession', 'signature_kind', 'signature_id']
    ).reset_index(drop=True)


def build_signature_summary_minimal(signature_summary):
    if signature_summary.empty:
        return pd.DataFrame(
            columns=['primary_accession', 'signature_id', 'signature_kind', 'mutations_present', 'mutations_missing']
        )
    return signature_summary[
        ['primary_accession', 'signature_id', 'signature_kind', 'mutations_present', 'mutations_missing']
    ].copy()


def build_completed_signatures_only(sequence_mutations, mutation_catalog):
    signature_summary = build_signature_summary(sequence_mutations, mutation_catalog)
    if signature_summary.empty:
        return pd.DataFrame(columns=['primary_accession', 'signature_id', 'signature_kind'])

    completed = signature_summary[signature_summary['combination_status'] == 'complete'].copy()
    if completed.empty:
        return pd.DataFrame(columns=['primary_accession', 'signature_id', 'signature_kind'])

    return completed[['primary_accession', 'signature_id', 'signature_kind']].drop_duplicates().reset_index(drop=True)


def build_sequence_relevant_mutation_summary(sequence_mutations):
    if sequence_mutations.empty:
        return pd.DataFrame(columns=['primary_accession', 'relevant_mutations_present', 'total_relevant_mutation_count'])

    deduped = sequence_mutations.drop_duplicates(subset=['primary_accession', 'mutation_id']).copy()
    deduped = deduped.sort_values(['primary_accession', 'mutation_id'])

    records = []
    for primary_accession, group in deduped.groupby('primary_accession', sort=False):
        mutation_ids = [mutation_id for mutation_id in group['mutation_id'].tolist() if mutation_id]
        records.append(
            {
                'primary_accession': primary_accession,
                'relevant_mutations_present': ';'.join(mutation_ids),
                'total_relevant_mutation_count': len(mutation_ids),
            }
        )

    return pd.DataFrame(records)


def build_layout_candidates(sequence_mutations, mutation_catalog):
    signature_summary = build_signature_summary(sequence_mutations, mutation_catalog)
    return {
        'baseline_sequence_mutations': sequence_mutations.drop_duplicates().reset_index(drop=True),
        'long_without_combination_components': build_long_without_combination_components(sequence_mutations),
        'combination_signature_summary': build_combination_signature_summary(sequence_mutations, mutation_catalog),
        'completed_signatures_only': build_completed_signatures_only(sequence_mutations, mutation_catalog),
        'sequence_relevant_mutation_summary': build_sequence_relevant_mutation_summary(sequence_mutations),
        'signature_summary': signature_summary,
        'signature_summary_minimal': build_signature_summary_minimal(signature_summary),
    }


def build_layout_report(candidates):
    report_rows = []
    for name, table in candidates.items():
        approx_bytes = int(table.memory_usage(deep=True).sum()) if not table.empty else 0
        report_rows.append(
            {
                'layout_name': name,
                'row_count': int(len(table)),
                'column_count': int(len(table.columns)),
                'populated_value_count': int(table.replace('', pd.NA).count().sum()) if not table.empty else 0,
                'approx_bytes': approx_bytes,
                'approx_size_human': _format_human_bytes(approx_bytes),
                'columns': ';'.join(table.columns.tolist()),
            }
        )
    return pd.DataFrame(report_rows).sort_values(['row_count', 'column_count', 'layout_name']).reset_index(drop=True)


def write_layout_outputs(candidates, report, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report.to_csv(output_dir / 'mutation_storage_layout_report.tsv', sep='\t', index=False)
    for name, table in candidates.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description='Explore alternative storage layouts for mutation tables.')
    parser.add_argument('--db', required=True, help='Path to SQLite database containing sequence_mutations and optionally mutation_catalog.')
    parser.add_argument('--mutation_catalog_tsv', help='Optional TSV to use when the DB does not contain a mutation_catalog table.')
    parser.add_argument('--output_dir', required=True, help='Directory where candidate TSVs and the summary report will be written.')
    args = parser.parse_args()

    sequence_mutations, mutation_catalog, catalog_source = load_mutation_tables(args.db, args.mutation_catalog_tsv)
    candidates = build_layout_candidates(sequence_mutations, mutation_catalog)
    report = build_layout_report(candidates)
    write_layout_outputs(candidates, report, args.output_dir)

    print(f'Catalog source: {catalog_source}')
    if catalog_source == 'inferred_from_sequence_mutations':
        print('Warning: mutation_catalog table was missing, so combination requirements were inferred from observed sequence_mutations rows.')

    print('Generated mutation storage layout candidates:')
    for _, row in report.iterrows():
        print(
            f"  - {row['layout_name']}: rows={row['row_count']:,}, columns={row['column_count']:,}, approx_size={row['approx_size_human']} ({row['approx_bytes']:,} bytes)"
        )


if __name__ == '__main__':
    main()