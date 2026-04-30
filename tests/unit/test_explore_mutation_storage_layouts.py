import sqlite3
import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / 'scripts' / 'ExploreMutationStorageLayouts.py'


def _seed_mutation_storage_db(db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            '''
            CREATE TABLE mutation_catalog (
                mutation_id TEXT,
                protein_name TEXT,
                segment TEXT,
                aa_position TEXT,
                alt_residue TEXT,
                reference_accession TEXT,
                mutation_type TEXT,
                signature_id TEXT,
                signature_kind TEXT,
                combination_id TEXT,
                combination_size TEXT,
                resistance_category TEXT,
                drug TEXT
            )
            '''
        )
        cur.executemany(
            'INSERT INTO mutation_catalog VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
            [
                ('NS3:107I', 'NS3', '1', '107', 'I', 'NC_004102', 'snp', 'NS3:107I', 'single', '', '', 'I', 'drugA'),
                ('NS5A:30Y', 'NS5A', '1', '30', 'Y', 'NC_004102', 'snp', 'NS5A:30Y', 'single', '', '', 'II', 'drugB'),
                ('NS3:107I', 'NS3', '1', '107', 'I', 'NC_004102', 'snp', 'comboA', 'combination', 'comboA', '2', 'III', 'drugA'),
                ('NS3:155K', 'NS3', '1', '155', 'K', 'NC_004102', 'snp', 'comboA', 'combination', 'comboA', '2', 'III', 'drugA'),
            ],
        )
        cur.execute(
            '''
            CREATE TABLE sequence_mutations (
                primary_accession TEXT,
                mutation_id TEXT,
                protein_name TEXT,
                segment TEXT,
                aa_position INTEGER,
                alt_residue TEXT,
                combination_id TEXT
            )
            '''
        )
        cur.executemany(
            'INSERT INTO sequence_mutations VALUES (?, ?, ?, ?, ?, ?, ?)',
            [
                ('SEQ1', 'NS3:107I', 'NS3', '1', 107, 'I', ''),
                ('SEQ1', 'NS3:107I', 'NS3', '1', 107, 'I', 'comboA'),
                ('SEQ1', 'NS3:155K', 'NS3', '1', 155, 'K', 'comboA'),
                ('SEQ1', 'NS5A:30Y', 'NS5A', '1', 30, 'Y', ''),
                ('SEQ2', 'NS3:107I', 'NS3', '1', 107, 'I', ''),
                ('SEQ2', 'NS3:107I', 'NS3', '1', 107, 'I', 'comboA'),
                ('SEQ2', 'NS5A:30Y', 'NS5A', '1', 30, 'Y', ''),
            ],
        )
        conn.commit()
    finally:
        conn.close()


def _seed_sequence_mutations_only_db(db_path: Path):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            '''
            CREATE TABLE sequence_mutations (
                primary_accession TEXT,
                mutation_id TEXT,
                protein_name TEXT,
                segment TEXT,
                aa_position INTEGER,
                alt_residue TEXT,
                combination_id TEXT
            )
            '''
        )
        cur.executemany(
            'INSERT INTO sequence_mutations VALUES (?, ?, ?, ?, ?, ?, ?)',
            [
                ('SEQ1', 'NS3:107I', 'NS3', '1', 107, 'I', ''),
                ('SEQ1', 'NS3:107I', 'NS3', '1', 107, 'I', 'comboA'),
                ('SEQ1', 'NS3:155K', 'NS3', '1', 155, 'K', 'comboA'),
                ('SEQ2', 'NS3:107I', 'NS3', '1', 107, 'I', 'comboA'),
            ],
        )
        conn.commit()
    finally:
        conn.close()


def test_explore_mutation_storage_layouts_generates_compressed_candidates(tmp_path: Path):
    db_path = tmp_path / 'mutation_storage.db'
    output_dir = tmp_path / 'layout_outputs'
    _seed_mutation_storage_db(db_path)

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), '--db', str(db_path), '--output_dir', str(output_dir)],
        check=True,
        capture_output=True,
        text=True,
    )

    assert 'Generated mutation storage layout candidates:' in result.stdout

    report = pd.read_csv(output_dir / 'mutation_storage_layout_report.tsv', sep='\t', keep_default_na=False)
    assert 'approx_size_human' in report.columns
    baseline_rows = int(report.loc[report['layout_name'] == 'baseline_sequence_mutations', 'row_count'].iloc[0])
    combo_summary_rows = int(report.loc[report['layout_name'] == 'combination_signature_summary', 'row_count'].iloc[0])
    completed_only_rows = int(report.loc[report['layout_name'] == 'completed_signatures_only', 'row_count'].iloc[0])
    per_sequence_rows = int(report.loc[report['layout_name'] == 'sequence_relevant_mutation_summary', 'row_count'].iloc[0])
    signature_summary_rows = int(report.loc[report['layout_name'] == 'signature_summary', 'row_count'].iloc[0])

    assert baseline_rows == 7
    assert combo_summary_rows == 2
    assert completed_only_rows == 5
    assert per_sequence_rows == 2
    assert signature_summary_rows == 6
    assert result.stdout.count('approx_size=') == len(report)
    assert ' bytes)' in result.stdout

    combo_summary = pd.read_csv(output_dir / 'combination_signature_summary.tsv', sep='\t', keep_default_na=False)
    combo_summary = combo_summary.sort_values(['primary_accession']).reset_index(drop=True)

    assert combo_summary['primary_accession'].tolist() == ['SEQ1', 'SEQ2']
    assert combo_summary['signature_id'].tolist() == ['comboA', 'comboA']
    assert combo_summary['mutations_present'].tolist() == ['NS3:107I;NS3:155K', 'NS3:107I']
    assert combo_summary['mutations_missing'].tolist() == ['', 'NS3:155K']
    assert combo_summary['combination_status'].tolist() == ['complete', 'partial']

    signature_summary = pd.read_csv(output_dir / 'signature_summary.tsv', sep='\t', keep_default_na=False)
    signature_summary = signature_summary.sort_values(['primary_accession', 'signature_kind', 'signature_id']).reset_index(drop=True)
    assert signature_summary['signature_id'].tolist() == ['comboA', 'NS3:107I', 'NS5A:30Y', 'comboA', 'NS3:107I', 'NS5A:30Y']

    completed_only = pd.read_csv(output_dir / 'completed_signatures_only.tsv', sep='\t', keep_default_na=False)
    completed_only = completed_only.sort_values(['primary_accession', 'signature_kind', 'signature_id']).reset_index(drop=True)
    assert completed_only.to_dict('records') == [
        {'primary_accession': 'SEQ1', 'signature_id': 'comboA', 'signature_kind': 'combination'},
        {'primary_accession': 'SEQ1', 'signature_id': 'NS3:107I', 'signature_kind': 'single'},
        {'primary_accession': 'SEQ1', 'signature_id': 'NS5A:30Y', 'signature_kind': 'single'},
        {'primary_accession': 'SEQ2', 'signature_id': 'NS3:107I', 'signature_kind': 'single'},
        {'primary_accession': 'SEQ2', 'signature_id': 'NS5A:30Y', 'signature_kind': 'single'},
    ]

    per_sequence = pd.read_csv(output_dir / 'sequence_relevant_mutation_summary.tsv', sep='\t', keep_default_na=False)
    per_sequence = per_sequence.sort_values(['primary_accession']).reset_index(drop=True)
    assert per_sequence.to_dict('records') == [
        {
            'primary_accession': 'SEQ1',
            'relevant_mutations_present': 'NS3:107I;NS3:155K;NS5A:30Y',
            'total_relevant_mutation_count': 3,
        },
        {
            'primary_accession': 'SEQ2',
            'relevant_mutations_present': 'NS3:107I;NS5A:30Y',
            'total_relevant_mutation_count': 2,
        },
    ]

    minimal = pd.read_csv(output_dir / 'signature_summary_minimal.tsv', sep='\t', keep_default_na=False)
    assert list(minimal.columns) == ['primary_accession', 'signature_id', 'signature_kind', 'mutations_present', 'mutations_missing']


def test_explore_mutation_storage_layouts_falls_back_when_mutation_catalog_missing(tmp_path: Path):
    db_path = tmp_path / 'sequence_mutations_only.db'
    output_dir = tmp_path / 'layout_outputs'
    _seed_sequence_mutations_only_db(db_path)

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), '--db', str(db_path), '--output_dir', str(output_dir)],
        check=True,
        capture_output=True,
        text=True,
    )

    assert 'Catalog source: inferred_from_sequence_mutations' in result.stdout
    assert 'Warning: mutation_catalog table was missing' in result.stdout

    report = pd.read_csv(output_dir / 'mutation_storage_layout_report.tsv', sep='\t', keep_default_na=False)
    combo_summary_rows = int(report.loc[report['layout_name'] == 'combination_signature_summary', 'row_count'].iloc[0])
    completed_only_rows = int(report.loc[report['layout_name'] == 'completed_signatures_only', 'row_count'].iloc[0])
    per_sequence_rows = int(report.loc[report['layout_name'] == 'sequence_relevant_mutation_summary', 'row_count'].iloc[0])
    assert combo_summary_rows == 2
    assert completed_only_rows == 2
    assert per_sequence_rows == 2
    assert 'approx_size=' in result.stdout

    combo_summary = pd.read_csv(output_dir / 'combination_signature_summary.tsv', sep='\t', keep_default_na=False)
    combo_summary = combo_summary.sort_values(['primary_accession']).reset_index(drop=True)
    assert combo_summary['mutations_present'].tolist() == ['NS3:107I;NS3:155K', 'NS3:107I']
    assert combo_summary['mutations_missing'].tolist() == ['', 'NS3:155K']

    completed_only = pd.read_csv(output_dir / 'completed_signatures_only.tsv', sep='\t', keep_default_na=False)
    completed_only = completed_only.sort_values(['primary_accession', 'signature_kind', 'signature_id']).reset_index(drop=True)
    assert completed_only.to_dict('records') == [
        {'primary_accession': 'SEQ1', 'signature_id': 'comboA', 'signature_kind': 'combination'},
        {'primary_accession': 'SEQ1', 'signature_id': 'NS3:107I', 'signature_kind': 'single'},
    ]

    per_sequence = pd.read_csv(output_dir / 'sequence_relevant_mutation_summary.tsv', sep='\t', keep_default_na=False)
    per_sequence = per_sequence.sort_values(['primary_accession']).reset_index(drop=True)
    assert per_sequence.to_dict('records') == [
        {
            'primary_accession': 'SEQ1',
            'relevant_mutations_present': 'NS3:107I;NS3:155K',
            'total_relevant_mutation_count': 2,
        },
        {
            'primary_accession': 'SEQ2',
            'relevant_mutations_present': 'NS3:107I',
            'total_relevant_mutation_count': 1,
        },
    ]