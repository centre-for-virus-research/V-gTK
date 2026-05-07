import csv
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEV_DIR = REPO_ROOT / 'dev'
SCRIPTS_DIR = REPO_ROOT / 'scripts'

if str(DEV_DIR) not in sys.path:
	sys.path.insert(0, str(DEV_DIR))
if str(SCRIPTS_DIR) not in sys.path:
	sys.path.insert(0, str(SCRIPTS_DIR))

from NormalizeHcvMutationCatalog import HcvMutationCatalogNormalizer  # type: ignore[reportMissingImports]
from probe_hcv_catalog_columns import probe_catalog_enrichment  # type: ignore[reportMissingImports]


def write_table(path: Path, rows, fieldnames, delimiter=','):
	path.parent.mkdir(parents=True, exist_ok=True)
	with path.open('w', encoding='utf-8', newline='') as handle:
		writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
		writer.writeheader()
		for row in rows:
			writer.writerow(row)


def result_for(results, column_name):
	match = [result for result in results if result['column'] == column_name]
	assert len(match) == 1
	return match[0]


def test_probe_catalog_enrichment_reports_overlaps_and_derived_column_mismatches(tmp_path: Path):
	table_dir = tmp_path / 'Tables'

	write_table(
		table_dir / 'gene_info.tsv',
		[
			{'description': 'Non-structural protein 3', 'display_name': 'NS3', 'name': 'NS3', 'parent_name': 'whole_genome'},
			{'description': 'Whole genome', 'display_name': 'Whole genome', 'name': 'whole_genome', 'parent_name': 'NULL'},
		],
		['description', 'display_name', 'name', 'parent_name'],
		delimiter='\t',
	)

	write_table(
		table_dir / 'variation.csv',
		[
			{
				'description': '',
				'display_name': '',
				'feature_name': 'NS3',
				'name': 'phdr_ras:NS3:107I',
				'ref_end': '10',
				'ref_seq_name': 'REF_MASTER',
				'ref_start': '8',
				'type': 'aminoAcidSimplePolymorphism',
				'phdr_ras_id': 'NS3:107I',
			},
			{
				'description': '',
				'display_name': '',
				'feature_name': 'NS3',
				'name': 'phdr_ras:NS3:155K',
				'ref_end': '20',
				'ref_seq_name': 'REF_MASTER',
				'ref_start': '18',
				'type': 'aminoAcidSimplePolymorphism',
				'phdr_ras_id': 'NS3:155K',
			},
		],
		['description', 'display_name', 'feature_name', 'name', 'ref_end', 'ref_seq_name', 'ref_start', 'type', 'phdr_ras_id'],
	)

	write_table(
		table_dir / 'variation_metatag.csv',
		[],
		['feature_name', 'metatag_name', 'metatag_value', 'ref_seq_name', 'variation_name'],
	)

	write_table(
		table_dir / 'phdr_alignment_ras.csv',
		[
			{'id': 'NS3:107I:AL_1a', 'display_structure': 'V107I', 'phdr_ras_id': 'NS3:107I', 'alignment_name': 'AL_1a'},
			{'id': 'NS3:155K:AL_1a', 'display_structure': 'R155K', 'phdr_ras_id': 'NS3:155K', 'alignment_name': 'AL_1a'},
		],
		['id', 'display_structure', 'phdr_ras_id', 'alignment_name'],
	)

	write_table(
		table_dir / 'phdr_alignment_ras_drug.csv',
		[
			{
				'id': 'row_1',
				'resistance_category': 'category_I',
				'display_resistance_category': 'I',
				'numeric_resistance_category': '1',
				'any_in_vitro_evidence': '1',
				'in_vitro_max_ec50_midpoint': '0.8',
				'any_in_vivo_evidence': '1',
				'in_vivo_baseline': '1',
				'in_vivo_treatment_emergent': '1',
				'phdr_alignment_ras_id': 'NS3:107I:AL_1a',
				'phdr_drug_id': 'drug_a',
			},
			{
				'id': 'row_2',
				'resistance_category': 'category_II',
				'display_resistance_category': 'II',
				'numeric_resistance_category': '2',
				'any_in_vitro_evidence': '1',
				'in_vitro_max_ec50_midpoint': '2.4',
				'any_in_vivo_evidence': '0',
				'in_vivo_baseline': '',
				'in_vivo_treatment_emergent': '',
				'phdr_alignment_ras_id': 'NS3:155K:AL_1a',
				'phdr_drug_id': 'drug_b',
			},
		],
		[
			'id',
			'resistance_category',
			'display_resistance_category',
			'numeric_resistance_category',
			'any_in_vitro_evidence',
			'in_vitro_max_ec50_midpoint',
			'any_in_vivo_evidence',
			'in_vivo_baseline',
			'in_vivo_treatment_emergent',
			'phdr_alignment_ras_id',
			'phdr_drug_id',
		],
	)

	write_table(
		table_dir / 'phdr_drug.csv',
		[
			{'id': 'drug_a', 'producer': 'Maker A', 'drug_category': 'Cat A'},
			{'id': 'drug_b', 'producer': 'Maker B', 'drug_category': 'Cat B'},
		],
		['id', 'producer', 'drug_category'],
	)

	write_table(
		table_dir / 'phdr_resistance_finding.csv',
		[
			{'phdr_alignment_ras_drug_id': 'row_1', 'phdr_publication_id': 'PUB_1'},
			{'phdr_alignment_ras_drug_id': 'row_1', 'phdr_publication_id': 'PUB_2'},
			{'phdr_alignment_ras_drug_id': 'row_2', 'phdr_publication_id': 'PUB_3'},
		],
		['phdr_alignment_ras_drug_id', 'phdr_publication_id'],
	)

	write_table(
		table_dir / 'phdr_publication.csv',
		[
			{'id': 'PUB_1', 'url': 'https://doi.org/1'},
			{'id': 'PUB_2', 'url': 'https://doi.org/2'},
			{'id': 'PUB_3', 'url': 'https://doi.org/3'},
		],
		['id', 'url'],
	)

	normalizer = HcvMutationCatalogNormalizer(
		variation_path=table_dir / 'variation.csv',
		variation_metatag_path=table_dir / 'variation_metatag.csv',
		phdr_alignment_ras_path=table_dir / 'phdr_alignment_ras.csv',
		phdr_alignment_ras_drug_path=table_dir / 'phdr_alignment_ras_drug.csv',
		gene_info_path=table_dir / 'gene_info.tsv',
		output_path=table_dir / 'ignored.tsv',
	)
	base_catalog_path = normalizer.normalize()['generalized_mutation_catalog']
	base_catalog = pd.read_csv(base_catalog_path, sep='\t')
	base_catalog['drug_producer'] = base_catalog['drug'].map({'drug_a': 'Maker A', 'drug_b': 'Maker B'})
	base_catalog['drug_category'] = base_catalog['drug'].map({'drug_a': 'Cat A', 'drug_b': 'Cat B'})
	base_catalog['pubmed_id'] = base_catalog['id'].map({'row_1': 'PUB_1;PUB_2', 'row_2': 'PUB_3'})
	base_catalog['DOI'] = base_catalog['id'].map({'row_1': 'https://doi.org/1;https://doi.org/2', 'row_2': 'https://doi.org/wrong'})
	base_catalog.to_csv(table_dir / 'generalized_mutation_catalog_with_extra_info.tsv', sep='\t', index=False)

	report = probe_catalog_enrichment(table_dir)

	assert report['row_count_matches_normalized_catalog'] is True
	assert any(overlap['table'] == 'phdr_drug.csv' and overlap['shared_columns'] == ['drug_category', 'id'] for overlap in report['column_overlaps'])
	assert any(overlap['table'] == 'phdr_resistance_finding.csv' and overlap['shared_columns'] == [] for overlap in report['column_overlaps'])
	assert any(
		candidate['table'] == 'phdr_resistance_finding.csv'
		and {'catalog_column': 'id', 'source_column': 'phdr_alignment_ras_drug_id', 'kind': 'direct value join'} in candidate['matches']
		for candidate in report['candidate_matches']
	)

	assert result_for(report['normalized_catalog_columns'], 'protein_name')['status'] == 'ok'
	assert result_for(report['normalized_catalog_columns'], 'mutation_id')['status'] == 'ok'
	assert result_for(report['derived_columns'], 'drug_producer')['status'] == 'ok'
	assert result_for(report['derived_columns'], 'drug_category')['status'] == 'ok'
	assert result_for(report['derived_columns'], 'pubmed_id')['status'] == 'ok'

	doi_result = result_for(report['derived_columns'], 'DOI')
	assert doi_result['status'] == 'mismatch'
	assert doi_result['mismatched_rows'] == 1
	assert doi_result['sample_mismatches'][0]['row_index'] == 1
	assert doi_result['sample_mismatches'][0]['expected'] == 'https://doi.org/3'