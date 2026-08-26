import pytest
import sqlite3
import pandas as pd
import os
import sys

# Assume scripts/ is in PYTHONPATH when running pytest, but adding it anyway
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../scripts')))

import AnnotateMutations


HCV_CATALOG_HEADER = [
    'protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession',
    'mutation_id', 'mutation_type', 'signature_id', 'signature_kind', 'combination_id',
    'combination_size', 'phenotype', 'resistance_category', 'drug'
]

HCV_DB_COLUMNS = [
    'mutation_id', 'protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession',
    'mutation_type', 'signature_id', 'signature_kind', 'combination_id', 'combination_size',
    'phenotype', 'resistance_category', 'drug', 'drug_category', 'drug_producer', 'pubmed_id',
    'DOI', 'any_in_vitro_evidence', 'in_vitro_max_ec50_midpoint', 'any_in_vivo_evidence',
    'in_vivo_baseline', 'in_vivo_treatment_emergent'
]


def _mutation_ids_by_accession(summary_df):
    result = {}
    for _, row in summary_df.iterrows():
        mutation_ids = [value for value in str(row['relevant_mutations_present'] or '').split(';') if value]
        result[row['primary_accession']] = mutation_ids
    return result


def _completed_signatures_by_accession(completed_df):
    result = {}
    for primary_accession, group in completed_df.groupby('primary_accession', sort=False):
        result[primary_accession] = set(group['signature_id'].tolist())
    return result

def test_translate_codon():
    assert AnnotateMutations.translate_codon('ATG') == 'M'
    assert AnnotateMutations.translate_codon('A-G') == 'X'
    assert AnnotateMutations.translate_codon('TAA') == '*'

def test_annotate_mutations_db_creation(tmp_path):
    db_path = tmp_path / "test.db"
    catalog_path = tmp_path / "catalog.tsv"
    
    # 1. Create a dummy catalog
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'I', 'REF1', 'NS3:2I', 'snp', 'sig1', 'single', 'C1', '2', '', 'I', 'DrugA'],
        ['NS3', '1', '3', 'K', 'REF1', 'NS3:3K', 'snp', 'sig2', 'combination', 'C1', '2', '', 'I', 'DrugA'],
        ['NS5A', '', '1', 'Y', 'REF1', 'NS5A:1Y', 'snp', 'sig3', 'single', '', '', '', '', '']
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    # 2. Setup SQLite DB
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    # Using 1-based indexing for cds_start as in the typical outputs
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF1', 'NS3', '1', 4, 12),  # Nucleotides 4-12 are NS3
        ('REF1', 'NS5A', '', 13, 15), # Nucleotides 13-15 are NS5A
        ('seq1', 'polyprotein', '', 4, 15),
        ('seq2', 'polyprotein', '', 4, 15)
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    # seq1 padded alignment:
    # 1-3: ABC
    # 4-6: ATA (AA 1: I), 7-9: ATC (AA 2: I), 10-12: AAA (AA 3: K)  => NS3 has 2I, 3K
    # 13-15: TAC (AA 1: Y) => NS5A has 1Y
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF1', 'REF1', 'REF1', 'ABCATAATCAAATAC'),
        ('seq1', 'seq1', 'REF1', 'ABCATAATCAAATAC'),
        ('seq2', 'seq2', 'REF1', 'ABCATAGGGAAA')
    ])
    
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF1', 'master'))
    conn.commit()
    conn.close()
    
    # 3. Call the script logic
    # Mocking sys.argv
    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', 'HCV']
    AnnotateMutations.main()
    
    # 4. Verify outcomes
    conn = sqlite3.connect(str(db_path))
    catalog_table = pd.read_sql_query("SELECT * FROM mutation_catalog ORDER BY mutation_id, drug", conn)
    mutation_summary = pd.read_sql_query(
        "SELECT * FROM sequence_relevant_mutation_summary ORDER BY primary_accession",
        conn,
    )
    completed = pd.read_sql_query(
        "SELECT * FROM completed_signatures_only ORDER BY primary_accession, signature_id",
        conn,
    )
    
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    assert mutation_ids['seq1'] == ['NS3:2I', 'NS3:3K', 'NS5A:1Y']
    assert mutation_ids['seq2'] == ['NS3:3K']
    assert mutation_summary.set_index('primary_accession').loc['seq1', 'total_relevant_mutation_count'] == 3
    assert mutation_summary.set_index('primary_accession').loc['seq2', 'total_relevant_mutation_count'] == 1

    completed_by_accession = _completed_signatures_by_accession(completed)
    assert completed_by_accession['seq1'] == {'sig2', 'sig3'}
    assert completed_by_accession['seq2'] == {'sig2'}

    assert set(catalog_table.columns) >= {
        'mutation_id', 'protein_name', 'segment', 'aa_position', 'alt_residue',
        'reference_accession', 'mutation_type', 'signature_id', 'signature_kind',
        'combination_id', 'combination_size', 'phenotype', 'resistance_category', 'drug'
    }
    assert sorted(catalog_table['mutation_id'].tolist()) == ['NS3:2I', 'NS3:3K', 'NS5A:1Y']
    
    conn.close()

def test_alignment_gap_proportions_in_db(tmp_path):
    """
    Check included sequences in the DB contain the right amount of gaps 
    and aren't too full of them.
    In a real dataset, we expect sequences not to be primarily gaps (e.g. < 50%).
    """
    db_path = tmp_path / "gap_test.db"
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, alignment TEXT)")
    
    # seq1: 30 bases, 3 gaps = 10% gaps (OK)
    seq1_aln = "ATG" * 9 + "---"
    # seq2: 30 bases, 27 gaps = 90% gaps (Too many)
    seq2_aln = "ATG" + "-" * 27
    
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?)", [
        ('seq1', seq1_aln),
        ('seq2', seq2_aln)
    ])
    conn.commit()
    
    # Fetch and check
    df = pd.read_sql_query("SELECT sequence_id, alignment FROM sequence_alignment", conn)
    gap_fractions = df['alignment'].apply(lambda x: x.count('-') / len(x) if len(x) > 0 else 1.0)
    
    assert gap_fractions.iloc[0] == 0.1
    assert gap_fractions.iloc[1] == 0.9
    conn.close()

def test_annotation_with_reference_indels(tmp_path):
    """
    Check that the padded alignment with indels has the expected reference amino acid 
    in the positions calculated by the mutation annotation.
    """
    db_path = tmp_path / "indel_test.db"
    catalog_path = tmp_path / "indel_catalog.tsv"
    
    # 1. Dummy catalog
    catalog_data = [
        HCV_CATALOG_HEADER,
        # Reference AA is M (ATG), we look for a mutation to V.
        ['NS3', '1', '2', 'V', 'REF_MASTER', 'NS3:2V', 'snp', 'sig1', 'single', '', '', '', '', ''],
        # Mutation to Q at position 3 (CAA).
        ['NS3', '1', '3', 'Q', 'REF_MASTER', 'NS3:3Q', 'snp', 'sig2', 'single', '', '', '', '', '']
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    # Features in the DB are query-sequence coordinates, not aligned-string offsets.
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_MASTER', 'NS3', '1', 4, 12),
        ('mutated_seq', 'polyprotein', '1', 4, 12)
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    
    # Padded Alignment layout: (1-based indices relative to alignment string)
    # 123 | 456 | 789 | 10,11,12 | 13,14,15
    # pos | AA1 | AA2 |   AA3    |    AA4
    # Pad | ATG | ATG |   ---    |    CAT
    ref_aln = "GGG" + "ATG" + "ATG" + "---" + "CAT"
    
    # Mutated has V at AA2 (GTA) and Q at AA4 (CAA). And still has deletion at AA3.
    mut_aln = "GGG" + "ATG" + "GTA" + "---" + "CAA"
    
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_MASTER', 'REF_MASTER', 'REF_MASTER', ref_aln),
        ('mutated_seq', 'mutated_seq', 'REF_MASTER', mut_aln)
    ])
    
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_MASTER', 'master'))
    conn.commit()
    conn.close()
    
    # 2. Run annotation
    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()
    
    # 3. Verify
    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    completed = pd.read_sql_query("SELECT * FROM completed_signatures_only", conn)
    
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    assert 'REF_MASTER' not in mutation_ids, "Reference should not trigger any mutation hits"
    assert mutation_ids['mutated_seq'] == ['NS3:2V', 'NS3:3Q']

    completed_by_accession = _completed_signatures_by_accession(completed)
    assert completed_by_accession['mutated_seq'] == {'sig1', 'sig2'}
    
    conn.close()


def test_extract_feature_codon_returns_none_when_query_coordinates_do_not_map():
    coord_map = AnnotateMutations.build_alignment_coordinate_map('GGG---ATA')

    codon = AnnotateMutations.extract_feature_codon('GGG---ATA', coord_map, cds_start=8, aa_pos=1)

    assert codon is None


def test_annotate_mutations_falls_back_to_db_gff_coordinate_mapping(tmp_path):
    db_path = tmp_path / "fallback_test.db"
    catalog_path = tmp_path / "fallback_catalog.tsv"

    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'I', 'REF_MASTER_NC_004102', 'NS3:2I', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS5A', '1', '1', 'Y', 'REF_MASTER_NC_004102', 'NS5A:1Y', 'snp', 'sig2', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?)", [
        ('REF_MASTER', 'polyprotein', 4, 12),
        ('seq1', 'polyprotein', 4, 15),
        ('seq2', 'polyprotein', 4, 15),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_MASTER', 'REF_MASTER', 'REF_MASTER', 'GGGCCGCCACAC'),
        ('seq1', 'seq1', 'REF_MASTER', 'GGGCCCATATAC'),
        ('seq2', 'seq2', 'REF_MASTER', 'GGGCCCAAATAC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_MASTER', 'master'))
    cursor.execute("CREATE TABLE gff_features (reference_accession TEXT, feature_type TEXT, product TEXT, start INTEGER, end INTEGER)")
    cursor.executemany("INSERT INTO gff_features VALUES (?, ?, ?, ?, ?)", [
        ('REF_MASTER', 'gene', 'NS3', 4, 9),
        ('REF_MASTER', 'mat_peptide', 'NS5A', 10, 12),
    ])
    cursor.execute("CREATE TABLE genes (description TEXT, display_name TEXT, name TEXT, parent_name TEXT)")
    cursor.executemany("INSERT INTO genes VALUES (?, ?, ?, ?)", [
        ('Non-structural protein 3', 'NS3', 'NS3', 'whole_genome'),
        ('Non-structural protein 5A', 'NS5A', 'NS5A', 'whole_genome'),
        ('Whole genome', 'Whole genome', 'whole_genome', 'NULL'),
    ])
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query(
        "SELECT * FROM sequence_relevant_mutation_summary ORDER BY primary_accession",
        conn,
    )
    completed = pd.read_sql_query(
        "SELECT * FROM completed_signatures_only ORDER BY primary_accession, signature_id",
        conn,
    )
    conn.close()

    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    assert mutation_ids == {
        'seq1': ['NS3:2I', 'NS5A:1Y'],
        'seq2': ['NS5A:1Y'],
    }
    completed_by_accession = _completed_signatures_by_accession(completed)
    assert completed_by_accession == {
        'seq1': {'sig1', 'sig2'},
        'seq2': {'sig2'},
    }


def test_annotate_mutations_can_build_reference_maps_from_features_table(tmp_path):
    db_path = tmp_path / "feature_fallback.db"

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute(
        "CREATE TABLE features (accession TEXT, master_ref_accession TEXT, reference_accession TEXT, cds_start INTEGER, cds_end INTEGER, product TEXT)"
    )
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?, ?)", [
        ('NC_004102', 'NC_004102', 'NC_004102', 4, 9, 'NS3'),
        ('NC_004102', 'NC_004102', 'NC_004102', 10, 12, 'NS5A'),
        ('EU781827', 'NC_004102', 'EU781827', 4, 9, 'NS3'),
        ('EU781827', 'NC_004102', 'EU781827', 10, 12, 'NS5A'),
        ('seq1', 'NC_004102', 'EU781827', 4, 12, 'polyprotein'),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('EU781827', 'EU781827', 'EU781827', 'GGGCCGCCACAC'),
        ('seq1', 'seq1', 'EU781827', 'GGGCCCATATAC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?)", [
        ('NC_004102', 'master'),
        ('EU781827', 'reference'),
    ])
    cursor.execute("CREATE TABLE genes (description TEXT, display_name TEXT, name TEXT, parent_name TEXT)")
    cursor.executemany("INSERT INTO genes VALUES (?, ?, ?, ?)", [
        ('Non-structural protein 3', 'NS3', 'NS3', 'whole_genome'),
        ('Non-structural protein 5A', 'NS5A', 'NS5A', 'whole_genome'),
        ('Whole genome', 'Whole genome', 'whole_genome', 'NULL'),
    ])
    conn.commit()

    alias_lookup = AnnotateMutations.load_gene_alias_lookup(conn)
    db_gff_maps = AnnotateMutations.load_db_gff_feature_maps(conn, alias_lookup)
    seq_aln = pd.read_sql_query("SELECT * FROM sequence_alignment", conn)
    meta_data = pd.read_sql_query("SELECT * FROM meta_data", conn)
    conn.close()

    catalog = pd.DataFrame([
        {
            'protein_name': 'NS3',
            'segment': '1',
            'aa_position': '2',
            'alt_residue': 'I',
            'reference_accession': 'REF_MASTER_NC_004102',
            'mutation_id': 'NS3:2I',
            'combination_id': '',
        },
        {
            'protein_name': 'NS5A',
            'segment': '1',
            'aa_position': '1',
            'alt_residue': 'Y',
            'reference_accession': 'REF_MASTER_NC_004102',
            'mutation_id': 'NS5A:1Y',
            'combination_id': '',
        },
    ])
    catalog, invalid_positions = AnnotateMutations.prepare_catalog(catalog, alias_lookup)

    assert invalid_positions == 0
    assert set(db_gff_maps) == {'NC_004102', 'EU781827'}

    mutations_found, diagnostics, resolved_maps = AnnotateMutations.annotate_from_reference_coordinates(
        catalog,
        seq_aln[['sequence_id', 'primary_accession', 'alignment', 'alignment_name']].copy(),
        meta_data,
        alias_lookup,
        db_gff_maps,
        False,
    )

    assert sorted(m['mutation_id'] for m in mutations_found) == ['NS3:2I', 'NS5A:1Y']
    assert all(m['primary_accession'] == 'seq1' for m in mutations_found)
    assert resolved_maps['EU781827']['resolved_accession'] == 'EU781827'
    assert diagnostics['mutation_hits'] == 2


def test_reference_coordinate_mapping_keeps_catalog_reference_hints_when_master_exists():
    catalog = pd.DataFrame(
        [
            {
                'protein_name': 'NS3',
                'segment': '1',
                'aa_position': '1',
                'alt_residue': 'I',
                'reference_accession': 'CATREF',
                'mutation_id': 'NS3:1I',
                'mutation_type': 'snp',
                'signature_id': 'sig1',
                'signature_kind': 'single',
                'combination_id': '',
                'combination_size': '',
                'phenotype': '',
            }
        ]
    )
    catalog, invalid_positions = AnnotateMutations.prepare_catalog(catalog, {})
    assert invalid_positions == 0

    seq_aln = pd.DataFrame(
        [
            {'sequence_id': 'CATREF', 'primary_accession': 'CATREF', 'alignment_name': 'ALIGNREF', 'alignment': 'ATG'},
            {'sequence_id': 'seq1', 'primary_accession': 'seq1', 'alignment_name': 'ALIGNREF', 'alignment': 'ATA'},
        ]
    )
    meta_data = pd.DataFrame(
        [
            {'primary_accession': 'MASTERREF', 'accession_type': 'master'},
        ]
    )
    db_gff_maps = {
        'CATREF': {
            'NS3': {
                'product': 'NS3',
                'raw_product': 'NS3',
                'cds_start': 1,
                'cds_end': 3,
                'reference_accession': 'CATREF',
                'feature_type': 'gene',
                'source': 'db_gff',
            }
        }
    }

    mutations_found, diagnostics, resolved_maps = AnnotateMutations.annotate_from_reference_coordinates(
        catalog,
        seq_aln,
        meta_data,
        {},
        db_gff_maps,
        False,
    )

    assert resolved_maps['CATREF']['resolved_accession'] == 'CATREF'
    assert [mutation['mutation_id'] for mutation in mutations_found] == ['NS3:1I']
    assert diagnostics['mutation_hits'] == 1


def test_annotate_mutations_uses_descriptive_hcv_feature_products_without_gene_aliases(tmp_path):
    db_path = tmp_path / "descriptive_hcv_features.db"
    catalog_path = tmp_path / "descriptive_hcv_catalog.tsv"

    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'I', 'REF_MASTER_NC_004102', 'NS3:2I', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS5A', '1', '1', 'Y', 'REF_MASTER_NC_004102', 'NS5A:1Y', 'snp', 'sig2', 'single', '', '', '', '', ''],
        ['NS5B', '1', '1', 'F', 'REF_MASTER_NC_004102', 'NS5B:1F', 'snp', 'sig3', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_MASTER', 'protease/helicase protein NS3', '1', 1, 6),
        ('REF_MASTER', 'nonstructural protein NS5A', '1', 7, 9),
        ('REF_MASTER', 'RNA-dependent RNA polymerase NS5B', '1', 10, 12),
        ('seq1', 'polyprotein', '1', 1, 12),
        ('seq2', 'polyprotein', '1', 1, 12),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_MASTER', 'REF_MASTER', 'REF_MASTER', 'ATAATCTACTTT'),
        ('seq1', 'seq1', 'REF_MASTER', 'ATAATCTACTTT'),
        ('seq2', 'seq2', 'REF_MASTER', 'ATAATCGGGTTC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_MASTER', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    catalog_table = pd.read_sql_query("SELECT * FROM mutation_catalog ORDER BY mutation_id", conn)
    mutation_summary = pd.read_sql_query(
        "SELECT * FROM sequence_relevant_mutation_summary ORDER BY primary_accession",
        conn,
    )
    completed = pd.read_sql_query(
        "SELECT * FROM completed_signatures_only ORDER BY primary_accession, signature_id",
        conn,
    )
    conn.close()

    assert catalog_table['mutation_id'].tolist() == ['NS3:2I', 'NS5A:1Y', 'NS5B:1F']
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    assert mutation_ids.get('seq1') == ['NS3:2I', 'NS5A:1Y', 'NS5B:1F']
    assert mutation_ids.get('seq2') == ['NS3:2I', 'NS5B:1F']
    completed_by_accession = _completed_signatures_by_accession(completed)
    assert completed_by_accession.get('seq1') == {'sig1', 'sig2', 'sig3'}
    assert completed_by_accession.get('seq2') == {'sig1', 'sig3'}


def test_annotate_mutations_can_fetch_genbank_gff_when_enabled(tmp_path, monkeypatch):
    db_path = tmp_path / "genbank_fallback.db"
    catalog_path = tmp_path / "genbank_catalog.tsv"

    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS5A', '1', '1', 'Y', 'REF_FETCH', 'NS5A:1Y', 'snp', 'sig1', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?)", [
        ('REF_FETCH', 'polyprotein', 4, 12),
        ('seq1', 'polyprotein', 4, 12),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_FETCH', 'REF_FETCH', 'REF_FETCH', 'GGGCCGCCACAC'),
        ('seq1', 'seq1', 'REF_FETCH', 'GGGCCCAAATAC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_FETCH', 'master'))
    cursor.execute("CREATE TABLE genes (description TEXT, display_name TEXT, name TEXT, parent_name TEXT)")
    cursor.execute("INSERT INTO genes VALUES (?, ?, ?, ?)", ('Non-structural protein 5A', 'NS5A', 'NS5A', 'whole_genome'))
    conn.commit()
    conn.close()

    def fake_fetch_reference_gff_map(accession, alias_lookup):
        assert accession == 'REF_FETCH'
        return {
            'NS5A': {
                'product': 'NS5A',
                'raw_product': 'NS5A',
                'cds_start': 10,
                'cds_end': 12,
                'reference_accession': accession,
                'feature_type': 'mat_peptide',
                'source': 'genbank_gff',
            }
        }

    monkeypatch.setattr(AnnotateMutations, 'fetch_reference_gff_map', fake_fetch_reference_gff_map)

    sys.argv = [
        'AnnotateMutations.py',
        '--db', str(db_path),
        '--mutation_catalog', str(catalog_path),
        '--catalog_column_profile', 'HCV',
        '--virus', '',
        '--allow_genbank_reference_gff',
    ]
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    completed = pd.read_sql_query("SELECT * FROM completed_signatures_only", conn)
    conn.close()

    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    assert mutation_ids == {'seq1': ['NS5A:1Y']}
    completed_by_accession = _completed_signatures_by_accession(completed)
    assert completed_by_accession == {'seq1': {'sig1'}}


def test_parse_gff_text_to_feature_map_accepts_mature_protein_regions():
    alias_lookup = {
        'nonstructuralproteinns5a': 'NS5A',
        'rnadependentrnapolymerasens5b': 'NS5B',
    }
    gff_text = (
        "##gff-version 3\n"
        "NC_004102.1\tRefSeq\tmature_protein_region_of_CDS\t6258\t7601\t.\t+\t.\tID=id-NS5A;product=nonstructural protein NS5A\n"
        "NC_004102.1\tRefSeq\tmature_protein_region_of_CDS\t7602\t9374\t.\t+\t.\tID=id-NS5B;product=RNA-dependent RNA polymerase NS5B\n"
    )

    feature_map = AnnotateMutations.parse_gff_text_to_feature_map(gff_text, 'NC_004102', alias_lookup, 'genbank_gff')

    assert feature_map['NS5A']['cds_start'] == 6258
    assert feature_map['NS5A']['cds_end'] == 7601
    assert feature_map['NS5B']['cds_start'] == 7602
    assert feature_map['NS5B']['cds_end'] == 9374


def test_annotate_mutations_crashes_when_reference_mapping_unavailable(tmp_path, capsys):
    db_path = tmp_path / "unmapped.db"
    catalog_path = tmp_path / "unmapped_catalog.tsv"

    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'I', 'REF_MASTER_NC_004102', 'NS3:2I', 'snp', 'sig1', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?)", [
        ('REF_MASTER', 'polyprotein', 4, 15),
        ('seq1', 'polyprotein', 4, 15),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_MASTER', 'REF_MASTER', 'REF_MASTER', 'GGGCCCAAATAC'),
        ('seq1', 'seq1', 'REF_MASTER', 'GGGCCCATAAAC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_MASTER', 'master'))
    cursor.execute("CREATE TABLE genes (description TEXT, display_name TEXT, name TEXT, parent_name TEXT)")
    cursor.execute("INSERT INTO genes VALUES (?, ?, ?, ?)", ('Non-structural protein 3', 'NS3', 'NS3', 'whole_genome'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    with pytest.raises(SystemExit) as excinfo:
        AnnotateMutations.main()

    captured = capsys.readouterr()
    assert excinfo.value.code == 2
    assert '--allow_genbank_reference_gff' in captured.err


def test_build_catalog_reference_table_hcv_profile_keeps_required_and_hcv_columns():
    catalog = pd.DataFrame([
        {
            'mutation_id': 'm1',
            'protein_name': 'NS3',
            'segment': '1',
            'aa_position': '2',
            'alt_residue': 'I',
            'reference_accession': 'REF1',
            'mutation_type': 'snp',
            'signature_id': 'sig1',
            'signature_kind': 'single',
            'combination_id': '',
            'combination_size': '',
            'phenotype': 'RAS',
            'resistance_category': 'III',
            'drug': 'drugA',
            'alignment_name': 'AL_1',
            '_canonical_protein': 'NS3',
        }
    ])

    result = AnnotateMutations.build_catalog_reference_table(catalog, 'HCV')

    assert result.columns.tolist() == HCV_DB_COLUMNS
    assert result.iloc[0]['phenotype'] == 'RAS'
    assert result.iloc[0]['drug'] == 'drugA'


def test_build_catalog_reference_table_all_columns_preserves_input_columns():
    catalog = pd.DataFrame([
        {
            'mutation_id': 'm1',
            'protein_name': 'HA',
            'segment': '4',
            'aa_position': '10',
            'alt_residue': 'N',
            'reference_accession': 'REF_FLU',
            'mutation_type': 'snp',
            'signature_id': 'sig1',
            'signature_kind': 'single',
            'combination_id': '',
            'combination_size': '',
            'phenotype': 'escape',
            'serotypes_tested': 'H1,H3',
            'alignment_name': 'flu_alignment',
            '_segment_norm': '4',
        }
    ])

    result = AnnotateMutations.build_catalog_reference_table(catalog, 'all_columns')

    assert result.columns.tolist() == [
        'mutation_id', 'protein_name', 'segment', 'aa_position', 'alt_residue',
        'reference_accession', 'mutation_type', 'signature_id', 'signature_kind',
        'combination_id', 'combination_size', 'phenotype', 'serotypes_tested', 'alignment_name'
    ]
    assert result.iloc[0]['serotypes_tested'] == 'H1,H3'


def test_write_mutation_tables_handles_empty_mutation_hits(tmp_path):
    conn = sqlite3.connect(str(tmp_path / 'empty_hits.db'))
    try:
        catalog = pd.DataFrame([
            {
                'mutation_id': 'm1',
                'protein_name': 'NS3',
                'segment': '1',
                'aa_position': '2',
                'alt_residue': 'I',
                'reference_accession': 'REF1',
                'mutation_type': 'snp',
                'signature_id': 'sig1',
                'signature_kind': 'single',
                'combination_id': '',
                'combination_size': '',
                'phenotype': '',
                'resistance_category': 'I',
                'drug': 'drugA',
            }
        ])

        AnnotateMutations.write_mutation_tables(conn, catalog, [], 'HCV')

        mutation_catalog = pd.read_sql_query('SELECT * FROM mutation_catalog', conn)
        relevant_summary = pd.read_sql_query('SELECT * FROM sequence_relevant_mutation_summary', conn)
        completed = pd.read_sql_query('SELECT * FROM completed_signatures_only', conn)

        assert mutation_catalog['mutation_id'].tolist() == ['m1']
        assert relevant_summary.empty
        assert relevant_summary.columns.tolist() == [
            'primary_accession',
            'relevant_mutations_present',
            'total_relevant_mutation_count',
        ]
        assert completed.empty
        assert completed.columns.tolist() == ['primary_accession', 'signature_id', 'signature_kind']
    finally:
        conn.close()


def test_stress_substitutions(tmp_path):
    db_path = tmp_path / "subst.db"
    catalog_path = tmp_path / "subst_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'V', 'REF_M', 'NS3:2V', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS3', '1', '2', 'I', 'REF_M', 'NS3:2I', 'snp', 'sig2', 'single', '', '', '', '', ''],
        ['NS3', '1', '3', 'K', 'REF_M', 'NS3:3K', 'snp', 'sig3', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_M', 'NS3', '1', 4, 12),
        ('seq1', 'polyprotein', '1', 4, 12),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_M', 'REF_M', 'REF_M', 'GGGATGATCAAG'),
        ('seq1', 'seq1', 'REF_M', 'GGGATGGTAaag'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_M', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    assert 'NS3:2V' in mutation_ids['seq1']
    assert 'NS3:3K' in mutation_ids['seq1']
    assert 'NS3:2I' not in mutation_ids['seq1']
    conn.close()


def test_stress_query_deletions(tmp_path):
    db_path = tmp_path / "deletions.db"
    catalog_path = tmp_path / "deletions_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '1', 'V', 'REF_M', 'NS3:1V', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS3', '1', '2', 'V', 'REF_M', 'NS3:2V', 'snp', 'sig2', 'single', '', '', '', '', ''],
        ['NS3', '1', '3', 'V', 'REF_M', 'NS3:3V', 'snp', 'sig3', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_M', 'NS3', '1', 1, 9),
        ('seq_del_mid', 'polyprotein', '1', 1, 9),
        ('seq_del_first', 'polyprotein', '1', 1, 9),
        ('seq_del_last', 'polyprotein', '1', 1, 9),
        ('seq_del_partial', 'polyprotein', '1', 1, 9),
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_M', 'REF_M', 'REF_M', 'ATGAAACCC'),
        ('seq_del_mid', 'seq_del_mid', 'REF_M', 'ATG---CCC'),
        ('seq_del_first', 'seq_del_first', 'REF_M', '---AAACCC'),
        ('seq_del_last', 'seq_del_last', 'REF_M', 'ATGAAA---'),
        ('seq_del_partial', 'seq_del_partial', 'REF_M', 'AT--AACCC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_M', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    assert 'seq_del_mid' not in mutation_ids or 'NS3:2V' not in mutation_ids['seq_del_mid']
    assert 'seq_del_first' not in mutation_ids or 'NS3:1V' not in mutation_ids['seq_del_first']
    assert 'seq_del_last' not in mutation_ids or 'NS3:3V' not in mutation_ids['seq_del_last']
    assert 'seq_del_partial' not in mutation_ids or 'NS3:1V' not in mutation_ids['seq_del_partial']
    conn.close()


def test_stress_query_insertions(tmp_path):
    db_path = tmp_path / "insertions.db"
    catalog_path = tmp_path / "insertions_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '1', 'M', 'REF_M', 'NS3:1M', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS3', '1', '2', 'V', 'REF_M', 'NS3:2V', 'snp', 'sig2', 'single', '', '', '', '', ''],
        ['NS3', '1', '3', 'P', 'REF_M', 'NS3:3P', 'snp', 'sig3', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_M', 'NS3', '1', 1, 9),
        ('seq_ins_3nt', 'polyprotein', '1', 1, 9),
        ('seq_ins_1nt', 'polyprotein', '1', 1, 9),
    ])
    
    ref_aln = 'ATG' + '---' + 'AAA' + 'CCC'
    seq_ins = 'ATG' + 'GGG' + 'GTA' + 'CCC'
    seq_ins_1nt = 'ATG' + 'G--' + 'GTA' + 'CCC'
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_M', 'REF_M', 'REF_M', ref_aln),
        ('seq_ins_3nt', 'seq_ins_3nt', 'REF_M', seq_ins),
        ('seq_ins_1nt', 'seq_ins_1nt', 'REF_M', seq_ins_1nt),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_M', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    assert 'NS3:2V' in mutation_ids['seq_ins_3nt']
    assert 'NS3:3P' in mutation_ids['seq_ins_3nt']
    conn.close()


def test_stress_boundary_coordinates(tmp_path):
    db_path = tmp_path / "boundary.db"
    catalog_path = tmp_path / "boundary_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '0', 'V', 'REF_M', 'NS3:0V', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS3', '1', '-1', 'V', 'REF_M', 'NS3:-1V', 'snp', 'sig2', 'single', '', '', '', '', ''],
        ['NS3', '1', '4', 'V', 'REF_M', 'NS3:4V', 'snp', 'sig3', 'single', '', '', '', '', ''],
        ['NS3', '1', 'abc', 'V', 'REF_M', 'NS3:abcV', 'snp', 'sig4', 'single', '', '', '', '', ''],
        ['NS3', '1', '2', 'V', 'REF_M', 'NS3:2V', 'snp', 'sig5', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_M', 'NS3', '1', 1, 9),
        ('seq1', 'polyprotein', '1', 1, 9),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_M', 'REF_M', 'REF_M', 'ATGAAACCC'),
        ('seq1', 'seq1', 'REF_M', 'ATGGTACCC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_M', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    assert 'NS3:2V' in mutation_ids['seq1']
    assert not any(m in mutation_ids['seq1'] for m in ['NS3:0V', 'NS3:-1V', 'NS3:4V', 'NS3:abcV'])
    conn.close()


@pytest.mark.skip(reason="Assumes multi-reference-coordinate resolution no longer supported under single-master system")
def test_stress_multiple_reference_groups(tmp_path):
    db_path = tmp_path / "multiref.db"
    catalog_path = tmp_path / "multiref_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'V', 'REF_A', 'NS3:2V_A', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS3', '1', '2', 'T', 'REF_B', 'NS3:2T_B', 'snp', 'sig2', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_A', 'NS3', '1', 4, 12),
        ('REF_B', 'NS3', '1', 1, 9),
        ('seq_a', 'polyprotein', '1', 4, 12),
        ('seq_b', 'polyprotein', '1', 1, 9),
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_A', 'REF_A', 'REF_A', 'GGGATGATGGGG'),
        ('seq_a', 'seq_a', 'REF_A', 'GGGATGGTAGGG'),
        ('REF_B', 'REF_B', 'REF_B', 'ATGATGGGG'),
        ('seq_b', 'seq_b', 'REF_B', 'ATGACGGGG'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?)", [
        ('REF_A', 'master'),
        ('REF_B', 'master'),
    ])
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    assert 'NS3:2V_A' in mutation_ids['seq_a']
    assert 'NS3:2T_B' in mutation_ids['seq_b']
    conn.close()


def test_stress_combination_signatures(tmp_path):
    db_path = tmp_path / "combo.db"
    catalog_path = tmp_path / "combo_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'V', 'REF_M', 'NS3:2V', 'snp', 'sig1', 'combination', 'C1', '2', '', 'I', 'DrugA'],
        ['NS3', '1', '3', 'T', 'REF_M', 'NS3:3T', 'snp', 'sig2', 'combination', 'C1', '2', '', 'I', 'DrugA'],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_M', 'NS3', '1', 1, 9),
        ('seq_full', 'polyprotein', '1', 1, 9),
        ('seq_partial', 'polyprotein', '1', 1, 9),
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_M', 'REF_M', 'REF_M', 'ATGAAACCC'),
        ('seq_full', 'seq_full', 'REF_M', 'ATGGTAACG'),
        ('seq_partial', 'seq_partial', 'REF_M', 'ATGGTACCC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_M', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    completed = pd.read_sql_query("SELECT * FROM completed_signatures_only", conn)
    completed_by_acc = _completed_signatures_by_accession(completed)
    
    assert 'seq_full' in completed_by_acc and 'sig1' in completed_by_acc['seq_full']
    assert 'seq_partial' not in completed_by_acc or 'sig1' not in completed_by_acc['seq_partial']
    conn.close()


@pytest.mark.skip(reason="Assumes multi-segment master coordinate resolution no longer supported under single-master system")
def test_stress_segmented_annotations(tmp_path):
    db_path = tmp_path / "segmented.db"
    catalog_path = tmp_path / "segmented_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'V', 'REF_M1', 'NS3:2V', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS5A', '4', '2', 'V', 'REF_M4', 'NS5A:2V', 'snp', 'sig2', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_M1', 'NS3', '1', 1, 9),
        ('REF_M4', 'NS5A', '4', 1, 9),
        ('seq_s1', 'polyprotein', '1', 1, 9),
        ('seq_s4', 'polyprotein', '4', 1, 9),
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_M1', 'REF_M1', 'REF_M1', 'ATGAAACCC'),
        ('seq_s1', 'seq_s1', 'REF_M1', 'ATGGTACCC'),
        ('REF_M4', 'REF_M4', 'REF_M4', 'ATGAAACCC'),
        ('seq_s4', 'seq_s4', 'REF_M4', 'ATGGTACCC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?)", [
        ('REF_M1', 'master'),
        ('REF_M4', 'master'),
    ])
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    assert 'NS3:2V' in mutation_ids['seq_s1']
    assert 'NS5A:2V' in mutation_ids['seq_s4']
    conn.close()


def test_coordinate_mapping_handles_deletions_indels(tmp_path):
    db_path = tmp_path / "indels.db"
    catalog_path = tmp_path / "indels_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', 'D', 'REF_INDEL', 'NS3:2D', 'snp', 'sig1', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    # In features, cds_start is sequence-specific start (1 in the ungapped sequence of REF_INDEL)
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_INDEL', 'NS3', '1', 1, 9),
        ('seq_q1', 'polyprotein', '1', 1, 12),
    ])
    
    # REF_INDEL has a deletion of 3 nt (gaps ---) at columns 3, 4, 5
    # Sequence-specific nucleotide indices map to columns:
    # 1->0, 2->1, 3->2 (ATG), 4->6, 5->7, 6->8 (GAC), 7->9, 8->10, 9->11 (CAA)
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_INDEL', 'REF_INDEL', 'REF_INDEL', 'ATG---GACCAA'),
        ('seq_q1', 'seq_q1', 'REF_INDEL', 'ATGTTTGACCAA'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_INDEL', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    # Codon 2 of NS3 in REF_INDEL starts at sequence index 4, which maps to alignment columns 6, 7, 8.
    # At columns 6, 7, 8, seq_q1 has 'GAC' -> translates to 'D'. 
    # Therefore, mutation 'NS3:2D' should be successfully annotated!
    assert 'NS3:2D' in mutation_ids['seq_q1']
    conn.close()


def test_coordinate_mapping_handles_truncations(tmp_path):
    db_path = tmp_path / "trunc.db"
    catalog_path = tmp_path / "trunc_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '1', 'K', 'REF_TRUNC', 'NS3:1K', 'snp', 'sig1', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    # REF_TRUNC is truncated at the start; its first base starts at master coordinate 4.
    # Its sequence-specific cds_start is 1 (starts at 'AAA').
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_TRUNC', 'NS3', '1', 1, 9),
        ('seq_q2', 'polyprotein', '1', 1, 12),
    ])
    
    # REF_TRUNC has gaps (---) at columns 0, 1, 2
    # Its sequence-specific nucleotide index 1 maps to column 3.
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_TRUNC', 'REF_TRUNC', 'REF_TRUNC', '---AAACCCGGG'),
        ('seq_q2', 'seq_q2', 'REF_TRUNC', 'ATGAAACCCGGG'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_TRUNC', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    # Codon 1 of NS3 in REF_TRUNC starts at sequence index 1, which maps to alignment columns 3, 4, 5.
    # At columns 3, 4, 5, seq_q2 has 'AAA' -> 'K'.
    # Therefore, mutation 'NS3:1K' should be successfully annotated!
    assert 'NS3:1K' in mutation_ids['seq_q2']
    conn.close()


def test_coordinate_mapping_detects_frame_shift_from_wrong_coordinates(tmp_path):
    db_path = tmp_path / "frameshift.db"
    catalog_path = tmp_path / "frameshift_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        # Catalog has NS3 codon 1 with expected residue K (lysine).
        ['NS3', '1', '1', 'K', 'REF_FS', 'NS3:1K', 'snp', 'sig1', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        # Set cds_start to 2 instead of 1, which is a 1 nt shift (frame shift!)
        ('REF_FS', 'NS3', '1', 2, 9),
        ('seq_q3', 'polyprotein', '1', 1, 9),
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_FS', 'REF_FS', 'REF_FS', 'ATGAAACCC'),
        ('seq_q3', 'seq_q3', 'REF_FS', 'ATGAAACCC'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_FS', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    # Because cds_start was set to 2 (offset), codon 1 starts at nucleotide position 2 (TGA instead of ATG).
    # 'TGA' translates to '*' (stop codon), not 'K'.
    # Therefore, mutation 'NS3:1K' should NOT be annotated for seq_q3 because of the frame shift.
    assert 'seq_q3' not in mutation_ids or 'NS3:1K' not in mutation_ids['seq_q3']
    conn.close()



def test_stress_translation_edge_cases(tmp_path):
    db_path = tmp_path / "translation.db"
    catalog_path = tmp_path / "translation_catalog.tsv"
    
    catalog_data = [
        HCV_CATALOG_HEADER,
        ['NS3', '1', '2', '_', 'REF_M', 'NS3:2Stop', 'snp', 'sig1', 'single', '', '', '', '', ''],
        ['NS3', '1', '3', 'X', 'REF_M', 'NS3:3X', 'snp', 'sig2', 'single', '', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_M', 'NS3', '1', 1, 9),
        ('seq_stop', 'polyprotein', '1', 1, 9),
        ('seq_invalid', 'polyprotein', '1', 1, 9),
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ('REF_M', 'REF_M', 'REF_M', 'ATGAAACCC'),
        ('seq_stop', 'seq_stop', 'REF_M', 'ATGTAACCC'),
        ('seq_invalid', 'seq_invalid', 'REF_M', 'ATGAAANNN'),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.execute("INSERT INTO meta_data VALUES (?, ?)", ('REF_M', 'master'))
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--catalog_column_profile', 'HCV', '--virus', '']
    AnnotateMutations.main()

    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    
    assert 'NS3:2Stop' in mutation_ids['seq_stop']
    assert 'NS3:3X' in mutation_ids['seq_invalid']
    conn.close()
