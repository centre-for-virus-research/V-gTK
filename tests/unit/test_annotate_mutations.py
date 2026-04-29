import pytest
import sqlite3
import pandas as pd
import os
import sys

# Assume scripts/ is in PYTHONPATH when running pytest, but adding it anyway
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../scripts')))

import AnnotateMutations


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
    assert AnnotateMutations.translate_codon('TAA') == '_'

def test_annotate_mutations_db_creation(tmp_path):
    db_path = tmp_path / "test.db"
    catalog_path = tmp_path / "catalog.tsv"
    
    # 1. Create a dummy catalog
    catalog_data = [
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id', 'mutation_type', 'signature_id', 'signature_kind', 'combination_id', 'combination_size', 'resistance_category', 'drug'],
        ['NS3', '1', '2', 'I', 'REF1', 'NS3:2I', 'snp', 'sig1', 'single', 'C1', '2', 'I', 'DrugA'],
        ['NS3', '1', '3', 'K', 'REF1', 'NS3:3K', 'snp', 'sig2', 'combination', 'C1', '2', 'I', 'DrugA'],
        ['NS5A', '', '1', 'Y', 'REF1', 'NS5A:1Y', 'snp', 'sig3', 'single', '', '', '', '']
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
        ('seq1', 'NS3', '1', 4, 12),  # Nucleotides 4-12 are NS3
        ('seq1', 'NS5A', '', 13, 15), # Nucleotides 13-15 are NS5A
        ('seq2', 'NS3', '1', 4, 12)
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, alignment TEXT)")
    # seq1 padded alignment:
    # 1-3: ABC
    # 4-6: ATA (AA 1: I), 7-9: ATC (AA 2: I), 10-12: AAA (AA 3: K)  => NS3 has 2I, 3K
    # 13-15: TAC (AA 1: Y) => NS5A has 1Y
    cursor.execute("INSERT INTO sequence_alignment VALUES (?, ?)", ('seq1', 'ABCATAATCAAATAC'))
    
    # seq2 padded alignment:
    # 1-3: ABC
    # 4-6: ATA (I), 7-9: GGG (G), 10-12: AAA (K) => NS3 has 3K but not 2I
    cursor.execute("INSERT INTO sequence_alignment VALUES (?, ?)", ('seq2', 'ABCATAGGGAAA'))
    conn.commit()
    conn.close()
    
    # 3. Call the script logic
    # Mocking sys.argv
    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--virus', 'HCV']
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
        'combination_id', 'combination_size', 'resistance_category', 'drug'
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
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id', 'mutation_type', 'signature_id', 'signature_kind', 'combination_id', 'combination_size', 'resistance_category', 'drug'],
        # Reference AA is M (ATG), we look for a mutation to V.
        ['NS3', '1', '2', 'V', 'REF1', 'NS3:2V', 'snp', 'sig1', 'single', '', '', '', ''],
        # Reference AA is H (CAT) at pos 4 (pos 3 is a gaping deletion). We look for Q.
        ['NS3', '1', '4', 'Q', 'REF1', 'NS3:4Q', 'snp', 'sig2', 'single', '', '', '', '']
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    # Let's say cds_start is 4
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_MASTER', 'NS3', '1', 4, 15),
        ('mutated_seq', 'NS3', '1', 4, 15)
    ])
    
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, alignment TEXT)")
    
    # Padded Alignment layout: (1-based indices relative to alignment string)
    # 123 | 456 | 789 | 10,11,12 | 13,14,15
    # pos | AA1 | AA2 |   AA3    |    AA4
    # Pad | ATG | ATG |   ---    |    CAT
    ref_aln = "GGG" + "ATG" + "ATG" + "---" + "CAT"
    
    # Mutated has V at AA2 (GTA) and Q at AA4 (CAA). And still has deletion at AA3.
    mut_aln = "GGG" + "ATG" + "GTA" + "---" + "CAA"
    
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?)", [
        ('REF_MASTER', ref_aln),
        ('mutated_seq', mut_aln)
    ])
    conn.commit()
    conn.close()
    
    # 2. Run annotation
    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--virus', '']
    AnnotateMutations.main()
    
    # 3. Verify
    conn = sqlite3.connect(str(db_path))
    mutation_summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    completed = pd.read_sql_query("SELECT * FROM completed_signatures_only", conn)
    
    mutation_ids = _mutation_ids_by_accession(mutation_summary)
    assert 'REF_MASTER' not in mutation_ids, "Reference should not trigger any mutation hits"
    assert mutation_ids['mutated_seq'] == ['NS3:2V', 'NS3:4Q']

    completed_by_accession = _completed_signatures_by_accession(completed)
    assert completed_by_accession['mutated_seq'] == {'sig1', 'sig2'}
    
    conn.close()


def test_annotate_mutations_falls_back_to_db_gff_coordinate_mapping(tmp_path):
    db_path = tmp_path / "fallback_test.db"
    catalog_path = tmp_path / "fallback_catalog.tsv"

    catalog_data = [
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id', 'mutation_type', 'signature_id', 'signature_kind', 'combination_id', 'combination_size', 'resistance_category', 'drug'],
        ['NS3', '1', '2', 'I', 'REF_MASTER_NC_004102', 'NS3:2I', 'snp', 'sig1', 'single', '', '', '', ''],
        ['NS5A', '1', '1', 'Y', 'REF_MASTER_NC_004102', 'NS5A:1Y', 'snp', 'sig2', 'single', '', '', '', ''],
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

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--virus', '']
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
    assert set(db_gff_maps) == {'NC_004102'}

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
    assert resolved_maps['EU781827']['resolved_accession'] == 'NC_004102'
    assert diagnostics['mutation_hits'] == 2


def test_annotate_mutations_uses_descriptive_hcv_feature_products_without_gene_aliases(tmp_path):
    db_path = tmp_path / "descriptive_hcv_features.db"
    catalog_path = tmp_path / "descriptive_hcv_catalog.tsv"

    catalog_data = [
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id', 'mutation_type', 'signature_id', 'signature_kind', 'combination_id', 'combination_size', 'resistance_category', 'drug'],
        ['NS3', '1', '2', 'I', 'REF_MASTER_NC_004102', 'NS3:2I', 'snp', 'sig1', 'single', '', '', '', ''],
        ['NS5A', '1', '1', 'Y', 'REF_MASTER_NC_004102', 'NS5A:1Y', 'snp', 'sig2', 'single', '', '', '', ''],
        ['NS5B', '1', '1', 'F', 'REF_MASTER_NC_004102', 'NS5B:1F', 'snp', 'sig3', 'single', '', '', '', ''],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('seq1', 'protease/helicase protein NS3', '1', 1, 6),
        ('seq1', 'nonstructural protein NS5A', '1', 7, 9),
        ('seq1', 'RNA-dependent RNA polymerase NS5B', '1', 10, 12),
        ('seq2', 'protease/helicase protein NS3', '1', 1, 6),
        ('seq2', 'nonstructural protein NS5A', '1', 7, 9),
        ('seq2', 'RNA-dependent RNA polymerase NS5B', '1', 10, 12),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?)", [
        ('seq1', 'ATAATCTACTTT'),
        ('seq2', 'ATAATCGGGTTC'),
    ])
    conn.commit()
    conn.close()

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--virus', '']
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
    assert mutation_ids == {
        'seq1': ['NS3:2I', 'NS5A:1Y', 'NS5B:1F'],
        'seq2': ['NS3:2I', 'NS5B:1F'],
    }
    completed_by_accession = _completed_signatures_by_accession(completed)
    assert completed_by_accession == {
        'seq1': {'sig1', 'sig2', 'sig3'},
        'seq2': {'sig1', 'sig3'},
    }


def test_annotate_mutations_can_fetch_genbank_gff_when_enabled(tmp_path, monkeypatch):
    db_path = tmp_path / "genbank_fallback.db"
    catalog_path = tmp_path / "genbank_catalog.tsv"

    catalog_data = [
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id', 'mutation_type', 'signature_id', 'signature_kind', 'combination_id', 'combination_size', 'resistance_category', 'drug'],
        ['NS5A', '1', '1', 'Y', 'REF_FETCH', 'NS5A:1Y', 'snp', 'sig1', 'single', '', '', '', ''],
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
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id', 'mutation_type', 'signature_id', 'signature_kind', 'combination_id', 'combination_size', 'resistance_category', 'drug'],
        ['NS3', '1', '2', 'I', 'REF_MASTER_NC_004102', 'NS3:2I', 'snp', 'sig1', 'single', '', '', '', ''],
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

    sys.argv = ['AnnotateMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--virus', '']
    with pytest.raises(SystemExit) as excinfo:
        AnnotateMutations.main()

    captured = capsys.readouterr()
    assert excinfo.value.code == 2
    assert '--allow_genbank_reference_gff' in captured.err
