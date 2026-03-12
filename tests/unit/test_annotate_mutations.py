import pytest
import sqlite3
import pandas as pd
import os
import sys

# Assume scripts/ is in PYTHONPATH when running pytest, but adding it anyway
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../scripts')))

import AnnotateMutations

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
    mut = pd.read_sql_query("SELECT * FROM sequence_mutations", conn)
    
    # seq1 should have NS3:2I, NS3:3K, NS5A:1Y
    seq1_muts = set(mut[mut['primary_accession']=='seq1']['mutation_id'])
    assert 'NS3:2I' in seq1_muts
    assert 'NS3:3K' in seq1_muts
    assert 'NS5A:1Y' in seq1_muts
    
    # seq2 should have NS3:3K
    seq2_muts = list(mut[mut['primary_accession']=='seq2']['mutation_id'])
    assert ['NS3:3K'] == seq2_muts
    assert mut[(mut['primary_accession']=='seq2') & (mut['mutation_id']=='NS3:3K')]['combination_id'].iloc[0] == 'C1'
    
    # check drug resistance logic
    dr = pd.read_sql_query("SELECT * FROM sequence_drug_resistance", conn)
    assert len(dr) == 2
    
    # seq1 has both 2I and 3K, so C1 is complete (req=2)
    s1_dr = dr[dr['primary_accession']=='seq1'].iloc[0]
    assert s1_dr['combination_id'] == 'C1'
    assert s1_dr['mutations_detected'] == 2
    assert s1_dr['combination_status'] == 'complete'
    
    # seq2 only has 3K, so C1 is partial
    s2_dr = dr[dr['primary_accession']=='seq2'].iloc[0]
    assert s2_dr['mutations_detected'] == 1
    assert s2_dr['combination_status'] == 'partial'
    
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
    mut = pd.read_sql_query("SELECT * FROM sequence_mutations", conn)
    
    # REF_MASTER should have NO mutations
    ref_muts = mut[mut['primary_accession']=='REF_MASTER']
    assert len(ref_muts) == 0, "Reference should not trigger any mutation hits"
    
    seq_muts = set(mut[mut['primary_accession']=='mutated_seq']['mutation_id'])
    assert 'NS3:2V' in seq_muts, "Failed to detect mutation at position 2"
    assert 'NS3:4Q' in seq_muts, "Failed to detect mutation at position 4 after an indel"
    
    conn.close()
