import pytest
import sqlite3
import pandas as pd
import os
import sys

# Ensure scripts/ is in PYTHONPATH when running pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../scripts')))

import VerifyMutations

def test_translate_codon():
    assert VerifyMutations.translate_codon('ATG') == 'M'
    assert VerifyMutations.translate_codon('TAA') == '_'
    assert VerifyMutations.translate_codon('TGA') == '_'
    assert VerifyMutations.translate_codon('TAG') == '_'
    assert VerifyMutations.translate_codon('A-G') == 'X'
    assert VerifyMutations.translate_codon('AT') == 'X'
    assert VerifyMutations.translate_codon('NNN') == 'X'

def test_infer_compact_product_name():
    assert VerifyMutations.infer_compact_product_name('protease/helicase protein NS3') == 'NS3'
    assert VerifyMutations.infer_compact_product_name('nonstructural protein NS5A') == 'NS5A'
    assert VerifyMutations.infer_compact_product_name('RNA-dependent RNA polymerase NS5B') == 'NS5B'
    assert VerifyMutations.infer_compact_product_name('envelope protein E2') == 'E2'
    assert VerifyMutations.infer_compact_product_name('p7 protein') == 'p7'
    assert VerifyMutations.infer_compact_product_name('Core protein') == 'Core'
    assert VerifyMutations.infer_compact_product_name('polyprotein') == 'Whole genome'
    assert VerifyMutations.infer_compact_product_name('wholegenome') == 'Whole genome'
    assert VerifyMutations.infer_compact_product_name('something_else') == ''

def test_canonicalize_product():
    alias_lookup = {
        'ns3protein': 'protease/helicase protein NS3',
        'nonstructuralprotein5a': 'nonstructural protein NS5A'
    }
    assert VerifyMutations.canonicalize_product('ns3protein', alias_lookup) == 'NS3'
    assert VerifyMutations.canonicalize_product('Non-structural protein 5A', alias_lookup) == 'NS5A'
    assert VerifyMutations.canonicalize_product('', {}) == ''

def test_extract_accession_tokens():
    assert VerifyMutations.extract_accession_tokens('NC_004102.1') == ['NC_004102.1', 'NC_004102']
    assert VerifyMutations.extract_accession_tokens('M58335') == ['M58335']
    assert VerifyMutations.extract_accession_tokens(None) == []

def test_build_alignment_coordinate_map():
    alignment = "A-T-G-C"
    coord_map = VerifyMutations.build_alignment_coordinate_map(alignment)
    # ungapped residues are: A (pos 1), T (pos 2), G (pos 3), C (pos 4)
    # alignment indices are: A -> 0, T -> 2, G -> 4, C -> 6
    assert coord_map == {1: 0, 2: 2, 3: 4, 4: 6}

def test_resolve_aligned_codon_indices():
    # ungapped sequence: ATGCCCGGG
    # alignment string: A-T-G-C-C-C-G-G-G
    # ungapped positions:
    # A(1)->0, T(2)->2, G(3)->4
    # C(4)->6, C(5)->8, C(6)->10
    # G(7)->12, G(8)->14, G(9)->16
    alignment = "A-T-G-C-C-C-G-G-G"
    coord_map = VerifyMutations.build_alignment_coordinate_map(alignment)
    
    # AA pos 1 = codons 1, 2, 3 -> align indices 0, 2, 4
    assert VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=1, aa_pos=1) == (0, 2, 4)
    # AA pos 2 = codons 4, 5, 6 -> align indices 6, 8, 10
    assert VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=1, aa_pos=2) == (6, 8, 10)
    # AA pos 3 = codons 7, 8, 9 -> align indices 12, 14, 16
    assert VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=1, aa_pos=3) == (12, 14, 16)
    
    # Boundary/invalid cases
    assert VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=-1, aa_pos=1) is None
    assert VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=1, aa_pos=0) is None
    assert VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start="abc", aa_pos=1) is None

def test_extract_aligned_codon():
    alignment = "ATGCCCGGG"
    assert VerifyMutations.extract_aligned_codon(alignment, (0, 1, 2)) == "ATG"
    assert VerifyMutations.extract_aligned_codon(alignment, (3, 4, 5)) == "CCC"
    assert VerifyMutations.extract_aligned_codon(alignment, (6, 7, 8)) == "GGG"
    
    # Invalid indices
    assert VerifyMutations.extract_aligned_codon(alignment, ()) is None
    assert VerifyMutations.extract_aligned_codon(alignment, (-1, 0, 1)) is None
    assert VerifyMutations.extract_aligned_codon(alignment, (7, 8, 9)) is None

def test_is_sequence_coverage_sufficient():
    # Sequence with no gaps
    ref_alignment = "ATGATGATGATGATGATGATGATGATGATG" # 10 codons
    ref_coord_map = VerifyMutations.build_alignment_coordinate_map(ref_alignment)
    gene_entry = {'cds_start': 1}
    
    # 1. Fully covered query
    query_alignment = "ATGATGATGATGATGATGATGATGATGATG"
    assert VerifyMutations.is_sequence_coverage_sufficient(
        query_alignment, ref_coord_map, gene_entry, start_aa=1, end_aa=4, target_pos=2
    ) is True
    
    # 2. Query with gap at target_pos (codon 2, cols 3,4,5)
    query_gap_target = "ATG---ATGATGATGATGATGATGATGATG"
    assert VerifyMutations.is_sequence_coverage_sufficient(
        query_gap_target, ref_coord_map, gene_entry, start_aa=1, end_aa=4, target_pos=2
    ) is False
    
    # 3. Query with too many gaps in window (6 gaps total out of 8 positions)
    # Window 1 to 8: (8 AAs). We have gaps at positions 3, 4, 5, 6, 7, 8
    query_too_many_gaps = "ATGATG" + "---" * 6 + "ATGATG"
    assert VerifyMutations.is_sequence_coverage_sufficient(
        query_too_many_gaps, ref_coord_map, gene_entry, start_aa=1, end_aa=8, target_pos=2
    ) is False
    
    # 4. Query with exactly 5 gaps in window (should pass)
    query_five_gaps = "ATGATG" + "---" * 5 + "ATGATGATG"
    assert VerifyMutations.is_sequence_coverage_sufficient(
        query_five_gaps, ref_coord_map, gene_entry, start_aa=1, end_aa=8, target_pos=2
    ) is True

def test_get_ref_alignment():
    ref_alignments_map = {
        ('REF1', 'REF1'): 'ATG',
        'REF2': 'CCC',
    }
    assert VerifyMutations.get_ref_alignment('REF1', ['REF2'], ref_alignments_map) == ('ATG', 'REF1')
    assert VerifyMutations.get_ref_alignment('REF2', ['REF1'], ref_alignments_map) == ('CCC', 'REF2')
    assert VerifyMutations.get_ref_alignment('REF_UNKNOWN', ['REF2'], ref_alignments_map) == ('CCC', 'REF2')
    assert VerifyMutations.get_ref_alignment('REF_UNKNOWN', ['REF_OTHER'], ref_alignments_map) == (None, None)

def test_resolve_feature_map():
    db_gff_maps = {
        'REF1': {
            'NS3': {'product': 'NS3'}
        },
        'REF2': {
            'NS3': {'product': 'NS3'}
        }
    }
    assert VerifyMutations.resolve_feature_map('REF1', set(), [], db_gff_maps) == ('REF1', db_gff_maps['REF1'])
    assert VerifyMutations.resolve_feature_map('REF_UNKNOWN', {'REF2'}, [], db_gff_maps) == ('REF2', db_gff_maps['REF2'])
    assert VerifyMutations.resolve_feature_map('REF_UNKNOWN', set(), ['REF1'], db_gff_maps) == ('REF1', db_gff_maps['REF1'])

def test_verify_mutations_end_to_end_success(tmp_path, capsys):
    db_path = tmp_path / "verify_success.db"
    catalog_path = tmp_path / "verify_catalog.tsv"
    
    # Catalog
    catalog_data = [
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id'],
        ['NS3', '1', '2', 'V', 'REF_MASTER', 'NS3:2V'],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?)", [
        ('REF_MASTER', 'master'),
        ('seq1', 'query'),
    ])
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_MASTER', 'NS3', '1', 4, 12),
        ('seq1', 'polyprotein', '1', 4, 12),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?)", [
        ('REF_MASTER', 'REF_MASTER', 'GGGATGATGGGG'),
        ('seq1', 'REF_MASTER', 'GGGATGGTAGGG'),
    ])
    cursor.execute("CREATE TABLE sequence_relevant_mutation_summary (primary_accession TEXT, relevant_mutations_present TEXT)")
    cursor.execute("INSERT INTO sequence_relevant_mutation_summary VALUES (?, ?)", ('seq1', 'NS3:2V;'))
    conn.commit()
    conn.close()
    
    # Mock sys.argv
    sys.argv = ['VerifyMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--sample_size', '1']
    
    with pytest.raises(SystemExit) as excinfo:
        VerifyMutations.main()
        
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "Result: VERIFICATION SUCCESSFUL" in captured.out
    assert "Mismatches:        0" in captured.out

def test_verify_mutations_end_to_end_mismatch(tmp_path, capsys):
    db_path = tmp_path / "verify_mismatch.db"
    catalog_path = tmp_path / "verify_catalog.tsv"
    
    # Catalog
    catalog_data = [
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id'],
        ['NS3', '1', '2', 'V', 'REF_MASTER', 'NS3:2V'],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?)", [
        ('REF_MASTER', 'master'),
        ('seq1', 'query'),
    ])
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ('REF_MASTER', 'NS3', '1', 4, 12),
        ('seq1', 'polyprotein', '1', 4, 12),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    # seq1 has ATG (M) instead of GTA (V) at codon 2
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?)", [
        ('REF_MASTER', 'REF_MASTER', 'GGGATGATGGGG'),
        ('seq1', 'REF_MASTER', 'GGGATGATGGGG'),
    ])
    cursor.execute("CREATE TABLE sequence_relevant_mutation_summary (primary_accession TEXT, relevant_mutations_present TEXT)")
    cursor.execute("INSERT INTO sequence_relevant_mutation_summary VALUES (?, ?)", ('seq1', 'NS3:2V;'))
    conn.commit()
    conn.close()
    
    sys.argv = ['VerifyMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--sample_size', '1']
    
    with pytest.raises(SystemExit) as excinfo:
        VerifyMutations.main()
        
    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "Result: VERIFICATION FAILED" in captured.err
    assert "Mismatches:        1" in captured.out

def test_verify_mutations_end_to_end_ref_mismatch(tmp_path, capsys):
    db_path = tmp_path / "verify_ref_mismatch.db"
    catalog_path = tmp_path / "verify_catalog.tsv"
    
    # Catalog reference is REF_OTHER, but in DB master is REF_MASTER
    catalog_data = [
        ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id'],
        ['NS3', '1', '2', 'V', 'REF_OTHER', 'NS3:2V'],
    ]
    with open(catalog_path, 'w') as f:
        for row in catalog_data:
            f.write('\t'.join(row) + '\n')
            
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?)", [
        ('REF_MASTER', 'master'),
    ])
    conn.commit()
    conn.close()
    
    sys.argv = ['VerifyMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--sample_size', '1']
    
    with pytest.raises(SystemExit) as excinfo:
        VerifyMutations.main()
        
    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "Result: Reference sequence flagging verification failed." in captured.err


def _make_catalog(catalog_path, rows):
    """Helper to write a catalog TSV file."""
    header = ['protein_name', 'segment', 'aa_position', 'alt_residue', 'reference_accession', 'mutation_id']
    with open(catalog_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')


def test_genotype_specific_reference_leading_gaps(tmp_path):
    """
    Test that a genotype-specific reference with leading gaps in the alignment
    correctly maps cds_start coordinates via the reference's own alignment row.
    
    Scenario: 
      - Master reference 'MASTER' has NS3 starting at genome position 100
      - Genotype reference 'GENO_REF' is in a separate alignment group,
        also has NS3 starting at position 100 of its own sequence,
        but in the padded alignment it starts at column 50 (50 leading dashes).
      - A query sequence 'SEQ1' is in GENO_REF's group and has mutation NS3:1M.
      - We verify that the codon lookup uses GENO_REF's alignment row (columns 100-102
        of the aligned string, accounting for leading gaps).
    """
    # GENO_REF alignment: 50 leading gaps, then sequence starting at ungapped pos 1
    # NS3 cds_start=100. Codon for AA pos 1 = ungapped positions 100,101,102
    # = alignment columns 100+50-1=149, 150, 151 (0-indexed: 149,150,151)
    leading_gaps = '-' * 50
    # Build a 200-nt sequence where positions 100-102 are ATG (M)
    genome = 'N' * 99 + 'ATG' + 'N' * 98   # 200 nt, ATG at pos 100-102 (1-indexed)
    geno_ref_aln = leading_gaps + genome  # total 250 chars
    # Query has ATG at same position -> M
    query_aln = leading_gaps + genome
    # Replace AA pos 1 codon in query with GTG -> V (mutation NS3:1V instead of M)
    query_list = list(query_aln)
    # ungapped pos 100 -> col 149, 101->150, 102->151
    query_list[149] = 'G'
    query_list[150] = 'T'
    query_list[151] = 'G'
    query_aln = ''.join(query_list)
    
    db_path = tmp_path / "geno_ref_test.db"
    catalog_path = tmp_path / "catalog.tsv"
    _make_catalog(catalog_path, [['NS3', '1', '1', 'V', 'MASTER', 'NS3:1V']])
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?)", [
        ('MASTER', 'master'),
        ('GENO_REF', 'query'),
        ('SEQ1', 'query'),
    ])
    # Features: GENO_REF has NS3 cds_start=100 (its own coordinate)
    cursor.execute("CREATE TABLE features (accession TEXT, reference_accession TEXT, master_ref_accession TEXT, product TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?, ?)", [
        ('MASTER', 'MASTER', 'MASTER', 'protease/helicase protein NS3', 100, 200),
        ('GENO_REF', 'MASTER', 'MASTER', 'protease/helicase protein NS3', 100, 200),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?)", [
        ('GENO_REF', 'GENO_REF', geno_ref_aln),
        ('SEQ1', 'GENO_REF', query_aln),
    ])
    cursor.execute("CREATE TABLE sequence_relevant_mutation_summary (primary_accession TEXT, relevant_mutations_present TEXT)")
    cursor.execute("INSERT INTO sequence_relevant_mutation_summary VALUES (?, ?)", ('SEQ1', 'NS3:1V'))
    conn.commit()
    conn.close()
    
    sys.argv = ['VerifyMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--sample_size', '1']
    
    with pytest.raises(SystemExit) as excinfo:
        VerifyMutations.main()
    
    assert excinfo.value.code == 0, "Should pass: SEQ1 has V at NS3:1 as expected"


def test_coordinate_map_uses_reference_alignment_row(tmp_path):
    """
    Regression test: verify that coordinate resolution correctly uses the
    alignment-based coordinate map (not raw genome coordinates), so that
    a reference with leading gaps maps cds_start correctly.
    
    This tests the fundamental correctness of build_alignment_coordinate_map +
    resolve_aligned_codon_indices when the reference has leading padding gaps.
    """
    # Reference has 10 leading gaps, then starts at ungapped pos 1
    ref_aln = '----------ATGCCCGGG'
    # NS3 cds_start = 1 (first nt of this reference), AA pos 1 = codons 1,2,3
    # ungapped pos 1 -> col 10, pos 2->11, pos 3->12
    coord_map = VerifyMutations.build_alignment_coordinate_map(ref_aln)
    
    # AA pos 1 with cds_start=1: nuc_start = 1, positions 1,2,3 -> cols 10,11,12
    indices = VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=1, aa_pos=1)
    assert indices == (10, 11, 12), f"Expected (10,11,12) but got {indices}"
    
    # AA pos 2: nuc_start = 4, positions 4,5,6 -> cols 13,14,15
    indices2 = VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=1, aa_pos=2)
    assert indices2 == (13, 14, 15), f"Expected (13,14,15) but got {indices2}"
    
    # Verify codon extraction works
    codon1 = VerifyMutations.extract_aligned_codon(ref_aln, indices)
    assert codon1 == 'ATG'
    assert VerifyMutations.translate_codon(codon1) == 'M'
    
    codon2 = VerifyMutations.extract_aligned_codon(ref_aln, indices2)
    assert codon2 == 'CCC'
    assert VerifyMutations.translate_codon(codon2) == 'P'


def test_truncated_sequences_excluded_from_visual_window():
    """
    Test that is_sequence_coverage_sufficient returns False for sequences
    that are all-gap at the target position, so they are excluded from
    the visual display (not shown as X-filled rows with 0% identity).
    """
    # Reference alignment: 30 codons (90 nt), no gaps
    ref_aln = 'ATG' * 30  # 90 nt, cds_start=1
    ref_coord_map = VerifyMutations.build_alignment_coordinate_map(ref_aln)
    gene_entry = {'cds_start': 1}
    
    # Truncated sequence: only has first 30 nt, rest are gaps
    truncated_aln = 'ATG' * 10 + '-' * 60
    
    # Target at position 15 (codon 15 = nt 43-45), which is in the gap region
    assert not VerifyMutations.is_sequence_coverage_sufficient(
        truncated_aln, ref_coord_map, gene_entry, start_aa=10, end_aa=20, target_pos=15
    ), "Truncated sequence should fail coverage check at gapped position"
    
    # But a sequence with a small number of gaps elsewhere (<=5) and real data at target should pass
    # Gaps only at positions 1-4, target at 15 (real data)
    mostly_covered = '-' * 12 + 'ATG' * 27  # first 4 codons gapped, rest real
    assert VerifyMutations.is_sequence_coverage_sufficient(
        mostly_covered, ref_coord_map, gene_entry, start_aa=5, end_aa=20, target_pos=15
    ), "Sequence with <=5 gaps in window should pass if target is covered"


def test_verify_mutations_genotype_specific_features(tmp_path, capsys):
    """
    End-to-end test: two alignment groups (MASTER and GENO_REF) with different
    cds_start values for the same protein. Verify that each group uses its own
    reference's cds_start coordinates, not the master's.
    
    Bug being tested: Previously, AnnotateMutations annotated using the correct
    genotype-specific cds_start, but if VerifyMutations used the wrong cds_start
    (e.g., master's), it would look up the wrong codon and report a mismatch.
    """
    # Master: NS3 cds_start=10, Query REF: NS3 cds_start=4
    # Master alignment: GGGGGGGGGGATGCCC (NS3 starts at pos 10, codon 1 = ATG = M)
    # GENO_REF alignment: GGGATGCCC (NS3 starts at pos 4, codon 1 = ATG = M)
    # Both references have M at NS3:1 in their OWN coordinate frame
    # SEQ1 is in GENO_REF group and has GTA (V) at NS3:1 position
    
    master_aln = 'GGGGGGGGG' + 'ATG' + 'CCC'   # 15 nt, NS3 cds_start=10
    geno_ref_aln = 'GGG' + 'ATG' + 'CCC'        # 9 nt, NS3 cds_start=4
    # SEQ1 in GENO_REF group: GTA at positions 4-6 -> V
    seq1_aln = 'GGG' + 'GTA' + 'CCC'
    
    db_path = tmp_path / "geno_specific.db"
    catalog_path = tmp_path / "catalog.tsv"
    _make_catalog(catalog_path, [['NS3', '1', '1', 'V', 'MASTER', 'NS3:1V']])
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?)", [
        ('MASTER', 'master'),
        ('GENO_REF', 'query'),
        ('SEQ1', 'query'),
    ])
    cursor.execute("CREATE TABLE features (accession TEXT, reference_accession TEXT, master_ref_accession TEXT, product TEXT, cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?, ?)", [
        ('MASTER', 'MASTER', 'MASTER', 'protease/helicase protein NS3', 10, 15),
        ('GENO_REF', 'MASTER', 'MASTER', 'protease/helicase protein NS3', 4, 9),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?)", [
        ('MASTER', 'MASTER', master_aln),
        ('GENO_REF', 'GENO_REF', geno_ref_aln),
        ('SEQ1', 'GENO_REF', seq1_aln),
    ])
    cursor.execute("CREATE TABLE sequence_relevant_mutation_summary (primary_accession TEXT, relevant_mutations_present TEXT)")
    cursor.execute("INSERT INTO sequence_relevant_mutation_summary VALUES (?, ?)", ('SEQ1', 'NS3:1V'))
    conn.commit()
    conn.close()
    
    sys.argv = ['VerifyMutations.py', '--db', str(db_path), '--mutation_catalog', str(catalog_path), '--sample_size', '1']
    
    with pytest.raises(SystemExit) as excinfo:
        VerifyMutations.main()
    
    assert excinfo.value.code == 0, "Should pass: SEQ1 has V at NS3:1 using GENO_REF's cds_start=4"
    captured = capsys.readouterr()
    assert 'Mismatches:        0' in captured.out


def test_resolve_aligned_codon_indices_with_cds_offset():
    """
    Test resolve_aligned_codon_indices with cds_start > 1 (non-trivial CDS offset).
    This simulates a gene that starts partway into the genome sequence.
    
    E.g.: genome = NNNNNNATGCCCGGG (6 leading Ns, then NS3 starts at pos 7)
    cds_start = 7, AA pos 1 -> nuc positions 7,8,9 -> alignment indices for those
    """
    genome = 'NNNNNN' + 'ATG' + 'CCC' + 'GGG'  # 15 nt
    coord_map = VerifyMutations.build_alignment_coordinate_map(genome)
    
    # cds_start=7, AA pos 1: nuc_start = 7 + 0 = 7 -> positions 7,8,9 -> cols 6,7,8 (0-indexed)
    indices1 = VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=7, aa_pos=1)
    assert indices1 == (6, 7, 8), f"Expected (6,7,8) but got {indices1}"
    codon1 = VerifyMutations.extract_aligned_codon(genome, indices1)
    assert codon1 == 'ATG'
    assert VerifyMutations.translate_codon(codon1) == 'M'
    
    # AA pos 2: nuc_start = 7 + 3 = 10 -> positions 10,11,12 -> cols 9,10,11
    indices2 = VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=7, aa_pos=2)
    assert indices2 == (9, 10, 11), f"Expected (9,10,11) but got {indices2}"
    codon2 = VerifyMutations.extract_aligned_codon(genome, indices2)
    assert codon2 == 'CCC'
    assert VerifyMutations.translate_codon(codon2) == 'P'


def test_annotate_mutations_matches_verify_mutations_coordinates(tmp_path):
    """
    Critical regression test: AnnotateMutations and VerifyMutations must use the 
    SAME alignment indices to find a mutation's codon. This tests the scenario
    where the genotype-specific reference has a different cds_start than the master.
    
    If AnnotateMutations finds NS3:1V using GENO_REF's cds_start=4 (cols 3,4,5),
    then VerifyMutations must also use cols 3,4,5, NOT the master's cds_start=10 (cols 9,10,11).
    """
    # Genotype reference: NS3 at cds_start=4, alignment = 'GGGATGCCC'
    geno_ref_aln = 'GGG' + 'ATG' + 'CCC'  # ungapped pos 4,5,6 = ATG = M at NS3:1
    coord_map = VerifyMutations.build_alignment_coordinate_map(geno_ref_aln)
    
    # Using GENO_REF's cds_start=4:
    indices_correct = VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=4, aa_pos=1)
    assert indices_correct == (3, 4, 5), f"Correct indices should be (3,4,5), got {indices_correct}"
    
    # If we erroneously used MASTER's cds_start=10:
    indices_wrong = VerifyMutations.resolve_aligned_codon_indices(coord_map, cds_start=10, aa_pos=1)
    assert indices_wrong is None or indices_wrong != (3, 4, 5), \
        "Using master's cds_start should give different/invalid result"
    
    # The correct codon in GENO_REF at NS3:1 = 'ATG' = M
    codon = VerifyMutations.extract_aligned_codon(geno_ref_aln, indices_correct)
    assert codon == 'ATG'
    assert VerifyMutations.translate_codon(codon) == 'M'
    
    # A query sequence with GTA (V) at NS3:1 (using same alignment frame as GENO_REF)
    query_aln = 'GGG' + 'GTA' + 'CCC'
    query_codon = VerifyMutations.extract_aligned_codon(query_aln, indices_correct)
    assert query_codon == 'GTA'
    assert VerifyMutations.translate_codon(query_codon) == 'V'
