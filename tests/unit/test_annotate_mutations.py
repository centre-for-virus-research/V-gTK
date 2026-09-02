"""Unit and adversarial coverage for scripts/AnnotateMutations.py.

Lens: the annotator turns a padded alignment plus a curated catalogue into
resistance calls that a clinician-facing database then reports as fact.  The
failure that matters is therefore never a traceback - it is a *plausible* row.

Section 0 holds the core behaviour: codon translation, database creation, gap
accounting, reference indels, GFF and features-table coordinate mapping, the
catalogue reference table, and the stress cases.

Sections 1 to 14 hunt three further shapes, each one an instance of "the safety
gate is keyed on something the data does not always spell the same way":

  * missing data promoted to evidence - a sequence with nothing sequenced at a
    position still produces a call;
  * a suppression gate that silently stops firing - the wild-type and genotype
    gates key on raw catalogue text, so a re-spelled position or a non-numeric
    genotype vocabulary disables them without any diagnostic;
  * a virus-specific assumption reached through a virus-agnostic door - the
    `influenza` and `all_columns` catalogue profiles, the `segment` column and
    non-HCV product names.

Each of those began as an xfail describing a defect found by reading and
running the annotator; the defects are fixed and the markers are gone, so a
regression now fails loudly instead of being absorbed by a tolerated xfail.
The docstrings are kept in the past tense on purpose - what the code used to
do, and what real input triggered it, is the reason the assertion is worth its
runtime.

Companion suites: test_hcv_mutation_edge_cases.py (catalogue and coordinate
arithmetic), test_hcv_genotype_gating.py (rules A and B), and the HCV
invariance check described in TESTING.md.

Everything here is synthetic and writes only into tmp_path.
"""

import contextlib
import io
import os
import sqlite3
import sys

import pandas as pd
import pytest

# Assume scripts/ is in PYTHONPATH when running pytest, but adding it anyway
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../scripts')))

import AnnotateMutations
import ExploreMutationStorageLayouts as EL


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


# The column order main() expects from a catalogue TSV, plus the two HCV
# profile columns, so the same header serves every end-to-end case.
CATALOG_HEADER = [
    "mutation_id", "protein_name", "segment", "aa_position", "alt_residue",
    "reference_accession", "mutation_type", "signature_id", "signature_kind",
    "combination_id", "combination_size", "phenotype", "resistance_category", "drug",
]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _catalog_row(**overrides):
    """A minimal catalogue row carrying every column AnnotateMutations requires."""
    row = {
        "mutation_id": "NS3:2K",
        "protein_name": "NS3",
        "segment": "1",
        "aa_position": "2",
        "alt_residue": "K",
        "reference_accession": "REF1",
        "mutation_type": "aminoAcidSimplePolymorphism",
        "signature_id": "NS3:2K",
        "signature_kind": "single",
        "combination_id": "",
        "combination_size": "",
        "phenotype": "",
        "resistance_category": "",
        "drug": "",
    }
    row.update(overrides)
    return row


def _annotate(catalog_rows, alignments, feature_map, meta_data=None, master="REF1"):
    """Run the annotator over a tiny synthetic alignment.

    `alignments` is a list of (accession, padded_alignment) pairs; the first is
    the master reference and supplies every row's alignment_name.  Returns
    (mutations_found, diagnostics-as-dict).
    """
    seq_aln = pd.DataFrame(
        [
            {
                "sequence_id": accession,
                "primary_accession": accession,
                "alignment": alignment,
                "alignment_name": alignments[0][0],
            }
            for accession, alignment in alignments
        ]
    )
    if meta_data is None:
        meta_data = pd.DataFrame([{"primary_accession": master, "accession_type": "master"}])
    prepared, _invalid = AnnotateMutations.prepare_catalog(pd.DataFrame(catalog_rows), {})
    with contextlib.redirect_stdout(io.StringIO()):
        found, diagnostics, _maps = AnnotateMutations.annotate_from_reference_coordinates(
            prepared, seq_aln, meta_data, {}, {master: feature_map}, False
        )
    return found, dict(diagnostics)


def _write_catalog(path, rows, header=CATALOG_HEADER):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row.get(column, "") for column in header) + "\n")
    return path


def _build_db(path, *, features, alignments, meta, extra_tables=()):
    """Materialise the three tables the annotator reads, plus any extras.

    `features`/`alignments`/`meta` are (create-sql, rows) pairs so a case can
    choose its own column set - schema drift is half of what is being tested.
    """
    conn = sqlite3.connect(str(path))
    try:
        cursor = conn.cursor()
        for create_sql, rows in (features, alignments, meta) + tuple(extra_tables):
            cursor.execute(create_sql)
            if rows:
                placeholders = ", ".join("?" * len(rows[0]))
                table = create_sql.split()[2]
                cursor.executemany(f"INSERT INTO {table} VALUES ({placeholders})", rows)
        conn.commit()
    finally:
        conn.close()
    return path


def _run_main(db_path, catalog_path, monkeypatch, profile="HCV", virus="", extra_args=()):
    """Invoke the CLI entry point with stdout captured."""
    monkeypatch.setattr(
        sys, "argv",
        ["AnnotateMutations.py", "--db", str(db_path), "--mutation_catalog", str(catalog_path),
         "--catalog_column_profile", profile, "--virus", virus, *extra_args],
    )
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        AnnotateMutations.main()
    return buffer.getvalue()


def _summary(db_path):
    """{accession: [mutation_id, ...]} from sequence_relevant_mutation_summary."""
    conn = sqlite3.connect(str(db_path))
    try:
        frame = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    finally:
        conn.close()
    return {
        row["primary_accession"]: [
            value for value in str(row["relevant_mutations_present"] or "").split(";") if value
        ]
        for _, row in frame.iterrows()
    }


NS3_FEATURE_MAP = {
    "NS3": {"product": "NS3", "cds_start": 1, "cds_end": 12, "feature_type": "mat_peptide"}
}


# --------------------------------------------------------------------------
# 0. Core behaviour: translation, database creation, coordinate mapping
# --------------------------------------------------------------------------

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
    # 'NNN' is an unsequenced codon, so the residue is unresolved and must not
    # satisfy the 'X' catalog row.  This assertion used to read the other way
    # round, pinning the behaviour that
    # test_hcv_mutation_edge_cases.test_ambiguous_or_missing_codon_does_not_satisfy_an_x_catalog_row
    # documents as a defect: an 'X' alt residue turned missing data into a
    # positive call and so flagged the worst-covered sequences hardest.
    assert 'seq_invalid' not in mutation_ids or 'NS3:3X' not in mutation_ids['seq_invalid']
    conn.close()


# --------------------------------------------------------------------------
# 1. Missing data promoted to a positive call
# --------------------------------------------------------------------------

def test_a_deletion_inside_the_covered_span_is_still_called():
    """Guard-rail for section 1: a genuine deletion must survive any fix here.

    Every test below asks the annotator to call *fewer* deletions.  This pins
    the case it must keep calling, so the cheap fix - stop reading gaps as
    deletions at all - fails loudly instead of passing quietly.
    """
    found, _diagnostics = _annotate(
        [_catalog_row(alt_residue="del", mutation_id="NS3:2-")],
        [("REF1", "ATGAAACCCGGG"), ("DELETED", "ATG---CCCGGG")],
        NS3_FEATURE_MAP,
    )

    assert [hit["primary_accession"] for hit in found] == ["DELETED"]
    assert found[0]["observed_residue"] == AnnotateMutations.DELETION_RESIDUE


def test_terminal_padding_is_not_called_as_a_deletion():
    """Pins the guard that already works: padding beyond the covered span.

    A partial GenBank record is squared up to the full genome width with gaps.
    Those columns are sequence the submitter never reported, and reading them
    as deletions would flag the worst-covered records the hardest.
    """
    found, _diagnostics = _annotate(
        [_catalog_row(alt_residue="del", mutation_id="NS3:2-")],
        [("REF1", "ATGAAACCCGGG"), ("PARTIAL", "---" * 2 + "CCCGGG")],
        NS3_FEATURE_MAP,
    )

    assert found == [], "gaps outside the covered span are missing data, not a deletion"


def test_a_sequence_with_nothing_sequenced_gets_no_deletion_calls():
    """A row with zero reported bases is the strongest possible missing-data case.

    Real-world trigger: a sequence that failed alignment, or was padded into
    the matrix before its bases were merged in, is stored as an all-gap row of
    full genome width.  alignment_covered_span() has no covered span to return
    and answers None; residue_from_aligned_codon() treats None as "a unit test
    handed me a bare codon" and returns '-'.  Every catalogued deletion then
    fires on the one sequence that carries no evidence whatsoever - the exact
    inversion of the terminal-padding rule immediately above, which exists to
    stop *partially* unsequenced records being called.
    """
    found, _diagnostics = _annotate(
        [_catalog_row(alt_residue="del", mutation_id="NS3:2-")],
        [("REF1", "ATGAAACCCGGG"), ("ALL_GAP", "-" * 12)],
        NS3_FEATURE_MAP,
    )

    assert found == [], (
        "a sequence with no sequenced bases reported "
        f"{[hit['mutation_id'] for hit in found]}"
    )


# --------------------------------------------------------------------------
# 2. Rule B (wild type) silently stops firing
# --------------------------------------------------------------------------

WILD_TYPE_META = pd.DataFrame(
    [
        {"primary_accession": "REF1", "accession_type": "master", "genotype": "1", "subtype": "a"},
        {"primary_accession": "SEQ1", "accession_type": "query", "genotype": "1", "subtype": "a"},
    ]
)


def _wild_type_case(aa_position, wild_type_residues="1a:K:99.0"):
    """Codon 2 of the query reads K, which is also the catalogued 1a wild type."""
    return _annotate(
        [
            _catalog_row(
                aa_position=aa_position,
                alt_residue="K",
                relevant_genotypes="1a",
                wild_type_residues=wild_type_residues,
            )
        ],
        [("REF1", "ATGAAACCCGGG"), ("SEQ1", "ATGAAACCCGGG")],
        NS3_FEATURE_MAP,
        meta_data=WILD_TYPE_META,
    )


def test_wild_type_anchor_is_suppressed_for_a_plainly_spelled_position():
    """Control for the two cases below: rule B works when the spelling agrees."""
    found, _diagnostics = _wild_type_case("2")

    assert found == [], "K is the 1a wild type at NS3:2; it is not a resistance call"


def test_wild_type_suppression_survives_a_float_spelled_position():
    """A position respelled as '2.0' must still suppress its own wild type.

    Real-world trigger: the catalogue is regenerated through pandas and
    spreadsheets, where an integer column comes back as '107.0'.
    prepare_catalog() already absorbs that for coordinates - int(float(value))
    - so the codon is still read at the right place and the call still fires.
    Only rule B breaks: build_wild_type_tables() keys the table on
    clean_cell(aa_position) == '2.0' while lookup_wild_type_residues() asks for
    str(2).  The lookup misses, the residue is labelled 'wt_unknown', and the
    genotype's own wild-type residue is emitted as a resistance call against
    every sequence that carries it - including the reference.  Nothing in the
    diagnostics distinguishes this from a genuinely uncurated position.
    """
    found, _diagnostics = _wild_type_case("2.0")

    assert found == [], (
        "wild-type suppression was disabled by the position spelling; emitted "
        f"{[(hit['primary_accession'], hit['residue_status']) for hit in found]}"
    )


def test_wild_type_suppression_survives_a_case_shifted_subtype_code():
    """meta_data spelling the subtype '1A' must still match a catalogue '1a'.

    Real-world trigger: subtype codes reach meta_data from several typing
    tools and hand curation, and build_subtype_code() concatenates them
    verbatim.  normalize_residue() folds case for residues precisely because
    "a curator typo cannot silently disable a catalog row", but the genotype
    vocabulary gets no such treatment.  Here the exact-subtype tier misses and
    only the leading-digit fallback saves the call - which means the wild type
    used is the union over every subtype of genotype 1, not 1a's own.
    """
    tables = ({("1a", "NS3", "2"): {"K"}}, {})

    residues, tier = AnnotateMutations.lookup_wild_type_residues(tables, "1A", "1", "NS3", 2)

    assert residues == {"K"} and tier == AnnotateMutations.SCOPE_TIER_SUBTYPE, (
        f"subtype '1A' resolved to {residues!r} at tier {tier!r}"
    )


# --------------------------------------------------------------------------
# 3. The genotype machinery is hard-wired to a numeric vocabulary
# --------------------------------------------------------------------------

def test_numeric_genotype_vocabulary_still_matches_at_genotype_level():
    """Control: HCV's own vocabulary works - '6xd' falls back onto genotype 6."""
    assert AnnotateMutations.classify_genotype_scope(frozenset({"6a"}), "6", "6xd") == AnnotateMutations.SCOPE_TIER_GENOTYPE


def test_a_non_numeric_genotype_vocabulary_is_not_globally_out_of_scope():
    """relevant_genotypes is documented as virus-agnostic; it is not.

    Real-world trigger: the column is described as "generic, virus-agnostic"
    and any virus is invited to supply it, but genotype_of_code() implements
    HCV's convention that a subtype code starts with its genotype's digits.
    An influenza catalogue scoped 'H1', a rabies one scoped 'Africa-2' or a
    SARS-CoV-2 one scoped 'BA.2' all yield '' from genotype_of_code(), so the
    genotype tier can never match.  A sequence typed 'H1' with no finer
    subtype, or typed 'H1N1' against a catalogue bucket 'H1', is then declared
    out_of_scope and *every* call for it is suppressed - a silent, build-wide
    false-negative with nothing in the output to distinguish it from a clean
    sequence.
    """
    scope = frozenset({"H1", "H3"})

    assert AnnotateMutations.classify_genotype_scope(scope, "H1", "") != AnnotateMutations.SCOPE_TIER_OUT_OF_SCOPE
    assert AnnotateMutations.classify_genotype_scope(scope, "H1", "H1N1") != AnnotateMutations.SCOPE_TIER_OUT_OF_SCOPE


def test_a_non_numeric_genotype_keeps_a_usable_wild_type_fallback():
    """The genotype-level wild type must survive a non-numeric vocabulary.

    Real-world trigger: the same HCV assumption, one layer down.  Every 'H1',
    'H3' entry is filed under genotype key '' - collapsing the two subtypes'
    wild types into one union - and lookup_wild_type_residues() skips the
    genotype tier whenever the sequence's own genotype is truthy, so that
    union is unreachable.  A sequence whose exact subtype code is absent from
    the catalogue therefore gets no wild type at all and its every catalogued
    residue is emitted as a change.
    """
    catalog = pd.DataFrame(
        [{"protein_name": "HA", "aa_position": "2", "wild_type_residues": "H1:K;H3:R"}]
    )
    tables = AnnotateMutations.build_wild_type_tables(catalog)

    residues, _tier = AnnotateMutations.lookup_wild_type_residues(tables, "H1N1", "H1", "HA", 2)

    assert residues == {"K"}, f"genotype H1 resolved to {residues!r}"


def test_signature_scope_is_the_union_over_the_whole_signature():
    """Characterisation: rule A is per-signature, so scope widens per row.

    Not a defect on its own - build_signature_genotype_scope() documents the
    union - but it is the mechanism by which a genotype-restricted catalogue
    row becomes callable outside its genotype: a row scored only in 1a, filed
    under a signature that also carries genotype-3 rows, passes the gate for a
    genotype-3 sequence.  Pinned here so that a change to per-row scoping is a
    deliberate edit rather than an accident.
    """
    catalog = pd.DataFrame(
        [
            {"signature_id": "sigX", "relevant_genotypes": "1a"},
            {"signature_id": "sigX", "relevant_genotypes": "3"},
        ]
    )

    scope = AnnotateMutations.build_signature_genotype_scope(catalog)

    assert scope == {"sigX": frozenset({"1a", "3"})}
    assert AnnotateMutations.classify_genotype_scope(scope["sigX"], "3", "3a") == AnnotateMutations.SCOPE_TIER_GENOTYPE


# --------------------------------------------------------------------------
# 4. The non-HCV catalogue profiles
# --------------------------------------------------------------------------

def _profile_case(tmp_path, name, header, catalog_row):
    catalog_path = _write_catalog(tmp_path / f"{name}.tsv", [catalog_row], header=header)
    db_path = _build_db(
        tmp_path / f"{name}.db",
        features=(
            "CREATE TABLE features (accession TEXT, product TEXT, cds_start INTEGER, cds_end INTEGER)",
            [("REF1", catalog_row["protein_name"], 1, 9)],
        ),
        alignments=(
            "CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
            "alignment_name TEXT, alignment TEXT)",
            [("REF1", "REF1", "REF1", "ATGAAACCC"), ("SEQ1", "SEQ1", "REF1", "ATGTTTCCC")],
        ),
        meta=(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)",
            [("REF1", "master")],
        ),
    )
    return db_path, catalog_path


def test_hcv_profile_completes_end_to_end(tmp_path, monkeypatch):
    """Control for the two cases below: the HCV profile writes every table."""
    db_path, catalog_path = _profile_case(
        tmp_path, "hcv_profile", CATALOG_HEADER,
        _catalog_row(mutation_id="NS3:2F", alt_residue="F", drug="DrugA"),
    )

    _run_main(db_path, catalog_path, monkeypatch, profile="HCV")

    assert _summary(db_path) == {"SEQ1": ["NS3:2F"]}


def test_influenza_profile_completes_end_to_end(tmp_path, monkeypatch):
    """--catalog_column_profile influenza is an advertised option; it cannot run.

    Real-world trigger: running the pipeline for any virus that is not HCV.
    argparse accepts `influenza` and `all_columns`, annotation succeeds and the
    hit count is printed, then write_mutation_tables() reaches
    build_completed_signatures_only(), whose helper projects the catalogue onto
    a column list containing the two HCV drug-resistance columns.  A catalogue
    without them raises a bare KeyError out of pandas - not an
    AnnotationMappingError, so main()'s handler does not catch it and the exit
    is a traceback with status 1.

    The damage is not just the crash: mutation_catalog has already been
    dropped, rewritten and committed by then, so the database is left carrying
    a fresh catalogue alongside the previous run's summary tables.  A rerun
    that fails this way looks, to any later query, like a successful one.
    """
    header = [column for column in CATALOG_HEADER
              if column not in {"resistance_category", "drug"}] + ["serotypes_tested"]
    db_path, catalog_path = _profile_case(
        tmp_path, "flu_profile", header,
        _catalog_row(mutation_id="HA:2F", protein_name="HA", segment="4",
                     alt_residue="F", serotypes_tested="H1"),
    )

    _run_main(db_path, catalog_path, monkeypatch, profile="influenza", virus="influenza")

    assert _summary(db_path) == {"SEQ1": ["HA:2F"]}


def test_all_columns_profile_completes_without_the_hcv_drug_columns(tmp_path, monkeypatch):
    """`all_columns` is the documented escape hatch for an unmodelled virus.

    It selects every non-private catalogue column, so it is the profile a new
    virus reaches for.  It fails in exactly the same place as `influenza`,
    which means no virus can be onboarded without first inventing empty
    `drug` and `resistance_category` columns.
    """
    header = [column for column in CATALOG_HEADER
              if column not in {"resistance_category", "drug"}]
    db_path, catalog_path = _profile_case(
        tmp_path, "all_columns_profile", header,
        _catalog_row(mutation_id="N:2F", protein_name="N", alt_residue="F"),
    )

    _run_main(db_path, catalog_path, monkeypatch, profile="all_columns")

    assert _summary(db_path) == {"SEQ1": ["N:2F"]}


def test_layout_builders_summarise_a_catalogue_without_the_drug_columns():
    """The failure above, isolated: drug resistance is an HCV-shaped question.

    `drug` and `resistance_category` describe a phenotype a rabies or influenza
    catalogue does not have and cannot invent a value for.  They used to be
    hard-projected, so their absence raised a bare KeyError out of pandas.  They
    are now filled empty instead of demanded, which keeps the output schema the
    same for every virus - a downstream query does not have to know which virus
    produced the table it is reading.
    """
    catalog = pd.DataFrame([{
        "mutation_id": "HA:2F", "signature_id": "sig1", "signature_kind": "single",
        "protein_name": "HA", "segment": "4", "aa_position": "2", "alt_residue": "F",
        "reference_accession": "REF1", "combination_id": "",
    }])
    hits = pd.DataFrame([{
        "primary_accession": "SEQ1", "mutation_id": "HA:2F", "protein_name": "HA",
        "segment": "4", "aa_position": 2, "alt_residue": "F", "combination_id": "",
    }])

    completed = EL.build_completed_signatures_only(hits, catalog)
    assert completed["signature_id"].tolist() == ["sig1"]

    summary = EL.build_signature_summary(hits, catalog)
    assert summary.loc[0, "drug"] == ""
    assert summary.loc[0, "resistance_category"] == ""


# --------------------------------------------------------------------------
# 5. Segments
# --------------------------------------------------------------------------

def test_a_segment_four_catalog_entry_is_not_called_on_a_segment_six_sequence(
    tmp_path, monkeypatch
):
    """Segment is present on both sides of the join and is thrown away.

    Real-world trigger: any segmented virus - which the `influenza` profile
    exists to serve.  meta_data, sequence_alignment, features and the
    catalogue all carry `segment`, but main() projects sequence_alignment down
    to four columns that exclude it, annotate_from_reference_coordinates()
    resolves one master coordinate map for the entire run, and then scans every
    row of the alignment table against it.

    Here an HA (segment 4) catalogue entry is read out of an NA (segment 6)
    sequence's alignment columns.  The emitted row is not merely a false
    positive: it carries segment '4' while the residue came from segment 6, so
    the output claims a provenance it does not have.  Nothing in the
    diagnostics counts a cross-segment read.
    """
    catalog_path = _write_catalog(
        tmp_path / "segmented.tsv",
        [_catalog_row(mutation_id="HA:2F", protein_name="HA", segment="4", alt_residue="F")],
    )
    db_path = _build_db(
        tmp_path / "segmented.db",
        features=(
            "CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, "
            "cds_start INTEGER, cds_end INTEGER)",
            [("REF_HA", "HA", "4", 1, 9), ("REF_NA", "NA", "6", 1, 9)],
        ),
        alignments=(
            "CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
            "alignment_name TEXT, segment TEXT, alignment TEXT)",
            [
                ("REF_HA", "REF_HA", "REF_HA", "4", "ATGAAACCC"),
                ("REF_NA", "REF_NA", "REF_NA", "6", "ATGAAACCC"),
                ("SEQ_NA", "SEQ_NA", "REF_NA", "6", "ATGTTTCCC"),
            ],
        ),
        meta=(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)",
            [("REF_HA", "master", "4"), ("REF_NA", "master", "6"), ("SEQ_NA", "query", "6")],
        ),
    )

    _run_main(db_path, catalog_path, monkeypatch, virus="influenza")

    assert _summary(db_path).get("SEQ_NA", []) == [], (
        "a segment 4 catalogue entry was read out of a segment 6 alignment"
    )


# --------------------------------------------------------------------------
# 6. One malformed annotation row aborts the whole run
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "start_value, why",
    [
        (None, "a NULL coordinate, e.g. a GFF feature whose start failed to parse"),
        ("1.0", "a coordinate written back through pandas as a float string"),
    ],
)
def test_one_unusable_feature_coordinate_does_not_discard_the_others(start_value, why):
    """A bad row must cost one feature, not the entire annotation source.

    Real-world trigger: `gff_features` and `features` are populated by upstream
    stages from GenBank records of wildly varying quality, and their coordinate
    columns are untyped in SQLite - a NULL start or a REAL written as '1.0'
    both reach make_feature_entry() as-is.  int() then raises ValueError, which
    nothing catches: load_db_gff_feature_maps() has no per-row handler and
    main() only handles AnnotationMappingError.  The whole build dies on a
    traceback, and the message names neither the table nor the accession, so
    the operator has one Python exception and several hundred thousand feature
    rows to search.  Every well-formed feature in the same table is lost with
    it.
    """
    conn = sqlite3.connect(":memory:")
    try:
        conn.execute(
            "CREATE TABLE gff_features (reference_accession TEXT, feature_type TEXT, "
            "product TEXT, start, end)"
        )
        conn.executemany(
            "INSERT INTO gff_features VALUES (?, ?, ?, ?, ?)",
            [("REF1", "gene", "NS3", start_value, 9), ("REF1", "gene", "NS5A", 10, 12)],
        )
        conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
        conn.execute("INSERT INTO meta_data VALUES ('REF1', 'master')")
        conn.commit()

        feature_maps = AnnotateMutations.load_db_gff_feature_maps(conn, {})
    finally:
        conn.close()

    assert "NS5A" in feature_maps.get("REF1", {}), f"the good feature was lost to {why}"


@pytest.mark.parametrize("aa_position", ["inf", "1e400", "-inf"])
def test_an_infinite_catalog_position_is_counted_invalid_rather_than_raising(aa_position):
    """The invalid-position counter exists precisely so a bad row is survivable.

    Real-world trigger: the same spreadsheet round-trip that turns '107' into
    '107.0' turns an overlong numeric token into '1e400', which float() reads
    as infinity.  prepare_catalog() already tolerates 'abc' and 'nan' by
    tallying them under invalid_positions and carrying on; an infinite value
    takes the identical path but escapes through the one exception type the
    handler does not name, so a single catalogue cell aborts the run before any
    sequence is annotated.
    """
    prepared, invalid = AnnotateMutations.prepare_catalog(
        pd.DataFrame([_catalog_row(aa_position=aa_position)]), {}
    )

    assert invalid == 1 and prepared["_aa_position_int"].tolist() == [None]


# --------------------------------------------------------------------------
# 7. Product-name canonicalisation crosses virus boundaries
# --------------------------------------------------------------------------

def test_hcv_product_heuristics_resolve_hcv_products():
    """Control: the compact-name inference still does its job, once asked for HCV."""
    hcv = AnnotateMutations.virus_profile("HCV")
    assert AnnotateMutations.canonicalize_product("protease/helicase protein NS3", {}, hcv) == "NS3"
    assert AnnotateMutations.canonicalize_product("nonstructural protein NS5A", {}, hcv) == "NS5A"

    # The polyprotein collapse is about the SHAPE of the annotation, not about
    # any one virus's proteins, so it applies without a virus profile.
    assert AnnotateMutations.canonicalize_product("polyprotein", {}) == "Whole genome"


@pytest.mark.parametrize(
    "product",
    ["nucleocapsid core protein", "core-like protein", "envelope protein E1"],
)
def test_hcv_product_names_are_not_forced_onto_another_virus(product):
    """A non-HCV product must not be renamed into HCV's protein vocabulary.

    Real-world trigger: this repository builds rabies alongside HCV, and
    canonicalize_product() is applied to both the catalogue's protein_name and
    every feature product with no virus context - `--virus` is declared,
    documented as "Virus context for specific logics (e.g. HCV)", and then
    never referenced again.  A GenBank product description containing the word
    'core', or an 'E1'/'E2' token, is rewritten to the HCV mature peptide of
    that name.  The two sides of the join then disagree - a catalogue saying
    'N' cannot meet a feature renamed 'Core' - and the rows silently drop out
    under proteins_missing_from_reference, which is printed as a warning and
    never counted in the run's diagnostics.
    """
    assert AnnotateMutations.canonicalize_product(product, {}) == product


# --------------------------------------------------------------------------
# 8. What the run summary does not say
# --------------------------------------------------------------------------

def test_catalog_rows_dropped_for_a_missing_protein_are_counted_in_the_summary():
    """The mapping summary must account for every catalogue row it discarded.

    Real-world trigger: a gene-alias change, a reference whose mature peptides
    were never annotated, or the product-renaming above leaves some proteins
    unresolvable.  Only a total failure raises - mappable_catalog_rows == 0 -
    so losing 9 catalogue proteins out of 10 is not an error.  The count lives
    in proteins_missing_from_reference, which is printed as a warning line and
    dropped on the floor, while '[AnnotateMutations] Mapping summary: ...' -
    the line an operator or a downstream check actually parses - reports a
    healthy-looking run.
    """
    _found, diagnostics = _annotate(
        [
            _catalog_row(protein_name="NS4B", mutation_id="NS4B:2K"),
            _catalog_row(protein_name="NS3", mutation_id="NS3:2K"),
        ],
        [("REF1", "ATGAAACCCGGG")],
        NS3_FEATURE_MAP,
    )

    assert any("protein" in key for key in diagnostics), (
        f"NS4B was dropped but the summary only reports {sorted(diagnostics)}"
    )


def test_blank_alignment_name_warning_reports_how_many_rows_were_filled():
    """A warning that overstates what it did is worse than no warning.

    Real-world trigger: any build where some rows reach sequence_alignment
    without an alignment_name.  The message is the operator's only signal that
    the coordinate reference was guessed rather than declared, and it reads
    'Filled 2 blank alignment_name values' for a table with one blank row and
    one that already named the master - so the number cannot be used to judge
    how much of the run rests on the guess.
    """
    seq_aln = pd.DataFrame(
        [
            {"sequence_id": "REF1", "primary_accession": "REF1",
             "alignment": "ATGAAACCCGGG", "alignment_name": "REF1"},
            {"sequence_id": "SEQ1", "primary_accession": "SEQ1",
             "alignment": "ATGAAACCCGGG", "alignment_name": ""},
        ]
    )
    prepared, _invalid = AnnotateMutations.prepare_catalog(pd.DataFrame([_catalog_row()]), {})
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        AnnotateMutations.annotate_from_reference_coordinates(
            prepared, seq_aln,
            pd.DataFrame([{"primary_accession": "REF1", "accession_type": "master"}]),
            {}, {"REF1": NS3_FEATURE_MAP}, False,
        )

    filled = [line for line in buffer.getvalue().splitlines() if "Filled" in line]
    assert filled and "Filled 1 " in filled[0], filled


def test_conflicting_genotypes_for_one_accession_are_not_resolved_by_row_order():
    """Two genotype calls for one accession is a data problem, not a coin toss.

    Real-world trigger: meta_data carries one row per accession *per segment*
    in a segmented build, and cross-run merges can leave an accession with a
    stale row beside its current one.  The genotype that survives is whichever
    pandas yields last, and it then drives both rule A and rule B for that
    sequence - so an identical database read in a different row order produces
    different resistance calls, with nothing recorded to say a conflict existed.
    """
    meta_data = pd.DataFrame(
        [
            {"primary_accession": "SEQ1", "accession_type": "query",
             "genotype": "1", "subtype": "a"},
            {"primary_accession": "SEQ1", "accession_type": "query",
             "genotype": "3", "subtype": "a"},
        ]
    )

    genotype_map = AnnotateMutations.build_sequence_genotype_map(meta_data)

    assert genotype_map.get("SEQ1") != ("3", "3a"), (
        "the second row silently overwrote the first"
    )


# --------------------------------------------------------------------------
# 9. Translation and codon assembly
# --------------------------------------------------------------------------

def test_translate_codon_does_not_translate_a_gapped_four_column_codon():
    """Removing a gap and translating what is left invents a residue.

    Real-world trigger: translate_codon() is the module's public translation
    entry point and is called with whatever a caller assembled.  Its guard
    rejects a *short* codon but runs `.replace('-', '')` first, so 'AT-G' -
    three bases spread over four alignment columns, which is what a codon
    spanning an insertion looks like - is silently read as ATG/Met.  The
    annotator's own path is protected by extract_aligned_codon()'s len == 3
    check, so this is a latent hazard rather than a live miscall today; the
    guard belongs in the function that claims to have one.
    """
    assert AnnotateMutations.translate_codon("AT-G") == "X"
    assert AnnotateMutations.translate_codon("A-T-G") == "X"


def test_translate_codon_rejects_a_short_or_unknown_codon():
    """Control: the cases the length guard and the table lookup do cover."""
    assert AnnotateMutations.translate_codon("AT") == "X"
    assert AnnotateMutations.translate_codon("---") == "X"
    assert AnnotateMutations.translate_codon("NNN") == "X"
    assert AnnotateMutations.translate_codon("atg") == "M"


# --------------------------------------------------------------------------
# 10. Degenerate inputs that are handled correctly - regression pins
# --------------------------------------------------------------------------

def test_an_empty_sequence_alignment_exits_with_the_mapping_error_code(tmp_path, monkeypatch):
    """No alignments must be a clean exit 2, not a traceback."""
    catalog_path = _write_catalog(tmp_path / "empty_aln.tsv", [_catalog_row()])
    db_path = _build_db(
        tmp_path / "empty_aln.db",
        features=("CREATE TABLE features (accession TEXT, product TEXT, cds_start INTEGER, "
                  "cds_end INTEGER)", [("REF1", "NS3", 1, 9)]),
        alignments=("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
                    "alignment_name TEXT, alignment TEXT)", []),
        meta=("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)",
              [("REF1", "master")]),
    )

    with pytest.raises(SystemExit) as excinfo:
        _run_main(db_path, catalog_path, monkeypatch)

    assert excinfo.value.code == 2


def test_a_header_only_catalog_exits_with_the_mapping_error_code(tmp_path, monkeypatch):
    """An empty catalogue must fail as a mapping error, not silently succeed."""
    catalog_path = _write_catalog(tmp_path / "empty_cat.tsv", [])
    db_path = _build_db(
        tmp_path / "empty_cat.db",
        features=("CREATE TABLE features (accession TEXT, product TEXT, cds_start INTEGER, "
                  "cds_end INTEGER)", [("REF1", "NS3", 1, 9)]),
        alignments=("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
                    "alignment_name TEXT, alignment TEXT)",
                    [("REF1", "REF1", "REF1", "ATGAAACCC")]),
        meta=("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)",
              [("REF1", "master")]),
    )

    with pytest.raises(SystemExit) as excinfo:
        _run_main(db_path, catalog_path, monkeypatch)

    assert excinfo.value.code == 2


def test_an_alignment_shorter_than_the_master_is_counted_not_crashed():
    """A ragged alignment row must be tallied out of bounds, never guessed at.

    The master is 12 columns wide and the truncated row only 6, so the codon
    for NS3:3 has no columns to read in it.  The alternative - indexing off the
    end and translating whatever comes back - is what the bounds check exists
    to prevent, and the tally is what tells an operator it happened.
    """
    found, diagnostics = _annotate(
        [_catalog_row(aa_position="3", alt_residue="W", mutation_id="NS3:3W")],
        [("REF1", "ATGAAACCCGGG"), ("TRUNCATED", "ATGAAA")],
        NS3_FEATURE_MAP,
    )

    assert found == []
    assert diagnostics.get("codon_out_of_bounds") == 1


def test_a_null_alignment_on_a_query_row_is_skipped_not_called():
    """A NULL alignment carries no evidence and must produce no call."""
    found, _diagnostics = _annotate(
        [_catalog_row(alt_residue="K")],
        [("REF1", "ATGAAACCCGGG"), ("NULL_ALN", None)],
        NS3_FEATURE_MAP,
    )

    assert [hit["primary_accession"] for hit in found] == ["REF1"], (
        "only the reference, which genuinely reads K, should be called"
    )


def test_features_rows_are_kept_regardless_of_their_reference_accession_columns():
    """Characterisation: the reference-agreement guard in the features branch is a no-op.

    load_db_gff_feature_maps() tests
    `reference_accession not in {ref_candidate, master_candidate, reference_accession}`,
    where the third member is the value being tested - so the set always
    contains it and the branch can never be taken.  Whatever the guard was
    meant to exclude is currently included, and only the allowed-accession
    filter above it does any filtering.  Pinned so that repairing the guard is
    a deliberate change with a visible blast radius rather than a silent one.
    """
    conn = sqlite3.connect(":memory:")
    try:
        conn.execute(
            "CREATE TABLE features (accession TEXT, master_ref_accession TEXT, "
            "reference_accession TEXT, cds_start INTEGER, cds_end INTEGER, product TEXT)"
        )
        conn.executemany(
            "INSERT INTO features VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("REF1", "REF1", "REF1", 1, 9, "NS3"),
                # accession agrees with neither of the two reference columns.
                ("GENREF", "REF1", "REF1", 1, 9, "NS5A"),
            ],
        )
        conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)")
        conn.execute("INSERT INTO meta_data VALUES ('REF1', 'master')")
        conn.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT)"
        )
        conn.execute("INSERT INTO sequence_alignment VALUES ('SEQ1', 'GENREF')")
        conn.commit()

        feature_maps = AnnotateMutations.load_db_gff_feature_maps(conn, {})
    finally:
        conn.close()

    assert "NS5A" in feature_maps.get("GENREF", {})


# --------------------------------------------------------------------------
# 11. The virus profile: --virus is read, and drives behaviour
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "virus, expected",
    [
        ("HCV", "hcv"), ("hcv", "hcv"), ("Hepatitis C virus", "hcv"),
        ("influenza", "influenza"), ("IAV", "influenza"), ("flu", "influenza"),
        ("all_columns", "generic"),
        ("rabies", "rabv"), ("RABV", "rabv"), ("torono", "generic"),
        ("", "generic"), (None, "generic"),
    ],
)
def test_virus_names_resolve_to_a_profile(virus, expected):
    """An unknown virus resolves to the generic profile, it does not fail.

    `--catalog_column_profile` used to be required with
    choices=['HCV','influenza','all_columns'], and BOTH pipeline call sites
    derive that argument from the virus name.  A build for any other virus -
    rabies, which this repository also produces - therefore died at argument
    parsing before a single line of the annotator ran.
    """
    assert AnnotateMutations.resolve_virus_name(virus) == expected


def test_virus_falls_back_to_the_column_profile_and_then_to_generic():
    """The two arguments carry the same information, so either may stand in.

    vgtk-init.nf and vgtk-rabv.sh both set the column profile FROM the virus
    name, so a caller that passes only one of them still gets the right
    vocabulary.
    """
    assert AnnotateMutations.resolve_virus_name("", "HCV") == "hcv"
    assert AnnotateMutations.resolve_virus_name("influenza", "all_columns") == "influenza"
    assert AnnotateMutations.resolve_virus_name("", "") == AnnotateMutations.GENERIC_VIRUS


@pytest.mark.parametrize(
    "profile_name, expected",
    [
        ("HCV", "hcv columns"),
        ("influenza", ["serotypes_tested"]),
        ("all_columns", None),
        ("generic", None),
        ("rabies", None),
        ("", None),
    ],
)
def test_unspecified_or_unknown_virus_keeps_every_catalog_column(profile_name, expected):
    """The default column set is 'whatever the curator supplied'.

    Keeping every column is the only default that cannot lose data, and it
    means onboarding a virus does not start by inventing empty `drug` and
    `resistance_category` columns for a pathogen that has no drugs.
    """
    columns = AnnotateMutations.catalog_columns_for_profile(profile_name)
    if expected == "hcv columns":
        assert "drug" in columns and "resistance_category" in columns
    else:
        assert columns == expected


def test_hcv_and_influenza_vocabularies_do_not_leak_into_each_other():
    """Each virus renames only into its own protein vocabulary.

    'core', 'E1' and 'E2' are ordinary words in non-HCV product descriptions,
    and 'NA'/'PA'/'M1' are ordinary tokens outside influenza.  Both sides of
    the protein join run through canonicalize_product(), so a rename that fits
    only one side makes the row vanish rather than mismatch loudly.
    """
    hcv = AnnotateMutations.virus_profile("HCV")
    flu = AnnotateMutations.virus_profile("influenza")
    generic = AnnotateMutations.virus_profile("rabies")

    assert AnnotateMutations.canonicalize_product("nonstructural protein NS5A", {}, hcv) == "NS5A"
    assert AnnotateMutations.canonicalize_product("nucleocapsid core protein", {}, hcv) == "Core"

    # The same description under another virus keeps its own name.
    assert AnnotateMutations.canonicalize_product("nucleocapsid core protein", {}, flu) == \
        "nucleocapsid core protein"
    assert AnnotateMutations.canonicalize_product("nucleocapsid core protein", {}, generic) == \
        "nucleocapsid core protein"

    # Influenza resolves its own polymerase products; HCV must not touch them.
    assert AnnotateMutations.canonicalize_product("polymerase PB2", {}, flu) == "PB2"
    assert AnnotateMutations.canonicalize_product("polymerase PB2", {}, hcv) == "polymerase PB2"
    # PA-X is a distinct protein from PA and must win the longer match.
    assert AnnotateMutations.canonicalize_product("PA-X protein", {}, flu) == "PA-X"


def test_gene_alias_lookup_is_a_plain_mapping_carrying_its_profile():
    """The profile rides on the alias lookup, which is threaded everywhere already.

    Every existing call site reads the lookup as a dict, so it stays one; the
    alternative was adding a virus parameter to the six functions between
    main() and canonicalize_product().
    """
    lookup = AnnotateMutations.GeneAliasLookup({"ns5a": "NS5A"}, profile=AnnotateMutations.virus_profile("HCV"))

    assert lookup["ns5a"] == "NS5A"
    assert lookup.get("missing", "fallback") == "fallback"
    assert AnnotateMutations.profile_of_lookup(lookup) is AnnotateMutations.virus_profile("HCV")
    # A bare dict yields the generic profile rather than raising.
    assert AnnotateMutations.profile_of_lookup({}) == AnnotateMutations.VIRUS_PROFILES[AnnotateMutations.GENERIC_VIRUS]


# --------------------------------------------------------------------------
# 12. 'NA' is neuraminidase, not a null
# --------------------------------------------------------------------------

def test_a_catalog_protein_named_na_survives_being_read(tmp_path, monkeypatch):
    """NA is influenza's neuraminidase and one of pandas' default null tokens.

    Real-world trigger: any influenza catalogue.  `pd.read_csv(dtype=str)`
    without keep_default_na=False turns the literal 'NA' into NaN, and
    canonicalize_product() then renders it back as the string 'nan' - so every
    neuraminidase row, which is what oseltamivir resistance is catalogued
    against, silently failed to resolve against the NA feature.  Found by
    running the annotator against the shipped 8-segment influenza build, where
    segment 6 produced zero calls and the warning line read 'Missing protein
    annotations in master feature map: nan'.

    This repository already knows the hazard: see test_na_is_not_a_null.py and
    the keep_default_na=False call sites in gisaid_tidy.py and
    merge_into_gB_matrix.py.
    """
    catalog_path = _write_catalog(
        tmp_path / "na.tsv",
        [_catalog_row(mutation_id="NA:2F", protein_name="NA", segment="6", alt_residue="F")],
    )
    db_path = _build_db(
        tmp_path / "na.db",
        features=("CREATE TABLE features (accession TEXT, product TEXT, cds_start INTEGER, "
                  "cds_end INTEGER)", [("REF1", "NA", 1, 9)]),
        alignments=("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
                    "alignment_name TEXT, alignment TEXT)",
                    [("REF1", "REF1", "REF1", "ATGAAACCC"), ("SEQ1", "SEQ1", "REF1", "ATGTTTCCC")]),
        meta=("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)",
              [("REF1", "master")]),
    )

    _run_main(db_path, catalog_path, monkeypatch, profile="influenza", virus="influenza")

    assert _summary(db_path) == {"SEQ1": ["NA:2F"]}


@pytest.mark.parametrize("protein", ["NA", "NULL", "None", "N/A", "nan"])
def test_catalog_cells_that_look_like_null_tokens_are_kept(protein, tmp_path):
    """A curated catalogue is entitled to use any of pandas' sentinels as a value."""
    path = _write_catalog(tmp_path / "sentinels.tsv", [_catalog_row(protein_name=protein)])
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    assert frame.loc[0, "protein_name"] == protein


# --------------------------------------------------------------------------
# 13. Segment-aware annotation, positively
# --------------------------------------------------------------------------

def _segmented_db(tmp_path, name="multiseg"):
    """Two segments of different widths, each with its own master and protein."""
    return _build_db(
        tmp_path / f"{name}.db",
        features=(
            "CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, "
            "cds_start INTEGER, cds_end INTEGER)",
            [("REF_HA", "HA", "4", 1, 9), ("REF_NA", "NA", "6", 1, 12)],
        ),
        alignments=(
            "CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
            "alignment_name TEXT, segment TEXT, alignment TEXT)",
            [
                ("REF_HA", "REF_HA", "REF_HA", "4", "ATGAAACCC"),
                ("SEQ_HA", "SEQ_HA", "REF_HA", "4", "ATGTTTCCC"),
                ("REF_NA", "REF_NA", "REF_NA", "6", "ATGAAACCCGGG"),
                ("SEQ_NA", "SEQ_NA", "REF_NA", "6", "ATGAAATGGGGG"),
            ],
        ),
        meta=(
            "CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)",
            [("REF_HA", "master", "4"), ("REF_NA", "master", "6"),
             ("SEQ_HA", "query", "4"), ("SEQ_NA", "query", "6")],
        ),
    )


def test_each_segment_is_annotated_against_its_own_reference(tmp_path, monkeypatch):
    """Both segments annotate, each from its own master's column space.

    The segments here are deliberately different widths, which is what a real
    influenza build looks like: on the shipped 8-segment database the alignment
    widths run from 844 to 2283 columns.  One coordinate space for all of them
    cannot be right for more than one.
    """
    catalog_path = _write_catalog(
        tmp_path / "multiseg.tsv",
        [
            _catalog_row(mutation_id="HA:2F", protein_name="HA", segment="4", alt_residue="F"),
            _catalog_row(mutation_id="NA:3W", protein_name="NA", segment="6",
                         aa_position="3", alt_residue="W"),
        ],
    )

    _run_main(_segmented_db(tmp_path), catalog_path, monkeypatch, virus="influenza")

    assert _summary(tmp_path / "multiseg.db") == {"SEQ_HA": ["HA:2F"], "SEQ_NA": ["NA:3W"]}


def test_a_catalog_segment_no_sequence_carries_is_reported(tmp_path, monkeypatch):
    """A catalogue naming a segment the build does not have must say so.

    Silently annotating nothing for it is indistinguishable from annotating it
    and finding nothing, which is the difference between a missing input and a
    clean result.
    """
    catalog_path = _write_catalog(
        tmp_path / "stranded.tsv",
        [
            _catalog_row(mutation_id="HA:2F", protein_name="HA", segment="4", alt_residue="F"),
            _catalog_row(mutation_id="PB2:2F", protein_name="PB2", segment="1", alt_residue="F"),
        ],
    )

    output = _run_main(
        _segmented_db(tmp_path, "stranded"), catalog_path, monkeypatch, virus="influenza"
    )

    assert "no sequence in this build carries" in output
    assert "1" in output  # the count of stranded catalog rows
    assert _summary(tmp_path / "stranded.db") == {"SEQ_HA": ["HA:2F"]}


def test_a_single_segment_build_is_not_treated_as_segmented(tmp_path, monkeypatch):
    """One segment means the catalogue's own segment column is not consulted.

    This is what keeps every non-segmented virus, and every caller that stores
    no segment at all, on exactly the coordinate resolution it has always had.
    The catalogue row below names segment '1' while the alignment says '1.0' -
    a spelling difference the shipped HCV database really does contain - and it
    must still be annotated.
    """
    catalog_path = _write_catalog(
        tmp_path / "single.tsv", [_catalog_row(mutation_id="NS3:2F", alt_residue="F")]
    )
    db_path = _build_db(
        tmp_path / "single.db",
        features=("CREATE TABLE features (accession TEXT, product TEXT, cds_start INTEGER, "
                  "cds_end INTEGER)", [("REF1", "NS3", 1, 9)]),
        alignments=("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
                    "alignment_name TEXT, segment TEXT, alignment TEXT)",
                    [("REF1", "REF1", "REF1", "1.0", "ATGAAACCC"),
                     ("SEQ1", "SEQ1", "REF1", "1.0", "ATGTTTCCC")]),
        meta=("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT)",
              [("REF1", "master")]),
    )

    _run_main(db_path, catalog_path, monkeypatch, profile="HCV", virus="HCV")

    assert _summary(db_path) == {"SEQ1": ["NS3:2F"]}


def test_segment_values_are_compared_through_normalize_segment(tmp_path, monkeypatch):
    """'4' and '4.0' are the same segment.

    Verified on a shipped build: HCV_full stores sequence_alignment.segment as
    both '1' (135544 rows) and '1.0' (2506 rows).  Comparing the raw text would
    split one segment into two buckets, give each its own coordinate reference,
    and silently halve the sequences any catalogue row could reach.
    """
    catalog_path = _write_catalog(
        tmp_path / "drift.tsv",
        [_catalog_row(mutation_id="HA:2F", protein_name="HA", segment="4.0", alt_residue="F")],
    )
    db_path = _build_db(
        tmp_path / "drift.db",
        features=("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, "
                  "cds_start INTEGER, cds_end INTEGER)",
                  [("REF_HA", "HA", "4", 1, 9), ("REF_NA", "NA", "6", 1, 9)]),
        alignments=("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
                    "alignment_name TEXT, segment TEXT, alignment TEXT)",
                    [("REF_HA", "REF_HA", "REF_HA", "4", "ATGAAACCC"),
                     ("SEQ_HA", "SEQ_HA", "REF_HA", "4.0", "ATGTTTCCC"),
                     ("REF_NA", "REF_NA", "REF_NA", "6", "ATGAAACCC")]),
        meta=("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, segment TEXT)",
              [("REF_HA", "master", "4"), ("REF_NA", "master", "6")]),
    )

    _run_main(db_path, catalog_path, monkeypatch, virus="influenza")

    assert _summary(db_path) == {"SEQ_HA": ["HA:2F"]}


# --------------------------------------------------------------------------
# 14. The annotator and the verifier cannot drift apart
# --------------------------------------------------------------------------

def test_verify_mutations_canonicalises_exactly_as_the_annotator_does():
    """One protein vocabulary, or the verification step reports its own artefacts.

    VerifyMutations re-derives what the annotator claimed, so it has to
    canonicalise product names by identical rules.  It used to carry verbatim
    private copies of infer_compact_product_name() and canonicalize_product()
    with HCV's five patterns hard-coded.  The moment the annotator learned
    which virus it was annotating, those copies would have gone on applying
    HCV's vocabulary to every virus - and the verifier would have reported
    mismatches caused by nothing but its own rules.  It imports them now, so
    this cannot regress without deleting the import.
    """
    import VerifyMutations as VM

    assert VM.canonicalize_product.__module__ == VM.__name__
    for product in ["nonstructural protein NS5A", "envelope protein E2",
                    "p7 protein", "Core protein", "polyprotein", "haemagglutinin"]:
        for virus in ("HCV", "influenza", "rabies"):
            profile = AnnotateMutations.virus_profile(virus)
            assert VM.canonicalize_product(product, {}, profile) == \
                AnnotateMutations.canonicalize_product(product, {}, profile), (
                    f"{product!r} under {virus} diverges between the two scripts"
                )


def test_a_catalog_row_with_no_segment_is_counted_when_it_fires(tmp_path, monkeypatch):
    """A wildcard row in a segmented build is a call resting on an assumption.

    It is tried against every segment because there is nothing better to do
    with it, and the call it produces may well be right - but which segment it
    came from is this code's choice rather than the catalogue's, so the run
    says how many such calls it made.
    """
    catalog_path = _write_catalog(
        tmp_path / "wildcard.tsv",
        [_catalog_row(mutation_id="HA:2F", protein_name="HA", segment="", alt_residue="F")],
    )

    output = _run_main(
        _segmented_db(tmp_path, "wildcard"), catalog_path, monkeypatch, virus="influenza"
    )

    assert "emitted_from_segment_unstated_catalog_row=1" in output
    assert _summary(tmp_path / "wildcard.db") == {"SEQ_HA": ["HA:2F"]}
