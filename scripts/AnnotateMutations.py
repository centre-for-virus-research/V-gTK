#!/usr/bin/env python3

import sqlite3
import pandas as pd
import argparse
import os
import sys

# Standard genetic code dictionary
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def translate_codon(codon):
    codon = codon.upper().replace('-', '')
    if len(codon) < 3:
        return 'X'
    return CODON_TABLE.get(codon, 'X')

def main():
    parser = argparse.ArgumentParser(description="Annotate mutations and drug resistance.")
    parser.add_argument("--db", required=True, help="Path to SQLite database.")
    parser.add_argument("--mutation_catalog", required=True, help="Path to mutation catalog TSV.")
    parser.add_argument("--virus", default="", help="Virus context for specific logics (e.g. HCV)")
    args = parser.parse_args()

    if not os.path.isfile(args.db):
        print(f"Error: Database {args.db} not found.", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.mutation_catalog):
        print(f"Error: Catalog {args.mutation_catalog} not found.", file=sys.stderr)
        sys.exit(1)

    print("Loading mutation catalog...")
    catalog = pd.read_csv(args.mutation_catalog, sep='\t', dtype=str)
    
    required_cols = ['protein_name', 'segment', 'aa_position', 'alt_residue', 
                     'reference_accession', 'mutation_id', 'mutation_type', 
                     'signature_id', 'signature_kind']
    for col in required_cols:
        if col not in catalog.columns:
            print(f"Error: Catalog missing required column '{col}'.", file=sys.stderr)
            sys.exit(1)

    conn = sqlite3.connect(args.db)
    
    print("Loading features and sequence alignments...")
    features = pd.read_sql_query("SELECT * FROM features", conn)
    seq_aln = pd.read_sql_query("SELECT sequence_id, alignment FROM sequence_alignment", conn)
    
    aln_dict = dict(zip(seq_aln['sequence_id'], seq_aln['alignment']))
        
    has_segment = 'segment' in features.columns
    
    mutations_found = []
    
    print("Extracting mutations...")
    for idx, row in catalog.iterrows():
        protein_name = row['protein_name']
        segment_val = row.get('segment', '')
        try:
            aa_pos = int(float(row['aa_position']))
        except ValueError:
            continue
            
        alt_residue = row['alt_residue']
        mutation_id = row['mutation_id']
        
        # Match feature
        if has_segment and pd.notna(segment_val) and str(segment_val).strip() != '':
            feat_match = features[(features['product'] == protein_name) & (features['segment'] == str(segment_val))]
        else:
            feat_match = features[features['product'] == protein_name]
            
        if feat_match.empty:
            continue
            
        # Group by accession
        for _, feat in feat_match.iterrows():
            acc = feat['accession']
            if acc not in aln_dict:
                continue
                
            cds_start = int(feat['cds_start'])
            # 1-based aa position to 0-based nucleotide index
            nuc_idx = (cds_start - 1) + (aa_pos - 1) * 3
            
            aln = aln_dict[acc]
            if nuc_idx + 3 <= len(aln):
                codon = aln[nuc_idx:nuc_idx+3]
                aa = translate_codon(codon)
                
                if aa == alt_residue:
                    mutations_found.append({
                        'primary_accession': acc,
                        'mutation_id': mutation_id,
                        'protein_name': protein_name,
                        'segment': feat['segment'] if has_segment else '',
                        'aa_position': aa_pos,
                        'alt_residue': alt_residue,
                        'combination_id': row.get('combination_id', '')
                    })

    df_mut = pd.DataFrame(mutations_found)
    cursor = conn.cursor()
    
    print("Writing sequence_mutations table...")
    cursor.execute("DROP TABLE IF EXISTS sequence_mutations")
    if not df_mut.empty:
        df_mut.to_sql('sequence_mutations', conn, if_exists='replace', index=False)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_seq_mut_acc_mutid ON sequence_mutations(primary_accession, mutation_id)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_seq_mut_seg_prot ON sequence_mutations(segment, protein_name)")
    else:
        cursor.execute("CREATE TABLE sequence_mutations (primary_accession TEXT, mutation_id TEXT, protein_name TEXT, segment TEXT, aa_position INTEGER, alt_residue TEXT, combination_id TEXT)")
        
    # HCV combination logic
    if args.virus.upper() == 'HCV' and 'combination_id' in catalog.columns:
        print("Evaluating HCV combination logics...")
        comb_found = []
        if not df_mut.empty:
            comb_catalog = catalog[catalog['signature_kind'] == 'combination'].dropna(subset=['combination_id'])
            
            # Map combinations to sizes
            comb_req = {}
            for _, row in comb_catalog.iterrows():
                cid = row['combination_id']
                if pd.notna(cid) and str(cid).strip() != '':
                    comb_req[cid] = {
                        'size': int(float(row['combination_size'])) if pd.notna(row['combination_size']) else 1,
                        'drug': row.get('drug', ''),
                        'resistance_category': row.get('resistance_category', '')
                    }
                    
            if comb_req:
                # Group by accession and combination_id
                df_comb = df_mut[df_mut['combination_id'] != '']
                if not df_comb.empty:
                    grouped = df_comb.groupby(['primary_accession', 'combination_id']).size().reset_index(name='count')
                    for _, row in grouped.iterrows():
                        acc = row['primary_accession']
                        cid = row['combination_id']
                        count = row['count']
                        
                        if cid in comb_req:
                            req_size = comb_req[cid]['size']
                            status = 'complete' if count >= req_size else 'partial'
                            
                            comb_found.append({
                                'primary_accession': acc,
                                'combination_id': cid,
                                'combination_status': status,
                                'mutations_detected': count,
                                'mutations_required': req_size,
                                'resistance_category': comb_req[cid]['resistance_category'],
                                'drug': comb_req[cid]['drug']
                            })
                            
        df_drug = pd.DataFrame(comb_found)
        cursor.execute("DROP TABLE IF EXISTS sequence_drug_resistance")
        if not df_drug.empty:
            df_drug.to_sql('sequence_drug_resistance', conn, if_exists='replace', index=False)
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_seq_drug_acc_comb ON sequence_drug_resistance(primary_accession, combination_id)")
        else:
            cursor.execute("CREATE TABLE sequence_drug_resistance (primary_accession TEXT, combination_id TEXT, combination_status TEXT, mutations_detected INTEGER, mutations_required INTEGER, resistance_category TEXT, drug TEXT)")
            
    conn.commit()
    conn.close()
    print("Done.")

if __name__ == "__main__":
    main()
