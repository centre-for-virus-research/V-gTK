#!/usr/bin/env python3

import sqlite3
import pandas as pd
import argparse
import os
import re
import sys
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

# The verifier's whole job is to re-derive what the annotator claimed, so it
# must canonicalise protein names by exactly the same rules. It used to carry
# verbatim private copies of these two functions with HCV's patterns hard-coded
# into them; once AnnotateMutations learned which virus it was annotating, those
# copies would have gone on applying HCV's vocabulary to every virus, and the
# verifier would have reported mismatches that were artefacts of its own rules.
# Importing removes the possibility by construction.
from AnnotateMutations import (
    canonicalize_product as _canonicalize_product,
    infer_compact_product_name as _infer_compact_product_name,
    virus_profile,
)

#: Which virus's protein vocabulary to canonicalise with. HCV is the default
#: because every other default in this script's argument parser is an HCV path;
#: main() replaces it from --virus.
ACTIVE_VIRUS_PROFILE = virus_profile('HCV')

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

def normalize_lookup_key(value):
    return re.sub(r'[^a-z0-9]+', '', str(value or '').strip().lower())

def infer_compact_product_name(product_name, profile=None):
    """Delegates to AnnotateMutations so the two scripts cannot drift apart."""
    return _infer_compact_product_name(product_name, profile or ACTIVE_VIRUS_PROFILE)


def canonicalize_product(product_name, alias_lookup, profile=None):
    """Delegates to AnnotateMutations so the two scripts cannot drift apart."""
    return _canonicalize_product(product_name, alias_lookup, profile or ACTIVE_VIRUS_PROFILE)


def load_gene_alias_lookup(conn):
    alias_lookup = {}
    try:
        genes = pd.read_sql_query('SELECT * FROM genes', conn)
    except Exception:
        return alias_lookup

    for _, row in genes.iterrows():
        canonical = ''
        for candidate in (row.get('display_name'), row.get('name'), row.get('description')):
            candidate = str(candidate or '').strip()
            if candidate:
                canonical = candidate
                break
        if not canonical:
            continue
        for candidate in (row.get('description'), row.get('display_name'), row.get('name')):
            candidate = str(candidate or '').strip()
            if candidate:
                alias_lookup[normalize_lookup_key(candidate)] = canonical
    return alias_lookup

def extract_accession_tokens(value):
    text = str(value or '').strip()
    if not text:
        return []
    seen = []
    for token in re.findall(r'[A-Z]{1,4}_?\d+(?:\.\d+)?', text.upper()):
        if token not in seen:
            seen.append(token)
        token_no_version = token.split('.')[0]
        if token_no_version not in seen:
            seen.append(token_no_version)
    if text not in seen:
        seen.insert(0, text)
    return seen

def build_alignment_coordinate_map(alignment):
    coord_map = {}
    ungapped_position = 0
    for alignment_index, base in enumerate(str(alignment or '')):
        if base != '-':
            ungapped_position += 1
            coord_map[ungapped_position] = alignment_index
    return coord_map

def resolve_aligned_codon_indices(coord_map, cds_start, aa_pos):
    try:
        cds_start = int(cds_start)
        aa_pos = int(aa_pos)
    except (TypeError, ValueError):
        return None

    if cds_start < 1 or aa_pos < 1:
        return None

    nuc_start = cds_start + (aa_pos - 1) * 3
    query_positions = [nuc_start + offset for offset in range(3)]
    alignment_indices = []
    for position in query_positions:
        alignment_index = coord_map.get(position)
        if alignment_index is None or alignment_index < 0:
            return None
        alignment_indices.append(alignment_index)
    return tuple(alignment_indices)

def extract_aligned_codon(alignment, alignment_indices):
    alignment = str(alignment or '')
    if not alignment_indices:
        return None
    if min(alignment_indices) < 0 or max(alignment_indices) >= len(alignment):
        return None
    codon = ''.join(alignment[index] for index in alignment_indices)
    if len(codon) != 3:
        return None
    return codon

def is_sequence_coverage_sufficient(alignment, ref_coord_map, gene_entry, start_aa, end_aa, target_pos):
    target_indices = resolve_aligned_codon_indices(ref_coord_map, gene_entry['cds_start'], target_pos)
    if target_indices is None:
        return False
    target_codon = extract_aligned_codon(alignment, target_indices)
    if target_codon is None or any(base in "-Nn?Xx" for base in target_codon):
        return False

    gap_count = 0
    for pos in range(start_aa, end_aa + 1):
        indices = resolve_aligned_codon_indices(ref_coord_map, gene_entry['cds_start'], pos)
        if indices is None:
            gap_count += 1
        else:
            codon = extract_aligned_codon(alignment, indices)
            if codon is None or any(base in "-Nn?Xx" for base in codon):
                gap_count += 1
    return gap_count <= 5


def load_master_accessions(conn):
    try:
        meta = pd.read_sql_query('SELECT primary_accession, accession_type FROM meta_data', conn)
    except Exception:
        return set()

    if not {'primary_accession', 'accession_type'}.issubset(meta.columns):
        return set()

    mask = meta['accession_type'].fillna('').astype(str).str.strip().str.lower() == 'master'
    return set(meta.loc[mask, 'primary_accession'].dropna().astype(str).str.strip())

def merge_feature_entry(feature_maps, reference_accession, product_name, start, end, feature_type, alias_lookup, raw_product=None, source='unknown'):
    entry = {
        'product': canonicalize_product(product_name, alias_lookup),
        'raw_product': raw_product or product_name,
        'cds_start': int(start),
        'cds_end': int(end),
        'reference_accession': reference_accession,
        'feature_type': feature_type,
        'source': source,
    }
    if not entry['product']:
        return
    reference_map = feature_maps.setdefault(reference_accession, {})
    existing_entry = reference_map.get(entry['product'])
    if existing_entry is None:
        reference_map[entry['product']] = entry
    else:
        def feature_priority(feature_type, product_name):
            canonical = normalize_lookup_key(product_name)
            if canonical in {'polyprotein', 'wholegenome'}:
                return 99
            priorities = {'mat_peptide': 0, 'mature_protein_region_of_CDS': 0, 'gene': 1, 'CDS': 2, 'cds': 2}
            return priorities.get(feature_type, 10)
        existing_priority = feature_priority(existing_entry.get('feature_type', ''), existing_entry.get('product', ''))
        new_priority = feature_priority(entry.get('feature_type', ''), entry.get('product', ''))
        if new_priority < existing_priority:
            reference_map[entry['product']] = entry
        elif new_priority == existing_priority:
            existing_span = int(existing_entry['cds_end']) - int(existing_entry['cds_start'])
            new_span = int(entry['cds_end']) - int(entry['cds_start'])
            if new_span < existing_span:
                reference_map[entry['product']] = entry

def load_db_gff_feature_maps(conn, alias_lookup):
    table_names = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", conn)['name'].tolist()
    candidate_tables = [name for name in ['gff_features', 'reference_gff_features', 'gff', 'reference_gff'] if name in table_names]
    feature_maps = {}
    master_accessions = load_master_accessions(conn)
    
    # Load all alignment names from sequence_alignment table to include genotype-specific references
    try:
        alignment_names = pd.read_sql_query('SELECT DISTINCT alignment_name FROM sequence_alignment', conn)['alignment_name'].dropna().tolist()
        allowed_accessions = master_accessions.union(alignment_names)
    except Exception:
        allowed_accessions = master_accessions

    for table_name in candidate_tables:
        df = pd.read_sql_query(f'SELECT * FROM {table_name}', conn)
        columns = set(df.columns)
        if {'reference_accession', 'start', 'end'}.issubset(columns):
            feature_type_col = 'feature_type' if 'feature_type' in columns else None
            product_col = next((name for name in ['product', 'gene_name', 'name', 'raw_product'] if name in columns), None)
            if product_col is None:
                continue
            for _, row in df.iterrows():
                reference_accession = str(row.get('reference_accession', '') or '').strip()
                if not reference_accession:
                    continue
                if allowed_accessions and reference_accession not in allowed_accessions:
                    continue
                feature_type = str(row.get(feature_type_col, 'gene') or 'gene') if feature_type_col else 'gene'
                merge_feature_entry(
                    feature_maps,
                    reference_accession,
                    row.get(product_col),
                    row.get('start'),
                    row.get('end'),
                    feature_type,
                    alias_lookup,
                    raw_product=row.get(product_col),
                    source=f'db_table:{table_name}',
                )

    if 'features' in table_names:
        df = pd.read_sql_query('SELECT * FROM features', conn)
        columns = set(df.columns)
        accession_col = 'accession' if 'accession' in columns else 'primary_accession' if 'primary_accession' in columns else None
        start_col = 'cds_start_OG_seq' if 'cds_start_OG_seq' in columns else 'cds_start' if 'cds_start' in columns else 'start' if 'start' in columns else None
        end_col = 'cds_end_OG_seq' if 'cds_end_OG_seq' in columns else 'cds_end' if 'cds_end' in columns else 'end' if 'end' in columns else None
        product_col = next((name for name in ['product', 'gene_name', 'name', 'raw_product'] if name in columns), None)
        feature_type_col = 'feature_type' if 'feature_type' in columns else None

        if accession_col and start_col and end_col and product_col:
            for _, row in df.iterrows():
                reference_accession = str(row.get(accession_col, '') or '').strip()
                if not reference_accession:
                    continue
                ref_candidate = str(row.get('reference_accession', '') or '').strip() if 'reference_accession' in columns else ''
                master_candidate = str(row.get('master_ref_accession', '') or '').strip() if 'master_ref_accession' in columns else ''
                is_allowed = reference_accession in allowed_accessions if allowed_accessions else True
                if not is_allowed:
                    continue
                if reference_accession not in {ref_candidate, master_candidate, reference_accession}:
                    continue
                feature_type = str(row.get(feature_type_col, 'gene') or 'gene') if feature_type_col else 'gene'
                merge_feature_entry(
                    feature_maps,
                    reference_accession,
                    row.get(product_col),
                    row.get(start_col),
                    row.get(end_col),
                    feature_type,
                    alias_lookup,
                    raw_product=row.get(product_col),
                    source='db_table:features',
                )

    return feature_maps

def get_ref_alignment(ref_acc, master_candidates, ref_alignments_map):
    if (ref_acc, ref_acc) in ref_alignments_map:
        return ref_alignments_map[(ref_acc, ref_acc)], ref_acc
    if ref_acc in ref_alignments_map:
        return ref_alignments_map[ref_acc], ref_acc

    for cand in master_candidates:
        if cand in ref_alignments_map:
            return ref_alignments_map[cand], cand

    return None, None

def resolve_feature_map(ref_acc, master_accessions, catalog_reference_hints, db_gff_maps):
    candidates = [ref_acc] + list(master_accessions) + list(catalog_reference_hints)
    seen = []
    for cand in candidates:
        for token in extract_accession_tokens(cand):
            if token and token not in seen:
                seen.append(token)
    for cand in seen:
        if cand in db_gff_maps:
            if any(normalize_lookup_key(p) not in {'polyprotein', 'wholegenome'} for p in db_gff_maps[cand]):
                return cand, db_gff_maps[cand]
    return None, None

def main():
    parser = argparse.ArgumentParser(description="Verify mutation annotations against database alignment columns.")
    parser.add_argument("--db", default="/home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/HCV_full.db", help="Path to SQLite database.")
    parser.add_argument("--mutation_catalog", default="/home3/oml4h/RABV-gTK/generic/hcv/Tables/generalized_mutation_catalog_with_extra_info.tsv", help="Path to mutation catalog TSV.")
    parser.add_argument("--sample_size", type=int, default=100, help="Number of annotated sequences to sample and verify.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling.")
    parser.add_argument("--hcv_test_ns3_36a", action="store_true", help="Run verification specifically for NS3:36A on the 171 query accessions.")
    parser.add_argument("--min_identity", type=float, default=0.65, help="Minimum nucleotide identity threshold for alignment validation.")
    parser.add_argument("--virus", default="HCV",
        help="Virus whose protein-name vocabulary to canonicalise with. Must match the "
             "--virus AnnotateMutations was run with, or this script will report "
             "mismatches that are artefacts of a different vocabulary.")
    args = parser.parse_args()

    global ACTIVE_VIRUS_PROFILE
    ACTIVE_VIRUS_PROFILE = virus_profile(args.virus)

    if not os.path.isfile(args.db):
        print(f"Error: Database {args.db} not found.", file=sys.stderr)
        sys.exit(1)
    if not args.hcv_test_ns3_36a and not os.path.isfile(args.mutation_catalog):
        print(f"Error: Catalog {args.mutation_catalog} not found.", file=sys.stderr)
        sys.exit(1)

    if args.hcv_test_ns3_36a:
        catalog_refs = []
    else:
        print("1. Loading mutation catalog...")
        catalog = pd.read_csv(args.mutation_catalog, sep='\t', dtype=str)
        catalog_refs = catalog['reference_accession'].dropna().unique()
        print(f"Reference accessions in catalog: {list(catalog_refs)}")
    
    conn = sqlite3.connect(args.db)
    
    try:
        master_accessions = load_master_accessions(conn)
        print(f"Master reference accessions in DB: {list(master_accessions)}")

        # 2. Check reference sequence flags
        if not args.hcv_test_ns3_36a:
            print("\n2. Checking reference sequence flagging...")
            reference_checking_failed = False
            for ref_raw in catalog_refs:
                tokens = extract_accession_tokens(ref_raw)
                found_master = False
                for token in tokens:
                    if token in master_accessions:
                        found_master = True
                        print(f"  [OK] Catalog reference '{ref_raw}' maps to DB master accession '{token}'.")
                        break
                if not found_master:
                    print(f"  [FAIL] Catalog reference '{ref_raw}' (tokens: {tokens}) does not map to any DB master accession.")
                    reference_checking_failed = True
                    for token in tokens:
                        cursor = conn.cursor()
                        cursor.execute("SELECT accession_type FROM meta_data WHERE primary_accession = ?", (token,))
                        row = cursor.fetchone()
                        if row:
                            print(f"    Token '{token}' exists in meta_data with accession_type = '{row[0]}'")
                        else:
                            print(f"    Token '{token}' does not exist in meta_data.")

            if reference_checking_failed:
                print("\nResult: Reference sequence flagging verification failed.", file=sys.stderr)
                sys.exit(2)

        # Load gene alias lookup and feature maps
        gene_alias_lookup = load_gene_alias_lookup(conn)
        db_gff_maps = load_db_gff_feature_maps(conn, gene_alias_lookup)

        catalog_reference_hints = []
        for ref_raw in catalog_refs:
            for token in extract_accession_tokens(ref_raw):
                if token and token not in catalog_reference_hints:
                    catalog_reference_hints.append(token)

        # 3. Load sequences with annotations
        if args.hcv_test_ns3_36a:
            test_accessions = [
                "KC770153", "KC770157", "HQ623402", "HM441243", "HM441242", "HM441244", "KP162539", "HM441239", 
                "HM441238", "HM441241", "HM441240", "KT162582", "KT162585", "MH714577", "MH714582", "MH714584", 
                "KT233299", "MH714576", "KP162569", "KP162370", "KX107853", "KX107860", "KT233290", "MH714589", 
                "MH714580", "MH714574", "MH714575", "MH714586", "MH714562", "PP856528", "HQ623298", "HQ623316", 
                "HQ623277", "HQ623296", "HQ623276", "HQ623305", "HQ623311", "KF489392", "KF489394", "KP162683", 
                "KF489391", "KF489395", "KF489396", "KP162575", "KF489390", "KY420587", "MG601771", "KY420586", 
                "OR917741", "KY420613", "MK962414", "KY420614", "KY420585", "KY420612", "HQ623361", "HQ623353", 
                "OM312335", "OM312307", "OM312295", "OM312283", "KT162436", "KT162426", "KT162409", "KT162373", 
                "KT162478", "KT162459", "KT162465", "KT162477", "KT162413", "KT162439", "KT162472", "KT162446", 
                "KT162448", "KP163378", "KT162458", "KT162378", "KP163131", "KT162438", "KT162367", "KT162375", 
                "KP163451", "KT162407", "KT162363", "KP163329", "KT162471", "KT162433", "KT162372", "KP163452", 
                "KT162418", "KP163279", "KT162371", "KT162362", "KT162432", "KT162361", "KP163175", "KP163186", 
                "KP163108", "KP162566", "KP163127", "KP162355", "KP163410", "KP163343", "KT162443", "KP162884", 
                "KM588737", "KM588733", "KM588739", "KT162467", "KP162691", "KP163542", "EU315121", "ON603924", 
                "ON603932", "MN650798", "MK609555", "MK609556", "MK609558", "MK609562", "MK609559", "MK609563", 
                "MK609561", "EU315118", "MK609560", "MK609557", "PQ799529", "EU709995", "EU710042", "DQ355728", 
                "EU710022", "EU710057", "EU710053", "DQ355722", "EU710051", "EU847448", "EU847451", "EU847454", 
                "DQ355652", "DQ355656", "EU709973", "DQ355719", "DQ355715", "PQ799566", "DQ355653", "EU709968", 
                "EU710014", "EU847453", "EU710034", "DQ355729", "EU710039", "DQ355714", "EU710029", "DQ355647", 
                "EU847455", "DQ355654", "DQ355717", "DQ355651", "DQ355648", "DQ355655", "EU710030", "DQ355727", 
                "DQ355649", "DQ355720", "DQ355723", "DQ355650", "DQ355730", "EU847450", "EU847449", "EU710074", 
                "EU847452", "MK962420", "EU312148"
            ]
            sampled_df = pd.DataFrame({
                'primary_accession': test_accessions,
                'relevant_mutations_present': ['NS3:36A'] * len(test_accessions)
            })
            print(f"Bypassed sampling. Running verification for NS3:36A on the {len(test_accessions)} specified accessions.")
        else:
            print("\n3. Loading annotated sequences from database...")
            annotated_df = pd.read_sql_query(
                "SELECT primary_accession, relevant_mutations_present FROM sequence_relevant_mutation_summary", conn
            )
            print(f"Total annotated sequences in DB: {len(annotated_df)}")
            if annotated_df.empty:
                print("Warning: sequence_relevant_mutation_summary table is empty. Nothing to verify.", file=sys.stderr)
                sys.exit(0)

            # Sample sequences
            sample_size = min(args.sample_size, len(annotated_df))
            sampled_df = annotated_df.sample(n=sample_size, random_state=args.seed).reset_index(drop=True)
            print(f"Sampled {sample_size} sequences (seed={args.seed}) for verification.")

        # Query alignments for the sampled sequences in chunks
        sampled_accessions = sampled_df['primary_accession'].tolist()
        chunk_size = 500
        alignments_df_list = []
        for i in range(0, len(sampled_accessions), chunk_size):
            chunk = sampled_accessions[i : i + chunk_size]
            placeholders = ','.join('?' for _ in chunk)
            chunk_df = pd.read_sql_query(
                f"SELECT primary_accession, alignment, alignment_name FROM sequence_alignment WHERE primary_accession IN ({placeholders})",
                conn,
                params=chunk
            )
            alignments_df_list.append(chunk_df)
        alignments_df = pd.concat(alignments_df_list, ignore_index=True)
        alignments_map = {row['primary_accession']: row for _, row in alignments_df.iterrows()}

        # Preload reference alignments in chunks
        print("Preloading reference alignments...")
        ref_accessions = set(alignments_df['alignment_name'].dropna().unique())
        for cand in list(master_accessions) + list(catalog_reference_hints):
            for token in extract_accession_tokens(cand):
                if token:
                    ref_accessions.add(token)

        ref_accessions_list = list(ref_accessions)
        ref_alignments_df_list = []
        for i in range(0, len(ref_accessions_list), chunk_size):
            chunk = ref_accessions_list[i : i + chunk_size]
            ref_placeholders = ','.join('?' for _ in chunk)
            chunk_df = pd.read_sql_query(
                f"SELECT primary_accession, alignment_name, alignment FROM sequence_alignment WHERE primary_accession IN ({ref_placeholders})",
                conn,
                params=chunk
            )
            ref_alignments_df_list.append(chunk_df)
        ref_alignments_df = pd.concat(ref_alignments_df_list, ignore_index=True)
        ref_alignments_map = {}
        for _, row in ref_alignments_df.iterrows():
            p_acc = row['primary_accession']
            aln_name = row['alignment_name']
            aln = row['alignment']
            ref_alignments_map[(p_acc, aln_name)] = aln
            ref_alignments_map[p_acc] = aln

        # -------------------------------------------------------------------
        # Resolve master reference coord map and feature map ONCE, before the
        # per-sequence loop.  Catalog positions (aa_position) are defined
        # relative to the master reference, and since ALL padded alignments
        # share the same column space, the master's column indices apply to
        # every sequence regardless of which genotype reference it was
        # locally aligned against.
        # -------------------------------------------------------------------
        master_ref_acc_for_coords = next(iter(master_accessions), None)
        master_coord_map = None
        master_feature_map = None
        master_feature_acc = None

        coord_search_order = (
            list(master_accessions) + catalog_reference_hints
        )
        for _cand in coord_search_order:
            for _token in extract_accession_tokens(_cand):
                _ref_aln = ref_alignments_map.get(_token)
                if _ref_aln is None:
                    continue
                _fmap_acc, _fmap = resolve_feature_map(
                    _token, master_accessions, catalog_reference_hints, db_gff_maps
                )
                if _fmap is not None:
                    master_ref_acc_for_coords = _token
                    master_coord_map = build_alignment_coordinate_map(_ref_aln)
                    master_feature_map = _fmap
                    master_feature_acc = _fmap_acc
                    break
            if master_coord_map is not None:
                break

        if master_coord_map is None:
            print(
                "[ERROR] Could not resolve a master reference alignment for coordinate resolution. "
                f"Candidates tried: {coord_search_order}",
                file=sys.stderr,
            )
            sys.exit(2)

        print(f"  [INFO] Coordinate reference: {master_ref_acc_for_coords} (used for all catalog position resolution)")

        print("\n4. Verifying mutations against alignments...")
        stats = Counter()
        failures = []
        checked_mutations = []

        for idx, row in sampled_df.iterrows():
            acc = row['primary_accession']
            muts_str = row['relevant_mutations_present']
            
            if not muts_str:
                continue

            mutations = [m.strip() for m in muts_str.split(';') if m.strip()]
            
            aln_row = alignments_map.get(acc)
            if aln_row is None:
                print(f"  [ERROR] No alignment found for sequence {acc} in database.")
                stats['sequence_errors'] += 1
                continue

            query_alignment = aln_row['alignment']
            ref_acc = aln_row['alignment_name']

            # Use per-sequence reference alignment only for identity scoring
            ref_alignment, resolved_ref_acc = get_ref_alignment(ref_acc, master_accessions, ref_alignments_map)
            if ref_alignment is None:
                print(f"  [ERROR] Could not resolve reference alignment for group {ref_acc} for sequence {acc}.")
                stats['reference_alignment_errors'] += 1
                continue

            # Use master coord_map and feature_map for all codon position resolution
            ref_coord_map = master_coord_map
            feature_map = master_feature_map
            resolved_ref_feature_acc = master_feature_acc

            # Calculate full nucleotide alignment identity directly from the padded alignment in the DB
            total_nuc_overlap = 0
            matching_nuc_overlap = 0
            for q_base, r_base in zip(query_alignment, ref_alignment):
                q_up = q_base.upper()
                r_up = r_base.upper()
                if q_up != '-' and q_up != 'N' and r_up != '-' and r_up != 'N':
                    total_nuc_overlap += 1
                    if q_up == r_up:
                        matching_nuc_overlap += 1
            
            overlap_identity = matching_nuc_overlap / total_nuc_overlap if total_nuc_overlap > 0 else 0.0
            is_alignment_valid = (overlap_identity >= args.min_identity)
            if not is_alignment_valid:
                print(f"  [WARNING] Sequence {acc} has low nucleotide identity to reference {resolved_ref_acc}: {overlap_identity:.1%} ({matching_nuc_overlap}/{total_nuc_overlap})")
                stats['low_identity_alignments'] += 1

            for mut in mutations:
                stats['mutations_checked'] += 1
                
                match = re.match(r'^([^:]+):(\d+)([A-Z_X])$', mut)
                if not match:
                    print(f"  [ERROR] Cannot parse mutation '{mut}' for sequence {acc}.")
                    stats['parse_errors'] += 1
                    continue

                protein_name = match.group(1)
                aa_pos = int(match.group(2))
                expected_aa = match.group(3)

                gene_entry = feature_map.get(protein_name)
                if gene_entry is None:
                    print(f"  [ERROR] Protein '{protein_name}' not found in master feature map of {resolved_ref_feature_acc} (mutation: {mut}, sequence: {acc}).")
                    stats['protein_not_found_errors'] += 1
                    continue

                # Always resolve column indices using master reference coordinates
                alignment_indices = resolve_aligned_codon_indices(
                    master_coord_map, gene_entry['cds_start'], aa_pos
                )
                if alignment_indices is None:
                    print(f"  [ERROR] Cannot resolve alignment columns for protein '{protein_name}' pos {aa_pos} (mutation: {mut}, sequence: {acc}).")
                    stats['column_resolution_errors'] += 1
                    continue

                codon = extract_aligned_codon(query_alignment, alignment_indices)
                if codon is None:
                    print(f"  [ERROR] Extracted codon is out of bounds for sequence {acc} at indices {alignment_indices} (mutation: {mut}).")
                    stats['codon_out_of_bounds_errors'] += 1
                    continue

                actual_aa = translate_codon(codon)
                
                checked_mutations.append({
                    'mutation': mut,
                    'protein': protein_name,
                    'position': aa_pos,
                    'expected_aa': expected_aa,
                    'sequence_id': acc,
                    'actual_aa': actual_aa if is_alignment_valid else f"Invalid (Low Identity {overlap_identity:.1%})",
                    'codon': codon,
                    'gene_entry': gene_entry,
                    'ref_coord_map': master_coord_map,
                    'resolved_ref_acc': master_ref_acc_for_coords,
                    'ref_alignment': ref_alignments_map.get(master_ref_acc_for_coords, ref_alignment),
                    'overlap_identity': overlap_identity,
                })

                if is_alignment_valid and actual_aa == expected_aa:
                    stats['matches'] += 1
                else:
                    stats['mismatches'] += 1
                    fail_info = {
                        'sequence_id': acc,
                        'mutation': mut,
                        'protein': protein_name,
                        'position': aa_pos,
                        'expected_aa': expected_aa,
                        'actual_aa': actual_aa if is_alignment_valid else f"Invalid (Low Identity {overlap_identity:.1%})",
                        'codon': codon,
                        'alignment_name': ref_acc,
                        'resolved_ref_acc': resolved_ref_acc,
                        'alignment_indices': alignment_indices
                    }
                    failures.append(fail_info)
                    if not is_alignment_valid:
                        print(
                            f"  [MISMATCH] Sequence {acc} fails alignment validation: low nucleotide identity ({overlap_identity:.1%}) to reference {resolved_ref_acc}."
                        )
                    else:
                        print(
                            f"  [MISMATCH] Sequence {acc} has mutation {mut} but alignment has codon {codon} -> {actual_aa} at columns {alignment_indices}."
                        )

            stats['sequences_checked'] += 1

        # Select 4 random mutation signatures and print alignment snippets
        if checked_mutations:
            import random
            rng = random.Random(args.seed)
            
            # Group by mutation signature and reference accession to keep coordinates aligned
            grouped_mutations = {}
            for item in checked_mutations:
                sig_ref_key = (item['protein'], item['position'], item['expected_aa'], item['resolved_ref_acc'])
                if sig_ref_key not in grouped_mutations:
                    grouped_mutations[sig_ref_key] = []
                grouped_mutations[sig_ref_key].append(item)
            
            # Select 4 mutation signature-reference combinations randomly
            if args.hcv_test_ns3_36a:
                selected_keys = list(grouped_mutations.keys())
            else:
                selected_keys = rng.sample(list(grouped_mutations.keys()), min(8, len(grouped_mutations)))
            
            print("\n====================================================================================================")
            print("5. VISUAL VERIFICATION: ALIGNMENT SNIPPETS FOR RANDOM MUTATION SIGNATURES")
            print("====================================================================================================")
            
            for sig_ref_key in selected_keys:
                protein, position, expected_aa, resolved_ref_acc = sig_ref_key
                items = grouped_mutations[sig_ref_key]
                
                # Use reference sequence and coordinate metadata from first matching item in this group
                first_item = items[0]
                gene_entry = first_item['gene_entry']
                ref_coord_map = first_item['ref_coord_map']
                ref_alignment = first_item['ref_alignment']
                
                # Identify mutated sequences in the sampled set for this specific reference
                mutated_seqs = list(set(item['sequence_id'] for item in items))
                
                # Find wild-type / other sequences that share the exact same reference sequence
                other_seqs = []
                for s_acc in sampled_df['primary_accession']:
                    if s_acc in mutated_seqs:
                        continue
                    s_aln_row = alignments_map.get(s_acc)
                    if s_aln_row is None:
                        continue
                    _, s_resolved_ref_acc = get_ref_alignment(s_aln_row['alignment_name'], master_accessions, ref_alignments_map)
                    if s_resolved_ref_acc == resolved_ref_acc:
                        other_seqs.append(s_acc)
                
                # Determine AA position window (20 AA window centered around mutation if possible)
                max_aa = (gene_entry['cds_end'] - gene_entry['cds_start'] + 1) // 3
                start_aa = max(1, position - 9)
                end_aa = start_aa + 19
                if end_aa > max_aa:
                    end_aa = max_aa
                    start_aa = max(1, end_aa - 19)

                # Filter sequences by coverage — only show sequences that have real data at the target position
                covered_mutated = [
                    s_acc for s_acc in mutated_seqs
                    if is_sequence_coverage_sufficient(alignments_map[s_acc]['alignment'], ref_coord_map, gene_entry, start_aa, end_aa, position)
                ]
                covered_others = [
                    s_acc for s_acc in other_seqs
                    if is_sequence_coverage_sufficient(alignments_map[s_acc]['alignment'], ref_coord_map, gene_entry, start_aa, end_aa, position)
                ]

                # Do NOT fall back to uncovered sequences — showing gap-only sequences creates
                # misleading 0% identity warnings. Only display sequences with real coverage.
                filtered_mutated = covered_mutated
                filtered_others = covered_others

                # Display up to 5 comparison query sequences (prioritizing mutated ones, and filling with others)
                selected_mutated = filtered_mutated[:5]
                num_needed = 5 - len(selected_mutated)
                selected_others = rng.sample(filtered_others, min(num_needed, len(filtered_others))) if num_needed > 0 and filtered_others else []
                display_seqs = selected_mutated + selected_others

                # Skip this signature if no sequences have coverage at the target position
                if not display_seqs:
                    print(f"\nMutation Signature: {protein}:{position}{expected_aa} [SKIPPED - no covered sequences in sample]")
                    continue
                
                # Extract reference alignment window
                ref_window_data = []
                for pos in range(start_aa, end_aa + 1):
                    indices = resolve_aligned_codon_indices(ref_coord_map, gene_entry['cds_start'], pos)
                    if indices is None:
                        ref_window_data.append(('X', '---'))
                    else:
                        codon = extract_aligned_codon(ref_alignment, indices)
                        if codon is None:
                            ref_window_data.append(('X', '---'))
                        else:
                            aa = translate_codon(codon)
                            ref_window_data.append((aa, codon))
                
                # Extract query alignment windows
                query_window_data = {}
                for s_acc in display_seqs:
                    s_aln_row = alignments_map[s_acc]
                    s_alignment = s_aln_row['alignment']
                    s_window = []
                    for pos in range(start_aa, end_aa + 1):
                        indices = resolve_aligned_codon_indices(ref_coord_map, gene_entry['cds_start'], pos)
                        if indices is None:
                            s_window.append(('X', '---'))
                        else:
                            codon = extract_aligned_codon(s_alignment, indices)
                            if codon is None:
                                s_window.append(('X', '---'))
                            else:
                                aa = translate_codon(codon)
                                s_window.append((aa, codon))
                    query_window_data[s_acc] = s_window
                
                print(f"\nMutation Signature: {protein}:{position}{expected_aa}")
                print(f"Reference: {resolved_ref_acc} ({protein} CDS: {gene_entry['cds_start']}-{gene_entry['cds_end']})")
                print(f"Window: AA positions {start_aa} to {end_aa} (Target: {position})")
                print("-" * 135)
                print(f"{'Sequence ID':<20} | {'Type':<9} | {'Win Identity':<15} | {'Full Nuc Identity':<18} | {'Amino Acid Window (20 AA)':<60}")
                print("-" * 135)
                
                # Print Reference sequence representation
                ref_aas = [item[0] for item in ref_window_data]
                ref_aa_parts = []
                for i, pos in enumerate(range(start_aa, end_aa + 1)):
                    aa, _ = ref_window_data[i]
                    is_target = (pos == position)
                    if is_target:
                        ref_aa_parts.append(f"[{aa}]")
                    else:
                        ref_aa_parts.append(f" {aa} ")
                print(f"{resolved_ref_acc:<20} | {'Reference':<9} | {'100.0% (Ref)':<15} | {'100.0% (Ref)':<18} | {''.join(ref_aa_parts)}")
                
                # Compute target-position index in the window list
                target_window_idx = position - start_aa

                # Print query sequences representations
                # Verification check: all mutated sequences must have the expected AA at the target position
                target_mismatches = []
                target_gaps = []
                for s_acc in display_seqs:
                    s_type = 'Mutated' if s_acc in mutated_seqs else 'WT/Other'
                    s_window = query_window_data[s_acc]
                    q_aas = [item[0] for item in s_window]

                    # Calculate window identity (for informational display)
                    matching = sum(1 for r_aa, q_aa in zip(ref_aas, q_aas) if r_aa == q_aa and r_aa != 'X')
                    total_non_gap = sum(1 for r_aa in ref_aas if r_aa != 'X')
                    identity_percent = (matching / total_non_gap) * 100 if total_non_gap else 0.0
                    identity_str = f"{identity_percent:.1f}% ({matching}/{total_non_gap})"

                    # Retrieve overlap identity from checked_mutations items
                    s_item = next((it for it in items if it['sequence_id'] == s_acc), None)
                    if s_item is not None:
                        full_identity_val = s_item['overlap_identity']
                        full_identity_str = f"{full_identity_val:.1%}"
                    else:
                        # Fallback calculation if not in checked_mutations (e.g. wild-type/other sequence)
                        s_aln_row = alignments_map.get(s_acc)
                        if s_aln_row is not None:
                            s_alignment = s_aln_row['alignment']
                            total_nuc_overlap = sum(1 for q in s_alignment if q != '-')
                            matching_nuc_overlap = sum(1 for q, r in zip(s_alignment, ref_alignment) if q != '-' and q.upper() == r.upper())
                            full_identity_val = matching_nuc_overlap / total_nuc_overlap if total_nuc_overlap > 0 else 0.0
                            full_identity_str = f"{full_identity_val:.1%}"
                        else:
                            full_identity_str = "N/A"

                    # Check target position specifically
                    target_aa = q_aas[target_window_idx] if 0 <= target_window_idx < len(q_aas) else 'X'
                    if s_type == 'Mutated':
                        if target_aa in ('X', '-'):
                            target_gaps.append(s_acc)
                        elif target_aa != expected_aa:
                            target_mismatches.append((s_acc, target_aa))

                    s_aa_parts = []
                    for i, pos in enumerate(range(start_aa, end_aa + 1)):
                        aa, _ = s_window[i]
                        is_target = (pos == position)
                        if is_target:
                            s_aa_parts.append(f"[{aa}]")
                        else:
                            s_aa_parts.append(f" {aa} ")
                    print(f"{s_acc:<20} | {s_type:<9} | {identity_str:<15} | {full_identity_str:<18} | {''.join(s_aa_parts)}")
                print("-" * 135)

                # Alignment verification based on target position, not whole-window identity.
                # Low window identity is EXPECTED for cross-genotype comparisons.
                if target_mismatches:
                    for seq_id, actual in target_mismatches:
                        print(f"Alignment verification: FAIL (sequence {seq_id} has {actual} at target {protein}:{position}, expected {expected_aa})")
                elif target_gaps:
                    for seq_id in target_gaps:
                        print(f"Alignment verification: WARNING (sequence {seq_id} has gap/ambiguous at target {protein}:{position})")
                else:
                    # Check if window identity is low (informational, not a failure)
                    low_window = [s for s in display_seqs if
                        (sum(1 for r, q in zip(ref_aas, [query_window_data[s][i][0] for i in range(len(ref_aas))]) if r == q and r != 'X') /
                          max(1, sum(1 for r in ref_aas if r != 'X'))) * 100 < 30.0
                    ]
                    if low_window:
                        print(f"Alignment verification: PASS (target codon correct; note: {len(low_window)} sequence(s) show low window identity, expected for cross-genotype comparisons)")
                    else:
                        print(f"Alignment verification: PASS (target codon correct for all displayed sequences)")
                print("-" * 135)
        else:
            print("\n[WARNING] No verified mutations available to print alignment snippets.")

        print("\nVerification Summary:")
        print(f"  Sequences Checked: {stats['sequences_checked']}")
        print(f"  Mutations Checked: {stats['mutations_checked']}")
        print(f"  Matches:           {stats['matches']}")
        print(f"  Mismatches:        {stats['mismatches']}")
        
        err_keys = [k for k in stats.keys() if k.endswith('errors')]
        if err_keys:
            print("  Errors:")
            for k in sorted(err_keys):
                print(f"    {k}: {stats[k]}")

        if stats['mismatches'] > 0 or any(stats[k] > 0 for k in err_keys):
            print("\nResult: VERIFICATION FAILED", file=sys.stderr)
            sys.exit(2)
        else:
            print("\nResult: VERIFICATION SUCCESSFUL")
            sys.exit(0)

    finally:
        conn.close()

if __name__ == "__main__":
    main()
