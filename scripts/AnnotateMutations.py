#!/usr/bin/env python3

import sqlite3
import pandas as pd
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

from Bio import Entrez, SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from ExploreMutationStorageLayouts import build_completed_signatures_only, build_sequence_relevant_mutation_summary

REQUIRED_MUTATION_CATALOG_COLUMNS = [
    'mutation_id',
    'protein_name',
    'segment',
    'aa_position',
    'alt_residue',
    'reference_accession',
    'mutation_type',
    'signature_id',
    'signature_kind',
    'combination_id',
    'combination_size',
    'phenotype',
]

CATALOG_COLUMN_PROFILES = {
    'HCV': ['resistance_category', 'drug','drug_category','drug_producer','pubmed_id','DOI','any_in_vitro_evidence','in_vitro_max_ec50_midpoint','any_in_vivo_evidence','in_vivo_baseline','in_vivo_treatment_emergent'],
    'influenza': ['serotypes_tested'],
    'all_columns': None,
}

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
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_codon(codon):
    codon = codon.upper().replace('-', '')
    if len(codon) < 3:
        return 'X'
    return CODON_TABLE.get(codon, 'X')


# ---------------------------------------------------------------------------
# Residue vocabulary.
#
# One spelling per concept, everywhere: a stop codon is '*' and a deleted
# residue is '-'.  The legacy spellings still have to be *readable* because
# catalogs shipped before this change spell stop '_' (the old CODON_TABLE
# value) and deletion 'del' (PHDR's own wording), but they are never written.
# ---------------------------------------------------------------------------
STOP_RESIDUE = '*'
DELETION_RESIDUE = '-'

def clean_cell(value):
    """Text of a possibly-missing dataframe cell, with NaN read as empty.

    ``str(value or '')`` is not enough: float('nan') is truthy, so a missing TSV
    cell becomes the literal string 'nan' and then looks like a real genotype
    code or residue.  Every parser below goes through here instead.
    """
    if value is None:
        return ''
    try:
        if pd.isna(value):
            return ''
    except (TypeError, ValueError):
        pass
    return str(value).strip()


LEGACY_RESIDUE_SPELLINGS = {
    '_': STOP_RESIDUE,
    'STOP': STOP_RESIDUE,
    'DEL': DELETION_RESIDUE,
    'DELETION': DELETION_RESIDUE,
}


def normalize_residue(value):
    """Fold a residue token from any source onto the standard vocabulary.

    Accepts the legacy spellings ('_' for stop, 'del' for a deletion) and
    tolerates the stray whitespace and lower case that hand-maintained TSVs
    accumulate, so a curator typo cannot silently disable a catalog row.
    """
    text = clean_cell(value)
    if not text:
        return ''
    upper = text.upper()
    if upper in LEGACY_RESIDUE_SPELLINGS:
        return LEGACY_RESIDUE_SPELLINGS[upper]
    if text in (STOP_RESIDUE, DELETION_RESIDUE):
        return text
    return upper


def alignment_covered_span(alignment):
    """First and last non-gap column of a padded alignment, or None if all gaps.

    Everything outside this span is padding added to square the alignment up -
    it is sequence the submitter never reported, not sequence that is missing
    from the virus.
    """
    text = str(alignment or '')
    first = None
    last = None
    for index, base in enumerate(text):
        if base != '-':
            if first is None:
                first = index
            last = index
    if first is None:
        return None
    return (first, last)


def residue_from_aligned_codon(codon, covered_span=None, alignment_indices=None):
    """Translate an aligned codon, or report a deletion when it is all gaps.

    A deletion has no codon to translate: the evidence for it is that every
    aligned column of the codon is a gap.  translate_codon() strips the gaps and
    then fails its length guard, so it can only ever answer 'X' here - which is
    indistinguishable from an unsequenced codon.  Reading the gap *before*
    translating is what makes catalogued deletions (NS5A 29/30/32del) findable.

    A gap only counts as a deletion when it sits *inside* the sequence's covered
    span.  A partial GenBank record is padded with gaps out to the full genome
    width, and in the shipped HCV build that padding covers NS5A 29/30/32 in 35
    sequences - reading it as a deletion would turn "we never sequenced this"
    into a resistance call on the worst-covered records.  Callers that have no
    alignment context (a unit test handing over a bare codon) get the plain
    deletion reading.
    """
    if codon is None:
        return None
    text = str(codon).strip()
    if text and set(text) == {'-'}:
        if covered_span is None or not alignment_indices:
            return DELETION_RESIDUE
        first_covered, last_covered = covered_span
        if min(alignment_indices) > first_covered and max(alignment_indices) < last_covered:
            return DELETION_RESIDUE
        return translate_codon(text)
    return translate_codon(text)


# ---------------------------------------------------------------------------
# Genotype scope (rule A) and per-genotype wild type (rule B).
# ---------------------------------------------------------------------------
ALIGNMENT_NAME_PREFIX = 'AL_'

# WT + position + ALT, where WT may be several residues separated by '/' and
# either side may be spelled 'del': R155C, D168A, K/Q80K, L/M31L, 32del.
DISPLAY_STRUCTURE_COMPONENT_RE = re.compile(
    r'^(?P<wild_type>(?:del|[A-Za-z*])(?:/(?:del|[A-Za-z*]))*)?'
    r'(?P<position>\d+)'
    r'(?P<alt_residue>del|[A-Za-z*])$'
)

SCOPE_TIER_SUBTYPE = 'subtype'
SCOPE_TIER_GENOTYPE = 'genotype'
SCOPE_TIER_UNSCOPED = 'unscoped'
SCOPE_TIER_OUT_OF_SCOPE = 'out_of_scope'

RESIDUE_STATUS_CHANGE = 'change'
RESIDUE_STATUS_ANCHOR = 'anchor'
RESIDUE_STATUS_WT_UNKNOWN = 'wt_unknown'

#: Generic, virus-agnostic catalog columns. Both are OPTIONAL: a catalog with
#: neither behaves exactly as it did before genotype gating existed, which is
#: what keeps non-HCV viruses - which have no subtypes at all - working
#: unchanged. Everything HCV-specific (alignment codes such as AL_1a) is
#: resolved into these at catalog build time by
#: scripts/BuildCatalogGenotypeColumns.py, so no pipeline code depends on it.
RELEVANT_GENOTYPES_COLUMN = 'relevant_genotypes'
WILD_TYPE_RESIDUES_COLUMN = 'wild_type_residues'

#: Semicolon between entries, colon between an entry's fields.
GENOTYPE_ENTRY_SEP = ';'
GENOTYPE_FIELD_SEP = ':'

CALL_STATUS_EMITTED = 'emitted'
CALL_STATUS_SUPPRESSED_OUT_OF_SCOPE = 'suppressed_out_of_scope'
CALL_STATUS_SUPPRESSED_WILD_TYPE = 'suppressed_wild_type'

SCOPE_TIER_RANK = {
    SCOPE_TIER_SUBTYPE: 0,
    SCOPE_TIER_GENOTYPE: 1,
    SCOPE_TIER_UNSCOPED: 2,
    SCOPE_TIER_OUT_OF_SCOPE: 3,
}



def genotype_of_code(value):
    """Leading digits of a genotype/subtype code: ``6xd`` -> ``6``, ``1a`` -> ``1``."""
    match = re.match(r'\d+', clean_cell(value))
    return match.group(0) if match else ''


def build_subtype_code(genotype, subtype):
    """Join a genotype and a subtype letter into a catalog-style code.

    meta_data stores them apart ('1' + 'a'); the catalog spells them together
    ('1a').  A subtype that already carries its genotype is left alone, and a
    value the catalog has no vocabulary for (e.g. 'NA') simply fails to match
    any subtype bucket and falls through to the genotype tier.
    """
    genotype_text = clean_cell(genotype)
    subtype_text = clean_cell(subtype)
    if not subtype_text:
        return genotype_text
    if not genotype_text:
        return subtype_text
    if subtype_text.lower().startswith(genotype_text.lower()):
        return subtype_text
    return f'{genotype_text}{subtype_text}'



def parse_genotype_entry_list(value):
    """Parse a semicolon-separated ``code[:field...[:frequency]]`` column.

    Returns ``{code: [fields...]}``. Used for both generic genotype columns:

        relevant_genotypes   1a:36.09;3:0.10        -> {'1a': ['36.09'], '3': ['0.10']}
        relevant_genotypes   1a;1b                  -> {'1a': [],        '1b': []}
        wild_type_residues   1a:Q:60.89;1b:R:92.26  -> {'1a': ['Q','60.89'], ...}
        wild_type_residues   1a:Q;1b:R              -> {'1a': ['Q'], '1b': ['R']}

    The trailing frequency is OPTIONAL throughout, so a virus that knows its
    wild types but has no frequency data can still supply the column.
    """
    parsed = {}
    for entry in clean_cell(value).split(GENOTYPE_ENTRY_SEP):
        entry = entry.strip()
        if not entry:
            continue
        fields = [part.strip() for part in entry.split(GENOTYPE_FIELD_SEP)]
        code = fields[0]
        if code:
            parsed[code] = fields[1:]
    return parsed


def build_wild_type_tables(catalog):
    """Read the per-genotype wild type from the generic ``wild_type_residues``
    column.

    Deliberately NOT derived from ``alignment_name`` / ``display_structure``.
    Those are PHDR/HCV artefacts - a rabies or influenza catalog will never
    have them - so all genotype resolution happens once at catalog build time
    (see scripts/BuildCatalogGenotypeColumns.py) and the pipeline reads only
    columns any virus can supply.

    The column is optional. Without it the tables are empty, nothing is ever
    suppressed, and behaviour is exactly what it was before genotype gating
    existed.

    Returns ``(by_subtype, by_genotype)`` keyed by ``(code, protein, position)``.
    A genotype-level entry is the union of its subtypes' wild types, which is
    the right fallback for a sequence whose own subtype has no data.
    """
    by_subtype = defaultdict(set)
    by_genotype = defaultdict(set)
    if WILD_TYPE_RESIDUES_COLUMN not in catalog.columns:
        return {}, {}
    protein_column = '_canonical_protein' if '_canonical_protein' in catalog.columns else 'protein_name'
    for _, row in catalog.iterrows():
        protein = clean_cell(row.get(protein_column, ''))
        position = clean_cell(row.get('aa_position', ''))
        if not protein or not position:
            continue
        for code, fields in parse_genotype_entry_list(row.get(WILD_TYPE_RESIDUES_COLUMN, '')).items():
            if not fields:
                continue
            residue = normalize_residue(fields[0])
            if not residue:
                continue
            by_subtype[(code, protein, str(position))].add(residue)
            by_genotype[(genotype_of_code(code), protein, str(position))].add(residue)
    return dict(by_subtype), dict(by_genotype)


def lookup_wild_type_residues(wild_type_tables, subtype_code, genotype, protein, position):
    """Wild type for this sequence's genotype: exact subtype first, then genotype.

    Returns ``(residues, tier)``; ``(None, '')`` when no wild type is known,
    which must NOT be read as 'the residue is a change' - it is 'we do not
    know', and the caller records it as such rather than suppressing.
    """
    by_subtype, by_genotype = wild_type_tables
    if subtype_code:
        residues = by_subtype.get((subtype_code, protein, str(position)))
        if residues:
            return residues, SCOPE_TIER_SUBTYPE
    if genotype:
        residues = by_genotype.get((genotype, protein, str(position)))
        if residues:
            return residues, SCOPE_TIER_GENOTYPE
    return None, ''


def parse_relevant_genotypes(value):
    """Genotype codes from the generic column.

    Entries are semicolon separated and each may carry an optional trailing
    frequency: ``1a:36.09;3:0.10`` and ``1a;3`` both yield ``{'1a', '3'}``.
    A comma is still accepted so a catalog written before the separator was
    settled keeps working.
    """
    text = clean_cell(value).replace(',', GENOTYPE_ENTRY_SEP)
    return set(parse_genotype_entry_list(text))


def build_signature_genotype_scope(catalog):
    """Genotype scope per signature: the union over every row of that signature.

    Read only from the generic ``relevant_genotypes`` column. Deriving scope
    from ``alignment_name`` would tie the pipeline to a PHDR/HCV artefact that
    no other virus has, so that resolution happens once at catalog build time
    (scripts/BuildCatalogGenotypeColumns.py) and never here.

    The column is optional. Without it no signature has a scope, the gate lets
    everything through, and behaviour is what it was before gating existed - so
    a virus with no genotypes, or no subtypes, needs to supply nothing.
    A signature whose entry is blank likewise gets an empty set, which the gate
    reads as 'applies to any genotype'.
    """
    if 'signature_id' not in catalog.columns:
        return {}
    if RELEVANT_GENOTYPES_COLUMN not in catalog.columns:
        return {}

    scope = defaultdict(set)
    for _, row in catalog.iterrows():
        signature_id = clean_cell(row.get('signature_id', ''))
        if not signature_id:
            continue
        scope[signature_id] |= parse_relevant_genotypes(row.get(RELEVANT_GENOTYPES_COLUMN, ''))
    return {signature_id: frozenset(codes) for signature_id, codes in scope.items()}


def classify_genotype_scope(scope_codes, genotype, subtype_code):
    """Rule A, matched at genotype level.

    Exact subtype wins; otherwise any bucket whose leading digits agree with the
    sequence's genotype counts, because observed subtypes such as 6n, 6xd, 6xc
    and 4v have no catalog bucket of their own and exact matching would throw
    away nearly every call.  An empty scope means the entry applies anywhere.
    """
    if not scope_codes:
        return SCOPE_TIER_UNSCOPED
    if subtype_code and subtype_code in scope_codes:
        return SCOPE_TIER_SUBTYPE
    if genotype and any(genotype_of_code(code) == genotype for code in scope_codes):
        return SCOPE_TIER_GENOTYPE
    return SCOPE_TIER_OUT_OF_SCOPE


def build_sequence_genotype_map(meta_data):
    """``{accession: (genotype, subtype_code)}`` from meta_data, tolerant of schema drift."""
    genotype_map = {}
    if meta_data is None or getattr(meta_data, 'empty', True):
        return genotype_map
    columns = set(meta_data.columns)
    if 'primary_accession' not in columns:
        return genotype_map
    genotype_column = next(
        (name for name in ['nearest_reference_genotype', 'genotype'] if name in columns), None
    )
    subtype_column = next(
        (name for name in ['nearest_reference_subtype', 'subtype'] if name in columns), None
    )
    if genotype_column is None and subtype_column is None:
        return genotype_map

    def clean(value):
        text = clean_cell(value)
        return '' if text.lower() in {'', 'nan', 'none', 'null'} else text

    for _, row in meta_data.iterrows():
        accession = clean(row.get('primary_accession'))
        if not accession:
            continue
        genotype = clean(row.get(genotype_column)) if genotype_column else ''
        subtype = clean(row.get(subtype_column)) if subtype_column else ''
        subtype_code = build_subtype_code(genotype, subtype)
        if not genotype:
            genotype = genotype_of_code(subtype_code)
        genotype_map[accession] = (genotype, subtype_code)
    return genotype_map


class AnnotationMappingError(RuntimeError):
    pass


def normalize_lookup_key(value):
    return re.sub(r'[^a-z0-9]+', '', str(value or '').strip().lower())


def normalize_segment(value):
    if pd.isna(value):
        return ''
    text = str(value or '').strip()
    if not text:
        return ''
    try:
        number = float(text)
        if number.is_integer():
            return str(int(number))
    except ValueError:
        pass
    return text


def infer_compact_product_name(product_name):
    text = str(product_name or '').strip()
    if not text:
        return ''

    ns_match = re.search(r'\b(NS[2-5](?:A|B)?)\b', text, flags=re.IGNORECASE)
    if ns_match:
        return ns_match.group(1).upper()

    envelope_match = re.search(r'\b(E[12])\b', text, flags=re.IGNORECASE)
    if envelope_match:
        return envelope_match.group(1).upper()

    if re.search(r'\bp7\b', text, flags=re.IGNORECASE):
        return 'p7'

    if re.search(r'\bcore\b', text, flags=re.IGNORECASE):
        return 'Core'

    normalized = normalize_lookup_key(text)
    if normalized in {'polyprotein', 'wholegenome'}:
        return 'Whole genome'

    return ''


def canonicalize_product(product_name, alias_lookup):
    product_name = str(product_name or '').strip()
    if not product_name:
        return ''

    canonical = alias_lookup.get(normalize_lookup_key(product_name), product_name)
    inferred = infer_compact_product_name(canonical)
    if inferred:
        return inferred
    return canonical


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


def feature_priority(feature_type, product_name):
    canonical = normalize_lookup_key(product_name)
    if canonical in {'polyprotein', 'wholegenome'}:
        return 99
    priorities = {
        'mat_peptide': 0,
        'mature_protein_region_of_CDS': 0,
        'gene': 1,
        'CDS': 2,
        'cds': 2,
    }
    return priorities.get(feature_type, 10)


def choose_feature_entry(existing_entry, new_entry):
    if existing_entry is None:
        return new_entry
    existing_priority = feature_priority(existing_entry.get('feature_type', ''), existing_entry.get('product', ''))
    new_priority = feature_priority(new_entry.get('feature_type', ''), new_entry.get('product', ''))
    if new_priority < existing_priority:
        return new_entry
    if new_priority == existing_priority:
        existing_span = int(existing_entry['cds_end']) - int(existing_entry['cds_start'])
        new_span = int(new_entry['cds_end']) - int(new_entry['cds_start'])
        if new_span < existing_span:
            return new_entry
    return existing_entry


def make_feature_entry(product_name, start, end, reference_accession, feature_type, alias_lookup, raw_product=None, source='unknown'):
    canonical_product = canonicalize_product(product_name, alias_lookup)
    if not canonical_product:
        return None
    return {
        'product': canonical_product,
        'raw_product': raw_product or product_name,
        'cds_start': int(start),
        'cds_end': int(end),
        'reference_accession': reference_accession,
        'feature_type': feature_type,
        'source': source,
    }


def parse_gff_attributes(attributes_text):
    parsed = {}
    for part in str(attributes_text or '').split(';'):
        if '=' not in part:
            continue
        key, value = part.split('=', 1)
        parsed[key.strip()] = value.strip()
    return parsed


def parse_gff_text_to_feature_map(gff_text, reference_accession, alias_lookup, source_label):
    feature_map = {}
    for line in str(gff_text or '').splitlines():
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) != 9:
            continue
        _, _, feature_type, start, end, _, _, _, attributes = parts
        if feature_type not in {'mat_peptide', 'gene', 'CDS', 'mature_protein_region_of_CDS'}:
            continue
        parsed_attributes = parse_gff_attributes(attributes)
        raw_product = parsed_attributes.get('gene') or parsed_attributes.get('product') or parsed_attributes.get('Name') or parsed_attributes.get('ID')
        entry = make_feature_entry(
            raw_product,
            start,
            end,
            reference_accession,
            feature_type,
            alias_lookup,
            raw_product=raw_product,
            source=source_label,
        )
        if entry is None:
            continue
        key = entry['product']
        feature_map[key] = choose_feature_entry(feature_map.get(key), entry)
    return feature_map


def fetch_reference_gff_map(accession, alias_lookup):
    Entrez.email = os.environ.get('ENTREZ_EMAIL', 'vgtk@example.com')
    handle = Entrez.efetch(db='nuccore', id=accession, rettype='gff3', retmode='text')
    try:
        gff_text = handle.read()
    finally:
        handle.close()
    return parse_gff_text_to_feature_map(gff_text, accession, alias_lookup, 'genbank_gff')


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
    entry = make_feature_entry(
        product_name,
        start,
        end,
        reference_accession,
        feature_type,
        alias_lookup,
        raw_product=raw_product,
        source=source,
    )
    if entry is None:
        return
    reference_map = feature_maps.setdefault(reference_accession, {})
    reference_map[entry['product']] = choose_feature_entry(reference_map.get(entry['product']), entry)


def feature_map_has_gene_information(feature_map):
    return any(
        normalize_lookup_key(product_name) not in {'polyprotein', 'wholegenome'}
        for product_name in (feature_map or {})
    )


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
        elif {'reference_accession', 'gff_text'}.issubset(columns) or {'reference_accession', 'content'}.issubset(columns):
            text_col = 'gff_text' if 'gff_text' in columns else 'content'
            for _, row in df.iterrows():
                reference_accession = str(row.get('reference_accession', '') or '').strip()
                if not reference_accession:
                    continue
                if allowed_accessions and reference_accession not in allowed_accessions:
                    continue
                parsed = parse_gff_text_to_feature_map(row.get(text_col, ''), reference_accession, alias_lookup, f'db_table:{table_name}')
                if not parsed:
                    continue
                reference_map = feature_maps.setdefault(reference_accession, {})
                for product_name, entry in parsed.items():
                    reference_map[product_name] = choose_feature_entry(reference_map.get(product_name), entry)

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


def resolve_reference_feature_map(reference_hint, fallback_candidates, db_gff_maps, alias_lookup, allow_genbank_gff):
    attempted = []
    candidates = []
    for value in [reference_hint] + list(fallback_candidates):
        for token in extract_accession_tokens(value):
            if token and token not in candidates:
                candidates.append(token)

    for candidate in candidates:
        attempted.append(candidate)
        if candidate in db_gff_maps and feature_map_has_gene_information(db_gff_maps[candidate]):
            return candidate, db_gff_maps[candidate], attempted, 'db_gff'

    if not allow_genbank_gff:
        raise AnnotationMappingError(
            f'No DB GFF annotations could be resolved for reference candidates: {attempted}. Re-run with --allow_genbank_reference_gff to fetch GFF from NCBI.'
        )

    last_error = None
    for candidate in candidates:
        try:
            feature_map = fetch_reference_gff_map(candidate, alias_lookup)
            if feature_map:
                return candidate, feature_map, attempted, 'genbank_gff'
        except Exception as exc:
            last_error = exc

    if last_error is not None:
        raise AnnotationMappingError(
            f'Unable to fetch GFF annotations for reference candidates: {attempted}. Last error: {last_error}'
        )
    raise AnnotationMappingError(f'Unable to resolve any GFF annotations from candidates: {attempted}')


def prepare_catalog(catalog, alias_lookup):
    prepared = catalog.copy()
    prepared['_canonical_protein'] = prepared['protein_name'].apply(lambda value: canonicalize_product(value, alias_lookup))
    prepared['_segment_norm'] = prepared['segment'].apply(normalize_segment)
    prepared['_alt_residue_norm'] = prepared['alt_residue'].apply(normalize_residue)

    aa_positions = []
    invalid_positions = 0
    for value in prepared['aa_position']:
        try:
            aa_positions.append(int(float(value)))
        except (TypeError, ValueError):
            aa_positions.append(None)
            invalid_positions += 1
    prepared['_aa_position_int'] = aa_positions

    return prepared, invalid_positions


def validate_required_catalog_columns(catalog):
    missing = [column for column in REQUIRED_MUTATION_CATALOG_COLUMNS if column not in catalog.columns]
    if missing:
        raise ValueError(f"Catalog missing required column(s): {', '.join(missing)}")


def build_catalog_reference_table(catalog, catalog_column_profile):
    if catalog_column_profile == 'all_columns':
        selected_columns = [column for column in catalog.columns if not column.startswith('_')]
    else:
        selected_columns = REQUIRED_MUTATION_CATALOG_COLUMNS + CATALOG_COLUMN_PROFILES[catalog_column_profile]

    catalog_reference = catalog.copy()
    for column in selected_columns:
        if column not in catalog_reference.columns:
            catalog_reference[column] = ''

    catalog_reference = catalog_reference[selected_columns].copy()
    catalog_reference = catalog_reference.fillna('')
    catalog_reference = catalog_reference.drop_duplicates().reset_index(drop=True)
    return catalog_reference




def build_alignment_coordinate_map(alignment):
    coord_map = {}
    ungapped_position = 0
    for alignment_index, base in enumerate(str(alignment or '')):
        if base != '-':
            ungapped_position += 1
            coord_map[ungapped_position] = alignment_index
    return coord_map


def prepare_sequence_alignment(seq_aln):
    seq_aln = seq_aln.copy()
    if 'sequence_id' not in seq_aln.columns:
        if 'primary_accession' in seq_aln.columns:
            seq_aln['sequence_id'] = seq_aln['primary_accession']
        else:
            raise AnnotationMappingError('sequence_alignment must contain sequence_id or primary_accession')
    if 'primary_accession' not in seq_aln.columns:
        seq_aln['primary_accession'] = seq_aln['sequence_id']
    if 'alignment' not in seq_aln.columns:
        raise AnnotationMappingError('sequence_alignment must contain an alignment column')
    if 'alignment_name' not in seq_aln.columns:
        seq_aln['alignment_name'] = ''
    return seq_aln


def find_reference_alignment_row(group_df, reference_id, master_candidates):
    candidates = extract_accession_tokens(reference_id) + list(master_candidates)
    seen = []
    for candidate in candidates:
        if candidate and candidate not in seen:
            seen.append(candidate)

    for candidate in seen:
        mask = (
            group_df['sequence_id'].astype(str).str.strip() == candidate
        ) | (
            group_df['primary_accession'].astype(str).str.strip() == candidate
        )
        if mask.any():
            return group_df.loc[mask].iloc[0], candidate
    return None, None


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


def extract_feature_codon(alignment, coord_map, cds_start, aa_pos):
    alignment_indices = resolve_aligned_codon_indices(coord_map, cds_start, aa_pos)
    if alignment_indices is None:
        return None
    return extract_aligned_codon(alignment, alignment_indices)




def annotate_from_reference_coordinates(catalog, seq_aln, meta_data, alias_lookup, db_gff_maps, allow_genbank_reference_gff, call_evidence=None):
    """Annotate every sequence, gated by genotype scope and per-genotype wild type.

    ``call_evidence``, when a list is supplied, is filled with one record per
    *evaluated* residue match - emitted and suppressed alike - each carrying the
    scope tier it matched on and whether the residue was a change, an anchor
    (equal to the wild type) or wild-type-unknown.  The returned mutation list
    contains only the emitted calls.
    """
    master_candidates = []
    if meta_data is not None and not meta_data.empty and {'primary_accession', 'accession_type'}.issubset(meta_data.columns):
        masters = meta_data[
            meta_data['accession_type'].fillna('').str.strip().str.lower() == 'master'
        ]['primary_accession'].dropna().astype(str).str.strip().unique().tolist()
        master_candidates = masters

    catalog_reference_hints = [value for value in catalog['reference_accession'].dropna().astype(str).str.strip().unique().tolist() if value]
    seq_aln = seq_aln.copy()
    seq_aln['_reference_id'] = seq_aln['alignment_name'].fillna('').astype(str).str.strip()

    if (seq_aln['_reference_id'] == '').any():
        if len(master_candidates) == 1:
            seq_aln.loc[seq_aln['_reference_id'] == '', '_reference_id'] = master_candidates[0]
            print(
                f"[AnnotateMutations][warn] Filled {int((seq_aln['_reference_id'] == master_candidates[0]).sum())} blank alignment_name values with master accession {master_candidates[0]}"
            )
        elif len(master_candidates) == 0 and not catalog_reference_hints:
            dummy_ref = seq_aln['sequence_id'].iloc[0]
            seq_aln.loc[seq_aln['_reference_id'] == '', '_reference_id'] = dummy_ref
            print(f"[AnnotateMutations][warn] Filled blank alignment_name with first sequence accession {dummy_ref}")
        else:
            matched_ref = None
            for hint in catalog_reference_hints:
                for token in extract_accession_tokens(hint):
                    if token in seq_aln['sequence_id'].values:
                        matched_ref = token
                        break
                if matched_ref:
                    break
            if matched_ref:
                seq_aln.loc[seq_aln['_reference_id'] == '', '_reference_id'] = matched_ref
                print(f"[AnnotateMutations][warn] Filled blank alignment_name with matched catalog reference accession {matched_ref}")
            elif len(seq_aln['sequence_id'].unique()) == 1:
                single_acc = seq_aln['sequence_id'].iloc[0]
                seq_aln.loc[seq_aln['_reference_id'] == '', '_reference_id'] = single_acc
                print(f"[AnnotateMutations][warn] Filled blank alignment_name with unique sequence accession {single_acc}")
            else:
                first_acc = seq_aln['sequence_id'].iloc[0]
                seq_aln.loc[seq_aln['_reference_id'] == '', '_reference_id'] = first_acc
                print(f"[AnnotateMutations][warn] Filled blank alignment_name with default first sequence accession {first_acc}")

    if seq_aln.empty:
        raise AnnotationMappingError('No usable reference identifiers were found in sequence_alignment.alignment_name.')

    diagnostics = Counter()
    mutations_found = []

    # -----------------------------------------------------------------------
    # Resolve the master reference's alignment and feature map.
    #
    # All catalog positions (aa_position) are defined relative to the master
    # reference (e.g. NC_004102 for HCV).  Because every sequence in the DB
    # is padded to the SAME alignment width as the master, the master's column
    # positions apply universally.  We must NOT use per-genotype-reference OG
    # coordinates to resolve catalog positions — those produce shifted column
    # indices for sequences aligned to divergent genotype references (e.g.
    # D10988 NS3 cds_start_OG=3345 vs NC_004102 cds_start=3420 shifts all
    # downstream column positions by ~71 bp relative to the master numbering).
    # -----------------------------------------------------------------------
    coord_candidates = []
    for cand in list(master_candidates) + catalog_reference_hints:
        for token in extract_accession_tokens(cand):
            if token and token not in coord_candidates:
                coord_candidates.append(token)

    master_coord_ref_acc = None
    master_coord_map = None
    master_feature_map = None

    for cand in coord_candidates:
        mask = (
            seq_aln['sequence_id'].astype(str).str.strip() == cand
        ) | (
            seq_aln['primary_accession'].astype(str).str.strip() == cand
        )
        if not mask.any():
            continue
        cand_coord_map = build_alignment_coordinate_map(seq_aln.loc[mask].iloc[0]['alignment'])
        for token in extract_accession_tokens(cand):
            if token in db_gff_maps and feature_map_has_gene_information(db_gff_maps[token]):
                master_coord_ref_acc = cand
                master_coord_map = cand_coord_map
                master_feature_map = db_gff_maps[token]
                break
        if master_coord_map is not None:
            break

    if master_coord_map is None or master_feature_map is None:
        # Fallback: try any reference group row with a usable feature map
        for ref_id in seq_aln['_reference_id'].dropna().unique():
            mask = seq_aln['sequence_id'].astype(str).str.strip() == str(ref_id).strip()
            if not mask.any():
                continue
            for token in extract_accession_tokens(str(ref_id)):
                if token in db_gff_maps and feature_map_has_gene_information(db_gff_maps[token]):
                    master_coord_ref_acc = str(ref_id)
                    master_coord_map = build_alignment_coordinate_map(seq_aln.loc[mask].iloc[0]['alignment'])
                    master_feature_map = db_gff_maps[token]
                    break
            if master_coord_map is not None:
                break

    if master_coord_map is None or master_feature_map is None:
        if not allow_genbank_reference_gff:
            raise AnnotationMappingError(
                'Could not resolve a master reference alignment + feature map. '
                f'Candidates tried: {coord_candidates}. Re-run with --allow_genbank_reference_gff to fetch from NCBI.'
            )
        # Try fetching from NCBI
        for cand in coord_candidates:
            mask = (
                seq_aln['sequence_id'].astype(str).str.strip() == cand
            ) | (
                seq_aln['primary_accession'].astype(str).str.strip() == cand
            )
            if not mask.any():
                continue
            try:
                fetched_map = fetch_reference_gff_map(cand, alias_lookup)
                if fetched_map:
                    master_coord_ref_acc = cand
                    master_coord_map = build_alignment_coordinate_map(seq_aln.loc[mask].iloc[0]['alignment'])
                    master_feature_map = fetched_map
                    break
            except Exception:
                continue
        if master_coord_map is None or master_feature_map is None:
            raise AnnotationMappingError(
                'Could not resolve a master reference alignment + feature map. '
                f'Candidates tried: {coord_candidates}'
            )

    print(
        f'[AnnotateMutations] Coordinate reference: {master_coord_ref_acc} '
        f"(all catalog positions resolved using this sequence's alignment column space)"
    )

    # -----------------------------------------------------------------------
    # Build a single catalog position lookup using universal alignment column
    # indices derived from the master reference coordinate system.
    # -----------------------------------------------------------------------
    grouped_positions = {}   # (protein_name, aa_pos, alignment_indices_tuple) -> {alt_residue: [catalog_rows]}
    proteins_missing_from_reference = Counter()
    mappable_catalog_rows = 0

    for _, row in catalog.iterrows():
        protein_name = row['_canonical_protein']
        aa_pos = row['_aa_position_int']
        if not protein_name or aa_pos is None:
            diagnostics['invalid_catalog_rows'] += 1
            continue

        gene_entry = master_feature_map.get(protein_name)
        if gene_entry is None:
            proteins_missing_from_reference[protein_name] += 1
            continue

        alignment_indices = resolve_aligned_codon_indices(
            master_coord_map,
            gene_entry['cds_start'],
            aa_pos,
        )
        if alignment_indices is None:
            diagnostics['reference_coordinate_missing_in_alignment'] += 1
            continue

        position_catalog = grouped_positions.setdefault(
            (protein_name, aa_pos, tuple(alignment_indices)), {}
        )
        position_catalog.setdefault(row['_alt_residue_norm'], []).append(row)
        mappable_catalog_rows += 1

    if mappable_catalog_rows == 0:
        missing_summary = ', '.join(
            protein for protein, _ in proteins_missing_from_reference.most_common(10)
        )
        raise AnnotationMappingError(
            'No catalog coordinates could be mapped onto master reference CDS annotations. '
            f'Proteins not found in master feature map: {missing_summary}'
        )

    # -----------------------------------------------------------------------
    # Annotate every sequence using the master-derived universal column
    # positions.  Because all padded alignments share the same column space,
    # these indices are valid for every row regardless of alignment_name.
    # -----------------------------------------------------------------------
    signature_scope = build_signature_genotype_scope(catalog)
    wild_type_tables = build_wild_type_tables(catalog)
    sequence_genotypes = build_sequence_genotype_map(meta_data)

    for _, seq_row in seq_aln.iterrows():
        alignment = seq_row['alignment']
        primary_accession = seq_row['sequence_id']
        genotype, subtype_code = sequence_genotypes.get(str(primary_accession).strip(), ('', ''))
        if not genotype:
            diagnostics['sequences_without_genotype'] += 1
        covered_span = alignment_covered_span(alignment)
        for (protein_name, aa_pos, alignment_indices), alt_lookup in grouped_positions.items():
            codon = extract_aligned_codon(alignment, alignment_indices)
            if codon is None:
                diagnostics['codon_out_of_bounds'] += 1
                continue
            aa = residue_from_aligned_codon(codon, covered_span, alignment_indices)
            matched_rows = alt_lookup.get(aa, [])
            if not matched_rows:
                continue

            wild_type_residues, wild_type_tier = lookup_wild_type_residues(
                wild_type_tables, subtype_code, genotype, protein_name, aa_pos
            )
            for row in matched_rows:
                signature_id = clean_cell(row.get('signature_id', ''))
                scope_codes = signature_scope.get(signature_id, frozenset())
                scope_tier = classify_genotype_scope(scope_codes, genotype, subtype_code)

                if wild_type_residues is None:
                    residue_status = RESIDUE_STATUS_WT_UNKNOWN
                elif aa in wild_type_residues:
                    residue_status = RESIDUE_STATUS_ANCHOR
                else:
                    residue_status = RESIDUE_STATUS_CHANGE

                if scope_tier == SCOPE_TIER_OUT_OF_SCOPE:
                    call_status = CALL_STATUS_SUPPRESSED_OUT_OF_SCOPE
                elif residue_status == RESIDUE_STATUS_ANCHOR:
                    call_status = CALL_STATUS_SUPPRESSED_WILD_TYPE
                else:
                    call_status = CALL_STATUS_EMITTED

                record = {
                    'primary_accession': primary_accession,
                    'mutation_id': row['mutation_id'],
                    'protein_name': protein_name,
                    'segment': row['_segment_norm'],
                    'aa_position': aa_pos,
                    'alt_residue': row['_alt_residue_norm'],
                    'combination_id': clean_cell(row.get('combination_id', '')),
                    'signature_id': signature_id,
                    'signature_kind': clean_cell(row.get('signature_kind', '')),
                    'observed_residue': aa,
                    'sequence_genotype': genotype,
                    'sequence_subtype': subtype_code,
                    'relevant_genotypes': ','.join(sorted(scope_codes)),
                    'scope_tier': scope_tier,
                    'wild_type_residues': ''.join(sorted(wild_type_residues)) if wild_type_residues else '',
                    'wild_type_scope_tier': wild_type_tier,
                    'residue_status': residue_status,
                    'call_status': call_status,
                }
                if call_evidence is not None:
                    call_evidence.append(record)

                if call_status != CALL_STATUS_EMITTED:
                    diagnostics[call_status] += 1
                    continue

                mutations_found.append(record)
                diagnostics['mutation_hits'] += 1
                diagnostics[f'emitted_{scope_tier}'] += 1
                diagnostics[f'emitted_{residue_status}'] += 1
                if aa == DELETION_RESIDUE:
                    diagnostics['emitted_deletions'] += 1

    if proteins_missing_from_reference:
        preview = ', '.join(
            protein for protein, _ in proteins_missing_from_reference.most_common(10)
        )
        print(f'[AnnotateMutations][warn] Missing protein annotations in master feature map: {preview}')

    resolved_maps = {
        master_coord_ref_acc: {
            'resolved_accession': master_coord_ref_acc,
            'feature_map': master_feature_map,
            'attempted': coord_candidates,
            'source': 'master_coord_ref',
            'aligned_reference_accession': master_coord_ref_acc,
            'aligned_reference_map': master_coord_map,
        }
    }

    return mutations_found, diagnostics, resolved_maps

MUTATION_CALL_IDENTITY_COLUMNS = [
    'primary_accession',
    'mutation_id',
    'protein_name',
    'segment',
    'aa_position',
    'alt_residue',
    'combination_id',
]

MUTATION_CALL_PROVENANCE_COLUMNS = [
    'signature_id',
    'signature_kind',
    'observed_residue',
    'sequence_genotype',
    'sequence_subtype',
    'relevant_genotypes',
    'scope_tier',
    'wild_type_residues',
    'wild_type_scope_tier',
    'residue_status',
    'call_status',
]

MUTATION_CALL_COLUMNS = MUTATION_CALL_IDENTITY_COLUMNS + MUTATION_CALL_PROVENANCE_COLUMNS


def build_mutation_call_table(records):
    """One row per evaluated residue match, best scope tier first.

    Keeping the *reason* next to the call is the point: a reader can tell a
    genuine change from an entry that only fired because no wild type was known
    for that genotype, and a suppressed anchor from one that was never in scope.
    """
    df_calls = pd.DataFrame(records, columns=MUTATION_CALL_COLUMNS)
    if df_calls.empty:
        return df_calls
    df_calls = df_calls.fillna('')
    df_calls['_scope_rank'] = df_calls['scope_tier'].map(SCOPE_TIER_RANK).fillna(len(SCOPE_TIER_RANK))
    df_calls = df_calls.sort_values(
        MUTATION_CALL_IDENTITY_COLUMNS + ['_scope_rank', 'signature_id'], kind='mergesort'
    )
    return df_calls.drop(columns=['_scope_rank']).drop_duplicates().reset_index(drop=True)


def write_mutation_tables(conn, catalog, mutations_found, catalog_column_profile, call_evidence=None,
                          publications=None, clinical_trials=None):
    # The compact layout tables are keyed on the identity columns alone, so the
    # provenance columns are kept out of them and land in
    # sequence_mutation_calls instead - joining on a wider frame would change
    # which columns the layout builders collide on.
    df_mut = pd.DataFrame(mutations_found, columns=MUTATION_CALL_IDENTITY_COLUMNS)
    if not df_mut.empty:
        df_mut = df_mut.drop_duplicates(subset=MUTATION_CALL_IDENTITY_COLUMNS)
        df_mut = df_mut.fillna('')

    cursor = conn.cursor()

    print('Writing mutation_catalog table...')
    df_catalog = build_catalog_reference_table(catalog, catalog_column_profile)
    cursor.execute('DROP TABLE IF EXISTS mutation_catalog')
    df_catalog.to_sql('mutation_catalog', conn, if_exists='replace', index=False)
    if 'mutation_id' in df_catalog.columns:
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_mut_catalog_mutid ON mutation_catalog(mutation_id)')
    if 'protein_name' in df_catalog.columns and 'segment' in df_catalog.columns:
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_mut_catalog_seg_prot ON mutation_catalog(segment, protein_name)')
    if 'combination_id' in df_catalog.columns:
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_mut_catalog_comb ON mutation_catalog(combination_id)')

    if publications is not None and not publications.empty:
        print(f'Writing publications table ({len(publications)} rows)...')
        cursor.execute('DROP TABLE IF EXISTS publications')
        publications.to_sql('publications', conn, if_exists='replace', index=False)
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_publications_pmid ON publications(pubmed_id)')

    if clinical_trials is not None and not clinical_trials.empty:
        print(f'Writing clinical_trials table ({len(clinical_trials)} rows)...')
        cursor.execute('DROP TABLE IF EXISTS clinical_trials')
        clinical_trials.to_sql('clinical_trials', conn, if_exists='replace', index=False)
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_clinical_trials_nct ON clinical_trials(nct_id)')

    print('Writing compact mutation summary tables...')
    df_catalog_for_layouts = catalog.fillna('').copy()
    df_relevant_summary = build_sequence_relevant_mutation_summary(df_mut)
    df_completed_signatures = build_completed_signatures_only(df_mut, df_catalog_for_layouts)

    cursor.execute('DROP TABLE IF EXISTS sequence_mutations')
    cursor.execute('DROP TABLE IF EXISTS sequence_drug_resistance')

    print('Writing per-call genotype scope / wild-type evidence...')
    df_calls = build_mutation_call_table(mutations_found if call_evidence is None else call_evidence)
    cursor.execute('DROP TABLE IF EXISTS sequence_mutation_calls')
    if not df_calls.empty:
        df_calls.to_sql('sequence_mutation_calls', conn, if_exists='replace', index=False)
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_seq_mut_calls_acc ON sequence_mutation_calls(primary_accession)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_seq_mut_calls_status ON sequence_mutation_calls(call_status)')
    else:
        cursor.execute(
            'CREATE TABLE sequence_mutation_calls ('
            + ', '.join(f'{column} TEXT' for column in MUTATION_CALL_COLUMNS)
            + ')'
        )

    cursor.execute('DROP TABLE IF EXISTS sequence_relevant_mutation_summary')
    if not df_relevant_summary.empty:
        df_relevant_summary.to_sql('sequence_relevant_mutation_summary', conn, if_exists='replace', index=False)
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_seq_rel_mut_acc ON sequence_relevant_mutation_summary(primary_accession)')
    else:
        cursor.execute('CREATE TABLE sequence_relevant_mutation_summary (primary_accession TEXT, relevant_mutations_present TEXT, total_relevant_mutation_count INTEGER)')

    cursor.execute('DROP TABLE IF EXISTS completed_signatures_only')
    if not df_completed_signatures.empty:
        df_completed_signatures.to_sql('completed_signatures_only', conn, if_exists='replace', index=False)
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_completed_sig_acc_sig ON completed_signatures_only(primary_accession, signature_id)')
    else:
        cursor.execute('CREATE TABLE completed_signatures_only (primary_accession TEXT, signature_id TEXT, signature_kind TEXT)')

    conn.commit()

def load_clinical_trials_table(clinical_trial_path):
    """The trial registry entries the catalogue's NCT identifiers refer to.

    ``mutation_catalog.clinical_trials`` holds semicolon-separated NCT numbers,
    scoped per (mutation, genotype, drug) - trial support genuinely varies by
    genotype, so NS5A:31M against daclatasvir cites one trial in genotype 1a and
    nine in 1b. This loads the registry rows so an NCT resolves to the trial's
    name rather than staying an opaque accession.

    Optional. Without it the NCT identifiers are still present and still
    correct, just not resolvable inside the database.
    """
    if not clinical_trial_path or not os.path.isfile(clinical_trial_path):
        return None
    rows = []
    with open(clinical_trial_path, newline='', encoding='utf-8', errors='replace') as handle:
        for row in csv.DictReader(handle):
            nct = (row.get('nct_id') or '').strip()
            if not nct:
                continue
            rows.append({
                'nct_id': nct,
                'trial_id': (row.get('id') or '').strip(),
                'trial_name': (row.get('display_name') or '').strip(),
            })
    return pd.DataFrame(rows) if rows else None


def load_publications_table(publications_path):
    """Read the publication metadata that PMIDs in the catalogue refer to.

    ``mutation_catalog.pubmed_id`` holds semicolon-separated PubMed IDs and is
    already genotype-scoped - the same signature and drug cites different
    publications in different genotypes, which survived the catalogue build
    intact. But a bare PMID is not usable evidence on its own: nothing in the
    database says what 27773808 is.

    This loads the 128 rows of title / authors / year / journal / url so a PMID
    resolves to something a reader can act on, and flags the ones whose titles
    identify them as clinical trials.

    Optional. Without it the database is exactly as it was - PMIDs present,
    unresolvable.
    """
    if not publications_path or not os.path.isfile(publications_path):
        return None
    rows = []
    with open(publications_path, newline='', encoding='utf-8', errors='replace') as handle:
        for row in csv.DictReader(handle):
            title = (row.get('title') or '').strip()
            rows.append({
                'pubmed_id': (row.get('id') or '').strip(),
                'title': title,
                'authors': (row.get('authors_short') or '').strip(),
                'year': (row.get('year') or '').strip(),
                'journal': (row.get('journal') or '').strip(),
                'url': (row.get('url') or '').strip(),
            })
    return pd.DataFrame(rows) if rows else None


def main():
    parser = argparse.ArgumentParser(description="Annotate mutations and drug resistance.")
    parser.add_argument("--db", required=True, help="Path to SQLite database.")
    parser.add_argument("--mutation_catalog", required=True, help="Path to mutation catalog TSV.")
    parser.add_argument("--virus", default="", help="Virus context for specific logics (e.g. HCV)")
    parser.add_argument("--publications", default=None,
        help="Optional publication metadata CSV (id/title/authors_short/year/journal/url). "
             "Loaded into a publications table so the PMIDs already in mutation_catalog.pubmed_id "
             "resolve to something readable. Without it those PMIDs stay bare numbers.")
    parser.add_argument("--clinical_trials", default=None,
        help="Optional clinical trial registry CSV (id/display_name/nct_id). Loaded into a "
             "clinical_trials table so the NCT ids in mutation_catalog.clinical_trials resolve "
             "to trial names.")
    parser.add_argument(
        "--catalog_column_profile",
        required=True,
        choices=['HCV', 'influenza', 'all_columns'],
        help="Select which mutation catalog columns are written to the database.",
    )
    parser.add_argument(
        "--allow_genbank_reference_gff",
        action="store_true",
        help="Allow fallback to fetch reference GFF annotations from NCBI when no usable annotated gene features or DB GFF annotations are available.",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.db):
        print(f"Error: Database {args.db} not found.", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.mutation_catalog):
        print(f"Error: Catalog {args.mutation_catalog} not found.", file=sys.stderr)
        sys.exit(1)

    print("Loading mutation catalog...")
    catalog = pd.read_csv(args.mutation_catalog, sep='\t', dtype=str)

    try:
        validate_required_catalog_columns(catalog)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    conn = sqlite3.connect(args.db)

    try:
        print('Loading sequence alignments...')
        seq_aln = prepare_sequence_alignment(pd.read_sql_query('SELECT * FROM sequence_alignment', conn))
        meta_data = None
        try:
            meta_data = pd.read_sql_query('SELECT * FROM meta_data', conn)
        except Exception:
            pass

        gene_alias_lookup = load_gene_alias_lookup(conn)
        catalog, invalid_positions = prepare_catalog(catalog, gene_alias_lookup)
        if invalid_positions:
            print(f'[AnnotateMutations][warn] Skipping {invalid_positions} catalog rows with invalid aa_position values')

        db_gff_maps = load_db_gff_feature_maps(conn, gene_alias_lookup)

        print('Extracting mutations...')
        call_evidence = []
        mutations_found, diagnostics, resolved_maps = annotate_from_reference_coordinates(
            catalog,
            seq_aln[['sequence_id', 'primary_accession', 'alignment', 'alignment_name']].copy(),
            meta_data,
            gene_alias_lookup,
            db_gff_maps,
            args.allow_genbank_reference_gff,
            call_evidence=call_evidence,
        )
        print(
            '[AnnotateMutations] Resolved reference coordinate maps for: '
            + ', '.join(
                f'{ref}->{info["resolved_accession"]} ({info["source"]}; aligned_ref={info["aligned_reference_accession"]})'
                for ref, info in resolved_maps.items()
            )
        )

        print('[AnnotateMutations] Mapping summary: ' + ', '.join(f'{key}={value}' for key, value in sorted(diagnostics.items())))

        if not mutations_found and diagnostics.get('mutation_hits', 0) == 0:
            if diagnostics.get('feature_rows_matched', 0) == 0 and diagnostics.get('reference_groups_without_positions', 0) > 0:
                raise AnnotationMappingError('No mutations were annotated because coordinate mapping failed for all available reference groups.')
            print('[AnnotateMutations][warn] No mutation hits were found after successful coordinate mapping')

        write_mutation_tables(conn, catalog, mutations_found, args.catalog_column_profile,
                              call_evidence=call_evidence,
                              publications=load_publications_table(args.publications),
                              clinical_trials=load_clinical_trials_table(args.clinical_trials))
    except AnnotationMappingError as exc:
        print(f'Error: {exc}', file=sys.stderr)
        conn.close()
        sys.exit(2)

    conn.close()
    print('Done.')

if __name__ == "__main__":
    main()
