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

#: ---------------------------------------------------------------------------
#: Virus profiles.
#:
#: Everything in this module that is not true of every virus lives here, keyed
#: by virus name.  ``--virus`` selects one.  It used to be declared, documented
#: as "Virus context for specific logics (e.g. HCV)", and then never read: the
#: HCV protein-name patterns below were applied to rabies and influenza alike,
#: so a GenBank product description containing the word 'core' or an 'E1' token
#: was rewritten into HCV's mature-peptide vocabulary.  Both sides of the
#: protein join go through canonicalize_product(), so the two sides then
#: disagreed and the rows dropped out silently.
#:
#: ``columns`` is the EXTRA catalogue columns written to the mutation_catalog
#: table on top of REQUIRED_MUTATION_CATALOG_COLUMNS.  None means "every column
#: the catalogue supplies", which is the default for a virus with no entry
#: here: keeping what the curator wrote is the only choice that cannot lose
#: data, and inventing empty HCV drug-resistance columns is not a prerequisite
#: for annotating a virus that has no drugs.
#:
#: ``product_patterns`` are tried in order against a feature product name or a
#: catalogue protein name; the first to match supplies the compact name.
#: Anything not listed is left exactly as written.
GENERIC_VIRUS = 'generic'

VIRUS_PROFILES = {
    GENERIC_VIRUS: {
        'columns': None,
        'product_patterns': [],
    },
    'hcv': {
        'columns': [
            'resistance_category', 'drug', 'drug_category', 'drug_producer', 'pubmed_id',
            'DOI', 'any_in_vitro_evidence', 'in_vitro_max_ec50_midpoint',
            'any_in_vivo_evidence', 'in_vivo_baseline', 'in_vivo_treatment_emergent',
        ],
        # The ten HCV mature peptides. 'core' and 'E1'/'E2' are the reason this
        # list must not be applied to another virus: they are ordinary English
        # in a great many non-HCV product descriptions.
        'product_patterns': [
            (r'\b(NS[2-5](?:A|B)?)\b', 'upper'),
            (r'\b(E[12])\b', 'upper'),
            (r'\b(p7)\b', 'p7'),
            (r'\b(core)\b', 'Core'),
        ],
    },
    # generic/rabv/ is the repository's own namesake build. It has no
    # drug-resistance vocabulary and no compact-name inference to do: its
    # proteins (N, P, M, G, L) are already the names GenBank uses. Registering
    # it means `--virus rabies` resolves instead of warning, which is the whole
    # difference between an onboarded virus and a tolerated one.
    'rabv': {
        'columns': None,
        'product_patterns': [],
    },
    'other': {
        'columns': None,
        'product_patterns': [],
    },
    'influenza': {
        'columns': ['serotypes_tested'],
        # Canonical IAV protein names, as spelled by the `genes` table of a
        # real influenza build.  PA-X is listed before PA so the longer name
        # wins; regex alternation is left-biased.
        'product_patterns': [
            (r'\b(PB1-F2|PB2|PB1|PA-X|PA|HA|NP|NA|M1|M2|NS1|NS2|NEP)\b', 'upper'),
        ],
    },
}

#: Aliases so a caller may name the virus the way its assets and profiles do.
VIRUS_ALIASES = {
    'hepatitis_c_virus': 'hcv',
    'hepatitis_c': 'hcv',
    'hepc': 'hcv',
    'flu': 'influenza',
    'iav': 'influenza',
    'influenza_a': 'influenza',
    'influenza_a_virus': 'influenza',
    'rabies': 'rabv',
    'rabies_virus': 'rabv',
    'all_columns': GENERIC_VIRUS,
}


def resolve_virus_name(*candidates):
    """The first candidate naming a virus profile, else the generic profile.

    Both pipeline call sites derive --catalog_column_profile FROM the virus
    name, so the two arguments carry the same information and either may stand
    in for the other.  An unrecognised name is not an error: it resolves to the
    generic profile, which keeps every catalogue column and applies no
    virus-specific protein-name inference.  Before this, a virus the argparse
    `choices` list had never heard of - 'rabies', say - failed at argument
    parsing, which is why the two shipped non-HCV paths were unreachable.
    """
    for candidate in candidates:
        key = clean_cell(candidate).lower().replace('-', '_').replace(' ', '_')
        if not key:
            continue
        key = VIRUS_ALIASES.get(key, key)
        if key in VIRUS_PROFILES:
            return key
    return GENERIC_VIRUS


def virus_profile(virus_name):
    """The profile dict for a virus name, generic when it is not one we know."""
    return VIRUS_PROFILES.get(resolve_virus_name(virus_name), VIRUS_PROFILES[GENERIC_VIRUS])


class GeneAliasLookup(dict):
    """The gene alias map, carrying the virus profile that goes with it.

    A plain dict everywhere it is read, so every existing ``.get()`` call site
    and every test that hands in ``{}`` keeps working.  Carrying the profile
    on the object is what lets canonicalize_product() know which virus it is
    canonicalising for without threading a parameter through the six functions
    between main() and it.
    """

    def __init__(self, aliases=None, profile=None):
        super().__init__(aliases or {})
        self.virus_profile = profile or VIRUS_PROFILES[GENERIC_VIRUS]


def profile_of_lookup(alias_lookup):
    """The virus profile carried by an alias lookup, generic if it carries none."""
    return getattr(alias_lookup, 'virus_profile', None) or VIRUS_PROFILES[GENERIC_VIRUS]


#: Retained so ``--catalog_column_profile`` keeps naming the same column sets
#: it always did; the values now come from the virus registry.
CATALOG_COLUMN_PROFILES = {
    'HCV': VIRUS_PROFILES['hcv']['columns'],
    'influenza': VIRUS_PROFILES['influenza']['columns'],
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
    """Translate a bare codon, answering 'X' for anything that is not one.

    The width is checked on the text as given, *before* the gaps come out.  A
    codon spread over four alignment columns ('AT-G') is three bases the
    aligner placed either side of an insertion; collapsing the gap and reading
    ATG would invent a residue the sequence does not carry.  Only a codon
    exactly three columns wide is translated.

    An all-gap codon still answers 'X' here.  The deletion reading belongs to
    residue_from_aligned_codon(), which can see the surrounding alignment and
    can therefore tell a deleted codon from an unsequenced one.
    """
    text = str(codon or '').upper()
    if len(text) != 3:
        return 'X'
    return CODON_TABLE.get(text.replace('-', ''), 'X')


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

#: What translate_codon() answers when it cannot name a residue.  It collapses
#: three very different situations - an unsequenced codon ('NNN'), a genuine
#: IUPAC ambiguity ('RGA') and a codon the alignment could not assemble - onto
#: one token, so it can only ever mean "we do not know", never "the residue is
#: X".  Nothing is allowed to satisfy a catalogue row with it.
UNKNOWN_RESIDUE = 'X'

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


def coerce_coordinate(value):
    """A 1-based nucleotide coordinate from an untyped cell, or None.

    ``features.cds_start`` and ``gff_features.start`` are TEXT in the shipped
    schema and are filled by upstream stages from GenBank records of every
    quality, so a NULL, a blank, or a REAL written back as '1.0' all arrive
    here verbatim.  int() raises on each of them and no caller has a per-row
    handler, so one unusable cell would abort the run and take every
    well-formed feature in the same table with it.  Answering None costs the
    caller a single feature, which it counts and reports.
    """
    text = clean_cell(value)
    if not text:
        return None
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    if number != number or number in (float('inf'), float('-inf')):
        return None
    return int(number)


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



def normalize_genotype_code(value):
    """Fold a genotype/subtype code onto one spelling before it is compared.

    Residues are already case-folded by normalize_residue() on the reasoning
    that "a curator typo cannot silently disable a catalog row"; the genotype
    vocabulary had no such treatment and is compared in three separate places.
    A meta_data row typed '1A' against a catalogue bucket '1a' therefore missed
    the exact-subtype tier, and rule B fell back to the union of every subtype
    of genotype 1 - or, for a vocabulary with no genotype tier at all, to no
    wild type whatsoever, which emits the residue as a change.

    Lower case is the canonical form because that is what every code in the
    shipped HCV catalogue and database already uses, so folding moves nothing.
    """
    return clean_cell(value).lower()


def genotype_of_code(value):
    """Leading digits of a genotype/subtype code: ``6xd`` -> ``6``, ``1a`` -> ``1``."""
    match = re.match(r'\d+', clean_cell(value))
    return match.group(0) if match else ''


def parent_genotype_of_code(value):
    """The genotype a subtype code belongs to, for any vocabulary.

    HCV spells a subtype as its genotype's digits plus a suffix, so the parent
    of '6xd' is '6' and matching at that level is not a nicety - the observed
    subtypes 6n, 6xd, 6xc and 4v have no catalogue bucket of their own, and
    exact matching would throw away nearly every call.

    A code with no leading digits IS its own parent.  The leading-digit rule
    was applied unconditionally even though both columns that use it are
    documented as generic and virus-agnostic, so every vocabulary that is not
    digit-led collapsed onto the empty string: influenza 'H1', rabies
    'Africa-2', SARS-CoV-2 'BA.2'.  classify_genotype_scope() then found no
    genotype tier and answered out_of_scope, suppressing every call for that
    sequence build-wide, and build_wild_type_tables() filed every code under
    one '' key that lookup_wild_type_residues() never asks for.  Falling back
    to the code itself costs HCV nothing - every code in the shipped catalogue
    and database is digit-led - and gives every other virus the tier back.
    """
    return genotype_of_code(value) or normalize_genotype_code(value)


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
        # Key on the PARSED position, exactly as the annotation loop will ask
        # for it.  Keying on the raw text let a catalogue that respelled '2' as
        # '2.0' - which prepare_catalog() absorbs for coordinates, so the call
        # still fires - build a table the lookup could never hit.  Rule B then
        # went quiet without a word and the genotype's own wild-type residue
        # was emitted as a resistance call, labelled 'wt_unknown' and so
        # indistinguishable from a genuinely uncurated position.
        position = parse_aa_position(row.get('aa_position', ''))
        if not protein or position is None:
            continue
        for code, fields in parse_genotype_entry_list(row.get(WILD_TYPE_RESIDUES_COLUMN, '')).items():
            if not fields:
                continue
            residue = normalize_residue(fields[0])
            if not residue:
                continue
            by_subtype[(normalize_genotype_code(code), protein, str(position))].add(residue)
            by_genotype[(parent_genotype_of_code(code), protein, str(position))].add(residue)
    return dict(by_subtype), dict(by_genotype)


def lookup_wild_type_residues(wild_type_tables, subtype_code, genotype, protein, position):
    """Wild type for this sequence's genotype: exact subtype first, then genotype.

    Returns ``(residues, tier)``; ``(None, '')`` when no wild type is known,
    which must NOT be read as 'the residue is a change' - it is 'we do not
    know', and the caller records it as such rather than suppressing.
    """
    by_subtype, by_genotype = wild_type_tables
    subtype_key = normalize_genotype_code(subtype_code)
    if subtype_key:
        residues = by_subtype.get((subtype_key, protein, str(position)))
        if residues:
            return residues, SCOPE_TIER_SUBTYPE
    if genotype:
        residues = by_genotype.get((parent_genotype_of_code(genotype), protein, str(position)))
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
        scope[signature_id] |= {
            normalize_genotype_code(code)
            for code in parse_relevant_genotypes(row.get(RELEVANT_GENOTYPES_COLUMN, ''))
        }
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
    subtype_key = normalize_genotype_code(subtype_code)
    if subtype_key and subtype_key in scope_codes:
        return SCOPE_TIER_SUBTYPE
    genotype_key = parent_genotype_of_code(genotype)
    if genotype_key and any(parent_genotype_of_code(code) == genotype_key for code in scope_codes):
        return SCOPE_TIER_GENOTYPE
    return SCOPE_TIER_OUT_OF_SCOPE


def build_sequence_genotype_map(meta_data, diagnostics=None):
    """``{accession: (genotype, subtype_code)}`` from meta_data, tolerant of schema drift.

    An accession that appears more than once with DIFFERENT typing keeps the
    first reading and the disagreement is counted and warned about.  It used to
    be resolved by row order alone: meta_data carries one row per accession per
    segment in a segmented build, and a cross-run merge can leave a stale row
    beside a current one, so whichever pandas yielded last drove both rule A
    and rule B for that sequence.  The same database read in a different order
    produced different resistance calls, with nothing recorded to say a
    conflict had ever existed.
    """
    genotype_map = {}
    conflicts = []
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
        resolved = (genotype, subtype_code)
        existing = genotype_map.get(accession)
        if existing is not None:
            if existing != resolved:
                conflicts.append((accession, existing, resolved))
            continue
        genotype_map[accession] = resolved

    if conflicts:
        if diagnostics is not None:
            diagnostics['meta_data_genotype_conflicts'] += len(conflicts)
        preview = ', '.join(
            f'{accession} kept {kept} ignored {dropped}'
            for accession, kept, dropped in conflicts[:5]
        )
        print(
            f'[AnnotateMutations][warn] {len(conflicts)} meta_data row(s) disagree with an '
            f'earlier row about the same accession; the first reading was kept: {preview}'
        )
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


def infer_compact_product_name(product_name, profile=None):
    """The compact protein name a product description implies, or ''.

    The patterns come from the virus profile, so a name is only rewritten into
    a vocabulary the virus actually has.  The polyprotein / whole-genome
    collapse below is applied for every virus because it is about the SHAPE of
    the annotation - a feature covering the whole coding region - rather than
    about any particular protein.
    """
    text = str(product_name or '').strip()
    if not text:
        return ''

    if profile is None:
        profile = VIRUS_PROFILES[GENERIC_VIRUS]

    for pattern, mode in profile.get('product_patterns', ()):
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return match.group(1).upper() if mode == 'upper' else mode

    normalized = normalize_lookup_key(text)
    if normalized in {'polyprotein', 'wholegenome'}:
        return 'Whole genome'

    return ''


def canonicalize_product(product_name, alias_lookup, profile=None):
    """One spelling for a protein, from a product description or a catalog cell.

    Applied to BOTH sides of the protein join - the catalogue's protein_name
    and every feature's product - so the two only meet if they canonicalise the
    same way.  `profile` selects the virus vocabulary; when it is not given the
    lookup is asked for the one it carries, and a bare dict yields the generic
    profile, which renames nothing.
    """
    product_name = str(product_name or '').strip()
    if not product_name:
        return ''

    if profile is None:
        profile = profile_of_lookup(alias_lookup)

    canonical = alias_lookup.get(normalize_lookup_key(product_name), product_name)
    inferred = infer_compact_product_name(canonical, profile)
    if inferred:
        return inferred
    return canonical


def load_gene_alias_lookup(conn, profile=None):
    alias_lookup = GeneAliasLookup(profile=profile)
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
        existing_span = existing_entry['cds_end'] - existing_entry['cds_start']
        new_span = new_entry['cds_end'] - new_entry['cds_start']
        if new_span < existing_span:
            return new_entry
    return existing_entry


def make_feature_entry(product_name, start, end, reference_accession, feature_type, alias_lookup, raw_product=None, source='unknown'):
    canonical_product = canonicalize_product(product_name, alias_lookup)
    if not canonical_product:
        return None
    cds_start = coerce_coordinate(start)
    cds_end = coerce_coordinate(end)
    if cds_start is None or cds_end is None:
        return None
    return {
        'product': canonical_product,
        'raw_product': raw_product or product_name,
        'cds_start': cds_start,
        'cds_end': cds_end,
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
    """Fold one feature row into the map. Returns False if the row was unusable.

    The return value exists so a caller can tally what it dropped: a feature
    silently missing from the map resurfaces much later as a whole protein
    'not found in master feature map', which is a far harder thing to trace
    back to the one malformed row that caused it.
    """
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
        return False
    reference_map = feature_maps.setdefault(reference_accession, {})
    reference_map[entry['product']] = choose_feature_entry(reference_map.get(entry['product']), entry)
    return True


def feature_map_has_gene_information(feature_map):
    return any(
        normalize_lookup_key(product_name) not in {'polyprotein', 'wholegenome'}
        for product_name in (feature_map or {})
    )


def load_db_gff_feature_maps(conn, alias_lookup):
    table_names = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", conn)['name'].tolist()
    candidate_tables = [name for name in ['gff_features', 'reference_gff_features', 'gff', 'reference_gff'] if name in table_names]
    feature_maps = {}
    unusable_rows = Counter()
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
                if not merge_feature_entry(
                    feature_maps,
                    reference_accession,
                    row.get(product_col),
                    row.get('start'),
                    row.get('end'),
                    feature_type,
                    alias_lookup,
                    raw_product=row.get(product_col),
                    source=f'db_table:{table_name}',
                ):
                    unusable_rows[table_name] += 1
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

                is_allowed = reference_accession in allowed_accessions if allowed_accessions else True
                if not is_allowed:
                    continue
                # NOTE: a guard used to sit here testing
                #   reference_accession not in {ref_candidate, master_candidate, reference_accession}
                # whose set always contained the value under test, so it could
                # never fire.  It is gone rather than repaired: on the shipped
                # HCV build 120 rows in allowed_accessions disagree with both
                # reference columns, and every one of them is a genotype
                # reference's own annotation.  Making the guard live would drop
                # them and change HCV output, so the permissive reading it has
                # always had is the one kept - now deliberately.

                feature_type = str(row.get(feature_type_col, 'gene') or 'gene') if feature_type_col else 'gene'
                if not merge_feature_entry(
                    feature_maps,
                    reference_accession,
                    row.get(product_col),
                    row.get(start_col),
                    row.get(end_col),
                    feature_type,
                    alias_lookup,
                    raw_product=row.get(product_col),
                    source='db_table:features',
                ):
                    unusable_rows['features'] += 1

    if unusable_rows:
        detail = ', '.join(f'{table}={count}' for table, count in sorted(unusable_rows.items()))
        print(
            '[AnnotateMutations][warn] Skipped feature rows with an unusable product name '
            f'or coordinate: {detail}'
        )

    return feature_maps


#: A catalogue position is written as digits.  A trailing '.0' is tolerated
#: because a whole number that made a round trip through a float column is the
#: one respelling that carries no ambiguity; everything else is rejected.
AA_POSITION_RE = re.compile(r'^(?P<digits>\d+)(?:\.0*)?$')


def parse_aa_position(value):
    """A catalogue amino-acid position as a positive integer, or None.

    Every rejection lands in the caller's ``invalid_positions`` tally, which is
    reported before annotation starts.  Letting a bad cell through is worse
    than dropping its row: a position that is merely *wrong* still resolves to
    a codon, and the call it produces is attributed to a residue nobody
    curated.

    Rejects, in the order they bite in practice:
      * text that is not a number at all ('abc', '', NaN);
      * '1e400', which float() reads as infinity - int() then raises
        OverflowError, which is not a ValueError and so escaped the previous
        handler and aborted the entire run before a single sequence was read;
      * scientific notation such as '1e1'.  int(float('1e1')) is 10, but a
        catalogue cell in that shape is far more often a mangled accession or
        identifier than a deliberate position, and reading it as residue 10 is
        unrecoverable once written;
      * a fractional position such as '6.7', which int(float(...)) truncated to
        6, silently relocating a mutation to a different residue;
      * zero and negatives, which used to survive this far and were then
        tallied under 'reference_coordinate_missing_in_alignment' - a catalogue
        defect hiding behind an alignment-coverage counter.
    """
    match = AA_POSITION_RE.match(clean_cell(value))
    if not match:
        return None
    position = int(match.group('digits'))
    return position if position >= 1 else None


# NOTE: resolve_reference_feature_map() and find_reference_alignment_row() used
# to live here.  Neither had a caller anywhere in scripts/ or tests/, and both
# reached canonicalize_product() through their own make_feature_entry chain -
# so once product canonicalisation became virus-aware they were a second,
# permanently-HCV code path that no test could ever have caught drifting.
# Deleted rather than threaded.


def prepare_catalog(catalog, alias_lookup):
    prepared = catalog.copy()
    prepared['_canonical_protein'] = prepared['protein_name'].apply(lambda value: canonicalize_product(value, alias_lookup))
    prepared['_segment_norm'] = prepared['segment'].apply(normalize_segment)
    prepared['_alt_residue_norm'] = prepared['alt_residue'].apply(normalize_residue)

    aa_positions = []
    invalid_positions = 0
    for value in prepared['aa_position']:
        position = parse_aa_position(value)
        aa_positions.append(position)
        if position is None:
            invalid_positions += 1
    prepared['_aa_position_int'] = aa_positions

    return prepared, invalid_positions


def validate_required_catalog_columns(catalog):
    missing = [column for column in REQUIRED_MUTATION_CATALOG_COLUMNS if column not in catalog.columns]
    if missing:
        raise ValueError(f"Catalog missing required column(s): {', '.join(missing)}")


def catalog_columns_for_profile(catalog_column_profile):
    """The extra catalogue columns a profile name selects, or None for all of them.

    Accepts a legacy profile name ('HCV', 'influenza', 'all_columns') or any
    virus name.  An unrecognised name yields None - every column the catalogue
    supplies - rather than raising: the profile argument reaches this from a
    shell variable at both call sites, and losing a curator's columns is a
    worse answer to a typo than keeping all of them.
    """
    if catalog_column_profile in CATALOG_COLUMN_PROFILES:
        return CATALOG_COLUMN_PROFILES[catalog_column_profile]
    return virus_profile(catalog_column_profile)['columns']


def build_catalog_reference_table(catalog, catalog_column_profile):
    extra_columns = catalog_columns_for_profile(catalog_column_profile)
    if extra_columns is None:
        selected_columns = [column for column in catalog.columns if not column.startswith('_')]
    else:
        selected_columns = REQUIRED_MUTATION_CATALOG_COLUMNS + extra_columns

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


#: The alignment columns the annotator reads. `segment` is not optional: it is
#: what keeps a segmented build from resolving one coordinate space for all of
#: its segments at once.
ANNOTATION_ALIGNMENT_COLUMNS = ['sequence_id', 'primary_accession', 'alignment', 'alignment_name', 'segment']


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
    if 'segment' not in seq_aln.columns:
        seq_aln['segment'] = ''
    return seq_aln


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




def resolve_master_coordinate_space(seq_aln, coord_candidates, db_gff_maps, alias_lookup,
                                    allow_genbank_reference_gff):
    """``(accession, coord_map, feature_map)`` for one set of alignment rows.

    Answers ``(None, None, None)`` rather than raising, because with more than
    one segment in play a segment that cannot be resolved must not take the
    other seven down with it.  The caller decides whether an empty result is
    fatal.

    All catalog positions are defined relative to the master reference, and
    every sequence sharing that master's alignment is padded to the same width,
    so the master's column positions apply to all of them.  Per-genotype
    reference OG coordinates must NOT be used: they produce shifted column
    indices for sequences aligned to a divergent genotype reference (D10988 NS3
    cds_start_OG=3345 against NC_004102's 3420 shifts everything downstream by
    ~71 bp relative to the master numbering).
    """
    def coord_map_for(candidate):
        mask = (
            seq_aln['sequence_id'].astype(str).str.strip() == candidate
        ) | (
            seq_aln['primary_accession'].astype(str).str.strip() == candidate
        )
        if not mask.any():
            return None
        return build_alignment_coordinate_map(seq_aln.loc[mask].iloc[0]['alignment'])

    for candidate in coord_candidates:
        cand_coord_map = coord_map_for(candidate)
        if cand_coord_map is None:
            continue
        for token in extract_accession_tokens(candidate):
            if token in db_gff_maps and feature_map_has_gene_information(db_gff_maps[token]):
                return candidate, cand_coord_map, db_gff_maps[token]

    # Fallback: any reference group in this set with a usable feature map.
    for ref_id in seq_aln['_reference_id'].dropna().unique():
        mask = seq_aln['sequence_id'].astype(str).str.strip() == str(ref_id).strip()
        if not mask.any():
            continue
        for token in extract_accession_tokens(str(ref_id)):
            if token in db_gff_maps and feature_map_has_gene_information(db_gff_maps[token]):
                return (
                    str(ref_id),
                    build_alignment_coordinate_map(seq_aln.loc[mask].iloc[0]['alignment']),
                    db_gff_maps[token],
                )

    if not allow_genbank_reference_gff:
        return None, None, None

    for candidate in coord_candidates:
        cand_coord_map = coord_map_for(candidate)
        if cand_coord_map is None:
            continue
        try:
            fetched_map = fetch_reference_gff_map(candidate, alias_lookup)
        except Exception:
            continue
        if fetched_map:
            return candidate, cand_coord_map, fetched_map

    return None, None, None


def build_segment_position_index(catalog, master_coord_map, master_feature_map, diagnostics,
                                 protein_miss_counts, proteins_resolved):
    """``{(protein, aa_pos, alignment_indices): {alt_residue: [rows]}}`` for one segment.

    Every catalog row that does not make it into the index is accounted for:
    an unusable protein name or position, a protein absent from this
    reference's features, a position past the protein's own end, or a
    coordinate the master's alignment does not carry.
    """
    grouped_positions = {}
    mappable_catalog_rows = 0

    for _, row in catalog.iterrows():
        protein_name = row['_canonical_protein']
        aa_pos = row['_aa_position_int']
        if not protein_name or aa_pos is None:
            diagnostics['invalid_catalog_rows'] += 1
            continue

        gene_entry = master_feature_map.get(protein_name)
        if gene_entry is None:
            protein_miss_counts[protein_name] += 1
            continue
        proteins_resolved.add(protein_name)

        # A position has to fall inside the protein it is catalogued against.
        # Only the alignment bounded it before, so a position past the feature
        # end resolved into the NEXT gene's reading frame and the residue found
        # there was reported under this protein's name - correct-looking output
        # attributed to the wrong protein.  Mature-peptide annotations are
        # per-reference and can be truncated, and choose_feature_entry() breaks
        # a priority tie by preferring the SHORTEST span, so this is reachable
        # without any curation error at all.
        feature_end = gene_entry.get('cds_end')
        if feature_end is not None:
            last_nucleotide = gene_entry['cds_start'] + aa_pos * 3 - 1
            if last_nucleotide > feature_end:
                diagnostics['catalog_position_past_feature_end'] += 1
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

    return grouped_positions, mappable_catalog_rows


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
        blank_reference_count = int((seq_aln['_reference_id'] == '').sum())
        if len(master_candidates) == 1:
            seq_aln.loc[seq_aln['_reference_id'] == '', '_reference_id'] = master_candidates[0]
            print(
                f'[AnnotateMutations][warn] Filled {blank_reference_count} blank alignment_name '
                f'values with master accession {master_candidates[0]}'
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
    # Segments.
    #
    # A segmented virus keeps every segment's alignment in the same table, each
    # padded to its OWN width, with its own master reference and its own
    # proteins.  Resolving one coordinate space for the whole run and scanning
    # every row against it produced two failures on real data: where the widths
    # happened to overlap, an HA (segment 4) catalogue entry was read out of an
    # NA (segment 6) sequence's columns and emitted carrying segment '4' - a
    # call asserting a provenance it did not have - and on the shipped 8-segment
    # influenza build it simply aborted, because the one master it picked
    # (segment 6, which annotates only NA) left every other segment's proteins
    # 'not found in master feature map'.
    #
    # So the run is bucketed by segment.  With a single bucket - every
    # non-segmented virus, and every caller that passes no segment column at
    # all - this is exactly the single coordinate space it has always been, and
    # catalogue segment values are not consulted.  Segment values are compared
    # through normalize_segment() because the shipped HCV build really does
    # store both '1' and '1.0' in this column.
    # -----------------------------------------------------------------------
    if 'segment' in seq_aln.columns:
        seq_aln['_segment_norm'] = seq_aln['segment'].apply(normalize_segment)
    else:
        seq_aln['_segment_norm'] = ''
    segment_keys = sorted(set(seq_aln['_segment_norm'].tolist()))
    segment_aware = len(segment_keys) > 1

    coord_candidates = []
    for cand in list(master_candidates) + catalog_reference_hints:
        for token in extract_accession_tokens(cand):
            if token and token not in coord_candidates:
                coord_candidates.append(token)

    signature_scope = build_signature_genotype_scope(catalog)
    wild_type_tables = build_wild_type_tables(catalog)
    sequence_genotypes = build_sequence_genotype_map(meta_data, diagnostics)

    resolved_maps = {}
    protein_miss_counts = Counter()
    proteins_resolved = set()
    mappable_catalog_rows = 0

    if segment_aware:
        catalog_segments = set(catalog['_segment_norm'].tolist())
        stranded = {
            value for value in catalog_segments
            if value and value not in set(segment_keys)
        }
        if stranded:
            stranded_rows = int(catalog['_segment_norm'].isin(stranded).sum())
            diagnostics['catalog_rows_without_matching_segment'] += stranded_rows
            print(
                '[AnnotateMutations][warn] '
                f'{stranded_rows} catalog row(s) name segment(s) '
                f"{', '.join(sorted(stranded))}, which no sequence in this build carries"
            )

    for segment_key in segment_keys:
        if segment_aware:
            segment_aln = seq_aln[seq_aln['_segment_norm'] == segment_key]
            # A catalogue row with no segment of its own applies to every
            # segment: it is the shape a non-segmented catalogue has, and there
            # is nothing better to do with it than try it everywhere.
            segment_catalog = catalog[catalog['_segment_norm'].isin({segment_key, ''})]
        else:
            segment_aln = seq_aln
            segment_catalog = catalog

        label = f'segment {segment_key}' if segment_aware else 'all sequences'
        if segment_catalog.empty:
            diagnostics['segments_without_catalog_rows'] += 1
            print(f'[AnnotateMutations][warn] No catalog rows apply to {label}')
            continue

        coord_ref_acc, coord_map, feature_map = resolve_master_coordinate_space(
            segment_aln, coord_candidates, db_gff_maps, alias_lookup, allow_genbank_reference_gff
        )
        if coord_map is None or feature_map is None:
            diagnostics['segments_without_reference_coordinates'] += 1
            print(
                f'[AnnotateMutations][warn] Could not resolve a reference coordinate space for '
                f'{label}; its sequences are not annotated'
            )
            continue

        if segment_aware:
            print(
                f'[AnnotateMutations] Coordinate reference for segment {segment_key}: '
                f'{coord_ref_acc} ({len(segment_aln)} sequences)'
            )
        else:
            print(
                f'[AnnotateMutations] Coordinate reference: {coord_ref_acc} '
                f"(all catalog positions resolved using this sequence's alignment column space)"
            )

        grouped_positions, segment_mappable = build_segment_position_index(
            segment_catalog, coord_map, feature_map, diagnostics,
            protein_miss_counts, proteins_resolved,
        )
        mappable_catalog_rows += segment_mappable
        resolved_maps[coord_ref_acc] = {
            'resolved_accession': coord_ref_acc,
            'feature_map': feature_map,
            'attempted': coord_candidates,
            'source': 'master_coord_ref',
            'segment': segment_key,
            'aligned_reference_accession': coord_ref_acc,
            'aligned_reference_map': coord_map,
        }
        if not grouped_positions:
            continue

        for _, seq_row in segment_aln.iterrows():
            alignment = seq_row['alignment']
            primary_accession = seq_row['sequence_id']
            genotype, subtype_code = sequence_genotypes.get(str(primary_accession).strip(), ('', ''))
            if not genotype:
                diagnostics['sequences_without_genotype'] += 1
            covered_span = alignment_covered_span(alignment)
            if covered_span is None:
                # No non-gap column anywhere: this row reports no bases at all,
                # so it is evidence for nothing.  It has to be stopped here
                # rather than in residue_from_aligned_codon(), which reads a
                # missing covered span as 'the caller gave me a bare codon' and
                # answers with the deletion residue - which would fire every
                # catalogued deletion against the one sequence carrying no data
                # whatsoever.
                diagnostics['sequences_without_sequenced_bases'] += 1
                continue
            for (protein_name, aa_pos, alignment_indices), alt_lookup in grouped_positions.items():
                codon = extract_aligned_codon(alignment, alignment_indices)
                if codon is None:
                    diagnostics['codon_out_of_bounds'] += 1
                    continue
                aa = residue_from_aligned_codon(codon, covered_span, alignment_indices)
                if aa == UNKNOWN_RESIDUE:
                    # 'X' is the absence of a reading, not a residue.  Letting it
                    # match meant any catalogue row spelled 'X' - a shape the
                    # normalizer's [A-Z*] token grammar admits - fired on every
                    # unsequenced or ambiguous codon, flagging the worst-covered
                    # records the hardest.
                    diagnostics['codon_residue_unresolved'] += 1
                    continue
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
                    if segment_aware and not row['_segment_norm']:
                        # This build has several segments and the catalogue row
                        # named none, so it was tried against every one of them.
                        # The call is not wrong, but which segment it came from
                        # rests on an assumption rather than on the catalogue.
                        diagnostics['emitted_from_segment_unstated_catalog_row'] += 1
                    diagnostics[f'emitted_{scope_tier}'] += 1
                    diagnostics[f'emitted_{residue_status}'] += 1
                    if aa == DELETION_RESIDUE:
                        diagnostics['emitted_deletions'] += 1

    if not resolved_maps:
        message = (
            'Could not resolve a master reference alignment + feature map. '
            f'Candidates tried: {coord_candidates}.'
        )
        if not allow_genbank_reference_gff:
            message += ' Re-run with --allow_genbank_reference_gff to fetch from NCBI.'
        raise AnnotationMappingError(message)

    proteins_missing_from_reference = Counter({
        protein: count for protein, count in protein_miss_counts.items()
        if protein not in proteins_resolved
    })

    if mappable_catalog_rows == 0:
        # Only a MAPPING failure is fatal - not one protein of the catalogue
        # could be found in any reference's features, which means the build is
        # misconfigured and every call would be missing.  If proteins did
        # resolve and their rows were simply all rejected (every position past
        # its protein's end, say), that is a catalogue problem already counted
        # in the diagnostics, and an empty result is the honest answer to it
        # rather than a crash.
        missing_summary = ', '.join(
            protein for protein, _ in proteins_missing_from_reference.most_common(10)
        )
        if not proteins_resolved:
            raise AnnotationMappingError(
                'No catalog coordinates could be mapped onto master reference CDS annotations. '
                f'Proteins not found in master feature map: {missing_summary}'
            )
        print(
            '[AnnotateMutations][warn] Every catalog row was rejected before annotation; '
            'see the mapping summary for the reason'
        )

    if proteins_missing_from_reference:
        # These counts have to reach the diagnostics Counter, not just stdout:
        # the 'Mapping summary' line built from it is what an operator and any
        # downstream check actually parse, and only a TOTAL mapping failure
        # raises.  Losing nine catalogue proteins out of ten was not an error
        # and left no trace in the summary at all.
        diagnostics['catalog_rows_without_reference_feature'] = sum(
            proteins_missing_from_reference.values()
        )
        diagnostics['proteins_without_reference_feature'] = len(proteins_missing_from_reference)
        preview = ', '.join(
            protein for protein, _ in proteins_missing_from_reference.most_common(10)
        )
        print(f'[AnnotateMutations][warn] Missing protein annotations in master feature map: {preview}')

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
    parser.add_argument(
        "--virus", default="",
        help="Virus this build is for (e.g. HCV, influenza). Selects the protein-name "
             "vocabulary used to canonicalise catalog protein names against reference "
             "features, and the default set of catalog columns written to the database. "
             f"Known: {', '.join(sorted(VIRUS_PROFILES))}. Anything else is annotated with "
             "the generic profile, which renames no proteins and keeps every catalog column.",
    )
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
        default=None,
        help="Which mutation catalog columns are written to the database. Defaults to the "
             "column set for --virus, and to every column the catalog supplies when the "
             "virus is unknown or unspecified. Accepts a virus name or the legacy values "
             "'HCV', 'influenza' and 'all_columns'. This used to be required with a fixed "
             "choices list, which made every virus outside that list fail at argument "
             "parsing.",
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

    virus_name = resolve_virus_name(args.virus, args.catalog_column_profile)
    profile = VIRUS_PROFILES[virus_name]
    catalog_column_profile = args.catalog_column_profile or args.virus or GENERIC_VIRUS
    requested_virus = clean_cell(args.virus) or clean_cell(args.catalog_column_profile)
    if requested_virus and virus_name == GENERIC_VIRUS and resolve_virus_name(requested_virus) == GENERIC_VIRUS:
        if clean_cell(args.catalog_column_profile) != 'all_columns':
            print(
                f"[AnnotateMutations][warn] No virus profile named {requested_virus!r}; using the "
                f"generic profile (no protein-name inference, every catalog column kept)"
            )
    print(
        f'[AnnotateMutations] Virus profile: {virus_name}; catalog column profile: '
        f'{catalog_column_profile}'
    )

    print("Loading mutation catalog...")
    # keep_default_na=False: pandas' default NA sentinels include the literal
    # 'NA', which is the standard name of influenza's neuraminidase - the
    # protein oseltamivir resistance is catalogued against.  Without it every
    # NA row loses its protein_name to NaN, canonicalize_product() renders that
    # back as the string 'nan', and the entire segment silently fails to
    # resolve.  'NULL', 'None', 'N/A' and 'nan' are sentinels too, and a
    # curated catalogue is entitled to use any of them as a real value.
    catalog = pd.read_csv(args.mutation_catalog, sep='\t', dtype=str, keep_default_na=False)

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

        gene_alias_lookup = load_gene_alias_lookup(conn, profile=profile)
        catalog, invalid_positions = prepare_catalog(catalog, gene_alias_lookup)
        if invalid_positions:
            print(f'[AnnotateMutations][warn] Skipping {invalid_positions} catalog rows with invalid aa_position values')

        db_gff_maps = load_db_gff_feature_maps(conn, gene_alias_lookup)

        print('Extracting mutations...')
        call_evidence = []
        mutations_found, diagnostics, resolved_maps = annotate_from_reference_coordinates(
            catalog,
            # `segment` travels with the alignment: dropping it here was what
            # made a segmented build resolve one coordinate space for every
            # segment at once.
            seq_aln[ANNOTATION_ALIGNMENT_COLUMNS].copy(),
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

        write_mutation_tables(conn, catalog, mutations_found, catalog_column_profile,
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
