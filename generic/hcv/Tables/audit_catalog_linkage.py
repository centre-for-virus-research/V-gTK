#!/usr/bin/env python3
"""Audit how the HCV mutation catalogue links back to the PHDR source tables.

Every value in the catalogue is re-derived from the phdr_*.csv tables and
compared. Works on both layouts:

  * wide  - generalized_mutation_catalog_with_extra_info.tsv
            (pubmed_id / DOI / clinical_trials as semicolon lists)
  * long  - generalized_mutation_catalog_evidence_linked.tsv
            (one evidence_id + data_source per row, genotype instead of
            alignment_name; built by scripts/BuildHcvEvidenceCatalog.py)

Checks are grouped:

  source    PHDR tables agree with each other (foreign keys, evidence flags)
  keys      catalogue ids, alignment codes and drugs resolve to PHDR
  values    resistance, drug and wild-type columns match their source
  genotype  the row's genotype survives into the database
  evidence  publications and trials attach to the right entry, and to each other

Each check reports PASS, WARN (a real but tolerated data shape) or FAIL
(the catalogue says something the source does not), with examples.

Usage:
    python audit_catalog_linkage.py                         # wide file
    python audit_catalog_linkage.py --catalog <long.tsv>    # long file
    python audit_catalog_linkage.py --json audit.json       # machine-readable
Exit code is 1 if any check FAILs.
"""

import argparse
import collections
import csv
import json
import math
import re
import sys
from pathlib import Path

csv.field_size_limit(min(sys.maxsize, 2 ** 31 - 1))

HERE = Path(__file__).resolve().parent
SEP = ';'
MAX_EXAMPLES = 5

#: The annotator's own code decides what reaches the database's mutation_catalog
#: table, so the genotype check below runs it rather than restating its column
#: list here.
SCRIPTS = HERE.parents[2] / 'scripts'

ARD_VALUE_COLUMNS = ['resistance_category', 'display_resistance_category', 'numeric_resistance_category',
                     'any_in_vitro_evidence', 'in_vitro_max_ec50_midpoint', 'any_in_vivo_evidence',
                     'in_vivo_baseline', 'in_vivo_treatment_emergent']
NUMERIC_COLUMNS = ARD_VALUE_COLUMNS[2:] + ['combination_size', 'component_order']

EVIDENCE_COLUMNS = ['evidence_id', 'data_source', 'evidence_url', 'evidence_label',
                    'evidence_type', 'regimens', 'linked_evidence_ids', 'finding_ids']


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def read_rows(path, delimiter=','):
    with open(path, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames), [{k: (v or '') for k, v in r.items()} for r in reader]


def split(value):
    return [x for x in (value or '').split(SEP) if x]


def genotype_of_alignment(name):
    name = (name or '').strip()
    return name[3:] if name.upper().startswith('AL_') else name


def same_number(a, b):
    if a == b:
        return True
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-9)
    except ValueError:
        return False


class Audit:
    def __init__(self):
        self.results = []

    def check(self, group, name, failures, *, level='FAIL', total=None, explain='', examples=None):
        """Record a check. ``failures`` is a count or a list of offending items."""
        items = failures if isinstance(failures, list) else None
        count = len(items) if items is not None else int(failures)
        status = 'PASS' if count == 0 else level
        self.results.append({
            'group': group, 'check': name, 'status': status, 'count': count, 'total': total,
            'explain': explain,
            'examples': (examples if examples is not None else (items or []))[:MAX_EXAMPLES],
        })

    def info(self, group, name, value, explain=''):
        self.results.append({'group': group, 'check': name, 'status': 'INFO', 'count': value,
                             'total': None, 'explain': explain, 'examples': []})

    def print(self, catalog, layout):
        print(f'Catalogue : {catalog}  ({layout} layout)\n')
        group = None
        for r in self.results:
            if r['group'] != group:
                group = r['group']
                print(f'== {group}')
            total = f"/{r['total']}" if r['total'] is not None else ''
            print(f"  [{r['status']:<4}] {r['check']}: {r['count']}{total}")
            if r['status'] != 'PASS' and r['explain']:
                print(f"         {r['explain']}")
            if r['status'] in ('FAIL', 'WARN'):
                for ex in r['examples']:
                    print(f'           e.g. {ex}')
        tally = collections.Counter(r['status'] for r in self.results)
        print('\nSummary: ' + ', '.join(f'{k} {tally[k]}' for k in ('PASS', 'WARN', 'FAIL', 'INFO') if tally[k]))
        return tally['FAIL']


# ---------------------------------------------------------------------------
# source tables
# ---------------------------------------------------------------------------

class Source:
    def __init__(self, tables):
        t = Path(tables)
        _, self.ras = read_rows(t / 'phdr_ras.csv')
        _, self.aln_ras = read_rows(t / 'phdr_alignment_ras.csv')
        _, self.ard = read_rows(t / 'phdr_alignment_ras_drug.csv')
        _, self.drugs = read_rows(t / 'phdr_drug.csv')
        _, self.findings = read_rows(t / 'phdr_resistance_finding.csv')
        _, self.pubs = read_rows(t / 'phdr_publication.csv')
        _, self.result_trial = read_rows(t / 'phdr_result_trial.csv')
        _, self.result_regimen = read_rows(t / 'phdr_result_regimen.csv')
        _, self.trials = read_rows(t / 'phdr_clinical_trial.csv')
        _, self.typical = read_rows(t / 'phdr_alignment_typical_aa.csv')

        self.ard_by_id = {r['id']: r for r in self.ard}
        self.aln_ras_by_id = {r['id']: r for r in self.aln_ras}
        self.drug_by_id = {r['id']: r for r in self.drugs}
        self.pub_by_id = {r['id']: r for r in self.pubs}
        self.registry_of = {r['id']: (r['nct_id'].strip() or r['id'].strip()) for r in self.trials}
        self.trials_by_result = collections.defaultdict(set)
        for r in self.result_trial:
            self.trials_by_result[r['phdr_in_vivo_result_id']].add(r['phdr_clinical_trial_id'])
        self.regimens_by_result = collections.defaultdict(set)
        for r in self.result_regimen:
            self.regimens_by_result[r['phdr_in_vivo_result_id']].add(r['phdr_regimen_id'])

        # Expected evidence per catalogue id, at finding resolution.
        self.pubs_for = collections.defaultdict(set)
        self.trials_for = collections.defaultdict(set)
        self.pub_types = collections.defaultdict(set)        # (id, pub) -> {in_vitro, in_vivo}
        self.pub_trials = collections.defaultdict(set)       # (id, pub) -> {registry}
        self.trial_pubs = collections.defaultdict(set)       # (id, registry) -> {pub}
        self.link_findings = collections.defaultdict(set)    # (id, evidence) -> {finding}
        self.link_regimens = collections.defaultdict(set)    # (id, evidence) -> {regimen}
        for f in self.findings:
            key, pub = f['phdr_alignment_ras_drug_id'], f['phdr_publication_id']
            self.pubs_for[key].add(pub)
            self.link_findings[(key, pub)].add(f['id'])
            if f['phdr_in_vitro_result_id']:
                self.pub_types[(key, pub)].add('in_vitro')
            if f['phdr_in_vivo_result_id']:
                self.pub_types[(key, pub)].add('in_vivo')
                regimens = self.regimens_by_result.get(f['phdr_in_vivo_result_id'], set())
                self.link_regimens[(key, pub)] |= regimens
                for curator in self.trials_by_result.get(f['phdr_in_vivo_result_id'], ()):
                    reg = self.registry_of.get(curator, curator)
                    self.trials_for[key].add(reg)
                    self.pub_trials[(key, pub)].add(reg)
                    self.trial_pubs[(key, reg)].add(pub)
                    self.link_findings[(key, reg)].add(f['id'])
                    self.link_regimens[(key, reg)] |= regimens

        # Dominant residue per (genotype code, feature, codon).
        self.dominant = {}
        self.typical_residues = collections.defaultdict(set)
        for r in self.typical:
            k = (genotype_of_alignment(r['alignment_name']), r['feature_name'], r['codon_label'])
            pct = float(r['pct_members'] or 0)
            self.typical_residues[k].add(r['aa_residue'])
            if k not in self.dominant or pct > self.dominant[k][1]:
                self.dominant[k] = (r['aa_residue'], pct)


def audit_source(a, s):
    g = 'source'
    ras_ids = {r['id'] for r in s.ras}
    a.check(g, 'alignment_ras -> ras', [r['id'] for r in s.aln_ras if r['phdr_ras_id'] not in ras_ids],
            total=len(s.aln_ras))
    a.check(g, 'alignment_ras_drug -> alignment_ras / drug',
            [r['id'] for r in s.ard if r['phdr_alignment_ras_id'] not in s.aln_ras_by_id
             or r['phdr_drug_id'] not in s.drug_by_id], total=len(s.ard))
    a.check(g, 'resistance_finding -> alignment_ras_drug / publication',
            [f['id'] for f in s.findings if f['phdr_alignment_ras_drug_id'] not in s.ard_by_id
             or f['phdr_publication_id'] not in s.pub_by_id], total=len(s.findings))
    a.check(g, 'finding has neither in-vitro nor in-vivo result',
            [f['id'] for f in s.findings if not (f['phdr_in_vitro_result_id'] or f['phdr_in_vivo_result_id'])],
            total=len(s.findings))
    in_vivo_ids = {f['phdr_in_vivo_result_id'] for f in s.findings if f['phdr_in_vivo_result_id']}
    trial_ids = {t['id'] for t in s.trials}
    a.check(g, 'result_trial -> in-vivo finding / clinical_trial',
            [r['id'] for r in s.result_trial if r['phdr_in_vivo_result_id'] not in in_vivo_ids
             or r['phdr_clinical_trial_id'] not in trial_ids], total=len(s.result_trial))

    vitro = {f['phdr_alignment_ras_drug_id'] for f in s.findings if f['phdr_in_vitro_result_id']}
    vivo = {f['phdr_alignment_ras_drug_id'] for f in s.findings if f['phdr_in_vivo_result_id']}
    a.check(g, 'any_in_vitro_evidence agrees with findings',
            [r['id'] for r in s.ard if (r['any_in_vitro_evidence'] == '1') != (r['id'] in vitro)], total=len(s.ard))
    a.check(g, 'any_in_vivo_evidence agrees with findings',
            [r['id'] for r in s.ard if (r['any_in_vivo_evidence'] == '1') != (r['id'] in vivo)], total=len(s.ard))

    nct_names = collections.defaultdict(list)
    for t in s.trials:
        if t['nct_id'].strip():
            nct_names[t['nct_id'].strip()].append(t['id'])
    a.check(g, 'NCT numbers shared by several curator trial ids',
            [f'{k}: {v}' for k, v in sorted(nct_names.items()) if len(v) > 1], level='WARN',
            explain='Trial arms / sponsor codes. Key on the registry id and keep the names, or the join double-counts.')
    a.check(g, 'trials with no NCT number', [t['id'] for t in s.trials if not t['nct_id'].strip()], level='WARN',
            explain='Non-ClinicalTrials.gov registrations (e.g. Japanese UMIN-CTR). Still real trials.')
    a.info(g, 'in-vivo results with no trial attached',
           len(in_vivo_ids - set(s.trials_by_result)),
           'In-vivo evidence whose only source is the publication (e.g. cohort studies).')
    kinds = collections.defaultdict(set)
    for (key, pub), types in s.pub_types.items():
        kinds[pub] |= types
    a.info(g, 'publications supplying both in-vitro and in-vivo findings',
           sum(1 for v in kinds.values() if len(v) == 2),
           'A publication id alone does not say which result columns it supports.')


# ---------------------------------------------------------------------------
# catalogue
# ---------------------------------------------------------------------------

def audit_keys(a, s, entries):
    g = 'keys'
    keyed = [e for e in entries if e['id']]
    ids = {e['id'] for e in keyed}
    a.check(g, 'PHDR alignment_ras_drug entries missing from catalogue', sorted(set(s.ard_by_id) - ids),
            total=len(s.ard_by_id))
    a.check(g, 'catalogue ids not in PHDR', sorted(ids - set(s.ard_by_id)), total=len(ids))
    a.check(g, 'id != source_phdr_ras_id : alignment (or AL_genotype) : drug',
            [e['id'] for e in keyed if e['id'] != f"{e['source_phdr_ras_id']}:{e['_alignment']}:{e['drug']}"],
            total=len(keyed))
    a.check(g, 'phdr_alignment_ras_id / phdr_drug_id disagree with id',
            [e['id'] for e in keyed if e['phdr_alignment_ras_id'] != s.ard_by_id.get(e['id'], {}).get('phdr_alignment_ras_id')
             or e['phdr_drug_id'] != e['drug']], total=len(keyed))
    a.check(g, 'display_structure differs from phdr_alignment_ras',
            [e['id'] for e in keyed
             if s.aln_ras_by_id.get(e['phdr_alignment_ras_id'], {}).get('display_structure') != e['display_structure']],
            total=len(keyed))
    a.check(g, 'source_phdr_ras_id != signature_id',
            [e['id'] for e in keyed if e['source_phdr_ras_id'] != e['signature_id']], total=len(keyed))
    unkeyed = [e for e in entries if not e['id']]
    a.check(g, 'rows with no PHDR id that still carry values',
            [e['mutation_id'] for e in unkeyed
             if any(e.get(c) for c in ('_alignment', 'drug', 'resistance_category'))
             or e['_pubs'] or e['_trials']], total=len(unkeyed),
            explain='Decomposed components should be pure anchors.')
    a.info(g, 'rows with no PHDR id (decomposed combination components)', len(unkeyed))


def audit_values(a, s, entries):
    g = 'values'
    keyed = [e for e in entries if e['id'] in s.ard_by_id]
    wrong, drift = [], []
    for e in keyed:
        src = s.ard_by_id[e['id']]
        for col in ARD_VALUE_COLUMNS:
            if not same_number(e[col], src[col]):
                wrong.append(f"{e['id']} {col}: {e[col]!r} vs {src[col]!r}")
            elif re.fullmatch(r'-?\d+\.0', e[col]) and e[col][:-2] == src[col]:
                drift.append(f"{e['id']} {col}: {e[col]!r} vs source {src[col]!r}")
    a.check(g, 'resistance / evidence values differ from alignment_ras_drug', wrong, total=len(keyed) * len(ARD_VALUE_COLUMNS))
    a.check(g, 'same value, different text (e.g. 1 -> 1.0)', drift, level='WARN',
            total=len(keyed) * len(ARD_VALUE_COLUMNS),
            explain='A pandas float round trip. Breaks string joins and exact de-duplication.')
    a.check(g, 'combination_size / component_order written as floats',
            [f"{e['mutation_id']} {c}={e[c]!r}" for e in entries for c in ('combination_size', 'component_order')
             if re.fullmatch(r'\d+\.0', e.get(c, ''))], level='WARN', total=len(entries) * 2)

    bad_drug = []
    for e in entries:
        if not e['drug']:
            continue
        d = s.drug_by_id.get(e['drug'])
        if d is None or e['drug_producer'].strip() != d['producer'].strip() \
                or e['drug_category'] != d['drug_category'].strip():
            bad_drug.append(f"{e['id']}: {e['drug_producer']!r}/{e['drug_category']!r}")
    a.check(g, 'drug_producer / drug_category differ from phdr_drug', bad_drug, total=sum(1 for e in entries if e['drug']))
    a.check(g, 'whitespace in phdr_drug values',
            [f"{d['id']}: {d['producer']!r}" for d in s.drugs if d['producer'] != d['producer'].strip()], level='WARN')

    wt_bad, blank_code = [], []
    own_missing = []
    for e in entries:
        key_feature, key_pos = e['protein_name'], e['aa_position']
        codes = set()
        for entry in split(e['wild_type_residues']):
            parts = entry.split(':')
            code, residue = parts[0], parts[1] if len(parts) > 1 else ''
            codes.add(code)
            if not code:
                blank_code.append(f"{e['mutation_id']}: {entry!r}")
                continue
            dom = s.dominant.get((code, key_feature, key_pos))
            if dom is None or dom[0] != residue or (len(parts) > 2 and abs(float(parts[2]) - dom[1]) > 0.006):
                wt_bad.append(f"{e['mutation_id']} {entry} (source dominant {dom})")
        own = genotype_of_alignment(e['_alignment'])
        if own and (own, key_feature, key_pos) in s.dominant and own not in codes:
            own_missing.append(f"{e['id']}: no {own} entry in wild_type_residues")
    a.check(g, 'wild_type_residues entry is not the dominant residue in typical_aa', wt_bad, total=len(entries))
    a.check(g, 'wild_type_residues entry with a blank genotype code', blank_code, total=len(entries),
            explain='BuildCatalogGenotypeColumns maps AL_MASTER and the 42 single-sequence AL_*_unassigned_* '
                    "alignments to '' and its `code is None` guard never fires, so the first 100% singleton wins. "
                    'Never looked up at runtime (lookups need a genotype), but wrong in the database.')
    a.check(g, "row's own genotype missing from wild_type_residues", own_missing, total=len(entries))


def audit_genotype(a, s, entries, rows, layout):
    g = 'genotype'
    scope = collections.defaultdict(set)
    for e in entries:
        if e['_alignment']:
            scope[e['signature_id']].add(genotype_of_alignment(e['_alignment']))
    a.check(g, 'relevant_genotypes != union of alignment codes for the signature',
            [f"{e['signature_id']}: {e['relevant_genotypes']!r}" for e in entries
             if set(split(e['relevant_genotypes'])) != scope.get(e['signature_id'], set())], total=len(entries))

    if layout == 'long':
        a.check(g, "genotype disagrees with the alignment inside the PHDR id",
                [e['id'] for e in entries if e['id'] and e['genotype'] != genotype_of_alignment(e['id'].rsplit(':', 2)[1])],
                total=sum(1 for e in entries if e['id']))
        a.check(g, 'unkeyed row with a genotype',
                [e['mutation_id'] for e in entries if not e['id'] and e['genotype']])

    genotypes = collections.defaultdict(set)
    for e in entries:
        if e['_alignment']:
            genotypes[(e['signature_id'], e['drug'], e['mutation_id'])].add(genotype_of_alignment(e['_alignment']))
    multi = {k for k, v in genotypes.items() if len(v) > 1}
    a.info(g, '(signature, drug, mutation) keys scored in more than one genotype', len(multi))

    # What the database holds: run the annotator's own catalogue writer and ask
    # whether each resistance category can still be pinned to one genotype.
    try:
        sys.path.insert(0, str(SCRIPTS))
        import pandas as pd
        import AnnotateMutations as AM
    except ImportError as exc:  # pragma: no cover - outside the vgtk env
        a.info(g, 'database check skipped', str(exc))
        return
    frame = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith('_')} for r in rows])
    table = AM.build_catalog_reference_table(frame, 'HCV')
    key = ['signature_id', 'drug', 'mutation_id'] + (['genotype'] if 'genotype' in table.columns else [])
    keyed_table = table[table['drug'] != '']
    categories = keyed_table.groupby(key)['resistance_category'].nunique()
    ambiguous = categories[categories > 1]
    a.check(g, 'database rows whose resistance category cannot be pinned to a genotype',
            [' / '.join(k) for k in ambiguous.index], total=len(categories),
            explain='The HCV column profile drops alignment_name and relevant_genotypes is signature-wide, so '
                    'the 1a and 1b rows of one finding are indistinguishable. The profile must keep a row-level '
                    'genotype column.')
    a.info(g, 'mutation_catalog rows as written by AnnotateMutations', len(table))


def audit_evidence_wide(a, s, entries):
    g = 'evidence'
    keyed = [e for e in entries if e['id']]
    a.check(g, 'pubmed_id set != publications in resistance_finding',
            [f"{e['id']}: {sorted(e['_pubs'])} vs {sorted(s.pubs_for.get(e['id'], ()))}" for e in keyed
             if set(e['_pubs']) != s.pubs_for.get(e['id'], set())], total=len(keyed))
    doi_bad = []
    for e in keyed:
        dois = split(e['DOI'])
        want = [s.pub_by_id.get(p, {}).get('url', '') for p in e['_pubs']]
        if dois != want:
            doi_bad.append(f"{e['id']}: {dois} vs {want}")
    a.check(g, 'DOI list not positionally paired with pubmed_id', doi_bad, total=len(keyed))
    a.check(g, 'clinical_trials set != trials via in-vivo findings',
            [f"{e['id']}: missing {sorted(s.trials_for.get(e['id'], set()) - set(e['_trials']))} "
             f"extra {sorted(set(e['_trials']) - s.trials_for.get(e['id'], set()))}"
             for e in keyed if set(e['_trials']) != s.trials_for.get(e['id'], set())], total=len(keyed))

    unpaired = []
    mixed = []
    for e in keyed:
        pubs = s.pubs_for.get(e['id'], set())
        trials = s.trials_for.get(e['id'], set())
        if len(pubs) > 1 and trials:
            # Can "which paper reported which trial" be recovered from the row?
            # Only if every paper reported every trial.
            if any(s.pub_trials.get((e['id'], p), set()) != trials for p in pubs):
                unpaired.append(e['id'])
        types = {frozenset(s.pub_types.get((e['id'], p), ())) for p in pubs}
        if len({t for ts in types for t in ts}) == 2:
            mixed.append(e['id'])
    a.check(g, 'rows where paper <-> trial pairing cannot be recovered', unpaired, level='FAIL', total=len(keyed),
            explain='pubmed_id and clinical_trials are independent lists; the finding that ties a trial '
                    'to the paper that reported it is not in the row.')
    a.check(g, 'rows mixing in-vitro and in-vivo sources in one pubmed_id list', mixed, level='FAIL', total=len(keyed),
            explain='Cannot tell which paper supports in_vitro_max_ec50_midpoint vs in_vivo_* columns.')
    a.check(g, 'pubmed_id holding non-PubMed references',
            sorted({p for e in keyed for p in e['_pubs'] if not p.isdigit()}), level='WARN',
            explain='Conference abstracts in a column named pubmed_id.')


def row_key(r):
    """One catalogue row: a PHDR id spans one row per combination component."""
    return (r['id'], r['signature_id'], r['mutation_id'], r['component_order'])


def audit_evidence_long(a, s, rows):
    g = 'evidence'
    by_key = collections.defaultdict(list)
    for r in rows:
        by_key[row_key(r)].append(r)

    want = set()
    for k in by_key:
        key = k[0]
        want |= {k + ('pub', p) for p in s.pubs_for.get(key, ())}
        want |= {k + ('trial', t) for t in s.trials_for.get(key, ())}
    got = collections.Counter()
    for r in rows:
        if r['id'] and r['evidence_id']:
            kind = 'trial' if r['data_source'] == 'clinical_trial' else 'pub'
            got[row_key(r) + (kind, r['evidence_id'])] += 1
    a.check(g, 'evidence links in PHDR missing from file', sorted(map(str, want - set(got))), total=len(want))
    a.check(g, 'evidence links in file not in PHDR', sorted(map(str, set(got) - want)), total=len(got))
    a.check(g, 'evidence link written more than once', [str(k) for k, v in got.items() if v > 1], total=len(got))
    a.check(g, 'catalogue entries with PHDR evidence but a blank evidence row',
            [r['id'] for r in rows if r['id'] and not r['evidence_id'] and s.pubs_for.get(r['id'])], total=len(rows))

    src_bad, type_bad, link_bad, find_bad, url_bad, regimen_bad = [], [], [], [], [], []
    for r in rows:
        if not r['evidence_id']:
            if any(r[c] for c in EVIDENCE_COLUMNS):
                src_bad.append(f"{r['mutation_id']}: evidence columns set without evidence_id")
            continue
        key, ev, source = r['id'], r['evidence_id'], r['data_source']
        if source == 'clinical_trial':
            ok = ev in s.trials_for.get(key, set())
            want_type, want_link = 'in_vivo', s.trial_pubs.get((key, ev), set())
            want_url = f'https://clinicaltrials.gov/study/{ev}' if ev.startswith('NCT') else ''
        else:
            ok = source == ('pubmed' if ev.isdigit() else 'conference_abstract')
            want_type = SEP.join(sorted(s.pub_types.get((key, ev), ())))
            want_link = s.pub_trials.get((key, ev), set())
            want_url = s.pub_by_id.get(ev, {}).get('url', '')
        if not ok:
            src_bad.append(f'{key} {ev}: data_source={source!r}')
        if r['evidence_type'] != want_type:
            type_bad.append(f"{key} {ev}: {r['evidence_type']!r} vs {want_type!r}")
        if set(split(r['linked_evidence_ids'])) != want_link:
            link_bad.append(f"{key} {ev}: {r['linked_evidence_ids']!r} vs {sorted(want_link)}")
        if set(split(r['finding_ids'])) != s.link_findings.get((key, ev), set()):
            find_bad.append(f'{key} {ev}')
        if r['evidence_url'] != want_url:
            url_bad.append(f"{key} {ev}: {r['evidence_url']!r} vs {want_url!r}")
        if set(split(r.get('regimens', ''))) != s.link_regimens.get((key, ev), set()):
            regimen_bad.append(f"{key} {ev}: {r.get('regimens', '')!r} vs "
                               f"{sorted(s.link_regimens.get((key, ev), set()))}")
    n = sum(1 for r in rows if r['evidence_id'])
    a.check(g, 'data_source does not match the evidence id', src_bad, total=n)
    a.check(g, 'evidence_type differs from the findings', type_bad, total=n)
    a.check(g, 'linked_evidence_ids differs from finding-level pairing', link_bad, total=n)
    a.check(g, 'finding_ids differ from resistance_finding', find_bad, total=n)
    a.check(g, 'evidence_url differs from publication / registry', url_bad, total=n)
    a.check(g, 'regimens differ from result_regimen', regimen_bad, total=n)

    # The fan-out must not change the entry's own values.
    fixed = [c for c in rows[0] if c not in EVIDENCE_COLUMNS]
    drift = []
    for key, group in by_key.items():
        if key[0] and len({tuple(r[c] for c in fixed) for r in group}) > 1:
            drift.append(key)
    a.check(g, 'entry values differ between its evidence rows', [k[0] + ' / ' + k[2] for k in drift], total=len(by_key))
    a.info(g, 'rows by data_source',
           dict(collections.Counter(r['data_source'] or '(none)' for r in rows)))


def collapse_entries(rows, layout):
    """One dict per catalogue entry, with parsed evidence lists attached."""
    if layout == 'wide':
        for r in rows:
            r['_pubs'] = split(r['pubmed_id'])
            r['_trials'] = split(r['clinical_trials'])
            r['_alignment'] = r['alignment_name']
        return rows
    entries, seen = [], {}
    for r in rows:
        # Unkeyed component rows are distinct catalogue rows even when id is blank.
        k = row_key(r)
        if k not in seen:
            e = dict(r)
            e['_pubs'], e['_trials'] = [], []
            e['_alignment'] = f"AL_{r['genotype']}" if r['genotype'] else ''
            seen[k] = e
            entries.append(e)
        if r['evidence_id']:
            (seen[k]['_trials'] if r['data_source'] == 'clinical_trial' else seen[k]['_pubs']).append(r['evidence_id'])
    return entries


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--catalog', default=str(HERE / 'generalized_mutation_catalog_with_extra_info.tsv'))
    parser.add_argument('--tables', default=str(HERE))
    parser.add_argument('--json', default=None, help='also write results as JSON')
    args = parser.parse_args(argv)

    s = Source(args.tables)
    _, rows = read_rows(args.catalog, delimiter='\t')
    layout = 'long' if 'evidence_id' in rows[0] else 'wide'
    entries = collapse_entries(rows, layout)

    a = Audit()
    audit_source(a, s)
    audit_keys(a, s, entries)
    audit_values(a, s, entries)
    audit_genotype(a, s, entries, rows, layout)
    if layout == 'wide':
        audit_evidence_wide(a, s, entries)
    else:
        audit_evidence_long(a, s, rows)

    failures = a.print(args.catalog, layout)
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as handle:
            json.dump({'catalog': str(args.catalog), 'layout': layout, 'rows': len(rows),
                       'entries': len(entries), 'results': a.results}, handle, indent=1, default=str)
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
