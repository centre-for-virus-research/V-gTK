#!/usr/bin/env python3
"""Expand the HCV catalogue to one evidence reference per row.

THE BUILD
---------
    NormaliseHcvMutationCatalog.py  PHDR tables -> generalized_mutation_catalog.tsv
                                    one row per (component, genotype, drug)
    BuildHcvEvidenceCatalog.py      -> generalized_mutation_catalog_evidence_linked.tsv
                                    (what the pipeline reads)

This script adds everything the normaliser does not: drug producer and
category (phdr_drug.csv), ``relevant_genotypes`` and ``wild_type_residues``
(computed with BuildCatalogGenotypeColumns.py), and the evidence block. Both
steps are deterministic, so the shipped file is always a rebuild of the PHDR
export and nothing is added by hand.

WHY ONE EVIDENCE REFERENCE PER ROW
----------------------------------
The earlier hand-assembled catalogue
(generalized_mutation_catalog_with_extra_info.tsv) kept its evidence as three
independent semicolon lists - ``pubmed_id``, ``DOI`` and ``clinical_trials`` -
on one row per (signature, genotype, drug). That loses two links PHDR records:

* which paper reported which trial. Every PHDR finding cites exactly one
  publication; a trial only appears behind an in-vivo result that publication
  reports. They are not alternative sources.
* which paper supports the in-vitro columns and which the in-vivo ones.

It also keeps the row's genotype only in ``alignment_name`` (``AL_1a``), a PHDR
join key the database load drops, so 113 genotype-specific resistance
categories reached ``mutation_catalog`` with nothing saying which genotype they
belonged to.

WHAT IT WRITES
--------------
One row per catalogue entry x evidence reference. Every input column is carried
unchanged, except:

    alignment_name       replaced by ``genotype`` (AL_1a -> 1a); the two are
                         redundant and only the generic one is kept
    pubmed_id, DOI,      replaced by the evidence block:
    clinical_trials
      evidence_id          PMID, conference abstract code, or trial registry id
      data_source          pubmed | conference_abstract | clinical_trial
      evidence_url         DOI for publications, registry URL for NCT trials
      evidence_label       "Komatsu et al. 2017, Gastroenterology" or trial name(s)
      evidence_type        in_vitro | in_vivo | in_vitro;in_vivo
      linked_evidence_ids  paper -> the trials it reported for this entry;
                           trial -> the papers that reported it
      finding_ids          the PHDR resistance_finding ids behind the link

Numbers are written in their shortest exact form (``1.0`` -> ``1``,
``1807.2999999999997`` -> ``1807.3``), and stray whitespace in drug fields is
stripped. The input may be the normaliser's output or the older wide catalogue;
both give the same file. An entry with no PHDR evidence
(the decomposed combination components) is written once, evidence blank.

THE PHDR CHAIN
--------------
    alignment_ras_drug.id  (RAS : alignment : drug)
      <- resistance_finding.phdr_alignment_ras_drug_id
           .phdr_publication_id                        -> publication
           .phdr_in_vitro_result_id                    (no trial behind it)
           .phdr_in_vivo_result_id
               <- result_trial.phdr_in_vivo_result_id
                    .phdr_clinical_trial_id            -> clinical_trial

Check the output with generic/hcv/Tables/audit_catalog_linkage.py.
"""

import argparse
import collections
import csv
import re
import sys
from pathlib import Path

from BuildCatalogGenotypeColumns import alignment_to_genotype_code, build_dominant_residues, compute_columns
from evidence_sources import (
    CLINICAL_TRIAL, EVIDENCE_COLUMNS, VALUE_SEP, publication_source, trial_registry_id, trial_url,
)

csv.field_size_limit(min(sys.maxsize, 2 ** 31 - 1))

#: Wide-catalogue columns the evidence block replaces.
WIDE_EVIDENCE_COLUMNS = ('pubmed_id', 'DOI', 'clinical_trials')

#: Written after ``phenotype``, from phdr_drug.csv.
DRUG_COLUMNS = ('drug_producer', 'drug_category')
#: Written last, computed by BuildCatalogGenotypeColumns.
GENOTYPE_COLUMNS = ('relevant_genotypes', 'wild_type_residues')

#: Numeric columns, written in shortest exact form.
NUMERIC_COLUMNS = ('combination_size', 'component_order', 'numeric_resistance_category',
                   'any_in_vitro_evidence', 'in_vitro_max_ec50_midpoint', 'any_in_vivo_evidence',
                   'in_vivo_baseline', 'in_vivo_treatment_emergent')


def read_rows(path, delimiter=','):
    with open(path, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or []), [{k: (v or '') for k, v in r.items()} for r in reader]


def tidy_number(value):
    """``2.0`` -> ``2``; ``1807.2999999999997`` -> ``1807.3``; text and blanks unchanged."""
    try:
        number = round(float(value), 9)
    except (TypeError, ValueError):
        return value
    text = f'{number:.9f}'.rstrip('0').rstrip('.')
    return '0' if text in ('-0', '') else text


def _natural(text):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', text)]


def load_evidence(tables):
    """``{alignment_ras_drug_id: [evidence record, ...]}``, publications first, then trials."""
    tables = Path(tables)
    _, publications = read_rows(tables / 'phdr_publication.csv')
    _, trials = read_rows(tables / 'phdr_clinical_trial.csv')
    _, result_trials = read_rows(tables / 'phdr_result_trial.csv')
    _, findings = read_rows(tables / 'phdr_resistance_finding.csv')

    pub_by_id = {p['id']: p for p in publications}
    registry_of = {t['id']: trial_registry_id(t['nct_id'], t['id']) for t in trials}
    name_of = {t['id']: t['display_name'].strip() for t in trials}
    trials_by_result = collections.defaultdict(set)
    for rt in result_trials:
        trials_by_result[rt['phdr_in_vivo_result_id']].add(rt['phdr_clinical_trial_id'])

    def blank(evidence_id, source):
        return {'evidence_id': evidence_id, 'data_source': source, 'types': set(),
                'linked': set(), 'findings': set(), 'labels': set()}

    per_key = collections.defaultdict(dict)
    for f in findings:
        key, pub_id = f['phdr_alignment_ras_drug_id'], f['phdr_publication_id']
        in_vivo = f['phdr_in_vivo_result_id']
        bucket = per_key[key]

        pub = bucket.setdefault(('pub', pub_id), blank(pub_id, publication_source(pub_id)))
        pub['findings'].add(f['id'])
        if f['phdr_in_vitro_result_id']:
            pub['types'].add('in_vitro')
        if not in_vivo:
            continue
        pub['types'].add('in_vivo')
        for curator_id in trials_by_result.get(in_vivo, ()):
            registry_id = registry_of.get(curator_id, curator_id)
            trial = bucket.setdefault(('trial', registry_id), blank(registry_id, CLINICAL_TRIAL))
            trial['findings'].add(f['id'])
            trial['types'].add('in_vivo')
            trial['labels'].add(name_of.get(curator_id, curator_id))
            trial['linked'].add(pub_id)
            pub['linked'].add(registry_id)

    out = {}
    for key, bucket in per_key.items():
        pubs = sorted((v for (kind, _), v in bucket.items() if kind == 'pub'),
                      key=lambda v: (pub_by_id.get(v['evidence_id'], {}).get('year', ''), v['evidence_id']))
        trials_ = sorted((v for (kind, _), v in bucket.items() if kind == 'trial'),
                         key=lambda v: v['evidence_id'])
        records = []
        for v in pubs:
            p = pub_by_id.get(v['evidence_id'], {})
            label = ' '.join(x for x in (p.get('authors_short', ''), p.get('year', '')) if x)
            if p.get('journal'):
                label = f"{label}, {p['journal']}"
            records.append(_record(v, p.get('url', ''), label))
        for v in trials_:
            records.append(_record(v, trial_url(v['evidence_id']), VALUE_SEP.join(sorted(v['labels']))))
        out[key] = records
    return out


def _record(v, url, label):
    return {
        'evidence_id': v['evidence_id'],
        'data_source': v['data_source'],
        'evidence_url': url,
        'evidence_label': label,
        'evidence_type': VALUE_SEP.join(sorted(v['types'])),
        'linked_evidence_ids': VALUE_SEP.join(sorted(v['linked'])),
        'finding_ids': VALUE_SEP.join(sorted(v['findings'], key=_natural)),
    }


def output_fields(fieldnames):
    """Input header -> output header.

    ``genotype`` takes ``alignment_name``'s place, drug details follow
    ``phenotype``, then the evidence block, then the genotype columns. Wide
    evidence lists and any earlier genotype columns are dropped and rebuilt, so
    the normaliser's output and the older wide catalogue give the same header.
    """
    rebuilt = set(WIDE_EVIDENCE_COLUMNS) | set(DRUG_COLUMNS) | set(GENOTYPE_COLUMNS)
    out = []
    for name in fieldnames:
        if name == 'alignment_name':
            out.append('genotype')
        elif name not in rebuilt:
            out.append(name)
        if name == 'phenotype':
            out.extend(DRUG_COLUMNS)
    for name in ('genotype', *DRUG_COLUMNS):
        if name not in out:
            out.append(name)
    return out + EVIDENCE_COLUMNS + list(GENOTYPE_COLUMNS)


def build(catalog_path, tables, output_path):
    tables = Path(tables)
    fieldnames, rows = read_rows(catalog_path, delimiter='\t')
    fields = output_fields(fieldnames)

    _, drugs = read_rows(tables / 'phdr_drug.csv')
    drug_by_id = {d['id'].strip(): d for d in drugs}
    relevant, wild_type, _ = compute_columns(
        rows, build_dominant_residues(str(tables / 'phdr_alignment_typical_aa.csv')), {}, True, {})
    evidence = load_evidence(tables)

    written = 0
    with open(output_path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n',
                                extrasaction='ignore')
        writer.writeheader()
        for row, rel, wt in zip(rows, relevant, wild_type):
            entry = dict(row)
            entry['genotype'] = alignment_to_genotype_code(row.get('alignment_name'))
            entry['relevant_genotypes'], entry['wild_type_residues'] = rel, wt
            for column in NUMERIC_COLUMNS:
                if column in entry:
                    entry[column] = tidy_number(entry[column])
            entry['drug'] = entry.get('drug', '').strip()
            drug = drug_by_id.get(entry['drug'], {})
            entry['drug_producer'] = (drug.get('producer') or '').strip()
            entry['drug_category'] = (drug.get('drug_category') or '').strip()
            for record in evidence.get(row.get('id', ''), []) or [dict.fromkeys(EVIDENCE_COLUMNS, '')]:
                writer.writerow({**entry, **record})
                written += 1

    print(f'Wrote {output_path}: {written} rows from {len(rows)} catalogue entries, {len(fields)} columns')
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--catalog', required=True,
                        help="NormaliseHcvMutationCatalog.py's output (generalized_mutation_catalog.tsv)")
    parser.add_argument('--tables', required=True, help='directory holding the phdr_*.csv tables')
    parser.add_argument('--output', required=True, help='long-format catalogue TSV to write')
    args = parser.parse_args(argv)
    return build(Path(args.catalog), Path(args.tables), Path(args.output))


if __name__ == '__main__':
    raise SystemExit(main())
