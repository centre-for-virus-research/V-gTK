#!/usr/bin/env python3
"""Fold per-genotype knowledge out of the PHDR source tables into the catalog TSV.

WHY THIS EXISTS
---------------
The pipeline must not depend on ``alignment_name``. Alignment codes (``AL_1a``,
``AL_6xd``) are a PHDR/HCV artefact; a rabies or influenza catalog will never
have them. So everything genotype-related is resolved *here*, at catalog build
time, into two generic columns that any virus can supply by any means:

``relevant_genotypes``
    Which genotypes/subtypes this signature was actually scored in, as a plain
    semicolon-separated list.

        1a;1b
        3                       a bare genotype
                                empty - applies to every genotype

``wild_type_residues``
    The dominant residue at this position for each genotype the entry can be
    applied to, again semicolon-separated, with an optional ``:frequency``.

        1a:Q:60.89;1b:R:92.26   with frequency
        1a:Q;1b:R               without

Both columns are optional. A catalog lacking them behaves exactly as before -
no genotype gating, no wild-type suppression - so non-HCV viruses need supply
nothing, and are never expected to have subtypes at all.

WHY TWO SOURCES
---------------
``phdr_alignment_typical_aa.csv`` gives the empirical dominant residue per
alignment, but it only records residues at >=10% frequency. 770 of the 807
catalogued (RAS, alignment) pairs sit below that, so it cannot say anything
about them, and it is strictly per-position so it says nothing about
combinations at all.

``var_almt_note.csv`` covers that regime: 71,016 rows of observed frequency for
every variation in every alignment, combinations included. Neither table alone
is enough; together they cover the catalog completely.

WHY THE DOMINANT RESIDUE AND NOT EVERY TYPICAL ONE
--------------------------------------------------
Suppressing every residue the alignment calls "typical" would discard real
findings. At NS3:80 in genotype 1a the residues are Q at 60.9% and K at 36.0%:
K is common, but Q80K is a genuine simeprevir resistance substitution. Marking
only the dominant residue as wild type suppresses Q80Q while still calling
Q80K. Recording the frequency alongside is what lets a reader see that the
"wild type" at that position is only a 61% majority.
"""

import argparse
import collections
import csv
import os
import re
import sys

csv.field_size_limit(min(sys.maxsize, 2 ** 31 - 1))

#: Separator between entries in both generated columns.
ENTRY_SEP = ';'
#: Separator between the fields of one entry.
FIELD_SEP = ':'

#: ``AL_1a`` -> ``1a``. Anything not shaped like an alignment code is ignored,
#: which is what keeps a non-PHDR catalog from acquiring nonsense genotypes.
_ALIGNMENT_CODE = re.compile(r'^(?:AL_)?(?P<code>[0-9][0-9A-Za-z]*)$', re.IGNORECASE)

#: The leading digits of a subtype are its genotype: 6n -> 6, 1a -> 1.
_GENOTYPE_DIGITS = re.compile(r'^(\d+)')


def alignment_to_genotype_code(alignment_name):
    """``AL_6xd`` -> ``6xd``; ``AL_6_unassigned_JX183558`` and ``AL_MASTER`` -> None.

    The ``_unassigned_`` pseudo-alignments each describe a single sequence, not
    a genotype, so they are deliberately excluded.
    """
    text = '' if alignment_name is None else str(alignment_name).strip()
    if text.lower() in ('', 'nan', 'none'):
        return ''
    match = _ALIGNMENT_CODE.match(text)
    return match.group('code').lower() if match else ''


def genotype_of(code):
    """Leading digits of a subtype code, or '' when it has none."""
    match = _GENOTYPE_DIGITS.match(code or '')
    return match.group(1) if match else ''


def read_csv(path):
    with open(path, newline='', encoding='utf-8', errors='replace') as handle:
        return list(csv.DictReader(handle))


def build_variant_frequencies(var_almt_note_path):
    """(variation_name, genotype_code) -> observed frequency percentage.

    ``var_almt_note`` reports ``ncbi_curated_frequency`` already as a
    percentage, alongside the present/absent counts that produced it. The
    counts are kept because a frequency computed from 7 sequences means
    something very different from one computed from 20,456, and any rule that
    uses the number needs to be able to tell those apart.
    """
    frequencies = {}
    counts = {}
    if not var_almt_note_path or not os.path.isfile(var_almt_note_path):
        return frequencies, counts
    for row in read_csv(var_almt_note_path):
        code = alignment_to_genotype_code(row.get('alignment_name'))
        if code is None:
            continue
        name = (row.get('variation_name') or '').strip()
        if not name:
            continue
        try:
            frequencies[(name, code)] = float(row.get('ncbi_curated_frequency') or 0.0)
        except ValueError:
            continue
        try:
            present = int(row.get('ncbi_curated_total_present') or 0)
            absent = int(row.get('ncbi_curated_total_absent') or 0)
            counts[(name, code)] = present + absent
        except ValueError:
            pass
    return frequencies, counts


def build_dominant_residues(typical_aa_path):
    """(genotype_code, feature, position) -> (residue, percentage).

    The dominant residue only - see the module docstring for why every typical
    residue would be the wrong thing to record.
    """
    best = {}
    if not typical_aa_path or not os.path.isfile(typical_aa_path):
        return best
    for row in read_csv(typical_aa_path):
        code = alignment_to_genotype_code(row.get('alignment_name'))
        if code is None:
            continue
        key = (code, (row.get('feature_name') or '').strip(), (row.get('codon_label') or '').strip())
        try:
            pct = float(row.get('pct_members') or 0.0)
        except ValueError:
            continue
        residue = (row.get('aa_residue') or '').strip()
        if not residue:
            continue
        if key not in best or pct > best[key][1]:
            best[key] = (residue, pct)
    return best


def format_entries(entries, with_frequency):
    """``[('1a', 'Q', 60.89)]`` -> ``1a:Q:60.89`` (or ``1a:Q`` without frequency)."""
    out = []
    for entry in entries:
        code = entry[0]
        fields = [code] + [str(f) for f in entry[1:-1]]
        if with_frequency and entry[-1] is not None:
            fields.append(f'{entry[-1]:.2f}')
        out.append(FIELD_SEP.join(fields))
    return ENTRY_SEP.join(out)


def compute_columns(catalog_rows, dominant, frequencies, with_frequency=True):
    """Return ``relevant_genotypes`` and ``wild_type_residues`` per catalog row."""
    # Scope is a property of the SIGNATURE, not of one row: a signature scored
    # in both AL_1a and AL_1b must carry both on every one of its rows.
    signature_codes = collections.defaultdict(set)
    for row in catalog_rows:
        code = alignment_to_genotype_code(row.get('alignment_name'))
        if code:
            signature_codes[row.get('signature_id')].add(code)

    # Which subtypes exist at all, so a genotype can be expanded to its siblings.
    subtypes_by_genotype = collections.defaultdict(set)
    for (code, _feature, _pos) in dominant:
        subtypes_by_genotype[genotype_of(code)].add(code)

    relevant_out, wild_type_out = [], []
    for row in catalog_rows:
        codes = sorted(signature_codes.get(row.get('signature_id'), set()))
        # A plain list of genotype codes. Deliberately NO frequency: the number
        # that would go here is the frequency of the *mutant*, while the number
        # in wild_type_residues is the frequency of the *wild type*. Two
        # different quantities in identically-formatted adjacent columns is a
        # trap, so the frequency lives only where the residue it describes does.
        relevant_out.append(ENTRY_SEP.join(codes))

        # Wild types for every subtype that the genotype gate can route here -
        # an entry scored in 1b is applied to genotype-1 sequences, which
        # include 1a, 1c, 1l..., so each needs its own wild type or nothing is
        # suppressed for them.
        feature = (row.get('protein_name') or '').strip()
        position = (row.get('aa_position') or '').strip()
        siblings = set()
        for code in codes:
            siblings |= subtypes_by_genotype.get(genotype_of(code), set())
        if not codes:
            # An entry with no scope applies to every genotype, so it needs a
            # wild type for every genotype or nothing can ever suppress it.
            # These are not an edge case: the unscoped rows are precisely the
            # no-change anchors (R155R, K/Q80K) that decompose out of composite
            # signatures, and leaving them without a wild type is what lets
            # them fire on every wild-type sequence.
            siblings = {code for (code, feat, pos) in dominant
                        if feat == feature and pos == position}
        entries = []
        for code in sorted(siblings):
            hit = dominant.get((code, feature, position))
            if hit:
                entries.append((code, hit[0], hit[1]))
        wild_type_out.append(format_entries(entries, with_frequency))

    return relevant_out, wild_type_out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--catalog', required=True, help='catalog TSV to rewrite in place')
    parser.add_argument('--typical_aa', required=True, help='phdr_alignment_typical_aa.csv')
    parser.add_argument('--var_almt_note', default=None, help='var_almt_note.csv (optional; supplies frequencies)')
    parser.add_argument('--output', default=None, help='write here instead of in place')
    parser.add_argument('--no_frequency', action='store_true',
                        help='emit residues and genotypes without the optional frequency field')
    args = parser.parse_args(argv)

    with open(args.catalog, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    dominant = build_dominant_residues(args.typical_aa)
    frequencies, counts = build_variant_frequencies(args.var_almt_note)
    relevant, wild_type = compute_columns(rows, dominant, frequencies, not args.no_frequency)

    for column in ('relevant_genotypes', 'wild_type_residues'):
        if column not in fieldnames:
            fieldnames.append(column)
    for row, rel, wt in zip(rows, relevant, wild_type):
        row['relevant_genotypes'] = rel
        row['wild_type_residues'] = wt

    destination = args.output or args.catalog
    with open(destination, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    scoped = sum(1 for value in relevant if value)
    with_wt = sum(1 for value in wild_type if value)
    print(f'Wrote {destination}: {len(rows)} rows, {len(fieldnames)} columns')
    print(f'  relevant_genotypes populated : {scoped}/{len(rows)}')
    print(f'  wild_type_residues populated : {with_wt}/{len(rows)}')
    print(f'  dominant-residue table       : {len(dominant)} (genotype, feature, position) entries')
    print(f'  variant frequency table      : {len(frequencies)} (variation, genotype) entries')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
