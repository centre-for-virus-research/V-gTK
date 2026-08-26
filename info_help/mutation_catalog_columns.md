# Mutation catalogue column contract

What the pipeline reads from a mutation catalogue, what is optional, and what a
non-HCV virus is expected to supply. The short answer for a new virus is:
**nothing beyond the columns you already have.**

---

## The generic columns

Two columns carry all genotype-related knowledge. **Both are optional.** A
catalogue with neither behaves exactly as it did before genotype gating existed:
every entry applies to every sequence and nothing is ever suppressed.

### `relevant_genotypes`

Which genotypes or subtypes this signature was actually curated in.

```
relevant_genotypes
1a;1b        two subtypes
1a;3;3a      a bare genotype alongside subtypes
3            a bare genotype, no subtype
             empty: applies to every genotype
```

- **Semicolon** between entries. Nothing else — this column is a plain list.
- **No frequency here.** A frequency attached to a genotype in this column would
  be the frequency of the *mutant*, whereas the number in `wild_type_residues`
  is the frequency of the *wild type*. Two different quantities in
  identically-formatted adjacent columns is a trap, so the frequency lives only
  alongside the residue it describes.
- Values are free text. `1a`, `3`, `L`, `Asian` are all acceptable; the only
  structure the pipeline imposes is that **leading digits are the genotype**, so
  `6n` and `6xd` both belong to genotype `6`. A label with no leading digits is
  matched only as an exact string.

### `wild_type_residues`

The dominant residue at this position for each genotype the entry can reach.

```
wild_type_residues
1a:Q;1b:R                       residues only
1a:Q:60.89;1b:R:92.26           with the optional frequency
                                empty: nothing is suppressed for this row
```

- Same separators. The trailing frequency is again **optional**.
- The residue is the **dominant** one, not every typical one. See
  [`hcv_mutation_annotation.md` §3](hcv_mutation_annotation.md) for why that
  distinction is clinically load-bearing.
- Coverage should include every subtype the genotype gate can route to this row —
  an entry scoped to `1b` is applied to genotype-1 sequences, so `1a`, `1c` and
  `1l` each need their own entry or nothing will be suppressed for them.

### `clinical_trials`

The registered clinical trials supporting this finding, **for this row's
genotype and drug**.

```
clinical_trials
NCT01717326;NCT02092350;NCT02105454
                        empty: no in-vivo trial evidence for this genotype/drug
```

- Semicolon separated, plain list, same as `relevant_genotypes`.
- Keyed per **(mutation, genotype, drug)** — not per mutation. Trial support
  genuinely varies by genotype: `NS5A:31M` against daclatasvir cites one trial in
  genotype 1a and **nine** in 1b. 49 of the 64 (RAS, drug) pairs curated in more
  than one genotype have different trial sets, so unioning them would attach 1b's
  evidence to a 1a call.
- Only in-vivo findings carry a trial. An in-vitro EC50 has no trial behind it
  and correctly gets nothing — 826 of 1,740 findings have an in-vivo result.
- Populated on 746 of 1,869 rows, drawing on 97 distinct trials.

Resolve an NCT identifier to a trial name through the `clinical_trials` table in
the database, loaded with `--clinical_trials`.

---

## What a non-HCV virus needs

| virus shape | `relevant_genotypes` | `wild_type_residues` |
|---|---|---|
| no genotypes at all (rabies, most) | omit the column | omit the column |
| genotypes but no subtypes | `1;2;3` | `1:Q;2:R` or `1:Q:98.1;2:R:95.4` |
| genotypes and subtypes (HCV, influenza) | `1a;1b` | `1a:Q;1b:R` or `1a:Q:60.89;1b:R:92.26` |

`clinical_trials` is optional everywhere and is omitted unless a virus has a
trial registry to link against.

**Subtypes are never required.** The genotype is the leading digits of whatever
you write; if you write `1` the genotype is `1` and there is no subtype. Nothing
in the pipeline demands a subtype column, and nothing infers one.

Omitting both columns is a supported, tested configuration — it is what every
non-HCV profile does today.

---

## Columns the pipeline does NOT read

| column | why it is not read |
|---|---|
| `alignment_name` | PHDR/HCV artefact (`AL_1a`, `AL_6xd`). No other virus has it. Resolved into `relevant_genotypes` at build time. |
| `display_structure` | Same. Its wild types are a curator's shorthand and, at `NS3:80`, spell `K/Q80K` — which would suppress the simeprevir RAS Q80K. Superseded by `wild_type_residues`. |

Both are still **carried** in the HCV catalogue as provenance, and both are
useful when auditing where a value came from. Neither is consulted at runtime.
This is deliberate: making the pipeline depend on an alignment column would tie
it to one virus's curation format.

---

## Residue spellings

| meaning | write | also accepted on read |
|---|---|---|
| amino acid | one uppercase letter | lower case, surrounding whitespace |
| stop codon | `*` | `_` (legacy) |
| deletion | `-` | `del` (legacy) |

Legacy spellings are read so that an older catalogue still loads; the standard
spellings are what gets written.

A deletion is **not** detected by translating a codon — there is nothing to
translate. It is detected as a gap in the aligned columns, and only when that gap
lies strictly inside the sequence's covered span. A gap in terminal padding is
missing data, not a deletion, and is read as `X`.

---

## Required columns

Unchanged by this work, listed for completeness:

```
protein_name  segment  aa_position  alt_residue  reference_accession
mutation_id   mutation_type  signature_id  signature_kind
combination_id  combination_size  phenotype
```

`phenotype` is worth a note: `_build_output_row` produced it but `output_fields`
omitted it, so every regeneration of the catalogue silently dropped a column the
annotator requires. Both it and `relevant_genotypes` are now declared.

---

## Regenerating the HCV columns

```
python scripts/BuildCatalogGenotypeColumns.py \
    --catalog       generic/hcv/Tables/generalized_mutation_catalog_with_extra_info.tsv \
    --typical_aa    generic/hcv/Tables/phdr_alignment_typical_aa.csv \
    --var_almt_note generic/hcv/Tables/var_almt_note.csv
```

Rewrites the catalogue in place, adding or refreshing the two columns and
touching nothing else — verified by comparing every cell of the 33 pre-existing
columns before and after.

See [`cli_reference.md`](cli_reference.md) for the full option list.
