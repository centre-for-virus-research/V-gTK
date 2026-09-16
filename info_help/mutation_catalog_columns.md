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

### `genotype`

The genotype **this row** was curated in: the one its resistance category and
evidence belong to.

```
genotype
1a
        empty: an unscoped entry (a component split out of a combination)
```

- Different from `relevant_genotypes`, which is the **signature's** whole scope and
  drives the genotype gate. The same finding scored in 1a and 1b has two rows,
  `genotype` `1a` and `1b`, and both carry `relevant_genotypes` `1a;1b`.
- Without it the database cannot tell those rows apart. NS3:V107I against
  grazoprevir is `category_I` in 1a and `category_II` in 1b. Before the
  column existed, 113 such categories reached `mutation_catalog` unattributable.
- Replaces `alignment_name` (`AL_1a`), which said the same thing in PHDR's
  vocabulary. The two were redundant, so only the generic column is kept.

### Evidence columns

One evidence reference per row. A catalogue entry with three papers and two
trials has five rows, identical in every other column.

| column | holds |
|---|---|
| `evidence_id` | the reference: a PMID, a conference abstract code, or a trial registry id |
| `data_source` | `pubmed`, `conference_abstract` or `clinical_trial` |
| `evidence_url` | DOI for a publication; registry URL for an NCT trial; blank otherwise |
| `evidence_label` | human-readable: `Komatsu et al. 2017, Gastroenterology`, or the trial name(s) |
| `evidence_type` | `in_vitro`, `in_vivo` or `in_vitro;in_vivo`: which result columns the reference supports for this row |
| `linked_evidence_ids` | a paper's trials for this row, or a trial's papers; semicolon separated |
| `finding_ids` | the source findings behind the link (PHDR `resistance_finding` ids for HCV) |

- **Why not lists.** `pubmed_id` and `clinical_trials` used to be independent
  semicolon lists. A trial only ever appears behind an in-vivo result that a paper
  reports, so two lists lost which paper reported which trial (206 HCV entries),
  and which papers backed the in-vitro versus the in-vivo columns (208).
- **Scoped per row.** Evidence is keyed on (mutation, genotype, drug). `NS5A:31M`
  against daclatasvir cites one trial in 1a and ten in 1b.
- A row with no evidence has every evidence column blank.
- The ids and sources come from `scripts/evidence_sources.py`, shared by the
  catalogue builder and the database writer so the two always agree.

The columns reach `mutation_catalog` intact, and each lookup table is one row
per `evidence_id`, so resolving a reference is a plain equality join:

```sql
SELECT mc.mutation_id, mc.genotype, mc.drug, mc.resistance_category,
       mc.evidence_type, p.title, ct.trial_name, mc.linked_evidence_ids
FROM mutation_catalog mc
LEFT JOIN publications    p  ON mc.data_source IN ('pubmed', 'conference_abstract')
                             AND p.evidence_id = mc.evidence_id
LEFT JOIN clinical_trials ct ON mc.data_source = 'clinical_trial'
                             AND ct.evidence_id = mc.evidence_id
WHERE mc.mutation_id = 'NS3:107I' AND mc.drug = 'grazoprevir';
```

A trial is keyed on its registry id: the NCT number, or the registry's own id
when there is none (`UMIN000015627`, a Japanese UMIN-CTR registration that two
daclatasvir entries cite). PHDR files five NCT numbers under two curator ids
each (`ALLY-2` / `NCT02032888`, `ASTRAL-1` / `GS-US-342-1138`, trial arms such
as `Magellan-1, Part 1` and `Part 2`). Keying on the registry id collapses those,
and every curator name is kept in `evidence_label` and `clinical_trials.trial_name`.

---

## What a non-HCV virus needs

| virus shape | `relevant_genotypes` | `wild_type_residues` |
|---|---|---|
| no genotypes at all (rabies, most) | omit the column | omit the column |
| genotypes but no subtypes | `1;2;3` | `1:Q;2:R` or `1:Q:98.1;2:R:95.4` |
| genotypes and subtypes (HCV, influenza) | `1a;1b` | `1a:Q;1b:R` or `1a:Q:60.89;1b:R:92.26` |

`genotype` and the evidence columns are optional everywhere. A virus with
publication or trial evidence should supply the evidence columns and pass
`--publications` / `--clinical_trials` so the ids resolve inside the database;
either without the other is inert.

**Subtypes are never required.** The genotype is the leading digits of whatever
you write; if you write `1` the genotype is `1` and there is no subtype. Nothing
in the pipeline demands a subtype column, and nothing infers one.

Omitting both columns is a supported, tested configuration — it is what every
non-HCV profile does today.

---

## Columns the pipeline does NOT read

| column | why it is not read |
|---|---|
| `display_structure` | PHDR/HCV artefact. Its wild types are a curator's shorthand and, at `NS3:80`, spell `K/Q80K`, which would suppress the simeprevir RAS Q80K. Superseded by `wild_type_residues`. |
| `id`, `phdr_alignment_ras_id`, `phdr_drug_id`, `source_*` | PHDR keys, carried so any row can be traced back to the source tables. |

These are **carried** in the HCV catalogue as provenance and are useful when
auditing where a value came from, but none is consulted at runtime. `alignment_name`
is no longer carried at all: `genotype` holds the same information in generic
form. This is deliberate. Making the pipeline depend on an alignment column
would tie it to one virus's curation format.

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

## Regenerating the HCV catalogue

Two deterministic steps from the PHDR export. Nothing is added by hand:

```
python scripts/NormaliseHcvMutationCatalog.py      # -> generic/hcv/Tables/generalized_mutation_catalog.tsv
python scripts/BuildHcvEvidenceCatalog.py \
    --catalog generic/hcv/Tables/generalized_mutation_catalog.tsv \
    --tables  generic/hcv/Tables \
    --output  generic/hcv/Tables/generalized_mutation_catalog_evidence_linked.tsv
python generic/hcv/Tables/audit_catalog_linkage.py \
    --catalog generic/hcv/Tables/generalized_mutation_catalog_evidence_linked.tsv
```

The second step adds the drug columns, `genotype`, the evidence block, and
`relevant_genotypes` / `wild_type_residues` (computed with
`BuildCatalogGenotypeColumns.py`). The audit re-derives every value and link from
the PHDR tables and exits 1 on any disagreement. A unit test asserts the shipped
file is byte-identical to a fresh build.

`generalized_mutation_catalog_with_extra_info.tsv` is the earlier hand-assembled
wide catalogue. It is no longer read by anything in the pipeline. Built through the
second step, it gives the same file as the normaliser's output.

See [`cli_reference.md`](cli_reference.md) for the full option list.
