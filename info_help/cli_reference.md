# CLI and pipeline options: the HCV mutation path

Every switch that affects mutation annotation, and what it does. Options are
listed as the scripts actually declare them.

---

## Pipeline parameters (`nextflow.config`)

| parameter | default | effect |
|---|---|---|
| `mutation_catalog` | `null` | Path to the catalogue TSV. **This is the master switch.** When null, `ANNOTATE_MUTATIONS` and `VERIFY_MUTATIONS` never run and no mutation tables are created. Set for `HCV_*` profiles only. |
| `mutation_virus` | `null` | Virus context passed to the annotator as `--virus`. Set to `HCV` on the HCV profiles. Selects virus-specific handling; anything else takes the generic path. |

Set together on the HCV profiles:

```groovy
mutation_catalog = "${projectDir}/generic/hcv/Tables/generalized_mutation_catalog_with_extra_info.tsv"
mutation_virus   = "HCV"
```

`VERIFY_MUTATIONS` is gated separately in the workflow body — it runs only when
`params.mutation_catalog` is a non-null, non-empty, existing path.

---

## `scripts/AnnotateMutations.py`

Annotates a finished database against a catalogue. Run by the
`ANNOTATE_MUTATIONS` step.

| option | required | default | effect |
|---|---|---|---|
| `--db` | yes | — | SQLite database to annotate, modified in place. |
| `--mutation_catalog` | yes | — | Catalogue TSV. |
| `--virus` | no | `""` | Virus context, e.g. `HCV`. Selects virus-specific handling. |

```
python scripts/AnnotateMutations.py \
    --db test_out/HCV_OM_test/HCV_OM_test.db \
    --mutation_catalog generic/hcv/Tables/generalized_mutation_catalog_with_extra_info.tsv \
    --virus HCV
```

**Behaviour controlled by the catalogue, not by flags.** Genotype gating and
wild-type suppression have no command-line switch. They activate when the
catalogue carries `relevant_genotypes` and `wild_type_residues`, and are inert
when it does not. This is deliberate: whether gating applies is a property of the
data, not of the invocation, so the same command produces the same result for a
given catalogue.

To turn gating off, use a catalogue without those columns.

**Writes:** `mutation_catalog`, `sequence_relevant_mutation_summary`,
`completed_signatures_only`, and `sequence_mutation_calls` (every evaluated call
with the reason it was emitted or suppressed).

---

## `scripts/BuildCatalogGenotypeColumns.py`

Folds per-genotype knowledge out of the PHDR source tables into the catalogue.
**HCV-specific** — it is the only place alignment codes are understood.

| option | required | default | effect |
|---|---|---|---|
| `--catalog` | yes | — | Catalogue TSV, rewritten in place unless `--output` is given. |
| `--typical_aa` | yes | — | `phdr_alignment_typical_aa.csv`. Supplies the dominant residue per alignment. |
| `--var_almt_note` | no | none | `var_almt_note.csv`. Supplies observed frequencies. Without it the columns are written without the optional frequency field. |
| `--output` | no | in place | Write elsewhere instead of rewriting the input. |
| `--no_frequency` | no | off | Emit `1a:Q;1b:R` rather than `1a:Q:60.89;1b:R:92.26`. Smaller; loses the ability to see that a wild type is only a bare majority. |

```
python scripts/BuildCatalogGenotypeColumns.py \
    --catalog       generic/hcv/Tables/generalized_mutation_catalog_with_extra_info.tsv \
    --typical_aa    generic/hcv/Tables/phdr_alignment_typical_aa.csv \
    --var_almt_note generic/hcv/Tables/var_almt_note.csv
```

Prints how many rows got each column, and the size of the two lookup tables it
built. It touches only the two generated columns.

---

## `scripts/NormalizeHcvMutationCatalog.py`

Builds the generalized catalogue from the PHDR source tables. Run by hand when
the PHDR export is refreshed, not by the pipeline.

| option | default | effect |
|---|---|---|
| `--variation` | `generic/hcv/Tables/variation.csv` | Variation definitions. |
| `--variation_metatag` | `generic/hcv/Tables/variation_metatag.csv` | Conjunction membership. |
| `--phdr_alignment_ras` | `generic/hcv/Tables/phdr_alignment_ras.csv` | RAS × alignment scope. |
| `--phdr_alignment_ras_drug` | `generic/hcv/Tables/phdr_alignment_ras_drug.csv` | Drug and resistance category. |
| `--gene_info` | `generic/hcv/Tables/gene_info.tsv` | Gene names. |

After regenerating, re-run `BuildCatalogGenotypeColumns.py` — the generated
columns are not reproduced by the normaliser's own inputs.

---

## `scripts/VerifyMutations.py`

Independent check of a database's annotations against the catalogue. Run by
`VERIFY_MUTATIONS` and useful standalone.

| option | default | effect |
|---|---|---|
| `--db` | a hard-coded HCV path | Database to verify. **Override this** — the default points at one developer's build. |
| `--mutation_catalog` | the shipped HCV catalogue | Catalogue to verify against. |
| `--sample_size` | `100` | Annotated sequences to sample. Raise for a fuller check; the cost is linear. |
| `--seed` | `42` | Sampling seed. Fixed so a run is reproducible. |
| `--min_identity` | `0.65` | Minimum nucleotide identity for a sequence's alignment to be trusted. Below this the sequence is reported as failing alignment validation rather than being scored. |
| `--hcv_test_ns3_36a` | off | Narrow diagnostic for NS3:36A on a specific accession set. Development aid. |

---

## Reading the results

Start with `sequence_mutation_calls`. A missing call is either
`suppressed_out_of_scope`, `suppressed_wild_type`, or absent because nothing
matched — and only that table distinguishes them.

```sql
-- why did this sequence lose a call?
SELECT signature_id, call_status, scope_tier, residue_status,
       observed_residue, wild_type_residues, sequence_genotype, sequence_subtype
FROM sequence_mutation_calls
WHERE primary_accession = 'NC_004102';

-- overall shape of a run
SELECT call_status, scope_tier, residue_status, COUNT(*)
FROM sequence_mutation_calls GROUP BY 1,2,3 ORDER BY 4 DESC;
```

---

## Turning behaviour off

| to disable | do this |
|---|---|
| mutation annotation entirely | leave `mutation_catalog` null |
| genotype gating | use a catalogue without `relevant_genotypes` |
| wild-type suppression | use a catalogue without `wild_type_residues` |
| the frequency field | rebuild with `--no_frequency` |

There is no flag that half-disables a rule. Either the catalogue carries the
column and the rule applies, or it does not and behaviour is what it was before
gating existed.
