# CLI and pipeline options: the HCV mutation path

Every switch that affects mutation annotation, and what it does. Options are
listed as the scripts actually declare them.

---

## Pipeline parameters (`nextflow.config`)

| parameter | default | effect |
|---|---|---|
| `mutation_catalog` | `null` | Path to the catalogue TSV. **This is the master switch.** When null, `ANNOTATE_MUTATIONS` and `VERIFY_MUTATIONS` never run and no mutation tables are created. Set for `HCV_*` profiles only. |
| `mutation_virus` | `null` | Virus context passed to the annotator as `--virus`. Set to `HCV` on the HCV profiles. Selects virus-specific handling; anything else takes the generic path. |
| `mutation_publications` | `null` | Optional publication metadata CSV. |
| `mutation_clinical_trials` | `null` | Optional trial registry CSV. |

Set together on the HCV profiles:

```groovy
mutation_catalog = "${projectDir}/generic/hcv/Tables/generalized_mutation_catalog_evidence_linked.tsv"
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
| `--publications` | no | none | Publication metadata CSV. Loads a `publications` table, **one row per `evidence_id`**, with `data_source` (`pubmed` or `conference_abstract`), title, authors, year, journal and url, so publication references in `mutation_catalog.evidence_id` resolve. Without it they stay bare references. |
| `--clinical_trials` | no | none | Trial registry CSV (`id,display_name,nct_id`). Loads a `clinical_trials` table, **one row per `evidence_id`** (the NCT number, or the registry's own id when there is none, e.g. `UMIN000015627`), so trial references in `mutation_catalog.evidence_id` resolve without the join fanning out. Curator rows sharing a registry id are merged, keeping every distinct id and name semicolon separated. A repeated `evidence_id` in either table stops the write. |

```
python scripts/AnnotateMutations.py \
    --db test_out/HCV_OM_test/HCV_OM_test.db \
    --mutation_catalog generic/hcv/Tables/generalized_mutation_catalog_evidence_linked.tsv \
    --virus HCV \
    --publications generic/hcv/Tables/phdr_publication.csv \
    --clinical_trials generic/hcv/Tables/phdr_clinical_trial.csv
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
with the reason it was emitted or suppressed). With `--publications` and
`--clinical_trials` it also writes the `publications` and `clinical_trials`
lookup tables that `mutation_catalog.evidence_id` joins to (see
[`mutation_catalog_columns.md`](mutation_catalog_columns.md#evidence-columns)).

**Catalogue rows and calls.** A catalogue entry can span many rows: one per
genotype, drug and evidence reference. A call is made from a row's identity
(mutation, segment, residue, combination, signature), so rows sharing it are
matched once. That keeps a long-format catalogue from multiplying the work, and
it is why the `Mapping summary` counters count distinct catalogue identities, not
raw rows.

---

## `scripts/NormaliseHcvMutationCatalog.py`

Step 1 of the HCV catalogue build. Flattens the PHDR variation, combination and
drug tables into `generalized_mutation_catalog.tsv`: one row per mutation
component × genotype (alignment) × drug. Combinations are split into their
components, and PHDR protein names are mapped to canonical names via `gene_info.tsv`.
Run by hand when the PHDR export is refreshed, not by the pipeline.

| option | default | effect |
|---|---|---|
| `--variation` | `generic/hcv/Tables/variation.csv` | Variation definitions. |
| `--variation_metatag` | `generic/hcv/Tables/variation_metatag.csv` | Conjunction membership. |
| `--phdr_alignment_ras` | `generic/hcv/Tables/phdr_alignment_ras.csv` | RAS × alignment scope. |
| `--phdr_alignment_ras_drug` | `generic/hcv/Tables/phdr_alignment_ras_drug.csv` | Drug and resistance category. |
| `--gene_info` | `generic/hcv/Tables/gene_info.tsv` | Gene names. |
| `--output_path` | `generic/hcv/Tables/generalized_mutation_catalog.tsv` | Where to write. |

Then run `BuildHcvEvidenceCatalog.py`. The normaliser's output is not what the
pipeline reads.

---

## `scripts/BuildHcvEvidenceCatalog.py`

Step 2 of the HCV catalogue build. Turns the normaliser's output into
`generalized_mutation_catalog_evidence_linked.tsv`, the file the pipeline reads:

- one row per catalogue entry × evidence reference (`evidence_id`, `data_source`,
  `evidence_url`, `evidence_label`, `evidence_type`, `linked_evidence_ids`,
  `finding_ids`), following PHDR's finding → publication / in-vivo result → trial
  chain;
- `genotype` in place of `alignment_name`;
- `drug_producer` and `drug_category` from `phdr_drug.csv`;
- `relevant_genotypes` and `wild_type_residues`, computed with
  `BuildCatalogGenotypeColumns.py`;
- numbers in shortest exact form.

| option | required | effect |
|---|---|---|
| `--catalog` | yes | The normaliser's output. The older wide catalogue is also accepted and gives the same file. |
| `--tables` | yes | Directory holding the `phdr_*.csv` tables. |
| `--output` | yes | Catalogue TSV to write. |

```
python scripts/BuildHcvEvidenceCatalog.py \
    --catalog generic/hcv/Tables/generalized_mutation_catalog.tsv \
    --tables  generic/hcv/Tables \
    --output  generic/hcv/Tables/generalized_mutation_catalog_evidence_linked.tsv
python generic/hcv/Tables/audit_catalog_linkage.py \
    --catalog generic/hcv/Tables/generalized_mutation_catalog_evidence_linked.tsv
```

The audit re-derives every value and link from the PHDR tables, reports PASS /
WARN / FAIL per check, and exits 1 on any FAIL. `--json` writes the results as
JSON.

---

## `scripts/BuildCatalogGenotypeColumns.py`

Folds per-genotype knowledge out of the PHDR source tables into a catalogue.
**HCV-specific**: it is the only place alignment codes are understood.
`BuildHcvEvidenceCatalog.py` imports it; the command line remains for rewriting
a wide catalogue that still carries `alignment_name`.

| option | required | default | effect |
|---|---|---|---|
| `--catalog` | yes | — | Catalogue TSV with `alignment_name`, rewritten in place unless `--output` is given. |
| `--typical_aa` | yes | — | `phdr_alignment_typical_aa.csv`. Supplies the dominant residue per alignment. `AL_MASTER` and the single-sequence `AL_*_unassigned_*` alignments are skipped. |
| `--var_almt_note` | no | none | `var_almt_note.csv`. Supplies observed frequencies. |
| `--clinical_trial` | no | none | `phdr_clinical_trial.csv`. Supplies trial registry ids (NCT, or the registry's own id when there is none). |
| `--result_trial` | no | none | `phdr_result_trial.csv`. Links in-vivo results to trials. |
| `--resistance_finding` | no | none | `phdr_resistance_finding.csv`. All three trial tables are needed together to populate the wide `clinical_trials` column. |
| `--output` | no | in place | Write elsewhere instead of rewriting the input. |
| `--no_frequency` | no | off | Emit `1a:Q;1b:R` rather than `1a:Q:60.89;1b:R:92.26`. |

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
| trial resolution | leave `mutation_clinical_trials` null (the NCT ids stay, unresolved) |

There is no flag that half-disables a rule. Either the catalogue carries the
column and the rule applies, or it does not and behaviour is what it was before
gating existed.
