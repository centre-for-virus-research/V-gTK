# info_help

Behaviour, design decisions and options, written so they can be checked without
reading the code.

Each document states not just *what* the pipeline does but *why it was decided
that way*, including the alternatives that were rejected and what they would have
cost. Where a number appears it was measured, and the document says what it was
measured against.

---

## Documents

| document | covers |
|---|---|
| [`hcv_mutation_annotation.md`](hcv_mutation_annotation.md) | How HCV resistance calls are made: the genotype scope rule, the wild-type suppression rule, why the reference used to be annotated against itself, and what the change measures out to. **Start here.** |
| [`mutation_catalog_columns.md`](mutation_catalog_columns.md) | The mutation catalogue column contract. What is required, what is optional, and what a non-HCV virus needs to supply (answer: nothing new). |
| [`genotype_provenance.md`](genotype_provenance.md) | Where a stored genotype or subtype came from — the `genotype_origin` vocabulary, its precedence, and why curation outranks inference. |
| [`guide_alignment_and_insertions.md`](guide_alignment_and_insertions.md) | What happens to a region present in a sub-reference but absent from the master, why insertions are lost without a guide alignment, and how to supply one for a segmented or non-segmented virus. |
| [`coordinate_spaces.md`](coordinate_spaces.md) | The three coordinate spaces a CDS passes through on its way into `features`, the 5' trim between a master's GFF and its alignment row, what each coordinate column actually holds, and the genome shapes the projection cannot yet handle. |
| [`hcv_tree_validation_ictv.md`](hcv_tree_validation_ictv.md) | Validation of the pipeline's HCV trees against the ICTV reference set: genotype and subtype monophyly, backbone agreement and where support runs out, and two metadata subtype errors the tree catches. |
| [`reference_recombination_screening.md`](reference_recombination_screening.md) | Screening the curated reference set for recombination with UShER/RIPPLES: why only the reference set is tractable, the 16-thread ceiling, chunking, and what gets stored. |
| [`tree_rooting.md`](tree_rooting.md) | How trees are rooted: midpoint everywhere, why HCV genotype 8 was evaluated as an outgroup and rejected, why single-tip and MMseqs-cluster outgroups do not work, and what happens to UShER node labels. |
| [`lineage_growth_rates.md`](lineage_growth_rates.md) | Estimating how fast a lineage or an allele is growing: which estimators are implemented and what each assumes, why node dating on an UShER tree is a bound rather than an estimate, what the denominator actually is, and how independent origins separate a mutation's effect from its lineage's. |
| [`cli_reference.md`](cli_reference.md) | Every command-line option on the scripts in this path, what each does, and how to turn each behaviour off. |

---

## The short version

**Mutation calling is gated on genotype and on wild type.** A catalogued
substitution is only called when the entry applies to the sequence's genotype,
and when the observed residue actually differs from what is normal for that
genotype. Before this, the reference sequence H77 — isolated in 1977, decades
before any antiviral it could be resistant to — carried 17 resistance calls
against itself. It now carries none, and no genuine finding was lost.

**The pipeline never reads an alignment column.** `AL_1a` and friends are a
PHDR/HCV curation artefact. All genotype knowledge is resolved into two generic,
optional TSV columns at catalogue build time, so the annotator works the same for
any virus and a virus without subtypes needs to supply nothing.

**Suppression is on the dominant residue, not on every typical one.** At
`NS3:80` in genotype 1a the residues are Q at 60.9% and K at 36.0%. Q80K is both
a common polymorphism and a real simeprevir resistance substitution. Suppressing
everything "typical" would have silenced it; suppressing only the dominant
residue keeps it. This is the single most consequential decision in the design
and there is a standing test on it.

**Evidence is attached per genotype, not per mutation.** Each catalogue row
carries the registered clinical trials (NCT identifiers) and publications
supporting it *in that genotype, for that drug*, and both reach the database on
`mutation_catalog`. `NS5A:31M` against daclatasvir has one trial in genotype 1a
and nine in 1b; flattening them would attach one genotype's evidence to
another's call. The `publications` and `clinical_trials` tables resolve those
identifiers to something readable, and each is keyed one row per identifier, so
joining to them cannot double the evidence a query reports.

**Nothing is suppressed silently.** Every evaluated call, emitted or not, lands
in `sequence_mutation_calls` with the reason. A call that was gated out, one
whose residue was wild type, and one that was never evaluated are three different
things and the table distinguishes them.

---

## Scope

These documents describe the HCV mutation path and the genotype provenance that
supports it. They do not describe the pipeline as a whole — see the repository
`README.md` for that, and `TESTING.md` for how to run the suite.
