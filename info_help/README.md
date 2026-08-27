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

**Evidence is attached per genotype, not per mutation.** Each call carries the
registered clinical trials (NCT identifiers) and publications supporting it *in
that genotype, for that drug*. `NS5A:31M` against daclatasvir has one trial in
genotype 1a and nine in 1b; flattening them would attach one genotype's evidence
to another's call.

**Nothing is suppressed silently.** Every evaluated call, emitted or not, lands
in `sequence_mutation_calls` with the reason. A call that was gated out, one
whose residue was wild type, and one that was never evaluated are three different
things and the table distinguishes them.

---

## Scope

These documents describe the HCV mutation path and the genotype provenance that
supports it. They do not describe the pipeline as a whole — see the repository
`README.md` for that, and `TESTING.md` for how to run the suite.
