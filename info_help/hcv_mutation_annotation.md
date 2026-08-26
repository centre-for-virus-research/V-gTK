# HCV mutation annotation: behaviour, design decisions, and options

Everything the HCV mutation path does, why it does it, and every switch that
changes it — written so the decisions are trackable without reading the code.

Companion documents in this folder:

- [`mutation_catalog_columns.md`](mutation_catalog_columns.md) — the exact column contract, including
  what a non-HCV virus must and must not supply.
- [`genotype_provenance.md`](genotype_provenance.md) — where a stored genotype came from.
- [`cli_reference.md`](cli_reference.md) — every command-line option on the scripts involved.

---

## 1. What the pipeline is doing

For each query sequence the annotator asks, position by position, whether the
residue it observes matches a catalogued resistance-associated substitution
(RAS). A match becomes a **call**. Calls roll up into two summary tables and one
per-call evidence table in the SQLite database.

Until this revision the question it asked was simply *"does this sequence have
residue X at position N?"* — with no reference to what residue is normal there,
and no reference to which genotype the catalogue entry was curated in. Both
omissions produced systematic false positives, described in §2.

The question it asks now is:

> Is this catalogue entry applicable to **this sequence's genotype**, and is the
> observed residue **different from the wild type for that genotype**?

Both halves must hold for a call to be emitted.

---

## 2. The two problems this fixes

### 2.1 The reference was annotated against itself

`NC_004102` is strain **H77**, HCV genotype 1a, isolated in **1977** — thirty-four
years before the first direct-acting antiviral was approved. It is the sequence
HCV positions are numbered against. It cannot carry treatment-emergent
resistance.

It was nonetheless annotated with **17 "relevant mutations present"** and 22
completed signatures, including `NS5A:30Q` and `NS5A:31L`, which are H77's own
residues.

Two mechanisms produced this.

**Wild-type anchors.** PHDR curates composite findings such as
`Y56H+K/Q80K+R155R+D168V`. Only `Y56H` and `D168V` are substitutions; `R155R` and
`K/Q80K` record the *background* they were observed on. The catalogue build
decomposes composites into standalone single-residue entries and keeps only the
alt residue, so `R155R` becomes a rule reading "R at 155" — which matches every
untreated genotype-1 sequence.

**Cross-genotype application.** All 1,827 catalogue rows are scored against the
single reference `REF_MASTER_NC_004102` (genotype 1a), but the catalogue's own
entries were curated in 23 different genotype buckets. `K30Q` is a genuine
substitution in a genotype where position 30 is K — but in genotype 1a the wild
type *is* Q, so the same rule fires on every untreated 1a sequence. Measured on
`HCV_OM_test.db`, 85% of calls landed on a genotype the entry was never scored in.

### 2.2 Catalogued deletions were undetectable

Twelve catalogue rows are NS5A deletions at positions 29, 30 and 32 — every one
classed `category_I` or `category_II` against daclatasvir, pibrentasvir,
elbasvir, ledipasvir or velpatasvir. An in-frame deletion in NS5A domain I
removes a residue the inhibitor binds and confers broad class-wide resistance.

They were spelled `del` in the catalogue and compared against a translated
codon, which is always a single character. The comparison could never succeed
and emitted no diagnostic, so a sequence carrying `NS5A:32del` was reported as
having no NS5A resistance at all.

---

## 3. The rules, exactly

### Rule A — genotype scope, matched at genotype level

A catalogue entry applies to a sequence when the sequence's **genotype** (the
leading digits of its subtype) matches the genotype of any code in the entry's
`relevant_genotypes`.

```
entry relevant_genotypes = 1b      sequence 1a   -> in scope   (both genotype 1)
entry relevant_genotypes = 1b      sequence 6n   -> out of scope
entry relevant_genotypes = (empty) sequence 6n   -> in scope   (unscoped = any)
```

**Why genotype level and not exact subtype.** Observed sequence subtypes include
`6n`, `6xd`, `6xc` and `4v`, none of which the RAS catalogue has an entry for.
Exact-subtype matching keeps only **198 of 5,817 calls (3.4%)**; genotype-level
keeps 803. The catalogue is simply not curated at the granularity the data is
classified at.

The matched tier is recorded per call as `scope_tier`: `subtype` when the
sequence's exact subtype appears in the entry, `genotype` when only the genotype
digits matched, `unscoped` when the entry has no scope at all.

### Rule B — wild-type suppression

A call is emitted only when the observed residue **differs from the dominant
residue for that sequence's genotype** at that position.

```
sequence 1a, observed Q at NS5A:30, wild type for 1a is Q  -> suppressed (anchor)
sequence 1b, observed Q at NS5A:30, wild type for 1b is R  -> emitted    (change)
sequence 8a, observed Q at NS5A:30, no wild type known     -> emitted    (wt_unknown)
```

An unknown wild type **never suppresses**. "We do not know" is not "the residue
is normal", and conflating them would silently hide findings in exactly the
genotypes with the least reference data.

Wild-type lookup follows the same order as scope: exact subtype first, then the
union across the genotype's subtypes.

### Why the *dominant* residue and not every typical one

This is the most consequential decision in the design, and it is worth stating
plainly because getting it wrong is clinically dangerous.

At `NS3:80` in genotype 1a the empirical residues are **Q at 60.9% and K at
36.0%**. Q80K is a common natural polymorphism — *and* a genuine simeprevir
resistance substitution.

- Suppressing every residue the alignment calls "typical" would suppress K, and
  **Q80K would stop being reported in every genotype-1a sequence**.
- Suppressing only the *dominant* residue suppresses Q80Q and still calls Q80K.

The catalogue's own `display_structure` spells the 1a finding `K/Q80K` — it
lists K itself as a wild type — so deriving wild types from `display_structure`
produces exactly the dangerous outcome. This is why the wild-type source is
`phdr_alignment_typical_aa.csv` and the rule is "dominant residue", not
"any typical residue".

`NS3:80` is the standing control for this behaviour, covered by a test.

---

## 4. Where the data comes from

The pipeline reads **only** the generalized catalogue TSV. It never opens the
PHDR CSVs, and it never reads `alignment_name`.

```
phdr_alignment_typical_aa.csv  ─┐
var_almt_note.csv              ─┼─► BuildCatalogGenotypeColumns.py ─► catalog TSV ─► AnnotateMutations.py
phdr_alignment_ras.csv         ─┘        (HCV-specific)              (generic)        (generic)
```

**Why the split.** `AL_1a`, `AL_6xd` and the rest are PHDR/HCV artefacts. A
rabies or influenza catalogue will never have them. Anything genotype-related is
therefore resolved once, at catalogue build time, into two generic columns that
any virus can supply by any means it likes. No pipeline code depends on an
alignment column.

### Why two source tables are needed

| | `phdr_alignment_typical_aa.csv` | `var_almt_note.csv` |
|---|---|---|
| rows | 12,113 | 71,016 |
| alignments | 138 | 139 |
| gives | dominant residue per position | observed frequency per variation |
| floor | records residues at **≥10%** only | no floor |
| combinations | no — strictly per position | yes |

770 of the 807 catalogued (RAS, alignment) pairs sit **below 10% frequency**, so
`typical_aa` cannot describe them at all, and 395 of the 807 are combinations,
which it does not model. `var_almt_note` covers both gaps. Neither table alone is
sufficient.

Coverage: 138 alignments means **97% of observed sequences (266/273) get an
exact-subtype wild type** rather than a genotype-level approximation. The
exceptions are genotype 8 (there is no `AL_8` at any granularity), `7b`, `6xf`
and `6xh` — 7 sequences, which fall through to `wt_unknown` and are therefore
still reported.

### What was checked for loss

All **807 (RAS, alignment) pairs** in `phdr_alignment_ras.csv` survive into the
generalized catalogue, across all 23 alignments. Nothing was lost to
deduplication. `variation.csv`'s omissions are all justified — 20 of its 26
columns are empty, including the entire epitope block.

---

## 4b. Clinical trial evidence

Each catalogue row carries the registered trials supporting it **for that row's
genotype and drug**, in a `clinical_trials` column of semicolon-separated NCT
identifiers.

The chain PHDR ships:

```
resistance_finding.phdr_alignment_ras_drug_id    NS3:107I:AL_1a:grazoprevir
    -> resistance_finding.phdr_in_vivo_result_id
    -> result_trial.phdr_clinical_trial_id
    -> clinical_trial.nct_id                     NCT01717326
```

Note that the first key already contains the RAS, the genotype **and** the drug,
so the linkage is genotype-scoped at source. That is not a formality: trial
support really does vary by genotype.

| `NS5A:31M` vs daclatasvir | trials |
|---|---|
| genotype 1a | 1 |
| genotype 1b | **9** |
| genotype 2a | 5 |
| genotype 3a | 0 |

49 of the 64 (RAS, drug) pairs curated in more than one genotype have different
trial sets, so a row takes **its own** genotype's trials rather than the union
across the signature's scope. Unioning would attach 1b's evidence to a 1a call.

Only in-vivo findings carry a trial — an in-vitro EC50 has none, and rows without
a drug get nothing. 746 of 1,869 rows are populated, drawing on 97 distinct
trials out of the registry's 102.

**A note on what was nearly built instead.** Before the trial tables arrived
there was no registry identifier anywhere in the PHDR export — `publication.id`
is a PubMed ID and `url` is a DOI — so trial status could only have been guessed
from whether a publication's title contained the word "trial". That heuristic was
written, then deleted when the real data landed. It is worth recording only
because a heuristic that looks like a curated fact is worse than no field at all.

**A related data-shape fact.** The column named `pubmed_id` is really
"publication reference": 8 of its 128 values are conference abstracts with no
PMID (`AASLD_2015_Abs_718`, `EASL_2017_Abs_THU-257`). Anything parsing that
column as an integer will break on them.

---

## 5. What ends up in the database

| table | contents |
|---|---|
| `mutation_catalog` | the catalogue as loaded |
| `sequence_relevant_mutation_summary` | per sequence, the mutations present and their count |
| `completed_signatures_only` | per sequence, the signatures fully satisfied |
| `sequence_mutation_calls` | **every evaluated call**, emitted or suppressed, with the reason |

`sequence_mutation_calls` exists because a suppressed call and a call that was
never evaluated are indistinguishable from the summary tables alone. Each row
carries `call_status` (`emitted` / `suppressed_out_of_scope` /
`suppressed_wild_type`), `scope_tier`, `residue_status` (`change` / `anchor` /
`wt_unknown`), the wild-type residues consulted, the observed residue, and the
sequence's genotype and subtype.

If a call you expected is missing, that table says why.

---

## 6. Measured effect

On `test_out/HCV_OM_test/HCV_OM_test.db` — 338 sequences, 238 annotated, 5,817
calls before the change:

| outcome | calls | share |
|---|---|---|
| dropped — wild type for its own genotype | 2,819 | 48.5% |
| dropped — out of genotype scope | 2,459 | 42.3% |
| kept — genuine change | 503 | 8.6% |
| kept — wild type unknown | 36 | 0.6% |

**H77 goes from 22 calls to 0**, with no genuine change lost.

Controls, both covered by tests:

- H77 is **not** flagged at `NS3:80` — it carries Q, the 1a wild type.
- A 1a sequence carrying K at `NS3:80` **is** still flagged — Q80K survives.

---

## 7. Open decisions and known limits

**Deletions are detected but rarely in scope.** With `-` standardised, 214
deletion matches are now evaluated where the old code found none. In this
dataset all 214 are out of genotype scope — NS5A 29/30/32del are scoped to 1a/1b
and the carriers are genotypes 2/3/5/6/8 — so none are emitted. The mechanism
works; the dataset simply has no in-scope carrier.

**A gap is only a deletion when it is internal.** Partial GenBank records padded
out to full genome width produce all-gap codons at NS5A 29/30/32 in 35
sequences. Those are missing data, not deletions, and are read as `X`. Only a gap
strictly inside the sequence's covered span counts. Without this the change would
have flagged the worst-covered records hardest.

**Stop codons are standardised to `*` but nothing uses them yet.** The codon
table previously emitted `_` while the catalogue grammar admits `*`. No
catalogue row currently uses `*`, so this was a latent trap rather than active
data loss — it is closed before the next catalogue revision can fall into it.

**Sequences with no genotype are out of scope** for any scoped entry; unscoped
entries still apply. Recorded, not silent: such calls land in
`sequence_mutation_calls` as `suppressed_out_of_scope` and `diagnostics` counts
`sequences_without_genotype`.

**Scope is per signature, not per mutation.** A signature scored in both `AL_1a`
and `AL_1b` carries both codes on every one of its rows. Per-mutation scope would
be a one-line change to `build_signature_genotype_scope` if that turns out to be
wanted.

**The frequency field is advisory.** `wild_type_residues` may carry a trailing
frequency, and nothing in the pipeline currently gates on it. It is recorded so a
reader can see that the "wild type" at `NS3:80` is only a 61% majority. If a
frequency-based rule is added later it must be gated on the denominator: 84 of
139 alignments have a maximum N below 10, and a percentage computed from 7
sequences means very little.

**Mutant frequency is not currently stored.** `var_almt_note.csv` gives the
observed frequency of each *variant* within each genotype — Q80K is present in
36.09% of curated genotype-1a sequences, for instance. That is a different
quantity from the wild-type frequency above and it is not carried in the
catalogue, because putting it in `relevant_genotypes` would have meant two
different numbers in identically-formatted adjacent columns.

It would be useful for one specific thing: catching common-but-not-majority
polymorphisms that the wild-type rule cannot suppress by construction. 27
catalogued pairs are ≥10% frequent within their own genotype and still emit —
`AL_5a NS3:168E` at 51.7%, `AL_1b NS5A:37L` at 44.5%, `AL_2b NS5A:31M` at 38.2%.
Those are baseline polymorphisms being reported as resistance events. If that
matters, the fix is a third optional column rather than overloading either
existing one.
