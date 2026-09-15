# Lineage growth rates

`scripts/lineage_growth_rates.py` is a bolt-on. It reads a finished database
and writes tables and a report; it touches the database only with `--write_db`.
It is not part of the Nextflow workflow and nothing in the pipeline calls it.

```bash
python scripts/lineage_growth_rates.py --db DB --list-proteins

python scripts/lineage_growth_rates.py \
    --db test_out/HCV_XML_full_plus_update/HCV_full_right_tree.db \
    --protein NS5A --out dev/growth/NS5A --fit-window 25
```

Everything runs single-threaded and read-only. A 137,000-tip UShER tree with a
60,000-sequence cohort takes about half a minute and 700 MB.

---

## What it measures, and why there are so many numbers

A lineage's growth rate is not one quantity with one estimator. The estimators
below make different assumptions, and they are all reported because the places
they disagree are the places the answer is not safe.

### Without a tree

You do not need a tree to estimate a change in allele or lineage frequency.
You need dates.

| Estimator | What it gives | Where it fails |
|---|---|---|
| Logistic regression of carrier frequency on time | The growth-rate advantage per year. Under a two-type exponential model the slope *is* the difference in growth rates, which is why it is the standard measure of a variant's fitness advantage | Assumes one constant advantage over the window |
| The same fit on binned counts, quasi-binomial | The same slope with honest intervals | Needs enough sequences per bin |
| Multinomial logistic regression across a site's residues | Mutually consistent advantages for every residue at once | Same assumptions, plus a pivot to read against |
| Poisson log-linear fit to counts | Absolute growth, and (with a sampling offset) relative growth | Counts follow sequencing effort; the two are reported separately for that reason |
| Frequency Increment Test | Whether increments exceed drift, without assuming logistic | Needs three or more polymorphic time points |
| Mann-Kendall with Sen's slope | A monotone trend, non-parametrically | Says nothing about rate |
| Early-versus-late Fisher exact | An odds ratio and an exact p-value | Throws away the trajectory |

The logistic slope is the headline. Frequencies are immune to sequencing effort
in a way counts are not: sequence ten times as much in 2020 and the counts
change, the frequency does not.

### With a tree

The tree adds ancestry, which is the one thing dates cannot supply.

| Estimator | What it gives |
|---|---|
| Weighted parsimony reconstruction | How many times an allele arose **independently**, as a range over all most-parsimonious histories |
| Per-origin logistic fits | Whether a frequency rise is one lucky lineage or the mutation repeatedly |
| Lineages-through-time slope | Net diversification rate, `lambda - mu` |
| Magallon-Sanderson method of moments | Net rate from a clade's size and age alone, over a range of assumed extinction fractions |
| Maximum-likelihood Yule rate | Pure-birth rate; an upper bound on the net rate |
| Pybus & Harvey gamma | Whether branching decelerated - i.e. whether one rate is a fair summary |
| Coalescent skyline | Ne through time |
| Exponential-growth coalescent (ML, profiled) | A growth rate read off the genealogy's *shape*, independent of how common the lineage is |
| Local branching index (Neher et al. 2014) | The tree-shape fitness proxy: where branching is dense right now |
| Expansion score | Recent sequences over the previous window, assumption-free |
| Whole-tree clade scan | Every clade's frequency growth, FDR corrected |

---

## The three things most likely to mislead you

### 1. Node dating

Every age-based estimator needs node *times*, and an UShER newick has none: its
branch lengths are mutation counts, and in a database built from a de-duplicated
alignment most of them are zero.

* `--dating mrca-bound` (default) dates each internal node to the earliest
  sample below it. That is a **bound**, not an estimate. Clades look younger
  than they are, so Magallon-Sanderson, Yule and the coalescent rates come out
  too fast, and the lineages-through-time curve starts at its maximum and only
  falls - which makes the LTT slope and gamma meaningless (gamma near -50 is
  the dating, not a biological slowdown). The report says so wherever it applies.
* `--dating root-to-tip` is the TempEst regression. It refuses when there are
  fewer than 20 dated tips or the clock explains less than 10% of the variance,
  and falls back with a message saying why. On HCV across all genotypes it
  refuses, correctly - there is no single clock across the genotypes.
* `--dating branch-lengths` takes the lengths to be years already, for a tree
  out of a dating program.

**The frequency-based estimators do not use node times at all.** They are
unaffected by any of this, which is why they are the headline.

### 2. The denominator

`--allele-source catalog` reads `sequence_mutations`, which records the residues
the annotator *called*, not every sequence it looked at. The denominator is
therefore the set of sequences that produced any call for that protein, and the
report states the number. On the shipped HCV build that is 259 sequences out of
135,544 - frequencies are relative to those 259, not to the database.

`--allele-source features` and `alignment` genotype every sequence in the cohort
from the alignment, so the denominator is the cohort. `features` is exact: the
protein is an annotated product and its coordinates come from the database.
`alignment` needs `--protein-start` in master coordinates, because a
sub-peptide of a polyprotein (HCV NS5A) has no feature row to take them from.

`--calibrate-start` infers that start by testing every frame against the
catalogue calls and refuses below 90% agreement. On HCV NS5A it finds master
coordinate 6389 at 85.8% agreement against a 46.6% runner-up and refuses. If
you override it with `--protein-start 6389` the reading-frame check then reports
that 75% of sequences carry a whole-codon deletion, which is what a wrong frame
looks like. Both guards exist because a wrong frame does not fail - it produces
hundreds of alleles with tiny p-values describing a reading frame.

### 3. Recency and the clade scan

A clade first sampled at the very end of the window always comes out near the
top of the scan with a tiny q-value, because over that window it did go from
absent to present. Check `first_date` and `n_time_bins_occupied`, and set
`--fit-window` so the question is asked over a period the cohort actually
covers. Without it, a database spanning 1905-2024 spends most of its time axis
on decades holding a handful of sequences.

---

## Independent origins

This is the part a tree is genuinely needed for, and it is what the frequency
fits cannot do. A mutation rising in frequency might be spreading because it
helps, or it might have arisen once inside a lineage spreading for its own
reasons. Splitting carriers by origin separates them.

* `ancestral_state` says whether the residue was reconstructed at the root. A
  wild-type residue is `present`, and its informative count is the **losses** -
  how many times something else arose in its place. Forcing the root to
  "absent", which an earlier version did, made every wild-type residue look like
  a mutation that arose independently in each sequence carrying it: on a
  polytomous tree, once per tip.
* `origins_min` and `origins_max` bracket the count over **all**
  most-parsimonious histories, including both root states when they tie.
  Parsimony often does not identify a unique history, and a single number is a
  statement about the tie-breaking rule.
* `origins_per_carrier` near 1 means every carrier is its own origin: the
  residue keeps arising and does not then spread - within-host selection, not
  transmission fitness. Near 0 means one origin did all the spreading, and the
  frequency rise belongs to that lineage rather than to the residue.

On the shipped HCV build, NS5A:31M has 73 carriers in 73 independent origins and
none of them expands: strongly recurrent, not transmitting.

---

## Genotypes

`--genotype` needs a genotype column. Databases built before genotype columns
existed have none, and the error says so. `--genotype-reflist <curated list>`
derives them from the tree using the pipeline's own `clade_from_tree`
assignment, so a genotype derived here is the genotype the pipeline would have
stored.

That assignment refuses to call a genotype when the enclosing clade spans more
than half the tree. On a heavily polytomous UShER tree that is the usual
outcome - most queries hang off the root polytomy and have no local
neighbourhood - and the run warns when most sequences come back unlabelled. The
IQ-TREE backbone is properly resolved and works much better for this, at the
cost of only holding cluster representatives.

---

## Statistics

There is no scipy or statsmodels in the `vgtk` environment, so `growth_stats.py`
implements what it needs: the normal, chi-squared and t tails, the normal
quantile, Fisher's exact test, Mann-Kendall, Benjamini-Hochberg, and a GLM
fitted by IRLS with binomial and Poisson families. Adding scipy to
`environment.yml` for a bolt-on would change the environment every pipeline run
has to solve.

Two choices are worth knowing about:

* **Firth penalisation** is applied automatically when either outcome class has
  fewer than 25 observations, or when the plain fit separates. Without it, an
  allele whose carriers are all in the last two years has an infinite
  maximum-likelihood slope, an infinite standard error, and a Wald p-value of
  1.0 for the strongest signal in the dataset. `logistic_method` records which
  was used. Significance comes from a likelihood-ratio test, not from Wald, so
  on a penalised fit the reported q-value and the Wald confidence interval can
  disagree - believe the likelihood ratio.
* **A weak Gaussian prior** (sd 5 per year) on the multinomial slopes. A residue
  seen twice, both times late, is perfectly separated and the unpenalised fit
  returns 1e31 with a p-value of exactly zero. A prior far wider than any real
  growth advantage leaves determined estimates alone and is decisive only where
  the data determine nothing.

---

## Output

One TSV per table plus `report.md`, and with `--plots` two PNGs.

| File | Contents |
|---|---|
| `alleles.tsv` | Every tree-free estimator per allele, FDR corrected |
| `trajectories.tsv` | Binned frequencies with Wilson intervals |
| `multinomial.tsv` | Jointly-fitted advantages per site |
| `origins.tsv` | Independent origins, losses, ancestral state per allele |
| `origin_lineages.tsv` | Each origin separately, with its own growth rate |
| `phylodynamics.tsv` | LTT, Magallon-Sanderson, Yule, gamma, coalescent per origin |
| `clades.tsv` | The whole-tree clade scan |

`--write_db` appends the same tables into the database under a `run_id`, with a
`lineage_growth_runs` row recording the parameters.

## Tests

`tests/unit/test_growth_stats.py` checks the distributions against published
table values and recovers known coefficients from simulated data.
`tests/unit/test_tree_growth_metrics.py` recovers a planted clock rate, a
simulated Yule rate and a simulated coalescent growth rate, and pins the
parsimony behaviour. `tests/unit/test_lineage_growth_rates.py` builds a small
database with a planted logistic sweep and runs the script end to end.

---

## Segmented builds (influenza)

A segmented database breaks four assumptions that a single-segment build hides.
All four were found by running the tool against `test_out/IAV_DB/iav-db.db` and
each is now handled, but they are worth knowing because they all fail *quietly*.

### One master per segment

An influenza build carries eight masters, one per segment. Resolving HA
coordinates against segment 2's master yields a residue for every sequence and
every one is wrong. `--segment` is therefore **required** whenever the build has
more than one master; the run refuses rather than guessing.

### One tree per segment

The trees table holds a tree per segment. Selecting "the longest UShER newick"
gave a neuraminidase cohort the haemagglutinin tree. Here that failed loudly
because the two share no tips — but on a build whose segments share accession
names it would have fitted the wrong tree silently. Tree choice now filters on
the segment first, and falls back only with a warning.

### `cds_start` is an alignment column, not an ungapped position

`features.cds_start` / `cds_end` in these builds are **alignment columns**. The
HA master AB573800 is 1782 columns wide with 132 gaps, so its ungapped length is
1650 — exactly its `aln_end` — while `cds_end` is 1782, exactly the width.
Reading those coordinates as ungapped positions puts the frame 18 codons out and
runs 44 codons off the end of the sequence.

On HCV the two readings are identical, because that master carries no gaps,
which is why this went unnoticed. `--coord-space auto` (the default) translates
the master's own CDS both ways and takes the reading with no internal stop
codons. On HA that is decisive: `column` gives 0 stops and 0 unresolved codons,
`ungapped` gives 44 unresolved. The reading used is printed and recorded in the
report.

### Padding is not deletion

A sequence that stops short of the CDS end is padded with gaps. Translating
those as whole-codon deletions gives every partial sequence "deletion" alleles
across the tail — and because sequencing completeness improved over time, those
correlate with date and come out as **the fastest-growing alleles in the run**.
On the HA test build the top four alleles were exactly that: `573-`, `577-`,
`574-`, `570-`.

Codons are now translated with the covered-span rule `AnnotateMutations`
already used: an all-gap codon is a deletion only when it sits strictly inside
the sequence's own first and last non-gap column, and padding outside that span
is unknown (`X`), not a deletion.

### Indel columns are not a broken frame

The frame check originally treated a high whole-codon-deletion rate as evidence
of a wrong start. Across subtypes that is wrong: HA subtypes differ by real
indels against any single master, so an H3 sequence read against an H1 master
genuinely lacks codons the master has — 74.6% deletion at one residue on the
test build, with a perfectly correct frame.

**Internal stop codons diagnose the frame; deletions do not.** High-deletion
sites are now reported as indel columns instead, with the caveat that alleles
there describe the alignment rather than a substitution.

### Site numbering

Positions are codon numbers counting from the CDS start in the master's own
frame. For influenza HA that is **neither H1 nor H3 numbering** — the master
frame is shared across subtypes that differ by indels — so numbers are
comparable within a run and not to the literature.

---

## Cost on a large database

Two things dominate, and both are handled:

* **Reading alignments.** The segment filter is pushed into SQL. Without it,
  every HA query drags the other seven segments off disk to discard them: eight
  times the I/O against a 30 GB file for the same answer.
* **Fitting.** Sequences sharing a sampling date carry identical information,
  and the binomial logit likelihood for grouped counts is the *same function* as
  the one for the individual rows. Collapsing onto distinct dates is therefore
  exact, not an approximation, and there is a standing test asserting the two
  agree to 1e-8. It is the difference between minutes and hours: 81,000
  sequences carry only a few thousand distinct dates.

Note that these databases carry **no indexes**, so every query is a full scan.
A cohort load on the 30 GB influenza build takes about 90 seconds on its own.

### Spliced products are refused

`matrix protein 2` and `nuclear export protein` are each annotated as **two**
feature rows, one per exon. A single `cds_start..cds_end` span cannot represent
a spliced CDS, and taking the first row translates 26 nucleotides of M2 exon 1
and reports it as the protein. Any product with more than one feature row on the
master is now refused, with the spans listed and a pointer to
`--allele-source alignment --protein-start --positions` for a specific exon.

### What the battery found on segments 7 and 8

Running every segment of `test_out/IAV_DB/iav-db.db` turned up something worth
knowing: **M1 and NS1 do not translate cleanly under either coordinate
reading.** M1 carries a stop codon at residue 38 in 100% of sequences, NS1 at
residue 176. Both readings give 10 and 19 internal stops respectively on the
master itself, so this is not a choice-of-convention problem.

These are the two spliced segments (M1/M2 on segment 7, NS1/NEP on segment 8),
and the stored coordinates do not yield a clean contiguous CDS on those masters.
The frame check catches it and both reports say so at the top; the allele tables
for those runs should not be used. Whether the upstream annotation can express
these correctly is a separate question — the point here is that the tool refuses
to present the numbers as sound rather than quietly ranking a reading frame.

HA, NA, NP and PB2 all translate cleanly, and HA's high-deletion site is
correctly reported as an indel column rather than a frame failure.

---

## Using a topology that is not in the database yet

`--tree-file <newick>` reads the topology from a file instead of the `trees`
table. This is the normal case while a rebuild is in flight: the database still
holds the old one, and growth rates computed against it describe a topology that
has been superseded.

On the influenza build the difference is not marginal:

| source | tips | in segment 4 | dated |
|---|---|---|---|
| `iqtree_refset_4...`, stored in the database | 9,470 | 9,470 | 4,853 |
| `graft/starter.rooted.nwk`, the pre-placement backbone | 207,040 | 207,026 | 84,089 |
| **newest finished chunk**, `chunk_0022/usher.pb` | **317,040** | **317,018** | **131,701** |

The stored one holds only cluster representatives, so a growth rate measured on
it is partly measuring the clustering.

**Use the newest finished chunk, not the starter backbone.** A chunk rebuild
places samples cumulatively, so each finished chunk is a strict superset of the
one before it and of the backbone it grew from: here 110,000 tips the starter
does not contain, and not one that it lacks. Extract a newick from the protobuf
with

```bash
matUtils extract -i treepatch/work/seg4/usher/chunk_0022/usher.pb \
    -T 4 --write-tree newest_chunk.nwk
```

Unlike the HCV UShER tree, whose branch lengths were almost all zero, these
carry real ones (55% non-zero), so `--dating root-to-tip` is worth attempting
instead of falling straight back to the mrca-bound.

A missing `--tree-file` is refused rather than quietly falling back to the
database, because that would answer a different question from the one asked.

**Caveat worth recording in any write-up:** the newest finished chunk is a
*partial* rebuild — 22 of 33 chunks placed at the time of writing, so roughly a
third of the samples are not on it yet. It is much the best topology available
and improves with every chunk, but it is not the finished tree. Clade-level
results — independent origins and the clade scan — should be re-run when the
rebuild completes. The frequency-based rates never touch the topology and will
not change.

`dev/growth/track_rebuild.sh` reports how far the rebuild has got and projects a
finish. It fits the *trend* rather than averaging a flat rate, because the cost
per chunk climbs as the tree grows — on segment 4 it rose from 5.9 h for the
first chunk to 14.4 h for the twenty-second, about +0.45 h per chunk with an R²
of 0.94. Averaging the recent rate instead puts the finish two days early. Both
numbers are printed. It is read-only and safe to run at any time.
