# Divergent run mode

Opt-in path for virus families too divergent for the standard reference-guided
pipeline. Testbed: *Cystoviridae* (NCBI TaxID **10877** — not 10880, which NCBI
now maps to the order Reovirales).

```
# Cystoviridae: a single-segment family, 33 records
nextflow run vgtk-init.nf -profile divergent_test

# Hantaviridae: the per-segment path, 45 records, three backbones
nextflow run vgtk-init.nf -profile hanta_divergent_test
```

or on your own data:

```
nextflow run vgtk-init.nf \
  --run_divergent_mode true \
  --divergent_genbank /path/to/records.gb \
  --publish_dir results/divergent
```

## When you need it

The standard path aligns every sequence to one master accession with
`nextalign`. That works while queries share exact k-mers with the master. On a
family where distinct species are 30–50% divergent it does not degrade — it
**stops the run**: references go through `nextalign_master`, which raises
`RuntimeError` once the relaxation cascade is exhausted
(`scripts/NextalignAlignment.py`).

Symptoms that point here:

- `CRITICAL ERROR: Fundamental alignment failure for reference sequence`
- references aligning at 60–80% identity with double-digit gap fractions
- UShER placements costing hundreds or thousands of substitutions

## Check your family is suitable before you start

Seed collapse is **not** the test. Across nine families screened, the median
inter-species seed match rate was `0.0000` for every family and every segment —
the standard path fails on all of them. What decides whether this mode can help
is whether **every segment carries a gene recoverable across the divergent
taxa**. Without that, there is no per-segment backbone to build and the mode
collapses to a single-segment analysis.

```bash
python dev/divergent_mode/screen_candidate_families.py \
  --families "MyFamily:<taxid>" \
  --cache-dir dev/divergent_mode/screen_cache \
  -o /tmp/screen.tsv --api-key "$(cat dev/ncbi_key)"
```

Then check the protein-level divergence is inside the usable band — nucleotides
too far apart but proteins still alignable. Below roughly 0.20 median
amino-acid identity, codon-awareness does not rescue anything:

```bash
python dev/divergent_mode/quantify_divergence.py \
  --genbank dev/divergent_mode/screen_cache/txid<taxid>.gb \
  --family MyFamily -o /tmp/div.tsv
```

### Where this mode applies

| group | per-segment core gene? | examples |
|---|---|---|
| Negative-sense segmented RNA | **yes, on every segment** | Hantaviridae, Peribunyaviridae, Nairoviridae, Tospoviridae, Arenaviridae |
| Segmented dsRNA | **RdRp segment only** | Cystoviridae, Chrysoviridae, Partitiviridae, Spinareoviridae |

On a segmented dsRNA family, expect a single-segment (RdRp) backbone and plan
for it — that is the correct answer for those families, not a failure. On
Cystoviridae, L_P2 and L_P4 cover 0.95 and 1.00 of the divergent taxa while M
and S peak at 0.45 and 0.60.

## What it does differently

| Problem | Standard path | Divergent mode |
|---|---|---|
| Seeding | nucleotide k-mers vs one master | protein homology (tblastn) per gene |
| Frame | gap penalties discourage sub-codon gaps | gaps are codon-wide **by construction** |
| Deep topology | parsimony onto one reference | ML backbone over cluster representatives |
| Dense clusters | same single parsimony run | UShER per cluster, grafted onto the backbone |

### Why there are no frameshifts

CDS are translated, aligned as protein with MAFFT, then the nucleotide
alignment is rebuilt one codon per residue and three gap characters per protein
gap. A gap in the output can only come from a gap in the protein alignment,
which is one residue wide by definition. A 1 or 2 nt gap is unreachable, not
merely penalised. `DM_CODON_ALIGN` re-checks this on its own output and aborts
if it ever fails.

## Gene selection is scoped to representatives, and that matters

A gene only helps the deep topology if it is recoverable across the *divergent*
taxa. Scoring coverage over all taxa gives the wrong answer whenever one dense
cluster dominates the dataset: all 48 Cystoviridae phi6 isolates carry every
gene, so the unscoped ranking picks five genes, three of which collapse outside
phi6.

So the pipeline partitions each segment on its own most complete gene, then
scores coverage over dense clusters rather than isolates. On Cystoviridae that
picks **L_P2 (RdRp) + L_P4 (packaging NTPase)** — the only two genes above 90%
coverage among the divergent taxa. M and S have none above 60%; P3, the
host-attachment spike, is at 10%.

If a segment reports `below_min_coverage`, that is the pipeline telling you
either that the segment has no family-wide gene, or that your reference set is
too narrow for it. Neither is a bug. Read that segment's
`backbone_gene_selection.tsv`, then widen the references or lower
`--divergent_min_gene_coverage` knowingly. The run only aborts if *every*
segment comes up empty.

## The cut-off

`--divergent_cutoff auto` (the default) finds the widest empty band in the
pairwise distance distribution and puts the cut-off at its midpoint. On
Cystoviridae 1,145 pairs sit at or below 6.0% divergence, 1,066 at or above
21.2%, and none in between, so the choice is unambiguous (it picks 13.6%).

A band qualifies on two tests: it must be at least 0.01 wide in absolute
terms, and it must span at least `--divergent_auto_min_band_frac` (default
0.15) of the observed distance range. The second is the guard against a
continuum, where no single gap ever holds a real share of the range —
measured, it separates cleanly: continuum distributions score 0.04–0.09,
genuinely clustered ones 0.25–0.61.

The band is judged against the **spread of the data, not against its own lower
edge**. That distinction is load-bearing. The original rule asked for
`gap >= 2 x band_lo`, i.e. `band_hi >= 3 x band_lo`; since candidate bands are
capped at `--divergent_auto_ceiling` (0.35), no band whose lower edge exceeded
0.35/3 = 0.117 could ever qualify, however empty it was. The rule was
unsatisfiable above that scale rather than merely strict, because it assumed
dense clusters are tight. Puumala orthohantavirus is not: its isolates sit
8–16% apart on segment L, with a clean empty band from 15.9% to 23.4% before
any other species. That band was rejected and segment L fell back to 1.5%,
which shattered all 15 taxa into singletons and skipped the graft step
entirely.

It still falls back to 1.5% when no band qualifies, and `partition_cutoff.json`
now records `reason`, `band_gap`, `band_frac` and `min_band_frac` so a fallback
can be read rather than guessed at. **A fallback still matters**: it usually
means the dataset is a divergence continuum rather than clusters, and the
backbone/subtree split is then arbitrary rather than natural. Check the
recorded reason before trusting the tree.

Pass a number (e.g. `--divergent_cutoff 0.015`) to force a specific cut-off.

## Branch-length units

UShER reports branch lengths as a **count** of substitutions; the ML backbone
reports **substitutions per site**. The graft step divides subtree lengths by
the alignment width so the unified tree carries one unit. Without that, a
within-cluster branch of 119 substitutions dwarfs a between-genus branch of 1.2
subs/site and every downstream distance reads the mixture as comparable.
`--no-rescale-subtree-branches` turns the conversion off; there is rarely a
reason to.

## One backbone per segment

Everything up to the codon alignment is family-wide and runs once. From gene
selection onwards the run fans out: **each segment gets its own seed gene, its
own distance partition, its own gene selection, its own supermatrix, its own
ML backbone and its own grafted tree**, published under
`divergent_mode/by_segment/<segment>/`.

A segment is its own analysis, not a column of a shared one. It has its own
accessions, its own dense clusters, and — in a family that reassorts, which
is most segmented families — its own topology. Concatenating L, M and S into
one supermatrix would impose a single tree on three histories, and on this
dataset there is no isolate key linking an L record to its M and S anyway.

A segment carrying no gene recoverable across the divergent taxa drops out of
the run. That is the right answer on a segmented dsRNA family. It is **not**
silent: every requested segment gets a row in
`divergent_mode/backbone_summary.tsv` saying whether a backbone was built and,
if not, which gene came closest and why it was not enough. The run fails
outright only if no segment produces a backbone.

```
segment  backbone_built  status              genes        best_gene  best_coverage
L        Y               ok                  L_P2,L_P4    L_P2       1.0000
M        N               below_min_coverage               M_P3       0.7143
S        N               below_min_coverage               S_P5       0.8571
```

That is the Cystoviridae fixture: one backbone, from L, with M and S named
and explained rather than absent. On Hantaviridae all three segments build a
backbone *and* graft a dense cluster — L 4 Puumala isolates, M 2, S 3 — so
every per-segment code path is exercised on all three.

### How coverage is counted, and the trap in it

Coverage is the fraction of a segment's **dense clusters** — not isolates —
in which a gene is recoverable. Clusters, because 48 near-identical phi6
isolates are one piece of evidence about the family, not 48.

The denominator is every genome kept for that segment, from
`segments/segment_manifest.tsv`. It is deliberately *not* the membership of
the partition. A partition is built from one seed gene's alignment, so a
genome missing that gene is missing from the partition too; scoring against
partition membership measures each gene over the taxa that already carry the
seed gene, and the seed gene then scores 1.00 by construction. That mistake
promoted Cystoviridae's M_P3 — the host-attachment spike, present in 10% of
the divergent taxa — to a backbone gene at full coverage.

Scoping one segment's selection against another segment's partition is worse
and harder to see: the two share no accessions, so every gene scores zero and
the segment reads as "no usable gene" when nothing about it was measured.
`SelectBackboneGenes.py` refuses both — one segment per call, and the
partition must contain that segment's accessions.

## What it does *not* do

**It stops before the database build.** The legacy path resolves every
coordinate against one master accession, and the premise here is that no such
master exists. Wiring these steps into `CreateSqliteDB` would change the
coordinate contract that `AnnotateMutations` and the insertions table depend
on, so that is separate work with its own schema review.

Outputs are published under `<publish_dir>/divergent_mode/`. The first four
are family-wide; everything else is per segment.

| Path | Contents |
|---|---|
| `segments/` | per-segment FASTA + manifest |
| `orthogroups/` | protein orthologue groups (the tblastn query set) |
| `cds/` | extracted CDS + `cds_extraction_report.tsv` |
| `codon_aln/` | codon alignments + frameshift report |
| `backbone_summary.tsv` | one row per segment: backbone built, or why not |
| `by_segment/<seg>/<seg>.backbone_gene_selection.tsv` | the gene ranking, written even when the segment is skipped |
| `by_segment/<seg>/seed/partition/` | the partition that scopes this segment's selection |
| `by_segment/<seg>/backbone/` | supermatrix, partitions, coordinate map |
| `by_segment/<seg>/final/partition/` | distance matrix, cluster assignment, chosen cut-off |
| `by_segment/<seg>/backbone_ml/` | IQ-TREE tree, ancestral sequences, parent–child pairs |
| `by_segment/<seg>/grafted/` | unified tree + graft report |

## Non-regression

With `--run_divergent_mode false` (the default) the branch returns before any
legacy process is instantiated. Verified by running the `test` profile from an
unmodified checkout and from the modified tree with the same pinned
`--xml_dir`: every data artefact and all 17 output DB tables are content-
identical. The few tables that differ also differ between two runs of the
*unmodified* baseline — the pipeline is not byte-reproducible run to run
(`GENBANK_PARSER` record order varies and `IQ_TREE` takes no `-seed`), which
predates this work.

Tests: `tests/unit/test_divergent_mode.py` (49 tests), plus the
`divergent_test` (Cystoviridae, one backbone) and `hanta_divergent_test`
(Hantaviridae, three backbones, three grafted clusters) Nextflow profiles.
