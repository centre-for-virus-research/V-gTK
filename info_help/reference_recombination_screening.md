# Screening the reference set for recombination

Recombinant sequences break the assumption every genotype call rests on: that a
genome descends from one parent. A recombinant hiding in the **reference set** is
the worst case, because every query that matches it inherits the wrong genotype.
This screen finds them with UShER/RIPPLES.

Script: [`scripts/ScreenReferenceRecombination.py`](../scripts/ScreenReferenceRecombination.py)
Tests: `tests/unit/test_reference_recombination.py`

---

## Why the reference set and not everything

RIPPLES detects recombination on a mutation-annotated tree by trying to explain a
branch's mutations with **two** parents instead of one. Its cost per branch grows
quadratically with the number of mutations on that branch, because it enumerates
breakpoint *pairs*.

For SARS-CoV-2, which RIPPLES was built for, branches carry a handful of
mutations. For a diverse virus they do not. Measured on this repository's HCV
data, 2026-08-27:

| tree | long branches | breakpoint pairs on the first branch | implied work |
|---|---|---|---|
| 33k-sample query MAT | 1,790 | **19,341,090** | ~3.5 × 10¹⁰ placements |
| 238-sequence reference MAT | 156 | 2,864,421 | ~3 min/branch |

The query tree is simply not screenable this way. The reference set is — and it
is the case that matters most, because it is the foundation. It is also small and
changes rarely, so a slow screen is acceptable there.

**Threads do not rescue it.** Going from 8 to 48 threads made no measurable
difference to the first branch: RIPPLES does not parallelise well *within* a
branch. The way to use a big machine is several processes over disjoint branch
ranges — see [Chunking](#chunking).

**For query sequences at scale, use a different method.** A genotype-support
matrix (precompute `SUP[position, base, genotype]` once from the reference panel;
scoring a query is then a single gather, independent of panel size) runs at
~700 sequences/second — about 24 minutes for a million, single-threaded. That is
a separate tool and is not what this document describes.

---

## Thread ceiling

`MAX_THREADS = 16`, enforced in code, not left to the caller.

This runs on a shared machine that reports 224 cores. Both UShER and RIPPLES
default to taking every core they can see, which starves other users. Any
`--threads` above 16 is **clamped, and the clamp is announced** so a
silently-reduced run is never mistaken for the one you asked for.

If you need more throughput, run more chunks — not more threads.

---

## Running it

```bash
python scripts/ScreenReferenceRecombination.py \
    --db test_out/HCV_full_15_jul26/HCV_full_aug/HCV_full.db \
    --outdir reference_recombination \
    --threads 8 \
    --write_db
```

That does the whole chain: pulls the reference set out of the database, writes a
FASTA, converts it with `faToVcf`, builds a MAT with `usher`, runs `ripples`, and
stores what it finds.

Requires `faToVcf`, `usher` and `ripples` on `PATH` (all ship with the UShER
suite: `conda install -c bioconda usher`). A missing binary is reported by name
with the install command rather than as a bare `FileNotFoundError`.

### Inputs

| option | default | effect |
|---|---|---|
| `--db` | — | Database holding the reference set. Everything is derived from it. |
| `--alignment` | — | Reference alignment FASTA instead of `--db`. Requires `--tree`. |
| `--tree` | derived | Starting tree. Without it, the best tree in the `trees` table is pruned to the reference set. |
| `--outdir` | `reference_recombination` | Working directory; all intermediates are kept. |
| `--skip_build` | off | Reuse an existing `reference.pb` — use when re-running with different RIPPLES parameters. |

Two details that are not cosmetic:

- **The master is written first in the FASTA.** `faToVcf` treats the first record
  as the VCF reference, so the order fixes the coordinate frame. Writing anything
  else first would express every variant against the wrong reference.
- **The alignment must be rectangular.** A ragged one is refused with an error
  rather than passed to `faToVcf`, which would produce silently wrong
  coordinates.

Where no `--tree` is given, the tree is chosen by **coverage first**: whichever
tree in the `trees` table contains most of the reference set wins, because a
missing reference is simply absent from the screen. Among trees with equal
coverage an IQ-TREE beats an UShER tree. On HCV the UShER tree usually wins — it
holds all 238 references where the IQ-TREE, built from cluster representatives,
holds 219.

Only the topology is used: branch lengths and support values are stripped,
because UShER re-derives lengths from the VCF.

**Prefer an explicitly curated tree where you have one.** For HCV that is the
ICTV reference tree:

```bash
--tree dev/HCV_extra/HCV_ref_midpointrooted.treefile
```

It covers all 238 references and its topology has been checked against ICTV
(see [`hcv_tree_validation_ictv.md`](hcv_tree_validation_ictv.md)), which neither
database tree has.

### RIPPLES parameters

| option | default | RIPPLES flag | notes |
|---|---|---|---|
| `--branch_length` | 3 | `-l` | Minimum mutations on a branch before it is tested. **Raise this on a diverse virus** — it is the main lever on runtime. |
| `--num_descendants` | **3** | `-n` | RIPPLES' own default is 10. A curated reference set has a handful of sequences per subtype, so 10 skips most of the tree. |
| `--min_range` | 1000 | `-r` | Minimum genomic span of the mutations on a candidate branch. |
| `--max_range` | 10000000 | `-R` | Maximum span. |
| `--parsimony_improvement` | 3 | `-p` | How much better the two-parent explanation must be. Raise to reduce false positives. |
| `--timeout` | none | — | Seconds before RIPPLES is stopped. **Partial results are kept**, not discarded. |

Every parameter is passed explicitly rather than relying on the tool's defaults,
which change between versions.

---

## Chunking

156 branches at ~3 minutes each is about eight hours serially. Split it:

```bash
BRANCHES=156       # from "Found N long branches" in the log of a first run
CHUNKS=8           # 8 chunks x 2 threads = 16, the ceiling

python - <<'EOF' > chunks.txt
import sys; sys.path.insert(0, "scripts")
from ScreenReferenceRecombination import chunk_bounds
for lo, hi in chunk_bounds(156, 8):
    print(lo, hi)
EOF

while read LO HI; do
    python scripts/ScreenReferenceRecombination.py \
        --db "$DB" --outdir "chunk_${LO}" --skip_build \
        --threads 2 --start_index "$LO" --end_index "$HI" &
done < chunks.txt
wait
```

`chunk_bounds` guarantees the ranges are contiguous, cover every branch exactly
once, and differ in size by at most one. Build the MAT once first, then point
every chunk at it with `--skip_build`.

Keep `chunks × threads ≤ 16`.

The script prints the split-it-up hint automatically whenever it sees more than
20 long branches and was not told `--chunks`.

---

## What gets stored

`--write_db` creates one table and one column.

**`reference_recombination`** — one row per (accession, event). RIPPLES reports
events against internal nodes, so `descendants.tsv` is what maps a node back to
accessions; each accession beneath a recombinant node gets a row.

| column | meaning |
|---|---|
| `primary_accession` | the sequence, or NULL if the node had no descendants listed |
| `recomb_node_id` | the node RIPPLES called recombinant |
| `breakpoint_1_interval`, `breakpoint_2_interval` | as RIPPLES writes them, e.g. `3200-3260` |
| `breakpoint_1_start/end`, `breakpoint_2_start/end` | the same, parsed to integers for querying |
| `donor_node_id`, `acceptor_node_id` | the two proposed parents |
| `donor_is_sibling`, `acceptor_is_sibling` | whether each parent is a sibling of the recombinant |
| `donor_parsimony`, `acceptor_parsimony` | parsimony of each partial placement |
| `original_parsimony`, `min_starting_parsimony`, `recomb_parsimony` | the scores the call rests on |
| `detected_by`, `detected_at` | provenance |

The column names come from the RIPPLES binary itself, not from documentation.
The two interval columns are renamed (`breakpoint-1_interval` →
`breakpoint_1_interval`) because a hyphen is not usable in an unquoted SQL
identifier. A unique index on
`(primary_accession, recomb_node_id, breakpoint_1_interval, breakpoint_2_interval)`
makes re-running idempotent.

**`meta_data.recombination_status`** — added if absent:

| value | meaning |
|---|---|
| `recombinant_reference` | RIPPLES found an event covering this sequence |
| `screened_no_evidence` | screened, nothing found |
| NULL | never screened |

The distinction between the last two matters. A clean screen and no screen are
different states, and a query that treats NULL as "clean" would be wrong. A
re-run that finds nothing **clears** a stale flag rather than leaving it.

---

## Reading the results

```sql
-- which references are flagged, and where are the breakpoints?
SELECT primary_accession, breakpoint_1_interval, breakpoint_2_interval,
       donor_node_id, acceptor_node_id, recomb_parsimony, original_parsimony
FROM reference_recombination
WHERE primary_accession IS NOT NULL
ORDER BY CAST(breakpoint_1_start AS INTEGER);

-- how much of the reference set has actually been screened?
SELECT recombination_status, COUNT(*) FROM meta_data
WHERE lower(COALESCE(accession_type,'')) IN ('reference','master')
GROUP BY 1;
```

An event is worth believing in proportion to how much `recomb_parsimony`
improves on `original_parsimony`. `--parsimony_improvement` sets the floor.

**These are candidates.** RIPPLES proposes a two-parent explanation that fits
better than one parent; that is evidence, not proof. Before acting on a flag —
especially before removing a sequence from the reference set — confirm it with an
independent method (3SEQ, RDP, or a windowed identity scan against genotype
prototypes) and check the breakpoint falls somewhere biologically plausible. For
HCV that is usually the NS2/NS3 junction.

---

## Failure modes handled

| situation | behaviour |
|---|---|
| no `accession_type='master'` row | error naming the problem |
| master absent from `sequence_alignment` | error |
| ragged reference alignments | error, before `faToVcf` sees them |
| no tree covers the reference set | error |
| RIPPLES finds nothing | empty table, everything marked `screened_no_evidence` |
| RIPPLES writes no header (it omits it when there are no events) | parsed anyway |
| a malformed breakpoint interval | that interval becomes NULL; the rest of the event is kept |
| an event whose node has no descendants listed | stored with NULL accession rather than dropped |
| `--timeout` expires | partial results kept and parsed |
| a required binary is missing | error naming the binary and the conda install |

The principle throughout: never let a partial or malformed result masquerade as
"no recombination found".
