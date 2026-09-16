# Guide alignments and insertions

What happens to a region that exists in a sub-reference but not in the master,
and how to make sure it survives.

---

## The failure mode this avoids

In the large SARS-CoV-2 alignments everything is projected onto Wuhan-Hu-1
coordinates. Any region inserted relative to that reference has nowhere to go, so
it is stripped. The alignment stays tidy and the evolutionary signal inside those
regions is gone — you cannot see a substitution within an insertion, because the
insertion has no columns to hold it.

The same thing happens here **if no guide alignment is supplied**. Nextalign
aligns to reference coordinates by construction; insertions are pulled out to
`*.insertions.csv`, which records that sequence X had N bases after position P.
That is honest bookkeeping and enough for presence/absence, but it is not an
alignment. Nothing downstream can compare two insertions to each other.

## What a guide alignment changes

`PadAlignment.process_master_alignment` iterates the records of the **guide
alignment**, not of the master. For each reference it takes that reference's own
gapped row and projects the queries nextalign assigned to it into that row's
coordinate space (`insert_gaps`, via a banded alignment between the guide
reference and the nextalign reference). The merged MSA therefore lives in the
guide alignment's column space.

The consequence: columns where the master is gapped and a sub-reference has bases
**exist in the output**, and queries carrying the insertion put their residues
there. Curating a sub-reference with an interesting insertion into the guide
alignment is how you make that region analysable.

Measured on a minimal case — master `ACGTACGT------ACGTACGT`, sub-reference
`ACGTACGTGGGGGGACGTACGT`:

| row | insertion columns |
|---|---|
| master | `------` |
| sub-reference | `GGGGGG` |
| query with the insertion | `GGG`**`T`**`GG` |
| query without it | `------` |

The substitution inside the insertion is retained and alignable against every
other sequence that has the region. Flanking bases do not shift. Every row is the
guide alignment's width.

Standing test: `tests/unit/test_guide_alignment_insertions.py`.

---

## How to supply one

| virus type | how |
|---|---|
| segmented | a directory of `refset_<segment>_aln.fasta`, one per segment |
| non-segmented | a directory containing a **single** FASTA, any name |

```groovy
ref_set_aligned = "${projectDir}/generic/<virus>/ref_set_aligned"
```

Nothing gates this on `is_segmented`. Where no segment value is available the
resolver uses the sole FASTA in the directory; if there is more than one it
**refuses rather than guessing**, because picking the wrong backbone would align
everything against the wrong coordinate frame and fail silently. The influenza
profiles point the parameter at a curated directory (`nextflow.config`).

## When none is supplied: the built backbone

Since September 2026 a fresh build with no `ref_set_aligned` **builds one**
(`BUILD_REFERENCE_ALIGNMENT`, `scripts/BuildReferenceAlignment.py`), so RABV,
HCV and any other simple genome get shared reference indels as columns without
curating anything. It runs after BLAST and uses only the reference FASTA, the
reference list and each master's GFF.

| step | why |
|---|---|
| MAFFT on **whole reference records** | the master's row is its entire record, so the 5' trim `CalcAlignmentCord` corrects for is zero |
| each master CDS **re-aligned as protein and back-translated** | a nucleotide aligner puts gaps anywhere: on HCV, plain MAFFT left 8 of 34 master gap runs inside genes out of frame. Gaps inside a gene are now whole codons |
| references with a broken frame in a gene added by nucleotide (`--add --keeplength`) | they still get a row without disturbing the codon alignment |
| insertion columns carried by fewer than `ref_aln_min_insertion_support` references (default 2) dropped | one unusual reference would otherwise widen every row in the database; on HCV 326 of 523 insertion columns were single-reference. The dropped bases go to the insertions table, not nowhere |
| the most accurate MAFFT strategy the size allows | L-INS-i within a DP budget, FFT-NS-i above it |

The step fails rather than writing a backbone that breaks a guarantee: the master
row must degap to the master record, no master gap run inside a codon-aligned
gene may be out of frame, and every reference base must be in its row or
recorded in `dropped_insertions.tsv`.

Outputs, published under `ref_set_aligned/`: `refset_<segment>_aln.fasta` (the
same naming the resolver expects; `<segment>` is the master's segment, `0` if it
has none), `dropped_insertions.tsv`, `build_report.tsv`.

| params | effect |
|---|---|
| `ref_set_aligned = <dir>` | use a curated backbone; nothing is built |
| `build_ref_alignment = false` | build nothing; references are projected onto the master and their insertions stripped, as before |
| `ref_aln_min_insertion_support = 1` | keep every reference insertion column |

**Update mode never builds one.** It keeps deriving the backbone from the
reference rows already in the database, so the column layout of an existing
build does not move. Adding references in update mode is not supported by this
path yet.

When no backbone is used at all, Nextflow passes the string `UNSET`, which is
normalised to "no guide alignment" rather than treated as a path.

### Caching a built backbone

The build is deterministic for a given reference list, insertion-support setting
and builder, and for HCV's 238 references it is the slowest step of a fresh build
(826 s on 8 threads). So a profile can point `ref_set_aligned` at a backbone built
once and kept under `generic/<virus>/ref_set_aligned/`. The HCV profiles
(`HCV_test`, `HCV_full`, `HCV_XML_full`) do. A cached backbone is used exactly like a
built one: its `dropped_insertions.tsv` still reaches the insertions table.

A cache goes stale silently, so `CHECK_REFERENCE_BACKBONE`
(`scripts/CheckReferenceBackbone.py`) runs before anything aligns against a
supplied backbone. It reads `backbone_manifest.tsv` beside the cache and:

- **fails** if the reference list's accessions differ from the backbone's rows,
  or if `ref_aln_min_insertion_support` differs from the cached setting;
- **warns** if `BuildReferenceAlignment.py` has changed since the cache was made;
- passes a directory with no manifest through untouched, as a hand-curated
  backbone (the influenza ones).

To rebuild a cache:

```bash
nextflow run vgtk-init.nf -profile HCV_test --ref_set_aligned null   # builds and publishes <publish_dir>/ref_set_aligned
cp <publish_dir>/ref_set_aligned/{refset_*_aln.fasta,dropped_insertions.tsv,build_report.tsv} generic/hcv/ref_set_aligned/
python scripts/CheckReferenceBackbone.py --backbone_dir generic/hcv/ref_set_aligned \
    --ref_list generic/hcv/ref_list_subtype_genotype.txt --min_insertion_support 2 \
    --write_manifest --source "<which run, when>"
```

`--ref_set_aligned null` means "no backbone supplied" and triggers a build, so
any profile's cache can be bypassed from the command line.

### Downstream effects of a wider backbone

- `MMseqsClustering --ref_list`: two-step completeness is measured against the
  master's length, not the alignment width, or complete genomes look incomplete.
- `UsherPlacement --ref_list`: the faToVcf reference is the master when it is a
  cluster representative. faToVcf drops columns where its reference is gapped,
  so UShER does not see insertion columns; IQ-TREE does.
- `GenerateTables --reference_insertions`: reference insertions come from the
  builder's `dropped_insertions.tsv`, not from Nextalign's reference-vs-master
  run, which would count bases already held as columns a second time.

## Choosing what goes in it

The guide alignment defines the coordinate system, so it is a curation decision,
not an automatic one. Include a reference or sub-reference when its insertion is
one you want to study; each one you add widens every row in the merged MSA by the
length of that insertion. The master should stay the master — the point is not to
replace it but to give the alignment somewhere to put what the master lacks.

## Known defect

A reference present in both the guide alignment and its own nextalign
subalignment is written to the merged MSA twice. It inflates sequence counts and
can confuse dedup by name. Tracked as an `xfail` in the test module above.

---

## Related

`ref_set_aligned` is separate from clustering. Whether a genome with an insertion
is *chosen* as a cluster representative is an MMseqs2 question
(`--cluster-mode 2`, longest member wins, so an insertion-bearing genome is
favoured but not guaranteed); whether its insertion *survives into the alignment*
is this question. A genome can be a representative and still lose its insertion
if there is no guide alignment, and a curated guide alignment entry keeps its
insertion whether or not it was ever a representative.
