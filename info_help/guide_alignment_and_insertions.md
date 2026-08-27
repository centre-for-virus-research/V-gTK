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
everything against the wrong coordinate frame and fail silently. Currently only
the influenza profiles set the parameter (`nextflow.config`); HCV and RABV can
adopt it by pointing it at a directory. When it is unset Nextflow passes the
string `UNSET`, which is normalised to "no guide alignment" rather than treated
as a path.

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
