# Missing H3 (and everything that isn't H1/N2) in the IAV database

> **Fixed.** The rebuild restored H3 retention to 99.89% (segment 4) and N1 to 99.82% (segment 6), and the projectability and segment-retention guards in `CollectFilteredSequences.py` and `CheckSegmentRetention.py` prevent recurrence.

**Build investigated:** `test_out/IAV_DB_genbank_n_gisaid` (pipeline run, 27 Jul 2026) and the
derived rebuild `test_out/IAV_DB_genbank_n_gisaid_rebuild/db_build/iav-db.db` (26 Aug 2026).
**Trigger:** report of "a cluster of H1 sequences placed inside the H5 clade" when colouring by
`subtype_blast` (subtree from `node_8`, cluster at `node_39151`).

---

## 1. Summary

The HA (segment 4) and NA (segment 6) alignments — and therefore the clustering, the base trees and
every tree placed on them — contain **only H1 and only N2**. Everything else was silently deleted
between the Nextalign step and the padded merged MSA.

| HA subtype | in `meta_data` (segment 4) | in `refset_4_aln_merged_MSA.fasta` | kept |
|---|---|---|---|
| **H3** | **242,732** | **46** | **0.02 %** |
| H1 | 201,575 | 201,259 | 99.8 % |
| H5 | 63,008 | 2,456 | 3.9 % |
| H9 | 26,616 | 8 | 0.03 % |
| H7 | 8,534 | 8 | 0.09 % |
| H6 | 6,677 | 11 | 0.16 % |
| H4 | 5,753 | 8 | 0.14 % |
| H2, H10–H19 | ~10,500 | ~45 | <0.5 % |
| **total** | **573,405** | **205,487** | **35.8 %** |

| NA subtype | in `meta_data` (segment 6) | in `refset_6_aln_merged_MSA.fasta` | kept |
|---|---|---|---|
| **N2** | 245,169 | 244,877 | **99.9 %** |
| N1 | 224,200 | 64 | 0.03 % |
| N8 | 12,266 | 11 | 0.09 % |
| N6 | 10,849 | 59 | 0.5 % |
| N3–N11 | ~16,600 | ~30 | <0.2 % |
| **total** | **516,559** | **248,414** | **48.1 %** |

Segments 1, 2, 3, 5, 7 and 8 are unaffected (each retains >99.8 %).

The cause is a **projectability test that asks the wrong question**. Whether a BLAST reference is
usable is decided by "does nextalign manage to align this reference against one of the eight
*master* sequences?" The segment-4 master is `AB573800` (**H1**N1) and the segment-6 master is
`AB472016` (H5**N2**). References from any other HA/NA subtype are too divergent from those masters
for nextalign's seed matcher, so they fail — and every query sequence assigned to them is added to
`filtered_sequences_ids.txt` and skipped by `PadAlignment`.

**633,988 sequences (16 % of the database) were dropped this way.** They are almost entirely HA and
NA of non-H1 / non-N2 subtype.

---

## 2. Evidence chain

Each step below was verified against the surviving Nextflow work directories of the 27 Jul run.

### 2.1 BLAST assignment is correct

`inputs/query_uniq_tophit_annotated.tsv` assigns 572,876 segment-4 queries across 113 of the 117 HA
references, with a sensible distribution (`MZ198323` H3N2 → 219,844; `MH559851` H1N1 → 170,711;
`DQ864721` H5N1 → 60,751; `JF795094` H9N2 → 21,427 …).

Cross-checking the BLAST-derived subtype against the GenBank/GISAID `serotype` field for all 565,081
segment-4 records where both are typed: **99.955 % agreement, 254 disagreements**. So
`subtype_blast` itself is sound, and the reference set is sound — a maximum-likelihood tree of the
117 aligned HA references recovers every subtype as a clean monophyletic group.

### 2.2 Nextalign produced the sequences

`work/e2/0965e0dd28b032b7cd95bf163c7b01/Nextalign/query_aln/` contains full per-reference
sub-alignments, including the H3 ones:

```
MZ198323.aligned.fasta   219,716 sequences   (394 MB)
MH559851.aligned.fasta   170,555
DQ864721.aligned.fasta    60,409
JF795094.aligned.fasta    21,426
```

Genuine nextalign failures across the whole run total only ~6,300 sequences (too short, ambiguous
indels, low seed match). The H3 data existed and was correctly aligned.

### 2.3 The reference-vs-master alignment is where the subtypes are lost

`NextalignAlignment.process()` (`scripts/NextalignAlignment.py:280`) also aligns the *entire*
418-sequence reference set against each of the 8 master sequences, into
`Nextalign/reference_aln/<master>/<master>.aligned.fasta`. Nextalign silently omits any reference it
cannot seed-match. The result:

| master | segment | master serotype | references aligned | references for that segment |
|---|---|---|---|---|
| KJ889313 | 1 | — | 28 | 32 |
| CY087822 | 2 | — | 29 | 32 |
| MW333949 | 3 | — | 34 | 37 |
| **AB573800** | **4 (HA)** | **H1N1** | **30** | **117** |
| NC_007369 | 5 | — | 52 | 54 |
| **AB472016** | **6 (NA)** | **H5N2** | **29** | **83** |
| NC_002016 | 7 | — | 23 | 28 |
| KY243296 | 8 | — | 30 | 35 |

The 30 HA references that survived against `AB573800` are:

```
AB573800 CY004504 CY085702 CY160213 CY195263 EF101749 GU183815 KC345638 KC859230 KM276920
KR700735 KR701057 KU976737 LC430960 MF370343 MG962504 MH559851 MK339605 MK340189 MW772910
MW848773 MZ088617 MZ477681 MZ945794 NC_002017 OK205843 OP613226 OR674119 OR861585 PP829968
```

— every one of them H1 (or `MZ945794`, "HxN2", which sits inside the H1 clade), plus two H5
references (`OK205843`, `PP829968`) that happened to clear the seed threshold. The 29 NA references
that survived against `AB472016` are **all N2** without exception.

### 2.4 That list becomes a deletion list

`CollectFilteredSequences.collect_unprojectable_queries()`
(`scripts/CollectFilteredSequences.py:27-79`) builds `aligned_reference_ids` **from
`reference_aln/*/*.aligned.fasta`** — i.e. from the table above — and then flags every query in any
`query_aln/<ref>/` directory whose `<ref>` is not in that set:

```python
aligned_reference_ids = set()
for ref_aln in nextalign_path.glob("reference_aln/*/*.aligned.fasta"):   # <-- wrong source
    ...
for query_ref_dir in nextalign_path.glob("query_aln/*"):
    ref_id = _normalize_accession(query_ref_dir.name)
    if ref_id in aligned_reference_ids:
        continue
    # ...every sequence in this group -> filtered
```

Output (`test_out/IAV_DB_genbank_n_gisaid/filtered_sequences_summary.txt`):

```
filtered_sequences_total          640286
unprojectable_queries_total       633988
unprojectable_reference_count     140
```

140 references = the 87 HA + 54 NA references excluded above (a couple of others). The remaining
6,298 entries are the real nextalign failures.

### 2.5 `PadAlignment` executes the deletion

`vgtk-init.nf:612` passes that list straight through:

```
python scripts/PadAlignment.py -nd Nextalign -m ref_list_refmast.txt -o . -d . \
  -i Nextalign/query_aln --keep_intermediate_files \
  --precomputed_ref_dir ".../generic/influenza/ref_set_aligned" \
  --skip_ids filtered_sequences_ids.txt
```

and `PadAlignment.process_master_alignment()` (`scripts/PadAlignment.py:440-445`) drops them:

```
[skip_ids] Loaded 640286 accessions to skip from filtered_sequences_ids.txt
...
[skip_ids] Skipped 219715 sequence(s) from MZ198323 subalignment    <-- all the human H3N2
[skip_ids] Skipped 157052 sequence(s) from CY053683 subalignment
[skip_ids] Skipped  60408 sequence(s) from DQ864721 subalignment    <-- all the H5N1
[skip_ids] Skipped  21425 sequence(s) from JF795094 subalignment    <-- all the H9N2
```

The resulting `refset_4_aln_merged_MSA.fasta` has 205,614 records. Walking it block by block, the 30
reference groups that contain any query sequence are **exactly** the 30 references from §2.3; the
other 87 blocks contain the reference sequence alone.

**This is the crux of the defect.** `PadAlignment` was given `--precomputed_ref_dir`, so it projects
each sub-alignment through `generic/influenza/ref_set_aligned/refset_4_aln.fasta`, which contains
**all 117** HA references and is a perfectly good MSA (1,782 columns, 123–177 gaps per reference).
Nothing about segment 4 was actually unprojectable. `CollectFilteredSequences` never looks at that
directory — it answers the projectability question from `reference_aln/`, a completely different and
far stricter artefact, and its (wrong) answer overrides the correct one.

### 2.6 Consequence for clustering and the trees

`refset_4_aln_merged_MSA_dedup.fasta` → 203,605 unique sequences → **9,470 MMseqs clusters**:

| representative subtype | clusters | members |
|---|---|---|
| H1 | 9,149 | 199,514 |
| H5 | 172 | 2,391 |
| **H3** | **19** | **46** |
| all others | ~130 | ~200 |

The IQ-TREE backbone (`segment_redo/IQTree_MMseqClusters_refset_4_aln_merged_MSA_dedup/iqtree.treefile`)
therefore has 9,470 tips of which 9,149 are H1. The only non-H1 structure of any size is a
**179-tip H5 island** (172 H5, plus 7 tips whose GenBank `serotype` is blank or `mixed` but whose
closest reference is the H5N2 `OK205843` — those are real H5). The backbone itself is internally
consistent: the H5 clade is 96 % pure, and the members of every cluster whose representative sits in
it are 2,390/2,486 H5 with `subtype_blast` = H5 for all but 9.

The DB build log records the damage plainly
(`test_out/IAV_DB_genbank_n_gisaid_rebuild/db_build/db_summary.txt`):

```
Segment | Passed QC | Failed QC
   1    |  462975   |    925
   2    |  454052   |   1033
   3    |  465133   |    956
   4    |  205397   | 367894      <-- HA
   5    |  472024   |    798
   6    |  248412   | 268067      <-- NA
   7    |  492067   |   1044
   8    |  474133   |    734

alignment_filtering: reference not present in master-projected reference_aln;
query cannot be projected into merged segment alignment | 633988
```

`CollectFilteredSequences.py:313-318` even prints an `[ALERT]` to stderr naming that 633,988 figure.
Nothing in the pipeline treats it as an error, so the run completes with exit code 0.

---

## 3. How this produces "H1 sequences inside the H5 clade"

The HA backbone that every placement is made against is, in effect, *an H1 tree with one H5 island
hanging in it*. There is no H3 backbone (19 tips), no H9, no H7, no H2, no H6 — 64 % of the HA data
never reached the alignment at all, and the subtypes that were removed are precisely the ones that
would have broken up the deep structure of the tree.

Two consequences follow, and between them they account for the reported cluster:

1. **UShER has nowhere correct to put a divergent sequence.** Parsimony placement onto this backbone
   can only choose between the dense H1 radiation and the single long-branch H5 island. Any sequence
   carrying many differences that no H1 backbone tip shares is most parsimoniously attached deep —
   which on this tree means inside or below the H5 clade. Because such sequences share the *same*
   systematic cause, they arrive together and appear as a coherent cluster rather than as scattered
   singletons. That is the shape of the artefact at `node_39151`.

2. **The `subtype_blast` label is right and the topology is wrong.** `subtype_blast` is the subtype
   of the BLAST closest reference, and §2.1 shows it is 99.96 % concordant with the submitted
   serotype. So a cluster labelled H1 sitting inside the H5 clade should be read as *"these really
   are H1 and the tree has put them in the wrong place"*, not as a labelling error. (For
   completeness: only 54 sequences in the whole segment-4 MSA have a metadata/BLAST subtype
   disagreement, and only 2 of those are H5-by-metadata called H1 by BLAST — far too few to be the
   cluster in the figure.)

I could not confirm this against the exact tree in the plot: `test_out/IAV_trees_31_8_26/` holds
UShER trees for segments 1, 2, 3, 5, 7 and 8 only, and no segment-4 tree with `node_*` labels exists
under `test_out/`. If you can point me at the seg4 `.tree` file behind that figure I can extract the
`node_39151` descendants and confirm directly which reference group and which alignment blocks they
came from. What is certain is that the backbone they were placed on is missing 64 % of HA, including
99.98 % of all H3.

---

## 4. Fix

### 4.1 The bug

`collect_unprojectable_queries()` must test projectability against the alignments `PadAlignment`
actually projects through, not against `reference_aln/`:

* when `--precomputed_ref_dir` is in use (the flu path), the authority is
  `generic/influenza/ref_set_aligned/refset_<segment>_aln.fasta` — a reference is projectable if it
  appears in the alignment **for its own segment**;
* only when there is no precomputed directory does `reference_aln/` become the right source.

`CollectFilteredSequences.py` currently takes only `-n <nextalign_dir>`; it needs the same
`--precomputed_ref_dir` / ref-list arguments `PadAlignment` gets, and the two must share one
projectability function so they cannot diverge again.

Note that `PadAlignment` already has the correct check —
`find_orphan_query_references()` (`PadAlignment.py:402-427`) tests against the per-segment
`reference_alignment_file`. Its per-segment warnings ("366 references … skipped") are themselves
noise, because on the segment-1 pass every segment-2..8 reference is trivially "absent"; that
warning should be scoped to references belonging to the segment being processed.

### 4.2 Guard rails

* **Fail the run, don't warn.** A `--max_filtered_fraction` (say 2 %) that aborts
  `COLLECT_FILTERED_SEQUENCES` would have stopped this build at the first attempt. 16 % of a database
  disappearing should not exit 0.
* **Per-segment retention check.** Assert `merged_MSA` count ≈ BLAST-assigned count per segment. The
  35.8 % / 48.1 % figures for segments 4 and 6 against >99.8 % everywhere else are unmissable.
* **Subtype-coverage assertion for flu.** After clustering, check that the representative-subtype
  distribution is not wildly different from the input-subtype distribution. 242,732 H3 sequences
  yielding 19 cluster representatives is a hard error, not a warning.

### 4.3 Reference-set robustness (secondary)

`nextalign_master()` (`NextalignAlignment.py:124`) has a five-step relaxation cascade, but it treats
"nextalign exited 0 and the file passes frame validation" as success even when nextalign silently
dropped 87 of 117 input references — a whole-file check cannot see per-sequence omissions. Even with
the cascade, aligning an H3/H7/H13 HA against an H1 master is not a well-posed problem. Two options:

* skip `nextalign_master` entirely on the flu path, since `--precomputed_ref_dir` already supplies
  hand-curated per-segment reference MSAs and nothing downstream needs `reference_aln/`; or
* if `reference_aln/` is needed, have `nextalign_master` assert that the output contains every input
  reference for that segment and fail loudly if not.

### 4.4 Rebuild scope

The fix invalidates segments 4 and 6 only. Segments 1, 2, 3, 5, 7 and 8 (alignments, clusters,
IQ-TREE backbones, and the UShER trees already in `test_out/IAV_trees_31_8_26/`) are unaffected and
do not need to be recomputed. For HA and NA the whole chain from `PAD_ALIGNMENT` onward has to be
redone — and the HA backbone will grow substantially, since 242,732 H3 sequences currently
contribute 19 tips.

---

## 5. Reproducing the key numbers

```bash
BASE=/home3/oml4h/RABV-gTK/test_out/IAV_DB_genbank_n_gisaid
REB=/home3/oml4h/RABV-gTK/test_out/IAV_DB_genbank_n_gisaid_rebuild
NA_WORK=/home3/oml4h/RABV-gTK/work/e2/0965e0dd28b032b7cd95bf163c7b01
PAD_WORK=/home3/oml4h/RABV-gTK/work/07/744c10f5555ccc1916c50eb1d298b5

# alignment size vs BLAST assignment, per segment
for i in 1 2 3 4 5 6 7 8; do
  echo "seg$i aln=$(grep -c '>' $BASE/refset_${i}_aln_merged_MSA.fasta)" \
       "blast=$(awk -F'\t' -v s=$i '$5==s' $REB/inputs/query_uniq_tophit_annotated.tsv | wc -l)"
done

# which references nextalign could align to each master
for m in $(ls $NA_WORK/Nextalign/reference_aln); do
  echo "$m -> $(grep -c '>' $NA_WORK/Nextalign/reference_aln/$m/$m.aligned.fasta)"
done

# the deletion, in the pipeline's own words
grep 'skip_ids] Skipped' $PAD_WORK/.command.log | sort -k3 -nr | head
cat $BASE/filtered_sequences_summary.txt
```
