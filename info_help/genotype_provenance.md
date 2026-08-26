# Genotype provenance: where a stored genotype came from

Every genotype and subtype in `meta_data` now records which source produced it.
This document explains the vocabulary, the precedence, and the problem it solves.

---

## The problem

A genotype in `meta_data.nearest_reference_genotype` could have come from a
curator, from a phylogenetic tree, from an EPA-ng placement, or from a BLAST top
hit. The column looked identical in every case, so a confidently-curated value
and a guess from a distant BLAST match were indistinguishable.

That is not hypothetical. In `generic/hcv/ref_list_subtype_genotype.txt`,
**69 of 238 references are curated with `subtype = NA`** — a curator stating that
the genotype is known but the subtype is not assignable. All 69 come out of the
database carrying an invented subtype letter, sourced from whatever sequence they
happened to align against:

```
JX183558   curated NA   ->   stored genotype 6, subtype 'a'
KC197239   curated NA   ->   stored genotype 2, subtype 'a'
```

The values are plausible — 6a and 2a are real subtypes — which is what makes them
dangerous.

**Current behaviour is deliberately unchanged**: a curated `NA` is still treated
as missing and may be filled by inference. What changes is that the fill is now
*visible*. `subtype_origin` will read `blast_tophit`, not `curated_reflist`, so a
consumer can tell the difference and filter on it.

---

## The vocabulary

`meta_data` gains two columns, `genotype_origin` and `subtype_origin`. Both take
one of:

| token | meaning |
|---|---|
| `curated_reflist` | The reference list stated it. A human decision. |
| `tree_usher` | Inherited from the UShER tree neighbourhood. |
| `tree_iqtree` | Inherited from the IQ-TREE backbone neighbourhood. |
| `epa_placement` | From an EPA-ng phylogenetic placement. |
| `blast_tophit` | From the genotype of the best BLAST hit. |
| `gisaid_declared` | Declared in GISAID metadata, not inferred. |
| `ncbi_declared` | Declared in GenBank metadata, not inferred. |
| `unresolved` | No source supplied a value. |

---

## Precedence

First source that supplies a value wins, and the value it supplied is what gets
stored:

```
curated_reflist  >  tree_usher  >  tree_iqtree  >  epa_placement
                 >  blast_tophit  >  gisaid_declared  >  ncbi_declared
                 >  unresolved
```

**Curation outranks inference.** A curated reference-list entry is never
overwritten by a tree or a BLAST hit. That is the direct fix for the class of
problem above — and it was a deliberate choice over the alternative of letting
phylogeny win, on the grounds that a curator who wrote a genotype down knew
something the tree does not.

Trees outrank placement, which outranks BLAST, because a tree uses the whole
alignment while a top hit uses one pairwise comparison. Vendor-declared values
rank below inference because a submitter's declared genotype is frequently
wrong, and unlike the curated list nobody has checked it.

---

## Genotype and subtype are resolved together

A previous defect resolved genotype and subtype by two independent climbs
through the fallback chain, so a sequence could end up labelled with a genotype
from one neighbourhood and a subtype from another — a combination that exists
nowhere in the data. Both are now taken from the same source wherever the source
supplies both.

Where a source supplies a genotype but no subtype, the subtype falls through to
the next source and `subtype_origin` records that lower source. The two origin
columns can therefore legitimately differ, and when they do it is because the
data genuinely came from two places — not because of an accident.

---

## Reading it

```sql
-- how were genotypes actually determined in this build?
SELECT genotype_origin, COUNT(*) FROM meta_data GROUP BY 1 ORDER BY 2 DESC;

-- everything whose subtype was inferred rather than curated
SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype,
       genotype_origin, subtype_origin
FROM meta_data
WHERE subtype_origin NOT IN ('curated_reflist', 'unresolved');

-- disagreement between the two
SELECT * FROM meta_data WHERE genotype_origin <> subtype_origin;
```

If a genotype looks wrong, the origin column tells you which mechanism to go and
examine.

---

## Effect on other viruses

Both columns are written for every virus, because every virus goes through the
same resolution chain. For a virus with no genotypes at all, every row reads
`unresolved`, which is accurate and costs two small text columns.

Nothing about this requires subtypes to exist. `subtype_origin` reads
`unresolved` when there is no subtype to resolve.

---

## Related: `info.creation_type`

Adjacent provenance problem, fixed at the same time. `create_db` only fell back
to `"last updated"` when `db_status` was empty, and the pipeline never passed
`-ds`, so **every `--update` run stamped itself `"new db"`**. The shipped
`test_out/update_test` database contains two `"new db"` rows — one for its
creation and one for its update, indistinguishable.

The pipeline now passes the status explicitly, so a database records whether it
was built fresh or updated. Given a database you can now answer both "where did
this genotype come from" and "was this thing built once or incrementally".
