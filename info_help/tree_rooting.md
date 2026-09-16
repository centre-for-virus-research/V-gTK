# Tree rooting

Script: [`scripts/TreeReRoot.py`](../scripts/TreeReRoot.py)
Tests: `tests/unit/test_tree_reroot.py`

**Everything is midpoint rooted.** No profile sets an outgroup. Midpoint needs no
prior knowledge and makes no claim about which lineage is basal, which is the
right default when nobody has a well-supported outgroup.

The mechanism for outgroup rooting exists and is tested, for a species that
later turns out to have a defensible one:

```groovy
root_outgroup = null                    // midpoint (every profile today)
root_outgroup = "ACC1,ACC2,ACC3"        // outgroup, if ever warranted
```

Comma-separated, because a Nextflow parameter is a single string. The script also
accepts the space-separated form for direct command-line use.

---

## HCV: genotype 8 was considered and rejected

Genotype 8 is the obvious candidate - the most divergent HCV lineage, and
monophyletic in this reference set (4/4 tips, checked against ICTV; see
[`hcv_tree_validation_ictv.md`](hcv_tree_validation_ictv.md)). It works
mechanically. Measured on the ICTV reference tree:

| rooting | one side | other side |
|---|---|---|
| midpoint | 43 tips: genotypes **2 + 7** | 195 tips |
| genotype 8 outgroup | 4 tips: genotype 8 | 234 tips |

The two conventions genuinely disagree about which lineage is basal, so the
choice is not cosmetic.

**It was rejected anyway.** Those four sequences (`MH590698`-`MH590701`) are all
Canadian, collected 2015-2017, from a lineage first described in 2018. Anchoring
the entire HCV tree on a group that geographically narrow and that recently
characterised asserts more than the data supports. Midpoint asserts nothing.

The supplied ICTV reference tree is itself already midpoint rooted; re-running
midpoint on it reproduces it exactly. Rerooting never changes the unrooted
topology - every rooting tried above gave RF = 0 against the input.

## Two things learned while evaluating it

**Do not root on a single sequence.** `MH590700` is the best-covered genotype 8
reference (9,547 nt of ~9,646, all 10 genes at >=95.2%). Rooting on it alone
still gives sides of **1 and 237 tips**, leaving the other three genotype 8
references in the ingroup. The script emits a note whenever the outgroup is a
single tip.

**Expanding one accession to its MMseqs cluster does not rescue that.** All four
genotype 8 references are *singleton* clusters at 95% identity, because HCV
genotypes are roughly 30% divergent and even within genotype 8 the four sequences
are not 95% identical. The cluster is the tip. If an outgroup is ever wanted,
name the whole clade.

The script refuses outright if the named tips are not an exclusive clade -
rooting on a paraphyletic group silently produces a different tree than intended.

## Segmented viruses

Midpoint, and the mechanism would need work before that could change: one
outgroup string cannot express a per-segment outgroup, and each segment has its
own tree.

## What happens to UShER node labels

They are dropped, and the script says so. This is correct, and preserving them
would be worse.

A support value is a property of a *bipartition*, so it survives rerooting - the
script records every split's support before rerooting and reattaches it after.
An internal node *name* identifies a node in one particular rooting. Two nodes
can collapse onto the same split: in a bifurcating tree both root children
describe the same bipartition, so preserving names by split duplicates one and
loses the other.

Measured on a 238-tip HCV MAT, UShER assigns `node_N` deterministically from the
topology:

| input tree | parsimony | labels regenerated |
|---|---|---|
| with `node_N` labels | 224,233 | 237 |
| labels stripped | 224,233 | 237, **identical clades (237/237)** |

So UShER neither needs the labels nor loses anything without them. A label
preserved through a reroot would describe the *old* rooting while a fresh UShER
run assigns its own, giving two numbering schemes that disagree.

They remain useful *within* one tree version: RIPPLES reports recombination
against `recomb_node_id`, `donor_node_id` and `acceptor_node_id`, all `node_N`.
Node IDs recorded against a previous rooting will not line up with a rerooted
tree, which is what the script's note says.

Update placement is unaffected: `UsherPlacement.py` reads the stored newick back
as a starting tree, and UShER regenerates labels from it either way.

## Where rooting happens

`CreateSqliteDB.py` midpoint roots every tree before storing it in the `trees`
table (IQ-TREE, UShER, VeryFastTree, per-segment manifest trees), so both
`vgtk-init.nf` and `vgtk-rabv.sh` get it without a separate step. Tree-based
genotype assignment (`clade_from_tree.py`) reads the same rooted trees, so a
stored tree and the genotypes called from it always agree.

A tree that cannot be rooted (one tip, missing branch lengths) is stored
unrooted with a `[warn]` line rather than stopping the build.

`root_outgroup` is still not passed through; the build is midpoint only.

One IQ-TREE per segment. Both callers pass the first treefile as `-it`/`-ut`
and list it in the tree manifest too; the manifest row carries the segment, so
the repeat is skipped. Older databases still hold an unlabelled `iqtree` copy,
which is deleted when it matches a segment-labelled tree exactly.

Update runs never rebuild IQ-TREE. They pass no IQ-TREE at all, so the seed
database's IQ-TREE rows (and any reference-only tree) come through untouched,
and only the UShER tree is replaced. A seed built before rooting was added
therefore keeps its unrooted IQ-TREE.

Rough cost, HCV full build: 11,825-tip IQ-TREE 2s, 138,095-tip UShER tree 24s.

## Reference-only tree

Each IQ-TREE also gets a copy pruned to the accessions in the run's `ref_list`
(matched without version suffixes), then midpoint rooted. Pruned first, rooted
second: the midpoint of the reference tips alone is not where the full tree's
midpoint falls. Pass-through nodes left by pruning are removed and their branch
lengths joined, so reference-to-reference distances are unchanged.

Stored as name `iqtree_reference_only` (or `iqtree_reference_only_<segment key>`)
under source `iqtree_reference_only`. The separate source matters: UsherPlacement
and ValidateDbTree pick a backbone by source, and must not take this one for
the full tree.

The IQ-TREE only holds cluster representatives, so references merged into
another cluster are absent. The build log prints how many made it (HCV: 226 of
238).

## Why not Biopython's root_at_midpoint

It is quadratic, and it is also wrong on some small trees. Its
`root_with_outgroup` does not split a branch that hangs directly off the root,
so `(A:0.1,(B:0.2,C:0.3):0.4)` came out as `(A:0.5,(B,C):0.1)`, 0.1 longer
than the input, and a two-tip tree grew to 1.5x its length. IQ-TREE always
writes its first taxon directly under the root, so outgroup rooting hit the
same fault. `TreeReRoot.root_on_branch` does the split itself, and the tests
now check that tip-to-tip distances survive rooting rather than comparing with
Biopython.
