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

## Not yet wired in

No pipeline process calls `TreeReRoot.py`. The parameter is declared and the
script is tested; nothing consumes it yet.
