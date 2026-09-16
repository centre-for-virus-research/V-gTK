#!/usr/bin/env python3

import argparse
import copy
import hashlib
import os
import sys
from io import StringIO

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

'''
python scripts/TreeReRoot.py --input_tree generic/rabv/tree/ref.treefile --output_tree ref_midpoint_rooted.treefile --order_node decrease
python scripts/TreeReRoot.py --input_tree generic/rabv/tree/ref_plus_outgroup.treefile --output_tree ref_plus_outgroup_rerooted.treefile --outgroup NC_009528.2 NC_009527.1 --order_node decrease 

'''
def _neighbours(clade, parent):
    """Every node adjacent to `clade`, treating the tree as unrooted.

    A phylogeny's stored root is an arbitrary artefact of how the file was
    written; the midpoint does not depend on it. Walking children *and* parent
    is what makes the search below see the real tree.
    """
    for child in clade.clades:
        yield child, (child.branch_length or 0.0)
    up = parent.get(id(clade))
    if up is not None:
        yield up, (clade.branch_length or 0.0)


def _farthest_from(start, parent):
    """Node farthest from `start`, its distance, and the predecessor map.

    Iterative rather than recursive: a 138,095-tip tree is far deeper than the
    default recursion limit.
    """
    dist = {id(start): 0.0}
    prev = {id(start): None}
    best, best_distance = start, 0.0
    stack = [start]
    while stack:
        current = stack.pop()
        base = dist[id(current)]
        for neighbour, weight in _neighbours(current, parent):
            if id(neighbour) in dist:
                continue
            reached = base + weight
            dist[id(neighbour)] = reached
            prev[id(neighbour)] = current
            # Only a tip may be a diameter endpoint. With zero-length branches
            # an internal node can tie a leaf, and returning it would place the
            # root off the taxa-to-taxa path.
            if reached > best_distance and not neighbour.clades:
                best, best_distance = neighbour, reached
            stack.append(neighbour)
    return best, best_distance, prev


def _add_lengths(first, second):
    if first is None and second is None:
        return None
    return (first or 0.0) + (second or 0.0)


def root_on_branch(tree, node, node_side_length):
    """Place a new root on the branch above `node`, `node_side_length` from it.

    Written out rather than using ``Tree.root_with_outgroup``, which does not
    split the branch when `node` hangs directly off the root: rooting
    (A:0.1,B:0.2,C:0.3) on A with 0.05 gives A:0.05 and leaves 0.1 on the other
    side, adding 0.05 of length the tree never had. IQ-TREE always writes its
    first taxon under the root, so that case is the common one, not a corner.

    Each node on the path from `node` up to the old root is flipped to become
    the child of the node it used to hang from, taking that branch's length
    with it. Linear in the tree size. Returns the tree.
    """
    if node is tree.root:
        raise ValueError("Cannot root on the branch above the current root.")

    parent = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parent[id(child)] = clade

    total = node.branch_length
    root_children = tree.root.clades
    if parent[id(node)] is tree.root and len(root_children) == 2:
        # Already a two-child root on this unrooted branch: just move the root
        # along it. Measured over both halves, which are one branch.
        sibling = root_children[0] if root_children[1] is node else root_children[1]
        if total is not None and node_side_length is not None:
            total += sibling.branch_length or 0.0
            node.branch_length = min(max(node_side_length, 0.0), total)
            sibling.branch_length = total - node.branch_length
        tree.root = Clade(clades=[node, sibling])
        return tree

    if total is None or node_side_length is None:
        node_length = other_length = None
    else:
        node_length = min(max(node_side_length, 0.0), total)
        other_length = total - node_length

    path = [parent[id(node)]]
    while id(path[-1]) in parent:
        path.append(parent[id(path[-1])])

    path[0].clades = [child for child in path[0].clades if child is not node]
    length = other_length
    for index, clade in enumerate(path):
        above_length = clade.branch_length
        clade.branch_length = length
        if index + 1 < len(path):
            above = path[index + 1]
            above.clades = [child for child in above.clades if child is not clade]
            clade.clades.append(above)
            length = above_length

    node.branch_length = node_length
    new_root = Clade(clades=[node, path[0]])

    # The old root now hangs off the path. With a single child left it is a
    # pass-through node, so splice it out and carry its length down.
    old_root = path[-1]
    if len(old_root.clades) == 1:
        holder = path[-2] if len(path) > 1 else new_root
        only = old_root.clades[0]
        only.branch_length = _add_lengths(only.branch_length, old_root.branch_length)
        holder.clades = [only if child is old_root else child for child in holder.clades]

    tree.root = new_root
    return tree


def midpoint_root(tree):
    """Root `tree` at the midpoint of its two most distant tips, in O(n).

    Biopython's ``Tree.root_at_midpoint`` re-roots the tree at every tip in
    turn and recomputes all depths each time, which is quadratic. Measured on
    this repository's HCV trees: 250 tips 0.07s, 500 0.28s, 1000 1.21s, 2000
    6.16s - a clean 4x per doubling. Extrapolated, the 11,825-tip IQ-TREE takes
    about 4 minutes and the 138,095-tip UShER tree about 8 hours. It is also
    wrong on some small trees: it roots (A:0.1,(B:0.2,C:0.3):0.4) as
    (A:0.5,(B,C):0.1), 0.1 longer than the input.

    Two farthest-point searches find the diameter instead, which is the
    standard linear method: the farthest node from any start is an end of some
    diameter, and the farthest node from *that* is the other end. Walking the
    path between them to its halfway point gives the branch and the offset,
    and root_on_branch splits that branch there.

    Returns the tree, rooted in place.
    """
    tips = tree.get_terminals()
    if len(tips) < 2:
        return tree

    parent = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parent[id(child)] = clade

    start, _, _ = _farthest_from(tips[0], parent)
    end, distance, previous = _farthest_from(start, parent)
    if distance <= 0:
        # Every branch is zero length; there is no midpoint to find.
        return tree

    # Walk back from `end` towards `start` until the halfway point is passed.
    half = distance / 2.0
    travelled = 0.0
    current = end
    while True:
        step_to = previous[id(current)]
        if parent.get(id(current)) is step_to:
            child, weight = current, current.branch_length or 0.0
            child_side = half - travelled
        else:
            child, weight = step_to, step_to.branch_length or 0.0
            child_side = weight - (half - travelled)
        if travelled + weight >= half:
            break
        travelled += weight
        current = step_to

    return root_on_branch(tree, child, child_side)


class Tree_Rerooter:

    def __init__(self, input_tree, output_tree, outgroup, allow_extra_outgroup_descendants, root_fraction, order_node):
        self.input_tree = input_tree
        self.output_tree = output_tree
        self.outgroup = self.split_outgroup(outgroup)
        self.allow_extra_outgroup_descendants = allow_extra_outgroup_descendants
        self.root_fraction = root_fraction
        self.order_node = order_node

    @staticmethod
    def split_outgroup(outgroup):
        """Accept space-separated, comma-separated, or both.

        Nextflow passes a parameter as one string, so `--outgroup A,B` has to
        work; the space-separated form is kept for direct command-line use.
        Order is preserved and duplicates are dropped, so a name repeated across
        both forms does not become a phantom extra tip.
        """
        names = []

        for chunk in outgroup or []:
            for name in str(chunk).split(","):
                name = name.strip()
                if name and name not in names:
                    names.append(name)

        return names

    def terminal_lookup(self, tree):
        lookup = {}
        duplicates = set()

        for tip in tree.get_terminals():
            if tip.name is None:
                raise ValueError("The tree contains an unnamed terminal.")

            if tip.name in lookup:
                duplicates.add(tip.name)
            else:
                lookup[tip.name] = tip

        if duplicates:
            raise ValueError(
                "Duplicate terminal names: " + ", ".join(sorted(duplicates))
            )

        return lookup

    def canonical_split(self, descendants, all_taxa):
        first = frozenset(descendants)
        second = frozenset(all_taxa - descendants)

        if len(first) < len(second):
            return first

        if len(second) < len(first):
            return second

        return min(first, second, key=lambda side: tuple(sorted(side)))

    @staticmethod
    def _tip_hash(name):
        return int.from_bytes(
            hashlib.blake2b(name.encode("utf-8"), digest_size=8).digest(), "big"
        )

    def split_keys(self, tree, all_taxa):
        """{id(node): (tip count, split key)} for every node, in one pass.

        Collecting each node's descendant set costs tips x depth, which on an
        UShER tree (138,095 tips, deep ladders) runs for hours. Instead each
        tip gets a fixed 64-bit hash and a node's side of the split is the XOR
        of its tips - the other side is the total XOR with that removed. The
        key is the smaller side's size and hash (both hashes on a tie), so the
        same bipartition gets the same key under any rooting, like
        canonical_split but linear.
        """
        total_hash = 0
        for name in all_taxa:
            total_hash ^= self._tip_hash(name)
        total = len(all_taxa)

        info = {}
        for clade in reversed(list(tree.find_clades(order="level"))):
            if not clade.clades:
                named = clade.name is not None
                count = 1 if named else 0
                value = self._tip_hash(clade.name) if named else 0
            else:
                count = 0
                value = 0
                for child in clade.clades:
                    child_count, child_value = info[id(child)][:2]
                    count += child_count
                    value ^= child_value
            info[id(clade)] = (count, value)

        keys = {}
        for node_id, (count, value) in info.items():
            other_count = total - count
            other_value = total_hash ^ value
            if count < other_count:
                key = (count, value)
            elif other_count < count:
                key = (other_count, other_value)
            else:
                key = (count, min(value, other_value))
            keys[node_id] = (count, key)
        return keys

    def get_supports(self, tree, all_taxa):
        supports = {}
        keys = self.split_keys(tree, all_taxa)

        for node in tree.get_nonterminals():
            if node is tree.root:
                continue

            count, split = keys[id(node)]

            if count <= 1:
                continue

            if count >= len(all_taxa) - 1:
                continue

            if split not in supports:
                supports[split] = node.confidence
            elif supports[split] is None and node.confidence is not None:
                supports[split] = node.confidence

        return supports

    def find_outgroup(self, tree):
        lookup = self.terminal_lookup(tree)
        missing = []

        for name in self.outgroup:
            if name not in lookup:
                missing.append(name)

        if missing:
            raise ValueError(
                "Outgroup tip(s) not found: " + ", ".join(missing)
            )

        outgroup_tips = []

        for name in self.outgroup:
            outgroup_tips.append(lookup[name])

        if len(outgroup_tips) == 1:
            outgroup_clade = outgroup_tips[0]
        else:
            outgroup_clade = tree.common_ancestor(outgroup_tips)

        descendants = set()

        for tip in outgroup_clade.get_terminals():
            if tip.name is not None:
                descendants.add(tip.name)

        extra_tips = descendants - set(self.outgroup)

        if extra_tips and not self.allow_extra_outgroup_descendants:
            raise ValueError(
                "The outgroup tips are not an exclusive clade. "
                "Their MRCA also contains: " + ", ".join(sorted(extra_tips))
            )

        if outgroup_clade is tree.root:
            raise ValueError(
                "The selected outgroup MRCA is the complete tree."
            )

        if len(self.outgroup) == 1:
            # Rooting on one tip puts the root on that tip's terminal branch, so
            # the position depends entirely on one sequence - its assembly, its
            # coverage, its long-branch attraction. If that tip has relatives in
            # the tree they end up in the ingroup, which is rarely what anyone
            # means by "root on the outgroup".
            #
            # Expanding to the sequence's MMseqs cluster does not rescue this for
            # a divergent genus: measured on HCV, all four genotype 8 references
            # are singleton clusters at 95% identity, so the cluster IS the tip.
            # Name the whole clade instead.
            sibling_names = sorted(
                tip.name for tip in tree.get_terminals()
                if tip.name is not None and tip.name != self.outgroup[0]
            )
            print(
                "NOTE: rooting on a single tip (" + self.outgroup[0] + "). The "
                "root lands on that one terminal branch, so its position rests "
                "on a single sequence, and any close relatives of it stay in the "
                "ingroup. Naming the whole outgroup clade is usually what is "
                "meant - e.g. all four HCV genotype 8 references rather than one "
                "of them. " + str(len(sibling_names)) + " other tips are present.",
                file=sys.stderr
            )

        return outgroup_clade

    def warn_if_internal_labels_will_be_lost(self, tree):
        """Internal node labels are dropped. That is correct, but say so.

        A support value is a property of a *bipartition*, so it survives a
        reroot - that is what get_supports/restore_supports do. An internal node
        *name* identifies a node in one particular rooting, and two nodes can
        collapse onto the same split (in a bifurcating tree both root children
        describe the same bipartition), so names cannot be carried across.

        Keeping them would be worse than dropping them. Measured on a 238-tip
        HCV MAT: UShER assigns node_N deterministically from the topology, and
        regenerates an identical set - 237/237 labels on identical clades -
        whether or not the input tree carried any, at identical parsimony
        (224,233 both ways). A preserved label would therefore describe the OLD
        rooting while a fresh UShER run on this tree assigns its own, giving two
        numbering schemes that disagree.

        For an IQ-TREE contree none of this applies: the internal label IS the
        support, already parsed into `confidence`, and there is nothing to lose.
        """
        labelled = 0

        for node in tree.get_nonterminals():
            if node.name and node.confidence is None:
                labelled += 1

        if labelled:
            print(
                "NOTE: dropping " + str(labelled) + " internal node label(s). "
                "They identify nodes in the input rooting and cannot carry over. "
                "UShER regenerates node_N deterministically from the topology, so "
                "re-feeding this tree gives a fresh self-consistent set; but node "
                "IDs recorded against the previous rooting (stored RIPPLES "
                "results, for example) will not line up with it.",
                file=sys.stderr
            )

    def restore_supports(self, tree, all_taxa, supports):
        keys = self.split_keys(tree, all_taxa)

        for node in tree.get_nonterminals():
            node.name = None

            if node is tree.root:
                node.confidence = None
                continue

            count, split = keys[id(node)]

            if count <= 1:
                node.confidence = None
                continue

            if count >= len(all_taxa) - 1:
                node.confidence = None
                continue

            node.confidence = supports.get(split)

        self.remove_duplicate_root_support(tree)

    def remove_duplicate_root_support(self, tree):
        if len(tree.root.clades) != 2:
            return

        first = tree.root.clades[0]
        second = tree.root.clades[1]

        first_size = len(first.get_terminals())
        second_size = len(second.get_terminals())

        if first_size <= second_size:
            second.confidence = None
        else:
            first.confidence = None

    def format_support(self, value):
        if value is None:
            return ""

        value = float(value)

        if value.is_integer():
            return str(int(value))

        return f"{value:.10f}".rstrip("0").rstrip(".")

    def format_branch_length(self, value):
        """Render a branch length without inventing or destroying precision.

        A fixed "%.10f" turns anything below ~5e-11 into 0.0000000000, and a
        zero-length branch is not the same tree - it merges two nodes as far as
        every downstream consumer is concerned. Repeated general format keeps
        short lengths short and falls back to scientific notation only when a
        value genuinely needs it.
        """
        value = float(value)

        if value == 0:
            return "0"

        for precision in range(10, 18):
            text = f"{value:.{precision}g}"
            if float(text) == value:
                return text

        return repr(value)

    def quote_name(self, name):
        if name is None:
            return ""

        special_characters = set("()[],:;' \t\r\n")

        for character in name:
            if character in special_characters:
                return "'" + name.replace("'", "''") + "'"

        return name

    def tree_to_newick(self, node, is_root=False):
        if node.is_terminal():
            text = self.quote_name(node.name)
        else:
            children = []

            for child in node.clades:
                children.append(self.tree_to_newick(child))

            text = "(" + ",".join(children) + ")"

            if node.confidence is not None:
                text += self.format_support(node.confidence)
            elif node.name:
                text += self.quote_name(node.name)

        if node.branch_length is not None and not is_root:
            text += ":" + self.format_branch_length(node.branch_length)

        return text


    def order_nodes(self, tree):
        """Ladderize by descendant tip count, as Tree.ladderize does.

        Tree.ladderize recounts every subtree at every node; counting once,
        children first, gives the same stable order in linear time.
        """
        if self.order_node not in ("increase", "decrease"):
            return

        counts = {}
        for clade in reversed(list(tree.find_clades(order="level"))):
            if clade.clades:
                counts[id(clade)] = sum(counts[id(child)] for child in clade.clades)
                clade.clades.sort(
                    key=lambda child: counts[id(child)],
                    reverse=self.order_node == "decrease"
                )
            else:
                counts[id(clade)] = 1

    def write_tree(self, tree):
        output_directory = os.path.dirname(self.output_tree)

        if output_directory:
            os.makedirs(output_directory, exist_ok=True)

        with open(self.output_tree, "w") as output:
            output.write(self.tree_to_newick(tree.root, is_root=True) + ";\n")

    def surviving_support_count(self, tree):
        count = 0

        for node in tree.get_nonterminals():
            if node is not tree.root and node.confidence is not None:
                count += 1

        return count

    def root(self, tree):
        """Root `tree` in place (outgroup if given, else midpoint) and return it.

        Supports are carried across by bipartition, internal node names are
        dropped, and nodes are ordered per `order_node`. No file I/O, so the
        database builder can root trees it already holds in memory.
        """
        if self.root_fraction < 0 or self.root_fraction > 1:
            raise ValueError("--root_fraction must be between 0 and 1.")

        # Bio.Phylo does not reject arbitrary text: "this is not newick" parses
        # as a single unnamed-structure tree with one tip. That used to be
        # caught only by accident, when root_at_midpoint fell over an unbound
        # variable on a one-tip tree; now that midpoint rooting handles the
        # degenerate case cleanly, the input has to be checked here instead.
        tip_count = len(tree.get_terminals())
        if tip_count < 2:
            raise ValueError(
                f"{self.input_tree or 'the input'} does not contain a usable "
                f"tree: parsed {tip_count} tip(s). Rerooting needs at least two."
            )

        terminal_names = self.terminal_lookup(tree)
        all_taxa = set(terminal_names.keys())
        supports = self.get_supports(tree, all_taxa)
        self.warn_if_internal_labels_will_be_lost(tree)

        if self.outgroup:
            outgroup_clade = self.find_outgroup(tree)
            branch_length = outgroup_clade.branch_length
            siblings = [
                child for child in tree.root.clades if child is not outgroup_clade
            ]

            if branch_length is not None and len(siblings) == 1 and len(tree.root.clades) == 2:
                # Under a two-child root the branch above the outgroup runs on
                # through the root into its sibling; that whole length is the
                # one --root_fraction is a fraction of.
                sibling_length = siblings[0].branch_length or 0.0
                branch_length += sibling_length
                root_on_branch(
                    tree, outgroup_clade, branch_length * self.root_fraction
                )
            elif branch_length is None:
                root_on_branch(tree, outgroup_clade, None)
            else:
                root_on_branch(
                    tree, outgroup_clade, branch_length * self.root_fraction
                )

            self.rooting_method = "outgroup: " + ", ".join(self.outgroup)
        else:
            missing = [
                node for node in tree.find_clades()
                if node is not tree.root and node.branch_length is None
            ]
            if missing:
                raise ValueError(
                    "Midpoint rooting needs branch lengths, and "
                    + str(len(missing))
                    + " branch(es) have none. Supply --outgroup to root a "
                    "cladogram, or use a tree with branch lengths."
                )

            if self.root_fraction != 0.5:
                print(
                    "WARNING: --root_fraction only applies to outgroup rooting; "
                    "it is ignored for midpoint rooting.",
                    file=sys.stderr
                )

            midpoint_root(tree)
            self.rooting_method = "midpoint"

        tree.rooted = True
        self.restore_supports(tree, all_taxa, supports)
        self.order_nodes(tree)
        self.all_taxa = all_taxa
        return tree

    def reroot_tree(self):
        original_tree = Phylo.read(self.input_tree, "newick")
        rerooted_tree = self.root(copy.deepcopy(original_tree))
        all_taxa = self.all_taxa
        self.write_tree(rerooted_tree)

        output_tree = Phylo.read(self.output_tree, "newick")
        output_taxa = set()

        for tip in output_tree.get_terminals():
            if tip.name is not None:
                output_taxa.add(tip.name)

        if output_taxa != all_taxa:
            raise RuntimeError(
                "The output tree does not contain the same tips as the input tree."
            )

        print("Rooting method: " + self.rooting_method)
        print("Input tips: " + str(len(all_taxa)))
        print(
            "Bootstrap labels retained: "
            + str(self.surviving_support_count(output_tree))
        )
        print("Node order: " + self.order_node)
        print("Output tree: " + self.output_tree)


def root_newick(newick, outgroup=None, order_node="increase"):
    """Root a newick string (or a parsed Bio.Phylo tree, rooted in place) and
    return the rooted newick string.

    Midpoint rooting unless `outgroup` names tips. Raises ValueError on a tree
    that cannot be rooted (fewer than two tips, missing branch lengths,
    duplicate tip names).
    """
    rerooter = Tree_Rerooter(None, None, outgroup or [], False, 0.5, order_node)
    if isinstance(newick, str):
        newick = Phylo.read(StringIO(newick), "newick")
    tree = rerooter.root(newick)
    return rerooter.tree_to_newick(tree.root, is_root=True) + ";"


def prune_to_tips(tree, keep):
    """Drop every tip whose name is not in `keep`, in place, in O(n).

    Biopython's ``Tree.prune`` walks from the root once per removed tip, which
    is quadratic when most of an 11,000-tip tree goes. Here each node is
    visited once, children first: a node left with no children is removed, and
    a node left with one child is spliced out, its branch length added to the
    child's so tip-to-tip distances are unchanged. The surviving child keeps
    its own support value, since the bipartition it describes is what remains.

    Returns the tree, or None when no tip survives.
    """
    order = list(tree.find_clades(order="level"))
    alive = {}
    replacement = {}

    for clade in reversed(order):
        if not clade.clades:
            alive[id(clade)] = clade.name in keep
            continue

        children = []
        for child in clade.clades:
            if not alive[id(child)]:
                continue
            children.append(replacement.get(id(child), child))
        clade.clades = children
        alive[id(clade)] = bool(children)

        if len(children) == 1 and clade is not tree.root:
            only = children[0]
            if clade.branch_length is not None or only.branch_length is not None:
                only.branch_length = (clade.branch_length or 0.0) + (only.branch_length or 0.0)
            replacement[id(clade)] = only

    if not alive[id(tree.root)]:
        return None

    while len(tree.root.clades) == 1:
        tree.root = tree.root.clades[0]
        tree.root.branch_length = None

    return tree


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Reroot a Newick or IQ-TREE contree using an outgroup. "
            "If no outgroup is provided, midpoint rooting is used."
        )
    )

    parser.add_argument(
        "-i",
        "--input_tree",
        required=True,
        help="Input Newick or IQ-TREE contree file"
    )

    parser.add_argument(
        "-o",
        "--output_tree",
        required=True,
        help="Output rerooted tree file"
    )

    parser.add_argument(
        "--outgroup",
        nargs="+",
        default=[],
        help=(
            "Optional outgroup tip name or names. Accepts either "
            "space-separated (--outgroup A B) or comma-separated "
            "(--outgroup A,B); the comma form is what a single Nextflow "
            "parameter string can carry."
        )
    )

    parser.add_argument(
        "--allow_extra_outgroup_descendants",
        action="store_true",
        help="Allow the outgroup MRCA to contain additional tips"
    )

    parser.add_argument(
        "--root_fraction",
        type=float,
        default=0.5,
        help="Root position on the outgroup branch (default: 0.5)"
    )

    parser.add_argument(
        "--order_node",
        choices=["increase", "decrease", "none"],
        default="increase",
        help=(
            "Order nodes by the number of descendant tips: increase, "
            "decrease, or none (default: increase)"
        )
    )

    args = parser.parse_args()

    rerooter = Tree_Rerooter(
        args.input_tree,
        args.output_tree,
        args.outgroup,
        args.allow_extra_outgroup_descendants,
        args.root_fraction,
        args.order_node
    )

    try:
        rerooter.reroot_tree()
    except Exception as error:
        print("ERROR: " + str(error), file=sys.stderr)
        raise SystemExit(1)
