#!/usr/bin/env python3

import argparse
import copy
import os
import sys

from Bio import Phylo

'''
python scripts/TreeReRoot.py --input_tree generic/rabv/tree/ref.treefile --output_tree ref_midpoint_rooted.treefile --order_node decrease
python scripts/TreeReRoot.py --input_tree generic/rabv/tree/ref_plus_outgroup.treefile --output_tree ref_plus_outgroup_rerooted.treefile --outgroup NC_009528.2 NC_009527.1 --order_node decrease 

'''
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

    def get_supports(self, tree, all_taxa):
        supports = {}

        for node in tree.get_nonterminals():
            if node is tree.root:
                continue

            descendants = set()

            for tip in node.get_terminals():
                if tip.name is not None:
                    descendants.add(tip.name)

            if len(descendants) <= 1:
                continue

            if len(descendants) >= len(all_taxa) - 1:
                continue

            split = self.canonical_split(descendants, all_taxa)

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
        for node in tree.get_nonterminals():
            node.name = None

            if node is tree.root:
                node.confidence = None
                continue

            descendants = set()

            for tip in node.get_terminals():
                if tip.name is not None:
                    descendants.add(tip.name)

            if len(descendants) <= 1:
                node.confidence = None
                continue

            if len(descendants) >= len(all_taxa) - 1:
                node.confidence = None
                continue

            split = self.canonical_split(descendants, all_taxa)
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
        if self.order_node == "increase":
            tree.ladderize(reverse=False)
        elif self.order_node == "decrease":
            tree.ladderize(reverse=True)

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

    def reroot_tree(self):
        if self.root_fraction < 0 or self.root_fraction > 1:
            raise ValueError("--root_fraction must be between 0 and 1.")

        original_tree = Phylo.read(self.input_tree, "newick")
        terminal_names = self.terminal_lookup(original_tree)
        all_taxa = set(terminal_names.keys())
        supports = self.get_supports(original_tree, all_taxa)
        self.warn_if_internal_labels_will_be_lost(original_tree)
        rerooted_tree = copy.deepcopy(original_tree)

        if self.outgroup:
            outgroup_clade = self.find_outgroup(rerooted_tree)
            branch_length = outgroup_clade.branch_length

            if branch_length is None:
                rerooted_tree.root_with_outgroup(outgroup_clade)
            else:
                rerooted_tree.root_with_outgroup(
                    outgroup_clade,
                    outgroup_branch_length=branch_length * self.root_fraction
                )

            rooting_method = "outgroup: " + ", ".join(self.outgroup)
        else:
            missing = [
                node for node in rerooted_tree.find_clades()
                if node is not rerooted_tree.root and node.branch_length is None
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

            rerooted_tree.root_at_midpoint()
            rooting_method = "midpoint"

        rerooted_tree.rooted = True
        self.restore_supports(rerooted_tree, all_taxa, supports)
        self.order_nodes(rerooted_tree)
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

        print("Rooting method: " + rooting_method)
        print("Input tips: " + str(len(all_taxa)))
        print(
            "Bootstrap labels retained: "
            + str(self.surviving_support_count(output_tree))
        )
        print("Node order: " + self.order_node)
        print("Output tree: " + self.output_tree)


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
