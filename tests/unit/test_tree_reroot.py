"""Tests for scripts/TreeReRoot.py.

The script reroots a tree by outgroup or, with no outgroup, at the midpoint. The
hard part is not the rerooting - Bio.Phylo does that - it is keeping bootstrap
support attached to the right branch afterwards. Rerooting moves branches around,
and a support value belongs to a *bipartition*, not to a node, so the script
records every split's support before rerooting and reattaches it after. These
tests pin that behaviour and the edge cases where it breaks.

Flaws found on the version imported from commit 0ca83df are marked `xfail` with
what they cost; each is fixed in the same change that adds these tests, so an
xpass here means a fix landed.
"""

import subprocess
import sys
from pathlib import Path

import pytest
from Bio import Phylo


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "TreeReRoot.py"

sys.path.insert(0, str(REPO_ROOT / "scripts"))


def run(tmp_path, newick, *args, name="in.nwk"):
    """Run the script on a newick string; return (returncode, stdout+stderr, out_text)."""
    src = tmp_path / name
    src.write_text(newick)
    dst = tmp_path / "out.nwk"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "-i", str(src), "-o", str(dst), *args],
        capture_output=True, text=True, timeout=300,
    )
    out = dst.read_text() if dst.exists() else ""
    return proc.returncode, proc.stdout + proc.stderr, out


def tips(newick_text):
    from io import StringIO
    return {t.name for t in Phylo.read(StringIO(newick_text), "newick").get_terminals()}


def supports(newick_text):
    """{canonical split: confidence} for informative splits.

    Mirrors the script's own first-wins rule so the helper cannot disagree with
    it about which value a split carries.
    """
    from io import StringIO
    tree = Phylo.read(StringIO(newick_text), "newick")
    allt = {t.name for t in tree.get_terminals()}
    out = {}
    for node in tree.get_nonterminals():
        desc = frozenset(t.name for t in node.get_terminals())
        if not (1 < len(desc) < len(allt) - 0):
            continue
        if len(desc) >= len(allt) - 1:
            continue
        side = min(desc, frozenset(allt - desc), key=lambda s: (len(s), sorted(s)))
        if node.confidence is None:
            continue
        out.setdefault(side, node.confidence)
    return out


# ---------------------------------------------------------------------------
# Core behaviour
# ---------------------------------------------------------------------------

#: An UNROOTED tree with a trifurcating root, which is what IQ-TREE emits and
#: what this script is for. A bifurcating root would make both root children
#: describe the same bipartition, so a support written on each would be
#: contradictory - real tools do not do that, and neither should a fixture.
BALANCED = "(A:1,B:1,((C:1,D:1)88:1,(E:1,F:1)72:1)60:2);"


class TestMidpointIsTheDefault:
    def test_no_outgroup_means_midpoint(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED)
        assert rc == 0, log
        assert "Rooting method: midpoint" in log

    def test_every_tip_survives(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED)
        assert tips(out) == {"A", "B", "C", "D", "E", "F"}

    def test_the_output_is_rooted_at_a_bifurcation(self, tmp_path):
        from io import StringIO
        rc, log, out = run(tmp_path, BALANCED)
        tree = Phylo.read(StringIO(out), "newick")
        assert len(tree.root.clades) == 2, "a rooted tree should have two root children"

    def test_it_refuses_an_unreadable_tree(self, tmp_path):
        rc, log, out = run(tmp_path, "this is not newick")
        assert rc != 0


class TestSupportSurvivesRerooting:
    """The whole point of the script."""

    def test_supports_move_with_their_split(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED)
        before, after = supports(BALANCED), supports(out)
        for split, value in after.items():
            assert split in before, f"invented a support for {sorted(split)}"
            assert before[split] == value, f"support changed on {sorted(split)}"

    def test_no_support_is_silently_dropped_except_the_root_split(self, tmp_path):
        """Rerooting collapses one bipartition into the root, so exactly one
        input support legitimately disappears - never more."""
        rc, log, out = run(tmp_path, BALANCED)
        lost = set(supports(BALANCED)) - set(supports(out))
        assert len(lost) <= 1, f"lost {len(lost)} supports: {[sorted(s) for s in lost]}"

    def test_the_root_bipartition_is_not_written_twice(self, tmp_path):
        """Both root children describe the same split; labelling both would
        double-count it in any downstream parse."""
        from io import StringIO
        rc, log, out = run(tmp_path, BALANCED)
        tree = Phylo.read(StringIO(out), "newick")
        labelled = [c for c in tree.root.clades if c.confidence is not None]
        assert len(labelled) <= 1

    def test_it_reports_how_many_it_kept(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED)
        assert "Bootstrap labels retained:" in log


class TestOutgroupRooting:
    def test_single_tip_outgroup(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "A")
        assert rc == 0, log
        assert "Rooting method: outgroup: A" in log
        assert tips(out) == {"A", "B", "C", "D", "E", "F"}

    def test_multi_tip_monophyletic_outgroup(self, tmp_path):
        """C and D form a real clade in the fixture; A and B are root children,
        so their MRCA is the whole tree and is correctly refused."""
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "C", "D")
        assert rc == 0, log
        assert tips(out) == {"A", "B", "C", "D", "E", "F"}

    def test_an_outgroup_whose_mrca_is_the_root_is_refused(self, tmp_path):
        """On a trifurcating (unrooted) tree the exclusive-clade check fires
        first, so the message names the intruders rather than saying "complete
        tree". Either refusal is correct; what matters is that it does not
        silently root on everything."""
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "A", "B")
        assert rc != 0
        assert ("complete tree" in log) or ("not an exclusive clade" in log)

    def test_a_missing_outgroup_tip_is_named(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "NOPE")
        assert rc != 0
        assert "NOPE" in log, "the error must name the tip that was not found"

    def test_a_non_monophyletic_outgroup_is_refused(self, tmp_path):
        """Rooting on a paraphyletic group silently produces a different tree
        than intended, so it must be an error by default."""
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "C", "E")
        assert rc != 0
        assert "not an exclusive clade" in log

    def test_the_refusal_lists_the_intruders(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "C", "E")
        assert "D" in log and "F" in log

    def test_outgroup_of_everything_is_refused(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup",
                           "A", "B", "C", "D", "E", "F")
        assert rc != 0

    def test_root_fraction_out_of_range_is_refused(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "A",
                           "--root_fraction", "1.5")
        assert rc != 0
        assert "root_fraction" in log

    @pytest.mark.parametrize("fraction", ["0", "0.25", "0.5", "1"])
    def test_root_fraction_endpoints_are_accepted(self, tmp_path, fraction):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "A",
                           "--root_fraction", fraction)
        assert rc == 0, log


class TestNodeOrdering:
    @pytest.mark.parametrize("order", ["increase", "decrease", "none"])
    def test_orders_are_accepted_and_preserve_tips(self, tmp_path, order):
        rc, log, out = run(tmp_path, BALANCED, "--order_node", order)
        assert rc == 0, log
        assert tips(out) == {"A", "B", "C", "D", "E", "F"}

    def test_an_unknown_order_is_refused(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--order_node", "sideways")
        assert rc != 0


class TestNameQuoting:
    def test_names_with_spaces_round_trip(self, tmp_path):
        rc, log, out = run(tmp_path, "(('tip one':1,'tip two':1):2,C:1);")
        assert rc == 0, log
        assert tips(out) == {"tip one", "tip two", "C"}

    def test_duplicate_tip_names_are_refused(self, tmp_path):
        """A duplicate name makes every split ambiguous."""
        rc, log, out = run(tmp_path, "((A:1,A:1):2,(C:1,D:1):2);")
        assert rc != 0
        assert "Duplicate" in log


# ---------------------------------------------------------------------------
# Flaws found in the imported version
# ---------------------------------------------------------------------------

class TestFlawsFoundOnImport:
    def test_internal_node_labels_are_dropped_but_not_silently(self, tmp_path):
        """Internal labels cannot survive a reroot, so the script warns.

        A support value is a property of a *bipartition* and carries across a
        reroot - that is what get_supports/restore_supports do. An internal node
        *name* identifies a node in one particular rooting; re-hanging the tree
        makes it meaningless, and two nodes can collapse onto the same split (in
        a bifurcating tree both root children describe the same bipartition, so
        preserving names by split would duplicate one and lose the other).

        Dropping them is therefore correct. Dropping them silently is not: for an
        UShER tree those `node_N` labels are how the trees table and RIPPLES
        refer to internal nodes.
        """
        rc, log, out = run(tmp_path, "((A:1,B:1)node_1:2,(C:1,D:1)node_2:2);")
        assert rc == 0, log
        assert "node_1" not in out and "node_2" not in out
        assert "NOTE" in log and "internal node label" in log
        assert "node_N" in log and "regenerates" in log, (
            "the note should say UShER regenerates them, not imply they are lost")

    def test_no_warning_when_there_are_no_internal_labels_to_lose(self, tmp_path):
        """An IQ-TREE contree carries support, not names - nothing is lost, so
        warning would be noise."""
        rc, log, out = run(tmp_path, BALANCED)
        assert "internal node label" not in log

    def test_small_branch_lengths_are_not_flattened_to_zero(self, tmp_path):
        """Lengths are written with %.10f, so anything below ~5e-11 becomes
        0.0000000000. A zero-length branch is not the same tree."""
        rc, log, out = run(tmp_path, "((A:0.00000000001,B:0.000000000012):2,(C:1,D:1):2);")
        assert rc == 0, log
        assert "0.0000000000" not in out, (
            "a sub-1e-10 branch length was rounded to zero")

    def test_a_tree_without_branch_lengths_fails_clearly(self, tmp_path):
        """Midpoint rooting needs lengths. Bio.Phylo raises
        `local variable 'tip1' referenced before assignment`, which tells the
        user nothing about what is wrong with their tree."""
        rc, log, out = run(tmp_path, "((A,B),(C,D));")
        assert rc != 0, "a cladogram cannot be midpoint rooted"
        assert "tip1" not in log, f"opaque Bio.Phylo internal leaked: {log.strip()[:120]}"
        assert "branch length" in log.lower()

    def test_comma_separated_outgroup_is_accepted(self, tmp_path):
        """Space separation alone is awkward to pass through Nextflow, where a
        single parameter string is the natural shape."""
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "C,D")
        assert rc == 0, log
        assert tips(out) == {"A", "B", "C", "D", "E", "F"}

    def test_root_fraction_with_midpoint_is_not_silently_ignored(self, tmp_path):
        """root_fraction only affects outgroup rooting. Accepting it silently
        for a midpoint run lets someone believe it did something."""
        rc, log, out = run(tmp_path, BALANCED, "--root_fraction", "0.9")
        assert "root_fraction" in log.lower(), (
            "midpoint rooting ignored --root_fraction without saying so")


# ---------------------------------------------------------------------------
# Round-trip guarantees that must hold for any input
# ---------------------------------------------------------------------------

TREES = [
    pytest.param(BALANCED, id="balanced-with-support"),
    pytest.param("((A:1,B:1):2,(C:1,D:1):2);", id="no-support"),
    pytest.param("(A:1,(B:1,(C:1,(D:1,E:1):1):1):1);", id="ladder"),
    pytest.param("((A:0.1,B:0.2)0.95:0.3,(C:0.4,D:0.5)0.80:0.6);", id="decimal-support"),
    pytest.param("((A:1,B:1,C:1):2,(D:1,E:1):2);", id="polytomy"),
]


@pytest.mark.parametrize("newick", TREES)
class TestRoundTrip:
    def test_tip_set_is_identical(self, tmp_path, newick):
        rc, log, out = run(tmp_path, newick)
        assert rc == 0, log
        assert tips(out) == tips(newick)

    def test_output_is_parseable(self, tmp_path, newick):
        from io import StringIO
        rc, log, out = run(tmp_path, newick)
        Phylo.read(StringIO(out), "newick")

    def test_total_branch_length_is_conserved(self, tmp_path, newick):
        """Rerooting moves the root along a branch; it never adds or removes
        evolutionary distance."""
        from io import StringIO
        rc, log, out = run(tmp_path, newick)
        def total(text):
            t = Phylo.read(StringIO(text), "newick")
            return sum(n.branch_length or 0 for n in t.find_clades())
        assert total(out) == pytest.approx(total(newick), abs=1e-6)


# ---------------------------------------------------------------------------
# Pipeline wiring: which profiles root how
# ---------------------------------------------------------------------------

CONFIG = REPO_ROOT / "nextflow.config"
HCV_REF_LIST = REPO_ROOT / "generic" / "hcv" / "ref_list_subtype_genotype.txt"
GENOTYPE_8 = ["MH590698", "MH590699", "MH590700", "MH590701"]


def profile_blocks():
    """{profile name: its config text}."""
    import re
    text = CONFIG.read_text()
    starts = [(m.group(1), m.start()) for m in re.finditer(r'^\s{4}(\w+)\s*\{', text, re.M)]
    out = {}
    for i, (name, pos) in enumerate(starts):
        end = starts[i + 1][1] if i + 1 < len(starts) else len(text)
        out[name] = text[pos:end]
    return out


class TestRootingDefaultsToMidpoint:
    def test_the_top_level_default_is_null(self):
        """Midpoint needs no prior knowledge, so it is the safe default for a
        species nobody has chosen an outgroup for."""
        import re
        text = CONFIG.read_text()
        block = text[re.search(r'^params\s*\{', text, re.M).start():]
        assert re.search(r'^\s+root_outgroup\s*=\s*null', block, re.M)

    def test_no_non_hcv_profile_sets_an_outgroup(self):
        """Notably the segmented profiles: a single outgroup string cannot
        express a per-segment outgroup, so influenza must midpoint root."""
        offenders = [n for n, b in profile_blocks().items()
                     if "root_outgroup" in b and not n.upper().startswith("HCV")]
        assert not offenders, f"non-HCV profiles setting an outgroup: {offenders}"


class TestEverythingMidpointRootsForNow:
    """No profile sets an outgroup, HCV included.

    Genotype 8 was evaluated as the HCV outgroup and rejected. It is the most
    divergent HCV lineage and monophyletic in the reference set (4/4), so it
    works mechanically - measured, it splits the tree 4 | 234 where midpoint
    splits it 43 | 195 on the genotype 2+7 branch. But those four sequences are
    all Canadian, collected 2015-2017, from a lineage first described in 2018.
    Anchoring the whole HCV tree on a group that narrow and that recently
    characterised is weaker than not choosing at all, so midpoint stands.

    The mechanism is kept for a species that has a defensible outgroup. These
    tests exist so that reintroducing one is a deliberate act with a reason
    attached, not a quiet edit.
    """

    @pytest.mark.parametrize("profile", ["HCV_test", "HCV_full", "HCV_XML_full", "HCV_update"])
    def test_hcv_profiles_do_not_set_an_outgroup(self, profile):
        blocks = profile_blocks()
        assert profile in blocks, f"{profile} missing from nextflow.config"
        assert "root_outgroup" not in blocks[profile], (
            f"{profile} sets an outgroup. Midpoint is the current decision for "
            f"HCV; if that changed, update this test and say why in "
            f"info_help/tree_rooting.md")

    def test_no_profile_at_all_sets_one(self):
        offenders = [n for n, b in profile_blocks().items() if "root_outgroup" in b]
        assert not offenders, f"profiles setting an outgroup: {offenders}"

    def test_the_mechanism_is_still_available(self):
        """Removing the HCV outgroup must not remove the ability to set one."""
        import re
        text = CONFIG.read_text()
        block = text[re.search(r'^params\s*\{', text, re.M).start():]
        assert re.search(r'^\s+root_outgroup\s*=\s*null', block, re.M),             "the root_outgroup parameter itself should remain declared"


class TestParamIsDeclared:
    def test_root_outgroup_is_in_the_allowlist(self):
        """vgtk-init.nf rejects params it does not know about."""
        text = (REPO_ROOT / "vgtk-init.nf").read_text()
        block = text[text.index("def scriptDefinedParams"):]
        block = block[:block.index("]")]
        assert '"root_outgroup"' in block or "'root_outgroup'" in block


class TestOutgroupSplitting:
    """The comma form is what a single Nextflow parameter string can carry."""

    @pytest.mark.parametrize("given,expected", [
        (["A,B"], ["A", "B"]),
        (["A", "B"], ["A", "B"]),
        (["A,B", "C"], ["A", "B", "C"]),
        (["A, B , C"], ["A", "B", "C"]),
        (["A,,B"], ["A", "B"]),
        ([""], []),
        ([], []),
        (None, []),
    ])
    def test_forms(self, given, expected):
        from TreeReRoot import Tree_Rerooter
        assert Tree_Rerooter.split_outgroup(given) == expected

    def test_duplicates_are_dropped_and_order_kept(self):
        """A name repeated across both forms must not become a phantom extra tip
        in the exclusive-clade check."""
        from TreeReRoot import Tree_Rerooter
        assert Tree_Rerooter.split_outgroup(["B,A", "A"]) == ["B", "A"]


class TestSingleTipRootingIsFlagged:
    """Rooting on one sequence puts the root on that tip's terminal branch.

    Measured on HCV: rooting on MH590700 alone gives sides of 1 and 237 tips,
    with the other three genotype 8 references left in the ingroup. Expanding to
    the MMseqs cluster does not help - all four genotype 8 references are
    singleton clusters at 95% identity, because HCV genotypes are ~30% divergent.
    The fix is to name the whole clade, so the script says so.
    """

    def test_a_single_tip_outgroup_is_noted(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "C")
        assert rc == 0, log
        assert "single tip" in log

    def test_the_note_suggests_the_clade(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "C")
        assert "clade" in log

    def test_a_multi_tip_outgroup_is_not_noted(self, tmp_path):
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "C,D")
        assert "single tip" not in log

    def test_it_is_a_note_not_a_refusal(self, tmp_path):
        """Sometimes one tip really is the whole outgroup; do not block it."""
        rc, log, out = run(tmp_path, BALANCED, "--outgroup", "C")
        assert rc == 0
        assert tips(out) == {"A", "B", "C", "D", "E", "F"}


# ---------------------------------------------------------------------------
# Linear-time midpoint rooting
#
# Biopython's Tree.root_at_midpoint re-roots the tree at every tip in turn and
# recomputes all depths each time. Measured on this repository's HCV trees:
# 250 tips 0.07s, 500 0.28s, 1000 1.21s, 2000 6.16s - a clean 4x per doubling.
# The 138,095-tip UShER tree extrapolates to about eight hours, and midpoint
# rooting is now part of the pipeline, so that is a per-run cost on every real
# dataset.
#
# midpoint_root finds the diameter with two farthest-point searches instead.
# The placement arithmetic afterwards is deliberately Biopython's own: its
# `root_with_outgroup` does not split a branch around the new root, and the
# `max_distance` its formula expects is measured AFTER rooting at the first
# tip, which is larger than the true diameter. Both of those are easy to get
# wrong in ways that still produce a plausible-looking tree, so the contract
# asserted here is equality with Biopython, not merely "a midpoint".
# ---------------------------------------------------------------------------

import copy as _copy
import random as _random

from Bio.Phylo.BaseTree import Clade as _Clade, Tree as _Tree

import TreeReRoot as _TRR


def _random_tree(n, seed, zero_length_fraction=0.0):
    rng = _random.Random(seed)

    def length():
        if zero_length_fraction and rng.random() < zero_length_fraction:
            return 0.0
        return rng.uniform(0.001, 2.0)

    clades = [_Clade(branch_length=length(), name=f"t{i}") for i in range(n)]
    while len(clades) > 1:
        a = clades.pop(rng.randrange(len(clades)))
        b = clades.pop(rng.randrange(len(clades)))
        clades.append(_Clade(branch_length=length(), clades=[a, b]))
    return _Tree(root=clades[0], rooted=False)


def _tip_depths(tree):
    return {t.name: round(tree.distance(tree.root, t), 9) for t in tree.get_terminals()}


def _pairwise(tree):
    lookup = {t.name: t for t in tree.get_terminals()}
    names = sorted(lookup)
    return {
        (a, b): round(tree.distance(lookup[a], lookup[b]), 9)
        for i, a in enumerate(names) for b in names[i + 1:]
    }


class TestMidpointRootIsExact:
    """Rooting may move the root, never change a tip-to-tip distance.

    Biopython's root_at_midpoint fails this on small trees: it roots
    (A:0.1,(B:0.2,C:0.3):0.4) as (A:0.5,(B,C):0.1) and inflates a two-tip tree
    to 1.5x its length, because root_with_outgroup does not split a branch that
    hangs directly off the root. Before the fix, 45 of 180 random 2-5 tip trees
    came out distorted.
    """

    @pytest.mark.parametrize("n", [2, 3, 4, 5, 8])
    @pytest.mark.parametrize("seed", range(12))
    def test_distances_kept_and_root_at_half_diameter(self, n, seed):
        original = _random_tree(n, seed)
        parent = {}
        for clade in original.find_clades(order="level"):
            for child in clade.clades:
                parent[id(child)] = clade
        end, _, _ = _TRR._farthest_from(original.get_terminals()[0], parent)
        _, diameter, _ = _TRR._farthest_from(end, parent)

        rooted = _TRR.midpoint_root(_copy.deepcopy(original))
        assert _pairwise(rooted) == _pairwise(original)
        deepest = max(rooted.distance(rooted.root, t) for t in rooted.get_terminals())
        assert deepest == pytest.approx(diameter / 2.0, abs=1e-9)

    def test_known_small_tree(self):
        from io import StringIO
        tree = Phylo.read(StringIO("(A:0.1,(B:0.2,C:0.3):0.4);"), "newick")
        _TRR.midpoint_root(tree)
        depths = _tip_depths(tree)
        assert depths == {"A": 0.4, "B": 0.3, "C": 0.4}


class TestRootOnBranch:
    @pytest.mark.parametrize("newick,path,offset", [
        ("(A:0.1,B:0.2,C:0.3);", [0], 0.05),
        ("(A:0.5,(B:0.2,C:0.3):0.4,D:1);", [1], 0.1),
        ("(A:0.5,((B:0.2,E:1):0.7,C:0.3):0.4,D:1);", [1, 0], 0.3),
        ("((A:0.1,B:0.2):0.3,(C:0.1,D:0.2):0.4);", [0], 0.5),
    ])
    def test_splits_the_branch_without_adding_length(self, newick, path, offset):
        from io import StringIO
        tree = Phylo.read(StringIO(newick), "newick")
        before = _pairwise(tree)
        node = tree.root
        for index in path:
            node = node.clades[index]
        _TRR.root_on_branch(tree, node, offset)
        assert _pairwise(tree) == before
        assert node in tree.root.clades
        assert node.branch_length == pytest.approx(offset)


class TestPruneToTips:
    def test_keeps_distances_between_survivors(self):
        tree = _random_tree(200, 77)
        keep = {f"t{i}" for i in range(0, 200, 7)}
        expected = {k: v for k, v in _pairwise(tree).items() if k[0] in keep and k[1] in keep}
        pruned = _TRR.prune_to_tips(tree, keep)
        assert {t.name for t in pruned.get_terminals()} == keep
        assert _pairwise(pruned) == expected
        assert all(len(c.clades) != 1 for c in pruned.find_clades())

    def test_nothing_kept_returns_none(self):
        assert _TRR.prune_to_tips(_random_tree(10, 1), set()) is None

    def test_root_newick_accepts_a_pruned_tree(self):
        pruned = _TRR.prune_to_tips(_random_tree(50, 3), {"t1", "t2", "t3"})
        assert tips(_TRR.root_newick(pruned)) == {"t1", "t2", "t3"}


class TestMidpointRootMatchesBiopython:
    # Biopython is only a trustworthy reference where it does not hit the
    # small-tree bug above.
    @pytest.mark.parametrize("n,seed", [
        (8, 1), (25, 2), (60, 3), (150, 4), (400, 5), (900, 6),
    ])
    def test_identical_rooting(self, n, seed):
        original = _random_tree(n, seed)
        reference = _copy.deepcopy(original)
        reference.root_at_midpoint()
        fast = _copy.deepcopy(original)
        _TRR.midpoint_root(fast)
        assert _tip_depths(fast) == _tip_depths(reference)

    @pytest.mark.parametrize("n,seed", [(30, 11), (120, 12), (500, 13)])
    def test_identical_with_zero_length_branches(self, n, seed):
        # Zero-length branches let an internal node tie a leaf for "farthest",
        # which would put the root off the tip-to-tip path.
        original = _random_tree(n, seed, zero_length_fraction=0.3)
        reference = _copy.deepcopy(original)
        reference.root_at_midpoint()
        fast = _copy.deepcopy(original)
        _TRR.midpoint_root(fast)
        assert _tip_depths(fast) == _tip_depths(reference)

    def test_root_sits_at_half_the_diameter(self):
        """The defining property, asserted independently of Biopython."""
        tree = _random_tree(300, 21)
        parent = {}
        for clade in tree.find_clades(order="level"):
            for child in clade.clades:
                parent[id(child)] = clade
        end, _, _ = _TRR._farthest_from(tree.get_terminals()[0], parent)
        _, diameter, _ = _TRR._farthest_from(end, parent)

        _TRR.midpoint_root(tree)
        deepest = max(tree.distance(tree.root, t) for t in tree.get_terminals())
        assert deepest == pytest.approx(diameter / 2.0, abs=1e-9)

    def test_single_tip_is_left_alone(self):
        tree = _Tree(root=_Clade(branch_length=1.0, name="only"))
        assert _TRR.midpoint_root(tree) is tree

    def test_all_zero_lengths_is_left_alone(self):
        # There is no midpoint to find, and re-rooting would move the root
        # arbitrarily while looking like it had done something meaningful.
        tree = _random_tree(20, 31)
        for clade in tree.find_clades():
            clade.branch_length = 0.0
        before = _tip_depths(tree)
        _TRR.midpoint_root(tree)
        assert _tip_depths(tree) == before


class TestFarthestFrom:
    def test_returns_a_tip_not_an_internal_node(self):
        tree = _random_tree(50, 41, zero_length_fraction=0.5)
        parent = {}
        for clade in tree.find_clades(order="level"):
            for child in clade.clades:
                parent[id(child)] = clade
        found, _, _ = _TRR._farthest_from(tree.get_terminals()[0], parent)
        assert not found.clades, "a diameter endpoint must be a tip"

    def test_finds_the_true_diameter_on_a_known_tree(self):
        # ((A:1,B:1):1,C:5)  ->  diameter is B..C (or A..C) = 1 + 1 + 5 = 7
        a = _Clade(branch_length=1.0, name="A")
        b = _Clade(branch_length=1.0, name="B")
        ab = _Clade(branch_length=1.0, clades=[a, b])
        c = _Clade(branch_length=5.0, name="C")
        tree = _Tree(root=_Clade(branch_length=0.0, clades=[ab, c]), rooted=False)
        parent = {}
        for clade in tree.find_clades(order="level"):
            for child in clade.clades:
                parent[id(child)] = clade
        end, _, _ = _TRR._farthest_from(tree.get_terminals()[0], parent)
        _, diameter, _ = _TRR._farthest_from(end, parent)
        assert diameter == pytest.approx(7.0)
