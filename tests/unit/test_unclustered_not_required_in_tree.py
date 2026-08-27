"""An unclustered accession cannot be in a cluster-derived tree.

The UShER tree is built from MMseqs cluster output. An accession that was never
assigned a cluster was never offered to clustering, so it can never appear as a
tree tip. Requiring it there is not a stricter check - it is a wrong one.

This is not hypothetical. `TEST_SUBSAMPLE_CLUSTER_INPUT` caps the clustering
input at `params.test_max_cluster_seqs` (120 for the HCV_OM_test profile), so a
test database has hundreds of accessions that were never clustered. Before this
fix every subsampled test run failed with "UShER tree does not match accessions"
- 218 of 338 HCV test accessions reported missing, every one of them unclustered.

Measured on two databases:

* `HCV_OM_test.db` - 338 accessions, 120 tree tips. All 218 absent from the tree
  were unclustered; all 116 clustered accessions were present.
* `HCV_full.db` - 274,606 accessions, 135,544 tree tips. All 139,062 absent were
  unclustered; **zero** clustered accessions were missing.

So the rule holds in production, not just under test subsampling.
"""

import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

import ValidateDbTree as V


REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATOR = REPO_ROOT / "scripts" / "ValidateDbTree.py"
HCV_TEST_DB = REPO_ROOT / "test_out" / "HCV_OM_test" / "HCV_OM_test.db"


class TestIsClusterPlaceholder:
    """The predicate the exclusion rests on."""

    @pytest.mark.parametrize("value", [
        None, "", "   ",
        "na-see-tree", "na- see tree", "na - see tree", "NA-SEE-TREE",
    ])
    def test_placeholders(self, value):
        """SQL NULL, blank, and the pipeline's own 'na-see-tree' sentinel."""
        assert V.is_cluster_placeholder(value)

    @pytest.mark.parametrize("value", ["AB123456", "PZ746522", "  AB123456  "])
    def test_real_cluster_ids_are_not_placeholders(self, value):
        assert not V.is_cluster_placeholder(value)

    @pytest.mark.parametrize("value", ["NA", "na", "N/A", "None", "null", "-"])
    def test_bare_na_is_deliberately_not_a_placeholder(self, value):
        """`NA` is a real value in this codebase - influenza's neuraminidase
        segment is literally named NA. Treating it as null here would silently
        exclude real clusters from the tree check. Only the explicit
        'na-see-tree' sentinel counts. If this ever starts failing, check that
        the change does not resurrect the NA-is-null bug."""
        assert not V.is_cluster_placeholder(value)


# ---------------------------------------------------------------------------
# The rule itself, expressed as the set arithmetic the validator performs
# ---------------------------------------------------------------------------

def _missing_in_tree(meta_to_cluster, tree_tips):
    """Reproduce the validator's rule: expected - tips - unclustered."""
    expected = set(meta_to_cluster)
    unclustered = {
        acc for acc, cluster in meta_to_cluster.items()
        if V.is_cluster_placeholder(cluster)
    }
    return sorted(expected - set(tree_tips) - unclustered)


class TestTheRule:
    def test_unclustered_absent_from_the_tree_is_not_a_failure(self):
        meta = {"A": "A", "B": "A", "C": None, "D": ""}
        assert _missing_in_tree(meta, {"A", "B"}) == []

    def test_a_clustered_accession_missing_from_the_tree_still_fails(self):
        """The check must not be weakened for accessions that *could* be there."""
        meta = {"A": "A", "B": "A", "C": None}
        assert _missing_in_tree(meta, {"A"}) == ["B"]

    def test_unclustered_present_in_the_tree_is_fine(self):
        """UShER places references and master without a cluster assignment."""
        meta = {"A": "A", "REF": None}
        assert _missing_in_tree(meta, {"A", "REF"}) == []

    def test_everything_unclustered_and_no_tips(self):
        meta = {"A": None, "B": None}
        assert _missing_in_tree(meta, set()) == []

    def test_nothing_unclustered_behaves_as_before(self):
        """With every accession clustered the rule is the original one."""
        meta = {"A": "A", "B": "A", "C": "C"}
        assert _missing_in_tree(meta, {"A"}) == ["B", "C"]


# ---------------------------------------------------------------------------
# The real databases this was measured on
# ---------------------------------------------------------------------------

def _cluster_column(conn):
    columns = [row[1] for row in conn.execute("PRAGMA table_info(meta_data)")]
    return V.find_cluster_column(columns)


@pytest.mark.skipif(not HCV_TEST_DB.exists(), reason=f"not built: {HCV_TEST_DB}")
class TestAgainstTheHcvTestDatabase:
    @pytest.fixture
    def conn(self):
        c = sqlite3.connect(f"file:{HCV_TEST_DB}?mode=ro", uri=True)
        yield c
        c.close()

    def test_the_database_really_does_have_unclustered_accessions(self, conn):
        """If this stops being true the profile's subsampling changed, and the
        rest of this class is no longer testing what it claims to."""
        column = _cluster_column(conn)
        assert column, "no cluster column in the test DB"
        rows = conn.execute(
            f"SELECT primary_accession, {column} FROM meta_data").fetchall()
        unclustered = [a for a, c in rows if V.is_cluster_placeholder(c)]
        assert len(unclustered) > 100, (
            f"expected the test profile's subsampling to leave many accessions "
            f"unclustered, found {len(unclustered)} of {len(rows)}")

    def test_no_clustered_accession_is_absent_from_the_usher_tree(self, conn):
        """The invariant that makes the exclusion safe. If this ever fails,
        something really is wrong with placement and the validator should say so.
        """
        import re
        column = _cluster_column(conn)
        row = conn.execute(
            "SELECT newick FROM trees WHERE lower(source)='usher' LIMIT 1").fetchone()
        if not row:
            pytest.skip("no usher tree in this DB")
        tips = set(re.findall(r'[(,]([^(),:;]+):', row[0]))
        rows = conn.execute(
            f"SELECT primary_accession, {column} FROM meta_data").fetchall()
        clustered = {a for a, c in rows if not V.is_cluster_placeholder(c)}
        assert not (clustered - tips), (
            f"clustered accessions absent from the tree: "
            f"{sorted(clustered - tips)[:10]}")


# ---------------------------------------------------------------------------
# End to end: the validator itself must pass on the test database
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not HCV_TEST_DB.exists(), reason=f"not built: {HCV_TEST_DB}")
def test_validator_passes_on_the_subsampled_hcv_test_db(tmp_path):
    """The regression this file exists for.

    Before the fix this exited 1 with
    'Validation failed: UShER tree does not match accessions'.
    """
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "--db", str(HCV_TEST_DB),
         "--outdir", str(tmp_path), "--test-mode", "--check-update-integrity"],
        capture_output=True, text=True, timeout=1800,
    )
    assert result.returncode == 0, (
        f"validator failed:\n{result.stdout[-3000:]}\n{result.stderr[-2000:]}")
    assert "Validation status: PASS" in result.stdout
    assert "does not match accessions" not in result.stdout


@pytest.mark.skipif(not HCV_TEST_DB.exists(), reason=f"not built: {HCV_TEST_DB}")
def test_the_exclusion_is_reported_not_silent(tmp_path):
    """Dropping 222 accessions from a check without saying so would be worse
    than the bug it fixes."""
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), "--db", str(HCV_TEST_DB),
         "--outdir", str(tmp_path), "--test-mode", "--check-update-integrity"],
        capture_output=True, text=True, timeout=1800,
    )
    assert "Unclustered accessions" in result.stdout
    report = (tmp_path / "db_tree_validation.txt").read_text()
    assert "Unclustered accessions" in report
    assert "Missing in tree (clustered meta_data -> tree)" in report
