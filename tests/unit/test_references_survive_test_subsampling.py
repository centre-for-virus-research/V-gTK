"""References must survive test-mode subsampling, and be in the tree.

`TEST_SUBSAMPLE_CLUSTER_INPUT` keeps CI fast by capping the clustering input.
It used to do that with a plain `seqkit sample` over the whole alignment, which
has no idea what a reference is. On the HCV test profile that dropped **153 of
237 references** and the master was only retained by luck of the random seed.

That is not a cosmetic loss. References seed the MMseqs clusters, define the
IQ-TREE topology, anchor UShER placement, and are what genotype calls are made
against. Losing them produced a database whose UShER tree held 85 of 238
references, which then failed DB validation with a misleading message about
accessions not matching.

Two things changed, and both are tested here:

1. The cap is a **query** cap. Every reference is kept unconditionally, and the
   process fails outright if any reference is missing from its own output. It
   does not warn - nobody reads Nextflow warnings.
2. `ValidateDbTree` fails when any reference or master is absent from every
   tree, so if this regresses the validator says so plainly instead of blaming
   accession coverage.
"""

import re
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

import ValidateDbTree as V


REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE = REPO_ROOT / "vgtk-init.nf"
VALIDATOR = REPO_ROOT / "scripts" / "ValidateDbTree.py"

requires_seqkit = pytest.mark.skipif(
    shutil.which("seqkit") is None, reason="seqkit not on PATH")


# ---------------------------------------------------------------------------
# The Nextflow process contract
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def subsample_process():
    text = PIPELINE.read_text()
    start = text.index("process TEST_SUBSAMPLE_CLUSTER_INPUT")
    end = text.index("process ", start + 10)
    return text[start:end]


class TestSubsampleProcess:
    def test_it_receives_the_reference_list(self, subsample_process):
        """Without the ref list it cannot tell a reference from a query."""
        assert "path ref_list" in subsample_process

    def test_the_workflow_actually_passes_the_reference_list(self):
        text = PIPELINE.read_text()
        call = re.search(r'TEST_SUBSAMPLE_CLUSTER_INPUT\(([^)]*)\)', text)
        assert call, "process is never called"
        assert "ref_list" in call.group(1), (
            f"ref list not passed at the call site: {call.group(1)}")

    def test_the_cap_is_applied_to_queries_not_the_whole_file(self, subsample_process):
        """The original bug: `seqkit sample` over the whole alignment."""
        assert "queries_all.fasta" in subsample_process
        assert re.search(r'seqkit sample[^\n]*queries_all\.fasta', subsample_process), \
            "sampling must be applied to the query subset"
        assert not re.search(r'seqkit sample -n "\$MAX_SEQS" -s 42 "!\{dedup_msa\}"',
                             subsample_process), \
            "sampling the whole alignment is what dropped 153 references"

    def test_reference_selection_matches_by_id_not_full_header(self, subsample_process):
        """`seqkit grep -n` would match the full header and miss any record with
        a description. The id list is produced by `seqkit seq -n -i`, so the
        grep must match by id."""
        for pattern in ("seqkit grep -f ref_ids.txt", "seqkit grep -f query_ids.txt"):
            assert pattern in subsample_process, f"missing: {pattern}"
        assert "seqkit grep -n -f ref_ids.txt" not in subsample_process
        assert "seqkit grep -n -f query_ids.txt" not in subsample_process

    def test_carriage_returns_are_stripped_from_the_reference_list(self, subsample_process):
        """The backslash must be DOUBLED in the .nf source.

        The shell block is a Groovy triple-quoted string, which processes
        escape sequences: a single backslash-r there becomes a real carriage
        return in the generated .command.sh, which breaks the line and fails
        the process with exit 127. Asserting the doubled form is the only way
        to catch that from a unit test."""
        assert r"tr -d '\\r'" in subsample_process, (
            "CRLF ref lists exist in the repo, so the CR must be stripped - "
            "and the backslash must be doubled to survive Groovy escaping")

    def test_the_shell_block_contains_no_stray_triple_quote(self, subsample_process):
        """A triple quote inside the block ends the Groovy string early and
        the whole pipeline fails to compile. Easy to introduce in a comment
        that quotes shell syntax."""
        marker = chr(39) * 3
        body = subsample_process[subsample_process.index('shell:'):]
        start = body.index(marker) + 3
        rest = body[start:]
        assert marker in rest, 'shell block should be terminated'
        inner = rest[:rest.index(marker)]
        assert marker not in inner, (
            'a triple quote inside the shell block would end the Groovy '
            'string early and break compilation')

    def test_a_dropped_reference_is_a_hard_failure_not_a_warning(self, subsample_process):
        """Nextflow warnings are invisible in practice, so this must exit 1."""
        assert "exit 1" in subsample_process
        assert "ERROR" in subsample_process
        post = subsample_process[subsample_process.index("KEPT="):]
        assert "comm -23 ref_ids.txt kept_ids.txt" in post, \
            "there must be a post-condition comparing kept ids against references"

    def test_the_process_is_test_mode_only(self):
        text = PIPELINE.read_text()
        assert re.search(
            r'if\(\s*params\.test\s*==\s*"1"\s*\)\s*\{\s*\n\s*TEST_SUBSAMPLE_CLUSTER_INPUT',
            text), "subsampling must be gated on test mode"


class TestTestModeIsToldTruthfully:
    def test_validate_db_tree_is_only_given_test_mode_in_test_runs(self):
        """`--test-mode` relaxes real checks. It used to be hardcoded at both
        call sites, so production validation was weaker than test validation."""
        text = PIPELINE.read_text()
        for match in re.finditer(r'ValidateDbTree\.py(.{0,400}?)\'\'\'', text, re.S):
            block = match.group(1)
            assert "--test-mode" not in block, (
                "--test-mode is passed unconditionally to ValidateDbTree; it must "
                "be added to EXTRA_ARGS under a params.test check instead")

    def test_test_mode_is_still_reachable(self):
        """Gating it must not mean never passing it."""
        text = PIPELINE.read_text()
        gated = re.findall(
            r'if \[ "!\{params\.test\}" = "1" \]; then\s*\n\s*EXTRA_ARGS="\$\{EXTRA_ARGS\} --test-mode"',
            text)
        assert len(gated) == 2, (
            f"expected both ValidateDbTree call sites to gate --test-mode, found {len(gated)}")


# ---------------------------------------------------------------------------
# The shell logic, executed
# ---------------------------------------------------------------------------

@requires_seqkit
class TestSubsamplingKeepsReferences:
    def _run(self, tmp_path, n_refs, n_queries, cap, description=""):
        msa = tmp_path / "msa.fasta"
        refs = [f"REF{i}" for i in range(n_refs)]
        queries = [f"Q{i}" for i in range(n_queries)]
        with open(msa, "w") as fh:
            for name in refs + queries:
                fh.write(f">{name}{description}\nACGTACGTAC\n")
        ref_list = tmp_path / "ref_list.txt"
        ref_list.write_text("".join(f"{r}\treference\t1\t1\tNA\n" for r in refs))

        script = f'''
set -e
cd {tmp_path}
MAX_SEQS={cap}
cut -f1 ref_list.txt | sed '/^$/d' | sort -u > ref_ids_all.txt
seqkit seq -n -i msa.fasta | sort -u > msa_ids.txt
comm -12 ref_ids_all.txt msa_ids.txt > ref_ids.txt
comm -23 msa_ids.txt ref_ids.txt > query_ids.txt
N_QUERY=$(wc -l < query_ids.txt)
if [ "$N_QUERY" -le "$MAX_SEQS" ]; then
    cp msa.fasta out.fasta
else
    seqkit grep -f query_ids.txt msa.fasta -o queries_all.fasta
    seqkit sample -n "$MAX_SEQS" -s 42 queries_all.fasta -o queries_kept.fasta
    seqkit grep -f ref_ids.txt msa.fasta -o refs_kept.fasta
    cat refs_kept.fasta queries_kept.fasta > out.fasta
fi
seqkit seq -n -i out.fasta | sort -u > kept_ids.txt
comm -23 ref_ids.txt kept_ids.txt > lost.txt
'''
        subprocess.run(["bash", "-c", script], check=True, capture_output=True)
        kept = {l.strip() for l in (tmp_path / "kept_ids.txt").read_text().splitlines() if l.strip()}
        lost = [l.strip() for l in (tmp_path / "lost.txt").read_text().splitlines() if l.strip()]
        return set(refs), kept, lost

    def test_no_reference_is_lost_when_queries_exceed_the_cap(self, tmp_path):
        refs, kept, lost = self._run(tmp_path, n_refs=238, n_queries=500, cap=120)
        assert lost == [], f"references dropped: {lost[:10]}"
        assert refs <= kept

    def test_references_survive_even_when_they_outnumber_the_cap(self, tmp_path):
        """The HCV case: 238 references, cap of 120. All must be kept."""
        refs, kept, lost = self._run(tmp_path, n_refs=238, n_queries=400, cap=120)
        assert lost == []
        assert len(refs & kept) == 238

    def test_queries_are_actually_capped(self, tmp_path):
        refs, kept, lost = self._run(tmp_path, n_refs=10, n_queries=500, cap=50)
        queries_kept = {k for k in kept if k.startswith("Q")}
        assert len(queries_kept) <= 50, "the cap must still bite on queries"
        assert len(queries_kept) > 0

    def test_nothing_is_dropped_when_queries_are_under_the_cap(self, tmp_path):
        refs, kept, lost = self._run(tmp_path, n_refs=10, n_queries=20, cap=120)
        assert lost == []
        assert len(kept) == 30

    def test_descriptive_headers_do_not_lose_references(self, tmp_path):
        """`seqkit grep -n` matches the FULL header, not the id. Since the id
        list comes from `seqkit seq -n -i` (ids only), using -n would fail to
        match any record carrying a description and drop it - references
        included. Influenza headers routinely carry descriptions; HCV's do not,
        which is why this was invisible on the HCV test."""
        refs, kept, lost = self._run(
            tmp_path, n_refs=20, n_queries=200, cap=50,
            description=" Influenza A virus segment 4 (HA)")
        assert lost == [], f"references lost with descriptive headers: {lost[:10]}"
        assert len(refs & kept) == 20

    def test_crlf_reference_list_does_not_lose_references(self, tmp_path):
        """generic/rabv/ref_list_clades.txt and test_data/rabv_test_ref_list.txt
        are CRLF. A trailing CR on the accession matches nothing, so every
        reference would look like a query."""
        msa = tmp_path / "msa.fasta"
        refs = [f"REF{i}" for i in range(5)]
        with open(msa, "w") as fh:
            for name in refs + [f"Q{i}" for i in range(5)]:
                fh.write(f">{name}\nACGTACGTAC\n")
        (tmp_path / "ref_list.txt").write_text(
            "".join(f"{r}\treference\r\n" for r in refs))
        script = f'''
set -e
cd {tmp_path}
cut -f1 ref_list.txt | tr -d '\r' | sed '/^$/d' | sort -u > ref_ids_all.txt
seqkit seq -n -i msa.fasta | sort -u > msa_ids.txt
comm -12 ref_ids_all.txt msa_ids.txt > ref_ids.txt
'''
        subprocess.run(["bash", "-c", script], check=True, capture_output=True)
        matched = {l.strip() for l in (tmp_path / "ref_ids.txt").read_text().splitlines() if l.strip()}
        assert matched == set(refs), f"CRLF broke reference matching: {matched}"


# ---------------------------------------------------------------------------
# The validator check
# ---------------------------------------------------------------------------

def _db_with(tmp_path, refs_in_tree, all_refs, queries=()):
    path = tmp_path / "t.db"
    conn = sqlite3.connect(str(path))
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, "
                 "segment TEXT, cluster_95pct TEXT)")
    conn.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
    conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment TEXT)")
    conn.execute("CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, "
                 "segment TEXT, newick TEXT, created_at TEXT)")
    for i, r in enumerate(all_refs):
        kind = "master" if i == 0 else "reference"
        conn.execute("INSERT INTO meta_data VALUES (?,?,?,?)", (r, kind, "1", r))
    for q in queries:
        conn.execute("INSERT INTO meta_data VALUES (?,?,?,?)", (q, "query", "1", None))
    tips = list(refs_in_tree) + list(queries)
    newick = "(" + ",".join(f"{t}:0.1" for t in tips) + ");"
    conn.execute("INSERT INTO trees VALUES ('u','usher',NULL,'1',?,'now')", (newick,))
    conn.commit()
    return conn


class TestValidatorCatchesMissingReferences:
    def test_all_references_present_is_fine(self, tmp_path):
        conn = _db_with(tmp_path, ["M", "R1", "R2"], ["M", "R1", "R2"], ["Q1"])
        refs = {a for (a,) in conn.execute(
            "SELECT primary_accession FROM meta_data "
            "WHERE lower(accession_type) IN ('reference','master')")}
        tips = set(re.findall(r'[(,]([^(),:;]+):', conn.execute(
            "SELECT newick FROM trees").fetchone()[0]))
        assert not (refs - tips)

    def test_a_missing_reference_is_detectable(self, tmp_path):
        """The exact shape of the HCV failure: references absent from the tree."""
        conn = _db_with(tmp_path, ["M", "R1"], ["M", "R1", "R2"], ["Q1"])
        refs = {a for (a,) in conn.execute(
            "SELECT primary_accession FROM meta_data "
            "WHERE lower(accession_type) IN ('reference','master')")}
        tips = set(re.findall(r'[(,]([^(),:;]+):', conn.execute(
            "SELECT newick FROM trees").fetchone()[0]))
        assert refs - tips == {"R2"}

    def test_the_validator_declares_the_check(self):
        """The failure message must name references, not blame accession
        coverage - the original message sent us looking in the wrong place."""
        text = VALIDATOR.read_text()
        assert "references_missing_from_tree" in text
        assert "reference/master" in text
        assert "absent from the UShER tree" in text

    def test_the_check_applies_only_to_the_usher_tree(self):
        """IQ-TREE is built from cluster representatives, so a non-centroid
        reference is legitimately absent from it - measured, rabv_update's
        IQ-TREE holds 13 of 28 references while its UShER tree holds 28/28.
        Applying the check to IQ-TREE would fail correct base_tree_only runs."""
        text = VALIDATOR.read_text()
        assert 'fetch_trees(conn, source="usher")' in text
        assert "usher_tree_present" in text

    def test_the_check_runs_before_the_branching(self):
        """It must apply on every path, not only the segment-tree branch."""
        text = VALIDATOR.read_text()
        check = text.index("if references_missing_from_tree:")
        branch = text.index("if tree_source == args.segment_tree_source:", check - 4000)
        assert check < branch, "the reference check must precede the branch"


# ---------------------------------------------------------------------------
# Against the real databases
# ---------------------------------------------------------------------------

REFERENCE_DBS = {
    "HCV_full": REPO_ROOT / "test_out" / "HCV_full_15_jul26" / "HCV_full_aug" / "HCV_full.db",
    "rabv_update": REPO_ROOT / "test_out" / "update_test" / "rabv-jul0425-update-test.db",
    "IAV": REPO_ROOT / "test_out" / "IAV_DB" / "iav-db.db",
}


@pytest.mark.parametrize("name", sorted(REFERENCE_DBS))
def test_every_reference_is_in_a_tree_in_the_reference_dbs(name):
    """Measured: HCV_full 238/238, rabv_update 28/28, IAV 396/396.

    This is what makes the validator check safe to enforce rather than warn
    about - it is a real invariant of a correctly built database.
    """
    path = REFERENCE_DBS[name]
    if not path.exists():
        pytest.skip(f"not built: {path}")
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        refs = {a for (a,) in conn.execute(
            "SELECT primary_accession FROM meta_data "
            "WHERE lower(coalesce(accession_type,'')) IN ('reference','master')")}
        if not refs:
            pytest.skip("no references in this DB")
        rows = conn.execute(
            "SELECT newick FROM trees WHERE newick IS NOT NULL AND length(newick) > 0"
        ).fetchall()
        if not rows:
            pytest.skip("no trees in this DB")
        tips = set()
        for (newick,) in rows:
            tips |= set(re.findall(r'[(,]([^(),:;]+)[:,)]', newick))
        missing = refs - tips
        assert not missing, (
            f"{name}: {len(missing)} references absent from every tree: "
            f"{sorted(missing)[:10]}")
    finally:
        conn.close()
