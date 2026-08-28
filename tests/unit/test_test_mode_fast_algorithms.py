"""Test runs use the cheapest version of every algorithm.

A smoke test only needs a tree and a clustering to exist and be self-consistent.
It does not need the best tree, so test mode swaps in the fastest available
variant of each step. Production must be untouched by any of this.

Measured on the HCV test data:

* MMseqs `linclust` (linear time) vs `cluster`: **29.3s -> 3.4s**, and on this
  input it produced the identical 276 representatives. linclust is documented as
  "less sensitive", so identical output is not guaranteed in general - only the
  speed is.
* IQ-TREE `--fast` ("fast search to resemble FastTree") on the 277-representative
  HCV alignment: **~13 min -> 14s**, at a slightly worse likelihood.

The risk this file guards against is the swap leaking into production, where a
rough topology or a less sensitive clustering would be silently wrong.
"""

import inspect
import re
from pathlib import Path

import pytest

import MMseqsClustering


REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE = REPO_ROOT / "vgtk-init.nf"


@pytest.fixture(scope="module")
def clustering_source():
    return inspect.getsource(MMseqsClustering.run_mmseqs_clustering)


class TestMmseqsFastPath:
    def test_fast_selects_linclust(self, clustering_source):
        assert 'subcommand = "linclust" if fast else "cluster"' in clustering_source

    def test_the_subcommand_is_actually_used_in_the_command(self, clustering_source):
        assert '"mmseqs", subcommand,' in clustering_source, \
            "the chosen subcommand must reach the command line"

    def test_max_seqs_is_suppressed_in_fast_mode(self, clustering_source):
        """linclust has no prefilter, so --max-seqs is rejected outright."""
        assert "if max_seqs is not None and not fast:" in clustering_source

    def test_fast_defaults_to_off(self):
        params = inspect.signature(MMseqsClustering.run_mmseqs_clustering).parameters
        assert params["fast"].default is False, \
            "production must get the sensitive algorithm without asking"

    def test_the_options_the_pipeline_relies_on_are_still_passed(self, clustering_source):
        """linclust keeps --cluster-mode/--cov-mode/--min-seq-id; the pipeline's
        choice of --cluster-mode 2 (greedy incremental, longest wins) is load
        bearing and must survive the swap."""
        for opt in ("--cluster-mode", "--cov-mode", "--min-seq-id", "--threads"):
            assert opt in clustering_source, f"{opt} dropped from the command"

    def test_the_cli_exposes_test_mode(self):
        source = (REPO_ROOT / "scripts" / "MMseqsClustering.py").read_text()
        assert '"--test_mode"' in source

    def test_both_call_sites_pass_it(self):
        source = (REPO_ROOT / "scripts" / "MMseqsClustering.py").read_text()
        assert source.count('fast=str(args.test_mode).strip() == "1"') == 2, \
            "every run_mmseqs_clustering call site must forward test mode"

    def test_only_the_literal_one_enables_it(self):
        """A truthy string like 'false' or '0' must not turn on the fast path."""
        source = (REPO_ROOT / "scripts" / "MMseqsClustering.py").read_text()
        assert 'strip() == "1"' in source


class TestIqTreeFastPath:
    @pytest.fixture(scope="class")
    def pipeline_text(self):
        return PIPELINE.read_text()

    def test_fast_is_added_only_under_test_mode(self, pipeline_text):
        block = pipeline_text[pipeline_text.index("process IQ_TREE"):]
        block = block[:block.index("\nprocess ", 10)]
        assert 'IQTREE_SPEED_ARGS="--fast"' in block
        gate = re.search(
            r'if \[ "!\{params\.test\}" = "1" \]; then\s*\n\s*IQTREE_SPEED_ARGS="--fast"',
            block)
        assert gate, "--fast must be gated on params.test"

    def test_the_variable_reaches_the_command_line(self, pipeline_text):
        assert "-mem !{params.iqtree_mem} ${IQTREE_SPEED_ARGS}" in pipeline_text

    def test_it_defaults_to_empty(self, pipeline_text):
        block = pipeline_text[pipeline_text.index("process IQ_TREE"):]
        block = block[:block.index("\nprocess ", 10)]
        assert 'IQTREE_SPEED_ARGS=""' in block, \
            "production must get no speed flag at all"

    def test_fast_is_not_applied_anywhere_unconditionally(self, pipeline_text):
        """A bare --fast on the iqtree line would degrade production silently."""
        for line in pipeline_text.splitlines():
            if "IQTREE_BIN" in line and "-s " in line:
                assert "--fast" not in line, f"unconditional --fast: {line.strip()[:90]}"


class TestMmseqsIsToldAboutTestMode:
    def test_the_pipeline_passes_test_mode_to_the_clustering_script(self):
        text = PIPELINE.read_text()
        call = text[text.index("MMseqsClustering.py"):]
        call = call[:call.index("'''")]
        assert "--test_mode !{params.test}" in call, \
            "MMseqsClustering must be told whether this is a test run"


class TestProductionIsUnaffected:
    """Every fast path must be reachable only through params.test."""

    def test_every_speed_switch_is_gated(self):
        text = PIPELINE.read_text()
        for marker in ('IQTREE_SPEED_ARGS="--fast"',):
            idx = text.index(marker)
            preceding = text[max(0, idx - 300):idx]
            assert 'params.test' in preceding, f"{marker} is not gated on test mode"

    def test_clustering_default_is_the_sensitive_algorithm(self, ):
        source = inspect.getsource(MMseqsClustering.run_mmseqs_clustering)
        # With fast=False the ternary must yield "cluster"
        assert 'if fast else "cluster"' in source
