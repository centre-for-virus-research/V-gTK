"""The frozen GenBank XML fixture for the segmented test.

`FETCH_GENBANK` was 37.8% of a fresh `segmented_test` (129s of 341s total task
time), and it is the only step whose output depends on what NCBI returns on the
day. `segmented_xml_test` reads a committed snapshot instead: same taxid, same
reference list, same publish_dir, 3m54s -> 1m48s.

These tests keep the fixture and the profile honest. They do NOT run the
pipeline - they check that the snapshot exists, matches its checksums, and that
the profile actually points at it.
"""

import hashlib
import re
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPO_ROOT / "test_data" / "iav_11320"
XML_DIR = FIXTURE / "GenBank-XML"
CONFIG = REPO_ROOT / "nextflow.config"
FETCH_SCRIPT = REPO_ROOT / "scripts" / "fetch_segmented_fixture.sh"

requires_fixture = pytest.mark.skipif(
    not XML_DIR.is_dir(), reason=f"fixture not present: {XML_DIR}")


class TestFixtureContents:
    @requires_fixture
    def test_it_holds_xml(self):
        files = sorted(XML_DIR.glob("*.xml"))
        assert files, "fixture directory has no XML"
        assert len(files) > 10, f"only {len(files)} XML files - snapshot looks truncated"

    @requires_fixture
    def test_records_are_present(self):
        """A snapshot of empty batches would let the pipeline run and produce
        nothing, which is worse than failing."""
        total = sum(f.read_text(errors="replace").count("<GBSeq>")
                    for f in XML_DIR.glob("*.xml"))
        assert total > 400, f"only {total} GenBank records in the fixture"

    @requires_fixture
    def test_checksums_match(self):
        """The fixture is the test's input. If it drifts silently the test is
        no longer reproducible, which is the whole reason it exists."""
        checksums = FIXTURE / "checksums.sha256"
        if not checksums.exists():
            pytest.skip("no checksums file")
        bad = []
        for line in checksums.read_text().splitlines():
            if not line.strip():
                continue
            digest, name = line.split(None, 1)
            path = FIXTURE / name.strip()
            if not path.exists():
                bad.append(f"{name.strip()} missing")
                continue
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != digest:
                bad.append(f"{name.strip()} changed")
        assert not bad, f"fixture drifted: {bad[:5]}"

    @requires_fixture
    def test_manifest_records_provenance(self):
        """Whoever inherits this needs to know what it is and how to refresh it."""
        manifest = FIXTURE / "manifest.txt"
        assert manifest.exists(), "fixture has no manifest"
        text = manifest.read_text()
        assert "11320" in text, "manifest does not record the taxid"
        assert "fetch_segmented_fixture.sh" in text, \
            "manifest does not say how to refresh the fixture"

    @requires_fixture
    def test_it_is_small_enough_to_live_in_the_repo(self):
        size = sum(f.stat().st_size for f in XML_DIR.glob("*.xml"))
        assert size < 50 * 1024 * 1024, (
            f"fixture is {size/1e6:.0f} MB - too large to commit; reduce the "
            f"record set rather than growing the repo")


class TestProfileWiring:
    @pytest.fixture(scope="class")
    def config(self):
        return CONFIG.read_text()

    def test_the_profile_exists(self, config):
        assert re.search(r'^\s{4}segmented_xml_test\s*\{', config, re.M)

    def test_it_points_at_the_fixture(self, config):
        block = config[config.index("segmented_xml_test"):]
        block = block[:block.index("\n    }")]
        assert 'xml_dir' in block
        assert "test_data/iav_11320/GenBank-XML" in block

    def test_it_keeps_the_same_taxid_as_the_live_profile(self, config):
        """It must be a drop-in, not a different test."""
        def taxid(name):
            b = config[config.index(name + " {"):]
            b = b[:b.index("\n    }")]
            m = re.search(r'tax_id\s*=\s*"(\d+)"', b)
            return m.group(1) if m else None
        assert taxid("segmented_xml_test") == taxid("segmented_test") == "11320"

    def test_the_live_profile_still_exists(self, config):
        """Keep a path that exercises the real fetch, so NCBI changing under us
        is discoverable rather than invisible."""
        assert re.search(r'^\s{4}segmented_test\s*\{', config, re.M)

    def test_the_live_profile_does_not_use_the_fixture(self, config):
        block = config[config.index("segmented_test {"):]
        block = block[:block.index("\n    }")]
        assert "xml_dir" not in block


class TestRefreshScript:
    def test_it_exists_and_is_executable(self):
        assert FETCH_SCRIPT.exists()
        assert FETCH_SCRIPT.stat().st_mode & 0o111, "script is not executable"

    def test_it_targets_the_right_taxid_and_output(self):
        text = FETCH_SCRIPT.read_text()
        assert 'TAX_ID="11320"' in text
        assert "test_data/iav_11320" in text

    def test_it_uses_the_same_reference_list_as_the_profile(self):
        text = FETCH_SCRIPT.read_text()
        assert "generic/influenza/ref_list_refmast.txt" in text
