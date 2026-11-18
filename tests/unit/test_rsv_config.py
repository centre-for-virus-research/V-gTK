"""Unit tests for Human RSV reference curation and Nextflow profiles, plus
divergent-mode reference list fixtures."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATE_REF_LIST = REPO_ROOT / "scripts" / "ValidateRefList.py"
RSV_REF_LIST = REPO_ROOT / "generic" / "rsv" / "ref_list.txt"
RSV_README = REPO_ROOT / "generic" / "rsv" / "README.md"
CYSTO_REF_LIST = REPO_ROOT / "test_data" / "cystoviridae" / "cystoviridae_ref_list.txt"
HANTA_REF_LIST = REPO_ROOT / "test_data" / "hantaviridae" / "hantaviridae_ref_list.txt"
NEXTFLOW_CONFIG = REPO_ROOT / "nextflow.config"


def _parse_ref_list(path: Path):
    rows = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        rows.append(line.split("\t"))
    return rows


class TestRsvReferenceSetAndProfiles:
    def test_rsv_assets_exist(self):
        assert RSV_REF_LIST.is_file(), f"Missing {RSV_REF_LIST}"
        assert RSV_README.is_file(), f"Missing {RSV_README}"

    def test_rsv_ref_list_passes_validate_ref_list_cli(self, tmp_path: Path):
        out = tmp_path / "rsv_valid.txt"
        proc = subprocess.run(
            [
                sys.executable,
                str(VALIDATE_REF_LIST),
                "-r",
                str(RSV_REF_LIST),
                "-s",
                "N",
                "-o",
                str(out),
            ],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"
        assert out.read_text(encoding="utf-8").strip() == "Reference list validation passed."

    def test_rsv_ref_list_covers_subgroups_a_and_b_with_bare_unique_accessions(self):
        rows = _parse_ref_list(RSV_REF_LIST)
        assert len(rows) == 30
        for row in rows:
            assert len(row) == 5, f"Expected 5 columns in RSV ref_list row: {row}"

        accs = [r[0] for r in rows]
        assert len(accs) == len(set(accs)), "Duplicate accessions in generic/rsv/ref_list.txt"
        for acc in accs:
            assert "." not in acc, f"Expected bare primary accession without version suffix: {acc}"

        masters = [r for r in rows if r[1] == "master"]
        assert len(masters) == 1
        assert masters[0][0] == "NC_038235"

        genotypes = {r[3] for r in rows}
        subtypes = {r[4] for r in rows}
        assert genotypes == {"A", "B"}
        # Both pre-duplication prototypes and modern G-duplication lineages (ON1, BA9) must be present
        for expected_subtype in ("GA1", "GA5", "NA1", "ON1", "GB1", "GB2", "GB3", "BA", "BA9"):
            assert expected_subtype in subtypes, f"Missing RSV lineage {expected_subtype}"

    def test_nextflow_config_defines_portable_rsv_profiles(self):
        cfg = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        for profile, expected_test in (("RSV_test", "1"), ("RSV_full", "0")):
            m = re.search(rf"\b{profile}\s*\{{.*?params\s*\{{(.*?)\}}\s*\}}", cfg, re.S)
            assert m is not None, f"Profile {profile} missing from nextflow.config"
            body = m.group(1)
            assert 'tax_id            = "3049954"' in body
            assert 'is_segmented      = "N"' in body
            assert f'test              = "{expected_test}"' in body
            assert 'generic/rsv/ref_list.txt' in body
            assert f'publish_dir       = "${{projectDir}}/test_out/{profile}"' in body
            assert "/dev/" not in body, f"{profile} publish_dir must not point into gitignored dev/"


@pytest.mark.parametrize(
    "ref_list_path",
    [CYSTO_REF_LIST, HANTA_REF_LIST],
    ids=["cystoviridae", "hantaviridae"],
)
def test_divergent_fixture_ref_lists_pass_validate_ref_list(ref_list_path: Path, tmp_path: Path):
    out = tmp_path / "valid.txt"
    proc = subprocess.run(
        [
            sys.executable,
            str(VALIDATE_REF_LIST),
            "-r",
            str(ref_list_path),
            "-s",
            "Y",
            "-o",
            str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"
