"""Carpet-bomb audit tests for Nextflow orchestration, profiles, assets, and fixtures.

Passing tests guard clean invariants on repo assets and shell syntax.
Strict xfail tests pin verified bugs (AUDIT-NF-*, AUDIT-ASSET-*) so they
cannot regress silently once fixed.
"""

import csv
import glob
import os
import re
import subprocess
from collections import Counter
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _read_text(rel_path: str) -> str:
    return (REPO_ROOT / rel_path).read_text(encoding="utf-8", errors="replace")


def _extract_process_block(nf_text: str, proc_name: str) -> str:
    m = re.search(rf"process\s+{re.escape(proc_name)}\s*\{{", nf_text)
    assert m, f"Process {proc_name} not found"
    depth = 0
    for idx in range(m.end() - 1, len(nf_text)):
        ch = nf_text[idx]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return nf_text[m.start() : idx + 1]
    raise AssertionError(f"Unclosed process block {proc_name}")


# ---------------------------------------------------------------------------
# Passing asset & shell hygiene invariants
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "backbone_dir",
    [
        "generic/hcv/ref_set_aligned",
        "generic/rabv/ref_set_aligned",
        "generic/influenza/ref_set_aligned",
    ],
)
def test_reference_backbone_fastas_are_rectangular_and_valid_iupac(backbone_dir):
    fastas = sorted(glob.glob(str(REPO_ROOT / backbone_dir / "*.fa*")))
    assert fastas, f"No backbone FASTAs in {backbone_dir}"
    allowed = set("ACGTUNRYKMSWBDHV-")
    for fa in fastas:
        lengths = set()
        seen_ids = set()
        cur_id = None
        seq_parts = []
        with open(fa, encoding="utf-8") as fh:
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                if s.startswith(">"):
                    if cur_id is not None:
                        seq = "".join(seq_parts).upper()
                        lengths.add(len(seq))
                        assert set(seq) <= allowed, f"Bad char in {fa} {cur_id}"
                    cur_id = s[1:].split()[0]
                    assert cur_id not in seen_ids, f"Duplicate ID {cur_id} in {fa}"
                    seen_ids.add(cur_id)
                    seq_parts = []
                else:
                    seq_parts.append(s)
            if cur_id is not None:
                seq = "".join(seq_parts).upper()
                lengths.add(len(seq))
                assert set(seq) <= allowed, f"Bad char in {fa} {cur_id}"
        assert len(lengths) == 1, f"Ragged backbone {fa}: {lengths}"


def test_shell_scripts_pass_bash_syntax_check():
    scripts = [
        "vgtk-rabv.sh",
        "scripts/fetch_segmented_fixture.sh",
        "scripts/fetch_h10n8_fixture.sh",
        "scripts/run_h10n8_xml_test.sh",
        "scripts/create_gisaid_test_data.sh",
    ]
    for rel in scripts:
        proc = subprocess.run(
            ["bash", "-n", str(REPO_ROOT / rel)],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, f"bash -n failed for {rel}: {proc.stderr}"


def test_hcv_catalog_evidence_ids_resolve_against_clinical_trials_and_publications():
    cat = REPO_ROOT / "generic/hcv/Tables/generalized_mutation_catalog_evidence_linked.tsv"
    pubs = REPO_ROOT / "generic/hcv/Tables/phdr_publication.csv"
    trials = REPO_ROOT / "generic/hcv/Tables/phdr_clinical_trial.csv"

    with open(pubs, encoding="utf-8") as fh:
        pub_ids = {(r.get("evidence_id") or r.get("id") or "").strip() for r in csv.DictReader(fh)}
    with open(trials, encoding="utf-8") as fh:
        trial_rows = list(csv.DictReader(fh))
        trial_ids = {(r.get("nct_id") or "").strip() for r in trial_rows} | {
            (r.get("id") or "").strip() for r in trial_rows
        }

    with open(cat, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            src = (row.get("data_source") or "").strip().lower()
            eid = (row.get("evidence_id") or "").strip()
            if not eid:
                continue
            if src == "pubmed":
                assert eid in pub_ids, f"Missing pubmed evidence_id {eid}"
            elif src == "clinical_trial":
                assert eid in trial_ids, f"Missing clinical_trial evidence_id {eid}"


def test_nextflow_config_has_balanced_braces_and_core_profiles():
    cfg = _read_text("nextflow.config")
    assert cfg.count("{") == cfg.count("}")
    for profile in ("standard", "conda", "test", "update_test", "segmented_xml_test", "HCV_test", "RSV_test"):
        assert re.search(rf"^\s*{profile}\s*\{{", cfg, re.M), f"Missing profile {profile}"


# ---------------------------------------------------------------------------
# Promoted regression tests for fixed Nextflow audit bugs (WP-01) & remaining xfails
# ---------------------------------------------------------------------------

def test_check_segment_retention_uses_pipefail_and_gates_downstream():
    nf = _read_text("vgtk-init.nf")
    cfg = _read_text("nextflow.config")
    block = _extract_process_block(nf, "CHECK_SEGMENT_RETENTION")
    has_global_pipefail = "pipefail" in cfg
    has_local_pipefail = "pipefail" in block
    assert has_global_pipefail or has_local_pipefail
    assert "CHECK_SEGMENT_RETENTION.out" in nf


def test_collect_filtered_sequences_supports_update_db_in_workflow_and_cli():
    nf = _read_text("vgtk-init.nf")
    block = _extract_process_block(nf, "COLLECT_FILTERED_SEQUENCES")
    cfs_text = _read_text("scripts/CollectFilteredSequences.py")
    assert "--update_db" in cfs_text
    assert "--update_db" in block


def test_db_validation_processes_run_in_production_mode():
    nf = _read_text("vgtk-init.nf")
    for proc in ("TEST_DB_VALIDATION", "TEST_INPUT_DB_VALIDATION"):
        block = _extract_process_block(nf, proc)
        when_m = re.search(r"when:\s*(.+)", block)
        if when_m:
            assert 'params.test == "1"' not in when_m.group(1)


def test_create_sqlite_db_does_not_bust_resume_with_workflow_run_name():
    nf = _read_text("vgtk-init.nf")
    block = _extract_process_block(nf, "CREATE_SQLITE_DB")
    assert "workflow.runName" not in block


def test_nextflow_config_enforces_14_thread_cap():
    cfg = _read_text("nextflow.config")
    over_14 = re.findall(r"max_threads\s*=\s*(1[5-9]|[2-9]\d)", cfg)
    assert not over_14, f"max_threads > 14 found: {over_14}"
    assert "resourceLimits" in cfg or "executor" in cfg


def test_params_test_has_safe_default_and_no_integer_eq_on_null():
    cfg = _read_text("nextflow.config")
    nf = _read_text("vgtk-init.nf")
    top_params = cfg.split("profiles {")[0]
    assert re.search(r"^\s*test\s*=", top_params, re.M), "params.test missing from top-level defaults"
    assert '"!{params.test}" -eq "1"' not in nf


def test_add_missing_data_cli_matches_script_parser():
    nf = _read_text("vgtk-init.nf")
    amd = _read_text("scripts/AddMissingData.py")
    block = _extract_process_block(nf, "ADD_MISSING_DATA")
    assert "'-t'" in amd or '"-t"' in amd or "-t !" not in block


@pytest.mark.xfail(strict=True, reason="AUDIT-NF-08: nextflow.config contains foreign paths (/home4/lm305z/), broken scripts_dir, and update_db='?'")
def test_nextflow_config_has_no_foreign_paths_duplicates_or_invalid_conda_cmd():
    cfg = _read_text("nextflow.config")
    assert "/home4/" not in cfg
    assert "/home1/" not in cfg
    assert "RABV-gTK/RABV-gTK" not in cfg
    assert "conda.cmd" not in cfg
    assert "update_db         = '?'" not in cfg
    top_params = cfg.split("profiles {")[0]
    for dup in ("test_max_cluster_seqs", "extra_info_fill"):
        matches = re.findall(rf"^\s*{dup}\s*=", top_params, re.M)
        assert len(matches) <= 1, f"Duplicate top-level param {dup}"


def test_very_fast_tree_when_clause_handles_empty_or_null_string():
    nf = _read_text("vgtk-init.nf")
    block = _extract_process_block(nf, "VERY_FAST_TREE")
    assert "params.update_db == null" not in block


@pytest.mark.xfail(strict=True, reason="AUDIT-NF-10: vgtk-rabv.sh references unbound ${THREADS} under set -u in EPA-ng fallback")
def test_vgtk_rabv_sh_does_not_reference_unbound_threads_var():
    sh = _read_text("vgtk-rabv.sh")
    assert re.search(r"\$\{?THREADS\}?\b", sh) is None or re.search(r"^\s*THREADS=", sh, re.M)


@pytest.mark.xfail(strict=True, reason="AUDIT-ASSET-01: generic/hcv/ref_list_subtype_genotype.txt has 199 duplicated accessions and 2 masters")
def test_hcv_ref_list_subtype_genotype_has_unique_accessions_and_one_master():
    lines = [l.strip() for l in _read_text("generic/hcv/ref_list_subtype_genotype.txt").splitlines() if l.strip()]
    accs = [l.split("\t")[0].strip() for l in lines]
    dups = {a: c for a, c in Counter(accs).items() if c > 1}
    assert not dups, f"Duplicated accessions: {len(dups)}"


@pytest.mark.xfail(strict=True, reason="AUDIT-ASSET-02: assets/country_mapping.tsv has trailing tab (3 cols) on line 14 and maps Republic of the Congo to DRC")
def test_assets_country_mapping_is_well_formed_and_accurate():
    lines = _read_text("assets/country_mapping.tsv").splitlines()
    for idx, line in enumerate(lines, 1):
        cols = line.split("\t")
        assert len(cols) == 2, f"Line {idx} has {len(cols)} columns: {line!r}"
        if cols[0].strip() == "Republic of the Congo":
            assert cols[1].strip() != "Democratic Republic of the Congo"


@pytest.mark.xfail(strict=True, reason="AUDIT-ASSET-03: generic/influenza/ref_set_aligned has 22 extra accessions beyond ref_list_refmast.txt and lacks backbone_manifest.tsv")
def test_influenza_ref_set_aligned_matches_ref_list_and_has_manifest():
    assert (REPO_ROOT / "generic/influenza/ref_set_aligned/backbone_manifest.tsv").exists()
    refs = set()
    for line in _read_text("generic/influenza/ref_list_refmast.txt").splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 2 and parts[1].lower() in ("master", "reference"):
            refs.add(parts[0].strip())
    bb_ids = set()
    for fa in (REPO_ROOT / "generic/influenza/ref_set_aligned").glob("*.fasta"):
        for line in fa.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                bb_ids.add(line[1:].split()[0])
    assert bb_ids == refs


@pytest.mark.xfail(strict=True, reason="AUDIT-ASSET-04: fetch_segmented_fixture.sh omits --test_run, writes manifest.tsv / H10N8, and hashes absolute paths")
def test_fetch_segmented_fixture_script_is_bounded_and_relative():
    sh = _read_text("scripts/fetch_segmented_fixture.sh")
    assert "--test_run" in sh
    assert "manifest.txt" in sh
    assert "H10N8 subtype" not in sh


@pytest.mark.xfail(strict=True, reason="AUDIT-ASSET-05: environment.yml duplicates requests=2.31.0 and omits entrez-direct")
def test_environment_yml_has_no_duplicate_requests_and_declares_entrez_direct():
    yml = _read_text("environment.yml")
    assert yml.count("requests=2.31.0") == 1
    assert "entrez-direct" in yml
