"""The evidence-linked HCV catalogue, its builder, and the audit that checks it.

scripts/BuildHcvEvidenceCatalog.py turns the normaliser's output into the
catalogue the pipeline reads, one evidence reference per row.
generic/hcv/Tables/audit_catalog_linkage.py re-derives every link from the PHDR
tables. These tests pin: a fresh build passes the audit, the shipped file IS a
fresh build, the older hand-assembled wide catalogue rebuilds to the same file,
and the audit still catches the defects it was written for in the wide layout.
"""

import sys
from pathlib import Path

import pytest

import BuildHcvEvidenceCatalog as builder
import evidence_sources as ES

REPO_ROOT = Path(__file__).resolve().parents[2]
TABLES = REPO_ROOT / "generic" / "hcv" / "Tables"
NORMALISED_CATALOG = TABLES / "generalized_mutation_catalog.tsv"
WIDE_CATALOG = TABLES / "generalized_mutation_catalog_with_extra_info.tsv"
LONG_CATALOG = TABLES / "generalized_mutation_catalog_evidence_linked.tsv"

requires_tables = pytest.mark.skipif(
    not (TABLES / "phdr_resistance_finding.csv").exists() or not NORMALISED_CATALOG.exists(),
    reason="HCV PHDR tables not present")

sys.path.insert(0, str(TABLES))
import audit_catalog_linkage as audit  # noqa: E402


@pytest.fixture(scope="module")
def source():
    return audit.Source(TABLES)


@pytest.fixture(scope="module")
def fresh_build(tmp_path_factory):
    out = tmp_path_factory.mktemp("hcv") / "evidence_linked.tsv"
    builder.build(NORMALISED_CATALOG, TABLES, out)
    return out


def run_audit(source, catalog):
    _, rows = audit.read_rows(catalog, delimiter="\t")
    layout = "long" if "evidence_id" in rows[0] else "wide"
    entries = audit.collapse_entries(rows, layout)
    a = audit.Audit()
    audit.audit_source(a, source)
    audit.audit_keys(a, source, entries)
    audit.audit_values(a, source, entries)
    audit.audit_genotype(a, source, entries, rows, layout)
    if layout == "wide":
        audit.audit_evidence_wide(a, source, entries)
    else:
        audit.audit_evidence_long(a, source, rows)
    return {r["check"]: r for r in a.results}


def _rows(path):
    return audit.read_rows(path, delimiter="\t")


@requires_tables
def test_fresh_build_passes_audit(source, fresh_build):
    failed = {k: v["examples"] for k, v in run_audit(source, fresh_build).items() if v["status"] == "FAIL"}
    assert failed == {}


@requires_tables
def test_shipped_catalogue_is_a_fresh_build(fresh_build):
    """The pipeline reads the shipped file; it must not drift from its inputs."""
    assert LONG_CATALOG.read_bytes() == fresh_build.read_bytes()


@requires_tables
@pytest.mark.skipif(not WIDE_CATALOG.exists(), reason="wide catalogue not present")
def test_the_hand_assembled_wide_catalogue_rebuilds_to_the_same_file(fresh_build, tmp_path):
    """Nothing in the old hand-made file is lost by building from the normaliser."""
    out = tmp_path / "from_wide.tsv"
    builder.build(WIDE_CATALOG, TABLES, out)
    assert out.read_bytes() == fresh_build.read_bytes()


@requires_tables
def test_every_normalised_row_survives(fresh_build):
    _, normalised = _rows(NORMALISED_CATALOG)
    _, long = _rows(fresh_build)

    def key(r):
        return (r["id"], r["signature_id"], r["mutation_id"], builder.tidy_number(r["component_order"]))

    assert len(normalised) == len({key(r) for r in normalised})
    assert {key(r) for r in normalised} == {key(r) for r in long}


@requires_tables
def test_genotype_replaces_alignment_name(fresh_build):
    fields, rows = _rows(fresh_build)
    assert "alignment_name" not in fields
    assert fields.index("genotype") == _rows(NORMALISED_CATALOG)[0].index("alignment_name")
    assert not set(fields) & set(builder.WIDE_EVIDENCE_COLUMNS)
    for r in rows:
        expected = r["id"].rsplit(":", 2)[1][len("AL_"):] if r["id"] else ""
        assert r["genotype"] == expected


@requires_tables
def test_paper_and_trial_are_paired(fresh_build):
    _, rows = _rows(fresh_build)
    grz = {r["evidence_id"]: r for r in rows
           if r["id"] == "NS3:107I:AL_1a:grazoprevir" and r["mutation_id"] == "NS3:107I"}
    # The JBC paper is in-vitro only and reported no trial.
    assert grz["28228479"]["evidence_type"] == "in_vitro"
    assert grz["28228479"]["linked_evidence_ids"] == ""
    # C-EDGE TE was reported by two papers; C-WORTHy by one.
    assert grz["NCT02105701"]["linked_evidence_ids"] == "27720838;27773808"
    assert grz["NCT01717326"]["linked_evidence_ids"] == "27773808"
    assert grz["NCT01717326"]["data_source"] == ES.CLINICAL_TRIAL
    assert grz["NCT01717326"]["evidence_url"] == "https://clinicaltrials.gov/study/NCT01717326"
    assert {r["genotype"] for r in grz.values()} == {"1a"}


@requires_tables
def test_trial_arms_collapse_to_one_row_keeping_every_name(fresh_build):
    _, rows = _rows(fresh_build)
    magellan = {r["evidence_label"] for r in rows if r["evidence_id"] == "NCT02446717"}
    assert magellan and all(";" not in label or label == "Magellan-1, Part 1;Magellan-1, Part 2"
                            for label in magellan)


@requires_tables
def test_non_nct_trial_is_kept(fresh_build):
    _, rows = _rows(fresh_build)
    umin = [r for r in rows if r["evidence_id"] == "UMIN000015627"]
    assert {r["id"] for r in umin} == {"NS5A:31M:AL_1b:daclatasvir", "NS5A:93H:AL_1b:daclatasvir"}
    assert all(r["data_source"] == ES.CLINICAL_TRIAL and r["evidence_url"] == "" for r in umin)


@requires_tables
@pytest.mark.skipif(not WIDE_CATALOG.exists(), reason="wide catalogue not present")
def test_audit_catches_wide_catalogue_defects(source):
    results = run_audit(source, WIDE_CATALOG)
    for check in ("database rows whose resistance category cannot be pinned to a genotype",
                  "rows where paper <-> trial pairing cannot be recovered",
                  "rows mixing in-vitro and in-vivo sources in one pubmed_id list"):
        assert results[check]["status"] == "FAIL", check
    # Fixed at source in BuildCatalogGenotypeColumns.py.
    for check in ("wild_type_residues entry with a blank genotype code",
                  "clinical_trials set != trials via in-vivo findings"):
        assert results[check]["status"] == "PASS", check


@pytest.mark.parametrize("value, expected", [
    ("2.0", "2"), ("2", "2"), ("0.8", "0.8"), ("", ""), ("1807.2999999999997", "1807.3"),
    ("7.800000000000001", "7.8"), ("0.0001", "0.0001"), ("-0.0", "0"), ("category_I", "category_I"),
])
def test_tidy_number(value, expected):
    assert builder.tidy_number(value) == expected


@pytest.mark.parametrize("value, expected", [
    ("27773808", ES.PUBMED),
    ("AASLD_2017_Abs_1176", ES.CONFERENCE_ABSTRACT),
    ("EASL_2017_Abs_THU-257", ES.CONFERENCE_ABSTRACT),
])
def test_publication_source(value, expected):
    assert ES.publication_source(value) == expected


def test_trial_registry_id_falls_back_to_the_curator_id():
    assert ES.trial_registry_id("NCT01717326", "C-WORTHY Part D") == "NCT01717326"
    assert ES.trial_registry_id("", "UMIN000015627") == "UMIN000015627"
    assert ES.trial_url("UMIN000015627") == ""


def test_output_fields_are_the_same_for_normalised_and_wide_input():
    normalised = ["mutation_id", "alignment_name", "id", "drug", "phenotype", "relevant_genotypes"]
    wide = ["mutation_id", "alignment_name", "id", "drug", "phenotype", "drug_producer", "drug_category",
            "pubmed_id", "DOI", "relevant_genotypes", "wild_type_residues", "clinical_trials"]
    expected = (["mutation_id", "genotype", "id", "drug", "phenotype", "drug_producer", "drug_category"]
                + ES.EVIDENCE_COLUMNS + ["relevant_genotypes", "wild_type_residues"])
    assert builder.output_fields(normalised) == expected
    assert builder.output_fields(wide) == expected
