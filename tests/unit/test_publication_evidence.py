"""Publication and clinical-trial evidence behind a genotype/mutation/drug call.

The catalogue carries **one evidence reference per row**: ``evidence_id`` with a
``data_source`` saying what it is (pubmed, conference_abstract, clinical_trial).
It used to hold three independent semicolon lists (pubmed_id, DOI,
clinical_trials), which lost which paper reported which trial and which paper
backed the in-vitro versus the in-vivo columns.

Both lookup tables are keyed on ``evidence_id``, so the database join is a plain
equality with one row per key.

Evidence is **genotype-scoped at source**: a PHDR entry id is (RAS, alignment,
drug), and trial support genuinely varies by genotype - NS5A:31M against
daclatasvir cites one trial in genotype 1a and nine (plus a UMIN registration)
in 1b. Flattening would attach one genotype's evidence to another's call.
"""

import collections
import csv
import sqlite3
from pathlib import Path

import pytest

import AnnotateMutations as AM
import evidence_sources as ES
from BuildCatalogGenotypeColumns import build_trial_links

csv.field_size_limit(2 ** 31 - 1)

REPO_ROOT = Path(__file__).resolve().parents[2]
TABLES = REPO_ROOT / "generic" / "hcv" / "Tables"
PUBLICATIONS = TABLES / "phdr_publication.csv"
CLINICAL_TRIAL = TABLES / "phdr_clinical_trial.csv"
RESULT_TRIAL = TABLES / "phdr_result_trial.csv"
RESISTANCE_FINDING = TABLES / "phdr_resistance_finding.csv"
CATALOG = TABLES / "generalized_mutation_catalog_evidence_linked.tsv"
HCV_DB = REPO_ROOT / "test_out" / "HCV_OM_test" / "HCV_OM_test.db"

requires_publications = pytest.mark.skipif(
    not PUBLICATIONS.exists(), reason=f"not present: {PUBLICATIONS}")
requires_trials = pytest.mark.skipif(
    not (CLINICAL_TRIAL.exists() and RESULT_TRIAL.exists() and RESISTANCE_FINDING.exists()),
    reason="trial registry tables not present")
requires_catalog = pytest.mark.skipif(
    not CATALOG.exists(), reason=f"not present: {CATALOG}")


def _catalog_rows():
    return list(csv.DictReader(open(CATALOG, encoding="utf-8"), delimiter="\t"))


def _source_trials():
    return list(csv.DictReader(open(CLINICAL_TRIAL, encoding="utf-8", errors="replace")))


# --------------------------------------------------------------------------
# The trial registry itself
# --------------------------------------------------------------------------

class TestClinicalTrialRegistry:
    @requires_trials
    def test_identifiers_are_real_registry_numbers(self):
        rows = _source_trials()
        with_nct = [r for r in rows if (r.get("nct_id") or "").strip().startswith("NCT")]
        assert len(with_nct) >= 100
        assert len(with_nct) / len(rows) > 0.95

    @requires_trials
    def test_table_loads_with_the_expected_shape(self):
        frame = AM.load_clinical_trials_table(str(CLINICAL_TRIAL))
        assert frame is not None and len(frame) == 98
        assert frame.columns.tolist() == ["evidence_id", "nct_id", "trial_id", "trial_name", "url"]
        ncts = frame[frame["nct_id"] != ""]
        assert len(ncts) == 97
        assert (ncts["evidence_id"] == ncts["nct_id"]).all()
        assert ncts["url"].str.startswith("https://clinicaltrials.gov/study/NCT").all()

    @requires_trials
    def test_one_row_per_registry_id(self):
        """PHDR keys trials on a curator label, so five NCT numbers appear twice
        (ALLY-2 / NCT02032888, ASTRAL-1 / GS-US-342-1138, trial arms). Loaded
        verbatim they fan the join out."""
        frame = AM.load_clinical_trials_table(str(CLINICAL_TRIAL))
        assert frame["evidence_id"].is_unique
        expected = {ES.trial_registry_id(r["nct_id"], r["id"]) for r in _source_trials()}
        assert set(frame["evidence_id"]) == expected

    @requires_trials
    def test_merging_keeps_every_distinct_value(self):
        frame = AM.load_clinical_trials_table(str(CLINICAL_TRIAL))
        by_id = {row.evidence_id: row for row in frame.itertuples()}

        assert by_id["NCT02446717"].trial_id == "Magellan-1_Part_1;Magellan-1_Part_2"
        assert by_id["NCT02446717"].trial_name == "Magellan-1, Part 1;Magellan-1, Part 2"
        assert by_id["NCT02201940"].trial_name == "ASTRAL-1;GS-US-342-1138 (ASTRAL-1)"
        assert by_id["NCT02032888"].trial_id == "ALLY-2;NCT02032888"
        assert ";" not in by_id["NCT02265237"].trial_id

        # Set equality both ways: no value invented, no value dropped.
        source = _source_trials()
        for column, field in (("trial_id", "id"), ("trial_name", "display_name")):
            assert {v for cell in frame[column] for v in cell.split(";")} == \
                   {(r[field] or "").strip() for r in source}

    @requires_trials
    def test_the_separator_could_not_have_been_a_comma(self):
        """'Magellan-1, Part 1' is one trial name, not two."""
        source = _source_trials()
        assert not any(";" in (v or "") for r in source for v in r.values())
        assert any("," in (r["display_name"] or "") for r in source)

    @requires_trials
    def test_a_trial_without_an_nct_number_is_loaded_under_its_own_id(self):
        """UMIN000015627 is a Japanese UMIN-CTR registration. Two daclatasvir
        entries cite it; it used to be dropped here and from the catalogue."""
        frame = AM.load_clinical_trials_table(str(CLINICAL_TRIAL))
        umin = frame[frame["evidence_id"] == "UMIN000015627"]
        assert len(umin) == 1
        assert umin.iloc[0]["nct_id"] == "" and umin.iloc[0]["url"] == ""

    def test_absent_file_is_not_an_error(self):
        assert AM.load_clinical_trials_table(None) is None
        assert AM.load_clinical_trials_table("/nonexistent/trials.csv") is None

    @requires_trials
    def test_every_referenced_trial_resolves(self):
        registry = {r["id"].strip() for r in _source_trials()}
        referenced = {r["phdr_clinical_trial_id"].strip() for r in
                      csv.DictReader(open(RESULT_TRIAL, encoding="utf-8", errors="replace"))}
        assert not (referenced - registry), f"unresolved trial ids: {sorted(referenced - registry)[:5]}"


# --------------------------------------------------------------------------
# The link from a trial to a (mutation, genotype, drug)
# --------------------------------------------------------------------------

class TestTrialLinkage:
    @requires_trials
    def test_links_are_keyed_on_mutation_genotype_and_drug(self):
        links = build_trial_links(str(CLINICAL_TRIAL), str(RESULT_TRIAL), str(RESISTANCE_FINDING))
        assert len(links) > 400
        registry = {ES.trial_registry_id(r["nct_id"], r["id"]) for r in _source_trials()}
        for (ras, genotype, drug), ids in links.items():
            assert ras and genotype and drug
            assert set(ids) <= registry

    @requires_trials
    def test_trial_support_differs_between_genotypes(self):
        links = build_trial_links(str(CLINICAL_TRIAL), str(RESULT_TRIAL), str(RESISTANCE_FINDING))
        by_pair = collections.defaultdict(dict)
        for (ras, genotype, drug), ids in links.items():
            by_pair[(ras, drug)][genotype] = tuple(ids)
        assert len(by_pair[("NS5A:31M", "daclatasvir")]["1a"]) == 1
        assert len(by_pair[("NS5A:31M", "daclatasvir")]["1b"]) == 10
        multi = {k: v for k, v in by_pair.items() if len(v) > 1}
        differing = [k for k, v in multi.items() if len(set(v.values())) > 1]
        assert len(differing) > 30

    @requires_trials
    def test_only_in_vivo_findings_carry_a_trial(self):
        findings = list(csv.DictReader(open(RESISTANCE_FINDING, encoding="utf-8", errors="replace")))
        with_in_vivo = [f for f in findings if (f.get("phdr_in_vivo_result_id") or "").strip()]
        assert 0 < len(with_in_vivo) < len(findings)


# --------------------------------------------------------------------------
# The catalogue's evidence columns
# --------------------------------------------------------------------------

class TestEvidenceColumns:
    @requires_catalog
    def test_old_list_columns_are_gone(self):
        rows = _catalog_rows()
        assert set(ES.EVIDENCE_COLUMNS) <= set(rows[0])
        assert not {"pubmed_id", "DOI", "clinical_trials", "alignment_name"} & set(rows[0])

    @requires_catalog
    def test_one_reference_per_row_with_a_known_source(self):
        for row in _catalog_rows():
            assert ";" not in row["evidence_id"]
            if row["evidence_id"]:
                assert row["data_source"] in (ES.PUBMED, ES.CONFERENCE_ABSTRACT, ES.CLINICAL_TRIAL)
            else:
                assert not any(row[c] for c in ES.EVIDENCE_COLUMNS)

    @requires_catalog
    def test_publication_sources_follow_the_id(self):
        for row in _catalog_rows():
            if row["data_source"] in ES.PUBLICATION_SOURCES:
                assert row["data_source"] == ES.publication_source(row["evidence_id"])

    @requires_catalog
    def test_trials_attach_only_to_in_vivo_evidence_with_a_drug(self):
        for row in _catalog_rows():
            if row["data_source"] == ES.CLINICAL_TRIAL:
                assert row["evidence_type"] == "in_vivo"
                assert row["drug"]
                assert row["linked_evidence_ids"], "a trial always comes through a publication"

    @requires_catalog
    def test_paper_trial_links_agree_in_both_directions(self):
        rows = _catalog_rows()
        entry = lambda r: (r["id"], r["mutation_id"], r["component_order"])  # noqa: E731
        forward = {(entry(r), r["evidence_id"], t) for r in rows
                   if r["data_source"] in ES.PUBLICATION_SOURCES for t in r["linked_evidence_ids"].split(";") if t}
        backward = {(entry(r), p, r["evidence_id"]) for r in rows
                    if r["data_source"] == ES.CLINICAL_TRIAL for p in r["linked_evidence_ids"].split(";") if p}
        assert forward and forward == backward

    @requires_catalog
    def test_evidence_is_scoped_to_the_row_s_own_genotype(self):
        """Per genotype, not the signature's union across genotypes."""
        by_pair = collections.defaultdict(lambda: collections.defaultdict(set))
        for row in _catalog_rows():
            if row["genotype"] and row["evidence_id"]:
                by_pair[(row["signature_id"], row["drug"])][row["genotype"]].add(row["evidence_id"])
        differing = [k for k, v in by_pair.items() if len(v) > 1 and len({frozenset(s) for s in v.values()}) > 1]
        assert len(differing) > 100


# --------------------------------------------------------------------------
# Publications
# --------------------------------------------------------------------------

class TestPublications:
    @requires_publications
    def test_loads_with_the_expected_shape(self):
        frame = AM.load_publications_table(str(PUBLICATIONS))
        assert frame is not None and len(frame) == 128
        assert frame.columns.tolist() == ["evidence_id", "data_source", "title", "authors", "year", "journal", "url"]
        assert frame["evidence_id"].is_unique

    @requires_publications
    def test_conference_abstracts_are_labelled_not_parsed(self):
        """8 of the 128 references are conference abstracts with no PMID."""
        frame = AM.load_publications_table(str(PUBLICATIONS))
        abstracts = frame[frame["data_source"] == ES.CONFERENCE_ABSTRACT]["evidence_id"].tolist()
        assert len(abstracts) == 8
        assert all(any(tag in a for tag in ("AASLD", "EASL")) for a in abstracts)
        assert frame[frame["data_source"] == ES.PUBMED]["evidence_id"].str.isdigit().all()

    def test_absent_file_is_not_an_error(self):
        assert AM.load_publications_table(None) is None
        assert AM.load_publications_table("/nonexistent/publications.csv") is None


# --------------------------------------------------------------------------
# The database, as AnnotateMutations writes it
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def catalog_db(tmp_path_factory):
    """mutation_catalog and both lookup tables exactly as AnnotateMutations writes them.

    Built from the shipped assets rather than read from test_out/, so the suite
    does not depend on when somebody last ran the pipeline.
    """
    if not (CATALOG.exists() and CLINICAL_TRIAL.exists() and PUBLICATIONS.exists()):
        pytest.skip("HCV catalogue assets not present")
    import pandas as pd
    catalog = pd.read_csv(CATALOG, sep="\t", dtype=str, keep_default_na=False)
    path = tmp_path_factory.mktemp("catalog_db") / "catalog.db"
    conn = sqlite3.connect(str(path))
    try:
        AM.write_mutation_tables(
            conn, catalog, [], "HCV",
            publications=AM.load_publications_table(str(PUBLICATIONS)),
            clinical_trials=AM.load_clinical_trials_table(str(CLINICAL_TRIAL)),
        )
    finally:
        conn.close()
    return path


def _query(db, sql):
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


class TestTheDatabaseJoin:
    def test_the_catalog_table_carries_genotype_and_evidence(self, catalog_db):
        columns = {r[1] for r in _query(catalog_db, "PRAGMA table_info(mutation_catalog)")}
        assert {"genotype", "relevant_genotypes", "wild_type_residues", *ES.EVIDENCE_COLUMNS} <= columns
        assert not {"pubmed_id", "DOI", "clinical_trials", "alignment_name"} & columns

    def test_every_row_survives_the_load(self, catalog_db):
        """Rows differ in genotype or evidence, so de-duplication removes nothing."""
        (count,) = _query(catalog_db, "SELECT COUNT(*) FROM mutation_catalog")[0]
        assert count == len(_catalog_rows())

    def test_a_resistance_category_belongs_to_one_genotype(self, catalog_db):
        ambiguous = _query(catalog_db,
            "SELECT signature_id, drug, mutation_id, genotype FROM mutation_catalog WHERE drug != '' "
            "GROUP BY 1, 2, 3, 4 HAVING COUNT(DISTINCT resistance_category) > 1")
        assert ambiguous == []

    def test_lookup_tables_are_one_row_per_evidence_id(self, catalog_db):
        for table in ("publications", "clinical_trials"):
            total, distinct = _query(catalog_db, f"SELECT COUNT(*), COUNT(DISTINCT evidence_id) FROM {table}")[0]
            assert total == distinct > 0, table

    def test_every_cited_reference_resolves_through_its_own_table(self, catalog_db):
        dangling = _query(catalog_db, """
            SELECT mc.evidence_id, mc.data_source FROM mutation_catalog mc
            LEFT JOIN publications p ON p.evidence_id = mc.evidence_id
            LEFT JOIN clinical_trials ct ON ct.evidence_id = mc.evidence_id
            WHERE mc.evidence_id != ''
              AND NOT ((mc.data_source = 'clinical_trial' AND ct.evidence_id IS NOT NULL)
                    OR (mc.data_source IN ('pubmed', 'conference_abstract') AND p.evidence_id IS NOT NULL))""")
        assert dangling == []
        (uncited,) = _query(catalog_db,
            "SELECT COUNT(*) FROM clinical_trials WHERE evidence_id NOT IN (SELECT evidence_id FROM mutation_catalog)")[0]
        assert uncited == 0, "registry entries nothing cites"

    def test_the_join_does_not_fan_out(self, catalog_db):
        (cited,) = _query(catalog_db, "SELECT COUNT(*) FROM mutation_catalog WHERE data_source = 'clinical_trial'")[0]
        (joined,) = _query(catalog_db,
            "SELECT COUNT(*) FROM mutation_catalog mc JOIN clinical_trials ct ON ct.evidence_id = mc.evidence_id")[0]
        assert cited == joined > 0


@pytest.mark.skipif(not HCV_DB.exists(), reason="HCV reference database not built here")
class TestShippedDatabase:
    def test_lookup_tables_are_present(self):
        """Column-level assertions live in TestTheDatabaseJoin, which builds its
        own database - a stale test_out/ should not fail the suite."""
        names = {r[0] for r in _query(HCV_DB, "SELECT name FROM sqlite_master WHERE type='table'")}
        assert {"mutation_catalog", "publications", "clinical_trials"} <= names
