"""Clinical-trial and publication evidence behind a genotype/mutation/drug call.

Two independent evidence chains reach a call, and **both are genotype-scoped at
source**. Flattening either would attach one genotype's evidence to another
genotype's call.

**Clinical trials.** `phdr_clinical_trial.csv` carries real registry
identifiers — 102 of its 103 rows have an NCT number. `phdr_result_trial.csv`
links them to the in-vivo results that `phdr_resistance_finding.csv` attaches to
a `(RAS, alignment, drug)` key — a key that already contains the genotype. Trial
support genuinely varies by genotype: `NS5A:31M` against daclatasvir cites one
trial in genotype 1a and **nine** in 1b, and 49 of the 64 (RAS, drug) pairs
curated in more than one genotype have different trial sets.

**Publications.** `mutation_catalog.pubmed_id` and `.DOI` were already populated
and already genotype-scoped — 130 of the 173 multi-genotype (signature, drug)
pairs cite different papers. They were bare numbers with nothing in the database
saying what they were; the `publications` table resolves them.
"""

import collections
import csv
import sqlite3
from pathlib import Path

import pytest

import AnnotateMutations as AM
from BuildCatalogGenotypeColumns import build_trial_links

csv.field_size_limit(2 ** 31 - 1)

REPO_ROOT = Path(__file__).resolve().parents[2]
TABLES = REPO_ROOT / "generic" / "hcv" / "Tables"
PUBLICATIONS = TABLES / "phdr_publication.csv"
CLINICAL_TRIAL = TABLES / "phdr_clinical_trial.csv"
RESULT_TRIAL = TABLES / "phdr_result_trial.csv"
RESISTANCE_FINDING = TABLES / "phdr_resistance_finding.csv"
CATALOG = TABLES / "generalized_mutation_catalog_with_extra_info.tsv"
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


# --------------------------------------------------------------------------
# The trial registry itself
# --------------------------------------------------------------------------

class TestClinicalTrialRegistry:
    @requires_trials
    def test_identifiers_are_real_nct_numbers(self):
        """Not a heuristic. Before these tables arrived, trial status could only
        have been guessed from a publication title."""
        rows = list(csv.DictReader(open(CLINICAL_TRIAL, encoding="utf-8", errors="replace")))
        with_nct = [r for r in rows if (r.get("nct_id") or "").strip().startswith("NCT")]
        assert len(with_nct) >= 100
        assert len(with_nct) / len(rows) > 0.95

    @requires_trials
    def test_table_loads_with_the_expected_shape(self):
        frame = AM.load_clinical_trials_table(str(CLINICAL_TRIAL))
        assert frame is not None and len(frame) >= 100
        for column in ("nct_id", "trial_id", "trial_name"):
            assert column in frame.columns
        assert frame["nct_id"].str.startswith("NCT").all()

    def test_absent_file_is_not_an_error(self):
        """Optional: without it the NCT ids are still correct, just unresolvable."""
        assert AM.load_clinical_trials_table(None) is None
        assert AM.load_clinical_trials_table("/nonexistent/trials.csv") is None

    @requires_trials
    def test_every_referenced_trial_resolves(self):
        """A dangling trial reference would silently drop evidence."""
        registry = {r["id"].strip() for r in
                    csv.DictReader(open(CLINICAL_TRIAL, encoding="utf-8", errors="replace"))}
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
        for (ras, genotype, drug), ncts in links.items():
            assert ras and genotype and drug
            assert all(n.startswith("NCT") for n in ncts)

    @requires_trials
    def test_trial_support_differs_between_genotypes(self):
        """The reason the linkage must not be flattened to one list per mutation.

        NS5A:31M against daclatasvir: one trial in 1a, nine in 1b.
        """
        links = build_trial_links(str(CLINICAL_TRIAL), str(RESULT_TRIAL), str(RESISTANCE_FINDING))
        by_pair = collections.defaultdict(dict)
        for (ras, genotype, drug), ncts in links.items():
            by_pair[(ras, drug)][genotype] = tuple(ncts)

        multi = {k: v for k, v in by_pair.items() if len(v) > 1}
        differing = [k for k, v in multi.items() if len(set(v.values())) > 1]
        assert len(differing) > 30, (
            f"only {len(differing)} of {len(multi)} multi-genotype pairs differ - "
            f"trial evidence may have been flattened across genotypes"
        )

    @requires_trials
    def test_only_in_vivo_findings_carry_a_trial(self):
        """An in-vitro EC50 has no trial behind it, and must not acquire one."""
        findings = list(csv.DictReader(open(RESISTANCE_FINDING, encoding="utf-8", errors="replace")))
        with_in_vivo = [f for f in findings if (f.get("phdr_in_vivo_result_id") or "").strip()]
        assert 0 < len(with_in_vivo) < len(findings)


# --------------------------------------------------------------------------
# The catalogue column
# --------------------------------------------------------------------------

class TestClinicalTrialsColumn:
    @requires_catalog
    def test_column_exists_and_holds_only_nct_identifiers(self):
        rows = _catalog_rows()
        assert "clinical_trials" in rows[0]
        values = [n for r in rows for n in (r["clinical_trials"] or "").split(";") if n]
        assert values, "no row carries a trial"
        assert all(n.startswith("NCT") for n in values)

    @requires_catalog
    def test_semicolon_separated_like_the_other_generic_columns(self):
        rows = _catalog_rows()
        assert not any("," in (r["clinical_trials"] or "") for r in rows)

    @requires_catalog
    def test_a_row_carries_only_its_own_genotype_s_trials(self):
        """Per-row, not the signature's union across genotypes."""
        rows = _catalog_rows()
        by_pair = collections.defaultdict(dict)
        for row in rows:
            alignment = (row.get("alignment_name") or "").strip()
            if alignment:
                by_pair[(row["signature_id"], (row.get("drug") or "").strip())][alignment] = \
                    row["clinical_trials"]
        differing = [k for k, v in by_pair.items() if len(v) > 1 and len(set(v.values())) > 1]
        assert len(differing) > 50, (
            "rows of one signature all carry the same trials - the per-genotype "
            "detail has been unioned away"
        )

    @requires_catalog
    def test_rows_without_a_drug_carry_no_trial(self):
        """Trials attach to a drug. A row with none cannot have trial support."""
        for row in _catalog_rows():
            if not (row.get("drug") or "").strip():
                assert not (row["clinical_trials"] or "").strip()


# --------------------------------------------------------------------------
# Publications
# --------------------------------------------------------------------------

class TestPublications:
    @requires_publications
    def test_loads_with_the_expected_shape(self):
        frame = AM.load_publications_table(str(PUBLICATIONS))
        assert frame is not None and len(frame) == 128
        for column in ("pubmed_id", "title", "authors", "year", "journal", "url"):
            assert column in frame.columns

    @requires_publications
    def test_most_ids_are_pubmed_ids_but_conference_abstracts_are_not(self):
        """`pubmed_id` is really 'publication reference'.

        8 of the 128 entries are conference proceedings with no PMID -
        AASLD_2015_Abs_718, EASL_2017_Abs_THU-257 and similar. That is real
        evidence, not malformed data, so anything parsing this column as an
        integer will break on it.
        """
        rows = list(csv.DictReader(open(PUBLICATIONS, encoding="utf-8", errors="replace")))
        numeric = [r for r in rows if r["id"].strip().isdigit()]
        conference = [r["id"].strip() for r in rows if not r["id"].strip().isdigit()]
        assert len(numeric) > 0.9 * len(rows)
        assert conference, "expected some conference abstracts"
        assert all(any(tag in c for tag in ("AASLD", "EASL")) for c in conference), conference

    def test_absent_file_is_not_an_error(self):
        assert AM.load_publications_table(None) is None
        assert AM.load_publications_table("/nonexistent/publications.csv") is None

    @requires_catalog
    def test_publications_are_genotype_scoped_too(self):
        rows = _catalog_rows()
        by_pair = collections.defaultdict(dict)
        for row in rows:
            alignment = (row.get("alignment_name") or "").strip()
            if alignment:
                by_pair[(row["signature_id"], (row.get("drug") or "").strip())][alignment] = \
                    (row.get("pubmed_id") or "").strip()
        multi = {k: v for k, v in by_pair.items() if len(v) > 1}
        differing = [k for k, v in multi.items() if len(set(v.values())) > 1]
        assert len(differing) > 100

    @requires_catalog
    def test_pubmed_ids_are_semicolon_separated(self):
        """Values are PMIDs or conference abstract references, never comma-joined."""
        multi = [r["pubmed_id"] for r in _catalog_rows() if ";" in (r.get("pubmed_id") or "")]
        assert multi
        for value in multi[:50]:
            assert "," not in value, "separator drifted from semicolon"
            for reference in value.split(";"):
                reference = reference.strip()
                assert reference
                assert reference.isdigit() or any(t in reference for t in ("AASLD", "EASL"))


@pytest.mark.skipif(not HCV_DB.exists(), reason="HCV reference database not built here")
class TestShippedDatabase:
    def test_catalog_carries_pubmed_ids(self):
        conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
        total, with_pmid = conn.execute(
            "SELECT COUNT(*), SUM(CASE WHEN TRIM(COALESCE(pubmed_id,'')) != '' THEN 1 ELSE 0 END) "
            "FROM mutation_catalog"
        ).fetchone()
        conn.close()
        assert total > 0 and with_pmid > 0.9 * total

    def test_lookup_tables_appear_once_the_database_is_rebuilt(self):
        """Documents current state: the shipped DB predates these options.

        When this starts failing the reference database has been rebuilt with
        --publications / --clinical_trials, which is the intended end state.
        """
        conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
        names = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        conn.close()
        assert "mutation_catalog" in names
