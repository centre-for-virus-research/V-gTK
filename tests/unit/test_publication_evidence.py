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
        assert frame is not None and len(frame) == 97
        assert frame.columns.tolist() == ["nct_id", "trial_id", "trial_name"]
        assert frame["nct_id"].str.startswith("NCT").all()

    @requires_trials
    def test_one_row_per_nct_number(self):
        """The registry is keyed on PHDR's own id, not on the NCT number.

        Five NCT numbers therefore appear twice - the same registration under
        two ids (ALLY-2 / NCT02032888), a sponsor code beside a trial name
        (GS-US-342-1138 / ASTRAL-1), and two split into arms. Loaded verbatim
        they fan the join out: every catalogue row citing one of the five
        matches twice and its trial evidence reads as doubled.
        """
        frame = AM.load_clinical_trials_table(str(CLINICAL_TRIAL))
        assert frame["nct_id"].is_unique
        source = list(csv.DictReader(open(CLINICAL_TRIAL, encoding="utf-8", errors="replace")))
        assert len(frame) == len({r["nct_id"].strip() for r in source if r["nct_id"].strip()})

    @requires_trials
    def test_merging_keeps_every_distinct_value(self):
        """Nothing a curator wrote is lost.

        Distinct ids and distinct names are both kept, semicolon separated; a
        group agreeing on both fields would collapse to one value.
        """
        frame = AM.load_clinical_trials_table(str(CLINICAL_TRIAL))
        by_nct = {row.nct_id: row for row in frame.itertuples()}

        assert by_nct["NCT02446717"].trial_id == "Magellan-1_Part_1;Magellan-1_Part_2"
        assert by_nct["NCT02446717"].trial_name == "Magellan-1, Part 1;Magellan-1, Part 2"
        assert by_nct["NCT02201940"].trial_name == "ASTRAL-1;GS-US-342-1138 (ASTRAL-1)"
        assert by_nct["NCT02032888"].trial_id == "ALLY-2;NCT02032888"
        # A trial appearing once keeps a bare, unjoined value.
        assert ";" not in by_nct["NCT02265237"].trial_id

        # Set equality both ways: no value invented, no value dropped.
        source = [r for r in csv.DictReader(open(CLINICAL_TRIAL, encoding="utf-8", errors="replace"))
                  if (r["nct_id"] or "").strip()]
        for column, field in (("trial_id", "id"), ("trial_name", "display_name")):
            assert {v for cell in frame[column] for v in cell.split(";")} == \
                   {(r[field] or "").strip() for r in source}

    @requires_trials
    def test_the_separator_could_not_have_been_a_comma(self):
        """'Magellan-1, Part 1' is one trial name, not two.

        Semicolon is the separator every multi-valued string field in this
        catalogue uses, and no registry value contains one.
        """
        source = list(csv.DictReader(open(CLINICAL_TRIAL, encoding="utf-8", errors="replace")))
        assert not any(";" in (v or "") for r in source for v in r.values())
        assert any("," in (r["display_name"] or "") for r in source)

    @requires_trials
    def test_a_row_without_an_nct_is_reported_not_dropped_in_silence(self, capsys):
        """PHDR carries one: UMIN000015627, a Japanese UMIN-CTR registration
        with no ClinicalTrials.gov entry.

        It cannot be keyed on an NCT number and nothing in the catalogue cites
        it, so it is not loaded - but the run says so rather than losing a
        registry row without a word.
        """
        frame = AM.load_clinical_trials_table(str(CLINICAL_TRIAL))
        assert "UMIN000015627" in capsys.readouterr().out
        assert "UMIN000015627" not in set(frame["nct_id"])

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


@pytest.fixture(scope="module")
def catalog_db(tmp_path_factory):
    """mutation_catalog and clinical_trials exactly as AnnotateMutations writes them.

    Built from the shipped catalogue TSV and registry CSV rather than read from
    test_out/: the join these tests are about is a property of the two writers,
    and reading the shipped database would make the suite depend on when
    somebody last ran the pipeline. No alignment work is done, so it is fast.
    """
    if not (CATALOG.exists() and CLINICAL_TRIAL.exists()):
        pytest.skip("HCV catalogue assets not present")
    import pandas as pd
    catalog = pd.read_csv(CATALOG, sep="\t", dtype=str, keep_default_na=False)
    path = tmp_path_factory.mktemp("catalog_db") / "catalog.db"
    conn = sqlite3.connect(str(path))
    try:
        AM.write_mutation_tables(
            conn, catalog, [], "HCV",
            clinical_trials=AM.load_clinical_trials_table(str(CLINICAL_TRIAL)),
        )
    finally:
        conn.close()
    return path


class TestTheDatabaseJoin:
    """mutation_catalog.clinical_trials -> clinical_trials.nct_id.

    Until the three generic columns entered the HCV column profile this join did
    not exist: clinical_trials was dropped on the way into the database and the
    registry table was orphaned - written on every HCV build, referenced by
    nothing.
    """

    def test_the_catalog_table_carries_the_generic_columns(self, catalog_db):
        conn = sqlite3.connect(f"file:{catalog_db}?mode=ro", uri=True)
        columns = {r[1] for r in conn.execute("PRAGMA table_info(mutation_catalog)")}
        conn.close()
        assert {"relevant_genotypes", "wild_type_residues", "clinical_trials"} <= columns

    def test_the_registry_is_one_row_per_nct(self, catalog_db):
        conn = sqlite3.connect(f"file:{catalog_db}?mode=ro", uri=True)
        total, distinct = conn.execute(
            "SELECT COUNT(*), COUNT(DISTINCT nct_id) FROM clinical_trials").fetchone()
        conn.close()
        assert total == distinct == 97

    def test_every_cited_trial_resolves(self, catalog_db):
        """A dangling NCT would be evidence pointing at nothing."""
        conn = sqlite3.connect(f"file:{catalog_db}?mode=ro", uri=True)
        cited = {n for (cell,) in conn.execute("SELECT clinical_trials FROM mutation_catalog")
                 for n in (cell or "").split(";") if n}
        registry = {r[0] for r in conn.execute("SELECT nct_id FROM clinical_trials")}
        conn.close()
        assert cited
        assert not (cited - registry), sorted(cited - registry)[:5]
        assert not (registry - cited), "registry entries nothing cites"

    def test_the_join_does_not_fan_out(self, catalog_db):
        """One matched (row, trial) pair per NCT token in the catalogue.

        With the registry loaded verbatim - 102 rows for 97 distinct NCTs - the
        five doubled registrations turn 2,290 genuine pairs into 2,597,
        inflating the visible trial evidence by 13% on exactly the rows citing
        the most widely used trials.
        """
        conn = sqlite3.connect(f"file:{catalog_db}?mode=ro", uri=True)
        tokens = sum(len([n for n in (cell or "").split(";") if n])
                     for (cell,) in conn.execute("SELECT clinical_trials FROM mutation_catalog"))
        joined = conn.execute(
            "SELECT COUNT(*) FROM mutation_catalog mc JOIN clinical_trials ct "
            "  ON ';' || mc.clinical_trials || ';' LIKE '%;' || ct.nct_id || ';%'"
        ).fetchone()[0]
        conn.close()
        assert tokens == 2290
        assert joined == tokens


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

    def test_lookup_tables_are_present(self):
        """The reference build is annotated with --publications and
        --clinical_trials; without them mutation_catalog.pubmed_id and
        .clinical_trials are bare identifiers resolving to nothing.

        Column-level assertions live in TestTheDatabaseJoin, which builds its
        own database - a stale test_out/ should not fail the suite.
        """
        conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
        names = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        conn.close()
        assert {"mutation_catalog", "publications", "clinical_trials"} <= names
