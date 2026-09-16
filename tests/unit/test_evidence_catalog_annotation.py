"""A long-format (one evidence reference per row) catalogue through the annotator and verifier.

A catalogue entry now spans one row per genotype, drug and evidence reference.
These tests pin what must NOT follow from that:

* the calls - sequence_mutation_calls and the two summary tables - are identical
  to those from one row per entry, because calls are made from a row's identity;
* rows for the same finding in two genotypes stay distinguishable in the
  database, and each keeps its own resistance category;
* a lookup table with a repeated evidence_id is refused rather than written,
  since it would double the evidence behind every row citing it;
* VerifyMutations, which reads the catalogue only for its reference accession,
  verifies a database annotated from the long catalogue unchanged.

Everything runs on a synthetic 30-base build in tmp_path.
"""

import sqlite3
import sys

import pandas as pd
import pytest

import AnnotateMutations as AM
import VerifyMutations

REFERENCE_30MER = "ATGAAACCCGGGTTTAGATAGCATGCATGC"   # codon 6 = AGA = R

BASE_COLUMNS = ["protein_name", "segment", "aa_position", "alt_residue", "reference_accession",
                "mutation_id", "mutation_type", "signature_id", "signature_kind", "combination_id",
                "combination_size", "phenotype", "resistance_category", "drug"]
LONG_COLUMNS = BASE_COLUMNS + ["genotype", "evidence_id", "data_source", "evidence_url", "evidence_label",
                               "evidence_type", "linked_evidence_ids", "finding_ids",
                               "relevant_genotypes", "wild_type_residues"]


def _entry(category, **extra):
    row = dict(zip(BASE_COLUMNS, ["NS3", "1", "6", "R", "REF1", "NS3:6R", "snp", "NS3:6R", "single",
                                  "", "", "drug_resistance", category, "drugA"]))
    row.update(relevant_genotypes="1a;1b", wild_type_residues="1a:K;1b:K")
    row.update(extra)
    return row


def _write_catalog(path, rows, columns):
    pd.DataFrame(rows, columns=columns).fillna("").to_csv(path, sep="\t", index=False)


def _build_db(path):
    conn = sqlite3.connect(str(path))
    cur = conn.cursor()
    cur.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, cds_start INTEGER, cds_end INTEGER)")
    cur.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ("REF1", "NS3", "1", 1, 30), ("Q1A", "polyprotein", "1", 1, 30),
        ("Q1B", "polyprotein", "1", 1, 30), ("Q3A", "polyprotein", "1", 1, 30)])
    cur.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, alignment_name TEXT, alignment TEXT)")
    cur.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)",
                    [(a, a, "REF1", REFERENCE_30MER) for a in ("REF1", "Q1A", "Q1B", "Q3A")])
    cur.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, "
                "nearest_reference_genotype TEXT, nearest_reference_subtype TEXT)")
    cur.executemany("INSERT INTO meta_data VALUES (?, ?, ?, ?)", [
        ("REF1", "master", "1", "a"), ("Q1A", "query", "1", "a"),
        ("Q1B", "query", "1", "b"), ("Q3A", "query", "3", "a")])
    conn.commit()
    conn.close()


def _annotate(db, catalog, **lookups):
    argv = sys.argv
    sys.argv = ["AnnotateMutations.py", "--db", str(db), "--mutation_catalog", str(catalog),
                "--catalog_column_profile", "HCV", "--virus", "HCV"]
    for flag, path in lookups.items():
        sys.argv += [f"--{flag}", str(path)]
    try:
        AM.main()
    finally:
        sys.argv = argv


def _table(db, name):
    conn = sqlite3.connect(str(db))
    try:
        frame = pd.read_sql_query(f"SELECT * FROM {name}", conn).fillna("")
    finally:
        conn.close()
    return frame.sort_values(list(frame.columns)).reset_index(drop=True)


@pytest.fixture
def two_builds(tmp_path):
    """The same two findings, one row each vs fanned out over genotype and evidence."""
    wide, long = tmp_path / "wide.tsv", tmp_path / "long.tsv"
    _write_catalog(wide, [_entry("category_I")], BASE_COLUMNS + ["relevant_genotypes", "wild_type_residues"])
    evidence = [
        dict(evidence_id="27773808", data_source="pubmed", evidence_type="in_vivo",
             linked_evidence_ids="NCT01717326", finding_ids="GZR_107"),
        dict(evidence_id="28228479", data_source="pubmed", evidence_type="in_vitro", finding_ids="GZR_166"),
        dict(evidence_id="NCT01717326", data_source="clinical_trial", evidence_type="in_vivo",
             linked_evidence_ids="27773808", finding_ids="GZR_107"),
    ]
    rows = [_entry("category_I", genotype="1a", **e) for e in evidence]
    rows += [_entry("category_II", genotype="1b", **evidence[0]), _entry("category_II", genotype="1b", **evidence[2])]
    _write_catalog(long, rows, LONG_COLUMNS)

    builds = {}
    for name, catalog in (("wide", wide), ("long", long)):
        db = tmp_path / f"{name}.db"
        _build_db(db)
        _annotate(db, catalog)
        builds[name] = db
    return builds


@pytest.mark.parametrize("table", ["sequence_mutation_calls", "sequence_relevant_mutation_summary",
                                   "completed_signatures_only"])
def test_evidence_fan_out_does_not_change_calls(two_builds, table):
    wide, long = _table(two_builds["wide"], table), _table(two_builds["long"], table)
    assert len(long) > 0 or table == "completed_signatures_only"
    pd.testing.assert_frame_equal(wide, long)


def test_calls_are_the_ones_expected(two_builds):
    calls = _table(two_builds["long"], "sequence_mutation_calls").set_index("primary_accession")["call_status"]
    # The master (typed 1a here) is annotated like any other sequence.
    assert calls.to_dict() == {"REF1": "emitted", "Q1A": "emitted", "Q1B": "emitted",
                               "Q3A": "suppressed_out_of_scope"}


def test_every_long_row_reaches_the_catalog_table_with_its_genotype(two_builds):
    catalog = _table(two_builds["long"], "mutation_catalog")
    assert len(catalog) == 5
    categories = catalog.groupby("genotype")["resistance_category"].agg(set).to_dict()
    assert categories == {"1a": {"category_I"}, "1b": {"category_II"}}
    assert set(catalog["evidence_id"]) == {"27773808", "28228479", "NCT01717326"}


def test_repeated_evidence_id_in_a_lookup_table_stops_the_write(tmp_path):
    conn = sqlite3.connect(str(tmp_path / "dup.db"))
    catalog = pd.DataFrame([_entry("category_I", genotype="1a", evidence_id="27773808", data_source="pubmed")])
    publications = pd.DataFrame([{"evidence_id": "27773808", "title": "a"},
                                 {"evidence_id": "27773808", "title": "b"}])
    try:
        with pytest.raises(ValueError, match="one row per evidence_id"):
            AM.write_mutation_tables(conn, catalog, [], "HCV", publications=publications)
    finally:
        conn.close()


def test_lookup_tables_resolve_every_reference(tmp_path):
    catalog = tmp_path / "long.tsv"
    _write_catalog(catalog, [
        _entry("category_I", genotype="1a", evidence_id="27773808", data_source="pubmed"),
        _entry("category_I", genotype="1a", evidence_id="AASLD_2017_Abs_1176", data_source="conference_abstract"),
        _entry("category_I", genotype="1a", evidence_id="NCT01717326", data_source="clinical_trial"),
        _entry("category_I", genotype="1a", evidence_id="UMIN000015627", data_source="clinical_trial"),
    ], LONG_COLUMNS)
    pubs, trials = tmp_path / "pubs.csv", tmp_path / "trials.csv"
    pubs.write_text("id,title,authors_short,year,journal,url\n"
                    "27773808,Paper,Komatsu et al.,2017,Gastroenterology,https://doi.org/x\n"
                    "AASLD_2017_Abs_1176,Abstract,Someone et al.,2017,Hepatology,https://doi.org/y\n")
    trials.write_text("id,display_name,nct_id\nC-WORTHY,C-WORTHy,NCT01717326\n"
                      "C-WORTHY Part D,C-WORTHy Part D,NCT01717326\nUMIN000015627,UMIN000015627,\n")
    db = tmp_path / "lookups.db"
    _build_db(db)
    _annotate(db, catalog, publications=pubs, clinical_trials=trials)

    conn = sqlite3.connect(str(db))
    try:
        unresolved = conn.execute("""
            SELECT mc.evidence_id FROM mutation_catalog mc
            LEFT JOIN publications p ON p.evidence_id = mc.evidence_id
            LEFT JOIN clinical_trials ct ON ct.evidence_id = mc.evidence_id
            WHERE NOT ((mc.data_source = 'clinical_trial' AND ct.evidence_id IS NOT NULL)
                    OR (mc.data_source != 'clinical_trial' AND p.evidence_id IS NOT NULL))""").fetchall()
        sources = dict(conn.execute("SELECT evidence_id, data_source FROM publications").fetchall())
        arms = conn.execute("SELECT trial_name FROM clinical_trials WHERE evidence_id = 'NCT01717326'").fetchone()
    finally:
        conn.close()
    assert unresolved == []
    assert sources == {"27773808": "pubmed", "AASLD_2017_Abs_1176": "conference_abstract"}
    assert arms == ("C-WORTHy;C-WORTHy Part D",)


def test_verify_mutations_passes_on_a_database_annotated_from_the_long_catalog(two_builds, tmp_path, capsys):
    long_catalog = tmp_path / "long.tsv"
    argv = sys.argv
    sys.argv = ["VerifyMutations.py", "--db", str(two_builds["long"]), "--mutation_catalog",
                str(long_catalog), "--sample_size", "10", "--virus", "HCV"]
    try:
        with pytest.raises(SystemExit) as excinfo:
            VerifyMutations.main()
    finally:
        sys.argv = argv
    out = capsys.readouterr().out
    assert excinfo.value.code == 0, out[-2000:]
    assert "Result: VERIFICATION SUCCESSFUL" in out
    assert "Mutations Checked: 3" in out and "Mismatches:        0" in out
