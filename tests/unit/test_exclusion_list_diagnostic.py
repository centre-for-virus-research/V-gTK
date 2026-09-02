"""A reference-list `exclusion_list` row must be marked excluded.

A CI run of the segmented profile failed with:

    sequence_alignment vs meta_data: ok=False missing=22
    missing: NC_002204, NC_002205, ... NC_006306, NC_006307

Those 22 are the influenza B references that `ref_list_refmast.txt` marks
`exclusion_list` - they are in the list precisely so they never take part, and
they cannot align to influenza A data. They were in meta_data *without*
`exclusion_status=1`, so `get_expected_accessions` counted them as rows that
ought to have alignments.

The chain that should prevent that was verified end to end locally and is
correct at every step:

* `load_reference_file_table` parses the list with keep_default_na=False and
  detects the header, giving 418 rows: 8 master, 388 reference, 22 exclusion_list
* `load_reference_file_dict` keeps every row with a non-empty accession_type,
  so all 22 reach `ref_seq_dict`
* `GenBankParser` sets `accession_type` from that dict
* `FilterAndExtractSequences` sets `exclusion_status = 1` when the type is
  `exclusion_list`, and nothing later resets it to 0

An earlier theory - that a versioned accession (`NC_002204.1`) broke a lookup -
is **disproved**: `GenBankParser` sets `gi_number = primary_accession`, both come
from `GBSeq_primary-accession`, and none of the 518 records in the committed
fixture carry a version in that field (the version lives in a separate
`GBSeq_accession-version` element).

So the cause is still unknown, and these tests do not pretend otherwise. What
they do is make a recurrence self-describing: the validator now reports what KIND
of row is missing, so the next occurrence names the bug instead of listing bare
accessions.
"""

import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATOR = REPO_ROOT / "scripts" / "ValidateDbTree.py"
REF_LIST = REPO_ROOT / "generic" / "influenza" / "ref_list_refmast.txt"

sys.path.insert(0, str(REPO_ROOT / "scripts"))


# --------------------------------------------------------------------------
# The chain that assigns the flag
# --------------------------------------------------------------------------

@pytest.mark.skipif(not REF_LIST.exists(), reason="influenza ref list absent")
class TestReferenceListParsing:
    def test_the_exclusion_rows_survive_parsing(self):
        from ExportRefListFromUpdateDb import load_reference_file_table
        table = load_reference_file_table(str(REF_LIST))
        counts = table["accession_type"].value_counts().to_dict()
        assert counts.get("exclusion_list", 0) == 22, (
            f"expected 22 exclusion_list rows, got {counts}")

    def test_they_reach_the_dict_the_parser_uses(self):
        """load_reference_file_dict drops rows with a blank accession_type; the
        exclusion rows must not be dropped with them."""
        from ExportRefListFromUpdateDb import load_reference_file_dict
        d = load_reference_file_dict(str(REF_LIST))
        excluded = [k for k, v in d.items() if v.strip().lower() == "exclusion_list"]
        assert len(excluded) == 22
        assert "NC_002204" in excluded

    def test_the_header_is_not_mistaken_for_data(self):
        """The list has no header row; consuming the first row as one would
        shift every column and silently mislabel an accession."""
        from ExportRefListFromUpdateDb import load_reference_file_table
        table = load_reference_file_table(str(REF_LIST))
        assert len(table) == 418
        assert table.iloc[0]["primary_accession"].startswith(("A", "C", "D", "E",
                                                              "F", "G", "H", "J",
                                                              "K", "L", "M", "N"))


class TestGiNumberIsNotVersioned:
    """Disproves the version-suffix theory, so nobody re-derives it."""

    def test_gi_number_is_assigned_from_the_bare_accession(self):
        source = (REPO_ROOT / "scripts" / "GenBankParser.py").read_text()
        assert "content['gi_number'] = content['primary_accession']" in source, (
            "gi_number is no longer the bare primary accession; the version-suffix "
            "theory may need revisiting")

    @pytest.mark.skipif(
        not (REPO_ROOT / "test_data" / "iav_11320" / "GenBank-XML").is_dir(),
        reason="XML fixture absent")
    def test_no_primary_accession_in_the_fixture_carries_a_version(self):
        import re
        xml_dir = REPO_ROOT / "test_data" / "iav_11320" / "GenBank-XML"
        versioned = 0
        total = 0
        for f in xml_dir.glob("*.xml"):
            for m in re.finditer(r"<GBSeq_primary-accession>([^<]*)<", f.read_text(errors="replace")):
                total += 1
                if "." in m.group(1):
                    versioned += 1
        assert total > 400, f"only {total} records scanned"
        assert versioned == 0, f"{versioned} of {total} primary accessions are versioned"


# --------------------------------------------------------------------------
# The diagnostic, exercised against CI's exact condition
# --------------------------------------------------------------------------

def _db_with_unexcluded_exclusion_row(path):
    """Build the shape CI produced: an exclusion_list row that is NOT excluded
    and has no alignment."""
    conn = sqlite3.connect(str(path))
    conn.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, "
                 "exclusion_status TEXT, exclusion_criteria TEXT, segment TEXT, cluster_95pct TEXT)")
    conn.executemany("INSERT INTO meta_data VALUES (?,?,?,?,?,?)", [
        ("MASTER1", "master", "0", "", "1", "MASTER1"),
        ("QUERY1", "query", "0", "", "1", "MASTER1"),
        # the bug: excluded by the list, but not marked as such
        ("NC_002204", "exclusion_list", "0", "", "1", None),
    ])
    conn.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
    conn.executemany("INSERT INTO sequences VALUES (?,?)",
                     [("MASTER1", "ACGT"), ("QUERY1", "ACGT"), ("NC_002204", "ACGT")])
    conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment TEXT, segment TEXT)")
    conn.executemany("INSERT INTO sequence_alignment VALUES (?,?,?)",
                     [("MASTER1", "ACGT", "1"), ("QUERY1", "ACGT", "1")])
    conn.execute("CREATE TABLE features (accession TEXT)")
    conn.executemany("INSERT INTO features VALUES (?)", [("MASTER1",), ("QUERY1",)])
    conn.execute("CREATE TABLE trees (name TEXT, source TEXT, newick TEXT, segment TEXT, segment_key TEXT)")
    conn.execute("INSERT INTO trees VALUES ('u','usher','(MASTER1:1,QUERY1:1);','1',NULL)")
    conn.execute("CREATE TABLE host_taxa (other_col TEXT)")
    conn.execute("INSERT INTO host_taxa VALUES ('x')")
    conn.commit()
    conn.close()


class TestTheDiagnosticNamesTheCause:
    def test_it_reports_the_kind_of_row_not_just_the_accession(self, tmp_path):
        """CI printed 22 bare accessions and nothing about what they were."""
        db = tmp_path / "ci_shape.db"
        _db_with_unexcluded_exclusion_row(db)
        result = subprocess.run(
            [sys.executable, str(VALIDATOR), "--db", str(db),
             "--outdir", str(tmp_path), "--test-mode"],
            capture_output=True, text=True, timeout=600)
        assert "accession_type='exclusion_list'" in result.stdout, (
            f"the breakdown did not appear:\n{result.stdout[-2000:]}")

    def test_it_explains_why_that_combination_is_a_bug(self, tmp_path):
        db = tmp_path / "ci_shape2.db"
        _db_with_unexcluded_exclusion_row(db)
        result = subprocess.run(
            [sys.executable, str(VALIDATOR), "--db", str(db),
             "--outdir", str(tmp_path), "--test-mode"],
            capture_output=True, text=True, timeout=600)
        assert "will never have an alignment" in result.stdout
        assert "FilterAndExtractSequences" in result.stdout, \
            "the message should point at the step that assigns the flag"


@pytest.mark.skipif(not (REPO_ROOT / "test_out" / "IAV_DB" / "iav-db.db").exists(),
                    reason="IAV DB not built")
class TestTheInvariantHoldsInTheBuiltDatabase:
    def test_every_exclusion_list_row_is_marked_excluded(self):
        from ExportRefListFromUpdateDb import load_reference_file_dict
        if not REF_LIST.exists():
            pytest.skip("ref list absent")
        d = load_reference_file_dict(str(REF_LIST))
        excluded = {k for k, v in d.items() if v.strip().lower() == "exclusion_list"}
        conn = sqlite3.connect(
            f"file:{REPO_ROOT / 'test_out' / 'IAV_DB' / 'iav-db.db'}?mode=ro", uri=True)
        try:
            rows = dict(conn.execute(
                "SELECT primary_accession, exclusion_status FROM meta_data"))
        finally:
            conn.close()
        present = excluded & set(rows)
        if not present:
            pytest.skip("no exclusion_list accessions in this build")
        bad = [a for a in present if str(rows[a]).strip() != "1"]
        assert not bad, (
            f"{len(bad)} exclusion_list accessions are in meta_data without "
            f"exclusion_status=1: {sorted(bad)[:10]}. This is the CI failure.")
