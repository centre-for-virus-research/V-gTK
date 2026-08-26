"""Extended invariant checks in scripts/ValidateDbTree.py.

ValidateDbTree.py is the only automated check that runs on a finished database,
so anything it does not look at is shipped unexamined. The checks exercised here
all target databases that open fine, load fine, and are quietly wrong:

  * a metadata column that stopped being written (the influenza neuraminidase
    gene name "NA" being read as a null erased a whole column exactly so),
  * a segment label that is blank or holds a gene name,
  * child rows addressed to accessions meta_data no longer describes,
  * an alignment that carries residues the submitted sequence never had,
  * a gene whose name was eaten by a wrong delimiter, breaking the hierarchy,
  * a tree leaf that resolves to nothing,
  * a collection date whose day and month were silently dropped.

Real databases are opened read-only through a file: URI and every real-DB test
skips cleanly when its database is absent.
"""

import datetime
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

import ValidateDbTree as vdt


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "ValidateDbTree.py"

HCV_DB = REPO_ROOT / "test_out" / "HCV_OM_test" / "HCV_OM_test.db"
IAV_DB = REPO_ROOT / "test_out" / "IAV_DB" / "iav-db.db"
RABV_UPDATE_DB = REPO_ROOT / "test_out" / "update_test" / "rabv-jul0425-update-test.db"
RABV_TREEFREE_DB = REPO_ROOT / "test_out" / "basic_test_treefree" / "rabv-jul0425.db"


def open_ro(db_path: Path) -> sqlite3.Connection:
    """Real databases are shared build artefacts; never open them writable."""
    return sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)


def finding_codes(result):
    return {finding["code"] for finding in result.get("findings", [])}


def finding(result, code):
    for item in result.get("findings", []):
        if item["code"] == code:
            return item
    return None


# ---------------------------------------------------------------------------
# synthetic DB helpers
# ---------------------------------------------------------------------------

META_COLUMNS = [
    "primary_accession", "organism", "accession_type", "exclusion_status",
    "segment", "length", "real_length", "a", "t", "g", "c", "n",
    "collection_date", "collection_year", "collection_mon", "collection_day",
    "country", "host", "host_taxa_id", "host_scientific_name", "country_validated",
    "taxonomy", "accession_version", "nearest_reference_genotype", "nearest_reference_subtype",
]


def _meta_row(accession, sequence, **overrides):
    row = {col: "" for col in META_COLUMNS}
    row.update(
        {
            "primary_accession": accession,
            "organism": "Test virus",
            "accession_type": "query",
            "exclusion_status": "0",
            "segment": "1",
            "length": str(len(sequence)),
            "real_length": str(sum(sequence.upper().count(b) for b in "ATGC")),
            "a": str(sequence.upper().count("A")),
            "t": str(sequence.upper().count("T")),
            "g": str(sequence.upper().count("G")),
            "c": str(sequence.upper().count("C")),
            "n": str(sequence.upper().count("N")),
            "collection_date": "12-Mar-2019",
            "collection_year": "2019",
            "collection_mon": "3",
            "collection_day": "12",
            "country": "Kenya",
            "host": "Homo sapiens",
            "host_taxa_id": "9606",
            "host_scientific_name": "Homo sapiens",
            "country_validated": "Kenya",
            "taxonomy": "Viruses",
            "accession_version": accession + ".1",
        }
    )
    row.update(overrides)
    return row


def build_db(
    tmp_path,
    name="synthetic.db",
    meta_rows=None,
    sequences=None,
    alignments=None,
    genes=None,
    trees=None,
    features=None,
    insertions=None,
):
    """Build a minimal but structurally faithful DB in tmp_path.

    Never writes into the repo: callers must pass pytest's tmp_path.
    """
    db_path = Path(tmp_path) / name
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            "CREATE TABLE meta_data (" + ", ".join(f'"{c}" TEXT' for c in META_COLUMNS) + ")"
        )
        rows = meta_rows if meta_rows is not None else [
            _meta_row("ACC1", "ACGTACGT"),
            _meta_row("ACC2", "ACGTACGA"),
        ]
        cur.executemany(
            "INSERT INTO meta_data VALUES (" + ", ".join("?" * len(META_COLUMNS)) + ")",
            [[row[c] for c in META_COLUMNS] for row in rows],
        )

        cur.execute("CREATE TABLE sequences (header TEXT, sequence TEXT, segment TEXT)")
        seqs = sequences if sequences is not None else [
            ("ACC1", "ACGTACGT", "1"),
            ("ACC2", "ACGTACGA", "1"),
        ]
        cur.executemany("INSERT INTO sequences VALUES (?, ?, ?)", seqs)

        cur.execute(
            "CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT, segment TEXT, sequence_id TEXT)"
        )
        alns = alignments if alignments is not None else [
            ("ACC1", "ACC1", "ACGTACGT", "1", "ACC1"),
            ("ACC2", "ACC1", "ACGTACGA", "1", "ACC2"),
        ]
        cur.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?, ?)", alns)

        cur.execute(
            "CREATE TABLE features (accession TEXT, master_ref_accession TEXT, reference_accession TEXT, aln_start TEXT, aln_end TEXT, cds_start TEXT, cds_end TEXT, product TEXT, segment TEXT)"
        )
        feats = features if features is not None else [
            ("ACC1", "ACC1", "ACC1", "1", "8", "1", "8", "N protein", "1"),
            ("ACC2", "ACC1", "ACC1", "1", "8", "1", "8", "N protein", "1"),
        ]
        cur.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", feats)

        cur.execute(
            "CREATE TABLE insertions (primary_accession TEXT, reference TEXT, insertion TEXT, segment TEXT)"
        )
        cur.executemany(
            "INSERT INTO insertions VALUES (?, ?, ?, ?)",
            insertions if insertions is not None else [],
        )

        cur.execute("CREATE TABLE genes (description TEXT, display_name TEXT, name TEXT, parent_name TEXT)")
        cur.executemany(
            "INSERT INTO genes VALUES (?, ?, ?, ?)",
            genes if genes is not None else [
                ("Nucleoprotein", "N protein", "N protein", "whole_genome"),
                ("Whole genome", "Whole genome", "whole_genome", "NULL"),
            ],
        )

        cur.execute("CREATE TABLE trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT)")
        cur.executemany(
            "INSERT INTO trees VALUES (?, ?, ?, ?, ?, ?)",
            trees if trees is not None else [
                ("usher_seg1", "usher", "seg_1", "1", "(ACC1:0.1,ACC2:0.1);", "2026-01-01"),
            ],
        )

        cur.execute("CREATE TABLE host_taxa (taxa_id TEXT, name TEXT, name_type TEXT, taxonomy_level TEXT)")
        cur.execute("INSERT INTO host_taxa VALUES ('9606', 'Homo sapiens', 'scientific name', 'species')")
        conn.commit()
    finally:
        conn.close()
    return db_path


# ---------------------------------------------------------------------------
# meta_data column population
# ---------------------------------------------------------------------------

def test_clean_db_reports_no_empty_core_columns(tmp_path):
    """A fully populated build must not raise a column-population warning, or the
    check is noise and will be ignored the day it matters."""
    db = build_db(tmp_path)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_meta_column_population(conn, meta_cols)
    assert "core_meta_column_entirely_empty" not in finding_codes(result)


def test_entirely_blank_core_column_is_flagged(tmp_path):
    """A column that quietly stopped being written stays a legal empty string in
    every row. Nothing else in the pipeline notices: the influenza gene name "NA"
    being coerced to a null erased a column exactly this way, and 65% of rows
    left VALIDATE_SEGMENT with an empty segment without a single error."""
    rows = [_meta_row("ACC1", "ACGTACGT", host=""), _meta_row("ACC2", "ACGTACGA", host="")]
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_meta_column_population(conn, meta_cols)
    item = finding(result, "core_meta_column_entirely_empty")
    assert item is not None and "host" in item["examples"]
    assert result["ok"] is False


def test_unpopulated_genotype_labels_are_reported_separately(tmp_path):
    """The influenza build carried nearest_reference_genotype on every one of 518
    rows and populated none of them, because the reference list was never seeded
    with labels. Nothing failed; the column simply stayed empty for months."""
    db = build_db(tmp_path)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_meta_column_population(conn, meta_cols)
    item = finding(result, "reference_label_columns_unpopulated")
    assert item is not None
    assert set(item["examples"]) == {"nearest_reference_genotype", "nearest_reference_subtype"}


def test_column_population_ignores_excluded_rows(tmp_path):
    """Excluded rows are frequently blank on purpose. If they counted as
    population the check would pass on a DB where every usable row is empty."""
    rows = [
        _meta_row("ACC1", "ACGTACGT", country=""),
        _meta_row("ACC2", "ACGTACGA", country="Kenya", exclusion_status="1"),
    ]
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        where_sql, params = vdt.get_meta_nonexcluded_filter(meta_cols, "exclusion_status", "1")
        result = vdt.validate_meta_column_population(conn, meta_cols, where_sql=where_sql, params=params)
    item = finding(result, "core_meta_column_entirely_empty")
    assert item is not None and "country" in item["examples"]


# ---------------------------------------------------------------------------
# segment labels
# ---------------------------------------------------------------------------

def test_blank_segment_on_non_excluded_rows_is_flagged(tmp_path):
    """A blank segment silently partitions the DB wrongly: per-segment
    alignments, trees and features all key off it. The tree-free RABV build ships
    with segment blank on all 78 rows while every other non-segmented build uses
    '1'."""
    rows = [_meta_row("ACC1", "ACGTACGT", segment=""), _meta_row("ACC2", "ACGTACGA", segment="")]
    db = build_db(
        tmp_path,
        meta_rows=rows,
        sequences=[("ACC1", "ACGTACGT", ""), ("ACC2", "ACGTACGA", "")],
        alignments=[("ACC1", "ACC1", "ACGTACGT", "", "ACC1"), ("ACC2", "ACC1", "ACGTACGA", "", "ACC2")],
        features=[("ACC1", "ACC1", "ACC1", "1", "8", "1", "8", "N protein", "")],
        trees=[],
    )
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_segment_labels(conn, meta_cols)
    assert finding(result, "blank_segment_in_meta_data")["count"] == 2


def test_segment_holding_a_gene_name_is_flagged(tmp_path):
    """Ten divergent segment normalisers inverted PB2/PB1 by scraping digits out
    of gene names. A segment that IS a gene name is the same accident one step
    earlier, and every downstream group-by silently splits on it."""
    rows = [_meta_row("ACC1", "ACGTACGT", segment="N protein"), _meta_row("ACC2", "ACGTACGA", segment="1")]
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_segment_labels(conn, meta_cols)
    item = finding(result, "segment_holds_a_gene_name")
    assert item is not None and item["examples"] == ["N protein"]


def test_child_table_segment_absent_from_meta_is_flagged(tmp_path):
    """If features carry segment '7' but no meta_data row does, those feature
    rows are unreachable through any per-segment query - the classic footprint of
    a segment normaliser applied in one place and not another."""
    db = build_db(
        tmp_path,
        features=[("ACC1", "ACC1", "ACC1", "1", "8", "1", "8", "N protein", "7")],
    )
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_segment_labels(conn, meta_cols)
    item = finding(result, "segment_absent_from_meta_data:features")
    assert item is not None and item["examples"] == ["7"]


def test_single_segment_db_must_label_segment_one(tmp_path):
    """Non-segmented viruses (HCV, RABV) label every row segment='1'. A lone
    segment with any other label means the label came from somewhere unexpected
    and will not match the reference list."""
    rows = [_meta_row("ACC1", "ACGTACGT", segment="L"), _meta_row("ACC2", "ACGTACGA", segment="L")]
    db = build_db(tmp_path, meta_rows=rows,
                  sequences=[("ACC1", "ACGTACGT", "L"), ("ACC2", "ACGTACGA", "L")],
                  alignments=[("ACC1", "ACC1", "ACGTACGT", "L", "ACC1"), ("ACC2", "ACC1", "ACGTACGA", "L", "ACC2")],
                  features=[("ACC1", "ACC1", "ACC1", "1", "8", "1", "8", "P", "L")],
                  genes=[("Whole genome", "Whole genome", "whole_genome", "NULL")],
                  trees=[("t", "usher", "k", "L", "(ACC1:0.1,ACC2:0.1);", "2026-01-01")])
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_segment_labels(conn, meta_cols)
    assert finding(result, "single_segment_not_labelled_1") is not None


# ---------------------------------------------------------------------------
# referential integrity
# ---------------------------------------------------------------------------

def test_orphan_child_rows_are_flagged(tmp_path):
    """An --update run that rewrites meta_data but not its children leaves rows
    addressed to accessions the DB no longer describes. They join to nothing, so
    they are invisible rather than loud."""
    db = build_db(
        tmp_path,
        sequences=[("ACC1", "ACGTACGT", "1"), ("GHOST", "ACGT", "1")],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_referential_integrity(conn, {"ACC1", "ACC2"})
    item = finding(result, "orphan_rows:sequences")
    assert item is not None and item["examples"] == ["GHOST"]


def test_non_excluded_accession_without_a_stored_sequence_is_flagged(tmp_path):
    """meta_data promises a sequence for every non-excluded accession. A missing
    sequences row makes the record appear in listings and exports and then serve
    nothing when opened."""
    db = build_db(tmp_path, sequences=[("ACC1", "ACGTACGT", "1")])
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_referential_integrity(conn, {"ACC1", "ACC2"})
    item = finding(result, "missing_from_child:sequences")
    assert item is not None and item["examples"] == ["ACC2"]


def test_unresolvable_alignment_name_is_flagged(tmp_path):
    """sequence_alignment.alignment_name names the reference a query was BLASTed
    against, and CreateSqliteDB._add_reference_columns uses exactly that value to
    look a query's genotype up in the reference table. When it names an accession
    meta_data does not carry, the lookup misses and nearest_reference_genotype
    comes back as an empty string - no exception, no log line, just an unlabelled
    query. All four shipped databases resolve every pointer, so this is a real
    invariant rather than an aspiration."""
    db = build_db(
        tmp_path,
        alignments=[
            ("ACC1", "MISSING_REF", "ACGTACGT", "1", "ACC1"),
            ("ACC2", "ACC1", "ACGTACGA", "1", "ACC2"),
        ],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_referential_integrity(conn, {"ACC1", "ACC2"})
    item = finding(result, "reference_pointer_unresolved:sequence_alignment.alignment_name")
    assert item is not None and item["examples"] == ["MISSING_REF"]


def test_clean_db_has_no_referential_findings(tmp_path):
    db = build_db(tmp_path)
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_referential_integrity(conn, {"ACC1", "ACC2"})
    assert result["ok"] is True, result["findings"]


# ---------------------------------------------------------------------------
# duplicate keys
# ---------------------------------------------------------------------------

def test_duplicate_accession_segment_pair_is_flagged(tmp_path):
    """An incremental update that inserts instead of upserting doubles rows.
    Joins still work and counts still look plausible, but every aggregate
    silently double-counts the duplicated accession."""
    rows = [_meta_row("ACC1", "ACGTACGT"), _meta_row("ACC1", "ACGTACGT"), _meta_row("ACC2", "ACGTACGA")]
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_duplicate_records(conn)
    item = finding(result, "duplicate_key:meta_data")
    assert item is not None and item["examples"] == ["ACC1 x2"]


def test_same_accession_on_two_segments_is_not_a_duplicate(tmp_path):
    """A segmented submission legitimately carries one accession per segment.
    Keying duplicates on the accession alone would condemn every influenza DB."""
    rows = [_meta_row("ACC1", "ACGTACGT", segment="4"), _meta_row("ACC1", "ACGTACGA", segment="6")]
    db = build_db(
        tmp_path,
        meta_rows=rows,
        sequences=[("ACC1", "ACGTACGT", "4")],
        alignments=[("ACC1", "ACC1", "ACGTACGT", "4", "ACC1")],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_duplicate_records(conn)
    assert result["ok"] is True, result["findings"]


def test_ambiguous_base_column_n_does_not_fake_duplicates(tmp_path):
    """meta_data has a column literally named "n" - the ambiguous-base count.
    A duplicate query written as `COUNT(*) AS n ... HAVING n > 1` silently binds
    to that column instead of the alias and reports every sequence with more than
    one N as a duplicated accession. Same shape as the segment/Segment collision:
    a legitimate column name colliding with a query's own vocabulary."""
    rows = [_meta_row("ACC1", "ACGTNNNN"), _meta_row("ACC2", "ACGTNNNN")]
    db = build_db(
        tmp_path,
        meta_rows=rows,
        sequences=[("ACC1", "ACGTNNNN", "1"), ("ACC2", "ACGTNNNN", "1")],
        alignments=[("ACC1", "ACC1", "ACGTNNNN", "1", "ACC1"), ("ACC2", "ACC1", "ACGTNNNN", "1", "ACC2")],
    )
    with sqlite3.connect(str(db)) as conn:
        n_values = [r[0] for r in conn.execute("SELECT n FROM meta_data")]
        assert n_values == ["4", "4"], "fixture must have >1 ambiguous base to be meaningful"
        result = vdt.validate_duplicate_records(conn)
    assert result["ok"] is True, result["findings"]


# ---------------------------------------------------------------------------
# sequence content / length coherence
# ---------------------------------------------------------------------------

def test_empty_and_non_nucleotide_sequences_are_flagged(tmp_path):
    """A zero-length sequence and a sequence full of protein residues are both
    legal TEXT. Only an alphabet check separates a nucleotide DB from one that
    quietly ingested the wrong FASTA."""
    rows = [_meta_row("ACC1", ""), _meta_row("ACC2", "ACGTACGA")]
    rows[0]["length"] = "0"
    db = build_db(
        tmp_path,
        meta_rows=rows,
        sequences=[("ACC1", "", "1"), ("ACC2", "EIQLLKKEKE", "1")],
        alignments=[("ACC1", "ACC1", "--------", "1", "ACC1"), ("ACC2", "ACC1", "ACGTACGA", "1", "ACC2")],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_sequence_content(conn)
    codes = finding_codes(result)
    assert "empty_sequence_rows" in codes
    assert "sequence_outside_nucleotide_alphabet" in codes


def test_iupac_ambiguity_codes_are_accepted(tmp_path):
    """R, Y, S, W, K, M, B, D, H, V and N are ordinary GenBank ambiguity codes;
    all four shipped databases contain them. Rejecting them would make the
    alphabet check fire on every real build."""
    seq = "ACGTRYSWKMBDHVN"
    rows = [_meta_row("ACC1", seq), _meta_row("ACC2", seq)]
    db = build_db(
        tmp_path,
        meta_rows=rows,
        sequences=[("ACC1", seq, "1"), ("ACC2", seq.lower(), "1")],
        alignments=[("ACC1", "ACC1", seq, "1", "ACC1"), ("ACC2", "ACC1", seq.lower(), "1", "ACC2")],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_sequence_content(conn)
    assert "sequence_outside_nucleotide_alphabet" not in finding_codes(result)


def test_meta_length_must_match_the_stored_sequence(tmp_path):
    """meta_data.length and the sequences row are produced by different pipeline
    steps. If either is written against the wrong accession both tables still
    look internally plausible - only the cross-check shows the swap. The relation
    holds exactly on all four shipped databases."""
    rows = [_meta_row("ACC1", "ACGTACGT"), _meta_row("ACC2", "ACGTACGA")]
    rows[0]["length"] = "9999"
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_sequence_content(conn)
    item = finding(result, "meta_length_disagrees_with_stored_sequence")
    assert item is not None and item["count"] == 1


def test_real_length_must_match_base_counts(tmp_path):
    """real_length is the unambiguous base count, i.e. a+t+g+c. When the base
    counter and the length calculator disagree, coverage and quality filters
    silently rank sequences on a number that describes a different sequence."""
    rows = [_meta_row("ACC1", "ACGTACGT"), _meta_row("ACC2", "ACGTACGA")]
    rows[1]["real_length"] = "3"
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_sequence_content(conn)
    assert finding(result, "real_length_disagrees_with_base_counts")["count"] == 1


# ---------------------------------------------------------------------------
# alignment geometry
# ---------------------------------------------------------------------------

def test_alignment_longer_than_raw_sequence_is_flagged(tmp_path):
    """An alignment row is the stored sequence with gaps inserted, so removing
    the gaps can never yield more residues than the raw sequence. When the
    padding step over-reaches it splices reference flanks onto the query and the
    DB serves an aligned sequence containing bases the submitter never sent -
    silently, because the raw sequence is stored in a different table."""
    db = build_db(
        tmp_path,
        alignments=[
            ("ACC1", "ACC1", "ACGTACGTAA", "1", "ACC1"),
            ("ACC2", "ACC1", "ACGTACGA--", "1", "ACC2"),
        ],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_alignment_geometry(conn)
    item = finding(result, "alignment_longer_than_raw_sequence")
    assert item is not None and item["examples"] == ["ACC1 ungapped=10 raw=8"]


def test_ragged_alignment_width_within_a_segment_is_flagged(tmp_path):
    """Every row of one segment's alignment must be the same width or column
    coordinates mean different things per row - feature projection then lands on
    the wrong bases without erroring."""
    db = build_db(
        tmp_path,
        alignments=[
            ("ACC1", "ACC1", "ACGTACGT", "1", "ACC1"),
            ("ACC2", "ACC1", "ACGTACGA----", "1", "ACC2"),
        ],
        sequences=[("ACC1", "ACGTACGT", "1"), ("ACC2", "ACGTACGA", "1")],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_alignment_geometry(conn)
    assert finding(result, "ragged_alignment_width") is not None


def test_all_gap_alignment_row_is_flagged(tmp_path):
    """A row of pure gaps means the sequence failed to align but was written
    anyway. It contributes a blank column to every downstream MSA view."""
    db = build_db(
        tmp_path,
        alignments=[
            ("ACC1", "ACC1", "--------", "1", "ACC1"),
            ("ACC2", "ACC1", "ACGTACGA", "1", "ACC2"),
        ],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_alignment_geometry(conn)
    assert finding(result, "all_gap_alignment_rows") is not None


def test_gapped_alignment_is_not_reported_as_a_bad_alphabet(tmp_path):
    """The alphabet predicate is a GLOB negated character class, and a literal
    '-' must sit LAST inside such a class or GLOB reads it as a range separator.
    Placed anywhere else the class stops matching gap characters and every
    gapped alignment in the DB - i.e. all of them - is reported as having an
    illegal alphabet. Same shape as the segment/Segment collision: a legitimate
    value ('-') that the surrounding syntax reinterprets."""
    seq = "ACGTACGT"
    db = build_db(
        tmp_path,
        alignments=[
            ("ACC1", "ACC1", "AC--GT--ACGT", "1", "ACC1"),
            ("ACC2", "ACC1", "AC..GT..ACGA", "1", "ACC2"),
        ],
        sequences=[("ACC1", seq, "1"), ("ACC2", "ACGTACGA", "1")],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_alignment_geometry(conn)
    assert "alignment_outside_nucleotide_alphabet" not in finding_codes(result)


def test_clean_alignment_passes(tmp_path):
    db = build_db(tmp_path)
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_alignment_geometry(conn)
    assert result["ok"] is True, result["findings"]


# ---------------------------------------------------------------------------
# genes table
# ---------------------------------------------------------------------------

def test_gene_with_a_null_looking_name_is_flagged(tmp_path):
    """Influenza's neuraminidase gene is called "NA". Read with pandas' default
    NA handling it becomes a null, and the shipped IAV database stores
    ('Neuraminidase', '', '', 'whole_genome') - a gene with no name, which no
    feature product can ever resolve to."""
    db = build_db(
        tmp_path,
        genes=[
            ("Neuraminidase", "", "", "whole_genome"),
            ("Whole genome", "Whole genome", "whole_genome", "NULL"),
        ],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_gene_table(conn)
    item = finding(result, "gene_with_no_name")
    assert item is not None and item["examples"] == ["Neuraminidase"]


def test_gene_name_that_swallowed_the_next_field_is_flagged(tmp_path):
    """generic/rabv/Tables/gene_info.csv ends its last row with four spaces
    instead of a tab, so the whole_genome row loads as a single field named
    'whole_genome    NULL'. generic/other/Tables/gene_info.csv uses a comma the
    same way. Nothing errors; the DB just gains a gene nobody can name."""
    db = build_db(
        tmp_path,
        genes=[
            ("Nucleoprotein", "N", "N", "whole_genome"),
            ("Whole genome", "Whole genome", "whole_genome    NULL", ""),
        ],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_gene_table(conn)
    codes = finding_codes(result)
    assert "gene_name_looks_like_a_merged_field" in codes
    assert "gene_parent_does_not_resolve" in codes


def test_dangling_gene_parent_is_flagged(tmp_path):
    """Every gene points at parent 'whole_genome'. If the whole_genome row is
    missing or misnamed the hierarchy is broken and gene selection falls back to
    an empty list rather than raising."""
    db = build_db(tmp_path, genes=[("Nucleoprotein", "N", "N", "whole_genome")])
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_gene_table(conn)
    assert finding(result, "gene_parent_does_not_resolve")["examples"] == ["N -> whole_genome"]


def test_literal_null_parent_sentinel_is_accepted(tmp_path):
    """The shipped gene_info files write the root's parent as the literal string
    'NULL'. Treating that as a dangling reference would make the check fire on
    every correct database."""
    db = build_db(tmp_path)
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_gene_table(conn)
    assert result["ok"] is True, result["findings"]


@pytest.mark.parametrize(
    "gene_info",
    [
        REPO_ROOT / "generic" / "rabv" / "Tables" / "gene_info.csv",
        REPO_ROOT / "generic" / "other" / "Tables" / "gene_info.csv",
    ],
    ids=["rabv", "other"],
)
@pytest.mark.xfail(
    reason="shipped gene_info files end their whole_genome row with spaces / a comma instead of a tab, so the row loads as 3 fields and the gene is named 'whole_genome    NULL'",
    strict=False,
)
def test_shipped_gene_info_rows_are_tab_separated(gene_info):
    """The gene_info file is loaded as TSV. A row with fewer than 4 tab-separated
    fields is not an error to any reader - the trailing values simply merge into
    the previous column, which is how a RABV build ends up with a gene called
    'whole_genome    NULL' and five children pointing at a parent that no longer
    exists."""
    if not gene_info.exists():
        pytest.skip(f"{gene_info} not present")
    lines = [line for line in gene_info.read_text().splitlines() if line.strip()]
    header = lines[0].split("\t")
    for line in lines[1:]:
        assert len(line.split("\t")) == len(header), f"ragged row in {gene_info.name}: {line!r}"


# ---------------------------------------------------------------------------
# trees
# ---------------------------------------------------------------------------

def test_tree_tip_not_in_meta_data_is_flagged(tmp_path):
    """Backbone genomes baked into the reference alignment become tree leaves
    without ever entering meta_data, sequences or sequence_alignment. Clicking
    one in the toolkit resolves to nothing. The existing validator only checks
    the single 'best' tree and, in segmented mode, returns before that check
    runs at all."""
    db = build_db(
        tmp_path,
        trees=[("usher_seg1", "usher", "seg_1", "1", "(ACC1:0.1,ACC2:0.1,NC_007357:0.1);", "2026-01-01")],
    )
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_tree_tip_resolution(conn)
    item = finding(result, "tree_tip_absent_from_meta_data")
    assert item is not None and item["examples"] == ["NC_007357"]


def test_tree_check_skips_cleanly_on_a_tree_free_build(tmp_path):
    """--tree_free is a supported run mode. The check must skip, not fail, or it
    would condemn every tree-free database."""
    db = build_db(tmp_path, trees=[])
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_tree_tip_resolution(conn)
    assert result.get("skipped") is True and result["ok"] is True


def test_unparsable_newick_is_flagged(tmp_path):
    db = build_db(tmp_path, trees=[("broken", "usher", "k", "1", "(ACC1:0.1,ACC2", "2026-01-01")])
    with sqlite3.connect(str(db)) as conn:
        result = vdt.validate_tree_tip_resolution(conn)
    assert finding(result, "unparsable_tree") is not None


def test_segment_without_a_tree_is_flagged(tmp_path):
    """A segmented build where one segment's tree job failed still writes a
    complete-looking DB; that segment just has no phylogeny."""
    rows = [_meta_row("ACC1", "ACGTACGT", segment="4"), _meta_row("ACC2", "ACGTACGA", segment="6")]
    db = build_db(
        tmp_path,
        meta_rows=rows,
        sequences=[("ACC1", "ACGTACGT", "4"), ("ACC2", "ACGTACGA", "6")],
        alignments=[("ACC1", "ACC1", "ACGTACGT", "4", "ACC1"), ("ACC2", "ACC2", "ACGTACGA", "6", "ACC2")],
        features=[("ACC1", "ACC1", "ACC1", "1", "8", "1", "8", "N protein", "4"),
                  ("ACC2", "ACC2", "ACC2", "1", "8", "1", "8", "N protein", "6")],
        trees=[("usher_seg4", "usher", "seg_4", "4", "(ACC1:0.1);", "2026-01-01")],
    )
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_segment_tree_coverage(conn, meta_cols)
    assert finding(result, "segment_without_tree")["examples"] == ["6"]


# ---------------------------------------------------------------------------
# collection dates
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "text,expected",
    [
        ("13-Feb-2006", (2006, 2, 13)),
        ("Feb-2006", (2006, 2, None)),
        ("2006", (2006, None, None)),
        ("2011-10-18", (2011, 10, 18)),
        ("2011-10", (2011, 10, None)),
        ("", (None, None, None)),
        ("unknown", (None, None, None)),
    ],
)
def test_collection_date_parser_handles_both_shipped_formats(text, expected):
    """GenBank records use 'DD-Mon-YYYY'; GISAID and recent INSDC submissions use
    ISO 'YYYY-MM-DD'. Both appear in every shipped database, so the validator has
    to read both or it will invent findings for whichever it cannot parse."""
    assert vdt._parse_collection_date(text) == expected


def test_iso_date_losing_day_and_month_is_flagged(tmp_path):
    """The day/month splitter historically only understood GenBank's
    'DD-Mon-YYYY'. Every ISO-dated record kept its year and silently lost day and
    month, so day-resolution filters dropped exactly the newest, best-annotated
    samples. All four shipped databases contain rows in this state: 8 in the IAV
    build, 35 in the RABV update build, 10 in the tree-free build."""
    rows = [
        _meta_row("ACC1", "ACGTACGT", collection_date="2011-10-18", collection_year="2011",
                  collection_mon="", collection_day=""),
        _meta_row("ACC2", "ACGTACGA"),
    ]
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_collection_dates(conn, meta_cols)
    item = finding(result, "date_precision_lost_in_split_columns")
    assert item is not None and item["count"] == 1
    assert "collection_mon+collection_day" in item["examples"][0]


def test_future_collection_date_is_flagged(tmp_path):
    """A sample cannot be collected tomorrow. A future date means a two-digit
    year was expanded into the wrong century, or a submission date was copied
    into the collection field - both produce a date that sorts and filters
    normally and is simply wrong."""
    future = datetime.date.today().replace(year=datetime.date.today().year + 3)
    rows = [
        _meta_row("ACC1", "ACGTACGT", collection_date=future.strftime("%d-%b-%Y"),
                  collection_year=str(future.year), collection_mon=str(future.month),
                  collection_day=str(future.day)),
        _meta_row("ACC2", "ACGTACGA"),
    ]
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_collection_dates(conn, meta_cols)
    assert finding(result, "collection_date_in_the_future")["count"] == 1


def test_collection_year_disagreeing_with_collection_date_is_flagged(tmp_path):
    """collection_year is the column every temporal plot groups on. When it
    disagrees with collection_date one of the two was written from a different
    record, and the plots are quietly wrong."""
    rows = [
        _meta_row("ACC1", "ACGTACGT", collection_date="13-Feb-2006", collection_year="2016",
                  collection_mon="2", collection_day="13"),
        _meta_row("ACC2", "ACGTACGA"),
    ]
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_collection_dates(conn, meta_cols)
    assert finding(result, "collection_year_disagrees_with_collection_date")["count"] == 1


def test_year_only_collection_date_does_not_trip_precision_check(tmp_path):
    """Most GenBank records carry only a year. Flagging their blank day/month
    would bury the real ISO-date losses in thousands of false positives."""
    rows = [_meta_row("ACC1", "ACGTACGT", collection_date="2006", collection_year="2006",
                      collection_mon="", collection_day=""),
            _meta_row("ACC2", "ACGTACGA")]
    db = build_db(tmp_path, meta_rows=rows)
    with sqlite3.connect(str(db)) as conn:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_collection_dates(conn, meta_cols)
    assert finding(result, "date_precision_lost_in_split_columns") is None


# ---------------------------------------------------------------------------
# CLI wiring
# ---------------------------------------------------------------------------

def _run_cli(db_path, outdir, *extra):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--db", str(db_path), "--outdir", str(outdir),
         "--allow-no-trees", *extra],
        capture_output=True,
        text=True,
    )


def test_invariant_findings_are_warnings_by_default(tmp_path):
    """A legitimate build can trip these checks (a reference list with no
    genotype labels, for instance), so a finding must not turn a good pipeline
    run red. It must still be written into the report."""
    db = build_db(tmp_path)
    outdir = tmp_path / "out"
    proc = _run_cli(db, outdir)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    report = (outdir / "db_tree_validation.txt").read_text()
    assert "Extended DB invariant checks" in report
    assert "reference_label_columns_unpopulated" in report
    assert "Validation status: PASS" in report


def test_strict_invariants_turns_findings_fatal(tmp_path):
    """CI for a production build wants the checks enforced; --strict-invariants
    is the switch that does it."""
    db = build_db(tmp_path)
    proc = _run_cli(db, tmp_path / "out_strict", "--strict-invariants")
    assert proc.returncode != 0
    assert "DB invariant checks failed" in (proc.stdout + proc.stderr)


def test_invariant_checks_can_be_skipped(tmp_path):
    db = build_db(tmp_path)
    outdir = tmp_path / "out_skip"
    proc = _run_cli(db, outdir, "--skip-invariant-checks")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "Extended DB invariant checks" not in (outdir / "db_tree_validation.txt").read_text()


# ---------------------------------------------------------------------------
# real databases (read-only)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not HCV_DB.exists(), reason=f"{HCV_DB} not present")
def test_hcv_db_satisfies_every_new_invariant():
    """The HCV build is the reference-quality artefact: non-segmented, with a
    genotype-labelled reference list and mutation tables. It must pass every new
    check, otherwise the checks are miscalibrated for a good DB."""
    conn = open_ro(HCV_DB)
    try:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        info = vdt.get_expected_accessions(conn)
        where_sql, params = vdt.get_meta_nonexcluded_filter(meta_cols, "exclusion_status", "1")
        results = vdt.run_extended_invariant_checks(
            conn, meta_cols, info["expected_accessions"], where_sql=where_sql, where_params=params
        )
    finally:
        conn.close()
    failures = {r["title"]: r["findings"] for r in results if vdt.result_failed(r)}
    assert failures == {}, failures


@pytest.mark.skipif(not IAV_DB.exists(), reason=f"{IAV_DB} not present")
@pytest.mark.xfail(
    reason="shipped IAV DB: 22 tree leaves (NC_0073xx reference genomes and 8 backbone accessions) appear in the segment trees but in no other table",
    strict=False,
)
def test_iav_tree_tips_all_resolve_to_meta_data():
    """Every segment tree in the influenza build carries 2-3 leaves that exist in
    no other table - not meta_data, not sequences, not sequence_alignment. They
    came from the reference alignment rather than the curated reference list. The
    build passed validation because segmented runs return before the existing
    tree/meta cross-check executes."""
    conn = open_ro(IAV_DB)
    try:
        result = vdt.validate_tree_tip_resolution(conn)
    finally:
        conn.close()
    assert result["ok"] is True, result["findings"]


@pytest.mark.skipif(not IAV_DB.exists(), reason=f"{IAV_DB} not present")
@pytest.mark.xfail(
    reason="shipped IAV DB: OP279655's ungapped alignment is 1410 residues against a 1408 bp raw sequence - the padding step appended two reference bases",
    strict=False,
)
def test_iav_alignments_are_never_longer_than_their_raw_sequence():
    """OP279655 (segment 6) is stored as 1408 bp but its alignment row ungaps to
    1410 residues, ending '...cctatataa' where the raw sequence ends
    '...cctatat'. The extra bases come from the reference tail, so the DB serves
    an aligned sequence the submitter never deposited."""
    conn = open_ro(IAV_DB)
    try:
        result = vdt.validate_alignment_geometry(conn)
    finally:
        conn.close()
    assert result["ok"] is True, result["findings"]


@pytest.mark.skipif(not IAV_DB.exists(), reason=f"{IAV_DB} not present")
@pytest.mark.xfail(
    reason="shipped IAV DB predates the keep_na_strings fix: the neuraminidase gene 'NA' was read as a null and stored with a blank name",
    strict=False,
)
def test_iav_genes_all_have_names():
    """CreateSqliteDB now passes keep_na_strings=True when loading gene_info, but
    the shipped database still carries ('Neuraminidase', '', '', 'whole_genome').
    The validator is what would have caught it at build time."""
    conn = open_ro(IAV_DB)
    try:
        result = vdt.validate_gene_table(conn)
    finally:
        conn.close()
    assert result["ok"] is True, result["findings"]


@pytest.mark.skipif(not RABV_TREEFREE_DB.exists(), reason=f"{RABV_TREEFREE_DB} not present")
@pytest.mark.xfail(
    reason="tree-free RABV build leaves meta_data.segment blank on all rows where HCV/RABV update builds use '1'",
    strict=False,
)
def test_treefree_rabv_labels_its_single_segment():
    """Non-segmented databases label every row segment='1' (HCV: 338/338, RABV
    update: 228/228). The tree-free build writes '' for all 78 rows and drops the
    segment column from sequences and features entirely, so any code that keys on
    segment sees a different shape depending on which run mode produced the DB."""
    conn = open_ro(RABV_TREEFREE_DB)
    try:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_segment_labels(conn, meta_cols)
    finally:
        conn.close()
    assert result["ok"] is True, result["findings"]


@pytest.mark.parametrize(
    "db_path",
    [IAV_DB, RABV_UPDATE_DB, RABV_TREEFREE_DB],
    ids=["iav", "rabv_update", "rabv_treefree"],
)
@pytest.mark.xfail(
    reason="shipped DBs predate the ISO-date fix in scripts/date_utils.py: every YYYY-MM-DD collection_date lost its day and month (8 rows IAV, 35 RABV update, 10 tree-free, 5945 in the 2.3 GB production HCV DB)",
    strict=False,
)
def test_real_dbs_kept_day_and_month_for_iso_collection_dates(db_path):
    """The day/month splitter only understood GenBank's 'DD-Mon-YYYY'. ISO dates
    kept their year and silently lost day and month, so the newest and
    best-annotated records - the ones INSDC is standardising on - are exactly the
    ones that drop out of any day-resolution filter. Fixed in date_utils.py for
    new builds; every shipped database still carries the damage, which is why the
    validator needs to see it rather than trusting the parser."""
    if not db_path.exists():
        pytest.skip(f"{db_path} not present")
    conn = open_ro(db_path)
    try:
        meta_cols = vdt.get_table_columns(conn, "meta_data")
        result = vdt.validate_collection_dates(conn, meta_cols)
    finally:
        conn.close()
    assert finding(result, "date_precision_lost_in_split_columns") is None, result["findings"]


@pytest.mark.parametrize(
    "db_path",
    [HCV_DB, IAV_DB, RABV_UPDATE_DB, RABV_TREEFREE_DB],
    ids=["hcv", "iav", "rabv_update", "rabv_treefree"],
)
def test_real_dbs_have_coherent_lengths_and_referential_integrity(db_path):
    """meta_data.length == LENGTH(sequences.sequence), real_length == a+t+g+c,
    and no orphan child rows: these hold exactly on every shipped database
    across every run mode (segmented and not, fresh and --update, tree-free).
    They are cheap, exact invariants, which makes them the right kind of tripwire
    for a build that silently writes a row against the wrong accession."""
    if not db_path.exists():
        pytest.skip(f"{db_path} not present")
    conn = open_ro(db_path)
    try:
        info = vdt.get_expected_accessions(conn)
        content = vdt.validate_sequence_content(conn)
        refs = vdt.validate_referential_integrity(conn, info["expected_accessions"])
        dupes = vdt.validate_duplicate_records(conn)
    finally:
        conn.close()
    assert content["ok"] is True, content["findings"]
    assert refs["ok"] is True, refs["findings"]
    assert dupes["ok"] is True, dupes["findings"]
