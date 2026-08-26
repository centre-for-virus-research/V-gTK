"""Edge cases in the non-segmented (HCV / RABV) SQLite build path.

Every case here is about *silent* wrongness in `CreateSqliteDB.create_db()`:
rows that vanish from `features`/`sequences`/`sequence_alignment`/`insertions`
without a single line of output, segment values invented or mangled, and curated
genotype/subtype labels overwritten by inferred ones.

Two mechanisms produce most of it and are worth stating once:

* A single-segment virus has no segment. The build therefore *infers* that it is
  unsegmented (`_should_force_unsegmented_segment_one`) and stamps segment ``'1'``
  everywhere. Anything that makes one table's segment string differ from
  meta_data's - a zero-padded ``'01'`` in the reference list, an explicit ``'1'``
  in features while the matrix is blank - both defeats that inference *and*
  makes the row fail the `valid_pairs` join, which drops it with no report.

* `valid_pairs` is an inner join on ``(accession, segment)`` applied to four
  child tables. It has no "dropped N rows" counter, so a join key that is merely
  spelled differently is indistinguishable from a genuinely excluded record.

The fixtures are deliberately HCV/RABV-shaped: one master/reference set with
genotype+subtype, a handful of queries, blank segment columns.
"""

import io
import contextlib
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from CreateSqliteDB import CreateSqliteDB


REPO_ROOT = Path(__file__).resolve().parents[2]
HCV_DB = REPO_ROOT / "test_out" / "HCV_OM_test" / "HCV_OM_test.db"
HCV_REF_LIST = REPO_ROOT / "generic" / "hcv" / "ref_list_subtype_genotype.txt"


# --------------------------------------------------------------------------
# fixture plumbing
# --------------------------------------------------------------------------

def _write_tsv(path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _write_csv(path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def _build_db(
    tmp_path,
    meta_rows,
    meta_cols=("primary_accession", "accession_type", "segment"),
    feat_rows=None,
    feat_cols=("accession", "product"),
    aln_rows=None,
    aln_cols=("sequence_id", "primary_accession", "aligned_seq"),
    fasta_ids=None,
    ins_rows=None,
    ins_cols=("primary_accession", "insertions"),
    ref_rows=None,
    db_name="db",
    **kwargs,
):
    """Run a full create_db() over a minimal non-segmented input set.

    Returns (sqlite3.Connection, captured_stdout). stdout is captured because
    several assertions here are about the *absence* of any warning.
    """
    work = tmp_path / db_name
    work.mkdir(parents=True, exist_ok=True)
    p = lambda name: str(work / name)

    accs = [row[0] for row in meta_rows]
    _write_tsv(p("meta.tsv"), meta_rows, list(meta_cols))
    _write_tsv(p("features.tsv"),
               feat_rows if feat_rows is not None else [[a, "polyprotein"] for a in accs],
               list(feat_cols))
    _write_tsv(p("aln.tsv"),
               aln_rows if aln_rows is not None else [[a, accs[0], "ATGC"] for a in accs],
               list(aln_cols))
    _write_tsv(p("gene.tsv"), [["polyprotein", "Polyprotein"]], ["name", "description"])
    _write_csv(p("m49_country.csv"), [["001", "World"]], ["m49_code", "name"])
    _write_csv(p("m49_inter.csv"), [["X", "Inter"]], ["code", "name"])
    _write_csv(p("m49_region.csv"), [["Y", "Region"]], ["code", "name"])
    _write_csv(p("m49_sub.csv"), [["Z", "Sub"]], ["code", "name"])
    _write_tsv(p("proj.tsv"), [["Python", "3"]], ["Software", "Version"])
    _write_tsv(p("ins.tsv"),
               ins_rows if ins_rows is not None else [[accs[0], "none"]],
               list(ins_cols))
    _write_tsv(p("host.tsv"), [[accs[0], "Homo sapiens"]], ["primary_accession", "host"])

    with open(p("seqs.fa"), "w", encoding="utf-8") as handle:
        for acc in (fasta_ids if fasta_ids is not None else accs):
            handle.write(f">{acc}\nATGCATGCATGC\n")

    reference_tsv = None
    if ref_rows is not None:
        reference_tsv = p("ref_list.txt")
        with open(reference_tsv, "w", encoding="utf-8") as handle:
            for row in ref_rows:
                handle.write("\t".join(row) + "\n")

    creator = CreateSqliteDB(
        meta_data=p("meta.tsv"),
        features=p("features.tsv"),
        pad_aln=p("aln.tsv"),
        gene_info=p("gene.tsv"),
        m49_countries=p("m49_country.csv"),
        m49_interm_region=p("m49_inter.csv"),
        m49_regions=p("m49_region.csv"),
        m49_sub_regions=p("m49_sub.csv"),
        proj_settings=p("proj.tsv"),
        fasta_sequence_file=p("seqs.fa"),
        insertions=p("ins.tsv"),
        host_taxa_file=p("host.tsv"),
        base_dir=str(work),
        output_dir="SqliteDB",
        db_name=db_name,
        db_status="new db",
        reference_tsv=reference_tsv,
        **kwargs,
    )

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        creator.create_db()

    conn = sqlite3.connect(str(work / "SqliteDB" / f"{db_name}.db"))
    return conn, buffer.getvalue()


def _creator_stub(tmp_path):
    """A CreateSqliteDB instance whose constructor touches no files.

    Used for unit-level calls to the segment-inference helpers. Built with
    keyword arguments only so it survives signature growth.
    """
    placeholder = str(tmp_path / "unused.tsv")
    return CreateSqliteDB(
        meta_data=placeholder,
        features=placeholder,
        pad_aln=placeholder,
        gene_info=None,
        m49_countries=placeholder,
        m49_interm_region=placeholder,
        m49_regions=placeholder,
        m49_sub_regions=placeholder,
        proj_settings=placeholder,
        fasta_sequence_file=placeholder,
        insertions=placeholder,
        host_taxa_file=placeholder,
        base_dir=str(tmp_path),
        output_dir="SqliteDB",
        db_name="stub",
        db_status="new db",
    )


def _rows(conn, sql):
    return conn.execute(sql).fetchall()


# --------------------------------------------------------------------------
# 1. reference-list segment spellings feed an un-normalised backfill
# --------------------------------------------------------------------------

def test_canonical_reference_segment_one_keeps_every_child_row(tmp_path):
    """Control: the shipped HCV/RABV reference lists write segment ``1``.

    generic/hcv/ref_list_subtype_genotype.txt and generic/rabv/ref_list_clades.txt
    both carry a literal ``1`` in the segment column, so this is the path every
    real non-segmented build takes. It must keep every child row and stamp
    segment '1' throughout - the rest of this section is what happens when the
    same value is spelled differently.
    """
    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["Q1", "query", ""], ["R1", "reference", ""]],
        ref_rows=[["R1", "reference", "1", "1", "a"]],
        db_name="canon",
    )
    assert _rows(conn, "SELECT primary_accession, segment FROM meta_data ORDER BY 1") == [
        ("Q1", "1"), ("R1", "1")]
    assert _rows(conn, "SELECT accession, segment FROM features ORDER BY 1") == [
        ("Q1", "1"), ("R1", "1")]
    assert _rows(conn, "SELECT header FROM sequences ORDER BY 1") == [("Q1",), ("R1",)]


def test_whitespace_padded_reference_segment_is_tolerated(tmp_path):
    """A hand-edited reference list with ``' 1 '`` must behave like ``'1'``.

    Reference lists are curated by hand and by spreadsheet export; stray spaces
    around the segment column are the most common cosmetic difference. This one
    IS handled (the backfill strips) and is pinned so the fix for the sibling
    cases below does not regress it.
    """
    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["Q1", "query", ""], ["R1", "reference", ""]],
        ref_rows=[["R1", "reference", " 1 ", "1", "a"]],
        db_name="ws",
    )
    assert _rows(conn, "SELECT primary_accession, segment FROM meta_data ORDER BY 1") == [
        ("Q1", "1"), ("R1", "1")]
    assert _rows(conn, "SELECT COUNT(*) FROM features")[0][0] == 2


# FIXED: create_db's reference-segment backfill now uses the shared
# normaliser, so this is a regression guard rather than a known bug.
def test_zero_padded_reference_segment_does_not_delete_the_reference(tmp_path):
    """A zero-padded ``01`` in the reference list erases a reference's data.

    Real-world trigger: a curator writes the segment column as ``01`` (or a tool
    zero-pads it). segment_utils exists precisely because '04' and '4' are the
    same segment, but CreateSqliteDB's reference backfill predates it and keeps
    the raw digits. The reference survives in meta_data, so nothing looks wrong -
    but its CDS features, its alignment row and its sequence are all gone, and a
    reference with no sequence silently degrades every downstream lookup.
    """
    conn, output = _build_db(
        tmp_path,
        meta_rows=[["Q1", "query", ""], ["R1", "reference", ""]],
        ref_rows=[["R1", "reference", "01", "1", "a"]],
        db_name="zeropad",
    )
    assert "R1" in {r[0] for r in _rows(conn, "SELECT accession FROM features")}, (
        f"reference features silently dropped; build said nothing:\n{output}")
    assert "R1" in {r[0] for r in _rows(conn, "SELECT header FROM sequences")}


# FIXED: create_db's reference-segment backfill now uses the shared
# normaliser, so this is a regression guard rather than a known bug.
def test_float_spelled_reference_segment_is_not_read_as_segment_ten(tmp_path):
    """``1.0`` in the reference list must mean segment 1, not segment 10.

    Real-world trigger: the reference list is opened in Excel or round-tripped
    through pandas without dtype=str, which renders an integer column with any
    blank cell as ``1.0``. segment_utils documents this exact corruption ('4.0'
    -> 40) as one of the bugs it was created to end; this call site never adopted
    it. The resulting DB has meta_data.segment '10' for references and '' for
    every query - both wrong, both silent.
    """
    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["Q1", "query", ""], ["R1", "reference", ""]],
        ref_rows=[["R1", "reference", "1.0", "1", "a"]],
        db_name="floatseg",
    )
    segments = dict(_rows(conn, "SELECT primary_accession, segment FROM meta_data"))
    assert segments == {"Q1": "1", "R1": "1"}, segments


# --------------------------------------------------------------------------
# 2. the unsegmented inference itself
# --------------------------------------------------------------------------

def test_empty_observed_segment_set_forces_segment_one(tmp_path):
    """No row anywhere carries a segment -> the build calls it unsegmented.

    This is load-bearing, not incidental: generic/hcv/ref_list.txt,
    generic/rabv/ref_list.txt and generic/other/ref_list-hpv.txt all leave the
    segment column empty, so a real HCV/RABV build reaches
    _should_force_unsegmented_segment_one with an *empty* observed set and relies
    on ``set().issubset({'1'})`` being True to get segment '1' into the DB.
    Pinned so nobody "fixes" the empty case without noticing who depends on it.
    """
    creator = _creator_stub(tmp_path)
    empty = pd.DataFrame({"primary_accession": ["A"], "segment": [""]})
    assert creator._should_force_unsegmented_segment_one(None, [empty]) is True
    assert creator._should_force_unsegmented_segment_one(None, [pd.DataFrame()]) is True

    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["Q1", "query", ""], ["R1", "reference", ""]],
        ref_rows=[["R1", "reference", "", "1", "a"]],
        db_name="allblank",
    )
    assert _rows(conn, "SELECT DISTINCT segment FROM meta_data") == [("1",)]
    assert _rows(conn, "SELECT DISTINCT segment FROM features") == [("1",)]


@pytest.mark.xfail(
    reason="When no explicit is_segmented flag is passed, _should_force_unsegmented_"
           "segment_one still decides 'this virus is unsegmented' purely from the segments "
           "observed in THIS batch, so a segmented build whose batch happens to hold only "
           "segment-1 rows invents segment '1' for rows that legitimately have none. The "
           "reference list is an unambiguous signal that the virus is segmented (here it "
           "declares segments 1,2,4,6,8) and is already loaded into ref_to_seg, but the "
           "inference fallback never consults it.",
    strict=False,
)
def test_segmented_reference_list_prevents_inventing_segment_one(tmp_path):
    """A batch of segment-1 influenza rows must not make blanks into PB2.

    Real-world trigger: a test=1 subset, a re-run restricted to one segment, or a
    batch in which segment validation failed for the rest. Segment 1 is PB2, so
    stamping '1' on a row whose segment is unknown asserts a specific, checkable,
    wrong fact about it - and it is written to meta_data, features, sequences and
    the alignment at once.

    An explicit is_segmented flag now settles this for callers that pass one;
    this case is the residual inference fallback, which is what any direct
    invocation of CreateSqliteDB.py without --is_segmented still gets. The
    reference list would answer the question here without any new flag.
    """
    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["S1", "reference", ""], ["Q1", "query", "1"], ["Q2", "query", ""]],
        aln_rows=[["Q1", "S1", "ATGC"], ["Q2", "S1", "ATGC"]],
        ref_rows=[["S1", "reference", "1", "", ""],
                  ["S2", "reference", "2", "", ""],
                  ["S4", "reference", "4", "", ""],
                  ["S6", "reference", "6", "", ""],
                  ["S8", "reference", "8", "", ""]],
        db_name="segmented_subset",
    )
    segments = dict(_rows(conn, "SELECT primary_accession, segment FROM meta_data"))
    assert segments["Q2"] == "", (
        f"Q2 had no segment and the reference list says this virus is segmented: {segments}")


# --------------------------------------------------------------------------
# 3. valid_pairs: the unreported inner join
# --------------------------------------------------------------------------

@pytest.mark.xfail(
    reason="valid_pairs is built from meta_data BEFORE _force_segment_one_df runs, so a "
           "features table that already carries the explicit segment '1' cannot match a "
           "meta_data matrix whose segment column is blank. Every features row is dropped "
           "and the resulting features table is empty, while the build summary still "
           "reports the sequences as passing QC.",
    strict=False,
)
def test_features_segment_one_survives_blank_metadata_segment(tmp_path):
    """features says segment 1, the matrix says nothing - and features empties out.

    Real-world trigger: for an unsegmented virus the gB matrix leaves segment
    blank (the reference list has no segment column), while the feature
    projection stamps segment 1 from the reference it projected against. The two
    spellings mean the same thing and the build reconciles them minutes later by
    forcing '1' everywhere - but by then valid_pairs has already discarded the
    rows. Nothing is printed; the operator sees "Query sequences passing QC: 2"
    over an empty features table.
    """
    conn, output = _build_db(
        tmp_path,
        meta_rows=[["Q1", "query", ""], ["Q2", "query", ""]],
        feat_rows=[["Q1", "1", "polyprotein"], ["Q2", "1", "polyprotein"]],
        feat_cols=("accession", "segment", "product"),
        db_name="featseg",
    )
    count = _rows(conn, "SELECT COUNT(*) FROM features")[0][0]
    assert count == 2, f"features table emptied silently; build output was:\n{output}"


@pytest.mark.xfail(
    reason="The child-table filters are guarded by `if valid_pairs:`. When every meta_data "
           "row is excluded the set is empty and the guard is read as 'nothing to filter' "
           "rather than 'nothing may pass', so features/sequence_alignment/sequences/"
           "insertions are written in full for accessions the same build marked excluded. "
           "One surviving accession changes the outcome for all the others.",
    strict=False,
)
def test_all_excluded_batch_does_not_publish_excluded_child_rows(tmp_path):
    """A batch where everything failed QC must not ship its sequences anyway.

    Real-world trigger: a small update batch, or a run against the wrong
    reference, in which every query trips alignment filtering. The exclusion
    bookkeeping in meta_data is correct - exclusion_status '1' on both rows - but
    the sequences, alignments and CDS features of those excluded records are
    published, which is the opposite of what the same code does when a single
    record passes.
    """
    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["Q1", "query", "", "missing fasta sequence"],
                   ["Q2", "query", "", "missing fasta sequence"]],
        meta_cols=("primary_accession", "accession_type", "segment", "exclusion"),
        feat_rows=[["Q1", "polyprotein"], ["Q2", "polyprotein"]],
        aln_rows=[["Q1", "Q1", "ATGC"], ["Q2", "Q1", "ATGC"]],
        ins_rows=[["Q1", "x"], ["Q2", "x"]],
        db_name="allexcl",
    )
    assert _rows(conn, "SELECT DISTINCT exclusion_status FROM meta_data") == [("1",)]
    assert _rows(conn, "SELECT COUNT(*) FROM features")[0][0] == 0
    assert _rows(conn, "SELECT COUNT(*) FROM sequences")[0][0] == 0


def test_partially_excluded_batch_drops_only_the_excluded_child_rows(tmp_path):
    """The non-degenerate case, pinned as the intended contract.

    With at least one accession passing QC, valid_pairs is non-empty and the
    excluded accession's features/sequences/alignment rows are correctly removed
    while the passing one's are kept. This is the behaviour the all-excluded case
    above should also produce.
    """
    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["Q1", "query", "", ""],
                   ["Q2", "query", "", "missing fasta sequence"]],
        meta_cols=("primary_accession", "accession_type", "segment", "exclusion"),
        feat_rows=[["Q1", "polyprotein"], ["Q2", "polyprotein"]],
        aln_rows=[["Q1", "Q1", "ATGC"], ["Q2", "Q1", "ATGC"]],
        ins_rows=[["Q1", "x"], ["Q2", "x"]],
        db_name="partexcl",
    )
    assert _rows(conn, "SELECT accession FROM features") == [("Q1",)]
    assert _rows(conn, "SELECT header FROM sequences") == [("Q1",)]


@pytest.mark.xfail(
    reason="valid_pairs matches accessions with a case-sensitive string comparison and "
           "reports nothing when a row fails to match, so a features file that spells an "
           "accession in lower case loses every feature for that accession without a "
           "single line of output.",
    strict=False,
)
def test_lowercase_feature_accession_is_not_silently_discarded(tmp_path):
    """Case-only differences in an accession should not delete a record's features.

    Real-world trigger: a GFF or supplier feature table whose seqid column is
    lower case, while the GenBank matrix uses the canonical upper case. INSDC
    accessions are conventionally upper case but not case-significant, so the two
    files describe the same sequence. If dropping is the intended answer it still
    needs to be counted and printed - a silent zero is indistinguishable from
    "this record genuinely has no CDS".
    """
    conn, output = _build_db(
        tmp_path,
        meta_rows=[["AB123", "query", ""], ["CD456", "query", ""]],
        feat_rows=[["ab123", "polyprotein"], ["CD456", "polyprotein"]],
        db_name="acccase",
    )
    kept = {r[0].upper() for r in _rows(conn, "SELECT accession FROM features")}
    assert kept == {"AB123", "CD456"}, (
        f"lost features for AB123 with no diagnostic; output was:\n{output}")


# --------------------------------------------------------------------------
# 4. HCV genotype / subtype assignment
# --------------------------------------------------------------------------

@pytest.mark.xfail(
    reason="_add_reference_columns treats a reference's own blank subtype as 'no direct "
           "information' rather than as the curated statement 'subtype unknown', so "
           "`fallback_subtype.where(direct_subtype == '', direct_subtype)` lets the BLAST/"
           "tree fallback overwrite it with the subtype of a DIFFERENT reference - while "
           "the genotype stays the reference's own. The stored pair is then a genotype "
           "from one sequence and a subtype from another.",
    strict=False,
)
def test_reference_with_curated_unknown_subtype_keeps_it_blank(tmp_path):
    """A reference curated as subtype NA must not be given someone else's subtype.

    Real-world trigger: generic/hcv/ref_list_subtype_genotype.txt records subtype
    'NA' for 69 of its 238 references - genotype known, subtype deliberately not
    assigned. Here R2 is curated genotype 6 / subtype NA and is BLAST-aligned
    against R1 (genotype 1, subtype a); it comes out of the build as genotype 6
    with subtype 'a', i.e. the specific claim "6a" about a sequence whose subtype
    the curators refused to call. Nothing in the DB distinguishes that fabricated
    letter from a curated one.
    """
    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["R1", "reference", ""], ["R2", "reference", ""], ["Q1", "query", ""]],
        aln_rows=[["R2", "R1", "ATGC"], ["Q1", "R1", "ATGC"]],
        ref_rows=[["R1", "reference", "1", "1", "a"],
                  ["R2", "reference", "1", "6", "NA"]],
        db_name="natype",
    )
    labels = {
        acc: (gt, st) for acc, gt, st in _rows(
            conn,
            "SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype "
            "FROM meta_data")
    }
    assert labels["R2"] == ("6", ""), labels


@pytest.mark.xfail(
    reason="_tree_based_reference_labels merges genotype and subtype independently across "
           "the candidate trees ('the first tree to label an accession wins' is applied per "
           "field), so a query whose UShER neighbourhood supplies a genotype but no subtype "
           "takes its subtype from a different tree's neighbourhood. The stored genotype/"
           "subtype pair can be one no reference in either tree carries.",
    strict=False,
)
def test_query_labels_are_not_assembled_from_two_different_trees(tmp_path):
    """genotype from the UShER tree + subtype from the IQ-TREE = a chimeric label.

    Real-world trigger: HCV builds carry both trees, over different tip sets (the
    IQ-TREE backbone holds cluster representatives only), and 69 of the HCV
    references have no subtype at all - so "UShER gives a genotype but no
    subtype" is the normal case, not a corner. Here Q sits beside genotype-1
    references with no subtype in the UShER tree and beside genotype-2c
    references in the IQ-TREE, and is stored as genotype 1, subtype c - "1c",
    which neither tree's neighbourhood supports.
    """
    work = tmp_path / "chimera"
    work.mkdir()
    usher = work / "usher.nh"
    usher.write_text("((Q:0.1,(R1:0.1,R2:0.1):0.1):0.1,(X1:0.5,X2:0.5):0.5);", encoding="utf-8")
    iqtree = work / "iqtree.treefile"
    iqtree.write_text("((Q:0.1,(C1:0.1,C2:0.1):0.1):0.1,(Y1:0.5,Y2:0.5):0.5);", encoding="utf-8")

    accs = ["Q", "R1", "R2", "C1", "C2"]
    conn, _ = _build_db(
        tmp_path,
        meta_rows=[[a, "query" if a == "Q" else "reference", ""] for a in accs],
        aln_rows=[["Q", "R1", "ATGC"]],
        ref_rows=[["R1", "reference", "1", "1", "NA"],
                  ["R2", "reference", "1", "1", "NA"],
                  ["C1", "reference", "1", "2", "c"],
                  ["C2", "reference", "1", "2", "c"]],
        db_name="chimera_db",
        usher_tree=str(usher),
        iqtree_file=str(iqtree),
    )
    genotype, subtype = _rows(
        conn,
        "SELECT nearest_reference_genotype, nearest_reference_subtype "
        "FROM meta_data WHERE primary_accession='Q'")[0]
    assert (genotype, subtype) in {("1", ""), ("2", "c")}, (
        f"Q was labelled genotype {genotype!r} + subtype {subtype!r}, a pair no "
        "reference in either tree carries")


def test_duplicate_reference_rows_resolve_to_the_first_occurrence(tmp_path):
    """Duplicated accessions in a reference list are first-wins, deterministically.

    Real-world trigger: generic/hcv/ref_list_subtype_genotype.txt has 437 lines
    for 238 accessions - every entry is listed twice. Those duplicates are
    byte-identical so the collapse is harmless, which is exactly why a
    *conflicting* duplicate would be easy to introduce and impossible to notice:
    _load_reference_lookup drops it with keep='first' and prints nothing. This
    pins the resolution order so it stays deterministic (file order), and
    documents that a conflict is not reported.
    """
    conn, output = _build_db(
        tmp_path,
        meta_rows=[["R1", "reference", ""], ["R9", "reference", ""], ["Q1", "query", ""]],
        aln_rows=[["Q1", "R1", "ATGC"]],
        ref_rows=[["R1", "reference", "1", "1", "a"],
                  ["R1", "reference", "1", "6", "z"],
                  ["R9", "reference", "1", "2", "b"]],
        db_name="dupref",
    )
    labels = {
        acc: (gt, st) for acc, gt, st in _rows(
            conn,
            "SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype "
            "FROM meta_data")
    }
    assert labels["R1"] == ("1", "a")
    assert "conflict" not in output.lower()


def test_reference_absent_from_every_tree_keeps_its_curated_labels(tmp_path):
    """A labelled reference missing from all trees still gets its own labels.

    Real-world trigger: references are routinely absent from the UShER/IQ-TREE
    tips (deduplication, cluster representatives, tree_free mode). The direct
    lookup must win regardless, or a reference's curated genotype would depend on
    whether it happened to be sampled into a tree.
    """
    work = tmp_path / "notree"
    work.mkdir()
    tree = work / "iqtree.treefile"
    tree.write_text("(Q1:0.1,R1:0.1);", encoding="utf-8")

    conn, _ = _build_db(
        tmp_path,
        meta_rows=[["R1", "reference", ""], ["R2", "reference", ""], ["Q1", "query", ""]],
        aln_rows=[["Q1", "R1", "ATGC"]],
        ref_rows=[["R1", "reference", "1", "1", "a"],
                  ["R2", "reference", "1", "4", "d"]],
        db_name="notree_db",
        iqtree_file=str(tree),
    )
    labels = {
        acc: (gt, st) for acc, gt, st in _rows(
            conn,
            "SELECT primary_accession, nearest_reference_genotype, nearest_reference_subtype "
            "FROM meta_data")
    }
    assert labels["R2"] == ("4", "d")


# --------------------------------------------------------------------------
# 5. the shipped HCV database
# --------------------------------------------------------------------------

@pytest.mark.skipif(not HCV_DB.exists(), reason=f"HCV test DB not present: {HCV_DB}")
def test_hcv_db_child_tables_agree_with_meta_data_on_segment(tmp_path):
    """Every child row in the shipped HCV DB joins back to meta_data.

    This is the invariant valid_pairs is supposed to guarantee. It holds in
    HCV_OM_test.db (all segments are '1'), and pinning it here means the
    zero-padded / '1.0' reference-list cases above are demonstrably departures
    from a real build rather than theoretical ones.
    """
    conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
    try:
        pairs = set(_rows(conn, "SELECT primary_accession, segment FROM meta_data"))
        for table, key in [("features", "accession"),
                           ("sequences", "header"),
                           ("sequence_alignment", "primary_accession")]:
            orphans = [
                row for row in _rows(conn, f"SELECT DISTINCT {key}, segment FROM {table}")
                if row not in pairs
            ]
            assert not orphans, f"{table} rows with no matching meta_data pair: {orphans[:5]}"
    finally:
        conn.close()


@pytest.mark.skipif(
    not HCV_DB.exists() or not HCV_REF_LIST.exists(),
    reason="HCV test DB or generic/hcv reference list not present",
)
@pytest.mark.xfail(
    reason="Shipped evidence for the curated-unknown-subtype bug: 69 references whose "
           "reference-list subtype is 'NA' are stored in HCV_OM_test.db with an inferred "
           "subtype letter taken from the reference they aligned against.",
    strict=False,
)
def test_hcv_db_preserves_curated_unknown_subtypes():
    """The real database already carries the fabricated subtypes.

    Real-world trigger: this is not a constructed input - it is the shipped
    HCV_OM_test.db built from generic/hcv/ref_list_subtype_genotype.txt. Every
    one of the 69 'NA'-subtype references comes out with a letter, so a consumer
    reading nearest_reference_subtype for a reference cannot tell a curated
    subtype from a guess.
    """
    curated_unknown = set()
    with open(HCV_REF_LIST, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[0] and parts[4].strip().upper() in {"NA", ""}:
                curated_unknown.add(parts[0].strip())

    conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
    try:
        rows = _rows(
            conn,
            "SELECT primary_accession, nearest_reference_subtype FROM meta_data "
            "WHERE accession_type IN ('reference', 'master')")
    finally:
        conn.close()

    invented = sorted(acc for acc, subtype in rows
                      if acc in curated_unknown and (subtype or "").strip() != "")
    assert not invented, (
        f"{len(invented)} references curated as subtype NA carry an inferred subtype, "
        f"e.g. {invented[:5]}")
