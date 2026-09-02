"""Genotype scoping, per-genotype wild-type suppression and residue spelling.

Lens: scripts/AnnotateMutations.py, scripts/NormalizeHcvMutationCatalog.py and
the curated PHDR catalogue at
generic/hcv/Tables/generalized_mutation_catalog_with_extra_info.tsv.

Three separate pieces of information exist in the catalogue and were being
thrown away by the annotator:

  * ``alignment_name`` - the genotype/subtype bucket PHDR scored a finding in;
  * ``display_structure`` - the wild-type residue(s) for that bucket, written
    out in full ("R155C", "K/Q80K", "32del");
  * the deletion findings, whose alt residue is spelled ``del`` and can
    therefore never equal a translated codon.

The tests below pin the rules that put them back to work: a call is emitted
only when the sequence's genotype is in scope for the entry AND the observed
residue differs from that genotype's wild type.

Real-DB assertions open the database read-only and skip when it is absent.
Every synthetic case writes only into ``tmp_path``.
"""

import csv
import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

import AnnotateMutations as AM
import BuildCatalogGenotypeColumns as BCG
import NormalizeHcvMutationCatalog as NM


REPO_ROOT = Path(__file__).resolve().parents[2]
HCV_TABLES = REPO_ROOT / "generic" / "hcv" / "Tables"
CATALOG_TSV = HCV_TABLES / "generalized_mutation_catalog_with_extra_info.tsv"
GENE_INFO_TSV = HCV_TABLES / "gene_info.tsv"
HCV_DB = REPO_ROOT / "test_out" / "HCV_OM_test" / "HCV_OM_test.db"

requires_hcv_db = pytest.mark.skipif(
    not HCV_DB.exists(), reason=f"HCV test database not present at {HCV_DB}"
)
requires_hcv_assets = pytest.mark.skipif(
    not CATALOG_TSV.exists(), reason=f"HCV mutation catalog not present at {CATALOG_TSV}"
)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _read_catalog():
    return pd.read_csv(CATALOG_TSV, sep="\t", dtype=str)


def _open_hcv_db():
    return sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)


def _catalog_row(**overrides):
    """A catalog row carrying the genotype scope and wild-type columns too."""
    row = {
        "mutation_id": "NS3:6R",
        "protein_name": "NS3",
        "segment": "1",
        "aa_position": "6",
        "alt_residue": "R",
        "reference_accession": "REF1",
        "mutation_type": "aminoAcidSimplePolymorphism",
        "signature_id": "NS3:6R",
        "signature_kind": "single",
        "combination_id": "",
        "combination_size": "",
        "phenotype": "",
        "alignment_name": "AL_1a",
        "display_structure": "K6R",
        "relevant_genotypes": "1a",
    }
    row.update(overrides)
    # The pipeline reads only the generic columns, so a test row must carry the
    # wild type the way a real catalog does. scripts/BuildCatalogGenotypeColumns.py
    # is what derives this from the PHDR tables at build time; here we do the
    # equivalent from display_structure so the fixtures stay readable.
    if "wild_type_residues" not in overrides:
        row["wild_type_residues"] = _wild_type_column(
            row.get("display_structure", ""), row.get("alignment_name", ""), row.get("aa_position", "")
        )
    return row


_COMPONENT = re.compile(r'^(?P<wt>[A-Za-z](?:/[A-Za-z])*)?(?P<pos>\d+)(?P<alt>[A-Za-z]|del|-)$')


def _wild_type_column(display_structure, alignment_name, position):
    """Build a ``wild_type_residues`` value from a display_structure component."""
    code = BCG.alignment_to_genotype_code(alignment_name)
    if not code or not display_structure:
        return ""
    for component in str(display_structure).split("+"):
        match = _COMPONENT.match(component.strip())
        if match and match.group("pos") == str(position) and match.group("wt"):
            return f"{code}:{match.group('wt').split('/')[0]}"
    return ""


def _annotate(catalog_rows, alignments, feature_map, meta_rows, master="REF1"):
    """Run the annotator over a tiny synthetic alignment.

    ``alignments`` is a list of (accession, padded_alignment) pairs; the first
    is treated as the master reference.  ``meta_rows`` supplies the per-sequence
    genotype the scope gate reads.  Returns (emitted, evidence, diagnostics).
    """
    seq_aln = pd.DataFrame(
        [
            {
                "sequence_id": accession,
                "primary_accession": accession,
                "alignment": alignment,
                "alignment_name": alignments[0][0],
            }
            for accession, alignment in alignments
        ]
    )
    meta = pd.DataFrame(meta_rows)
    prepared, _invalid = AM.prepare_catalog(pd.DataFrame(catalog_rows), {})
    evidence = []
    found, diagnostics, _maps = AM.annotate_from_reference_coordinates(
        prepared, seq_aln, meta, {}, {master: feature_map}, False, call_evidence=evidence
    )
    return found, evidence, dict(diagnostics)


def _meta(*rows):
    """meta_data rows: (accession, genotype, subtype); the master is first."""
    records = []
    for index, (accession, genotype, subtype) in enumerate(rows):
        records.append(
            {
                "primary_accession": accession,
                "accession_type": "master" if index == 0 else "",
                "nearest_reference_genotype": genotype,
                "nearest_reference_subtype": subtype,
            }
        )
    return records


NS3_FEATURE_MAP = {
    "NS3": {"product": "NS3", "cds_start": 1, "cds_end": 30, "feature_type": "mat_peptide"}
}
# codon 1 ATG=M, 2 AAA=K, 3 CCC=P, 4 GGG=G, 5 TTT=F, 6 AGA=R, 7 TAG=stop
REFERENCE_30MER = "ATGAAACCCGGGTTTAGATAGCATGCATGC"


# --------------------------------------------------------------------------
# 1. Residue vocabulary: one spelling per concept
# --------------------------------------------------------------------------

def test_codon_table_spells_stop_as_a_star():
    """Stop is '*' in the genetic code table, not the private '_' spelling.

    The normalizer's token grammar already admits '*' (``[A-Z*]|del``), so the
    day PHDR publishes a nonsense variant the two halves of the pipeline have to
    agree on how it is written.  '_' met nothing on the catalog side, so such a
    row could only ever match zero sequences, silently.
    """
    assert AM.CODON_TABLE["TAA"] == "*"
    assert AM.CODON_TABLE["TAG"] == "*"
    assert AM.CODON_TABLE["TGA"] == "*"
    assert AM.translate_codon("TAA") == "*"
    assert "_" not in set(AM.CODON_TABLE.values())


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("del", "-"),
        ("DEL", "-"),
        (" del ", "-"),
        ("_", "*"),
        ("*", "*"),
        ("r", "R"),
        ("R ", "R"),
        (" R", "R"),
        ("", ""),
        (None, ""),
        (float("nan"), ""),
    ],
)
def test_normalize_residue_accepts_legacy_spellings_and_writes_the_standard_one(raw, expected):
    """Legacy 'del'/'_' must still load; the standard '-'/'*' is what comes out.

    The catalogue is a hand-maintained TSV regenerated from curator
    spreadsheets, so lower case and a stray space have to survive too - the old
    raw comparison turned either into a row that matched nothing at all with no
    diagnostic recorded.  A missing cell arrives as float NaN, which is truthy:
    ``str(value or '')`` would turn it into the literal residue 'NAN'.
    """
    assert AM.normalize_residue(raw) == expected


def test_deletion_is_detected_from_the_gap_not_from_a_translated_codon():
    """A deleted codon has no codon to translate; the gap itself is the evidence.

    translate_codon strips '-' and then fails its length guard, so a fully
    gapped codon can only come back 'X'.  The catalogue spells its twelve
    deletion findings 'del', which no translation can ever equal - so
    NS5A 29/30/32del were undetectable by construction.
    """
    assert AM.residue_from_aligned_codon("---") == "-"
    assert AM.residue_from_aligned_codon("NNN") == "X"
    assert AM.residue_from_aligned_codon("A-G") == "X"
    assert AM.residue_from_aligned_codon("ATG") == "M"
    assert AM.residue_from_aligned_codon(None) is None


def test_terminal_padding_is_missing_data_not_a_deletion():
    """Gaps outside the sequence's covered span are padding, never a deletion.

    Partial GenBank records are padded out to the full genome width.  In the
    shipped HCV build that padding covers NS5A 29/30/32 in 35 sequences; reading
    it as a deletion would convert "never sequenced" into a resistance call on
    exactly the worst-covered records.
    """
    alignment = "---ATGAAA---CCCGGG---"
    span = AM.alignment_covered_span(alignment)
    assert span == (3, 17)
    # columns 9-11 are an internal gap between covered bases -> a deletion
    assert AM.residue_from_aligned_codon("---", span, (9, 10, 11)) == "-"
    # columns 0-2 and 18-20 are terminal padding -> undecidable, not a deletion
    assert AM.residue_from_aligned_codon("---", span, (0, 1, 2)) == "X"
    assert AM.residue_from_aligned_codon("---", span, (18, 19, 20)) == "X"
    assert AM.alignment_covered_span("-----") is None


def test_catalogued_deletion_is_annotated_when_the_codon_is_deleted():
    """An end-to-end deletion call, the case the 'del' spelling made impossible."""
    deleted = "ATGAAACCCGGGTTT---TAGCATGCATGC"   # codon 6 removed, internal
    catalog = [_catalog_row(alt_residue="del", mutation_id="NS3:6del",
                            signature_id="NS3:6del", display_structure="6del")]
    found, _evidence, _diag = _annotate(
        catalog,
        [("REF1", REFERENCE_30MER), ("DEL1", deleted)],
        NS3_FEATURE_MAP,
        _meta(("REF1", "1", "a"), ("DEL1", "1", "a")),
    )
    assert {hit["primary_accession"] for hit in found} == {"DEL1"}
    assert found[0]["alt_residue"] == "-"
    assert found[0]["observed_residue"] == "-"


# --------------------------------------------------------------------------
# 2. Genotype codes and the scope gate (rule A)
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "raw,expected",
    [("AL_1a", "1a"), ("AL_6", "6"), ("al_3a", "3a"), ("1b", "1b"), ("", ""), (None, ""),
     (float("nan"), "")],
)
def test_alignment_name_to_genotype_code(raw, expected):
    assert BCG.alignment_to_genotype_code(raw) == expected


@pytest.mark.parametrize(
    "code,expected",
    [("1a", "1"), ("6xd", "6"), ("3", "3"), ("10a", "10"), ("", ""), (float("nan"), "")],
)
def test_genotype_of_code_takes_the_leading_digits(code, expected):
    assert AM.genotype_of_code(code) == expected


@pytest.mark.parametrize(
    "genotype,subtype,expected",
    [("1", "a", "1a"), ("6", "xd", "6xd"), ("1", "1a", "1a"), ("1", "", "1"),
     ("", "1a", "1a"), ("1", "NA", "1NA")],
)
def test_build_subtype_code_joins_meta_data_halves(genotype, subtype, expected):
    """meta_data stores genotype and subtype apart; the catalog spells them together."""
    assert AM.build_subtype_code(genotype, subtype) == expected


def test_scope_gate_matches_at_genotype_level_not_subtype_level():
    """6n must match a finding scored in AL_6a, because 6n has no bucket of its own.

    Observed subtypes in the shipped build include 6n, 6xd, 6xc and 4v, none of
    which PHDR scores anything in.  Requiring an exact subtype match would throw
    away 605 of the 803 in-scope calls and leave those sequences unannotated for
    reasons that have nothing to do with their biology.
    """
    scope = frozenset({"6a", "6l"})
    assert AM.classify_genotype_scope(scope, "6", "6a") == AM.SCOPE_TIER_SUBTYPE
    assert AM.classify_genotype_scope(scope, "6", "6n") == AM.SCOPE_TIER_GENOTYPE
    assert AM.classify_genotype_scope(scope, "6", "6xd") == AM.SCOPE_TIER_GENOTYPE
    assert AM.classify_genotype_scope(scope, "1", "1a") == AM.SCOPE_TIER_OUT_OF_SCOPE


def test_an_entry_with_no_genotype_scope_applies_to_any_genotype():
    """The 63 signatures with no alignment bucket are anchors, not genotype-1 findings.

    They exist because the normalizer promotes conjunction components with an
    empty phdr_ras_id to standalone signatures.  Gating them to nothing would
    hide them; gating them to a genotype would invent one.  Rule B is what
    filters them, not the scope gate.
    """
    assert AM.classify_genotype_scope(frozenset(), "3", "3a") == AM.SCOPE_TIER_UNSCOPED
    assert AM.classify_genotype_scope(frozenset(), "", "") == AM.SCOPE_TIER_UNSCOPED


def test_a_sequence_with_no_genotype_cannot_satisfy_a_scoped_entry():
    """An unclassified sequence is out of scope for a genotype-scoped finding.

    Recorded rather than silently dropped: the call lands in
    sequence_mutation_calls with call_status 'suppressed_out_of_scope' so the
    reason is visible, and diagnostics counts the sequences it happened to.
    """
    assert AM.classify_genotype_scope(frozenset({"1a"}), "", "") == AM.SCOPE_TIER_OUT_OF_SCOPE

    found, evidence, diagnostics = _annotate(
        [_catalog_row()],
        [("REF1", REFERENCE_30MER), ("Q1", REFERENCE_30MER)],
        NS3_FEATURE_MAP,
        _meta(("REF1", "1", "a"), ("Q1", "", "")),
    )
    assert "Q1" not in {hit["primary_accession"] for hit in found}
    suppressed = [row for row in evidence if row["primary_accession"] == "Q1"]
    assert suppressed and all(
        row["call_status"] == AM.CALL_STATUS_SUPPRESSED_OUT_OF_SCOPE for row in suppressed
    )
    assert diagnostics["sequences_without_genotype"] == 1


def test_out_of_scope_call_is_suppressed_and_in_scope_call_survives():
    """A 1a-only finding must not be scored against a genotype 3 sequence."""
    found, evidence, _diag = _annotate(
        [_catalog_row(display_structure="K6R")],
        [("REF1", REFERENCE_30MER), ("G1A", REFERENCE_30MER), ("G3A", REFERENCE_30MER)],
        NS3_FEATURE_MAP,
        _meta(("REF1", "1", "a"), ("G1A", "1", "a"), ("G3A", "3", "a")),
    )
    assert {hit["primary_accession"] for hit in found} == {"REF1", "G1A"}
    by_accession = {row["primary_accession"]: row for row in evidence}
    assert by_accession["G1A"]["scope_tier"] == AM.SCOPE_TIER_SUBTYPE
    assert by_accession["G3A"]["scope_tier"] == AM.SCOPE_TIER_OUT_OF_SCOPE
    assert by_accession["G3A"]["call_status"] == AM.CALL_STATUS_SUPPRESSED_OUT_OF_SCOPE


# --------------------------------------------------------------------------
# 3. Wild-type suppression (rule B)
# --------------------------------------------------------------------------

def test_wild_type_anchor_is_suppressed_and_a_real_change_is_emitted():
    """The residue that IS the wild type for this genotype must not be a call.

    PHDR ships wild-type conjunction anchors (R155R, A156A) so a combination can
    say "155 unchanged AND 168V".  Promoted to standalone signatures they used
    to fire on every sequence carrying the ordinary wild-type residue.
    """
    # codon 6 is AGA=R in the reference; the catalog says wild type for 1a is R.
    anchor_catalog = [_catalog_row(alt_residue="R", display_structure="R6R")]
    found, evidence, _diag = _annotate(
        anchor_catalog,
        [("REF1", REFERENCE_30MER), ("Q1", REFERENCE_30MER)],
        NS3_FEATURE_MAP,
        _meta(("REF1", "1", "a"), ("Q1", "1", "a")),
    )
    assert found == []
    assert {row["residue_status"] for row in evidence} == {AM.RESIDUE_STATUS_ANCHOR}
    assert {row["call_status"] for row in evidence} == {AM.CALL_STATUS_SUPPRESSED_WILD_TYPE}

    # the same position with a wild type of K makes the observed R a change
    change_catalog = [_catalog_row(alt_residue="R", display_structure="K6R")]
    found, evidence, _diag = _annotate(
        change_catalog,
        [("REF1", REFERENCE_30MER), ("Q1", REFERENCE_30MER)],
        NS3_FEATURE_MAP,
        _meta(("REF1", "1", "a"), ("Q1", "1", "a")),
    )
    assert {hit["primary_accession"] for hit in found} == {"REF1", "Q1"}
    assert {row["residue_status"] for row in evidence} == {AM.RESIDUE_STATUS_CHANGE}


def test_unknown_wild_type_emits_the_call_and_says_so():
    """No wild type known is 'we do not know', not 'the residue is unchanged'.

    Suppressing on an absent wild type would hide genuine resistance for every
    genotype PHDR has not scored that position in.  The call is emitted and
    labelled wt_unknown so a reader can tell the two apart.
    """
    catalog = [_catalog_row(alt_residue="R", display_structure="", alignment_name="",
                            relevant_genotypes="")]
    found, evidence, _diag = _annotate(
        catalog,
        [("REF1", REFERENCE_30MER)],
        NS3_FEATURE_MAP,
        _meta(("REF1", "1", "a")),
    )
    assert len(found) == 1
    assert found[0]["residue_status"] == AM.RESIDUE_STATUS_WT_UNKNOWN
    assert found[0]["scope_tier"] == AM.SCOPE_TIER_UNSCOPED
    assert found[0]["wild_type_residues"] == ""
    assert evidence[0]["call_status"] == AM.CALL_STATUS_EMITTED


def test_wild_type_falls_back_from_subtype_to_genotype():
    """Exact subtype first, then any bucket sharing the leading genotype digits."""
    catalog = [
        _catalog_row(alt_residue="R", alignment_name="AL_1a", display_structure="K6R",
                     relevant_genotypes="1a,1b"),
        _catalog_row(alt_residue="R", alignment_name="AL_1b", display_structure="R6R",
                     relevant_genotypes="1a,1b"),
    ]
    prepared, _invalid = AM.prepare_catalog(pd.DataFrame(catalog), {})
    tables = AM.build_wild_type_tables(prepared)

    residues, tier = AM.lookup_wild_type_residues(tables, "1a", "1", "NS3", 6)
    assert residues == {"K"} and tier == AM.SCOPE_TIER_SUBTYPE
    residues, tier = AM.lookup_wild_type_residues(tables, "1b", "1", "NS3", 6)
    assert residues == {"R"} and tier == AM.SCOPE_TIER_SUBTYPE
    # 1c has no bucket of its own, so the genotype-1 union answers
    residues, tier = AM.lookup_wild_type_residues(tables, "1c", "1", "NS3", 6)
    assert residues == {"K", "R"} and tier == AM.SCOPE_TIER_GENOTYPE
    # genotype 3 was never scored here at all
    assert AM.lookup_wild_type_residues(tables, "3a", "3", "NS3", 6) == (None, "")


# --------------------------------------------------------------------------
# 4. The catalogue asset itself
# --------------------------------------------------------------------------

@requires_hcv_assets
def test_catalog_carries_relevant_genotypes_as_a_per_signature_union():
    """relevant_genotypes is the signature's whole scope, not the row's bucket.

    A signature is scored once per alignment bucket, so scoping a call by the
    bucket of the single row it came from would drop the other buckets the same
    finding was curated in - 112 signatures span more than one.
    """
    catalog = _read_catalog().fillna("")
    assert "relevant_genotypes" in catalog.columns
    # The generic columns are appended; order between them is not contractual.
    assert {"relevant_genotypes", "wild_type_residues", "clinical_trials"} <= set(catalog.columns)

    expected = {}
    for signature_id, group in catalog.groupby("signature_id"):
        codes = {
            BCG.alignment_to_genotype_code(value)
            for value in group["alignment_name"]
            if BCG.alignment_to_genotype_code(value)
        }
        expected[signature_id] = set(codes)

    for _, row in catalog.iterrows():
        # Entries are semicolon separated, each optionally carrying ':frequency'.
        got = {e.split(":")[0] for e in str(row["relevant_genotypes"]).split(";") if e}
        assert got == expected[row["signature_id"]], row["signature_id"]

    # expected is now a set of codes per signature, not a joined string.
    multi = [codes for codes in expected.values() if len(codes) > 1]
    empty = [codes for codes in expected.values() if not codes]
    assert len(expected) == 596
    assert len(multi) == 112
    assert len(empty) == 63


@requires_hcv_assets
def test_normalizer_emits_relevant_genotypes_and_phenotype(tmp_path):
    """A regenerated catalogue must not lose the two columns downstream requires.

    output_fields is the only thing write_tsv looks at, so a key
    _build_output_row computes but output_fields omits vanishes with no warning
    - which is how 'phenotype' came to be dropped on regeneration while
    AnnotateMutations lists it as required.
    """
    def dump(name, rows, fields):
        path = tmp_path / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})
        return path

    variation = dump(
        "variation.csv",
        [{"feature_name": "NS3", "name": "phdr_ras:NS3:80K", "type": "aminoAcidSimplePolymorphism",
          "ref_seq_name": "REF_MASTER_NC_004102", "phdr_ras_id": "phdr_ras:NS3:80K"}],
        ["description", "display_name", "feature_name", "name", "ref_end", "ref_seq_name",
         "ref_start", "type", "phdr_ras_id"],
    )
    metatag = dump("variation_metatag.csv", [],
                   ["feature_name", "metatag_name", "metatag_value", "ref_seq_name",
                    "variation_name"])
    alignment = dump(
        "phdr_alignment_ras.csv",
        [{"id": "NS3:80K:AL_1a", "alignment_name": "AL_1a", "display_structure": "K/Q80K",
          "phdr_ras_id": "phdr_ras:NS3:80K"},
         {"id": "NS3:80K:AL_3", "alignment_name": "AL_3", "display_structure": "Q80K",
          "phdr_ras_id": "phdr_ras:NS3:80K"}],
        ["id", "alignment_name", "display_structure", "phdr_ras_id"],
    )
    drug = dump("phdr_alignment_ras_drug.csv", [],
                ["id", "phdr_alignment_ras_id", "phdr_drug_id", "resistance_category",
                 "display_resistance_category", "numeric_resistance_category",
                 "any_in_vitro_evidence", "in_vitro_max_ec50_midpoint", "any_in_vivo_evidence",
                 "in_vivo_baseline", "in_vivo_treatment_emergent"])

    output = tmp_path / "regenerated.tsv"
    NM.HcvMutationCatalogNormalizer(
        variation_path=variation,
        variation_metatag_path=metatag,
        phdr_alignment_ras_path=alignment,
        phdr_alignment_ras_drug_path=drug,
        gene_info_path=GENE_INFO_TSV,
        output_path=output,
    ).normalize()

    regenerated = pd.read_csv(output, sep="\t", dtype=str).fillna("")
    assert "relevant_genotypes" in regenerated.columns
    assert set(regenerated["relevant_genotypes"]) == {"1a,3"}
    missing = [column for column in AM.REQUIRED_MUTATION_CATALOG_COLUMNS
               if column not in regenerated.columns]
    assert missing == []


# --------------------------------------------------------------------------
# 5. Real data: the reference genome must carry no resistance
# --------------------------------------------------------------------------

def _annotate_real_sequences(accessions):
    """Run the annotator over a handful of real sequences from the shipped build."""
    conn = _open_hcv_db()
    try:
        wanted = ["NC_004102"] + [accession for accession in accessions
                                  if accession != "NC_004102"]
        placeholders = ",".join("?" for _ in wanted)
        seq_aln = pd.read_sql_query(
            "SELECT sequence_id, primary_accession, alignment, alignment_name "
            f"FROM sequence_alignment WHERE sequence_id IN ({placeholders})",
            conn, params=wanted,
        )
        meta = pd.read_sql_query(
            "SELECT primary_accession, accession_type, nearest_reference_genotype, "
            "nearest_reference_subtype FROM meta_data", conn,
        )
        # The HCV profile has to be named now: canonicalize_product() only
        # applies a virus's protein-name patterns when it has been told which
        # virus, so that 'core' and 'E1' in a rabies or influenza product
        # description are not rewritten into HCV's mature peptides.
        alias_lookup = AM.load_gene_alias_lookup(conn, profile=AM.virus_profile('HCV'))
        db_gff_maps = AM.load_db_gff_feature_maps(conn, alias_lookup)
    finally:
        conn.close()

    if seq_aln.empty:
        pytest.skip("requested sequences are not present in this build")

    catalog, _invalid = AM.prepare_catalog(_read_catalog(), alias_lookup)
    evidence = []
    found, diagnostics, _maps = AM.annotate_from_reference_coordinates(
        catalog, seq_aln, meta, alias_lookup, db_gff_maps, False, call_evidence=evidence
    )
    return found, evidence, dict(diagnostics)


@requires_hcv_db
@requires_hcv_assets
def test_master_reference_h77_keeps_no_genuine_resistance_call():
    """H77 was isolated in 1977; it cannot carry a direct-acting-antiviral RAS.

    NC_004102 is HCV-1a strain H77, the residue-numbering standard the whole
    catalogue is written against, and the shipped build reports it with 17
    relevant mutations.  Every one of them is either H77's own wild-type residue
    for genotype 1 or a finding scored only in another genotype's bucket, so
    after rules A and B not one survives as a change.
    """
    found, evidence, _diag = _annotate_real_sequences([])
    h77_evidence = [row for row in evidence if row["primary_accession"] == "NC_004102"]
    assert h77_evidence, "expected the reference to be evaluated at all"

    changes = [row for row in found
               if row["primary_accession"] == "NC_004102"
               and row["residue_status"] == AM.RESIDUE_STATUS_CHANGE]
    assert changes == [], f"reference flagged with genuine changes: {changes}"

    mutations = {row["mutation_id"] for row in h77_evidence}
    suppressed = {row["mutation_id"] for row in h77_evidence
                  if row["call_status"] != AM.CALL_STATUS_EMITTED}
    assert len(mutations) == 17
    assert len(suppressed) >= 15


@requires_hcv_db
@requires_hcv_assets
def test_ns3_80_polymorphism_control_h77_is_not_flagged_but_others_are():
    """NS3:80 is a genuine natural polymorphism site; H77 is the negative control.

    HQ537007 is genotype 1a and carries L at NS3 80, which is not a wild type in
    any genotype-1 bucket, so it must stay flagged.  H77 carries Q - the 1a wild
    type - and must not be flagged there whatever else changes.
    """
    found, evidence, _diag = _annotate_real_sequences(["HQ537007"])
    at_80 = [row for row in evidence
             if row["protein_name"] == "NS3" and row["aa_position"] == 80]
    assert at_80, "NS3:80 was never evaluated"

    flagged = {row["primary_accession"] for row in found
               if row["protein_name"] == "NS3" and row["aa_position"] == 80}
    assert "NC_004102" not in flagged
    if "HQ537007" in {row["primary_accession"] for row in at_80}:
        assert "HQ537007" in flagged


@requires_hcv_db
@requires_hcv_assets
@pytest.mark.xfail(
    reason="the 63 alignment-less anchor signatures apply to every genotype, so "
           "one whose position no genotype-1 bucket scores (NS5B 293, NS5B 479) "
           "still reaches the 1977 reference as a wt_unknown call",
    strict=False,
)
def test_master_reference_h77_keeps_no_calls_at_all():
    """Zero, not near-zero: H77 should come out of annotation completely clean.

    NS5B:293L and NS5B:479P are manufactured anchors - conjunction components
    with an empty phdr_ras_id that the normalizer promoted to standalone
    signatures with no drug, no resistance category and no alignment bucket.
    Rule A lets them through because an empty scope means "any genotype", and
    rule B cannot suppress them because the catalogue only ever scores NS5B 293
    and 479 inside AL_2a, so no genotype-1 wild type is known for either.  The
    fix is upstream, in how those anchors are promoted.
    """
    found, _evidence, _diag = _annotate_real_sequences([])
    h77 = [row for row in found if row["primary_accession"] == "NC_004102"]
    assert h77 == [], f"reference still carries calls: {[row['mutation_id'] for row in h77]}"


@requires_hcv_db
@requires_hcv_assets
@pytest.mark.xfail(
    reason="PHDR spells the 1a NS3:80 finding 'K/Q80K', i.e. it lists K as a "
           "typical 1a residue, so rule B classes the simeprevir RAS Q80K as a "
           "wild-type anchor and suppresses it in genotype 1a and 1d",
    strict=False,
)
def test_q80k_is_still_callable_in_genotype_1a():
    """Q80K is the clinically important simeprevir RAS in genotype 1a.

    Rule B is right in general and right here by the letter of the catalogue -
    ``K/Q80K`` says K is one of the typical residues in the AL_1a alignment - but
    the consequence is that NS3:80K is now suppressed for every 1a and 1d
    sequence (and, via ``K80K``, for 6a).  Deciding whether a residue that is
    both typical and a RAS should be called needs curator input, not a code
    change: the catalogue as shipped cannot distinguish the two.
    """
    found, evidence, _diag = _annotate_real_sequences(["KJ439780"])
    at_80 = [row for row in evidence
             if row["primary_accession"] == "KJ439780"
             and row["protein_name"] == "NS3" and row["aa_position"] == 80]
    if not at_80:
        pytest.skip("KJ439780 not present in this build")
    assert any(row["observed_residue"] == "K" for row in at_80)
    assert any(row["primary_accession"] == "KJ439780" and row["mutation_id"] == "NS3:80K"
               for row in found)


# --------------------------------------------------------------------------
# 6. The evidence table records why every call was emitted or suppressed
# --------------------------------------------------------------------------

def test_sequence_mutation_calls_table_records_the_reason_for_every_call(tmp_path):
    """Suppression must leave a trace, or it is indistinguishable from a bug.

    A call that vanished because it was out of genotype scope, one that vanished
    because the residue is the wild type, and one that was never evaluated look
    identical from the compact summary tables alone.
    """
    db_path = tmp_path / "calls.db"
    catalog_path = tmp_path / "catalog.tsv"
    columns = ["protein_name", "segment", "aa_position", "alt_residue", "reference_accession",
               "mutation_id", "mutation_type", "signature_id", "signature_kind",
               "combination_id", "combination_size", "phenotype", "resistance_category",
               "drug", "alignment_name", "display_structure", "relevant_genotypes",
               "wild_type_residues"]
    # wild_type_residues is what the pipeline actually reads; alignment_name and
    # display_structure are carried only as provenance and are never consulted.
    rows = [
        # observed R at codon 6; wild type for 1a is K -> a change
        ["NS3", "1", "6", "R", "REF1", "NS3:6R", "snp", "NS3:6R", "single", "", "", "", "I",
         "drugA", "AL_1a", "K6R", "1a", "1a:K"],
        # observed R at codon 6; wild type for 3a is R -> an anchor
        ["NS3", "1", "6", "R", "REF1", "NS3:6Ranchor", "snp", "NS3:6Ranchor", "single", "", "",
         "", "I", "drugA", "AL_3a", "R6R", "3a", "3a:R"],
    ]
    with catalog_path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE features (accession TEXT, product TEXT, segment TEXT, "
                   "cds_start INTEGER, cds_end INTEGER)")
    cursor.executemany("INSERT INTO features VALUES (?, ?, ?, ?, ?)", [
        ("REF1", "NS3", "1", 1, 30),
        ("Q1A", "polyprotein", "1", 1, 30),
        ("Q3A", "polyprotein", "1", 1, 30),
    ])
    cursor.execute("CREATE TABLE sequence_alignment (sequence_id TEXT, primary_accession TEXT, "
                   "alignment_name TEXT, alignment TEXT)")
    cursor.executemany("INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)", [
        ("REF1", "REF1", "REF1", REFERENCE_30MER),
        ("Q1A", "Q1A", "REF1", REFERENCE_30MER),
        ("Q3A", "Q3A", "REF1", REFERENCE_30MER),
    ])
    cursor.execute("CREATE TABLE meta_data (primary_accession TEXT, accession_type TEXT, "
                   "nearest_reference_genotype TEXT, nearest_reference_subtype TEXT)")
    cursor.executemany("INSERT INTO meta_data VALUES (?, ?, ?, ?)", [
        ("REF1", "master", "1", "a"),
        ("Q1A", "", "1", "a"),
        ("Q3A", "", "3", "a"),
    ])
    conn.commit()
    conn.close()

    argv = sys.argv
    sys.argv = ["AnnotateMutations.py", "--db", str(db_path), "--mutation_catalog",
                str(catalog_path), "--catalog_column_profile", "HCV", "--virus", "HCV"]
    try:
        AM.main()
    finally:
        sys.argv = argv

    conn = sqlite3.connect(str(db_path))
    try:
        calls = pd.read_sql_query("SELECT * FROM sequence_mutation_calls", conn)
        summary = pd.read_sql_query("SELECT * FROM sequence_relevant_mutation_summary", conn)
    finally:
        conn.close()

    for column in ("scope_tier", "residue_status", "call_status", "wild_type_residues",
                   "relevant_genotypes", "sequence_genotype", "sequence_subtype",
                   "observed_residue"):
        assert column in calls.columns

    by_key = {(row.primary_accession, row.mutation_id): row for row in calls.itertuples()}
    assert by_key[("Q1A", "NS3:6R")].call_status == AM.CALL_STATUS_EMITTED
    assert by_key[("Q1A", "NS3:6R")].scope_tier == AM.SCOPE_TIER_SUBTYPE
    assert by_key[("Q1A", "NS3:6R")].residue_status == AM.RESIDUE_STATUS_CHANGE
    assert by_key[("Q1A", "NS3:6R")].wild_type_residues == "K"

    assert by_key[("Q3A", "NS3:6R")].call_status == AM.CALL_STATUS_SUPPRESSED_OUT_OF_SCOPE
    assert by_key[("Q3A", "NS3:6Ranchor")].call_status == AM.CALL_STATUS_SUPPRESSED_WILD_TYPE
    assert by_key[("Q3A", "NS3:6Ranchor")].residue_status == AM.RESIDUE_STATUS_ANCHOR

    present = dict(zip(summary["primary_accession"], summary["relevant_mutations_present"]))
    assert present.get("Q1A") == "NS3:6R"
    assert "Q3A" not in present


def test_empty_mutation_list_still_creates_the_calls_table(tmp_path):
    """An empty run must leave a queryable table, not a missing one."""
    conn = sqlite3.connect(str(tmp_path / "empty.db"))
    try:
        catalog = pd.DataFrame([_catalog_row(drug="drugA", resistance_category="I")])
        AM.write_mutation_tables(conn, catalog, [], "HCV")
        calls = pd.read_sql_query("SELECT * FROM sequence_mutation_calls", conn)
        assert calls.empty
        assert calls.columns.tolist() == AM.MUTATION_CALL_COLUMNS
    finally:
        conn.close()
