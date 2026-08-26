"""Edge cases in the HCV mutation pipeline where data is silently mishandled.

Lens: scripts/AnnotateMutations.py, scripts/VerifyMutations.py,
scripts/NormalizeHcvMutationCatalog.py and the curated PHDR assets under
generic/hcv/Tables/.

The recurring shape hunted here is the one that already bit the influenza
build: information that exists on both sides of a join being dropped because a
per-row/per-column meaning is applied file-wide, and legitimate values that
look like sentinels.  In this pipeline the discarded information is the
*wild-type residue* and the *genotype scope* of every catalogued mutation, and
the sentinel-shaped values are ``del``/``*``/``X`` alt residues.

Real-DB assertions open the database read-only and skip when it is absent.
Every synthetic case writes only into ``tmp_path``.
"""

import csv
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import AnnotateMutations as AM
import NormalizeHcvMutationCatalog as NM


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
HCV_TABLES = REPO_ROOT / "generic" / "hcv" / "Tables"
CATALOG_TSV = HCV_TABLES / "generalized_mutation_catalog_with_extra_info.tsv"
NORMALIZED_TSV = HCV_TABLES / "generalized_mutation_catalog.tsv"
VARIATION_CSV = HCV_TABLES / "variation.csv"
VARIATION_METATAG_CSV = HCV_TABLES / "variation_metatag.csv"
PHDR_ALIGNMENT_RAS_CSV = HCV_TABLES / "phdr_alignment_ras.csv"
PHDR_ALIGNMENT_RAS_DRUG_CSV = HCV_TABLES / "phdr_alignment_ras_drug.csv"
GENE_INFO_TSV = HCV_TABLES / "gene_info.tsv"

HCV_DB = REPO_ROOT / "test_out" / "HCV_OM_test" / "HCV_OM_test.db"

# NC_004102 (HCV-1a strain H77) mature-peptide starts, 1-based inclusive, as
# stored in the `features` table of a real HCV build.
NC_004102_CDS_STARTS = {"NS3": 3420, "NS5A": 6258, "NS5B": 7602}

requires_hcv_db = pytest.mark.skipif(
    not HCV_DB.exists(), reason=f"HCV test database not present at {HCV_DB}"
)
requires_hcv_assets = pytest.mark.skipif(
    not CATALOG_TSV.exists(), reason=f"HCV mutation catalog not present at {CATALOG_TSV}"
)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _read_hcv_db(query):
    """Read a query from the shared HCV build, read-only so an active build is safe."""
    conn = sqlite3.connect(f"file:{HCV_DB}?mode=ro", uri=True)
    try:
        return pd.read_sql_query(query, conn)
    finally:
        conn.close()


def _catalog_row(**overrides):
    """A minimal catalog row carrying every column AnnotateMutations requires."""
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
    }
    row.update(overrides)
    return row


def _annotate(catalog_rows, alignments, feature_map, master="REF1"):
    """Run annotate_from_reference_coordinates over a tiny synthetic alignment.

    `alignments` is a list of (accession, padded_alignment) pairs; the first is
    treated as the master reference.  Returns (mutations_found, diagnostics).
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
    meta = pd.DataFrame([{"primary_accession": master, "accession_type": "master"}])
    prepared, _invalid = AM.prepare_catalog(pd.DataFrame(catalog_rows), {})
    found, diagnostics, _maps = AM.annotate_from_reference_coordinates(
        prepared, seq_aln, meta, {}, {master: feature_map}, False
    )
    return found, dict(diagnostics)


def _write_normalizer_inputs(tmp_path, variation_rows, metatag_rows,
                             alignment_rows=None, drug_rows=None):
    """Materialise the four PHDR source tables the normalizer consumes."""
    def dump(name, rows, fields):
        path = tmp_path / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})
        return path

    variation = dump(
        "variation.csv", variation_rows,
        ["description", "display_name", "feature_name", "name", "ref_end",
         "ref_seq_name", "ref_start", "type", "phdr_ras_id"],
    )
    metatag = dump(
        "variation_metatag.csv", metatag_rows,
        ["feature_name", "metatag_name", "metatag_value", "ref_seq_name", "variation_name"],
    )
    alignment = dump(
        "phdr_alignment_ras.csv", alignment_rows or [],
        ["id", "alignment_name", "display_structure", "phdr_ras_id"],
    )
    drug = dump(
        "phdr_alignment_ras_drug.csv", drug_rows or [],
        ["id", "phdr_alignment_ras_id", "phdr_drug_id", "resistance_category",
         "display_resistance_category", "numeric_resistance_category",
         "any_in_vitro_evidence", "in_vitro_max_ec50_midpoint",
         "any_in_vivo_evidence", "in_vivo_baseline", "in_vivo_treatment_emergent"],
    )
    output = tmp_path / "generalized_mutation_catalog.tsv"
    NM.HcvMutationCatalogNormalizer(
        variation_path=variation,
        variation_metatag_path=metatag,
        phdr_alignment_ras_path=alignment,
        phdr_alignment_ras_drug_path=drug,
        gene_info_path=GENE_INFO_TSV,
        output_path=output,
    ).normalize()
    return pd.read_csv(output, sep="\t", dtype=str).fillna("")


# --------------------------------------------------------------------------
# 1. Coordinate arithmetic - correct behaviour, pinned against curated data
# --------------------------------------------------------------------------

@requires_hcv_assets
def test_codon_arithmetic_matches_curated_phdr_nucleotide_coordinates():
    """Pin the 1-based inclusive codon arithmetic against PHDR's own coordinates.

    variation.csv ships curated ref_start/ref_end nucleotide coordinates for
    every catalogued residue, relative to REF_MASTER_NC_004102.  Independently,
    AnnotateMutations computes the codon as cds_start + 3 * (aa_position - 1)
    over a coordinate map whose keys are 1-based ungapped positions.  If anyone
    ever slides a cds_start to 0-based, or makes cds_end exclusive, or drops the
    (aa_pos - 1), the two disagree and every residue call shifts frame - which
    produces plausible-looking wrong amino acids rather than an error.  This
    test is the tripwire for that.
    """
    variations = pd.read_csv(VARIATION_CSV, dtype=str).fillna("")
    identity_coord_map = {position: position - 1 for position in range(1, 10_000)}

    checked = 0
    for _, row in variations.iterrows():
        feature = row["feature_name"]
        if row["type"] == "conjunction" or feature not in NC_004102_CDS_STARTS:
            continue
        match = re.match(r"^phdr_ras:\w+:(\d+)([A-Z*]|del)$", row["name"])
        if not match:
            continue
        aa_position = int(match.group(1))
        indices = AM.resolve_aligned_codon_indices(
            identity_coord_map, NC_004102_CDS_STARTS[feature], aa_position
        )
        assert indices is not None, row["name"]
        # +1 converts the 0-based alignment column back to a 1-based nucleotide.
        assert indices[0] + 1 == int(row["ref_start"]), row["name"]
        assert indices[2] + 1 == int(row["ref_end"]), row["name"]
        checked += 1

    assert checked > 200, f"expected the full curated set, only checked {checked}"


def test_translate_codon_never_invents_a_residue_from_gaps_or_ambiguity():
    """Gapped and ambiguous codons must degrade to X, never to a real residue.

    translate_codon strips '-' before length-checking, so a codon such as 'A-G'
    becomes the 2-mer 'AG'.  If the length guard were ever relaxed (or the strip
    moved after a table lookup) a deletion would translate to whatever the
    remaining bases spell, and a deleted residue would be reported as a
    substitution.  Locking this down also documents that missing data and
    genuine ambiguity are indistinguishable downstream - both are 'X'.
    """
    assert AM.translate_codon("ATG") == "M"
    assert AM.translate_codon("atg") == "M"
    assert AM.translate_codon("TAA") == "*"
    for undecidable in ("A-G", "---", "AT-", "-TG", "NNN", "RGA", "AAN", "AT"):
        assert AM.translate_codon(undecidable) == "X", undecidable


# --------------------------------------------------------------------------
# 2. Wild-type residues reported as mutations
# --------------------------------------------------------------------------

@pytest.mark.xfail(
    reason="AnnotateMutations matches alt_residue alone and never compares the "
           "reference/wild-type residue, so a catalog row whose alt residue IS "
           "the wild type flags every sequence including the reference itself",
    strict=False,
)
def test_wild_type_residue_is_not_reported_as_a_mutation():
    """A catalog row whose alt residue equals the reference residue must not hit.

    Real-world trigger: PHDR's variation.csv contains wild-type conjunction
    anchors such as NS3 R155R and NS3 A156A (they exist so a combination can
    say 'R155 unchanged AND D168V').  Normalisation promotes each of them to a
    standalone catalogue entry, and annotation then reports them on every
    sequence that carries the ordinary wild-type residue - silently converting
    'nothing to see here' into a drug-resistance call.
    """
    reference = "ATGAAACCCGGGTTTAGATAGCATGCATGC"   # codon 6 = AGA = R
    query = "ATGAAACCCGGGTTTAGATAGCATGCATGC"       # identical to the reference
    feature_map = {"NS3": {"product": "NS3", "cds_start": 1, "cds_end": 30,
                           "feature_type": "mat_peptide"}}
    found, _diag = _annotate(
        [_catalog_row(alt_residue="R", mutation_id="NS3:6R")],
        [("REF1", reference), ("Q1", query)],
        feature_map,
    )
    flagged = {hit["primary_accession"] for hit in found}
    assert "REF1" not in flagged, "the reference sequence carries no mutations by definition"
    assert flagged == set()


@requires_hcv_db
@pytest.mark.xfail(
    reason="the master reference NC_004102 is itself annotated with wild-type "
           "residues as 'relevant mutations present' and 'completed signatures'",
    strict=False,
)
def test_master_reference_is_not_annotated_with_its_own_residues():
    """The reference genome must have zero mutations relative to itself.

    NC_004102 is HCV-1a strain H77, the residue-numbering standard the whole
    catalogue is written against.  In the shipped build it is reported with 17
    relevant mutations and 22 completed signatures - among them the singles
    NS5A:30Q and NS5A:31L, which are precisely H77's wild-type residues.  Any
    clinician-facing view of `sequence_relevant_mutation_summary` therefore
    reads wild type as resistance.
    """
    summary = _read_hcv_db(
        "SELECT * FROM sequence_relevant_mutation_summary "
        "WHERE primary_accession = 'NC_004102'"
    )
    if summary.empty:
        pytest.skip("master reference NC_004102 not present in this build")
    present = [
        value
        for value in str(summary.iloc[0]["relevant_mutations_present"] or "").split(";")
        if value
    ]
    completed = _read_hcv_db(
        "SELECT * FROM completed_signatures_only WHERE primary_accession = 'NC_004102'"
    )
    assert present == [], f"reference flagged with its own residues: {present}"
    assert completed.empty


@requires_hcv_db
@requires_hcv_assets
@pytest.mark.xfail(
    reason="~50% of relevant_mutations_present entries are drug-less 'single' "
           "signatures created from non-RAS conjunction anchors",
    strict=False,
)
def test_relevant_mutations_summary_only_contains_drug_annotated_signatures():
    """Every entry surfaced as a 'relevant mutation' should carry a drug annotation.

    Real-world trigger: 55 of the 232 non-conjunction rows in variation.csv have
    an empty `phdr_ras_id`, i.e. PHDR is saying 'this residue is a combination
    component, not a resistance-associated substitution'.  The normalizer's
    `row.get("phdr_ras_id") or row.get("name")` treats empty as missing and
    substitutes the variation name, so they become standalone signatures with no
    drug, no resistance category and no genotype - and they then dominate the
    per-sequence summary.
    """
    catalog = pd.read_csv(CATALOG_TSV, sep="\t", dtype=str).fillna("")
    singles = catalog[catalog["signature_kind"] == "single"]
    drugless = set(singles.loc[singles["drug"] == "", "mutation_id"]) - set(
        singles.loc[singles["drug"] != "", "mutation_id"]
    )

    summary = _read_hcv_db("SELECT * FROM sequence_relevant_mutation_summary")
    if summary.empty:
        pytest.skip("no mutation annotations in this build")

    total = 0
    drugless_calls = 0
    for value in summary["relevant_mutations_present"]:
        for mutation_id in str(value or "").split(";"):
            if not mutation_id:
                continue
            total += 1
            if mutation_id in drugless:
                drugless_calls += 1

    assert total > 0
    assert drugless_calls == 0, (
        f"{drugless_calls}/{total} summary entries are drug-less non-RAS anchors"
    )


# --------------------------------------------------------------------------
# 3. Genotype / reference resolution
# --------------------------------------------------------------------------

@requires_hcv_db
@requires_hcv_assets
@pytest.mark.xfail(
    reason="the catalog's alignment_name (AL_1a, AL_3a, ...) is not in the HCV "
           "column profile, so the DB mutation_catalog keeps no genotype scope",
    strict=False,
)
def test_db_mutation_catalog_retains_the_catalog_genotype_scope():
    """The genotype a mutation was curated for must survive into the database.

    Every catalogue row carries `alignment_name` - AL_1a, AL_1b, AL_3a and so on
    - which is the genotype/subtype bucket the PHDR curators scored the mutation
    in, and 110 of the 232 mutations differ in scope between buckets.  The HCV
    column profile in AnnotateMutations does not list `alignment_name` (nor
    `display_structure`, which holds the wild-type letter), so
    build_catalog_reference_table drops both and then de-duplicates the
    survivors.  Downstream nobody can tell that NS5A:30R is a 1a finding, and
    nobody can re-scope the calls afterwards.
    """
    source = pd.read_csv(CATALOG_TSV, sep="\t", dtype=str)
    assert source["alignment_name"].notna().any()
    columns = set(_read_hcv_db("SELECT * FROM mutation_catalog LIMIT 1").columns)
    assert "alignment_name" in columns
    assert "display_structure" in columns


@requires_hcv_db
@requires_hcv_assets
@pytest.mark.xfail(
    reason="catalog rows are applied to every sequence regardless of genotype; "
           "~90% of calls in the shipped HCV DB are out of the curated scope",
    strict=False,
)
def test_mutation_calls_respect_the_genotype_they_were_curated_for():
    """A 1a-only resistance mutation must not be called on a genotype 3 sequence.

    Both halves of the join exist: meta_data carries
    `nearest_reference_genotype` / `nearest_reference_subtype` per sequence, and
    the catalogue carries `alignment_name` per row.  AnnotateMutations reads
    neither, so a mutation curated only in AL_1a is scored against genotype 2,
    3, 4, 5, 6 and 8 sequences alike.  Because many of those residues are simply
    the wild type of the other genotype (e.g. NS5A 30R is 100% conserved in
    4r), this manufactures resistance calls that look entirely plausible.
    """
    catalog = pd.read_csv(CATALOG_TSV, sep="\t", dtype=str)
    scope = catalog.groupby("mutation_id")["alignment_name"].apply(
        lambda values: {value for value in values.dropna()}
    )
    meta = _read_hcv_db(
        "SELECT primary_accession, nearest_reference_genotype, "
        "nearest_reference_subtype FROM meta_data"
    )
    summary = _read_hcv_db("SELECT * FROM sequence_relevant_mutation_summary")
    if summary.empty:
        pytest.skip("no mutation annotations in this build")
    merged = summary.merge(meta, on="primary_accession", how="left")

    total = 0
    out_of_scope = 0
    for _, row in merged.iterrows():
        genotype = str(row["nearest_reference_genotype"] or "").strip()
        subtype = str(row["nearest_reference_subtype"] or "").strip()
        if not genotype:
            continue
        labels = {f"AL_{genotype}", f"AL_{genotype}{subtype}"}
        for mutation_id in str(row["relevant_mutations_present"] or "").split(";"):
            if not mutation_id:
                continue
            curated = scope.get(mutation_id, set())
            if not curated:
                continue
            total += 1
            if not (curated & labels):
                out_of_scope += 1

    assert total > 0
    assert out_of_scope == 0, (
        f"{out_of_scope}/{total} calls applied outside their curated genotype"
    )


# --------------------------------------------------------------------------
# 4. Alt-residue vocabulary: del, stop codons, ambiguity, whitespace, case
# --------------------------------------------------------------------------

@pytest.mark.xfail(
    reason="alt_residue 'del' can never equal a translate_codon result, so "
           "catalogued deletions are silently never annotated",
    strict=False,
)
def test_catalogued_deletions_are_annotated_when_the_codon_is_deleted():
    """NS5A deletion RASs must be detectable, or at least loudly unsupported.

    Real-world trigger: the shipped catalogue contains NS5A:29del, NS5A:30del
    and NS5A:32del (mutation_type aminoAcidDeletion) - NS5A N-terminal deletions
    are a documented DAA-resistance finding.  Annotation compares
    translate_codon(codon) against alt_residue, and translate_codon returns 'X'
    for a fully gapped codon, never the literal 'del'.  The rows therefore
    contribute zero hits with no diagnostic counter incremented at all: a
    deletion carrier looks clean.
    """
    reference = "ATGAAACCCGGGTTTAGATAGCATGCATGC"
    deleted = "ATGAAACCCGGGTTT---TAGCATGCATGC"   # codon 6 removed
    feature_map = {"NS3": {"product": "NS3", "cds_start": 1, "cds_end": 30,
                           "feature_type": "mat_peptide"}}
    found, diagnostics = _annotate(
        [_catalog_row(alt_residue="del", mutation_id="NS3:6del")],
        [("REF1", reference), ("DEL1", deleted)],
        feature_map,
    )
    assert {hit["primary_accession"] for hit in found} == {"DEL1"}, diagnostics


@pytest.mark.xfail(
    reason="CODON_TABLE encodes stop as '_' while the catalog token grammar "
           "accepts '*', so a nonsense mutation can never match",
    strict=False,
)
def test_stop_codon_catalog_entry_matches_a_real_stop_codon():
    """A '*' stop entry must match a TAG/TAA/TGA codon.

    Real-world trigger: NormalizeHcvMutationCatalog's MUTATION_TOKEN_RE is
    ``^(\\d+)([A-Z*]|del)$`` - it explicitly admits '*' - so the day PHDR adds a
    nonsense variant the normalizer will happily emit it.  AnnotateMutations'
    CODON_TABLE spells stop as '_', so the row silently matches nothing forever.
    The two spellings of the same concept never meet.
    """
    reference = "ATGAAACCCGGGTTTAGATAGCATGCATGC"   # codon 7 = TAG = stop
    feature_map = {"NS3": {"product": "NS3", "cds_start": 1, "cds_end": 30,
                           "feature_type": "mat_peptide"}}
    found, _diag = _annotate(
        [_catalog_row(aa_position="7", alt_residue="*", mutation_id="NS3:7*")],
        [("REF1", reference)],
        feature_map,
    )
    assert {hit["primary_accession"] for hit in found} == {"REF1"}


@pytest.mark.xfail(
    reason="alt_residue is compared raw: stray whitespace or lower case makes a "
           "catalog row unmatchable and no diagnostic is recorded",
    strict=False,
)
@pytest.mark.parametrize("alt_residue", ["R ", " R", "r"])
def test_whitespace_or_case_in_alt_residue_does_not_silently_disable_a_row(alt_residue):
    """A residue with a stray space or lower case must still match, or be reported.

    Real-world trigger: the catalogue is a hand-maintained TSV regenerated from
    curator spreadsheets; a trailing space in one cell is invisible in every
    viewer.  prepare_catalog canonicalises protein_name and coerces aa_position
    but leaves alt_residue untouched, and the annotation loop's
    ``alt_lookup.get(aa, [])`` simply finds nothing.  The run still exits 0, the
    mapping summary still reports mutation_hits for the other rows, and the
    disabled mutation is indistinguishable from one that is genuinely absent.
    """
    reference = "ATGAAACCCGGGTTTAGATAGCATGCATGC"   # codon 6 = AGA = R
    feature_map = {"NS3": {"product": "NS3", "cds_start": 1, "cds_end": 30,
                           "feature_type": "mat_peptide"}}
    found, diagnostics = _annotate(
        [_catalog_row(alt_residue=alt_residue)],
        [("REF1", reference)],
        feature_map,
    )
    assert found or diagnostics.get("unmatchable_catalog_rows"), (
        f"alt_residue {alt_residue!r} silently matched nothing; diagnostics={diagnostics}"
    )


@pytest.mark.xfail(
    reason="an 'X' alt residue matches gapped and ambiguous codons, turning "
           "missing data into a positive mutation call",
    strict=False,
)
def test_ambiguous_or_missing_codon_does_not_satisfy_an_x_catalog_row():
    """Absence of evidence must not be reported as evidence of a mutation.

    translate_codon collapses three very different situations onto 'X': a
    deleted codon ('---'), an unsequenced codon ('NNN') and a genuine IUPAC
    ambiguity ('RGA').  Any catalogue row with alt_residue 'X' - a shape the
    normalizer's [A-Z*] token grammar admits - therefore fires on every
    truncated or low-coverage sequence.  GenBank HCV records are routinely
    partial, so this would flag the worst-covered sequences the hardest.
    """
    reference = "ATGAAACCCGGGTTTAGATAGCATGCATGC"
    gapped = "ATGAAACCCGGGTTT---TAGCATGCATGC"
    ambiguous = "ATGAAACCCGGGTTTNNNTAGCATGCATGC"
    feature_map = {"NS3": {"product": "NS3", "cds_start": 1, "cds_end": 30,
                           "feature_type": "mat_peptide"}}
    found, _diag = _annotate(
        [_catalog_row(alt_residue="X", mutation_id="NS3:6X")],
        [("REF1", reference), ("GAP", gapped), ("AMB", ambiguous)],
        feature_map,
    )
    assert {hit["primary_accession"] for hit in found} == set()


# --------------------------------------------------------------------------
# 5. Catalog parsing: positions outside the feature, non-integer positions
# --------------------------------------------------------------------------

@pytest.mark.xfail(
    reason="resolve_aligned_codon_indices bounds the position against the "
           "alignment but never against gene_entry['cds_end'], so an "
           "out-of-range position reads a codon from the next gene",
    strict=False,
)
def test_position_past_the_feature_end_does_not_read_the_downstream_gene():
    """A position beyond a protein's length must not silently borrow the next one.

    Real-world trigger: mature-peptide annotations in the `features` table are
    per-reference and can be truncated (partial GenBank records, or
    choose_feature_entry preferring the *shortest* span on a priority tie).  A
    catalogued NS3 position past the truncated NS3 end then resolves into the
    NS4A/NS5A reading frame and the resulting residue is reported as an NS3
    mutation - correct-looking output attributed to the wrong protein.
    """
    reference = "ATGAAACCCGGGTTTAGATAGCATGCATGC"
    feature_map = {
        "NS3": {"product": "NS3", "cds_start": 1, "cds_end": 9,
                "feature_type": "mat_peptide"},          # only 3 codons long
        "NS5A": {"product": "NS5A", "cds_start": 10, "cds_end": 30,
                 "feature_type": "mat_peptide"},
    }
    found, _diag = _annotate(
        [_catalog_row(aa_position="6", alt_residue="R")],
        [("REF1", reference)],
        feature_map,
    )
    assert found == [], "NS3 has 3 codons; position 6 must not resolve into NS5A"


@pytest.mark.xfail(
    reason="prepare_catalog uses int(float(value)), silently accepting '6.7', "
           "'1e1', '0' and '-3' as positions",
    strict=False,
)
@pytest.mark.parametrize("aa_position", ["6.7", "1e1", "0", "-3"])
def test_non_integer_or_non_positive_positions_are_flagged_as_invalid(aa_position):
    """Positions that are not whole positive integers must be counted as invalid.

    Real-world trigger: the catalogue is regenerated by a chain that passes
    positions through pandas and spreadsheets, where '107' can come back as
    '107.0' and a large accession-like token as '1e5'.  int(float(value))
    swallows all of them, so '6.7' becomes 6 and '1e1' becomes 10 - a mutation
    silently relocated to a different residue.  Zero and negative positions do
    get skipped, but they are tallied under
    'reference_coordinate_missing_in_alignment', which hides a catalog defect
    behind an alignment-coverage counter.
    """
    prepared, invalid = AM.prepare_catalog(
        pd.DataFrame([_catalog_row(aa_position=aa_position)]), {}
    )
    assert invalid == 1, (
        f"aa_position {aa_position!r} coerced to "
        f"{prepared['_aa_position_int'].tolist()} instead of being rejected"
    )


# --------------------------------------------------------------------------
# 6. Normalizer: catalogue generation
# --------------------------------------------------------------------------

@requires_hcv_assets
def test_normalizer_regenerates_the_shipped_catalog_byte_for_byte(tmp_path):
    """Regenerating the catalogue from the PHDR tables must be reproducible.

    The catalogue is a checked-in build product; if regeneration is not stable
    (dict iteration order, unsorted groupby) then a rebuild silently reshuffles
    rows and every downstream diff becomes unreadable, hiding real curation
    changes.  This pins the whole normalizer end to end against the committed
    output for the columns it declares.
    """
    output = tmp_path / "regenerated.tsv"
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / "NormalizeHcvMutationCatalog.py"),
         "--variation", str(VARIATION_CSV),
         "--variation_metatag", str(VARIATION_METATAG_CSV),
         "--phdr_alignment_ras", str(PHDR_ALIGNMENT_RAS_CSV),
         "--phdr_alignment_ras_drug", str(PHDR_ALIGNMENT_RAS_DRUG_CSV),
         "--gene_info", str(GENE_INFO_TSV),
         "--output_path", str(output)],
        capture_output=True, text=True, check=True,
    )
    assert output.exists(), result.stderr

    regenerated = pd.read_csv(output, sep="\t", dtype=str).fillna("")
    shipped = pd.read_csv(NORMALIZED_TSV, sep="\t", dtype=str).fillna("")
    assert len(regenerated) == len(shipped)
    # The regenerated file may legitimately declare *more* columns than the
    # committed one (a newly added column has not been checked in yet); it must
    # never drop or reorder the ones already there, nor change a single value.
    missing = [column for column in shipped.columns if column not in regenerated.columns]
    assert missing == [], f"regeneration dropped columns: {missing}"
    shared = list(shipped.columns)
    assert regenerated[shared].equals(shipped[shared])


@requires_hcv_assets
@pytest.mark.xfail(
    reason="output_fields omits 'phenotype' although _build_output_row produces "
           "it, so regeneration drops a column AnnotateMutations requires",
    strict=False,
)
def test_normalizer_writes_every_column_annotate_mutations_requires(tmp_path):
    """A regenerated catalogue must still satisfy AnnotateMutations' schema check.

    Real-world trigger: someone re-runs NormalizeHcvMutationCatalog.py after a
    PHDR data refresh.  `_build_output_row` computes a "phenotype" key, but
    `output_fields` does not list it and write_tsv rebuilds each row from
    fieldnames only - so the column vanishes with no warning.  The committed
    file has 29 columns, a fresh one has 28, and the next pipeline run dies with
    "Catalog missing required column(s): phenotype" long after the catalogue was
    regenerated.
    """
    output = tmp_path / "regenerated.tsv"
    subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / "NormalizeHcvMutationCatalog.py"),
         "--variation", str(VARIATION_CSV),
         "--variation_metatag", str(VARIATION_METATAG_CSV),
         "--phdr_alignment_ras", str(PHDR_ALIGNMENT_RAS_CSV),
         "--phdr_alignment_ras_drug", str(PHDR_ALIGNMENT_RAS_DRUG_CSV),
         "--gene_info", str(GENE_INFO_TSV),
         "--output_path", str(output)],
        capture_output=True, text=True, check=True,
    )
    regenerated = pd.read_csv(output, sep="\t", dtype=str)
    missing = [
        column
        for column in AM.REQUIRED_MUTATION_CATALOG_COLUMNS
        if column not in regenerated.columns
    ]
    assert missing == [], f"regenerated catalog is missing {missing}"


@pytest.mark.xfail(
    reason="_signature_id_from_variation falls back to `name` when phdr_ras_id "
           "is blank, promoting non-RAS conjunction anchors to standalone "
           "'single' signatures",
    strict=False,
)
def test_variation_without_a_phdr_ras_id_is_not_promoted_to_a_single_signature(tmp_path):
    """An empty phdr_ras_id means 'not a RAS' and must not become a signature.

    Real-world trigger: 55 of the 232 non-conjunction rows in the shipped
    variation.csv have a blank phdr_ras_id - NS3:155R, NS3:156A, NS3:168D and
    friends.  They exist purely so a conjunction can express 'R155 unchanged AND
    D168V'.  ``row.get("phdr_ras_id") or row.get("name", "")`` reads blank as
    absent and substitutes the variation name, so each becomes its own
    ``signature_kind="single"`` row with empty drug, empty resistance_category
    and empty alignment_name - and then gets called on nearly every sequence.
    """
    rows = _write_normalizer_inputs(
        tmp_path,
        variation_rows=[
            # conjunction anchor only: PHDR left phdr_ras_id blank on purpose
            {"feature_name": "NS3", "name": "phdr_ras:NS3:155R",
             "ref_seq_name": "REF_MASTER_NC_004102", "ref_start": "3882",
             "ref_end": "3884", "type": "aminoAcidSimplePolymorphism",
             "phdr_ras_id": ""},
            {"feature_name": "NS3", "name": "phdr_ras:NS3:168V",
             "ref_seq_name": "REF_MASTER_NC_004102", "ref_start": "3921",
             "ref_end": "3923", "type": "aminoAcidSimplePolymorphism",
             "phdr_ras_id": "NS3:168V"},
            {"feature_name": "NS3", "name": "phdr_ras:NS3:155R+168V",
             "ref_seq_name": "REF_MASTER_NC_004102", "ref_start": "3882",
             "ref_end": "3923", "type": "conjunction",
             "phdr_ras_id": "NS3:155R+168V"},
        ],
        metatag_rows=[
            {"feature_name": "NS3", "metatag_name": "CONJUNCT_NAME_1",
             "metatag_value": "phdr_ras:NS3:155R",
             "variation_name": "phdr_ras:NS3:155R+168V"},
            {"feature_name": "NS3", "metatag_name": "CONJUNCT_NAME_2",
             "metatag_value": "phdr_ras:NS3:168V",
             "variation_name": "phdr_ras:NS3:155R+168V"},
        ],
    )
    singles = rows[rows["signature_kind"] == "single"]
    assert set(singles["mutation_id"]) == {"NS3:168V"}, (
        "NS3:155R has no phdr_ras_id and must stay a combination component only"
    )


@pytest.mark.xfail(
    reason="_build_combination_memberships consumes every metatag on a "
           "conjunction, not just CONJUNCT_NAME_*, so an unrelated metatag "
           "becomes a phantom required component",
    strict=False,
)
def test_only_conjunct_name_metatags_become_combination_components(tmp_path):
    """Combination membership must be driven by CONJUNCT_NAME_* metatags alone.

    Real-world trigger: variation_metatag.csv already carries other metatag
    kinds (SIMPLE_AA_PATTERN, MIN_COMBINED_TRIPLET_FRACTION); today they only
    hang off single variations, so nothing breaks.  The membership builder
    filters on the *parent* being a conjunction, never on the metatag name, and
    assigns component_order 999 to anything it cannot index.  The first
    conjunction to acquire a non-conjunct metatag whose value happens to parse
    as a mutation token silently gains a required component - inflating
    combination_size and making the signature impossible to complete.
    """
    rows = _write_normalizer_inputs(
        tmp_path,
        variation_rows=[
            {"feature_name": "NS3", "name": "phdr_ras:NS3:155R",
             "ref_seq_name": "REF_MASTER_NC_004102", "ref_start": "3882",
             "ref_end": "3884", "type": "aminoAcidSimplePolymorphism",
             "phdr_ras_id": "NS3:155R"},
            {"feature_name": "NS3", "name": "phdr_ras:NS3:168V",
             "ref_seq_name": "REF_MASTER_NC_004102", "ref_start": "3921",
             "ref_end": "3923", "type": "aminoAcidSimplePolymorphism",
             "phdr_ras_id": "NS3:168V"},
            {"feature_name": "NS3", "name": "phdr_ras:NS3:155R+168V",
             "ref_seq_name": "REF_MASTER_NC_004102", "ref_start": "3882",
             "ref_end": "3923", "type": "conjunction",
             "phdr_ras_id": "NS3:155R+168V"},
        ],
        metatag_rows=[
            {"feature_name": "NS3", "metatag_name": "CONJUNCT_NAME_1",
             "metatag_value": "phdr_ras:NS3:155R",
             "variation_name": "phdr_ras:NS3:155R+168V"},
            {"feature_name": "NS3", "metatag_name": "CONJUNCT_NAME_2",
             "metatag_value": "phdr_ras:NS3:168V",
             "variation_name": "phdr_ras:NS3:155R+168V"},
            # not a component - a display/reporting metatag that happens to
            # carry a mutation token as its value
            {"feature_name": "NS3", "metatag_name": "DISPLAY_ANCHOR",
             "metatag_value": "phdr_ras:NS3:36A",
             "variation_name": "phdr_ras:NS3:155R+168V"},
        ],
    )
    combination = rows[rows["combination_id"] == "NS3:155R+168V"]
    assert set(combination["combination_size"]) == {"2"}
    assert set(combination["mutation_id"]) == {"NS3:155R", "NS3:168V"}


@pytest.mark.xfail(
    reason="_build_output_row coerces a blank segment to the literal '1', so an "
           "unknown segment is indistinguishable from segment 1",
    strict=False,
)
def test_blank_segment_is_not_defaulted_to_segment_one():
    """An unresolvable segment must stay blank rather than become segment '1'.

    Real-world trigger: this is the same shape as the influenza segment bug -
    a per-column default applied file-wide.  ProteinNameMapper hands back '' for
    anything it cannot place (whole_genome, or a gene_info row whose parent is
    NULL), and ``mutation_row.get("segment", "") or "1"`` turns that into the
    literal segment '1'.  For HCV every value is genuinely 1 so nothing shows;
    the moment this normalizer is pointed at a segmented virus, every
    unresolved mutation is filed under segment 1 (PB2 for influenza) and the
    error is invisible because '1' is a legitimate value.
    """
    mutation_row = {
        "canonical_protein_name": "whole_genome",
        "segment": "",
        "aa_position": "5",
        "alt_residue": "A",
        "reference_accession": "REF_MASTER_NC_004102",
        "mutation_id": "whole_genome:5A",
        "mutation_type": "aminoAcidSimplePolymorphism",
    }
    row = NM.HcvMutationCatalogNormalizer._build_output_row(
        mutation_row=mutation_row,
        signature_id="whole_genome:5A",
        signature_kind="single",
        combination_id="",
        combination_size="",
        component_order="",
        source_variation_name="",
        source_phdr_ras_id="",
        alignment_context={},
        drug_row={},
    )
    assert row["segment"] == "", f"blank segment silently became {row['segment']!r}"
