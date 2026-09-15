"""Coordinate-space integrity when a master CDS is projected onto an alignment.

Lens: the three coordinate spaces described at the top of
``scripts/CalcAlignmentCord.py``. Space 1 is the master's own ungapped genome
(what a GFF holds), space 2 is alignment columns, space 3 is one row's own
ungapped numbering. The failure this suite exists for is never a traceback -
it is a coordinate that is merely *wrong*, which produces a CDS that still
translates, still looks like a protein, and is attributed to the right gene.

The specific defect that motivated it: gap counts measured in space 2 were
subtracted from space-1 coordinates, which is a no-op while the master's own
alignment row is ungapped (RABV, HCV) and shifts every ``og`` coordinate
upstream by the master's gap count as soon as a guide alignment introduces
columns the master lacks (influenza segments 2, 4, 6, 7 and 8). Those gap
counts are multiples of three, so the reading frame survived and the only
symptom was a translation starting on a UTR codon.

Everything here is synthetic. Where a test encodes a limitation rather than a
guarantee it is an xfail carrying the reason, so that fixing the limitation
fails the suite loudly instead of going unnoticed.
"""

import random

import pytest
from Bio.Seq import Seq

from CalcAlignmentCord import CalculateAlignmentCoordinates


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def processor():
    """A bare instance; none of the coordinate maths touches instance state."""
    return CalculateAlignmentCoordinates(None, None, None, None, None, None, None)


def coord_map(master_alignment):
    """space 1 -> space 2, exactly as find_gaps_in_fasta builds it."""
    mapping = {}
    residue = 0
    for aln_pos, base in enumerate(master_alignment, start=1):
        if base != "-":
            residue += 1
            mapping[residue] = aln_pos
    return mapping


def project(master_alignment, query_alignment, cds_list, span=None):
    """Project every CDS onto one row, returning the raw adjusted entries."""
    proc = processor()
    span_start, span_end = span if span else (None, None)
    return proc.recalculate_cds_coordinates_with_span(
        "Q",
        proc.get_gap_ranges(query_alignment),
        cds_list,
        start_offset=1,
        genome_cord_start=span_start,
        genome_cord_end=span_end,
        master_coord_to_aln_pos=coord_map(master_alignment),
    )


def covered_span(master_alignment, query_alignment):
    """The space-1 span a row covers, as CalculateGenomeCoordinates reports it."""
    first = last = None
    residue = 0
    for aln_pos, base in enumerate(master_alignment, start=1):
        if base == "-":
            continue
        residue += 1
        if query_alignment[aln_pos - 1] != "-":
            if first is None:
                first = residue
            last = residue
    return first, last


def degap(text):
    return text.replace("-", "")


def translate(nt):
    nt = degap(nt)
    nt = nt[: len(nt) - len(nt) % 3]
    return str(Seq(nt).translate()) if nt else ""


def has_internal_stop(protein):
    body = protein[:-1] if protein.endswith("*") else protein
    return "*" in body


# A 12-codon ORF: ATG, ten sense codons, TAA.
ORF = "ATG" + "AAACCCGGGTTTAAGCCGGTCACGTGCAGT" + "TAA"
assert len(ORF) % 3 == 0


def build_master(insertion_columns, at_master_coord, genome=ORF):
    """Master row carrying ``insertion_columns`` gaps before ``at_master_coord``.

    Returns ``(master_alignment, cds_list)`` for a CDS spanning the whole genome.
    """
    cut = at_master_coord - 1
    master = genome[:cut] + "-" * insertion_columns + genome[cut:]
    return master, [{"start": "1", "end": str(len(genome)), "product": "ORF"}]


# ---------------------------------------------------------------------------
# 1. the master's own row: space 3 must come back out as space 1
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("gaps", [0, 3, 6, 9, 102, 129])
def test_master_own_og_coords_equal_its_genome_coordinates(gaps):
    """The invariant AnnotateMutations resolves every catalogue position through.

    A master's ungapped numbering *is* the master genome, so its own
    ``og`` coordinates must be its plain GFF coordinates no matter how many
    gaps a guide alignment puts in its row. Counting the row's gaps against a
    genome coordinate returned ``1 - gaps`` here, and every catalogue position
    downstream then resolved to a codon that many bases early.
    """
    master, cds_list = build_master(gaps, at_master_coord=1)
    entry = project(master, master, cds_list)[0]

    assert (entry["og_start"], entry["og_end"]) == (1, len(ORF))
    assert degap(master)[entry["og_start"] - 1:entry["og_end"]] == ORF


@pytest.mark.parametrize("gaps", [3, 6, 9, 129])
@pytest.mark.parametrize("at", [1, 4, 16, 31])
def test_master_og_coords_survive_an_insertion_anywhere_in_the_genome(gaps, at):
    """Upstream, inside and downstream of the CDS all have to behave."""
    master, cds_list = build_master(gaps, at_master_coord=at)
    entry = project(master, master, cds_list)[0]

    assert (entry["og_start"], entry["og_end"]) == (1, len(ORF))
    assert not has_internal_stop(translate(ORF))


# ---------------------------------------------------------------------------
# 2. the reported field failure, both ways round
# ---------------------------------------------------------------------------


def test_query_carrying_an_insertion_keeps_its_start_codon():
    """The influenza M1/NS1 signature: no ATG at the recorded start.

    The master is gapped where this query has bases, so a gap count taken
    against the master coordinate is zero and the CDS window sits three bases
    early - beginning on the UTR codon that precedes it.
    """
    master, cds_list = build_master(3, at_master_coord=1)
    query = "TGA" + ORF  # carries the insertion, and it happens to be a stop

    entry = project(master, query, cds_list)[0]
    raw = degap(query)
    cds = raw[entry["og_start"] - 1:entry["og_end"]]

    assert cds.startswith("ATG")
    assert cds == ORF
    assert not has_internal_stop(translate(cds))


def test_query_lacking_the_insertion_keeps_coordinates_in_range():
    """The same master, for a row gapped exactly where the master is.

    Counting from the master coordinate put ``og_start`` at 0 here, because
    space-1 position 1 falls inside this row's own leading gap.
    """
    master, cds_list = build_master(3, at_master_coord=1)
    query = "---" + ORF

    entry = project(master, query, cds_list)[0]

    assert entry["og_start"] >= 1
    assert (entry["og_start"], entry["og_end"]) == (1, len(ORF))
    assert degap(query)[entry["og_start"] - 1:entry["og_end"]] == ORF


def test_gap_count_that_is_not_a_multiple_of_three_still_resolves():
    """Nothing enforces codon-sized insertions in a guide alignment.

    Every count in generic/influenza/ref_set_aligned happens to be a multiple
    of three, which is why the defect stayed in frame and therefore quiet. A
    curated insertion of any other length would have produced an obvious
    frameshift instead - the projection must be right for both.
    """
    for gaps in (1, 2, 4, 5, 7):
        master, cds_list = build_master(gaps, at_master_coord=1)
        query = "N" * gaps + ORF
        entry = project(master, query, cds_list)[0]
        cds = degap(query)[entry["og_start"] - 1:entry["og_end"]]
        assert cds == ORF, f"gaps={gaps}"


# ---------------------------------------------------------------------------
# 3. the general round-trip property
# ---------------------------------------------------------------------------


def assert_round_trip(master, query, cds_list, span=None):
    """og coordinates must cut from raw what the aln coordinates cut from the row.

    Only asserted where both boundary columns carry a residue in this row;
    a boundary inside a gap has no residue to name and is covered separately.
    """
    for entry in project(master, query, cds_list, span=span):
        start_col, end_col = entry["start"], entry["end"]
        if query[start_col - 1] == "-" or query[end_col - 1] == "-":
            continue
        from_alignment = degap(query[start_col - 1:end_col])
        from_raw = degap(query)[entry["og_start"] - 1:entry["og_end"]]
        assert from_raw == from_alignment, entry


@pytest.mark.parametrize("gaps", [0, 1, 3, 12])
@pytest.mark.parametrize(
    "query_pattern",
    [
        "full",          # no gaps at all
        "internal",      # a deletion inside the CDS
        "leading",       # 5' padding
        "trailing",      # 3' padding
        "both",          # padded at both ends
    ],
)
def test_og_and_alignment_coordinates_agree(gaps, query_pattern):
    master, cds_list = build_master(gaps, at_master_coord=7)
    width = len(master)
    row = list(master.replace("-", "N"))  # a row with a base in every column

    if query_pattern == "internal":
        row[15:18] = "---"
    elif query_pattern == "leading":
        row[:6] = "-" * 6
    elif query_pattern == "trailing":
        row[-6:] = "-" * 6
    elif query_pattern == "both":
        row[:4] = "-" * 4
        row[-4:] = "-" * 4
    query = "".join(row)
    assert len(query) == width

    assert_round_trip(master, query, cds_list, span=covered_span(master, query))


# ---------------------------------------------------------------------------
# 4. degenerate rows
# ---------------------------------------------------------------------------


def test_row_of_pure_gaps_yields_no_features():
    """An all-gap row covers nothing, so it must claim no coordinates."""
    master, cds_list = build_master(3, at_master_coord=1)
    query = "-" * len(master)

    first, last = covered_span(master, query)
    assert (first, last) == (None, None)

    # With no covered span there is nothing to clamp against; the entries that
    # come back must still be in range rather than negative or zero.
    for entry in project(master, query, cds_list):
        assert entry["og_start"] >= 0
        assert entry["og_end"] >= entry["og_start"] - 1


def test_feature_outside_the_covered_span_is_dropped_not_invented():
    """A partial row must not be credited with a gene it never reached."""
    genome = ORF + ORF
    master = genome
    cds_list = [
        {"start": "1", "end": str(len(ORF)), "product": "FIRST"},
        {"start": str(len(ORF) + 1), "end": str(len(genome)), "product": "SECOND"},
    ]
    query = ORF + "-" * len(ORF)  # only the first gene is present

    entries = project(master, query, cds_list, span=covered_span(master, query))
    assert [e["product"] for e in entries] == ["FIRST"]


def test_overlapping_features_each_keep_their_own_frame():
    """PA / PA-X style overlap: two products sharing columns must not merge."""
    master, _ = build_master(6, at_master_coord=10)
    cds_list = [
        {"start": "1", "end": "18", "product": "LONG"},
        {"start": "4", "end": "15", "product": "SHORT"},
    ]
    entries = {e["product"]: e for e in project(master, master, cds_list)}

    assert entries["LONG"]["og_start"] == 1 and entries["LONG"]["og_end"] == 18
    assert entries["SHORT"]["og_start"] == 4 and entries["SHORT"]["og_end"] == 15
    assert entries["LONG"]["start"] != entries["SHORT"]["start"]


def test_feature_at_the_extreme_first_and_last_columns():
    master, _ = build_master(3, at_master_coord=1)
    genome_len = len(ORF)
    cds_list = [
        {"start": "1", "end": "3", "product": "HEAD"},
        {"start": str(genome_len - 2), "end": str(genome_len), "product": "TAIL"},
    ]
    entries = {e["product"]: e for e in project(master, master, cds_list)}

    assert (entries["HEAD"]["og_start"], entries["HEAD"]["og_end"]) == (1, 3)
    assert (entries["TAIL"]["og_start"], entries["TAIL"]["og_end"]) == (
        genome_len - 2,
        genome_len,
    )


# ---------------------------------------------------------------------------
# 5. gap-boundary clamping
# ---------------------------------------------------------------------------


def test_feature_end_inside_a_gap_clamps_to_the_last_real_base():
    """Backwards clamping is correct at the END of a feature."""
    proc = processor()
    gaps = proc.get_gap_ranges("ATG-CGTAA")
    assert gaps == [[4, 4]]
    assert proc.count_gaps_before_position(gaps, 3) == 0
    assert proc.count_gaps_before_position(gaps, 4) == 1
    assert proc.count_gaps_before_position(gaps, 9) == 1


def test_feature_start_inside_an_internal_gap_clamps_forwards():
    """og_start clamps FORWARDS, to the first residue at or after the column.

    It used to share count_gaps_before_position with og_end, which clamps
    backwards. That is right for an end and wrong for a start: a CDS beginning
    inside one of this row's gaps resolved to the last residue BEFORE the
    feature, so the raw slice began early and the reading frame shifted.
    """
    master, _ = build_master(0, at_master_coord=1)
    cds_list = [{"start": "4", "end": "18", "product": "ORF"}]
    query = list(master)
    query[3:6] = "---"  # an internal deletion covering the CDS start
    query = "".join(query)

    entry = project(master, query, cds_list)[0]
    from_raw = degap(query)[entry["og_start"] - 1:entry["og_end"]]
    from_alignment = degap(query[entry["start"] - 1:entry["end"]])
    assert from_raw == from_alignment


def test_leading_padding_over_the_cds_start_never_produces_zero():
    """og_start cannot go below 1, unclamped by any span.

    It could: a row whose leading padding covered the CDS start column produced
    0, which is not a valid 1-based coordinate. The span clamp hid it in the
    pipeline, because the span begins at the first covered master coordinate,
    but the arithmetic itself was unguarded. It is now guarded by construction -
    only whole gap columns PRECEDING the position are subtracted.
    """
    master, cds_list = build_master(0, at_master_coord=1)
    query = "-" * 6 + degap(master)[6:]

    entry = project(master, query, cds_list)[0]  # deliberately unclamped
    assert entry["og_start"] >= 1


# The 5' trim is exercised in test_alignment_coordinate_contract.py, because
# correcting it needs the master's full record and therefore the whole script.


# ---------------------------------------------------------------------------
# 6. fuzz
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("seed", range(40))
def test_fuzz_projection_invariants_hold(seed):
    """Random masters and rows; the invariants must not depend on the shape.

    Checks the three things that are true by construction for every row:
    coordinates are ordered, they are inside the row's own sequence, and the
    raw slice matches the alignment slice wherever both boundaries carry a base.
    """
    rng = random.Random(seed)
    genome = "".join(rng.choice("ACGT") for _ in range(rng.randrange(30, 90)))

    # master row: scatter insertion columns through the genome
    master_chars = []
    for base in genome:
        if rng.random() < 0.12:
            master_chars.append("-" * rng.randrange(1, 4))
        master_chars.append(base)
    master = "".join(master_chars)
    width = len(master)

    # query row: a base wherever it likes, gaps elsewhere
    query = "".join(
        "-" if rng.random() < 0.2 else rng.choice("ACGT") for _ in range(width)
    )

    start = rng.randrange(1, len(genome))
    end = rng.randrange(start, len(genome) + 1)
    cds_list = [{"start": str(start), "end": str(end), "product": "P"}]
    span = covered_span(master, query)
    if span == (None, None):
        pytest.skip("row covers nothing")

    raw = degap(query)
    for entry in project(master, query, cds_list, span=span):
        assert entry["og_start"] <= entry["og_end"] + 1
        assert entry["og_end"] <= len(raw)
        assert 1 <= entry["start"] <= entry["end"] <= width
        if query[entry["start"] - 1] != "-" and query[entry["end"] - 1] != "-":
            assert (
                raw[entry["og_start"] - 1:entry["og_end"]]
                == degap(query[entry["start"] - 1:entry["end"]])
            )


@pytest.mark.parametrize("seed", range(20))
def test_fuzz_master_row_is_always_its_own_genome(seed):
    """However the master is gapped, its own og coordinates are its genome."""
    rng = random.Random(1000 + seed)
    genome = "".join(rng.choice("ACGT") for _ in range(rng.randrange(30, 120)))

    master_chars = []
    for base in genome:
        if rng.random() < 0.15:
            master_chars.append("-" * rng.randrange(1, 5))
        master_chars.append(base)
    master = "".join(master_chars)

    start = rng.randrange(1, len(genome))
    end = rng.randrange(start, len(genome) + 1)
    cds_list = [{"start": str(start), "end": str(end), "product": "P"}]

    entry = project(master, master, cds_list)[0]
    assert (entry["og_start"], entry["og_end"]) == (start, end)
    assert genome[entry["og_start"] - 1:entry["og_end"]] == genome[start - 1:end]
