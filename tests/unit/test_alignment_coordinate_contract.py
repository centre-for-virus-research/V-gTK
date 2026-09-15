"""The contract between CalcAlignmentCord's output and the things that read it.

test_alignment_coordinate_integrity.py checks the projection arithmetic in
isolation. This module checks the two things that arithmetic exists to serve:

  * the whole script, run end to end over real files with a **gapped master** -
    every fixture in test_calc_alignment_cord.py uses an ungapped master, where
    a master genome coordinate and an alignment column are the same number, so
    none of them can see a space-1/space-2 confusion at all;

  * the cross-script contract with ``AnnotateMutations``, which takes
    ``cds_start_OG_seq`` from the features table and resolves catalogue residues
    through ``build_alignment_coordinate_map``. That map is keyed by the
    master's *ungapped* position, so the two scripts have to agree on what the
    OG columns mean. Nothing else in the suite pins that agreement down, and it
    is exactly what drifted.

Plus a set of deliberately awkward simulated rows - ambiguity codes, lower
case, runs of N, insertions that are not codon-sized - because the projection
is pure column arithmetic and must not care what the residues actually are.
"""

import csv
from pathlib import Path

import pytest

import AnnotateMutations
from CalcAlignmentCord import CalculateAlignmentCoordinates


# A 12-codon ORF, as in the sibling module.
ORF = "ATG" + "AAACCCGGGTTTAAGCCGGTCACGTGCAGT" + "TAA"
PROTEIN = "MKPGFKPVTCS*"


def degap(text):
    return text.replace("-", "")


def coord_map(master_alignment):
    mapping, residue = {}, 0
    for aln_pos, base in enumerate(master_alignment, start=1):
        if base != "-":
            residue += 1
            mapping[residue] = aln_pos
    return mapping


def build_case(tmp_path: Path, master_row: str, rows: dict, cds, master_record=None,
               region_length=None):
    """Write a padded alignment + GFF + ref list, and run the real entry point.

    ``cds`` is a list of (start, end, product) in MASTER GENOME coordinates.
    ``master_record``, when given, is the master's full submitted sequence and is
    written to a master_seq directory, which is how the script locates any 5'
    trim between the GFF's record coordinates and the master's alignment row.
    Returns the parsed features.tsv rows.
    """
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir(exist_ok=True)
    fasta = [">MASTERG", master_row]
    for name, row in rows.items():
        fasta += [f">{name}", row]
    (alignment_dir / "MASTERG.aligned_merged_MSA.fasta").write_text(
        "\n".join(fasta) + "\n", encoding="utf-8"
    )

    gff = ["##gff-version 3"]
    if region_length is None and master_record is not None:
        region_length = len(master_record)
    if region_length is not None:
        gff.append(
            f"MASTERG\tRefSeq\tregion\t1\t{region_length}\t.\t+\t.\tID=MASTERG:1..{region_length}"
        )
    for index, (start, end, product) in enumerate(cds, start=1):
        gff.append(
            f"MASTERG\tRefSeq\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds{index};product={product}"
        )
    gff_path = tmp_path / "MASTERG.gff3"
    gff_path.write_text("\n".join(gff) + "\n", encoding="utf-8")

    master_seq_dir = None
    if master_record is not None:
        seq_dir = tmp_path / "master_seq"
        seq_dir.mkdir(exist_ok=True)
        (seq_dir / "MASTERG.fasta").write_text(
            f">MASTERG\n{master_record}\n", encoding="utf-8"
        )
        master_seq_dir = str(seq_dir)

    master_list = tmp_path / "master_list.tsv"
    master_list.write_text("MASTERG\n", encoding="utf-8")

    hits = tmp_path / "query_uniq_tophits.tsv"
    hits.write_text(
        "".join(f"{name}\tMASTERG\t99.0\tplus\n" for name in rows), encoding="utf-8"
    )

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir),
        master_gff=[str(gff_path)],
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(master_list),
        blast_uniq_hits=str(hits),
        master_seq_dir=master_seq_dir,
    )
    processor.find_gaps_in_fasta()

    with (tmp_path / "Tables" / "features.tsv").open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---------------------------------------------------------------------------
# 1. the whole script, with a gapped master
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("gaps", [3, 6, 9, 30])
def test_end_to_end_gapped_master_writes_its_own_genome_coordinates(tmp_path, gaps):
    """The master's row of features.tsv must carry its plain GFF coordinates."""
    master_row = "-" * gaps + ORF
    rows = {"Q_INS": "N" * gaps + ORF, "Q_NOINS": "-" * gaps + ORF}

    features = build_case(tmp_path, master_row, rows, [(1, len(ORF), "ORF")])
    master_feature = next(r for r in features if r["accession"] == "MASTERG")

    assert master_feature["cds_start_OG_seq"] == "1"
    assert master_feature["cds_end_OG_seq"] == str(len(ORF))
    # cds_start/cds_end are alignment columns, so they DO shift with the gaps.
    assert master_feature["cds_start"] == str(gaps + 1)
    assert master_feature["cds_end"] == str(gaps + len(ORF))


def test_end_to_end_every_row_translates_without_an_internal_stop(tmp_path):
    """The symptom that started this: CDS slices that would not translate."""
    from Bio.Seq import Seq

    gaps = 3
    master_row = "-" * gaps + ORF
    rows = {
        "Q_INS": "TGA" + ORF,      # carries the insertion, which is itself a stop
        "Q_NOINS": "-" * gaps + ORF,
    }
    raw = {name: degap(row) for name, row in rows.items()}
    raw["MASTERG"] = degap(master_row)

    features = build_case(tmp_path, master_row, rows, [(1, len(ORF), "ORF")])

    for row in features:
        start = int(row["cds_start_OG_seq"])
        end = int(row["cds_end_OG_seq"])
        assert start >= 1, row
        cds = raw[row["accession"]][start - 1:end]
        assert cds.startswith("ATG"), row
        protein = str(Seq(cds).translate())
        assert protein == PROTEIN, row
        assert "*" not in protein[:-1], row


def test_end_to_end_partial_row_is_clamped_and_stays_in_range(tmp_path):
    """A row covering only the back half must not claim the whole gene."""
    gaps = 6
    master_row = "-" * gaps + ORF + ORF
    half = len(ORF)
    rows = {"Q_BACK": "-" * (gaps + half) + ORF}

    features = build_case(
        tmp_path,
        master_row,
        rows,
        [(1, half, "FIRST"), (half + 1, 2 * half, "SECOND")],
    )
    back = [r for r in features if r["accession"] == "Q_BACK"]

    assert [r["product"] for r in back] == ["SECOND"]
    entry = back[0]
    assert int(entry["cds_start_OG_seq"]) >= 1
    assert degap(rows["Q_BACK"])[
        int(entry["cds_start_OG_seq"]) - 1:int(entry["cds_end_OG_seq"])
    ] == ORF


# ---------------------------------------------------------------------------
# 2. the cross-script contract with AnnotateMutations
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("gaps", [0, 3, 6, 9, 129])
def test_og_coords_resolve_correct_codons_through_annotate_mutations(tmp_path, gaps):
    """CalcAlignmentCord's OG output must be usable by AnnotateMutations as-is.

    This is the join the defect broke. AnnotateMutations reads
    ``cds_start_OG_seq`` and looks up ``cds_start + (aa_pos-1)*3`` in
    ``build_alignment_coordinate_map(master_alignment)``, which is keyed by the
    master's ungapped position. If the OG value is short by the master's gap
    count every catalogue position lands on the wrong codon, and the lookup
    still succeeds - so nothing raises and the residues are merely wrong.
    """
    master_row = "-" * gaps + ORF
    rows = {"Q": "N" * gaps + ORF}
    features = build_case(tmp_path, master_row, rows, [(1, len(ORF), "ORF")])

    master_feature = next(r for r in features if r["accession"] == "MASTERG")
    cds_start = int(master_feature["cds_start_OG_seq"])

    master_coord_map = AnnotateMutations.build_alignment_coordinate_map(master_row)

    for aa_pos, expected_residue in enumerate(PROTEIN, start=1):
        indices = AnnotateMutations.resolve_aligned_codon_indices(
            master_coord_map, cds_start, aa_pos
        )
        assert indices is not None, (aa_pos, gaps)
        codon = AnnotateMutations.extract_aligned_codon(master_row, indices)
        assert AnnotateMutations.translate_codon(codon) == expected_residue, (
            aa_pos,
            gaps,
            codon,
        )


def test_the_previous_arithmetic_would_fail_that_contract(tmp_path):
    """Guard the guard: show the contract test above is actually load-bearing.

    Recomputes the OG start the way the defect did - master genome coordinate
    minus the row's gap count in alignment columns - and asserts it resolves the
    WRONG residue. If this ever starts agreeing, the contract test above has
    stopped discriminating and needs rewriting.
    """
    gaps = 3
    master_row = "-" * gaps + ORF
    processor = CalculateAlignmentCoordinates(None, None, None, None, None, None, None)

    gap_ranges = processor.get_gap_ranges(master_row)
    buggy_start = 1 - processor.count_gaps_before_position(gap_ranges, 1)
    assert buggy_start == 0  # not even a valid 1-based coordinate

    # Take the shape the defect actually expressed in the field: a CDS starting
    # after the insertion, which is the influenza NS1/M1 case. Master genome
    # coordinate 7 sits at alignment column 10, because columns 7-9 are the
    # inserted ones.
    master_row2 = ORF[:6] + "-" * gaps + ORF[6:]
    gap_ranges2 = processor.get_gap_ranges(master_row2)
    aln_column = coord_map(master_row2)[7]
    assert aln_column == 10

    correct = aln_column - processor.count_gaps_before_position(gap_ranges2, aln_column)
    buggy = 7 - processor.count_gaps_before_position(gap_ranges2, 7)
    assert correct == 7, "the master's own row must report its genome coordinate"
    assert buggy == 6, "the previous arithmetic lost one base per gap it crossed"

    master_coord_map = AnnotateMutations.build_alignment_coordinate_map(master_row2)

    good = AnnotateMutations.resolve_aligned_codon_indices(master_coord_map, correct, 1)
    assert AnnotateMutations.extract_aligned_codon(master_row2, good) == ORF[6:9]

    # The same lookup with the old value succeeds - it just answers with the
    # wrong codon, which is why this never surfaced as an error.
    bad = AnnotateMutations.resolve_aligned_codon_indices(master_coord_map, buggy, 1)
    assert bad is not None
    assert AnnotateMutations.extract_aligned_codon(master_row2, bad) != ORF[6:9]


# ---------------------------------------------------------------------------
# 2b. the 5' trim between a master's GFF and its alignment row
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("utr_len", [0, 1, 9, 26, 69])
def test_five_prime_trim_between_gff_and_alignment_row_is_corrected(tmp_path, utr_len):
    """A master's row can start partway into its record; the GFF cannot.

    The guide alignment's copy of a reference is frequently trimmed at the 5'
    end - by 27, 10, 24, 69, 45, 25 and 26 bases on influenza segments 1, 2, 3,
    4, 5, 7 and 8, with segment 6 alone flush. The GFF is in RECORD coordinates,
    so counting residues down the trimmed row numbers them differently and every
    annotated coordinate lands that many bases out.

    Note what the OG columns must contain here: the row's OWN numbering, which
    for a trimmed master is NOT the record coordinate. AnnotateMutations indexes
    the alignment as stored, and what is stored is the trimmed row.
    """
    record = "G" * utr_len + ORF  # the submitted record: UTR, then the ORF
    master_row = ORF[:6] + "---" + ORF[6:]  # the row: UTR trimmed, 3 gap columns
    assert degap(master_row) == ORF

    cds_list = [(utr_len + 1, utr_len + len(ORF), "ORF")]
    features = build_case(tmp_path, master_row, {}, cds_list, master_record=record)
    master_feature = next(r for r in features if r["accession"] == "MASTERG")

    # OG is the row's own numbering: the whole row, 1..33.
    assert (master_feature["cds_start_OG_seq"], master_feature["cds_end_OG_seq"]) == (
        "1",
        str(len(ORF)),
    )
    assert degap(master_row)[0:len(ORF)] == ORF

    # The alignment columns must span the row, start to finish.
    assert (master_feature["cds_start"], master_feature["cds_end"]) == (
        "1",
        str(len(master_row)),
    )

    # And the mapped start column must actually hold the start codon.
    start_column = int(master_feature["cds_start"])
    assert master_row[start_column - 1:start_column + 2] == "ATG"


def test_trimmed_master_translates_cleanly_from_the_stored_row(tmp_path):
    """The end-to-end symptom: a CDS truncated by exactly the trim.

    Before the trim was located, og_start was the record coordinate while og_end
    was clamped to the row's length, so the slice ran short by the trim and lost
    the stop codon - which is why influenza HA translated as a correct-looking
    protein missing its C-terminus.
    """
    from Bio.Seq import Seq

    utr_len = 26
    record = "G" * utr_len + ORF
    master_row = ORF[:6] + "---" + ORF[6:]
    cds_list = [(utr_len + 1, utr_len + len(ORF), "ORF")]

    features = build_case(tmp_path, master_row, {}, cds_list, master_record=record)
    entry = next(r for r in features if r["accession"] == "MASTERG")

    row_sequence = degap(master_row)
    cds = row_sequence[int(entry["cds_start_OG_seq"]) - 1:int(entry["cds_end_OG_seq"])]
    assert len(cds) == len(ORF), "the CDS must not be truncated by the trim"
    assert str(Seq(cds).translate()) == PROTEIN
    assert cds.startswith("ATG") and cds.endswith("TAA")


def test_trim_is_reported_when_no_master_sequence_is_supplied(tmp_path, capfd):
    """Without the record we cannot locate a trim - so say so, don't assume 0."""
    utr_len = 9
    master_row = ORF[:6] + "---" + ORF[6:]
    cds_list = [(utr_len + 1, utr_len + len(ORF), "ORF")]

    # region says the record is longer than the row, but no master_seq_dir.
    build_case(tmp_path, master_row, {}, cds_list,
               region_length=utr_len + len(ORF))

    captured = capfd.readouterr()
    assert "trim" in (captured.err + captured.out).lower()


# ---------------------------------------------------------------------------
# 3. simulated awkward rows - the arithmetic must not care about the residues
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "label,row_builder",
    [
        ("lowercase", lambda gaps: ("n" * gaps + ORF.lower())),
        ("mixed_case", lambda gaps: ("N" * gaps + "".join(
            ch.lower() if i % 2 else ch for i, ch in enumerate(ORF)))),
        ("iupac", lambda gaps: ("N" * gaps + "RYS" + ORF[3:])),
        ("n_run", lambda gaps: ("N" * gaps + "N" * 9 + ORF[9:])),
    ],
)
def test_projection_is_indifferent_to_residue_content(tmp_path, label, row_builder):
    """Column arithmetic must not depend on what letters sit in the columns."""
    gaps = 3
    master_row = "-" * gaps + ORF
    rows = {"Q": row_builder(gaps)}
    assert len(rows["Q"]) == len(master_row), label

    features = build_case(tmp_path, master_row, rows, [(1, len(ORF), "ORF")])
    query = next(r for r in features if r["accession"] == "Q")

    # This row carries bases in the inserted columns, so in its OWN numbering the
    # CDS begins after them - og_start is gaps+1, not 1. The master's row, which
    # is gapped there, reports 1. Both are correct; that is the whole point of
    # the OG columns being per-row.
    assert (query["cds_start_OG_seq"], query["cds_end_OG_seq"]) == (
        str(gaps + 1),
        str(gaps + len(ORF)),
    )
    assert (query["cds_start"], query["cds_end"]) == (
        str(gaps + 1),
        str(gaps + len(ORF)),
    )
    master_feature = next(r for r in features if r["accession"] == "MASTERG")
    assert (master_feature["cds_start_OG_seq"], master_feature["cds_end_OG_seq"]) == (
        "1",
        str(len(ORF)),
    )


@pytest.mark.parametrize("gaps", [1, 2, 4, 5, 7, 11])
def test_insertions_that_are_not_codon_sized(tmp_path, gaps):
    """Nothing enforces codon-sized insertions in a curated guide alignment.

    Every count in generic/influenza/ref_set_aligned happens to be a multiple of
    three, which is precisely why the defect stayed in frame and went unseen. A
    non-multiple would have frameshifted loudly; the projection must be correct
    for both so that the next curated insertion cannot reintroduce it quietly.
    """
    master_row = "-" * gaps + ORF
    rows = {"Q": "N" * gaps + ORF}
    features = build_case(tmp_path, master_row, rows, [(1, len(ORF), "ORF")])
    raw = {"MASTERG": degap(master_row), "Q": degap(rows["Q"])}

    # Assert the semantics, not a fixed number: each row's OG coordinates must
    # cut the ORF out of that row's own sequence, wherever it happens to sit.
    for row in features:
        start = int(row["cds_start_OG_seq"])
        end = int(row["cds_end_OG_seq"])
        assert start >= 1, (gaps, row)
        assert raw[row["accession"]][start - 1:end] == ORF, (gaps, row)


def test_long_gap_run_does_not_overflow_the_row(tmp_path):
    """A pathologically long insertion still yields in-range coordinates."""
    gaps = 300
    master_row = "-" * gaps + ORF
    rows = {"Q_INS": "A" * gaps + ORF, "Q_NOINS": "-" * gaps + ORF}

    features = build_case(tmp_path, master_row, rows, [(1, len(ORF), "ORF")])
    raw = {"MASTERG": degap(master_row), **{k: degap(v) for k, v in rows.items()}}

    for row in features:
        start, end = int(row["cds_start_OG_seq"]), int(row["cds_end_OG_seq"])
        assert 1 <= start <= end <= len(raw[row["accession"]]), row
        assert raw[row["accession"]][start - 1:end] == ORF


def test_cds_spanning_an_insertion_boundary(tmp_path):
    """A gene straddling the inserted columns is the realistic influenza shape."""
    gaps = 6
    master_row = ORF[:15] + "-" * gaps + ORF[15:]
    rows = {"Q_INS": ORF[:15] + "CCCGGG" + ORF[15:], "Q_NOINS": master_row}

    features = build_case(tmp_path, master_row, rows, [(1, len(ORF), "ORF")])
    master_feature = next(r for r in features if r["accession"] == "MASTERG")

    assert master_feature["cds_start_OG_seq"] == "1"
    assert master_feature["cds_end_OG_seq"] == str(len(ORF))

    q = next(r for r in features if r["accession"] == "Q_INS")
    raw = degap(rows["Q_INS"])
    assert raw[int(q["cds_start_OG_seq"]) - 1:int(q["cds_end_OG_seq"])] == raw


def test_duplicate_products_do_not_collapse(tmp_path):
    """Two features sharing a product name must both survive with their spans."""
    master_row = "-" * 3 + ORF
    rows = {"Q": "N" * 3 + ORF}
    features = build_case(
        tmp_path, master_row, rows, [(1, 18, "SAME"), (19, len(ORF), "SAME")]
    )
    master_rows = [r for r in features if r["accession"] == "MASTERG"]

    spans = sorted((r["cds_start_OG_seq"], r["cds_end_OG_seq"]) for r in master_rows)
    assert spans == [("1", "18"), ("19", str(len(ORF)))]
