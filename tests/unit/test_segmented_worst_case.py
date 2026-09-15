"""The worst segmented build the projection has to survive.

An influenza build is eight independent coordinate systems processed in one
pass: eight masters, eight alignment widths, eight 5' trims, eight gap counts.
Every one of those is a chance for one segment's numbers to be applied to
another's sequences, and the result of that is not a crash - it is a database
of plausible coordinates.

The real asset measured in generic/influenza/ref_set_aligned is deliberately
reproduced here rather than simplified: masters gapped by 3, 3, 6, 132, 0, 102,
9 and 6 columns, and trimmed at the 5' end by 27, 10, 24, 69, 45, 0, 25 and 26
bases. Segment 5 is flush and ungapped, segment 4 is the pathological one, and
the suite asserts that each segment is resolved against its own master with its
own trim - which is the thing that silently was not happening.
"""

import csv
from pathlib import Path

import pytest
from Bio.Seq import Seq

from CalcAlignmentCord import CalculateAlignmentCoordinates


ORF = "ATG" + "AAACCCGGGTTTAAGCCGGTCACGTGCAGT" + "TAA"
PROTEIN = "MKPGFKPVTCS*"

# (segment, master accession, gap columns inserted, 5' bases trimmed)
SEGMENTS = [
    ("1", "KJ889313", 3, 27),
    ("2", "CY087822", 3, 10),
    ("3", "MW333949", 6, 24),
    ("4", "AB573800", 132, 69),
    ("5", "NC_007369", 0, 45),
    ("6", "AB472016", 102, 0),
    ("7", "NC_002016", 9, 25),
    ("8", "KY243296", 6, 26),
]


def degap(text):
    return text.replace("-", "")


def master_row_for(gaps):
    """The master's row: the ORF with ``gaps`` inserted columns after base 6."""
    return ORF[:6] + "-" * gaps + ORF[6:]


def record_for(trim):
    """The master's submitted record: ``trim`` bases of UTR, then the ORF."""
    return "G" * trim + ORF


def build_segmented_case(tmp_path: Path, segments=SEGMENTS, queries_for=None,
                         segment_map=None, with_master_seqs=True):
    """Write one alignment file per segment and run the real entry point once."""
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir(exist_ok=True)
    seq_dir = tmp_path / "master_seq"
    seq_dir.mkdir(exist_ok=True)
    gff_paths = []
    ref_rows = []
    hit_rows = []

    for segment, master, gaps, trim in segments:
        row = master_row_for(gaps)
        record = record_for(trim)
        fasta = [f">{master}", row]
        for name, query_row in (queries_for(segment) if queries_for else {}).items():
            fasta += [f">{name}", query_row]
            hit_rows.append(f"{name}\t{master}\t99.0\tplus\n")
        # refset_<segment>_ naming is how a segmented build resolves its master
        (alignment_dir / f"refset_{segment}_aln_merged_MSA.fasta").write_text(
            "\n".join(fasta) + "\n", encoding="utf-8"
        )

        gff_path = tmp_path / f"{master}.gff3"
        gff_path.write_text(
            "##gff-version 3\n"
            f"{master}\tRefSeq\tregion\t1\t{len(record)}\t.\t+\t.\tID=r\n"
            f"{master}\tRefSeq\tCDS\t{trim + 1}\t{trim + len(ORF)}\t.\t+\t0\tID=c;product=PROD{segment}\n",
            encoding="utf-8",
        )
        gff_paths.append(str(gff_path))
        if with_master_seqs:
            (seq_dir / f"{master}.fasta").write_text(f">{master}\n{record}\n", encoding="utf-8")
        ref_rows.append(f"{master}\tmaster\t{segment}\n")

    ref_list = tmp_path / "ref_list.tsv"
    ref_list.write_text("".join(ref_rows), encoding="utf-8")
    hits = tmp_path / "hits.tsv"
    hits.write_text("".join(hit_rows), encoding="utf-8")

    segment_map_path = None
    if segment_map is not None:
        segment_map_path = tmp_path / "segment_map.tsv"
        segment_map_path.write_text(
            "primary_accession\tsegment\n"
            + "".join(f"{acc}\t{seg}\n" for acc, seg in segment_map.items()),
            encoding="utf-8",
        )

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir),
        master_gff=gff_paths,
        tmp_dir=str(tmp_path),
        output_dir="Tables",
        output_file="features.tsv",
        master_accession=str(ref_list),
        blast_uniq_hits=str(hits),
        segment_map_tsv=str(segment_map_path) if segment_map_path else None,
        master_seq_dir=str(seq_dir) if with_master_seqs else None,
    )
    processor.find_gaps_in_fasta()
    with (tmp_path / "Tables" / "features.tsv").open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---------------------------------------------------------------------------
# 1. eight segments at once
# ---------------------------------------------------------------------------


def test_every_segment_resolves_against_its_own_master_and_trim(tmp_path):
    """Eight masters, eight gap counts, eight trims, one pass."""
    features = build_segmented_case(tmp_path)

    by_product = {r["product"]: r for r in features}
    assert len(by_product) == len(SEGMENTS), "every segment must contribute its product"

    for segment, master, gaps, trim in SEGMENTS:
        entry = by_product[f"PROD{segment}"]
        assert entry["accession"] == master
        # OG is the row's own numbering: the whole ORF, 1..33, whatever the trim.
        assert (entry["cds_start_OG_seq"], entry["cds_end_OG_seq"]) == ("1", str(len(ORF))), segment
        # The alignment columns span the row, whose width is 33 + that segment's gaps.
        assert (entry["cds_start"], entry["cds_end"]) == ("1", str(len(ORF) + gaps)), segment


def test_every_segment_translates_cleanly(tmp_path):
    """The end-to-end property: each master's CDS is a complete protein."""
    features = build_segmented_case(tmp_path)

    for segment, master, gaps, trim in SEGMENTS:
        entry = next(r for r in features if r["product"] == f"PROD{segment}")
        row = degap(master_row_for(gaps))
        cds = row[int(entry["cds_start_OG_seq"]) - 1:int(entry["cds_end_OG_seq"])]
        assert len(cds) == len(ORF), segment
        assert str(Seq(cds).translate()) == PROTEIN, segment


def test_one_segments_coordinates_never_leak_into_another(tmp_path):
    """Segment 4's 132 gaps must not widen segment 5's coordinates."""
    features = build_segmented_case(tmp_path)
    seg4 = next(r for r in features if r["product"] == "PROD4")
    seg5 = next(r for r in features if r["product"] == "PROD5")

    assert seg4["cds_end"] == str(len(ORF) + 132)
    assert seg5["cds_end"] == str(len(ORF) + 0)
    assert seg4["cds_end_OG_seq"] == seg5["cds_end_OG_seq"] == str(len(ORF))


def test_queries_in_every_segment_round_trip(tmp_path):
    """A query per segment, some carrying the insertion and some not."""
    def queries_for(segment):
        gaps = dict((s, g) for s, _m, g, _t in SEGMENTS)[segment]
        return {
            f"Q{segment}_INS": ORF[:6] + "A" * gaps + ORF[6:],   # carries it
            f"Q{segment}_DEL": master_row_for(gaps),             # lacks it
        }

    features = build_segmented_case(tmp_path, queries_for=queries_for)
    rows = {(r["accession"], r["product"]): r for r in features}

    for segment, master, gaps, trim in SEGMENTS:
        for suffix, row_seq in (("INS", ORF[:6] + "A" * gaps + ORF[6:]),
                                ("DEL", master_row_for(gaps))):
            entry = rows[(f"Q{segment}_{suffix}", f"PROD{segment}")]
            raw = degap(row_seq)
            start, end = int(entry["cds_start_OG_seq"]), int(entry["cds_end_OG_seq"])
            assert 1 <= start <= end <= len(raw), (segment, suffix, entry)
            sliced = raw[start - 1:end]
            assert sliced.startswith("ATG"), (segment, suffix)
            assert sliced.endswith("TAA"), (segment, suffix)


# ---------------------------------------------------------------------------
# 2. things that go wrong across segments
# ---------------------------------------------------------------------------


def test_segment_mismatch_between_record_and_master_is_refused(tmp_path):
    """A query mapped to the wrong segment's master must not be projected."""
    def queries_for(segment):
        return {f"Q{segment}": master_row_for(dict((s, g) for s, _m, g, _t in SEGMENTS)[segment])}

    bad_map = {master: segment for segment, master, _g, _t in SEGMENTS}
    bad_map["Q4"] = "7"          # a segment-4 row declared as segment 7
    bad_map.update({f"Q{s}": s for s, _m, _g, _t in SEGMENTS if s != "4"})

    with pytest.raises(ValueError, match="Segment mismatch"):
        build_segmented_case(tmp_path, queries_for=queries_for, segment_map=bad_map)


def test_the_same_accession_in_two_segments_keeps_separate_coordinates(tmp_path):
    """Accession reuse across segments must not collapse to one row.

    Influenza submissions do reuse an isolate name across segments, and the
    projection keys on the record id within one alignment file.
    """
    shared = "SHARED_ACC"

    def queries_for(segment):
        if segment in ("4", "7"):
            gaps = dict((s, g) for s, _m, g, _t in SEGMENTS)[segment]
            return {shared: master_row_for(gaps)}
        return {}

    features = build_segmented_case(tmp_path, queries_for=queries_for)
    shared_rows = [r for r in features if r["accession"] == shared]

    assert len(shared_rows) == 2
    products = sorted(r["product"] for r in shared_rows)
    assert products == ["PROD4", "PROD7"]
    # Different segments, different alignment widths - the columns must differ.
    assert len({r["cds_end"] for r in shared_rows}) == 2


def test_a_segment_whose_master_is_absent_from_its_alignment(tmp_path):
    """A master missing from its own alignment file must fail loudly.

    find_gaps_in_fasta falls back to ``fasta_records[0]`` when the master is not
    present, which on its own would silently project one segment's annotation
    onto whatever sequence happened to sort first. It does not get that far:
    CalculateGenomeCoordinates is called first and refuses, naming the master.
    That refusal is the safety property worth pinning - if the fallback ever
    starts being reached, this test fails.
    """
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir()
    (alignment_dir / "refset_7_aln_merged_MSA.fasta").write_text(
        f">SOMEONE_ELSE\n{master_row_for(9)}\n", encoding="utf-8"
    )
    gff = tmp_path / "NC_002016.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "NC_002016\tRefSeq\tCDS\t26\t58\t.\t+\t0\tID=c;product=M1\n",
        encoding="utf-8",
    )
    ref_list = tmp_path / "ref_list.tsv"
    ref_list.write_text("NC_002016\tmaster\t7\n", encoding="utf-8")
    hits = tmp_path / "hits.tsv"
    hits.write_text("", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir), master_gff=[str(gff)],
        tmp_dir=str(tmp_path), output_dir="Tables", output_file="features.tsv",
        master_accession=str(ref_list), blast_uniq_hits=str(hits),
    )

    with pytest.raises(ValueError, match="NC_002016"):
        processor.find_gaps_in_fasta()


def test_segment_with_no_queries_still_emits_its_master(tmp_path):
    features = build_segmented_case(tmp_path)
    assert {r["accession"] for r in features} == {m for _s, m, _g, _t in SEGMENTS}


def test_without_master_sequences_every_trimmed_segment_is_reported(tmp_path, capfd):
    """No master_seq_dir: the trim cannot be corrected, so it must be reported."""
    build_segmented_case(tmp_path, with_master_seqs=False)
    captured = capfd.readouterr()
    text = (captured.out + captured.err).lower()

    assert "trim" in text
    # Segment 6 is flush, so it has nothing to report; the other seven do.
    for segment, master, gaps, trim in SEGMENTS:
        if trim:
            assert master.lower() in text, master


def test_alignment_widths_differ_per_segment_and_stay_within_their_own(tmp_path):
    """No coordinate may exceed the width of the file it came from."""
    features = build_segmented_case(tmp_path)
    width_of = {f"PROD{s}": len(ORF) + g for s, _m, g, _t in SEGMENTS}

    for entry in features:
        assert int(entry["cds_end"]) <= width_of[entry["product"]], entry
        assert int(entry["cds_start"]) >= 1, entry
