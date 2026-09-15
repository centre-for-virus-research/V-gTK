"""Genome shapes this pipeline has never been run against.

RABV, HCV and influenza are all linear, plus-strand-annotated records whose
genes are single, non-overlapping, non-wrapping CDS features. Almost every
coordinate assumption in the projection holds *because* of that, silently. This
module runs the real projection against the shapes a next virus would bring:

  * circular genomes, where a CDS wraps the origin and its GFF start is GREATER
    than its end (HBV, and every other circular genome with a gene over the
    join);
  * spliced CDSs, where one product arrives as two GFF lines - which influenza
    already has, in M2 and NEP;
  * ambisense genomes, where genes sit on BOTH strands (arenaviruses,
    bunyaviruses) and a minus-strand gene must be reverse-complemented before it
    means anything;
  * records annotated with no CDS at all;
  * overlapping reading frames.

Where the pipeline cannot yet handle a shape the test is an xfail carrying what
actually happens, so the limitation is recorded rather than discovered later by
someone's database. These are not hypothetical: an ambisense or circular virus
is a plausible next target and nothing in the code would report the problem.
"""

import csv
from pathlib import Path

import pytest
from Bio.Seq import Seq

from CalcAlignmentCord import CalculateAlignmentCoordinates
from GffToDictionary import GffDictionary


ORF = "ATG" + "AAACCCGGGTTTAAGCCGGTCACGTGCAGT" + "TAA"
PROTEIN = "MKPGFKPVTCS*"


def degap(text):
    return text.replace("-", "")


def run_projection(tmp_path: Path, master_row: str, rows: dict, gff_lines,
                   master_record=None):
    """Run the real entry point over one alignment file and return features.tsv."""
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir(exist_ok=True)
    fasta = [">MASTERX", master_row]
    for name, row in rows.items():
        fasta += [f">{name}", row]
    (alignment_dir / "MASTERX.aligned_merged_MSA.fasta").write_text(
        "\n".join(fasta) + "\n", encoding="utf-8"
    )

    gff_path = tmp_path / "MASTERX.gff3"
    gff_path.write_text("##gff-version 3\n" + "\n".join(gff_lines) + "\n", encoding="utf-8")

    master_list = tmp_path / "master_list.tsv"
    master_list.write_text("MASTERX\n", encoding="utf-8")
    hits = tmp_path / "hits.tsv"
    hits.write_text("".join(f"{n}\tMASTERX\t99.0\tplus\n" for n in rows), encoding="utf-8")

    master_seq_dir = None
    if master_record is not None:
        seq_dir = tmp_path / "master_seq"
        seq_dir.mkdir(exist_ok=True)
        (seq_dir / "MASTERX.fasta").write_text(f">MASTERX\n{master_record}\n", encoding="utf-8")
        master_seq_dir = str(seq_dir)

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


def cds_line(start, end, product, strand="+", ident=None):
    ident = ident or f"cds-{product}-{start}"
    return f"MASTERX\tRefSeq\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={ident};product={product}"


# ---------------------------------------------------------------------------
# 1. circular genomes
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    reason="A CDS that wraps the origin of a circular genome has GFF start > end. "
           "recalculate_cds_coordinates_with_span computes overlap_start/overlap_end "
           "and drops the feature on 'overlap_start > overlap_end' - the same branch "
           "that legitimately discards a feature lying outside a partial sequence's "
           "span. So a wrapping gene is silently absent from features.tsv with no "
           "diagnostic. Affects any circular virus (HBV, and circular DNA phages).",
)
def test_circular_genome_cds_wrapping_the_origin_is_not_silently_dropped(tmp_path):
    genome = ORF + "TTTTTTTTT"          # 42 nt, circular
    master_row = genome
    # A gene running from position 37, over the origin, to position 6.
    gff = [cds_line(37, 6, "WRAPPING")]

    features = run_projection(tmp_path, master_row, {}, gff, master_record=genome)
    products = [r["product"] for r in features]
    assert "WRAPPING" in products, "a wrapping CDS must be represented, or rejected loudly"


def test_circular_genome_non_wrapping_genes_are_unaffected(tmp_path):
    """The ordinary genes of a circular genome must still project normally."""
    genome = ORF + "TTTTTTTTT"
    gff = [cds_line(1, len(ORF), "NORMAL")]

    features = run_projection(tmp_path, genome, {}, gff, master_record=genome)
    entry = next(r for r in features if r["product"] == "NORMAL")

    assert (entry["cds_start_OG_seq"], entry["cds_end_OG_seq"]) == ("1", str(len(ORF)))
    assert genome[0:len(ORF)] == ORF


@pytest.mark.xfail(
    strict=True,
    reason="Nothing in the pipeline records that a genome is circular. meta_data "
           "carries a 'topology' column parsed by GenBankParser, but the coordinate "
           "projection never consults it, so it cannot tell a wrapping gene from a "
           "malformed annotation and treats both as the latter.",
)
def test_projection_knows_about_topology(tmp_path):
    genome = ORF + "TTTTTTTTT"
    gff = [cds_line(1, len(ORF), "NORMAL")]
    features = run_projection(tmp_path, genome, {}, gff, master_record=genome)
    assert any("topology" in key.lower() for key in features[0])


# ---------------------------------------------------------------------------
# 2. spliced CDS - influenza already has this
# ---------------------------------------------------------------------------


def test_spliced_cds_currently_becomes_two_independent_features(tmp_path):
    """Two exons of one product arrive as two rows, not one joined feature.

    This is the real influenza M2/NEP shape: NC_002016 annotates M2 as
    ``26..51`` plus ``740..1007`` under a single ID. The projection has no
    concept of a joined feature, so it emits both spans under the same product
    name. Documented rather than asserted as correct - see the xfail below for
    what it costs.
    """
    genome = ORF + "GGGGGG" + ORF
    gff = [
        cds_line(1, 18, "SPLICED", ident="cds-SPLICED"),
        cds_line(40, 72, "SPLICED", ident="cds-SPLICED"),
    ]

    features = run_projection(tmp_path, genome, {}, gff, master_record=genome)
    spliced = [r for r in features if r["product"] == "SPLICED"]

    assert len(spliced) == 2
    spans = sorted((r["cds_start_OG_seq"], r["cds_end_OG_seq"]) for r in spliced)
    assert spans == [("1", "18"), ("40", "72")]


@pytest.mark.xfail(
    strict=True,
    reason="Two exons of one product are emitted as two features with the same "
           "name. AnnotateMutations.choose_feature_entry breaks a same-priority tie "
           "by preferring the SHORTEST span, so a spliced product resolves to its "
           "smaller exon: influenza M2 would resolve to 26..51, 26 nt, and every "
           "catalogue position past codon 8 is then discarded as "
           "'catalog_position_past_feature_end'. The joined coding length is never "
           "represented anywhere.",
)
def test_spliced_cds_should_expose_the_joined_coding_length(tmp_path):
    genome = ORF + "GGGGGG" + ORF
    gff = [
        cds_line(1, 18, "SPLICED", ident="cds-SPLICED"),
        cds_line(40, 72, "SPLICED", ident="cds-SPLICED"),
    ]
    features = run_projection(tmp_path, genome, {}, gff, master_record=genome)
    spliced = [r for r in features if r["product"] == "SPLICED"]

    coding_length = sum(
        int(r["cds_end_OG_seq"]) - int(r["cds_start_OG_seq"]) + 1 for r in spliced
    )
    # 18 + 33 = 51 nt of coding sequence, which no single row reports.
    assert any(
        int(r["cds_end_OG_seq"]) - int(r["cds_start_OG_seq"]) + 1 == coding_length
        for r in spliced
    )


# ---------------------------------------------------------------------------
# 3. ambisense / minus-strand genes
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    reason="GffToDictionary unpacks the strand column and then throws it away: the "
           "entry it builds is {start, end, product}. An ambisense genome "
           "(arenavirus, bunyavirus, and the plant geminiviruses) has genes on both "
           "strands, and a minus-strand gene's residues are the reverse complement "
           "of the columns its coordinates name. Nothing downstream can discover "
           "this, because the strand never leaves the GFF parser.",
)
def test_gff_dictionary_preserves_strand():
    import tempfile, os
    handle, path = tempfile.mkstemp(suffix=".gff3")
    os.close(handle)
    try:
        with open(path, "w", encoding="utf-8") as f:
            f.write(
                "##gff-version 3\n"
                "SEQ\tRefSeq\tCDS\t1\t33\t.\t-\t0\tID=cds1;product=MINUS\n"
            )
        entry = GffDictionary(path).gff_dict["CDS"][0]
        assert entry.get("strand") == "-"
    finally:
        os.unlink(path)


@pytest.mark.xfail(
    strict=True,
    reason="A minus-strand CDS is projected exactly like a plus-strand one, so "
           "slicing the row between its coordinates yields the reverse complement "
           "of the gene. It still translates - into a different, wrong protein - so "
           "nothing raises and no QC catches it.",
)
def test_minus_strand_cds_translates_to_its_protein(tmp_path):
    minus_gene = str(Seq(ORF).reverse_complement())
    genome = minus_gene
    gff = [cds_line(1, len(genome), "MINUS", strand="-")]

    features = run_projection(tmp_path, genome, {}, gff, master_record=genome)
    entry = next(r for r in features if r["product"] == "MINUS")
    cds = genome[int(entry["cds_start_OG_seq"]) - 1:int(entry["cds_end_OG_seq"])]

    assert str(Seq(cds).translate()) == PROTEIN


# ---------------------------------------------------------------------------
# 4. degenerate annotations
# ---------------------------------------------------------------------------


def test_gff_with_no_cds_at_all_raises_rather_than_producing_empty_output(tmp_path):
    """A GFF with only gene/region lines. Currently a bare KeyError."""
    genome = ORF
    gff = ["MASTERX\tRefSeq\tregion\t1\t33\t.\t+\t.\tID=r1",
           "MASTERX\tRefSeq\tgene\t1\t33\t.\t+\t.\tID=g1;gene=G"]

    with pytest.raises(KeyError):
        run_projection(tmp_path, genome, {}, gff, master_record=genome)


@pytest.mark.xfail(
    strict=True,
    reason="A GFF carrying no CDS fails with a bare KeyError('CDS') out of "
           "find_gaps_in_fasta, naming neither the file nor the master. Every other "
           "input error in this script raises a message that says what was wrong "
           "with which input.",
)
def test_gff_with_no_cds_names_the_offending_file(tmp_path):
    genome = ORF
    gff = ["MASTERX\tRefSeq\tgene\t1\t33\t.\t+\t.\tID=g1;gene=G"]
    with pytest.raises(Exception, match="CDS.*MASTERX|MASTERX.*CDS"):
        run_projection(tmp_path, genome, {}, gff, master_record=genome)


def test_single_codon_cds(tmp_path):
    genome = ORF
    features = run_projection(tmp_path, genome, {}, [cds_line(1, 3, "TINY")],
                              master_record=genome)
    entry = next(r for r in features if r["product"] == "TINY")
    assert (entry["cds_start_OG_seq"], entry["cds_end_OG_seq"]) == ("1", "3")


def test_cds_running_past_the_end_of_the_record(tmp_path):
    """A mis-annotated GFF must not produce coordinates off the end of the row."""
    genome = ORF
    features = run_projection(tmp_path, genome, {}, [cds_line(1, 900, "TOOLONG")],
                              master_record=genome)
    for entry in features:
        assert int(entry["cds_end_OG_seq"]) <= len(genome), entry
        assert int(entry["cds_start_OG_seq"]) >= 1, entry


def test_overlapping_reading_frames_stay_separate(tmp_path):
    """HBV-style: four genes overlapping in different frames, all preserved."""
    genome = ORF + ORF
    gff = [
        cds_line(1, 33, "FRAME0"),
        cds_line(2, 34, "FRAME1"),
        cds_line(3, 35, "FRAME2"),
        cds_line(10, 60, "LONG"),
    ]
    features = run_projection(tmp_path, genome, {}, gff, master_record=genome)
    spans = {r["product"]: (r["cds_start_OG_seq"], r["cds_end_OG_seq"])
             for r in features if r["accession"] == "MASTERX"}

    assert spans["FRAME0"] == ("1", "33")
    assert spans["FRAME1"] == ("2", "34")
    assert spans["FRAME2"] == ("3", "35")
    assert spans["LONG"] == ("10", "60")


def test_master_row_that_is_entirely_gaps(tmp_path):
    """A master contributing no residues must not yield invalid coordinates."""
    genome = ORF
    master_row = "-" * len(ORF)
    features = run_projection(tmp_path, master_row, {"Q": ORF},
                              [cds_line(1, 33, "P")], master_record=genome)
    for entry in features:
        assert int(entry["cds_start_OG_seq"]) >= 0, entry


def test_duplicate_accession_in_one_alignment_file(tmp_path):
    """The same id twice in one FASTA must not silently double every feature."""
    genome = ORF
    alignment_dir = tmp_path / "padded_alignment"
    alignment_dir.mkdir()
    (alignment_dir / "MASTERX.aligned_merged_MSA.fasta").write_text(
        f">MASTERX\n{genome}\n>Q_DUP\n{genome}\n>Q_DUP\n{genome}\n", encoding="utf-8"
    )
    gff_path = tmp_path / "MASTERX.gff3"
    gff_path.write_text(
        "##gff-version 3\n" + cds_line(1, 33, "P") + "\n", encoding="utf-8"
    )
    master_list = tmp_path / "m.tsv"
    master_list.write_text("MASTERX\n", encoding="utf-8")
    hits = tmp_path / "h.tsv"
    hits.write_text("Q_DUP\tMASTERX\t99.0\tplus\n", encoding="utf-8")

    processor = CalculateAlignmentCoordinates(
        paded_alignment=str(alignment_dir), master_gff=[str(gff_path)],
        tmp_dir=str(tmp_path), output_dir="Tables", output_file="features.tsv",
        master_accession=str(master_list), blast_uniq_hits=str(hits),
    )
    processor.find_gaps_in_fasta()
    with (tmp_path / "Tables" / "features.tsv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    dup_rows = [r for r in rows if r["accession"] == "Q_DUP"]
    # Recorded so a change in dedup behaviour is visible; both copies are identical.
    assert len({(r["cds_start_OG_seq"], r["cds_end_OG_seq"]) for r in dup_rows}) == 1
