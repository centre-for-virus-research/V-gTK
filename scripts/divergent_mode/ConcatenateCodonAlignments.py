#!/usr/bin/env python3
"""Concatenate per-gene codon alignments into a partitioned supermatrix.

Taxa missing a gene are padded with gaps for that block. A partition file is
written alongside so IQ-TREE can fit per-gene rates, and gene boundaries are
recorded so a column in the supermatrix can be mapped back to a gene and a
codon index for mutation reporting.
"""
import argparse
import csv
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--codon-dir", required=True)
    ap.add_argument("--genes", required=True, help="comma-separated group names, in order")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--prefix", default="backbone")
    ap.add_argument("--min-gene-coverage", type=float, default=0.0,
                    help="drop a taxon present in fewer than this fraction of genes")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]

    blocks = {}
    for gene in genes:
        path = os.path.join(args.codon_dir, f"{gene}.codon.fasta")
        recs = {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}
        widths = {len(s) for s in recs.values()}
        if len(widths) != 1:
            raise SystemExit(f"{gene}: ragged alignment, widths {sorted(widths)}")
        blocks[gene] = (recs, widths.pop())

    taxa = sorted(set().union(*[set(r.keys()) for r, _ in blocks.values()]))
    if args.min_gene_coverage > 0:
        keep = [t for t in taxa
                if sum(1 for g in genes if t in blocks[g][0]) / len(genes)
                >= args.min_gene_coverage]
        dropped = sorted(set(taxa) - set(keep))
        taxa = keep
        if dropped:
            print(f"dropped {len(dropped)} taxa below gene coverage: "
                  f"{', '.join(dropped[:8])}{'...' if len(dropped) > 8 else ''}")

    records = []
    for taxon in taxa:
        pieces = []
        for gene in genes:
            recs, width = blocks[gene]
            pieces.append(recs.get(taxon, "-" * width))
        records.append(SeqRecord(Seq("".join(pieces)), id=taxon, description=""))

    aln_path = os.path.join(args.outdir, f"{args.prefix}.codon.fasta")
    SeqIO.write(records, aln_path, "fasta")

    part_path = os.path.join(args.outdir, f"{args.prefix}.partitions.nex")
    coord_rows = []
    offset = 0
    with open(part_path, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for gene in genes:
            width = blocks[gene][1]
            fh.write(f"  charset {gene} = {offset + 1}-{offset + width};\n")
            coord_rows.append({"gene": gene, "nt_start": offset + 1,
                               "nt_end": offset + width,
                               "codons": width // 3,
                               "n_taxa_present": len(blocks[gene][0])})
            offset += width
        fh.write("end;\n")

    map_path = os.path.join(args.outdir, f"{args.prefix}.gene_coordinates.tsv")
    with open(map_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(coord_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(coord_rows)

    print(f"{len(records)} taxa x {offset} nt ({offset // 3} codons) over {len(genes)} genes")
    for r in coord_rows:
        print(f"  {r['gene']:<10} {r['nt_start']}-{r['nt_end']} "
              f"({r['codons']} codons, {r['n_taxa_present']} taxa)")
    print(f"alignment  -> {aln_path}")
    print(f"partitions -> {part_path}")
    print(f"coords     -> {map_path}")


if __name__ == "__main__":
    main()
