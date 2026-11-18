#!/usr/bin/env python3
"""Codon-aware nucleotide alignment by protein-guided back-projection.

Translate each CDS, align the proteins with MAFFT, then rebuild the nucleotide
alignment by walking the protein alignment and emitting one codon per residue
and three gap characters per protein gap.

That construction is the frameshift guarantee. A gap in the output is always a
whole codon wide and always codon-phased, because the only thing that can
produce a gap is a gap in the *protein* alignment, which is one residue wide by
definition. No nucleotide aligner is involved, so none of its gap penalties can
place a 1 or 2 nt gap inside a coding sequence.

Why MAFFT rather than MACSE: MACSE is not in the `vgtk` environment, and the
codon integrity it provides is obtained here by construction instead. MACSE's
real advantage is over sequences with genuine frameshifts and sequencing
errors, where it can insert compensating single-base gaps; those sequences are
rejected upstream by ExtractDivergentCds.py instead, which is the behaviour the
brief asks for (zero frameshifts in coding regions).
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_mafft(in_faa, out_faa, threads, accurate=True):
    """MAFFT on amino acids. L-INS-i for the small, divergent sets here."""
    cmd = ["mafft", "--thread", str(threads), "--anysymbol", "--quiet"]
    if accurate:
        cmd += ["--maxiterate", "1000", "--localpair"]
    else:
        cmd += ["--auto"]
    cmd.append(in_faa)
    with open(out_faa, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"mafft failed: {proc.stderr[-400:]}")
    return out_faa


def back_project(nt_records, aln_proteins):
    """Rebuild a nucleotide alignment from an aligned protein alignment.

    Every column of the output is a codon position; every gap is 3 nt wide.
    """
    nt_by_id = {r.id: str(r.seq) for r in nt_records}
    out = []
    for prot in aln_proteins:
        nt = nt_by_id.get(prot.id)
        if nt is None:
            raise KeyError(f"no nucleotide record for {prot.id}")
        pieces = []
        cursor = 0
        for residue in str(prot.seq):
            if residue == "-":
                pieces.append("---")
                continue
            codon = nt[cursor:cursor + 3]
            if len(codon) < 3:
                codon = codon + "-" * (3 - len(codon))
            pieces.append(codon)
            cursor += 3
        out.append(SeqRecord(Seq("".join(pieces)), id=prot.id, description=""))
    return out


def check_codon_integrity(records):
    """Every sequence a codon multiple, every gap run a codon-phased triple."""
    problems = []
    width = len(records[0].seq) if records else 0
    for rec in records:
        s = str(rec.seq)
        if len(s) != width:
            problems.append((rec.id, "ragged_width", len(s)))
        if len(s) % 3 != 0:
            problems.append((rec.id, "not_triplet", len(s)))
        run = 0
        start = 0
        for i, ch in enumerate(s + "X"):
            if ch == "-":
                if run == 0:
                    start = i
                run += 1
            else:
                if run and (run % 3 != 0 or start % 3 != 0):
                    problems.append((rec.id, f"frameshift_gap@{start}", run))
                run = 0
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cds-dir", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-taxa", type=int, default=4)
    ap.add_argument("--fast", action="store_true",
                    help="MAFFT --auto instead of L-INS-i")
    ap.add_argument("--groups", default="", help="comma-separated group names")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    wanted = {g.strip() for g in args.groups.split(",") if g.strip()}

    rows = []
    workdir = tempfile.mkdtemp(prefix="codon_align_")
    try:
        for path in sorted(glob.glob(os.path.join(args.cds_dir, "*.cds.fasta"))):
            label = os.path.basename(path)[: -len(".cds.fasta")]
            if wanted and label not in wanted:
                continue
            nt_records = list(SeqIO.parse(path, "fasta"))
            if len(nt_records) < args.min_taxa:
                continue

            prot_path = os.path.join(workdir, f"{label}.faa")
            SeqIO.write(
                [
                    SeqRecord(r.seq.translate(table=11).rstrip("*"),
                              id=r.id, description="")
                    for r in nt_records
                ],
                prot_path, "fasta",
            )
            aln_path = os.path.join(workdir, f"{label}.aln.faa")
            run_mafft(prot_path, aln_path, args.threads, accurate=not args.fast)
            aln_prot = list(SeqIO.parse(aln_path, "fasta"))

            # Translation drops the terminal stop, so trim the nucleotide
            # record to match the residues being projected.
            trimmed = []
            for r in nt_records:
                aa = str(r.seq.translate(table=11))
                n_res = len(aa.rstrip("*"))
                trimmed.append(SeqRecord(r.seq[: n_res * 3], id=r.id, description=""))

            codon_aln = back_project(trimmed, aln_prot)
            problems = check_codon_integrity(codon_aln)

            out_nt = os.path.join(args.outdir, f"{label}.codon.fasta")
            SeqIO.write(codon_aln, out_nt, "fasta")
            out_aa = os.path.join(args.outdir, f"{label}.aa.fasta")
            SeqIO.write(aln_prot, out_aa, "fasta")

            rows.append(
                {
                    "group": label,
                    "n_taxa": len(codon_aln),
                    "aa_columns": len(aln_prot[0].seq),
                    "nt_columns": len(codon_aln[0].seq),
                    "frameshift_gaps": sum(1 for p in problems if "frameshift" in p[1]),
                    "other_problems": sum(1 for p in problems if "frameshift" not in p[1]),
                }
            )
            print(f"{label:<12} taxa={len(codon_aln):<3} aa={len(aln_prot[0].seq):<5} "
                  f"nt={len(codon_aln[0].seq):<6} frameshift_gaps="
                  f"{rows[-1]['frameshift_gaps']}", flush=True)
            for p in problems[:5]:
                print(f"    PROBLEM {p}")
    finally:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)

    if not rows:
        sys.exit("no alignments produced")
    man = os.path.join(args.outdir, "codon_alignment_report.tsv")
    with open(man, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    total_fs = sum(r["frameshift_gaps"] for r in rows)
    print(f"\n{len(rows)} alignments, total frameshift gaps in coding regions: {total_fs}")
    print(f"report -> {man}")


if __name__ == "__main__":
    main()
