#!/usr/bin/env python3
"""
Normla mode:
    python scripts/CladeAssignment.py   --ref-aln generic/rabv/tree/ref_plus_am3ca_am5.fa   --ref-tree generic/rabv/tree/ref_tree_am3c_am5.treefile   --query tmp/Pad-alignment/alUnc509RefseqsMafftHandModified.fa   --taxon-major generic/rabv/ref_major_clades.tsv   --taxon-minor generic/rabv/ref_minor_clades.tsv   --threads 6   --meta-data tmp/GenBank-matrix/gB_matrix_raw.tsv   --steps all   --skip-mafft

Update mode:
    python scripts/CladeAssignment.py   --ref-aln generic/rabv/tree/ref_plus_am3ca_am5.fa   --ref-tree generic/rabv/tree/ref_tree_am3c_am5.treefile   --query tmp/Update/Pad-alignment/alUnc509RefseqsMafftHandModified.fa   --aligned-out tmp/Update/Pad-alignment/alUnc509RefseqsMafftHandModified.fa   --taxon-major generic/rabv/ref_major_clades.tsv   --taxon-minor generic/rabv/ref_minor_clades.tsv   --threads 6   --meta-data tmp/Update/GenBank-matrix/gB_matrix_raw.tsv   --steps all   --skip-mafft --update
"""

import os
import sys
import shutil
import argparse
import subprocess
import pandas as pd
from os.path import join


class CladeAssignment:
    def __init__(
        self,
        ref_aln,
        ref_tree,
        query_fa,
        taxon_major,
        taxon_minor,
        base_dir="tmp",
        output_dir="CladeAssignment",
        steps="all",
        threads=6,
        dry_run=False,
        # executables
        mafft_exe="mafft",
        epa_exe="epa-ng",
        gappa_exe="gappa",
        # mafft opts
        mafft_anysymbol=True,
        mafft_extra="",
        #aligned_out="",
        # epa opts
        epa_model="TVM+F+R5",
        epa_redo=True,
        epa_workdir="",
        epa_extra="",
        jplace="",
        # gappa opts
        gappa_major_outdir="",
        gappa_minor_outdir="",
        gappa_extra="",
        # metadata / matrix
        meta_data="tmp/Genbank-matrix/gB_matrix_raw.tsv",
        meta_id_col="primary_accession",
        out_matrix="",
        strip_version=True,
        skip_mafft=False,
        
    ):
        if pd is None:
            self._die("pandas is required for the matrix update step. Install: pip install pandas")

        # --------------------
        # inputs (store first)
        # --------------------
        self.ref_aln = ref_aln
        self.ref_tree = ref_tree
        self.query_fa = query_fa
        self.taxon_major = taxon_major
        self.taxon_minor = taxon_minor

        # flags
        self.skip_mafft = skip_mafft
        # dirs
        self.base_dir = base_dir
        self.output_dir = output_dir

        # tools
        self.mafft_exe = mafft_exe
        self.epa_exe = epa_exe
        self.gappa_exe = gappa_exe

        # options
        self.steps = steps
        self.threads = threads
        self.dry_run = dry_run

        self.mafft_anysymbol = mafft_anysymbol
        self.mafft_extra = mafft_extra

        self.epa_model = epa_model
        self.epa_redo = epa_redo
        self.epa_extra = epa_extra

        self.gappa_extra = gappa_extra

        # metadata
        self.meta_data = meta_data
        self.meta_id_col = meta_id_col
        self.strip_version = strip_version

        # --------------------
        # FIX #1: make INPUT paths absolute
        # This prevents failures when we run tools with cwd=self.epa_workdir
        # --------------------
        self.ref_aln = os.path.abspath(self.ref_aln)
        self.ref_tree = os.path.abspath(self.ref_tree)
        self.query_fa = os.path.abspath(self.query_fa)
        self.taxon_major = os.path.abspath(self.taxon_major)
        self.taxon_minor = os.path.abspath(self.taxon_minor)
        self.meta_data = os.path.abspath(self.meta_data)

        # --------------------
        # output root
        # normal: <base_dir>/<output_dir>
        # update : <base_dir>/update/<output_dir>
        # --------------------
        #outdir_rel = join(self.base_dir, "update", self.output_dir) if self.update else join(self.base_dir, self.output_dir)
        #self.outdir = os.path.abspath(outdir_rel)

        outdir_rel = join(self.base_dir, self.output_dir)
        self.outdir = os.path.abspath(outdir_rel)

        # derived outputs (inside outdir by default)
        #self.aligned_out = os.path.abspath(aligned_out.strip() or join(self.outdir, "input_seqs_with_ref_alignment.fa"))
        if self.skip_mafft:
            self.aligned_out = os.path.abspath(self.query_fa)
        else:
            self.aligned_out = os.path.abspath(join(self.outdir, "input_seqs_with_ref_alignment.fa"))

        self.epa_workdir = os.path.abspath(epa_workdir.strip() or join(self.outdir, "epa-ng"))
        self.jplace = os.path.abspath(jplace.strip() or join(self.epa_workdir, "epa_result.jplace"))

        self.gappa_major_outdir = os.path.abspath(gappa_major_outdir.strip() or join(self.outdir, "gappa_major_clades_assigned"))
        self.gappa_minor_outdir = os.path.abspath(gappa_minor_outdir.strip() or join(self.outdir, "gappa_minor_clades_assigned"))

        self.major_per_query = os.path.abspath(join(self.gappa_major_outdir, "per_query.tsv"))
        self.minor_per_query = os.path.abspath(join(self.gappa_minor_outdir, "per_query.tsv"))


        if out_matrix.strip():
            self.out_matrix = os.path.abspath(out_matrix.strip())
        else:
            self.out_matrix = os.path.abspath(self.meta_data)

    def _die(self, msg):
        print("[error]", msg, file=sys.stderr)
        sys.exit(1)

    def _check_file(self, path, label):
        if not os.path.exists(path):
            self._die(label + " not found: " + path)
        if os.path.isdir(path):
            self._die(label + " is a directory (expected a file): " + path)

    def _ensure_dir(self, path):
        os.makedirs(path, exist_ok=True)

    def _check_exe(self, name):
        if shutil.which(name) is None:
            self._die("Executable not found on PATH: " + name)

    def _run(self, cmd, cwd=None):
        print("[cmd]", " ".join(cmd), file=sys.stderr)
        if cwd:
            print("[cwd]", cwd, file=sys.stderr)
        if self.dry_run:
            return
        subprocess.run(cmd, cwd=cwd, check=True)

    def _norm_id(self, s):
        if s is None:
            return ""
        s = str(s).strip()
        if self.strip_version and "." in s:
            s = s.split(".", 1)[0]
        return s

    def validate_inputs(self):
        self._check_file(self.ref_aln, "Reference alignment (--ref-aln)")
        self._check_file(self.ref_tree, "Reference tree (--ref-tree)")
        self._check_file(self.query_fa, "Query FASTA (--query)")
        self._check_file(self.taxon_major, "Major clade taxon file (--taxon-major)")
        self._check_file(self.taxon_minor, "Minor clade taxon file (--taxon-minor)")
        self._check_file(self.meta_data, "Meta/matrix file (--meta-data)")
        self._ensure_dir(self.outdir)

    # ----------------------
    # Steps 1..4
    # ----------------------
    def run_mafft(self):
        self._check_exe(self.mafft_exe)
        self._ensure_dir(self.outdir)

        cmd = [self.mafft_exe, "--add", self.query_fa, "--keeplength"]
        if self.mafft_anysymbol:
            cmd.append("--anysymbol")
        cmd += ["--thread", str(self.threads), self.ref_aln]

        if self.mafft_extra.strip():
            cmd += self.mafft_extra.strip().split()

        print("[mafft] writing:", self.aligned_out, file=sys.stderr)

        if self.dry_run:
            print("[cmd]", " ".join(cmd), ">", self.aligned_out, file=sys.stderr)
            return

        with open(self.aligned_out, "w") as fh:
            subprocess.run(cmd, stdout=fh, check=True)

    def run_epa_ng(self):
        self._check_exe(self.epa_exe)
        self._ensure_dir(self.epa_workdir)
        self._check_file(self.aligned_out, "MAFFT output (aligned query+ref)")

        cmd = [self.epa_exe]
        if self.epa_redo:
            cmd.append("--redo")

        cmd += [
            "-m",
            self.epa_model,
            "-t",
            self.ref_tree,   # now absolute
            "-s",
            self.ref_aln,    # now absolute
            "-q",
            self.aligned_out,
            "-T",
            str(self.threads),
        ]

        if self.epa_extra.strip():
            cmd += self.epa_extra.strip().split()

        self._run(cmd, cwd=self.epa_workdir)

        if (not self.dry_run) and (not os.path.exists(self.jplace)):
            self._die("EPA-ng finished but jplace not found where expected: " + self.jplace)

    def run_gappa_major(self):
        self._check_exe(self.gappa_exe)
        self._ensure_dir(self.gappa_major_outdir)

        cmd = [
            self.gappa_exe,
            "examine",
            "assign",
            "--jplace-path",
            self.jplace,
            "--taxon-file",
            self.taxon_major,  # absolute
            "--out-dir",
            self.gappa_major_outdir,
            "--per-query-results",
        ]
        if self.gappa_extra.strip():
            cmd += self.gappa_extra.strip().split()

        self._run(cmd)

        if (not self.dry_run) and (not os.path.exists(self.major_per_query)):
            self._die("Expected gappa per_query.tsv not found: " + self.major_per_query)

    def run_gappa_minor(self):
        self._check_exe(self.gappa_exe)
        self._ensure_dir(self.gappa_minor_outdir)

        cmd = [
            self.gappa_exe,
            "examine",
            "assign",
            "--jplace-path",
            self.jplace,
            "--taxon-file",
            self.taxon_minor,  # absolute
            "--out-dir",
            self.gappa_minor_outdir,
            "--per-query-results",
        ]
        if self.gappa_extra.strip():
            cmd += self.gappa_extra.strip().split()

        self._run(cmd)

        if (not self.dry_run) and (not os.path.exists(self.minor_per_query)):
            self._die("Expected gappa per_query.tsv not found: " + self.minor_per_query)

    # ----------------------
    # Step 5: Parse per_query.tsv
    # ----------------------
    def parse_per_query(self, per_query_path):
        df = pd.read_csv(per_query_path, sep="\t", dtype=str).fillna("")
        needed = ["name", "LWR", "taxopath"]
        for c in needed:
            if c not in df.columns:
                self._die(f"{per_query_path} missing required column {c!r}. Found: {list(df.columns)}")

        def to_float(x):
            try:
                return float(str(x))
            except Exception:
                return -1.0

        out = {}
        df["__nid__"] = df["name"].map(self._norm_id)

        for nid, g in df.groupby("__nid__", sort=False):
            rows = []
            for _, r in g.iterrows():
                clade = str(r.get("taxopath", "")).strip()
                score_str = str(r.get("LWR", "")).strip()
                score_f = to_float(score_str)
                rows.append((score_f, score_str, clade))

            rows.sort(key=lambda x: x[0], reverse=True)

            best_clade = rows[0][2] if rows else ""
            all_clades = ";".join([x[2] if x[2] != "" else "NULL" for x in rows])
            all_scores = ";".join([x[1] for x in rows])

            out[nid] = {
                "best_clade": best_clade if best_clade != "" else "",
                "all_clades": all_clades,
                "all_scores": all_scores,
            }
        return out

    # ----------------------
    # Step 6: Update gB_matrix
    # ----------------------
    def update_matrix(self):
        major_map = self.parse_per_query(self.major_per_query)
        minor_map = self.parse_per_query(self.minor_per_query)

        mat = pd.read_csv(self.meta_data, sep="\t", dtype=str).fillna("")
        if self.meta_id_col not in mat.columns:
            self._die(
                f"meta-data file missing id column {self.meta_id_col!r}. "
                f"Columns: {list(mat.columns)}"
            )

        cols = [
            "EPA_major_clade",
            "EPA_minor_clade",
            "EPA_major_all",
            "EPA_minor_all",
            "major_LWR_score",
            "minor_LWR_score",
        ]
        for c in cols:
            if c not in mat.columns:
                mat[c] = ""

        ids = mat[self.meta_id_col].map(self._norm_id)

        def get_major(nid, key, default=""):
            return major_map.get(nid, {}).get(key, default)

        def get_minor(nid, key, default=""):
            return minor_map.get(nid, {}).get(key, default)

        mat["EPA_major_clade"] = [get_major(nid, "best_clade", "") for nid in ids]
        mat["EPA_minor_clade"] = [get_minor(nid, "best_clade", "") for nid in ids]

        mat["EPA_major_all"] = [get_major(nid, "all_clades", "") for nid in ids]
        mat["EPA_minor_all"] = [get_minor(nid, "all_clades", "") for nid in ids]
        mat["major_LWR_score"] = [get_major(nid, "all_scores", "") for nid in ids]
        mat["minor_LWR_score"] = [get_minor(nid, "all_scores", "") for nid in ids]

        self._ensure_dir(os.path.dirname(self.out_matrix))
        mat.to_csv(self.out_matrix, sep="\t", index=False)
        print("[matrix] written:", self.out_matrix, file=sys.stderr)

    # ----------------------
    # Orchestrate
    # ----------------------
    def run_all(self):
        print(f"[outdir] {self.outdir}", file=sys.stderr)

        '''
        if self.steps in ("all", "mafft"):
            self.run_mafft()
        '''
        if self.steps in ("all", "mafft"):
            if self.steps == "mafft":
                self.run_mafft()
            else:
                if self.skip_mafft:
                    print("[mafft] skipped (--skip-mafft enabled)", file=sys.stderr)
                else:
                    self.run_mafft()

        if self.steps in ("all", "epa"):
            if (not self.dry_run) and (not os.path.exists(self.query_fa)):
                self._die("Missing MAFFT output. Run --steps mafft first, or use --steps all.")
            self.run_epa_ng()

        if self.steps in ("all", "gappa-major"):
            if (not self.dry_run) and (not os.path.exists(self.jplace)):
                self._die("Missing jplace. Run --steps epa first, or use --steps all.")
            self.run_gappa_major()

        if self.steps in ("all", "gappa-minor"):
            if (not self.dry_run) and (not os.path.exists(self.jplace)):
                self._die("Missing jplace. Run --steps epa first, or use --steps all.")
            self.run_gappa_minor()

        if self.steps in ("all", "update-matrix"):
            if not self.dry_run:
                self._check_file(self.major_per_query, "Major per_query.tsv")
                self._check_file(self.minor_per_query, "Minor per_query.tsv")
            self.update_matrix()
        
        print("[done] ok", file=sys.stderr)
        print("[outdir]", self.outdir, file=sys.stderr)


def build_arg_parser():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="EPA-ng placement + GAPPA assignment + update gB_matrix",
    )

    p.add_argument("--ref-aln", required=True, help="Reference alignment FASTA")
    p.add_argument("--ref-tree", required=True, help="Reference tree (.treefile etc.)")
    p.add_argument("--query", required=True, help="Query sequences FASTA")
    p.add_argument("--taxon-major", required=True, help="GAPPA taxon TSV for major clades")
    p.add_argument("--taxon-minor", required=True, help="GAPPA taxon TSV for minor clades")

    p.add_argument("--base-dir", default="tmp", help="Base directory")
    p.add_argument("--output-dir", default="CladeAssignment", help="Output directory name under base dir")

    p.add_argument(
        "--steps",
        default="all",
        choices=["all", "mafft", "epa", "gappa-major", "gappa-minor", "update-matrix"],
        help="Which step to run",
    )
    p.add_argument("--dry-run", action="store_true", help="Print commands only")
    p.add_argument("--threads", type=int, default=6, help="Threads for MAFFT and EPA-ng")

    # tools
    p.add_argument("--mafft-exe", default="mafft", help="MAFFT executable")
    p.add_argument("--epa-exe", default="epa-ng", help="EPA-ng executable")
    p.add_argument("--gappa-exe", default="gappa", help="GAPPA executable")

    # MAFFT options
    p.add_argument("--mafft-anysymbol", action="store_true", default=True, help="Use --anysymbol")
    p.add_argument("--mafft-extra", default="", help="Extra MAFFT args (space-separated)")
    #p.add_argument("--aligned-out", default="", help="Override MAFFT output FASTA path")

    # EPA-ng options
    p.add_argument("--epa-model", default="TVM+F+R5", help="EPA-ng model (-m)")
    p.add_argument("--epa-redo", action="store_true", default=True, help="Use --redo")
    p.add_argument("--epa-workdir", default="", help="Override EPA-ng working directory")
    p.add_argument("--epa-extra", default="", help="Extra EPA-ng args (space-separated)")
    p.add_argument("--jplace", default="", help="Override expected epa_result.jplace path")

    # GAPPA options
    p.add_argument("--gappa-major-outdir", default="", help="Override major output dir")
    p.add_argument("--gappa-minor-outdir", default="", help="Override minor output dir")
    p.add_argument("--gappa-extra", default="", help="Extra GAPPA args (space-separated)")

    # metadata / matrix options
    p.add_argument(
        "--meta-data",
        default="tmp/Genbank-matrix/gB_matrix_raw.tsv",
        help="Input gB_matrix file (TSV)",
    )
    p.add_argument(
        "--meta-id-col",
        default="primary_accession",
        help="Column to match per_query names against",
    )
    p.add_argument(
        "--out-matrix",
        default="",
        help="Output matrix path. If empty: normal mode overwrites --meta-data; update mode writes to update outdir.",
    )
    p.add_argument(
        "--no-strip-version",
        action="store_true",
        help="Do NOT strip .1/.2 version when matching accessions",
    )

    p.add_argument(
        "--skip-mafft",
        action="store_true",
        help="Skip MAFFT step even when --steps is 'all' (expects --aligned-out to already exist, or default aligned_out path to exist).",
    )
    
    return p


def main():
    args = build_arg_parser().parse_args()

    runner = CladeAssignment(
        ref_aln=args.ref_aln,
        ref_tree=args.ref_tree,
        query_fa=args.query,
        taxon_major=args.taxon_major,
        taxon_minor=args.taxon_minor,
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        steps=args.steps,
        threads=args.threads,
        dry_run=args.dry_run,
        mafft_exe=args.mafft_exe,
        epa_exe=args.epa_exe,
        gappa_exe=args.gappa_exe,
        mafft_anysymbol=args.mafft_anysymbol,
        mafft_extra=args.mafft_extra,
        #aligned_out=args.aligned_out,
        epa_model=args.epa_model,
        epa_redo=args.epa_redo,
        epa_workdir=args.epa_workdir,
        epa_extra=args.epa_extra,
        jplace=args.jplace,
        gappa_major_outdir=args.gappa_major_outdir,
        gappa_minor_outdir=args.gappa_minor_outdir,
        gappa_extra=args.gappa_extra,
        meta_data=args.meta_data,
        meta_id_col=args.meta_id_col,
        out_matrix=args.out_matrix,
        strip_version=(not args.no_strip_version),
        skip_mafft=args.skip_mafft,
    )
    runner.validate_inputs()

    try:
        runner.run_all()
    except subprocess.CalledProcessError as e:
        runner._die("Command failed with exit code " + str(e.returncode))

if __name__ == "__main__":
    main()