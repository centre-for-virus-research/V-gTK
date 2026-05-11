import argparse
import subprocess
import os

class GenotypePipeline:
    def __init__(self, reference_fasta, major_clade, minor_clade, tree, ref_accession):
        self.reference_fasta = reference_fasta
        self.major_clade = major_clade
        self.minor_clade = minor_clade
        self.tree = tree
        self.ref_accession = ref_accession

        # Derived filenames
        self.base_prefix = os.path.splitext(os.path.basename(self.reference_fasta))[0]
        self.vcf_output = f"{self.base_prefix}.vcf"
        self.pb_backbone = f"{self.base_prefix}_backbone.pb"
        self.pb_genotype = f"{self.base_prefix}_genotype.pb"
        self.pb_subgenotype = f"{self.base_prefix}_subgenotype.pb"

    def run_command(self, command, description):
        print(f"\n[INFO] {description}")
        print(f"[CMD] {command}")
        result = subprocess.run(command, shell=True)
        if result.returncode != 0:
            raise RuntimeError(f"[ERROR] Command failed: {command}")
        print("[INFO] Step completed.")

    def convert_fasta_to_vcf(self):
        cmd = f"faToVcf -ref={self.ref_accession} {self.reference_fasta} {self.vcf_output}"
        self.run_command(cmd, "Converting aligned FASTA to VCF")

    def build_usher_tree(self):
        cmd = f"usher --vcf {self.vcf_output} --tree {self.tree} -o {self.pb_backbone}"
        self.run_command(cmd, "Building UShER protobuf tree")

    def annotate_major_clades(self):
        cmd = f"matUtils annotate -i {self.pb_backbone} -c {self.major_clade} -o {self.pb_genotype}"
        self.run_command(cmd, "Annotating major clades")

    def annotate_minor_clades(self):
        cmd = f"matUtils annotate -i {self.pb_genotype} -c {self.minor_clade} -o {self.pb_subgenotype}"
        self.run_command(cmd, "Annotating minor clades")

    def run(self):
        self.convert_fasta_to_vcf()
        self.build_usher_tree()
        self.annotate_major_clades()
        self.annotate_minor_clades()
        print("\n[INFO] All steps completed successfully.")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Genotyping using UShER")
    parser.add_argument("--reference_fasta", required=True, help="Aligned reference FASTA file")
    parser.add_argument("--major_clade", required=True, help="Major clade TSV file")
    parser.add_argument("--minor_clade", required=True, help="Minor clade TSV file")
    parser.add_argument("--tree", required=True, help="Reference Newick tree file")
    parser.add_argument("--ref_accession", required=True, help="Reference accession name (e.g., NC_001542)")
    return parser.parse_args()

def main():
    args = parse_arguments()
    pipeline = GenotypePipeline(
        reference_fasta=args.reference_fasta,
        major_clade=args.major_clade,
        minor_clade=args.minor_clade,
        tree=args.tree,
        ref_accession=args.ref_accession
    )
    pipeline.run()

if __name__ == "__main__":
    main()

