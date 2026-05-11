import os
import argparse
import subprocess
import pandas as pd
from os.path import join

class NCBI_GFF_Downloader:
    def __init__(self, accession_ids, base_dir, output_dir):
        self.accession_ids = accession_ids
        self.base_dir = base_dir
        self.output_dir = output_dir

    def get_accession_list(self):
        """Parse accession_ids argument. Can be a comma-separated string or a file path."""
        if os.path.isfile(self.accession_ids):
            try:
                # Try reading as TSV without header
                df = pd.read_csv(self.accession_ids, sep='\t', header=None, dtype=str)
                
                # If 2nd column exists and contains 'master', filter by it
                if df.shape[1] >= 2:
                    # Check if any value in column 1 is 'master' (case insensitive)
                    if df[1].str.lower().eq('master').any():
                        print(f"Found 'master' indicators in {self.accession_ids}. Filtering...")
                        masters = df[df[1].str.lower() == 'master']
                        return masters[0].tolist()
                
                # If no 'master' found or only 1 column, return all from first column
                print(f"No 'master' indicators found in {self.accession_ids}. Using all accessions from first column.")
                return df[0].tolist()
                
            except Exception as e:
                print(f"Error reading accession file: {e}")
                return []
        else:
            # Assume comma-separated string
            return [x.strip() for x in self.accession_ids.split(',') if x.strip()]

    def download_gff(self):
        os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)
        
        accessions = self.get_accession_list()
        if not accessions:
            print("No accessions found to download.")
            return

        print(f"Downloading GFFs for {len(accessions)} accessions...")

        for each_accession in accessions:
            try:
                output_file = join(self.base_dir, self.output_dir, each_accession + ".gff3")
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                     print(f"File {output_file} already exists. Skipping.")
                     continue

                command = ["efetch", "-db", "nuccore", "-id", each_accession, "-format", "gff3"]
                
                with open(output_file, "w") as output:
                    subprocess.run(command, check=True, stdout=output, stderr=subprocess.PIPE)
                    print(f"Successfully downloaded GFF3 file: {each_accession}")
            except subprocess.CalledProcessError as e:
                print(f"Error downloading GFF3 file for {each_accession}: {e.stderr.decode().strip()}")
            except Exception as e:
                print(f"An error occurred for {each_accession}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download GFF3 file from NCBI using efetch")
    parser.add_argument("-id", "--accession_ids", required=True, help="NCBI accession ID (comma separated) or path to reference list file (TSV)")
    parser.add_argument("-b", "--base_dir", help="Base directory", default="tmp")
    parser.add_argument("-o", "--output_dir", help="Directory where the GFF files are saved", default="Gff")
        
    args = parser.parse_args()
    downloader = NCBI_GFF_Downloader(args.accession_ids, args.base_dir, args.output_dir)
    downloader.download_gff()
