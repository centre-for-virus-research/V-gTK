import os
import argparse
import subprocess
from os.path import join
from ExportRefListFromUpdateDb import load_master_accessions, load_master_accessions_from_file, load_reference_file_accessions

class NCBI_GFF_Downloader:
    def __init__(self, accession_ids, base_dir, output_dir, update_db=None):
        self.accession_ids = accession_ids
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.update_db = update_db

    def get_accession_list(self):
        """Parse accession_ids argument. Can be a comma-separated string or a file path."""
        if self.update_db:
            masters = load_master_accessions(self.update_db)
            if masters:
                return masters
        if os.path.isfile(self.accession_ids):
            try:
                masters = load_master_accessions_from_file(self.accession_ids)
                if masters:
                    return masters
                return load_reference_file_accessions(self.accession_ids)
            except Exception as e:
                print(f"Error reading accession file: {e}")
                return []
        else:
            # Assume comma-separated string
            return [x.strip() for x in self.accession_ids.split(',') if x.strip()]

    @staticmethod
    def looks_like_gff3(path):
        """True if the file's first non-empty line is a GFF3 version pragma.

        efetch writes its own error text to stdout and can still exit 0, so a
        failed fetch leaves a non-empty file that the size-only skip check
        below then trusts forever. generic/influenza/ref_list.txt.gff3 was a
        committed error log for that reason.
        """
        try:
            with open(path, encoding='utf-8', errors='replace') as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    return line.startswith('##gff-version')
        except OSError:
            return False
        return False

    @staticmethod
    def validate_accession(accession):
        """Reject anything that is plainly a path rather than an accession.

        A path reaching efetch as -id is split on '/' and '-' and every
        fragment is reported as an unknown id, so the failure is only visible
        in the output file. Caught here instead.
        """
        if os.sep in accession or '/' in accession or any(c.isspace() for c in accession):
            raise ValueError(
                f"'{accession}' is not an NCBI accession. A file path was most likely "
                f"passed to --accession_ids where the file could not be read."
            )

    def download_gff(self):
        os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)

        accessions = self.get_accession_list()
        if not accessions:
            print("No accessions found to download.")
            return

        print(f"Downloading GFFs for {len(accessions)} accessions...")

        for each_accession in accessions:
            # Validated before an output path is derived from it, so a bad
            # accession can never name a file that _discard then deletes.
            try:
                self.validate_accession(each_accession)
            except ValueError as e:
                print(f"An error occurred for {each_accession}: {e}")
                continue

            output_file = join(self.base_dir, self.output_dir, each_accession + ".gff3")
            try:
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    if self.looks_like_gff3(output_file):
                        print(f"File {output_file} already exists. Skipping.")
                        continue
                    print(f"File {output_file} exists but is not GFF3; re-downloading.")

                command = ["efetch", "-db", "nuccore", "-id", each_accession, "-format", "gff3"]

                with open(output_file, "w") as output:
                    subprocess.run(command, check=True, stdout=output, stderr=subprocess.PIPE)

                if not self.looks_like_gff3(output_file):
                    raise ValueError(
                        f"efetch returned no GFF3 for '{each_accession}'. First line: "
                        f"{self._first_line(output_file)!r}"
                    )
                print(f"Successfully downloaded GFF3 file: {each_accession}")
            except subprocess.CalledProcessError as e:
                self._discard(output_file)
                print(f"Error downloading GFF3 file for {each_accession}: {e.stderr.decode().strip()}")
            except Exception as e:
                self._discard(output_file)
                print(f"An error occurred for {each_accession}: {e}")

    @staticmethod
    def _first_line(path):
        try:
            with open(path, encoding='utf-8', errors='replace') as handle:
                return handle.readline().strip()
        except OSError:
            return ''

    @staticmethod
    def _discard(path):
        """Remove a failed download so it is not mistaken for a real GFF3."""
        try:
            if os.path.exists(path):
                os.remove(path)
                print(f"Removed incomplete download: {path}")
        except OSError as e:
            print(f"Could not remove incomplete download {path}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download GFF3 file from NCBI using efetch")
    parser.add_argument("-id", "--accession_ids", required=True, help="NCBI accession ID (comma separated) or path to reference list file (TSV)")
    parser.add_argument("-b", "--base_dir", help="Base directory", default="tmp")
    parser.add_argument("-o", "--output_dir", help="Directory where the GFF files are saved", default="Gff")
    parser.add_argument("--update_db", help="Existing update DB; when set, derive master accession IDs from DB", default=None)
        
    args = parser.parse_args()
    downloader = NCBI_GFF_Downloader(args.accession_ids, args.base_dir, args.output_dir, args.update_db)
    downloader.download_gff()
