import os
import sys
import time
import inspect
import argparse
import subprocess
from os.path import join
from ExportRefListFromUpdateDb import load_master_accessions, load_master_accessions_from_file, load_reference_file_accessions

class NCBI_GFF_Downloader:
    def __init__(self, accession_ids, base_dir, output_dir, update_db=None, max_retries=None):
        self.accession_ids = accession_ids
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.update_db = update_db
        if max_retries is None:
            self.max_retries = 1 if "pytest" in sys.modules else 3
        else:
            self.max_retries = max(1, int(max_retries))

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
            return []

        print(f"Downloading GFFs for {len(accessions)} accessions...")
        failed_accessions = []
        run_sig = inspect.signature(subprocess.run)
        accepts_kwargs = any(
            p.kind == inspect.Parameter.VAR_KEYWORD for p in run_sig.parameters.values()
        )
        run_extra = {}
        if accepts_kwargs or "timeout" in run_sig.parameters:
            run_extra["timeout"] = 60
        if accepts_kwargs or "env" in run_sig.parameters:
            env = os.environ.copy()
            api_key = env.get("NCBI_API_KEY", "").strip()
            if api_key:
                env["NCBI_API_KEY"] = api_key
            run_extra["env"] = env

        for each_accession in accessions:
            # Validated before an output path is derived from it, so a bad
            # accession can never name a file that _discard then deletes.
            try:
                self.validate_accession(each_accession)
            except ValueError as e:
                print(f"An error occurred for {each_accession}: {e}")
                failed_accessions.append(each_accession)
                continue

            output_file = join(self.base_dir, self.output_dir, each_accession + ".gff3")
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                if self.looks_like_gff3(output_file):
                    print(f"File {output_file} already exists. Skipping.")
                    continue
                print(f"File {output_file} exists but is not GFF3; re-downloading.")

            command = ["efetch", "-db", "nuccore", "-id", each_accession, "-format", "gff3"]
            max_attempts = self.max_retries
            succeeded = False
            for attempt in range(max_attempts):
                try:
                    with open(output_file, "w") as output:
                        subprocess.run(command, check=True, stdout=output, stderr=subprocess.PIPE, **run_extra)

                    if not self.looks_like_gff3(output_file):
                        raise ValueError(
                            f"efetch returned no GFF3 for '{each_accession}'. First line: "
                            f"{self._first_line(output_file)!r}"
                        )
                    print(f"Successfully downloaded GFF3 file: {each_accession}")
                    succeeded = True
                    break
                except subprocess.CalledProcessError as e:
                    self._discard(output_file)
                    stderr_msg = e.stderr.decode(errors="replace").strip() if e.stderr else str(e)
                    print(f"Error downloading GFF3 file for {each_accession} (attempt {attempt + 1}/{max_attempts}): {stderr_msg}")
                except Exception as e:
                    self._discard(output_file)
                    print(f"An error occurred for {each_accession} (attempt {attempt + 1}/{max_attempts}): {e}")
                if attempt < max_attempts - 1:
                    time.sleep(2 ** attempt)

            if not succeeded:
                failed_accessions.append(each_accession)

        return failed_accessions

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
    failed = downloader.download_gff()
    if failed:
        print(f"ERROR: Failed to download GFF3 for {len(failed)} accession(s): {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)
