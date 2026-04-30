#!/usr/bin/env python3
import os
import re
import csv
import gzip
import requests
import tarfile
import hashlib
from os.path import join
from itertools import islice
from datetime import datetime
from argparse import ArgumentParser


class ValidateMatrix:
    def __init__(self, url, taxa_path, base_dir, output_dir, gb_matrix, country_file, assets, host_map, country_map):
        self.url = url
        self.taxa_path = taxa_path
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.gb_matrix = gb_matrix
        self.country_file = country_file
        self.assets = assets
        self.host_map_path = host_map
        self.country_map_path = country_map

        self.dump_file = "taxadump.tar.gz"
        self.md5_url = f"{self.url}.md5"

        os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)
        os.makedirs(join(self.base_dir, self.taxa_path), exist_ok=True)

        # mapping stats
        self.host_mapped_count = 0
        self.country_mapped_count = 0

    # ---------------- Download & verify taxdump ----------------
    def get_remote_md5(self):
        print(f"Fetching MD5 checksum from {self.md5_url}")
        response = requests.get(self.md5_url, timeout=60)
        response.raise_for_status()
        return response.text.split()[0]

    def get_local_md5(self, file_path):
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def download(self):
        print(f"Downloading {self.url}")
        response = requests.get(self.url, stream=True, timeout=300)
        response.raise_for_status()
        outp = join(self.base_dir, self.taxa_path, self.dump_file)
        with open(outp, 'wb') as file:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    file.write(chunk)
        print("Download complete.")

    def verify_and_download(self):
        local_file = join(self.base_dir, self.taxa_path, self.dump_file)
        if os.path.exists(local_file):
            try:
                remote_md5 = self.get_remote_md5()
                local_md5 = self.get_local_md5(local_file)
                if local_md5 == remote_md5:
                    print("Checksum matches. Using existing taxadump file.")
                    return
                else:
                    print("Checksum mismatch. Downloading new taxadump file.")
            except Exception as e:
                print(f"Warning: Could not verify MD5 ({e}). Downloading fresh copy.")
        else:
            print("Local taxadump file not found. Downloading...")
        self.download()

    # ---------------- Extract names.dmp ----------------
    def read_tar(self):
        print(f"Extracting {self.dump_file}")
        tar_path = join(self.base_dir, self.taxa_path, self.dump_file)
        with tarfile.open(tar_path, 'r:*') as tar:
            # names.dmp → write tabbed form: tax_id \t name_txt \t unique_name \t name_class
            names_member = tar.getmember('names.dmp')
            with tar.extractfile(names_member) as fh, \
                open(join(self.base_dir, self.taxa_path, "names.dmp"), 'wb') as write_taxa:
                write_taxa.write(fh.read())
                #for line in fh:
                #    split_line = [p.strip() for p in line.decode('utf-8').split('|') if p is not None]
                #    write_taxa.write("\t".join(split_line) + '\n')

            nodes_member = tar.getmember('nodes.dmp')
            with tar.extractfile(nodes_member) as fh, \
                open(join(self.base_dir, self.taxa_path, "nodes.dmp"), "wb") as f:
                f.write(fh.read())
        print("Extraction complete.")

    # ---------------- Build taxa dict ----------------
    """
    def taxa_name_dump_to_dict(self):
        taxa_dump = {}
        scientific_name = {}
        file_name = join(self.base_dir, self.taxa_path, "names.dmp")
        with open(file_name) as f:
            for each_line in f:
                split_line = each_line.strip().split('\t')
                if 'genbank common name' in split_line or 'common name' in split_line or 'scientific name' in split_line or 'synonym' in split_line:
                    taxa_id = split_line[0]
                    name = split_line[1]
                    taxa_dump[name] = taxa_id
                    if 'scientific name' in split_line:
                        scientific_name[taxa_id] = name
        return taxa_dump, scientific_name

    """
    def taxa_name_dump_to_dict(self, *, allowed_classes=None, case_insensitive=True):
        if allowed_classes is None:
            allowed_classes = {
                "scientific name",
                "synonym",
                "common name",
                "genbank common name",
                "blast name",
                "equivalent name",
             }

        file_path = join(self.base_dir, self.taxa_path, "names.dmp")

        # autodetect gzip
        opener = gzip.open if file_path.endswith(".gz") else open
        taxa_dump = {}
        scientific = {}

        with opener(file_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # names.dmp is pipe-delimited with trailing pipe
                # e.g. '2\t|\tBacteria\t|\tBacteria <bacteria>\t|\tscientific name\t|'
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 4:
                    continue  # malformed line

                tax_id, name_txt, unique_name, name_class = parts[0], parts[1], parts[2], parts[3]

                if name_class == "scientific name":
                    # last write wins if duplicated; that’s fine for NCBI
                    scientific[tax_id] = name_txt

                if name_class in allowed_classes:
                    key = name_txt.lower() if case_insensitive else name_txt
                    # If a name maps to multiple tax_ids, last one wins.
                    # If you’d rather keep *all* matches, store a set instead.
                    taxa_dump[key] = tax_id

        return taxa_dump, scientific

    # ---------------- Generic: read CSV/TSV with header ----------------
    def _open_dict_reader(self, path):
        f = open(path, encoding="utf-8")
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=[',', '\t', ';', '|'])
            delim = dialect.delimiter
        except Exception:
            delim = '\t' if '\t' in sample and ',' not in sample else ','
        reader = csv.DictReader(f, delimiter=delim)
        return reader, delim, f

    # ---------------- Read mapping files ----------------
    def read_mapping_file(self, path, key_col, value_col):
        mapping = {}
        with open(path) as f:
            for each_line in islice(f, 1, None):
                source, replaced_by = each_line.strip().split("\t")
                mapping[source] = replaced_by
        return mapping

    # ---------------- Country dictionary ----------------
    def country_to_dict(self, infile):
        if not infile or not os.path.exists(infile):
            raise FileNotFoundError(f"M49 country file not found: {infile}")

        reader, delim, f = self._open_dict_reader(infile)
        out = {}
        with f:
            if not reader.fieldnames:
                raise ValueError(f"M49 file '{infile}' is empty or missing headers.")
            headers = [h.strip().lower() for h in reader.fieldnames]
            if 'display_name' not in headers or 'm49_code' not in headers:
                raise ValueError(f"M49 file '{infile}' must contain headers 'display_name' and 'm49_code'.")

            display_key = reader.fieldnames[headers.index('display_name')]
            m49_key = reader.fieldnames[headers.index('m49_code')]

            for row in reader:
                name = (row.get(display_key) or '').strip()
                code = (row.get(m49_key) or '').strip()
                if name and code:
                    out[name] = code
        print(f"Loaded {len(out)} M49 countries from {infile} (delimiter='{delim}')")
        return out

    # ---------------- Helpers ----------------
    @staticmethod
    def _is_na(val):
        return val is None or str(val).strip() == ""

    @staticmethod
    def _validate_date_str(date_str):
        if date_str is None:
            return "NA"
        s = str(date_str).strip()
        if not s:
            return "NA"
        for fmt in ('%d-%b-%Y', '%b-%Y', '%Y'):
            try:
                datetime.strptime(s, fmt)
                return "Yes"
            except ValueError:
                pass
        return "NA"

    def validate_country(self, country, country_dict, country_map):
        if self._is_na(country):
            return "NA", ""
        original = str(country).split(':')[0].strip()
        mapped = country_map.get(original, original)
        comment = ""
        if mapped != original:
            self.country_mapped_count += 1
            comment = f"Country mapped from mapping file: '{original}' -> '{mapped}'"
        validated = country_dict.get(mapped, "NA")
        return validated, comment

    def resolve_host(self, host_value, dict_tuple, host_map):
        comment = ""
        taxa_dict, scientific_dict = dict_tuple
        if self._is_na(host_value):
            return ("NA", "NA", "NA", comment)

        host = re.sub(r'\s*\([^)]*\)', '', str(host_value).strip())
        host_key = host.lower()

        if host_key in taxa_dict:
            taxaid = taxa_dict[host_key]
            sci_name = scientific_dict.get(taxaid, "NA")
            return ("Yes", taxaid, sci_name, comment)

        repl = host_map.get(host) or host_map.get(host_key)
        if repl:
            self.host_mapped_count += 1
            comment = f"Host mapped from mapping file: '{host}' -> '{repl}'"
            repl_key = repl.lower()
            if repl_key in taxa_dict:
                taxaid = taxa_dict[repl_key]
                sci_name = scientific_dict.get(taxaid, "NA")
                return ("Yes", taxaid, sci_name, comment)
            else:
                return ("NA", "NA", repl, comment)

        return ("NA", "NA", "NA", comment)

    # ---------------- Process & Write ----------------
    def read_meta_sheet(self, country_dict, taxa_dict, host_map, country_map):
        print(f"Reading meta data {self.gb_matrix}")

        with open(self.gb_matrix, encoding="utf-8", newline='') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            rows = [row for row in reader]
            in_headers = reader.fieldnames

        def is_processed(row):
            val = str(row.get('exclusion_status', '')).strip()
            return val != '1'

        validated_cols = ['collection_date_validated', 'country_validated',
                          'host_validated', 'host_taxa_id', 'host_scientific_name']
        out_headers = [h for h in in_headers if h not in validated_cols] + validated_cols

        total = sum(1 for r in rows if is_processed(r))
        skipped = len(rows) - total
        print(f"Processing {total} accessions, skipping {skipped} excluded")

        failed_rows = []
        for row in rows:
            row_comment = row.get('comment', '').strip()
            if is_processed(row):
                # Clean up NA
                if row_comment.upper() == "NA":
                    row_comment = ""

                # Helper to avoid duplicates
                def add_unique_comment(existing, new_info):
                    if not new_info:
                        return existing.strip(' ;|')
                    existing_norm = re.sub(r'\s+', ' ', existing.lower())
                    new_norm = re.sub(r'\s+', ' ', new_info.lower())
                    if new_norm in existing_norm:
                        return existing.strip(' ;|')
                    if existing.strip() == "":
                        return new_info.strip()
                    else:
                        return f"{existing.strip(' ;|')}; {new_info.strip()}"

                # Date
                row['collection_date_validated'] = self._validate_date_str(row.get('collection_date'))

                # Country
                country_validated, country_comment = self.validate_country(row.get('country'), country_dict, country_map)
                row['country_validated'] = country_validated
                row_comment = add_unique_comment(row_comment, country_comment)

                # Host
                host_validated, host_taxid, host_sci, host_comment = self.resolve_host(row.get('host'), taxa_dict, host_map)
                row['host_validated'] = host_validated
                row['host_taxa_id'] = host_taxid
                row['host_scientific_name'] = host_sci
                row_comment = add_unique_comment(row_comment, host_comment)

                # Cleanup
                row_comment = re.sub(r'(^NA[;| ]*|[;| ]*NA$)', '', row_comment, flags=re.IGNORECASE).strip(' ;|')
                row['comment'] = row_comment
            else:
                for c in validated_cols:
                    row[c] = ""

            if is_processed(row):
                if (row['collection_date_validated'] != "Yes") or \
                   (row['country_validated'] == "NA") or \
                   (row['host_validated'] != "Yes"):
                    failed_rows.append(row)

        outdir = join(self.base_dir, self.output_dir)
        os.makedirs(outdir, exist_ok=True)
        validated_path = join(outdir, 'gB_matrix_validated.tsv')
        with open(validated_path, 'w', encoding='utf-8', newline='') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=out_headers, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

        failed_fields = ['primary_accession', 'collection_date_validated', 'host_validated', 'country_validated',
                         'host', 'host_scientific_name', 'country', 'collection_date']
        failed_path = join(outdir, 'gB_matrix_failed_validation.tsv')
        with open(failed_path, 'w', encoding='utf-8', newline='') as f_fail:
            writer = csv.DictWriter(f_fail, fieldnames=failed_fields, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in failed_rows:
                writer.writerow(row)

        print("\n######---Validation summary---######")
        print(f"Total accessions in file: {len(rows)}")
        print(f"Processed (not excluded): {total}")
        print(f"Validated (all three passed): {int(total - len(failed_rows))}")
        print(f"Missing information (any failed): {len(failed_rows)}")
        print(f"Host names corrected via mapping file: {self.host_mapped_count}")
        print(f"Country names corrected via mapping file: {self.country_mapped_count}")
        print(f"Results saved to {outdir}\n")

        with open(self.gb_matrix, 'w', encoding='utf-8', newline='') as f_overwrite:
            writer = csv.DictWriter(f_overwrite, fieldnames=out_headers, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    def process(self):
        country_dict = self.country_to_dict(self.country_file)
        host_map = self.read_mapping_file(self.host_map_path, "host", "replaced_by")
        country_map = self.read_mapping_file(self.country_map_path, "country", "replaced_by")
        self.verify_and_download()
        self.read_tar()
        taxa_dict = self.taxa_name_dump_to_dict()
        self.read_meta_sheet(country_dict, taxa_dict, host_map, country_map)


def main():
    parser = ArgumentParser(description='Validate the gB_matrix file based on country, date, and host columns.')
    parser.add_argument('-u', '--url', default="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
    parser.add_argument('-t', '--taxa_path', default='Taxa')
    parser.add_argument('-b', '--base_dir', default='tmp')
    parser.add_argument('-o', '--output_dir', default='Validate-matrix')
    parser.add_argument('-g', '--gb_matrix', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
    parser.add_argument('-c', '--country', default='assets/m49_country.csv')
    parser.add_argument('-a', '--assets', default='assets/')
    parser.add_argument('-m', '--host_map', default="generic/rabv/host_mapping.tsv")
    parser.add_argument('-n', '--country_map', default="generic/rabv/country_mapping.tsv")
    args = parser.parse_args()

    validator = ValidateMatrix(
        args.url, args.taxa_path, args.base_dir, args.output_dir,
        args.gb_matrix, args.country, args.assets, args.host_map, args.country_map
    )
    validator.process()


if __name__ == "__main__":
    main()

