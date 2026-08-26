import os
import re
import csv
import glob
import shutil
import sqlite3
import numpy as np
import pandas as pd
from datetime import datetime
from Bio import SeqIO
from argparse import ArgumentParser

""" TO-DO: 
parse country and geo location from location column
clean Subtype to match "serotype" pattern in genbank (now is done during the  serotype clean step downstream)
"""



def _restore_blank_as_missing(df):
    """Make blanks missing again after a lossless read.

    Read with ``keep_default_na=False`` nothing is nulled, which is what lets
    influenza's ``NA`` segment survive - pandas' default sentinel list contains
    the literal string "NA", and neuraminidase is genuinely named NA. But the
    rest of this module legitimately relies on blank meaning missing
    (``dropna(how="all")``, ``pd.notnull(row["Segment_Id"])``), so the empty
    string - and only the empty string - is converted back to NaN.

    Net effect versus the old behaviour: '' still nulls; 'NA', 'NULL', 'None'
    and 'nan' no longer do.
    """
    return df.replace('', np.nan)


SEGS = ["PB1", "PB2", "PA", "HA", "NA", "NP", "NS", "MP"]
ACCEPTABLE_NUCS = set('ACGTNRYKMWSBDHVXZ')

def extract_segment_epi(header, epi_ids_set, epi_dict, segment_names=None, log_handle=None):
    """Extract Segment_Id (EPI#) from a fasta header"""
    if segment_names is None:
        segment_names = set(SEGS)

    clean_header = re.sub(r'\s*\|\s*', '|', header.replace('>', '').strip())
    fields = clean_header.split('|')

    # 1. Look for EPI#######
    for field in fields:
        m = re.fullmatch(r"(EPI\d+)", field.strip())
        if m:
            epi = m.group(1)
            if epi in epi_ids_set:
                return epi

    # 2. Look for a numeric field; if in metadata, return as 'EPI'+number
    for field in fields:
        n = field.strip()
        if n.isdigit() and len(n) >= 3:  # check length threshold
            epi_candidate = f"EPI{n}"
            if epi_candidate in epi_ids_set:
                return epi_candidate

    # 3. Look for EPI_ISL_#### and segment, cross-check metadata dict
    isolate_id = None
    segment = None
    for field in fields:
        m = re.fullmatch(r"EPI_ISL_\d+", field.strip())
        if m:
            isolate_id = m.group(0)
    for field in fields:
        if field.strip() in segment_names:
            segment = field.strip()
    if isolate_id and segment:
        epi = epi_dict.get((isolate_id, segment))
        if epi:
            return epi

    # 4. No match, log failed header 
    if log_handle:
        if any(re.fullmatch(r"EPI\d+", field.strip()) for field in fields):
            error_type = "Segment_ID (EPI) ambiguous"
        else:
            error_type = "Segment_ID (EPI) not_found"
        log_handle.write(f"{header}\t{error_type}\n")

    return None


class GISAIDTidy:
    def __init__(self, data_dir, output_dir, filetype="xls", log_file=None, db_file=None, test_mode=False):
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.filetype = filetype
        self.log_file = log_file or os.path.abspath("gisaid_processing.log")
        self.db_file = db_file
        self.test_mode = test_mode
        if self.test_mode:
            self.log_file = log_file or os.path.abspath("test_gisaid_processing.log")
        os.makedirs(self.output_dir, exist_ok=True)
        self.seen_files = self.load_seen_files()
        self.processed_files = []
        self.log(f"Pipeline started at {self.timestamp()}. Filetype: {self.filetype}. Test mode: {self.test_mode}")

    def timestamp(self):
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def log(self, message):
        with open(self.log_file, "a") as logf:
            logf.write(f"[{self.timestamp()}] {message}\n")
        print(f"[{self.timestamp()}] {message}")

    def load_seen_files(self):
        """Load list of previously processed files from log."""
        if os.path.exists(self.log_file):
            with open(self.log_file) as logf:
                seen = [line.strip().split(": ")[-1] for line in logf if line.startswith("[") and "Processed file:" in line]
            return set(seen)
        return set()

    def update_log_with_file(self, filename):
        self.log(f"Processed file: {filename}")

    def clean_metadata(self, metadata_files, db_file=None):
        """Parse and tidy metadata files based on filetype. Output only keeps one row per unique Segment_Id."""
        
        existing_accessions = set()
        if db_file and os.path.exists(db_file):
            try:
                conn = sqlite3.connect(db_file)
                cursor = conn.cursor()
                try:
                    cursor.execute("SELECT primary_accession FROM meta_data")
                    existing_accessions.update(row[0] for row in cursor.fetchall() if row[0])
                except Exception:
                    pass
                try:
                    cursor.execute(
                        "SELECT primary_accession FROM meta_data WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) <> '' "
                        "AND LOWER(TRIM(COALESCE(CAST(exclusion_status AS TEXT), ''))) NOT IN ('', '0', 'false', 'no', 'na', 'none', 'nan')"
                    )
                    existing_accessions.update(row[0] for row in cursor.fetchall() if row[0])
                except Exception:
                    pass
                try:
                    cursor.execute("SELECT primary_accession FROM excluded_accessions")
                    existing_accessions.update(row[0] for row in cursor.fetchall() if row[0])
                except Exception:
                    pass
                conn.close()
                self.log(f"Loaded {len(existing_accessions)} existing accessions from DB to exclude.")
            except Exception as e:
                self.log(f"Error reading DB: {e}")

        tidy_rows = []

        # List of columns on original files to remove
        remove_cols = [
            "PB2 Segment_Id", "PB1 Segment_Id", "PA Segment_Id", "HA Segment_Id", "NP Segment_Id",
            "NA Segment_Id", "MP Segment_Id", "NS Segment_Id", "HE Segment_Id", "P3 Segment_Id",
            "Animal_Vaccin_Product", "Adamantanes_Resistance_geno", "Oseltamivir_Resistance_geno",
            "Zanamivir_Resistance_geno", "Peramivir_Resistance_geno", "Other_Resistance_geno",
            "Adamantanes_Resistance_pheno", "Oseltamivir_Resistance_pheno", "Zanamivir_Resistance_pheno",
            "Peramivir_Resistance_pheno", "Other_Resistance_pheno"
        ]

        for file in metadata_files:
            if file.endswith('.tsv') and self.filetype == "tsv":
                if self.test_mode:
                    df = _restore_blank_as_missing(pd.read_csv(file, sep='\t', dtype=str, nrows=100, keep_default_na=False))
                    self.log(f"Test mode active: loaded first 100 rows from {file}")
                else:
                    df = _restore_blank_as_missing(pd.read_csv(file, sep='\t', dtype=str, keep_default_na=False))
            elif (file.endswith('.xls') or file.endswith('.xlsx')) and self.filetype == "xls":
                df = _restore_blank_as_missing(pd.read_excel(file, dtype=str, keep_default_na=False))
                if self.test_mode:
                    df = df.head(100)
                    self.log(f"Test mode active: keeping first 100 rows from {file}")
            else:
                self.log(f"Skipped unsupported metadata file: {file}")
                continue
            
            for col in df.select_dtypes(include=['object']):
                df[col] = df[col].str.replace(r"[\r\n\f\v]", "", regex=True)
            df = df.dropna(how="all")

            for idx, row in df.iterrows():
                # Handle Long Format (Segment_Id + Segment columns)
                if "Segment_Id" in row and pd.notnull(row["Segment_Id"]):
                    item = str(row["Segment_Id"])
                    m = re.search(r"(EPI\d+)", item)
                    if m:
                        epi = m.group(1)
                        if epi not in existing_accessions:
                            tidy_row = row.to_dict()
                            tidy_row['Segment_Id'] = epi
                            
                            # Ensure Segment is set if available
                            if "Segment" in tidy_row and pd.notnull(tidy_row["Segment"]):
                                tidy_row["Segment"] = str(tidy_row["Segment"]).strip()
                            
                            # Add division if missing
                            div = tidy_row.get("division")
                            if not div or pd.isna(div) or str(div).strip() == "":
                                tidy_row["division"] = "VRL"
                            tidy_rows.append(tidy_row)
                            continue # Skip the wide format check for this row

                # Handle Wide Format (PB1 Segment_Id etc)
                for seg in SEGS:
                    seg_col = f"{seg} Segment_Id"
                    if seg_col in row and pd.notnull(row[seg_col]):
                        for item in re.split(r'[,;]', str(row[seg_col])):
                            m = re.search(r"(EPI\d+)", item)
                            if m:
                                epi = m.group(1)
                                if epi in existing_accessions:
                                    continue

                                tidy_row = row.to_dict()
                                tidy_row['Segment_Id'] = epi
                                tidy_row['Segment'] = seg
                                
                                # Add division if missing
                                div = tidy_row.get("division")
                                if not div or pd.isna(div) or str(div).strip() == "":
                                    tidy_row["division"] = "VRL"

                                tidy_rows.append(tidy_row)
            self.update_log_with_file(file)

        if tidy_rows:
            tidy_df = pd.DataFrame(tidy_rows)
            tidy_df = tidy_df.drop(columns=[col for col in remove_cols if col in tidy_df.columns], errors='ignore')
            
            for col in tidy_df.columns:
                if tidy_df[col].dtype == "object":
                    tidy_df[col] = tidy_df[col].str.strip()
            tidy_df = tidy_df.replace(to_replace=r"[\r\n\f\v]", value="", regex=True)
            tidy_df = tidy_df.dropna(how="all")
            out_path = os.path.join(self.output_dir, "metadata.tsv")
            removed_dups_path = os.path.join(self.output_dir, "metadata_removed_duplicates.tsv")

            if os.path.exists(out_path):
                old_df = _restore_blank_as_missing(pd.read_csv(out_path, sep='\t', dtype=str, keep_default_na=False))
                
                for col in old_df.columns:
                    if old_df[col].dtype == "object":
                        old_df[col] = old_df[col].str.strip()
                old_df = old_df.replace(to_replace=r"[\r\n\f\v]", value="", regex=True)
                old_df = old_df.dropna(how="all")
                common_cols = [col for col in tidy_df.columns if col in old_df.columns]
                tidy_df = tidy_df[common_cols]
                old_df = old_df[common_cols]
                full_df = pd.concat([old_df, tidy_df], ignore_index=True)
            else:
                full_df = tidy_df.copy()

            dup_mask = full_df.duplicated(subset=["Segment_Id"], keep="first")
            if dup_mask.any():
                removed = full_df[dup_mask]
                removed.to_csv(removed_dups_path, sep="\t", index=False)
                self.log(f"{dup_mask.sum()} duplicate rows by Segment_Id removed. Logged in {removed_dups_path}")

            n_before = len(full_df)
            full_df = full_df.drop_duplicates(subset=["Segment_Id"], keep="first")
            n_after = len(full_df)
            n_new = max(0, n_after - (len(old_df) if 'old_df' in locals() else 0))
            n_prev_duplicate = len(tidy_df) - n_new if n_new < len(tidy_df) else 0

            full_df.to_csv(out_path, sep='\t', index=False)
            self.log(f"Metadata cleaned. Output: {out_path}")
            self.log(f"New unique rows added in this update: {n_new}")
            self.log(f"Rows skipped as already present in previous metadata: {n_prev_duplicate}")
        else:
            self.log("No valid metadata rows with EPI# found. Please check file format and column names.")



    def clean_fasta(self, fasta_files):
        """Incremental append: add new sequences (by Segment_Id)."""
        metadata_path = os.path.join(self.output_dir, "metadata.tsv")
        if not os.path.exists(metadata_path):
            self.log("Metadata file not found. Cannot parse FASTA headers robustly.")
            return

        meta_df = _restore_blank_as_missing(pd.read_csv(metadata_path, sep='\t', dtype=str, keep_default_na=False))
        epi_ids_set = set(meta_df['Segment_Id'].dropna()) if 'Segment_Id' in meta_df.columns else set()
        epi_dict = {(row["Isolate_Id"], row["Segment"]): row["Segment_Id"]
                    for _, row in meta_df.iterrows()
                    if "Isolate_Id" in row and "Segment" in row and "Segment_Id" in row}

        out_fasta = os.path.join(self.output_dir, "all_nuc.fas")
        orphan_path = os.path.join(self.output_dir, "orphan_fasta_headers.tsv")
        non_iupac_log_path = os.path.join(self.output_dir, "non_IUPAC_fasta.tsv")

        previous_ids = set()
        if os.path.exists(out_fasta):
            for rec in SeqIO.parse(out_fasta, "fasta"):
                previous_ids.add(rec.id)

        n_new = 0
        n_prev_duplicate = 0
        seen_ids = set(previous_ids)

        mode = "a" if os.path.exists(out_fasta) else "w"
        with open(out_fasta, mode) as out_f, \
            open(orphan_path, "a") as orphan_log, \
            open(non_iupac_log_path, "a") as non_iupac_log:

            for fasta in fasta_files:
                for record in SeqIO.parse(fasta, "fasta"):
                    epi = extract_segment_epi(
                        record.description, epi_ids_set, epi_dict,
                        segment_names=set(SEGS), log_handle=orphan_log
                    )
                    if not epi:
                        continue

                    if epi in seen_ids:
                        if epi in previous_ids:
                            n_prev_duplicate += 1
                        continue

                    seq = str(record.seq).upper()
                    total = len(seq)
                    n_valid = sum((b in ACCEPTABLE_NUCS) for b in seq)
                    ratio_valid = (n_valid / total) if total > 0 else 0.0

                    if ratio_valid >= 0.99: # define threshold
                        record.id = epi
                        record.description = ""
                        SeqIO.write(record, out_f, "fasta")
                        seen_ids.add(epi)
                        n_new += 1
                    else:
                        invalid = "".join(sorted(set(seq) - ACCEPTABLE_NUCS))
                        non_iupac_log.write(
                            f"{record.description}\tratio_valid:{ratio_valid:.4f}\tinvalid_bases:{invalid}\n"
                        )
                self.update_log_with_file(fasta)

        self.log(f"FASTA cleaned. New sequences written: {n_new}. Output: {out_fasta}")
        self.log(f"Orphan headers written to: {orphan_path}")
        self.log(f"Non-IUPAC headers written to: {non_iupac_log_path}")
        self.log(f"New unique sequences added in this update: {n_new}")
        self.log(f"Sequences skipped as already present in previous fasta: {n_prev_duplicate}")


    def validate(self, metadata_file, fasta_file):
        """Cross-validate Segment_Ids between metadata and fasta, report missing in either direction and save lists to txt."""
        meta_df = _restore_blank_as_missing(pd.read_csv(metadata_file, sep='\t', dtype=str, keep_default_na=False))
        meta_epis = set(meta_df['Segment_Id'].dropna())
        fasta_epis = set()
        for record in SeqIO.parse(fasta_file, "fasta"):
            fasta_epis.add(record.id)
        missing_in_fasta = meta_epis - fasta_epis
        missing_in_meta = fasta_epis - meta_epis

        missing_fasta_path = os.path.join(self.output_dir, "missing_in_fasta.txt")
        missing_meta_path = os.path.join(self.output_dir, "missing_in_metadata.txt")
        with open(missing_fasta_path, "w") as f:
            for x in sorted(missing_in_fasta):
                f.write(f"{x}\n")
        with open(missing_meta_path, "w") as f:
            for x in sorted(missing_in_meta):
                f.write(f"{x}\n")

        self.log(f"{len(missing_in_fasta)} IDs in metadata but missing in fasta: {sorted(list(missing_in_fasta))[:3]} (check non_IUPAC_fasta.log for possible excluded sequences)")
        self.log(f"{len(missing_in_meta)} IDs in fasta but missing in metadata: {sorted(list(missing_in_meta))[:3]}")
        self.log(f"Missing IDs written to {missing_fasta_path} and {missing_meta_path}")
        return missing_in_fasta, missing_in_meta


    def update(self):
        """Check for new files to process since last run."""
        if self.filetype == "tsv":
            all_metadata = glob.glob(os.path.join(self.data_dir, "*.tsv"))
        else:
            all_metadata = glob.glob(os.path.join(self.data_dir, "*.xls")) + \
                           glob.glob(os.path.join(self.data_dir, "*.xlsx"))
        all_fasta = glob.glob(os.path.join(self.data_dir, "*.fa*"))

        new_metadata = [f for f in all_metadata if f not in self.seen_files]
        new_fasta = [f for f in all_fasta if f not in self.seen_files]

        if self.test_mode:
            new_metadata = new_metadata[:1]
            new_fasta = new_fasta[:1]
            if new_metadata:
                self.log(f"Test mode active: processing only the first metadata file: {new_metadata[0]}")
            if new_fasta:
                self.log(f"Test mode active: processing only the first fasta file: {new_fasta[0]}")

        if not new_metadata and not new_fasta:
            self.log("No new files to process.")
            return

        if new_metadata:
            self.clean_metadata(new_metadata, db_file=self.db_file)
        if new_fasta:
            self.clean_fasta(new_fasta)

        tidy_meta = os.path.join(self.output_dir, "metadata.tsv")
        tidy_fasta = os.path.join(self.output_dir, "all_nuc.fas")
        if os.path.exists(tidy_meta) and os.path.exists(tidy_fasta):
            self.validate(tidy_meta, tidy_fasta)

    def process(self):
        self.update()

if __name__ == "__main__":
    parser = ArgumentParser(description="Clean and tidy GISAID metadata and fasta files.")
    parser.add_argument("--data_dir", required=True, help="Input directory with raw files (.tsv, .xls, .fas)")
    parser.add_argument("--output_dir", help="Output directory for tidy files", default='tmp/gisaid-data')
    parser.add_argument("--filetype", choices=["xls", "tsv"], default="xls", help="Filetype of metadata to process: 'tsv' or 'xls' (both xls and xlsx)")
    parser.add_argument("--db_file", help="Path to existing SQLite DB to exclude accessions from", default=None)
    parser.add_argument("--test", action="store_true", help="Run in test mode (process limited files only)")
    args = parser.parse_args()

    tidy = GISAIDTidy(data_dir=args.data_dir, output_dir=args.output_dir, filetype=args.filetype, db_file=args.db_file, test_mode=args.test)
    tidy.process()
