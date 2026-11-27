import os
import io
import re
import csv
import glob
import gzip
import lzma
import shutil
import sqlite3
import tarfile
from contextlib import contextmanager
import numpy as np
import pandas as pd
from datetime import datetime
from Bio import SeqIO
from argparse import ArgumentParser
from ExportRefListFromUpdateDb import load_reference_file_accessions

"""
Clean and tidy GISAID metadata and FASTA exports for both:
  - EpiFlu (segmented influenza, keyed on `Segment_Id` = `EPI\\d+`)
  - EpiCoV (unsegmented SARS-CoV-2, keyed on `Accession_Id` = `EPI_ISL_\\d+`)
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
_ACCEPTABLE_DELETE_TABLE = str.maketrans("", "", "".join(sorted(ACCEPTABLE_NUCS)))

EPICOV_ACC_COLUMNS = ("Accession_Id", "Accession ID", "accession_id", "covv_accession_id", "Isolate_Id")
EPICOV_NAME_COLUMNS = ("Virus name", "Virus_name", "virus_name", "covv_virus_name", "Isolate_Name")
EPICOV_COLL_DATE_COLUMNS = ("Collection date", "Collection_Date", "collection_date", "covv_collection_date")
EPICOV_SUBM_DATE_COLUMNS = ("Submission date", "Submission_Date", "submission_date", "covv_subm_date")
EPICOV_LOCATION_COLUMNS = ("Location", "location", "covv_location")


def normalize_epicov_virus_name(name):
    """Canonical stripped virus name for matching EpiCoV metadata <-> FASTA headers."""
    if name is None or (isinstance(name, float) and np.isnan(name)):
        return ""
    text = re.sub(r"\s+", " ", str(name).strip())
    if text.lower().startswith("hcov-19/"):
        text = text[8:].strip()
    return text


def parse_gisaid_location(location_str):
    """Split GISAID 'Continent / Country / Region ...' into (country, geo_loc)."""
    if not location_str or (isinstance(location_str, float) and np.isnan(location_str)):
        return "", ""
    parts = [p.strip() for p in str(location_str).split("/") if p.strip()]
    if len(parts) >= 2:
        country = parts[1]
        geo_loc = " / ".join(parts[2:]) if len(parts) > 2 else ""
        return country, geo_loc
    if len(parts) == 1:
        return parts[0], ""
    return "", ""


def extract_segment_epi(header, epi_ids_set, epi_dict, segment_names=None, log_handle=None):
    """Extract Segment_Id (EPI#) from an EpiFlu fasta header."""
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


def extract_epicov_isl(header, isl_ids_set, triple_map, unique_name_map, log_handle=None):
    """Extract Accession_Id (EPI_ISL_#) from an EpiCoV FASTA header.

    Supports all three standard GISAID EpiCoV FASTA header layouts:
      1. Explicit EPI_ISL_# token: '>hCoV-19/...|EPI_ISL_402124|2020-01-01' or '>EPI_ISL_402124'
      2. Bulk GISAID tar.xz format: '>hCoV-19/Country/Isolate/Year|Collection_date|Submission_date'
      3. Unambiguous virus/isolate name match against metadata.tsv
    """
    clean_header = re.sub(r'\s*\|\s*', '|', header.replace('>', '').strip())
    fields = [f.strip() for f in clean_header.split('|')]

    # 1. Explicit EPI_ISL_#### token anywhere in the header fields
    for field in fields:
        m = re.search(r"\b(EPI_ISL_\d+)\b", field, flags=re.IGNORECASE)
        if m:
            isl = m.group(1).upper()
            if isl in isl_ids_set:
                return isl

    # 2. Bulk GISAID format: Virus name | Collection date | Submission date
    norm_name = normalize_epicov_virus_name(fields[0]) if fields else ""
    if norm_name and len(fields) >= 3:
        triple = (norm_name.casefold(), fields[1], fields[2])
        isl = triple_map.get(triple)
        if isl:
            return isl

    # 3. Unambiguous virus name match (or 2-field Virus name | Collection date)
    if norm_name:
        isl = unique_name_map.get(norm_name.casefold())
        if isl:
            return isl

    if log_handle:
        log_handle.write(f"{header}\tAccession_Id (EPI_ISL) not_found\n")
    return None


def _iter_fasta_blocks_fast(fasta_path, chunk_bytes=16 * 1024 * 1024):
    """Stream (header_without_gt, raw_body_str) in 16 MiB chunks so skipped records avoid per-line Python overhead."""
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
        carry = ""
        first_chunk = True
        while True:
            chunk = fh.read(chunk_bytes)
            if not chunk:
                break
            data = carry + chunk
            if first_chunk:
                first_gt = data.find(">")
                if first_gt == -1:
                    carry = ""
                    continue
                data = data[first_gt + 1:]
                first_chunk = False
            parts = data.split("\n>")
            carry = parts[-1]
            for rec_block in parts[:-1]:
                hdr, sep_nl, body = rec_block.partition("\n")
                if sep_nl or hdr:
                    yield hdr.strip(), body
        if carry and not first_chunk:
            hdr, sep_nl, body = carry.partition("\n")
            if sep_nl or hdr:
                yield hdr.strip(), body


def _iter_fasta_fast(fasta_path):
    """Stream (header_without_gt, uppercase_seq) from a FASTA file."""
    for hdr, body in _iter_fasta_blocks_fast(fasta_path):
        yield hdr, body.replace("\n", "").replace("\r", "").replace(" ", "").upper()


class GISAIDTidy:
    # Columns safe to drop before reading (wide-format `* Segment_Id` columns must
    # stay until after wide-to-long expansion in clean_metadata).
    READ_DROP_COLS = frozenset({
        "Animal_Vaccin_Product", "Adamantanes_Resistance_geno", "Oseltamivir_Resistance_geno",
        "Zanamivir_Resistance_geno", "Peramivir_Resistance_geno", "Other_Resistance_geno",
        "Adamantanes_Resistance_pheno", "Oseltamivir_Resistance_pheno", "Zanamivir_Resistance_pheno",
        "Peramivir_Resistance_pheno", "Other_Resistance_pheno",
        "AA Substitutions", "AA_Substitutions",
    })

    def __init__(self, data_dir, output_dir, filetype="xls", log_file=None, db_file=None, test_mode=False, mode="auto", ref_list=None, min_seq_length=None, max_n_count=None):
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.filetype = filetype
        self.log_file = log_file or os.path.abspath("gisaid_processing.log")
        self.db_file = db_file
        self.test_mode = test_mode
        self.mode = (mode or "auto").lower()
        self.ref_list = ref_list
        self.min_seq_length = int(min_seq_length) if min_seq_length is not None and str(min_seq_length).strip().lower() not in ("", "null", "none", "unset") else None
        self.max_n_count = int(max_n_count) if max_n_count is not None and str(max_n_count).strip().lower() not in ("", "null", "none", "unset") else None
        self._resolved_mode = None
        if self.test_mode:
            self.log_file = log_file or os.path.abspath("test_gisaid_processing.log")
        os.makedirs(self.output_dir, exist_ok=True)
        self.seen_files = self.load_seen_files()
        self.processed_files = []
        self.log(f"Pipeline started at {self.timestamp()}. Filetype: {self.filetype}. Mode: {self.mode}. Test mode: {self.test_mode}")

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

    def _load_required_ref_ids(self):
        if not self.ref_list or not os.path.isfile(self.ref_list):
            return set()
        try:
            return {str(a).strip() for a in load_reference_file_accessions(self.ref_list) if str(a).strip()}
        except Exception as exc:
            self.log(f"Warning: could not read ref_list {self.ref_list}: {exc}")
            return set()

    def _detect_table_mode(self, columns):
        if self.mode in ("epiflu", "epicov"):
            return self.mode
        col_set = set(columns)
        has_flu_seg = "Segment_Id" in col_set or any(f"{seg} Segment_Id" in col_set for seg in SEGS)
        if has_flu_seg:
            return "epiflu"
        if any(c in col_set for c in EPICOV_ACC_COLUMNS):
            return "epicov"
        return "epiflu"

    @staticmethod
    def _first_present_col(columns, candidates):
        col_set = set(columns)
        for c in candidates:
            if c in col_set:
                return c
        return None

    def _load_existing_db_accessions(self, db_file):
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
        return existing_accessions

    @staticmethod
    @contextmanager
    def _open_delimited_stream(file):
        lower = file.lower()
        if lower.endswith((".tar.xz", ".txz", ".tar.gz", ".tgz")):
            tar = tarfile.open(file, mode="r:*")
            try:
                yielded = False
                for member in tar:
                    if member.isfile() and member.name.lower().endswith((".tsv", ".csv")):
                        raw = tar.extractfile(member)
                        if raw is not None:
                            sep = "," if member.name.lower().endswith(".csv") else "\t"
                            yield io.TextIOWrapper(raw, encoding="utf-8", errors="replace"), sep
                            yielded = True
                            break
                if not yielded:
                    raise ValueError(f"No .tsv/.csv member found inside archive: {file}")
            finally:
                tar.close()
            return
        if lower.endswith(".gz"):
            sep = "," if lower.endswith(".csv.gz") else "\t"
            with gzip.open(file, "rb") as raw:
                yield io.TextIOWrapper(raw, encoding="utf-8", errors="replace"), sep
            return
        if lower.endswith(".xz"):
            sep = "," if lower.endswith(".csv.xz") else "\t"
            with lzma.open(file, "rb") as raw:
                yield io.TextIOWrapper(raw, encoding="utf-8", errors="replace"), sep
            return
        sep = "," if lower.endswith(".csv") else "\t"
        with open(file, "r", encoding="utf-8", errors="replace") as fh:
            yield fh, sep

    @staticmethod
    def _has_supported_gisaid_header(columns):
        col_set = {str(c).strip() for c in columns}
        if "Segment_Id" in col_set or any(f"{seg} Segment_Id" in col_set for seg in SEGS):
            return True
        if any(c in col_set for c in EPICOV_ACC_COLUMNS):
            return True
        return False

    def _apply_metadata_pre_qc(self, chunk):
        """Filter query chunk by optional min_seq_length and max_n_count when GISAID length/N-Content columns exist."""
        if self.min_seq_length is None and self.max_n_count is None:
            return chunk
        if chunk.empty:
            return chunk
        len_col = self._first_present_col(chunk.columns, ("Sequence length", "sequence_length", "covv_length"))
        n_col = self._first_present_col(chunk.columns, ("N-Content", "n-content", "n_content"))
        mask = pd.Series(True, index=chunk.index)
        seq_lens = None
        if len_col:
            seq_lens = pd.to_numeric(chunk[len_col], errors="coerce")
            if self.min_seq_length is not None:
                # Keep rows where length is missing/unreported OR >= min_seq_length
                has_len = chunk[len_col].fillna("").astype(str).str.strip().ne("")
                mask &= (~has_len) | (seq_lens >= self.min_seq_length)
        if n_col and self.max_n_count is not None:
            has_n = chunk[n_col].fillna("").astype(str).str.strip().ne("")
            n_frac = pd.to_numeric(chunk[n_col], errors="coerce").fillna(0.0)
            if seq_lens is not None:
                est_n = n_frac * seq_lens.fillna(29903.0)
                mask &= (~has_n) | (est_n <= (self.max_n_count + 0.5))
            else:
                mask &= (~has_n) | (n_frac <= (self.max_n_count / 29903.0 + 1e-6))
        return chunk.loc[mask]

    def _read_metadata_df(self, file, required_refs):
        lower = file.lower()
        is_delimited = (
            lower.endswith((".tsv", ".csv", ".tsv.gz", ".csv.gz", ".tsv.xz", ".csv.xz", ".tar.xz", ".txz", ".tar.gz", ".tgz"))
            and self.filetype in ("tsv", "csv")
        )
        if is_delimited:
            with self._open_delimited_stream(file) as (stream, sep):
                header_df = pd.read_csv(stream, sep=sep, dtype=str, nrows=0, keep_default_na=False)
                if not self._has_supported_gisaid_header(header_df.columns):
                    return None
                keep_cols = [c for c in header_df.columns if c not in self.READ_DROP_COLS]

            has_pre_qc = (self.min_seq_length is not None) or (self.max_n_count is not None)
            if self.test_mode and not required_refs and not has_pre_qc:
                with self._open_delimited_stream(file) as (stream, sep):
                    df = _restore_blank_as_missing(
                        pd.read_csv(stream, sep=sep, dtype=str, nrows=100, usecols=keep_cols, keep_default_na=False)
                    )
                self.log(f"Test mode active: loaded first 100 rows from {file}")
                return df
            if self.test_mode or has_pre_qc:
                kept_chunks = []
                query_budget = 100 if self.test_mode else None
                found_refs = set()
                with self._open_delimited_stream(file) as (stream, sep):
                    for chunk in pd.read_csv(stream, sep=sep, dtype=str, chunksize=250000, usecols=keep_cols, keep_default_na=False):
                        acc_col = self._first_present_col(chunk.columns, ("Segment_Id",) + EPICOV_ACC_COLUMNS)
                        if acc_col and required_refs:
                            is_ref = chunk[acc_col].isin(required_refs)
                            ref_rows = chunk[is_ref]
                            non_ref = chunk[~is_ref]
                            if len(ref_rows) > 0:
                                kept_chunks.append(ref_rows)
                                found_refs.update(ref_rows[acc_col].tolist())
                        else:
                            non_ref = chunk
                        non_ref = self._apply_metadata_pre_qc(non_ref)
                        if query_budget is None:
                            if len(non_ref) > 0:
                                kept_chunks.append(non_ref)
                        elif query_budget > 0 and len(non_ref) > 0:
                            take = non_ref.head(query_budget)
                            kept_chunks.append(take)
                            query_budget -= len(take)
                        if self.test_mode and query_budget == 0 and (not required_refs or found_refs >= {r for r in required_refs if r.startswith("EPI")}):
                            break
                df = pd.concat(kept_chunks, ignore_index=True) if kept_chunks else pd.DataFrame(columns=keep_cols)
                df = _restore_blank_as_missing(df)
                self.log(f"Loaded {len(df)} metadata rows from {file} (test={self.test_mode}, min_len={self.min_seq_length}, max_n={self.max_n_count})")
                return df
            with self._open_delimited_stream(file) as (stream, sep):
                return _restore_blank_as_missing(
                    pd.read_csv(stream, sep=sep, dtype=str, usecols=keep_cols, keep_default_na=False)
                )
        if (lower.endswith('.xls') or lower.endswith('.xlsx')) and self.filetype == "xls":
            df = _restore_blank_as_missing(pd.read_excel(file, dtype=str, keep_default_na=False))
            if self.test_mode:
                df = df.head(100)
                self.log(f"Test mode active: keeping first 100 rows from {file}")
            return df
        return None

    def clean_metadata(self, metadata_files, db_file=None):
        """Parse and tidy metadata files based on filetype and EpiFlu/EpiCoV schema."""
        existing_accessions = self._load_existing_db_accessions(db_file)
        required_refs = self._load_required_ref_ids()
        tidy_dfs = []
        tidy_rows = []
        key_col = "Segment_Id"

        remove_cols = [
            "PB2 Segment_Id", "PB1 Segment_Id", "PA Segment_Id", "HA Segment_Id", "NP Segment_Id",
            "NA Segment_Id", "MP Segment_Id", "NS Segment_Id", "HE Segment_Id", "P3 Segment_Id",
            "Animal_Vaccin_Product", "Adamantanes_Resistance_geno", "Oseltamivir_Resistance_geno",
            "Zanamivir_Resistance_geno", "Peramivir_Resistance_geno", "Other_Resistance_geno",
            "Adamantanes_Resistance_pheno", "Oseltamivir_Resistance_pheno", "Zanamivir_Resistance_pheno",
            "Peramivir_Resistance_pheno", "Other_Resistance_pheno"
        ]

        for file in metadata_files:
            df = self._read_metadata_df(file, required_refs)
            if df is None:
                self.log(f"Skipped unsupported metadata file: {file}")
                continue
            if df.empty:
                self.update_log_with_file(file)
                continue

            for col in df.select_dtypes(include=['object']):
                df[col] = df[col].str.replace(r"[\r\n\f\v]", "", regex=True)
            df = df.dropna(how="all")

            table_mode = self._detect_table_mode(df.columns)
            self._resolved_mode = table_mode

            if table_mode == "epicov":
                key_col = "Accession_Id"
                acc_src = self._first_present_col(df.columns, EPICOV_ACC_COLUMNS)
                loc_src = self._first_present_col(df.columns, EPICOV_LOCATION_COLUMNS)
                name_src = self._first_present_col(df.columns, EPICOV_NAME_COLUMNS)
                if not acc_src:
                    self.update_log_with_file(file)
                    continue

                extracted = df[acc_src].fillna("").astype(str).str.extract(r"(EPI_ISL_\d+)", flags=re.IGNORECASE, expand=False).str.upper()
                valid_mask = extracted.notna()
                if existing_accessions:
                    valid_mask &= ~extracted.isin(existing_accessions)
                sub = df.loc[valid_mask].copy()
                if sub.empty:
                    self.update_log_with_file(file)
                    continue
                sub["Accession_Id"] = extracted.loc[valid_mask]
                if "division" not in sub.columns:
                    sub["division"] = "VRL"
                else:
                    div_s = sub["division"].fillna("").astype(str).str.strip()
                    sub["division"] = div_s.where(div_s != "", "VRL")

                if loc_src:
                    parsed_locs = [parse_gisaid_location(v) for v in sub[loc_src].fillna("").tolist()]
                    countries = pd.Series([p[0] for p in parsed_locs], index=sub.index, dtype=object)
                    geo_locs = pd.Series([p[1] for p in parsed_locs], index=sub.index, dtype=object)
                    if "Country" not in sub.columns:
                        sub["Country"] = countries.replace("", np.nan)
                    else:
                        cur_c = sub["Country"].fillna("").astype(str).str.strip()
                        sub["Country"] = cur_c.where(cur_c != "", countries).replace("", np.nan)
                    if "Geo_Loc" not in sub.columns:
                        sub["Geo_Loc"] = geo_locs.replace("", np.nan)
                    else:
                        cur_g = sub["Geo_Loc"].fillna("").astype(str).str.strip()
                        sub["Geo_Loc"] = cur_g.where(cur_g != "", geo_locs).replace("", np.nan)
                elif name_src:
                    derived_c = []
                    for v in sub[name_src].fillna("").tolist():
                        norm_vname = normalize_epicov_virus_name(v)
                        vparts = [p.strip() for p in norm_vname.split("/") if p.strip()]
                        derived_c.append(vparts[0] if len(vparts) >= 2 else "")
                    countries = pd.Series(derived_c, index=sub.index, dtype=object)
                    if "Country" not in sub.columns:
                        sub["Country"] = countries.replace("", np.nan)
                    else:
                        cur_c = sub["Country"].fillna("").astype(str).str.strip()
                        sub["Country"] = cur_c.where(cur_c != "", countries).replace("", np.nan)

                tidy_dfs.append(sub)
            else:
                key_col = "Segment_Id"
                for _, row in df.iterrows():
                    # Handle Long Format (Segment_Id + Segment columns)
                    if "Segment_Id" in row and pd.notnull(row["Segment_Id"]):
                        item = str(row["Segment_Id"])
                        m = re.search(r"(EPI\d+)", item)
                        if m:
                            epi = m.group(1)
                            if epi not in existing_accessions:
                                tidy_row = row.to_dict()
                                tidy_row['Segment_Id'] = epi
                                if "Segment" in tidy_row and pd.notnull(tidy_row["Segment"]):
                                    tidy_row["Segment"] = str(tidy_row["Segment"]).strip()
                                div = tidy_row.get("division")
                                if not div or pd.isna(div) or str(div).strip() == "":
                                    tidy_row["division"] = "VRL"
                                tidy_rows.append(tidy_row)
                                continue

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
                                    div = tidy_row.get("division")
                                    if not div or pd.isna(div) or str(div).strip() == "":
                                        tidy_row["division"] = "VRL"
                                    tidy_rows.append(tidy_row)
            self.update_log_with_file(file)

        if tidy_rows:
            tidy_dfs.append(pd.DataFrame(tidy_rows))

        if tidy_dfs:
            tidy_df = pd.concat(tidy_dfs, ignore_index=True)
            tidy_df = tidy_df.drop(columns=[col for col in remove_cols if col in tidy_df.columns], errors='ignore')

            for col in tidy_df.columns:
                if tidy_df[col].dtype == "object":
                    tidy_df[col] = tidy_df[col].str.strip()
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

            dup_mask = full_df.duplicated(subset=[key_col], keep="first")
            if dup_mask.any():
                removed = full_df[dup_mask]
                removed.to_csv(removed_dups_path, sep="\t", index=False)
                self.log(f"{dup_mask.sum()} duplicate rows by {key_col} removed. Logged in {removed_dups_path}")

            full_df = full_df.drop_duplicates(subset=[key_col], keep="first")
            n_after = len(full_df)
            n_new = max(0, n_after - (len(old_df) if 'old_df' in locals() else 0))
            n_prev_duplicate = len(tidy_df) - n_new if n_new < len(tidy_df) else 0

            full_df.to_csv(out_path, sep='\t', index=False)
            self.log(f"Metadata cleaned ({self._resolved_mode}, key={key_col}). Output: {out_path}")
            self.log(f"New unique rows added in this update: {n_new}")
            self.log(f"Rows skipped as already present in previous metadata: {n_prev_duplicate}")
        else:
            self.log("No valid metadata rows with EPI# / EPI_ISL_# found. Please check file format and column names.")

    def _build_epicov_lookup(self, meta_df, required_refs=None):
        required_refs = set(required_refs or ())
        isl_ids_set = set(meta_df["Accession_Id"].dropna()) if "Accession_Id" in meta_df.columns else set()
        name_col = self._first_present_col(meta_df.columns, EPICOV_NAME_COLUMNS)
        coll_col = self._first_present_col(meta_df.columns, EPICOV_COLL_DATE_COLUMNS)
        subm_col = self._first_present_col(meta_df.columns, EPICOV_SUBM_DATE_COLUMNS)

        triple_map = {}
        name_candidates = {}
        ref_name_map = {}
        if name_col and "Accession_Id" in meta_df.columns:
            isl_vals = meta_df["Accession_Id"].fillna("").tolist()
            name_vals = meta_df[name_col].fillna("").tolist()
            coll_vals = meta_df[coll_col].fillna("").tolist() if coll_col else [""] * len(isl_vals)
            subm_vals = meta_df[subm_col].fillna("").tolist() if subm_col else [""] * len(isl_vals)
            for isl, raw_name, raw_cdate, raw_sdate in zip(isl_vals, name_vals, coll_vals, subm_vals):
                if not isl:
                    continue
                norm_name = normalize_epicov_virus_name(raw_name)
                if not norm_name:
                    continue
                key_name = norm_name.casefold()
                name_candidates.setdefault(key_name, []).append(isl)
                cdate = str(raw_cdate).strip()
                sdate = str(raw_sdate).strip()
                is_req_ref = isl in required_refs
                if is_req_ref:
                    ref_name_map[key_name] = isl
                if cdate or sdate:
                    tkey = (key_name, cdate, sdate)
                    if is_req_ref or tkey not in triple_map:
                        triple_map[tkey] = isl

        unique_name_map = {k: v[0] for k, v in name_candidates.items() if len(v) == 1}
        unique_name_map.update(ref_name_map)
        return isl_ids_set, triple_map, unique_name_map

    def clean_fasta(self, fasta_files):
        """Incremental append: add new sequences keyed by Segment_Id (EpiFlu) or Accession_Id (EpiCoV)."""
        metadata_path = os.path.join(self.output_dir, "metadata.tsv")
        if not os.path.exists(metadata_path):
            self.log("Metadata file not found. Cannot parse FASTA headers robustly.")
            return

        required_refs = self._load_required_ref_ids()
        meta_df = _restore_blank_as_missing(pd.read_csv(metadata_path, sep='\t', dtype=str, keep_default_na=False))
        is_epicov = "Accession_Id" in meta_df.columns and "Segment_Id" not in meta_df.columns

        if is_epicov:
            target_ids_set, triple_map, unique_name_map = self._build_epicov_lookup(meta_df, required_refs=required_refs)
        else:
            target_ids_set = set(meta_df['Segment_Id'].dropna()) if 'Segment_Id' in meta_df.columns else set()
            epi_dict = {
                (row["Isolate_Id"], row["Segment"]): row["Segment_Id"]
                for _, row in meta_df.iterrows()
                if "Isolate_Id" in row and "Segment" in row and "Segment_Id" in row
            }

        out_fasta = os.path.join(self.output_dir, "all_nuc.fas")
        orphan_path = os.path.join(self.output_dir, "orphan_fasta_headers.tsv")
        non_iupac_log_path = os.path.join(self.output_dir, "non_IUPAC_fasta.tsv")

        required_refs = self._load_required_ref_ids()
        has_pre_qc = (self.min_seq_length is not None) or (self.max_n_count is not None)
        orphan_handle_for_lookup = None if has_pre_qc else True

        previous_ids = set()
        if os.path.exists(out_fasta):
            for header, _ in _iter_fasta_blocks_fast(out_fasta):
                previous_ids.add(header.split()[0])

        n_new = 0
        n_prev_duplicate = 0
        n_pre_qc_dropped = 0
        seen_ids = set(previous_ids)

        mode = "a" if os.path.exists(out_fasta) else "w"
        with open(out_fasta, mode, encoding="utf-8") as out_f, \
             open(orphan_path, "a", encoding="utf-8") as orphan_log, \
             open(non_iupac_log_path, "a", encoding="utf-8") as non_iupac_log:

            eff_orphan_log = orphan_log if orphan_handle_for_lookup else None
            for fasta in fasta_files:
                for raw_header, raw_body in _iter_fasta_blocks_fast(fasta):
                    if is_epicov:
                        seq_id = extract_epicov_isl(
                            raw_header, target_ids_set, triple_map, unique_name_map, log_handle=eff_orphan_log
                        )
                    else:
                        seq_id = extract_segment_epi(
                            raw_header, target_ids_set, epi_dict,
                            segment_names=set(SEGS), log_handle=eff_orphan_log
                        )
                    if not seq_id:
                        continue

                    if seq_id in seen_ids:
                        if seq_id in previous_ids:
                            n_prev_duplicate += 1
                        continue

                    seq = raw_body.replace("\n", "").replace("\r", "").replace(" ", "").upper()
                    total = len(seq)
                    if has_pre_qc and seq_id not in required_refs:
                        if self.min_seq_length is not None and total < self.min_seq_length:
                            n_pre_qc_dropped += 1
                            continue
                        if self.max_n_count is not None and seq.count("N") > self.max_n_count:
                            n_pre_qc_dropped += 1
                            continue

                    n_invalid = len(seq.translate(_ACCEPTABLE_DELETE_TABLE))
                    n_valid = total - n_invalid
                    ratio_valid = (n_valid / total) if total > 0 else 0.0

                    if ratio_valid >= 0.99:
                        out_f.write(f">{seq_id}\n{seq}\n")
                        seen_ids.add(seq_id)
                        n_new += 1
                    else:
                        invalid = "".join(sorted(set(seq) - ACCEPTABLE_NUCS))
                        non_iupac_log.write(
                            f"{raw_header}\tratio_valid:{ratio_valid:.4f}\tinvalid_bases:{invalid}\n"
                        )
                    if self.test_mode and len(seen_ids) >= len(target_ids_set) and len(target_ids_set) > 0:
                        # All metadata targets found; avoid scanning multi-hundred-GB FASTA further in test mode
                        break
                self.update_log_with_file(fasta)

        if n_pre_qc_dropped:
            self.log(f"Sequences skipped by FASTA pre-QC (min_seq_length={self.min_seq_length}, max_n_count={self.max_n_count}): {n_pre_qc_dropped}")

        self.log(f"FASTA cleaned. New sequences written: {n_new}. Output: {out_fasta}")
        self.log(f"Orphan headers written to: {orphan_path}")
        self.log(f"Non-IUPAC headers written to: {non_iupac_log_path}")
        self.log(f"New unique sequences added in this update: {n_new}")
        self.log(f"Sequences skipped as already present in previous fasta: {n_prev_duplicate}")

    def validate(self, metadata_file, fasta_file):
        """Cross-validate IDs between metadata and fasta, report missing in either direction."""
        meta_df = _restore_blank_as_missing(pd.read_csv(metadata_file, sep='\t', dtype=str, keep_default_na=False))
        key_col = "Accession_Id" if ("Accession_Id" in meta_df.columns and "Segment_Id" not in meta_df.columns) else "Segment_Id"
        meta_epis = set(meta_df[key_col].dropna()) if key_col in meta_df.columns else set()
        fasta_epis = {header.split()[0] for header, _ in _iter_fasta_blocks_fast(fasta_file)}
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
        out_meta_abs = os.path.abspath(os.path.join(self.output_dir, "metadata.tsv"))
        if self.filetype == "tsv":
            tsv_like = sorted(
                f for f in (
                    glob.glob(os.path.join(self.data_dir, "*.tsv")) +
                    glob.glob(os.path.join(self.data_dir, "*.tsv.gz")) +
                    glob.glob(os.path.join(self.data_dir, "*.tsv.xz")) +
                    glob.glob(os.path.join(self.data_dir, "*metadata*.tar.xz")) +
                    glob.glob(os.path.join(self.data_dir, "*metadata*.tar.gz"))
                )
                if os.path.abspath(f) != out_meta_abs
            )
            if tsv_like:
                plain_tsv = [f for f in tsv_like if not f.lower().endswith((".tar.xz", ".tar.gz"))]
                all_metadata = plain_tsv if plain_tsv else tsv_like
            else:
                all_metadata = sorted(
                    f for f in glob.glob(os.path.join(self.data_dir, "*.csv"))
                    if os.path.basename(f).lower() != "lineages.csv"
                )
        elif self.filetype == "csv":
            all_metadata = sorted(
                f for f in (
                    glob.glob(os.path.join(self.data_dir, "*.csv")) +
                    glob.glob(os.path.join(self.data_dir, "*.csv.gz")) +
                    glob.glob(os.path.join(self.data_dir, "*.tsv"))
                )
                if os.path.abspath(f) != out_meta_abs and os.path.basename(f).lower() != "lineages.csv"
            )
        else:
            all_metadata = sorted(
                glob.glob(os.path.join(self.data_dir, "*.xls")) +
                glob.glob(os.path.join(self.data_dir, "*.xlsx"))
            )
        all_fasta = sorted(
            f for f in glob.glob(os.path.join(self.data_dir, "*.fa*"))
            if f.lower().endswith((".fa", ".fas", ".fasta", ".fna"))
        )

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
    parser = ArgumentParser(description="Clean and tidy GISAID metadata and fasta files (EpiFlu & EpiCoV).")
    parser.add_argument("--data_dir", required=True, help="Input directory with raw files (.tsv, .xls, .fas)")
    parser.add_argument("--output_dir", help="Output directory for tidy files", default='tmp/gisaid-data')
    parser.add_argument("--filetype", choices=["xls", "tsv", "csv"], default="xls", help="Filetype of metadata to process: 'tsv', 'csv', or 'xls'")
    parser.add_argument("--db_file", help="Path to existing SQLite DB to exclude accessions from", default=None)
    parser.add_argument("--mode", choices=["auto", "epiflu", "epicov"], default="auto", help="GISAID schema mode (default: auto)")
    parser.add_argument("--ref_list", help="Optional reference list TSV so test mode retains reference accessions", default=None)
    parser.add_argument("--min_seq_length", type=int, default=None, help="Optional minimum sequence length pre-filter on query rows")
    parser.add_argument("--max_n_count", type=int, default=None, help="Optional maximum N count pre-filter on query rows")
    parser.add_argument("--test", action="store_true", help="Run in test mode (process limited rows/files only)")
    args = parser.parse_args()

    tidy = GISAIDTidy(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        filetype=args.filetype,
        db_file=args.db_file,
        test_mode=args.test,
        mode=args.mode,
        ref_list=args.ref_list,
        min_seq_length=args.min_seq_length,
        max_n_count=args.max_n_count,
    )
    tidy.process()
