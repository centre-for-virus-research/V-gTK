import os
import csv
import re
from os.path import join, dirname, abspath
from argparse import ArgumentParser
from datetime import datetime
from collections import OrderedDict
import pandas as pd
from date_utils import split_date_components
from ExportRefListFromUpdateDb import load_reference_file_dict

"""
Transform an external TSV + FASTA into the gB_matrix-compatible format and append rows to the base gB_matrix.
- Key in input TSV: configurable with --key, default 'Segment_Id'
- Copies the key into: primary_accession, accession_version, locus columns
- Adds sequences from FASTA (header must match the key exactly)
- If a TSV row has no matching FASTA, sets exclusion = "missing fasta sequence"
- Duplicates in input by key: keeps the first, writes the others to merge_matrix_duplicates.tsv and log count
- Column mapping file: two-column TSV (source_col<TAB>target_col); mapped columns are copied/renamed
- Merge strategy: append-only (no existence check in base); writes a new merged file
"""

class NormalizeAndMerge:
    def __init__(self, gb_matrix, input_tsv, input_fasta, mapping_file, output, key="Segment_Id", dataset_source="external", prefer="input", drop_unmapped=False, log_file=None, ref_list=None, require_refs=False):
        self.gb_matrix = gb_matrix
        self.input_tsv = input_tsv
        self.input_fasta = input_fasta
        self.mapping_file = mapping_file
        self.output = output
        self.key = key
        self.dataset_source = dataset_source
        self.prefer = prefer                 # unused in append-only
        self.drop_unmapped = drop_unmapped
        self.ref_list = ref_list
        self.require_refs = bool(require_refs)

        out_dir = dirname(abspath(self.output)) or "."
        os.makedirs(out_dir, exist_ok=True)

        self.log_file = log_file or abspath("gisaid_processing.log")

    def timestamp(self):
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def log(self, message):
        with open(self.log_file, "a") as logf:
            logf.write(f"[{self.timestamp()}] {message}\n")
        print(f"[{self.timestamp()}] {message}")

    def read_mapping(self):
        """Read source_col<TAB>target_col pairs."""
        mapping = OrderedDict()
        with open(self.mapping_file, mode='r') as mf:
            reader = csv.reader(mf, delimiter='\t')
            for row in reader:
                if not row or len(row) < 2:
                    continue
                src, dst = row[0].strip(), row[1].strip()
                if src and dst:
                    mapping[src] = dst
        return mapping

    def parse_fasta(self):
        """Parse FASTA into sequence dict"""
        seqs = {}
        if not self.input_fasta:
            return seqs

        current_header = None
        chunks = []

        with open(self.input_fasta, 'r') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header is not None:
                        seqs[current_header] = ''.join(chunks)
                    current_header = line[1:].strip()
                    chunks = []
                else:
                    chunks.append(line.strip())
            if current_header is not None:
                seqs[current_header] = ''.join(chunks)

        return seqs

    def parse_fasta_stats(self):
        """Stream FASTA and return lightweight {header: (total_len, real_len, a, t, g, c, n)} without retaining sequence text."""
        stats = {}
        if not self.input_fasta:
            return stats

        def _summarize(seq_str):
            up = seq_str.upper()
            a = up.count("A")
            t = up.count("T")
            g = up.count("G")
            c = up.count("C")
            n = up.count("N")
            return (len(up), a + t + g + c, a, t, g, c, n)

        current_header = None
        chunks = []
        with open(self.input_fasta, 'r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                line = line.rstrip('\r\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header is not None:
                        stats[current_header] = _summarize(''.join(chunks))
                    current_header = line[1:].strip()
                    chunks = []
                else:
                    chunks.append(line.strip())
            if current_header is not None:
                stats[current_header] = _summarize(''.join(chunks))
        return stats

    def count_real_length(self, seq):
        if not isinstance(seq, str):
            return 0
        seq_lower = seq.lower()
        return sum(seq_lower.count(n) for n in ("a", "t", "g", "c"))
    
    def read_base_matrix_df(self) -> pd.DataFrame:
        # keep_default_na=False: pandas' default NA sentinels include the literal
        # string "NA", which is influenza's neuraminidase segment/gene name. Without
        # this, every NA-segment row silently loses its label on read.
        base_df = pd.read_csv(self.gb_matrix, sep="\t", dtype=str, keep_default_na=False)
        return base_df

    def read_input_tsv_df(self) -> pd.DataFrame:
        in_df = pd.read_csv(self.input_tsv, sep="\t", dtype=str, keep_default_na=False)
        return in_df

    def write_duplicates_report_df(self, dup_df: pd.DataFrame):
        if dup_df is None or dup_df.empty:
            return None
        dup_file = join(dirname(abspath(self.output)) or ".", 'merge_matrix_duplicates.tsv')
        dup_df.to_csv(dup_file, sep="\t", index=False)
        return dup_file

    #: Columns this class writes itself. They are canonical by construction and
    #: must never be namespaced away.
    PIPELINE_OWNED_COLUMNS = frozenset({
        "primary_accession", "accession_version", "locus", "gi_number",
        "sequence", "real_length", "dataset_source", "exclusion",
    })

    @staticmethod
    def slugify_column(name: str) -> str:
        """Reduce an arbitrary vendor header to a legal, lowercase SQL identifier.

        ``"HA INSDC_Upload"`` -> ``"ha_insdc_upload"``. Unquoted identifiers with
        spaces are not merely ugly: ``ALTER TABLE t ADD COLUMN HA INSDC_Upload TEXT``
        parses as a column ``HA`` of type ``INSDC_Upload TEXT``, silently creating
        the wrong column instead of raising.
        """
        slug = re.sub(r'[^0-9a-zA-Z]+', '_', str(name)).strip('_').lower()
        if not slug:
            slug = "column"
        if slug[0].isdigit():
            slug = f"col_{slug}"
        return slug

    @staticmethod
    def case_insensitive_duplicates(columns):
        seen, dupes = set(), []
        for col in columns:
            key = str(col).casefold()
            if key in seen:
                dupes.append(col)
            else:
                seen.add(key)
        return dupes

    @classmethod
    def columns_collide_ignoring_case(cls, columns) -> bool:
        return bool(cls.case_insensitive_duplicates(columns))

    def namespace_external_columns(self, df: pd.DataFrame, mapped_targets, canonical=()) -> pd.DataFrame:
        """Prefix un-mapped vendor columns with the dataset source.

        Exempt: deliberate mapping targets, the columns this class writes itself,
        and the join key. Everything else becomes ``<source>_<slug>`` so it cannot
        collide with a canonical column under SQLite's case-insensitive
        identifier rules.

        Idempotent: a column that already carries the prefix is left alone, so
        re-merging an already-merged matrix does not stack prefixes.
        """
        prefix = self.slugify_column(self.dataset_source or "external")
        # An external column spelled EXACTLY like a canonical one is not a
        # collision - it is the same field, and renaming it would silently move
        # the vendor's data out of the column the pipeline reads. Only differing
        # spellings (including case-only differences) get namespaced.
        protected = set(mapped_targets) | set(self.PIPELINE_OWNED_COLUMNS) | {self.key} | set(canonical)

        renames, taken = {}, {str(c).casefold() for c in protected}
        for col in df.columns:
            if col in protected:
                continue
            slug = self.slugify_column(col)
            new = slug if slug.startswith(f"{prefix}_") else f"{prefix}_{slug}"
            candidate, n = new, 2
            while candidate.casefold() in taken:
                candidate, n = f"{new}_{n}", n + 1
            taken.add(candidate.casefold())
            if candidate != col:
                renames[col] = candidate

        if renames:
            shown = {k: v for k, v in list(renames.items())[:8]}
            self.log(
                f"[INFO] Namespaced {len(renames)} un-mapped '{self.dataset_source}' "
                f"column(s) to avoid collision with canonical names, e.g. {shown}"
            )
            df = df.rename(columns=renames)
        return df

    def collapse_duplicate_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Pandas `reindex(columns=...)` fails when the source frame has duplicate
        column labels. This can happen after mapping/rename if a target column
        already exists in the input TSV.

        Strategy: keep one column per label and coalesce values left-to-right,
        filling empty/NaN entries with later duplicates.

        Matching is case-INSENSITIVE, because the destination is SQLite and its
        column identifiers are. Detecting a `segment`/`Segment` pair and then
        handing it to an exact-match collapse would log "Collapsing duplicates"
        and change nothing, leaving the original `duplicate column name` crash
        one step further downstream. The first spelling seen wins the name.
        """
        keys = [str(c).casefold() for c in df.columns]
        if len(keys) == len(set(keys)):
            return df

        collapsed = pd.DataFrame(index=df.index)
        seen_keys = []
        for key in keys:
            if key not in seen_keys:
                seen_keys.append(key)
        for key in seen_keys:
            mask = [k == key for k in keys]
            same_cols = df.loc[:, mask]
            col = df.columns[mask.index(True)]
            if same_cols.shape[1] == 1:
                collapsed[col] = same_cols.iloc[:, 0]
                continue

            merged = same_cols.iloc[:, 0].copy()
            for i in range(1, same_cols.shape[1]):
                candidate = same_cols.iloc[:, i]
                empty_mask = merged.isna() | merged.astype(str).str.strip().eq("")
                merged = merged.where(~empty_mask, candidate)
            collapsed[col] = merged

        return collapsed



    def process(self):
        base_df = self.read_base_matrix_df()
        self.log(f"[INFO] Loaded base matrix: {self.gb_matrix} (rows={len(base_df)})")
        base_df["dataset_source"] = "genbank" # define source name
        base_col_set = set(base_df.columns)
        mapping = self.read_mapping()
        in_df = self.read_input_tsv_df()
        self.log(f"[INFO] Loaded input TSV: {self.input_tsv} (rows={len(in_df)})")
        self.log(f"[INFO] Loaded mapping file: {self.mapping_file} (pairs={len(mapping)})")
        stream_stats_only = ("sequence" not in base_col_set) and bool(self.ref_list)
        if stream_stats_only:
            fasta_stats = self.parse_fasta_stats()
            fasta_map = None
            self.log(f"[INFO] Streamed FASTA base stats: {self.input_fasta} (entries={len(fasta_stats)})")
        else:
            fasta_map = self.parse_fasta()
            fasta_stats = None
            self.log(f"[INFO] Loaded FASTA: {self.input_fasta} (entries={len(fasta_map)})")

        if self.key not in in_df.columns:
            self.log(f"[WARN] Key column '{self.key}' not found in input TSV. No rows will be appended.")
            base_df.to_csv(self.output, sep="\t", index=False)
            self.log(f"[INFO] Merged matrix written: {self.output}")
            return

        # Duplicates handling in input tsv (keep first)
        key_series = in_df[self.key].fillna("").astype(str)
        missing_key_mask = key_series.eq("")
        dup_mask = in_df.duplicated(subset=[self.key], keep="first")

        dup_report = pd.DataFrame()
        if missing_key_mask.any():
            r1 = in_df[missing_key_mask].copy()
            r1.insert(0, "row_key", "")
            r1.insert(0, "reason", "missing-key")
            dup_report = pd.concat([dup_report, r1], ignore_index=True)
        if dup_mask.any():
            r2 = in_df[dup_mask].copy()
            r2.insert(0, "row_key", r2[self.key])
            r2.insert(0, "reason", "duplicate")
            dup_report = pd.concat([dup_report, r2], ignore_index=True)

        if not dup_report.empty:
            dup_file = self.write_duplicates_report_df(dup_report)
            self.log(f"[WARN] Duplicates/missing-key in input: {len(dup_report)} (report: {dup_file})")
        else:
            self.log("[INFO] No duplicates in input.")

        keep_mask = (~missing_key_mask) & (~dup_mask)
        norm_df = in_df.loc[keep_mask].copy()

        if mapping:
            norm_df = norm_df.rename(columns=mapping)

        # Namespace every column that arrived raw from the external source and was
        # not deliberately mapped onto a canonical name. Without this, a vendor
        # field that merely differs in case from a canonical one - GISAID ships
        # `Segment` against our `segment` - reaches SQLite, which treats column
        # identifiers case-insensitively and aborts the build with
        # "duplicate column name: Segment". Prefixing makes that unrepresentable
        # and, as a side effect, guarantees every identifier is legal SQL.
        norm_df = self.namespace_external_columns(
            norm_df, set(mapping.values()), canonical=set(base_df.columns))

        if self.columns_collide_ignoring_case(norm_df.columns):
            dup_cols = self.case_insensitive_duplicates(norm_df.columns)
            self.log(f"[WARN] Duplicate column labels after mapping: {dup_cols}. Collapsing duplicates.")
            norm_df = self.collapse_duplicate_columns(norm_df)

        # copy key to gB fields
        norm_df["primary_accession"] = norm_df[self.key]
        norm_df["accession_version"] = norm_df[self.key]
        norm_df["locus"] = norm_df[self.key]
        norm_df["gi_number"] = norm_df[self.key]

        # import sequences or lightweight base counts from FASTA
        if stream_stats_only:
            key_list = norm_df[self.key].tolist()
            zero_stat = (0, 0, 0, 0, 0, 0, 0)
            matched_stats = [fasta_stats.get(k) for k in key_list]
            missing_seq = pd.Series([st is None or st[0] == 0 for st in matched_stats], index=norm_df.index)
            norm_df.loc[missing_seq, "exclusion"] = "missing fasta sequence"
            norm_df["real_length"] = [st[1] if st else 0 for st in matched_stats]
        else:
            norm_df["sequence"] = norm_df[self.key].map(fasta_map).fillna("")
            norm_df.loc[norm_df["sequence"] == "", "exclusion"] = "missing fasta sequence"
            norm_df["real_length"] = norm_df["sequence"].apply(self.count_real_length)
            matched_stats = None
        
        norm_df["dataset_source"] = self.dataset_source

        ref_dict = load_reference_file_dict(self.ref_list) if (self.ref_list and os.path.isfile(self.ref_list)) else {}

        # Populate canonical gB_matrix columns when present in base_df so external/GISAID
        # references and queries behave identically to GenBank rows downstream.
        if "accession_type" in base_col_set:
            existing_type = norm_df["accession_type"].fillna("").astype(str).str.strip() if "accession_type" in norm_df.columns else pd.Series("", index=norm_df.index)
            inferred_type = norm_df["primary_accession"].map(lambda acc: ref_dict.get(str(acc).strip(), "query"))
            norm_df["accession_type"] = existing_type.where(existing_type != "", inferred_type)

        if "exclusion_status" in base_col_set:
            excl_text = norm_df["exclusion"].fillna("").astype(str).str.strip() if "exclusion" in norm_df.columns else pd.Series("", index=norm_df.index)
            is_excl_ref = (
                norm_df["accession_type"].fillna("").astype(str).str.strip().str.lower().eq("exclusion_list")
                if "accession_type" in norm_df.columns else pd.Series(False, index=norm_df.index)
            )
            norm_df["exclusion_status"] = (excl_text.ne("") | is_excl_ref).map(lambda flag: "1" if flag else "0")
            if "exclusion_criteria" in base_col_set:
                existing_crit = norm_df["exclusion_criteria"].fillna("").astype(str) if "exclusion_criteria" in norm_df.columns else pd.Series("", index=norm_df.index)
                norm_df["exclusion_criteria"] = existing_crit.where(existing_crit != "", excl_text)

        if "division" in base_col_set:
            existing_div = norm_df["division"].fillna("").astype(str).str.strip() if "division" in norm_df.columns else pd.Series("", index=norm_df.index)
            norm_df["division"] = existing_div.where(existing_div != "", "VRL")

        if stream_stats_only:
            if "length" in base_col_set:
                norm_df["length"] = [str(st[0]) if st else "0" for st in matched_stats]
            for idx_nt, base_nt in ((2, "a"), (3, "t"), (4, "g"), (5, "c"), (6, "n")):
                if base_nt in base_col_set:
                    norm_df[base_nt] = [str(st[idx_nt]) if st else "0" for st in matched_stats]
        else:
            if "length" in base_col_set:
                norm_df["length"] = norm_df["sequence"].map(lambda s: str(len(s)) if isinstance(s, str) else "0")
            for base_nt in ("a", "t", "g", "c", "n"):
                if base_nt in base_col_set:
                    norm_df[base_nt] = norm_df["sequence"].map(lambda s, ch=base_nt: str(s.lower().count(ch)) if isinstance(s, str) else "0")

        if "segment" in base_col_set and "segment_declared" not in norm_df.columns and self.key != "Segment_Id":
            base_segs = {str(v).strip() for v in base_df["segment"].tolist() if str(v).strip()}
            if base_segs == {"1"}:
                existing_seg = norm_df["segment"].fillna("").astype(str).str.strip() if "segment" in norm_df.columns else pd.Series("", index=norm_df.index)
                norm_df["segment"] = existing_seg.where(existing_seg != "", "1")

        if "collection_date" in norm_df.columns and {"collection_day", "collection_mon", "collection_year"} <= base_col_set:
            parsed_dates = norm_df["collection_date"].map(lambda d: split_date_components(str(d or "")))
            for part in ("day", "mon", "year"):
                col_name = f"collection_{part}"
                norm_df[col_name] = parsed_dates.map(lambda parts_dict, k=part: str(parts_dict.get(k, "")))

        base_cols = list(base_df.columns)

        if self.drop_unmapped:
            mapping_targets = set(mapping.values())
            required_cols = {'primary_accession', 'accession_version', 'locus', 'sequence', 'dataset_source', 'exclusion'}
            keep_cols = list(set(base_cols) | mapping_targets | required_cols | {self.key})
            ordered_keep = [c for c in base_cols if c in keep_cols] + [c for c in norm_df.columns if c in keep_cols and c not in base_cols]
            norm_df = norm_df[[c for c in ordered_keep if c in norm_df.columns]]
            extras = [c for c in norm_df.columns if c not in base_cols]
            final_cols = base_cols + extras
        else:
            extras = [c for c in norm_df.columns if c not in base_cols]
            final_cols = base_cols + extras

        norm_df = norm_df.reindex(columns=final_cols, fill_value="")
        merged_df = pd.concat([base_df.reindex(columns=final_cols, fill_value=""), norm_df], ignore_index=True)

        if self.require_refs and ref_dict:
            present_valid = set(
                merged_df.loc[
                    merged_df.get("exclusion_status", pd.Series("0", index=merged_df.index)).astype(str).str.strip() != "1",
                    "primary_accession",
                ].astype(str).str.strip()
            ) if "exclusion_status" in merged_df.columns else set(merged_df["primary_accession"].astype(str).str.strip())
            missing_critical = sorted(
                acc for acc, rtype in ref_dict.items()
                if str(rtype).strip().lower() != "exclusion_list" and acc not in present_valid
            )
            if missing_critical:
                raise ValueError(
                    f"Missing CRITICAL reference accessions after merging GenBank + {self.dataset_source}: "
                    f"{', '.join(missing_critical[:20])}"
                    + (f" (+{len(missing_critical) - 20} more)" if len(missing_critical) > 20 else "")
                )

        # Drop raw sequence column if base_df did not carry it (sequences are combined via CAT_FASTA)
        if "sequence" not in base_col_set and "sequence" in merged_df.columns and self.ref_list:
            merged_df = merged_df.drop(columns=["sequence"])

        merged_df.to_csv(self.output, sep="\t", index=False)
        self.log(f"[INFO] Merged matrix written: {self.output}")


if __name__ == "__main__":
    parser = ArgumentParser(description='Normalize external TSV+FASTA to gB_matrix format and append to base matrix.')
    parser.add_argument('-g', '--gb_matrix', help='Base gB_matrix file', default='gB_matrix.tsv')
    parser.add_argument('-t', '--input_tsv', help='Input TSV to transform', default='./tmp/gisaid-data/metadata.tsv')
    parser.add_argument('-f', '--input_fasta', help='Input FASTA with sequences', default='./tmp/gisaid-data/all_nuc.fas')
    parser.add_argument('-m', '--mapping_file', help='Column mapping file (source_col<TAB>target_col)',
                        default='./generic-influenza/column_mapping.tsv')
    parser.add_argument('-k', '--key', help='Key column name in input TSV', default='Segment_Id')
    parser.add_argument('-s', '--dataset_source', 
                        help='Source name value to annotate in the "dataset_source" column in "gB_matrix_merged.tsv"', default='external')
    parser.add_argument('-o', '--output', help='Full output path including filename',
                        default='./tmp/gisaid-data/gB_matrix_merged.tsv')
    parser.add_argument('--log_file', help='Path to main process log (default: gisaid_processing.log in CWD)', default=None)
    parser.add_argument('--prefer', choices=['input', 'base'], default='input',
                        help='(accepted for compatibility) preferred side when upsert mode is enabled')
    parser.add_argument('--drop_unmapped', action='store_true',
                        help='If set, drop columns from the new dataset that are not in base matrix, mapping targets, or required columns')
    parser.add_argument('-r', '--ref_list', default=None,
                        help='Optional reference list TSV to classify merged master/reference rows')
    parser.add_argument('--require_refs', action='store_true',
                        help='Fail if any master/reference accession in --ref_list is missing after merge')

    args = parser.parse_args()

    job = NormalizeAndMerge(
        gb_matrix=args.gb_matrix,
        input_tsv=args.input_tsv,
        input_fasta=args.input_fasta,
        mapping_file=args.mapping_file,
        output=args.output,
        key=args.key,
        dataset_source=args.dataset_source,
        prefer=args.prefer,
        drop_unmapped=args.drop_unmapped,
        log_file=args.log_file,
        ref_list=args.ref_list,
        require_refs=args.require_refs,
    )
    job.process()
