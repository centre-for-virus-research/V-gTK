import os
import csv
import re
import sqlite3
import sys
from os.path import join, dirname, abspath
from argparse import ArgumentParser
from datetime import datetime
from collections import OrderedDict, defaultdict
import pandas as pd

_SCRIPTS_DIR = dirname(abspath(__file__))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

from accession_utils import normalise_accession
from cross_db_dedup import (
    candidate_winner_sort_key,
    coalesce_winner_metadata,
    count_acgt,
    extract_canonical_isolate,
    extract_insdc_links,
    extract_passage_class,
    has_coverage_parity,
    parse_fasta_file,
    passages_compatible,
    resolve_row_segment_hint,
    unaligned_sequences_agree,
)

"""
Transform an external TSV + FASTA into the gB_matrix-compatible format and merge with base gB_matrix.
- Key in input TSV: configurable with --key, default 'Segment_Id'
- Copies the key into: primary_accession, accession_version, locus columns
- Adds sequences from FASTA (header must match the key exactly)
- If a TSV row has no matching FASTA, sets exclusion = "missing fasta sequence"
- Duplicates in input by key: keeps the first, writes the others to merge_matrix_duplicates.tsv and log count
- Column mapping file: two-column TSV (source_col<TAB>target_col); mapped columns are copied/renamed
- Cross-source deduplication (enabled by default, --no_dedup_cross_source to disable):
  detects dual uploads deposited in both GenBank (base / --update_db) and the external
  source (GISAID / Pathoplexus) via explicit INSDC cross-reference columns or matching
  canonical isolate + compatible passage + 0-mismatch sequence agreement, coalesces
  richer metadata onto the preferred record (default: base/GenBank), and retains the
  retired external row as an exclusion_status='1' tombstone so incremental updates never
  re-ingest it while downstream BLAST/alignment/tree stages skip it automatically.
"""


class NormalizeAndMerge:
    def __init__(
        self,
        gb_matrix,
        input_tsv,
        input_fasta,
        mapping_file,
        output,
        key="Segment_Id",
        dataset_source="external",
        prefer="base",
        drop_unmapped=False,
        log_file=None,
        base_fasta=None,
        update_db=None,
        dedup_cross_source=True,
        dedup_fasta_out=None,
        dedup_report=None,
    ):
        self.gb_matrix = gb_matrix
        self.input_tsv = input_tsv
        self.input_fasta = input_fasta
        self.mapping_file = mapping_file
        self.output = output
        self.key = key
        self.dataset_source = dataset_source
        self.prefer = prefer
        self.drop_unmapped = drop_unmapped
        self.base_fasta = base_fasta
        self.update_db = update_db
        self.dedup_cross_source = dedup_cross_source
        self.dedup_fasta_out = dedup_fasta_out
        self.dedup_report = dedup_report

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

    CROSS_DB_REPORT_COLUMNS = [
        "winner_accession",
        "retired_accession",
        "match_tier",
        "mismatches",
        "segment",
        "canonical_isolate",
        "winner_origin",
    ]

    def write_cross_db_duplicates_report_df(self, dup_df: pd.DataFrame):
        if self.dedup_report:
            dup_file = self.dedup_report
            out_dir = dirname(abspath(dup_file)) or "."
            os.makedirs(out_dir, exist_ok=True)
            if dup_df is None or dup_df.empty:
                pd.DataFrame(columns=self.CROSS_DB_REPORT_COLUMNS).to_csv(dup_file, sep="\t", index=False)
            else:
                dup_df.to_csv(dup_file, sep="\t", index=False)
            return dup_file
        if dup_df is None or dup_df.empty:
            return None
        dup_file = join(dirname(abspath(self.output)) or ".", 'cross_db_duplicates.tsv')
        dup_df.to_csv(dup_file, sep="\t", index=False)
        return dup_file

    def write_dedup_fasta(self, fasta_map: dict, retired_keys: set):
        if not self.dedup_fasta_out:
            return None
        out_dir = dirname(abspath(self.dedup_fasta_out)) or "."
        os.makedirs(out_dir, exist_ok=True)
        written = 0
        with open(self.dedup_fasta_out, "w", encoding="utf-8") as fh:
            for hdr, seq in fasta_map.items():
                if hdr in retired_keys:
                    continue
                fh.write(f">{hdr}\n{seq}\n")
                written += 1
        self.log(f"[INFO] Wrote deduplicated external FASTA: {self.dedup_fasta_out} (records={written}, retired_omitted={len(retired_keys)})")
        return self.dedup_fasta_out

    #: Columns this class writes itself. They are canonical by construction and
    #: must never be namespaced away.
    PIPELINE_OWNED_COLUMNS = frozenset({
        "primary_accession", "accession_version", "locus", "gi_number",
        "sequence", "real_length", "dataset_source", "exclusion",
        "exclusion_status", "exclusion_criteria", "canonical_accession", "alt_accessions",
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

    def _load_update_db_genbank_index(self):
        """Load kept rows from ``self.update_db`` (both GenBank/INSDC and external portals)."""
        records = []
        if not self.update_db or not os.path.isfile(self.update_db):
            return records

        conn = sqlite3.connect(f"file:{self.update_db}?mode=ro", uri=True)
        try:
            tables = {
                r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
            }
            if "meta_data" not in tables:
                return records
            meta_cols = [r[1] for r in conn.execute("PRAGMA table_info(meta_data)")]
            if "primary_accession" not in meta_cols:
                return records

            wanted = [
                c for c in (
                    "primary_accession", "segment", "segment_validated", "strain",
                    "Parsed_strain", "isolate", "definition", "collection_date",
                    "collection_year", "Passage_History", "isolation_source",
                    "real_length", "dataset_source", "exclusion_status", "accession_type",
                )
                if c in meta_cols
            ]
            seq_by_acc = {}
            if "sequences" in tables:
                seq_cols = [r[1] for r in conn.execute("PRAGMA table_info(sequences)")]
                hdr_col = "header" if "header" in seq_cols else ("primary_accession" if "primary_accession" in seq_cols else None)
                if hdr_col and "sequence" in seq_cols:
                    for hdr, seq in conn.execute(f'SELECT "{hdr_col}", sequence FROM sequences WHERE sequence IS NOT NULL'):
                        if hdr:
                            bare = normalise_accession(str(hdr).strip())
                            if bare:
                                seq_by_acc[bare.upper()] = str(seq)

            quoted = ", ".join(f'"{c}"' for c in wanted)
            for row in conn.execute(f"SELECT {quoted} FROM meta_data"):
                rec = {col: ("" if val is None else str(val)) for col, val in zip(wanted, row)}
                if str(rec.get("exclusion_status") or "").strip() == "1":
                    continue
                src = str(rec.get("dataset_source") or "genbank").strip().lower()
                acc = normalise_accession(rec.get("primary_accession"))
                if not acc:
                    continue
                acc_up = acc.upper()
                seq = seq_by_acc.get(acc_up, "")
                try:
                    rlen = int(float(rec.get("real_length") or 0))
                except ValueError:
                    rlen = 0
                if rlen <= 0 and seq:
                    rlen = count_acgt(seq)
                records.append({
                    "df_idx": None,
                    "accession": acc,
                    "accession_upper": acc_up,
                    "segment": resolve_row_segment_hint(rec),
                    "isolate": extract_canonical_isolate(rec),
                    "passage": extract_passage_class(rec),
                    "seq": seq,
                    "real_length": rlen,
                    "protected": True,  # already committed in update_db
                    "origin": "update_db",
                    "source": src or "genbank",
                    "collection_date": rec.get("collection_date", ""),
                })
        finally:
            conn.close()
        return records

    def deduplicate_cross_source(self, base_df: pd.DataFrame, norm_df: pd.DataFrame):
        """Reconcile cross-source AND intra-source 0-mm same-strain duplicates across
        ``base_df``, ``norm_df``, and ``update_db``.

        Returns ``(base_df, norm_df, retired_input_keys, report_df)``.
        """
        if not self.dedup_cross_source or (norm_df.empty and base_df.empty):
            return base_df, norm_df, set(), pd.DataFrame()

        base_fasta_map = parse_fasta_file(self.base_fasta) if self.base_fasta else {}
        if self.base_fasta:
            self.log(f"[INFO] Loaded base FASTA for cross-source deduplication: {self.base_fasta} (entries={len(base_fasta_map)})")

        by_acc = defaultdict(list)
        by_isolate = defaultdict(list)
        ext_update_by_isolate = defaultdict(list)

        update_db_entries = self._load_update_db_genbank_index()
        if update_db_entries:
            n_gb_upd = 0
            n_ext_upd = 0
            for entry in update_db_entries:
                src = str(entry.get("source") or "genbank").lower()
                if src in ("", "genbank", "ncbi", "insdc", "ena", "ddbj"):
                    n_gb_upd += 1
                    by_acc[entry["accession_upper"]].append(entry)
                    if entry["isolate"]:
                        by_isolate[entry["isolate"]].append(entry)
                else:
                    n_ext_upd += 1
                    if entry["isolate"]:
                        ext_update_by_isolate[entry["isolate"]].append(entry)
            self.log(
                f"[INFO] Indexed {n_gb_upd} GenBank + {n_ext_upd} external rows from update_db: {self.update_db}"
            )

        retired_input_keys = set()
        report_rows = []
        base_updates_by_idx = defaultdict(dict)
        norm_updates_by_pos = defaultdict(dict)

        # Build base (GenBank) candidate list and perform intra-GenBank 0-mm strain dedup
        base_records = base_df.fillna("").to_dict(orient="records")
        raw_base_entries = []
        for idx, rec in enumerate(base_records):
            if str(rec.get("exclusion_status") or "").strip() == "1":
                continue
            if str(rec.get("exclusion") or "").strip():
                continue
            raw_acc = rec.get("primary_accession") or rec.get("gi_number") or rec.get("accession_version")
            acc = normalise_accession(raw_acc)
            if not acc:
                continue
            acc_up = acc.upper()
            seq = (
                base_fasta_map.get(str(raw_acc).strip())
                or base_fasta_map.get(acc)
                or base_fasta_map.get(acc_up)
                or str(rec.get("sequence") or "")
            )
            try:
                rlen = int(float(rec.get("real_length") or 0))
            except ValueError:
                rlen = 0
            if rlen <= 0 and seq:
                rlen = count_acgt(seq)
            acc_type = str(rec.get("accession_type") or "").strip().lower()
            raw_base_entries.append({
                "df_idx": idx,
                "accession": str(rec.get("primary_accession") or acc).strip(),
                "accession_upper": acc_up,
                "segment": resolve_row_segment_hint(rec),
                "isolate": extract_canonical_isolate(rec),
                "passage": extract_passage_class(rec),
                "seq": seq,
                "real_length": rlen,
                "protected": acc_type in ("master", "reference"),
                "origin": "base_matrix",
                "source": "genbank",
                "collection_date": str(rec.get("collection_date") or ""),
            })

        raw_base_entries.sort(key=lambda e: candidate_winner_sort_key(e, prefer_source="genbank"))
        base_redirect_by_acc = {}
        for entry in raw_base_entries:
            idx = entry["df_idx"]
            matched_gb = None
            if not entry["protected"] and entry["isolate"] and entry["seq"]:
                for cand in by_isolate.get(entry["isolate"], ()):
                    if not cand["seq"]:
                        continue
                    if entry["segment"] and cand["segment"] and entry["segment"] != cand["segment"]:
                        continue
                    if not passages_compatible(cand["passage"], entry["passage"]):
                        continue
                    ok, _, _ = unaligned_sequences_agree(cand["seq"], entry["seq"], max_mm=0)
                    if ok:
                        matched_gb = cand
                        break
            if matched_gb is not None:
                winner_acc = matched_gb["accession"]
                base_redirect_by_acc[entry["accession_upper"]] = matched_gb
                if matched_gb["df_idx"] is not None:
                    w_idx = matched_gb["df_idx"]
                    w_row = dict(base_records[w_idx])
                    w_row.update(base_updates_by_idx[w_idx])
                    w_row.setdefault("alt_accessions", "")
                    upd = coalesce_winner_metadata(
                        w_row, base_records[idx], retired_source="genbank", update_split_date_cols=True
                    )
                    base_updates_by_idx[w_idx].update(upd)
                reason = f"intra_db_duplicate_of:{winner_acc}"
                base_updates_by_idx[idx].update({
                    "exclusion_status": "1",
                    "exclusion_criteria": reason,
                    "exclusion": reason,
                    "canonical_accession": winner_acc,
                })
                report_rows.append({
                    "winner_accession": winner_acc,
                    "retired_accession": entry["accession"],
                    "match_tier": "intra_db_same_strain",
                    "mismatches": "0",
                    "segment": entry["segment"] or matched_gb["segment"],
                    "canonical_isolate": entry["isolate"],
                    "winner_origin": matched_gb["origin"],
                })
            else:
                by_acc[entry["accession_upper"]].append(entry)
                if entry["isolate"]:
                    by_isolate[entry["isolate"]].append(entry)
        for retired_up, canon_entry in base_redirect_by_acc.items():
            if retired_up not in by_acc:
                by_acc[retired_up].append(canon_entry)

        norm_records = norm_df.fillna("").to_dict(orient="records")
        norm_indices = list(norm_df.index)
        kept_ext_by_isolate = defaultdict(list)
        for iso, entries in ext_update_by_isolate.items():
            kept_ext_by_isolate[iso].extend(entries)

        # Evaluate incoming external records in deterministic quality order so longest wins intra-GISAID ties
        ext_order = []
        for pos, ext_rec in enumerate(norm_records):
            ext_acc = str(ext_rec.get("primary_accession") or "").strip()
            ext_seq = str(ext_rec.get("sequence") or "")
            try:
                ext_len = int(float(ext_rec.get("real_length") or 0))
            except ValueError:
                ext_len = count_acgt(ext_seq)
            ext_acc_type = str(ext_rec.get("accession_type") or "").strip().lower()
            sort_stub = {
                "protected": ext_acc_type in ("master", "reference"),
                "source": self.dataset_source,
                "origin": "input",
                "real_length": ext_len,
                "collection_date": str(ext_rec.get("collection_date") or ""),
                "accession": ext_acc,
            }
            ext_order.append((candidate_winner_sort_key(sort_stub, prefer_source="genbank"), pos, ext_rec, ext_acc, ext_seq, ext_len, ext_acc_type))
        ext_order.sort(key=lambda t: t[0])

        for _, pos, ext_rec, ext_acc, ext_seq, ext_len, ext_acc_type in ext_order:
            if str(ext_rec.get("exclusion_status") or "").strip() == "1":
                continue
            if str(ext_rec.get("exclusion") or "").strip():
                continue
            if not ext_acc or ext_len <= 0:
                continue

            ext_seg = resolve_row_segment_hint(ext_rec)
            ext_isolate = extract_canonical_isolate(ext_rec)
            ext_passage = extract_passage_class(ext_rec)
            if ext_acc_type in ("master", "reference"):
                if ext_isolate:
                    kept_ext_by_isolate[ext_isolate].append({
                        "pos": pos, "accession": ext_acc, "segment": ext_seg,
                        "isolate": ext_isolate, "passage": ext_passage,
                        "seq": ext_seq, "real_length": ext_len,
                        "protected": True, "origin": "input",
                    })
                continue

            insdc_links = extract_insdc_links(ext_rec, segment_hint=ext_seg)

            matched_cand = None
            match_tier = ""
            match_mm = 0

            # Tier 1: Explicit INSDC upload cross-reference verified by sequence (or strain+length if base_fasta absent)
            for link_acc in insdc_links:
                for cand in by_acc.get(link_acc, ()):
                    if ext_seg and cand["segment"] and ext_seg != cand["segment"]:
                        continue
                    prefer_base = (self.prefer != "input") or cand["protected"]
                    w_len = cand["real_length"] if prefer_base else ext_len
                    r_len = ext_len if prefer_base else cand["real_length"]
                    if cand["seq"] and ext_seq:
                        allow_mm = (
                            2
                            if (ext_isolate and cand["isolate"] and ext_isolate == cand["isolate"] and passages_compatible(cand["passage"], ext_passage))
                            else 0
                        )
                        w_seq = cand["seq"] if prefer_base else ext_seq
                        r_seq = ext_seq if prefer_base else cand["seq"]
                        ok, mm, _ = unaligned_sequences_agree(w_seq, r_seq, max_mm=allow_mm)
                        if ok:
                            matched_cand = cand
                            match_tier = "insdc_link_seq_verified"
                            match_mm = mm
                            break
                    else:
                        if (
                            ext_isolate
                            and cand["isolate"]
                            and ext_isolate == cand["isolate"]
                            and passages_compatible(cand["passage"], ext_passage)
                            and has_coverage_parity(w_len, r_len)
                        ):
                            matched_cand = cand
                            match_tier = "insdc_link_strain_length_verified"
                            match_mm = 0
                            break
                if matched_cand is not None:
                    break

            # Tier 2: Unlinked dual upload — same canonical isolate + compatible passage + 0-mismatch sequence
            if matched_cand is None and ext_isolate and ext_seq:
                for cand in by_isolate.get(ext_isolate, ()):
                    if not cand["seq"]:
                        continue
                    if ext_seg and cand["segment"] and ext_seg != cand["segment"]:
                        continue
                    if not passages_compatible(cand["passage"], ext_passage):
                        continue
                    prefer_base = (self.prefer != "input") or cand["protected"]
                    w_seq = cand["seq"] if prefer_base else ext_seq
                    r_seq = ext_seq if prefer_base else cand["seq"]
                    ok, mm, _ = unaligned_sequences_agree(w_seq, r_seq, max_mm=0)
                    if ok:
                        matched_cand = cand
                        match_tier = "strain_and_zero_mm_sequence"
                        match_mm = mm
                        break

            # Tier 3: Intra-source duplicate (same external source + same isolate + compatible passage + 0-mm sequence)
            matched_ext_peer = None
            if matched_cand is None and ext_isolate and ext_seq:
                for peer in kept_ext_by_isolate.get(ext_isolate, ()):
                    if not peer.get("seq"):
                        continue
                    if ext_seg and peer.get("segment") and ext_seg != peer.get("segment"):
                        continue
                    if not passages_compatible(peer.get("passage", "unknown"), ext_passage):
                        continue
                    ok, _, _ = unaligned_sequences_agree(peer["seq"], ext_seq, max_mm=0)
                    if ok:
                        matched_ext_peer = peer
                        break

            if matched_cand is None and matched_ext_peer is None:
                if ext_isolate:
                    kept_ext_by_isolate[ext_isolate].append({
                        "pos": pos, "accession": ext_acc, "segment": ext_seg,
                        "isolate": ext_isolate, "passage": ext_passage,
                        "seq": ext_seq, "real_length": ext_len,
                        "protected": False, "origin": "input",
                    })
                continue

            if matched_ext_peer is not None:
                winner_acc = matched_ext_peer["accession"]
                if matched_ext_peer.get("pos") is not None:
                    w_pos = matched_ext_peer["pos"]
                    w_row = dict(norm_records[w_pos])
                    w_row.update(norm_updates_by_pos[w_pos])
                    w_row.setdefault("alt_accessions", "")
                    upd = coalesce_winner_metadata(
                        w_row, ext_rec, retired_source=self.dataset_source, update_split_date_cols=True
                    )
                    norm_updates_by_pos[w_pos].update(upd)
                reason = f"intra_db_duplicate_of:{winner_acc}"
                norm_updates_by_pos[pos].update({
                    "exclusion_status": "1",
                    "exclusion_criteria": reason,
                    "exclusion": reason,
                    "canonical_accession": winner_acc,
                })
                retired_input_keys.add(str(ext_rec.get(self.key) or ext_acc))
                report_rows.append({
                    "winner_accession": winner_acc,
                    "retired_accession": ext_acc,
                    "match_tier": "intra_db_same_strain",
                    "mismatches": "0",
                    "segment": ext_seg or matched_ext_peer.get("segment", ""),
                    "canonical_isolate": ext_isolate,
                    "winner_origin": matched_ext_peer.get("origin", self.dataset_source),
                })
                continue

            prefer_base = (self.prefer != "input") or matched_cand["protected"] or matched_cand["df_idx"] is None
            if prefer_base:
                winner_acc = matched_cand["accession"]
                retired_acc = ext_acc
                if matched_cand["df_idx"] is not None:
                    b_idx = matched_cand["df_idx"]
                    b_row = dict(base_records[b_idx])
                    b_row.update(base_updates_by_idx[b_idx])
                    b_row.setdefault("alt_accessions", "")
                    updates = coalesce_winner_metadata(
                        b_row, ext_rec, retired_source=self.dataset_source, update_split_date_cols=True
                    )
                    base_updates_by_idx[b_idx].update(updates)
                reason = f"cross_db_duplicate_of:{winner_acc}"
                norm_updates_by_pos[pos].update({
                    "exclusion_status": "1",
                    "exclusion_criteria": reason,
                    "exclusion": reason,
                    "canonical_accession": winner_acc,
                })
                retired_input_keys.add(str(ext_rec.get(self.key) or ext_acc))
            else:
                winner_acc = ext_acc
                retired_acc = matched_cand["accession"]
                b_idx = matched_cand["df_idx"]
                b_row = dict(base_records[b_idx])
                e_row = dict(ext_rec)
                e_row.update(norm_updates_by_pos[pos])
                e_row.setdefault("alt_accessions", "")
                updates = coalesce_winner_metadata(
                    e_row, b_row, retired_source="genbank", update_split_date_cols=True
                )
                norm_updates_by_pos[pos].update(updates)
                reason = f"cross_db_duplicate_of:{winner_acc}"
                base_updates_by_idx[b_idx].update({
                    "exclusion_status": "1",
                    "exclusion_criteria": reason,
                    "exclusion": reason,
                    "canonical_accession": winner_acc,
                })
                if ext_isolate:
                    kept_ext_by_isolate[ext_isolate].append({
                        "pos": pos, "accession": ext_acc, "segment": ext_seg,
                        "isolate": ext_isolate, "passage": ext_passage,
                        "seq": ext_seq, "real_length": ext_len,
                        "protected": False, "origin": "input",
                    })

            report_rows.append({
                "winner_accession": winner_acc,
                "retired_accession": retired_acc,
                "match_tier": match_tier,
                "mismatches": str(match_mm),
                "segment": ext_seg or matched_cand["segment"],
                "canonical_isolate": ext_isolate or matched_cand["isolate"],
                "winner_origin": matched_cand["origin"] if prefer_base else self.dataset_source,
            })

        if not report_rows:
            rep_path = self.write_cross_db_duplicates_report_df(pd.DataFrame())
            self.log(
                f"[INFO] Cross-source deduplication: 0 cross-database duplicates found"
                + (f" (report={rep_path})." if rep_path else ".")
            )
            return base_df, norm_df, retired_input_keys, pd.DataFrame()

        base_df = base_df.copy()
        norm_df = norm_df.copy()

        for b_idx, upd in base_updates_by_idx.items():
            row_label = base_df.index[b_idx]
            for col, val in upd.items():
                if col not in base_df.columns:
                    base_df[col] = ""
                base_df.at[row_label, col] = val

        for pos, upd in norm_updates_by_pos.items():
            row_label = norm_indices[pos]
            for col, val in upd.items():
                if col not in norm_df.columns:
                    norm_df[col] = ""
                norm_df.at[row_label, col] = val

        report_df = pd.DataFrame(report_rows)
        rep_path = self.write_cross_db_duplicates_report_df(report_df)
        self.log(
            f"[INFO] Cross-source deduplication retired {len(report_rows)} dual-uploaded record(s) "
            f"(prefer={self.prefer}, report={rep_path})"
        )
        return base_df, norm_df, retired_input_keys, report_df

    def process(self):
        base_df = self.read_base_matrix_df()
        orig_base_cols = list(base_df.columns)
        self.log(f"[INFO] Loaded base matrix: {self.gb_matrix} (rows={len(base_df)})")
        base_df["dataset_source"] = "genbank" # define source name
        mapping = self.read_mapping()
        in_df = self.read_input_tsv_df()
        self.log(f"[INFO] Loaded input TSV: {self.input_tsv} (rows={len(in_df)})")
        self.log(f"[INFO] Loaded mapping file: {self.mapping_file} (pairs={len(mapping)})")
        fasta_map = self.parse_fasta()
        self.log(f"[INFO] Loaded FASTA: {self.input_fasta} (entries={len(fasta_map)})")

        if self.key not in in_df.columns:
            self.log(f"[WARN] Key column '{self.key}' not found in input TSV. No rows will be appended.")
            self.write_dedup_fasta(fasta_map, set())
            self.write_cross_db_duplicates_report_df(pd.DataFrame())
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
        key_series = norm_df[self.key].copy()

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
        norm_df["primary_accession"] = key_series
        norm_df["accession_version"] = key_series
        norm_df["locus"] = key_series
        norm_df["gi_number"] = key_series

        # import sequences from FASTA to column 'sequence'
        norm_df["sequence"] = key_series.map(fasta_map).fillna("")
        if "exclusion" not in norm_df.columns:
            norm_df["exclusion"] = ""
        else:
            norm_df["exclusion"] = norm_df["exclusion"].fillna("")
        norm_df.loc[norm_df["sequence"] == "", "exclusion"] = "missing fasta sequence"
        norm_df["real_length"] = norm_df["sequence"].apply(self.count_real_length)
        
        norm_df["dataset_source"] = self.dataset_source

        base_df, norm_df, retired_input_keys, cross_dup_df = self.deduplicate_cross_source(base_df, norm_df)
        if cross_dup_df.empty and self.dedup_report and not os.path.exists(self.dedup_report):
            self.write_cross_db_duplicates_report_df(cross_dup_df)
        self.write_dedup_fasta(fasta_map, retired_input_keys)

        # Preserve original base column ordering at the front of the output table.
        base_cols = orig_base_cols + [c for c in base_df.columns if c not in orig_base_cols]

        if self.drop_unmapped:
            mapping_targets = set(mapping.values())
            required_cols = {
                'primary_accession', 'accession_version', 'locus', 'sequence',
                'dataset_source', 'exclusion', 'exclusion_status', 'exclusion_criteria',
                'canonical_accession', 'alt_accessions',
            }
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

        merged_df.to_csv(self.output, sep="\t", index=False)
        self.log(f"[INFO] Merged matrix written: {self.output}")


if __name__ == "__main__":
    parser = ArgumentParser(description='Normalize external TSV+FASTA to gB_matrix format and merge with base matrix.')
    parser.add_argument('-g', '--gb_matrix', help='Base gB_matrix file', default='gB_matrix.tsv')
    parser.add_argument('-t', '--input_tsv', help='Input TSV to transform', default='./tmp/gisaid-data/metadata.tsv')
    parser.add_argument('-f', '--input_fasta', help='Input FASTA with sequences', default='./tmp/gisaid-data/all_nuc.fas')
    parser.add_argument('--base_fasta', help='Optional FASTA for base gB_matrix records (e.g. GenBankParser sequences.fa)', default=None)
    parser.add_argument('--update_db', help='Optional existing SQLite DB for cross-batch duplicate detection in --update_db mode', default=None)
    parser.add_argument('--dedup_fasta_out', help='Optional output FASTA omitting retired external duplicate sequences', default=None)
    parser.add_argument('--dedup_report', help='Optional output TSV path for cross-source duplicate audit report', default=None)
    parser.add_argument('-m', '--mapping_file', help='Column mapping file (source_col<TAB>target_col)',
                        default='./generic-influenza/column_mapping.tsv')
    parser.add_argument('-k', '--key', help='Key column name in input TSV', default='Segment_Id')
    parser.add_argument('-s', '--dataset_source', 
                        help='Source name value to annotate in the "dataset_source" column in "gB_matrix_merged.tsv"', default='external')
    parser.add_argument('-o', '--output', help='Full output path including filename',
                        default='./tmp/gisaid-data/gB_matrix_merged.tsv')
    parser.add_argument('--log_file', help='Path to main process log (default: gisaid_processing.log in CWD)', default=None)
    parser.add_argument('--prefer', choices=['input', 'base'], default='base',
                        help='Preferred record to keep when a cross-database duplicate is detected (default: base/GenBank)')
    parser.add_argument('--drop_unmapped', action='store_true',
                        help='If set, drop columns from the new dataset that are not in base matrix, mapping targets, or required columns')
    dedup_grp = parser.add_mutually_exclusive_group()
    dedup_grp.add_argument('--dedup_cross_source', dest='dedup_cross_source', action='store_true', default=True,
                           help='Enable cross-source duplicate detection against base matrix / update_db (default)')
    dedup_grp.add_argument('--no_dedup_cross_source', dest='dedup_cross_source', action='store_false',
                           help='Disable cross-source duplicate detection (legacy append-only behaviour)')

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
        base_fasta=args.base_fasta,
        update_db=args.update_db,
        dedup_cross_source=args.dedup_cross_source,
        dedup_fasta_out=args.dedup_fasta_out,
        dedup_report=args.dedup_report,
    )
    job.process()
