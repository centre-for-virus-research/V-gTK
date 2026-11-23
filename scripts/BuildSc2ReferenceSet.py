#!/usr/bin/env python3
"""Build SARS-CoV-2 Pango reference lists for RABV-gTK (GenBank, GISAID-hybrid, and diversity-tier).

Sources combined:
  1. cov-lineages/pango-designation:
       - lineage_notes.txt (active vs withdrawn '*' lineages)
       - pango_designation/alias_key.json (Pango alias root / major parent)
       - lineages.csv (taxon -> designated Pango lineage)
  2. Optional GISAID metadata (.tsv, .csv, .gz, .xz, or metadata_tsv_*.tar.xz):
       - Resolves exact GISAID `EPI_ISL_*` accessions for designated taxa in lineages.csv,
         enforces strict length / N-content QC, optionally verifies sequences against
         GISAID `sequences.fasta`, and optionally extracts a compact GISAID subset directory.
  3. UCSC public SARS-CoV-2 UShER metadata (`public-latest.metadata.tsv.gz`):
       - Resolves high-completeness GenBank accessions (preferring exact `lineages.csv`
         designated taxa first, then clean UShER/Pangolin-concordant GenBank genomes),
         with strict hard QC (`length >= min_length`, `ATGC >= min_real_length`,
         `N <= max_n`, and `(length - ATGC)/length <= max_ambig_fraction`).

Outputs headered RABV-gTK `--ref_list` TSVs:
    primary_accession  accession_type  segment  genotype  subtype
with `NC_045512` (Wuhan-Hu-1 RefSeq) as the single `master` row (`genotype=B`, `subtype=B`).
"""

import argparse
from contextlib import contextmanager
import csv
import gzip
import io
import json
import lzma
import os
import re
import sqlite3
import sys
import tarfile
import time
import urllib.request

PANGO_NOTES_URL = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
PANGO_ALIAS_URL = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
PANGO_LINEAGES_URL = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv"
UCSC_METADATA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.metadata.tsv.gz"

MASTER_ACCESSION = "NC_045512"

# Compact panel of landmark Pango lineages spanning root -> VOCs -> Omicron -> recombinants
TEST_PANEL_LINEAGES = [
    "A", "A.2", "B.1", "B.1.1.7", "B.1.351", "P.1", "B.1.525", "B.1.617.2",
    "AY.4", "AY.43", "BA.1", "BA.1.1", "BA.2", "BA.2.12.1", "BA.2.75", "BA.4",
    "BA.5", "BA.5.2", "BQ.1", "BQ.1.1", "CH.1.1", "XBB", "XBB.1.5", "XBB.1.16",
    "EG.5.1", "HK.3", "BA.2.86", "JN.1", "JN.1.4", "KP.2", "KP.3", "KP.3.1.1",
    "LB.1", "XEC", "LP.8.1", "NB.1.8.1", "XFG",
]


def normalize_taxon(name):
    if not name:
        return ""
    text = str(name).strip()
    if text.lower().startswith("hcov-19/"):
        text = text[8:].strip()
    return text


def derive_pango_major_parent(lineage, alias_key):
    """Derive Pango alias root / major parent (genotype) for a full Pango lineage (subtype).

    Examples:
      - 'A'        -> 'A'
      - 'B'        -> 'B'
      - 'B.1.1.7'  -> 'B.1'
      - 'BA.5.2.1' -> 'BA.5'
      - 'JN.1.4'   -> 'JN.1'
      - 'KP.3.1.1' -> 'KP.3'
      - 'XBB.1.5'  -> 'XBB'
      - 'XFG.3.1'  -> 'XFG'
    """
    lin = (lineage or "").strip()
    if not lin:
        return "NA"
    parts = lin.split(".")
    prefix = parts[0]
    target = alias_key.get(prefix)
    # Recombinant roots (alias value is a list, or X-prefixed root) stay at the recombinant root
    if isinstance(target, list) or prefix.startswith("X"):
        return prefix
    if len(parts) == 1:
        return prefix
    return f"{prefix}.{parts[1]}"


@contextmanager
def _open_text_or_url(source, timeout=60):
    """Open local plain/.gz/.xz/.tar.* or remote URL as a UTF-8 text stream."""
    lower = source.lower()
    if os.path.isfile(source):
        if lower.endswith((".tar.xz", ".txz", ".tar.gz", ".tgz", ".tar")):
            tar = tarfile.open(source, mode="r:*")
            try:
                chosen_stream = None
                for member in tar:
                    mlow = os.path.basename(member.name).lower()
                    if member.isfile() and mlow.endswith((".tsv", ".csv")) and mlow != "readme.txt":
                        raw = tar.extractfile(member)
                        if raw is not None:
                            chosen_stream = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                            yield chosen_stream
                            break
                if chosen_stream is None:
                    raise ValueError(f"No .tsv/.csv/.txt member found inside archive: {source}")
            finally:
                tar.close()
            return
        if lower.endswith(".gz"):
            with gzip.open(source, "rb") as raw:
                yield io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            return
        if lower.endswith(".xz"):
            with lzma.open(source, "rb") as raw:
                yield io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            return
        with open(source, "r", encoding="utf-8", errors="replace") as fh:
            yield fh
        return

    req = urllib.request.urlopen(source, timeout=timeout)
    try:
        if lower.endswith(".gz"):
            with gzip.GzipFile(fileobj=req) as gz:
                yield io.TextIOWrapper(gz, encoding="utf-8", errors="replace")
        else:
            yield io.TextIOWrapper(req, encoding="utf-8", errors="replace")
    finally:
        req.close()


def load_alias_key(source=PANGO_ALIAS_URL):
    with _open_text_or_url(source) as fh:
        return json.load(fh)


def load_active_lineages(source=PANGO_NOTES_URL):
    active = set()
    withdrawn = set()
    with _open_text_or_url(source) as fh:
        first = True
        for line in fh:
            if first:
                first = False
                continue
            line = line.strip()
            if not line:
                continue
            lin = line.split("\t", 1)[0].strip()
            if lin.startswith("*"):
                withdrawn.add(lin.lstrip("*"))
            elif lin:
                active.add(lin)
    return active, withdrawn


def load_pango_designations(source=PANGO_LINEAGES_URL):
    """Return ({taxon: lineage}, set_of_all_lineages)."""
    taxon_to_lin = {}
    all_lineages = set()
    with _open_text_or_url(source) as fh:
        reader = csv.reader(fh)
        for row in reader:
            if not row or row[0].strip().lower() == "taxon" or len(row) < 2:
                continue
            taxon = normalize_taxon(row[0])
            lin = row[1].strip()
            if taxon and lin:
                taxon_to_lin[taxon] = lin
                all_lineages.add(lin)
    return taxon_to_lin, all_lineages


def _infer_delimiter_from_path_or_header(path, header_line):
    lower = path.lower()
    for suffix in (".gz", ".xz"):
        if lower.endswith(suffix):
            lower = lower[:-len(suffix)]
    if lower.endswith(".csv"):
        return ","
    if lower.endswith(".tsv"):
        return "\t"
    return "\t" if "\t" in header_line else ","


def resolve_from_gisaid_metadata(
    metadata_path,
    taxon_to_lin,
    min_length=29500,
    max_n=150,
    max_ambig_fraction=0.005,
    top_k=4,
):
    """Scan GISAID metadata (.tsv, .csv, .gz, .xz, or .tar.xz) and return top EPI_ISL_* candidates per lineage."""
    if not metadata_path or not os.path.isfile(metadata_path):
        return {}
    candidates_by_lin = {}

    def _push(lin, cand):
        lst = candidates_by_lin.get(lin)
        if lst is None:
            candidates_by_lin[lin] = [cand]
            return
        if any(x["accession"] == cand["accession"] for x in lst):
            return
        if len(lst) < top_k:
            lst.append(cand)
            lst.sort(key=lambda x: x["score"], reverse=True)
        elif cand["score"] > lst[-1]["score"]:
            lst[-1] = cand
            lst.sort(key=lambda x: x["score"], reverse=True)

    with _open_text_or_url(metadata_path) as fh:
        first_line = fh.readline()
        if not first_line:
            return {}
        sep = _infer_delimiter_from_path_or_header(metadata_path, first_line)
        header = next(csv.reader([first_line], delimiter=sep))
        cols = [c.strip() for c in header]
        col_map = {c.lower(): idx for idx, c in enumerate(cols)}
        idx_acc = next((col_map[c] for c in ("accession_id", "accession id", "covv_accession_id", "isolate_id") if c in col_map), None)
        idx_name = next((col_map[c] for c in ("virus_name", "virus name", "covv_virus_name", "isolate_name") if c in col_map), None)
        idx_date = next((col_map[c] for c in ("collection_date", "collection date", "covv_collection_date") if c in col_map), None)
        idx_subm = next((col_map[c] for c in ("submission_date", "submission date", "covv_subm_date") if c in col_map), None)
        idx_len = next((col_map[c] for c in ("sequence length", "sequence_length", "covv_length") if c in col_map), None)
        idx_n = next((col_map[c] for c in ("n-content", "n_content") if c in col_map), None)
        idx_complete = next((col_map[c] for c in ("is complete?", "is_complete") if c in col_map), None)
        if idx_acc is None or idx_name is None:
            raise ValueError(f"Could not locate accession/virus_name columns in {metadata_path}: {cols}")

        row_iter = (line.rstrip("\r\n").split("\t") for line in fh) if sep == "\t" else csv.reader(fh, delimiter=sep)
        req_idx = max(idx_acc, idx_name)
        for row in row_iter:
            if len(row) <= req_idx:
                continue
            taxon = normalize_taxon(row[idx_name])
            lin = taxon_to_lin.get(taxon)
            if not lin:
                continue
            acc = row[idx_acc].strip().upper()
            if not acc.startswith("EPI_ISL_"):
                continue
            cdate = row[idx_date].strip() if idx_date is not None and len(row) > idx_date else ""
            sdate = row[idx_subm].strip() if idx_subm is not None and len(row) > idx_subm else ""
            seq_len = None
            if idx_len is not None and len(row) > idx_len and row[idx_len].strip():
                try:
                    seq_len = int(row[idx_len])
                except ValueError:
                    seq_len = 0
                if seq_len < min_length or seq_len > 30500:
                    continue

            n_frac = None
            if idx_n is not None and len(row) > idx_n and row[idx_n].strip():
                try:
                    n_frac = float(row[idx_n])
                except ValueError:
                    n_frac = 1.0
                if max_ambig_fraction is not None and n_frac > max_ambig_fraction:
                    continue
                if max_n is not None and seq_len is not None and (n_frac * seq_len) > (max_n + 0.5):
                    continue

            is_comp = 1
            if idx_complete is not None:
                comp_val = row[idx_complete].strip().lower() if len(row) > idx_complete else ""
                is_comp = 1 if comp_val in ("true", "1", "yes") else 0

            eff_len = seq_len if seq_len is not None else 0
            eff_n_frac = n_frac if n_frac is not None else 0.0
            has_full_date = 1 if re.fullmatch(r"\d{4}-\d{2}-\d{2}", cdate) else 0
            score = (is_comp, 1 if eff_len >= min_length else 0, -eff_n_frac, has_full_date, eff_len)
            _push(lin, {
                "accession": acc,
                "taxon": taxon,
                "date": cdate,
                "subm_date": sdate,
                "length": eff_len,
                "source": "gisaid_designated",
                "score": score,
            })

    return candidates_by_lin


def verify_gisaid_candidates_from_fasta(
    gisaid_fasta_path,
    gisaid_candidates_by_lin,
    min_length=29500,
    min_real_length=29350,
    max_n=150,
    max_ambig_fraction=0.005,
):
    """Stream GISAID sequences.fasta once in 16 MiB blocks and return {lineage: best_verified_candidate}."""
    triple_map = {}
    taxon_map = {}
    acc_map = {}
    for cands in gisaid_candidates_by_lin.values():
        for c in cands:
            tkey = c["taxon"].casefold()
            taxon_map.setdefault(tkey, []).append(c)
            acc_map[c["accession"].upper()] = c
            cdate = c.get("date", "")
            sdate = c.get("subm_date", "")
            if cdate or sdate:
                triple_map.setdefault((tkey, cdate, sdate), []).append(c)

    if not taxon_map and not acc_map:
        return {}

    qc_by_acc = {}
    target_accs = set(acc_map.keys())

    def _process_record(raw_block):
        hdr, sep_nl, body = raw_block.partition("\n")
        if not sep_nl:
            return
        fields = [f.strip() for f in hdr.rstrip("\r").split("|")]
        matched = []
        for f in fields:
            if "EPI_ISL_" in f or "epi_isl_" in f:
                m = re.search(r"\b(EPI_ISL_\d+)\b", f, flags=re.IGNORECASE)
                if m:
                    cand = acc_map.get(m.group(1).upper())
                    if cand and cand["accession"] not in qc_by_acc:
                        matched.append(cand)
        if not matched and fields:
            tkey = normalize_taxon(fields[0]).casefold()
            if len(fields) >= 3:
                for cand in triple_map.get((tkey, fields[1], fields[2]), ()):
                    if cand["accession"] not in qc_by_acc:
                        matched.append(cand)
            if not matched:
                t_cands = taxon_map.get(tkey)
                if t_cands and len(t_cands) == 1 and t_cands[0]["accession"] not in qc_by_acc:
                    matched.append(t_cands[0])
                elif t_cands and len(fields) < 3:
                    for cand in t_cands:
                        if cand["accession"] not in qc_by_acc:
                            matched.append(cand)
        if not matched:
            return
        seq_upper = body.replace("\n", "").replace("\r", "").replace(" ", "").upper()
        total_len = len(seq_upper)
        atgc = seq_upper.count("A") + seq_upper.count("T") + seq_upper.count("G") + seq_upper.count("C")
        n_cnt = seq_upper.count("N")
        for cand in matched:
            qc_by_acc[cand["accession"]] = (total_len, atgc, n_cnt)

    chunk_bytes = 16 * 1024 * 1024
    with _open_text_or_url(gisaid_fasta_path) as fh:
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
                _process_record(rec_block)
                if len(qc_by_acc) >= len(target_accs):
                    carry = ""
                    break
            if len(qc_by_acc) >= len(target_accs):
                break
        if carry and not first_chunk and len(qc_by_acc) < len(target_accs):
            _process_record(carry)

    verified_by_lin = {}
    for lin, cands in gisaid_candidates_by_lin.items():
        scored = []
        for idx, cand in enumerate(cands):
            qc = qc_by_acc.get(cand["accession"])
            if qc is None:
                continue
            total_len, atgc_len, n_cnt = qc
            if total_len < min_length or total_len > 30500:
                continue
            if atgc_len < min_real_length:
                continue
            if max_n is not None and n_cnt > max_n:
                continue
            if max_ambig_fraction is not None and total_len > 0 and ((total_len - atgc_len) / total_len) > max_ambig_fraction:
                continue
            enriched = dict(cand)
            enriched["length"] = total_len
            enriched["real_length"] = atgc_len
            enriched["n_count"] = n_cnt
            scored.append(((-n_cnt, atgc_len, -idx), enriched))
        if scored:
            scored.sort(key=lambda x: x[0], reverse=True)
            verified_by_lin[lin] = scored[0][1]
    return verified_by_lin


def _is_ena_padded_prefix(acc):
    """ENA COG-UK style accessions (OA/OD/OU/OV/OW/OX/OY/OZ) often pad consensuses to 29,903 with Ns."""
    return bool(re.match(r"^O[ADUVWXYZ]\d+$", acc, flags=re.IGNORECASE))


def resolve_from_ucsc_genbank(ucsc_source, taxon_to_lin, target_lineages, min_length=29500, top_k=8):
    """Stream UCSC public metadata and retain top_k GenBank candidates per lineage."""
    best_exact = {}
    best_called = {}

    def _push_candidate(store, lin, cand):
        lst = store.get(lin)
        if lst is None:
            store[lin] = [cand]
            return
        if any(x["accession"] == cand["accession"] for x in lst):
            return
        if len(lst) < top_k:
            lst.append(cand)
            lst.sort(key=lambda x: x["score"], reverse=True)
        elif cand["score"] > lst[-1]["score"]:
            lst[-1] = cand
            lst.sort(key=lambda x: x["score"], reverse=True)

    with _open_text_or_url(ucsc_source, timeout=120) as fh:
        header_line = fh.readline()
        if not header_line:
            return {}, {}
        header = header_line.rstrip("\r\n").split("\t")
        idx_strain = header.index("strain")
        idx_gb = header.index("genbank_accession")
        idx_date = header.index("date")
        idx_country = header.index("country")
        idx_len = header.index("length")
        idx_pango = header.index("pangolin_lineage")
        idx_usher = header.index("pango_lineage_usher")

        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) <= idx_usher:
                continue
            gb_raw = parts[idx_gb].strip()
            if not gb_raw or gb_raw == "?":
                continue
            bare_gb = gb_raw.split(".", 1)[0].strip()
            if not bare_gb or bare_gb == MASTER_ACCESSION:
                continue
            try:
                length = int(parts[idx_len])
            except ValueError:
                length = 0
            if length < min_length or length > 30500:
                continue

            strain_field = parts[idx_strain].split("|", 1)[0].strip()
            taxon = normalize_taxon(strain_field)
            date_str = parts[idx_date].strip()
            country = parts[idx_country].strip()
            has_full_date = 1 if re.fullmatch(r"\d{4}-\d{2}-\d{2}", date_str) else 0
            not_ena_padded = 0 if (_is_ena_padded_prefix(bare_gb) and length == 29903) else 1
            is_refseq_like = 1 if 29750 <= length <= 29910 else 0
            len_closeness = -abs(29860 - length)
            cand_score = (not_ena_padded, is_refseq_like, has_full_date, len_closeness, length)

            desig_lin = taxon_to_lin.get(taxon)
            if desig_lin and desig_lin in target_lineages:
                pango_agree = 1 if (parts[idx_pango].strip() == desig_lin or parts[idx_usher].strip() == desig_lin) else 0
                exact_score = (pango_agree,) + cand_score
                _push_candidate(best_exact, desig_lin, {
                    "accession": bare_gb,
                    "taxon": taxon,
                    "date": date_str,
                    "country": country,
                    "length": length,
                    "source": "genbank_designated",
                    "score": exact_score,
                })

            pl = parts[idx_pango].strip()
            ul = parts[idx_usher].strip()
            called_lin = pl if (pl and pl in target_lineages) else (ul if (ul and ul in target_lineages) else None)
            if called_lin:
                concordant = 1 if (pl == ul == called_lin) else 0
                called_score = (concordant,) + cand_score
                _push_candidate(best_called, called_lin, {
                    "accession": bare_gb,
                    "taxon": taxon,
                    "date": date_str,
                    "country": country,
                    "length": length,
                    "source": "genbank_concordant" if concordant else "genbank_usher_fallback",
                    "score": called_score,
                })

    return best_exact, best_called


def load_qc_from_sqlite(db_path):
    """Load cached {bare_acc: (total_len, atgc_len, n_count)} from an existing RABV-gTK SQLite DB."""
    if not db_path or not os.path.isfile(db_path):
        return {}
    stats = {}
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()
        rows = cur.execute(
            "SELECT primary_accession, length, real_length, n FROM meta_data "
            "WHERE primary_accession IS NOT NULL AND TRIM(primary_accession) != ''"
        ).fetchall()
        for acc, length, real_len, n_cnt in rows:
            try:
                stats[str(acc).split(".", 1)[0].strip()] = (int(length), int(real_len), int(n_cnt))
            except (TypeError, ValueError):
                continue
    finally:
        conn.close()
    return stats


def _load_ncbi_api_key():
    env_key = os.environ.get("NCBI_API_KEY", "").strip()
    if len(env_key) == 36:
        return env_key
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for cand in (os.path.join(repo_root, "dev", "ncbi_key"), os.path.join(repo_root, "ncbi_key")):
        if os.path.isfile(cand):
            key = open(cand, encoding="utf-8").read().strip()
            if len(key) == 36:
                return key
    return None


def fetch_ncbi_fasta_qc(accessions, batch_size=200, existing_qc=None):
    """Fetch FASTA from NCBI for candidate accessions and return {bare_acc: (total_len, atgc_len, n_count)}."""
    stats = dict(existing_qc or {})
    unique_accs = [a for a in dict.fromkeys(a for a in accessions if a) if a not in stats]
    if not unique_accs:
        return stats

    api_key = _load_ncbi_api_key()
    interval = 0.12 if api_key else 0.36
    last_req = 0.0
    for i in range(0, len(unique_accs), batch_size):
        chunk = unique_accs[i:i + batch_size]
        text = None
        for attempt in range(1, 5):
            elapsed = time.time() - last_req
            if elapsed < interval:
                time.sleep(interval - elapsed)
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={','.join(chunk)}&rettype=fasta&retmode=text"
            if api_key:
                url += f"&api_key={api_key}"
            last_req = time.time()
            try:
                text = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", errors="replace")
                break
            except Exception as exc:
                if attempt == 4:
                    print(f"[warn] efetch QC batch {i} failed after 4 attempts: {exc}")
                else:
                    time.sleep(1.5 * attempt)
        if not text:
            continue

        cur_acc = None
        seq_parts = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_acc and seq_parts:
                    s = "".join(seq_parts).upper()
                    atgc = sum(s.count(b) for b in "ATGC")
                    stats[cur_acc] = (len(s), atgc, s.count("N"))
                token = line[1:].split()[0].split("|")[-1]
                cur_acc = token.split(".", 1)[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if cur_acc and seq_parts:
            s = "".join(seq_parts).upper()
            atgc = sum(s.count(b) for b in "ATGC")
            stats[cur_acc] = (len(s), atgc, s.count("N"))
    return stats


def extract_gisaid_subset(gisaid_fasta_path, selected_records, output_dir):
    """Stream a GISAID sequences.fasta once and write metadata.tsv + sequences.fasta for selected EPI_ISL_* rows."""
    os.makedirs(output_dir, exist_ok=True)
    by_taxon = {rec["taxon"].casefold(): rec for rec in selected_records.values() if rec["source"].startswith("gisaid")}
    if not by_taxon:
        print("[extract] No GISAID-sourced references selected; nothing to extract.")
        return

    meta_out = os.path.join(output_dir, "metadata.tsv")
    fasta_out = os.path.join(output_dir, "sequences.fasta")
    found_taxa = set()

    with open(gisaid_fasta_path, "r", encoding="utf-8", errors="replace") as fh, \
         open(fasta_out, "w", encoding="utf-8") as out_fa:
        keep_current = None
        chunks = []
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if keep_current is not None and chunks:
                    seq = "".join(chunks).upper()
                    out_fa.write(f">{keep_current['accession']}\n{seq}\n")
                    found_taxa.add(keep_current["taxon"].casefold())
                    if len(found_taxa) == len(by_taxon):
                        keep_current = None
                        break
                fields = [f.strip() for f in line[1:].split("|")]
                taxon = normalize_taxon(fields[0]).casefold() if fields else ""
                rec = by_taxon.get(taxon)
                if rec and taxon not in found_taxa:
                    keep_current = rec
                    chunks = []
                else:
                    keep_current = None
                    chunks = []
            elif keep_current is not None:
                chunks.append(line)
        if keep_current is not None and chunks and keep_current["taxon"].casefold() not in found_taxa:
            seq = "".join(chunks).upper()
            out_fa.write(f">{keep_current['accession']}\n{seq}\n")
            found_taxa.add(keep_current["taxon"].casefold())

    with open(meta_out, "w", encoding="utf-8", newline="") as out_tsv:
        writer = csv.writer(out_tsv, delimiter="\t", lineterminator="\n")
        writer.writerow(["Accession_Id", "Virus name", "Collection date", "Country", "Host"])
        for taxon_key in sorted(found_taxa):
            rec = by_taxon[taxon_key]
            tparts = [p.strip() for p in rec["taxon"].split("/") if p.strip()]
            country = tparts[0] if len(tparts) >= 2 else ""
            writer.writerow([rec["accession"], f"hCoV-19/{rec['taxon']}", rec.get("date", ""), country, "Human"])

    print(f"[extract] Wrote {len(found_taxa)}/{len(by_taxon)} GISAID reference sequences to {output_dir}")


def write_ref_list(path, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["primary_accession", "accession_type", "segment", "genotype", "subtype"])
        seen_acc = {MASTER_ACCESSION}
        writer.writerow([MASTER_ACCESSION, "master", "1", "B", "B"])
        for acc, gtype, stype in rows:
            if acc in seen_acc:
                continue
            seen_acc.add(acc)
            writer.writerow([acc, "reference", "1", gtype, stype])


def build_diversity_subset(lineages_sorted, chosen_by_lin, alias_key):
    """Select 1 highest-scoring GenBank representative per unique Pango major parent (`genotype`)."""
    best_per_major = {}
    for lin in lineages_sorted:
        rec = chosen_by_lin.get(lin)
        if not rec or rec["accession"].startswith("EPI_ISL_"):
            continue
        major = derive_pango_major_parent(lin, alias_key)
        is_exact_major = 1 if lin == major else 0
        n_cnt = rec.get("n_count", 0)
        rank = (is_exact_major, 1 if rec["source"] == "genbank_designated" else 0, -n_cnt, rec["length"])
        prev = best_per_major.get(major)
        if prev is None or rank > prev[0]:
            best_per_major[major] = (rank, lin, rec)
    out = []
    for major in sorted(best_per_major):
        _, lin, rec = best_per_major[major]
        out.append((rec["accession"], major, lin))
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--lineages-csv", default=PANGO_LINEAGES_URL, help="Path or URL to pango-designation lineages.csv")
    parser.add_argument("--lineage-notes", default=PANGO_NOTES_URL, help="Path or URL to pango-designation lineage_notes.txt")
    parser.add_argument("--alias-key", default=PANGO_ALIAS_URL, help="Path or URL to pango-designation alias_key.json")
    parser.add_argument("--ucsc-metadata", default=UCSC_METADATA_URL, help="Path or URL to UCSC public-latest.metadata.tsv.gz")
    parser.add_argument("--gisaid-metadata", default=None,
                        help="Optional GISAID metadata (.tsv, .csv, .gz, .xz, or metadata_tsv_*.tar.xz)")
    parser.add_argument("--gisaid-fasta", default=None, help="Optional GISAID sequences.fasta for sequence QC verification / subset extraction")
    parser.add_argument("--verify-gisaid-fasta", action="store_true",
                        help="Verify GISAID reference candidates directly against --gisaid-fasta (length, ATGC, N<=max-n)")
    parser.add_argument("--extract-gisaid-dir", default=None, help="Output directory for extracted GISAID reference subset")
    parser.add_argument("--prefer", choices=["gisaid", "genbank"], default="genbank",
                        help="Primary source when a lineage is available in both GISAID and GenBank (default: genbank)")
    parser.add_argument("--include-withdrawn", action="store_true", help="Include withdrawn '*' lineages from lineage_notes.txt")
    parser.add_argument("--min-length", type=int, default=29500, help="Minimum total sequence length (default: 29500)")
    parser.add_argument("--min-real-length", type=int, default=29350, help="Minimum unambiguous ATGC count (default: 29350)")
    parser.add_argument("--max-n", type=int, default=150, help="Maximum allowed N bases in any reference sequence (default: 150)")
    parser.add_argument("--max-ambig-fraction", type=float, default=0.005,
                        help="Maximum allowed non-ATGC fraction in any reference sequence (default: 0.005)")
    parser.add_argument("--qc-cache-db", default=None,
                        help="Optional existing RABV-gTK SQLite DB to seed NCBI length/ATGC/N QC before efetch")
    parser.add_argument("--output-all", required=True, help="Output TSV path for full Pango reference list")
    parser.add_argument("--output-diversity", default=None, help="Optional output TSV path for 1-per-major-clade diversity reference list")
    parser.add_argument("--output-gisaid-hybrid", default=None, help="Optional output TSV path preferring exact GISAID EPI_ISL_* with GenBank fallback")
    parser.add_argument("--output-test", default=None, help="Optional output TSV path for compact test panel (~37 landmark lineages)")
    parser.add_argument("--manifest", default=None, help="Optional provenance manifest TSV")
    parser.add_argument("--verify-ncbi-diversity", action="store_true",
                        help="Verify FASTA from NCBI for diversity/test candidates and enforce strict QC")
    parser.add_argument("--verify-ncbi-all", action="store_true",
                        help="Verify FASTA from NCBI for ALL candidate lineages and enforce strict QC")
    args = parser.parse_args(argv)

    print("[1/4] Loading Pango alias_key.json, lineage_notes.txt, and lineages.csv...")
    alias_key = load_alias_key(args.alias_key)
    active_lins, withdrawn_lins = load_active_lineages(args.lineage_notes)
    taxon_to_lin, desig_lins = load_pango_designations(args.lineages_csv)
    target_lins = desig_lins if args.include_withdrawn else (desig_lins - withdrawn_lins)
    print(f"      Designated lineages: {len(desig_lins):,} total, {len(target_lins):,} active target lineages ({len(taxon_to_lin):,} designated taxa)")

    gisaid_by_lin = {}
    if args.gisaid_metadata:
        print(f"[2/4] Resolving exact GISAID EPI_ISL_* candidates from {args.gisaid_metadata}...")
        gisaid_cands = resolve_from_gisaid_metadata(
            args.gisaid_metadata,
            taxon_to_lin,
            min_length=args.min_length,
            max_n=args.max_n,
            max_ambig_fraction=args.max_ambig_fraction,
        )
        print(f"      Matched {len(gisaid_cands):,} lineages in GISAID metadata")
        if args.verify_gisaid_fasta and args.gisaid_fasta:
            print(f"      Verifying GISAID candidates against {args.gisaid_fasta} "
                  f"(len>={args.min_length}, ATGC>={args.min_real_length}, N<={args.max_n}, ambig<={args.max_ambig_fraction})...")
            gisaid_by_lin = verify_gisaid_candidates_from_fasta(
                args.gisaid_fasta,
                gisaid_cands,
                min_length=args.min_length,
                min_real_length=args.min_real_length,
                max_n=args.max_n,
                max_ambig_fraction=args.max_ambig_fraction,
            )
            print(f"      Verified {len(gisaid_by_lin):,} clean GISAID reference lineages")
        else:
            gisaid_by_lin = {lin: lst[0] for lin, lst in gisaid_cands.items() if lst}
    else:
        print("[2/4] No --gisaid-metadata provided; skipping GISAID pass")

    print(f"[3/4] Scanning UCSC public GenBank metadata ({args.ucsc_metadata}) with min_length={args.min_length}...")
    gb_exact_lists, gb_called_lists = resolve_from_ucsc_genbank(
        args.ucsc_metadata, taxon_to_lin, target_lins, min_length=args.min_length, top_k=8
    )
    print(f"      GenBank exact designated matches: {len(gb_exact_lists):,}; GenBank called matches: {len(gb_called_lists):,}")

    lineages_sorted = sorted(target_lins)
    ncbi_qc = load_qc_from_sqlite(args.qc_cache_db) if args.qc_cache_db else {}
    if args.manifest and os.path.isfile(args.manifest):
        with open(args.manifest, "r", encoding="utf-8", errors="replace") as mfh:
            for mrow in csv.DictReader(mfh, delimiter="\t"):
                gba = (mrow.get("genbank_accession") or "").strip()
                pacc = (mrow.get("primary_accession") or "").strip()
                if gba and gba == pacc and gba not in ncbi_qc:
                    try:
                        ncbi_qc[gba] = (int(mrow["length"]), int(mrow["real_length"]), int(mrow["n_count"]))
                    except (KeyError, TypeError, ValueError):
                        pass
    if ncbi_qc:
        print(f"      Loaded {len(ncbi_qc):,} cached QC entries")

    def _passes_strict_qc(qc):
        if qc is None:
            return False
        total_len, atgc_len, n_cnt = qc
        if total_len < args.min_length or total_len > 30500:
            return False
        if atgc_len < args.min_real_length:
            return False
        if args.max_n is not None and n_cnt > args.max_n:
            return False
        if args.max_ambig_fraction is not None and total_len > 0:
            if ((total_len - atgc_len) / total_len) > args.max_ambig_fraction:
                return False
        return True

    strict_ncbi = bool(args.verify_ncbi_all or args.verify_ncbi_diversity or args.qc_cache_db)
    if args.verify_ncbi_all or args.verify_ncbi_diversity:
        if args.verify_ncbi_all:
            priority_lins = set(target_lins)
            per_lin_cap = 6
        else:
            major_roots = {derive_pango_major_parent(lin, alias_key) for lin in lineages_sorted}
            priority_lins = (major_roots | set(TEST_PANEL_LINEAGES)) & target_lins
            per_lin_cap = 6
        cands_to_check = []
        for lin in sorted(priority_lins):
            pool = (gb_exact_lists.get(lin, []) + gb_called_lists.get(lin, []))[:per_lin_cap]
            # If we already have a cached candidate for this lineage that passes strict QC, skip extra network fetches
            if any(_passes_strict_qc(ncbi_qc.get(c["accession"])) for c in pool):
                continue
            for c in pool:
                if c["accession"] not in ncbi_qc:
                    cands_to_check.append(c["accession"])
        unique_missing = list(dict.fromkeys(cands_to_check))
        print(f"      Fetching NCBI FASTA QC for {len(unique_missing):,} uncached candidate accessions...")
        ncbi_qc = fetch_ncbi_fasta_qc(unique_missing, existing_qc=ncbi_qc)
        print(f"      Total verified NCBI candidate sequences in QC map: {len(ncbi_qc):,}")

    def _pick_best_gb(lin):
        pool = gb_exact_lists.get(lin, []) + gb_called_lists.get(lin, [])
        if not pool:
            return None
        if strict_ncbi:
            scored = []
            for idx, cand in enumerate(pool):
                qc = ncbi_qc.get(cand["accession"])
                if not _passes_strict_qc(qc):
                    continue
                total_len, atgc_len, n_cnt = qc
                enriched = dict(cand)
                enriched["length"] = total_len
                enriched["real_length"] = atgc_len
                enriched["n_count"] = n_cnt
                scored.append((
                    (
                        1 if cand["source"] == "genbank_designated" else 0,
                        1 if atgc_len >= 29500 else 0,
                        -n_cnt,
                        atgc_len,
                        -idx,
                    ),
                    enriched,
                ))
            if scored:
                scored.sort(key=lambda x: x[0], reverse=True)
                return scored[0][1]
            return None
        return pool[0]

    gb_chosen = {}
    hybrid_chosen = {}
    manifest_rows = []

    for lin in lineages_sorted:
        major = derive_pango_major_parent(lin, alias_key)
        gb_rec = _pick_best_gb(lin)
        if gb_rec:
            gb_chosen[lin] = gb_rec

        gis_rec = gisaid_by_lin.get(lin)
        if args.prefer == "gisaid":
            hyb_rec = gis_rec or gb_rec
        else:
            hyb_rec = gb_rec or gis_rec
        if hyb_rec:
            hybrid_chosen[lin] = hyb_rec

        primary_rec = hyb_rec if args.gisaid_metadata else gb_rec
        manifest_rows.append({
            "lineage": lin,
            "genotype_major": major,
            "primary_accession": primary_rec["accession"] if primary_rec else "",
            "source": primary_rec["source"] if primary_rec else "unresolved",
            "gisaid_epi_isl": gis_rec["accession"] if gis_rec else "",
            "genbank_accession": gb_rec["accession"] if gb_rec else "",
            "genbank_source": gb_rec["source"] if gb_rec else "",
            "taxon": primary_rec["taxon"] if primary_rec else "",
            "length": primary_rec.get("length", "") if primary_rec else "",
            "real_length": primary_rec.get("real_length", "") if primary_rec else "",
            "n_count": primary_rec.get("n_count", "") if primary_rec else "",
            "date": primary_rec.get("date", "") if primary_rec else "",
        })

    print("[4/4] Writing reference lists...")
    primary_map = hybrid_chosen if (args.gisaid_metadata and args.prefer == "gisaid") else gb_chosen
    all_rows = [(primary_map[lin]["accession"], derive_pango_major_parent(lin, alias_key), lin)
                for lin in lineages_sorted if lin in primary_map]
    write_ref_list(args.output_all, all_rows)
    print(f"      Wrote {args.output_all} ({len(all_rows) + 1:,} rows incl. master {MASTER_ACCESSION})")

    if args.output_gisaid_hybrid and gisaid_by_lin:
        gisaid_first_rows = []
        gisaid_selected_map = {}
        for lin in lineages_sorted:
            rec = gisaid_by_lin.get(lin) or gb_chosen.get(lin)
            if rec:
                gisaid_first_rows.append((rec["accession"], derive_pango_major_parent(lin, alias_key), lin))
                gisaid_selected_map[lin] = rec
        write_ref_list(args.output_gisaid_hybrid, gisaid_first_rows)
        n_gis = sum(1 for r in gisaid_selected_map.values() if r["source"].startswith("gisaid"))
        n_gbf = len(gisaid_selected_map) - n_gis
        print(f"      Wrote {args.output_gisaid_hybrid} ({len(gisaid_first_rows) + 1:,} rows: {n_gis:,} exact GISAID EPI_ISL_* + {n_gbf:,} GenBank fallback + master)")
        if args.gisaid_fasta and args.extract_gisaid_dir:
            extract_gisaid_subset(args.gisaid_fasta, gisaid_selected_map, args.extract_gisaid_dir)

    if args.output_diversity:
        div_rows = build_diversity_subset(lineages_sorted, gb_chosen, alias_key)
        write_ref_list(args.output_diversity, div_rows)
        print(f"      Wrote {args.output_diversity} ({len(div_rows) + 1:,} diversity-spanning GenBank rows incl. master)")

    if args.output_test:
        test_rows = []
        for lin in TEST_PANEL_LINEAGES:
            rec = gb_chosen.get(lin)
            if rec:
                test_rows.append((rec["accession"], derive_pango_major_parent(lin, alias_key), lin))
        write_ref_list(args.output_test, test_rows)
        print(f"      Wrote {args.output_test} ({len(test_rows) + 1:,} test panel rows incl. master)")

    if args.manifest:
        os.makedirs(os.path.dirname(os.path.abspath(args.manifest)), exist_ok=True)
        with open(args.manifest, "w", encoding="utf-8", newline="") as fh:
            fieldnames = [
                "lineage", "genotype_major", "primary_accession", "source",
                "gisaid_epi_isl", "genbank_accession", "genbank_source",
                "taxon", "length", "real_length", "n_count", "date",
            ]
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest_rows)
        print(f"      Wrote provenance manifest {args.manifest} ({len(manifest_rows):,} lineages)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
