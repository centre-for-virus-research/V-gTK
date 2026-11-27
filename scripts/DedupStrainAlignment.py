#!/usr/bin/env python3
"""Post-alignment strain-aware sequence deduplication gate for vgtk-init.nf.

Runs after ``PAD_ALIGNMENT`` (when ``Parsed_strain``, ``segment``, and padded
coordinate alignments exist for every batch sequence) and before ``DEDUP_ALIGNMENT``,
``MMSEQS_CLUSTERING``, ``USHER_PLACEMENT``, ``GENERATE_TABLES``, and ``CREATE_SQLITE_DB``.

Guarantees across both fresh builds and ``--update_db`` runs (in ANY arrival order):
1. Intra-batch duplicates (GenBank+GenBank, GISAID+GISAID, and cross-DB GenBank+GISAID)
   with compatible passage and 0 unambiguous ACGT mismatches (or <=2 mismatches when
   explicitly INSDC-linked) are collapsed to a single canonical winner per
   ``(segment, canonical_isolate)``.
2. Incoming duplicates of an already-committed ``--update_db`` winner are retired from
   the incoming batch (marked ``exclusion_status='1'`` in ``meta_data`` and omitted
   from the clean padded MSA so UShER never places a 0-length duplicate cherry).
3. Reverse arrival order in ``--update_db`` mode: when an incoming preferred-source
   (GenBank/INSDC) record matches an existing external (GISAID/Pathoplexus) record
   already in ``--update_db``, the incoming GenBank record wins, inherits coalesced
   metadata from the ``--update_db`` record, stays in the clean MSA for placement,
   and records the superseded ``--update_db`` accession in ``--out_retired_update_ids``
   so ``UsherPlacement.py`` prunes it from the starter tree and ``CreateSqliteDB.py``
   retires it cleanly in SQLite.
4. Protected ``master`` and ``reference`` accessions are never retired or omitted.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from accession_utils import normalise_accession  # noqa: E402
import segment_utils  # noqa: E402
from cross_db_dedup import (  # noqa: E402
    aligned_sequences_agree,
    candidate_winner_sort_key,
    coalesce_winner_metadata,
    count_acgt,
    extract_canonical_isolate,
    extract_insdc_links,
    extract_passage_class,
    is_external_quality_upgrade,
    is_preferred_source,
    passages_compatible,
    primary_sequence_source,
    sequence_quality_stats,
)


def _normalize_seg(val: object) -> str:
    norm = segment_utils.normalise_segment(val)
    if norm is None or norm.casefold() in segment_utils.PANDAS_NULL_TOKENS:
        return ""
    return norm


def _is_excluded_status(val: object) -> bool:
    if val is None:
        return False
    text = str(val).strip().lower()
    return text not in ("", "0", "false", "no", "na", "none", "nan")


def _load_protected_accessions(ref_list_path: Optional[str]) -> Set[str]:
    protected: Set[str] = set()
    if not ref_list_path or not os.path.isfile(ref_list_path):
        return protected
    with open(ref_list_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            first = raw.split("\t")[0].split()[0].strip()
            if first.lower() in ("primary_accession", "accession", "id"):
                continue
            protected.add(first)
            bare = normalise_accession(first)
            if bare:
                protected.add(bare)
    return protected


_SEGMENT_IN_FILENAME_RE = re.compile(
    r"(?:^|[_.-])(?:segment|seg|refset)[_-]?(\d+)(?:$|[_.-])",
    re.IGNORECASE,
)


def _infer_msa_segment(msa_path: Path, record_ids: Iterable[str], acc_to_segs: Dict[str, Set[str]]) -> str:
    found = {
        segment_utils.normalise_segment(tok)
        for tok in _SEGMENT_IN_FILENAME_RE.findall(msa_path.name)
    }
    found.discard(None)
    if len(found) == 1:
        return next(iter(found))

    seg_counts: Dict[str, int] = defaultdict(int)
    for rid in record_ids:
        bare = normalise_accession(rid) or rid
        for seg in acc_to_segs.get(bare, ()):
            if seg:
                seg_counts[seg] += 1
    if len(seg_counts) == 1:
        return next(iter(seg_counts.keys()))
    if seg_counts:
        return sorted(seg_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    return ""


def _load_update_db_active_candidates(
    update_db_path: str,
    prefer_source: str = "genbank",
) -> Tuple[Dict[Tuple[str, str], Dict[str, object]], Dict[Tuple[str, str], Dict[str, str]], Dict[Tuple[str, str], str]]:
    """Load active ``meta_data`` + ``sequence_alignment`` rows from ``update_db`` (read-only)."""
    cand_by_key: Dict[Tuple[str, str], Dict[str, object]] = {}
    raw_by_key: Dict[Tuple[str, str], Dict[str, str]] = {}
    aln_by_key: Dict[Tuple[str, str], str] = {}

    if not update_db_path or str(update_db_path).strip().lower() in ("", "null", "none", "unset"):
        return cand_by_key, raw_by_key, aln_by_key
    if not os.path.isfile(update_db_path):
        return cand_by_key, raw_by_key, aln_by_key

    conn = sqlite3.connect(f"file:{update_db_path}?mode=ro", uri=True)
    try:
        tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "meta_data" not in tables or "sequence_alignment" not in tables:
            return cand_by_key, raw_by_key, aln_by_key

        meta_cols = [r[1] for r in conn.execute("PRAGMA table_info(meta_data)")]
        aln_cols = [r[1] for r in conn.execute("PRAGMA table_info(sequence_alignment)")]
        aln_acc_col = next((c for c in ("primary_accession", "accession", "sequence_id") if c in aln_cols), None)
        aln_seq_col = next((c for c in ("alignment", "aligned_seq", "sequence") if c in aln_cols), None)
        if not aln_acc_col or not aln_seq_col:
            return cand_by_key, raw_by_key, aln_by_key

        has_aln_seg = "segment" in aln_cols
        aln_sql = (
            f'SELECT "{aln_acc_col}", "{aln_seq_col}"'
            + (', "segment"' if has_aln_seg else "")
            + " FROM sequence_alignment"
        )
        raw_aln_by_acc: Dict[str, str] = {}
        for row in conn.execute(aln_sql):
            acc = str(row[0] or "").strip()
            seq = str(row[1] or "").strip()
            if not acc or not seq:
                continue
            bare = normalise_accession(acc) or acc
            seg = _normalize_seg(row[2]) if has_aln_seg else ""
            aln_by_key[(bare, seg)] = seq
            raw_aln_by_acc.setdefault(bare, seq)

        qcols = ", ".join(f'"{c}"' for c in meta_cols)
        for raw_tuple in conn.execute(f"SELECT {qcols} FROM meta_data"):
            row = {col: ("" if val is None else str(val)) for col, val in zip(meta_cols, raw_tuple)}
            if _is_excluded_status(row.get("exclusion_status", "")):
                continue
            acc_raw = row.get("primary_accession", "").strip()
            if not acc_raw:
                continue
            acc = normalise_accession(acc_raw) or acc_raw
            seg = _normalize_seg(row.get("segment", ""))
            seq = aln_by_key.get((acc, seg)) or raw_aln_by_acc.get(acc, "")
            if not seq:
                continue
            aln_by_key[(acc, seg)] = seq
            aln_acgt, aln_amb = sequence_quality_stats(seq)
            acc_type = row.get("accession_type", "").strip().lower()
            protected = acc_type in ("master", "reference")
            src = (row.get("dataset_source") or "genbank").strip().lower() or "genbank"
            seq_src = (row.get("sequence_source") or "").strip().lower() or primary_sequence_source(src)
            key = (acc, seg)
            raw_by_key[key] = row
            cand_by_key[key] = {
                "accession": acc,
                "bare_upper": acc.upper(),
                "segment": seg,
                "isolate": extract_canonical_isolate(row),
                "passage": extract_passage_class(row),
                "real_length": aln_acgt,
                "aligned_acgt": aln_acgt,
                "aligned_ambig": aln_amb,
                "protected": protected,
                "source": src,
                "sequence_source": seq_src,
                "collection_date": row.get("collection_date", ""),
                "origin": "update_db",
                "preferred": is_preferred_source(src, prefer_source) or protected,
            }
    finally:
        conn.close()

    return cand_by_key, raw_by_key, aln_by_key


def deduplicate_strain_alignments(
    df_meta: pd.DataFrame,
    msa_files: Sequence[Path],
    out_msa_dir: Path,
    ref_list_path: Optional[str] = None,
    update_db_path: Optional[str] = None,
    prefer_source: str = "genbank",
) -> Tuple[pd.DataFrame, Set[str], List[Dict[str, str]]]:
    """Deduplicate batch metadata + padded MSA files against themselves and ``update_db``.

    Returns ``(updated_df_meta, retired_update_db_accessions, audit_pairs)``.
    """
    df = df_meta.copy()
    if "primary_accession" not in df.columns:
        os.makedirs(out_msa_dir, exist_ok=True)
        for msa in msa_files:
            target = out_msa_dir / msa.name
            target.write_text(msa.read_text(encoding="utf-8", errors="replace"), encoding="utf-8")
        return df, set(), []

    if "segment" not in df.columns:
        df["segment"] = ""
    if "exclusion_status" not in df.columns:
        df["exclusion_status"] = ""
    if "exclusion_criteria" not in df.columns:
        df["exclusion_criteria"] = ""

    protected_ids = _load_protected_accessions(ref_list_path)

    # Normalise segment & accession lookup for batch rows
    acc_to_segs: Dict[str, Set[str]] = defaultdict(set)
    batch_rows_by_key: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for idx, row in df.iterrows():
        raw_acc = str(row.get("primary_accession") or "").strip()
        if not raw_acc:
            continue
        acc = normalise_accession(raw_acc) or raw_acc
        seg = _normalize_seg(row.get("segment", ""))
        acc_to_segs[acc].add(seg)
        batch_rows_by_key[(acc, seg)].append(idx)
        acc_type = str(row.get("accession_type") or "").strip().lower()
        if acc_type in ("master", "reference"):
            protected_ids.add(acc)
            protected_ids.add(raw_acc)

    # Read MSA sequences
    msa_records_per_file: List[Tuple[Path, str, List[object]]] = []
    batch_aln_by_key: Dict[Tuple[str, str], str] = {}
    batch_aln_fallback: Dict[str, str] = {}
    for msa_path in msa_files:
        with open(msa_path, "r", encoding="utf-8", errors="replace") as fh:
            records = list(SeqIO.parse(fh, "fasta"))
        rec_ids = [r.id.strip() for r in records if r.id.strip()]
        file_seg = _infer_msa_segment(msa_path, rec_ids, acc_to_segs)
        msa_records_per_file.append((msa_path, file_seg, records))
        for rec in records:
            rid = rec.id.strip()
            if not rid:
                continue
            acc = normalise_accession(rid) or rid
            seq_str = str(rec.seq)
            segs_for_acc = acc_to_segs.get(acc)
            if file_seg and (not segs_for_acc or file_seg in segs_for_acc):
                batch_aln_by_key[(acc, file_seg)] = seq_str
            elif segs_for_acc and len(segs_for_acc) == 1:
                batch_aln_by_key[(acc, next(iter(segs_for_acc)))] = seq_str
            else:
                batch_aln_by_key[(acc, file_seg)] = seq_str
            batch_aln_fallback.setdefault(acc, seq_str)

    # Load update_db candidates (if any)
    upd_cands, upd_raw, upd_aln = _load_update_db_active_candidates(update_db_path, prefer_source=prefer_source)

    # Build combined candidate pool per segment
    cands_by_seg: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    raw_by_key: Dict[Tuple[str, str], Dict[str, str]] = dict(upd_raw)
    aln_by_key: Dict[Tuple[str, str], str] = dict(upd_aln)

    # Determine whether segment normalization collapses "" and "1" for unsegmented builds
    all_segs = {seg for (_, seg) in batch_rows_by_key.keys()} | {seg for (_, seg) in upd_cands.keys()}
    nonblank_segs = {s for s in all_segs if s}
    unify_single_seg = len(nonblank_segs) <= 1 and ("" in all_segs)
    canonical_single_seg = next(iter(nonblank_segs), "")

    def _canon_seg(s: str) -> str:
        return canonical_single_seg if (unify_single_seg and not s) else s

    seen_batch_keys: Set[Tuple[str, str]] = set()
    for (acc, orig_seg), indices in batch_rows_by_key.items():
        idx = indices[0]
        row_series = df.loc[idx]
        row_dict = {col: ("" if pd.isna(val) else str(val)) for col, val in row_series.items()}
        if _is_excluded_status(row_dict.get("exclusion_status", "")):
            continue
        if str(row_dict.get("exclusion_criteria") or "").strip():
            continue
        seq = batch_aln_by_key.get((acc, orig_seg)) or batch_aln_fallback.get(acc, "")
        if not seq:
            continue
        cseg = _canon_seg(orig_seg)
        key = (acc, cseg)
        seen_batch_keys.add(key)
        raw_by_key[key] = row_dict
        aln_by_key[key] = seq
        protected = acc in protected_ids or str(row_dict.get("accession_type") or "").strip().lower() in ("master", "reference")
        src = (row_dict.get("dataset_source") or "genbank").strip().lower() or "genbank"
        seq_src = (row_dict.get("sequence_source") or "").strip().lower() or primary_sequence_source(src)
        aln_acgt, aln_amb = sequence_quality_stats(seq)
        cands_by_seg[cseg].append({
            "accession": acc,
            "bare_upper": acc.upper(),
            "segment": cseg,
            "orig_segment": orig_seg,
            "isolate": extract_canonical_isolate(row_dict),
            "passage": extract_passage_class(row_dict),
            "real_length": aln_acgt,
            "aligned_acgt": aln_acgt,
            "aligned_ambig": aln_amb,
            "protected": protected,
            "source": src,
            "sequence_source": seq_src,
            "collection_date": row_dict.get("collection_date", ""),
            "origin": "batch",
            "preferred": is_preferred_source(src, prefer_source) or protected,
        })

    for (acc, orig_seg), ucand in upd_cands.items():
        cseg = _canon_seg(orig_seg)
        key = (acc, cseg)
        if key in seen_batch_keys:
            # Re-supplied same accession in update batch: incoming row upserts it directly
            continue
        ucopy = dict(ucand)
        ucopy["segment"] = cseg
        ucopy["orig_segment"] = orig_seg
        raw_by_key[key] = upd_raw[(acc, orig_seg)]
        aln_by_key[key] = upd_aln[(acc, orig_seg)]
        cands_by_seg[cseg].append(ucopy)

    retired_batch_keys: Dict[Tuple[str, str], Tuple[str, str]] = {}
    retired_update_ids: Set[str] = set()
    winner_meta_updates: Dict[Tuple[str, str], Dict[str, str]] = defaultdict(dict)
    audit_pairs: List[Dict[str, str]] = []

    for cseg, seg_cands in sorted(cands_by_seg.items()):
        pref_by_iso: Dict[str, List[Dict[str, object]]] = defaultdict(list)
        pref_by_bare: Dict[str, Dict[str, object]] = {}
        externals: List[Dict[str, object]] = []
        ext_by_iso: Dict[str, List[Dict[str, object]]] = defaultdict(list)

        for c in seg_cands:
            if c["preferred"]:
                pref_by_bare[str(c["bare_upper"])] = c
                if c["isolate"]:
                    pref_by_iso[str(c["isolate"])].append(c)
            else:
                externals.append(c)
                if c["isolate"]:
                    ext_by_iso[str(c["isolate"])].append(c)

        # Helper for intra-pool 0-mm same-strain deduplication
        def _dedup_pool(by_iso: Dict[str, List[Dict[str, object]]]) -> Set[str]:
            retired_here: Set[str] = set()
            for iso, grp in sorted(by_iso.items()):
                if len(grp) < 2:
                    continue
                ordered = sorted(grp, key=lambda item: candidate_winner_sort_key(item, prefer_source))
                kept_winners: List[Dict[str, object]] = []
                for cand in ordered:
                    c_acc = str(cand["accession"])
                    c_key = (c_acc, cseg)
                    c_seq = aln_by_key.get(c_key, "")
                    matched_w = None
                    if not cand["protected"] and c_seq:
                        for w in kept_winners:
                            if not passages_compatible(str(w["passage"]), str(cand["passage"])):
                                continue
                            w_key = (str(w["accession"]), cseg)
                            ok, _, _ = aligned_sequences_agree(
                                aln_by_key[w_key],
                                c_seq,
                                int(w["real_length"]),
                                int(cand["real_length"]),
                                max_mm=0,
                            )
                            if ok:
                                matched_w = w
                                break
                    if matched_w is None:
                        kept_winners.append(cand)
                        continue

                    w_acc = str(matched_w["accession"])
                    w_key = (w_acc, cseg)
                    upd = coalesce_winner_metadata(raw_by_key[w_key], raw_by_key[c_key])
                    if upd:
                        raw_by_key[w_key].update(upd)
                        if matched_w["origin"] == "batch":
                            winner_meta_updates[w_key].update(upd)
                    retired_here.add(c_acc)
                    reason = f"intra_db_duplicate_of:{w_acc}"
                    if cand["origin"] == "batch":
                        retired_batch_keys[(c_acc, str(cand["orig_segment"]))] = (w_acc, reason)
                    else:
                        retired_update_ids.add(c_acc)
                    audit_pairs.append({
                        "segment": cseg or "1",
                        "retired_accession": c_acc,
                        "retired_origin": str(cand["origin"]),
                        "retired_source": str(cand["source"]),
                        "winner_accession": w_acc,
                        "winner_origin": str(matched_w["origin"]),
                        "winner_source": str(matched_w["source"]),
                        "match_tier": "intra_db_same_strain",
                        "canonical_isolate": iso,
                    })
            return retired_here

        # Step 1: Intra-preferred pool (GenBank+GenBank)
        retired_pref = _dedup_pool(pref_by_iso)
        if retired_pref:
            redirect: Dict[str, str] = {p["retired_accession"]: p["winner_accession"] for p in audit_pairs if p["segment"] == (cseg or "1")}
            for bare_k, w_info in list(pref_by_bare.items()):
                curr = str(w_info["accession"])
                while curr in redirect:
                    curr = redirect[curr]
                if curr != str(w_info["accession"]):
                    for cand in seg_cands:
                        if str(cand["accession"]) == curr:
                            pref_by_bare[bare_k] = cand
                            break
            for iso in list(pref_by_iso.keys()):
                pref_by_iso[iso] = [w for w in pref_by_iso[iso] if str(w["accession"]) not in retired_pref]

        # Step 2: Cross-DB deduplication (Option B symmetric quality competition, BOTH arrival orders!)
        retired_ext: Set[str] = set()
        for ext in sorted(externals, key=lambda item: candidate_winner_sort_key(item, prefer_source)):
            if ext["protected"]:
                continue
            e_acc = str(ext["accession"])
            e_key = (e_acc, cseg)
            e_seq = aln_by_key.get(e_key, "")
            if not e_seq:
                continue
            e_row = raw_by_key[e_key]
            links = extract_insdc_links(e_row, segment_hint=cseg)

            matched_w = None
            ext_wins = False
            match_tier = ""
            for lk in links:
                w = pref_by_bare.get(lk)
                if w is None or str(w["accession"]) in retired_pref:
                    continue
                if not passages_compatible(str(w["passage"]), str(ext["passage"])):
                    continue
                w_key = (str(w["accession"]), cseg)
                w_seq = aln_by_key[w_key]
                cand_ext_wins = (not w["protected"]) and is_external_quality_upgrade(w_seq, e_seq, is_aligned=True)
                if cand_ext_wins:
                    ok, _, _ = aligned_sequences_agree(
                        e_seq, w_seq, int(ext["real_length"]), int(w["real_length"]), max_mm=2
                    )
                else:
                    ok, _, _ = aligned_sequences_agree(
                        w_seq, e_seq, int(w["real_length"]), int(ext["real_length"]), max_mm=2
                    )
                if ok:
                    matched_w = w
                    ext_wins = cand_ext_wins
                    match_tier = "tier1_insdc_link"
                    break

            if matched_w is None and ext["isolate"] and ext["isolate"] in pref_by_iso:
                candidates = sorted(
                    pref_by_iso[str(ext["isolate"])],
                    key=lambda item: candidate_winner_sort_key(item, prefer_source),
                )
                for w in candidates:
                    if str(w["accession"]) in retired_pref:
                        continue
                    if not passages_compatible(str(w["passage"]), str(ext["passage"])):
                        continue
                    w_key = (str(w["accession"]), cseg)
                    w_seq = aln_by_key[w_key]
                    cand_ext_wins = (not w["protected"]) and is_external_quality_upgrade(w_seq, e_seq, is_aligned=True)
                    if cand_ext_wins:
                        ok, _, _ = aligned_sequences_agree(
                            e_seq, w_seq, int(ext["real_length"]), int(w["real_length"]), max_mm=0
                        )
                    else:
                        ok, _, _ = aligned_sequences_agree(
                            w_seq, e_seq, int(w["real_length"]), int(ext["real_length"]), max_mm=0
                        )
                    if ok:
                        matched_w = w
                        ext_wins = cand_ext_wins
                        match_tier = "tier2_exact_same_strain"
                        break

            if matched_w is None:
                continue

            if ext_wins:
                winner_cand, loser_cand = ext, matched_w
                retired_pref.add(str(matched_w["accession"]))
            else:
                winner_cand, loser_cand = matched_w, ext
                retired_ext.add(e_acc)

            win_acc = str(winner_cand["accession"])
            win_key = (win_acc, cseg)
            los_acc = str(loser_cand["accession"])
            los_key = (los_acc, cseg)
            upd = coalesce_winner_metadata(
                raw_by_key[win_key],
                raw_by_key[los_key],
                retired_source=str(loser_cand["source"]),
            )
            upd.setdefault("sequence_source", str(winner_cand.get("sequence_source") or primary_sequence_source(str(winner_cand["source"]))))
            if upd:
                raw_by_key[win_key].update(upd)
                if "dataset_source" in upd:
                    winner_cand["source"] = upd["dataset_source"]
                if winner_cand["origin"] == "batch":
                    winner_meta_updates[win_key].update(upd)
            reason = f"cross_db_duplicate_of:{win_acc}"
            if loser_cand["origin"] == "batch":
                retired_batch_keys[(los_acc, str(loser_cand["orig_segment"]))] = (win_acc, reason)
            else:
                retired_update_ids.add(los_acc)
            audit_pairs.append({
                "segment": cseg or "1",
                "retired_accession": los_acc,
                "retired_origin": str(loser_cand["origin"]),
                "retired_source": str(loser_cand["source"]),
                "winner_accession": win_acc,
                "winner_origin": str(winner_cand["origin"]),
                "winner_source": str(winner_cand["source"]),
                "sequence_source": str(winner_cand.get("sequence_source") or ""),
                "match_tier": match_tier,
                "canonical_isolate": str(ext["isolate"] or matched_w["isolate"]),
            })

        # Step 3: Intra-external pool (e.g. GISAID+GISAID across update_db + batch)
        surviving_ext_by_iso = {
            iso: [e for e in grp if str(e["accession"]) not in retired_ext]
            for iso, grp in ext_by_iso.items()
        }
        _dedup_pool(surviving_ext_by_iso)

    # Apply metadata coalescing to batch winners
    for (w_acc, cseg), upd in winner_meta_updates.items():
        for (b_acc, orig_seg), indices in batch_rows_by_key.items():
            if b_acc == w_acc and _canon_seg(orig_seg) == cseg:
                for col, val in upd.items():
                    if col not in df.columns:
                        df[col] = ""
                    for idx in indices:
                        df.at[idx, col] = val

    # Mark retired batch rows excluded in meta_data
    retired_msa_pairs: Set[Tuple[str, str]] = set()
    if retired_batch_keys and "canonical_accession" not in df.columns:
        df["canonical_accession"] = ""
    for (r_acc, orig_seg), (w_acc, reason) in retired_batch_keys.items():
        retired_msa_pairs.add((r_acc, _canon_seg(orig_seg)))
        for idx in batch_rows_by_key.get((r_acc, orig_seg), ()):
            df.at[idx, "exclusion_status"] = "1"
            existing_crit = str(df.at[idx, "exclusion_criteria"] or "").strip()
            if not existing_crit:
                df.at[idx, "exclusion_criteria"] = reason
            df.at[idx, "canonical_accession"] = w_acc

    # Write clean MSA files omitting retired batch query records (never omitting protected refs)
    os.makedirs(out_msa_dir, exist_ok=True)
    for msa_path, file_seg, records in msa_records_per_file:
        c_file_seg = _canon_seg(file_seg)
        kept_records = []
        for rec in records:
            rid = rec.id.strip()
            acc = normalise_accession(rid) or rid
            if acc in protected_ids or rid in protected_ids:
                kept_records.append(rec)
                continue
            segs_for_acc = acc_to_segs.get(acc, set())
            rec_seg = c_file_seg
            if not rec_seg and len(segs_for_acc) == 1:
                rec_seg = _canon_seg(next(iter(segs_for_acc)))
            if (acc, rec_seg) in retired_msa_pairs:
                continue
            kept_records.append(rec)
        out_path = out_msa_dir / msa_path.name
        with open(out_path, "w", encoding="utf-8") as out_fh:
            SeqIO.write(kept_records, out_fh, "fasta")

    return df, retired_update_ids, audit_pairs


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Strain-aware sequence deduplication on padded MSAs.")
    parser.add_argument("--meta_data", required=True, help="Input metadata TSV")
    parser.add_argument("--msa", nargs="+", required=True, help="Input padded alignment FASTA file(s)")
    parser.add_argument("--ref_list", default=None, help="Optional reference list TSV/TXT")
    parser.add_argument("--update_db", default=None, help="Optional existing SQLite DB for update mode")
    parser.add_argument("--prefer_source", default="genbank", help="Preferred repository source (default: genbank)")
    parser.add_argument("--out_meta", required=True, help="Output deduplicated metadata TSV")
    parser.add_argument("--out_msa_dir", required=True, help="Output directory for clean padded alignment FASTAs")
    parser.add_argument(
        "--out_retired_update_ids",
        required=True,
        help="Output file listing update_db accessions superseded by incoming preferred-source sequences",
    )
    parser.add_argument("--report_tsv", default=None, help="Optional TSV audit report of retired duplicates")
    args = parser.parse_args(argv)

    df_meta = pd.read_csv(args.meta_data, sep="\t", dtype=str, keep_default_na=False)
    msa_paths = [Path(p) for p in args.msa]
    out_msa_dir = Path(args.out_msa_dir)

    dedup_df, retired_upd_ids, pairs = deduplicate_strain_alignments(
        df_meta=df_meta,
        msa_files=msa_paths,
        out_msa_dir=out_msa_dir,
        ref_list_path=args.ref_list,
        update_db_path=args.update_db,
        prefer_source=args.prefer_source,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.out_meta)), exist_ok=True)
    dedup_df.to_csv(args.out_meta, sep="\t", index=False)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_retired_update_ids)), exist_ok=True)
    with open(args.out_retired_update_ids, "w", encoding="utf-8") as fh:
        for acc in sorted(retired_upd_ids):
            fh.write(f"{acc}\n")

    if args.report_tsv:
        os.makedirs(os.path.dirname(os.path.abspath(args.report_tsv)), exist_ok=True)
        fieldnames = [
            "segment", "retired_accession", "retired_origin", "retired_source",
            "winner_accession", "winner_origin", "winner_source", "sequence_source",
            "match_tier", "canonical_isolate",
        ]
        with open(args.report_tsv, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(pairs)

    print(
        f"[DedupStrainAlignment] Retired {len(pairs)} duplicate strain segment(s) "
        f"({len(pairs) - len(retired_upd_ids)} incoming batch, {len(retired_upd_ids)} superseded in update_db)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
