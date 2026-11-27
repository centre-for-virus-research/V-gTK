#!/usr/bin/env python3
"""Segmented virus strain-name linkage table builder for V-gTK.

Constructs a one-row-per-unique-strain SQLite table (`strain_segment_linkage`)
with one SQL-safe column per expected segment (`segment_1`..`segment_8` for IAV,
or `segment_L`, `segment_M`, `segment_S` for Bunyavirales / Arenaviridae, etc.)
pointing at a single canonical active sequence accession per `(strain, segment)`.

Handles real-world AIV/IAV strain naming heterogeneity:
  - Case differences, 2- vs 4-digit years, US state & country abbreviations
    (`NY` <-> `New_York`, `Viet_Nam` <-> `Vietnam`), double slashes (`//`),
    leading wrappers (`_A/`, `(A/`, `partial/A/`, `Recombinant_A/`, `RG-A/`),
    trailing subtype/segment/passage tags (`(H5N1)`, `_H7N2NSB`, `-HA`,
    `-original`), and omitted `A/` prefixes (`Guinea_Fowl/...`).
  - Reunites `Isolate_Id` (`EPI_ISL_*`) splits caused by per-segment intra-DB
    dedup via union-find across explicit `alt_accessions` `EPI_ISL_*` links
    (e.g. `A/ABUDHABI/2373/2025`), while keeping biologically distinct passaged
    genomes (`A/AICHI/2/2018` `MDCK 1` vs `Original`) as separate rows.
  - Attaches `Isolate_Id`-less GenBank rows onto matching canonical strain +
    compatible passage groups when no segment conflict exists.
  - Selects one deterministic canonical accession per `(strain_accession, segment)`
    cell and records any multi-candidate cells in `ambiguous_segments` plus the
    companion `strain_segment_alternates` table so no active sequence is invisible.
"""

from __future__ import annotations

import csv
import os
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import segment_utils
from cross_db_dedup import (
    candidate_winner_sort_key,
    classify_passage,
    extract_passage_class,
    merge_dataset_sources,
    passages_compatible,
)

EPI_ISL_RE = re.compile(r"EPI_ISL_\d+", re.IGNORECASE)
BIOSAMPLE_RE = re.compile(r"(?:BioSample:|SAM[NED][A-Z]?\d+)", re.IGNORECASE)

_PLACEHOLDERS = frozenset(
    {"", "na", "n/a", "none", "null", "nan", "unknown", "unk", "unassigned", "missing", "-", "."}
)
_LEADING_WRAPPER_RE = re.compile(
    r"^(?:[\(\[_]+|partial/|recombinant[_\s/-]+|reassortant[_\s/-]+|rg[-_\s]+)+",
    re.IGNORECASE,
)
_TRAILING_SUBTYPE_PAREN_RE = re.compile(
    r"\s*\(\s*(?:H\d{1,2}(?:N\d{1,2}|\?|x)?|H\?N\d{1,2}|N\d{1,2}|mixed)\s*\)\s*$",
    re.IGNORECASE,
)
_TRAILING_TAG_RE = re.compile(
    r"(?:[_-](?:H\d{1,2}N\d{1,2}[A-Z0-9]*|PB2|PB1|PA|HA|NP|NA|MP|M1|M2|NS|NS1|NEP|original|orig|egg|e\d+|mdck\d*|siat\d*|clinical|swab))+$",
    re.IGNORECASE,
)
_FLU_TYPE_PREFIX_RE = re.compile(r"^[abcd]/", re.IGNORECASE)


def load_location_aliases(path: Optional[str]) -> Dict[str, str]:
    """Load optional 2-column TSV (`alias<TAB>canonical`) for geographic normalisation."""
    aliases: Dict[str, str] = {}
    if not path or not os.path.isfile(path):
        return aliases
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return aliases
        # Allow headerless or titled files
        if header and header[0].strip().lower() not in ("alias", "raw", "from", "variant"):
            if len(header) >= 2 and header[0].strip() and header[1].strip():
                aliases[re.sub(r"[\s_-]+", "_", header[0].strip().casefold())] = re.sub(
                    r"[\s_-]+", "_", header[1].strip().casefold()
                )
        for row in reader:
            if len(row) < 2:
                continue
            k = re.sub(r"[\s_-]+", "_", row[0].strip().casefold())
            v = re.sub(r"[\s_-]+", "_", row[1].strip().casefold())
            if k and v:
                aliases[k] = v
    return aliases


def normalise_segmented_strain_key(
    raw: Optional[str],
    collection_year: Optional[str] = None,
    location_aliases: Optional[Mapping[str, str]] = None,
) -> str:
    """Canonicalise a segmented-virus strain name (AIV/IAV and general segmented viruses)."""
    if raw is None:
        return ""
    text = str(raw).strip().strip("'\"")
    if not text:
        return ""
    if text.casefold() in _PLACEHOLDERS:
        return ""

    # Strip trailing parenthesised subtypes FIRST before trimming wrapping parens/brackets,
    # and only strip leading wrappers like _A/, partial/A/, Recombinant_A/, RG-A/.
    prev = None
    while prev != text:
        prev = text
        text = _TRAILING_SUBTYPE_PAREN_RE.sub("", text).strip()
        if (text.startswith("(") and text.endswith(")")) or (text.startswith("[") and text.endswith("]")):
            text = text[1:-1].strip()
        text = _LEADING_WRAPPER_RE.sub("", text).strip()

    # Collapse multiple slashes and surrounding whitespace
    text = re.sub(r"\s*/+\s*", "/", text)
    text = text.strip(" _/;:,.")
    if not text or text.casefold() in _PLACEHOLDERS:
        return ""

    # Strip trailing segment / subtype / passage suffixes attached with '_' or '-' (e.g. _H7N2NSB, -HA, -original)
    parts = text.split("/")
    if parts:
        parts[-1] = _TRAILING_TAG_RE.sub("", parts[-1]).strip()
        # Re-strip parenthesised subtype if revealed
        parts[-1] = _TRAILING_SUBTYPE_PAREN_RE.sub("", parts[-1]).strip()

    # Normalise spaces/underscores per slash-separated field
    norm_parts: List[str] = []
    for p in parts:
        tok = re.sub(r"[\s_]+", "_", p.strip()).strip("_")
        if tok:
            norm_parts.append(tok.casefold())
    if not norm_parts:
        return ""

    lowered = "/".join(norm_parts)
    if lowered in _PLACEHOLDERS or len(lowered) < 3:
        return ""

    # Apply location alias normalisation on location field(s)
    if location_aliases and len(norm_parts) >= 3:
        # In flu names: 4-field = A/Loc/ID/Year (loc idx 1); 5+-field = A/Host/Loc/ID/Year (loc idx 2)
        # Missing A/ prefix: 4-field = Host/Loc/ID/Year (loc idx 1)
        loc_indices = []
        if _FLU_TYPE_PREFIX_RE.match(lowered):
            if len(norm_parts) == 4:
                loc_indices = [1]
            elif len(norm_parts) >= 5:
                loc_indices = [2, 1]
        elif len(norm_parts) >= 4:
            loc_indices = [1]
        for idx in loc_indices:
            if idx < len(norm_parts) and norm_parts[idx] in location_aliases:
                norm_parts[idx] = location_aliases[norm_parts[idx]]
                break

    # Normalise trailing 2-digit year to 4-digit year for flu-like slash-delimited strains
    if len(norm_parts) >= 4:
        tail = norm_parts[-1]
        if len(tail) == 2 and tail.isdigit():
            yr = str(collection_year or "").strip()
            if len(yr) == 4 and yr.isdigit() and yr.endswith(tail):
                norm_parts[-1] = yr
            else:
                two = int(tail)
                norm_parts[-1] = f"20{tail}" if two <= 35 else f"19{tail}"

    lowered = "/".join(norm_parts)
    if lowered in _PLACEHOLDERS or len(lowered) < 3:
        return ""
    return lowered


def segment_column_name(seg_label: str) -> str:
    """Return a clean SQL-safe column identifier for a segment label (e.g. ``segment_1``, ``segment_L``)."""
    clean = re.sub(r"[^A-Za-z0-9]+", "_", str(seg_label).strip()).strip("_")
    return f"segment_{clean or 'unknown'}"


def natural_segment_order(segments: Iterable[str]) -> List[str]:
    """Sort segments naturally: numeric ascending first, then alphabetic."""
    unique = []
    seen = set()
    for s in segments:
        norm = segment_utils.normalise_segment(s)
        if not norm or norm.casefold() in segment_utils.PANDAS_NULL_TOKENS:
            continue
        if norm not in seen:
            seen.add(norm)
            unique.append(norm)

    def _key(s: str):
        return (0, int(s), s) if s.isdigit() else (1, 0, s)

    return sorted(unique, key=_key)


def _is_excluded_row(row: Mapping[str, object]) -> bool:
    es = str(row.get("exclusion_status") or "").strip().lower()
    if es not in ("", "0", "false", "no", "na", "none", "nan"):
        return True
    ec = str(row.get("exclusion_criteria") or "").strip()
    if ec and ec.lower() not in ("na", "none", "null", "nan"):
        return True
    return False


def _extract_primary_strain_display(row: Mapping[str, object]) -> str:
    for col in ("Parsed_strain", "strain", "isolate", "Isolate_Name", "gisaid_isolate_name"):
        val = str(row.get(col) or "").strip()
        if val and val.casefold() not in _PLACEHOLDERS:
            return val
    return ""


def _extract_biosample(row: Mapping[str, object]) -> str:
    for col in ("biosample", "BioSample", "db_xref"):
        val = str(row.get(col) or "").strip()
        if not val:
            continue
        m = re.search(r"SAM[NED][A-Z]?\d+", val, re.IGNORECASE)
        if m:
            return m.group(0).upper()
    return ""


def _extract_year_hint(row: Mapping[str, object]) -> str:
    yr = str(row.get("collection_year") or "").strip()
    if len(yr) == 4 and yr.isdigit():
        return yr
    cdate = str(row.get("collection_date") or "").strip()
    if len(cdate) >= 4 and cdate[:4].isdigit():
        return cdate[:4]
    return ""


def _row_winner_key(row: Mapping[str, object]):
    acc = str(row.get("primary_accession") or "").strip()
    acc_type = str(row.get("accession_type") or "").strip().lower()
    try:
        rlen = int(float(row.get("real_length") or 0))
    except (ValueError, TypeError):
        rlen = 0
    return candidate_winner_sort_key(
        {
            "accession": acc,
            "protected": acc_type in ("master", "reference"),
            "source": str(row.get("dataset_source") or "genbank"),
            "origin": "batch",
            "real_length": rlen,
            "collection_date": str(row.get("collection_date") or ""),
        },
        prefer_source="genbank",
        quality_first=True,
    )


def build_strain_segment_linkage(
    rows: Iterable[Mapping[str, object]],
    expected_segments: Optional[Sequence[str]] = None,
    location_aliases: Optional[Mapping[str, str]] = None,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[str]]:
    """Build deterministic one-row-per-strain linkage records from active metadata rows.

    Returns ``(linkage_rows, alternate_rows, segment_order)``.
    """
    active: List[Dict[str, object]] = []
    observed_segs: Set[str] = set()
    for raw in rows:
        if _is_excluded_row(raw):
            continue
        acc = str(raw.get("primary_accession") or "").strip()
        if not acc:
            continue
        seg = segment_utils.normalise_segment(raw.get("segment") or raw.get("segment_validated"))
        if not seg or seg.casefold() in segment_utils.PANDAS_NULL_TOKENS:
            continue
        observed_segs.add(seg)
        display_strain = _extract_primary_strain_display(raw)
        year_hint = _extract_year_hint(raw)
        canon_key = normalise_segmented_strain_key(
            display_strain, collection_year=year_hint, location_aliases=location_aliases
        )
        iso_id_raw = str(raw.get("Isolate_Id") or raw.get("gisaid_isolate_id") or "").strip()
        iso_match = EPI_ISL_RE.search(iso_id_raw)
        iso_id = iso_match.group(0).upper() if iso_match else (iso_id_raw if iso_id_raw and iso_id_raw.casefold() not in _PLACEHOLDERS else "")
        alt_text = str(raw.get("alt_accessions") or "")
        alt_isos = sorted({m.upper() for m in EPI_ISL_RE.findall(alt_text)})
        biosample = _extract_biosample(raw)
        passage = extract_passage_class(raw)
        active.append(
            {
                "raw": raw,
                "acc": acc,
                "seg": seg,
                "display_strain": display_strain,
                "canon_key": canon_key,
                "iso_id": iso_id,
                "alt_isos": alt_isos,
                "biosample": biosample,
                "passage": passage,
                "source": str(raw.get("dataset_source") or "genbank"),
            }
        )

    seg_order = natural_segment_order(expected_segments if expected_segments else observed_segs)
    if expected_segments:
        # Ensure any extra observed segments aren't silently lost if unexpected
        seg_order = natural_segment_order(list(expected_segments))
    seg_set = set(seg_order)

    # Pass 0: Recover omitted 'A/' prefix when 'host/loc/id/year' matches an existing 'a/host/loc/id/year'
    flu_a_keys = {r["canon_key"] for r in active if r["canon_key"].startswith("a/")}
    for r in active:
        ck = r["canon_key"]
        if ck and not _FLU_TYPE_PREFIX_RE.match(ck) and f"a/{ck}" in flu_a_keys:
            r["canon_key"] = f"a/{ck}"

    # Pass 1: Union-find across explicit Isolate_Id <-> alt_accessions EPI_ISL_* equivalence links
    # (Reunites intra-DB split isolates like A/ABUDHABI/2373/2025 without merging non-identical passages like A/AICHI/2/2018)
    uf_parent: Dict[str, str] = {}

    def _uf_find(x: str) -> str:
        while uf_parent.get(x, x) != x:
            uf_parent[x] = uf_parent.get(uf_parent[x], uf_parent[x])
            x = uf_parent[x]
        return x

    def _uf_union(a: str, b: str) -> None:
        ra, rb = _uf_find(a), _uf_find(b)
        if ra != rb:
            if ra < rb:
                uf_parent[rb] = ra
            else:
                uf_parent[ra] = rb

    for r in active:
        iso = r["iso_id"]
        if iso:
            uf_parent.setdefault(iso, iso)
            for alt_iso in r["alt_isos"]:
                uf_parent.setdefault(alt_iso, alt_iso)
                _uf_union(iso, alt_iso)

    # Group rows that have an explicit Isolate_Id or BioSample first
    groups: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    group_method: Dict[str, str] = {}
    group_iso_members: Dict[str, Set[str]] = defaultdict(set)
    unlinked: List[Dict[str, object]] = []

    for r in active:
        if r["seg"] not in seg_set:
            continue
        if r["iso_id"]:
            root = _uf_find(r["iso_id"])
            gid = f"iso:{root}"
            groups[gid].append(r)
            group_iso_members[gid].add(r["iso_id"])
            for alt_iso in r["alt_isos"]:
                if _uf_find(alt_iso) == root:
                    group_iso_members[gid].add(alt_iso)
            group_method[gid] = "isolate_id_merged" if len(group_iso_members[gid]) > 1 else "isolate_id"
        elif r["biosample"]:
            gid = f"biosample:{r['biosample']}"
            groups[gid].append(r)
            group_method[gid] = "biosample"
        else:
            unlinked.append(r)

    # Index isolate_id / biosample groups by canonical strain key so Isolate_Id-less GenBank rows
    # can attach cleanly when there is a single compatible group with no segment conflict.
    groups_by_ck: Dict[str, List[str]] = defaultdict(list)
    group_segs: Dict[str, Set[str]] = defaultdict(set)
    group_passages: Dict[str, Set[str]] = defaultdict(set)
    for gid, members in groups.items():
        cks = {m["canon_key"] for m in members if m["canon_key"]}
        for ck in cks:
            groups_by_ck[ck].append(gid)
        for m in members:
            group_segs[gid].add(m["seg"])
            if m["passage"] != "unknown":
                group_passages[gid].add(m["passage"])

    # Pass 2: Attach Isolate_Id-less rows to matching isolate_id group if unique & non-conflicting,
    # else group by canonical_strain_key (+ passage class when multiple distinct passages share a strain).
    strain_only_buckets: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for r in unlinked:
        ck = r["canon_key"]
        if not ck:
            gid = f"acc:{r['acc']}"
            groups[gid].append(r)
            group_method[gid] = "singleton_fallback"
            continue

        compatible_gids = []
        for gid in groups_by_ck.get(ck, ()):
            if r["seg"] in group_segs[gid]:
                continue
            g_pas = group_passages[gid]
            if not g_pas or r["passage"] == "unknown" or any(passages_compatible(gp, r["passage"]) for gp in g_pas):
                compatible_gids.append(gid)

        if len(compatible_gids) == 1:
            target_gid = compatible_gids[0]
            groups[target_gid].append(r)
            group_segs[target_gid].add(r["seg"])
            if r["passage"] != "unknown":
                group_passages[target_gid].add(r["passage"])
        else:
            strain_only_buckets[ck].append(r)

    for ck, members in strain_only_buckets.items():
        distinct_known_passages = sorted({m["passage"] for m in members if m["passage"] != "unknown"})
        if len(distinct_known_passages) > 1:
            # Split by passage tier if non-compatible passages exist under the same strain name
            for m in members:
                p = m["passage"]
                gid = f"strain:{ck}|{p}"
                groups[gid].append(m)
                group_method[gid] = "canonical_strain_passage"
        else:
            gid = f"strain:{ck}"
            groups[gid].extend(members)
            group_method[gid] = "canonical_strain"

    # Pass 3: Construct final linkage rows + deterministic single-winner selection per segment cell
    linkage_rows: List[Dict[str, object]] = []
    alternate_rows: List[Dict[str, object]] = []

    for gid in sorted(groups.keys()):
        members = groups[gid]
        iso_members = sorted(group_iso_members.get(gid, ()))
        # Choose primary isolate_id: most frequent among members, tie-broken lexicographically
        if iso_members:
            member_iso_counts = Counter(m["iso_id"] for m in members if m["iso_id"])
            primary_iso = sorted(
                member_iso_counts.keys() or iso_members,
                key=lambda x: (-member_iso_counts.get(x, 0), x),
            )[0]
            alt_isos = [x for x in iso_members if x != primary_iso]
            strain_accession = primary_iso
        else:
            primary_iso = ""
            alt_isos = []
            strain_accession = gid

        # Pick best human-readable strain_name and canonical_strain_key
        ck_counts = Counter(m["canon_key"] for m in members if m["canon_key"])
        canonical_strain_key = sorted(ck_counts.keys(), key=lambda x: (-ck_counts[x], x))[0] if ck_counts else ""
        disp_counts = Counter(m["display_strain"] for m in members if m["display_strain"])
        strain_name = (
            sorted(disp_counts.keys(), key=lambda x: (-disp_counts[x], -len(x), x))[0]
            if disp_counts
            else ""
        )

        passages = sorted({m["passage"] for m in members if m["passage"] != "unknown"})
        if not passages:
            passage_class = "unknown"
        elif len(passages) == 1:
            passage_class = passages[0]
        else:
            passage_class = "mixed"

        sources_union = merge_dataset_sources(*(m["source"] for m in members))

        by_seg: Dict[str, List[Dict[str, object]]] = defaultdict(list)
        for m in members:
            by_seg[m["seg"]].append(m)

        rec: Dict[str, object] = {
            "strain_accession": strain_accession,
            "strain_name": strain_name,
            "canonical_strain_key": canonical_strain_key,
            "isolate_id": primary_iso,
            "alt_isolate_ids": ";".join(alt_isos),
            "passage_class": passage_class,
        }

        present_count = 0
        ambiguous_segs: List[str] = []
        for seg in seg_order:
            col = segment_column_name(seg)
            cands = by_seg.get(seg, [])
            if not cands:
                rec[col] = None
                continue
            present_count += 1
            if len(cands) == 1:
                rec[col] = cands[0]["acc"]
            else:
                ordered = sorted(cands, key=lambda item: _row_winner_key(item["raw"]))
                winner = ordered[0]["acc"]
                rec[col] = winner
                ambiguous_segs.append(seg)
                for alt in ordered[1:]:
                    alternate_rows.append(
                        {
                            "strain_accession": strain_accession,
                            "segment": seg,
                            "selected_accession": winner,
                            "alternate_accession": alt["acc"],
                            "dataset_source": alt["source"],
                        }
                    )

        rec["segments_present"] = present_count
        rec["complete_status"] = "Complete" if (seg_order and present_count == len(seg_order)) else "Incomplete"
        rec["dataset_sources"] = sources_union
        rec["linkage_method"] = group_method.get(gid, "canonical_strain")
        rec["ambiguous_segments"] = ";".join(ambiguous_segs)
        linkage_rows.append(rec)

    return linkage_rows, alternate_rows, seg_order


def write_strain_linkage_tables(
    conn: sqlite3.Connection,
    expected_segments: Optional[Sequence[str]] = None,
    location_aliases: Optional[Mapping[str, str]] = None,
    summary_tsv_path: Optional[str] = None,
) -> Dict[str, object]:
    """Read active kept rows from ``meta_data``, rebuild ``strain_segment_linkage`` and
    ``strain_segment_alternates`` inside ``conn``, and return summary metrics."""
    meta_cols = [r[1] for r in conn.execute("PRAGMA table_info(meta_data)")]
    wanted = [
        c
        for c in (
            "primary_accession",
            "accession_type",
            "segment",
            "segment_validated",
            "Parsed_strain",
            "strain",
            "isolate",
            "Isolate_Name",
            "gisaid_isolate_name",
            "Isolate_Id",
            "gisaid_isolate_id",
            "alt_accessions",
            "dataset_source",
            "Passage_History",
            "passage_history",
            "gisaid_passage_history",
            "isolation_source",
            "collection_year",
            "collection_date",
            "real_length",
            "db_xref",
            "exclusion_status",
            "exclusion_criteria",
        )
        if c in meta_cols
    ]
    if "primary_accession" not in wanted:
        raise ValueError("meta_data lacks primary_accession column")

    qcols = ", ".join(f'"{c}"' for c in wanted)
    where_parts = []
    if "exclusion_status" in wanted:
        where_parts.append("LOWER(TRIM(COALESCE(CAST(exclusion_status AS TEXT), ''))) IN ('', '0', 'false', 'no', 'na', 'none', 'nan')")
    if "exclusion_criteria" in wanted:
        where_parts.append("TRIM(COALESCE(CAST(exclusion_criteria AS TEXT), '')) = ''")
    sql = f"SELECT {qcols} FROM meta_data"
    if where_parts:
        sql += " WHERE " + " AND ".join(where_parts)

    def _row_iter():
        cur = conn.execute(sql)
        try:
            for tup in cur:
                yield dict(zip(wanted, tup))
        finally:
            cur.close()

    linkage_rows, alternate_rows, seg_order = build_strain_segment_linkage(
        _row_iter(),
        expected_segments=expected_segments,
        location_aliases=location_aliases,
    )

    seg_cols = [segment_column_name(s) for s in seg_order]
    all_cols = (
        [
            "strain_accession",
            "strain_name",
            "canonical_strain_key",
            "isolate_id",
            "alt_isolate_ids",
            "passage_class",
        ]
        + seg_cols
        + [
            "segments_present",
            "complete_status",
            "dataset_sources",
            "linkage_method",
            "ambiguous_segments",
        ]
    )

    conn.execute("DROP TABLE IF EXISTS strain_segment_linkage")
    conn.execute("DROP TABLE IF EXISTS strain_segment_alternates")

    col_defs = []
    for c in all_cols:
        if c == "strain_accession":
            col_defs.append('"strain_accession" TEXT PRIMARY KEY')
        elif c == "segments_present":
            col_defs.append('"segments_present" INTEGER')
        else:
            col_defs.append(f'"{c}" TEXT')
    conn.execute(f"CREATE TABLE strain_segment_linkage ({', '.join(col_defs)})")
    conn.execute(
        "CREATE TABLE strain_segment_alternates ("
        '"strain_accession" TEXT, "segment" TEXT, "selected_accession" TEXT, '
        '"alternate_accession" TEXT, "dataset_source" TEXT)'
    )

    if linkage_rows:
        ph = ", ".join("?" for _ in all_cols)
        qcol_list = ", ".join(f'"{c}"' for c in all_cols)
        batch = []
        for row in linkage_rows:
            batch.append(tuple(row.get(c) for c in all_cols))
            if len(batch) >= 20000:
                conn.executemany(f"INSERT INTO strain_segment_linkage ({qcol_list}) VALUES ({ph})", batch)
                batch.clear()
        if batch:
            conn.executemany(f"INSERT INTO strain_segment_linkage ({qcol_list}) VALUES ({ph})", batch)

    if alternate_rows:
        conn.executemany(
            "INSERT INTO strain_segment_alternates "
            "(strain_accession, segment, selected_accession, alternate_accession, dataset_source) "
            "VALUES (?, ?, ?, ?, ?)",
            [
                (
                    r["strain_accession"],
                    r["segment"],
                    r["selected_accession"],
                    r["alternate_accession"],
                    r["dataset_source"],
                )
                for r in alternate_rows
            ],
        )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_ssl_canonical_key ON strain_segment_linkage (canonical_strain_key)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_ssl_isolate_id ON strain_segment_linkage (isolate_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_ssl_complete ON strain_segment_linkage (complete_status)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_ssa_strain ON strain_segment_alternates (strain_accession)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_ssa_alt_acc ON strain_segment_alternates (alternate_accession)")

    complete_count = sum(1 for r in linkage_rows if r["complete_status"] == "Complete")
    ambig_rows = sum(1 for r in linkage_rows if r["ambiguous_segments"])
    summary = {
        "strain_rows": len(linkage_rows),
        "complete_strains": complete_count,
        "incomplete_strains": len(linkage_rows) - complete_count,
        "ambiguous_strains": ambig_rows,
        "alternate_accessions": len(alternate_rows),
        "segments": seg_order,
    }

    if summary_tsv_path:
        os.makedirs(os.path.dirname(os.path.abspath(summary_tsv_path)), exist_ok=True)
        with open(summary_tsv_path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["key", "value"])
            for k, v in summary.items():
                w.writerow([k, ",".join(v) if isinstance(v, list) else v])

    return summary
