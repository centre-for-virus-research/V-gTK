"""Cross-database dual-upload deduplication helpers (GenBank + GISAID / Pathoplexus).

Submitting laboratories routinely deposit the same clinical isolate into both
an INSDC member database (GenBank / ENA / DDBJ) and a restricted or specialized
portal (GISAID EpiFlu/EpiCoV/EpiRSV/EpiPox/EpiArbo, or Pathoplexus). When both
feeds are merged without cross-source reconciliation, every dual upload enters
the database twice as two kept rows and a 2-tip UShER polytomy cherry.

This module provides species-agnostic (and influenza-segment-aware) matching,
sequence-verification, passage-compatibility, and metadata-coalescing helpers
shared by:
  1. ``scripts/merge_into_gB_matrix.py`` (build-time and ``--update_db`` merges)
  2. ``dev/deduplicate_cross_db_uploads.py`` (surgical cleanup of extant DBs)
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from accession_utils import normalise_accession
from date_utils import split_date_components
from ValidateStrain import parse_strain_from_definition


FLU_SEG_NUM_TO_GENE = {
    "1": "pb2",
    "2": "pb1",
    "3": "pa",
    "4": "ha",
    "5": "np",
    "6": "na",
    "7": "mp",
    "8": "ns",
}

FLU_GENE_TO_SEG_NUM = {
    "pb2": "1",
    "pb1": "2",
    "pa": "3",
    "ha": "4",
    "np": "5",
    "na": "6",
    "mp": "7",
    "m": "7",
    "ns": "8",
    "he": "he",
    "p3": "p3",
}

GENERIC_INSDC_COLUMN_SUFFIXES = (
    "insdc_accession",
    "insdcaccessionbase",
    "insdc_accession_base",
    "insdc_accession_full",
    "insdcaccessionfull",
    "genbank_accession",
    "insdc_upload",
)

_PLACEHOLDER_STRAINS = frozenset({
    "", "na", "n/a", "none", "null", "unknown", "missing",
    "not applicable", "not collected", "not provided", "unidentified",
})

_SUBTYPE_TRAIL_RE = re.compile(
    r"\s*\(\s*(?:H\d{1,2}(?:N\d{1,2}|\?|x)?|H\?N\d{1,2}|N\d{1,2}|mixed)\s*\)\s*$",
    re.IGNORECASE,
)
_FLU_PREFIX_RE = re.compile(r"^[abcd]/", re.IGNORECASE)
_INSDC_TOKEN_SPLIT_RE = re.compile(r"[,;|/\s]+")
_INSDC_SHAPE_RE = re.compile(r"^[A-Z]{1,3}_?\d{5,9}$", re.IGNORECASE)

_ORIG_PASSAGE_RE = re.compile(
    r"\b(?:original|orig|clinical(?:\s+specimen)?|direct|primary|unpassaged|swab|nasopharyngeal|or|p0)\b",
    re.IGNORECASE,
)
_EGG_PASSAGE_RE = re.compile(
    r"(?:\b(?:egg|eggs|allantoic|amniotic|spf)\b|(?:^|[\s/;,-])e\d+(?:$|[\s/;,-]))",
    re.IGNORECASE,
)
_CELL_PASSAGE_RE = re.compile(
    r"(?:\b(?:mdck|siat\d*|hck|vero|cell|cells|llc-mk2|caco-?2|rhmk)\b|(?:^|[\s/;,-])[csm]\d+(?:$|[\s/;,-]))",
    re.IGNORECASE,
)

_UNAMBIGUOUS_ACGT = frozenset("ACGT")


def parse_fasta_file(fasta_path: Optional[str]) -> Dict[str, str]:
    """Load a FASTA file into ``{bare_or_header_id: sequence}``."""
    seqs: Dict[str, str] = {}
    if not fasta_path:
        return seqs
    header: Optional[str] = None
    chunks: List[str] = []
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(chunks)
                    seqs[header] = seq
                    bare = normalise_accession(header.split()[0])
                    if bare and bare not in seqs:
                        seqs[bare] = seq
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            seq = "".join(chunks)
            seqs[header] = seq
            bare = normalise_accession(header.split()[0])
            if bare and bare not in seqs:
                seqs[bare] = seq
    return seqs


def count_acgt(seq: Optional[str]) -> int:
    """Count unambiguous A/C/G/T bases (case-insensitive)."""
    if not seq or not isinstance(seq, str):
        return 0
    upper = seq.upper()
    return upper.count("A") + upper.count("C") + upper.count("G") + upper.count("T")


def normalize_isolate_key(raw: Optional[str], collection_year: Optional[str] = None) -> str:
    """Canonicalize a viral isolate / strain identifier for cross-DB matching.

    Preserves species generality (works for influenza ``A/host/loc/id/year`` as
    well as unsegmented isolate labels in RABV, HCV, Dengue, Mpox, etc.) while
    refusing generic placeholders so two ``unknown`` rows never match.
    """
    if raw is None:
        return ""
    text = str(raw).strip()
    if not text:
        return ""
    # Repeatedly strip trailing influenza subtype suffixes such as "(H3N2)".
    prev = None
    while prev != text:
        prev = text
        text = _SUBTYPE_TRAIL_RE.sub("", text).strip()

    text = re.sub(r"\s*/\s*", "/", text)
    text = re.sub(r"[\s_]+", "_", text)
    text = text.strip(" _/;:,.")
    lowered = text.casefold()
    if lowered in _PLACEHOLDER_STRAINS or len(lowered) < 4:
        return ""

    # Normalize 2-digit trailing influenza year when 4-digit year is unambiguous.
    if _FLU_PREFIX_RE.match(lowered) and lowered.count("/") >= 3:
        parts = lowered.split("/")
        tail = parts[-1]
        if len(tail) == 2 and tail.isdigit():
            yr = str(collection_year or "").strip()
            if len(yr) == 4 and yr.isdigit() and yr.endswith(tail):
                parts[-1] = yr
            else:
                two = int(tail)
                parts[-1] = f"20{tail}" if two <= 35 else f"19{tail}"
            lowered = "/".join(parts)
    return lowered


def extract_canonical_isolate(row: Dict[str, str]) -> str:
    """Extract and canonicalize the best isolate/strain key from a metadata row."""
    year_hint = str(row.get("collection_year") or "").strip()
    if not year_hint:
        cdate = str(row.get("collection_date") or "").strip()
        if len(cdate) >= 4 and cdate[:4].isdigit():
            year_hint = cdate[:4]

    for col in (
        "Parsed_strain", "strain", "isolate", "Isolate_Name",
        "gisaid_isolate_name", "virus_name", "Virus_Name",
    ):
        val = row.get(col)
        if val is not None and str(val).strip():
            norm = normalize_isolate_key(str(val), collection_year=year_hint)
            if norm:
                return norm

    definition = str(row.get("definition") or "").strip()
    if definition:
        parsed = parse_strain_from_definition(definition)
        if parsed:
            return normalize_isolate_key(parsed, collection_year=year_hint)
    return ""


def classify_passage(text: Optional[str]) -> str:
    """Classify passage history into ``'orig'``, ``'egg'``, ``'cell'``, or ``'unknown'``."""
    if not text:
        return "unknown"
    cleaned = str(text).strip()
    if not cleaned or cleaned.casefold() in _PLACEHOLDER_STRAINS:
        return "unknown"
    has_egg = bool(_EGG_PASSAGE_RE.search(cleaned))
    has_cell = bool(_CELL_PASSAGE_RE.search(cleaned))
    has_orig = bool(_ORIG_PASSAGE_RE.search(cleaned))
    if has_egg and not has_cell:
        return "egg"
    if has_cell and not has_egg:
        return "cell"
    if has_orig and not has_egg and not has_cell:
        return "orig"
    return "unknown"


def extract_passage_class(row: Dict[str, str]) -> str:
    """Extract passage class from standard or vendor-namespaced passage columns."""
    for col in (
        "Passage_History", "passage_history", "gisaid_passage_history",
        "gisaid_passage", "passage", "isolation_source",
    ):
        val = row.get(col)
        if val is not None and str(val).strip():
            cls = classify_passage(str(val))
            if cls != "unknown":
                return cls
    return "unknown"


def passages_compatible(cls_a: str, cls_b: str) -> bool:
    """Return False when two rows represent incompatible laboratory passage tiers
    (e.g. clinical ``orig`` specimen vs ``egg`` or ``cell`` isolate)."""
    if cls_a == "unknown" or cls_b == "unknown":
        return True
    return cls_a == cls_b


def resolve_row_segment_hint(row: Dict[str, str]) -> str:
    """Resolve canonical segment number/token from row fields when available."""
    for col in ("segment", "segment_validated", "segment_declared", "gisaid_segment", "Segment"):
        raw = str(row.get(col) or "").strip()
        if not raw:
            continue
        if raw.isdigit():
            return raw
        mapped = FLU_GENE_TO_SEG_NUM.get(raw.casefold())
        if mapped:
            return mapped
    return ""


def extract_insdc_links(row: Dict[str, str], segment_hint: Optional[str] = None) -> List[str]:
    """Extract bare INSDC accessions referenced by an external (GISAID/Pathoplexus) row.

    When a row belongs to a known influenza segment (e.g. segment ``4`` / ``HA``),
    only that segment's ``*_<gene>_insdc_upload`` column (plus generic single-genome
    INSDC columns) is consulted so multi-segment upload sheets never cross-link
    one segment to another segment's GenBank accession.
    """
    seg = (segment_hint if segment_hint is not None else resolve_row_segment_hint(row)).strip()
    target_gene = FLU_SEG_NUM_TO_GENE.get(seg, seg.casefold() if seg else "")

    candidate_cols: List[str] = []
    for col in row.keys():
        low = col.casefold()
        if low.endswith("_insdc_upload") and low != "insdc_upload":
            # Check if this is a segment-specific GISAID column like gisaid_ha_insdc_upload
            matched_gene = None
            for gene in FLU_GENE_TO_SEG_NUM:
                if low.endswith(f"_{gene}_insdc_upload") or low == f"{gene}_insdc_upload":
                    matched_gene = gene
                    break
            if matched_gene:
                if not target_gene or matched_gene == target_gene:
                    candidate_cols.append(col)
            else:
                candidate_cols.append(col)
        elif any(low == suf or low.endswith(f"_{suf}") for suf in GENERIC_INSDC_COLUMN_SUFFIXES):
            candidate_cols.append(col)

    found: List[str] = []
    seen: Set[str] = set()
    for col in candidate_cols:
        raw = row.get(col)
        if raw is None:
            continue
        text = str(raw).strip()
        if not text or text.casefold() in _PLACEHOLDER_STRAINS:
            continue
        for tok in _INSDC_TOKEN_SPLIT_RE.split(text):
            bare = normalise_accession(tok)
            if not bare:
                continue
            bare_up = bare.upper()
            if _INSDC_SHAPE_RE.match(bare_up) and bare_up not in seen:
                seen.add(bare_up)
                found.append(bare_up)
    return found


def has_coverage_parity(
    winner_real_len: int,
    retired_real_len: int,
    min_ratio: float = 0.94,
    max_shortfall_bp: int = 90,
) -> bool:
    """Ensure the kept sequence is not a short fragment relative to the retired copy.

    Allows up to 90 bp / 6% terminal UTR overhang (e.g. influenza full vRNA 1,778 bp
    vs complete CDS 1,701 bp = 77 bp UTR difference) while preventing partial amplicons
    from retiring full-length segments.
    """
    if winner_real_len <= 0:
        return False
    if retired_real_len <= 0:
        return True
    if winner_real_len >= retired_real_len:
        return True
    shortfall = retired_real_len - winner_real_len
    ratio = winner_real_len / float(retired_real_len)
    return ratio >= min_ratio and shortfall <= max_shortfall_bp


def count_aligned_acgt_mismatches(
    aln_a: str,
    aln_b: str,
    max_mm: int = 2,
) -> Tuple[int, int]:
    """Count unambiguous A/C/G/T mismatches between two coordinate-aligned strings.

    Returns ``(mismatches, acgt_overlap)``, stopping early once ``mismatches > max_mm``.
    """
    a_up = aln_a.upper()
    b_up = aln_b.upper()
    limit = min(len(a_up), len(b_up))
    mm = 0
    overlap = 0
    acgt = _UNAMBIGUOUS_ACGT
    for i in range(limit):
        ca = a_up[i]
        if ca not in acgt:
            continue
        cb = b_up[i]
        if cb not in acgt:
            continue
        overlap += 1
        if ca != cb:
            mm += 1
            if mm > max_mm:
                return mm, overlap
    return mm, overlap


def aligned_sequences_agree(
    aln_winner: str,
    aln_retired: str,
    winner_real_len: int,
    retired_real_len: int,
    max_mm: int = 0,
    min_ratio: float = 0.94,
    max_shortfall_bp: int = 90,
) -> Tuple[bool, int, int]:
    """Verify coverage parity and aligned sequence agreement."""
    if not has_coverage_parity(winner_real_len, retired_real_len, min_ratio=min_ratio, max_shortfall_bp=max_shortfall_bp):
        return False, -1, 0
    mm, overlap = count_aligned_acgt_mismatches(aln_winner, aln_retired, max_mm=max_mm)
    min_len = min(winner_real_len, retired_real_len)
    required_overlap = min(min_len, max(1, int(0.80 * min_len)))
    if overlap < required_overlap or mm > max_mm:
        return False, mm, overlap
    return True, mm, overlap


def _strip_alignment_gaps(seq: Optional[str]) -> str:
    if not seq or not isinstance(seq, str):
        return ""
    return seq.upper().replace("-", "").replace(".", "").replace(" ", "")


def unaligned_sequences_agree(
    seq_winner: Optional[str],
    seq_retired: Optional[str],
    max_mm: int = 0,
    min_ratio: float = 0.94,
    max_shortfall_bp: int = 90,
) -> Tuple[bool, int, int]:
    """Compare two raw (pre-alignment) sequences with terminal-UTR overhang tolerance.

    Handles exact matches, substring containment (where one repository trimmed 5'/3'
    UTR primers), equal-length mismatch counting, and bounded terminal offset shifts
    when one sequence overhangs at 5' and the other overhangs at 3'.
    """
    uw = _strip_alignment_gaps(seq_winner)
    ur = _strip_alignment_gaps(seq_retired)
    len_w = count_acgt(uw)
    len_r = count_acgt(ur)
    if not has_coverage_parity(len_w, len_r, min_ratio=min_ratio, max_shortfall_bp=max_shortfall_bp):
        return False, -1, 0

    if uw == ur:
        return True, 0, len_w
    if uw in ur or ur in uw:
        return True, 0, min(len_w, len_r)

    min_raw = min(len(uw), len(ur))
    required_overlap = min(min(len_w, len_r), max(1, int(0.85 * min(len_w, len_r))))

    # Try zero-shift first, then bounded terminal shifts (up to +/-60 bp) for 5'/3' UTR trims.
    max_shift = min(60, abs(len(uw) - len(ur)) + 15)
    best_mm = max_mm + 1
    best_overlap = 0
    acgt = _UNAMBIGUOUS_ACGT

    shifts = [0] + [d for step in range(1, max_shift + 1) for d in (step, -step)]
    for shift in shifts:
        start_w = max(0, shift)
        start_r = max(0, -shift)
        span = min(len(uw) - start_w, len(ur) - start_r)
        if span < required_overlap:
            continue
        mm = 0
        overlap = 0
        for i in range(span):
            cw = uw[start_w + i]
            if cw not in acgt:
                continue
            cr = ur[start_r + i]
            if cr not in acgt:
                continue
            overlap += 1
            if cw != cr:
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm and overlap >= required_overlap:
            if mm == 0:
                return True, 0, overlap
            if mm < best_mm:
                best_mm = mm
                best_overlap = overlap

    if best_mm <= max_mm and best_overlap >= required_overlap:
        return True, best_mm, best_overlap
    return False, best_mm, best_overlap


def date_precision_rank(date_str: Optional[str]) -> Tuple[int, Dict[str, object]]:
    """Return ``(rank, split_components)`` where rank is 0..3 (none, Y, Y-M, Y-M-D)."""
    if not date_str or not str(date_str).strip():
        return 0, {"day": "", "month": "", "year": ""}
    parts = split_date_components(str(date_str).strip())
    has_y = parts.get("year") not in ("", None)
    has_m = parts.get("month") not in ("", None)
    has_d = parts.get("day") not in ("", None)
    if has_y and has_m and has_d:
        return 3, parts
    if has_y and has_m:
        return 2, parts
    if has_y:
        return 1, parts
    return 0, parts


def can_upgrade_collection_date(winner_date: Optional[str], retired_date: Optional[str]) -> Tuple[bool, Dict[str, object]]:
    """Return ``(True, retired_parts)`` when ``retired_date`` strictly refines ``winner_date``
    without contradicting any date component already present on ``winner_date``."""
    w_rank, w_parts = date_precision_rank(winner_date)
    r_rank, r_parts = date_precision_rank(retired_date)
    if r_rank <= w_rank:
        return False, w_parts
    # Never upgrade to a future year or pre-1900 year.
    r_year = int(r_parts["year"])
    if r_year < 1900 or r_year > 2100:
        return False, w_parts
    if w_rank >= 1 and int(w_parts["year"]) != r_year:
        return False, w_parts
    if w_rank >= 2 and int(w_parts["month"]) != int(r_parts["month"]):
        return False, w_parts
    return True, r_parts


COALESCE_TEXT_FIELDS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("Originating_Lab", ("Originating_Lab", "gisaid_originating_lab", "originating_lab")),
    ("Submitting_Lab", ("Submitting_Lab", "gisaid_submitting_lab", "submitting_lab")),
    ("Originating_Sample_Id", ("Originating_Sample_Id", "gisaid_originating_sample_id", "originating_sample_id")),
    ("Submitting_Sample_Id", ("Submitting_Sample_Id", "gisaid_submitting_sample_id", "submitting_sample_id")),
    ("Passage_History", ("Passage_History", "gisaid_passage_history", "passage_history")),
    ("Isolate_Id", ("Isolate_Id", "gisaid_isolate_id", "isolate_id")),
    ("Lineage", ("Lineage", "gisaid_lineage", "lineage")),
    ("Clade", ("Clade", "gisaid_clade", "clade")),
    ("host", ("host",)),
    ("country", ("country",)),
    ("serotype", ("serotype",)),
)


def coalesce_winner_metadata(
    winner: Dict[str, str],
    retired: Dict[str, str],
    retired_source: str = "external",
    update_split_date_cols: bool = True,
) -> Dict[str, str]:
    """Return a dict of ``{column: new_value}`` updates to enrich ``winner`` from ``retired``.

    Preserves ``ValidateDbTree.py``'s ``date_precision_lost_in_split_columns`` invariant
    by updating ``collection_year``, ``collection_mon``, and ``collection_day`` in
    lockstep whenever ``collection_date`` is upgraded to higher precision.
    """
    updates: Dict[str, str] = {}

    w_date = str(winner.get("collection_date") or "").strip()
    r_date = str(retired.get("collection_date") or "").strip()
    upgrade_ok, r_parts = can_upgrade_collection_date(w_date, r_date)
    if upgrade_ok:
        updates["collection_date"] = r_date
        if update_split_date_cols:
            for col, key in (("collection_year", "year"), ("collection_mon", "month"), ("collection_day", "day")):
                if col in winner:
                    val = r_parts.get(key, "")
                    updates[col] = str(val) if val != "" else ""

    for target_col, source_candidates in COALESCE_TEXT_FIELDS:
        cur_val = str(winner.get(target_col) or "").strip()
        if cur_val and cur_val.casefold() not in _PLACEHOLDER_STRAINS:
            continue
        for src_col in source_candidates:
            cand = str(retired.get(src_col) or "").strip()
            if cand and cand.casefold() not in _PLACEHOLDER_STRAINS:
                updates[target_col] = cand
                break

    retired_acc = str(retired.get("primary_accession") or "").strip()
    isolate_id = str(
        retired.get("Isolate_Id") or retired.get("gisaid_isolate_id") or ""
    ).strip()
    tokens: List[str] = []
    existing_alt = str(winner.get("alt_accessions") or "").strip()
    if existing_alt:
        tokens.extend(t.strip() for t in existing_alt.split(";") if t.strip())
    if retired_acc:
        tag = f"{retired_source}:{retired_acc}" if retired_source else retired_acc
        if tag not in tokens:
            tokens.append(tag)
    if isolate_id and isolate_id not in tokens:
        tokens.append(isolate_id)
    if tokens:
        new_alt = ";".join(tokens)
        if new_alt != existing_alt:
            updates["alt_accessions"] = new_alt

    return updates


INSDC_SOURCE_NAMES = frozenset({"", "genbank", "ncbi", "insdc", "ena", "ddbj"})


def is_preferred_source(source: Optional[str], prefer_source: str = "genbank") -> bool:
    """Return True if ``source`` belongs to the preferred repository class."""
    src = (source or "genbank").strip().lower()
    pref = (prefer_source or "genbank").strip().lower()
    if pref in ("genbank", "ncbi", "insdc", "base"):
        return src in INSDC_SOURCE_NAMES
    return src == pref


def candidate_winner_sort_key(
    cand: Dict[str, object],
    prefer_source: str = "genbank",
) -> Tuple[int, int, int, int, int, str]:
    """Deterministic sort key (smallest wins) for choosing 1 record per strain.

    Priority:
      1. Protected (``master`` / ``reference``) always beats query records.
      2. Preferred repository (GenBank / INSDC by default) beats external portals.
      3. Already-committed ``update_db`` record beats a new re-upload within the
         same source tier so existing database accessions stay stable across updates.
      4. Longest unambiguous ``real_length`` wins.
      5. Higher collection-date precision (``YYYY-MM-DD`` > ``YYYY-MM`` > ``YYYY``) wins.
      6. Lexicographically earliest accession breaks remaining ties deterministically.
    """
    protected_rank = 0 if bool(cand.get("protected")) else 1
    src_rank = 0 if is_preferred_source(str(cand.get("source") or ""), prefer_source) else 1
    origin_rank = 0 if str(cand.get("origin") or "") == "update_db" else 1
    real_len = int(cand.get("real_length") or 0)
    d_rank, _ = date_precision_rank(str(cand.get("collection_date") or ""))
    acc = str(cand.get("accession") or "")
    return (protected_rank, src_rank, origin_rank, -real_len, -d_rank, acc)


# ---------------------------------------------------------------------------
# Fast iterative O(N) Newick parser, leaf renamer & unary-collapsing pruner
# ---------------------------------------------------------------------------

_NEWICK_TOKEN_RE = re.compile(
    r"""
    (\s+)
    | (\[[^\]]*\])
    | ([(),:;])
    | ('(?:''|[^'])*'|"(?:""|[^"])*")
    | ([^(),:;\[\]\s]+)
    """,
    re.VERBOSE,
)


def _format_branch_len(val: float) -> str:
    if val == 0.0:
        return "0"
    text = f"{val:.10f}".rstrip("0").rstrip(".")
    return text if text and text != "-0" else "0"


def _add_branch_lengths(len_a: Optional[str], len_b: Optional[str]) -> Optional[str]:
    if len_a is None and len_b is None:
        return None
    va = float(len_a) if len_a is not None and len_a != "" else 0.0
    vb = float(len_b) if len_b is not None and len_b != "" else 0.0
    return _format_branch_len(va + vb)


def _unquote_label(tok: str) -> str:
    if len(tok) >= 2 and tok[0] == tok[-1] and tok[0] in ("'", '"'):
        q = tok[0]
        return tok[1:-1].replace(q + q, q)
    return tok


def prune_and_rename_newick(
    newick: str,
    prune_leaves: Set[str],
    rename_leaves: Optional[Dict[str, str]] = None,
) -> Tuple[str, Dict[str, int]]:
    """Single-pass iterative ``O(N)`` Newick leaf renamer + pruner.

    * Avoids recursion completely so 500k+-tip UShER ladder trees never hit
      Python's recursion limit.
    * Collapses degree-1 internal nodes created when one member of a 2-tip
      cherry (e.g. UShER ``_make_polytomy``) is pruned, adding the collapsed
      parent's branch length onto the surviving child.
    """
    rename_map = rename_leaves or {}
    if not prune_leaves and not rename_map:
        return newick, {"orig_tips": 0, "kept_tips": 0, "pruned_tips": 0, "renamed_tips": 0}

    parents: List[int] = [-1]
    children: List[List[int]] = [[]]
    labels: List[str] = [""]
    lengths: List[Optional[str]] = [None]
    comments: List[str] = [""]

    cur = 0
    expect_length = False

    for m in _NEWICK_TOKEN_RE.finditer(newick):
        ws, comment, punct, quoted, bare = m.groups()
        if ws:
            continue
        if comment:
            comments[cur] += comment
            continue
        if punct:
            expect_length = False
            if punct == "(":
                nid = len(parents)
                parents.append(cur)
                children.append([])
                labels.append("")
                lengths.append(None)
                comments.append("")
                children[cur].append(nid)
                cur = nid
            elif punct == ",":
                p = parents[cur]
                if p == -1:
                    raise ValueError("Unexpected ',' at Newick root")
                nid = len(parents)
                parents.append(p)
                children.append([])
                labels.append("")
                lengths.append(None)
                comments.append("")
                children[p].append(nid)
                cur = nid
            elif punct == ")":
                p = parents[cur]
                if p == -1:
                    raise ValueError("Unbalanced ')' in Newick")
                cur = p
            elif punct == ":":
                expect_length = True
            elif punct == ";":
                break
        else:
            tok = quoted if quoted is not None else bare
            assert tok is not None
            if expect_length:
                lengths[cur] = tok
                expect_length = False
            else:
                labels[cur] = tok

    num_nodes = len(parents)
    orig_is_leaf = [len(children[u]) == 0 for u in range(num_nodes)]
    orig_leaves: Set[str] = {
        _unquote_label(labels[u]).strip()
        for u in range(num_nodes)
        if orig_is_leaf[u] and labels[u]
    }
    orig_tip_count = len(orig_leaves)

    alive = [True] * num_nodes
    pruned_count = 0
    renamed_count = 0

    for u in range(num_nodes):
        if not orig_is_leaf[u]:
            continue
        clean = _unquote_label(labels[u]).strip()
        if clean in rename_map:
            target = rename_map[clean]
            if target not in orig_leaves:
                labels[u] = target
                orig_leaves.add(target)
                renamed_count += 1
            else:
                alive[u] = False
                pruned_count += 1
        elif clean in prune_leaves:
            alive[u] = False
            pruned_count += 1

    if pruned_count == 0 and renamed_count == 0:
        return newick, {
            "orig_tips": orig_tip_count,
            "kept_tips": orig_tip_count,
            "pruned_tips": 0,
            "renamed_tips": 0,
        }

    rep = list(range(num_nodes))
    for u in range(num_nodes - 1, -1, -1):
        if orig_is_leaf[u]:
            continue
        surv = [rep[c] for c in children[u] if alive[c]]
        if not surv:
            alive[u] = False
        elif len(surv) == 1:
            sole = surv[0]
            if parents[u] != -1:
                lengths[sole] = _add_branch_lengths(lengths[sole], lengths[u])
            else:
                lengths[sole] = None
            rep[u] = sole
            alive[u] = True
        else:
            children[u] = surv
            rep[u] = u
            alive[u] = True

    if not alive[0]:
        raise RuntimeError("Pruning removed every leaf from Newick tree")

    root = rep[0]
    lengths[root] = None

    out: List[str] = []
    stack: List[List[int]] = [[root, 0]]
    kept_tips = 0

    while stack:
        frame = stack[-1]
        u, idx = frame[0], frame[1]
        ch = children[u]
        if not ch:
            kept_tips += 1
            out.append(labels[u])
            if comments[u]:
                out.append(comments[u])
            if lengths[u] is not None:
                out.append(f":{lengths[u]}")
            stack.pop()
        else:
            if idx == 0:
                out.append("(")
            else:
                out.append(",")
            frame[1] = idx + 1
            stack.append([ch[idx], 0])
        while stack:
            top_u, top_idx = stack[-1]
            top_ch = children[top_u]
            if top_ch and top_idx == len(top_ch):
                stack.pop()
                out.append(")")
                if labels[top_u]:
                    out.append(labels[top_u])
                if comments[top_u]:
                    out.append(comments[top_u])
                if lengths[top_u] is not None:
                    out.append(f":{lengths[top_u]}")
            else:
                break

    out.append(";")
    return "".join(out), {
        "orig_tips": orig_tip_count,
        "kept_tips": kept_tips,
        "pruned_tips": pruned_count,
        "renamed_tips": renamed_count,
    }
