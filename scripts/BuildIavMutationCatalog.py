#!/usr/bin/env python3
"""Convert curated raw influenza mutation CSVs (`dev/flu_raw_tables/Clean *.csv`)
into the HCV-compatible evidence-linked TSV + publications CSV consumed by
`scripts/AnnotateMutations.py`.

Key invariants enforced:
  - Uses `csv.reader` (never pandas default NA parsing) so `Segment == 'NA'`
    (Neuraminidase!) is NEVER dropped as NULL.
  - Maps `Clean HA.csv` H3 mature numbering (11..550) onto the pipeline's Seg 4
    master `AB573800` aligned AA coordinates (1..550) via the verified H3
    (`CY121117` A/Aichi/2/1968) <-> `AB573800` pan-HA coordinate alignment,
    preserving both `h3_numbering` and `h5_numbering` labels alongside the
    mapped master coordinate.
  - Explodes semicolon-separated papers and DOIs (including `No DOI - PMID - ...`)
    into one evidence row per reference in `mutation_catalog`, with a 1-row-per-
    `evidence_id` `publications.csv` matching HCV's exact equality-join contract.
  - Preserves HCV wild-type suppression semantics (`wild_type_residues`) for
    markers fixed in seasonal human H3N2/H1N1 (e.g. `PB2:627K` in H3N2,
    `HA:228S` in H3N2) so avian/zoonotic markers fire cleanly on H5N1/H7N9/H9N2
    (including `PB2:M631L` in dairy-cattle H5N1) without flooding seasonal flu.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

PROTEIN_TO_SEGMENT: Dict[str, str] = {
    "PB2": "1",
    "PB1": "2",
    "PB1-F2": "2",
    "PA": "3",
    "PA-X": "3",
    "HA": "4",
    "NP": "5",
    "NA": "6",
    "M1": "7",
    "M2": "7",
    "NS1": "8",
    "NEP": "8",
}

SEGMENT_MASTER_ACCESSION: Dict[str, str] = {
    "1": "KJ889313",
    "2": "CY087822",
    "3": "MW333949",
    "4": "AB573800",
    "5": "NC_007369",
    "6": "AB472016",
    "7": "NC_002016",
    "8": "KY243296",
}

# Verified H3 mature HA1+HA2 continuous numbering (A/Aichi/2/1968 CY121117,
# QDLPG = pos 1, GLFGAIAGF HA2 = pos 330) -> Seg 4 master AB573800 ungapped AA
# position (1..550, DTLCV = pos 1). Derived from refset_4_aln pan-HA guide MSA
# where every tested H1/H3/H5/H7/H9 reference aligns codon-for-codon.
# Breakpoints: (h3_start, h3_end, offset_subtracted_from_h3)
# Verified landmarks:
#   H3 11 -> 1 (offset -10)
#   H3 17 -> 7 (offset -10)
#   H3 158 -> 155 (offset -3)
#   H3 160 -> 157 (offset -3)
#   H3 180 -> 177 (offset -3)
#   H3 190 -> 187 (offset -3)
#   H3 225 -> 222 (offset -3)
#   H3 226 -> 223 (offset -3)
#   H3 228 -> 225 (offset -3)
#   H3 329 -> 327 (offset -2)
#   H3 330 -> 328 (offset -2)
# Exact per-position map built dynamically when a Seg 4 alignment (DB or FASTA)
# is supplied, with an exact fallback derived from CY121117 vs AB573800.
_DEFAULT_H3_TO_AB573800_EXACT: Optional[Dict[int, Optional[int]]] = None


def _compute_h3_to_master_map_from_alignments(master_aln: str, h3_aln: str) -> Dict[int, Optional[int]]:
    codon_map = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }

    def _tc(aln: str, col: int) -> str:
        c = aln[col : col + 3].upper().replace("U", "T")
        if len(c) < 3 or "-" in c:
            return "-"
        return codon_map.get(c, "X")

    h3_ungap = [_tc(h3_aln, c) for c in range(0, len(master_aln) - 2, 3) if _tc(h3_aln, c) != "-"]
    qdlpg_idx = "".join(h3_ungap).find("QDLPG")
    if qdlpg_idx < 0:
        qdlpg_idx = 0

    mapping: Dict[int, Optional[int]] = {}
    m_aa = 0
    h3_idx = -1
    for col in range(0, len(master_aln) - 2, 3):
        mc = _tc(master_aln, col)
        c3 = _tc(h3_aln, col)
        if mc != "-":
            m_aa += 1
        if c3 != "-":
            h3_idx += 1
            h3_pos = h3_idx - qdlpg_idx + 1
            mapping[h3_pos] = m_aa if mc != "-" else None
    return mapping


def _load_h3_to_master_map(ref_aln_fasta: Optional[str] = None, db_path: Optional[str] = None) -> Dict[int, Optional[int]]:
    """Return exact `{h3_mature_pos: ab573800_master_aa_pos}` map."""
    global _DEFAULT_H3_TO_AB573800_EXACT

    # 1. Explicit FASTA if provided
    if ref_aln_fasta and os.path.isfile(ref_aln_fasta):
        seqs: Dict[str, str] = {}
        cur_id = None
        chunks: List[str] = []
        with open(ref_aln_fasta, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if cur_id is not None:
                        seqs[cur_id] = "".join(chunks)
                    cur_id = line[1:].split()[0].split(".")[0]
                    chunks = []
                else:
                    chunks.append(line)
            if cur_id is not None:
                seqs[cur_id] = "".join(chunks)
        if "AB573800" in seqs and "CY121117" in seqs:
            return _compute_h3_to_master_map_from_alignments(seqs["AB573800"], seqs["CY121117"])

    # 2. Live or supplied SQLite DB if available
    candidate_dbs = [db_path] if db_path else [
        str(Path(__file__).resolve().parents[1] / "test_out" / "IAV_DB_genbank_n_gisaid_rebuild" / "db_build" / "iav-db.db"),
    ]
    for cdb in candidate_dbs:
        if cdb and os.path.isfile(cdb):
            try:
                import sqlite3
                conn = sqlite3.connect(f"file:{cdb}?mode=ro", uri=True)
                try:
                    m4 = conn.execute("SELECT alignment FROM sequence_alignment WHERE primary_accession='AB573800' LIMIT 1").fetchone()
                    h3 = conn.execute("SELECT alignment FROM sequence_alignment WHERE primary_accession='CY121117' LIMIT 1").fetchone()
                    if m4 and h3 and m4[0] and h3[0]:
                        return _compute_h3_to_master_map_from_alignments(m4[0], h3[0])
                finally:
                    conn.close()
            except Exception:
                pass

    # 3. Deterministic piecewise fallback (exact for AB573800 vs CY121117 pan-HA alignment)
    if _DEFAULT_H3_TO_AB573800_EXACT is None:
        m: Dict[int, Optional[int]] = {}
        for p in range(1, 562):
            if p < 11:
                m[p] = None
            elif p <= 47:
                m[p] = p - 10
            elif p <= 130:
                m[p] = p - 9
            elif p <= 145:
                m[p] = p - 6
            elif p <= 275:
                m[p] = p - 3
            else:
                m[p] = p - 2
        _DEFAULT_H3_TO_AB573800_EXACT = m
    return dict(_DEFAULT_H3_TO_AB573800_EXACT)


# Subtype-specific wild-type residues for classic mammalian-adaptation sites that
# are naturally fixed in seasonal human H3N2 or H1N1 (using HCV's exact
# wild_type_residues mechanism so avian/zoonotic H5N1/H7N9/H9N2 fire cleanly
# without flooding 1.17M normal seasonal human H3N2 sequences).
KNOWN_SUBTYPE_WILD_TYPES: Dict[Tuple[str, int], str] = {
    ("PB2", 627): "H3N2:K:99.0;H5N1:E:99.0;H1N1:E:95.0;H7N9:E:95.0;H9N2:E:99.0",
    ("PB2", 271): "H3N2:A:99.0;H1N1:A:99.0;H5N1:T:99.0;H7N9:T:99.0;H9N2:T:99.0",
    ("PB2", 591): "H1N1:R:95.0;H3N2:Q:99.0;H5N1:Q:99.0;H7N9:Q:99.0;H9N2:Q:99.0",
    # HA master pos 223 = H3 226 (Q226L); master pos 225 = H3 228 (G228S)
    ("HA", 223): "H3N2:I:90.0;H3N2:L:90.0;H5N1:Q:99.0;H1N1:Q:99.0;H7N9:Q:50.0;H9N2:Q:50.0",
    ("HA", 225): "H3N2:S:99.0;H5N1:G:99.0;H1N1:G:99.0;H7N9:G:99.0;H9N2:G:99.0",
}

_MUT_TOKEN_RE = re.compile(r"^([A-Z\*])(-?\d+)([A-Z\*∆Δ-])$")
_PMID_RE = re.compile(r"PMID\s*[-:]?\s*(\d+)", re.IGNORECASE)
_YEAR_RE = re.compile(r"\b(19\d{2}|20\d{2})\b")

CATALOG_COLUMNS: List[str] = [
    "protein_name",
    "segment",
    "aa_position",
    "alt_residue",
    "reference_accession",
    "mutation_id",
    "mutation_type",
    "signature_id",
    "signature_kind",
    "combination_id",
    "combination_size",
    "phenotype",
    "genotype",
    "resistance_category",
    "drug",
    "drug_category",
    "drug_producer",
    "any_in_vitro_evidence",
    "in_vitro_max_ec50_midpoint",
    "any_in_vivo_evidence",
    "in_vivo_baseline",
    "in_vivo_treatment_emergent",
    "evidence_id",
    "data_source",
    "evidence_url",
    "evidence_label",
    "evidence_type",
    "regimens",
    "linked_evidence_ids",
    "finding_ids",
    "relevant_genotypes",
    "wild_type_residues",
    "wt_residue",
    "serotypes_tested",
    "h3_numbering",
    "h5_numbering",
    "experimentally_verified",
    "discovery_method",
    "notes",
]


def _clean_val(val: Optional[str], *, preserve_na: bool = False) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    if not preserve_na and s.upper() == "NA":
        return ""
    return s


def _normalise_alt(alt_raw: str) -> str:
    a = alt_raw.strip()
    if a in ("∆", "Δ", "del", "DEL", "-"):
        return "-"
    return a.upper()


def _parse_evidence_pairs(papers_raw: str, dois_raw: str) -> List[Tuple[str, str, str, str]]:
    """Return list of ``(evidence_id, data_source, evidence_url, evidence_label)``."""
    papers = [p.strip() for p in papers_raw.split(";") if p.strip() and p.strip().upper() != "NA"]
    dois = [d.strip() for d in dois_raw.split(";") if d.strip() and d.strip().upper() != "NA"]

    def _parse_single_doi(tok: str) -> Tuple[str, str, str]:
        m_pmid = _PMID_RE.search(tok)
        if m_pmid:
            pmid = m_pmid.group(1)
            return f"pmid:{pmid}", "pubmed", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        cleaned = tok.strip().rstrip(".")
        if cleaned.lower().startswith("https://doi.org/"):
            doi_core = cleaned[len("https://doi.org/") :]
            return f"doi:{doi_core}", "publication", f"https://doi.org/{doi_core}"
        if cleaned.lower().startswith("http://doi.org/"):
            doi_core = cleaned[len("http://doi.org/") :]
            return f"doi:{doi_core}", "publication", f"https://doi.org/{doi_core}"
        if cleaned.startswith("10."):
            return f"doi:{cleaned}", "publication", f"https://doi.org/{cleaned}"
        if cleaned.lower().startswith("http"):
            slug = re.sub(r"[^a-zA-Z0-9._-]+", "_", cleaned)[:80]
            return f"url:{slug}", "publication", cleaned
        if cleaned.lower() in ("missing", "none", "unknown", ""):
            return "", "", ""
        slug = re.sub(r"[^a-zA-Z0-9._-]+", "_", cleaned)[:80]
        return f"ref:{slug}", "publication", ""

    results: List[Tuple[str, str, str, str]] = []
    n = max(len(papers), len(dois))
    if n == 0:
        return [("", "", "", "")]
    for i in range(n):
        label = papers[i] if i < len(papers) else (papers[-1] if papers else "")
        doi_tok = dois[i] if i < len(dois) else ""
        eid, dsrc, url = _parse_single_doi(doi_tok) if doi_tok else ("", "", "")
        if not eid and label:
            slug = re.sub(r"[^a-zA-Z0-9]+", "_", label.lower()).strip("_")[:64]
            eid = f"lit:{slug}"
            dsrc = "publication"
        results.append((eid, dsrc, url, label))
    return results


def _infer_phenotype(notes: str, protein: str) -> str:
    low = notes.lower()
    if "ph stability" in low or "ph " in low:
        return "ha_ph_stability"
    if "alpha 2,6" in low or "2,6 binding" in low or "receptor" in low:
        return "receptor_binding_alpha2_6"
    if "btn3a3" in low:
        return "btn3a3_evasion"
    if "ferret" in low or "transmission" in low or "guinea pig" in low:
        return "mammalian_transmission"
    if "mouse" in low or "mice" in low or "horse" in low or "dog" in low or "mammal" in low:
        return "mammalian_adaptation"
    return "host_adaptation"


def _parse_citation_fields(label: str) -> Tuple[str, str, str]:
    """Extract rough `(authors_short, journal, year)` from a short citation label like `Herfst et al. Science 2012`."""
    yr_m = _YEAR_RE.search(label)
    year = yr_m.group(1) if yr_m else ""
    authors = label
    journal = ""
    if "et al." in label:
        head, tail = label.split("et al.", 1)
        authors = f"{head.strip()} et al."
        journal = _YEAR_RE.sub("", tail).strip(" .,;")
    return authors, journal, year


def build_iav_catalog(
    raw_dir: str,
    out_catalog_tsv: str,
    out_publications_csv: str,
    ref_aln_fasta: Optional[str] = None,
    db_path: Optional[str] = None,
) -> Dict[str, int]:
    """Convert all `Clean *.csv` files in `raw_dir` into HCV-compatible TSV + publications CSV."""
    h3_to_master = _load_h3_to_master_map(ref_aln_fasta=ref_aln_fasta, db_path=db_path)
    csv_files = sorted(glob.glob(os.path.join(raw_dir, "Clean *.csv")))
    if not csv_files:
        raise FileNotFoundError(f"No 'Clean *.csv' files found in {raw_dir}")

    out_rows: List[Dict[str, str]] = []
    publications: "OrderedDict[str, Dict[str, str]]" = OrderedDict()
    skipped_empty = 0
    skipped_unmapped_ha = 0
    seen_row_keys = set()

    for path in csv_files:
        fn = os.path.basename(path)
        is_na_file = fn.lower() == "clean na.csv"
        is_ha_file = fn.lower() == "clean ha.csv"

        with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as fh:
            reader = csv.reader(fh)
            hdr = next(reader, None)
            if not hdr:
                continue
            for raw_row in reader:
                row = {hdr[i]: (raw_row[i].strip() if i < len(raw_row) else "") for i in range(len(hdr))}
                protein = _clean_val(row.get("Segment"), preserve_na=is_na_file)
                if not protein or (not is_na_file and protein.upper() == "NA"):
                    skipped_empty += 1
                    continue
                if protein not in PROTEIN_TO_SEGMENT:
                    skipped_empty += 1
                    continue

                seg = PROTEIN_TO_SEGMENT[protein]
                ref_acc = SEGMENT_MASTER_ACCESSION[seg]

                wt = _clean_val(row.get("WT"))
                pos_raw = _clean_val(row.get("Position"))
                if is_ha_file:
                    alt_raw = _clean_val(row.get("Mutation"))
                    h3_label = _clean_val(row.get("Mutation.H3.numbering"))
                    h5_label = _clean_val(row.get("Mutation.H5.numbering"))
                else:
                    alt_raw = _clean_val(row.get("Mutant"))
                    h3_label = ""
                    h5_label = ""

                if not alt_raw:
                    skipped_empty += 1
                    continue
                alt_norm = _normalise_alt(alt_raw)
                if not alt_norm or len(alt_norm) != 1 or alt_norm not in "ACDEFGHIKLMNPQRSTVWY*-":
                    skipped_empty += 1
                    continue

                # Determine exact coordinate in master space
                if is_ha_file:
                    m_h3 = _MUT_TOKEN_RE.match(h3_label)
                    try:
                        h3_pos = int(m_h3.group(2)) if m_h3 else int(pos_raw)
                    except ValueError:
                        skipped_unmapped_ha += 1
                        continue
                    master_pos = h3_to_master.get(h3_pos)
                    if master_pos is None or master_pos < 1:
                        skipped_unmapped_ha += 1
                        continue
                    aa_pos = int(master_pos)
                    mut_id = f"HA:{aa_pos}{alt_norm}"
                else:
                    try:
                        aa_pos = int(pos_raw)
                    except ValueError:
                        skipped_empty += 1
                        continue
                    if aa_pos < 1:
                        skipped_empty += 1
                        continue
                    mut_id = f"{protein}:{aa_pos}{alt_norm}"

                mut_type = "aminoAcidDeletion" if alt_norm == "-" else "aminoAcidSimplePolymorphism"
                notes_parts = [
                    _clean_val(row.get("Notes")),
                    _clean_val(row.get("Extra.Notes")),
                    _clean_val(row.get("Further.Notes")),
                ]
                notes = "; ".join(p.strip(" ;") for p in notes_parts if p.strip(" ;"))
                phenotype = _infer_phenotype(notes, protein)
                serotypes_tested = _clean_val(row.get("Which.virus.")).strip(" ;")
                exp_ver = _clean_val(row.get("Experimentally.Verified"))
                disc = _clean_val(row.get("How.was.the.mutation.discovered"))

                wt_residues_col = KNOWN_SUBTYPE_WILD_TYPES.get((protein, aa_pos), "")
                ev_pairs = _parse_evidence_pairs(
                    _clean_val(row.get("Paper..s.")),
                    _clean_val(row.get("DOI")),
                )

                for eid, dsrc, url, label in ev_pairs:
                    dedup_key = (mut_id, eid, h3_label, alt_norm)
                    if dedup_key in seen_row_keys:
                        continue
                    seen_row_keys.add(dedup_key)

                    if eid and eid not in publications:
                        authors, journal, year = _parse_citation_fields(label)
                        publications[eid] = {
                            "evidence_id": eid,
                            "id": eid,
                            "data_source": dsrc or "publication",
                            "title": label,
                            "authors_short": authors,
                            "authors": authors,
                            "year": year,
                            "journal": journal,
                            "url": url,
                        }

                    out_rows.append(
                        {
                            "protein_name": protein,
                            "segment": seg,
                            "aa_position": str(aa_pos),
                            "alt_residue": alt_norm,
                            "reference_accession": ref_acc,
                            "mutation_id": mut_id,
                            "mutation_type": mut_type,
                            "signature_id": mut_id,
                            "signature_kind": "single",
                            "combination_id": "",
                            "combination_size": "",
                            "phenotype": phenotype,
                            "genotype": "all",
                            "resistance_category": "",
                            "drug": "",
                            "drug_category": "",
                            "drug_producer": "",
                            "any_in_vitro_evidence": "1" if exp_ver.lower() == "yes" else "",
                            "in_vitro_max_ec50_midpoint": "",
                            "any_in_vivo_evidence": "1" if any(k in notes.lower() for k in ("mice", "mouse", "ferret", "horse", "dog", "guinea")) else "",
                            "in_vivo_baseline": "",
                            "in_vivo_treatment_emergent": "",
                            "evidence_id": eid,
                            "data_source": dsrc,
                            "evidence_url": url,
                            "evidence_label": label,
                            "evidence_type": "literature" if eid else "",
                            "regimens": "",
                            "linked_evidence_ids": "",
                            "finding_ids": "",
                            "relevant_genotypes": "",
                            "wild_type_residues": wt_residues_col,
                            "wt_residue": wt,
                            "serotypes_tested": serotypes_tested,
                            "h3_numbering": h3_label,
                            "h5_numbering": h5_label,
                            "experimentally_verified": exp_ver,
                            "discovery_method": disc,
                            "notes": notes,
                        }
                    )

    os.makedirs(os.path.dirname(os.path.abspath(out_catalog_tsv)), exist_ok=True)
    with open(out_catalog_tsv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CATALOG_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    os.makedirs(os.path.dirname(os.path.abspath(out_publications_csv)), exist_ok=True)
    pub_cols = ["evidence_id", "id", "data_source", "title", "authors_short", "authors", "year", "journal", "url"]
    with open(out_publications_csv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=pub_cols)
        writer.writeheader()
        writer.writerows(publications.values())

    stats = {
        "catalog_rows": len(out_rows),
        "unique_mutations": len({r["mutation_id"] for r in out_rows}),
        "publications": len(publications),
        "skipped_empty": skipped_empty,
        "skipped_unmapped_ha": skipped_unmapped_ha,
    }
    print(f"[BuildIavMutationCatalog] Wrote {out_catalog_tsv} & {out_publications_csv}: {stats}")
    return stats


def main() -> int:
    parser = argparse.ArgumentParser(description="Build HCV-compatible IAV mutation catalog from raw flu CSVs.")
    parser.add_argument("--raw_dir", default="dev/flu_raw_tables", help="Directory containing 'Clean *.csv' files")
    parser.add_argument("--out_catalog", default="generic/influenza/Tables/generalized_mutation_catalog_evidence_linked.tsv")
    parser.add_argument("--out_publications", default="generic/influenza/Tables/publications.csv")
    parser.add_argument("--ref_aln_fasta", default=None, help="Optional Seg 4 alignment FASTA with AB573800 and CY121117")
    parser.add_argument("--db", default=None, help="Optional IAV SQLite DB for exact AB573800/CY121117 alignment map")
    args = parser.parse_args()
    build_iav_catalog(
        raw_dir=args.raw_dir,
        out_catalog_tsv=args.out_catalog,
        out_publications_csv=args.out_publications,
        ref_aln_fasta=args.ref_aln_fasta,
        db_path=args.db,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
