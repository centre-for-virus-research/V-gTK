import argparse
import re
import sqlite3
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from AnnotateMutations import (  # noqa: E402
    AnnotationMappingError,
    build_alignment_coordinate_map,
    canonicalize_product,
    extract_accession_tokens,
    load_db_gff_feature_maps,
    load_gene_alias_lookup,
    normalize_segment,
    prepare_catalog,
    prepare_sequence_alignment,
    resolve_reference_feature_map,
)


DEFAULT_MIN_YEARLY_SEQUENCES = 25
DEFAULT_PRIORITY_DRUG_COUNT = 2


TARGET_MUTATION_PLOT_DEFINITIONS = [
    {
        "mutation_id": "NS5A:93H",
        "plot_label": "NS5A Y93H",
        "slug": "ns5a_y93h",
        "lineages": [
            {"label": "Genotype 1", "genotype": "1"},
            {"label": "Genotype 3", "genotype": "3"},
        ],
    },
    {
        "mutation_id": "NS5A:30K",
        "plot_label": "NS5A A30K",
        "slug": "ns5a_a30k",
        "lineages": [
            {"label": "Genotype 3", "genotype": "3"},
        ],
    },
    {
        "mutation_id": "NS3:80K",
        "plot_label": "NS3/4A Q80K",
        "slug": "ns3_4a_q80k",
        "lineages": [
            {"label": "Genotype 1a", "genotype": "1", "subtype": "a"},
        ],
    },
    {
        "mutation_id": "NS5B:282T",
        "plot_label": "NS5B S282T",
        "slug": "ns5b_s282t",
        "lineages": "all_genotypes",
    },
]


SIGNATURE_PLOT_COLUMNS = [
    "collection_year",
    "drug",
    "signature_id",
    "signature_kind",
    "protein_name",
    "resistance_category",
    "combination_size",
    "callable_sequences",
    "sequences_with_signature",
    "proportion",
    "ci_lower",
    "ci_upper",
    "smoothed_proportion",
]


TARGET_MUTATION_PLOT_COLUMNS = [
    "collection_year",
    "mutation_id",
    "plot_label",
    "genotype_label",
    "callable_sequences",
    "sequences_with_mutation",
    "proportion",
    "ci_lower",
    "ci_upper",
    "smoothed_proportion",
]


def _load_yearly_cohort(conn):
    return pd.read_sql_query(
        """
        SELECT DISTINCT
            primary_accession,
            CAST(collection_year AS INTEGER) AS collection_year
        FROM meta_data
        WHERE TRIM(COALESCE(primary_accession, '')) != ''
          AND COALESCE(TRIM(exclusion_status), '') != '1'
          AND TRIM(COALESCE(collection_year, '')) GLOB '[0-9][0-9][0-9][0-9]'
        """,
        conn,
    )


def _load_relevant_mutation_counts(conn):
    return pd.read_sql_query(
        """
        SELECT
            primary_accession,
            relevant_mutations_present,
            total_relevant_mutation_count
        FROM sequence_relevant_mutation_summary
        """,
        conn,
    )


def _load_completed_signature_counts(conn):
    return pd.read_sql_query(
        """
        SELECT
            primary_accession,
            COUNT(DISTINCT signature_id) AS completed_signature_count
        FROM completed_signatures_only
        GROUP BY primary_accession
        """,
        conn,
    )


def _load_signature_drug_hits(conn):
    return pd.read_sql_query(
        """
        SELECT DISTINCT
            completed_signatures_only.primary_accession,
            completed_signatures_only.signature_id,
            completed_signatures_only.signature_kind,
            TRIM(COALESCE(mutation_catalog.protein_name, '')) AS protein_name,
            TRIM(COALESCE(mutation_catalog.drug, '')) AS drug
        FROM completed_signatures_only
        INNER JOIN mutation_catalog
            ON mutation_catalog.signature_id = completed_signatures_only.signature_id
        """,
        conn,
    )


def _load_sequence_alignment(conn):
    return prepare_sequence_alignment(pd.read_sql_query("SELECT * FROM sequence_alignment", conn))


def _load_signature_catalog(conn):
    return pd.read_sql_query(
        """
        SELECT DISTINCT
            signature_id,
            signature_kind,
            protein_name,
            segment,
            aa_position,
            reference_accession,
            TRIM(COALESCE(resistance_category, '')) AS resistance_category,
            TRIM(COALESCE(combination_size, '')) AS combination_size,
            TRIM(COALESCE(drug, '')) AS drug
        FROM mutation_catalog
        WHERE TRIM(COALESCE(signature_id, '')) != ''
        """,
        conn,
    )


def _load_target_mutation_catalog(conn, target_mutation_ids=None):
    target_mutation_ids = [str(value).strip() for value in (target_mutation_ids or []) if str(value).strip()]
    if not target_mutation_ids:
        return pd.DataFrame(columns=["mutation_id", "protein_name", "segment", "aa_position", "reference_accession"])

    placeholders = ", ".join(["?"] * len(target_mutation_ids))
    return pd.read_sql_query(
        f"""
        SELECT DISTINCT
            TRIM(COALESCE(mutation_id, '')) AS mutation_id,
            TRIM(COALESCE(protein_name, '')) AS protein_name,
            TRIM(COALESCE(segment, '')) AS segment,
            TRIM(COALESCE(aa_position, '')) AS aa_position,
            TRIM(COALESCE(reference_accession, '')) AS reference_accession
        FROM mutation_catalog
        WHERE TRIM(COALESCE(mutation_id, '')) IN ({placeholders})
        """,
        conn,
        params=tuple(target_mutation_ids),
    )


def _load_reference_metadata(conn):
    try:
        return pd.read_sql_query("SELECT * FROM meta_data", conn)
    except Exception:
        return pd.DataFrame()


def _choose_reference_row(group_df, reference_id, resolved_accession):
    candidates = []
    for value in [reference_id, resolved_accession]:
        for token in extract_accession_tokens(value):
            if token and token not in candidates:
                candidates.append(token)

    for candidate in candidates:
        mask = (
            group_df["sequence_id"].astype(str).str.strip() == candidate
        ) | (
            group_df["primary_accession"].astype(str).str.strip() == candidate
        )
        if mask.any():
            return group_df.loc[mask].iloc[0]
    return None


def _resolve_reference_contexts(seq_aln, signature_catalog, meta_data, alias_lookup, db_gff_maps):
    seq_aln = seq_aln.copy()
    seq_aln["_reference_id"] = seq_aln["alignment_name"].fillna("").astype(str).str.strip()

    if (seq_aln["_reference_id"] == "").any():
        reference_hints = [
            value
            for value in signature_catalog["reference_accession"].dropna().astype(str).str.strip().unique().tolist()
            if value
        ]
        if len(reference_hints) == 1:
            seq_aln.loc[seq_aln["_reference_id"] == "", "_reference_id"] = reference_hints[0]
        else:
            raise ValueError(
                "sequence_alignment contains blank alignment_name values and signature catalog does not resolve to one reference_accession"
            )

    contexts = {}
    unresolved = {}
    meta_data = meta_data if isinstance(meta_data, pd.DataFrame) else pd.DataFrame()
    master_candidates = []
    if not meta_data.empty and {"primary_accession", "accession_type"}.issubset(meta_data.columns):
        master_candidates = (
            meta_data.loc[
                meta_data["accession_type"].fillna("").astype(str).str.strip().str.lower() == "master",
                "primary_accession",
            ]
            .dropna()
            .astype(str)
            .str.strip()
            .unique()
            .tolist()
        )

    for reference_id, group_df in seq_aln.groupby("_reference_id", dropna=False):
        reference_id = str(reference_id or "").strip()
        if not reference_id:
            continue
        reference_candidates = master_candidates if master_candidates else [reference_id]
        try:
            resolved_accession, feature_map, attempted, source = resolve_reference_feature_map(
                reference_candidates[0],
                reference_candidates[1:],
                db_gff_maps,
                alias_lookup,
                allow_genbank_gff=False,
            )
            reference_row = _choose_reference_row(group_df, reference_id, resolved_accession)
            if reference_row is None:
                raise AnnotationMappingError(
                    f"No aligned reference sequence row was found for reference group {reference_id}. Candidates tried: {attempted}"
                )
            contexts[reference_id] = {
                "group_df": group_df.copy(),
                "feature_map": feature_map,
                "aligned_reference_map": build_alignment_coordinate_map(reference_row["alignment"]),
                "resolved_accession": resolved_accession,
                "source": source,
            }
        except AnnotationMappingError as exc:
            unresolved[reference_id] = str(exc)

    if not contexts:
        detail = "; ".join(f"{ref}: {msg}" for ref, msg in unresolved.items())
        raise ValueError(
            "Unable to resolve reference coordinate context for signature coverage calculation. " + detail
        )
    return contexts


def _build_signature_position_index(signature_catalog, alias_lookup, reference_contexts):
    prepared_catalog, _ = prepare_catalog(signature_catalog, alias_lookup)
    records = []

    for reference_id, context in reference_contexts.items():
        feature_map = context["feature_map"]
        aligned_reference_map = context["aligned_reference_map"]
        reference_catalog = prepared_catalog[
            prepared_catalog["reference_accession"].astype(str).str.strip() == str(reference_id).strip()
        ]
        if reference_catalog.empty:
            reference_catalog = prepared_catalog.copy()

        for _, row in reference_catalog.iterrows():
            protein_name = row["_canonical_protein"]
            aa_pos = row["_aa_position_int"]
            if not protein_name or aa_pos is None:
                continue

            gene_entry = feature_map.get(protein_name)
            if gene_entry is None:
                continue

            gene_start = int(gene_entry["cds_start"])
            query_positions = [gene_start + (aa_pos - 1) * 3 + offset for offset in range(3)]
            alignment_indices = []
            missing_coordinate = False
            for position in query_positions:
                alignment_index = aligned_reference_map.get(position)
                if alignment_index is None:
                    missing_coordinate = True
                    break
                alignment_indices.append(alignment_index)
            if missing_coordinate:
                continue

            records.append(
                {
                    "reference_id": reference_id,
                    "signature_id": row["signature_id"],
                    "signature_kind": row["signature_kind"],
                    "drug": str(row.get("drug", "") or "").strip(),
                    "protein_name": protein_name,
                    "resistance_category": str(row.get("resistance_category", "") or "").strip(),
                    "combination_size": str(row.get("combination_size", "") or "").strip(),
                    "alignment_indices": tuple(alignment_indices),
                }
            )

    if not records:
        return pd.DataFrame(
            columns=[
                "reference_id",
                "signature_id",
                "signature_kind",
                "drug",
                "protein_name",
                "resistance_category",
                "combination_size",
                "alignment_indices_list",
            ]
        )

    signature_positions = pd.DataFrame(records).drop_duplicates()
    return signature_positions.groupby(
        [
            "reference_id",
            "signature_id",
            "signature_kind",
            "drug",
            "protein_name",
            "resistance_category",
            "combination_size",
        ],
        as_index=False,
    ).agg(alignment_indices_list=("alignment_indices", lambda values: list(dict.fromkeys(values))))


def _sequence_has_callable_signature_coverage(alignment, alignment_indices_list):
    alignment = str(alignment or "")
    for alignment_indices in alignment_indices_list:
        if max(alignment_indices) >= len(alignment):
            return False
        codon = "".join(alignment[index] for index in alignment_indices)
        if len(codon) != 3:
            return False
        if any(base in "-Nn?Xx" for base in codon):
            return False
    return True


def _sequence_has_full_region_coverage(alignment, alignment_indices):
    alignment = str(alignment or "")
    if not alignment_indices:
        return False
    if max(alignment_indices) >= len(alignment):
        return False
    return not any(alignment[index] in "-Nn?Xx" for index in alignment_indices)


def _resistance_category_score(value):
    text = str(value or "").strip().upper()
    if not text:
        return 0

    numeric_match = re.search(r"\d+", text)
    if numeric_match:
        return int(numeric_match.group(0))

    roman_scores = {
        "VI": 6,
        "V": 5,
        "IV": 4,
        "III": 3,
        "II": 2,
        "I": 1,
    }
    for token, score in roman_scores.items():
        if re.search(rf"(^|[^A-Z]){token}([^A-Z]|$)", text):
            return score

    if "HIGH" in text:
        return 3
    if "INTERMEDIATE" in text or "MEDIUM" in text:
        return 2
    if "LOW" in text:
        return 1
    return 0


def _slugify_label(value):
    text = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "").strip().lower()).strip("_")
    return text or "unknown"


def _build_signature_metadata(signature_catalog, alias_lookup, reference_contexts):
    prepared_catalog, _ = prepare_catalog(signature_catalog, alias_lookup)
    records = []

    for reference_id in reference_contexts:
        reference_catalog = prepared_catalog[
            prepared_catalog["reference_accession"].astype(str).str.strip() == str(reference_id).strip()
        ]
        if reference_catalog.empty:
            reference_catalog = prepared_catalog.copy()

        for _, row in reference_catalog.iterrows():
            protein_name = row["_canonical_protein"]
            if not protein_name:
                continue
            records.append(
                {
                    "reference_id": reference_id,
                    "signature_id": row["signature_id"],
                    "signature_kind": row["signature_kind"],
                    "drug": str(row.get("drug", "") or "").strip(),
                    "protein_name": protein_name,
                    "resistance_category": str(row.get("resistance_category", "") or "").strip(),
                    "combination_size": str(row.get("combination_size", "") or "").strip(),
                }
            )

    if not records:
        return pd.DataFrame(
            columns=[
                "reference_id",
                "signature_id",
                "signature_kind",
                "drug",
                "protein_name",
                "resistance_category",
                "combination_size",
            ]
        )
    return pd.DataFrame(records).drop_duplicates().reset_index(drop=True)


def _build_gene_region_index(reference_contexts, protein_names=None):
    records = []
    allowed_proteins = set(protein_names.tolist() if isinstance(protein_names, np.ndarray) else (protein_names or []))

    for reference_id, context in reference_contexts.items():
        for protein_name, gene_entry in context["feature_map"].items():
            if allowed_proteins and protein_name not in allowed_proteins:
                continue

            alignment_indices = []
            missing_coordinate = False
            for position in range(int(gene_entry["cds_start"]), int(gene_entry["cds_end"]) + 1):
                alignment_index = context["aligned_reference_map"].get(position)
                if alignment_index is None:
                    missing_coordinate = True
                    break
                alignment_indices.append(alignment_index)

            if missing_coordinate:
                continue

            records.append(
                {
                    "reference_id": reference_id,
                    "protein_name": protein_name,
                    "alignment_indices": tuple(alignment_indices),
                }
            )

    if not records:
        return pd.DataFrame(columns=["reference_id", "protein_name", "alignment_indices"])
    return pd.DataFrame(records).drop_duplicates().reset_index(drop=True)


def _build_mutation_position_index(target_mutation_catalog, alias_lookup, reference_contexts):
    if target_mutation_catalog.empty:
        return pd.DataFrame(columns=["reference_id", "mutation_id", "alignment_indices_list"])

    prepared_catalog_input = target_mutation_catalog.copy()
    prepared_catalog_input["signature_id"] = prepared_catalog_input["mutation_id"]
    prepared_catalog_input["signature_kind"] = "single"
    prepared_catalog_input["drug"] = ""
    prepared_catalog_input["resistance_category"] = ""
    prepared_catalog_input["combination_size"] = ""
    prepared_catalog, _ = prepare_catalog(prepared_catalog_input, alias_lookup)
    records = []

    for reference_id, context in reference_contexts.items():
        feature_map = context["feature_map"]
        aligned_reference_map = context["aligned_reference_map"]
        reference_catalog = prepared_catalog[
            prepared_catalog["reference_accession"].astype(str).str.strip() == str(reference_id).strip()
        ]
        if reference_catalog.empty:
            reference_catalog = prepared_catalog.copy()

        for _, row in reference_catalog.iterrows():
            protein_name = row["_canonical_protein"]
            aa_pos = row["_aa_position_int"]
            if not protein_name or aa_pos is None:
                continue

            gene_entry = feature_map.get(protein_name)
            if gene_entry is None:
                continue

            gene_start = int(gene_entry["cds_start"])
            query_positions = [gene_start + (aa_pos - 1) * 3 + offset for offset in range(3)]
            alignment_indices = []
            missing_coordinate = False
            for position in query_positions:
                alignment_index = aligned_reference_map.get(position)
                if alignment_index is None:
                    missing_coordinate = True
                    break
                alignment_indices.append(alignment_index)
            if missing_coordinate:
                continue

            records.append(
                {
                    "reference_id": reference_id,
                    "mutation_id": str(row.get("mutation_id", "") or "").strip(),
                    "alignment_indices": tuple(alignment_indices),
                }
            )

    if not records:
        return pd.DataFrame(columns=["reference_id", "mutation_id", "alignment_indices_list"])

    mutation_positions = pd.DataFrame(records).drop_duplicates()
    return mutation_positions.groupby(["reference_id", "mutation_id"], as_index=False).agg(
        alignment_indices_list=("alignment_indices", lambda values: list(dict.fromkeys(values)))
    )


def _prepare_genotype_metadata(meta_data):
    if meta_data.empty or "primary_accession" not in meta_data.columns:
        return pd.DataFrame(columns=["primary_accession", "genotype", "subtype"])

    genotype_column = next(
        (column for column in ["nearest_reference_genotype", "genotype"] if column in meta_data.columns),
        None,
    )
    subtype_column = next(
        (column for column in ["nearest_reference_subtype", "subtype"] if column in meta_data.columns),
        None,
    )
    if genotype_column is None and subtype_column is None:
        return pd.DataFrame(columns=["primary_accession", "genotype", "subtype"])

    metadata = meta_data[["primary_accession"]].copy()
    metadata["genotype"] = (
        meta_data[genotype_column].fillna("").astype(str).str.strip() if genotype_column is not None else ""
    )
    metadata["subtype"] = (
        meta_data[subtype_column].fillna("").astype(str).str.strip().str.lower() if subtype_column is not None else ""
    )
    return metadata.drop_duplicates(subset=["primary_accession"]).reset_index(drop=True)


def _finalize_signature_frequency_grid(signature_grid):
    if signature_grid.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    signature_grid = signature_grid.copy()
    for column in ["resistance_category", "combination_size"]:
        if column not in signature_grid.columns:
            signature_grid[column] = ""
    signature_grid["callable_sequences"] = pd.to_numeric(
        signature_grid["callable_sequences"], errors="coerce"
    ).fillna(0)
    signature_grid["sequences_with_signature"] = pd.to_numeric(
        signature_grid["sequences_with_signature"], errors="coerce"
    ).fillna(0)
    signature_grid["proportion"] = np.where(
        signature_grid["callable_sequences"] > 0,
        signature_grid["sequences_with_signature"] / signature_grid["callable_sequences"],
        np.nan,
    )

    ci_lower = pd.Series(np.nan, index=signature_grid.index, dtype=float)
    ci_upper = pd.Series(np.nan, index=signature_grid.index, dtype=float)
    callable_mask = signature_grid["callable_sequences"] > 0
    if callable_mask.any():
        masked_lower, masked_upper = _wilson_95_ci(
            signature_grid.loc[callable_mask, "sequences_with_signature"],
            signature_grid.loc[callable_mask, "callable_sequences"],
        )
        ci_lower.loc[callable_mask] = masked_lower.to_numpy()
        ci_upper.loc[callable_mask] = masked_upper.to_numpy()
    signature_grid["ci_lower"] = ci_lower
    signature_grid["ci_upper"] = ci_upper
    signature_grid["smoothed_proportion"] = np.nan

    for (_, signature_id), group in signature_grid.groupby(["drug", "signature_id"], sort=False):
        callable_group = group[group["callable_sequences"] > 0].sort_values("collection_year")
        if callable_group.empty:
            continue
        smoothed = _smooth_weighted_proportions(
            callable_group["collection_year"].to_numpy(),
            callable_group["proportion"].to_numpy(),
            callable_group["callable_sequences"].to_numpy(),
        )
        signature_grid.loc[callable_group.index, "smoothed_proportion"] = smoothed

    return signature_grid[SIGNATURE_PLOT_COLUMNS].sort_values(
        ["drug", "signature_id", "collection_year"]
    ).reset_index(drop=True)


def _finalize_target_mutation_frequency_grid(target_grid):
    if target_grid.empty:
        return pd.DataFrame(columns=TARGET_MUTATION_PLOT_COLUMNS)

    target_grid = target_grid.copy()
    target_grid["callable_sequences"] = pd.to_numeric(target_grid["callable_sequences"], errors="coerce").fillna(0)
    target_grid["sequences_with_mutation"] = pd.to_numeric(
        target_grid["sequences_with_mutation"], errors="coerce"
    ).fillna(0)
    target_grid["proportion"] = np.where(
        target_grid["callable_sequences"] > 0,
        target_grid["sequences_with_mutation"] / target_grid["callable_sequences"],
        np.nan,
    )

    ci_lower = pd.Series(np.nan, index=target_grid.index, dtype=float)
    ci_upper = pd.Series(np.nan, index=target_grid.index, dtype=float)
    callable_mask = target_grid["callable_sequences"] > 0
    if callable_mask.any():
        masked_lower, masked_upper = _wilson_95_ci(
            target_grid.loc[callable_mask, "sequences_with_mutation"],
            target_grid.loc[callable_mask, "callable_sequences"],
        )
        ci_lower.loc[callable_mask] = masked_lower.to_numpy()
        ci_upper.loc[callable_mask] = masked_upper.to_numpy()
    target_grid["ci_lower"] = ci_lower
    target_grid["ci_upper"] = ci_upper
    target_grid["smoothed_proportion"] = np.nan

    for (_, genotype_label), group in target_grid.groupby(["mutation_id", "genotype_label"], sort=False):
        callable_group = group[group["callable_sequences"] > 0].sort_values("collection_year")
        if callable_group.empty:
            continue
        smoothed = _smooth_weighted_proportions(
            callable_group["collection_year"].to_numpy(),
            callable_group["proportion"].to_numpy(),
            callable_group["callable_sequences"].to_numpy(),
        )
        target_grid.loc[callable_group.index, "smoothed_proportion"] = smoothed

    return target_grid[TARGET_MUTATION_PLOT_COLUMNS].sort_values(
        ["mutation_id", "genotype_label", "collection_year"]
    ).reset_index(drop=True)


def _mutation_token_present(series, mutation_id):
    pattern = rf"(?:^|;)\s*{re.escape(str(mutation_id).strip())}\s*(?:;|$)"
    return series.fillna("").astype(str).str.contains(pattern, regex=True)


def _build_target_mutation_frequency(
    cohort,
    mutation_counts,
    target_mutation_catalog,
    seq_aln,
    meta_data,
    alias_lookup,
    db_gff_maps,
    target_mutation_definitions=None,
):
    target_mutation_definitions = target_mutation_definitions or TARGET_MUTATION_PLOT_DEFINITIONS
    if not target_mutation_definitions or target_mutation_catalog.empty or seq_aln.empty:
        return pd.DataFrame(columns=TARGET_MUTATION_PLOT_COLUMNS)

    genotype_metadata = _prepare_genotype_metadata(meta_data)
    if genotype_metadata.empty:
        return pd.DataFrame(columns=TARGET_MUTATION_PLOT_COLUMNS)

    cohort = cohort.merge(genotype_metadata, on="primary_accession", how="left")
    cohort["genotype"] = cohort["genotype"].fillna("").astype(str).str.strip()
    cohort["subtype"] = cohort["subtype"].fillna("").astype(str).str.strip().str.lower()

    reference_contexts = _resolve_reference_contexts(seq_aln, target_mutation_catalog, meta_data, alias_lookup, db_gff_maps)
    mutation_positions = _build_mutation_position_index(target_mutation_catalog, alias_lookup, reference_contexts)
    if mutation_positions.empty:
        return pd.DataFrame(columns=TARGET_MUTATION_PLOT_COLUMNS)

    callable_records = []
    for reference_id, context in reference_contexts.items():
        group_sequences = context["group_df"][["primary_accession", "alignment"]].drop_duplicates("primary_accession")
        group_sequences = cohort.merge(group_sequences, on="primary_accession", how="inner")
        if group_sequences.empty:
            continue

        group_mutation_positions = mutation_positions[
            mutation_positions["reference_id"].astype(str).str.strip() == str(reference_id).strip()
        ]
        for _, mutation_row in group_mutation_positions.iterrows():
            callable_mask = group_sequences["alignment"].map(
                lambda alignment: _sequence_has_callable_signature_coverage(
                    alignment,
                    mutation_row["alignment_indices_list"],
                )
            )
            callable_sequences = group_sequences.loc[
                callable_mask,
                ["primary_accession", "collection_year", "genotype", "subtype"],
            ].drop_duplicates()
            if callable_sequences.empty:
                continue
            callable_sequences["mutation_id"] = mutation_row["mutation_id"]
            callable_records.append(callable_sequences)

    if not callable_records:
        return pd.DataFrame(columns=TARGET_MUTATION_PLOT_COLUMNS)

    callable_detail = pd.concat(callable_records, ignore_index=True).drop_duplicates(
        subset=["primary_accession", "collection_year", "mutation_id", "genotype", "subtype"]
    )
    mutation_presence = mutation_counts[["primary_accession", "relevant_mutations_present"]].drop_duplicates(
        subset=["primary_accession"]
    )
    mutation_hits = []
    for target_definition in target_mutation_definitions:
        mutation_id = str(target_definition.get("mutation_id", "") or "").strip()
        if not mutation_id:
            continue
        matched = mutation_presence.loc[
            _mutation_token_present(mutation_presence["relevant_mutations_present"], mutation_id),
            ["primary_accession"],
        ].copy()
        if matched.empty:
            continue
        matched["mutation_id"] = mutation_id
        mutation_hits.append(matched)
    mutation_hit_detail = (
        pd.concat(mutation_hits, ignore_index=True).drop_duplicates()
        if mutation_hits
        else pd.DataFrame(columns=["primary_accession", "mutation_id"])
    )

    target_records = []
    years = cohort[["collection_year"]].drop_duplicates().sort_values("collection_year")
    for target_definition in target_mutation_definitions:
        mutation_id = str(target_definition.get("mutation_id", "") or "").strip()
        plot_label = str(target_definition.get("plot_label", mutation_id) or mutation_id).strip()
        mutation_callable = callable_detail[callable_detail["mutation_id"] == mutation_id].copy()
        if mutation_callable.empty:
            continue

        lineage_definitions = target_definition.get("lineages", [])
        if lineage_definitions == "all_genotypes":
            lineage_definitions = [
                {"label": f"Genotype {genotype}", "genotype": genotype}
                for genotype in sorted(
                    {str(value).strip() for value in mutation_callable["genotype"].dropna().tolist() if str(value).strip()},
                    key=lambda value: (0, int(value)) if str(value).isdigit() else (1, str(value)),
                )
            ]

        for lineage_definition in lineage_definitions:
            genotype = str(lineage_definition.get("genotype", "") or "").strip()
            subtype = str(lineage_definition.get("subtype", "") or "").strip().lower()
            genotype_label = str(lineage_definition.get("label", genotype or subtype) or genotype or subtype).strip()
            lineage_callable = mutation_callable.copy()
            if genotype:
                lineage_callable = lineage_callable[lineage_callable["genotype"] == genotype]
            if subtype:
                lineage_callable = lineage_callable[lineage_callable["subtype"] == subtype]

            callable_summary = (
                lineage_callable.groupby("collection_year", as_index=False)
                .agg(callable_sequences=("primary_accession", "nunique"))
                if not lineage_callable.empty
                else pd.DataFrame(columns=["collection_year", "callable_sequences"])
            )

            lineage_hits = lineage_callable.merge(
                mutation_hit_detail,
                on=["primary_accession", "mutation_id"],
                how="inner",
            )
            numerator_summary = (
                lineage_hits.groupby("collection_year", as_index=False)
                .agg(sequences_with_mutation=("primary_accession", "nunique"))
                if not lineage_hits.empty
                else pd.DataFrame(columns=["collection_year", "sequences_with_mutation"])
            )

            target_grid = years.copy()
            target_grid["mutation_id"] = mutation_id
            target_grid["plot_label"] = plot_label
            target_grid["genotype_label"] = genotype_label
            target_grid = target_grid.merge(callable_summary, on="collection_year", how="left").merge(
                numerator_summary,
                on="collection_year",
                how="left",
            )
            target_records.append(target_grid)

    if not target_records:
        return pd.DataFrame(columns=TARGET_MUTATION_PLOT_COLUMNS)
    return _finalize_target_mutation_frequency_grid(pd.concat(target_records, ignore_index=True))


def _build_yearly_signature_frequency(
    cohort,
    signature_hits,
    signature_catalog,
    seq_aln,
    meta_data,
    alias_lookup,
    db_gff_maps,
):
    if signature_catalog.empty or seq_aln.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    reference_contexts = _resolve_reference_contexts(seq_aln, signature_catalog, meta_data, alias_lookup, db_gff_maps)
    signature_positions = _build_signature_position_index(signature_catalog, alias_lookup, reference_contexts)
    if signature_positions.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    signature_positions = signature_positions.groupby(
        ["reference_id", "signature_id", "signature_kind", "drug"], as_index=False
    ).agg(
        protein_name=(
            "protein_name",
            lambda values: ";".join(sorted({str(value).strip() for value in values if str(value).strip()})),
        ),
        resistance_category=(
            "resistance_category",
            lambda values: next((str(value).strip() for value in values if str(value).strip()), ""),
        ),
        combination_size=(
            "combination_size",
            lambda values: next((str(value).strip() for value in values if str(value).strip()), ""),
        ),
        alignment_indices_list=(
            "alignment_indices_list",
            lambda groups: [indices for group in groups for indices in group],
        ),
    )

    callable_records = []
    for reference_id, context in reference_contexts.items():
        group_sequences = context["group_df"][["primary_accession", "alignment"]].drop_duplicates("primary_accession")
        group_sequences = cohort.merge(group_sequences, on="primary_accession", how="inner")
        if group_sequences.empty:
            continue

        group_signature_positions = signature_positions[
            signature_positions["reference_id"].astype(str).str.strip() == str(reference_id).strip()
        ]
        for _, signature_row in group_signature_positions.iterrows():
            callable_mask = group_sequences["alignment"].map(
                lambda alignment: _sequence_has_callable_signature_coverage(
                    alignment,
                    signature_row["alignment_indices_list"],
                )
            )
            callable_sequences = group_sequences.loc[
                callable_mask,
                ["primary_accession", "collection_year"],
            ].drop_duplicates()
            if callable_sequences.empty:
                continue
            callable_sequences["signature_id"] = signature_row["signature_id"]
            callable_sequences["signature_kind"] = signature_row["signature_kind"]
            callable_sequences["drug"] = signature_row["drug"]
            callable_sequences["protein_name"] = signature_row["protein_name"]
            callable_sequences["resistance_category"] = signature_row["resistance_category"]
            callable_sequences["combination_size"] = signature_row["combination_size"]
            callable_records.append(callable_sequences)

    signatures = signature_positions[
        [
            "signature_id",
            "signature_kind",
            "drug",
            "protein_name",
            "resistance_category",
            "combination_size",
        ]
    ].drop_duplicates().reset_index(drop=True)
    year_signature_grid = cohort[["collection_year"]].drop_duplicates().sort_values("collection_year")
    year_signature_grid = year_signature_grid.assign(key=1).merge(
        signatures.assign(key=1), on="key", how="inner"
    ).drop(columns="key")

    if not callable_records:
        year_signature_grid["callable_sequences"] = 0
        year_signature_grid["sequences_with_signature"] = 0
        return _finalize_signature_frequency_grid(year_signature_grid)

    callable_detail = pd.concat(callable_records, ignore_index=True).drop_duplicates(
        subset=["primary_accession", "collection_year", "signature_id"]
    )
    callable_summary = (
        callable_detail.groupby(
            [
                "collection_year",
                "signature_id",
                "signature_kind",
                "drug",
                "protein_name",
                "resistance_category",
                "combination_size",
            ],
            as_index=False,
        )
        .agg(callable_sequences=("primary_accession", "nunique"))
    )

    callable_hits = callable_detail.merge(
        signature_hits[["primary_accession", "signature_id"]].drop_duplicates(),
        on=["primary_accession", "signature_id"],
        how="inner",
    )
    numerator_summary = (
        callable_hits.groupby(
            [
                "collection_year",
                "signature_id",
                "signature_kind",
                "drug",
                "protein_name",
                "resistance_category",
                "combination_size",
            ],
            as_index=False,
        )
        .agg(sequences_with_signature=("primary_accession", "nunique"))
    )

    signature_grid = year_signature_grid.merge(
        callable_summary,
        on=[
            "collection_year",
            "signature_id",
            "signature_kind",
            "drug",
            "protein_name",
            "resistance_category",
            "combination_size",
        ],
        how="left",
    ).merge(
        numerator_summary,
        on=[
            "collection_year",
            "signature_id",
            "signature_kind",
            "drug",
            "protein_name",
            "resistance_category",
            "combination_size",
        ],
        how="left",
    )

    return _finalize_signature_frequency_grid(signature_grid)


def _build_yearly_gene_signature_frequency(
    cohort,
    signature_hits,
    signature_catalog,
    seq_aln,
    meta_data,
    alias_lookup,
    db_gff_maps,
):
    if signature_catalog.empty or seq_aln.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    reference_contexts = _resolve_reference_contexts(seq_aln, signature_catalog, meta_data, alias_lookup, db_gff_maps)
    signature_metadata = _build_signature_metadata(signature_catalog, alias_lookup, reference_contexts)
    if signature_metadata.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    gene_regions = _build_gene_region_index(reference_contexts, signature_metadata["protein_name"].unique())
    if gene_regions.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    signature_metadata = signature_metadata.merge(
        gene_regions,
        on=["reference_id", "protein_name"],
        how="inner",
    )
    if signature_metadata.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    callable_records = []
    for reference_id, context in reference_contexts.items():
        group_sequences = context["group_df"][["primary_accession", "alignment"]].drop_duplicates("primary_accession")
        group_sequences = cohort.merge(group_sequences, on="primary_accession", how="inner")
        if group_sequences.empty:
            continue

        reference_signatures = signature_metadata[
            signature_metadata["reference_id"].astype(str).str.strip() == str(reference_id).strip()
        ]
        for protein_name, gene_df in reference_signatures.groupby("protein_name", sort=False):
            alignment_indices = gene_df["alignment_indices"].iloc[0]
            callable_mask = group_sequences["alignment"].map(
                lambda alignment: _sequence_has_full_region_coverage(alignment, alignment_indices)
            )
            full_gene_sequences = group_sequences.loc[
                callable_mask,
                ["primary_accession", "collection_year"],
            ].drop_duplicates()
            if full_gene_sequences.empty:
                continue

            for _, signature_row in gene_df.drop(columns="alignment_indices").drop_duplicates().iterrows():
                gene_callable = full_gene_sequences.copy()
                gene_callable["signature_id"] = signature_row["signature_id"]
                gene_callable["signature_kind"] = signature_row["signature_kind"]
                gene_callable["drug"] = signature_row["drug"]
                gene_callable["protein_name"] = protein_name
                gene_callable["resistance_category"] = signature_row["resistance_category"]
                gene_callable["combination_size"] = signature_row["combination_size"]
                callable_records.append(gene_callable)

    signatures = signature_metadata[
        [
            "signature_id",
            "signature_kind",
            "drug",
            "protein_name",
            "resistance_category",
            "combination_size",
        ]
    ].drop_duplicates().reset_index(drop=True)
    year_signature_grid = cohort[["collection_year"]].drop_duplicates().sort_values("collection_year")
    year_signature_grid = year_signature_grid.assign(key=1).merge(
        signatures.assign(key=1), on="key", how="inner"
    ).drop(columns="key")

    if not callable_records:
        year_signature_grid["callable_sequences"] = 0
        year_signature_grid["sequences_with_signature"] = 0
        return _finalize_signature_frequency_grid(year_signature_grid)

    callable_detail = pd.concat(callable_records, ignore_index=True).drop_duplicates(
        subset=["primary_accession", "collection_year", "signature_id", "protein_name"]
    )
    callable_summary = (
        callable_detail.groupby(
            [
                "collection_year",
                "signature_id",
                "signature_kind",
                "drug",
                "protein_name",
                "resistance_category",
                "combination_size",
            ],
            as_index=False,
        )
        .agg(callable_sequences=("primary_accession", "nunique"))
    )

    callable_hits = callable_detail.merge(
        signature_hits[["primary_accession", "signature_id"]].drop_duplicates(),
        on=["primary_accession", "signature_id"],
        how="inner",
    )
    numerator_summary = (
        callable_hits.groupby(
            [
                "collection_year",
                "signature_id",
                "signature_kind",
                "drug",
                "protein_name",
                "resistance_category",
                "combination_size",
            ],
            as_index=False,
        )
        .agg(sequences_with_signature=("primary_accession", "nunique"))
    )

    signature_grid = year_signature_grid.merge(
        callable_summary,
        on=[
            "collection_year",
            "signature_id",
            "signature_kind",
            "drug",
            "protein_name",
            "resistance_category",
            "combination_size",
        ],
        how="left",
    ).merge(
        numerator_summary,
        on=[
            "collection_year",
            "signature_id",
            "signature_kind",
            "drug",
            "protein_name",
            "resistance_category",
            "combination_size",
        ],
        how="left",
    )

    return _finalize_signature_frequency_grid(signature_grid)


def _select_priority_drugs(gene_signature_frequency, max_drugs=DEFAULT_PRIORITY_DRUG_COUNT):
    if gene_signature_frequency.empty:
        return pd.DataFrame(columns=["drug", "resistance_score", "combination_signature_count", "max_proportion", "max_hits"])

    ranking_df = gene_signature_frequency.copy()
    ranking_df = ranking_df[
        ranking_df["drug"].fillna("").astype(str).str.strip() != ""
    ]
    ranking_df = ranking_df[ranking_df["callable_sequences"] > 0]
    if ranking_df.empty:
        return pd.DataFrame(columns=["drug", "resistance_score", "combination_signature_count", "max_proportion", "max_hits"])

    ranking_df["resistance_score"] = ranking_df["resistance_category"].map(_resistance_category_score)
    ranking_df["is_combination"] = ranking_df["signature_kind"].fillna("").astype(str).str.strip().str.lower() == "combination"

    ranked = (
        ranking_df.groupby("drug", as_index=False)
        .agg(
            resistance_score=("resistance_score", "max"),
            combination_signature_count=("is_combination", "sum"),
            max_proportion=("proportion", lambda series: float(np.nanmax(series.to_numpy())) if len(series) else 0.0),
            max_hits=("sequences_with_signature", "max"),
        )
        .sort_values(
            ["resistance_score", "combination_signature_count", "max_proportion", "max_hits", "drug"],
            ascending=[False, False, False, False, True],
        )
        .head(int(max_drugs))
        .reset_index(drop=True)
    )
    return ranked


def _select_priority_signature_subset(gene_signature_frequency, drug):
    if gene_signature_frequency.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    drug_df = gene_signature_frequency[
        gene_signature_frequency["drug"].fillna("").astype(str).str.strip() == str(drug).strip()
    ].copy()
    if drug_df.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    candidate_df = drug_df[drug_df["callable_sequences"] > 0].copy()
    if candidate_df.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    candidate_df["resistance_score"] = candidate_df["resistance_category"].map(_resistance_category_score)
    selected_signature_ids = []
    for _, gene_df in candidate_df.groupby("protein_name", sort=False):
        highest_score = int(gene_df["resistance_score"].max()) if not gene_df.empty else 0
        priority = gene_df[
            (gene_df["signature_kind"].fillna("").astype(str).str.strip().str.lower() == "combination")
            | (gene_df["resistance_score"] == highest_score)
        ]
        if priority.empty:
            priority = gene_df
        selected_signature_ids.extend(priority["signature_id"].tolist())

    selected_signature_ids = sorted(set(selected_signature_ids))
    return drug_df[drug_df["signature_id"].isin(selected_signature_ids)].copy()


def _select_category_i_signature_subset(gene_signature_frequency, drug):
    if gene_signature_frequency.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    drug_df = gene_signature_frequency[
        gene_signature_frequency["drug"].fillna("").astype(str).str.strip() == str(drug).strip()
    ].copy()
    if drug_df.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    candidate_df = drug_df[drug_df["callable_sequences"] > 0].copy()
    if candidate_df.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    candidate_df["resistance_score"] = candidate_df["resistance_category"].map(_resistance_category_score)
    category_i_df = candidate_df[candidate_df["resistance_score"] == 1].copy()
    if category_i_df.empty:
        return pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)

    return drug_df[drug_df["signature_id"].isin(sorted(set(category_i_df["signature_id"].tolist())))].copy()


def _smooth_weighted_proportions(years, proportions, sample_sizes, span_fraction=0.22):
    years = np.asarray(years, dtype=float)
    proportions = np.asarray(proportions, dtype=float)
    sample_sizes = np.asarray(sample_sizes, dtype=float)

    if len(years) == 0:
        return np.array([])
    if len(years) == 1:
        return proportions.copy()

    year_range = max(years.max() - years.min(), 1.0)
    bandwidth = max(year_range * span_fraction, 2.0)
    smoothed = np.empty_like(proportions, dtype=float)

    for index, target_year in enumerate(years):
        distances = np.abs(years - target_year) / bandwidth
        kernel_weights = np.where(distances < 1, (1 - distances**3) ** 3, 0.0)
        combined_weights = kernel_weights * sample_sizes
        if not np.any(combined_weights > 0):
            smoothed[index] = proportions[index]
            continue
        smoothed[index] = np.average(proportions, weights=combined_weights)

    return np.clip(smoothed, 0.0, 1.0)


def _wilson_95_ci(successes, totals):
    successes = pd.Series(successes, dtype=float)
    totals = pd.Series(totals, dtype=float).clip(lower=1)
    z = 1.96
    proportions = successes / totals
    denominator = 1 + (z**2 / totals)
    centre = (proportions + (z**2 / (2 * totals))) / denominator
    margin = (
        (z / denominator)
        * np.sqrt((proportions * (1 - proportions) / totals) + (z**2 / (4 * totals**2)))
    )
    lower = (centre - margin).clip(lower=0.0)
    upper = (centre + margin).clip(upper=1.0)
    return lower, upper


def build_yearly_summaries(db_path, min_yearly_sequences=DEFAULT_MIN_YEARLY_SEQUENCES):
    conn = sqlite3.connect(str(db_path))
    try:
        cohort = _load_yearly_cohort(conn)
        mutation_counts = _load_relevant_mutation_counts(conn)
        signature_counts = _load_completed_signature_counts(conn)
        signature_hits = _load_signature_drug_hits(conn)
        signature_catalog = _load_signature_catalog(conn)
        target_mutation_catalog = _load_target_mutation_catalog(
            conn,
            [definition["mutation_id"] for definition in TARGET_MUTATION_PLOT_DEFINITIONS],
        )
        meta_data = _load_reference_metadata(conn)
        alias_lookup = load_gene_alias_lookup(conn)
        db_gff_maps = load_db_gff_feature_maps(conn, alias_lookup)
        try:
            seq_aln = _load_sequence_alignment(conn)
        except Exception:
            seq_aln = pd.DataFrame()
    finally:
        conn.close()

    if cohort.empty:
        raise ValueError("No non-excluded sequences with a valid four-digit collection_year were found in meta_data")

    cohort = cohort.copy()
    cohort["collection_year"] = cohort["collection_year"].astype(int)

    mutation_counts = mutation_counts.rename(columns={"total_relevant_mutation_count": "relevant_mutation_count"})
    mutation_counts["relevant_mutation_count"] = pd.to_numeric(
        mutation_counts["relevant_mutation_count"], errors="coerce"
    ).fillna(0)

    signature_counts["completed_signature_count"] = pd.to_numeric(
        signature_counts["completed_signature_count"], errors="coerce"
    ).fillna(0)

    mutation_cohort = cohort.merge(mutation_counts, on="primary_accession", how="left")
    mutation_cohort["relevant_mutation_count"] = mutation_cohort["relevant_mutation_count"].fillna(0)

    signature_cohort = cohort.merge(signature_counts, on="primary_accession", how="left")
    signature_cohort["completed_signature_count"] = signature_cohort["completed_signature_count"].fillna(0)

    yearly_denominator = (
        cohort.groupby("collection_year", as_index=False)
        .agg(total_sequences=("primary_accession", "nunique"))
        .sort_values("collection_year")
    )
    eligible_years = yearly_denominator.loc[
        yearly_denominator["total_sequences"] >= int(min_yearly_sequences),
        "collection_year",
    ]

    if eligible_years.empty:
        raise ValueError(
            f"No collection_year values met the minimum yearly sample size threshold of {int(min_yearly_sequences)}"
        )

    mutation_cohort = mutation_cohort[mutation_cohort["collection_year"].isin(eligible_years)]
    signature_cohort = signature_cohort[signature_cohort["collection_year"].isin(eligible_years)]
    cohort = cohort[cohort["collection_year"].isin(eligible_years)]

    yearly_mutation_summary = (
        mutation_cohort.groupby("collection_year", as_index=False)
        .agg(
            sequence_count=("primary_accession", "nunique"),
            mean_relevant_mutation_count=("relevant_mutation_count", "mean"),
        )
        .sort_values("collection_year")
    )

    yearly_signature_summary = (
        signature_cohort.groupby("collection_year", as_index=False)
        .agg(
            sequence_count=("primary_accession", "nunique"),
            mean_completed_signature_count=("completed_signature_count", "mean"),
        )
        .sort_values("collection_year")
    )
    yearly_denominator = yearly_denominator[yearly_denominator["collection_year"].isin(eligible_years)]

    yearly_signature_frequency = _build_yearly_signature_frequency(
        cohort,
        signature_hits,
        signature_catalog,
        seq_aln,
        meta_data,
        alias_lookup,
        db_gff_maps,
    )
    yearly_gene_signature_frequency = _build_yearly_gene_signature_frequency(
        cohort,
        signature_hits,
        signature_catalog,
        seq_aln,
        meta_data,
        alias_lookup,
        db_gff_maps,
    )
    targeted_mutation_frequency = _build_target_mutation_frequency(
        cohort,
        mutation_counts,
        target_mutation_catalog,
        seq_aln,
        meta_data,
        alias_lookup,
        db_gff_maps,
    )
    priority_drugs = _select_priority_drugs(yearly_gene_signature_frequency)

    yearly_sample_size = yearly_denominator.rename(columns={"total_sequences": "sequence_count"})

    return {
        "yearly_mutation_summary": yearly_mutation_summary,
        "yearly_signature_summary": yearly_signature_summary,
        "yearly_signature_frequency": yearly_signature_frequency,
        "yearly_gene_signature_frequency": yearly_gene_signature_frequency,
        "targeted_mutation_frequency": targeted_mutation_frequency,
        "priority_drugs": priority_drugs,
        "yearly_sample_size": yearly_sample_size,
    }


def _plot_mean_series(ax, df, value_column, title, ylabel, color):
    ax.plot(df["collection_year"], df[value_column], color=color, linewidth=2.4, marker="o", markersize=4.2)
    ax.set_title(title, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.tick_params(axis="both", labelsize=15)


def _apply_year_axis(ax, yearly_sample_size):
    if yearly_sample_size.empty:
        return
    min_year = int(yearly_sample_size["collection_year"].min())
    max_year = int(yearly_sample_size["collection_year"].max())
    ax.set_xlim(min_year - 0.5, max_year + 0.5)


def _plot_signature_proportions(ax, df, show_ci_cloud=False):
    if df.empty:
        ax.text(0.5, 0.5, "No callable signature coverage could be resolved", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    cmap = getattr(cm, "tab20")
    plot_df = df[df["callable_sequences"] > 0].copy()
    if plot_df.empty:
        ax.text(0.5, 0.5, "No callable signature coverage could be resolved", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    max_callable = float(plot_df["callable_sequences"].max())
    for index, (_, signature_df) in enumerate(plot_df.groupby(["drug", "signature_id"], sort=False)):
        signature_df = signature_df.sort_values("collection_year")
        color = cmap(index % 20)
        label_bits = [value for value in [signature_df["drug"].iloc[0], signature_df["signature_id"].iloc[0]] if value]
        label = " | ".join(label_bits) if label_bits else signature_df["signature_id"].iloc[0]
        marker_sizes = 18 + 40 * np.sqrt(signature_df["callable_sequences"] / max_callable)
        if show_ci_cloud:
            ax.fill_between(
                signature_df["collection_year"],
                signature_df["ci_lower"],
                signature_df["ci_upper"],
                color=color,
                alpha=0.12,
                linewidth=0,
                zorder=1,
            )
        ax.scatter(
            signature_df["collection_year"],
            signature_df["proportion"],
            color=color,
            s=marker_sizes,
            alpha=0.28,
            linewidths=0,
            zorder=3,
        )
        ax.plot(
            signature_df["collection_year"],
            signature_df["smoothed_proportion"],
            label=label,
            color=color,
            linewidth=2.7,
            zorder=4,
        )
    ax.set_title("Coverage-aware completed signature frequency by year", fontsize=24)
    ax.set_ylabel("Completed signature frequency", fontsize=18)
    ax.set_xlabel("Collection year", fontsize=18)
    ax.set_ylim(0, 1)
    ax.set_xlim(int(plot_df["collection_year"].min()) - 0.5, int(plot_df["collection_year"].max()) + 0.5)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.tick_params(axis="both", labelsize=16)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.18, 1.02),
        frameon=False,
        title="Drug | Signature",
        fontsize=16,
        title_fontsize=18,
        ncol=1,
    )


def _plot_subset_signature_proportions(ax, df, title, legend_title, label_columns, show_ci_cloud=False):
    if df.empty:
        ax.text(0.5, 0.5, "No signatures met the plotting criteria", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    cmap = getattr(cm, "tab20")
    plot_df = df[df["callable_sequences"] > 0].copy()
    if plot_df.empty:
        ax.text(0.5, 0.5, "No signatures met the plotting criteria", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    max_callable = float(plot_df["callable_sequences"].max())
    for index, (_, signature_df) in enumerate(plot_df.groupby(label_columns, sort=False)):
        signature_df = signature_df.sort_values("collection_year")
        color = cmap(index % 20)
        label = " | ".join(
            [str(signature_df[column].iloc[0]).strip() for column in label_columns if str(signature_df[column].iloc[0]).strip()]
        )
        marker_sizes = 18 + 40 * np.sqrt(signature_df["callable_sequences"] / max_callable)
        if show_ci_cloud:
            ax.fill_between(
                signature_df["collection_year"],
                signature_df["ci_lower"],
                signature_df["ci_upper"],
                color=color,
                alpha=0.12,
                linewidth=0,
                zorder=1,
            )
        ax.scatter(
            signature_df["collection_year"],
            signature_df["proportion"],
            color=color,
            s=marker_sizes,
            alpha=0.28,
            linewidths=0,
            zorder=3,
        )
        ax.plot(
            signature_df["collection_year"],
            signature_df["smoothed_proportion"],
            label=label,
            color=color,
            linewidth=2.7,
            zorder=4,
        )
    ax.set_title(title, fontsize=24)
    ax.set_ylabel("Completed signature frequency", fontsize=18)
    ax.set_xlabel("Collection year", fontsize=18)
    ax.set_ylim(0, 1)
    ax.set_xlim(int(plot_df["collection_year"].min()) - 0.5, int(plot_df["collection_year"].max()) + 0.5)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.tick_params(axis="both", labelsize=16)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.18, 1.02),
        frameon=False,
        title=legend_title,
        fontsize=16,
        title_fontsize=18,
        ncol=1,
    )


def _plot_target_mutation_proportions(ax, df, title, show_ci_cloud=False):
    if df.empty:
        ax.text(0.5, 0.5, "No genotype-specific mutation frequencies could be resolved", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    cmap = getattr(cm, "tab20")
    plot_df = df[df["callable_sequences"] > 0].copy()
    if plot_df.empty:
        ax.text(0.5, 0.5, "No genotype-specific mutation frequencies could be resolved", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    max_callable = float(plot_df["callable_sequences"].max())
    for index, (_, genotype_df) in enumerate(plot_df.groupby("genotype_label", sort=False)):
        genotype_df = genotype_df.sort_values("collection_year")
        color = cmap(index % 20)
        marker_sizes = 18 + 40 * np.sqrt(genotype_df["callable_sequences"] / max_callable)
        if show_ci_cloud:
            ax.fill_between(
                genotype_df["collection_year"],
                genotype_df["ci_lower"],
                genotype_df["ci_upper"],
                color=color,
                alpha=0.12,
                linewidth=0,
                zorder=1,
            )
        ax.scatter(
            genotype_df["collection_year"],
            genotype_df["proportion"],
            color=color,
            s=marker_sizes,
            alpha=0.28,
            linewidths=0,
            zorder=3,
        )
        ax.plot(
            genotype_df["collection_year"],
            genotype_df["smoothed_proportion"],
            label=genotype_df["genotype_label"].iloc[0],
            color=color,
            linewidth=2.7,
            zorder=4,
        )

    ax.set_title(title, fontsize=24)
    ax.set_ylabel("Mutation frequency", fontsize=18)
    ax.set_xlabel("Collection year", fontsize=18)
    ax.set_ylim(0, 1)
    ax.set_xlim(int(plot_df["collection_year"].min()) - 0.5, int(plot_df["collection_year"].max()) + 0.5)
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.tick_params(axis="both", labelsize=16)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.18, 1.02),
        frameon=False,
        title="Genotype",
        fontsize=16,
        title_fontsize=18,
        ncol=1,
    )


def _add_callable_sequence_overlay(ax, df, ylabel):
    plot_df = df[df["callable_sequences"] > 0].copy()
    if plot_df.empty:
        return None

    yearly_callable = (
        plot_df.groupby("collection_year", as_index=False)
        .agg(callable_sequences=("callable_sequences", "max"))
        .sort_values("collection_year")
    )
    if yearly_callable.empty:
        return None

    count_ax = ax.twinx()
    (count_line,) = count_ax.plot(
        yearly_callable["collection_year"],
        yearly_callable["callable_sequences"],
        color="#222222",
        linestyle="--",
        linewidth=2.0,
        marker="o",
        markersize=4.5,
        alpha=0.85,
        label=ylabel,
        zorder=2,
    )
    count_ax.set_ylabel(ylabel, fontsize=16, color="#222222")
    count_ax.tick_params(axis="y", labelsize=14, colors="#222222")
    count_ax.set_ylim(0, max(1.0, float(yearly_callable["callable_sequences"].max()) * 1.1))
    return count_line


def generate_plots(db_path, output_dir, min_yearly_sequences=DEFAULT_MIN_YEARLY_SEQUENCES):
    db_path = Path(db_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries = build_yearly_summaries(db_path, min_yearly_sequences=min_yearly_sequences)
    mutation_summary = summaries["yearly_mutation_summary"]
    signature_summary = summaries["yearly_signature_summary"]
    signature_frequency = summaries["yearly_signature_frequency"]
    gene_signature_frequency = summaries["yearly_gene_signature_frequency"]
    targeted_mutation_frequency = summaries["targeted_mutation_frequency"]
    priority_drugs = summaries["priority_drugs"]
    yearly_sample_size = summaries["yearly_sample_size"]

    figure, axes = plt.subplots(3, 1, figsize=(16, 20), constrained_layout=False)

    _plot_mean_series(
        axes[0],
        mutation_summary,
        "mean_relevant_mutation_count",
        "Average relevant mutation count per year",
        "Mean relevant mutation count",
        "#1f4e79",
    )
    _apply_year_axis(axes[0], yearly_sample_size)

    _plot_mean_series(
        axes[1],
        signature_summary,
        "mean_completed_signature_count",
        "Average completed signature count per year",
        "Mean completed signature count",
        "#8c2d04",
    )
    axes[1].set_xlabel("Collection year", fontsize=16)
    _apply_year_axis(axes[1], yearly_sample_size)

    _plot_signature_proportions(axes[2], signature_frequency)
    figure.subplots_adjust(right=0.76, hspace=0.34)

    combined_plot_path = output_dir / "hcv_mutation_trends.png"
    figure.savefig(str(combined_plot_path), dpi=260, bbox_inches="tight")
    plt.close(figure)

    mutation_plot_path = output_dir / "mean_relevant_mutation_count_by_year.png"
    signature_plot_path = output_dir / "mean_completed_signature_count_by_year.png"
    signature_frequency_plot_path = output_dir / "signature_frequency_by_year.png"
    signature_frequency_plot_with_ci_path = output_dir / "signature_frequency_by_year_with_ci.png"
    gene_plot_dir = output_dir / "plots_by_gene"
    drug_plot_dir = output_dir / "plots_by_drug"
    target_mutation_plot_dir = output_dir / "plots_by_target_mutation"
    gene_plot_dir.mkdir(parents=True, exist_ok=True)
    drug_plot_dir.mkdir(parents=True, exist_ok=True)
    target_mutation_plot_dir.mkdir(parents=True, exist_ok=True)

    for data_frame, value_column, title, ylabel, color, plot_path in [
        (
            mutation_summary,
            "mean_relevant_mutation_count",
            "Average relevant mutation count per year",
            "Mean relevant mutation count",
            "#1f4e79",
            mutation_plot_path,
        ),
        (
            signature_summary,
            "mean_completed_signature_count",
            "Average completed signature count per year",
            "Mean completed signature count",
            "#8c2d04",
            signature_plot_path,
        ),
    ]:
        figure, ax = plt.subplots(figsize=(11, 5.5))
        _plot_mean_series(ax, data_frame, value_column, title, ylabel, color)
        ax.set_xlabel("Collection year", fontsize=16)
        _apply_year_axis(ax, yearly_sample_size)
        figure.savefig(str(plot_path), dpi=240, bbox_inches="tight")
        plt.close(figure)

    figure, ax = plt.subplots(figsize=(18, 11))
    _plot_signature_proportions(ax, signature_frequency)
    figure.subplots_adjust(right=0.69)
    figure.savefig(str(signature_frequency_plot_path), dpi=260, bbox_inches="tight")
    plt.close(figure)

    figure, ax = plt.subplots(figsize=(18, 11))
    _plot_signature_proportions(ax, signature_frequency, show_ci_cloud=True)
    figure.subplots_adjust(right=0.69)
    figure.savefig(str(signature_frequency_plot_with_ci_path), dpi=260, bbox_inches="tight")
    plt.close(figure)

    gene_plot_paths = {}
    for gene in sorted(
        gene_signature_frequency["protein_name"].dropna().astype(str).str.strip().replace("", np.nan).dropna().unique().tolist()
    ):
        gene_df = gene_signature_frequency[gene_signature_frequency["protein_name"] == gene].copy()
        figure, ax = plt.subplots(figsize=(18, 11))
        _plot_subset_signature_proportions(
            ax,
            gene_df,
            f"{gene} signatures with full-gene coverage by year",
            "Drug | Signature",
            ["drug", "signature_id"],
        )
        callable_line = _add_callable_sequence_overlay(ax, gene_df, "Sequences with full-gene coverage")
        if callable_line is not None:
            handles, labels = ax.get_legend_handles_labels()
            handles.append(callable_line)
            labels.append(callable_line.get_label())
            ax.legend(
                handles,
                labels,
                loc="upper left",
                bbox_to_anchor=(1.18, 1.02),
                frameon=False,
                title="Drug | Signature",
                fontsize=16,
                title_fontsize=18,
                ncol=1,
            )
        figure.subplots_adjust(right=0.69)
        plot_path = gene_plot_dir / f"{_slugify_label(gene)}_signature_frequency_by_year.png"
        figure.savefig(str(plot_path), dpi=260, bbox_inches="tight")
        plt.close(figure)
        gene_plot_paths[gene] = plot_path

    target_mutation_plot_paths = {}
    for target_definition in TARGET_MUTATION_PLOT_DEFINITIONS:
        mutation_id = str(target_definition.get("mutation_id", "") or "").strip()
        plot_label = str(target_definition.get("plot_label", mutation_id) or mutation_id).strip()
        mutation_df = targeted_mutation_frequency[
            targeted_mutation_frequency["mutation_id"].fillna("").astype(str).str.strip() == mutation_id
        ].copy()
        if mutation_df.empty:
            continue
        figure, ax = plt.subplots(figsize=(18, 11))
        _plot_target_mutation_proportions(ax, mutation_df, f"{plot_label} frequency by year")
        figure.subplots_adjust(right=0.69)
        plot_path = target_mutation_plot_dir / f"{_slugify_label(target_definition.get('slug', mutation_id))}_frequency_by_year.png"
        figure.savefig(str(plot_path), dpi=260, bbox_inches="tight")
        plt.close(figure)
        target_mutation_plot_paths[mutation_id] = plot_path

    priority_drug_plot_paths = {}
    priority_drug_category_i_plot_paths = {}
    for drug in priority_drugs["drug"].tolist():
        drug_df = _select_priority_signature_subset(gene_signature_frequency, drug)
        if drug_df.empty:
            drug_df = pd.DataFrame(columns=SIGNATURE_PLOT_COLUMNS)
        else:
            figure, ax = plt.subplots(figsize=(18, 11))
            _plot_subset_signature_proportions(
                ax,
                drug_df,
                f"{drug} priority signatures with full-gene coverage by year",
                "Gene | Signature",
                ["protein_name", "signature_id"],
            )
            figure.subplots_adjust(right=0.69)
            plot_path = drug_plot_dir / f"{_slugify_label(drug)}_priority_signature_frequency_by_year.png"
            figure.savefig(str(plot_path), dpi=260, bbox_inches="tight")
            plt.close(figure)
            priority_drug_plot_paths[drug] = plot_path

        category_i_df = _select_category_i_signature_subset(gene_signature_frequency, drug)
        if category_i_df.empty:
            continue
        figure, ax = plt.subplots(figsize=(18, 11))
        _plot_subset_signature_proportions(
            ax,
            category_i_df,
            f"{drug} Category I signatures with full-gene coverage by year",
            "Gene | Signature",
            ["protein_name", "signature_id"],
        )
        figure.subplots_adjust(right=0.69)
        plot_path = drug_plot_dir / f"{_slugify_label(drug)}_category_i_signature_frequency_by_year.png"
        figure.savefig(str(plot_path), dpi=260, bbox_inches="tight")
        plt.close(figure)
        priority_drug_category_i_plot_paths[drug] = plot_path

    mutation_summary.to_csv(output_dir / "yearly_relevant_mutation_summary.tsv", sep="\t", index=False)
    signature_summary.to_csv(output_dir / "yearly_completed_signature_summary.tsv", sep="\t", index=False)
    signature_frequency.to_csv(output_dir / "yearly_signature_frequency.tsv", sep="\t", index=False)
    gene_signature_frequency.to_csv(output_dir / "yearly_gene_signature_frequency.tsv", sep="\t", index=False)
    targeted_mutation_frequency.to_csv(output_dir / "targeted_mutation_frequency.tsv", sep="\t", index=False)
    priority_drugs.to_csv(output_dir / "priority_antiviral_selection.tsv", sep="\t", index=False)
    yearly_sample_size.to_csv(output_dir / "yearly_sample_size.tsv", sep="\t", index=False)

    summaries["plot_paths"] = {
        "combined": combined_plot_path,
        "mutations": mutation_plot_path,
        "completed_signatures": signature_plot_path,
        "signature_frequency": signature_frequency_plot_path,
        "signature_frequency_with_ci": signature_frequency_plot_with_ci_path,
        "gene_signature_frequency": gene_plot_paths,
        "targeted_mutation_frequency": target_mutation_plot_paths,
        "priority_drug_signature_frequency": priority_drug_plot_paths,
        "priority_drug_category_i_signature_frequency": priority_drug_category_i_plot_paths,
    }
    return summaries


def main():
    parser = argparse.ArgumentParser(
        description="Plot HCV mutation and completed-signature trends by collection year from a SQLite DB"
    )
    parser.add_argument("--db", required=True, help="Path to the SQLite DB")
    parser.add_argument("--output_dir", required=True, help="Directory to write plot files into")
    parser.add_argument(
        "--min_yearly_sequences",
        type=int,
        default=DEFAULT_MIN_YEARLY_SEQUENCES,
        help="Only plot years with at least this many non-excluded sequences",
    )
    args = parser.parse_args()

    summaries = generate_plots(args.db, args.output_dir, min_yearly_sequences=args.min_yearly_sequences)
    print(f"Saved plots to {Path(args.output_dir).resolve()}")
    print(f"Years plotted: {int(summaries['yearly_mutation_summary']['collection_year'].min())} - {int(summaries['yearly_mutation_summary']['collection_year'].max())}")
    print(
        "Sample size range: "
        f"{int(summaries['yearly_sample_size']['sequence_count'].min())} - "
        f"{int(summaries['yearly_sample_size']['sequence_count'].max())} sequences per year"
    )


if __name__ == "__main__":
    # Tiny test DB example:
    # python scripts/PlotHcvMutationTrends.py --db /tmp/hcv_mutation_plot_test.db --output_dir /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/plots_mutations --min_yearly_sequences 2
    # Full DB example:
    # python scripts/PlotHcvMutationTrends.py --db /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/HCV_full.db --output_dir /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/plots_mutations --min_yearly_sequences 25
    main()