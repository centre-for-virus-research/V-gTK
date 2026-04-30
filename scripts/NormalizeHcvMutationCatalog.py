#!/usr/bin/env python3
"""Normalize HCV mutation assets into a single general catalog handoff file."""

import csv
import re
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


PHDR_PREFIX = "phdr_ras:"
MUTATION_TOKEN_RE = re.compile(r"^(?P<position>\d+)(?P<alt_residue>[A-Z*]|del)$")
CONJUNCT_INDEX_RE = re.compile(r"CONJUNCT_NAME_(?P<index>\d+)$")


def read_delimited(path: Path, delimiter: str = ",") -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [
            {key: "" if value is None else str(value).strip() for key, value in row.items()}
            for row in csv.DictReader(handle, delimiter=delimiter)
        ]


def write_tsv(path: Path, rows: Sequence[Dict[str, str]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def strip_phdr_prefix(value: str) -> str:
    value = (value or "").strip()
    return value[len(PHDR_PREFIX):] if value.startswith(PHDR_PREFIX) else value


def normalize_lookup_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (value or "").strip().lower())


class ProteinNameMapper:
    def __init__(self, gene_info_path: Path):
        self.gene_info_path = gene_info_path
        self.rows = read_delimited(gene_info_path, delimiter="\t")
        self.lookup: Dict[str, str] = {}
        self.segment_by_protein: Dict[str, str] = {}
        self.alias_rows: List[Dict[str, str]] = []
        self._build()

    def _build(self) -> None:
        for row in self.rows:
            canonical_name = (row.get("name") or "").strip()
            if not canonical_name:
                continue
            parent_name = (row.get("parent_name") or "").strip()
            segment = ""
            if canonical_name != "whole_genome":
                segment = "1"
            if parent_name and parent_name.upper() != "NULL" and parent_name != "whole_genome":
                segment = parent_name
            self.segment_by_protein[canonical_name] = segment
            alias_sources = {
                "name": canonical_name,
                "display_name": row.get("display_name", ""),
                "description": row.get("description", ""),
            }
            for source_field, source_value in alias_sources.items():
                source_value = (source_value or "").strip()
                if not source_value:
                    continue
                lookup_key = normalize_lookup_key(source_value)
                if not lookup_key:
                    continue
                self.lookup[lookup_key] = canonical_name
                self.alias_rows.append(
                    {
                        "source_value": source_value,
                        "source_field": source_field,
                        "normalized_lookup_key": lookup_key,
                        "canonical_protein_name": canonical_name,
                        "segment": segment,
                    }
                )

    def canonicalize(self, product_name: str) -> str:
        product_name = (product_name or "").strip()
        if not product_name:
            return ""
        lookup_key = normalize_lookup_key(product_name)
        return self.lookup.get(lookup_key, product_name)

    def segment_for(self, product_name: str) -> str:
        canonical_name = self.canonicalize(product_name)
        return self.segment_by_protein.get(canonical_name, "")


class HcvMutationCatalogNormalizer:
    output_fields = [
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
        "component_order",
        "source_variation_name",
        "source_phdr_ras_id",
        "alignment_name",
        "display_structure",
        "id",
        "resistance_category",
        "display_resistance_category",
        "numeric_resistance_category",
        "any_in_vitro_evidence",
        "in_vitro_max_ec50_midpoint",
        "any_in_vivo_evidence",
        "in_vivo_baseline",
        "in_vivo_treatment_emergent",
        "phdr_alignment_ras_id",
        "phdr_drug_id",
        "drug",
    ]

    def __init__(
        self,
        variation_path: Path,
        variation_metatag_path: Path,
        phdr_alignment_ras_path: Path,
        phdr_alignment_ras_drug_path: Path,
        gene_info_path: Path,
        output_path: Path,
    ):
        self.variation_path = variation_path
        self.variation_metatag_path = variation_metatag_path
        self.phdr_alignment_ras_path = phdr_alignment_ras_path
        self.phdr_alignment_ras_drug_path = phdr_alignment_ras_drug_path
        self.gene_info_path = gene_info_path
        self.output_path = output_path
        self.protein_mapper = ProteinNameMapper(gene_info_path)

    def normalize(self) -> Dict[str, Path]:
        variation_rows = read_delimited(self.variation_path)
        variation_metatag_rows = read_delimited(self.variation_metatag_path)
        alignment_rows = read_delimited(self.phdr_alignment_ras_path)
        drug_rows = read_delimited(self.phdr_alignment_ras_drug_path)

        mutation_catalog = self._build_mutation_catalog(variation_rows)
        combination_catalog = self._build_combination_catalog(variation_rows)
        combination_memberships = self._build_combination_memberships(
            variation_metatag_rows,
            mutation_catalog,
            combination_catalog,
        )
        rows = self._build_generalized_rows(
            mutation_catalog=mutation_catalog,
            combination_catalog=combination_catalog,
            combination_memberships=combination_memberships,
            alignment_contexts=self._build_alignment_context_map(alignment_rows),
            drug_rows_by_alignment=self._build_drug_row_map(drug_rows),
        )
        output_path = self.output_path
        write_tsv(output_path, rows, self.output_fields)
        return {"generalized_mutation_catalog": output_path}

    def _build_generalized_rows(
        self,
        mutation_catalog: Dict[str, Dict[str, str]],
        combination_catalog: Dict[str, Dict[str, str]],
        combination_memberships: Sequence[Dict[str, str]],
        alignment_contexts: Dict[str, List[Dict[str, str]]],
        drug_rows_by_alignment: Dict[str, List[Dict[str, str]]],
    ) -> List[Dict[str, str]]:
        rows: List[Dict[str, str]] = []
        membership_by_combo: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        for membership in combination_memberships:
            membership_by_combo[membership["combination_id"]].append(membership)
        for combo_id, memberships in membership_by_combo.items():
            membership_by_combo[combo_id] = self._sorted_rows(
                memberships,
                ("component_order", "component_mutation_id"),
            )

        for mutation in self._sorted_rows(
            mutation_catalog.values(),
            ("canonical_protein_name", "aa_position", "alt_residue", "mutation_id"),
        ):
            for context in self._contexts_for_signature(alignment_contexts, mutation["mutation_id"]):
                for drug_row in self._drug_rows_for_context(drug_rows_by_alignment, context):
                    rows.append(
                        self._build_output_row(
                            mutation_row=mutation,
                            signature_id=mutation["mutation_id"],
                            signature_kind="single",
                            combination_id="",
                            combination_size="",
                            component_order="",
                            source_variation_name=mutation["source_variation_name"],
                            source_phdr_ras_id=mutation["source_phdr_ras_id"],
                            alignment_context=context,
                            drug_row=drug_row,
                        )
                    )

        for combination in self._sorted_rows(
            combination_catalog.values(),
            ("canonical_protein_name", "combination_id"),
        ):
            memberships = membership_by_combo.get(combination["combination_id"], [])
            for context in self._contexts_for_signature(alignment_contexts, combination["combination_id"]):
                for drug_row in self._drug_rows_for_context(drug_rows_by_alignment, context):
                    for membership in memberships:
                        mutation = mutation_catalog[membership["component_mutation_id"]]
                        rows.append(
                            self._build_output_row(
                                mutation_row=mutation,
                                signature_id=combination["combination_id"],
                                signature_kind="combination",
                                combination_id=combination["combination_id"],
                                combination_size=str(len(memberships)),
                                component_order=membership["component_order"],
                                source_variation_name=combination["source_variation_name"],
                                source_phdr_ras_id=combination["source_phdr_ras_id"],
                                alignment_context=context,
                                drug_row=drug_row,
                            )
                        )

        return self._sorted_rows(
            rows,
            (
                "protein_name",
                "aa_position",
                "alt_residue",
                "mutation_id",
                "signature_id",
                "combination_id",
                "component_order",
                "alignment_name",
                "phdr_drug_id",
                "id",
            ),
        )

    @staticmethod
    def _build_alignment_context_map(alignment_rows: Iterable[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
        contexts_by_signature: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        for row in alignment_rows:
            signature_id = strip_phdr_prefix(row.get("phdr_ras_id") or row.get("id", "").rsplit(":", 1)[0])
            if not signature_id:
                continue
            contexts_by_signature[signature_id].append(
                {
                    "alignment_name": row.get("alignment_name", ""),
                    "display_structure": row.get("display_structure", ""),
                    "phdr_alignment_ras_id": row.get("id", ""),
                }
            )
        return contexts_by_signature

    @staticmethod
    def _build_drug_row_map(drug_rows: Iterable[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
        drug_rows_by_alignment: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        for row in drug_rows:
            key = row.get("phdr_alignment_ras_id", "")
            if key:
                drug_rows_by_alignment[key].append(row)
        return drug_rows_by_alignment

    @staticmethod
    def _contexts_for_signature(
        alignment_contexts: Dict[str, List[Dict[str, str]]],
        signature_id: str,
    ) -> List[Dict[str, str]]:
        contexts = alignment_contexts.get(signature_id)
        if contexts:
            return contexts
        return [{"alignment_name": "", "display_structure": "", "phdr_alignment_ras_id": ""}]

    @staticmethod
    def _drug_rows_for_context(
        drug_rows_by_alignment: Dict[str, List[Dict[str, str]]],
        context: Dict[str, str],
    ) -> List[Dict[str, str]]:
        phdr_alignment_ras_id = context.get("phdr_alignment_ras_id", "")
        if phdr_alignment_ras_id and phdr_alignment_ras_id in drug_rows_by_alignment:
            return drug_rows_by_alignment[phdr_alignment_ras_id]
        return [{}]

    @staticmethod
    def _build_output_row(
        mutation_row: Dict[str, str],
        signature_id: str,
        signature_kind: str,
        combination_id: str,
        combination_size: str,
        component_order: str,
        source_variation_name: str,
        source_phdr_ras_id: str,
        alignment_context: Dict[str, str],
        drug_row: Dict[str, str],
    ) -> Dict[str, str]:
        return {
            "protein_name": mutation_row["canonical_protein_name"],
            "segment": mutation_row.get("segment", "") or "1",
            "aa_position": mutation_row["aa_position"],
            "alt_residue": mutation_row["alt_residue"],
            "reference_accession": mutation_row["reference_accession"],
            "mutation_id": mutation_row["mutation_id"],
            "mutation_type": mutation_row["mutation_type"],
            "signature_id": signature_id,
            "signature_kind": signature_kind,
            "combination_id": combination_id,
            "combination_size": combination_size,
            "phenotype": "",
            "component_order": component_order,
            "source_variation_name": source_variation_name,
            "source_phdr_ras_id": source_phdr_ras_id,
            "alignment_name": alignment_context.get("alignment_name", ""),
            "display_structure": alignment_context.get("display_structure", ""),
            "id": drug_row.get("id", ""),
            "resistance_category": drug_row.get("resistance_category", ""),
            "display_resistance_category": drug_row.get("display_resistance_category", ""),
            "numeric_resistance_category": drug_row.get("numeric_resistance_category", ""),
            "any_in_vitro_evidence": drug_row.get("any_in_vitro_evidence", ""),
            "in_vitro_max_ec50_midpoint": drug_row.get("in_vitro_max_ec50_midpoint", ""),
            "any_in_vivo_evidence": drug_row.get("any_in_vivo_evidence", ""),
            "in_vivo_baseline": drug_row.get("in_vivo_baseline", ""),
            "in_vivo_treatment_emergent": drug_row.get("in_vivo_treatment_emergent", ""),
            "phdr_alignment_ras_id": drug_row.get("phdr_alignment_ras_id", alignment_context.get("phdr_alignment_ras_id", "")),
            "phdr_drug_id": drug_row.get("phdr_drug_id", ""),
            "drug": drug_row.get("phdr_drug_id", ""),
        }

    def _build_mutation_catalog(self, variation_rows: Iterable[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
        mutation_catalog: Dict[str, Dict[str, str]] = {}
        for row in variation_rows:
            mutation_type = row.get("type", "")
            if mutation_type == "conjunction":
                continue
            mutation_id = self._signature_id_from_variation(row)
            if not mutation_id:
                continue
            canonical_protein = self.protein_mapper.canonicalize(row.get("feature_name", ""))
            segment = self.protein_mapper.segment_for(canonical_protein)
            aa_position, alt_residue = self._parse_single_mutation_id(mutation_id)
            mutation_catalog[mutation_id] = {
                "mutation_id": mutation_id,
                "canonical_protein_name": canonical_protein,
                "segment": segment,
                "aa_position": str(aa_position),
                "alt_residue": alt_residue,
                "reference_accession": row.get("ref_seq_name", ""),
                "mutation_type": mutation_type,
                "source_variation_name": row.get("name", ""),
                "source_phdr_ras_id": mutation_id,
                "hcv_alignment_buckets": "",
                "drug_annotation_ids": "",
                "curated_combination_ids": "",
            }
        return mutation_catalog

    def _build_combination_catalog(self, variation_rows: Iterable[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
        combination_catalog: Dict[str, Dict[str, str]] = {}
        for row in variation_rows:
            mutation_type = row.get("type", "")
            if mutation_type != "conjunction":
                continue
            combination_id = self._signature_id_from_variation(row)
            if not combination_id:
                continue
            canonical_protein = self.protein_mapper.canonicalize(row.get("feature_name", ""))
            segment = self.protein_mapper.segment_for(canonical_protein)
            combination_catalog[combination_id] = {
                "combination_id": combination_id,
                "canonical_protein_name": canonical_protein,
                "segment": segment,
                "reference_accession": row.get("ref_seq_name", ""),
                "mutation_type": mutation_type,
                "source_variation_name": row.get("name", ""),
                "source_phdr_ras_id": combination_id,
                "hcv_alignment_buckets": "",
                "drug_annotation_ids": "",
                "component_mutation_ids": "",
            }
        return combination_catalog

    def _build_combination_memberships(
        self,
        variation_metatag_rows: Iterable[Dict[str, str]],
        mutation_catalog: Dict[str, Dict[str, str]],
        combination_catalog: Dict[str, Dict[str, str]],
    ) -> List[Dict[str, str]]:
        memberships: List[Dict[str, str]] = []
        for row in variation_metatag_rows:
            combination_id = strip_phdr_prefix(row.get("variation_name", ""))
            component_id = strip_phdr_prefix(row.get("metatag_value", ""))
            if not combination_id or combination_id not in combination_catalog:
                continue
            match = CONJUNCT_INDEX_RE.match(row.get("metatag_name", ""))
            component_order = int(match.group("index")) if match else 999
            component_details = mutation_catalog.get(component_id)
            if component_details is None:
                protein_name, aa_position, alt_residue = self._parse_signature_component(component_id)
                canonical_protein = self.protein_mapper.canonicalize(protein_name)
                component_details = {
                    "mutation_id": component_id,
                    "canonical_protein_name": canonical_protein,
                    "segment": self.protein_mapper.segment_for(canonical_protein),
                    "aa_position": str(aa_position),
                    "alt_residue": alt_residue,
                    "reference_accession": row.get("ref_seq_name", ""),
                }
                mutation_catalog[component_id] = {
                    "mutation_id": component_id,
                    "canonical_protein_name": canonical_protein,
                    "segment": component_details["segment"],
                    "aa_position": component_details["aa_position"],
                    "alt_residue": component_details["alt_residue"],
                    "reference_accession": row.get("ref_seq_name", ""),
                    "mutation_type": "derived_component",
                    "source_variation_name": component_id,
                    "source_phdr_ras_id": component_id,
                    "hcv_alignment_buckets": "",
                    "drug_annotation_ids": "",
                    "curated_combination_ids": "",
                }
            memberships.append(
                {
                    "combination_id": combination_id,
                    "component_order": str(component_order),
                    "component_mutation_id": component_id,
                    "canonical_protein_name": component_details["canonical_protein_name"],
                    "segment": component_details["segment"],
                    "aa_position": component_details["aa_position"],
                    "alt_residue": component_details["alt_residue"],
                    "reference_accession": component_details.get("reference_accession", ""),
                }
            )
        return memberships

    def _build_drug_interpretations(
        self,
        drug_rows: Iterable[Dict[str, str]],
        mutation_catalog: Dict[str, Dict[str, str]],
        combination_catalog: Dict[str, Dict[str, str]],
    ) -> List[Dict[str, str]]:
        interpretations: List[Dict[str, str]] = []
        for row in drug_rows:
            source_alignment_ras_id = row.get("phdr_alignment_ras_id", "")
            signature_id, hcv_alignment_bucket = self._parse_alignment_ras_id(source_alignment_ras_id)
            if signature_id in combination_catalog:
                signature_kind = "combination"
                canonical_protein = combination_catalog[signature_id]["canonical_protein_name"]
                segment = combination_catalog[signature_id]["segment"]
                reference_accession = combination_catalog[signature_id]["reference_accession"]
            else:
                signature_kind = "single"
                mutation = mutation_catalog.get(signature_id)
                if mutation is None:
                    protein_name, aa_position, alt_residue = self._parse_signature_component(signature_id)
                    canonical_protein = self.protein_mapper.canonicalize(protein_name)
                    mutation = {
                        "mutation_id": signature_id,
                        "canonical_protein_name": canonical_protein,
                        "segment": self.protein_mapper.segment_for(canonical_protein),
                        "aa_position": str(aa_position),
                        "alt_residue": alt_residue,
                        "reference_accession": "",
                        "mutation_type": "derived_from_drug_table",
                        "source_variation_name": signature_id,
                        "source_phdr_ras_id": signature_id,
                        "hcv_alignment_buckets": "",
                        "drug_annotation_ids": "",
                        "curated_combination_ids": "",
                    }
                    mutation_catalog[signature_id] = mutation
                canonical_protein = mutation["canonical_protein_name"]
                segment = mutation["segment"]
                reference_accession = mutation.get("reference_accession", "")
            interpretations.append(
                {
                    "drug_annotation_id": row.get("id", ""),
                    "signature_id": signature_id,
                    "signature_kind": signature_kind,
                    "canonical_protein_name": canonical_protein,
                    "segment": segment,
                    "reference_accession": reference_accession,
                    "hcv_alignment_bucket": hcv_alignment_bucket,
                    "drug_name": row.get("phdr_drug_id", ""),
                    "resistance_category": row.get("resistance_category", ""),
                    "display_resistance_category": row.get("display_resistance_category", ""),
                    "numeric_resistance_category": row.get("numeric_resistance_category", ""),
                    "any_in_vitro_evidence": row.get("any_in_vitro_evidence", ""),
                    "in_vitro_max_ec50_midpoint": row.get("in_vitro_max_ec50_midpoint", ""),
                    "any_in_vivo_evidence": row.get("any_in_vivo_evidence", ""),
                    "in_vivo_baseline": row.get("in_vivo_baseline", ""),
                    "in_vivo_treatment_emergent": row.get("in_vivo_treatment_emergent", ""),
                    "source_phdr_alignment_ras_id": source_alignment_ras_id,
                }
            )
        return interpretations

    @staticmethod
    def _sorted_rows(rows: Iterable[Dict[str, str]], sort_keys: Tuple[str, ...]) -> List[Dict[str, str]]:
        return sorted(
            rows,
            key=lambda row: tuple((row.get(key, "") or "") for key in sort_keys),
        )

    @staticmethod
    def _signature_id_from_variation(row: Dict[str, str]) -> str:
        return strip_phdr_prefix(row.get("phdr_ras_id") or row.get("name", ""))

    @staticmethod
    def _parse_alignment_ras_id(value: str) -> Tuple[str, str]:
        value = (value or "").strip()
        if not value:
            return "", ""
        signature_id, alignment_bucket = value.rsplit(":", 1)
        return strip_phdr_prefix(signature_id), alignment_bucket

    def _parse_single_mutation_id(self, mutation_id: str) -> Tuple[int, str]:
        protein_name, aa_position, alt_residue = self._parse_signature_component(mutation_id)
        if "+" in mutation_id:
            raise ValueError(f"Expected a single mutation id but found combination: {mutation_id}")
        return aa_position, alt_residue

    @staticmethod
    def _parse_signature_component(component_id: str) -> Tuple[str, int, str]:
        protein_name, token = component_id.split(":", 1)
        match = MUTATION_TOKEN_RE.match(token)
        if not match:
            raise ValueError(f"Unsupported mutation token format: {component_id}")
        return protein_name, int(match.group("position")), match.group("alt_residue")


def build_parser() -> ArgumentParser:
    repo_root = Path(__file__).resolve().parents[1]
    hcv_tables_dir = repo_root / "generic" / "hcv" / "Tables"

    parser = ArgumentParser(description="Normalize HCV mutation assets into a single general catalog TSV")
    parser.add_argument("--variation", default=str(hcv_tables_dir / "variation.csv"))
    parser.add_argument("--variation_metatag", default=str(hcv_tables_dir / "variation_metatag.csv"))
    parser.add_argument("--phdr_alignment_ras", default=str(hcv_tables_dir / "phdr_alignment_ras.csv"))
    parser.add_argument("--phdr_alignment_ras_drug", default=str(hcv_tables_dir / "phdr_alignment_ras_drug.csv"))
    parser.add_argument("--gene_info", default=str(hcv_tables_dir / "gene_info.tsv"))
    parser.add_argument(
        "--output_path",
        default=str(hcv_tables_dir / "generalized_mutation_catalog.tsv"),
        help="Path to the normalized TSV handoff file",
    )
    return parser


def main(args=None) -> Dict[str, Path]:
    parser = build_parser()
    parsed_args = parser.parse_args(args=args)
    normalizer = HcvMutationCatalogNormalizer(
        variation_path=Path(parsed_args.variation),
        variation_metatag_path=Path(parsed_args.variation_metatag),
        phdr_alignment_ras_path=Path(parsed_args.phdr_alignment_ras),
        phdr_alignment_ras_drug_path=Path(parsed_args.phdr_alignment_ras_drug),
        gene_info_path=Path(parsed_args.gene_info),
        output_path=Path(parsed_args.output_path),
    )
    outputs = normalizer.normalize()
    for name, path in outputs.items():
        print(f"{name}\t{path}")
    return outputs


if __name__ == "__main__":
    main()
