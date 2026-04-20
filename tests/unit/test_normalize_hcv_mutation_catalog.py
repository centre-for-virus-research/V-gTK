import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from NormalizeHcvMutationCatalog import HcvMutationCatalogNormalizer, main  # type: ignore[reportMissingImports]


DRUG_SOURCE_COLUMNS = [
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
]

KEY_COLUMNS = [
    "protein_name",
    "segment",
    "aa_position",
    "alt_residue",
    "reference_accession",
]


def write_csv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def row_for(rows, **criteria):
    matches = [row for row in rows if all(row.get(key, "") == value for key, value in criteria.items())]
    assert len(matches) == 1, f"Expected exactly one row for {criteria}, found {len(matches)}"
    return matches[0]


def test_normalize_hcv_mutation_catalog_flattens_to_keyed_rows(tmp_path: Path):
    tables_dir = tmp_path / "Tables"
    output_path = tables_dir / "generalized_mutation_catalog.tsv"

    write_tsv(
        tables_dir / "gene_info.tsv",
        [
            {
                "description": "Non-structural protein 3",
                "display_name": "NS3",
                "name": "NS3",
                "parent_name": "whole_genome",
            },
            {
                "description": "Whole genome",
                "display_name": "Whole genome",
                "name": "whole_genome",
                "parent_name": "NULL",
            },
        ],
        ["description", "display_name", "name", "parent_name"],
    )
    write_csv(
        tables_dir / "variation.csv",
        [
            {
                "description": "",
                "display_name": "",
                "feature_name": "NS3",
                "name": "phdr_ras:NS3:107I",
                "ref_end": "3740",
                "ref_seq_name": "REF_MASTER_NC_004102",
                "ref_start": "3738",
                "type": "aminoAcidSimplePolymorphism",
                "phdr_ras_id": "NS3:107I",
            },
            {
                "description": "",
                "display_name": "",
                "feature_name": "NS3",
                "name": "phdr_ras:NS3:155K",
                "ref_end": "3884",
                "ref_seq_name": "REF_MASTER_NC_004102",
                "ref_start": "3882",
                "type": "aminoAcidSimplePolymorphism",
                "phdr_ras_id": "NS3:155K",
            },
            {
                "description": "",
                "display_name": "",
                "feature_name": "NS3",
                "name": "phdr_ras:NS3:107I+155K",
                "ref_end": "",
                "ref_seq_name": "REF_MASTER_NC_004102",
                "ref_start": "",
                "type": "conjunction",
                "phdr_ras_id": "NS3:107I+155K",
            },
        ],
        [
            "description",
            "display_name",
            "feature_name",
            "name",
            "ref_end",
            "ref_seq_name",
            "ref_start",
            "type",
            "phdr_ras_id",
        ],
    )
    write_csv(
        tables_dir / "variation_metatag.csv",
        [
            {
                "feature_name": "NS3",
                "metatag_name": "CONJUNCT_NAME_1",
                "metatag_value": "phdr_ras:NS3:107I",
                "ref_seq_name": "REF_MASTER_NC_004102",
                "variation_name": "phdr_ras:NS3:107I+155K",
            },
            {
                "feature_name": "NS3",
                "metatag_name": "CONJUNCT_NAME_2",
                "metatag_value": "phdr_ras:NS3:155K",
                "ref_seq_name": "REF_MASTER_NC_004102",
                "variation_name": "phdr_ras:NS3:107I+155K",
            },
        ],
        ["feature_name", "metatag_name", "metatag_value", "ref_seq_name", "variation_name"],
    )
    write_csv(
        tables_dir / "phdr_alignment_ras.csv",
        [
            {
                "id": "NS3:107I:AL_1a",
                "display_structure": "V107I",
                "phdr_ras_id": "NS3:107I",
                "alignment_name": "AL_1a",
            },
            {
                "id": "NS3:107I+155K:AL_1a",
                "display_structure": "V107I+R155K",
                "phdr_ras_id": "NS3:107I+155K",
                "alignment_name": "AL_1a",
            },
        ],
        ["id", "display_structure", "phdr_ras_id", "alignment_name"],
    )
    write_csv(
        tables_dir / "phdr_alignment_ras_drug.csv",
        [
            {
                "id": "NS3:107I:AL_1a:grazoprevir",
                "resistance_category": "category_I",
                "display_resistance_category": "I",
                "numeric_resistance_category": "1",
                "any_in_vitro_evidence": "1",
                "in_vitro_max_ec50_midpoint": "0.8",
                "any_in_vivo_evidence": "1",
                "in_vivo_baseline": "1",
                "in_vivo_treatment_emergent": "1",
                "phdr_alignment_ras_id": "NS3:107I:AL_1a",
                "phdr_drug_id": "grazoprevir",
            },
            {
                "id": "NS3:107I+155K:AL_1a:voxilaprevir",
                "resistance_category": "category_III",
                "display_resistance_category": "III",
                "numeric_resistance_category": "3",
                "any_in_vitro_evidence": "1",
                "in_vitro_max_ec50_midpoint": "5.0",
                "any_in_vivo_evidence": "0",
                "in_vivo_baseline": "",
                "in_vivo_treatment_emergent": "",
                "phdr_alignment_ras_id": "NS3:107I+155K:AL_1a",
                "phdr_drug_id": "voxilaprevir",
            },
        ],
        DRUG_SOURCE_COLUMNS,
    )

    normalizer = HcvMutationCatalogNormalizer(
        variation_path=tables_dir / "variation.csv",
        variation_metatag_path=tables_dir / "variation_metatag.csv",
        phdr_alignment_ras_path=tables_dir / "phdr_alignment_ras.csv",
        phdr_alignment_ras_drug_path=tables_dir / "phdr_alignment_ras_drug.csv",
        gene_info_path=tables_dir / "gene_info.tsv",
        output_path=output_path,
    )

    outputs = normalizer.normalize()
    rows = read_tsv(outputs["generalized_mutation_catalog"])

    assert rows
    assert all(column in rows[0] for column in DRUG_SOURCE_COLUMNS)
    assert "drug" in rows[0]
    assert all(row["segment"] == "1" for row in rows)
    assert all(row[column] for row in rows for column in KEY_COLUMNS)
    assert all(row["segment"] != "whole_genome" for row in rows)

    single_row = row_for(rows, signature_id="NS3:107I", mutation_id="NS3:107I")
    assert single_row["protein_name"] == "NS3"
    assert single_row["alignment_name"] == "AL_1a"
    assert single_row["display_structure"] == "V107I"
    assert single_row["id"] == "NS3:107I:AL_1a:grazoprevir"
    assert single_row["phdr_alignment_ras_id"] == "NS3:107I:AL_1a"
    assert single_row["phdr_drug_id"] == "grazoprevir"
    assert single_row["drug"] == "grazoprevir"
    assert single_row["resistance_category"] == "category_I"

    combination_component_1 = row_for(
        rows,
        signature_id="NS3:107I+155K",
        mutation_id="NS3:107I",
        component_order="1",
    )
    assert combination_component_1["signature_kind"] == "combination"
    assert combination_component_1["combination_id"] == "NS3:107I+155K"
    assert combination_component_1["combination_size"] == "2"
    assert combination_component_1["phdr_drug_id"] == "voxilaprevir"
    assert combination_component_1["drug"] == "voxilaprevir"
    assert combination_component_1["resistance_category"] == "category_III"

    combination_component_2 = row_for(
        rows,
        signature_id="NS3:107I+155K",
        mutation_id="NS3:155K",
        component_order="2",
    )
    assert combination_component_2["signature_kind"] == "combination"
    assert combination_component_2["phdr_drug_id"] == "voxilaprevir"
    assert combination_component_2["drug"] == "voxilaprevir"

    undrugged_single = row_for(rows, signature_id="NS3:155K", mutation_id="NS3:155K")
    assert undrugged_single["alignment_name"] == ""
    assert undrugged_single["phdr_drug_id"] == ""
    assert undrugged_single["drug"] == ""
    assert undrugged_single["resistance_category"] == ""


def test_normalize_hcv_mutation_catalog_cli_sets_segment_to_one_and_maps_alias(tmp_path: Path):
    tables_dir = tmp_path / "Tables"
    output_path = tables_dir / "out.tsv"

    write_tsv(
        tables_dir / "gene_info.tsv",
        [
            {
                "description": "Non-structural protein 5A",
                "display_name": "NS5A",
                "name": "NS5A",
                "parent_name": "whole_genome",
            }
        ],
        ["description", "display_name", "name", "parent_name"],
    )
    write_csv(
        tables_dir / "variation.csv",
        [
            {
                "description": "",
                "display_name": "",
                "feature_name": "Non-structural protein 5A",
                "name": "phdr_ras:NS5A:29del",
                "ref_end": "6344",
                "ref_seq_name": "REF_MASTER_NC_004102",
                "ref_start": "6342",
                "type": "aminoAcidDeletion",
                "phdr_ras_id": "NS5A:29del",
            }
        ],
        [
            "description",
            "display_name",
            "feature_name",
            "name",
            "ref_end",
            "ref_seq_name",
            "ref_start",
            "type",
            "phdr_ras_id",
        ],
    )
    write_csv(tables_dir / "variation_metatag.csv", [], ["feature_name", "metatag_name", "metatag_value", "ref_seq_name", "variation_name"])
    write_csv(tables_dir / "phdr_alignment_ras.csv", [], ["id", "display_structure", "phdr_ras_id", "alignment_name"])
    write_csv(tables_dir / "phdr_alignment_ras_drug.csv", [], DRUG_SOURCE_COLUMNS)

    main(
        [
            "--variation",
            str(tables_dir / "variation.csv"),
            "--variation_metatag",
            str(tables_dir / "variation_metatag.csv"),
            "--phdr_alignment_ras",
            str(tables_dir / "phdr_alignment_ras.csv"),
            "--phdr_alignment_ras_drug",
            str(tables_dir / "phdr_alignment_ras_drug.csv"),
            "--gene_info",
            str(tables_dir / "gene_info.tsv"),
            "--output_path",
            str(output_path),
        ]
    )

    rows = read_tsv(output_path)
    assert len(rows) == 1
    assert rows[0]["protein_name"] == "NS5A"
    assert rows[0]["segment"] == "1"
    assert rows[0]["aa_position"] == "29"
    assert rows[0]["alt_residue"] == "del"
    assert rows[0]["reference_accession"] == "REF_MASTER_NC_004102"
    assert rows[0]["phdr_drug_id"] == ""
