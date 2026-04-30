from ValidateMatrix import ValidateMatrix


def test_resolve_host_matches_taxonomy_case_insensitively(tmp_path):
    validator = ValidateMatrix(
        url="https://example.invalid/taxdump.tar.gz",
        taxa_path="Taxa",
        base_dir=str(tmp_path),
        output_dir="Validate-matrix",
        gb_matrix="gB_matrix_raw.tsv",
        country_file="m49_country.csv",
        assets="assets",
        host_map="host_mapping.tsv",
        country_map="country_mapping.tsv",
    )

    taxa_dict = ({"homo sapiens": "9606"}, {"9606": "Homo sapiens"})

    assert validator.resolve_host("Homo sapiens", taxa_dict, {}) == (
        "Yes",
        "9606",
        "Homo sapiens",
        "",
    )


def test_resolve_host_uses_case_insensitive_mapping_then_taxonomy(tmp_path):
    validator = ValidateMatrix(
        url="https://example.invalid/taxdump.tar.gz",
        taxa_path="Taxa",
        base_dir=str(tmp_path),
        output_dir="Validate-matrix",
        gb_matrix="gB_matrix_raw.tsv",
        country_file="m49_country.csv",
        assets="assets",
        host_map="host_mapping.tsv",
        country_map="country_mapping.tsv",
    )

    taxa_dict = ({"canis lupus familiaris": "9615"}, {"9615": "Canis lupus familiaris"})
    host_map = {"dog": "Canis lupus familiaris"}

    assert validator.resolve_host("Dog", taxa_dict, host_map) == (
        "Yes",
        "9615",
        "Canis lupus familiaris",
        "Host mapped from mapping file: 'Dog' -> 'Canis lupus familiaris'",
    )