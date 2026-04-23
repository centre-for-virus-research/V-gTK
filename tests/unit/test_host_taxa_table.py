from pathlib import Path

from HostTaxaTable import HostTaxaTable


def test_load_taxa_ids_skips_rows_missing_host_taxa_id(tmp_path: Path):
    gb_matrix = tmp_path / "gB_matrix_validated.tsv"
    gb_matrix.write_text(
        "primary_accession\thost_taxa_id\n"
        "ACC1\t9606\n"
        "ACC2\n"
        "ACC3\tNA\n"
        "ACC4\t10090\n"
        "ACC5\t9606\n",
        encoding="utf-8",
    )

    table = HostTaxaTable(
        str(gb_matrix),
        "HostTaxa",
        "names.dmp",
        "nodes.dmp",
        str(tmp_path),
        "Host_taxa.tsv",
        "Host_taxa_children.tsv",
        "Host_taxa_lineage.tsv",
    )

    assert table.load_taxa_ids_from_tsv() == [9606, 10090]