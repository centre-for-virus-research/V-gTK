import os
import sqlite3
import sys

import pytest

#/home3/oml4h/miniconda3/etc/profile.d/conda.sh && conda activate vgtk && python scripts/PlotHcvMutationTrends.py --db /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/HCV_full.db --output_dir /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/plots_mutations --min_yearly_sequences 25

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../scripts")))

import PlotHcvMutationTrends


@pytest.fixture
def hcv_plot_db(tmp_path):
    db_path = tmp_path / "hcv_plot_test.db"
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE meta_data (
                primary_accession TEXT,
                collection_year TEXT,
                exclusion_status TEXT
            )
            """
        )
        cur.executemany(
            "INSERT INTO meta_data VALUES (?, ?, ?)",
            [
                ("A_2020_1", "2020", ""),
                ("A_2020_2", "2020", ""),
                ("A_2021_1", "2021", ""),
                ("A_2021_2", "2021", ""),
                ("A_2021_excluded", "2021", "1"),
                ("A_missing_year", "", ""),
            ],
        )

        cur.execute(
            """
            CREATE TABLE sequence_relevant_mutation_summary (
                primary_accession TEXT,
                relevant_mutations_present TEXT,
                total_relevant_mutation_count INTEGER
            )
            """
        )
        cur.executemany(
            "INSERT INTO sequence_relevant_mutation_summary VALUES (?, ?, ?)",
            [
                ("A_2020_1", "mut1;mut2", 2),
                ("A_2021_1", "mut3", 1),
                ("A_2021_excluded", "mut4;mut5;mut6", 3),
            ],
        )

        cur.execute(
            """
            CREATE TABLE completed_signatures_only (
                primary_accession TEXT,
                signature_id TEXT,
                signature_kind TEXT
            )
            """
        )
        cur.executemany(
            "INSERT INTO completed_signatures_only VALUES (?, ?, ?)",
            [
                ("A_2020_1", "sig_drug_a", "single"),
                ("A_2020_1", "sig_drug_b", "single"),
                ("A_2021_1", "sig_drug_a", "single"),
                ("A_2021_2", "sig_drug_a_combo", "combination"),
                ("A_2021_excluded", "sig_drug_b", "single"),
            ],
        )

        cur.execute(
            """
            CREATE TABLE mutation_catalog (
                mutation_id TEXT,
                protein_name TEXT,
                segment TEXT,
                aa_position TEXT,
                alt_residue TEXT,
                reference_accession TEXT,
                mutation_type TEXT,
                signature_id TEXT,
                signature_kind TEXT,
                combination_id TEXT,
                combination_size TEXT,
                phenotype TEXT,
                resistance_category TEXT,
                drug TEXT
            )
            """
        )
        cur.executemany(
            "INSERT INTO mutation_catalog VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("mut1", "NS3", "1", "1", "A", "REF", "snp", "sig_drug_a", "single", "", "", "drug_resistance", "I", "DrugA"),
                ("mut2", "NS3", "1", "2", "T", "REF", "snp", "sig_drug_b", "single", "", "", "drug_resistance", "I", "DrugB"),
                ("mut3", "NS5A", "", "3", "G", "REF", "snp", "sig_drug_a_combo", "combination", "combo_a", "2", "drug_resistance", "II", "DrugA"),
            ],
        )
        conn.commit()
    finally:
        conn.close()
    return db_path


def test_generate_plots_builds_expected_yearly_summaries(hcv_plot_db, tmp_path):
    output_dir = tmp_path / "plots_mutations"
    results = PlotHcvMutationTrends.generate_plots(hcv_plot_db, output_dir, min_yearly_sequences=2)

    mutation_summary = results["yearly_mutation_summary"].set_index("collection_year")
    assert mutation_summary.loc[2020, "sequence_count"] == 2
    assert mutation_summary.loc[2020, "mean_relevant_mutation_count"] == 1.0
    assert mutation_summary.loc[2021, "sequence_count"] == 2
    assert mutation_summary.loc[2021, "mean_relevant_mutation_count"] == 0.5

    signature_summary = results["yearly_signature_summary"].set_index("collection_year")
    assert signature_summary.loc[2020, "mean_completed_signature_count"] == 1.0
    assert signature_summary.loc[2021, "mean_completed_signature_count"] == 1.0

    drug_summary = results["yearly_drug_proportion"].sort_values(["drug", "collection_year"]).reset_index(drop=True)
    assert set(drug_summary["drug"].tolist()) == {"DrugA", "DrugB"}
    assert len(drug_summary) == 4

    drug_a_2020 = drug_summary[(drug_summary["drug"] == "DrugA") & (drug_summary["collection_year"] == 2020)].iloc[0]
    assert drug_a_2020["sequences_with_signature"] == 1
    assert drug_a_2020["total_sequences"] == 2
    assert drug_a_2020["proportion"] == 0.5
    assert 0.0 <= drug_a_2020["smoothed_proportion"] <= 1.0

    drug_a_2021 = drug_summary[(drug_summary["drug"] == "DrugA") & (drug_summary["collection_year"] == 2021)].iloc[0]
    assert drug_a_2021["sequences_with_signature"] == 2
    assert drug_a_2021["total_sequences"] == 2
    assert drug_a_2021["proportion"] == 1.0

    drug_b_2021 = drug_summary[(drug_summary["drug"] == "DrugB") & (drug_summary["collection_year"] == 2021)].iloc[0]
    assert drug_b_2021["sequences_with_signature"] == 0
    assert drug_b_2021["total_sequences"] == 2
    assert drug_b_2021["proportion"] == 0.0
    assert 0.0 <= drug_b_2021["smoothed_proportion"] <= 1.0

    yearly_sample_size = results["yearly_sample_size"].set_index("collection_year")
    assert yearly_sample_size.loc[2020, "sequence_count"] == 2
    assert yearly_sample_size.loc[2021, "sequence_count"] == 2

    for expected_path in [
        output_dir / "hcv_mutation_trends.png",
        output_dir / "mean_relevant_mutation_count_by_year.png",
        output_dir / "mean_completed_signature_count_by_year.png",
        output_dir / "antiviral_signature_proportion_by_year.png",
        output_dir / "antiviral_signature_proportion_by_year_with_ci.png",
        output_dir / "yearly_relevant_mutation_summary.tsv",
        output_dir / "yearly_completed_signature_summary.tsv",
        output_dir / "yearly_antiviral_signature_proportion.tsv",
        output_dir / "yearly_sample_size.tsv",
    ]:
        assert expected_path.exists()