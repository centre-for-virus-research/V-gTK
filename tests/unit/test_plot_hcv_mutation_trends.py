import importlib
import os
import sqlite3
import sys

import matplotlib.pyplot as plt
import pandas as pd
import pytest

#/home3/oml4h/miniconda3/etc/profile.d/conda.sh && conda activate vgtk && python scripts/PlotHcvMutationTrends.py --db /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/HCV_full.db --output_dir /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless/plots_mutations --min_yearly_sequences 25

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../scripts")))

PlotHcvMutationTrends = importlib.import_module("PlotHcvMutationTrends")


def _create_base_plot_db(db_path, *, include_years=True, include_alignment=True, blank_alignment_name=False):
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE meta_data (
                primary_accession TEXT,
                collection_year TEXT,
                exclusion_status TEXT,
                accession_type TEXT
            )
            """
        )
        meta_rows = [
            ("A_2020_1", "2020" if include_years else "", "", ""),
            ("A_2020_2", "2020" if include_years else "", "", ""),
            ("A_2021_1", "2021" if include_years else "", "", ""),
            ("A_2021_2", "2021" if include_years else "", "", ""),
            ("A_2021_excluded", "2021" if include_years else "", "1", ""),
            ("A_missing_year", "", "", ""),
            ("REF1", "", "", "master"),
        ]
        cur.executemany("INSERT INTO meta_data VALUES (?, ?, ?, ?)", meta_rows)

        cur.execute(
            """
            CREATE TABLE features (
                accession TEXT,
                product TEXT,
                segment TEXT,
                cds_start INTEGER,
                cds_end INTEGER
            )
            """
        )
        cur.executemany(
            "INSERT INTO features VALUES (?, ?, ?, ?, ?)",
            [
                ("REF1", "NS3", "1", 1, 6),
                ("REF1", "NS5A", "1", 7, 9),
            ],
        )

        if include_alignment:
            cur.execute(
                """
                CREATE TABLE sequence_alignment (
                    sequence_id TEXT,
                    primary_accession TEXT,
                    alignment_name TEXT,
                    alignment TEXT
                )
                """
            )
            alignment_name = "" if blank_alignment_name else "REF1"
            cur.executemany(
                "INSERT INTO sequence_alignment VALUES (?, ?, ?, ?)",
                [
                    ("REF1", "REF1", "REF1", "ATAATCTAC"),
                    ("A_2020_1", "A_2020_1", alignment_name, "ATAATCTAC"),
                    ("A_2020_2", "A_2020_2", alignment_name, "ATAATC---"),
                    ("A_2021_1", "A_2021_1", alignment_name, "ATAATCTAC"),
                    ("A_2021_2", "A_2021_2", alignment_name, "ATAATCTAC"),
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
                ("mut1", "NS3", "1", "1", "A", "REF1", "snp", "sig_drug_a", "single", "", "", "drug_resistance", "I", "DrugA"),
                ("mut2", "NS3", "1", "2", "T", "REF1", "snp", "sig_drug_b", "single", "", "", "drug_resistance", "I", "DrugB"),
                ("mut3", "NS5A", "1", "1", "G", "REF1", "snp", "sig_drug_a_combo", "combination", "combo_a", "2", "drug_resistance", "II", "DrugA"),
            ],
        )
        conn.commit()
    finally:
        conn.close()


@pytest.fixture
def hcv_plot_db(tmp_path):
    db_path = tmp_path / "hcv_plot_test.db"
    _create_base_plot_db(db_path)
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

    signature_frequency = results["yearly_signature_frequency"].set_index(["signature_id", "collection_year"])

    sig_drug_a_2020 = signature_frequency.loc[("sig_drug_a", 2020)]
    assert sig_drug_a_2020["drug"] == "DrugA"
    assert sig_drug_a_2020["callable_sequences"] == 2
    assert sig_drug_a_2020["sequences_with_signature"] == 1
    assert sig_drug_a_2020["proportion"] == 0.5
    assert 0.0 <= sig_drug_a_2020["ci_lower"] <= sig_drug_a_2020["ci_upper"] <= 1.0

    sig_drug_b_2021 = signature_frequency.loc[("sig_drug_b", 2021)]
    assert sig_drug_b_2021["callable_sequences"] == 2
    assert sig_drug_b_2021["sequences_with_signature"] == 0
    assert sig_drug_b_2021["proportion"] == 0.0

    sig_combo_2020 = signature_frequency.loc[("sig_drug_a_combo", 2020)]
    assert sig_combo_2020["callable_sequences"] == 1
    assert sig_combo_2020["sequences_with_signature"] == 0
    assert sig_combo_2020["proportion"] == 0.0

    sig_combo_2021 = signature_frequency.loc[("sig_drug_a_combo", 2021)]
    assert sig_combo_2021["protein_name"] == "NS5A"
    assert sig_combo_2021["callable_sequences"] == 2
    assert sig_combo_2021["sequences_with_signature"] == 1
    assert sig_combo_2021["proportion"] == 0.5
    assert 0.0 <= sig_combo_2021["smoothed_proportion"] <= 1.0

    gene_signature_frequency = results["yearly_gene_signature_frequency"].set_index(["signature_id", "collection_year"])
    assert gene_signature_frequency.loc[("sig_drug_a", 2020), "callable_sequences"] == 2
    assert gene_signature_frequency.loc[("sig_drug_a_combo", 2020), "callable_sequences"] == 1
    assert gene_signature_frequency.loc[("sig_drug_a_combo", 2021), "resistance_category"] == "II"

    priority_drugs = results["priority_drugs"]
    assert priority_drugs["drug"].tolist() == ["DrugA", "DrugB"]
    assert priority_drugs.iloc[0]["resistance_score"] >= priority_drugs.iloc[1]["resistance_score"]

    yearly_sample_size = results["yearly_sample_size"].set_index("collection_year")
    assert yearly_sample_size.loc[2020, "sequence_count"] == 2
    assert yearly_sample_size.loc[2021, "sequence_count"] == 2

    for expected_path in [
        output_dir / "hcv_mutation_trends.png",
        output_dir / "mean_relevant_mutation_count_by_year.png",
        output_dir / "mean_completed_signature_count_by_year.png",
        output_dir / "signature_frequency_by_year.png",
        output_dir / "signature_frequency_by_year_with_ci.png",
        output_dir / "yearly_relevant_mutation_summary.tsv",
        output_dir / "yearly_completed_signature_summary.tsv",
        output_dir / "yearly_signature_frequency.tsv",
        output_dir / "yearly_gene_signature_frequency.tsv",
        output_dir / "priority_antiviral_selection.tsv",
        output_dir / "yearly_sample_size.tsv",
        output_dir / "plots_by_gene" / "ns3_signature_frequency_by_year.png",
        output_dir / "plots_by_gene" / "ns5a_signature_frequency_by_year.png",
        output_dir / "plots_by_drug" / "druga_priority_signature_frequency_by_year.png",
        output_dir / "plots_by_drug" / "druga_category_i_signature_frequency_by_year.png",
        output_dir / "plots_by_drug" / "drugb_priority_signature_frequency_by_year.png",
        output_dir / "plots_by_drug" / "drugb_category_i_signature_frequency_by_year.png",
    ]:
        assert expected_path.exists()

    category_i_plot_paths = results["plot_paths"]["priority_drug_category_i_signature_frequency"]
    assert set(category_i_plot_paths) == {"DrugA", "DrugB"}


def test_smooth_weighted_proportions_handles_empty_and_singleton_inputs():
    assert PlotHcvMutationTrends._smooth_weighted_proportions([], [], []).size == 0

    smoothed = PlotHcvMutationTrends._smooth_weighted_proportions([2020], [0.75], [10])
    assert smoothed.tolist() == [0.75]


def test_sequence_has_callable_signature_coverage_rejects_gaps_and_ambiguity():
    assert PlotHcvMutationTrends._sequence_has_callable_signature_coverage("ATAATCTAC", [(0, 1, 2), (6, 7, 8)])
    assert not PlotHcvMutationTrends._sequence_has_callable_signature_coverage("ATAATC---", [(6, 7, 8)])
    assert not PlotHcvMutationTrends._sequence_has_callable_signature_coverage("ATNATCTAC", [(0, 1, 2)])
    assert not PlotHcvMutationTrends._sequence_has_callable_signature_coverage("ATA", [(0, 1, 3)])


def test_finalize_signature_frequency_grid_populates_ci_only_for_callable_rows():
    signature_grid = pd.DataFrame(
        [
            {
                "collection_year": 2020,
                "drug": "DrugA",
                "signature_id": "sig1",
                "signature_kind": "single",
                "protein_name": "NS3",
                "callable_sequences": 2,
                "sequences_with_signature": 1,
            },
            {
                "collection_year": 2021,
                "drug": "DrugA",
                "signature_id": "sig1",
                "signature_kind": "single",
                "protein_name": "NS3",
                "callable_sequences": 0,
                "sequences_with_signature": 0,
            },
        ]
    )

    finalized = PlotHcvMutationTrends._finalize_signature_frequency_grid(signature_grid)
    callable_row = finalized[finalized["collection_year"] == 2020].iloc[0]
    non_callable_row = finalized[finalized["collection_year"] == 2021].iloc[0]

    assert callable_row["proportion"] == 0.5
    assert 0.0 <= callable_row["ci_lower"] <= callable_row["ci_upper"] <= 1.0
    assert pd.isna(non_callable_row["proportion"])
    assert pd.isna(non_callable_row["ci_lower"])
    assert pd.isna(non_callable_row["ci_upper"])


def test_build_yearly_summaries_raises_when_no_valid_collection_years(tmp_path):
    db_path = tmp_path / "missing_years.db"
    _create_base_plot_db(db_path, include_years=False)

    with pytest.raises(ValueError, match="No non-excluded sequences with a valid four-digit collection_year"):
        PlotHcvMutationTrends.build_yearly_summaries(db_path, min_yearly_sequences=1)


def test_build_yearly_summaries_raises_when_threshold_filters_all_years(hcv_plot_db):
    with pytest.raises(ValueError, match="minimum yearly sample size threshold of 3"):
        PlotHcvMutationTrends.build_yearly_summaries(hcv_plot_db, min_yearly_sequences=3)


def test_build_yearly_summaries_returns_empty_signature_frequency_without_alignment(tmp_path):
    db_path = tmp_path / "no_alignment.db"
    _create_base_plot_db(db_path, include_alignment=False)

    results = PlotHcvMutationTrends.build_yearly_summaries(db_path, min_yearly_sequences=2)

    assert results["yearly_mutation_summary"].shape[0] == 2
    assert results["yearly_signature_frequency"].empty
    assert results["yearly_gene_signature_frequency"].empty
    assert results["priority_drugs"].empty


def test_build_yearly_summaries_raises_for_blank_alignment_names_with_multiple_reference_hints(tmp_path):
    db_path = tmp_path / "blank_alignment_name.db"
    _create_base_plot_db(db_path, blank_alignment_name=True)

    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute(
            "INSERT INTO mutation_catalog VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ("mut4", "NS3", "1", "1", "A", "REF2", "snp", "sig_ref2", "single", "", "", "drug_resistance", "I", "DrugC"),
        )
        conn.commit()
    finally:
        conn.close()

    with pytest.raises(ValueError, match="blank alignment_name values"):
        PlotHcvMutationTrends.build_yearly_summaries(db_path, min_yearly_sequences=2)


def test_plot_signature_proportions_handles_empty_and_zero_callable_frames():
    figure, ax = plt.subplots()
    PlotHcvMutationTrends._plot_signature_proportions(ax, pd.DataFrame(columns=PlotHcvMutationTrends.SIGNATURE_PLOT_COLUMNS))
    assert not getattr(ax, "axison")
    plt.close(figure)

    figure, ax = plt.subplots()
    zero_callable = pd.DataFrame(
        [
            {
                "collection_year": 2020,
                "drug": "DrugA",
                "signature_id": "sig1",
                "signature_kind": "single",
                "protein_name": "NS3",
                "callable_sequences": 0,
                "sequences_with_signature": 0,
                "proportion": float("nan"),
                "ci_lower": float("nan"),
                "ci_upper": float("nan"),
                "smoothed_proportion": float("nan"),
            }
        ]
    )
    PlotHcvMutationTrends._plot_signature_proportions(ax, zero_callable)
    assert not getattr(ax, "axison")
    plt.close(figure)


def test_add_callable_sequence_overlay_uses_yearly_gene_denominator():
    figure, ax = plt.subplots()
    gene_df = pd.DataFrame(
        [
            {
                "collection_year": 2020,
                "drug": "DrugA",
                "signature_id": "sig_a",
                "signature_kind": "single",
                "protein_name": "NS3",
                "callable_sequences": 5,
                "sequences_with_signature": 2,
                "proportion": 0.4,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.4,
            },
            {
                "collection_year": 2020,
                "drug": "DrugB",
                "signature_id": "sig_b",
                "signature_kind": "single",
                "protein_name": "NS3",
                "callable_sequences": 5,
                "sequences_with_signature": 1,
                "proportion": 0.2,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.2,
            },
            {
                "collection_year": 2021,
                "drug": "DrugA",
                "signature_id": "sig_a",
                "signature_kind": "single",
                "protein_name": "NS3",
                "callable_sequences": 4,
                "sequences_with_signature": 1,
                "proportion": 0.25,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.25,
            },
        ]
    )

    line = PlotHcvMutationTrends._add_callable_sequence_overlay(ax, gene_df, "Sequences with full-gene coverage")

    assert line is not None
    assert line.get_label() == "Sequences with full-gene coverage"
    assert line.get_xdata().tolist() == [2020, 2021]
    assert line.get_ydata().tolist() == [5, 4]
    assert len(figure.axes) == 2
    assert figure.axes[1].get_ylabel() == "Sequences with full-gene coverage"
    plt.close(figure)


def test_build_target_mutation_frequency_groups_requested_genotypes(monkeypatch):
    cohort = pd.DataFrame(
        [
            {"primary_accession": "A1", "collection_year": 2020},
            {"primary_accession": "A2", "collection_year": 2020},
            {"primary_accession": "A3", "collection_year": 2021},
            {"primary_accession": "A4", "collection_year": 2021},
        ]
    )
    mutation_counts = pd.DataFrame(
        [
            {"primary_accession": "A1", "relevant_mutations_present": "NS5A:93H", "total_relevant_mutation_count": 1},
            {"primary_accession": "A2", "relevant_mutations_present": "", "total_relevant_mutation_count": 0},
            {"primary_accession": "A3", "relevant_mutations_present": "NS5A:93H", "total_relevant_mutation_count": 1},
            {"primary_accession": "A4", "relevant_mutations_present": "NS5B:282T", "total_relevant_mutation_count": 1},
        ]
    )
    target_mutation_catalog = pd.DataFrame(
        [
            {"mutation_id": "NS5A:93H", "protein_name": "NS5A", "segment": "1", "aa_position": "1", "reference_accession": "REF1"},
            {"mutation_id": "NS5B:282T", "protein_name": "NS5B", "segment": "1", "aa_position": "2", "reference_accession": "REF1"},
        ]
    )
    seq_aln = pd.DataFrame([{"primary_accession": "REF1", "alignment": "AAAAAA"}])
    meta_data = pd.DataFrame(
        [
            {"primary_accession": "A1", "nearest_reference_genotype": "1", "nearest_reference_subtype": "a"},
            {"primary_accession": "A2", "nearest_reference_genotype": "3", "nearest_reference_subtype": "a"},
            {"primary_accession": "A3", "nearest_reference_genotype": "3", "nearest_reference_subtype": "a"},
            {"primary_accession": "A4", "nearest_reference_genotype": "2", "nearest_reference_subtype": "b"},
        ]
    )
    target_defs = [
        {
            "mutation_id": "NS5A:93H",
            "plot_label": "NS5A Y93H",
            "lineages": [
                {"label": "Genotype 1", "genotype": "1"},
                {"label": "Genotype 3", "genotype": "3"},
            ],
        },
        {
            "mutation_id": "NS5B:282T",
            "plot_label": "NS5B S282T",
            "lineages": "all_genotypes",
        },
    ]

    def fake_resolve_reference_contexts(*args, **kwargs):
        return {
            "REF1": {
                "group_df": pd.DataFrame(
                    [
                        {"primary_accession": "A1", "alignment": "AAAAAA"},
                        {"primary_accession": "A2", "alignment": "NNNAAA"},
                        {"primary_accession": "A3", "alignment": "AAAAAA"},
                        {"primary_accession": "A4", "alignment": "AAAAAA"},
                    ]
                )
            }
        }

    def fake_build_mutation_position_index(*args, **kwargs):
        return pd.DataFrame(
            [
                {"reference_id": "REF1", "mutation_id": "NS5A:93H", "alignment_indices_list": [(0, 1, 2)]},
                {"reference_id": "REF1", "mutation_id": "NS5B:282T", "alignment_indices_list": [(3, 4, 5)]},
            ]
        )

    monkeypatch.setattr(PlotHcvMutationTrends, "_resolve_reference_contexts", fake_resolve_reference_contexts)
    monkeypatch.setattr(PlotHcvMutationTrends, "_build_mutation_position_index", fake_build_mutation_position_index)

    summary = PlotHcvMutationTrends._build_target_mutation_frequency(
        cohort,
        mutation_counts,
        target_mutation_catalog,
        seq_aln,
        meta_data,
        alias_lookup={},
        db_gff_maps={},
        target_mutation_definitions=target_defs,
    )

    summary = summary.set_index(["mutation_id", "genotype_label", "collection_year"])
    assert summary.loc[("NS5A:93H", "Genotype 1", 2020), "callable_sequences"] == 1
    assert summary.loc[("NS5A:93H", "Genotype 1", 2020), "sequences_with_mutation"] == 1
    assert summary.loc[("NS5A:93H", "Genotype 3", 2020), "callable_sequences"] == 0
    assert summary.loc[("NS5A:93H", "Genotype 3", 2021), "callable_sequences"] == 1
    assert summary.loc[("NS5A:93H", "Genotype 3", 2021), "sequences_with_mutation"] == 1
    assert summary.loc[("NS5B:282T", "Genotype 2", 2021), "callable_sequences"] == 1
    assert summary.loc[("NS5B:282T", "Genotype 2", 2021), "sequences_with_mutation"] == 1


def test_plot_target_mutation_proportions_handles_empty_and_nonempty_frames():
    figure, ax = plt.subplots()
    PlotHcvMutationTrends._plot_target_mutation_proportions(ax, pd.DataFrame(columns=PlotHcvMutationTrends.TARGET_MUTATION_PLOT_COLUMNS), "Target")
    assert not getattr(ax, "axison")
    plt.close(figure)

    figure, ax = plt.subplots()
    target_df = pd.DataFrame(
        [
            {
                "collection_year": 2020,
                "mutation_id": "NS5A:93H",
                "plot_label": "NS5A Y93H",
                "genotype_label": "Genotype 1",
                "callable_sequences": 2,
                "sequences_with_mutation": 1,
                "proportion": 0.5,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.5,
            },
            {
                "collection_year": 2021,
                "mutation_id": "NS5A:93H",
                "plot_label": "NS5A Y93H",
                "genotype_label": "Genotype 3",
                "callable_sequences": 3,
                "sequences_with_mutation": 1,
                "proportion": 1 / 3,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 1 / 3,
            },
        ]
    )
    PlotHcvMutationTrends._plot_target_mutation_proportions(ax, target_df, "NS5A Y93H frequency by year")
    assert getattr(ax, "axison")
    assert ax.get_ylabel() == "Mutation frequency"
    assert ax.get_legend() is not None
    plt.close(figure)


def test_resistance_category_scoring_and_priority_selection_prefers_severe_combination_signatures():
    assert PlotHcvMutationTrends._resistance_category_score("III") > PlotHcvMutationTrends._resistance_category_score("II")
    assert PlotHcvMutationTrends._resistance_category_score("5") > PlotHcvMutationTrends._resistance_category_score("III")

    gene_signature_frequency = pd.DataFrame(
        [
            {
                "collection_year": 2020,
                "drug": "DrugB",
                "signature_id": "sig_b",
                "signature_kind": "single",
                "protein_name": "NS3",
                "resistance_category": "I",
                "combination_size": "",
                "callable_sequences": 5,
                "sequences_with_signature": 2,
                "proportion": 0.4,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.4,
            },
            {
                "collection_year": 2020,
                "drug": "DrugA",
                "signature_id": "sig_a_combo",
                "signature_kind": "combination",
                "protein_name": "NS5A",
                "resistance_category": "III",
                "combination_size": "2",
                "callable_sequences": 4,
                "sequences_with_signature": 1,
                "proportion": 0.25,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.25,
            },
            {
                "collection_year": 2020,
                "drug": "DrugC",
                "signature_id": "sig_c",
                "signature_kind": "single",
                "protein_name": "NS5B",
                "resistance_category": "II",
                "combination_size": "",
                "callable_sequences": 4,
                "sequences_with_signature": 1,
                "proportion": 0.25,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.25,
            },
        ]
    )

    priority_drugs = PlotHcvMutationTrends._select_priority_drugs(gene_signature_frequency, max_drugs=2)
    assert priority_drugs["drug"].tolist() == ["DrugA", "DrugC"]

    priority_subset = PlotHcvMutationTrends._select_priority_signature_subset(gene_signature_frequency, "DrugA")
    assert priority_subset["signature_id"].tolist() == ["sig_a_combo"]

    category_i_subset = PlotHcvMutationTrends._select_category_i_signature_subset(gene_signature_frequency, "DrugB")
    assert category_i_subset["signature_id"].tolist() == ["sig_b"]


def test_select_category_i_signature_subset_filters_to_category_i_rows():
    gene_signature_frequency = pd.DataFrame(
        [
            {
                "collection_year": 2020,
                "drug": "DrugA",
                "signature_id": "sig_a_i",
                "signature_kind": "single",
                "protein_name": "NS3",
                "resistance_category": "I",
                "combination_size": "",
                "callable_sequences": 5,
                "sequences_with_signature": 2,
                "proportion": 0.4,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.4,
            },
            {
                "collection_year": 2021,
                "drug": "DrugA",
                "signature_id": "sig_a_i",
                "signature_kind": "single",
                "protein_name": "NS3",
                "resistance_category": "I",
                "combination_size": "",
                "callable_sequences": 4,
                "sequences_with_signature": 1,
                "proportion": 0.25,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.25,
            },
            {
                "collection_year": 2020,
                "drug": "DrugA",
                "signature_id": "sig_a_ii",
                "signature_kind": "combination",
                "protein_name": "NS5A",
                "resistance_category": "II",
                "combination_size": "2",
                "callable_sequences": 5,
                "sequences_with_signature": 1,
                "proportion": 0.2,
                "ci_lower": 0.0,
                "ci_upper": 1.0,
                "smoothed_proportion": 0.2,
            },
        ]
    )

    subset = PlotHcvMutationTrends._select_category_i_signature_subset(gene_signature_frequency, "DrugA")

    assert subset["signature_id"].unique().tolist() == ["sig_a_i"]
    assert subset["collection_year"].tolist() == [2020, 2021]