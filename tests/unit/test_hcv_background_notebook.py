import importlib
import os
import sys

import pandas as pd


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../scripts")))

HcvBackgroundNotebook = importlib.import_module("HcvBackgroundNotebook")


def test_normalize_focus_country_maps_expected_variants():
    assert HcvBackgroundNotebook.normalize_focus_country("United Kingdom") == "United Kingdom"
    assert HcvBackgroundNotebook.normalize_focus_country("United Kingdom: Scotland") == "United Kingdom"
    assert HcvBackgroundNotebook.normalize_focus_country("China") == "China"
    assert HcvBackgroundNotebook.normalize_focus_country("USA") == "USA"
    assert HcvBackgroundNotebook.normalize_focus_country("United States of America") == "USA"
    assert HcvBackgroundNotebook.normalize_focus_country("France") == ""


def test_build_focus_country_mutation_frequency_counts_hits_and_denominators():
    focus_df = pd.DataFrame(
        [
            {
                "primary_accession": "A1",
                "collection_year": 2020,
                "focus_country": "United Kingdom",
                "relevant_mutations_present": "NS5A:93H;NS5B:282T",
            },
            {
                "primary_accession": "A2",
                "collection_year": 2020,
                "focus_country": "United Kingdom",
                "relevant_mutations_present": "",
            },
            {
                "primary_accession": "A3",
                "collection_year": 2021,
                "focus_country": "China",
                "relevant_mutations_present": "NS5A:30K",
            },
        ]
    )

    summary = HcvBackgroundNotebook.build_focus_country_mutation_frequency(
        focus_df,
        mutation_definitions=[
            {"mutation_id": "NS5A:93H", "plot_label": "NS5A Y93H"},
            {"mutation_id": "NS5A:30K", "plot_label": "NS5A A30K"},
        ],
    )

    summary = summary.set_index(["focus_country", "mutation_id", "collection_year"])
    assert summary.loc[("United Kingdom", "NS5A:93H", 2020), "sequence_count"] == 2
    assert summary.loc[("United Kingdom", "NS5A:93H", 2020), "mutation_count"] == 1
    assert summary.loc[("United Kingdom", "NS5A:93H", 2020), "proportion"] == 0.5
    assert summary.loc[("China", "NS5A:30K", 2021), "sequence_count"] == 1
    assert summary.loc[("China", "NS5A:30K", 2021), "mutation_count"] == 1
    assert summary.loc[("China", "NS5A:30K", 2021), "proportion"] == 1.0


def test_build_yearly_genotype_pivot_returns_wide_year_by_genotype_table():
    genotype_counts = pd.DataFrame(
        [
            {"collection_year": 2020, "genotype_label": "1", "sequence_count": 3},
            {"collection_year": 2020, "genotype_label": "3", "sequence_count": 2},
            {"collection_year": 2021, "genotype_label": "1", "sequence_count": 4},
        ]
    )

    pivot = HcvBackgroundNotebook.build_yearly_genotype_pivot(genotype_counts)

    assert pivot.loc[2020, "1"] == 3
    assert pivot.loc[2020, "3"] == 2
    assert pivot.loc[2021, "1"] == 4
    assert pivot.loc[2021, "3"] == 0


def test_filter_country_mutation_summary_filters_country_and_label():
    country_mutation_df = pd.DataFrame(
        [
            {"focus_country": "United Kingdom", "plot_label": "NS5A Y93H", "collection_year": 2020},
            {"focus_country": "China", "plot_label": "NS5A Y93H", "collection_year": 2020},
            {"focus_country": "United Kingdom", "plot_label": "NS5A A30K", "collection_year": 2021},
        ]
    )

    filtered = HcvBackgroundNotebook.filter_country_mutation_summary(
        country_mutation_df,
        country="United Kingdom",
        plot_label="NS5A Y93H",
    )

    assert filtered.shape[0] == 1
    assert filtered.iloc[0]["focus_country"] == "United Kingdom"
    assert filtered.iloc[0]["plot_label"] == "NS5A Y93H"


def test_build_yearly_sequence_counts_fills_missing_years_and_excludes_1905():
    background_df = pd.DataFrame(
        [
            {"primary_accession": "A1", "collection_year": 1905},
            {"primary_accession": "A2", "collection_year": 2020},
            {"primary_accession": "A3", "collection_year": 2022},
        ]
    )

    yearly_counts = HcvBackgroundNotebook.build_yearly_sequence_counts(background_df)

    assert yearly_counts["collection_year"].tolist() == [2020, 2021, 2022]
    assert yearly_counts["sequence_count"].tolist() == [1, 0, 1]


def test_determine_frequency_ylim_uses_data_max_instead_of_defaulting_to_one():
    lower, upper = HcvBackgroundNotebook.determine_frequency_ylim(pd.Series([0.01, 0.08, 0.10]))

    assert lower == 0.0
    assert upper == 0.12