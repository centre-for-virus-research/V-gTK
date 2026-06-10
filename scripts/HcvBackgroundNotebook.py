"""Notebook-style background plots for the 12 May HCV database."""

# %% [markdown]
# # HCV Background Notebook
#
# This script is written in a notebook-friendly `.py` format for VS Code and
# Python Interactive. It reuses the cohort and mutation-loading helpers from the
# main HCV trend plotting workflow, then builds a lighter set of background
# figures:
#
# - yearly sequence counts
# - stacked genotype counts by year
# - focus-country sequence counts for the UK, China, and USA
# - separate key-mutation plots by country and by mutation

# %%
from __future__ import annotations

import sqlite3
from collections.abc import Sequence
from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure

from PlotHcvMutationTrends import (  # noqa: E402
    TARGET_MUTATION_PLOT_DEFINITIONS,
    _load_reference_metadata,
    _load_relevant_mutation_counts,
    _load_yearly_cohort,
    _mutation_token_present,
    _slugify_label,
)


DB_PATH = Path("/home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless_update/HCV_full.db")
OUTPUT_DIR = Path("/home3/oml4h/RABV-gTK/test_out/hcv_background_notebook")

#python /home3/oml4h/RABV-gTK/scripts/PlotHcvMutationTrends.py --db /home3/oml4h/RABV-gTK/test_out/HCV_full_XML_treeless_update/HCV_full.db --output_dir /home3/oml4h/RABV-gTK/test_out/hcv_background_notebook/extra
FOCUS_COUNTRIES = ("United Kingdom", "China", "USA")
EXCLUDED_COLLECTION_YEARS = {1905}
AXIS_LABEL_SIZE = 18
TICK_LABEL_SIZE = 15
TITLE_SIZE = 20
LEGEND_SIZE = 14
MAX_YEAR_TICKS = 8

KEY_MUTATION_DEFINITIONS = [
    {
        "mutation_id": definition["mutation_id"],
        "plot_label": definition["plot_label"],
    }
    for definition in TARGET_MUTATION_PLOT_DEFINITIONS
]


def load_background_dataset(db_path: Path) -> pd.DataFrame:
    conn = sqlite3.connect(str(db_path))
    try:
        cohort = _load_yearly_cohort(conn)
        meta_data = _load_reference_metadata(conn)
        mutation_counts = _load_relevant_mutation_counts(conn)
    finally:
        conn.close()

    metadata_columns = [
        column
        for column in [
            "primary_accession",
            "country",
            "country_validated",
            "nearest_reference_genotype",
            "nearest_reference_subtype",
            "host_scientific_name",
        ]
        if column in meta_data.columns
    ]
    metadata = meta_data[metadata_columns].drop_duplicates(subset=["primary_accession"])

    background_df = cohort.merge(metadata, on="primary_accession", how="left")
    background_df = background_df.merge(
        mutation_counts[["primary_accession", "relevant_mutations_present"]].drop_duplicates(subset=["primary_accession"]),
        on="primary_accession",
        how="left",
    )

    background_df["collection_year"] = pd.to_numeric(background_df["collection_year"], errors="coerce")
    background_df = background_df.dropna(subset=["collection_year"]).copy()
    background_df["collection_year"] = background_df["collection_year"].astype(int)
    background_df["genotype"] = background_df.get("nearest_reference_genotype", "").fillna("").astype(str).str.strip()
    background_df["subtype"] = background_df.get("nearest_reference_subtype", "").fillna("").astype(str).str.strip().str.lower()
    background_df["country_display"] = background_df.get("country", "").fillna("").astype(str).str.strip()
    if "country_validated" in background_df.columns:
        validated = background_df["country_validated"].fillna("").astype(str).str.strip()
        background_df.loc[background_df["country_display"] == "", "country_display"] = validated
    background_df["relevant_mutations_present"] = background_df["relevant_mutations_present"].fillna("").astype(str)
    return background_df


def normalize_focus_country(country_value: str) -> str:
    text = str(country_value or "").strip()
    if not text:
        return ""
    lowered = text.lower()
    if lowered.startswith("united kingdom") or lowered == "uk":
        return "United Kingdom"
    if lowered.startswith("china"):
        return "China"
    if lowered.startswith("usa") or lowered.startswith("united states"):
        return "USA"
    return ""


def build_yearly_sequence_counts(background_df: pd.DataFrame) -> pd.DataFrame:
    yearly_counts = (
        background_df.groupby("collection_year", as_index=False)
        .agg(sequence_count=("primary_accession", "nunique"))
        .sort_values("collection_year")
    )
    years = build_complete_year_index(background_df["collection_year"])
    if not years:
        return yearly_counts.iloc[0:0].copy()
    return pd.DataFrame({"collection_year": years}).merge(yearly_counts, on="collection_year", how="left").fillna({"sequence_count": 0})


def build_yearly_genotype_counts(background_df: pd.DataFrame) -> pd.DataFrame:
    genotype_df = background_df.copy()
    genotype_df["genotype_label"] = genotype_df["genotype"].replace("", "Unknown")
    genotype_counts = (
        genotype_df.groupby(["collection_year", "genotype_label"], as_index=False)
        .agg(sequence_count=("primary_accession", "nunique"))
        .sort_values(["collection_year", "genotype_label"])
    )
    years = build_complete_year_index(genotype_df["collection_year"])
    genotype_labels = sorted(genotype_counts["genotype_label"].dropna().astype(str).unique().tolist())
    if not years or not genotype_labels:
        return genotype_counts.iloc[0:0].copy()
    year_label_grid = pd.MultiIndex.from_product(
        [years, genotype_labels],
        names=["collection_year", "genotype_label"],
    ).to_frame(index=False)
    return year_label_grid.merge(genotype_counts, on=["collection_year", "genotype_label"], how="left").fillna({"sequence_count": 0})


def build_focus_country_dataset(background_df: pd.DataFrame) -> pd.DataFrame:
    focus_df = background_df.copy()
    focus_df["focus_country"] = focus_df["country_display"].map(normalize_focus_country)
    focus_df = focus_df[focus_df["focus_country"].isin(FOCUS_COUNTRIES)].copy()
    return focus_df


def build_focus_country_sequence_counts(focus_df: pd.DataFrame) -> pd.DataFrame:
    country_counts = (
        focus_df.groupby(["focus_country", "collection_year"], as_index=False)
        .agg(sequence_count=("primary_accession", "nunique"))
        .sort_values(["focus_country", "collection_year"])
    )
    years = build_complete_year_index(focus_df["collection_year"])
    if not years:
        return country_counts.iloc[0:0].copy()
    year_country_grid = pd.MultiIndex.from_product(
        [FOCUS_COUNTRIES, years],
        names=["focus_country", "collection_year"],
    ).to_frame(index=False)
    return year_country_grid.merge(country_counts, on=["focus_country", "collection_year"], how="left").fillna({"sequence_count": 0})


def build_focus_country_mutation_frequency(
    focus_df: pd.DataFrame,
    mutation_definitions: list[dict[str, str]] | None = None,
) -> pd.DataFrame:
    mutation_definitions = mutation_definitions or KEY_MUTATION_DEFINITIONS
    years = build_complete_year_index(focus_df["collection_year"])
    if not years:
        return pd.DataFrame(
            columns=[
                "focus_country",
                "collection_year",
                "mutation_id",
                "plot_label",
                "sequence_count",
                "mutation_count",
                "proportion",
            ]
        )

    records = []
    for country in FOCUS_COUNTRIES:
        country_df = focus_df[focus_df["focus_country"] == country].copy()
        yearly_counts = (
            country_df.groupby("collection_year", as_index=False)
            .agg(sequence_count=("primary_accession", "nunique"))
            .sort_values("collection_year")
        )

        for mutation_definition in mutation_definitions:
            mutation_id = mutation_definition["mutation_id"]
            plot_label = mutation_definition["plot_label"]
            mutation_mask = _mutation_token_present(country_df["relevant_mutations_present"], mutation_id)
            mutation_counts = (
                country_df.loc[mutation_mask]
                .groupby("collection_year", as_index=False)
                .agg(mutation_count=("primary_accession", "nunique"))
                .sort_values("collection_year")
            )

            year_grid = pd.DataFrame({"collection_year": years})
            year_grid["focus_country"] = country
            year_grid["mutation_id"] = mutation_id
            year_grid["plot_label"] = plot_label
            year_grid = year_grid.merge(yearly_counts, on="collection_year", how="left")
            year_grid = year_grid.merge(mutation_counts, on="collection_year", how="left")
            year_grid["sequence_count"] = year_grid["sequence_count"].fillna(0)
            year_grid["mutation_count"] = year_grid["mutation_count"].fillna(0)
            year_grid["proportion"] = year_grid["mutation_count"] / year_grid["sequence_count"].where(
                year_grid["sequence_count"] > 0
            )
            records.append(year_grid)

    return pd.concat(records, ignore_index=True)


def build_complete_year_index(year_values: pd.Series) -> list[int]:
    years = sorted({int(year) for year in pd.Series(year_values).dropna().astype(int).tolist() if int(year) not in EXCLUDED_COLLECTION_YEARS})
    if not years:
        return []
    return [year for year in range(years[0], years[-1] + 1) if year not in EXCLUDED_COLLECTION_YEARS]


def select_sparse_year_ticks(years: list[int], max_ticks: int = MAX_YEAR_TICKS) -> list[int]:
    if not years:
        return []
    if len(years) <= max_ticks:
        return years
    step = max(1, ceil(len(years) / max_ticks))
    tick_years = years[::step]
    if tick_years[-1] != years[-1]:
        tick_years.append(years[-1])
    return tick_years


def style_axis_text(ax: plt.Axes, *, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, fontsize=TITLE_SIZE)
    ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)


def style_year_axis(ax: plt.Axes, years: Sequence[int | float], *, categorical: bool = False) -> None:
    normalized_years = [int(year) for year in years]
    tick_years = select_sparse_year_ticks(normalized_years)
    if categorical:
        tick_positions = [normalized_years.index(year) for year in tick_years]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([str(year) for year in tick_years], rotation=0, ha="center")
    else:
        ax.set_xticks(tick_years)
        ax.set_xlim(normalized_years[0] - 0.5, normalized_years[-1] + 0.5)
    ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE)
    ax.tick_params(axis="y", labelsize=TICK_LABEL_SIZE)


def build_yearly_genotype_pivot(genotype_counts: pd.DataFrame) -> pd.DataFrame:
    return (
        genotype_counts.pivot(index="collection_year", columns="genotype_label", values="sequence_count")
        .fillna(0)
        .sort_index()
    )


def filter_country_mutation_summary(country_mutation_df: pd.DataFrame, *, country: str = "", plot_label: str = "") -> pd.DataFrame:
    filtered_df = country_mutation_df.copy()
    if country:
        filtered_df = filtered_df[filtered_df["focus_country"] == country]
    if plot_label:
        filtered_df = filtered_df[filtered_df["plot_label"] == plot_label]
    return filtered_df.sort_values(["focus_country", "plot_label", "collection_year"]).reset_index(drop=True)


def determine_frequency_ylim(proportions: pd.Series) -> tuple[float, float]:
    values = pd.to_numeric(proportions, errors="coerce").dropna()
    positive_values = values[values > 0]
    if positive_values.empty:
        return (0.0, 1.0)

    max_value = float(positive_values.max())
    upper = min(1.0, max_value * 1.2)
    if upper <= max_value:
        upper = min(1.0, max_value + 0.02)
    return (0.0, upper)


def apply_frequency_axis(ax: plt.Axes, proportions: pd.Series, *, log_scale: bool = False) -> None:
    values = pd.to_numeric(proportions, errors="coerce").dropna()
    positive_values = values[values > 0]
    if log_scale and not positive_values.empty:
        min_positive = float(positive_values.min())
        max_value = float(positive_values.max())
        lower = max(1e-4, min_positive / 1.8)
        upper = min(1.0, max_value * 1.35)
        if upper <= lower:
            upper = min(1.0, max_value * 2)
        ax.set_yscale("log")
        ax.set_ylim(lower, upper)
        return

    ax.set_ylim(*determine_frequency_ylim(proportions))


def save_figure(fig: Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=220, bbox_inches="tight")


def create_sequences_per_year_figure(yearly_counts: pd.DataFrame) -> Figure:
    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.bar(yearly_counts["collection_year"], yearly_counts["sequence_count"], color="#1f4e79", alpha=0.85)
    ax.plot(yearly_counts["collection_year"], yearly_counts["sequence_count"], color="#0b2239", linewidth=1.8)
    years = yearly_counts["collection_year"].astype(int).tolist()
    style_axis_text(ax, title="HCV sequences per year", xlabel="Collection year", ylabel="Sequences")
    style_year_axis(ax, years)
    ax.grid(axis="y", alpha=0.25)
    return fig


def create_genotype_stacked_bar_figure(genotype_counts: pd.DataFrame) -> Figure:
    pivot = build_yearly_genotype_pivot(genotype_counts)
    fig, ax = plt.subplots(figsize=(14, 7))
    pivot.plot(kind="bar", stacked=True, ax=ax, width=0.9, colormap="tab20")
    years = pivot.index.astype(int).tolist()
    style_axis_text(ax, title="Genotype composition by year", xlabel="Collection year", ylabel="Sequences")
    style_year_axis(ax, years, categorical=True)
    ax.legend(title="Genotype", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=LEGEND_SIZE, title_fontsize=LEGEND_SIZE)
    return fig


def create_focus_country_sequences_figure(country_counts: pd.DataFrame) -> Figure:
    fig, ax = plt.subplots(figsize=(12, 6))
    palette = {"United Kingdom": "#7a0019", "China": "#005f73", "USA": "#3a5a40"}
    years = sorted(country_counts["collection_year"].dropna().astype(int).unique().tolist())
    for country, country_df in country_counts.groupby("focus_country", sort=False):
        country_df = country_df.sort_values("collection_year")
        ax.plot(
            country_df["collection_year"],
            country_df["sequence_count"],
            linewidth=2.4,
            marker="o",
            label=country,
            color=palette.get(str(country)),
        )
    style_axis_text(ax, title="Focus-country sequence counts by year", xlabel="Collection year", ylabel="Sequences")
    style_year_axis(ax, years)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=LEGEND_SIZE)
    return fig


def create_country_mutation_frequency_figure(country_mutation_df: pd.DataFrame, country: str, *, log_scale: bool = False) -> Figure:
    country_df = filter_country_mutation_summary(country_mutation_df, country=country)
    fig, ax = plt.subplots(figsize=(13, 6.5))
    years = sorted(country_df["collection_year"].dropna().astype(int).unique().tolist())
    for plot_label, mutation_df in country_df.groupby("plot_label", sort=False):
        mutation_df = mutation_df.sort_values("collection_year")
        y_values = mutation_df["proportion"].where(mutation_df["proportion"] > 0) if log_scale else mutation_df["proportion"]
        ax.plot(
            mutation_df["collection_year"],
            y_values,
            linewidth=2.3,
            marker="o",
            label=plot_label,
        )
    ylabel = "Mutation frequency (log scale)" if log_scale else "Mutation frequency"
    style_axis_text(ax, title=f"Key mutation frequencies by year in {country}", xlabel="Collection year", ylabel=ylabel)
    style_year_axis(ax, years)
    apply_frequency_axis(ax, country_df["proportion"], log_scale=log_scale)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=LEGEND_SIZE)
    return fig


def create_mutation_country_frequency_figure(country_mutation_df: pd.DataFrame, plot_label: str, *, log_scale: bool = False) -> Figure:
    palette = {"United Kingdom": "#7a0019", "China": "#005f73", "USA": "#3a5a40"}
    mutation_df = filter_country_mutation_summary(country_mutation_df, plot_label=plot_label)
    fig, ax = plt.subplots(figsize=(13, 6.5))
    years = sorted(mutation_df["collection_year"].dropna().astype(int).unique().tolist())
    for country, country_df in mutation_df.groupby("focus_country", sort=False):
        country_df = country_df.sort_values("collection_year")
        y_values = country_df["proportion"].where(country_df["proportion"] > 0) if log_scale else country_df["proportion"]
        ax.plot(
            country_df["collection_year"],
            y_values,
            linewidth=2.3,
            marker="o",
            label=country,
            color=palette.get(str(country)),
        )
    ylabel = "Mutation frequency (log scale)" if log_scale else "Mutation frequency"
    style_axis_text(ax, title=f"{plot_label} frequency by year in focus countries", xlabel="Collection year", ylabel=ylabel)
    style_year_axis(ax, years)
    apply_frequency_axis(ax, mutation_df["proportion"], log_scale=log_scale)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=LEGEND_SIZE)
    return fig


def create_country_mutation_count_figure(country_mutation_df: pd.DataFrame, country: str) -> Figure:
    country_df = filter_country_mutation_summary(country_mutation_df, country=country)
    fig, ax = plt.subplots(figsize=(13, 6.5))
    years = sorted(country_df["collection_year"].dropna().astype(int).unique().tolist())
    for plot_label, mutation_df in country_df.groupby("plot_label", sort=False):
        mutation_df = mutation_df.sort_values("collection_year")
        ax.plot(
            mutation_df["collection_year"],
            mutation_df["mutation_count"],
            linewidth=2.3,
            marker="o",
            label=plot_label,
        )
    style_axis_text(ax, title=f"Key mutation counts by year in {country}", xlabel="Collection year", ylabel="Sequences with mutation")
    style_year_axis(ax, years)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=LEGEND_SIZE)
    return fig

# %% [markdown]
# ## Notebook Setup
#
# These values are easy to tweak while exploring. Change `db_path` or
# `output_dir` here without touching the helper functions above.

# %%
db_path = DB_PATH
output_dir = OUTPUT_DIR
output_dir.mkdir(parents=True, exist_ok=True)


# %% [markdown]
# ## Load The Background Dataset
#
# This gives one row per sequence with year, genotype, country, and the raw
# mutation string from the SQLite database.

# %%
background_df = load_background_dataset(db_path)
background_df.head()


# %%
background_df[["collection_year", "genotype", "subtype", "country_display"]].describe(include="all")


# %% [markdown]
# ## Quick Metadata Checks

# %%
background_df["collection_year"].agg(["min", "max", "nunique"])


# %%
background_df["country_display"].value_counts().head(15)


# %%
background_df["genotype"].replace("", "Unknown").value_counts().sort_index()


# %% [markdown]
# ## Sequences Per Year

# %%
yearly_counts = build_yearly_sequence_counts(background_df)
yearly_counts.tail(15)


# %%
fig = create_sequences_per_year_figure(yearly_counts)
plt.show()


# %% [markdown]
# ## Genotypes Per Year

# %%
genotype_counts = build_yearly_genotype_counts(background_df)
genotype_counts.head(20)


# %%
genotype_pivot = build_yearly_genotype_pivot(genotype_counts)
genotype_pivot.tail(15)


# %%
fig = create_genotype_stacked_bar_figure(genotype_counts)
plt.show()


# %% [markdown]
# ## Focus Countries: UK, China, USA

# %%
focus_df = build_focus_country_dataset(background_df)
focus_df.head()


# %%
focus_df["focus_country"].value_counts()


# %%
focus_country_counts = build_focus_country_sequence_counts(focus_df)
focus_country_counts.tail(20)


# %%
fig = create_focus_country_sequences_figure(focus_country_counts)
plt.show()


# %% [markdown]
# ## Key Mutations In Focus Countries

# %%
country_mutation_df = build_focus_country_mutation_frequency(focus_df)
country_mutation_df.head(20)


# %%
filter_country_mutation_summary(country_mutation_df, country="United Kingdom").head(20)


# %%
filter_country_mutation_summary(country_mutation_df, plot_label="NS5A Y93H").head(20)


# %% [markdown]
# ## Country-Specific Mutation Frequency Plots

# %%
fig = create_country_mutation_frequency_figure(country_mutation_df, "United Kingdom")
plt.show()


# %%
fig = create_country_mutation_frequency_figure(country_mutation_df, "China")
plt.show()


# %%
fig = create_country_mutation_frequency_figure(country_mutation_df, "USA")
plt.show()


# %% [markdown]
# ## Mutation-Specific Plots Across Countries

# %%
for mutation_definition in KEY_MUTATION_DEFINITIONS:
    fig = create_mutation_country_frequency_figure(country_mutation_df, mutation_definition["plot_label"])
    plt.show()


# %% [markdown]
# ## Country-Specific Mutation Counts

# %%
fig = create_country_mutation_count_figure(country_mutation_df, "United Kingdom")
plt.show()


# %%
fig = create_country_mutation_count_figure(country_mutation_df, "China")
plt.show()


# %%
fig = create_country_mutation_count_figure(country_mutation_df, "USA")
plt.show()


# %% [markdown]
# ## Optional Save Block
#
# Run this only when you want to write the current summaries and plots to disk.

# %%
yearly_counts.to_csv(output_dir / "yearly_sequence_count.tsv", sep="\t", index=False)
genotype_counts.to_csv(output_dir / "yearly_genotype_count.tsv", sep="\t", index=False)
focus_country_counts.to_csv(output_dir / "focus_country_sequence_count.tsv", sep="\t", index=False)
country_mutation_df.to_csv(output_dir / "focus_country_key_mutation_frequency.tsv", sep="\t", index=False)

save_figure(create_sequences_per_year_figure(yearly_counts), output_dir / "yearly_sequence_count.png")
plt.close("all")
save_figure(create_genotype_stacked_bar_figure(genotype_counts), output_dir / "genotype_stacked_bar_by_year.png")
plt.close("all")
save_figure(create_focus_country_sequences_figure(focus_country_counts), output_dir / "focus_country_sequence_count_by_year.png")
plt.close("all")

for country in FOCUS_COUNTRIES:
    save_figure(
        create_country_mutation_frequency_figure(country_mutation_df, country),
        output_dir / f"{_slugify_label(country)}_key_mutation_frequency_by_year.png",
    )
    plt.close("all")
    save_figure(
        create_country_mutation_count_figure(country_mutation_df, country),
        output_dir / f"{_slugify_label(country)}_key_mutation_count_by_year.png",
    )
    plt.close("all")

for mutation_definition in KEY_MUTATION_DEFINITIONS:
    save_figure(
        create_mutation_country_frequency_figure(country_mutation_df, mutation_definition["plot_label"]),
        output_dir / f"{_slugify_label(mutation_definition['plot_label'])}_focus_country_frequency_by_year.png",
    )
    plt.close("all")
    save_figure(
        create_mutation_country_frequency_figure(
            country_mutation_df,
            mutation_definition["plot_label"],
            log_scale=True,
        ),
        output_dir / f"{_slugify_label(mutation_definition['plot_label'])}_focus_country_frequency_by_year_log.png",
    )
    plt.close("all")
