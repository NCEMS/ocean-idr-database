# python-scripts/idr_analysis_setup.py

import os
import sqlite3
from itertools import combinations
from typing import List

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display, clear_output
import ipywidgets as widgets
from scipy.stats import mannwhitneyu

# --------------------
# Config
# --------------------

DB_PATH = "minimal_noenv.db"

DATASET_OPTIONS = [
    "omrgcv2",
    "3300003177",
    "3300003178",
    "3300003872",
    "3300003873",
]

DEFAULT_DATASETS = list(DATASET_OPTIONS)

MASTER_IDR_PROPS = [
    'Asphericity','Radius_of_gyration','Radius_of_gyration_scaled',
    'End_to_end_distance','End_to_end_distance_scaled',
    'Scaling_exponent','Prefactor','Kappa','Length','FCR',
    'NCPR','SHD','SCD','Molecular_weight','F_Neg','F_Pos',
    'Hydrophobicity','Fraction_aromatic','Fraction_aliphatic',
    'Fraction_polar','Complexity'
]

DEFAULT_SELECTED = {
    'FCR', 'NCPR', 'SHD', 'SCD',
    'F_Neg', 'F_Pos', 'Hydrophobicity',
}

LAST_RESULTS = None

# --------------------
# Core helpers
# --------------------

def run_query(query: str, params: tuple = ()) -> pd.DataFrame:
    with sqlite3.connect(DB_PATH) as conn:
        return pd.read_sql_query(query, conn, params=params)

def query_by_cog_root(cog_root_value: str) -> pd.DataFrame:
    query = """
    SELECT
        a.*,
        i.*,
        c.cog_root_ids
    FROM architecture AS a
    JOIN idr_props AS i
        ON a.Gene = i.Gene
    JOIN cog_root AS c
        ON a.Gene = c.Gene
    WHERE c.cog_root_ids = ?
    """
    return run_query(query, (cog_root_value,))

def query_omrgcv2_with_temps() -> pd.DataFrame:
    """
    Join architecture, idr_props, cog_root, and tara_temperatures
    for entries from the omrgcv2 dataset only.
    """
    query = """
    SELECT
        a.*,
        i.*,
        c.cog_root_ids,
        t.*
    FROM architecture AS a
    JOIN idr_props AS i
        ON a.Gene = i.Gene
    JOIN cog_root AS c
        ON a.Gene = c.Gene
    JOIN tara_temperatures AS t
        ON a.Gene = t.Gene
    WHERE a.dataset = 'omrgcv2'
    """
    return run_query(query)


def check_idr_props_list(user_props: List[str]) -> None:
    invalid = [p for p in user_props if p not in MASTER_IDR_PROPS]
    if invalid:
        allowed_str = ", ".join(MASTER_IDR_PROPS)
        invalid_str = ", ".join(invalid)
        raise ValueError(
            f"Invalid properties: {invalid_str}. "
            f"Allowed properties are: {allowed_str}."
        )

def plot_idr_property_boxplots2(df: pd.DataFrame,
                                properties: list[str],
                                group_col: str = "dataset",
                                ncols: int = 4,
                                figsize: tuple = (16, 10)):
    if not properties:
        raise ValueError("No properties provided to plot.")

    nplots = len(properties)
    nrows = (nplots + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    if group_col not in df.columns:
        raise ValueError(f"Expected column '{group_col}' not found in results.")

    datasets = df[group_col].dropna().unique()
    pairs = list(combinations(datasets, 2)) if len(datasets) >= 2 else []

    for i, prop in enumerate(properties):
        ax = axes[i]

        if prop not in df.columns:
            ax.set_visible(False)
            continue

        sns.boxplot(data=df, x=group_col, y=prop, ax=ax, fliersize=1)
        ax.set_xlabel("")
        ax.set_ylabel(prop)
        ax.tick_params(axis="x", rotation=45)

        pvals = {}
        significant = False
        for a, b in pairs:
            group_a = df.loc[df[group_col] == a, prop].dropna()
            group_b = df.loc[df[group_col] == b, prop].dropna()
            if len(group_a) > 0 and len(group_b) > 0:
                _, p = mannwhitneyu(group_a, group_b, alternative="two-sided")
                p_formatted = "<0.001" if p < 0.001 else f"{p:.3f}"
                pvals[f"{a} vs {b}"] = p_formatted
                if p < 0.05:
                    significant = True

        title = f"{{{prop}}}" if significant else prop
        ax.set_title(title)

        if 0 < len(pvals) <= 3:
            text = "\n".join([f"{k}: p={v}" for k, v in pvals.items()])
            ax.text(0.05, 0.95, text, transform=ax.transAxes,
                    fontsize=8, va="top", ha="left",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.6))
        elif len(pvals) > 3:
            print(f"\n{prop}:")
            for k, v in pvals.items():
                print(f"  {k}: p={v}")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.tight_layout()
    plt.show()

# --------------------
# UI factory
# --------------------

def create_idr_ui():
    """Build and return the full interactive UI as a widget container."""

    instructions = widgets.HTML(
        "<b>Instructions:</b> "
        "1) Enter a COG ID. "
        "2) Select dataset(s). "
        "3) Select IDR properties. "
        "4) Click <i>'Run query and plot'</i>."
    )

    # ---- COG ID input ----
    cog_input = widgets.Text(
        description="COG ID:",
        placeholder="e.g. COG0513",
        layout=widgets.Layout(width="250px")
    )

    # ---- Dataset selection ----
    dataset_checkboxes = [
        widgets.Checkbox(
            value=(ds in DEFAULT_DATASETS),
            description=ds,
            indent=False
        )
        for ds in DATASET_OPTIONS
    ]
    dataset_box = widgets.VBox(dataset_checkboxes)

    # ---- IDR property selection ----
    prop_checkboxes = [
        widgets.Checkbox(
            value=(prop in DEFAULT_SELECTED),
            description=prop,
            indent=False
        )
        for prop in MASTER_IDR_PROPS
    ]
    props_box = widgets.VBox(prop_checkboxes)

    # ---- Run button + output ----
    run_button = widgets.Button(
        description="▶ Run query and plot",
        button_style="primary",
        layout=widgets.Layout(width="250px")
    )

    output = widgets.Output()

    # ---- Callback ----
    def on_run_clicked(b):
        global LAST_RESULTS
        with output:
            clear_output()

            # 1. Get COG ID
            cog_id = cog_input.value.strip()
            if not cog_id:
                print("Please enter a COG ID.")
                return

            # 2. Get selected datasets
            selected_datasets = [
                cb.description for cb in dataset_checkboxes if cb.value
            ]
            if not selected_datasets:
                print("Please select at least one dataset.")
                return

            # 3. Get selected IDR properties
            selected_props = [
                cb.description for cb in prop_checkboxes if cb.value
            ]
            if not selected_props:
                print("Please select at least one IDR property.")
                return

            # 4. Validate properties
            try:
                check_idr_props_list(selected_props)
            except ValueError as e:
                print(e)
                return

            # 5. Run query
            results = query_by_cog_root(cog_id)

            if results is None or results.empty:
                print(f"No results found for COG ID: {cog_id}")
                return

            results = results.loc[:, ~results.columns.duplicated()].copy()

            # 6. Filter by selected datasets (if 'dataset' column exists)
            if "dataset" not in results.columns:
                print("Warning: 'dataset' column not found; cannot filter by dataset.")
            else:
                results = results[results["dataset"].isin(selected_datasets)].copy()
                if results.empty:
                    print(
                        "No results found for COG ID "
                        f"{cog_id} in selected dataset(s): "
                        f"{', '.join(selected_datasets)}"
                    )
                    return

            # 7. Coerce selected properties to numeric
            for col in selected_props:
                if col in results.columns:
                    results[col] = pd.to_numeric(results[col], errors='coerce')

            # 8. Summaries
            print(f"{len(results)} unique IDRs retrieved after filtering.\n")

            if "dataset" in results.columns:
                print("Number of IDRs per dataset:")
                print(results["dataset"].value_counts(), "\n")

            # Preview
            display(results.head())

            # 9. Plots
            try:
                plot_idr_property_boxplots2(results, properties=selected_props)
            except Exception as e:
                print("\nPlotting error:", e)

            # 10. Store results for export
            LAST_RESULTS = results

    run_button.on_click(on_run_clicked)

    # ---- Layout ----
    ui = widgets.VBox([
        instructions,
        cog_input,
        widgets.HTML("<b>Select dataset(s):</b>"),
        dataset_box,
        widgets.HTML("<b>Select IDR properties:</b>"),
        props_box,
        run_button,
        output
    ])

    return ui

def create_omrgcv2_temp_ui():
    """
    UI for:
      - omrgcv2 only
      - select IDR properties
      - select min/max temperature thresholds
      - compare proteins uniquely below min vs uniquely above max
    """
    instructions = widgets.HTML(
        "<b>Instructions:</b> "
        "1) Select IDR properties.<br>"
        "2) Enter a minimum and maximum temperature (°C).<br>"
        "3) Click <i>'Run temperature-based query'</i> to compare IDRs from cold vs warm bins."
    )

    # ---- IDR property selection ----
    prop_checkboxes = [
        widgets.Checkbox(
            value=(prop in DEFAULT_SELECTED),
            description=prop,
            indent=False
        )
        for prop in MASTER_IDR_PROPS
    ]
    props_box = widgets.VBox(prop_checkboxes)

    # ---- Temperature inputs ----
    min_temp_input = widgets.FloatText(
        description="Min T (°C):",
        value=5.0,
        layout=widgets.Layout(width="220px")
    )

    max_temp_input = widgets.FloatText(
        description="Max T (°C):",
        value=25.0,
        layout=widgets.Layout(width="220px")
    )

    # ---- Run button + output ----
    run_button = widgets.Button(
        description="▶ Run temperature-based query",
        button_style="primary",
        layout=widgets.Layout(width="300px")
    )

    output = widgets.Output()

    def on_run_clicked(b):
        global LAST_RESULTS
        with output:
            clear_output()

            # 1) Read selected properties
            selected_props = [
                cb.description for cb in prop_checkboxes if cb.value
            ]
            if not selected_props:
                print("Please select at least one IDR property.")
                return

            # 2) Validate properties
            try:
                check_idr_props_list(selected_props)
            except ValueError as e:
                print(e)
                return

            # 3) Get temperature thresholds
            try:
                min_T = float(min_temp_input.value)
                max_T = float(max_temp_input.value)
            except Exception:
                print("Please enter numeric values for temperatures.")
                return

            if min_T >= max_T:
                print("Minimum temperature must be less than maximum temperature.")
                return

            # 4) Run base query (omrgcv2 only, already joined)
            df = query_omrgcv2_with_temps()

            if df is None or df.empty:
                print("No records found for omrgcv2 with temperatures.")
                return

            # 5) Clean duplicated columns
            df = df.loc[:, ~df.columns.duplicated()].copy()

            # 6) Ensure required temperature columns exist
            # Adjust names here if your tara_temperatures schema differs
            temp_cols = ["min_temperature", "max_temperature"]
            for col in temp_cols:
                if col not in df.columns:
                    print(f"Expected temperature column '{col}' not found in joined table.")
                    return

            # 7) Coerce temps to numeric
            for col in temp_cols:
                df[col] = pd.to_numeric(df[col], errors="coerce")

            # Drop rows with missing temps
            df = df.dropna(subset=temp_cols)
            if df.empty:
                print("No entries with valid min_T and max_T after cleaning.")
                return

            # 8) Define the two disjoint groups:
            #    - "below_min": proteins whose entire observed range is < min_T  (max_T < min_T)
            #    - "above_max": proteins whose entire observed range is > max_T (min_T > max_T)
            below_min_df = df[df["max_temperature"] < min_T].copy()
            above_max_df = df[df["min_temperature"] > max_T].copy()

            # Keep only proteins uniquely in one bin (disjoint sets)
            # (If the same Gene somehow satisfies both, this removes overlaps.)
            if "Gene" in df.columns:
                below_genes = set(below_min_df["Gene"])
                above_genes = set(above_max_df["Gene"])
                overlap = below_genes & above_genes
                if overlap:
                    below_min_df = below_min_df[~below_min_df["Gene"].isin(overlap)]
                    above_max_df = above_max_df[~above_max_df["Gene"].isin(overlap)]

            if below_min_df.empty and above_max_df.empty:
                print(
                    "No proteins found with max_T < Min T or min_T > Max T.\n"
                    "Try adjusting your temperature thresholds."
                )
                return

            # 9) Tag groups
            below_min_df["temp_bin"] = f"< {min_T}°C"
            above_max_df["temp_bin"] = f"> {max_T}°C"

            combined = pd.concat(
                [df_ for df_ in [below_min_df, above_max_df] if not df_.empty],
                ignore_index=True
            )

            # 10) Coerce selected IDR properties to numeric
            for col in selected_props:
                if col in combined.columns:
                    combined[col] = pd.to_numeric(combined[col], errors="coerce")

            # 11) Report and preview
            print("Temperature-based groups (omrgcv2 only):\n")
            if not below_min_df.empty:
                print(f"  Proteins only in cold bin (< {min_T}°C): {len(below_min_df)}")
            if not above_max_df.empty:
                print(f"  Proteins only in warm bin (> {max_T}°C): {len(above_max_df)}")
            print()
            display(combined[["Gene", "temp_bin"] + selected_props].head())

            # 12) Plot comparison of selected properties between temp bins
            if "temp_bin" in combined.columns:
                try:
                    plot_idr_property_boxplots2(
                        combined,
                        properties=selected_props,
                        group_col="temp_bin"
                    )
                except Exception as e:
                    print("\nPlotting error:", e)

            # 13) Store for export
            LAST_RESULTS = combined

    run_button.on_click(on_run_clicked)

    ui = widgets.VBox([
        instructions,
        widgets.HTML("<b>Select IDR properties:</b>"),
        props_box,
        widgets.HBox([min_temp_input, max_temp_input]),
        run_button,
        output
    ])

    return ui
