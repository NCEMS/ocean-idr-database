# python-scripts/idr_analysis_setup.py

import os
import sqlite3
from itertools import combinations
from typing import List
from datetime import datetime

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from IPython.display import display, clear_output
import ipywidgets as widgets

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# --------------------
# Config
# --------------------

# The database must contain:
# architecture, cog_root, idr_props, environment_ranges
DB_PATH = "minimal_noenv-12May2026.db"

DATASET_OPTIONS = [
    "omrgcv2",
    "matou15",
    "3300003177",
    "3300003178",
    "3300003872",
    "3300003873",
]

DEFAULT_DATASETS = list(DATASET_OPTIONS)

ENV_DATASETS = ["omrgcv2", "matou15"]

MASTER_IDR_PROPS = [
    "Asphericity",
    "Radius_of_gyration",
    "Radius_of_gyration_scaled",
    "End_to_end_distance",
    "End_to_end_distance_scaled",
    "Scaling_exponent",
    "Prefactor",
    "Kappa",
    "Length",
    "FCR",
    "NCPR",
    "SHD",
    "SCD",
    "Molecular_weight",
    "F_Neg",
    "F_Pos",
    "Hydrophobicity",
    "Fraction_aromatic",
    "Fraction_aliphatic",
    "Fraction_polar",
    "Complexity",
]

DEFAULT_SELECTED = {
    "FCR",
    "NCPR",
    "SHD",
    "SCD",
    "F_Neg",
    "F_Pos",
    "Hydrophobicity",
}

LAST_RESULTS = None
LAST_PARAMS = None


# --------------------
# Core helpers
# --------------------

def get_last_params():
    return LAST_PARAMS or {}


def _safe_concat(dfs, **kwargs):
    """Concat after dropping empty/None and aligning columns to avoid FutureWarning."""
    dfs = [d for d in dfs if d is not None and hasattr(d, "empty") and not d.empty]
    if not dfs:
        return pd.DataFrame()

    cols = sorted(set().union(*(d.columns for d in dfs)))
    dfs = [d.reindex(columns=cols) for d in dfs]

    return pd.concat(dfs, **kwargs)


def run_query(query: str, params: tuple = ()) -> pd.DataFrame:
    with sqlite3.connect(DB_PATH) as conn:
        return pd.read_sql_query(query, conn, params=params)


def query_by_cog_root(cog_root_value: str) -> pd.DataFrame:
    """
    Return architecture + idr_props + cog_root annotations for a given COG root.

    Environmental ranges are joined when available. Currently, environment_ranges
    is expected to be populated only for dataset == 'omrgcv2' or 'matou15'.
    """
    query = """
    SELECT
        a.*,
        i.*,
        c.cog_root_ids,
        e.min_temperature,
        e.max_temperature,
        e.min_depth,
        e.max_depth
    FROM architecture AS a
    JOIN idr_props AS i
        ON a.dataset = i.dataset
       AND a.Gene = i.Gene
    JOIN cog_root AS c
        ON a.dataset = c.dataset
       AND a.Gene = c.Gene
    LEFT JOIN environment_ranges AS e
        ON a.dataset = e.dataset
       AND a.Gene = e.Gene
    WHERE c.cog_root_ids = ?
    """
    return run_query(query, (cog_root_value,))


def check_idr_props_list(user_props: List[str]) -> None:
    invalid = [p for p in user_props if p not in MASTER_IDR_PROPS]

    if invalid:
        allowed_str = ", ".join(MASTER_IDR_PROPS)
        invalid_str = ", ".join(invalid)
        raise ValueError(
            f"Invalid properties: {invalid_str}. "
            f"Allowed properties are: {allowed_str}."
        )


def plot_idr_property_boxplots2(
    df: pd.DataFrame,
    properties: list[str],
    group_col: str = "dataset",
    ncols: int = 4,
    figsize: tuple = (16, 10),
):
    """
    Create subplots of boxplots for IDR properties grouped by `group_col`.

    For each property, perform pairwise Mann–Whitney U tests between groups.
    Apply Benjamini–Hochberg FDR correction across all pairwise comparisons
    for all properties together.

    If any BH-adjusted p-value for a property is < 0.05, the title is {prop}.
    """

    if not properties:
        raise ValueError("No properties provided to plot.")

    if group_col not in df.columns:
        raise ValueError(f"Expected column '{group_col}' not found in results.")

    groups = df[group_col].dropna().unique()

    if len(groups) < 2:
        raise ValueError(
            f"Need at least 2 groups in '{group_col}' for statistical comparison; "
            f"found {len(groups)}."
        )

    pairs = list(combinations(groups, 2))

    # --------------------
    # First pass: compute all raw p-values
    # --------------------
    all_tests = []

    for prop in properties:
        if prop not in df.columns:
            continue

        for a, b in pairs:
            group_a = df.loc[df[group_col] == a, prop].dropna()
            group_b = df.loc[df[group_col] == b, prop].dropna()

            if len(group_a) > 0 and len(group_b) > 0:
                try:
                    _, p = mannwhitneyu(
                        group_a,
                        group_b,
                        alternative="two-sided",
                    )
                except ValueError:
                    continue

                all_tests.append((prop, f"{a} vs {b}", float(p)))

    # --------------------
    # If no tests, plot without stats
    # --------------------
    if not all_tests:
        nplots = len(properties)
        nrows = (nplots + ncols - 1) // ncols

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=figsize,
            squeeze=False,
        )
        axes = axes.flatten()

        last_ax_index = -1

        for i, prop in enumerate(properties):
            ax = axes[i]
            last_ax_index = i

            if prop not in df.columns:
                ax.set_visible(False)
                continue

            sns.boxplot(data=df, x=group_col, y=prop, ax=ax, fliersize=1)
            ax.set_xlabel("")
            ax.set_ylabel(prop)
            ax.tick_params(axis="x", rotation=45)
            ax.set_title(prop)

        for j in range(last_ax_index + 1, len(axes)):
            axes[j].set_visible(False)

        fig.tight_layout()
        plt.show()
        return

    # --------------------
    # Benjamini–Hochberg FDR across all tests
    # --------------------
    raw_pvals = [t[2] for t in all_tests]
    _, qvals, _, _ = multipletests(raw_pvals, alpha=0.05, method="fdr_bh")

    qmap = {}

    for (prop, comp, _p), q in zip(all_tests, qvals):
        qmap[(prop, comp)] = q

    # --------------------
    # Plot with q-values
    # --------------------
    nplots = len(properties)
    nrows = (nplots + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize,
        squeeze=False,
    )
    axes = axes.flatten()

    last_ax_index = -1

    for i, prop in enumerate(properties):
        ax = axes[i]
        last_ax_index = i

        if prop not in df.columns:
            ax.set_visible(False)
            continue

        sns.boxplot(data=df, x=group_col, y=prop, ax=ax, fliersize=1)
        ax.set_xlabel("")
        ax.set_ylabel(prop)
        ax.tick_params(axis="x", rotation=45)

        prop_qvals = {
            comp: q
            for (p, comp), q in qmap.items()
            if p == prop
        }

        significant = any(q < 0.05 for q in prop_qvals.values())
        title = f"{{{prop}}}" if significant else prop
        ax.set_title(title)

        if prop_qvals:
            if len(prop_qvals) <= 3:
                lines = []

                for comp, q in sorted(prop_qvals.items()):
                    q_str = "<0.001" if q < 0.001 else f"{q:.3f}"
                    lines.append(f"{comp}: q={q_str}")

                ax.text(
                    0.05,
                    0.95,
                    "\n".join(lines),
                    transform=ax.transAxes,
                    fontsize=8,
                    va="top",
                    ha="left",
                    bbox=dict(
                        boxstyle="round",
                        facecolor="white",
                        alpha=0.6,
                    ),
                )
            else:
                print(f"\nBH-FDR q-values for {prop}:")

                for comp, q in sorted(prop_qvals.items()):
                    q_str = "<0.001" if q < 0.001 else f"{q:.3f}"
                    print(f"  {comp}: q={q_str}")

    for j in range(last_ax_index + 1, len(axes)):
        axes[j].set_visible(False)

    fig.tight_layout()
    plt.show()


# --------------------
# UI factory
# --------------------

def create_idr_ui():
    """Build and return the full interactive UI as a widget container."""

    instructions = widgets.HTML(
        "<b>Instructions:</b><br>"
        "1) Enter a COG ID.<br>"
        "2) Select dataset(s).<br>"
        "3) Optionally set a minimum IDR length.<br>"
        "4) Optionally filter by architecture.<br>"
        "5) Optionally define independent temperature bins for OM-RGC.v2 "
        "(omrgcv2) and MATOU (matou15).<br>"
        "&nbsp;&nbsp;&nbsp;&nbsp;Cold bin: max_temperature &lt; cold threshold.<br>"
        "&nbsp;&nbsp;&nbsp;&nbsp;Hot bin: min_temperature &gt; hot threshold.<br>"
        "&nbsp;&nbsp;&nbsp;&nbsp;Rows between thresholds are dropped for that dataset when binning is active.<br>"
        "6) Select IDR properties.<br>"
        "7) Click <i>Run query and plot</i>."
    )

    # ---- COG ID input ----
    cog_input = widgets.Text(
        description="COG ID:",
        placeholder="e.g. COG0513",
        layout=widgets.Layout(width="250px"),
    )

    # ---- Dataset selection ----
    dataset_checkboxes = [
        widgets.Checkbox(
            value=(ds in DEFAULT_DATASETS),
            description=ds,
            indent=False,
        )
        for ds in DATASET_OPTIONS
    ]

    dataset_box = widgets.VBox(dataset_checkboxes)

    # ---- Min IDR length filter ----
    idr_length_input = widgets.IntText(
        description="Min IDR len:",
        value=0,
        layout=widgets.Layout(width="200px"),
    )

    # ---- Architecture filter ----
    arch_input = widgets.Text(
        description="Architectures:",
        placeholder="e.g. IDR-ORDERED, ORDERED-IDR",
        layout=widgets.Layout(width="400px"),
    )

    # ---- Temperature filter inputs for omrgcv2 ----
    omrg_min_temp_input = widgets.Text(
        description="cold <",
        placeholder="optional",
        layout=widgets.Layout(width="200px"),
    )

    omrg_max_temp_input = widgets.Text(
        description="hot >",
        placeholder="optional",
        layout=widgets.Layout(width="200px"),
    )

    # ---- Temperature filter inputs for matou15 ----
    matou_min_temp_input = widgets.Text(
        description="cold <",
        placeholder="optional",
        layout=widgets.Layout(width="200px"),
    )

    matou_max_temp_input = widgets.Text(
        description="hot >",
        placeholder="optional",
        layout=widgets.Layout(width="200px"),
    )

    # ---- IDR property selection ----
    prop_checkboxes = [
        widgets.Checkbox(
            value=(prop in DEFAULT_SELECTED),
            description=prop,
            indent=False,
        )
        for prop in MASTER_IDR_PROPS
    ]

    props_box = widgets.VBox(prop_checkboxes)

    # ---- Run button + output ----
    run_button = widgets.Button(
        description="▶ Run query and plot",
        button_style="primary",
        layout=widgets.Layout(width="250px"),
    )

    output = widgets.Output()

    # --------------------
    # Helper functions inside UI factory
    # --------------------

    def parse_and_validate_arch_list(raw: str):
        text = raw.strip()

        if not text:
            return [], None

        raw_items = [item.strip() for item in text.split(",")]
        archs = [a for a in raw_items if a]

        if not archs:
            return [], None

        valid_archs = []
        invalid_archs = []

        for arch in archs:
            parts = arch.upper().split("-")

            if parts and all(p in ("IDR", "ORDERED") for p in parts):
                valid_archs.append("-".join(parts))
            else:
                invalid_archs.append(arch)

        if invalid_archs:
            msg = (
                "Invalid architecture value(s): "
                + ", ".join(invalid_archs)
                + ". Each must be combinations of 'IDR' and 'ORDERED' joined by '-', "
                "e.g. IDR, ORDERED, IDR-ORDERED-IDR."
            )
            return None, msg

        seen = set()
        uniq = []

        for a in valid_archs:
            if a not in seen:
                seen.add(a)
                uniq.append(a)

        return uniq, None

    def parse_float_or_none(s: str):
        s = str(s).strip()

        if not s:
            return None

        try:
            return float(s)
        except ValueError:
            return None

    def _collect_current_params_for_save():
        selected_datasets = [
            cb.description for cb in dataset_checkboxes if cb.value
        ]

        selected_props = [
            cb.description for cb in prop_checkboxes if cb.value
        ]

        return {
            "query_datetime": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "cog_id": cog_input.value.strip(),
            "datasets": ", ".join(selected_datasets),
            "min_idr_length": str(idr_length_input.value),
            "architectures": arch_input.value.strip(),
            "omrgcv2_cold_temperature": omrg_min_temp_input.value.strip(),
            "omrgcv2_hot_temperature": omrg_max_temp_input.value.strip(),
            "matou15_cold_temperature": matou_min_temp_input.value.strip(),
            "matou15_hot_temperature": matou_max_temp_input.value.strip(),
            "selected_properties": ", ".join(selected_props),
        }

    def apply_temperature_binning(results: pd.DataFrame, selected_datasets: list[str]):
        """
        Apply independent hot/cold temperature binning for omrgcv2 and matou15.

        For each environmental dataset:
        - If both thresholds are blank, retain the dataset unbinned.
        - If both thresholds are provided and cold < hot, split into:
            cold: max_temperature < cold threshold
            hot:  min_temperature > hot threshold
        - Rows between thresholds are dropped for that dataset.
        - Other selected datasets are retained unbinned.
        """
        temp_filter_active = False

        bin_specs = {
            "omrgcv2": {
                "label": "Tara Oceans",
                "cold": parse_float_or_none(omrg_min_temp_input.value),
                "hot": parse_float_or_none(omrg_max_temp_input.value),
            },
            "matou15": {
                "label": "MATOU",
                "cold": parse_float_or_none(matou_min_temp_input.value),
                "hot": parse_float_or_none(matou_max_temp_input.value),
            },
        }

        binned_dfs = []
        unbinned_dfs = []

        # Process omrgcv2 and matou15 independently.
        for ds, spec in bin_specs.items():
            ds_df = results[results["dataset"] == ds].copy()

            if ds_df.empty:
                continue

            cold_T = spec["cold"]
            hot_T = spec["hot"]

            # No thresholds for this dataset: keep this dataset unbinned.
            if cold_T is None and hot_T is None:
                unbinned_dfs.append(ds_df)
                print(f"No temperature binning applied for {spec['label']} ({ds}).")
                continue

            # Partial or invalid thresholds: keep this dataset unbinned.
            if cold_T is None or hot_T is None or cold_T >= hot_T:
                unbinned_dfs.append(ds_df)
                print(
                    f"Temperature binning not applied for {spec['label']} ({ds}): "
                    "provide both cold and hot thresholds with cold < hot. "
                    "Rows for this dataset are retained unbinned."
                )
                continue

            required_temp_cols = {"min_temperature", "max_temperature"}

            if not required_temp_cols.issubset(ds_df.columns):
                unbinned_dfs.append(ds_df)
                print(
                    f"Temperature binning requested for {spec['label']} ({ds}), "
                    "but min_temperature / max_temperature columns were not found. "
                    "Rows for this dataset are retained unbinned."
                )
                continue

            ds_df["min_temperature"] = pd.to_numeric(
                ds_df["min_temperature"],
                errors="coerce",
            )

            ds_df["max_temperature"] = pd.to_numeric(
                ds_df["max_temperature"],
                errors="coerce",
            )

            before_dropna = len(ds_df)

            ds_df = ds_df.dropna(
                subset=["min_temperature", "max_temperature"]
            ).copy()

            dropped_missing = before_dropna - len(ds_df)

            if dropped_missing > 0:
                print(
                    f"Dropped {dropped_missing} {spec['label']} ({ds}) rows "
                    "with missing temperature metadata."
                )

            if ds_df.empty:
                print(
                    f"No usable temperature metadata for {spec['label']} ({ds}); "
                    "dropping this dataset from the binned comparison."
                )
                continue

            print(
                f"Applying temperature bins for {spec['label']} ({ds}): "
                f"cold = max_temperature < {cold_T}; "
                f"hot = min_temperature > {hot_T}"
            )

            cold_df = ds_df[ds_df["max_temperature"] < cold_T].copy()
            hot_df = ds_df[ds_df["min_temperature"] > hot_T].copy()

            # Ensure no gene appears in both bins for this dataset.
            if "Gene" in ds_df.columns:
                cold_genes = set(cold_df["Gene"])
                hot_genes = set(hot_df["Gene"])
                overlap = cold_genes & hot_genes

                if overlap:
                    cold_df = cold_df[~cold_df["Gene"].isin(overlap)].copy()
                    hot_df = hot_df[~hot_df["Gene"].isin(overlap)].copy()

                    print(
                        f"Removed {len(overlap)} {spec['label']} ({ds}) genes "
                        "that overlapped between cold and hot bins."
                    )

            if cold_df.empty and hot_df.empty:
                print(
                    f"No {spec['label']} ({ds}) rows fall into the requested "
                    "cold or hot bins. Rows for this dataset are dropped."
                )
                temp_filter_active = True
                continue

            if not cold_df.empty:
                cold_df["temp_bin"] = f"{ds}: cold < {cold_T}°C"
                binned_dfs.append(cold_df)
                print(f"{spec['label']} ({ds}) cold-bin rows: {len(cold_df)}")

            if not hot_df.empty:
                hot_df["temp_bin"] = f"{ds}: hot > {hot_T}°C"
                binned_dfs.append(hot_df)
                print(f"{spec['label']} ({ds}) hot-bin rows: {len(hot_df)}")

            temp_filter_active = True

        # Retain non-environment datasets unbinned.
        other_df = results[~results["dataset"].isin(bin_specs.keys())].copy()

        if not other_df.empty:
            unbinned_dfs.append(other_df)

        updated_results = _safe_concat(
            unbinned_dfs + binned_dfs,
            ignore_index=True,
        )

        if temp_filter_active:
            print(f"After temperature binning: total rows = {len(updated_results)}")
        else:
            print("No temperature binning applied.")

        return updated_results, temp_filter_active

    # --------------------
    # Callback
    # --------------------

    def on_run_clicked(b):
        global LAST_RESULTS, LAST_PARAMS

        with output:
            clear_output()

            # 1. COG ID
            cog_id = cog_input.value.strip()

            if not cog_id:
                print("Please enter a COG ID.")
                return

            # 2. Datasets
            selected_datasets = [
                cb.description for cb in dataset_checkboxes if cb.value
            ]

            if not selected_datasets:
                print("Please select at least one dataset.")
                return

            # 3. IDR properties
            selected_props = [
                cb.description for cb in prop_checkboxes if cb.value
            ]

            if not selected_props:
                print("Please select at least one IDR property.")
                return

            try:
                check_idr_props_list(selected_props)
            except ValueError as e:
                print(e)
                return

            # 4. Architectures
            arch_list, arch_err = parse_and_validate_arch_list(arch_input.value)

            if arch_err:
                print(arch_err)
                return

            if arch_list:
                print("Architecture filter requested for:")

                for a in arch_list:
                    print(f"  - {a}")

            # 5. Base query
            print("Running query...")

            try:
                results = query_by_cog_root(cog_id)
            except Exception as e:
                print("Database query failed:")
                print(e)
                return

            if results is None or results.empty:
                print(f"No results found for COG ID: {cog_id}")
                return

            # Drop duplicated column names created by SELECT a.*, i.*.
            # The first occurrence is retained.
            results = results.loc[:, ~results.columns.duplicated()].copy()

            print(f"Retrieved {len(results)} rows after base join.")

            # 6. Dataset filter
            if "dataset" not in results.columns:
                print("Warning: 'dataset' column not found; cannot filter by dataset.")
                return

            results = results[results["dataset"].isin(selected_datasets)].copy()

            if results.empty:
                print(
                    "No results found for COG ID "
                    f"{cog_id} in selected dataset(s): "
                    f"{', '.join(selected_datasets)}"
                )
                return

            print(f"{len(results)} rows remain after dataset filter.")

            # 7. Minimum IDR length filter
            min_len = idr_length_input.value

            if isinstance(min_len, int) and min_len > 0:
                if "Length" not in results.columns:
                    print(
                        "Min IDR length specified, but 'Length' column was not found."
                    )
                    return

                results["Length"] = pd.to_numeric(
                    results["Length"],
                    errors="coerce",
                )

                before = len(results)
                results = results[results["Length"] >= min_len].copy()

                if results.empty:
                    print(f"No IDRs with Length ≥ {min_len} after filtering.")
                    return

                print(
                    f"Applied Length ≥ {min_len} filter: "
                    f"{before} → {len(results)} rows."
                )
            else:
                print("No minimum IDR length filter applied.")

            # 8. Architecture filter
            if arch_list:
                if "architecture" not in results.columns:
                    print(
                        "Architecture filter requested, but 'architecture' column "
                        "not found."
                    )
                    return

                before = len(results)
                results = results[results["architecture"].isin(arch_list)].copy()

                if results.empty:
                    print(
                        "No entries matched the requested architecture(s): "
                        + ", ".join(arch_list)
                    )
                    return

                print(
                    f"Applied architecture filter: {before} → {len(results)} rows "
                    f"matching {', '.join(arch_list)}."
                )
            else:
                print("No architecture filter applied.")

            # 9. Independent temperature binning for omrgcv2 and matou15
            if "dataset" not in results.columns:
                print("Temperature binning error: 'dataset' column missing.")
                return

            results, temp_filter_active = apply_temperature_binning(
                results,
                selected_datasets,
            )

            if results.empty:
                print("No results left after applying all filters.")
                return

            # 10. Coerce selected properties to numeric
            for col in selected_props:
                if col in results.columns:
                    results[col] = pd.to_numeric(results[col], errors="coerce")
                else:
                    print(f"Warning: selected property '{col}' not found in results.")

            # 11. Summary + preview
            print(f"\nFinal row count: {len(results)}")

            if "dataset" in results.columns:
                print("\nNumber of IDRs per dataset:")
                print(results["dataset"].value_counts())
                print("")

            if temp_filter_active and "temp_bin" in results.columns:
                print("Number of IDRs per temperature bin:")
                print(results["temp_bin"].dropna().value_counts())
                print("")

            print("Preview of filtered results:")
            display(results.head())

            # 12. Plots
            if temp_filter_active and "temp_bin" in results.columns:
                plot_df = results.copy()

                if "dataset" not in plot_df.columns:
                    print("Plotting error: 'dataset' column missing.")
                else:
                    plot_df["plot_group"] = plot_df["dataset"].astype(str)

                    mask_binned = plot_df["temp_bin"].notna()

                    plot_df.loc[mask_binned, "plot_group"] = (
                        plot_df.loc[mask_binned, "temp_bin"].astype(str)
                    )

                    print(
                        "\nPlotting IDR properties by group, where:\n"
                        "  - omrgcv2 and/or matou15 may be split into temperature bins\n"
                        "  - other datasets are shown by dataset name\n"
                    )

                    try:
                        plot_idr_property_boxplots2(
                            plot_df,
                            properties=selected_props,
                            group_col="plot_group",
                        )
                    except Exception as e:
                        print("\nPlotting error using plot_group:", e)

            else:
                if "dataset" not in results.columns:
                    print("Plotting error: 'dataset' column missing.")
                else:
                    try:
                        plot_idr_property_boxplots2(
                            results,
                            properties=selected_props,
                            group_col="dataset",
                        )
                    except Exception as e:
                        print("\nPlotting error using dataset:", e)

            # 13. Store results
            LAST_RESULTS = results
            LAST_PARAMS = _collect_current_params_for_save()

            try:
                LAST_RESULTS.attrs["idr_params"] = LAST_PARAMS
            except Exception:
                pass

            print("\nResults stored in LAST_RESULTS and parameters captured for saving.")

    run_button.on_click(on_run_clicked)

    # --------------------
    # Layout
    # --------------------

    ui = widgets.VBox(
        [
            instructions,
            cog_input,
            widgets.HTML("<b>Select dataset(s):</b>"),
            dataset_box,
            widgets.HTML("<b>Optional: set minimum IDR length:</b>"),
            idr_length_input,
            widgets.HTML("<b>Optional: filter by architecture, comma-separated:</b>"),
            arch_input,
            widgets.HTML("<b>Optional: OM-RGC.v2 (omrgcv2) temperature bins:</b>"),
            widgets.HBox([omrg_min_temp_input, omrg_max_temp_input]),
            widgets.HTML("<b>Optional: MATOU (matou15) temperature bins:</b>"),
            widgets.HBox([matou_min_temp_input, matou_max_temp_input]),
            widgets.HTML("<b>Select IDR properties:</b>"),
            props_box,
            run_button,
            output,
        ]
    )

    return ui