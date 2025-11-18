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
LAST_PARAMS  = None

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
    # align columns to the union so no all-NA columns sneak in differently-typed
    cols = sorted(set().union(*(d.columns for d in dfs)))
    dfs = [d.reindex(columns=cols) for d in dfs]
    return pd.concat(dfs, **kwargs)

def run_query(query: str, params: tuple = ()) -> pd.DataFrame:
    with sqlite3.connect(DB_PATH) as conn:
        return pd.read_sql_query(query, conn, params=params)

def query_by_cog_root(cog_root_value: str) -> pd.DataFrame:
    """
    Return architecture + idr_props + cog_root for a given COG,
    with optional Tara Oceans temperature data joined when available.

    Temperature info (min_temperature, max_temperature) is only present
    for omrgcv2 entries in the tara_temperature table.
    """
    query = """
    SELECT
        a.*,
        i.*,
        c.cog_root_ids,
        t.min_temperature,
        t.max_temperature
    FROM architecture AS a
    JOIN idr_props AS i
        ON a.Gene = i.Gene
    JOIN cog_root AS c
        ON a.Gene = c.Gene
    LEFT JOIN tara_temperatures AS t
        ON a.Gene = t.Gene
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

def plot_idr_property_boxplots2(df: pd.DataFrame,
                                properties: list[str],
                                group_col: str = "dataset",
                                ncols: int = 4,
                                figsize: tuple = (16, 10)):
    """
    Create subplots of boxplots for IDR properties grouped by `group_col`.
    For each property, perform pairwise Mann–Whitney U tests between groups.

    Apply Benjamini–Hochberg FDR correction across ALL pairwise comparisons
    for ALL properties together.

    - If any BH-adjusted p-value (q) for a property is < 0.05, its title is {prop}.
    - Annotations show BH-adjusted p-values (q-values).
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

    from itertools import combinations
    pairs = list(combinations(groups, 2))

    # --------------------
    # First pass: compute all raw p-values
    # --------------------
    all_tests = []  # (prop, "A vs B", p)

    for prop in properties:
        if prop not in df.columns:
            continue

        for a, b in pairs:
            group_a = df.loc[df[group_col] == a, prop].dropna()
            group_b = df.loc[df[group_col] == b, prop].dropna()

            if len(group_a) > 0 and len(group_b) > 0:
                try:
                    _, p = mannwhitneyu(group_a, group_b, alternative="two-sided")
                except ValueError:
                    # e.g., identical values or other numerical quirks
                    continue
                all_tests.append((prop, f"{a} vs {b}", float(p)))

    # If no tests, just plot without stats
    if not all_tests:
        nplots = len(properties)
        nrows = (nplots + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        axes = axes.flatten()

        for i, prop in enumerate(properties):
            ax = axes[i]
            if prop not in df.columns:
                ax.set_visible(False)
                continue

            sns.boxplot(data=df, x=group_col, y=prop, ax=ax, fliersize=1)
            ax.set_xlabel("")
            ax.set_ylabel(prop)
            ax.tick_params(axis="x", rotation=45)
            ax.set_title(prop)

        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        fig.tight_layout()
        plt.show()
        return

    # --------------------
    # Benjamini–Hochberg FDR across all tests
    # --------------------
    raw_pvals = [t[2] for t in all_tests]
    # multipletests returns: reject, pvals_corrected, alphacSidak, alphacBonf
    _, qvals, _, _ = multipletests(raw_pvals, alpha=0.05, method="fdr_bh")

    # Map (prop, comparison) -> q
    qmap = {}
    for (prop, comp, _p), q in zip(all_tests, qvals):
        qmap[(prop, comp)] = q

    # --------------------
    # Plot with q-values
    # --------------------
    nplots = len(properties)
    nrows = (nplots + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
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

        # Get this property's q-values
        prop_qvals = {
            comp: q
            for (p, comp), q in qmap.items()
            if p == prop
        }

        # Mark title if any q < 0.05
        significant = any(q < 0.05 for q in prop_qvals.values())
        title = f"{{{prop}}}" if significant else prop
        ax.set_title(title)

        # Annotate if reasonable number of comparisons
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
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.6),
                )
            else:
                # Too many comparisons; print to stdout to avoid clutter
                print(f"\nBH-FDR q-values for {prop}:")
                for comp, q in sorted(prop_qvals.items()):
                    q_str = "<0.001" if q < 0.001 else f"{q:.3f}"
                    print(f"  {comp}: q={q_str}")

    # Hide unused axes
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
        "1) Enter a COG ID. <br>"
        "2) Select dataset(s). <br>"
        "3) (Optional) Set a minimum IDR length. <br>"
        "4) (Optional) Filter by architecture (comma-separated; each must be 'IDR'/'ORDERED' combos, e.g. 'IDR-ORDERED-IDR').<br>"
        "5) (Optional) For Tara Oceans (omrgcv2), set Min/Max temperature to define two bins:<br>"
        "&nbsp;&nbsp;&nbsp;&nbsp;Cold: max_temperature &lt; Min T; Warm: min_temperature &gt; Max T (others dropped for omrgcv2).<br>"
        "6) Select IDR properties. <br>"
        "7) Click <i>'Run query and plot'</i>."
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

    # ---- Min IDR length filter (optional) ----
    idr_length_input = widgets.IntText(
        description="Min IDR len:",
        value=0,
        layout=widgets.Layout(width="200px")
    )

    # ---- Architecture filter (optional, comma-separated) ----
    arch_input = widgets.Text(
        description="Architectures:",
        placeholder="e.g. IDR-ORDERED, ORDERED-IDR",
        layout=widgets.Layout(width="400px")
    )

    # ---- Temperature filter inputs for omrgcv2 (optional) ----
    min_temp_input = widgets.Text(
        description="Min T (°C):",
        placeholder="optional",
        layout=widgets.Layout(width="200px")
    )
    max_temp_input = widgets.Text(
        description="Max T (°C):",
        placeholder="optional",
        layout=widgets.Layout(width="200px")
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

    # ---- Run button + output ----
    run_button = widgets.Button(
        description="▶ Run query and plot",
        button_style="primary",
        layout=widgets.Layout(width="250px")
    )

    output = widgets.Output()

    # ---- Helpers ----

    def parse_and_validate_arch_list(raw: str):
        text = raw.strip()
        if not text:
            return [], None  # no filter

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

        # Deduplicate, preserve order
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

    # ---- Callback ----
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

            # 4. Architectures (optional)
            arch_list, arch_err = parse_and_validate_arch_list(arch_input.value)
            if arch_err:
                print(arch_err)
                return
            if arch_list:
                print("Architecture filter requested for:")
                for a in arch_list:
                    print(f"  - {a}")

            # 5. Base query (single join, includes temps via LEFT JOIN)
            print("Running query...")
            results = query_by_cog_root(cog_id)
            if results is None or results.empty:
                print(f"No results found for COG ID: {cog_id}")
                return

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

            # 7. Min IDR length filter
            min_len = idr_length_input.value
            if isinstance(min_len, int) and min_len > 0:
                if "Length" not in results.columns:
                    print(
                        "Min IDR length specified, but 'Length' column was not found."
                    )
                    return

                results["Length"] = pd.to_numeric(results["Length"], errors="coerce")
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
                        "Architecture filter requested, but 'architecture' column not found."
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

            # 9. Temperature-based Tara Oceans filter (omrgcv2 only)
            min_T = parse_float_or_none(min_temp_input.value)
            max_T = parse_float_or_none(max_temp_input.value)

            temp_filter_active = False

            if (
                min_T is not None
                and max_T is not None
                and min_T < max_T
                and "omrgcv2" in selected_datasets
            ):
                print(
                    f"Applying Tara Oceans temperature bins for omrgcv2:"
                    f" cold: max_temperature < {min_T},"
                    f" warm: min_temperature > {max_T}"
                )

                if "min_temperature" not in results.columns or "max_temperature" not in results.columns:
                    print(
                        "Temperature filter requested, but 'min_temperature' / "
                        "'max_temperature' columns not found in results."
                    )
                else:
                    # Split omrgcv2 vs others
                    omrg_mask = results["dataset"] == "omrgcv2"
                    omrg_df = results[omrg_mask].copy()
                    other_df = results[~omrg_mask].copy()

                    if not omrg_df.empty:
                        # Ensure numeric & drop missing
                        omrg_df["min_temperature"] = pd.to_numeric(
                            omrg_df["min_temperature"], errors="coerce"
                        )
                        omrg_df["max_temperature"] = pd.to_numeric(
                            omrg_df["max_temperature"], errors="coerce"
                        )
                        omrg_df = omrg_df.dropna(
                            subset=["min_temperature", "max_temperature"]
                        )

                        cold_df = omrg_df[
                            omrg_df["max_temperature"] < min_T
                        ].copy()
                        warm_df = omrg_df[
                            omrg_df["min_temperature"] > max_T
                        ].copy()

                        # (Optional) ensure unique by Gene across bins
                        if "Gene" in omrg_df.columns:
                            cold_genes = set(cold_df["Gene"])
                            warm_genes = set(warm_df["Gene"])
                            overlap = cold_genes & warm_genes
                            if overlap:
                                cold_df = cold_df[~cold_df["Gene"].isin(overlap)]
                                warm_df = warm_df[~warm_df["Gene"].isin(overlap)]

                        if cold_df.empty and warm_df.empty:
                            print(
                                "No omrgcv2 proteins fall into cold or warm bins; "
                                "omrgcv2 entries will be dropped, others retained."
                            )
                            results = other_df
                        else:
                            if not cold_df.empty:
                                cold_df["temp_bin"] = f"< {min_T}°C"
                                print(f"Cold-bin omrgcv2 rows: {len(cold_df)}")
                            if not warm_df.empty:
                                warm_df["temp_bin"] = f"> {max_T}°C"
                                print(f"Warm-bin omrgcv2 rows: {len(warm_df)}")

                            binned = _safe_concat(
                                [d for d in (cold_df, warm_df) if not d.empty],
                                ignore_index=True
                            )
                            results = _safe_concat(
                                [other_df, binned],
                                ignore_index=True
                            )
                            temp_filter_active = True
                            print(
                                f"After temperature binning: total rows = {len(results)}"
                            )
                    else:
                        print(
                            "No omrgcv2 rows present after previous filters; "
                            "temperature filter has no effect."
                        )
            else:
                if min_T is None and max_T is None:
                    print("No temperature filter applied.")
                elif "omrgcv2" not in selected_datasets:
                    print("Temperature values provided but omrgcv2 not selected; ignoring.")
                else:
                    print(
                        "Temperature filter not applied: please provide both Min and Max with Min < Max."
                    )

            if results.empty:
                print("No results left after applying all filters.")
                return

            # 10. Coerce selected props numeric
            for col in selected_props:
                if col in results.columns:
                    results[col] = pd.to_numeric(results[col], errors="coerce")

            # 11. Summary + preview
            print(f"\nFinal row count: {len(results)}")
            if "dataset" in results.columns:
                print("\nNumber of IDRs per dataset:")
                print(results["dataset"].value_counts())
                print("")
            print("Preview of filtered results:")
            display(results.head())

            # 12. Plots
            if temp_filter_active and "temp_bin" in results.columns:
                # Build a unified grouping column:
                # - binned omrgcv2 use temp_bin labels
                # - all other rows use their dataset label
                plot_df = results.copy()

                if "dataset" not in plot_df.columns:
                    print("Plotting error: 'dataset' column missing.")
                else:
                    # Initialize with dataset
                    plot_df["plot_group"] = plot_df["dataset"].astype(str)

                    # Overwrite with temp_bin where available (for binned omrgcv2 rows)
                    mask_binned = plot_df["temp_bin"].notna()
                    plot_df.loc[mask_binned, "plot_group"] = plot_df.loc[mask_binned, "temp_bin"].astype(str)

                    print(
                        "\nPlotting IDR properties by group, where:\n"
                        f"  - Tara Oceans (omrgcv2) entries are split into bins via temp_bin\n"
                        f"  - Other datasets are shown by their dataset name\n"
                    )

                    try:
                        plot_idr_property_boxplots2(
                            plot_df,
                            properties=selected_props,
                            group_col="plot_group"
                        )
                    except Exception as e:
                        print("\nPlotting error (plot_group):", e)
            else:
                # Standard dataset-based plots
                if "dataset" not in results.columns:
                    print("Plotting error: 'dataset' column missing.")
                else:
                    try:
                        plot_idr_property_boxplots2(
                            results,
                            properties=selected_props,
                            group_col="dataset"
                        )
                    except Exception as e:
                        print("\nPlotting error (dataset):", e)

            # 13. Store results
            LAST_RESULTS = results
            LAST_PARAMS = _collect_current_params_for_save()
            try:
                LAST_RESULTS.attrs["idr_params"] = LAST_PARAMS
            except Exception:
                pass
            print("\nResults stored in LAST_RESULTS (and parameters captured for saving).")

    def _collect_current_params_for_save():
        selected_datasets = [cb.description for cb in dataset_checkboxes if cb.value]
        selected_props    = [cb.description for cb in prop_checkboxes if cb.value]
        return {
            "query_datetime": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "cog_id": cog_input.value.strip(),
            "datasets": ", ".join(selected_datasets),
            "min_idr_length": str(idr_length_input.value),
            "architectures": arch_input.value.strip(),
            "min_temperature": min_temp_input.value.strip(),
            "max_temperature": max_temp_input.value.strip(),
            "selected_properties": ", ".join(selected_props),
        }

    run_button.on_click(on_run_clicked)

    # ---- Layout ----
    ui = widgets.VBox([
        instructions,
        cog_input,
        widgets.HTML("<b>Select dataset(s):</b>"),
        dataset_box,
        widgets.HTML("<b>Optional: set minimum IDR length:</b>"),
        idr_length_input,
        widgets.HTML("<b>Optional: filter by architecture (comma-separated):</b>"),
        arch_input,
        widgets.HTML("<b>Optional: Tara Oceans (omrgcv2) temperature bins:</b>"),
        widgets.HBox([min_temp_input, max_temp_input]),
        widgets.HTML("<b>Select IDR properties:</b>"),
        props_box,
        run_button,
        output
    ])

    return ui
