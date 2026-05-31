#!/usr/bin/env python3
"""
Step 06: Plot Plasmid Library Representation and Quality Metrics

Generates visualizations for plasmid library quality assessment:
- GC content vs counts (bias detection)
- Count distribution density plots by SPS
- Gini index comparison (library uniformity)
- Category breakdown by perfect/imperfect matches

Input: Purity-filtered guide-donor-bc0 counts from Step 04
Output: Publication-ready plots in plots/ directory

Usage:
    python 06_plot_plasmid_library_representation.py
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# ============================================================================
# Configuration - Read from environment variables set by SLURM script
# ============================================================================

import os

def get_config_from_env():
    """Get project configuration from environment variables."""
    project_dir_str = os.environ.get("PROJECT_DIR")
    if not project_dir_str:
        raise ValueError(
            "PROJECT_DIR environment variable not set. "
            "Run via SLURM script or set: export PROJECT_DIR=/path/to/project"
        )

    project_dir = Path(project_dir_str)
    datatype = os.environ.get("DATATYPE", "guide_donor_bc0")

    processed_data_dir = project_dir / "processed_data" / datatype
    keyfile_dir = project_dir / "keyfiles" / datatype

    # Oligo design file from env or search keyfiles directory
    oligo_design_file = os.environ.get("OLIGO_DESIGN_FILE")
    if oligo_design_file:
        oligo_pool_files = [Path(oligo_design_file)]
    else:
        # Search for oligo design files in keyfiles directory
        oligo_pool_files = list(keyfile_dir.glob("*Twist*.tsv"))
        if not oligo_pool_files:
            raise ValueError(f"No oligo design files found in {keyfile_dir}")

    purity_dir = processed_data_dir / "guide_donor_bc0_purity_filtered"

    return {
        "project_dir": project_dir,
        "processed_data_dir": processed_data_dir,
        "keyfile_dir": keyfile_dir,
        "purity_dir": purity_dir,
        "plots_dir": purity_dir / "plots",
        "oligo_pool_files": oligo_pool_files,
    }

# File patterns
INPUT_SUFFIX = "_guide_donor_bc0_purity_filtered.tsv"


@dataclass
class OligoConstants:
    """Constants for oligo sequence parsing."""
    DONOR_START: int = 51
    DONOR_END: int = 180
    SPS_START: int = 180


@dataclass
class LibraryConfig:
    """Configuration for library processing."""
    # Library to SPS mapping
    library_sps_map: Dict[str, str] = field(default_factory=lambda: {
        "V629": "CTTGACTCGACAGTGAGAGC",
        "V631": "CGTATGGAGATGACGTCAGG",
        "V634": "CGATAGACGACTGGACAGCA",
    })

    # Display order for plots
    library_order: List[str] = field(default_factory=lambda: ["V629", "V631", "V634"])


# Global constants
OLIGO = OligoConstants()
LIB_CONFIG = LibraryConfig()


def setup_directories(env_config: dict) -> None:
    """Create necessary output directories."""
    env_config["plots_dir"].mkdir(parents=True, exist_ok=True)


def gc_content(seq: str) -> float:
    """
    Calculate the GC content of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC content as fraction (0.0 to 1.0)
    """
    if not isinstance(seq, str) or len(seq) == 0:
        return 0.0
    seq_upper = seq.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    return gc_count / len(seq_upper)


def gini_index(counts: np.ndarray) -> float:
    """
    Calculate the Gini index for a distribution of counts.

    The Gini index measures inequality in a distribution:
    - 0 = perfect equality (all oligos have same counts)
    - 1 = maximum inequality (one oligo has all counts)

    Args:
        counts: Array of count values

    Returns:
        Gini index (0.0 to 1.0) or NaN if invalid
    """
    if len(counts) == 0:
        return np.nan

    counts = np.array(counts)
    n = len(counts)

    if n == 0 or np.sum(counts) == 0:
        return np.nan

    sorted_counts = np.sort(counts)
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * sorted_counts) / np.sum(counts) - (n + 1)) / n


def load_data_files(input_dir: Path, suffix: str) -> pd.DataFrame:
    """
    Load all data files matching the suffix and combine into single DataFrame.

    Args:
        input_dir: Directory containing input files
        suffix: File suffix to match

    Returns:
        Combined DataFrame with library column
    """
    files = sorted(input_dir.glob(f"*{suffix}"))
    df_list = []

    for filepath in files:
        library = filepath.stem.replace(suffix.replace(".tsv", ""), "")
        print(f"  Loading: {filepath.name} (library: {library})")
        df = pd.read_csv(filepath, sep="\t")
        df["library"] = library
        df_list.append(df)

    if not df_list:
        raise FileNotFoundError(f"No files found matching pattern *{suffix}")

    return pd.concat(df_list, ignore_index=True)


def load_designed_oligos(env_config: dict) -> pd.DataFrame:
    """
    Load and process designed oligo sequences from all pool files.

    Args:
        env_config: Configuration dictionary with oligo_pool_files

    Returns:
        Combined DataFrame with extracted sequence features
    """
    dfs = []

    for filepath in env_config["oligo_pool_files"]:
        if not filepath.exists():
            print(f"  WARNING: Oligo pool not found: {filepath}")
            continue

        print(f"  Loading: {filepath.name}")
        df = pd.read_csv(filepath, sep="\t")
        dfs.append(df)

    if not dfs:
        raise FileNotFoundError("No oligo pool files found")

    combined = pd.concat(dfs, ignore_index=True)

    # Extract sequence features
    combined["designed_donor"] = (
        combined["oligo_seq"]
        .str[OLIGO.DONOR_START:OLIGO.DONOR_END]
        .str.upper()
    )
    combined["SPS"] = combined["oligo_seq"].str[OLIGO.SPS_START:].str.upper()
    combined["gc_content"] = combined["designed_donor"].apply(gc_content)
    combined["oligo_seq"] = combined["oligo_seq"].str.upper()

    return combined


def check_non_string_donors(df: pd.DataFrame, env_config: dict) -> None:
    """
    Identify and report non-string donor sequences.

    Args:
        df: DataFrame with donor column
        env_config: Configuration dictionary with plots_dir
    """
    if "donor" not in df.columns:
        return

    non_string_mask = ~df["donor"].apply(lambda x: isinstance(x, str))
    non_string_count = non_string_mask.sum()

    if non_string_count > 0:
        print(f"  WARNING: Found {non_string_count} non-string donor entries")
        non_string_df = df[non_string_mask]
        non_string_file = env_config["plots_dir"] / "non_string_donors.tsv"
        non_string_df.to_csv(non_string_file, sep="\t", index=False)


def add_missing_oligos(
    grouped_df: pd.DataFrame,
    designed_oligos: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add missing oligos with zero counts for each library.

    Args:
        grouped_df: Grouped oligo counts DataFrame
        designed_oligos: DataFrame with all designed oligos

    Returns:
        DataFrame with missing oligos added
    """
    result = grouped_df.copy()

    for library, sps in LIB_CONFIG.library_sps_map.items():
        expected = designed_oligos[designed_oligos["SPS"] == sps]
        existing = result[result["library"] == library]["oligo_name"]
        missing = expected[~expected["oligo_name"].isin(existing)]

        if missing.empty:
            continue

        # Check for sequence overlap
        if "oligo_seq" in result.columns:
            existing_seqs = result[result["library"] == library]["oligo_seq"]
            overlap = missing[missing["oligo_seq"].isin(existing_seqs)]

            if not overlap.empty:
                print(f"  WARNING: Library {library} has {len(overlap)} missing oligos "
                      f"with duplicate oligo_seq")
                missing = missing[~missing.index.isin(overlap.index)]

        if not missing.empty:
            print(f"  Library {library}: Adding {len(missing)} missing oligos")
            to_append = missing.assign(
                library=library,
                grouped_guide_bc0_counts=0,
            )[["library", "oligo_name", "grouped_guide_bc0_counts",
               "oligo_seq", "designed_donor", "gc_content", "SPS"]]
            result = pd.concat([result, to_append], ignore_index=True)

    return result


def plot_gc_vs_counts(df: pd.DataFrame, env_config: dict) -> None:
    """
    Plot GC content vs counts faceted by library.

    Args:
        df: DataFrame with gc_content and grouped_oligo_name_counts columns
        env_config: Configuration dictionary with plots_dir
    """
    sns.set_style("whitegrid")

    # Filter to libraries in our order
    plot_df = df[df["library"].isin(LIB_CONFIG.library_order)]

    g = sns.relplot(
        data=plot_df,
        x="gc_content",
        y="grouped_oligo_name_counts",
        col="library",
        kind="scatter",
        col_wrap=3,
        height=4,
        aspect=1.2,
        alpha=0.5,
        facet_kws={"sharey": False, "sharex": False},
        col_order=LIB_CONFIG.library_order,
    )

    g.set_axis_labels("Designed GC Content", "Total Counts")
    g.set(yscale="log")
    g.fig.suptitle("GC Content vs Total Counts by Library (log scale)", y=1.02)

    # Add R² and p-value annotations
    for library in LIB_CONFIG.library_order:
        lib_df = plot_df[plot_df["library"] == library]
        valid_df = lib_df[
            (lib_df["grouped_oligo_name_counts"] > 0) &
            lib_df["gc_content"].notnull()
        ]

        if len(valid_df) > 2:
            r, p = stats.pearsonr(
                valid_df["gc_content"],
                np.log10(valid_df["grouped_oligo_name_counts"])
            )
            r2 = r ** 2
            ax_idx = LIB_CONFIG.library_order.index(library)
            ax = g.axes.flat[ax_idx]
            ax.text(
                0.05, 0.95,
                f"$R^2$={r2:.2f}\np={p:.2g}",
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=10,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
            )

    plt.tight_layout()
    outfile = env_config["plots_dir"] / "all_libraries_gc_content_vs_counts_faceted_logy.png"
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outfile.name}")


def plot_density_by_sps(df: pd.DataFrame, env_config: dict) -> None:
    """
    Plot density of normalized counts separated by SPS.

    Args:
        df: DataFrame with normalized_counts, library, SPS columns
        env_config: Configuration dictionary with plots_dir
    """
    unique_sps = df["SPS"].dropna().unique()
    palette = sns.color_palette("tab10", n_colors=len(LIB_CONFIG.library_order))
    library_palette = dict(zip(LIB_CONFIG.library_order, palette))

    fig, axes = plt.subplots(1, len(unique_sps), figsize=(5 * len(unique_sps), 6), sharey=True)
    if len(unique_sps) == 1:
        axes = [axes]

    for idx, sps in enumerate(unique_sps):
        sps_df = df[df["SPS"] == sps]
        ax = axes[idx]

        # Filter to valid normalized counts
        valid_df = sps_df[sps_df["normalized_counts"] > 0]

        if valid_df.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"SPS: {sps[:10]}...")
            continue

        sns.kdeplot(
            data=valid_df,
            x="normalized_counts",
            hue="library",
            common_norm=False,
            fill=True,
            alpha=0.5,
            palette=library_palette,
            bw_adjust=0.5,
            hue_order=LIB_CONFIG.library_order,
            ax=ax,
        )

        ax.set_xscale("log")
        ax.set_xlabel("Normalized Counts (log scale)")
        ax.set_ylabel("Density" if idx == 0 else "")
        ax.set_title(f"SPS: {sps[:10]}...")

        # Custom legend
        present_libs = [lib for lib in LIB_CONFIG.library_order
                        if lib in valid_df["library"].unique()]
        handles = [plt.Line2D([0], [0], color=library_palette[lib], lw=4)
                   for lib in present_libs]
        ax.legend(handles, present_libs, title="Library",
                  bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.grid(axis="y", linestyle=":", linewidth=0.8)

    plt.suptitle("Density Plot of Normalized Counts by Library (by SPS)", y=1.02)
    plt.tight_layout()
    outfile = env_config["plots_dir"] / "normalized_counts_density_plot_by_SPS_side_by_side.png"
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outfile.name}")


def plot_gini_indices(df: pd.DataFrame, env_config: dict) -> None:
    """
    Calculate and plot Gini indices for each library.

    Args:
        df: DataFrame with library and grouped_oligo_name_counts columns
        env_config: Configuration dictionary with plots_dir
    """
    gini_data = (
        df.groupby("library")["grouped_oligo_name_counts"]
        .apply(lambda x: gini_index(x[x > 0].values))
        .reset_index(name="gini_index")
    )

    print("\n  Gini indices by library:")
    for _, row in gini_data.iterrows():
        print(f"    {row['library']}: {row['gini_index']:.3f}")

    plt.figure(figsize=(8, 5))
    sns.barplot(
        data=gini_data,
        x="library",
        y="gini_index",
        palette="viridis",
        order=LIB_CONFIG.library_order,
    )
    plt.title("Gini Index of Oligo Counts by Library\n(0 = equal, 1 = unequal)")
    plt.xlabel("Library")
    plt.ylabel("Gini Index")
    plt.ylim(0, 1)
    plt.tight_layout()
    outfile = env_config["plots_dir"] / "gini_index_by_library.png"
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"  Saved: {outfile.name}")


def plot_category_percentages(df: pd.DataFrame, env_config: dict) -> None:
    """
    Plot percentage of total counts by perfect/imperfect category for each library.

    Args:
        df: DataFrame with match quality columns
        env_config: Configuration dictionary with plots_dir
    """
    required_cols = {"counts", "library_ID", "perfect_guide", "perfect_donor",
                     "perfect_middle_donor_region"}
    if not required_cols.issubset(df.columns):
        print("  Skipping category plot: missing required columns")
        return

    # Filter low-count entries
    filtered_df = df[df["counts"] > 2]

    if filtered_df.empty:
        print("  Skipping category plot: no data after filtering")
        return

    # Aggregate by category
    summary = (
        filtered_df
        .groupby(["library_ID", "perfect_guide", "perfect_donor", "perfect_middle_donor_region"],
                 as_index=False)
        .agg({"counts": "sum"})
    )

    # Calculate percentages
    summary["percent_of_total"] = (
        summary.groupby("library_ID")["counts"]
        .transform(lambda x: 100 * x / x.sum())
    )

    # Create label
    summary["label"] = (
        summary["perfect_guide"].astype(str) + "|" +
        summary["perfect_donor"].astype(str) + "|" +
        summary["perfect_middle_donor_region"].astype(str)
    )

    plt.figure(figsize=(12, 6))
    sns.barplot(data=summary, x="label", y="percent_of_total", hue="library_ID")
    plt.ylabel("Percent of Total Counts")
    plt.xlabel("Guide|Donor|MiddleDonorRegion (True=Perfect)")
    plt.title("Percent of Total Counts by Match Category")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(ticks=range(0, 105, 10))
    plt.grid(axis="y", linestyle=":", linewidth=0.8)
    plt.legend(title="Library", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    outfile = env_config["plots_dir"] / "percent_of_total_counts_by_category_all_libraries.png"
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outfile.name}")


def main():
    """Main analysis pipeline."""
    # Load configuration from environment
    env_config = get_config_from_env()
    purity_dir = env_config["purity_dir"]
    plots_dir = env_config["plots_dir"]

    start_time = datetime.now()
    print("=" * 70)
    print("Step 06: Plot Plasmid Library Representation")
    print(f"Started: {start_time}")
    print(f"Project: {env_config['project_dir']}")
    print("=" * 70)

    # Setup
    setup_directories(env_config)

    # Check input directory
    if not purity_dir.exists():
        print(f"ERROR: Input directory not found: {purity_dir}")
        sys.exit(1)

    # Load and combine data files
    print("\nLoading purity-filtered files...")
    combined_df = load_data_files(purity_dir, INPUT_SUFFIX)
    print(f"  Total rows: {len(combined_df):,}")

    # Calculate GC content
    combined_df["gc_content"] = combined_df["donor"].apply(gc_content)

    # Check for non-string donors
    check_non_string_donors(combined_df, env_config)

    # Load designed oligos
    print("\nLoading designed oligos...")
    designed_oligos = load_designed_oligos(env_config)
    print(f"  Total designed oligos: {len(designed_oligos):,}")

    # Group by library and oligo_name
    print("\nAggregating by oligo name...")
    grouped_df = (
        combined_df
        .groupby(["library", "oligo_name"], as_index=False)["grouped_guide_bc0_counts"]
        .sum()
    )

    # Merge with designed oligos
    grouped_df = grouped_df.merge(
        designed_oligos[["oligo_name", "oligo_seq", "designed_donor", "gc_content", "SPS"]],
        on="oligo_name",
        how="left",
    )

    # Add missing oligos
    print("\nAdding missing oligos...")
    grouped_df = add_missing_oligos(grouped_df, designed_oligos)

    # Save missing oligos report
    missing_df = grouped_df[grouped_df["grouped_guide_bc0_counts"] == 0]
    if not missing_df.empty:
        missing_file = plots_dir / "missing_oligos_with_zero_counts.tsv"
        missing_df.to_csv(missing_file, sep="\t", index=False)
        print(f"  Saved: {missing_file.name} ({len(missing_df):,} oligos)")

    # Rename and process
    grouped_df.rename(
        columns={"grouped_guide_bc0_counts": "grouped_oligo_name_counts"},
        inplace=True,
    )

    # Extract SPS and filter to top SPS per library
    grouped_df["SPS"] = grouped_df["oligo_seq"].str[-20:]

    sps_totals = (
        grouped_df
        .groupby(["library", "SPS"], as_index=False)["grouped_oligo_name_counts"]
        .sum()
        .sort_values(["library", "grouped_oligo_name_counts"], ascending=[True, False])
        .drop_duplicates(subset=["library"])
    )

    filtered_df = grouped_df.merge(
        sps_totals[["library", "SPS"]],
        on=["library", "SPS"],
        how="inner",
    )

    # Generate plots
    print("\nGenerating plots...")

    # GC vs counts
    print("\n  Creating GC content vs counts plot...")
    plot_gc_vs_counts(filtered_df, env_config)

    # Normalize counts for density plot
    filtered_df["normalized_counts"] = (
        filtered_df
        .groupby("library")["grouped_oligo_name_counts"]
        .transform(lambda x: x / x.sum() * x.count())
    )

    # Density plot
    print("  Creating density plots...")
    plot_density_by_sps(filtered_df, env_config)

    # Gini indices
    print("  Creating Gini index plot...")
    plot_gini_indices(filtered_df, env_config)

    # Category percentages (if columns available)
    print("  Creating category breakdown plot...")
    plot_category_percentages(combined_df, env_config)

    # Save summary files
    print("\nSaving summary files...")
    summary_file = plots_dir / "top_1000_rows_combined_df.tsv"
    combined_df.head(1000).to_csv(summary_file, sep="\t", index=False)
    print(f"  Saved: {summary_file.name}")

    if "perfect_guide" in combined_df.columns:
        not_perfect_file = plots_dir / "top_1000_rows_combined_df_not_perfect_guide.tsv"
        combined_df.query("not perfect_guide").head(1000).to_csv(
            not_perfect_file, sep="\t", index=False
        )
        print(f"  Saved: {not_perfect_file.name}")

    # Final summary
    elapsed = datetime.now() - start_time
    print(f"\n{'=' * 70}")
    print("Step F Complete")
    print(f"  Output directory: {plots_dir}")
    print(f"  Elapsed time: {elapsed}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
