"""
Steps 05 and 06 core logic: Library representation assessment and plotting.

Aggregates purity-filtered guide-donor-BC0 counts to the oligo-name level,
identifies missing oligos, filters to the dominant SPS subpool per library,
and produces publication-ready QC plots.

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


# ============================================================================
# Oligo structure constants (200mer design)
# ============================================================================

@dataclass(frozen=True)
class OligoConstants:
    """Positional constants for extracting features from designed 200mer oligos."""
    DONOR_START: int = 51
    DONOR_END: int = 180
    SPS_START: int = 180
    SPS_LENGTH: int = 20


OLIGO = OligoConstants()

PURITY_INPUT_SUFFIX = "_guide_donor_bc0_purity_filtered.tsv"


# ============================================================================
# Utility functions
# ============================================================================

def gc_content(seq: str) -> float:
    """
    Return the GC fraction of a DNA sequence (0.0–1.0).

    Returns 0.0 for empty or non-string inputs.
    """
    if not isinstance(seq, str) or len(seq) == 0:
        return 0.0
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)


def gini_index(counts: np.ndarray) -> float:
    """
    Return the Gini index for an array of counts (0 = equal, 1 = maximally unequal).

    Zero-count entries are included in the calculation so that missing oligos
    drag the index toward 1 (unequal).

    Parameters
    ----------
    counts : array-like
        Non-negative count values.

    Returns
    -------
    float
        Gini index, or ``np.nan`` if input is empty or all zeros.
    """
    counts = np.array(counts)
    n = len(counts)
    if n == 0 or np.sum(counts) == 0:
        return np.nan
    sorted_counts = np.sort(counts)
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * sorted_counts) / np.sum(counts) - (n + 1)) / n


# ============================================================================
# Data loading
# ============================================================================

def load_purity_files(
    input_dir: Path,
    suffix: str = PURITY_INPUT_SUFFIX,
) -> List[pd.DataFrame]:
    """
    Load all purity-filtered TSVs from *input_dir* and attach a ``library`` column.

    Parameters
    ----------
    input_dir : Path
        Directory containing purity-filtered files from Step 04.
    suffix : str
        File suffix pattern.

    Returns
    -------
    list[pd.DataFrame]
        One DataFrame per file, with a ``library`` column added.
    """
    df_list = []
    for filepath in sorted(input_dir.glob(f"*{suffix}")):
        library = filepath.stem.replace(suffix.replace(".tsv", ""), "")
        df = pd.read_csv(filepath, sep="\t")
        df["library"] = library
        df_list.append(df)
    return df_list


def load_designed_oligos(oligo_pool_files: List[Path]) -> pd.DataFrame:
    """
    Load oligo pool TSVs and extract donor, SPS, and GC content features.

    Parameters
    ----------
    oligo_pool_files : list[Path]
        Paths to Twist oligo pool TSV files (must have ``oligo_seq`` and
        ``oligo_name`` columns).

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with additional columns: ``designed_donor``,
        ``SPS``, ``gc_content``.

    Raises
    ------
    FileNotFoundError
        If none of the provided files exist.
    """
    dfs = []
    for filepath in oligo_pool_files:
        if filepath.exists():
            dfs.append(pd.read_csv(filepath, sep="\t"))

    if not dfs:
        raise FileNotFoundError(f"No oligo pool files found in: {oligo_pool_files}")

    combined = pd.concat(dfs, ignore_index=True)
    combined["oligo_seq"] = combined["oligo_seq"].str.upper()
    combined["designed_donor"] = combined["oligo_seq"].str[OLIGO.DONOR_START:OLIGO.DONOR_END]
    combined["SPS"] = combined["oligo_seq"].str[OLIGO.SPS_START:]
    combined["gc_content"] = combined["designed_donor"].apply(gc_content)
    return combined


# ============================================================================
# Step 05: Representation analysis
# ============================================================================

def aggregate_by_oligo_name(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sum guide-BC0 counts and count unique BC0 barcodes at the oligo-name level.

    Parameters
    ----------
    df : pd.DataFrame
        Must have columns: ``library``, ``oligo_name``,
        ``grouped_guide_bc0_counts``, ``bc0``.

    Returns
    -------
    pd.DataFrame
        Aggregated with columns: ``library``, ``oligo_name``,
        ``grouped_guide_bc0_counts``, ``unique_bc0``.
    """
    required = {"library", "oligo_name", "grouped_guide_bc0_counts", "bc0"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    return df.groupby(["library", "oligo_name"], as_index=False).agg(
        grouped_guide_bc0_counts=("grouped_guide_bc0_counts", "sum"),
        unique_bc0=("bc0", "nunique"),
    )


def add_missing_oligos(
    df: pd.DataFrame,
    designed_oligos: pd.DataFrame,
    library_sps_map: Dict[str, str],
) -> pd.DataFrame:
    """
    Append zero-count rows for every designed oligo absent from *df*.

    Only oligos whose SPS matches the expected SPS for a given library are
    considered "missing" — cross-pool oligos are not flagged.

    Parameters
    ----------
    df : pd.DataFrame
        Aggregated oligo-level counts (output of :func:`aggregate_by_oligo_name`
        after merging with designed_oligos for SPS/seq columns).
    designed_oligos : pd.DataFrame
        Full designed oligo pool with ``oligo_name``, ``SPS``, ``oligo_seq``,
        ``designed_donor``, ``gc_content`` columns.
    library_sps_map : dict[str, str]
        Maps each library ID to its expected SPS sequence.

    Returns
    -------
    pd.DataFrame
        Extended DataFrame with missing oligos appended as zero-count rows.
    """
    parts = [df]

    for library, sps in library_sps_map.items():
        expected = designed_oligos[designed_oligos["SPS"] == sps]
        existing_names = df.loc[df["library"] == library, "oligo_name"]
        missing = expected[~expected["oligo_name"].isin(existing_names)]

        if missing.empty:
            continue

        # Skip oligos that share a sequence with an existing row (renaming artifact)
        if "oligo_seq" in df.columns:
            existing_seqs = df.loc[df["library"] == library, "oligo_seq"]
            missing = missing[~missing["oligo_seq"].isin(existing_seqs)]

        if missing.empty:
            continue

        to_append = missing.assign(
            library=library,
            grouped_guide_bc0_counts=0,
            unique_bc0=0,
        )[["library", "oligo_name", "grouped_guide_bc0_counts", "unique_bc0",
           "oligo_seq", "designed_donor", "gc_content", "SPS"]]
        parts.append(to_append)

    return pd.concat(parts, ignore_index=True)


def filter_by_top_sps(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only the dominant SPS subpool per library (highest total count).

    Handles libraries that received off-target synthesis pool contamination by
    retaining only the most abundant SPS.

    The input DataFrame must have columns ``library``, ``SPS``, and
    ``grouped_oligo_name_counts``.
    """
    sps_totals = (
        df.groupby(["library", "SPS"], as_index=False)
        .agg(total_counts=("grouped_oligo_name_counts", "sum"))
        .sort_values(["library", "total_counts"], ascending=[True, False])
    )
    top_sps = sps_totals.drop_duplicates(subset=["library"], keep="first")

    return df.merge(top_sps[["library", "SPS"]], on=["library", "SPS"], how="inner")


def save_representation_outputs(df: pd.DataFrame, summary_dir: Path) -> None:
    """
    Write per-library and combined oligo-representation TSVs.

    Files written:
    - ``summary_dir/all_libraries_oligo_name_filtered.tsv``
    - ``summary_dir/{library}_oligo_name_filtered.tsv``  (one per library)
    - ``summary_dir/missing_oligos_with_zero_counts.tsv``  (if any)

    Parameters
    ----------
    df : pd.DataFrame
        Filtered oligo-level representation DataFrame.
    summary_dir : Path
        Output directory (created if absent).
    """
    summary_dir.mkdir(parents=True, exist_ok=True)

    df.to_csv(summary_dir / "all_libraries_oligo_name_filtered.tsv", sep="\t", index=False)

    for library in df["library"].unique():
        lib_df = df[df["library"] == library]
        lib_df.to_csv(
            summary_dir / f"{library}_oligo_name_filtered.tsv",
            sep="\t", index=False,
        )

    missing = df[df["grouped_oligo_name_counts"] == 0]
    if not missing.empty:
        missing.to_csv(
            summary_dir / "missing_oligos_with_zero_counts.tsv",
            sep="\t", index=False,
        )


def run_representation_analysis(
    purity_dir: Path,
    oligo_pool_files: List[Path],
    summary_dir: Path,
    library_sps_map: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    End-to-end Step 05 analysis: load → aggregate → fill missing → filter → save.

    Parameters
    ----------
    purity_dir : Path
        Directory containing Step 04 purity-filtered TSV files.
    oligo_pool_files : list[Path]
        Paths to Twist oligo pool design files.
    summary_dir : Path
        Output directory for TSV files.
    library_sps_map : dict[str, str], optional
        Maps library IDs to expected SPS sequences.  If ``None``, missing-oligo
        filling and SPS filtering are skipped.

    Returns
    -------
    pd.DataFrame
        The final filtered representation DataFrame (also written to disk).
    """
    # Load and combine purity-filtered files
    df_list = load_purity_files(purity_dir)
    if not df_list:
        raise FileNotFoundError(f"No purity-filtered files found in {purity_dir}")

    combined = pd.concat(df_list, ignore_index=True)
    combined["gc_content"] = combined["donor"].apply(gc_content) if "donor" in combined.columns else 0.0

    # Load designed oligos
    designed_oligos = load_designed_oligos(oligo_pool_files)

    # Aggregate to oligo-name level
    oligo_df = aggregate_by_oligo_name(combined)

    # Merge with designed oligo annotations
    oligo_df = oligo_df.merge(
        designed_oligos[["oligo_name", "oligo_seq", "designed_donor", "gc_content", "SPS"]],
        on="oligo_name",
        how="left",
    )
    for col in ("SPS", "oligo_seq"):
        if col in oligo_df.columns:
            oligo_df[col] = oligo_df[col].str.upper()

    # Fill missing oligos
    if library_sps_map:
        oligo_df = add_missing_oligos(oligo_df, designed_oligos, library_sps_map)

    # Rename for downstream consistency
    oligo_df.rename(
        columns={"grouped_guide_bc0_counts": "grouped_oligo_name_counts"},
        inplace=True,
    )

    # Filter to dominant SPS per library (only when SPS map is known)
    if library_sps_map and "grouped_oligo_name_counts" in oligo_df.columns and "SPS" in oligo_df.columns:
        oligo_df = filter_by_top_sps(oligo_df)

    save_representation_outputs(oligo_df, summary_dir)
    return oligo_df


# ============================================================================
# Step 06: Plotting
# ============================================================================

def plot_gc_vs_counts(
    df: pd.DataFrame,
    plots_dir: Path,
    library_order: Optional[List[str]] = None,
) -> None:
    """
    Scatter plot of designed GC content vs total oligo counts, faceted by library.

    Parameters
    ----------
    df : pd.DataFrame
        Must have columns: ``gc_content``, ``grouped_oligo_name_counts``, ``library``.
    plots_dir : Path
        Directory for output PNG files.
    library_order : list[str], optional
        Display order of libraries.
    """
    import matplotlib.pyplot as plt
    import scipy.stats as stats
    import seaborn as sns

    sns.set_style("whitegrid")
    if library_order is None:
        library_order = sorted(df["library"].unique())

    plot_df = df[df["library"].isin(library_order)]
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
        col_order=library_order,
    )
    g.set_axis_labels("Designed GC Content", "Total Counts")
    g.set(yscale="log")
    g.fig.suptitle("GC Content vs Total Counts by Library (log scale)", y=1.02)

    for library in library_order:
        lib_df = plot_df[plot_df["library"] == library]
        valid = lib_df[
            (lib_df["grouped_oligo_name_counts"] > 0) & lib_df["gc_content"].notnull()
        ]
        if len(valid) > 2:
            r, p = stats.pearsonr(
                valid["gc_content"], np.log10(valid["grouped_oligo_name_counts"])
            )
            ax_idx = library_order.index(library)
            ax = g.axes.flat[ax_idx]
            ax.text(
                0.05, 0.95,
                f"$R^2$={r**2:.2f}\np={p:.2g}",
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=10,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
            )

    plt.tight_layout()
    outfile = plots_dir / "all_libraries_gc_content_vs_counts_faceted_logy.png"
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()


def plot_density_by_sps(
    df: pd.DataFrame,
    plots_dir: Path,
    library_order: Optional[List[str]] = None,
) -> None:
    """
    Kernel density plot of normalised oligo counts, split by SPS subpool.

    Parameters
    ----------
    df : pd.DataFrame
        Must have columns: ``normalized_counts``, ``library``, ``SPS``.
    plots_dir : Path
        Directory for output PNG files.
    library_order : list[str], optional
        Display order of libraries.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    if library_order is None:
        library_order = sorted(df["library"].unique())

    unique_sps = df["SPS"].dropna().unique()
    palette = sns.color_palette("tab10", n_colors=len(library_order))
    library_palette = dict(zip(library_order, palette))

    fig, axes = plt.subplots(
        1, len(unique_sps),
        figsize=(5 * len(unique_sps), 6),
        sharey=True,
    )
    if len(unique_sps) == 1:
        axes = [axes]

    for idx, sps in enumerate(unique_sps):
        sps_df = df[df["SPS"] == sps]
        ax = axes[idx]
        valid = sps_df[sps_df["normalized_counts"] > 0]

        if valid.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"SPS: {sps[:10]}...")
            continue

        sns.kdeplot(
            data=valid,
            x="normalized_counts",
            hue="library",
            common_norm=False,
            fill=True,
            alpha=0.5,
            palette=library_palette,
            bw_adjust=0.5,
            hue_order=library_order,
            ax=ax,
        )
        ax.set_xscale("log")
        ax.set_xlabel("Normalized Counts (log scale)")
        ax.set_ylabel("Density" if idx == 0 else "")
        ax.set_title(f"SPS: {sps[:10]}...")

        present = [lib for lib in library_order if lib in valid["library"].unique()]
        handles = [
            plt.Line2D([0], [0], color=library_palette[lib], lw=4)
            for lib in present
        ]
        ax.legend(handles, present, title="Library",
                  bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.grid(axis="y", linestyle=":", linewidth=0.8)

    plt.suptitle("Density Plot of Normalized Counts by Library (by SPS)", y=1.02)
    plt.tight_layout()
    outfile = plots_dir / "normalized_counts_density_plot_by_SPS_side_by_side.png"
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()


def plot_gini_indices(
    df: pd.DataFrame,
    plots_dir: Path,
    library_order: Optional[List[str]] = None,
) -> None:
    """
    Bar plot of Gini index per library (lower = more uniform representation).

    Parameters
    ----------
    df : pd.DataFrame
        Must have columns: ``library``, ``grouped_oligo_name_counts``.
    plots_dir : Path
        Directory for output PNG files.
    library_order : list[str], optional
        Display order of libraries.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    if library_order is None:
        library_order = sorted(df["library"].unique())

    gini_data = (
        df.groupby("library")["grouped_oligo_name_counts"]
        .apply(lambda x: gini_index(x[x > 0].values))
        .reset_index(name="gini_index")
    )

    plt.figure(figsize=(8, 5))
    sns.barplot(
        data=gini_data,
        x="library",
        y="gini_index",
        palette="viridis",
        order=library_order,
    )
    plt.title("Gini Index of Oligo Counts by Library\n(0 = equal, 1 = unequal)")
    plt.xlabel("Library")
    plt.ylabel("Gini Index")
    plt.ylim(0, 1)
    plt.tight_layout()
    outfile = plots_dir / "gini_index_by_library.png"
    plt.savefig(outfile, dpi=150)
    plt.close()


def plot_category_percentages(
    df: pd.DataFrame,
    plots_dir: Path,
) -> None:
    """
    Stacked bar plot: percent of total counts by perfect/imperfect match category.

    Skipped silently if the required columns are absent.

    Parameters
    ----------
    df : pd.DataFrame
        Must have columns: ``counts``, ``library_ID``, ``perfect_guide``,
        ``perfect_donor``, ``perfect_middle_donor_region``.
    plots_dir : Path
        Directory for output PNG files.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    required = {"counts", "library_ID", "perfect_guide", "perfect_donor",
                "perfect_middle_donor_region"}
    if not required.issubset(df.columns):
        return

    filtered = df[df["counts"] > 2]
    if filtered.empty:
        return

    summary = (
        filtered
        .groupby(
            ["library_ID", "perfect_guide", "perfect_donor", "perfect_middle_donor_region"],
            as_index=False,
        )
        .agg({"counts": "sum"})
    )
    summary["percent_of_total"] = (
        summary.groupby("library_ID")["counts"]
        .transform(lambda x: 100 * x / x.sum())
    )
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
    outfile = plots_dir / "percent_of_total_counts_by_category_all_libraries.png"
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()


def run_library_plots(
    purity_dir: Path,
    oligo_pool_files: List[Path],
    plots_dir: Path,
    library_sps_map: Optional[Dict[str, str]] = None,
    library_order: Optional[List[str]] = None,
) -> None:
    """
    End-to-end Step 06: load purity files → group → plot.

    Parameters
    ----------
    purity_dir : Path
        Directory containing Step 04 purity-filtered TSV files.
    oligo_pool_files : list[Path]
        Paths to Twist oligo pool design files.
    plots_dir : Path
        Directory for output PNG files.
    library_sps_map : dict[str, str], optional
        Maps library IDs to expected SPS sequences for missing-oligo filling.
    library_order : list[str], optional
        Display order of libraries in plots.
    """
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Load purity files
    df_list = load_purity_files(purity_dir)
    if not df_list:
        raise FileNotFoundError(f"No purity-filtered files found in {purity_dir}")

    combined = pd.concat(df_list, ignore_index=True)

    if "donor" in combined.columns:
        combined["gc_content"] = combined["donor"].apply(gc_content)

    # Load designed oligos
    designed_oligos = load_designed_oligos(oligo_pool_files)

    # Aggregate to oligo name level
    grouped = (
        combined
        .groupby(["library", "oligo_name"], as_index=False)["grouped_guide_bc0_counts"]
        .sum()
    )
    grouped = grouped.merge(
        designed_oligos[["oligo_name", "oligo_seq", "designed_donor", "gc_content", "SPS"]],
        on="oligo_name",
        how="left",
    )

    # Add missing oligos
    if library_sps_map:
        grouped = add_missing_oligos(grouped, designed_oligos, library_sps_map)

    grouped.rename(
        columns={"grouped_guide_bc0_counts": "grouped_oligo_name_counts"},
        inplace=True,
    )

    # Save missing oligos report
    missing = grouped[grouped["grouped_oligo_name_counts"] == 0]
    if not missing.empty:
        missing.to_csv(plots_dir / "missing_oligos_with_zero_counts.tsv", sep="\t", index=False)

    # Filter to top SPS per library
    if "SPS" in grouped.columns:
        grouped = filter_by_top_sps(grouped)

    if library_order is None:
        library_order = sorted(grouped["library"].unique())

    # Plot: GC vs counts
    plot_gc_vs_counts(grouped, plots_dir, library_order)

    # Normalize for density plot
    grouped["normalized_counts"] = (
        grouped
        .groupby("library")["grouped_oligo_name_counts"]
        .transform(lambda x: x / x.sum() * x.count())
    )

    # Plot: density by SPS
    plot_density_by_sps(grouped, plots_dir, library_order)

    # Plot: Gini indices
    plot_gini_indices(grouped, plots_dir, library_order)

    # Plot: category percentages (only if match quality columns present)
    plot_category_percentages(combined, plots_dir)
