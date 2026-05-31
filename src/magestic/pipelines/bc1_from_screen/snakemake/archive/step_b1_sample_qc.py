#!/usr/bin/env python3
"""
Step B1: Sample Quality Control for BC1 Screen Data

Run AFTER Step B (bc1 counting) and BEFORE Step C (tier classification).

This script identifies problematic samples (failed PCR, poor gDNA extraction,
insufficient reads) before they propagate errors into DESeq2 analysis.

QC Metrics computed:
- total_reads: Sum of all read counts
- n_unique_bc1s: Number of unique BC1s detected
- shannon_entropy: Diversity metric (bits)
- normalized_entropy: Entropy / log2(n_bc1s)
- match_rate: Fraction of reads with BC1 match
- plate_type: "FWD" or "REV" (based on inner primer pair)
- mean_replicate_correlation: Mean Pearson r with same-plate-type replicates

QC Flags:
- LOW_READ_COUNT: total_reads < 50,000
- LOW_BC1_DIVERSITY: normalized_entropy < 0.5
- POOR_REPLICATE_CONCORDANCE: mean_correlation < 0.85 (within plate type)
- PCA_OUTLIER: >3 SD from condition cluster center

Usage:
    python step_b1_sample_qc.py

Output:
    qc/sample_qc_summary.tsv      - All metrics per sample
    qc/flagged_samples.tsv        - Samples failing QC
    qc_plots/*.png                - QC visualizations
"""

from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
import logging
import math
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)

# ============================================================================
# Configuration
# ============================================================================

PROJECT = "QTL"
SCREEN_NAME = "20250912_SpG_QTL_screen"

BASE_DIR = Path("/path/to")
PROJECT_DIR = BASE_DIR / "projects" / PROJECT / SCREEN_NAME
PROCESSED_DATA_DIR = (
    BASE_DIR / "processed_data/by_project" / PROJECT / SCREEN_NAME / "bc1_from_screen"
)

# Input files
LONG_TABLE_FILE = PROCESSED_DATA_DIR / "bc1_counts" / f"{SCREEN_NAME}_bc1_counts_long_format.tsv"
MATCH_SUMMARY_FILE = PROCESSED_DATA_DIR / "bc1_counts" / f"{SCREEN_NAME}_match_counts_summary.tsv"
KEYFILE = PROJECT_DIR / "keyfiles/bc1_from_screen/combined_key.tsv"

# Output directories
QC_DIR = PROCESSED_DATA_DIR / "qc"
QC_PLOTS_DIR = PROCESSED_DATA_DIR / "qc_plots"

# Primer pair to plate type mapping
# KR1964-KR1965 = plate_1B = REV (~98% REV reads)
# KR1966-KR1967 = plate_2C = FWD (~99% FWD reads)
PRIMER_TO_PLATE_TYPE = {
    "KR1964-KR1965": "REV",
    "KR1966-KR1967": "FWD",
    "KR2469-KR2472": "REV",
    "KR2470-KR2473": "FWD",
}


@dataclass
class SampleQCConfig:
    """Configuration for sample QC thresholds."""

    # Read count thresholds
    min_total_reads: int = 50_000           # Fail threshold
    low_read_warning: int = 100_000         # Warning threshold

    # Diversity thresholds
    min_unique_bc1s: int = 1000
    min_normalized_entropy: float = 0.5

    # Match rate
    min_match_rate: float = 0.95

    # Replicate concordance
    min_replicate_correlation: float = 0.85

    # PCA outlier
    pca_outlier_threshold: float = 3.0      # Z-score

    # Processing
    chunksize: int = 1_000_000              # Rows per chunk for large file
    n_top_bc1s_for_pca: int = 5000          # Most variable BC1s for PCA

    # Library-specific thresholds
    library_min_bc1s: Dict[str, int] = field(default_factory=lambda: {
        "yL437": 30_000,
        "yL442": 100_000,
    })


@dataclass
class SampleMetrics:
    """QC metrics for a single sample."""
    sample_name: str
    total_reads: int = 0
    n_unique_bc1s: int = 0
    shannon_entropy: float = 0.0
    normalized_entropy: float = 0.0
    simpson_index: float = 0.0
    plate_type: str = ""
    mean_replicate_correlation: float = 0.0
    min_replicate_correlation: float = 0.0
    pca_distance: float = 0.0

    # Metadata (populated from keyfile)
    library_ID: str = ""
    condition: str = ""
    generations: int = 0
    screen_date: str = ""
    screen_replicate: str = ""
    inner_primers: str = ""


# ============================================================================
# Metric Computation Functions
# ============================================================================


def compute_diversity_metrics(bc1_counts: Dict[str, int]) -> Tuple[float, float, float]:
    """
    Compute diversity metrics for a sample's BC1 counts.

    Args:
        bc1_counts: Dict mapping bc1 -> count

    Returns:
        Tuple of (shannon_entropy, normalized_entropy, simpson_index)
    """
    if not bc1_counts:
        return 0.0, 0.0, 0.0

    total = sum(bc1_counts.values())
    if total == 0:
        return 0.0, 0.0, 0.0

    n_bc1s = len(bc1_counts)

    # Calculate probabilities
    probs = [count / total for count in bc1_counts.values()]

    # Shannon entropy: H = -sum(p * log2(p))
    shannon_entropy = 0.0
    for p in probs:
        if p > 0:
            shannon_entropy -= p * math.log2(p)

    # Normalized entropy: H / H_max where H_max = log2(n)
    max_entropy = math.log2(n_bc1s) if n_bc1s > 1 else 1.0
    normalized_entropy = shannon_entropy / max_entropy if max_entropy > 0 else 0.0

    # Simpson index: 1 - sum(p^2)
    simpson_index = 1.0 - sum(p * p for p in probs)

    return shannon_entropy, normalized_entropy, simpson_index


def compute_sample_metrics_chunked(
    long_table_path: Path,
    keyfile: pd.DataFrame,
    config: SampleQCConfig,
) -> pd.DataFrame:
    """
    Compute QC metrics for all samples using chunked file reading.

    Args:
        long_table_path: Path to BC1 counts long format TSV
        keyfile: DataFrame with sample metadata
        config: QC configuration

    Returns:
        DataFrame with one row per sample containing all QC metrics
    """
    print(f"Computing sample metrics from {long_table_path}")
    print(f"  Chunk size: {config.chunksize:,} rows")

    # Build sample metadata lookup
    sample_metadata = {}
    for _, row in keyfile.iterrows():
        sample_name = row["sample_name"]
        sample_metadata[sample_name] = {
            "library_ID": row.get("library_ID", ""),
            "condition": row.get("condition", ""),
            "generations": row.get("generations", 0),
            "screen_date": str(row.get("screen_date", "")),
            "screen_replicate": row.get("screen_replicate", ""),
            "inner_primers": row.get("inner_primers", ""),
        }

    # Initialize per-sample aggregates
    sample_stats = defaultdict(lambda: {
        "total_reads": 0,
        "bc1_counts": defaultdict(int),
    })

    # Process file in chunks
    total_rows = 0
    for chunk_num, chunk in enumerate(pd.read_csv(
        long_table_path, sep="\t", chunksize=config.chunksize
    )):
        total_rows += len(chunk)
        if chunk_num % 10 == 0:
            print(f"  Processed {total_rows:,} rows...")

        # Aggregate counts per sample
        for _, row in chunk.iterrows():
            sample = row["sample_name"]
            bc1 = row["bc1"]
            count = row["sample_counts"]

            sample_stats[sample]["total_reads"] += count
            sample_stats[sample]["bc1_counts"][bc1] += count

    print(f"  Total rows processed: {total_rows:,}")
    print(f"  Unique samples found: {len(sample_stats):,}")

    # Compute final metrics for each sample
    metrics_list = []
    for sample_name, stats in sample_stats.items():
        # Diversity metrics
        shannon, normalized, simpson = compute_diversity_metrics(stats["bc1_counts"])

        # Get metadata
        meta = sample_metadata.get(sample_name, {})
        inner_primers = meta.get("inner_primers", "")
        plate_type = PRIMER_TO_PLATE_TYPE.get(inner_primers, "UNKNOWN")

        metrics = SampleMetrics(
            sample_name=sample_name,
            total_reads=stats["total_reads"],
            n_unique_bc1s=len(stats["bc1_counts"]),
            shannon_entropy=shannon,
            normalized_entropy=normalized,
            simpson_index=simpson,
            plate_type=plate_type,
            library_ID=meta.get("library_ID", ""),
            condition=meta.get("condition", ""),
            generations=meta.get("generations", 0),
            screen_date=meta.get("screen_date", ""),
            screen_replicate=meta.get("screen_replicate", ""),
            inner_primers=inner_primers,
        )
        metrics_list.append(metrics)

    # Convert to DataFrame
    df = pd.DataFrame([vars(m) for m in metrics_list])
    return df


def load_match_summary(match_summary_path: Path) -> pd.DataFrame:
    """Load match summary to get per-sample match rates."""
    if not match_summary_path.exists():
        print(f"  Warning: Match summary file not found: {match_summary_path}")
        return None
    return pd.read_csv(match_summary_path, sep="\t")


# ============================================================================
# Replicate Concordance
# ============================================================================


def compute_replicate_concordance(
    sample_metrics: pd.DataFrame,
    long_table_path: Path,
    keyfile: pd.DataFrame,
    config: SampleQCConfig,
) -> pd.DataFrame:
    """
    Compute pairwise correlations between replicates of same condition.

    Correlations are computed WITHIN plate type (FWD vs FWD, REV vs REV)
    since FWD and REV plates have different technical characteristics.

    Returns:
        DataFrame with sample_name, mean_replicate_correlation, min_replicate_correlation
    """
    print("Computing replicate concordance...")

    # Create condition+generations+library+plate_type grouping
    sample_metrics = sample_metrics.copy()
    sample_metrics["condition_group"] = (
        sample_metrics["condition"] + "_" +
        sample_metrics["generations"].astype(str) + "g_" +
        sample_metrics["library_ID"] + "_" +
        sample_metrics["plate_type"]
    )

    # Get unique conditions with multiple replicates
    group_counts = sample_metrics.groupby("condition_group").size()
    groups_with_replicates = group_counts[group_counts > 1].index.tolist()

    print(f"  Found {len(groups_with_replicates)} condition groups with replicates")

    if len(groups_with_replicates) == 0:
        # No replicates to compare
        sample_metrics["mean_replicate_correlation"] = np.nan
        sample_metrics["min_replicate_correlation"] = np.nan
        return sample_metrics[["sample_name", "mean_replicate_correlation", "min_replicate_correlation"]]

    # Build a wide matrix for correlation computation
    # For efficiency, we'll use a subset of top BC1s
    print("  Loading BC1 counts for correlation computation...")

    # First pass: identify most abundant BC1s across all samples
    bc1_totals = defaultdict(int)
    for chunk in pd.read_csv(long_table_path, sep="\t", chunksize=config.chunksize):
        for bc1, count in chunk.groupby("bc1")["sample_counts"].sum().items():
            bc1_totals[bc1] += count

    # Select top N BC1s
    top_bc1s = sorted(bc1_totals.keys(), key=lambda x: bc1_totals[x], reverse=True)
    top_bc1s = set(top_bc1s[:config.n_top_bc1s_for_pca])
    print(f"  Using top {len(top_bc1s):,} BC1s for correlation")

    # Second pass: build wide matrix for these BC1s
    sample_bc1_counts = defaultdict(lambda: defaultdict(int))
    for chunk in pd.read_csv(long_table_path, sep="\t", chunksize=config.chunksize):
        chunk_filtered = chunk[chunk["bc1"].isin(top_bc1s)]
        for _, row in chunk_filtered.iterrows():
            sample_bc1_counts[row["sample_name"]][row["bc1"]] += row["sample_counts"]

    # Convert to matrix
    sample_names = list(sample_bc1_counts.keys())
    bc1_list = list(top_bc1s)
    matrix = np.zeros((len(sample_names), len(bc1_list)))
    for i, sample in enumerate(sample_names):
        for j, bc1 in enumerate(bc1_list):
            matrix[i, j] = sample_bc1_counts[sample][bc1]

    # Log transform (add pseudocount)
    matrix_log = np.log2(matrix + 1)

    # Compute pairwise correlations
    print("  Computing pairwise correlations...")
    sample_to_idx = {name: i for i, name in enumerate(sample_names)}

    correlation_results = {}
    for sample_name in sample_metrics["sample_name"]:
        if sample_name not in sample_to_idx:
            correlation_results[sample_name] = (np.nan, np.nan)
            continue

        # Get this sample's condition group
        sample_row = sample_metrics[sample_metrics["sample_name"] == sample_name].iloc[0]
        group = sample_row["condition_group"]

        # Get other samples in same group
        group_samples = sample_metrics[
            (sample_metrics["condition_group"] == group) &
            (sample_metrics["sample_name"] != sample_name)
        ]["sample_name"].tolist()

        if not group_samples:
            correlation_results[sample_name] = (np.nan, np.nan)
            continue

        # Compute correlations with each replicate
        sample_idx = sample_to_idx[sample_name]
        sample_vec = matrix_log[sample_idx]

        correlations = []
        for other_name in group_samples:
            if other_name in sample_to_idx:
                other_idx = sample_to_idx[other_name]
                other_vec = matrix_log[other_idx]
                r, _ = stats.pearsonr(sample_vec, other_vec)
                correlations.append(r)

        if correlations:
            correlation_results[sample_name] = (np.mean(correlations), np.min(correlations))
        else:
            correlation_results[sample_name] = (np.nan, np.nan)

    # Add to dataframe
    sample_metrics["mean_replicate_correlation"] = sample_metrics["sample_name"].map(
        lambda x: correlation_results.get(x, (np.nan, np.nan))[0]
    )
    sample_metrics["min_replicate_correlation"] = sample_metrics["sample_name"].map(
        lambda x: correlation_results.get(x, (np.nan, np.nan))[1]
    )

    return sample_metrics[["sample_name", "mean_replicate_correlation", "min_replicate_correlation"]]


# ============================================================================
# PCA Analysis
# ============================================================================


def compute_pca_and_outliers(
    sample_metrics: pd.DataFrame,
    long_table_path: Path,
    config: SampleQCConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute PCA on log-transformed counts and detect outliers.

    Returns:
        Tuple of (PCA coordinates DataFrame, DataFrame with pca_distance per sample)
    """
    print("Computing PCA...")

    # Build wide matrix (same as correlation computation)
    bc1_totals = defaultdict(int)
    for chunk in pd.read_csv(long_table_path, sep="\t", chunksize=config.chunksize):
        for bc1, count in chunk.groupby("bc1")["sample_counts"].sum().items():
            bc1_totals[bc1] += count

    # Select top N BC1s by variance (total counts as proxy)
    top_bc1s = sorted(bc1_totals.keys(), key=lambda x: bc1_totals[x], reverse=True)
    top_bc1s = set(top_bc1s[:config.n_top_bc1s_for_pca])

    # Build matrix
    sample_bc1_counts = defaultdict(lambda: defaultdict(int))
    for chunk in pd.read_csv(long_table_path, sep="\t", chunksize=config.chunksize):
        chunk_filtered = chunk[chunk["bc1"].isin(top_bc1s)]
        for _, row in chunk_filtered.iterrows():
            sample_bc1_counts[row["sample_name"]][row["bc1"]] += row["sample_counts"]

    sample_names = list(sample_bc1_counts.keys())
    bc1_list = list(top_bc1s)
    matrix = np.zeros((len(sample_names), len(bc1_list)))
    for i, sample in enumerate(sample_names):
        for j, bc1 in enumerate(bc1_list):
            matrix[i, j] = sample_bc1_counts[sample][bc1]

    # Log transform
    matrix_log = np.log2(matrix + 1)

    # Standardize
    scaler = StandardScaler()
    matrix_scaled = scaler.fit_transform(matrix_log)

    # PCA
    n_components = min(10, len(sample_names) - 1, len(bc1_list) - 1)
    pca = PCA(n_components=n_components)
    pca_coords = pca.fit_transform(matrix_scaled)

    # Create PCA DataFrame
    pca_df = pd.DataFrame(
        pca_coords,
        columns=[f"PC{i+1}" for i in range(n_components)],
        index=sample_names
    )
    pca_df["sample_name"] = sample_names

    print(f"  Explained variance: {pca.explained_variance_ratio_[:3].sum():.1%} (first 3 PCs)")

    # Compute outlier distances within condition groups
    pca_df = pca_df.merge(
        sample_metrics[["sample_name", "condition", "generations", "library_ID"]],
        on="sample_name",
        how="left"
    )
    pca_df["condition_group"] = (
        pca_df["condition"].fillna("") + "_" +
        pca_df["generations"].fillna(0).astype(int).astype(str) + "g_" +
        pca_df["library_ID"].fillna("")
    )

    # Compute distance from group centroid
    outlier_distances = {}
    for group, group_df in pca_df.groupby("condition_group"):
        if len(group_df) < 3:
            # Not enough samples to detect outliers
            for sample in group_df["sample_name"]:
                outlier_distances[sample] = 0.0
            continue

        # Use first 5 PCs for distance
        coords = group_df[[f"PC{i+1}" for i in range(min(5, n_components))]].values
        centroid = coords.mean(axis=0)

        # Compute distances
        for i, sample in enumerate(group_df["sample_name"]):
            dist = np.linalg.norm(coords[i] - centroid)
            # Normalize by group std
            group_std = np.std([np.linalg.norm(c - centroid) for c in coords])
            z_score = dist / group_std if group_std > 0 else 0.0
            outlier_distances[sample] = z_score

    # Map back to sample_metrics
    distance_df = pd.DataFrame([
        {"sample_name": name, "pca_distance": dist}
        for name, dist in outlier_distances.items()
    ])

    return pca_df, distance_df


# ============================================================================
# Outlier Detection
# ============================================================================


def detect_outliers(
    sample_metrics: pd.DataFrame,
    config: SampleQCConfig,
) -> pd.DataFrame:
    """
    Flag samples that fail QC thresholds.

    Returns:
        DataFrame with columns: sample_name, flag, reason, metric_value, threshold, severity
    """
    print("Detecting outliers...")

    flags = []

    for _, row in sample_metrics.iterrows():
        sample = row["sample_name"]

        # Check read count
        if row["total_reads"] < config.min_total_reads:
            severity = "fail" if row["total_reads"] < 10_000 else "warning"
            flags.append({
                "sample_name": sample,
                "flag": "LOW_READ_COUNT",
                "reason": f"total_reads={row['total_reads']:,} < {config.min_total_reads:,}",
                "metric_value": row["total_reads"],
                "threshold": config.min_total_reads,
                "severity": severity,
            })

        # Check BC1 diversity
        if row["normalized_entropy"] < config.min_normalized_entropy:
            flags.append({
                "sample_name": sample,
                "flag": "LOW_BC1_DIVERSITY",
                "reason": f"normalized_entropy={row['normalized_entropy']:.3f} < {config.min_normalized_entropy}",
                "metric_value": row["normalized_entropy"],
                "threshold": config.min_normalized_entropy,
                "severity": "warning",
            })

        # Check library-specific BC1 count
        library = row.get("library_ID", "")
        min_bc1s = config.library_min_bc1s.get(library, config.min_unique_bc1s)
        if row["n_unique_bc1s"] < min_bc1s:
            flags.append({
                "sample_name": sample,
                "flag": "LOW_BC1_COUNT",
                "reason": f"n_unique_bc1s={row['n_unique_bc1s']:,} < {min_bc1s:,} (library {library})",
                "metric_value": row["n_unique_bc1s"],
                "threshold": min_bc1s,
                "severity": "warning",
            })

        # Check replicate concordance
        mean_corr = row.get("mean_replicate_correlation", np.nan)
        if pd.notna(mean_corr) and mean_corr < config.min_replicate_correlation:
            flags.append({
                "sample_name": sample,
                "flag": "POOR_REPLICATE_CONCORDANCE",
                "reason": f"mean_correlation={mean_corr:.3f} < {config.min_replicate_correlation}",
                "metric_value": mean_corr,
                "threshold": config.min_replicate_correlation,
                "severity": "warning",
            })

        # Check PCA outlier
        pca_dist = row.get("pca_distance", 0.0)
        if pca_dist > config.pca_outlier_threshold:
            flags.append({
                "sample_name": sample,
                "flag": "PCA_OUTLIER",
                "reason": f"pca_distance={pca_dist:.2f} > {config.pca_outlier_threshold}",
                "metric_value": pca_dist,
                "threshold": config.pca_outlier_threshold,
                "severity": "warning",
            })

    flags_df = pd.DataFrame(flags)
    print(f"  Flagged {len(flags_df)} issues across {flags_df['sample_name'].nunique() if len(flags_df) > 0 else 0} samples")

    return flags_df


# ============================================================================
# Visualization
# ============================================================================


def plot_pca(
    pca_df: pd.DataFrame,
    sample_metrics: pd.DataFrame,
    flagged_samples: List[str],
    output_dir: Path,
) -> Dict[str, Path]:
    """
    Create PCA plots colored by condition, with shapes for library.
    """
    print("Generating PCA plots...")

    output_dir.mkdir(parents=True, exist_ok=True)
    plot_paths = {}

    # Merge metadata
    plot_df = pca_df.merge(
        sample_metrics[["sample_name", "condition", "library_ID", "plate_type", "generations"]],
        on="sample_name",
        how="left"
    )

    # Get unique conditions and assign colors
    conditions = plot_df["condition"].fillna("Unknown").unique()
    n_conditions = len(conditions)
    cmap = plt.cm.get_cmap("tab20", n_conditions)
    condition_colors = {cond: cmap(i) for i, cond in enumerate(conditions)}

    # Library markers
    library_markers = {"yL437": "o", "yL442": "s"}

    # Plate type fill
    plate_fills = {"FWD": True, "REV": False}

    # Figure 1: PCA by condition
    fig, ax = plt.subplots(figsize=(12, 10))

    for library in plot_df["library_ID"].unique():
        for plate_type in plot_df["plate_type"].unique():
            subset = plot_df[
                (plot_df["library_ID"] == library) &
                (plot_df["plate_type"] == plate_type)
            ]
            if len(subset) == 0:
                continue

            marker = library_markers.get(library, "o")
            fill = plate_fills.get(plate_type, True)

            for condition in subset["condition"].unique():
                cond_subset = subset[subset["condition"] == condition]
                color = condition_colors.get(condition, "gray")

                # Check if any are flagged
                mask_flagged = cond_subset["sample_name"].isin(flagged_samples)

                # Plot non-flagged
                if (~mask_flagged).any():
                    ax.scatter(
                        cond_subset.loc[~mask_flagged, "PC1"],
                        cond_subset.loc[~mask_flagged, "PC2"],
                        c=[color],
                        marker=marker,
                        s=60 if fill else 80,
                        facecolors=color if fill else "none",
                        edgecolors=color,
                        linewidth=1.5,
                        label=f"{condition} ({library}, {plate_type})" if fill else None,
                        alpha=0.8,
                    )

                # Plot flagged with X marker
                if mask_flagged.any():
                    ax.scatter(
                        cond_subset.loc[mask_flagged, "PC1"],
                        cond_subset.loc[mask_flagged, "PC2"],
                        c=["red"],
                        marker="X",
                        s=100,
                        alpha=0.9,
                        zorder=5,
                    )

    ax.set_xlabel("PC1", fontsize=12)
    ax.set_ylabel("PC2", fontsize=12)
    ax.set_title("Sample PCA - Colored by Condition\n(Shape=Library, Fill=FWD/Open=REV, X=Flagged)", fontsize=14)

    # Legend - just conditions
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8, label=cond)
               for cond, c in condition_colors.items()]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1, 0.5), fontsize=8)

    plt.tight_layout()
    pca_path = output_dir / "pca_by_condition.png"
    fig.savefig(pca_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    plot_paths["pca_by_condition"] = pca_path

    # Figure 2: PCA by library
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for ax, library in zip(axes, ["yL437", "yL442"]):
        subset = plot_df[plot_df["library_ID"] == library]
        if len(subset) == 0:
            ax.set_title(f"{library} - No data")
            continue

        for condition in subset["condition"].unique():
            cond_subset = subset[subset["condition"] == condition]
            color = condition_colors.get(condition, "gray")

            # Separate by plate type
            for plate_type in ["FWD", "REV"]:
                pt_subset = cond_subset[cond_subset["plate_type"] == plate_type]
                if len(pt_subset) == 0:
                    continue

                fill = plate_fills.get(plate_type, True)
                ax.scatter(
                    pt_subset["PC1"],
                    pt_subset["PC2"],
                    c=[color],
                    marker="o",
                    s=60,
                    facecolors=color if fill else "none",
                    edgecolors=color,
                    linewidth=1.5,
                    alpha=0.8,
                )

        ax.set_xlabel("PC1", fontsize=11)
        ax.set_ylabel("PC2", fontsize=11)
        ax.set_title(f"{library}", fontsize=13)

    # Shared legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8, label=cond)
               for cond, c in condition_colors.items()]
    fig.legend(handles=handles, loc="center right", fontsize=8)

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    pca_lib_path = output_dir / "pca_by_library.png"
    fig.savefig(pca_lib_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    plot_paths["pca_by_library"] = pca_lib_path

    return plot_paths


def plot_read_count_distribution(
    sample_metrics: pd.DataFrame,
    config: SampleQCConfig,
    output_dir: Path,
) -> Path:
    """Create boxplot of read counts by condition."""
    print("Generating read count distribution plot...")

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(14, 6))

    # Group by condition
    conditions = sorted(sample_metrics["condition"].unique())
    data = [sample_metrics[sample_metrics["condition"] == c]["total_reads"].values for c in conditions]

    bp = ax.boxplot(data, labels=conditions, patch_artist=True)

    # Color by library distribution
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor("lightblue")
        patch.set_alpha(0.7)

    # Add threshold line
    ax.axhline(y=config.min_total_reads, color="red", linestyle="--",
               label=f"Min threshold ({config.min_total_reads:,})")
    ax.axhline(y=config.low_read_warning, color="orange", linestyle=":",
               label=f"Warning threshold ({config.low_read_warning:,})")

    ax.set_ylabel("Total Reads", fontsize=12)
    ax.set_xlabel("Condition", fontsize=12)
    ax.set_title("Read Count Distribution by Condition", fontsize=14)
    ax.tick_params(axis='x', rotation=45)
    ax.legend(loc="upper right")

    # Format y-axis with commas
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: f"{int(x):,}"))

    plt.tight_layout()
    output_path = output_dir / "read_count_distribution.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path


def plot_clustering_dendrogram(
    pca_df: pd.DataFrame,
    sample_metrics: pd.DataFrame,
    output_dir: Path,
) -> Path:
    """Create hierarchical clustering dendrogram."""
    print("Generating clustering dendrogram...")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Use first 5 PCs
    n_pcs = min(5, len([c for c in pca_df.columns if c.startswith("PC")]))
    coords = pca_df[[f"PC{i+1}" for i in range(n_pcs)]].values

    # Compute linkage
    linkage_matrix = hierarchy.linkage(coords, method="ward")

    # Get condition colors
    conditions = sample_metrics.set_index("sample_name")["condition"]
    unique_conditions = conditions.unique()
    cmap = plt.cm.get_cmap("tab20", len(unique_conditions))
    condition_colors = {cond: cmap(i) for i, cond in enumerate(unique_conditions)}

    # Map to colors
    sample_names = pca_df["sample_name"].tolist()
    leaf_colors = [condition_colors.get(conditions.get(s, ""), "gray") for s in sample_names]

    fig, ax = plt.subplots(figsize=(20, 8))

    # Dendrogram
    dendro = hierarchy.dendrogram(
        linkage_matrix,
        labels=sample_names,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=4,
    )

    ax.set_title("Hierarchical Clustering of Samples", fontsize=14)
    ax.set_ylabel("Distance", fontsize=12)

    plt.tight_layout()
    output_path = output_dir / "clustering_dendrogram.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path


def plot_correlation_heatmap(
    sample_metrics: pd.DataFrame,
    output_dir: Path,
) -> Path:
    """Create heatmap of replicate correlations."""
    print("Generating correlation heatmap...")

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Pivot to show mean correlation by condition and plate type
    pivot_data = sample_metrics.pivot_table(
        values="mean_replicate_correlation",
        index="condition",
        columns=["library_ID", "plate_type"],
        aggfunc="mean"
    )

    if pivot_data.empty:
        ax.text(0.5, 0.5, "No correlation data available", ha="center", va="center")
    else:
        im = ax.imshow(pivot_data.values, cmap="RdYlGn", vmin=0.5, vmax=1.0, aspect="auto")

        ax.set_xticks(range(len(pivot_data.columns)))
        ax.set_xticklabels([f"{c[0]}\n{c[1]}" for c in pivot_data.columns], fontsize=8)
        ax.set_yticks(range(len(pivot_data.index)))
        ax.set_yticklabels(pivot_data.index, fontsize=9)

        plt.colorbar(im, ax=ax, label="Mean Replicate Correlation")

    ax.set_title("Mean Replicate Correlation by Condition", fontsize=14)

    plt.tight_layout()
    output_path = output_dir / "replicate_correlation_heatmap.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path


# ============================================================================
# Main Pipeline
# ============================================================================


def setup_logging(log_dir: Path) -> logging.Logger:
    """Setup logging to file and console."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{datetime.today():%Y-%m-%d}_sample_qc.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ]
    )
    return logging.getLogger(__name__)


def main():
    """Main QC pipeline."""
    start_time = datetime.now()
    print(f"=" * 70)
    print(f"Step B1: Sample Quality Control")
    print(f"Started: {start_time}")
    print(f"=" * 70)

    # Parse arguments
    parser = argparse.ArgumentParser(description="Sample QC for BC1 screen data")
    parser.add_argument("--config", type=str, help="Path to config YAML (optional)")
    args = parser.parse_args()

    # Load configuration
    config = SampleQCConfig()

    # Create output directories
    QC_DIR.mkdir(parents=True, exist_ok=True)
    QC_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # Check inputs
    if not LONG_TABLE_FILE.exists():
        print(f"ERROR: Long table not found: {LONG_TABLE_FILE}")
        sys.exit(1)
    if not KEYFILE.exists():
        print(f"ERROR: Keyfile not found: {KEYFILE}")
        sys.exit(1)

    # ========================================================================
    # Step 1: Load keyfile
    # ========================================================================
    print(f"\n[1/6] Loading keyfile: {KEYFILE}")
    keyfile = pd.read_csv(KEYFILE, sep="\t")
    print(f"  Loaded {len(keyfile)} samples")
    print(f"  Libraries: {keyfile['library_ID'].unique().tolist()}")
    print(f"  Conditions: {keyfile['condition'].unique().tolist()}")

    # ========================================================================
    # Step 2: Compute per-sample metrics
    # ========================================================================
    print(f"\n[2/6] Computing per-sample QC metrics")
    sample_metrics = compute_sample_metrics_chunked(LONG_TABLE_FILE, keyfile, config)
    print(f"  Computed metrics for {len(sample_metrics)} samples")

    # ========================================================================
    # Step 3: Compute replicate concordance
    # ========================================================================
    print(f"\n[3/6] Computing replicate concordance")
    concordance_df = compute_replicate_concordance(
        sample_metrics, LONG_TABLE_FILE, keyfile, config
    )
    sample_metrics = sample_metrics.merge(concordance_df, on="sample_name", how="left")

    # ========================================================================
    # Step 4: PCA and outlier detection
    # ========================================================================
    print(f"\n[4/6] Computing PCA and outlier scores")
    pca_df, distance_df = compute_pca_and_outliers(sample_metrics, LONG_TABLE_FILE, config)
    sample_metrics = sample_metrics.merge(distance_df, on="sample_name", how="left")

    # ========================================================================
    # Step 5: Detect outliers
    # ========================================================================
    print(f"\n[5/6] Detecting outliers")

    # Add qc_pass column
    flagged_df = detect_outliers(sample_metrics, config)
    flagged_samples = set(flagged_df["sample_name"].unique()) if len(flagged_df) > 0 else set()
    sample_metrics["qc_pass"] = ~sample_metrics["sample_name"].isin(flagged_samples)

    # Print summary
    n_pass = sample_metrics["qc_pass"].sum()
    n_fail = len(sample_metrics) - n_pass
    print(f"  QC PASS: {n_pass} samples ({100*n_pass/len(sample_metrics):.1f}%)")
    print(f"  QC FAIL: {n_fail} samples ({100*n_fail/len(sample_metrics):.1f}%)")

    if len(flagged_df) > 0:
        print("\n  Flag summary:")
        for flag, count in flagged_df.groupby("flag").size().items():
            print(f"    {flag}: {count}")

    # ========================================================================
    # Step 6: Generate visualizations
    # ========================================================================
    print(f"\n[6/6] Generating QC plots")

    plot_pca(pca_df, sample_metrics, list(flagged_samples), QC_PLOTS_DIR)
    plot_read_count_distribution(sample_metrics, config, QC_PLOTS_DIR)
    plot_clustering_dendrogram(pca_df, sample_metrics, QC_PLOTS_DIR)
    plot_correlation_heatmap(sample_metrics, QC_PLOTS_DIR)

    # ========================================================================
    # Save outputs
    # ========================================================================
    print(f"\nSaving outputs to {QC_DIR}")

    # Sample QC summary
    summary_path = QC_DIR / "sample_qc_summary.tsv"
    sample_metrics.to_csv(summary_path, sep="\t", index=False)
    print(f"  Saved: {summary_path}")

    # Flagged samples
    flagged_path = QC_DIR / "flagged_samples.tsv"
    flagged_df.to_csv(flagged_path, sep="\t", index=False)
    print(f"  Saved: {flagged_path}")

    # PCA coordinates
    pca_path = QC_DIR / "pca_coordinates.tsv"
    pca_df.to_csv(pca_path, sep="\t", index=False)
    print(f"  Saved: {pca_path}")

    # Summary
    elapsed = datetime.now() - start_time
    print(f"\n{'=' * 70}")
    print(f"Step B1 Complete")
    print(f"Elapsed time: {elapsed}")
    print(f"{'=' * 70}")

    return sample_metrics, flagged_df


if __name__ == "__main__":
    main()
