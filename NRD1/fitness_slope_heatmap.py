#!/usr/bin/env python3
"""
Fitness calling for saturation mutagenesis screen.

Performs per-barcode normalization, series trimming, slope calculation,
and hit calling with barcode-level support analysis. Generates heatmaps,
dispersion plots, and comprehensive output tables.
"""

import argparse
import os
from pathlib import Path
from typing import Tuple, List, Optional
import warnings

import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================================
# CONFIGURATION
# ============================================================================

DEFAULT_MIN_TOTAL_COUNT = 10          # Keep barcode if avg G0 >= this
DEFAULT_EPS = 1e-8                    # Pseudocount to avoid log(0)
DEFAULT_USE_WEIGHTED = True           # Weight by G0 when collapsing across barcodes
DEFAULT_FC_TRIM_CUTOFF = 0.01         # Drop series after it hits k consecutive points below this FC
DEFAULT_SUSTAINED_K = 3               # Number of consecutive points below cutoff (post-G0) to trigger a trim

# Hit calling parameters
DEFAULT_TREAT_YPD = "YPD+IAA"
DEFAULT_TREAT_YPD_NOIAA = "YPD"
DEFAULT_HIT_TAIL_PROB = 0.02          # 2nd percentile cutoff for YPD adjusted slopes
DEFAULT_DIFF_FACTOR_RULE = 3          # Require slope_adj_YPD / 3 < slope_adj_YPD_noIAA
DEFAULT_DEPLETION_FC_CUT = 1e-4       # Depletion in last FC for YPD


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def trim_series(
    generations: np.ndarray,
    values: np.ndarray,
    cutoff_fc: float = 0.5,
    k: int = 2
) -> pd.DataFrame:
    """
    Single trimming rule: find first run of k consecutive post-G0 points with log2FC < log2(cutoff_fc);
    keep up to and including that generation. NA does NOT count toward the run.

    Args:
        generations: Array of generation numbers.
        values: Array of log2FC values.
        cutoff_fc: Fold change cutoff (default 0.5).
        k: Number of consecutive points below cutoff to trigger trim (default 2).

    Returns:
        DataFrame with columns 'generation' and 'y' (values).
    """
    cutoff_l2 = np.log2(cutoff_fc)

    # Sort by generation
    order = np.argsort(generations)
    g = generations[order]
    y = values[order]

    # Keep only finite values
    finite_mask = np.isfinite(y) & np.isfinite(g)
    g = g[finite_mask]
    y = y[finite_mask]

    if len(g) == 0:
        return pd.DataFrame({"generation": np.array([]), "y": np.array([])})

    # If no G0, no trim
    if not np.any(g == 0):
        return pd.DataFrame({"generation": g, "y": y})

    # Find post-G0 points
    post_mask = g > 0
    gens_post = g[post_mask]
    y_post = y[post_mask]

    if len(y_post) == 0:
        return pd.DataFrame({"generation": g, "y": y})

    # Find first run of k consecutive points below cutoff
    below = y_post < cutoff_l2
    cnt = 0
    first_idx = None

    for i, is_below in enumerate(below):
        cnt = cnt + 1 if is_below else 0
        if cnt >= k:
            first_idx = i
            break

    # If found a run, trim at that generation
    if first_idx is not None:
        cut_gen = gens_post[first_idx]
        keep = g <= cut_gen
        return pd.DataFrame({"generation": g[keep], "y": y[keep]})

    return pd.DataFrame({"generation": g, "y": y})


def safe_lm_slope(
    x: np.ndarray,
    y: np.ndarray,
    min_points: int = 3
) -> Optional[float]:
    """
    Fit a simple linear regression (y ~ x) and return the slope.
    Returns NA if fewer than min_points finite pairs.

    Args:
        x: X values (e.g., generations).
        y: Y values (e.g., log2FC).
        min_points: Minimum number of finite pairs required (default 3).

    Returns:
        Slope of the linear regression, or np.nan if insufficient data.
    """
    # Create dataframe with finite x and y
    df = pd.DataFrame({"x": np.asarray(x, dtype=float), "y": np.asarray(y, dtype=float)})
    df = df[df[["x", "y"]].notna().all(axis=1)]
    df = df.drop_duplicates(subset=["x"], keep="first").sort_values("x")

    if len(df) < min_points:
        return np.nan

    # Fit linear regression
    slope, intercept, r_value, p_value, std_err = linregress(df["x"], df["y"])
    return slope


def weighted_mean(
    values: np.ndarray,
    weights: np.ndarray,
    na_rm: bool = True
) -> float:
    """
    Compute weighted mean.

    Args:
        values: Array of values.
        weights: Array of weights.
        na_rm: If True, ignore NaN values.

    Returns:
        Weighted mean.
    """
    if na_rm:
        mask = np.isfinite(values) & np.isfinite(weights)
        values = values[mask]
        weights = weights[mask]

    return np.average(values, weights=weights)


# ============================================================================
# LOAD AND PROCESS DATA
# ============================================================================

def load_and_prepare_data(
    in_file: str,
    min_total_count: int = DEFAULT_MIN_TOTAL_COUNT,
    eps: float = DEFAULT_EPS,
    use_weighted: bool = DEFAULT_USE_WEIGHTED
) -> pd.DataFrame:
    """
    Load input CSV, add mutation_type column, convert to long format,
    and perform per-barcode normalization.

    Args:
        in_file: Path to input CSV file.
        min_total_count: Minimum average G0 count to keep a barcode.
        eps: Pseudocount to avoid log(0).
        use_weighted: Whether to weight by G0.

    Returns:
        DataFrame with normalized fitness data (long format).
    """
    # Load data
    wide = pd.read_csv(in_file)

    # Ensure bc1 is string and amino_acid_position is int
    wide["bc1"] = wide["bc1"].astype(str)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wide["amino_acid_position"] = pd.to_numeric(
            wide["amino_acid_position"], errors="coerce"
        ).astype("Int64")

    # Filter for passing bc0 fragments if column exists
    if "multiple_donor_bc0_fragments_passing_filter" in wide.columns:
        col = wide["multiple_donor_bc0_fragments_passing_filter"]
        wide = wide[col.isin([True, "TRUE", "True"])]

    # Classify mutations
    ref = wide.get("ref_amino_acid", pd.Series([np.nan] * len(wide))).str.strip()
    mut = wide.get("mutant_amino_acid", pd.Series([np.nan] * len(wide))).str.strip()

    mutation_type = np.where(
        ref.isna() | mut.isna(),
        np.nan,
        np.where(
            mut == "*",
            "PTC",
            np.where(ref == mut, "synonymous", "missense")
        )
    )

    # Insert mutation_type before first *_Normalized_Count_* column
    norm_cols = [col for col in wide.columns if "_Normalized_Count_" in col]
    if norm_cols:
        insert_pos = wide.columns.get_loc(norm_cols[0])
    else:
        insert_pos = len(wide.columns)

    wide.insert(insert_pos, "mutation_type", mutation_type)

    # Convert to long format
    id_vars = [col for col in wide.columns if col != "mutation_type" and not "_Normalized_Count_" in col]
    id_vars.insert(
        [col for col in wide.columns if col != "mutation_type"].index(
            [col for col in wide.columns if "_Normalized_Count_" in col][0]
            if [col for col in wide.columns if "_Normalized_Count_" in col] else wide.columns[0]
        ),
        "mutation_type"
    )

    # Actually do the melt more cleanly
    norm_cols = [col for col in wide.columns if "_Normalized_Count_" in col]
    id_cols = [col for col in wide.columns if col not in norm_cols]

    long = pd.melt(
        wide,
        id_vars=id_cols,
        value_vars=norm_cols,
        var_name="sample",
        value_name="fitness"
    )

    # Parse sample column to extract treatment, replicate, set_flag, generation
    long["treatment_raw"] = long["sample"].str.replace("_Normalized_Count_.*", "", regex=True)
    long["replicate"] = long["treatment_raw"].str.extract(r"_R(\d+)$", expand=False).astype("Int64")
    long["set_flag"] = long["treatment_raw"].str.contains("_set$", regex=True)
    long["treatment"] = (
        long["treatment_raw"]
        .str.replace(r"_R\d+$", "", regex=True)
        .str.replace(r"_set$", "", regex=True)
        .str.replace(r"\s+", "_", regex=True)
    )
    long["generation"] = long["sample"].str.extract(r"Normalized_Count_(\d+)", expand=False).astype("Int64")

    # Keep only relevant columns
    long = long[[
        "bc1", "gene_name", "amino_acid_position", "mutation_type",
        "sample", "treatment", "set_flag", "replicate", "generation", "fitness"
    ]]

    # Remove rows with missing amino_acid_position
    long = long[long["amino_acid_position"].notna()]

    # Calculate G0 means per barcode
    gen0_means = long[long["generation"] == 0].groupby("bc1")["fitness"].mean().reset_index()
    gen0_means.columns = ["bc1", "g0"]

    # Filter for valid barcodes
    gen0_means = gen0_means[
        (gen0_means["g0"].notna()) &
        (np.isfinite(gen0_means["g0"])) &
        (gen0_means["g0"] >= min_total_count)
    ]
    valid_bc = set(gen0_means["bc1"])

    # Join with means and normalize
    df = long[long["bc1"].isin(valid_bc)].copy()
    df = df.merge(gen0_means, on="bc1", how="left")

    # Normalize fitness to G0, compute FC and log2FC
    df["fitness"] = df.apply(
        lambda row: row["g0"] if row["generation"] == 0 else row["fitness"],
        axis=1
    )
    df["fc"] = np.maximum(df["fitness"], eps) / np.maximum(df["g0"], eps)
    df["log2fc"] = np.log2(df["fc"])
    df["w"] = np.maximum(df["g0"], 1) if use_weighted else 1.0

    return df


# ============================================================================
# COLLAPSE, TRIM, AND FIT SLOPES
# ============================================================================

def collapse_and_trim(
    df: pd.DataFrame,
    fc_trim_cutoff: float = DEFAULT_FC_TRIM_CUTOFF,
    sustained_k: int = DEFAULT_SUSTAINED_K,
    use_weighted: bool = DEFAULT_USE_WEIGHTED
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Collapse to mean log2FC by variant x treatment x generation,
    determine trimming point for each variant x treatment, and fit slopes on trimmed data.

    Args:
        df: Normalized fitness dataframe.
        fc_trim_cutoff: Fold change cutoff for trimming.
        sustained_k: Number of consecutive points below cutoff.
        use_weighted: Whether to use weighted mean.

    Returns:
        Tuple of (all_summaries_fc, collapsed_slopes) dataframes.
    """
    # Collapse to mean log2FC by variant x treatment x generation
    if use_weighted:
        all_summaries_fc = df.groupby(
            ["amino_acid_position", "mutation_type", "treatment", "generation"]
        ).apply(
            lambda g: pd.Series({
                "mean_log2fc": weighted_mean(g["log2fc"].values, g["w"].values, na_rm=True)
            }),
            include_groups=False
        ).reset_index()
    else:
        all_summaries_fc = df.groupby(
            ["amino_acid_position", "mutation_type", "treatment", "generation"]
        )["log2fc"].mean().reset_index()
        all_summaries_fc.columns = ["amino_acid_position", "mutation_type", "treatment", "generation", "mean_log2fc"]

    # Determine cut_gen from collapsed series
    trim_results = []
    for (pos, mut_type, treatment), group in all_summaries_fc.groupby(
        ["amino_acid_position", "mutation_type", "treatment"]
    ):
        ser = trim_series(
            group["generation"].values,
            group["mean_log2fc"].values,
            cutoff_fc=fc_trim_cutoff,
            k=sustained_k
        )

        if len(ser) == 0:
            cut_gen = group[group["mean_log2fc"].notna()]["generation"].max()
        else:
            cut_gen = ser["generation"].max()

        trim_results.append({
            "amino_acid_position": pos,
            "mutation_type": mut_type,
            "treatment": treatment,
            "cut_gen": cut_gen
        })

    trim_info = pd.DataFrame(trim_results)

    # Fit slopes on collapsed, trimmed series
    slope_results = []
    for (pos, mut_type, treatment), group in all_summaries_fc.groupby(
        ["amino_acid_position", "mutation_type", "treatment"]
    ):
        cut_gen = trim_info[
            (trim_info["amino_acid_position"] == pos) &
            (trim_info["mutation_type"] == mut_type) &
            (trim_info["treatment"] == treatment)
        ]["cut_gen"].values[0]

        filtered = group[
            (group["mean_log2fc"].notna()) &
            (group["generation"].notna()) &
            (group["generation"] <= cut_gen)
        ]

        if len(filtered) > 0:
            slope = safe_lm_slope(filtered["generation"].values, filtered["mean_log2fc"].values)
        else:
            slope = np.nan

        slope_results.append({
            "amino_acid_position": pos,
            "mutation_type": mut_type,
            "treatment": treatment,
            "slope": slope
        })

    collapsed_slopes = pd.DataFrame(slope_results)

    return all_summaries_fc, collapsed_slopes


def median_center_slopes(collapsed_slopes: pd.DataFrame) -> pd.DataFrame:
    """
    Median-center slopes within each treatment.

    Args:
        collapsed_slopes: DataFrame with slopes by variant and treatment.

    Returns:
        DataFrame with added slope_adj column.
    """
    slope_table = collapsed_slopes.copy()

    # Median-center within treatment
    for treatment in slope_table["treatment"].unique():
        mask = slope_table["treatment"] == treatment
        median_slope = slope_table[mask]["slope"].median()
        slope_table.loc[mask, "median_slope"] = median_slope
        slope_table.loc[mask, "slope_adj"] = slope_table.loc[mask, "slope"] - median_slope

    return slope_table


# ============================================================================
# DEPLETION ANALYSIS
# ============================================================================

def calculate_depletion(
    all_summaries_fc: pd.DataFrame,
    treat_ypd: str = DEFAULT_TREAT_YPD,
    depletion_fc_cut: float = DEFAULT_DEPLETION_FC_CUT
) -> pd.DataFrame:
    """
    Compare G0 vs last finite mean point for each variant x treatment.

    Args:
        all_summaries_fc: Collapsed mean log2FC dataframe.
        treat_ypd: YPD treatment name.
        depletion_fc_cut: FC cutoff for depletion.

    Returns:
        DataFrame with depletion flags.
    """
    # Get last generation for each variant x treatment
    last_by_group = all_summaries_fc[
        all_summaries_fc["mean_log2fc"].notna()
    ].groupby(
        ["amino_acid_position", "mutation_type", "treatment"]
    )["generation"].max().reset_index()
    last_by_group.columns = ["amino_acid_position", "mutation_type", "treatment", "last_gen"]

    # Get G0 values
    gen0 = all_summaries_fc[all_summaries_fc["generation"] == 0][
        ["amino_acid_position", "mutation_type", "treatment", "mean_log2fc"]
    ].copy()
    gen0.columns = ["amino_acid_position", "mutation_type", "treatment", "log2fc_gen0"]

    # Get last values
    last_vals = all_summaries_fc.merge(
        last_by_group,
        on=["amino_acid_position", "mutation_type", "treatment"]
    )
    last_vals = last_vals[last_vals["generation"] == last_vals["last_gen"]][
        ["amino_acid_position", "mutation_type", "treatment", "mean_log2fc"]
    ].copy()
    last_vals.columns = ["amino_acid_position", "mutation_type", "treatment", "log2fc_last"]

    # Merge and compute depletion
    fitness_extremes = gen0.merge(
        last_vals,
        on=["amino_acid_position", "mutation_type", "treatment"],
        how="outer"
    )

    fitness_extremes["depletion"] = (
        fitness_extremes["log2fc_last"].notna() &
        (fitness_extremes["log2fc_last"] < np.log2(depletion_fc_cut))
    )

    return fitness_extremes


# ============================================================================
# HIT CALLING
# ============================================================================

def call_hits(
    slope_table: pd.DataFrame,
    fitness_extremes: pd.DataFrame,
    treat_ypd: str = DEFAULT_TREAT_YPD,
    treat_ypd_noiaa: str = DEFAULT_TREAT_YPD_NOIAA,
    hit_tail_prob: float = DEFAULT_HIT_TAIL_PROB,
    diff_factor_rule: float = DEFAULT_DIFF_FACTOR_RULE,
    depletion_fc_cut: float = DEFAULT_DEPLETION_FC_CUT
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Call hits based on YPD vs YPD_noIAA slope comparison and depletion.

    Args:
        slope_table: DataFrame with adjusted slopes.
        fitness_extremes: DataFrame with depletion information.
        treat_ypd: YPD treatment name.
        treat_ypd_noiaa: YPD_noIAA treatment name.
        hit_tail_prob: Percentile cutoff for hit threshold.
        diff_factor_rule: Ratio cutoff for differential requirement.
        depletion_fc_cut: FC cutoff for depletion.

    Returns:
        Tuple of (hits, slope_table) dataframes.
    """
    # Calculate hit threshold from YPD distribution
    ypd_slopes = slope_table[slope_table["treatment"] == treat_ypd]["slope_adj"]
    hit_threshold_ypd = ypd_slopes.quantile(hit_tail_prob)

    # Pivot slope_table for YPD vs YPD_noIAA comparison
    slope_wide = slope_table[
        slope_table["treatment"].isin([treat_ypd, treat_ypd_noiaa])
    ].pivot_table(
        index=["amino_acid_position", "mutation_type"],
        columns="treatment",
        values=["slope", "median_slope", "slope_adj"],
        aggfunc="first"
    )
    slope_wide.columns = ["_".join(col).strip() for col in slope_wide.columns.values]
    slope_wide = slope_wide.reset_index()

    # Add depletion info
    depl_by_posmut = fitness_extremes[
        fitness_extremes["treatment"] == treat_ypd
    ][["amino_acid_position", "mutation_type", "depletion"]].copy()

    slope_wide = slope_wide.merge(
        depl_by_posmut,
        on=["amino_acid_position", "mutation_type"],
        how="left"
    )

    # Define column names (handle potential naming variations)
    col_ypd = f"slope_adj_{treat_ypd}"
    col_noiaa = f"slope_adj_{treat_ypd_noiaa}"

    # Call hits
    hits = slope_wide[
        (slope_wide[col_ypd].notna()) &
        (slope_wide[col_noiaa].notna()) &
        (slope_wide["depletion"] == True)
    ].copy()

    hits = hits[
        (hits[col_ypd] < hit_threshold_ypd) &
        (hits[col_ypd] / diff_factor_rule < hits[col_noiaa])
    ]

    hits["hit"] = True

    # Add hit flag to slope_table
    slope_table = slope_table.merge(
        hits[["amino_acid_position", "mutation_type", "hit"]],
        on=["amino_acid_position", "mutation_type"],
        how="left"
    )
    slope_table["hit"] = slope_table["hit"].fillna(False)

    return hits, slope_table


# ============================================================================
# BARCODE-LEVEL SUPPORT
# ============================================================================

def analyze_barcode_support(
    df: pd.DataFrame,
    trim_info: pd.DataFrame,
    slope_table: pd.DataFrame,
    fitness_extremes: pd.DataFrame,
    treat_ypd: str = DEFAULT_TREAT_YPD,
    treat_ypd_noiaa: str = DEFAULT_TREAT_YPD_NOIAA,
    hit_threshold_ypd: float = None,
    diff_factor_rule: float = DEFAULT_DIFF_FACTOR_RULE,
    depletion_fc_cut: float = DEFAULT_DEPLETION_FC_CUT,
    use_weighted: bool = DEFAULT_USE_WEIGHTED
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Analyze per-barcode slopes and support for hit calling.

    Args:
        df: Original normalized fitness dataframe.
        trim_info: Trimming information from collapse_and_trim.
        slope_table: Slope table with hit information.
        fitness_extremes: Depletion information.
        treat_ypd: YPD treatment name.
        treat_ypd_noiaa: YPD_noIAA treatment name.
        hit_threshold_ypd: Hit threshold (percentile from distribution).
        diff_factor_rule: Ratio cutoff.
        depletion_fc_cut: FC cutoff.
        use_weighted: Whether slopes were weighted.

    Returns:
        Tuple of (bc_hit_flags, support_counts) dataframes.
    """
    # Fit per-barcode slopes using same cut_gen window
    df_trimmed = df.merge(
        trim_info,
        on=["amino_acid_position", "mutation_type", "treatment"],
        how="inner"
    )

    df_trimmed = df_trimmed[
        (df_trimmed["log2fc"].notna()) &
        (df_trimmed["generation"].notna()) &
        (df_trimmed["generation"] <= df_trimmed["cut_gen"])
    ]

    per_bc_slopes = []
    for (pos, mut_type, treatment, bc), group in df_trimmed.groupby(
        ["amino_acid_position", "mutation_type", "treatment", "bc1"]
    ):
        slope = safe_lm_slope(group["generation"].values, group["log2fc"].values)
        w = group["w"].iloc[0] if len(group) > 0 else 1.0
        per_bc_slopes.append({
            "amino_acid_position": pos,
            "mutation_type": mut_type,
            "treatment": treatment,
            "bc1": bc,
            "slope": slope,
            "w": w
        })

    per_bc_slopes = pd.DataFrame(per_bc_slopes)

    # Median-center per-barcode slopes by treatment
    per_bc_slopes_adj = per_bc_slopes.copy()
    for treatment in per_bc_slopes_adj["treatment"].unique():
        mask = per_bc_slopes_adj["treatment"] == treatment
        median = per_bc_slopes_adj[mask]["slope"].median()
        per_bc_slopes_adj.loc[mask, "slope_adj"] = per_bc_slopes_adj.loc[mask, "slope"] - median

    # Per-barcode depletion (YPD only)
    df_ypd = df[df["treatment"] == treat_ypd].copy()
    df_ypd = df_ypd.sort_values(
        ["amino_acid_position", "mutation_type", "bc1", "generation"]
    )

    bc_depl_list = []
    for (pos, mut_type, bc), group in df_ypd.groupby(
        ["amino_acid_position", "mutation_type", "bc1"]
    ):
        finite_mask = group["log2fc"].notna()
        if finite_mask.sum() > 0:
            log2fc_vals = group[finite_mask]["log2fc"].values
            gen_vals = group[finite_mask]["generation"].values
            log2fc_last = log2fc_vals[np.argmax(gen_vals)]
        else:
            log2fc_last = np.nan

        depletion_bc = (
            np.isfinite(log2fc_last) and
            (log2fc_last < np.log2(depletion_fc_cut))
        )

        bc_depl_list.append({
            "amino_acid_position": pos,
            "mutation_type": mut_type,
            "bc1": bc,
            "log2fc_last": log2fc_last,
            "depletion_bc": depletion_bc
        })

    bc_depl_ypd = pd.DataFrame(bc_depl_list)

    # Wide per-barcode slopes for YPD vs YPD_noIAA
    bc_slope_wide = per_bc_slopes_adj[
        per_bc_slopes_adj["treatment"].isin([treat_ypd, treat_ypd_noiaa])
    ].pivot_table(
        index=["amino_acid_position", "mutation_type", "bc1"],
        columns="treatment",
        values="slope_adj",
        aggfunc="first"
    ).reset_index()

    bc_slope_wide = bc_slope_wide.merge(
        bc_depl_ypd,
        on=["amino_acid_position", "mutation_type", "bc1"],
        how="left"
    )

    # Dynamic column names
    col_ypd = treat_ypd if treat_ypd in bc_slope_wide.columns else f"slope_adj_{treat_ypd}"
    col_noiaa = treat_ypd_noiaa if treat_ypd_noiaa in bc_slope_wide.columns else f"slope_adj_{treat_ypd_noiaa}"

    if col_ypd not in bc_slope_wide.columns:
        col_ypd = [c for c in bc_slope_wide.columns if treat_ypd in c][0] if any(
            treat_ypd in c for c in bc_slope_wide.columns
        ) else None
    if col_noiaa not in bc_slope_wide.columns:
        col_noiaa = [c for c in bc_slope_wide.columns if treat_ypd_noiaa in c][0] if any(
            treat_ypd_noiaa in c for c in bc_slope_wide.columns
        ) else None

    # Hit calling thresholds
    if hit_threshold_ypd is None:
        hit_threshold_ypd = slope_table[
            slope_table["treatment"] == treat_ypd
        ]["slope_adj"].quantile(0.02)

    # Barcode-level hit flags
    bc_hit_flags = bc_slope_wide.copy()
    bc_hit_flags["supports_hit"] = (
        (bc_hit_flags[col_ypd].notna()) &
        (bc_hit_flags[col_noiaa].notna()) &
        (bc_hit_flags["depletion_bc"] == True) &
        (bc_hit_flags[col_ypd] < hit_threshold_ypd) &
        (bc_hit_flags[col_ypd] / diff_factor_rule < bc_hit_flags[col_noiaa])
    )

    bc_hit_flags["agrees_direction"] = (
        (bc_hit_flags[col_ypd].notna()) &
        (bc_hit_flags[col_noiaa].notna()) &
        (bc_hit_flags[col_ypd] < bc_hit_flags[col_noiaa]) &
        (bc_hit_flags[col_ypd] < 0)
    )

    # Support counts
    support_counts = bc_hit_flags.groupby(
        ["amino_acid_position", "mutation_type"]
    ).agg(
        n_barcodes=("bc1", "count"),
        n_support=("supports_hit", "sum"),
        n_agree=("agrees_direction", "sum")
    ).reset_index()

    support_counts["pct_support"] = support_counts.apply(
        lambda r: r["n_support"] / r["n_barcodes"] if r["n_barcodes"] > 0 else np.nan,
        axis=1
    )
    support_counts["pct_agree"] = support_counts.apply(
        lambda r: r["n_agree"] / r["n_barcodes"] if r["n_barcodes"] > 0 else np.nan,
        axis=1
    )

    return bc_hit_flags, support_counts


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def plot_heatmap_slopes(
    slope_table: pd.DataFrame,
    output_path: str,
    treat_ypd: str = DEFAULT_TREAT_YPD,
    treat_ypd_noiaa: str = DEFAULT_TREAT_YPD_NOIAA
) -> None:
    """
    Create heatmap of adjusted slopes (excluding PTC/Deletion).

    Args:
        slope_table: DataFrame with slope_adj values.
        output_path: Path to save PDF.
        treat_ypd: YPD treatment name.
        treat_ypd_noiaa: YPD_noIAA treatment name.
    """
    # Prepare data for heatmap
    heatmap_table = slope_table[
        ~slope_table["mutation_type"].isin(["PTC", "Deletion"])
    ].copy()

    # Pivot for heatmap
    heatmap_data = heatmap_table.pivot_table(
        index="treatment",
        columns="amino_acid_position",
        values="slope_adj",
        aggfunc="first"
    )

    # Sort treatments
    treatment_order = ["HIT", treat_ypd, treat_ypd_noiaa]
    treatment_order = [t for t in treatment_order if t in heatmap_data.index]
    heatmap_data = heatmap_data.loc[treatment_order]

    # Plot
    fig, axes = plt.subplots(
        len(heatmap_table["mutation_type"].unique()),
        1,
        figsize=(12, 8),
        squeeze=False
    )

    for idx, (mut_type, ax) in enumerate(
        zip(sorted(heatmap_table["mutation_type"].unique()), axes.flatten())
    ):
        data = heatmap_table[heatmap_table["mutation_type"] == mut_type].pivot_table(
            index="treatment",
            columns="amino_acid_position",
            values="slope_adj",
            aggfunc="first"
        )
        data = data.loc[[t for t in treatment_order if t in data.index]]

        sns.heatmap(
            data,
            ax=ax,
            cmap="RdBu_r",
            center=0,
            cbar_kws={"label": "Normalized Slope"},
            vmin=-data.abs().max().max() if data.abs().max().max() > 0 else -1,
            vmax=data.abs().max().max() if data.abs().max().max() > 0 else 1,
            linewidths=0.5,
            linecolor="grey"
        )
        ax.set_title(f"Mutation Type: {mut_type}", fontsize=12, fontweight="bold")
        ax.set_xlabel("Amino Acid Position")
        ax.set_ylabel("Treatment")

    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches="tight")
    plt.close()


def plot_dispersion(
    slope_table: pd.DataFrame,
    output_path: str,
    treat_ypd: str = DEFAULT_TREAT_YPD,
    treat_ypd_noiaa: str = DEFAULT_TREAT_YPD_NOIAA
) -> None:
    """
    Create dispersion plot: YPD vs YPD_noIAA adjusted slopes.

    Args:
        slope_table: DataFrame with slope_adj values and hit flags.
        output_path: Path to save PDF.
        treat_ypd: YPD treatment name.
        treat_ypd_noiaa: YPD_noIAA treatment name.
    """
    # Pivot for dispersion plot
    plot_data = slope_table[
        slope_table["treatment"].isin([treat_ypd, treat_ypd_noiaa])
    ].pivot_table(
        index=["amino_acid_position", "mutation_type"],
        columns="treatment",
        values=["slope_adj", "hit"],
        aggfunc="first"
    )

    plot_data.columns = ["_".join(col).strip() for col in plot_data.columns.values]
    plot_data = plot_data.reset_index()

    # Drop NA rows
    col_ypd = f"slope_adj_{treat_ypd}"
    col_noiaa = f"slope_adj_{treat_ypd_noiaa}"
    plot_data = plot_data.dropna(subset=[col_ypd, col_noiaa])

    # Color mapping
    color_map = {
        "missense": "red",
        "synonymous": "blue",
        "PTC": "grey",
        "Deletion": "black"
    }
    colors = plot_data["mutation_type"].map(color_map).fillna("grey")
    sizes = plot_data[f"hit_{treat_ypd}"].fillna(False).map({True: 100, False: 20})

    fig, ax = plt.subplots(figsize=(12, 8))

    for mut_type in sorted(plot_data["mutation_type"].unique()):
        subset = plot_data[plot_data["mutation_type"] == mut_type]
        subset_hit = subset[subset[f"hit_{treat_ypd}"] == True]
        subset_not_hit = subset[subset[f"hit_{treat_ypd}"] != True]

        if len(subset_not_hit) > 0:
            ax.scatter(
                subset_not_hit[col_noiaa],
                subset_not_hit[col_ypd],
                c=color_map.get(mut_type, "grey"),
                s=20,
                alpha=0.85,
                label=mut_type if mut_type in ["missense", "synonymous", "PTC", "Deletion"] else None
            )
        if len(subset_hit) > 0:
            ax.scatter(
                subset_hit[col_noiaa],
                subset_hit[col_ypd],
                c=color_map.get(mut_type, "grey"),
                s=100,
                alpha=0.85,
                marker="o"
            )

    ax.axhline(y=0, linestyle="--", color="grey", alpha=0.5)
    ax.axvline(x=0, linestyle="--", color="grey", alpha=0.5)
    ax.set_xlabel(f"Adjusted slope ({treat_ypd_noiaa})", fontsize=13)
    ax.set_ylabel(f"Adjusted slope ({treat_ypd})", fontsize=13)
    ax.set_title("Adjusted log2FC Slopes by Mutation Type", fontsize=13, fontweight="bold")
    ax.legend(title="Mutation Type", loc="best")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches="tight")
    plt.close()


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def main(
    in_file: str,
    out_dir: str,
    min_total_count: int = DEFAULT_MIN_TOTAL_COUNT,
    eps: float = DEFAULT_EPS,
    use_weighted: bool = DEFAULT_USE_WEIGHTED,
    fc_trim_cutoff: float = DEFAULT_FC_TRIM_CUTOFF,
    sustained_k: int = DEFAULT_SUSTAINED_K,
    treat_ypd: str = DEFAULT_TREAT_YPD,
    treat_ypd_noiaa: str = DEFAULT_TREAT_YPD_NOIAA,
    hit_tail_prob: float = DEFAULT_HIT_TAIL_PROB,
    diff_factor_rule: float = DEFAULT_DIFF_FACTOR_RULE,
    depletion_fc_cut: float = DEFAULT_DEPLETION_FC_CUT
) -> None:
    """
    Main analysis pipeline for fitness calling and hit identification.

    Args:
        in_file: Path to input CSV.
        out_dir: Output directory for results.
        min_total_count: Minimum average G0 count.
        eps: Pseudocount.
        use_weighted: Whether to weight by G0.
        fc_trim_cutoff: FC cutoff for trimming.
        sustained_k: Number of consecutive points for trim.
        treat_ypd: YPD treatment name.
        treat_ypd_noiaa: YPD_noIAA treatment name.
        hit_tail_prob: Percentile for hit threshold.
        diff_factor_rule: Ratio cutoff.
        depletion_fc_cut: FC cutoff for depletion.
    """
    os.makedirs(out_dir, exist_ok=True)

    print("Loading and preparing data...")
    df = load_and_prepare_data(
        in_file,
        min_total_count=min_total_count,
        eps=eps,
        use_weighted=use_weighted
    )

    print("Collapsing and trimming series...")
    all_summaries_fc, collapsed_slopes = collapse_and_trim(
        df,
        fc_trim_cutoff=fc_trim_cutoff,
        sustained_k=sustained_k,
        use_weighted=use_weighted
    )

    # Get trim_info for barcode analysis
    trim_results = []
    for (pos, mut_type, treatment), group in all_summaries_fc.groupby(
        ["amino_acid_position", "mutation_type", "treatment"]
    ):
        ser = trim_series(
            group["generation"].values,
            group["mean_log2fc"].values,
            cutoff_fc=fc_trim_cutoff,
            k=sustained_k
        )
        if len(ser) == 0:
            cut_gen = group[group["mean_log2fc"].notna()]["generation"].max()
        else:
            cut_gen = ser["generation"].max()
        trim_results.append({
            "amino_acid_position": pos,
            "mutation_type": mut_type,
            "treatment": treatment,
            "cut_gen": cut_gen
        })
    trim_info = pd.DataFrame(trim_results)

    # Median center slopes
    slope_table = median_center_slopes(collapsed_slopes)

    # Save mean_log2fc by group
    print(f"Saving mean log2FC by group to {os.path.join(out_dir, 'mean_log2FC_by_group.csv')}")
    all_summaries_fc.to_csv(
        os.path.join(out_dir, "mean_log2FC_by_group.csv"),
        index=False
    )

    # Calculate depletion
    print("Calculating depletion...")
    fitness_extremes = calculate_depletion(
        all_summaries_fc,
        treat_ypd=treat_ypd,
        depletion_fc_cut=depletion_fc_cut
    )

    print(f"Saving depletion data to {os.path.join(out_dir, 'depletion_log2FC.csv')}")
    fitness_extremes.to_csv(
        os.path.join(out_dir, "depletion_log2FC.csv"),
        index=False
    )

    # Call hits
    print("Calling hits...")
    hits, slope_table = call_hits(
        slope_table,
        fitness_extremes,
        treat_ypd=treat_ypd,
        treat_ypd_noiaa=treat_ypd_noiaa,
        hit_tail_prob=hit_tail_prob,
        diff_factor_rule=diff_factor_rule,
        depletion_fc_cut=depletion_fc_cut
    )

    # Calculate hit threshold for barcode analysis
    hit_threshold_ypd = slope_table[
        slope_table["treatment"] == treat_ypd
    ]["slope_adj"].quantile(hit_tail_prob)

    # Barcode support analysis
    print("Analyzing barcode-level support...")
    bc_hit_flags, support_counts = analyze_barcode_support(
        df,
        trim_info,
        slope_table,
        fitness_extremes,
        treat_ypd=treat_ypd,
        treat_ypd_noiaa=treat_ypd_noiaa,
        hit_threshold_ypd=hit_threshold_ypd,
        diff_factor_rule=diff_factor_rule,
        depletion_fc_cut=depletion_fc_cut,
        use_weighted=use_weighted
    )

    print(f"Saving barcode flags to {os.path.join(out_dir, 'barcode_support_flags.csv')}")
    bc_hit_flags.to_csv(
        os.path.join(out_dir, "barcode_support_flags.csv"),
        index=False
    )

    print(f"Saving support counts to {os.path.join(out_dir, 'hit_barcode_support_summary.csv')}")
    support_counts.to_csv(
        os.path.join(out_dir, "hit_barcode_support_summary.csv"),
        index=False
    )

    # Merge support counts into slope_table
    slope_table = slope_table.merge(
        support_counts,
        on=["amino_acid_position", "mutation_type"],
        how="left"
    )

    hits = hits.merge(
        support_counts,
        on=["amino_acid_position", "mutation_type"],
        how="left"
    )

    # Add HIT rows to slope_table for visualization
    hit_rows = slope_table[slope_table["hit"] == True][
        ["amino_acid_position", "mutation_type"]
    ].drop_duplicates().copy()
    hit_rows["treatment"] = "HIT"
    hit_rows["slope"] = np.nan
    hit_rows["median_slope"] = np.nan
    hit_rows["slope_adj"] = np.nan
    hit_rows["hit"] = True

    # Copy over support columns if they exist
    for col in ["n_barcodes", "n_support", "n_agree", "pct_support", "pct_agree"]:
        if col in slope_table.columns:
            hit_rows[col] = np.nan

    slope_table_out = pd.concat([slope_table, hit_rows], ignore_index=True)

    # Sort treatments
    treatment_order = ["HIT", treat_ypd, treat_ypd_noiaa]
    slope_table_out["treatment"] = pd.Categorical(
        slope_table_out["treatment"],
        categories=treatment_order,
        ordered=True
    )
    slope_table_out = slope_table_out.sort_values("treatment")

    # Save tables
    print(f"Saving slope table to {os.path.join(out_dir, 'slope_table_log2FC.csv')}")
    slope_table_out.to_csv(
        os.path.join(out_dir, "slope_table_log2FC.csv"),
        index=False
    )

    print(f"Saving hits to {os.path.join(out_dir, 'hits_summary.csv')}")
    hits.to_csv(
        os.path.join(out_dir, "hits_summary.csv"),
        index=False
    )

    # Generate plots
    print(f"Generating heatmap plot...")
    plot_heatmap_slopes(
        slope_table_out,
        os.path.join(out_dir, "heatmap_log2FC_slopes.pdf"),
        treat_ypd=treat_ypd,
        treat_ypd_noiaa=treat_ypd_noiaa
    )

    print(f"Generating dispersion plot...")
    plot_dispersion(
        slope_table,
        os.path.join(out_dir, "dispersion_log2FC_slopes.pdf"),
        treat_ypd=treat_ypd,
        treat_ypd_noiaa=treat_ypd_noiaa
    )

    print("Analysis complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fitness calling for saturation mutagenesis screen"
    )

    # Required arguments
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input CSV file with normalized counts"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for results"
    )

    # Optional analysis parameters
    parser.add_argument(
        "--min-total-count",
        type=int,
        default=DEFAULT_MIN_TOTAL_COUNT,
        help=f"Minimum average G0 count to keep barcode (default: {DEFAULT_MIN_TOTAL_COUNT})"
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=DEFAULT_EPS,
        help=f"Pseudocount to avoid log(0) (default: {DEFAULT_EPS})"
    )
    parser.add_argument(
        "--use-weighted",
        action="store_true",
        default=DEFAULT_USE_WEIGHTED,
        help="Weight by G0 when collapsing across barcodes"
    )
    parser.add_argument(
        "--fc-trim-cutoff",
        type=float,
        default=DEFAULT_FC_TRIM_CUTOFF,
        help=f"FC cutoff for series trimming (default: {DEFAULT_FC_TRIM_CUTOFF})"
    )
    parser.add_argument(
        "--sustained-k",
        type=int,
        default=DEFAULT_SUSTAINED_K,
        help=f"Number of consecutive points below cutoff to trim (default: {DEFAULT_SUSTAINED_K})"
    )

    # Hit calling parameters
    parser.add_argument(
        "--treat-ypd",
        default=DEFAULT_TREAT_YPD,
        help=f"YPD treatment name (default: {DEFAULT_TREAT_YPD})"
    )
    parser.add_argument(
        "--treat-ypd-noiaa",
        default=DEFAULT_TREAT_YPD_NOIAA,
        help=f"YPD_noIAA treatment name (default: {DEFAULT_TREAT_YPD_NOIAA})"
    )
    parser.add_argument(
        "--hit-tail-prob",
        type=float,
        default=DEFAULT_HIT_TAIL_PROB,
        help=f"Percentile cutoff for hit threshold (default: {DEFAULT_HIT_TAIL_PROB})"
    )
    parser.add_argument(
        "--diff-factor-rule",
        type=float,
        default=DEFAULT_DIFF_FACTOR_RULE,
        help=f"Ratio cutoff for differential requirement (default: {DEFAULT_DIFF_FACTOR_RULE})"
    )
    parser.add_argument(
        "--depletion-fc-cut",
        type=float,
        default=DEFAULT_DEPLETION_FC_CUT,
        help=f"FC cutoff for depletion (default: {DEFAULT_DEPLETION_FC_CUT})"
    )

    args = parser.parse_args()

    main(
        in_file=args.input,
        out_dir=args.output,
        min_total_count=args.min_total_count,
        eps=args.eps,
        use_weighted=args.use_weighted,
        fc_trim_cutoff=args.fc_trim_cutoff,
        sustained_k=args.sustained_k,
        treat_ypd=args.treat_ypd,
        treat_ypd_noiaa=args.treat_ypd_noiaa,
        hit_tail_prob=args.hit_tail_prob,
        diff_factor_rule=args.diff_factor_rule,
        depletion_fc_cut=args.depletion_fc_cut
    )
