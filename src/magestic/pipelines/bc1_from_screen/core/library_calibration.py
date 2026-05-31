"""
Library and yeast_sublibrary-specific bias correction for hit calling.

This module addresses systematic biases in log2 fold changes that arise from:

1. **Library-level bias** (yL437 vs yL442):
   - yL437 (bottlenecked): Fitness-correlated dropout → negative LFC bias
   - yL442 (non-bottlenecked): Different dynamics → positive LFC bias

2. **Yeast_sublibrary-level bias** (nuclease/strain background effects):
   - Each yeast_sublibrary corresponds to a trafo_ID (transformation)
   - trafo_1 (SpG_NGNG), trafo_2 (SpG_NGNH), trafo_3 (SpG_NGG): Same nuclease
   - trafo_4 (SpCas9_NGG): Different nuclease/strain background
   - Example: In Cobalt_chloride, trafo_4 has median LFC=-0.62 while trafo_3=-0.04
   - CRITICAL: trafo_3 and trafo_4 have SAME oligos but 0.58 LFC difference!

TERMINOLOGY NOTE (2026-02-14):
- "yeast_sublibrary" = trafo_ID at the level of trafo_1, trafo_2, trafo_3, trafo_4
- Replicate plates (trafo_1-1, trafo_1-2, etc.) can be combined
- This is distinct from "oligo sublibrary" (V library) or "PAM sublibrary"
- Use yeast_sublibrary (not sublibrary) to avoid ambiguity

The solution is to center LFCs by synonymous median per (library, yeast_sublibrary, condition),
which removes both library-level and yeast_sublibrary-level systematic bias.

IMPORTANT WORKFLOW NOTE (2026-02-13):
For yL437/yL442, even though yeast_sublibraries were pooled early (at plate wash stage),
yeast_sublibrary-specific mean fitness differences persist through outgrowth. This is
due to strain/nuclease background effects that don't equilibrate during outgrowth.

Recommended approach: Analyze yeast_sublibraries separately, then merge hits.
This is safer than pooled analysis even when libraries are physically combined.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Union
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation


@dataclass
class LibraryBiasStats:
    """Statistics for library/yeast_sublibrary-specific bias correction."""
    library: str
    condition: str
    yeast_sublibrary: str = ""  # Optional yeast_sublibrary (trafo_ID) for stratified analysis
    syn_n: int = 0
    syn_median_lfc: float = 0.0
    syn_mad: float = 0.0
    syn_mad_scaled: float = 0.0  # MAD * 1.4826 (robust std estimate; see variance_calibration.py)
    centering_applied: bool = False


# Type alias for bias keys - can be 2-tuple or 3-tuple
BiasKey = tuple  # Either (library, condition) or (library, yeast_sublibrary, condition)


def calculate_library_bias(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    lfc_col: str = "log2FoldChange",
    variant_type_col: str = "codon_variant_type",
    yeast_sublibrary_col: Optional[str] = None,
    min_synonymous: int = 50,
) -> Dict[BiasKey, LibraryBiasStats]:
    """
    Calculate library/yeast_sublibrary-specific bias using synonymous variants.

    For each grouping (library, [yeast_sublibrary], condition):
    1. Extract synonymous variants
    2. Calculate median LFC (this is the bias to correct)
    3. Calculate MAD for z-score denominator

    Args:
        df: DataFrame with DESeq2 results
        library_col: Column with library identifier
        condition_col: Column with condition/comparison name
        lfc_col: Column with log2 fold change values
        variant_type_col: Column with variant type (synonymous/missense/etc)
        yeast_sublibrary_col: Optional column for yeast_sublibrary (trafo_ID) stratification.
            If provided, bias is calculated per (library, yeast_sublibrary, condition).
            This is RECOMMENDED for yL437/yL442 due to nuclease background effects.
        min_synonymous: Minimum synonymous variants required for estimation

    Returns:
        Dictionary mapping bias key to LibraryBiasStats.
        Key is (library, condition) if yeast_sublibrary_col is None,
        or (library, yeast_sublibrary, condition) if yeast_sublibrary_col is provided.
    """
    bias_stats = {}

    # Create synonymous mask
    syn_mask = df[variant_type_col].str.lower().str.contains("synonymous", na=False)
    syn_df = df[syn_mask].copy()

    if len(syn_df) == 0:
        print("  WARNING: No synonymous variants found for bias calculation")
        return bias_stats

    # Determine grouping columns
    if yeast_sublibrary_col and yeast_sublibrary_col in syn_df.columns:
        group_cols = [library_col, yeast_sublibrary_col, condition_col]
        use_yeast_sublibrary = True
    else:
        group_cols = [library_col, condition_col]
        use_yeast_sublibrary = False
        if yeast_sublibrary_col:
            print(f"  WARNING: yeast_sublibrary_col '{yeast_sublibrary_col}' not found, using library-level only")

    # Group and calculate bias
    for group_key, group in syn_df.groupby(group_cols):
        lfc_values = group[lfc_col].dropna()

        if len(lfc_values) < min_synonymous:
            continue

        syn_median = float(lfc_values.median())
        syn_mad = float(median_abs_deviation(lfc_values, scale='normal'))

        if use_yeast_sublibrary:
            library, yeast_sublibrary, condition = group_key
            bias_key = (library, yeast_sublibrary, condition)
        else:
            library, condition = group_key
            yeast_sublibrary = ""
            bias_key = (library, condition)

        bias_stats[bias_key] = LibraryBiasStats(
            library=str(library),
            yeast_sublibrary=str(yeast_sublibrary),
            condition=str(condition),
            syn_n=len(lfc_values),
            syn_median_lfc=syn_median,
            syn_mad=float(median_abs_deviation(lfc_values)),  # Unscaled
            syn_mad_scaled=syn_mad,  # Scaled to approx std
        )

    return bias_stats


def center_by_library(
    df: pd.DataFrame,
    bias_stats: Dict[BiasKey, LibraryBiasStats],
    library_col: str = "library",
    condition_col: str = "comparison",
    yeast_sublibrary_col: Optional[str] = None,
    lfc_col: str = "log2FoldChange",
    output_col: str = "lfc_centered",
) -> pd.DataFrame:
    """
    Apply bias correction by centering LFCs to synonymous median.

    After centering:
    - lfc_centered = lfc - synonymous_median_for_this_group
    - Group is (library, yeast_sublibrary, condition) if yeast_sublibrary_col provided,
      otherwise (library, condition)
    - Synonymous variants should have median lfc_centered ≈ 0
    - Z-scores calculated from lfc_centered are properly calibrated

    Args:
        df: DataFrame with DESeq2 results
        bias_stats: Dictionary from calculate_library_bias()
        library_col: Column with library identifier
        condition_col: Column with condition/comparison name
        yeast_sublibrary_col: Optional column for yeast_sublibrary (must match what was
            used in calculate_library_bias)
        lfc_col: Column with log2 fold change values
        output_col: Name for the centered LFC column

    Returns:
        DataFrame with added lfc_centered column
    """
    df = df.copy()

    # Determine if bias_stats uses sublibrary (check key length)
    sample_key = next(iter(bias_stats.keys())) if bias_stats else None
    use_yeast_sublibrary = sample_key is not None and len(sample_key) == 3

    # Create lookup for bias values
    if use_yeast_sublibrary and yeast_sublibrary_col and yeast_sublibrary_col in df.columns:
        df["_bias_key"] = list(zip(df[library_col], df[yeast_sublibrary_col], df[condition_col]))
    else:
        df["_bias_key"] = list(zip(df[library_col], df[condition_col]))

    # Map bias values
    df["_syn_median"] = df["_bias_key"].map(
        lambda k: bias_stats[k].syn_median_lfc if k in bias_stats else 0.0
    )
    df["_syn_mad_scaled"] = df["_bias_key"].map(
        lambda k: bias_stats[k].syn_mad_scaled if k in bias_stats else 1.0
    )

    # Apply centering
    df[output_col] = df[lfc_col] - df["_syn_median"]

    # Store scaling factor for z-score calculation
    df["syn_mad_scaled"] = df["_syn_mad_scaled"]

    # Clean up temporary columns
    df = df.drop(columns=["_bias_key", "_syn_median", "_syn_mad_scaled"])

    return df


def validate_centering(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    yeast_sublibrary_col: Optional[str] = None,
    lfc_centered_col: str = "lfc_centered",
    variant_type_col: str = "codon_variant_type",
    max_acceptable_median: float = 0.1,
) -> pd.DataFrame:
    """
    Validate that centering worked correctly.

    After centering, synonymous median should be ~0 for each group.

    Args:
        df: DataFrame after centering
        library_col: Column with library identifier
        condition_col: Column with condition/comparison name
        yeast_sublibrary_col: Optional column for yeast_sublibrary stratification
        lfc_centered_col: Column with centered LFC values
        variant_type_col: Column with variant type
        max_acceptable_median: Maximum acceptable absolute median (flag if exceeded)

    Returns:
        DataFrame with validation results per group
    """
    # Filter to synonymous
    syn_mask = df[variant_type_col].str.lower().str.contains("synonymous", na=False)
    syn_df = df[syn_mask]

    # Determine grouping columns
    if yeast_sublibrary_col and yeast_sublibrary_col in syn_df.columns:
        group_cols = [library_col, yeast_sublibrary_col, condition_col]
    else:
        group_cols = [library_col, condition_col]

    # Calculate post-centering statistics
    validation_records = []

    for group_key, group in syn_df.groupby(group_cols):
        centered_values = group[lfc_centered_col].dropna()

        if len(centered_values) == 0:
            continue

        post_median = float(centered_values.median())
        post_mad = float(median_abs_deviation(centered_values, scale='normal'))

        if len(group_cols) == 3:
            library, yeast_sublibrary, condition = group_key
        else:
            library, condition = group_key
            yeast_sublibrary = ""

        validation_records.append({
            "library": library,
            "yeast_sublibrary": yeast_sublibrary,
            "condition": condition,
            "n_synonymous": len(centered_values),
            "post_centering_median": post_median,
            "post_centering_mad": post_mad,
            "centering_valid": abs(post_median) <= max_acceptable_median,
        })

    validation_df = pd.DataFrame(validation_records)

    # Report any failures
    if len(validation_df) > 0:
        n_invalid = (~validation_df["centering_valid"]).sum()
        if n_invalid > 0:
            group_name = "(library, yeast_sublibrary, condition)" if yeast_sublibrary_col else "(library, condition)"
            print(f"  WARNING: {n_invalid} {group_name} groups have |median| > {max_acceptable_median}")
            invalid = validation_df[~validation_df["centering_valid"]]
            for _, row in invalid.iterrows():
                if row.get('yeast_sublibrary'):
                    print(f"    {row['library']} × {row['yeast_sublibrary']} × {row['condition']}: median = {row['post_centering_median']:.3f}")
                else:
                    print(f"    {row['library']} × {row['condition']}: median = {row['post_centering_median']:.3f}")

    return validation_df


def calculate_z_scores_empirical(
    df: pd.DataFrame,
    lfc_centered_col: str = "lfc_centered",
    syn_mad_col: str = "syn_mad_scaled",
    output_col: str = "z_score_empirical",
) -> pd.DataFrame:
    """
    Calculate z-scores using empirical (synonymous-based) scaling.

    z_score_empirical = lfc_centered / syn_mad_scaled

    This gives properly calibrated z-scores where:
    - Synonymous variants have z ≈ N(0, 1) distribution
    - |z| > 3 corresponds to ~0.3% of synonymous (if normal)

    Args:
        df: DataFrame with centered LFC and MAD values
        lfc_centered_col: Column with centered LFC values
        syn_mad_col: Column with synonymous MAD (scaled)
        output_col: Name for z-score column

    Returns:
        DataFrame with added z_score_empirical column
    """
    df = df.copy()

    # Calculate z-scores, handling division by zero
    mad_values = df[syn_mad_col].replace(0, np.nan)
    df[output_col] = df[lfc_centered_col] / mad_values

    return df


def generate_bias_report(
    bias_stats: Dict[BiasKey, LibraryBiasStats],
) -> pd.DataFrame:
    """
    Generate a report of bias statistics.

    Args:
        bias_stats: Dictionary from calculate_library_bias()

    Returns:
        DataFrame with bias statistics per group
    """
    records = []
    for key, stats in bias_stats.items():
        records.append({
            "library": stats.library,
            "yeast_sublibrary": stats.yeast_sublibrary,
            "condition": stats.condition,
            "syn_n": stats.syn_n,
            "syn_median_lfc": stats.syn_median_lfc,
            "syn_mad": stats.syn_mad,
            "syn_mad_scaled": stats.syn_mad_scaled,
            "bias_magnitude": abs(stats.syn_median_lfc),
        })

    report_df = pd.DataFrame(records)

    if len(report_df) > 0:
        sort_cols = ["library", "yeast_sublibrary", "bias_magnitude"] if report_df["yeast_sublibrary"].any() else ["library", "bias_magnitude"]
        report_df = report_df.sort_values(
            sort_cols,
            ascending=[True] * (len(sort_cols) - 1) + [False]
        )

    return report_df


def run_bias_correction_pipeline(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    yeast_sublibrary_col: Optional[str] = None,
    lfc_col: str = "log2FoldChange",
    variant_type_col: str = "codon_variant_type",
    min_synonymous: int = 50,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run the complete bias correction pipeline.

    Steps:
    1. Calculate synonymous-based bias per group
    2. Center LFCs by subtracting synonymous median
    3. Calculate empirical z-scores
    4. Validate centering worked

    Args:
        df: DataFrame with DESeq2 results
        library_col: Column with library identifier
        condition_col: Column with condition/comparison name
        yeast_sublibrary_col: Optional column for yeast_sublibrary (trafo_ID) stratification.
            RECOMMENDED for yL437/yL442 due to nuclease background effects.
            If provided, bias is calculated per (library, yeast_sublibrary, condition).
        lfc_col: Column with log2 fold change values
        variant_type_col: Column with variant type
        min_synonymous: Minimum synonymous variants for estimation
        verbose: Whether to print progress

    Returns:
        Tuple of (corrected_df, bias_report, validation_report)
    """
    stratification_mode = "yeast_sublibrary-stratified" if yeast_sublibrary_col else "library-level"

    if verbose:
        print("\n" + "=" * 80)
        print(f"BIAS CORRECTION ({stratification_mode.upper()})")
        print("=" * 80)

    # Step 1: Calculate bias
    if verbose:
        grouping = f"(library, yeast_sublibrary, condition)" if yeast_sublibrary_col else "(library, condition)"
        print(f"\n1. Calculating synonymous-based bias per {grouping}...")

    bias_stats = calculate_library_bias(
        df, library_col, condition_col, lfc_col, variant_type_col,
        yeast_sublibrary_col=yeast_sublibrary_col, min_synonymous=min_synonymous
    )

    if verbose:
        print(f"   Calculated bias for {len(bias_stats)} groups")

        # Show summary by library (and sublibrary if applicable)
        if bias_stats:
            if yeast_sublibrary_col:
                by_yeast_sublibrary = {}
                for key, stats in bias_stats.items():
                    sublib = stats.yeast_sublibrary
                    if sublib not in by_yeast_sublibrary:
                        by_yeast_sublibrary[sublib] = []
                    by_yeast_sublibrary[sublib].append(stats.syn_median_lfc)

                for sublib, medians in sorted(by_yeast_sublibrary.items()):
                    mean_bias = np.mean(medians)
                    print(f"   {sublib}: mean bias = {mean_bias:+.3f}")
            else:
                by_library = {}
                for key, stats in bias_stats.items():
                    lib = stats.library
                    if lib not in by_library:
                        by_library[lib] = []
                    by_library[lib].append(stats.syn_median_lfc)

                for lib, medians in sorted(by_library.items()):
                    mean_bias = np.mean(medians)
                    print(f"   {lib}: mean bias = {mean_bias:+.3f}")

    # Step 2: Center LFCs
    if verbose:
        print("\n2. Centering LFCs by synonymous median...")
    df_corrected = center_by_library(
        df, bias_stats, library_col, condition_col,
        yeast_sublibrary_col=yeast_sublibrary_col, lfc_col=lfc_col
    )

    # Step 3: Calculate z-scores
    if verbose:
        print("\n3. Calculating empirical z-scores...")
    df_corrected = calculate_z_scores_empirical(df_corrected)

    # Step 4: Validate centering
    if verbose:
        print("\n4. Validating centering...")
    validation_df = validate_centering(
        df_corrected, library_col, condition_col,
        yeast_sublibrary_col=yeast_sublibrary_col,
        lfc_centered_col="lfc_centered",
        variant_type_col=variant_type_col
    )

    if verbose and len(validation_df) > 0:
        n_valid = validation_df["centering_valid"].sum()
        n_total = len(validation_df)
        print(f"   Centering valid for {n_valid}/{n_total} groups")

    # Generate bias report
    bias_report = generate_bias_report(bias_stats)

    if verbose:
        print("\n" + "=" * 80)
        print("BIAS CORRECTION COMPLETE")
        print("=" * 80)

    return df_corrected, bias_report, validation_df
