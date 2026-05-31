"""
Empirical FDR-based threshold derivation for hit calling.

This module derives per-condition z-score thresholds that achieve a target
false discovery rate (FDR) using synonymous variants as the null distribution.

Key insight: Fixed thresholds (e.g., z ≥ 3.0) don't achieve uniform FDR across
conditions because conditions have different noise profiles. By calibrating
thresholds empirically per condition, we ensure consistent FDR control.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np


@dataclass
class ThresholdCalibration:
    """Results of threshold calibration for a condition."""
    condition: str
    library: str
    target_fdr: float
    achieved_fdr: float
    z_threshold: float
    syn_n: int
    syn_n_above_thresh: int
    missense_n: int
    missense_n_above_thresh: int
    missense_syn_ratio: float
    calibration_status: str  # 'OK', 'WARN_HIGH_FDR', 'WARN_LOW_N'


def derive_fdr_threshold(
    z_scores: pd.Series,
    target_fdr: float = 0.01,
    z_min: float = 2.0,
    z_max: float = 6.0,
    z_step: float = 0.1,
) -> Tuple[float, float]:
    """
    Derive z-score threshold achieving target FDR on control population.

    Scans z-thresholds from z_min to z_max and finds the lowest threshold
    that achieves ≤ target_fdr hit rate on the control population.

    Args:
        z_scores: Series of z-scores from control population (synonymous)
        target_fdr: Target false discovery rate (default 1%)
        z_min: Minimum z-score threshold to consider
        z_max: Maximum z-score threshold to consider
        z_step: Step size for threshold search

    Returns:
        Tuple of (z_threshold, achieved_fdr)
    """
    if len(z_scores) == 0:
        return z_max, 0.0

    z_scores = z_scores.dropna()
    abs_z = z_scores.abs()

    # Scan thresholds from low to high
    for z_thresh in np.arange(z_min, z_max + z_step, z_step):
        fdr = (abs_z >= z_thresh).mean()
        if fdr <= target_fdr:
            return float(z_thresh), float(fdr)

    # If we can't achieve target FDR, return max threshold
    final_fdr = (abs_z >= z_max).mean()
    return z_max, float(final_fdr)


def calibrate_thresholds_per_condition(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    z_score_col: str = "z_score_empirical",
    variant_type_col: str = "codon_variant_type",
    target_fdr: float = 0.01,
    z_min: float = 2.0,
    z_max: float = 6.0,
    min_synonymous: int = 50,
) -> Dict[Tuple[str, str], ThresholdCalibration]:
    """
    Calibrate z-score thresholds for each (library, condition) to achieve target FDR.

    Uses synonymous variants as the null distribution. For each condition,
    finds the z-threshold where synonymous hit rate ≤ target_fdr.

    Args:
        df: DataFrame with z-scores and variant annotations
        library_col: Column with library identifier
        condition_col: Column with condition name
        z_score_col: Column with z-scores
        variant_type_col: Column with variant type
        target_fdr: Target FDR on synonymous (default 1%)
        z_min: Minimum z-threshold to consider
        z_max: Maximum z-threshold to consider
        min_synonymous: Minimum synonymous variants required

    Returns:
        Dictionary mapping (library, condition) to ThresholdCalibration
    """
    calibrations = {}

    # Create variant type masks
    syn_mask = df[variant_type_col].str.lower().str.contains("synonymous", na=False)
    mis_mask = df[variant_type_col].str.lower().str.contains("missense", na=False)

    for (library, condition), group in df.groupby([library_col, condition_col]):
        syn_group = group[syn_mask]
        mis_group = group[mis_mask]

        syn_z = syn_group[z_score_col].dropna()
        mis_z = mis_group[z_score_col].dropna()

        # Check for sufficient synonymous
        if len(syn_z) < min_synonymous:
            calibrations[(library, condition)] = ThresholdCalibration(
                condition=str(condition),
                library=str(library),
                target_fdr=target_fdr,
                achieved_fdr=0.0,
                z_threshold=z_max,  # Use conservative threshold
                syn_n=len(syn_z),
                syn_n_above_thresh=0,
                missense_n=len(mis_z),
                missense_n_above_thresh=0,
                missense_syn_ratio=0.0,
                calibration_status="WARN_LOW_N",
            )
            continue

        # Derive threshold
        z_threshold, achieved_fdr = derive_fdr_threshold(
            syn_z, target_fdr, z_min, z_max
        )

        # Calculate hit rates at this threshold
        syn_n_above = (syn_z.abs() >= z_threshold).sum()
        mis_n_above = (mis_z.abs() >= z_threshold).sum() if len(mis_z) > 0 else 0

        # Calculate missense/synonymous ratio
        if syn_n_above > 0:
            # Normalize by total counts
            syn_rate = syn_n_above / len(syn_z)
            mis_rate = mis_n_above / len(mis_z) if len(mis_z) > 0 else 0
            missense_syn_ratio = mis_rate / syn_rate if syn_rate > 0 else 0
        else:
            missense_syn_ratio = 0.0

        # Determine status
        if achieved_fdr > target_fdr * 1.5:
            status = "WARN_HIGH_FDR"
        elif achieved_fdr <= target_fdr:
            status = "OK"
        else:
            status = "OK"  # Within acceptable range

        calibrations[(library, condition)] = ThresholdCalibration(
            condition=str(condition),
            library=str(library),
            target_fdr=target_fdr,
            achieved_fdr=achieved_fdr,
            z_threshold=z_threshold,
            syn_n=len(syn_z),
            syn_n_above_thresh=int(syn_n_above),
            missense_n=len(mis_z),
            missense_n_above_thresh=int(mis_n_above),
            missense_syn_ratio=missense_syn_ratio,
            calibration_status=status,
        )

    return calibrations


def validate_ko_enrichment(
    df: pd.DataFrame,
    z_score_col: str = "z_score_empirical",
    variant_type_col: str = "codon_variant_type",
    top_percentile: float = 0.05,
) -> Dict[str, float]:
    """
    Validate that top hits are enriched for KO (knockout) variants.

    KO variants = stop_gain, frameshift, deletions
    These should be overrepresented in top hits if selection is working.

    Args:
        df: DataFrame with z-scores and variant annotations
        z_score_col: Column with z-scores
        variant_type_col: Column with variant type
        top_percentile: Fraction of top hits to consider (default 5%)

    Returns:
        Dictionary with enrichment statistics
    """
    z_scores = df[z_score_col].dropna()
    abs_z = z_scores.abs()

    if len(abs_z) == 0:
        return {"ko_enrichment": 0.0, "status": "NO_DATA"}

    # Define KO variants
    ko_types = ["stop_gain", "frameshift", "nonsense"]
    ko_mask = df[variant_type_col].str.lower().str.contains(
        "|".join(ko_types), na=False, regex=True
    )

    # Calculate KO fraction in all variants
    ko_fraction_all = ko_mask.sum() / len(df)

    if ko_fraction_all == 0:
        return {"ko_enrichment": 0.0, "status": "NO_KO_VARIANTS"}

    # Get top hits by |z|
    threshold = abs_z.quantile(1 - top_percentile)
    top_mask = abs_z >= threshold

    # Calculate KO fraction in top hits
    top_df = df[top_mask]
    ko_in_top = ko_mask[top_mask].sum()
    ko_fraction_top = ko_in_top / len(top_df) if len(top_df) > 0 else 0

    # Calculate enrichment
    enrichment = ko_fraction_top / ko_fraction_all if ko_fraction_all > 0 else 0

    # Determine status
    if enrichment >= 2.0:
        status = "GOOD"  # Strong enrichment
    elif enrichment >= 1.5:
        status = "OK"  # Modest enrichment
    elif enrichment >= 1.0:
        status = "WEAK"  # No depletion, but no enrichment
    else:
        status = "WARN"  # KO variants depleted in top hits

    return {
        "ko_fraction_all": ko_fraction_all,
        "ko_fraction_top": ko_fraction_top,
        "ko_enrichment": enrichment,
        "n_ko_all": int(ko_mask.sum()),
        "n_ko_top": int(ko_in_top),
        "n_top": len(top_df),
        "top_percentile": top_percentile,
        "status": status,
    }


def validate_missense_enrichment(
    df: pd.DataFrame,
    z_score_col: str = "z_score_empirical",
    variant_type_col: str = "codon_variant_type",
    z_thresholds: List[float] = [2.0, 3.0, 4.0, 5.0],
) -> pd.DataFrame:
    """
    Validate that missense/synonymous ratio increases with effect size.

    At higher z-score thresholds, we expect more missense relative to synonymous
    because synonymous are mostly null and missense include functional variants.

    Args:
        df: DataFrame with z-scores and variant annotations
        z_score_col: Column with z-scores
        variant_type_col: Column with variant type
        z_thresholds: Z-score thresholds to evaluate

    Returns:
        DataFrame with missense/synonymous ratios at each threshold
    """
    syn_mask = df[variant_type_col].str.lower().str.contains("synonymous", na=False)
    mis_mask = df[variant_type_col].str.lower().str.contains("missense", na=False)

    syn_z = df.loc[syn_mask, z_score_col].dropna()
    mis_z = df.loc[mis_mask, z_score_col].dropna()

    records = []
    for z_thresh in z_thresholds:
        syn_above = (syn_z.abs() >= z_thresh).sum()
        mis_above = (mis_z.abs() >= z_thresh).sum()

        syn_rate = syn_above / len(syn_z) if len(syn_z) > 0 else 0
        mis_rate = mis_above / len(mis_z) if len(mis_z) > 0 else 0

        ratio = mis_rate / syn_rate if syn_rate > 0 else float("inf")

        records.append({
            "z_threshold": z_thresh,
            "syn_n": len(syn_z),
            "syn_n_above": syn_above,
            "syn_rate": syn_rate,
            "missense_n": len(mis_z),
            "missense_n_above": mis_above,
            "missense_rate": mis_rate,
            "missense_syn_ratio": ratio,
        })

    result_df = pd.DataFrame(records)

    # Check if ratio increases with threshold (expected behavior)
    if len(result_df) >= 2:
        ratios = result_df["missense_syn_ratio"].values
        is_increasing = all(ratios[i] <= ratios[i+1] for i in range(len(ratios)-1))
        result_df["trend_valid"] = is_increasing
    else:
        result_df["trend_valid"] = True

    return result_df


def apply_calibrated_thresholds(
    df: pd.DataFrame,
    calibrations: Dict[Tuple[str, str], ThresholdCalibration],
    library_col: str = "library",
    condition_col: str = "comparison",
    z_score_col: str = "z_score_empirical",
    output_col: str = "is_significant_empirical",
) -> pd.DataFrame:
    """
    Apply calibrated z-score thresholds to call hits.

    Uses per-(library, condition) thresholds from calibration.

    Args:
        df: DataFrame with z-scores
        calibrations: Dictionary from calibrate_thresholds_per_condition()
        library_col: Column with library identifier
        condition_col: Column with condition name
        z_score_col: Column with z-scores
        output_col: Name for significance column

    Returns:
        DataFrame with added significance column
    """
    df = df.copy()

    # Create lookup key
    df["_cal_key"] = list(zip(df[library_col], df[condition_col]))

    # Map thresholds
    df["_z_thresh"] = df["_cal_key"].map(
        lambda k: calibrations[k].z_threshold if k in calibrations else 6.0
    )

    # Apply thresholds
    df[output_col] = df[z_score_col].abs() >= df["_z_thresh"]

    # Store the threshold used
    df["empirical_z_threshold"] = df["_z_thresh"]

    # Clean up
    df = df.drop(columns=["_cal_key", "_z_thresh"])

    return df


def calculate_syn_percentile(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    z_score_col: str = "z_score_empirical",
    variant_type_col: str = "codon_variant_type",
    output_col: str = "syn_percentile",
) -> pd.DataFrame:
    """
    Calculate percentile rank of each variant vs synonymous distribution.

    A variant at 99th syn_percentile means it has |z| greater than 99%
    of synonymous variants in that (library, condition).

    Args:
        df: DataFrame with z-scores
        library_col: Column with library identifier
        condition_col: Column with condition name
        z_score_col: Column with z-scores
        variant_type_col: Column with variant type
        output_col: Name for percentile column

    Returns:
        DataFrame with added syn_percentile column
    """
    df = df.copy()
    df[output_col] = np.nan

    syn_mask = df[variant_type_col].str.lower().str.contains("synonymous", na=False)

    for (library, condition), group in df.groupby([library_col, condition_col]):
        # Get synonymous distribution
        syn_group = group[syn_mask]
        syn_abs_z = syn_group[z_score_col].abs().dropna()

        if len(syn_abs_z) == 0:
            continue

        # For each variant, calculate percentile vs synonymous
        group_mask = (df[library_col] == library) & (df[condition_col] == condition)
        variant_abs_z = df.loc[group_mask, z_score_col].abs()

        # Calculate percentile (what fraction of synonymous is below this value)
        def get_percentile(z):
            if pd.isna(z):
                return np.nan
            return (syn_abs_z < z).mean() * 100

        df.loc[group_mask, output_col] = variant_abs_z.apply(get_percentile)

    return df


def generate_calibration_report(
    calibrations: Dict[Tuple[str, str], ThresholdCalibration],
) -> pd.DataFrame:
    """
    Generate a report of threshold calibrations.

    Args:
        calibrations: Dictionary from calibrate_thresholds_per_condition()

    Returns:
        DataFrame with calibration results
    """
    records = []
    for (library, condition), cal in calibrations.items():
        records.append({
            "library": cal.library,
            "condition": cal.condition,
            "target_fdr": cal.target_fdr,
            "achieved_fdr": cal.achieved_fdr,
            "z_threshold": cal.z_threshold,
            "syn_n": cal.syn_n,
            "syn_n_above_thresh": cal.syn_n_above_thresh,
            "missense_n": cal.missense_n,
            "missense_n_above_thresh": cal.missense_n_above_thresh,
            "missense_syn_ratio": cal.missense_syn_ratio,
            "calibration_status": cal.calibration_status,
        })

    report_df = pd.DataFrame(records)

    if len(report_df) > 0:
        report_df = report_df.sort_values(["library", "condition"])

    return report_df


def run_threshold_calibration_pipeline(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    z_score_col: str = "z_score_empirical",
    variant_type_col: str = "codon_variant_type",
    target_fdr: float = 0.01,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    Run the complete threshold calibration pipeline.

    Steps:
    1. Calibrate z-thresholds per (library, condition) for target FDR
    2. Apply calibrated thresholds
    3. Calculate synonymous percentiles
    4. Validate KO enrichment
    5. Validate missense/synonymous ratio trend

    Args:
        df: DataFrame with z-scores (from bias correction)
        library_col: Column with library identifier
        condition_col: Column with condition name
        z_score_col: Column with z-scores
        variant_type_col: Column with variant type
        target_fdr: Target FDR on synonymous
        verbose: Whether to print progress

    Returns:
        Tuple of (df_with_calls, calibration_report, validation_stats)
    """
    if verbose:
        print("\n" + "=" * 80)
        print("EMPIRICAL THRESHOLD CALIBRATION")
        print("=" * 80)

    # Step 1: Calibrate thresholds
    if verbose:
        print(f"\n1. Calibrating thresholds for {target_fdr*100:.1f}% FDR...")
    calibrations = calibrate_thresholds_per_condition(
        df, library_col, condition_col, z_score_col, variant_type_col, target_fdr
    )

    if verbose:
        n_cal = len(calibrations)
        n_ok = sum(1 for c in calibrations.values() if c.calibration_status == "OK")
        n_warn = n_cal - n_ok
        print(f"   Calibrated {n_cal} (library, condition) pairs")
        print(f"   {n_ok} OK, {n_warn} with warnings")

        # Show threshold distribution
        thresholds = [c.z_threshold for c in calibrations.values()]
        print(f"   Threshold range: {min(thresholds):.2f} - {max(thresholds):.2f}")
        print(f"   Threshold median: {np.median(thresholds):.2f}")

    # Step 2: Apply thresholds
    if verbose:
        print("\n2. Applying calibrated thresholds...")
    df = apply_calibrated_thresholds(
        df, calibrations, library_col, condition_col, z_score_col
    )

    n_hits = df["is_significant_empirical"].sum()
    if verbose:
        print(f"   {n_hits:,} hits called")

    # Step 3: Calculate synonymous percentiles
    if verbose:
        print("\n3. Calculating synonymous percentiles...")
    df = calculate_syn_percentile(
        df, library_col, condition_col, z_score_col, variant_type_col
    )

    # Step 4: Validate KO enrichment
    if verbose:
        print("\n4. Validating KO enrichment in top hits...")
    ko_validation = validate_ko_enrichment(df, z_score_col, variant_type_col)
    if verbose:
        print(f"   KO enrichment: {ko_validation['ko_enrichment']:.2f}x ({ko_validation['status']})")

    # Step 5: Validate missense/synonymous trend
    if verbose:
        print("\n5. Validating missense/synonymous ratio trend...")
    mis_syn_validation = validate_missense_enrichment(df, z_score_col, variant_type_col)
    if verbose:
        trend_ok = mis_syn_validation["trend_valid"].iloc[0] if len(mis_syn_validation) > 0 else False
        print(f"   Ratio increases with threshold: {trend_ok}")

    # Generate reports
    calibration_report = generate_calibration_report(calibrations)

    validation_stats = {
        "ko_enrichment": ko_validation,
        "missense_syn_validation": mis_syn_validation.to_dict("records"),
        "n_hits_total": int(n_hits),
        "target_fdr": target_fdr,
    }

    if verbose:
        print("\n" + "=" * 80)
        print("THRESHOLD CALIBRATION COMPLETE")
        print("=" * 80)

    return df, calibration_report, validation_stats
