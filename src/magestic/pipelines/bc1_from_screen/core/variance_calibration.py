#!/usr/bin/env python3
"""
Variance Calibration Module for BC1-From-Screen Pipeline

Calibrates expected variance using synonymous variants as the "neutral" reference set.
Based on the CRISPEY-BAR analysis approach (Chen et al. 2023).

Key Features:
- MAD-based outlier removal to exclude synonymous variants with real effects
- Heteroscedasticity correction: variance as function of read depth
- Per-BC1 inverse-variance weights for WLS aggregation

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import logging

logger = logging.getLogger(__name__)


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class VarianceCalibrationResult:
    """Results from variance calibration."""
    condition: str
    n_synonymous: int
    n_outliers: int
    outlier_pct: float
    fitted_params: Optional[Tuple[float, float]] = None  # (a, b) for Var = a + b/sqrt(count)
    r_squared: Optional[float] = None
    median_variance: float = 0.0
    outlier_ids: List[str] = field(default_factory=list)
    weights: Optional[np.ndarray] = None
    variance_by_variant: Optional[pd.DataFrame] = None
    calibration_method: str = "synonymous"  # "synonymous", "fallback_median", "none"


@dataclass
class WLSAggregationResult:
    """Results from WLS aggregation."""
    aggregated_df: pd.DataFrame
    n_bc1s_input: int
    n_variants_output: int
    median_effective_n: float


# =============================================================================
# Core Functions
# =============================================================================

def mad(arr: np.ndarray, scale: float = 1.4826) -> float:
    """
    Compute Median Absolute Deviation (MAD).

    MAD is a robust measure of variability less sensitive to outliers than std dev.
    The scale factor 1.4826 makes MAD consistent with std dev for normal data.

    Mathematical Derivation of 1.4826
    ---------------------------------
    For a standard normal distribution N(0,1):
    - The median of |X - median(X)| = Φ⁻¹(0.75) ≈ 0.6745
    - To convert MAD to standard deviation: σ = MAD / 0.6745 = MAD × 1.4826

    This is the "consistency constant" k = 1 / Φ⁻¹(0.75) where Φ⁻¹ is the
    inverse CDF of the standard normal distribution.

    Reference: Rousseeuw & Croux (1993), "Alternatives to the Median Absolute Deviation"

    Parameters
    ----------
    arr : array-like
        Input data
    scale : float
        Scale factor for consistency with std dev (default: 1.4826 for normal distribution)

    Returns
    -------
    float
        Scaled MAD value (comparable to standard deviation for normal data)
    """
    med = np.median(arr)
    return scale * np.median(np.abs(arr - med))


def identify_outliers_mad(
    values: np.ndarray,
    threshold: float = 3.5
) -> np.ndarray:
    """
    Identify outliers using MAD-based method.

    More robust than Z-score for distributions with outliers.
    A threshold of 3.5 MADs corresponds roughly to 3 std devs for normal data.

    Parameters
    ----------
    values : array-like
        Input values (e.g., log2 fold changes)
    threshold : float
        Number of MADs from median to consider outlier (default: 3.5)

    Returns
    -------
    np.ndarray
        Boolean mask where True indicates outlier
    """
    med = np.median(values)
    mad_value = mad(values)

    if mad_value == 0:
        # All values are identical
        return np.zeros(len(values), dtype=bool)

    mad_scores = np.abs(values - med) / mad_value
    return mad_scores > threshold


def variance_model(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """
    Model variance as a function of mean count.

    Var(log2FC) ~ a + b/sqrt(mean_count)

    This captures heteroscedasticity where low-count observations
    have higher variance due to Poisson sampling noise.

    Parameters
    ----------
    x : array-like
        Mean read counts
    a : float
        Baseline variance component
    b : float
        Count-dependent variance component

    Returns
    -------
    np.ndarray
        Predicted variance
    """
    return a + b / np.sqrt(x + 1)


def fit_variance_model(
    mean_counts: np.ndarray,
    observed_variances: np.ndarray,
    min_bin_size: int = 20
) -> Tuple[Optional[Tuple[float, float]], Optional[float]]:
    """
    Fit variance model to binned data.

    Parameters
    ----------
    mean_counts : array-like
        Mean read counts for each observation
    observed_variances : array-like
        Observed variance of log2FC for each observation
    min_bin_size : int
        Minimum number of observations per bin

    Returns
    -------
    tuple
        (fitted_params, r_squared) or (None, None) if fitting fails
    """
    # Bin by count ranges
    log_counts = np.log10(mean_counts + 1)
    n_bins = max(10, len(mean_counts) // min_bin_size)

    try:
        bins = pd.qcut(log_counts, q=n_bins, duplicates='drop')
    except ValueError:
        # Not enough unique values for binning
        return None, None

    binned_data = []
    for bin_label in bins.unique():
        mask = bins == bin_label
        if mask.sum() >= min_bin_size:
            binned_data.append({
                'mean_count': mean_counts[mask].mean(),
                'var': observed_variances[mask].mean(),
                'n': mask.sum()
            })

    binned_df = pd.DataFrame(binned_data)

    if len(binned_df) < 3:
        return None, None

    # Fit model
    try:
        popt, _ = curve_fit(
            variance_model,
            binned_df['mean_count'].values,
            binned_df['var'].values,
            p0=[0.01, 1.0],
            bounds=([0, 0], [np.inf, np.inf]),
            maxfev=10000
        )

        # Calculate R-squared
        y_pred = variance_model(binned_df['mean_count'].values, *popt)
        ss_res = np.sum((binned_df['var'].values - y_pred) ** 2)
        ss_tot = np.sum((binned_df['var'].values - binned_df['var'].mean()) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        return tuple(popt), r_squared
    except Exception as e:
        logger.warning(f"Variance model fitting failed: {e}")
        return None, None


def calibrate_variance_from_synonymous(
    bc1_data: pd.DataFrame,
    condition: str,
    lfc_col: str = "log2FoldChange",
    count_col: str = "mean_count",
    variant_col: str = "variant_id",
    variant_type_col: str = "codon_variant_type",
    mad_threshold: float = 3.5,
    min_synonymous: int = 50,
    verbose: bool = True
) -> VarianceCalibrationResult:
    """
    Calibrate variance using synonymous variants with outlier exclusion.

    Parameters
    ----------
    bc1_data : pd.DataFrame
        BC1-level data with columns for LFC, counts, variant IDs, and variant types
    condition : str
        Condition name for logging
    lfc_col : str
        Column name for log2 fold changes
    count_col : str
        Column name for mean read counts
    variant_col : str
        Column name for variant identifiers
    variant_type_col : str
        Column name for variant type classification
    mad_threshold : float
        MAD threshold for outlier detection (default: 3.5)
    min_synonymous : int
        Minimum number of synonymous variants required

    Returns
    -------
    VarianceCalibrationResult
        Calibration results including fitted parameters and weights
    """
    if verbose:
        logger.info(f"Calibrating variance for condition: {condition}")

    # Check required columns
    required_cols = [lfc_col, count_col, variant_col]
    missing = [c for c in required_cols if c not in bc1_data.columns]
    if missing:
        logger.warning(f"Missing required columns: {missing}")
        return VarianceCalibrationResult(
            condition=condition,
            n_synonymous=0,
            n_outliers=0,
            outlier_pct=0.0,
            calibration_method="none"
        )

    # Filter to synonymous variants
    if variant_type_col in bc1_data.columns:
        syn_mask = bc1_data[variant_type_col].str.contains('synonymous', case=False, na=False)
        syn_data = bc1_data[syn_mask].copy()
    else:
        logger.warning(f"Variant type column '{variant_type_col}' not found, skipping calibration")
        return VarianceCalibrationResult(
            condition=condition,
            n_synonymous=0,
            n_outliers=0,
            outlier_pct=0.0,
            calibration_method="none"
        )

    if verbose:
        logger.info(f"  Total BC1s: {len(bc1_data)}, Synonymous BC1s: {len(syn_data)}")

    if len(syn_data) < min_synonymous:
        logger.warning(f"Too few synonymous variants ({len(syn_data)} < {min_synonymous})")
        return VarianceCalibrationResult(
            condition=condition,
            n_synonymous=len(syn_data),
            n_outliers=0,
            outlier_pct=0.0,
            calibration_method="none"
        )

    # Step 1: Identify outliers among synonymous variants
    log2fc_values = syn_data[lfc_col].values
    outlier_mask = identify_outliers_mad(log2fc_values, threshold=mad_threshold)

    n_outliers = outlier_mask.sum()
    outlier_pct = 100 * n_outliers / len(syn_data)

    if verbose:
        logger.info(f"  Outliers (MAD>{mad_threshold}): {n_outliers} ({outlier_pct:.1f}%)")

    # Get outlier IDs
    outlier_ids = syn_data[outlier_mask][variant_col].tolist() if variant_col in syn_data.columns else []

    # Step 2: Use non-outlier synonymous for variance calibration
    neutral_data = syn_data[~outlier_mask]

    if len(neutral_data) < 20:
        logger.warning(f"Too few neutral synonymous after outlier removal ({len(neutral_data)})")
        median_var = syn_data[lfc_col].var() if len(syn_data) > 1 else 0.1
        return VarianceCalibrationResult(
            condition=condition,
            n_synonymous=len(syn_data),
            n_outliers=n_outliers,
            outlier_pct=outlier_pct,
            median_variance=median_var,
            outlier_ids=outlier_ids,
            calibration_method="fallback_median"
        )

    # Step 3: Calculate variance by variant (need multiple BC1s per variant)
    variance_by_variant = neutral_data.groupby(variant_col).agg({
        lfc_col: ['var', 'count', 'mean'],
        count_col: 'mean'
    }).reset_index()
    variance_by_variant.columns = [variant_col, 'var', 'n_bc1', 'mean_log2fc', 'mean_count']

    # Need at least 2 BC1s per variant for variance
    variance_by_variant = variance_by_variant[variance_by_variant['n_bc1'] >= 2]

    if len(variance_by_variant) < 20:
        logger.warning(f"Too few variants with multiple BC1s ({len(variance_by_variant)})")
        median_var = variance_by_variant['var'].median() if len(variance_by_variant) > 0 else 0.1
        return VarianceCalibrationResult(
            condition=condition,
            n_synonymous=len(syn_data),
            n_outliers=n_outliers,
            outlier_pct=outlier_pct,
            median_variance=median_var,
            outlier_ids=outlier_ids,
            variance_by_variant=variance_by_variant,
            calibration_method="fallback_median"
        )

    median_var = variance_by_variant['var'].median()

    if verbose:
        logger.info(f"  Variants with >=2 BC1s: {len(variance_by_variant)}")
        logger.info(f"  Median variance: {median_var:.4f}")

    # Step 4: Fit variance model
    popt, r_squared = fit_variance_model(
        variance_by_variant['mean_count'].values,
        variance_by_variant['var'].values
    )

    if popt is None:
        logger.warning("Could not fit variance model, using median variance")
        return VarianceCalibrationResult(
            condition=condition,
            n_synonymous=len(syn_data),
            n_outliers=n_outliers,
            outlier_pct=outlier_pct,
            median_variance=median_var,
            outlier_ids=outlier_ids,
            variance_by_variant=variance_by_variant,
            calibration_method="fallback_median"
        )

    if verbose:
        logger.info(f"  Variance model: Var = {popt[0]:.4f} + {popt[1]:.4f}/sqrt(count)")
        logger.info(f"  R-squared: {r_squared:.3f}")

    # Step 5: Calculate weights for all BC1s based on predicted variance
    all_predicted_var = variance_model(bc1_data[count_col].values, *popt)
    weights = 1 / np.maximum(all_predicted_var, 1e-6)  # Inverse variance weighting

    # Normalize weights to sum to n (preserves effective sample size interpretation)
    weights = weights * len(weights) / weights.sum()

    return VarianceCalibrationResult(
        condition=condition,
        n_synonymous=len(syn_data),
        n_outliers=n_outliers,
        outlier_pct=outlier_pct,
        fitted_params=popt,
        r_squared=r_squared,
        median_variance=median_var,
        outlier_ids=outlier_ids,
        weights=weights,
        variance_by_variant=variance_by_variant,
        calibration_method="synonymous"
    )


def aggregate_bc1s_with_wls(
    bc1_data: pd.DataFrame,
    weights: np.ndarray,
    variant_col: str = "variant_id",
    lfc_col: str = "log2FoldChange",
    additional_cols: Optional[List[str]] = None
) -> WLSAggregationResult:
    """
    Aggregate BC1-level estimates to variant-level using Weighted Least Squares.

    Parameters
    ----------
    bc1_data : pd.DataFrame
        BC1-level data with log2fc values
    weights : array-like
        Per-BC1 weights (inverse variance)
    variant_col : str
        Column to group by
    lfc_col : str
        Column with log2 fold changes
    additional_cols : list, optional
        Additional columns to include (first value per group)

    Returns
    -------
    WLSAggregationResult
        Aggregated results
    """
    bc1_data = bc1_data.copy()
    bc1_data['_weight'] = weights

    def weighted_stats(group):
        w = group['_weight'].values
        x = group[lfc_col].values

        # Handle edge cases
        if len(x) == 0 or np.sum(w) == 0:
            return pd.Series({
                'log2fc_wls': np.nan,
                'se_wls': np.nan,
                'n_bc1': 0,
                'n_eff': 0.0,
                'total_weight': 0.0
            })

        # Weighted mean
        weighted_mean = np.sum(w * x) / np.sum(w)

        # Effective sample size
        n_eff = np.sum(w)**2 / np.sum(w**2) if np.sum(w**2) > 0 else 1.0

        # Weighted variance and SE
        if len(x) > 1:
            weighted_var = np.sum(w * (x - weighted_mean)**2) / np.sum(w)
            weighted_se = np.sqrt(weighted_var / n_eff)
        else:
            weighted_var = 0.0
            weighted_se = np.nan

        result = {
            'log2fc_wls': weighted_mean,
            'se_wls': weighted_se,
            'n_bc1': len(group),
            'n_eff': n_eff,
            'total_weight': np.sum(w)
        }

        return pd.Series(result)

    aggregated = bc1_data.groupby(variant_col).apply(weighted_stats).reset_index()

    # Add additional columns (first value per group)
    if additional_cols:
        for col in additional_cols:
            if col in bc1_data.columns:
                first_vals = bc1_data.groupby(variant_col)[col].first().reset_index()
                aggregated = aggregated.merge(first_vals, on=variant_col, how='left')

    return WLSAggregationResult(
        aggregated_df=aggregated,
        n_bc1s_input=len(bc1_data),
        n_variants_output=len(aggregated),
        median_effective_n=aggregated['n_eff'].median()
    )


def run_variance_calibration_for_comparison(
    bc1_data: pd.DataFrame,
    comparison: str,
    lfc_col: str = "log2FoldChange",
    count_col: str = "baseMean",
    variant_col: str = "oligo_name",
    variant_type_col: str = "codon_variant_type",
    mad_threshold: float = 3.5,
    min_synonymous: int = 50,
    apply_wls: bool = True,
    verbose: bool = True
) -> Tuple[VarianceCalibrationResult, Optional[WLSAggregationResult]]:
    """
    Run full variance calibration workflow for a single comparison.

    This is the main entry point for pipeline integration.

    Parameters
    ----------
    bc1_data : pd.DataFrame
        BC1-level DESeq2 results
    comparison : str
        Comparison name (e.g., "Caffeine_vs_t0")
    lfc_col : str
        Column for log2 fold changes
    count_col : str
        Column for mean counts (baseMean from DESeq2)
    variant_col : str
        Column for variant/oligo identifier
    variant_type_col : str
        Column for variant type classification
    mad_threshold : float
        MAD threshold for synonymous outlier detection
    min_synonymous : int
        Minimum synonymous variants required
    apply_wls : bool
        Whether to apply WLS aggregation
    verbose : bool
        Whether to log progress

    Returns
    -------
    tuple
        (VarianceCalibrationResult, WLSAggregationResult or None)
    """
    # Run calibration
    cal_result = calibrate_variance_from_synonymous(
        bc1_data=bc1_data,
        condition=comparison,
        lfc_col=lfc_col,
        count_col=count_col,
        variant_col=variant_col,
        variant_type_col=variant_type_col,
        mad_threshold=mad_threshold,
        min_synonymous=min_synonymous,
        verbose=verbose
    )

    # Apply WLS if calibration succeeded and requested
    wls_result = None
    if apply_wls and cal_result.weights is not None:
        wls_result = aggregate_bc1s_with_wls(
            bc1_data=bc1_data,
            weights=cal_result.weights,
            variant_col=variant_col,
            lfc_col=lfc_col,
            additional_cols=[variant_type_col, 'gene', 'common_gene_name_locus']
        )
        if verbose:
            logger.info(f"  WLS: {wls_result.n_bc1s_input} BC1s -> {wls_result.n_variants_output} variants")
            logger.info(f"  Median effective n: {wls_result.median_effective_n:.1f}")

    return cal_result, wls_result


# =============================================================================
# Plotting Functions (optional, for diagnostics)
# =============================================================================

def plot_variance_calibration(
    cal_result: VarianceCalibrationResult,
    output_path: Optional[str] = None
) -> None:
    """
    Create diagnostic plot for variance calibration.

    Parameters
    ----------
    cal_result : VarianceCalibrationResult
        Calibration results
    output_path : str, optional
        Path to save plot (if None, displays interactively)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return

    if cal_result.variance_by_variant is None or cal_result.fitted_params is None:
        logger.warning("Cannot plot: missing variance data or fitted parameters")
        return

    var_df = cal_result.variance_by_variant
    popt = cal_result.fitted_params

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Variance vs mean count
    ax1 = axes[0]
    ax1.scatter(var_df['mean_count'], var_df['var'],
               alpha=0.5, s=20, label='Observed')

    # Fitted curve
    x_fit = np.logspace(0, np.log10(var_df['mean_count'].max() + 1), 100)
    y_fit = variance_model(x_fit, *popt)
    ax1.plot(x_fit, y_fit, 'r-', linewidth=2,
            label=f'Fitted: {popt[0]:.3f} + {popt[1]:.3f}/sqrt(n)')

    ax1.set_xlabel('Mean Read Count')
    ax1.set_ylabel('Variance of log2FC')
    ax1.set_title(f'{cal_result.condition}: Variance vs Count')
    ax1.set_xscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Weight distribution
    ax2 = axes[1]
    if cal_result.weights is not None:
        ax2.hist(cal_result.weights, bins=50, alpha=0.7, edgecolor='black')
        ax2.axvline(x=1.0, color='red', linestyle='--', label='Mean weight')
        ax2.set_xlabel('Weight')
        ax2.set_ylabel('Count')
        ax2.set_title(f'{cal_result.condition}: Weight Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        plt.close()
        logger.info(f"Saved plot: {output_path}")
    else:
        plt.show()
