#!/usr/bin/env python3
"""
Independent Hypothesis Weighting (IHW) Integration Module

This module provides FDR control using IHW (Ignatiadis et al. 2016, Nature Methods).
IHW increases discovery power by using a covariate (e.g., mean count) to learn
optimal hypothesis weights.

Since IHW is implemented in R (Bioconductor), this module uses rpy2 to interface
with the R package. If rpy2/IHW is not available, falls back to standard BH.

Prerequisites:
    pip install rpy2
    R -e 'BiocManager::install("IHW")'

Author: Kevin R. Roy

References:
    Ignatiadis N, Klaus B, Zaugg JB, Huber W (2016). Data-driven hypothesis
    weighting increases detection power in genome-scale multiple testing.
    Nature Methods 13:577-580. doi:10.1038/nmeth.3885
"""

from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)

# Check for rpy2 and IHW availability
HAS_RPY2 = False
HAS_IHW = False

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, FloatVector
    from rpy2.robjects.packages import importr
    HAS_RPY2 = True

    # Try to import IHW
    try:
        IHW = importr('IHW')
        HAS_IHW = True
        logger.debug("IHW R package loaded successfully")
    except Exception as e:
        logger.warning(f"IHW R package not available: {e}")
        logger.warning("Install with: R -e 'BiocManager::install(\"IHW\")'")
except ImportError as e:
    logger.warning(f"rpy2 not available: {e}")
    logger.warning("Install with: pip install rpy2")

# Fallback: statsmodels for BH correction
try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


@dataclass
class IHWResult:
    """Results from IHW multiple testing correction."""
    adjusted_pvalues: np.ndarray
    weights: np.ndarray
    rejections: np.ndarray  # Boolean mask of rejected hypotheses
    n_rejections: int
    alpha: float
    method: str  # "ihw", "bh_fallback", or "none"
    covariate_used: bool


def run_ihw(
    pvalues: np.ndarray,
    covariates: np.ndarray,
    alpha: float = 0.1,
    covariate_type: str = "ordinal",
    nbins: int = 10,
) -> IHWResult:
    """
    Run Independent Hypothesis Weighting for FDR control.

    Parameters
    ----------
    pvalues : array-like
        Raw p-values from hypothesis tests
    covariates : array-like
        Covariate values (e.g., mean counts) - must be independent of p-values
        under the null hypothesis
    alpha : float
        Target FDR level (default: 0.1)
    covariate_type : str
        Type of covariate: "ordinal" (continuous) or "nominal" (categorical)
    nbins : int
        Number of bins for covariate stratification (default: 10)

    Returns
    -------
    IHWResult
        Results including adjusted p-values, weights, and rejections
    """
    pvalues = np.asarray(pvalues)
    covariates = np.asarray(covariates)

    # Handle missing values
    valid_mask = ~(np.isnan(pvalues) | np.isnan(covariates))
    n_valid = valid_mask.sum()
    n_total = len(pvalues)

    if n_valid < n_total:
        logger.warning(f"Dropping {n_total - n_valid} tests with missing p-values or covariates")

    if n_valid == 0:
        logger.error("No valid p-values to test")
        return IHWResult(
            adjusted_pvalues=np.full(n_total, np.nan),
            weights=np.ones(n_total),
            rejections=np.zeros(n_total, dtype=bool),
            n_rejections=0,
            alpha=alpha,
            method="none",
            covariate_used=False
        )

    # Try IHW via rpy2
    if HAS_RPY2 and HAS_IHW:
        try:
            return _run_ihw_r(pvalues, covariates, valid_mask, alpha, covariate_type, nbins)
        except Exception as e:
            logger.warning(f"IHW failed: {e}, falling back to BH")

    # Fallback to standard BH
    if HAS_STATSMODELS:
        return _run_bh_fallback(pvalues, valid_mask, alpha)

    # No correction available
    logger.error("No multiple testing correction available (install statsmodels or rpy2+IHW)")
    return IHWResult(
        adjusted_pvalues=pvalues,
        weights=np.ones(len(pvalues)),
        rejections=pvalues < alpha,
        n_rejections=(pvalues < alpha).sum(),
        alpha=alpha,
        method="none",
        covariate_used=False
    )


def _run_ihw_r(
    pvalues: np.ndarray,
    covariates: np.ndarray,
    valid_mask: np.ndarray,
    alpha: float,
    covariate_type: str,
    nbins: int,
) -> IHWResult:
    """Run IHW using R via rpy2."""
    # Activate pandas/R conversion
    pandas2ri.activate()

    # Extract valid values
    valid_pvals = pvalues[valid_mask]
    valid_covs = covariates[valid_mask]

    # Convert to R vectors
    pvals_r = FloatVector(valid_pvals)
    covs_r = FloatVector(valid_covs)

    # Run IHW
    ihw_result = IHW.ihw(
        pvals_r,
        covs_r,
        alpha=alpha,
        covariate_type=covariate_type,
        nbins=nbins
    )

    # Extract results using IHW accessor functions
    adj_pvalues_valid = np.array(IHW.adj_pvalues(ihw_result))
    weights_valid = np.array(IHW.weights(ihw_result, levels_only=False))
    rejections_valid = np.array(IHW.rejected_hypotheses(ihw_result))

    # Map back to full arrays
    n_total = len(pvalues)
    adjusted_pvalues = np.full(n_total, np.nan)
    weights = np.ones(n_total)
    rejections = np.zeros(n_total, dtype=bool)

    adjusted_pvalues[valid_mask] = adj_pvalues_valid
    weights[valid_mask] = weights_valid
    rejections[valid_mask] = rejections_valid

    n_rejections = rejections.sum()

    logger.info(f"IHW: {n_rejections} rejections at alpha={alpha} (covariate-weighted)")

    return IHWResult(
        adjusted_pvalues=adjusted_pvalues,
        weights=weights,
        rejections=rejections,
        n_rejections=n_rejections,
        alpha=alpha,
        method="ihw",
        covariate_used=True
    )


def _run_bh_fallback(
    pvalues: np.ndarray,
    valid_mask: np.ndarray,
    alpha: float,
) -> IHWResult:
    """Run standard Benjamini-Hochberg as fallback."""
    valid_pvals = pvalues[valid_mask]

    # Run BH correction
    rejected, adj_pvals, _, _ = multipletests(valid_pvals, alpha=alpha, method='fdr_bh')

    # Map back to full arrays
    n_total = len(pvalues)
    adjusted_pvalues = np.full(n_total, np.nan)
    rejections = np.zeros(n_total, dtype=bool)

    adjusted_pvalues[valid_mask] = adj_pvals
    rejections[valid_mask] = rejected

    n_rejections = rejections.sum()

    logger.info(f"BH (fallback): {n_rejections} rejections at alpha={alpha}")

    return IHWResult(
        adjusted_pvalues=adjusted_pvalues,
        weights=np.ones(n_total),  # BH has uniform weights
        rejections=rejections,
        n_rejections=n_rejections,
        alpha=alpha,
        method="bh_fallback",
        covariate_used=False
    )


def apply_ihw_to_robust_wls_results(
    robust_wls_df: pd.DataFrame,
    pvalue_col: str = "p_value",
    covariate_col: str = "covariate",
    alpha: float = 0.1,
) -> Tuple[pd.DataFrame, IHWResult]:
    """
    Apply IHW to results from robust WLS aggregation.

    This is a convenience function for the bc1_from_screen pipeline.

    Parameters
    ----------
    robust_wls_df : pd.DataFrame
        DataFrame from aggregate_bc1s_robust_wls() containing p_value and covariate columns
    pvalue_col : str
        Column name for p-values
    covariate_col : str
        Column name for covariate (typically mean baseMean)
    alpha : float
        Target FDR level

    Returns
    -------
    tuple
        (DataFrame with adjusted p-values, IHWResult)
    """
    if pvalue_col not in robust_wls_df.columns:
        raise ValueError(f"p-value column '{pvalue_col}' not found in DataFrame")

    if covariate_col not in robust_wls_df.columns:
        logger.warning(f"Covariate column '{covariate_col}' not found, using n_bc1 as fallback")
        covariate_col = "n_bc1" if "n_bc1" in robust_wls_df.columns else None

    pvalues = robust_wls_df[pvalue_col].values

    if covariate_col:
        covariates = robust_wls_df[covariate_col].values
    else:
        # Use uniform covariate (equivalent to BH)
        covariates = np.ones(len(pvalues))
        logger.warning("No covariate available, IHW will be equivalent to BH")

    # Run IHW
    ihw_result = run_ihw(pvalues, covariates, alpha=alpha)

    # Add results to DataFrame
    result_df = robust_wls_df.copy()
    result_df['adj_pvalue_ihw'] = ihw_result.adjusted_pvalues
    result_df['ihw_weight'] = ihw_result.weights
    result_df['significant_ihw'] = ihw_result.rejections

    return result_df, ihw_result


def compare_ihw_vs_bh(
    pvalues: np.ndarray,
    covariates: np.ndarray,
    alpha: float = 0.1,
) -> dict:
    """
    Compare IHW vs standard BH for the same data.

    Useful for understanding the power gain from IHW.

    Parameters
    ----------
    pvalues : array-like
        Raw p-values
    covariates : array-like
        Covariate values
    alpha : float
        Target FDR level

    Returns
    -------
    dict
        Comparison statistics
    """
    # Run IHW
    ihw_result = run_ihw(pvalues, covariates, alpha=alpha)

    # Run BH
    valid_mask = ~np.isnan(pvalues)
    bh_result = _run_bh_fallback(pvalues, valid_mask, alpha)

    # Compare
    ihw_only = ihw_result.rejections & ~bh_result.rejections
    bh_only = bh_result.rejections & ~ihw_result.rejections
    both = ihw_result.rejections & bh_result.rejections

    comparison = {
        'n_tests': len(pvalues),
        'n_valid': valid_mask.sum(),
        'alpha': alpha,
        'ihw_rejections': ihw_result.n_rejections,
        'bh_rejections': bh_result.n_rejections,
        'ihw_only': ihw_only.sum(),
        'bh_only': bh_only.sum(),
        'both': both.sum(),
        'power_gain': ihw_result.n_rejections - bh_result.n_rejections,
        'power_gain_pct': 100 * (ihw_result.n_rejections - bh_result.n_rejections) / max(bh_result.n_rejections, 1),
        'ihw_method': ihw_result.method,
    }

    return comparison


# Availability check function
def check_ihw_availability() -> dict:
    """Check availability of IHW and its dependencies."""
    return {
        'rpy2_available': HAS_RPY2,
        'ihw_available': HAS_IHW,
        'statsmodels_available': HAS_STATSMODELS,
        'recommended_method': 'ihw' if HAS_IHW else ('bh' if HAS_STATSMODELS else 'none'),
    }
