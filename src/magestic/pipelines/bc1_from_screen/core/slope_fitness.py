"""
Slope-Based Multi-Timepoint Fitness Estimation
===============================================

This module implements fitness estimation from barcode counts across multiple
timepoints using linear regression, inspired by FitSeq and CRISPEY methodologies.

Core Concept:
    log2(relative_abundance) = β₀ + β₁ × generation + ε

Where β₁ (slope) represents the fitness (log2 fold change per generation).

Key Features:
    1. Weighted regression using inverse variance weights
    2. Linearity assessment via R² and residual analysis
    3. Detection of non-linear effects (curvature)
    4. Robust aggregation from BC1-level to variant-level

References:
    - FitSeq: Li et al. (2018) Nat Methods 15(4):267-270
    - PyFitSeq: https://github.com/FangfeiLi05/PyFitSeq
    - CRISPEY epistasis: Mullis et al. (2023) Cell 186(6):1307

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.robust.norms import HuberT


@dataclass
class SlopeFitnessResult:
    """Result of slope fitness calculation for a single BC1 or variant."""
    slope: float                    # β₁: log2 fold change per generation
    slope_se: float                 # Standard error of slope
    intercept: float                # β₀: initial abundance offset
    r_squared: float                # Model fit quality (0-1)
    z_score: float                  # slope / slope_se
    p_value: float                  # Two-sided p-value for slope ≠ 0
    n_timepoints: int               # Number of timepoints used
    generations: List[int]          # Generation values used
    residual_std: float             # Std of residuals (measures noise)
    linearity_score: float          # 1 - (curvature_term / linear_term)
    is_monotonic: bool              # Whether abundance changes monotonically
    endpoint_lfc: float             # Traditional t0 vs endpoint log2FC for comparison


@dataclass
class VariantSlopeFitness:
    """Aggregated slope fitness for a variant across multiple BC1s."""
    variant_id: str
    oligo_id: str
    slope: float                    # Weighted mean slope across BC1s
    slope_se: float                 # SE of the weighted mean
    z_score: float                  # slope / slope_se
    p_value: float                  # p-value for slope ≠ 0
    n_bc1s: int                     # Number of BC1s used
    mean_r_squared: float           # Mean R² across BC1s
    mean_linearity: float           # Mean linearity score
    concordance: float              # Fraction of BC1s with same sign
    bc1_slopes: List[float]         # Individual BC1 slopes
    bc1_weights: List[float]        # Weights used for each BC1
    endpoint_lfc: float             # Weighted mean endpoint LFC


def calculate_bc1_slope_fitness(
    counts: Dict[int, float],
    total_counts: Dict[int, float],
    variance_model: Optional[Tuple[float, float]] = None,
    min_timepoints: int = 3
) -> Optional[SlopeFitnessResult]:
    """
    Calculate slope fitness for a single BC1 across timepoints.

    Parameters
    ----------
    counts : Dict[int, float]
        BC1 counts at each generation, e.g., {0: 100, 4: 80, 8: 60, 12: 50, 16: 40, 20: 30}
    total_counts : Dict[int, float]
        Total library counts at each generation (for normalization)
    variance_model : Optional[Tuple[float, float]]
        Variance model (a, b) where Var = a + b/sqrt(count). If None, use unweighted regression.
    min_timepoints : int
        Minimum number of timepoints required for fitting

    Returns
    -------
    SlopeFitnessResult or None if insufficient data
    """
    # Filter to timepoints with counts
    valid_gens = sorted([g for g, c in counts.items() if c is not None and c > 0])

    if len(valid_gens) < min_timepoints:
        return None

    # Calculate log2 relative abundance at each timepoint
    # relative_abundance = (bc1_count / total_count) normalized to t0
    t0 = min(valid_gens)
    t0_freq = counts[t0] / total_counts[t0] if total_counts[t0] > 0 else 0

    if t0_freq == 0:
        return None

    generations = []
    log2_rel_abund = []
    weights = []

    for gen in valid_gens:
        freq = counts[gen] / total_counts[gen] if total_counts[gen] > 0 else 0
        if freq > 0:
            rel_abund = freq / t0_freq  # Normalize to t0
            generations.append(gen)
            log2_rel_abund.append(np.log2(rel_abund))

            # Calculate weight from variance model
            if variance_model is not None:
                a, b = variance_model
                var = a + b / np.sqrt(counts[gen]) if counts[gen] > 0 else a
                weights.append(1 / var if var > 0 else 1)
            else:
                weights.append(1)

    if len(generations) < min_timepoints:
        return None

    generations = np.array(generations)
    log2_rel_abund = np.array(log2_rel_abund)
    weights = np.array(weights)

    # Weighted linear regression: log2(rel_abund) ~ generation
    X = sm.add_constant(generations)

    try:
        if variance_model is not None:
            # WLS regression
            model = sm.WLS(log2_rel_abund, X, weights=weights)
        else:
            # OLS regression
            model = sm.OLS(log2_rel_abund, X)

        results = model.fit()

        intercept = results.params[0]
        slope = results.params[1]
        slope_se = results.bse[1]
        r_squared = results.rsquared
        p_value = results.pvalues[1]

        # Z-score
        z_score = slope / slope_se if slope_se > 0 else 0

        # Residual analysis
        residuals = results.resid
        residual_std = np.std(residuals)

        # Linearity assessment: fit quadratic and compare
        # linearity_score = 1 means perfectly linear, <1 means curvature present
        linearity_score = _calculate_linearity_score(generations, log2_rel_abund, weights)

        # Check monotonicity
        is_monotonic = _check_monotonicity(log2_rel_abund)

        # Calculate endpoint LFC for comparison
        endpoint_lfc = log2_rel_abund[-1] if len(log2_rel_abund) > 0 else 0

        return SlopeFitnessResult(
            slope=slope,
            slope_se=slope_se,
            intercept=intercept,
            r_squared=r_squared,
            z_score=z_score,
            p_value=p_value,
            n_timepoints=len(generations),
            generations=list(generations),
            residual_std=residual_std,
            linearity_score=linearity_score,
            is_monotonic=is_monotonic,
            endpoint_lfc=endpoint_lfc
        )

    except Exception as e:
        # Regression failed (e.g., singular matrix)
        return None


def _calculate_linearity_score(
    generations: np.ndarray,
    log2_rel_abund: np.ndarray,
    weights: np.ndarray
) -> float:
    """
    Calculate linearity score by comparing linear and quadratic fits.

    Returns 1.0 for perfectly linear, <1.0 if curvature is significant.
    """
    if len(generations) < 4:
        return 1.0  # Can't assess curvature with few points

    try:
        # Linear fit
        X_linear = sm.add_constant(generations)
        linear_model = sm.WLS(log2_rel_abund, X_linear, weights=weights)
        linear_results = linear_model.fit()
        linear_sse = np.sum(weights * linear_results.resid**2)

        # Quadratic fit
        X_quad = sm.add_constant(np.column_stack([generations, generations**2]))
        quad_model = sm.WLS(log2_rel_abund, X_quad, weights=weights)
        quad_results = quad_model.fit()
        quad_sse = np.sum(weights * quad_results.resid**2)

        # If quadratic doesn't improve much, data is linear
        if linear_sse == 0:
            return 1.0

        improvement = (linear_sse - quad_sse) / linear_sse
        linearity_score = 1.0 - improvement

        return max(0.0, min(1.0, linearity_score))

    except Exception:
        return 1.0


def _check_monotonicity(values: np.ndarray) -> bool:
    """Check if values are monotonically increasing or decreasing."""
    if len(values) < 2:
        return True

    diffs = np.diff(values)

    # All non-negative or all non-positive
    return np.all(diffs >= 0) or np.all(diffs <= 0)


def aggregate_bc1_slopes_to_variant(
    bc1_results: List[SlopeFitnessResult],
    variant_id: str,
    oligo_id: str,
    variance_model: Optional[Tuple[float, float]] = None,
    min_bc1s: int = 2
) -> Optional[VariantSlopeFitness]:
    """
    Aggregate BC1-level slope fitness to variant-level using robust WLS.

    Uses Huber M-estimation to downweight outlier BC1s (e.g., double-edits).

    Parameters
    ----------
    bc1_results : List[SlopeFitnessResult]
        Slope fitness results for each BC1 of this variant
    variant_id : str
        Variant identifier
    oligo_id : str
        Oligo identifier
    variance_model : Optional[Tuple[float, float]]
        Variance model (a, b) for weighting
    min_bc1s : int
        Minimum BC1s required for aggregation

    Returns
    -------
    VariantSlopeFitness or None if insufficient BC1s
    """
    # Filter to valid results
    valid_results = [r for r in bc1_results if r is not None and np.isfinite(r.slope)]

    if len(valid_results) < min_bc1s:
        return None

    slopes = np.array([r.slope for r in valid_results])
    slope_ses = np.array([r.slope_se for r in valid_results])
    r_squareds = np.array([r.r_squared for r in valid_results])
    linearities = np.array([r.linearity_score for r in valid_results])
    endpoint_lfcs = np.array([r.endpoint_lfc for r in valid_results])

    # Calculate weights based on precision (inverse variance)
    weights = 1 / (slope_ses**2 + 1e-10)  # Add small constant to prevent division by zero
    weights = weights / weights.sum()

    if len(valid_results) >= 3:
        try:
            # Robust WLS with Huber M-estimation
            X = np.ones((len(slopes), 1))
            # Pre-whiten
            sqrt_weights = np.sqrt(weights * len(slopes))
            y_weighted = slopes * sqrt_weights
            X_weighted = X * sqrt_weights[:, None]

            rlm_model = sm.RLM(y_weighted, X_weighted, M=HuberT())
            rlm_results = rlm_model.fit()

            aggregated_slope = rlm_results.params[0]
            aggregated_se = rlm_results.bse[0]

            # Get weights from M-estimation for output
            m_weights = rlm_results.weights

        except Exception:
            # Fall back to weighted mean
            aggregated_slope = np.average(slopes, weights=weights)
            aggregated_se = np.sqrt(1 / np.sum(1 / (slope_ses**2 + 1e-10)))
            m_weights = weights
    else:
        # Simple weighted mean for 2 BC1s
        aggregated_slope = np.average(slopes, weights=weights)
        aggregated_se = np.sqrt(1 / np.sum(1 / (slope_ses**2 + 1e-10)))
        m_weights = weights

    # Calculate z-score and p-value
    z_score = aggregated_slope / aggregated_se if aggregated_se > 0 else 0
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

    # Concordance: fraction of BC1s with same sign as aggregate
    if aggregated_slope != 0:
        same_sign = np.sum(np.sign(slopes) == np.sign(aggregated_slope))
        concordance = same_sign / len(slopes)
    else:
        concordance = 0.5

    # Aggregate endpoint LFC
    aggregated_endpoint_lfc = np.average(endpoint_lfcs, weights=weights)

    return VariantSlopeFitness(
        variant_id=variant_id,
        oligo_id=oligo_id,
        slope=aggregated_slope,
        slope_se=aggregated_se,
        z_score=z_score,
        p_value=p_value,
        n_bc1s=len(valid_results),
        mean_r_squared=np.mean(r_squareds),
        mean_linearity=np.mean(linearities),
        concordance=concordance,
        bc1_slopes=slopes.tolist(),
        bc1_weights=list(m_weights),
        endpoint_lfc=aggregated_endpoint_lfc
    )


def calculate_slope_fitness_dataframe(
    count_matrix: pd.DataFrame,
    generation_mapping: Dict[str, int],
    variance_model: Optional[Tuple[float, float]] = None,
    bc1_to_oligo: Optional[Dict[str, str]] = None,
    oligo_to_variant: Optional[Dict[str, str]] = None,
    min_timepoints: int = 3,
    min_bc1s: int = 2
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate slope fitness for all BC1s and aggregate to variants.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Count matrix with BC1s as rows and samples as columns.
        Sample names should map to generations via generation_mapping.
    generation_mapping : Dict[str, int]
        Maps sample names to generation numbers, e.g.,
        {"YPD_t0_rep1": 0, "YPD_t4_rep1": 4, "YPD_t8_rep1": 8, ...}
    variance_model : Optional[Tuple[float, float]]
        Variance model (a, b) for weighting
    bc1_to_oligo : Optional[Dict[str, str]]
        Maps BC1 to oligo ID
    oligo_to_variant : Optional[Dict[str, str]]
        Maps oligo to variant ID
    min_timepoints : int
        Minimum timepoints for slope calculation
    min_bc1s : int
        Minimum BC1s for variant aggregation

    Returns
    -------
    bc1_fitness : pd.DataFrame
        BC1-level slope fitness results
    variant_fitness : pd.DataFrame
        Variant-level aggregated slope fitness results
    """
    # Calculate total counts per sample
    sample_totals = count_matrix.sum(axis=0)

    # Identify unique generations and their samples
    gen_samples = {}  # generation -> list of sample names
    for sample, gen in generation_mapping.items():
        if sample in count_matrix.columns:
            if gen not in gen_samples:
                gen_samples[gen] = []
            gen_samples[gen].append(sample)

    bc1_results = []

    for bc1 in count_matrix.index:
        # Aggregate counts across replicates for each generation
        bc1_counts = {}
        total_counts = {}

        for gen, samples in sorted(gen_samples.items()):
            bc1_count = count_matrix.loc[bc1, samples].sum()
            total_count = sample_totals[samples].sum()

            if bc1_count > 0:  # Only include non-zero generations
                bc1_counts[gen] = bc1_count
                total_counts[gen] = total_count

        # Calculate slope fitness
        result = calculate_bc1_slope_fitness(
            counts=bc1_counts,
            total_counts=total_counts,
            variance_model=variance_model,
            min_timepoints=min_timepoints
        )

        if result is not None:
            bc1_results.append({
                'bc1': bc1,
                'oligo_id': bc1_to_oligo.get(bc1, bc1) if bc1_to_oligo else bc1,
                'variant_id': oligo_to_variant.get(bc1_to_oligo.get(bc1, bc1), 'unknown') if oligo_to_variant and bc1_to_oligo else 'unknown',
                'slope': result.slope,
                'slope_se': result.slope_se,
                'intercept': result.intercept,
                'r_squared': result.r_squared,
                'z_score': result.z_score,
                'p_value': result.p_value,
                'n_timepoints': result.n_timepoints,
                'residual_std': result.residual_std,
                'linearity_score': result.linearity_score,
                'is_monotonic': result.is_monotonic,
                'endpoint_lfc': result.endpoint_lfc
            })

    bc1_fitness = pd.DataFrame(bc1_results)

    # Aggregate to variant level
    if len(bc1_fitness) == 0:
        return bc1_fitness, pd.DataFrame()

    variant_results = []

    for (variant_id, oligo_id), group in bc1_fitness.groupby(['variant_id', 'oligo_id']):
        # Convert group rows back to SlopeFitnessResult objects using itertuples (5x faster)
        bc1_slope_results = []
        for row in group.itertuples(index=False):
            bc1_slope_results.append(SlopeFitnessResult(
                slope=row.slope,
                slope_se=row.slope_se,
                intercept=row.intercept,
                r_squared=row.r_squared,
                z_score=row.z_score,
                p_value=row.p_value,
                n_timepoints=row.n_timepoints,
                generations=[],  # Not stored in DataFrame
                residual_std=row.residual_std,
                linearity_score=row.linearity_score,
                is_monotonic=row.is_monotonic,
                endpoint_lfc=row.endpoint_lfc
            ))

        variant_result = aggregate_bc1_slopes_to_variant(
            bc1_results=bc1_slope_results,
            variant_id=variant_id,
            oligo_id=oligo_id,
            variance_model=variance_model,
            min_bc1s=min_bc1s
        )

        if variant_result is not None:
            variant_results.append({
                'variant_id': variant_result.variant_id,
                'oligo_id': variant_result.oligo_id,
                'slope_fitness': variant_result.slope,
                'slope_se': variant_result.slope_se,
                'z_score': variant_result.z_score,
                'p_value': variant_result.p_value,
                'n_bc1s': variant_result.n_bc1s,
                'mean_r_squared': variant_result.mean_r_squared,
                'mean_linearity': variant_result.mean_linearity,
                'concordance': variant_result.concordance,
                'endpoint_lfc': variant_result.endpoint_lfc
            })

    variant_fitness = pd.DataFrame(variant_results)

    return bc1_fitness, variant_fitness


def compare_slope_vs_endpoint_fitness(
    slope_results: pd.DataFrame,
    endpoint_results: pd.DataFrame,
    merge_on: str = 'variant_id'
) -> pd.DataFrame:
    """
    Compare slope-based fitness with traditional endpoint-based fitness.

    Returns merged DataFrame with correlation and concordance metrics.
    """
    merged = slope_results.merge(
        endpoint_results,
        on=merge_on,
        suffixes=('_slope', '_endpoint')
    )

    if len(merged) == 0:
        return merged

    # Calculate correlation
    corr = merged['slope_fitness'].corr(merged.get('log2FoldChange', merged.get('lfc', pd.Series())))

    # Calculate concordance (same sign)
    if 'log2FoldChange' in merged.columns or 'lfc' in merged.columns:
        endpoint_col = 'log2FoldChange' if 'log2FoldChange' in merged.columns else 'lfc'
        merged['sign_concordant'] = np.sign(merged['slope_fitness']) == np.sign(merged[endpoint_col])

    return merged
