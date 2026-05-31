"""
Hit calling functions for the BC1-From-Screen pipeline.

Implements abundance-tier-aware adaptive hit calling with concordance filtering.

Abundance tiers (for DESeq2 analysis strategy):
- abundance_high (≥100 reads): Individual BC1 analysis, z>=3.0, concordance>=0.75
- abundance_medium (20-99 reads): Oligo-aggregated, z>=3.75, concordance>=0.85

NOTE: These are distinct from confidence tiers in bc1_donor_bc0_pipeline.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
from scipy import stats

from ..config import HitCallingConfig


@dataclass
class HitCallResult:
    """Result of hit calling for a single comparison."""
    comparison: str
    library: str
    tier: str
    n_tested: int
    n_significant: int
    n_up: int
    n_down: int
    thresholds_used: Dict = field(default_factory=dict)

    @property
    def hit_rate(self) -> float:
        if self.n_tested == 0:
            return 0.0
        return self.n_significant / self.n_tested

    def to_dict(self) -> Dict:
        return {
            "comparison": self.comparison,
            "library": self.library,
            "tier": self.tier,
            "n_tested": self.n_tested,
            "n_significant": self.n_significant,
            "n_up": self.n_up,
            "n_down": self.n_down,
            "hit_rate": self.hit_rate,
            **{f"threshold_{k}": v for k, v in self.thresholds_used.items()},
        }


def calculate_z_scores(
    log2fc: pd.Series,
    method: str = "robust",
) -> pd.Series:
    """
    Calculate z-scores from log2 fold changes.

    Args:
        log2fc: Series of log2 fold changes
        method: 'robust' (MAD-based) or 'standard' (SD-based)

    Returns:
        Series of z-scores
    """
    if method == "robust":
        # MAD-based z-scores (more robust to outliers)
        median = log2fc.median()
        mad = np.median(np.abs(log2fc - median))
        # Scale MAD to be consistent with SD for normal distribution
        # 1.4826 = 1/Φ⁻¹(0.75), the consistency constant for normal data
        # See variance_calibration.py for full derivation
        mad_scaled = mad * 1.4826
        if mad_scaled == 0:
            return pd.Series(0, index=log2fc.index)
        return (log2fc - median) / mad_scaled
    else:
        # Standard z-scores
        mean = log2fc.mean()
        std = log2fc.std()
        if std == 0:
            return pd.Series(0, index=log2fc.index)
        return (log2fc - mean) / std


def calculate_concordance(
    replicate_effects: pd.DataFrame,
    direction_col: str = "direction",
) -> pd.Series:
    """
    Calculate concordance fraction (fraction of replicates agreeing on direction).

    Args:
        replicate_effects: DataFrame with replicate-level effects
        direction_col: Column indicating direction (positive/negative/neutral)

    Returns:
        Series of concordance fractions indexed by BC1
    """
    # Count direction votes per BC1
    direction_counts = replicate_effects.groupby("bc1")[direction_col].value_counts()

    # Get the dominant direction count
    total_reps = replicate_effects.groupby("bc1").size()
    max_counts = direction_counts.groupby("bc1").max()

    return max_counts / total_reps


def calculate_adaptive_base_mean_threshold(
    deseq_results: pd.DataFrame,
    percentile: float = 5.0,
    base_mean_col: str = "baseMean",
) -> float:
    """
    Calculate data-driven baseMean threshold.

    Args:
        deseq_results: DESeq2 results DataFrame
        percentile: Percentile to use for threshold
        base_mean_col: Name of baseMean column

    Returns:
        Threshold value
    """
    # Filter to significant results first
    significant = deseq_results[deseq_results["padj"] < 0.1]

    if len(significant) == 0:
        return 0.0

    return np.percentile(significant[base_mean_col], percentile)


def call_hits_single_comparison(
    deseq_results: pd.DataFrame,
    tier_mapping: pd.DataFrame,
    tier: str,
    config: HitCallingConfig,
    comparison_name: str,
    library_name: str,
    bc1_col: str = "bc1",
    lfc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    base_mean_col: str = "baseMean",
) -> Tuple[pd.DataFrame, HitCallResult]:
    """
    Call hits for a single comparison using abundance-tier-specific thresholds.

    Args:
        deseq_results: DESeq2 results DataFrame
        tier_mapping: BC1 to tier mapping
        tier: Target abundance tier ('abundance_high' or 'abundance_medium')
        config: Hit calling configuration
        comparison_name: Name of the comparison
        library_name: Name of the library
        bc1_col: Name of BC1 column
        lfc_col: Name of log2 fold change column
        padj_col: Name of adjusted p-value column
        base_mean_col: Name of baseMean column

    Returns:
        Tuple of (hits DataFrame, HitCallResult summary)
    """
    # Get tier-specific thresholds
    thresholds = config.get_thresholds(tier)
    z_threshold = thresholds["z_threshold"]
    concordance_threshold = thresholds["concordance_threshold"]

    # Filter to target tier
    tier_bc1s = tier_mapping[tier_mapping["tier"] == tier][bc1_col].tolist()
    results_filtered = deseq_results[deseq_results[bc1_col].isin(tier_bc1s)].copy()

    if len(results_filtered) == 0:
        return pd.DataFrame(), HitCallResult(
            comparison=comparison_name,
            library=library_name,
            tier=tier,
            n_tested=0,
            n_significant=0,
            n_up=0,
            n_down=0,
            thresholds_used=thresholds,
        )

    # Calculate z-scores
    results_filtered["z_score"] = calculate_z_scores(results_filtered[lfc_col])

    # Calculate adaptive baseMean threshold if enabled
    if config.use_adaptive_base_mean:
        base_mean_threshold = calculate_adaptive_base_mean_threshold(
            results_filtered,
            config.base_mean_percentile,
            base_mean_col,
        )
    else:
        base_mean_threshold = 0.0

    # Apply significance criteria
    results_filtered["is_significant"] = (
        (results_filtered[padj_col] < 0.1) &
        (np.abs(results_filtered["z_score"]) >= z_threshold) &
        (results_filtered[base_mean_col] >= base_mean_threshold)
    )

    # Determine direction
    results_filtered["direction"] = np.where(
        results_filtered[lfc_col] > 0, "up",
        np.where(results_filtered[lfc_col] < 0, "down", "neutral")
    )

    # Filter to significant hits
    hits = results_filtered[results_filtered["is_significant"]].copy()

    # Create summary
    summary = HitCallResult(
        comparison=comparison_name,
        library=library_name,
        tier=tier,
        n_tested=len(results_filtered),
        n_significant=len(hits),
        n_up=len(hits[hits["direction"] == "up"]),
        n_down=len(hits[hits["direction"] == "down"]),
        thresholds_used={
            **thresholds,
            "base_mean_threshold": base_mean_threshold,
        },
    )

    return hits, summary


def aggregate_hit_calls(
    all_hits: List[pd.DataFrame],
    all_summaries: List[HitCallResult],
    bc1_col: str = "bc1",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate hit calls across multiple comparisons.

    Args:
        all_hits: List of hit DataFrames from call_hits_single_comparison()
        all_summaries: List of HitCallResult objects
        bc1_col: Name of BC1 column

    Returns:
        Tuple of (aggregated hits, summary statistics)
    """
    if not all_hits:
        return pd.DataFrame(), pd.DataFrame()

    # Combine all hits
    combined_hits = pd.concat(all_hits, ignore_index=True)

    # Count hits per BC1
    hit_counts = combined_hits.groupby(bc1_col).agg({
        "comparison": "count",
        "direction": lambda x: (x == "up").sum(),
    }).rename(columns={
        "comparison": "n_comparisons_hit",
        "direction": "n_up",
    })
    hit_counts["n_down"] = hit_counts["n_comparisons_hit"] - hit_counts["n_up"]

    # Convert summaries to DataFrame
    summary_df = pd.DataFrame([s.to_dict() for s in all_summaries])

    return hit_counts.reset_index(), summary_df


def identify_consistent_hits(
    aggregated_hits: pd.DataFrame,
    min_comparisons: int = 2,
    min_direction_consistency: float = 0.75,
    bc1_col: str = "bc1",
) -> pd.DataFrame:
    """
    Identify BC1s with consistent hit calls across comparisons.

    Args:
        aggregated_hits: Aggregated hits from aggregate_hit_calls()
        min_comparisons: Minimum number of comparisons to be called a hit
        min_direction_consistency: Minimum fraction of comparisons with same direction
        bc1_col: Name of BC1 column

    Returns:
        DataFrame of consistent hits
    """
    if len(aggregated_hits) == 0:
        return pd.DataFrame()

    # Calculate direction consistency (protect against division by zero)
    aggregated_hits = aggregated_hits.copy()
    aggregated_hits["direction_consistency"] = np.where(
        aggregated_hits["n_comparisons_hit"] > 0,
        aggregated_hits[["n_up", "n_down"]].max(axis=1) / aggregated_hits["n_comparisons_hit"],
        0.0
    )

    # Filter to consistent hits
    consistent = aggregated_hits[
        (aggregated_hits["n_comparisons_hit"] >= min_comparisons) &
        (aggregated_hits["direction_consistency"] >= min_direction_consistency)
    ]

    # Add dominant direction
    consistent["dominant_direction"] = np.where(
        consistent["n_up"] > consistent["n_down"], "up", "down"
    )

    return consistent
