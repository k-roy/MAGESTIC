"""
Noise assessment and adaptive threshold handling for hit calling.

This module handles three distinct noise sources:
1. Noisy conditions: High synonymous MAD (condition-level noise)
2. Noisy loci: High within-locus variance (gene×condition noise)
3. SV-prone regions: Structural variant artifacts (region-level noise)

The framework applies adaptive penalties based on noise flags:
- No flags: Standard thresholds
- 1 flag: Moderate penalty (stricter thresholds)
- 2+ flags: Severe penalty (very strict thresholds)

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation


# ============================================================================
# Noise Flag Dataclass
# ============================================================================

@dataclass
class NoiseFlags:
    """Combined noise assessment for a variant."""

    # Condition-level
    noisy_condition: bool = False
    condition_mad: float = 0.0
    condition_mad_percentile: float = 50.0

    # Locus-level
    high_variance_locus: bool = False
    locus_middle50_std: float = 0.0
    locus_variance_percentile: float = 50.0

    # Region-level (SV prediction)
    is_high_sv: bool = False
    sv_prediction_score: float = 0.0

    # Variant-level (observed)
    within_variant_cv: float = 0.0

    @property
    def n_flags(self) -> int:
        """Count number of noise flags set."""
        return sum([self.noisy_condition, self.high_variance_locus, self.is_high_sv])

    @property
    def any_noise_flag(self) -> bool:
        """Check if any noise flag is set."""
        return self.n_flags > 0

    @property
    def noise_penalty(self) -> str:
        """Determine penalty level based on noise flags."""
        if self.n_flags >= 2:
            return "severe"
        elif self.n_flags == 1:
            return "moderate"
        else:
            return "none"


@dataclass
class NoisePenaltyThresholds:
    """Thresholds for different penalty levels."""
    min_bc1: int
    min_concordance: float
    max_cv: float
    confidence: str


# Default penalty configurations
PENALTY_THRESHOLDS = {
    "none": NoisePenaltyThresholds(
        min_bc1=2,
        min_concordance=0.75,
        max_cv=1.0,  # No CV restriction
        confidence="high",
    ),
    "moderate": NoisePenaltyThresholds(
        min_bc1=3,
        min_concordance=0.85,
        max_cv=0.5,
        confidence="medium",
    ),
    "severe": NoisePenaltyThresholds(
        min_bc1=4,
        min_concordance=0.90,
        max_cv=0.4,
        confidence="low",
    ),
}


# ============================================================================
# Noisy Condition Detection
# ============================================================================

def identify_noisy_conditions(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    lfc_col: str = "lfc_centered",
    variant_type_col: str = "codon_variant_type",
    mad_multiplier: float = 1.5,
) -> pd.DataFrame:
    """
    Identify conditions with high synonymous MAD (noisy conditions).

    A condition is noisy if its synonymous MAD > mad_multiplier × median_MAD.

    Args:
        df: DataFrame with LFCs
        library_col: Column with library identifier
        condition_col: Column with condition name
        lfc_col: Column with (centered) LFC values
        variant_type_col: Column with variant type
        mad_multiplier: Flag if MAD > multiplier × median

    Returns:
        DataFrame with condition noise statistics
    """
    # Filter to synonymous
    syn_mask = df[variant_type_col].str.lower().str.contains("synonymous", na=False)
    syn_df = df[syn_mask]

    # Calculate MAD per (library, condition)
    records = []
    for (library, condition), group in syn_df.groupby([library_col, condition_col]):
        lfc_values = group[lfc_col].dropna()

        if len(lfc_values) < 10:
            continue

        mad = median_abs_deviation(lfc_values, scale='normal')

        records.append({
            "library": library,
            "condition": condition,
            "syn_n": len(lfc_values),
            "syn_mad": float(mad),
        })

    result_df = pd.DataFrame(records)

    if len(result_df) == 0:
        return result_df

    # Calculate percentile and flag noisy conditions
    median_mad = result_df["syn_mad"].median()
    result_df["mad_percentile"] = result_df["syn_mad"].rank(pct=True) * 100
    result_df["is_noisy_condition"] = result_df["syn_mad"] > (mad_multiplier * median_mad)

    # Add threshold used
    result_df["mad_threshold"] = mad_multiplier * median_mad

    return result_df


# ============================================================================
# Noisy Locus Detection (Variance-Based)
# ============================================================================

def identify_noisy_loci(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    gene_col: str = "common_gene_name_locus",
    lfc_col: str = "lfc_centered",
    variance_percentile_threshold: float = 90.0,
    min_variants_per_locus: int = 4,
) -> pd.DataFrame:
    """
    Identify gene×condition pairs with high within-locus variance.

    Uses the middle 50% variance approach (from noisy_locus_v2.py):
    1. For each (gene, condition), calculate middle 50% variance
    2. Compare to genome-wide distribution
    3. Flag if variance > percentile_threshold

    Args:
        df: DataFrame with LFCs
        library_col: Column with library identifier
        condition_col: Column with condition name
        gene_col: Column with gene name
        lfc_col: Column with (centered) LFC values
        variance_percentile_threshold: Flag if variance > this percentile
        min_variants_per_locus: Minimum variants required

    Returns:
        DataFrame with locus variance statistics
    """
    records = []

    for (library, condition, gene), group in df.groupby([library_col, condition_col, gene_col]):
        if pd.isna(gene):
            continue

        lfc_values = group[lfc_col].dropna().values

        if len(lfc_values) < min_variants_per_locus:
            continue

        # Calculate middle 50%
        q25, q75 = np.percentile(lfc_values, [25, 75])
        middle_50 = lfc_values[(lfc_values >= q25) & (lfc_values <= q75)]

        if len(middle_50) < 2:
            continue

        middle_50_std = float(np.std(middle_50))

        records.append({
            "library": library,
            "condition": condition,
            "gene": gene,
            "n_variants": len(lfc_values),
            "middle_50_std": middle_50_std,
        })

    result_df = pd.DataFrame(records)

    if len(result_df) == 0:
        return result_df

    # Calculate percentile and flag high-variance loci
    result_df["variance_percentile"] = result_df["middle_50_std"].rank(pct=True) * 100
    result_df["is_high_variance_locus"] = result_df["variance_percentile"] >= variance_percentile_threshold

    return result_df


# ============================================================================
# SV-Prone Region Flagging
# ============================================================================

def flag_sv_prone_variants(
    df: pd.DataFrame,
    sv_score_col: str = "SV_prediction_score",
    sv_threshold: float = 0.15,
    output_col: str = "is_high_sv",
) -> pd.DataFrame:
    """
    Flag variants in SV-prone regions.

    Args:
        df: DataFrame with SV scores
        sv_score_col: Column with SV prediction scores
        sv_threshold: Threshold for flagging
        output_col: Name for output column

    Returns:
        DataFrame with SV flag added
    """
    df = df.copy()

    if sv_score_col in df.columns:
        df[output_col] = df[sv_score_col].fillna(0) > sv_threshold
    else:
        df[output_col] = False

    return df


# ============================================================================
# Combined Noise Assessment
# ============================================================================

def assess_variant_noise(
    df: pd.DataFrame,
    condition_noise_df: pd.DataFrame,
    locus_noise_df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    gene_col: str = "common_gene_name_locus",
    sv_score_col: str = "SV_prediction_score",
    cv_col: str = "within_variant_cv",
    sv_threshold: float = 0.15,
) -> pd.DataFrame:
    """
    Combine all noise assessments for each variant.

    Merges condition-level, locus-level, and region-level noise flags.

    Args:
        df: DataFrame with variants
        condition_noise_df: From identify_noisy_conditions()
        locus_noise_df: From identify_noisy_loci()
        library_col: Column with library identifier
        condition_col: Column with condition name
        gene_col: Column with gene name
        sv_score_col: Column with SV prediction scores
        cv_col: Column with within-variant CV (if available)
        sv_threshold: Threshold for SV flagging

    Returns:
        DataFrame with all noise flags added
    """
    df = df.copy()

    # Initialize columns
    df["noisy_condition"] = False
    df["condition_mad"] = 0.0
    df["condition_mad_percentile"] = 50.0

    df["high_variance_locus"] = False
    df["locus_middle50_std"] = 0.0
    df["locus_variance_percentile"] = 50.0

    df["is_high_sv"] = False
    df["sv_score"] = 0.0

    # Merge condition noise
    if len(condition_noise_df) > 0:
        cond_cols = ["library", "condition", "is_noisy_condition", "syn_mad", "mad_percentile"]
        cond_subset = condition_noise_df[cond_cols].rename(columns={
            "is_noisy_condition": "noisy_condition",
            "syn_mad": "condition_mad",
            "mad_percentile": "condition_mad_percentile",
        })
        df = df.merge(
            cond_subset,
            left_on=[library_col, condition_col],
            right_on=["library", "condition"],
            how="left",
            suffixes=("", "_cond"),
        )
        # Update columns from merge
        for col in ["noisy_condition", "condition_mad", "condition_mad_percentile"]:
            if f"{col}_cond" in df.columns:
                df[col] = df[f"{col}_cond"].fillna(df[col])
                df = df.drop(columns=[f"{col}_cond"])

        # Clean up duplicate key columns
        if "library_cond" in df.columns:
            df = df.drop(columns=["library_cond"])
        if "condition_cond" in df.columns:
            df = df.drop(columns=["condition_cond"])

    # Merge locus noise
    if len(locus_noise_df) > 0:
        locus_cols = ["library", "condition", "gene", "is_high_variance_locus", "middle_50_std", "variance_percentile"]
        locus_subset = locus_noise_df[locus_cols].rename(columns={
            "is_high_variance_locus": "high_variance_locus",
            "middle_50_std": "locus_middle50_std",
            "variance_percentile": "locus_variance_percentile",
        })
        df = df.merge(
            locus_subset,
            left_on=[library_col, condition_col, gene_col],
            right_on=["library", "condition", "gene"],
            how="left",
            suffixes=("", "_locus"),
        )
        # Update columns from merge
        for col in ["high_variance_locus", "locus_middle50_std", "locus_variance_percentile"]:
            if f"{col}_locus" in df.columns:
                df[col] = df[f"{col}_locus"].fillna(df[col])
                df = df.drop(columns=[f"{col}_locus"])

        # Clean up duplicate key columns
        for dup_col in ["library_locus", "condition_locus", "gene_locus"]:
            if dup_col in df.columns:
                df = df.drop(columns=[dup_col])

    # Flag SV-prone variants
    if sv_score_col in df.columns:
        df["sv_score"] = df[sv_score_col].fillna(0)
        df["is_high_sv"] = df["sv_score"] > sv_threshold
    else:
        df["is_high_sv"] = False
        df["sv_score"] = 0.0

    # Fill NaN noise flags with False/0
    df["noisy_condition"] = df["noisy_condition"].fillna(False).astype(bool)
    df["high_variance_locus"] = df["high_variance_locus"].fillna(False).astype(bool)
    df["is_high_sv"] = df["is_high_sv"].fillna(False).astype(bool)
    df["condition_mad"] = df["condition_mad"].fillna(0.0)
    df["condition_mad_percentile"] = df["condition_mad_percentile"].fillna(50.0)
    df["locus_middle50_std"] = df["locus_middle50_std"].fillna(0.0)
    df["locus_variance_percentile"] = df["locus_variance_percentile"].fillna(50.0)

    # Count noise flags
    df["n_noise_flags"] = (
        df["noisy_condition"].astype(int) +
        df["high_variance_locus"].astype(int) +
        df["is_high_sv"].astype(int)
    )

    # Determine penalty level
    df["noise_penalty"] = df["n_noise_flags"].map({
        0: "none",
        1: "moderate",
        2: "severe",
        3: "severe",
    })

    return df


# ============================================================================
# Apply Noise Penalties
# ============================================================================

def apply_noise_penalties(
    df: pd.DataFrame,
    n_bc1_col: str = "n_bc1",
    concordance_col: str = "concordance_fraction",
    cv_col: str = "within_variant_cv",
    is_significant_col: str = "is_significant_empirical",
    penalty_thresholds: Dict[str, NoisePenaltyThresholds] = None,
) -> pd.DataFrame:
    """
    Apply adaptive thresholds based on noise assessment.

    For variants with noise flags, applies stricter thresholds.
    Original significance is preserved in is_significant_pre_penalty.

    Args:
        df: DataFrame with noise flags
        n_bc1_col: Column with BC1 count
        concordance_col: Column with concordance fraction
        cv_col: Column with within-variant CV
        is_significant_col: Column with current significance
        penalty_thresholds: Custom penalty thresholds (optional)

    Returns:
        DataFrame with adjusted significance
    """
    if penalty_thresholds is None:
        penalty_thresholds = PENALTY_THRESHOLDS

    df = df.copy()

    # Save pre-penalty significance
    df["is_significant_pre_penalty"] = df[is_significant_col]

    # Initialize post-penalty significance
    df["is_significant_penalized"] = df[is_significant_col]

    # Initialize confidence
    df["confidence"] = "high"

    for penalty_level, thresholds in penalty_thresholds.items():
        mask = df["noise_penalty"] == penalty_level

        if mask.sum() == 0:
            continue

        # Check if variant passes stricter thresholds
        passes_bc1 = df.loc[mask, n_bc1_col] >= thresholds.min_bc1
        passes_concordance = df.loc[mask, concordance_col] >= thresholds.min_concordance

        if cv_col in df.columns:
            cv_values = df.loc[mask, cv_col]
            # For missing CV: pass the check (we can't assess CV, so don't penalize)
            # For present CV: check against threshold
            # Note: fillna(0) was wrong because it implied perfect precision for missing data
            passes_cv = cv_values.isna() | (cv_values <= thresholds.max_cv)
        else:
            passes_cv = True

        # Update significance
        df.loc[mask, "is_significant_penalized"] = (
            df.loc[mask, is_significant_col] &
            passes_bc1 &
            passes_concordance &
            passes_cv
        )

        # Update confidence
        df.loc[mask, "confidence"] = thresholds.confidence

    # Count penalized hits
    n_penalized = (
        df["is_significant_pre_penalty"] &
        ~df["is_significant_penalized"]
    ).sum()

    return df


# ============================================================================
# Noise QC Report
# ============================================================================

def generate_noise_qc_report(
    df: pd.DataFrame,
    condition_noise_df: pd.DataFrame,
    locus_noise_df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
) -> pd.DataFrame:
    """
    Generate QC report summarizing noise assessment.

    Args:
        df: DataFrame with noise flags
        condition_noise_df: From identify_noisy_conditions()
        locus_noise_df: From identify_noisy_loci()
        library_col: Column with library identifier
        condition_col: Column with condition name

    Returns:
        DataFrame with noise QC metrics per (library, condition)
    """
    records = []

    for (library, condition), group in df.groupby([library_col, condition_col]):
        n_total = len(group)
        n_hits = group["is_significant_penalized"].sum() if "is_significant_penalized" in group.columns else 0

        # Condition noise
        cond_match = condition_noise_df[
            (condition_noise_df["library"] == library) &
            (condition_noise_df["condition"] == condition)
        ]
        is_noisy_cond = cond_match["is_noisy_condition"].iloc[0] if len(cond_match) > 0 else False
        cond_mad = cond_match["syn_mad"].iloc[0] if len(cond_match) > 0 else 0

        # Locus noise
        locus_match = locus_noise_df[
            (locus_noise_df["library"] == library) &
            (locus_noise_df["condition"] == condition)
        ]
        n_high_var_loci = locus_match["is_high_variance_locus"].sum() if len(locus_match) > 0 else 0

        # SV counts
        n_high_sv = group["is_high_sv"].sum() if "is_high_sv" in group.columns else 0

        # Noise flag distribution
        n_with_any_flag = (group["n_noise_flags"] > 0).sum() if "n_noise_flags" in group.columns else 0
        n_with_2plus_flags = (group["n_noise_flags"] >= 2).sum() if "n_noise_flags" in group.columns else 0

        # Penalty distribution
        n_hits_pre = group["is_significant_pre_penalty"].sum() if "is_significant_pre_penalty" in group.columns else 0
        n_hits_post = group["is_significant_penalized"].sum() if "is_significant_penalized" in group.columns else 0
        n_penalized = n_hits_pre - n_hits_post

        records.append({
            "library": library,
            "condition": condition,
            "n_variants": n_total,
            "n_hits": n_hits_post,
            "syn_mad": cond_mad,
            "is_noisy_condition": is_noisy_cond,
            "n_high_variance_loci": n_high_var_loci,
            "n_high_sv_variants": n_high_sv,
            "n_with_noise_flag": n_with_any_flag,
            "n_with_2plus_flags": n_with_2plus_flags,
            "pct_with_noise_flag": 100 * n_with_any_flag / n_total if n_total > 0 else 0,
            "n_hits_penalized": n_penalized,
            "pct_hits_penalized": 100 * n_penalized / n_hits_pre if n_hits_pre > 0 else 0,
        })

    report_df = pd.DataFrame(records)

    if len(report_df) > 0:
        report_df = report_df.sort_values(["library", "condition"])

    return report_df


# ============================================================================
# Main Pipeline
# ============================================================================

def run_noise_assessment_pipeline(
    df: pd.DataFrame,
    library_col: str = "library",
    condition_col: str = "comparison",
    gene_col: str = "common_gene_name_locus",
    lfc_col: str = "lfc_centered",
    variant_type_col: str = "codon_variant_type",
    sv_score_col: str = "SV_prediction_score",
    n_bc1_col: str = "n_bc1",
    concordance_col: str = "concordance_fraction",
    cv_col: str = "within_variant_cv",
    is_significant_col: str = "is_significant_empirical",
    mad_multiplier: float = 1.5,
    variance_percentile_threshold: float = 90.0,
    sv_threshold: float = 0.15,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run the complete noise assessment pipeline.

    Steps:
    1. Identify noisy conditions (high synonymous MAD)
    2. Identify noisy loci (high within-locus variance)
    3. Flag SV-prone regions
    4. Combine noise assessments
    5. Apply adaptive penalties
    6. Generate QC report

    Args:
        df: DataFrame with variants
        library_col: Column with library identifier
        condition_col: Column with condition name
        gene_col: Column with gene name
        lfc_col: Column with (centered) LFC values
        variant_type_col: Column with variant type
        sv_score_col: Column with SV prediction scores
        n_bc1_col: Column with BC1 count
        concordance_col: Column with concordance fraction
        cv_col: Column with within-variant CV
        is_significant_col: Column with current significance
        mad_multiplier: Flag conditions with MAD > multiplier × median
        variance_percentile_threshold: Flag loci above this percentile
        sv_threshold: Threshold for SV flagging
        verbose: Whether to print progress

    Returns:
        Tuple of (df_with_noise, condition_noise_df, locus_noise_df, qc_report_df)
    """
    if verbose:
        print("\n" + "=" * 80)
        print("NOISE ASSESSMENT PIPELINE")
        print("=" * 80)

    # Step 1: Identify noisy conditions
    if verbose:
        print("\n1. Identifying noisy conditions...")
    condition_noise_df = identify_noisy_conditions(
        df, library_col, condition_col, lfc_col, variant_type_col, mad_multiplier
    )
    n_noisy_cond = condition_noise_df["is_noisy_condition"].sum() if len(condition_noise_df) > 0 else 0
    if verbose:
        print(f"   {n_noisy_cond} noisy conditions flagged")

    # Step 2: Identify noisy loci
    if verbose:
        print("\n2. Identifying noisy loci (high-variance gene×condition)...")
    locus_noise_df = identify_noisy_loci(
        df, library_col, condition_col, gene_col, lfc_col,
        variance_percentile_threshold
    )
    n_noisy_loci = locus_noise_df["is_high_variance_locus"].sum() if len(locus_noise_df) > 0 else 0
    if verbose:
        print(f"   {n_noisy_loci} high-variance loci flagged")

    # Step 3: Combine noise assessments
    if verbose:
        print("\n3. Combining noise assessments...")
    df_noise = assess_variant_noise(
        df, condition_noise_df, locus_noise_df,
        library_col, condition_col, gene_col,
        sv_score_col, cv_col, sv_threshold
    )

    n_high_sv = df_noise["is_high_sv"].sum()
    if verbose:
        print(f"   {n_high_sv} high-SV variants flagged")

    n_any_flag = (df_noise["n_noise_flags"] > 0).sum()
    if verbose:
        print(f"   {n_any_flag} variants with ≥1 noise flag ({100*n_any_flag/len(df_noise):.1f}%)")

    # Step 4: Apply penalties
    if verbose:
        print("\n4. Applying noise penalties...")

    # Check if required columns exist
    has_bc1 = n_bc1_col in df_noise.columns
    has_concordance = concordance_col in df_noise.columns
    has_significance = is_significant_col in df_noise.columns

    if has_bc1 and has_concordance and has_significance:
        df_noise = apply_noise_penalties(
            df_noise, n_bc1_col, concordance_col, cv_col, is_significant_col
        )

        n_pre = df_noise["is_significant_pre_penalty"].sum()
        n_post = df_noise["is_significant_penalized"].sum()
        n_penalized = n_pre - n_post
        if verbose:
            print(f"   {n_penalized} hits downgraded due to noise ({100*n_penalized/n_pre:.1f}% of pre-penalty hits)")
    else:
        if verbose:
            print("   Skipping penalty application (missing required columns)")
        # Use proper column check instead of DataFrame.get() which doesn't work like dict.get()
        if is_significant_col in df_noise.columns:
            df_noise["is_significant_pre_penalty"] = df_noise[is_significant_col]
            df_noise["is_significant_penalized"] = df_noise[is_significant_col]
        else:
            df_noise["is_significant_pre_penalty"] = False
            df_noise["is_significant_penalized"] = False
        df_noise["confidence"] = "high"

    # Step 5: Generate QC report
    if verbose:
        print("\n5. Generating noise QC report...")
    qc_report_df = generate_noise_qc_report(
        df_noise, condition_noise_df, locus_noise_df,
        library_col, condition_col
    )

    if verbose:
        print("\n" + "=" * 80)
        print("NOISE ASSESSMENT COMPLETE")
        print("=" * 80)

    return df_noise, condition_noise_df, locus_noise_df, qc_report_df
