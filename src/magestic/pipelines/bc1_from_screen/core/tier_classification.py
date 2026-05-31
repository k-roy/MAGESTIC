# -*- coding: utf-8 -*-
"""
Tier classification functions for the BC1-From-Screen pipeline.

Classifies BC1 barcodes into abundance tiers based on t0 read counts within each library.

Abundance tiers (for DESeq2 analysis strategy):
- abundance_high (>=100 reads): Individual BC1 analysis in DESeq2
- abundance_medium (20-99 reads): Aggregate by oligo before DESeq2
- abundance_low (<20 reads): Excluded from analysis (low confidence)

NOTE: These "abundance tiers" are distinct from the "confidence_tier" column
in bc1_donor_bc0_pipeline which classifies oligo attribution confidence.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np

from ..config import TierConfig, LibraryConfig


@dataclass
class TierSummary:
    """Summary statistics for abundance tier classification."""
    library: str
    abundance_high_count: int
    abundance_medium_count: int
    abundance_low_count: int
    abundance_high_reads: int
    abundance_medium_reads: int
    abundance_low_reads: int

    @property
    def total_bc1s(self) -> int:
        return self.abundance_high_count + self.abundance_medium_count + self.abundance_low_count

    @property
    def total_reads(self) -> int:
        return self.abundance_high_reads + self.abundance_medium_reads + self.abundance_low_reads

    def to_dict(self) -> Dict:
        return {
            "library": self.library,
            "abundance_high_count": self.abundance_high_count,
            "abundance_medium_count": self.abundance_medium_count,
            "abundance_low_count": self.abundance_low_count,
            "abundance_high_reads": self.abundance_high_reads,
            "abundance_medium_reads": self.abundance_medium_reads,
            "abundance_low_reads": self.abundance_low_reads,
            "total_bc1s": self.total_bc1s,
            "total_reads": self.total_reads,
        }


def calculate_t0_counts(
    long_table: pd.DataFrame,
    keyfile: pd.DataFrame,
    library_config: LibraryConfig,
    bc1_col: str = "bc1",
    sample_col: str = "sample_name",
    count_col: str = "sample_counts",
) -> pd.DataFrame:
    """
    Calculate total t0 counts per BC1 for a library.

    Args:
        long_table: Long-format BC1 counts table
        keyfile: Sample metadata with condition and generation info
        library_config: Library configuration
        bc1_col: Name of BC1 column
        sample_col: Name of sample column
        count_col: Name of count column

    Returns:
        DataFrame with bc1 and t0_total columns
    """
    # Identify t0 samples for this library
    # Handle generations as numeric for robust comparison
    generations_numeric = pd.to_numeric(keyfile["generations"], errors="coerce").fillna(-1).astype(int)
    t0_samples = keyfile[
        (keyfile["condition"] == "YPD") &
        (generations_numeric == 0) &
        (keyfile["screen_date"].isin(library_config.screen_dates))
    ][sample_col].tolist()

    if not t0_samples:
        raise ValueError(f"No t0 samples found for library {library_config.name}")

    # Filter to t0 samples
    t0_counts = long_table[long_table[sample_col].isin(t0_samples)]

    # Sum counts per BC1
    bc1_t0_totals = (
        t0_counts
        .groupby(bc1_col)[count_col]
        .sum()
        .reset_index()
        .rename(columns={count_col: "t0_total"})
    )

    return bc1_t0_totals


def classify_bc1_tiers(
    bc1_t0_counts: pd.DataFrame,
    tier_config: TierConfig,
    bc1_col: str = "bc1",
    t0_col: str = "t0_total",
) -> pd.DataFrame:
    """
    Classify BC1s into tiers based on t0 counts.

    Args:
        bc1_t0_counts: DataFrame with bc1 and t0_total columns
        tier_config: Tier threshold configuration
        bc1_col: Name of BC1 column
        t0_col: Name of t0 total column

    Returns:
        DataFrame with bc1, t0_total, and tier columns
    """
    result = bc1_t0_counts.copy()

    # Apply tier classification
    result["tier"] = result[t0_col].apply(tier_config.get_tier)

    return result


def create_tier_summary(
    tier_assignments: pd.DataFrame,
    library_name: str,
    t0_col: str = "t0_total",
    tier_col: str = "tier",
) -> TierSummary:
    """
    Create summary statistics for abundance tier assignments.

    Args:
        tier_assignments: DataFrame from classify_bc1_tiers()
        library_name: Name of the library
        t0_col: Name of t0 total column
        tier_col: Name of tier column

    Returns:
        TierSummary with counts and read totals per abundance tier
    """
    tier_counts = tier_assignments[tier_col].value_counts()
    tier_reads = tier_assignments.groupby(tier_col)[t0_col].sum()

    return TierSummary(
        library=library_name,
        abundance_high_count=tier_counts.get("abundance_high", 0),
        abundance_medium_count=tier_counts.get("abundance_medium", 0),
        abundance_low_count=tier_counts.get("abundance_low", 0),
        abundance_high_reads=tier_reads.get("abundance_high", 0),
        abundance_medium_reads=tier_reads.get("abundance_medium", 0),
        abundance_low_reads=tier_reads.get("abundance_low", 0),
    )


def filter_by_tier(
    long_table: pd.DataFrame,
    tier_assignments: pd.DataFrame,
    target_tiers: List[str],
    bc1_col: str = "bc1",
    tier_col: str = "tier",
) -> pd.DataFrame:
    """
    Filter long table to include only BC1s in specified abundance tiers.

    Args:
        long_table: Long-format BC1 counts table
        tier_assignments: DataFrame with bc1 and tier columns
        target_tiers: List of tiers to include (e.g., ["abundance_high", "abundance_medium"])
        bc1_col: Name of BC1 column
        tier_col: Name of tier column

    Returns:
        Filtered long table
    """
    # Get BC1s in target tiers
    target_bc1s = tier_assignments[
        tier_assignments[tier_col].isin(target_tiers)
    ][bc1_col].tolist()

    # Filter long table
    return long_table[long_table[bc1_col].isin(target_bc1s)]


def aggregate_abundance_medium_by_oligo(
    long_table: pd.DataFrame,
    reference_table: pd.DataFrame,
    bc1_col: str = "bc1",
    sample_col: str = "sample_name",
    count_col: str = "sample_counts",
    oligo_col: str = "oligo_name",
) -> pd.DataFrame:
    """
    Aggregate abundance_medium BC1 counts by oligo to create pseudo-BC1s.

    BC1s in the abundance_medium tier (20-99 t0 reads) are aggregated by their
    assigned oligo before DESeq2 analysis. This increases statistical power
    while preserving oligo-level resolution.

    Args:
        long_table: Long-format BC1 counts for abundance_medium tier only
        reference_table: BC1 to oligo mapping
        bc1_col: Name of BC1 column
        sample_col: Name of sample column
        count_col: Name of count column
        oligo_col: Name of oligo column

    Returns:
        Aggregated long table with oligo as the identifier
    """
    # Merge with reference to get oligo assignments
    merged = long_table.merge(
        reference_table[[bc1_col, oligo_col]],
        on=bc1_col,
        how="left"
    )

    # Aggregate by oligo and sample
    aggregated = (
        merged
        .groupby([oligo_col, sample_col])
        .agg({
            count_col: "sum",
            bc1_col: "count"  # Track number of BC1s aggregated
        })
        .reset_index()
        .rename(columns={bc1_col: "n_bc1s_aggregated"})
    )

    # Rename oligo column to bc1 for compatibility with downstream steps
    aggregated = aggregated.rename(columns={oligo_col: bc1_col})

    return aggregated
