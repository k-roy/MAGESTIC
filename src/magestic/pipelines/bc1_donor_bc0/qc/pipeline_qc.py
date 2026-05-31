"""
Quality Control functions for the BC1-Donor-BC0 pipeline.

QC checks include:
- Input file validation (presence, format, columns)
- Data source consistency (bc1 length, bc0 length, donor patterns)
- Merge quality (conflict rate, coverage across sources)
- Oligo matching (match rate by method, confidence distribution)
- BC1 purity (fraction of bc1 counts from dominant donor)

Usage:
    from bc1_donor_bc0_pipeline.qc import run_pipeline_qc, QCConfig

    qc_config = QCConfig()
    metrics = run_pipeline_qc(reference_table, qc_config)

    if metrics.has_warnings:
        print(metrics.get_warning_summary())
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class QCConfig:
    """Configuration for QC thresholds."""

    # Sequence validation
    expected_bc1_length: int = 20
    expected_bc0_length: int = 10
    expected_donor_length_range: Tuple[int, int] = (100, 200)

    # Match rate thresholds
    min_match_rate: float = 0.80          # At least 80% should match
    match_rate_warning: float = 0.90      # Warn if below 90%

    # Purity thresholds
    min_bc1_purity: float = 0.70          # Flag if bc1_purity < 0.7
    purity_warning: float = 0.85          # Warn if bc1_purity < 0.85

    # Conflict thresholds
    max_conflict_rate: float = 0.10       # Fail if >10% have conflicts
    conflict_warning: float = 0.05        # Warn if >5% have conflicts

    # Coverage thresholds
    min_pcr_replicates: int = 2           # Minimum PCR replicates expected
    min_total_reads: int = 10             # Minimum reads per bc1


@dataclass
class QCMetrics:
    """Container for QC metrics from pipeline validation."""

    # Input metrics
    n_input_rows: int = 0
    n_unique_bc1s: int = 0
    n_unique_bc0s: int = 0

    # Sequence validation
    bc1_length_issues: int = 0
    bc0_length_issues: int = 0
    donor_length_issues: int = 0

    # Data source coverage
    source_coverage: Dict[str, int] = field(default_factory=dict)
    has_2x150: bool = False
    has_2x300: bool = False
    has_gdb: bool = False
    multi_source_bc1s: int = 0

    # Merge quality
    conflict_count: int = 0
    conflict_rate: float = 0.0
    conflict_details: List[str] = field(default_factory=list)

    # Matching quality
    total_matched: int = 0
    total_unmatched: int = 0
    match_rate: float = 0.0
    match_method_distribution: Dict[str, int] = field(default_factory=dict)
    unambiguous_count: int = 0
    ambiguous_count: int = 0

    # Purity metrics
    mean_bc1_purity: float = 0.0
    median_bc1_purity: float = 0.0
    low_purity_count: int = 0
    purity_distribution: Dict[str, int] = field(default_factory=dict)

    # Read depth
    mean_reads_per_bc1: float = 0.0
    median_reads_per_bc1: float = 0.0
    low_read_count: int = 0

    # Warnings and flags
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    @property
    def has_warnings(self) -> bool:
        return len(self.warnings) > 0

    @property
    def has_errors(self) -> bool:
        return len(self.errors) > 0

    @property
    def passed(self) -> bool:
        return not self.has_errors

    def get_summary(self) -> str:
        """Return a human-readable summary of QC metrics."""
        lines = [
            "=" * 60,
            "BC1-Donor-BC0 Pipeline QC Summary",
            "=" * 60,
            "",
            f"Total BC1s: {self.n_unique_bc1s:,}",
            f"Total BC0s: {self.n_unique_bc0s:,}",
            "",
            "Data Sources:",
            f"  2x150 present: {self.has_2x150}",
            f"  2x300 present: {self.has_2x300}",
            f"  GDB present: {self.has_gdb}",
            f"  BC1s with multiple sources: {self.multi_source_bc1s:,}",
            "",
            "Matching Quality:",
            f"  Match rate: {self.match_rate:.1%}",
            f"  Matched: {self.total_matched:,}",
            f"  Unmatched: {self.total_unmatched:,}",
            f"  Unambiguous matches: {self.unambiguous_count:,}",
            "",
            "Match Method Distribution:",
        ]
        for method, count in sorted(self.match_method_distribution.items()):
            pct = 100 * count / max(self.total_matched, 1)
            lines.append(f"  {method}: {count:,} ({pct:.1f}%)")

        lines.extend([
            "",
            "BC1 Purity:",
            f"  Mean purity: {self.mean_bc1_purity:.3f}",
            f"  Median purity: {self.median_bc1_purity:.3f}",
            f"  Low purity BC1s: {self.low_purity_count:,}",
            "",
            "Conflicts:",
            f"  Conflict rate: {self.conflict_rate:.2%}",
            f"  BC1s with conflicts: {self.conflict_count:,}",
        ])

        if self.warnings:
            lines.extend(["", "WARNINGS:"])
            for w in self.warnings:
                lines.append(f"  - {w}")

        if self.errors:
            lines.extend(["", "ERRORS:"])
            for e in self.errors:
                lines.append(f"  - {e}")

        lines.append("=" * 60)
        return "\n".join(lines)


def validate_input_files(
    data_sources: Dict[str, Tuple[pd.DataFrame, 'DataSourceInfo']],
    config: QCConfig,
) -> QCMetrics:
    """
    Validate input data files before processing.

    Args:
        data_sources: Dict of data source name -> (DataFrame, DataSourceInfo)
        config: QC configuration

    Returns:
        QCMetrics with input validation results
    """
    metrics = QCMetrics()

    for source_name, (df, info) in data_sources.items():
        if df is None:
            continue

        metrics.source_coverage[source_name] = len(df)

        if source_name == "2x150":
            metrics.has_2x150 = True
        elif source_name == "2x300":
            metrics.has_2x300 = True
        elif source_name == "gdb":
            metrics.has_gdb = True

        # Check bc1 lengths
        if "bc1" in df.columns:
            bad_bc1 = df["bc1"].str.len() != config.expected_bc1_length
            metrics.bc1_length_issues += bad_bc1.sum()

        # Check bc0 lengths
        if "bc0" in df.columns:
            bad_bc0 = df["bc0"].str.len() != config.expected_bc0_length
            metrics.bc0_length_issues += bad_bc0.sum()

        # Check donor lengths
        if "donor" in df.columns:
            donor_lens = df["donor"].str.len()
            bad_donor = (
                (donor_lens < config.expected_donor_length_range[0]) |
                (donor_lens > config.expected_donor_length_range[1])
            )
            metrics.donor_length_issues += bad_donor.sum()

    # Add warnings
    if metrics.bc1_length_issues > 0:
        metrics.warnings.append(
            f"{metrics.bc1_length_issues} BC1s have unexpected length "
            f"(expected {config.expected_bc1_length}bp)"
        )

    if metrics.bc0_length_issues > 0:
        metrics.warnings.append(
            f"{metrics.bc0_length_issues} BC0s have unexpected length "
            f"(expected {config.expected_bc0_length}bp)"
        )

    return metrics


def validate_merge_quality(
    merged_df: pd.DataFrame,
    config: QCConfig,
) -> QCMetrics:
    """
    Validate merge quality after combining data sources.

    Args:
        merged_df: DataFrame after merging data sources
        config: QC configuration

    Returns:
        QCMetrics with merge validation results
    """
    metrics = QCMetrics()

    metrics.n_input_rows = len(merged_df)
    metrics.n_unique_bc1s = merged_df["bc1"].nunique() if "bc1" in merged_df.columns else 0
    metrics.n_unique_bc0s = merged_df["bc0"].nunique() if "bc0" in merged_df.columns else 0

    # Check for donor conflicts
    if "has_donor_conflict" in merged_df.columns:
        metrics.conflict_count = merged_df["has_donor_conflict"].sum()
        metrics.conflict_rate = metrics.conflict_count / max(len(merged_df), 1)

        if metrics.conflict_rate > config.max_conflict_rate:
            metrics.errors.append(
                f"Conflict rate {metrics.conflict_rate:.1%} exceeds maximum "
                f"({config.max_conflict_rate:.1%})"
            )
        elif metrics.conflict_rate > config.conflict_warning:
            metrics.warnings.append(
                f"Conflict rate {metrics.conflict_rate:.1%} is above warning threshold "
                f"({config.conflict_warning:.1%})"
            )

    # Check data source coverage
    source_cols = ["has_2x150", "has_2x300", "has_gdb"]
    for col in source_cols:
        if col in merged_df.columns:
            metrics.source_coverage[col] = merged_df[col].sum()

    # Count BC1s with multiple sources
    if all(col in merged_df.columns for col in source_cols):
        source_count = (
            merged_df["has_2x150"].astype(int) +
            merged_df["has_2x300"].astype(int) +
            merged_df["has_gdb"].astype(int)
        )
        metrics.multi_source_bc1s = (source_count > 1).sum()

    return metrics


def validate_matching_quality(
    matched_df: pd.DataFrame,
    config: QCConfig,
) -> QCMetrics:
    """
    Validate oligo matching quality.

    Args:
        matched_df: DataFrame after oligo matching
        config: QC configuration

    Returns:
        QCMetrics with matching validation results
    """
    metrics = QCMetrics()

    metrics.n_unique_bc1s = matched_df["bc1"].nunique() if "bc1" in matched_df.columns else 0

    # Match rate
    if "oligo_name" in matched_df.columns:
        matched_mask = matched_df["oligo_name"].notna() & (matched_df["oligo_name"] != "")
        metrics.total_matched = matched_mask.sum()
        metrics.total_unmatched = (~matched_mask).sum()
        metrics.match_rate = metrics.total_matched / max(len(matched_df), 1)

        if metrics.match_rate < config.min_match_rate:
            metrics.errors.append(
                f"Match rate {metrics.match_rate:.1%} is below minimum "
                f"({config.min_match_rate:.1%})"
            )
        elif metrics.match_rate < config.match_rate_warning:
            metrics.warnings.append(
                f"Match rate {metrics.match_rate:.1%} is below warning threshold "
                f"({config.match_rate_warning:.1%})"
            )

    # Match method distribution
    if "match_method" in matched_df.columns:
        method_counts = matched_df["match_method"].value_counts()
        metrics.match_method_distribution = method_counts.to_dict()

    # Unambiguous vs ambiguous
    if "is_unambiguous" in matched_df.columns:
        metrics.unambiguous_count = matched_df["is_unambiguous"].sum()
        metrics.ambiguous_count = len(matched_df) - metrics.unambiguous_count

    # BC1 purity
    if "bc1_purity" in matched_df.columns:
        purity = matched_df["bc1_purity"].dropna()
        metrics.mean_bc1_purity = purity.mean()
        metrics.median_bc1_purity = purity.median()
        metrics.low_purity_count = (purity < config.min_bc1_purity).sum()

        # Purity distribution bins
        bins = [0, 0.5, 0.7, 0.85, 0.95, 1.0]
        labels = ["<0.5", "0.5-0.7", "0.7-0.85", "0.85-0.95", "0.95-1.0"]
        purity_binned = pd.cut(purity, bins=bins, labels=labels, include_lowest=True)
        metrics.purity_distribution = purity_binned.value_counts().to_dict()

        if metrics.mean_bc1_purity < config.min_bc1_purity:
            metrics.warnings.append(
                f"Mean BC1 purity {metrics.mean_bc1_purity:.3f} is below threshold "
                f"({config.min_bc1_purity})"
            )

    # Read depth
    if "counts" in matched_df.columns:
        counts = matched_df["counts"]
        metrics.mean_reads_per_bc1 = counts.mean()
        metrics.median_reads_per_bc1 = counts.median()
        metrics.low_read_count = (counts < config.min_total_reads).sum()

    return metrics


def run_pipeline_qc(
    reference_table: pd.DataFrame,
    config: Optional[QCConfig] = None,
) -> QCMetrics:
    """
    Run comprehensive QC on the final reference table.

    This is the main entry point for QC validation after the pipeline
    has completed. It validates the final output.

    Args:
        reference_table: Final BC1 reference table from pipeline
        config: QC configuration (uses defaults if not provided)

    Returns:
        QCMetrics with all validation results
    """
    if config is None:
        config = QCConfig()

    logger.info("Running BC1-Donor-BC0 pipeline QC...")

    # Validate matching quality (this covers most metrics)
    metrics = validate_matching_quality(reference_table, config)

    # Add merge quality metrics if columns exist
    if "has_donor_conflict" in reference_table.columns:
        merge_metrics = validate_merge_quality(reference_table, config)
        metrics.conflict_count = merge_metrics.conflict_count
        metrics.conflict_rate = merge_metrics.conflict_rate
        metrics.multi_source_bc1s = merge_metrics.multi_source_bc1s

        # Merge warnings/errors
        metrics.warnings.extend(merge_metrics.warnings)
        metrics.errors.extend(merge_metrics.errors)

    # Check data source coverage
    for col in ["has_2x150", "has_2x300", "has_gdb"]:
        if col in reference_table.columns:
            setattr(metrics, col.replace("has_", "has_"), reference_table[col].any())

    logger.info(f"QC complete: {metrics.n_unique_bc1s:,} BC1s, "
                f"match rate {metrics.match_rate:.1%}")

    if metrics.has_errors:
        logger.error(f"QC FAILED with {len(metrics.errors)} error(s)")
    elif metrics.has_warnings:
        logger.warning(f"QC passed with {len(metrics.warnings)} warning(s)")
    else:
        logger.info("QC PASSED")

    return metrics


def generate_qc_report(
    metrics: QCMetrics,
    output_path: Path,
    include_timestamp: bool = True,
) -> Path:
    """
    Generate a QC report file.

    Args:
        metrics: QC metrics to report
        output_path: Path for the report file
        include_timestamp: Whether to include timestamp in report

    Returns:
        Path to the generated report
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        if include_timestamp:
            f.write(f"Generated: {datetime.now().isoformat()}\n\n")

        f.write(metrics.get_summary())

    logger.info(f"QC report saved to: {output_path}")
    return output_path
