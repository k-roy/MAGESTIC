"""
Quality Control functions for the Guide-Donor-BC0 pipeline.

QC checks include:
- Oligo design file validation (columns, lengths, uniqueness)
- Sequence matching quality (match rate, confidence distribution)
- BC0 purity (fraction of reads from dominant oligo per bc0)
- Library representation (oligo coverage across pool)

Usage:
    from guide_donor_bc0_pipeline.qc import run_pipeline_qc, QCConfig

    qc_config = QCConfig()
    metrics = run_pipeline_qc(
        oligo_counts_df=counts,
        bc0_purity_df=purity,
        oligo_design_df=design,
        config=qc_config
    )

    if metrics.has_warnings:
        print(metrics.get_warning_summary())
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


# Oligo structure constants (200mer design)
EXPECTED_OLIGO_LENGTH = 200
GUIDE_START, GUIDE_END = 20, 40
DONOR_START, DONOR_END = 51, 180
EXPECTED_GUIDE_LENGTH = GUIDE_END - GUIDE_START  # 20bp
EXPECTED_DONOR_LENGTH = DONOR_END - DONOR_START  # 129bp


@dataclass
class QCConfig:
    """Configuration for QC thresholds."""

    # Oligo design validation
    expected_oligo_length: int = EXPECTED_OLIGO_LENGTH
    expected_guide_length: int = EXPECTED_GUIDE_LENGTH
    expected_donor_length: int = EXPECTED_DONOR_LENGTH

    # Match rate thresholds
    min_match_rate: float = 0.80          # At least 80% should match
    match_rate_warning: float = 0.90      # Warn if below 90%

    # BC0 purity thresholds
    min_bc0_purity: float = 0.70          # Flag if purity < 0.7
    purity_warning: float = 0.85          # Warn if purity < 0.85
    min_reads_per_bc0: int = 10           # Minimum reads to consider purity valid

    # Library representation
    min_coverage: float = 0.80            # At least 80% of oligos covered
    coverage_warning: float = 0.95        # Warn if below 95%
    min_reads_for_present: int = 1        # Minimum reads to count as "present"


@dataclass
class QCMetrics:
    """Container for QC metrics from pipeline validation."""

    # Oligo design metrics
    n_designed_oligos: int = 0
    oligo_length_issues: int = 0
    guide_length_issues: int = 0
    donor_length_issues: int = 0
    duplicate_oligo_names: int = 0
    invalid_dna_sequences: int = 0

    # Matching quality metrics
    total_sequences: int = 0
    total_matched: int = 0
    total_unmatched: int = 0
    match_rate: float = 0.0
    match_type_distribution: Dict[str, int] = field(default_factory=dict)
    confidence_distribution: Dict[str, int] = field(default_factory=dict)
    mean_confidence: float = 0.0

    # BC0 purity metrics
    n_unique_bc0s: int = 0
    mean_bc0_purity: float = 0.0
    median_bc0_purity: float = 0.0
    low_purity_bc0s: int = 0
    purity_distribution: Dict[str, int] = field(default_factory=dict)

    # Library representation
    n_oligos_covered: int = 0
    n_oligos_missing: int = 0
    coverage_rate: float = 0.0
    coverage_depth_distribution: Dict[str, int] = field(default_factory=dict)
    min_reads_per_oligo: int = 0
    max_reads_per_oligo: int = 0
    median_reads_per_oligo: float = 0.0

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
            "Guide-Donor-BC0 Pipeline QC Summary",
            "=" * 60,
            "",
            "Oligo Design:",
            f"  Total designed oligos: {self.n_designed_oligos:,}",
            f"  Length issues: {self.oligo_length_issues}",
            f"  Duplicate names: {self.duplicate_oligo_names}",
            "",
            "Matching Quality:",
            f"  Match rate: {self.match_rate:.1%}",
            f"  Matched sequences: {self.total_matched:,}",
            f"  Unmatched sequences: {self.total_unmatched:,}",
            f"  Mean confidence: {self.mean_confidence:.3f}",
            "",
            "Match Type Distribution:",
        ]
        for match_type, count in sorted(self.match_type_distribution.items()):
            pct = 100 * count / max(self.total_matched, 1)
            lines.append(f"  {match_type}: {count:,} ({pct:.1f}%)")

        lines.extend([
            "",
            "BC0 Purity:",
            f"  Unique BC0s: {self.n_unique_bc0s:,}",
            f"  Mean purity: {self.mean_bc0_purity:.3f}",
            f"  Median purity: {self.median_bc0_purity:.3f}",
            f"  Low purity BC0s: {self.low_purity_bc0s:,}",
            "",
            "Library Representation:",
            f"  Coverage: {self.coverage_rate:.1%}",
            f"  Oligos covered: {self.n_oligos_covered:,}",
            f"  Oligos missing: {self.n_oligos_missing:,}",
            f"  Median reads/oligo: {self.median_reads_per_oligo:.0f}",
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


def validate_oligo_design(
    oligo_design_df: pd.DataFrame,
    config: QCConfig,
) -> QCMetrics:
    """
    Validate the oligo design file.

    Args:
        oligo_design_df: DataFrame with oligo_name and oligo_seq columns
        config: QC configuration

    Returns:
        QCMetrics with design validation results
    """
    metrics = QCMetrics()

    metrics.n_designed_oligos = len(oligo_design_df)

    # Check for required columns
    if "oligo_name" not in oligo_design_df.columns:
        metrics.errors.append("Missing required column: oligo_name")
        return metrics

    if "oligo_seq" not in oligo_design_df.columns:
        metrics.errors.append("Missing required column: oligo_seq")
        return metrics

    # Check oligo lengths
    oligo_lens = oligo_design_df["oligo_seq"].str.len()
    metrics.oligo_length_issues = (oligo_lens != config.expected_oligo_length).sum()

    if metrics.oligo_length_issues > 0:
        metrics.warnings.append(
            f"{metrics.oligo_length_issues} oligos have unexpected length "
            f"(expected {config.expected_oligo_length}bp)"
        )

    # Check for duplicate names
    name_counts = oligo_design_df["oligo_name"].value_counts()
    metrics.duplicate_oligo_names = (name_counts > 1).sum()

    if metrics.duplicate_oligo_names > 0:
        metrics.errors.append(
            f"{metrics.duplicate_oligo_names} duplicate oligo names found"
        )

    # Validate DNA sequences (only ACGT)
    valid_dna = oligo_design_df["oligo_seq"].str.match(r"^[ACGTacgt]+$", na=False)
    metrics.invalid_dna_sequences = (~valid_dna).sum()

    if metrics.invalid_dna_sequences > 0:
        metrics.warnings.append(
            f"{metrics.invalid_dna_sequences} oligos contain non-ACGT characters"
        )

    # Extract and validate guide lengths
    if config.expected_oligo_length == EXPECTED_OLIGO_LENGTH:
        guides = oligo_design_df["oligo_seq"].str[GUIDE_START:GUIDE_END]
        guide_lens = guides.str.len()
        metrics.guide_length_issues = (guide_lens != config.expected_guide_length).sum()

        # Extract and validate donor lengths
        donors = oligo_design_df["oligo_seq"].str[DONOR_START:DONOR_END]
        donor_lens = donors.str.len()
        metrics.donor_length_issues = (donor_lens != config.expected_donor_length).sum()

    return metrics


def validate_matching_quality(
    oligo_counts_df: pd.DataFrame,
    config: QCConfig,
) -> QCMetrics:
    """
    Validate sequence matching quality.

    Args:
        oligo_counts_df: DataFrame with matching results
        config: QC configuration

    Returns:
        QCMetrics with matching validation results
    """
    metrics = QCMetrics()

    metrics.total_sequences = len(oligo_counts_df)

    # Check match rate
    if "oligo_name" in oligo_counts_df.columns:
        matched_mask = oligo_counts_df["oligo_name"].notna() & (oligo_counts_df["oligo_name"] != "")
        metrics.total_matched = matched_mask.sum()
        metrics.total_unmatched = (~matched_mask).sum()
        metrics.match_rate = metrics.total_matched / max(metrics.total_sequences, 1)

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

    # Match type distribution
    if "match_type" in oligo_counts_df.columns:
        type_counts = oligo_counts_df["match_type"].value_counts()
        metrics.match_type_distribution = type_counts.to_dict()

    # Confidence distribution
    if "confidence" in oligo_counts_df.columns:
        conf = oligo_counts_df["confidence"].dropna()
        metrics.mean_confidence = conf.mean() if len(conf) > 0 else 0.0

        # Bin confidence scores
        bins = [0, 0.5, 0.7, 0.85, 0.95, 1.0]
        labels = ["<0.5", "0.5-0.7", "0.7-0.85", "0.85-0.95", "0.95-1.0"]
        conf_binned = pd.cut(conf, bins=bins, labels=labels, include_lowest=True)
        metrics.confidence_distribution = conf_binned.value_counts().to_dict()

    return metrics


def validate_bc0_purity(
    bc0_purity_df: pd.DataFrame,
    config: QCConfig,
) -> QCMetrics:
    """
    Validate BC0 purity metrics.

    Args:
        bc0_purity_df: DataFrame with bc0 purity calculations
        config: QC configuration

    Returns:
        QCMetrics with purity validation results
    """
    metrics = QCMetrics()

    if bc0_purity_df is None or len(bc0_purity_df) == 0:
        metrics.warnings.append("No BC0 purity data available")
        return metrics

    metrics.n_unique_bc0s = len(bc0_purity_df)

    # Purity statistics
    if "purity" in bc0_purity_df.columns:
        purity = bc0_purity_df["purity"].dropna()
        metrics.mean_bc0_purity = purity.mean()
        metrics.median_bc0_purity = purity.median()
        metrics.low_purity_bc0s = (purity < config.min_bc0_purity).sum()

        # Purity distribution
        bins = [0, 0.5, 0.7, 0.85, 0.95, 1.0]
        labels = ["<0.5", "0.5-0.7", "0.7-0.85", "0.85-0.95", "0.95-1.0"]
        purity_binned = pd.cut(purity, bins=bins, labels=labels, include_lowest=True)
        metrics.purity_distribution = purity_binned.value_counts().to_dict()

        if metrics.mean_bc0_purity < config.min_bc0_purity:
            metrics.warnings.append(
                f"Mean BC0 purity {metrics.mean_bc0_purity:.3f} is below threshold "
                f"({config.min_bc0_purity})"
            )

    return metrics


def validate_library_representation(
    representation_df: pd.DataFrame,
    n_designed_oligos: int,
    config: QCConfig,
) -> QCMetrics:
    """
    Validate library representation/coverage.

    Args:
        representation_df: DataFrame with per-oligo read counts
        n_designed_oligos: Total number of designed oligos
        config: QC configuration

    Returns:
        QCMetrics with representation validation results
    """
    metrics = QCMetrics()

    if representation_df is None or len(representation_df) == 0:
        metrics.warnings.append("No library representation data available")
        return metrics

    # Count covered oligos
    if "reads" in representation_df.columns:
        covered = representation_df["reads"] >= config.min_reads_for_present
        metrics.n_oligos_covered = covered.sum()
        metrics.n_oligos_missing = n_designed_oligos - metrics.n_oligos_covered
        metrics.coverage_rate = metrics.n_oligos_covered / max(n_designed_oligos, 1)

        if metrics.coverage_rate < config.min_coverage:
            metrics.errors.append(
                f"Library coverage {metrics.coverage_rate:.1%} is below minimum "
                f"({config.min_coverage:.1%})"
            )
        elif metrics.coverage_rate < config.coverage_warning:
            metrics.warnings.append(
                f"Library coverage {metrics.coverage_rate:.1%} is below warning "
                f"({config.coverage_warning:.1%})"
            )

        # Read depth statistics
        reads = representation_df["reads"]
        metrics.min_reads_per_oligo = reads.min()
        metrics.max_reads_per_oligo = reads.max()
        metrics.median_reads_per_oligo = reads.median()

        # Coverage depth distribution
        bins = [0, 1, 10, 100, 1000, float('inf')]
        labels = ["0", "1-9", "10-99", "100-999", "1000+"]
        depth_binned = pd.cut(reads, bins=bins, labels=labels, include_lowest=True)
        metrics.coverage_depth_distribution = depth_binned.value_counts().to_dict()

    return metrics


def run_pipeline_qc(
    oligo_counts_df: Optional[pd.DataFrame] = None,
    bc0_purity_df: Optional[pd.DataFrame] = None,
    oligo_design_df: Optional[pd.DataFrame] = None,
    representation_df: Optional[pd.DataFrame] = None,
    config: Optional[QCConfig] = None,
) -> QCMetrics:
    """
    Run comprehensive QC on pipeline outputs.

    This is the main entry point for QC validation after the pipeline
    has completed.

    Args:
        oligo_counts_df: DataFrame with oligo match counts
        bc0_purity_df: DataFrame with BC0 purity calculations
        oligo_design_df: DataFrame with oligo design
        representation_df: DataFrame with per-oligo coverage
        config: QC configuration (uses defaults if not provided)

    Returns:
        QCMetrics with all validation results
    """
    if config is None:
        config = QCConfig()

    logger.info("Running Guide-Donor-BC0 pipeline QC...")

    # Combined metrics
    combined = QCMetrics()

    # Validate oligo design
    if oligo_design_df is not None:
        design_metrics = validate_oligo_design(oligo_design_df, config)
        combined.n_designed_oligos = design_metrics.n_designed_oligos
        combined.oligo_length_issues = design_metrics.oligo_length_issues
        combined.duplicate_oligo_names = design_metrics.duplicate_oligo_names
        combined.warnings.extend(design_metrics.warnings)
        combined.errors.extend(design_metrics.errors)

    # Validate matching quality
    if oligo_counts_df is not None:
        match_metrics = validate_matching_quality(oligo_counts_df, config)
        combined.total_sequences = match_metrics.total_sequences
        combined.total_matched = match_metrics.total_matched
        combined.total_unmatched = match_metrics.total_unmatched
        combined.match_rate = match_metrics.match_rate
        combined.match_type_distribution = match_metrics.match_type_distribution
        combined.confidence_distribution = match_metrics.confidence_distribution
        combined.mean_confidence = match_metrics.mean_confidence
        combined.warnings.extend(match_metrics.warnings)
        combined.errors.extend(match_metrics.errors)

    # Validate BC0 purity
    if bc0_purity_df is not None:
        purity_metrics = validate_bc0_purity(bc0_purity_df, config)
        combined.n_unique_bc0s = purity_metrics.n_unique_bc0s
        combined.mean_bc0_purity = purity_metrics.mean_bc0_purity
        combined.median_bc0_purity = purity_metrics.median_bc0_purity
        combined.low_purity_bc0s = purity_metrics.low_purity_bc0s
        combined.purity_distribution = purity_metrics.purity_distribution
        combined.warnings.extend(purity_metrics.warnings)
        combined.errors.extend(purity_metrics.errors)

    # Validate library representation
    if representation_df is not None:
        rep_metrics = validate_library_representation(
            representation_df, combined.n_designed_oligos, config
        )
        combined.n_oligos_covered = rep_metrics.n_oligos_covered
        combined.n_oligos_missing = rep_metrics.n_oligos_missing
        combined.coverage_rate = rep_metrics.coverage_rate
        combined.coverage_depth_distribution = rep_metrics.coverage_depth_distribution
        combined.min_reads_per_oligo = rep_metrics.min_reads_per_oligo
        combined.max_reads_per_oligo = rep_metrics.max_reads_per_oligo
        combined.median_reads_per_oligo = rep_metrics.median_reads_per_oligo
        combined.warnings.extend(rep_metrics.warnings)
        combined.errors.extend(rep_metrics.errors)

    # Log summary
    logger.info(f"QC complete: {combined.n_designed_oligos:,} oligos, "
                f"match rate {combined.match_rate:.1%}, "
                f"coverage {combined.coverage_rate:.1%}")

    if combined.has_errors:
        logger.error(f"QC FAILED with {len(combined.errors)} error(s)")
    elif combined.has_warnings:
        logger.warning(f"QC passed with {len(combined.warnings)} warning(s)")
    else:
        logger.info("QC PASSED")

    return combined


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
