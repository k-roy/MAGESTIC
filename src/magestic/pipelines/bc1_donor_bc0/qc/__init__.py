"""
Quality Control module for the BC1-Donor-BC0 pipeline.

This module provides QC functions for validating:
- Input data integrity (files, columns, sequences)
- Merge quality (conflicts, coverage)
- Oligo matching (match rates, confidence)
- BC1 purity and consistency
"""

from .pipeline_qc import (
    run_pipeline_qc,
    validate_input_files,
    validate_merge_quality,
    validate_matching_quality,
    generate_qc_report,
    QCMetrics,
    QCConfig,
)

__all__ = [
    "run_pipeline_qc",
    "validate_input_files",
    "validate_merge_quality",
    "validate_matching_quality",
    "generate_qc_report",
    "QCMetrics",
    "QCConfig",
]
