"""
Core processing functions for the BC1-From-Screen pipeline.
"""

from .tier_classification import (
    classify_bc1_tiers,
    calculate_t0_counts,
    create_tier_summary,
)

from .count_matrix import (
    create_wide_count_matrix,
    consolidate_t0_samples,
    split_by_comparison,
)

from .hit_calling import (
    call_hits_single_comparison,
    aggregate_hit_calls,
    HitCallResult,
)

# Empirical hit calling framework (v2)
from .library_calibration import (
    calculate_library_bias,
    center_by_library,
    run_bias_correction_pipeline,
    LibraryBiasStats,
)

from .empirical_thresholds import (
    calibrate_thresholds_per_condition,
    validate_ko_enrichment,
    validate_missense_enrichment,
    run_threshold_calibration_pipeline,
    ThresholdCalibration,
)

from .noise_assessment import (
    identify_noisy_conditions,
    identify_noisy_loci,
    assess_variant_noise,
    apply_noise_penalties,
    run_noise_assessment_pipeline,
    NoiseFlags,
    NoisePenaltyThresholds,
)

from .causality import (
    CausalityEvidence,
    CausalityEvaluator,
    aggregate_to_aa_level,
    generate_causal_variants_output,
    generate_control_qc_output,
)

__all__ = [
    # Tier classification
    "classify_bc1_tiers",
    "calculate_t0_counts",
    "create_tier_summary",
    # Count matrix
    "create_wide_count_matrix",
    "consolidate_t0_samples",
    "split_by_comparison",
    # Hit calling (original)
    "call_hits_single_comparison",
    "aggregate_hit_calls",
    "HitCallResult",
    # Library calibration (v2)
    "calculate_library_bias",
    "center_by_library",
    "run_bias_correction_pipeline",
    "LibraryBiasStats",
    # Empirical thresholds (v2)
    "calibrate_thresholds_per_condition",
    "validate_ko_enrichment",
    "validate_missense_enrichment",
    "run_threshold_calibration_pipeline",
    "ThresholdCalibration",
    # Noise assessment (v2)
    "identify_noisy_conditions",
    "identify_noisy_loci",
    "assess_variant_noise",
    "apply_noise_penalties",
    "run_noise_assessment_pipeline",
    "NoiseFlags",
    "NoisePenaltyThresholds",
    # Causality framework (v2)
    "CausalityEvidence",
    "CausalityEvaluator",
    "aggregate_to_aa_level",
    "generate_causal_variants_output",
    "generate_control_qc_output",
]
