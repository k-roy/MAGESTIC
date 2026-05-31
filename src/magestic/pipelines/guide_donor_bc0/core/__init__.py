"""
Core processing functions for the Guide-Donor-BC0 pipeline.
"""

from .oligo_matching import (
    load_oligo_pool,
    build_lookup_dictionaries,
    resolve_oligo_name,
    process_counts_file,
    OligoMatchResult,
)

from .structure_inference import (
    NucleaseAnchors,
    NUCLEASE_ANCHORS,
    SPS_TO_LIBRARY,
    StructureInferenceResult,
    infer_amplicon_structure,
    infer_structure_from_multiple_fastqs,
    get_structure_for_nuclease,
    get_nuclease_from_sps,
    rev_comp,
    detect_nuclease_from_oligo,
)

from .fastq_processing import (
    SequenceConstants,
    SampleInfo,
    load_sample_keyfile,
    parse_sequence,
    process_collapsed_fastq,
    calculate_totals,
    write_counts_file,
    process_library,
)

from .bc0_purity import (
    PurityConfig,
    find_input_files,
    calculate_bc0_purity,
    INPUT_SUFFIX as PURITY_INPUT_SUFFIX,
    OUTPUT_SUFFIX as PURITY_OUTPUT_SUFFIX,
)

from .library_representation import (
    OligoConstants,
    gc_content,
    gini_index,
    load_purity_files,
    load_designed_oligos,
    aggregate_by_oligo_name,
    add_missing_oligos,
    filter_by_top_sps,
    save_representation_outputs,
    run_representation_analysis,
    plot_gc_vs_counts,
    plot_density_by_sps,
    plot_gini_indices,
    plot_category_percentages,
    run_library_plots,
)

__all__ = [
    # Oligo matching
    "load_oligo_pool",
    "build_lookup_dictionaries",
    "resolve_oligo_name",
    "process_counts_file",
    "OligoMatchResult",
    # Structure inference
    "NucleaseAnchors",
    "NUCLEASE_ANCHORS",
    "SPS_TO_LIBRARY",
    "StructureInferenceResult",
    "infer_amplicon_structure",
    "infer_structure_from_multiple_fastqs",
    "get_structure_for_nuclease",
    "get_nuclease_from_sps",
    "rev_comp",
    "detect_nuclease_from_oligo",
    # FASTQ processing (Step 02)
    "SequenceConstants",
    "SampleInfo",
    "load_sample_keyfile",
    "parse_sequence",
    "process_collapsed_fastq",
    "calculate_totals",
    "write_counts_file",
    "process_library",
    # BC0 purity (Step 04)
    "PurityConfig",
    "find_input_files",
    "calculate_bc0_purity",
    "PURITY_INPUT_SUFFIX",
    "PURITY_OUTPUT_SUFFIX",
    # Library representation (Steps 05–06)
    "OligoConstants",
    "gc_content",
    "gini_index",
    "load_purity_files",
    "load_designed_oligos",
    "aggregate_by_oligo_name",
    "add_missing_oligos",
    "filter_by_top_sps",
    "save_representation_outputs",
    "run_representation_analysis",
    "plot_gc_vs_counts",
    "plot_density_by_sps",
    "plot_gini_indices",
    "plot_category_percentages",
    "run_library_plots",
]
