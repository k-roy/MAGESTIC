"""
Core modules for the BC1-Donor-BC0 pipeline.

Contains:
- merge: Early merge of 2x150 and 2x300 data
- matching: Unified sliding window donor matching
- guide_disambiguation: Multi-PAM oligo disambiguation
- structure_inference: De novo inference of amplicon structure (NEW)
"""

from .merge import (
    BC1Record,
    merge_data_sources,
    resolve_donor_conflicts,
)

from .matching import (
    OligoMatch,
    OligoLookupTables,
    build_lookup_tables,
    match_sequence,
    match_dataframe,
    load_oligo_design,
    load_and_build_lookup,
    # New: Load from harmonized table (preferred)
    load_oligos_from_harmonized_table,
    build_lookup_from_harmonized,
)

from .guide_disambiguation import (
    ReassignmentResult,
    MultiPamOligoLookup,
    hamming_distance,
    load_oligo_designs_for_disambiguation,
    find_best_oligo_by_guide as find_best_guide_match,  # Alias for top-level import
    reassign_oligos_by_guide,
    run_guide_disambiguation,
)

from .structure_inference import (
    InferredConstantRegion,
    BC1StructureInferenceResult,
    MSD_SEQUENCE,
    MSD_SEQUENCE_RC,
    KNOWN_SPS_SEQUENCES,
    LBCAS12A_SCAFFOLD_VARIANTS,
    infer_bc1_amplicon_structure,
    infer_structure_from_multiple_fastqs,
    get_extraction_parameters,
    print_inference_report,
)
