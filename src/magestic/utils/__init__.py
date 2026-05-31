"""
magestic.utils — Shared utilities for MAGESTIC design and analysis pipelines.

Migrated from common/scripts/utils/ (16 modules).

Modules:
- sequence_utils: DNA sequence operations (reverse_complement, parse_cigar, etc.)
- alignment_utils: SAM/BAM parsing and alignment comparison
- oligo_annotations: Pre-parsed QTL oligo annotations from harmonized table
- gff_utils: GFF3 parsing with GeneFeature dataclass
- vcf_parsing: VCF file reading and FDR filtering
- variant_annotation: Variant consequence annotation
- strain_utils: Yeast strain constants and utilities
- coordinate_utils: Genomic coordinate conversion utilities
- fasta_io: FASTA file I/O
- sequence_anchors: PCR primer / sequencing anchor utilities
- translation: Codon table loading and translation
- collapse_utils: Read collapsing and deduplication
- config_utils: YAML config loading utilities
- path_utils: Path resolution helpers
- thread_limiting: SLURM thread-count environment setup
"""

from .sequence_utils import (
    reverse_complement,
    rev_comp,
    parse_cigar,
    get_alignment_length,
    hamming_distance,
    gc_content,
    find_kmer_positions,
    translate_codon,
)

from .alignment_utils import (
    AlignmentScore,
    AlignmentInfo,
    extract_alignment_info,
    parse_sam_file,
    select_best_alignment,
    select_best_alignment_paired,
)

from .oligo_annotations import (
    load_oligo_annotations,
    annotate_dataframe,
    get_oligo_annotation,
    get_annotation_columns_for_bc1_from_screen,
    get_annotation_columns_for_bc1_donor_bc0,
    HARMONIZED_OLIGO_TABLE,
)

__all__ = [
    # sequence_utils
    'reverse_complement',
    'rev_comp',
    'parse_cigar',
    'get_alignment_length',
    'hamming_distance',
    'gc_content',
    'find_kmer_positions',
    'translate_codon',
    # alignment_utils
    'AlignmentScore',
    'AlignmentInfo',
    'extract_alignment_info',
    'parse_sam_file',
    'select_best_alignment',
    'select_best_alignment_paired',
    # oligo_annotations
    'load_oligo_annotations',
    'annotate_dataframe',
    'get_oligo_annotation',
    'get_annotation_columns_for_bc1_from_screen',
    'get_annotation_columns_for_bc1_donor_bc0',
    'HARMONIZED_OLIGO_TABLE',
]
