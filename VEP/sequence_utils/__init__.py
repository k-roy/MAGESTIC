"""
Sequence Utilities Module

Common sequence operations used across the VEP pipeline:
- basics: reverse_complement, pad_sequence, gc_content
- translation: DNA to protein translation (CODON_TABLE, translate)
- one_hot: One-hot encoding for neural models
- fasta_io: Read/write FASTA files
- vcf_parsing: Parse VCF files and variants
- variant_ops: Apply variants to sequences
- pairwise_variants: Extract and apply pairwise variants between strains (VCF-guided)
- alignment_variants: Extract variants by direct sequence alignment (handles complex cases)

Author: Kevin R. Roy
"""

# Basic sequence operations
from .basics import (
    reverse_complement,
    pad_sequence,
    gc_content,
)

# Translation
from .translation import (
    CODON_TABLE,
    translate,
    get_codon_at_position,
)

# One-hot encoding
from .one_hot import (
    one_hot_encode,
    one_hot_decode,
)

# FASTA I/O
from .fasta_io import (
    read_fasta,
    write_fasta,
    iter_fasta,
)

# VCF parsing
from .vcf_parsing import (
    Variant,
    parse_vcf_genotype,
    load_variants_for_region,
    get_strain_variants,
)

# Variant operations
from .variant_ops import (
    apply_variants_to_sequence,
    insert_variant,
    filter_overlapping_variants,
    get_offset_at_position,
)

# Pairwise variant operations (added 2026-02-17, updated 2026-02-23)
from .pairwise_variants import (
    PairwiseVariant,
    LinkedVariantGroup,
    compute_strain_seq_position,
    extract_pairwise_variants_vcf_guided,
    extract_pairwise_variants_alignment,
    extract_pairwise_variants_auto,
    apply_pairwise_variant,
    apply_multiple_pairwise_variants,
    detect_linked_variants,
    detect_linked_variants_cds_aware,
    classify_linked_effect,
)

# Alignment-based variant extraction (added 2026-02-22)
# Use these functions when VCF-guided approach fails (e.g., both strains have
# different variants at the same S288C position)
# NOTE: Requires mappy (minimap2) - make import optional
try:
    from .alignment_variants import (
        AlignmentVariant,
        align_and_extract_variants,
        extract_bidirectional_variants,
        apply_alignment_variant,
        load_offset_table,
        s288c_to_strain_pos,
        reverse_offset_lookup,
        build_reverse_mapping,
    )
    _ALIGNMENT_VARIANTS_AVAILABLE = True
except ImportError:
    _ALIGNMENT_VARIANTS_AVAILABLE = False
    # Define placeholders that raise helpful errors
    def _missing_mappy(*args, **kwargs):
        raise ImportError(
            "alignment_variants requires mappy (minimap2). "
            "Install with: pip install mappy"
        )
    AlignmentVariant = None
    align_and_extract_variants = _missing_mappy
    extract_bidirectional_variants = _missing_mappy
    apply_alignment_variant = _missing_mappy
    load_offset_table = _missing_mappy
    s288c_to_strain_pos = _missing_mappy
    reverse_offset_lookup = _missing_mappy
    build_reverse_mapping = _missing_mappy

__all__ = [
    # basics
    'reverse_complement',
    'pad_sequence',
    'gc_content',
    # translation
    'CODON_TABLE',
    'translate',
    'get_codon_at_position',
    # one_hot
    'one_hot_encode',
    'one_hot_decode',
    # fasta_io
    'read_fasta',
    'write_fasta',
    'iter_fasta',
    # vcf_parsing
    'Variant',
    'parse_vcf_genotype',
    'load_variants_for_region',
    'get_strain_variants',
    # variant_ops
    'apply_variants_to_sequence',
    'insert_variant',
    'filter_overlapping_variants',
    'get_offset_at_position',
    # pairwise_variants
    'PairwiseVariant',
    'LinkedVariantGroup',
    'compute_strain_seq_position',
    'extract_pairwise_variants_vcf_guided',
    'extract_pairwise_variants_alignment',
    'extract_pairwise_variants_auto',
    'apply_pairwise_variant',
    'apply_multiple_pairwise_variants',
    'detect_linked_variants',
    'detect_linked_variants_cds_aware',
    'classify_linked_effect',
    # alignment_variants
    'AlignmentVariant',
    'align_and_extract_variants',
    'extract_bidirectional_variants',
    'apply_alignment_variant',
    'load_offset_table',
    's288c_to_strain_pos',
    'reverse_offset_lookup',
    'build_reverse_mapping',
]

__version__ = '0.3.0'  # Bumped for alignment_variants addition
