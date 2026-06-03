"""
Variant operation utilities.

Functions:
- apply_variants_to_sequence: Apply a list of variants to a reference sequence
- insert_variant: Insert a single variant into a sequence
- filter_overlapping_variants: Remove overlapping variants

Author: Kevin R. Roy
"""

from typing import Dict, List, Optional, Tuple

# Import Variant from vcf_parsing to avoid circular imports
from .vcf_parsing import Variant


def filter_overlapping_variants(
    variants: List[Variant],
    ref_start: int,
    ref_length: int
) -> List[Variant]:
    """Filter variants to remove overlapping ones (keep first).

    Args:
        variants: List of variants (will be sorted by position)
        ref_start: Start position of reference sequence
        ref_length: Length of reference sequence

    Returns:
        Filtered list with non-overlapping variants
    """
    filtered = []
    covered_until = -1

    for var in sorted(variants, key=lambda x: x.pos):
        seq_pos = var.pos - ref_start

        # Skip if outside sequence
        if seq_pos < 0 or seq_pos >= ref_length:
            continue

        # Skip if overlaps with previous variant
        if var.pos <= covered_until:
            continue

        filtered.append(var)
        var_end = var.pos + len(var.ref) - 1
        covered_until = var_end

    return filtered


def apply_variants_to_sequence(
    ref_seq: str,
    ref_start: int,
    variants: List[Variant],
    strain: Optional[str] = None
) -> Tuple[str, List[Variant], Dict[int, int]]:
    """Apply variants to a reference sequence.

    Variants are applied right-to-left to maintain coordinate consistency.
    Builds an offset table tracking cumulative coordinate changes.

    Args:
        ref_seq: Reference DNA sequence
        ref_start: 1-based genomic position of first base in ref_seq
        variants: List of variants to apply
        strain: Optional strain name - if provided, uses strain-specific allele

    Returns:
        Tuple of (modified_sequence, applied_variants, offset_table)
        - offset_table maps original positions to cumulative offset at that position
    """
    # Collect variants to apply
    strain_variants = []

    for var in sorted(variants, key=lambda x: x.pos):
        if strain is not None:
            if strain not in var.samples:
                continue
            strain_allele = var.samples[strain]

            # Skip symbolic alleles
            if strain_allele in ('DEL', 'INS', 'DUP', 'INV', 'CNV', '*'):
                continue
            if strain_allele.startswith('<'):
                continue

            # Only apply if strain differs from reference
            if strain_allele != var.ref:
                strain_variants.append(Variant(
                    chrom=var.chrom,
                    pos=var.pos,
                    ref=var.ref,
                    alt=strain_allele,
                    samples=var.samples
                ))
        else:
            # No strain filtering - apply all variants
            strain_variants.append(var)

    if not strain_variants:
        return ref_seq, [], {}

    # Filter overlapping variants
    filtered_variants = filter_overlapping_variants(
        strain_variants, ref_start, len(ref_seq)
    )

    # Build offset table
    applied = []
    offset_table = {}
    cumulative_offset = 0

    for var in sorted(filtered_variants, key=lambda x: x.pos):
        offset_change = len(var.alt) - len(var.ref)
        cumulative_offset += offset_change
        offset_table[var.pos] = cumulative_offset
        applied.append(var)

    # Apply variants from right to left
    seq_list = list(ref_seq)
    failed_variants = []

    for var in sorted(applied, key=lambda x: -x.pos):
        seq_pos = var.pos - ref_start

        # Check that sequence matches expected reference
        current_ref = ''.join(seq_list[seq_pos:seq_pos + len(var.ref)])
        if current_ref.upper() != var.ref.upper():
            failed_variants.append(var.pos)
            continue

        seq_list[seq_pos:seq_pos + len(var.ref)] = list(var.alt)

    # Rebuild offset table if any variants failed
    if failed_variants:
        applied = [v for v in applied if v.pos not in failed_variants]
        offset_table = {}
        cumulative_offset = 0
        for var in sorted(applied, key=lambda x: x.pos):
            cumulative_offset += len(var.alt) - len(var.ref)
            offset_table[var.pos] = cumulative_offset

    return ''.join(seq_list), applied, offset_table


def insert_variant(
    sequence: str,
    position: int,
    ref: str,
    alt: str,
    seq_start: int = 1
) -> Optional[str]:
    """Insert a single variant into a sequence.

    Args:
        sequence: DNA sequence
        position: 1-based genomic position of variant
        ref: Reference allele
        alt: Alternate allele
        seq_start: 1-based genomic position of first base in sequence

    Returns:
        Modified sequence, or None if variant couldn't be applied
    """
    seq_pos = position - seq_start

    if seq_pos < 0 or seq_pos >= len(sequence):
        return None

    # Check reference matches
    current_ref = sequence[seq_pos:seq_pos + len(ref)]
    if current_ref.upper() != ref.upper():
        return None

    return sequence[:seq_pos] + alt + sequence[seq_pos + len(ref):]


def get_offset_at_position(
    ref_pos: int,
    offset_table: Dict[int, int]
) -> int:
    """Get cumulative offset at a reference position.

    The offset at position P is the cumulative size change from all
    variants at positions < P.

    Args:
        ref_pos: Reference position to query
        offset_table: Dict mapping variant positions to cumulative offsets

    Returns:
        Cumulative offset at the given position
    """
    offset = 0
    for pos, off in sorted(offset_table.items()):
        if pos < ref_pos:
            offset = off
        else:
            break
    return offset


if __name__ == "__main__":
    # Test apply_variants_to_sequence with simple SNP
    ref = "ATGCATGC"
    v = Variant(chrom="chrI", pos=5, ref="A", alt="T", samples={})
    mod, applied, offsets = apply_variants_to_sequence(ref, 1, [v])
    assert mod == "ATGCTTGC", f"Expected ATGCTTGC, got {mod}"
    assert len(applied) == 1
    print("apply_variants SNP: PASS")

    # Test with insertion
    ref = "ATGCATGC"
    v_ins = Variant(chrom="chrI", pos=4, ref="C", alt="CCC", samples={})
    mod, applied, offsets = apply_variants_to_sequence(ref, 1, [v_ins])
    assert mod == "ATGCCCATGC", f"Expected ATGCCCATGC, got {mod}"
    assert offsets[4] == 2  # +2 bases
    print("apply_variants insertion: PASS")

    # Test with deletion
    ref = "ATGCATGC"
    v_del = Variant(chrom="chrI", pos=4, ref="CAT", alt="C", samples={})
    mod, applied, offsets = apply_variants_to_sequence(ref, 1, [v_del])
    assert mod == "ATGCGC", f"Expected ATGCGC, got {mod}"
    assert offsets[4] == -2  # -2 bases
    print("apply_variants deletion: PASS")

    # Test insert_variant
    result = insert_variant("ATGCATGC", 5, "A", "T", seq_start=1)
    assert result == "ATGCTTGC"
    print("insert_variant: PASS")

    # Test insert_variant with wrong ref
    result = insert_variant("ATGCATGC", 5, "G", "T", seq_start=1)
    assert result is None
    print("insert_variant wrong ref: PASS")

    # Test get_offset_at_position
    offsets = {10: 2, 20: 5, 30: 3}
    assert get_offset_at_position(5, offsets) == 0   # before any variant
    assert get_offset_at_position(15, offsets) == 2  # after pos 10
    assert get_offset_at_position(25, offsets) == 5  # after pos 20
    print("get_offset_at_position: PASS")

    print("\nAll variant_ops.py tests passed!")
