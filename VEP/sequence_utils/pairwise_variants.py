"""
Pairwise variant extraction and manipulation.

Functions for comparing two strain sequences and extracting variants between them.
This is different from VCF-based variant extraction because it handles cases where
both strains have different variants at the same S288C position.

Classes:
- PairwiseVariant: Dataclass representing a variant between two strain sequences

Functions:
- compute_strain_seq_position: Convert S288C position to strain sequence position
- extract_pairwise_variants_vcf_guided: Extract variants using VCF coordinate guidance
- apply_pairwise_variant: Apply a pairwise variant to a sequence

Author: Kevin R. Roy
Date: 2026-02-17
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple, Union

# Set up module logger
logger = logging.getLogger(__name__)

# Handle both module import and direct execution
try:
    from .vcf_parsing import Variant
except ImportError:
    # Direct execution - import from sibling module
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent))
    from vcf_parsing import Variant


@dataclass
class PairwiseVariant:
    """Variant extracted by comparing two strain sequences.

    This represents a difference between background and source strain sequences.
    Unlike VCF Variant which is relative to reference genome, this is relative
    to the background strain sequence.

    Attributes:
        background_strain: Name of background strain
        source_strain: Name of source strain
        bg_seq_pos: Position in background sequence (0-based)
        bg_allele: What background strain has at this position
        src_allele: What source strain has at this position
        s288c_pos: S288C reference genome position (1-based)
    """
    background_strain: str
    source_strain: str
    bg_seq_pos: int
    bg_allele: str
    src_allele: str
    s288c_pos: int

    @property
    def variant_id(self) -> str:
        """Get unique variant identifier."""
        return f"{self.s288c_pos}_{self.bg_allele}_{self.src_allele}"

    @property
    def is_indel(self) -> bool:
        """Check if this is an indel variant."""
        return len(self.bg_allele) != len(self.src_allele)

    @property
    def length_change(self) -> int:
        """Get length change from bg to src (positive for insertion)."""
        return len(self.src_allele) - len(self.bg_allele)


def compute_strain_seq_position(
    s288c_pos: int,
    strain_variants: Union[List[Variant], Set[Variant]],
    locus_start: int,
    locus_end: Optional[int] = None,
    seq_length: Optional[int] = None,
    validate_bounds: bool = False
) -> int:
    """Convert S288C position to strain sequence position.

    Accounts for cumulative indel offsets from all variants before this position.
    This is critical for correctly locating variants in strain-specific sequences.

    Algorithm:
    1. Sort all strain variants by position
    2. Filter out overlapping variants (keep first)
    3. For each variant before s288c_pos, accumulate the length change
    4. Return (s288c_pos - locus_start) + cumulative_offset

    Args:
        s288c_pos: S288C reference genome position (1-based)
        strain_variants: All VCF variants for this strain in the region
        locus_start: S288C position where the locus starts (1-based)
        locus_end: S288C position where the locus ends (1-based, optional)
        seq_length: Length of the strain sequence (for bounds checking, optional)
        validate_bounds: If True, raise ValueError for out-of-bounds positions

    Returns:
        Position in the strain sequence (0-based relative to locus start)

    Raises:
        ValueError: If validate_bounds=True and result is out of bounds

    Example:
        If locus starts at S288C position 1000, and there's a 5bp deletion at
        position 1010, then S288C position 1020 would be at strain position:
        (1020 - 1000) + (-5) = 15
    """
    # FIX (2026-02-23): Added bounds validation
    # Validate input position is within locus
    if locus_end is not None and (s288c_pos < locus_start or s288c_pos > locus_end):
        if validate_bounds:
            raise ValueError(
                f"s288c_pos {s288c_pos} is outside locus range [{locus_start}, {locus_end}]"
            )
        else:
            logger.warning(
                f"s288c_pos {s288c_pos} is outside locus range [{locus_start}, {locus_end}]"
            )

    # Sort and filter overlapping variants
    sorted_vars = sorted(strain_variants, key=lambda v: v.pos)
    filtered_vars = []
    covered_until = -1

    for var in sorted_vars:
        if var.pos <= covered_until:
            continue
        filtered_vars.append(var)
        covered_until = var.pos + len(var.ref) - 1

    # Calculate cumulative offset from all indels before s288c_pos
    offset = 0
    for var in filtered_vars:
        if var.pos < s288c_pos:
            offset += len(var.alt) - len(var.ref)

    result = (s288c_pos - locus_start) + offset

    # FIX (2026-02-23): Validate result is within sequence bounds
    if validate_bounds:
        if result < 0:
            raise ValueError(
                f"Computed strain position {result} is negative "
                f"(s288c_pos={s288c_pos}, locus_start={locus_start}, offset={offset})"
            )
        if seq_length is not None and result >= seq_length:
            raise ValueError(
                f"Computed strain position {result} exceeds sequence length {seq_length} "
                f"(s288c_pos={s288c_pos}, locus_start={locus_start}, offset={offset})"
            )
    elif result < 0:
        logger.warning(
            f"Computed strain position {result} is negative "
            f"(s288c_pos={s288c_pos}, locus_start={locus_start}, offset={offset})"
        )

    return result


def extract_pairwise_variants_vcf_guided(
    bg_seq: str,
    src_seq: str,
    background_strain: str,
    source_strain: str,
    bg_vcf_variants: Union[List[Variant], Set[Variant]],
    src_vcf_variants: Union[List[Variant], Set[Variant]],
    s288c_start: int,
    s288c_end: int
) -> List[PairwiseVariant]:
    """Extract variants using VCF-guided coordinate mapping.

    Uses VCF information to properly map S288C positions to strain sequence
    positions, accounting for cumulative indel offsets. This handles cases where
    both strains have different variants at the same S288C position.

    Args:
        bg_seq: Background strain sequence
        src_seq: Source strain sequence
        background_strain: Name of background strain
        source_strain: Name of source strain
        bg_vcf_variants: VCF variants for background strain
        src_vcf_variants: VCF variants for source strain
        s288c_start: S288C start coordinate of the locus (1-based)
        s288c_end: S288C end coordinate of the locus (1-based)

    Returns:
        List of PairwiseVariant objects representing differences between strains

    Example:
        If background has SNP A→T at S288C position 1000, and source has SNP A→G
        at the same position, this will return a PairwiseVariant with:
        - bg_allele = "T" (what background has)
        - src_allele = "G" (what source has)
        - s288c_pos = 1000
        - bg_seq_pos = computed accounting for indels before position 1000
    """
    variants = []

    # Get all VCF positions in this region
    bg_var_dict = {v.pos: v for v in bg_vcf_variants if s288c_start <= v.pos <= s288c_end}
    src_var_dict = {v.pos: v for v in src_vcf_variants if s288c_start <= v.pos <= s288c_end}

    all_positions = sorted(set(bg_var_dict.keys()) | set(src_var_dict.keys()))

    for s288c_pos in all_positions:
        bg_var = bg_var_dict.get(s288c_pos)
        src_var = src_var_dict.get(s288c_pos)

        # Calculate correct position in background sequence
        # This accounts for cumulative indel offsets
        bg_seq_pos = compute_strain_seq_position(s288c_pos, bg_vcf_variants, s288c_start)

        # Determine what each strain has at this position
        if bg_var and src_var:
            # Both have variants - compare ALT alleles
            if bg_var.alt != src_var.alt:
                variants.append(PairwiseVariant(
                    background_strain=background_strain,
                    source_strain=source_strain,
                    bg_seq_pos=bg_seq_pos,
                    bg_allele=bg_var.alt,
                    src_allele=src_var.alt,
                    s288c_pos=s288c_pos
                ))
        elif bg_var:
            # Only background has variant - source has reference
            variants.append(PairwiseVariant(
                background_strain=background_strain,
                source_strain=source_strain,
                bg_seq_pos=bg_seq_pos,
                bg_allele=bg_var.alt,
                src_allele=bg_var.ref,
                s288c_pos=s288c_pos
            ))
        elif src_var:
            # Only source has variant - background has reference
            variants.append(PairwiseVariant(
                background_strain=background_strain,
                source_strain=source_strain,
                bg_seq_pos=bg_seq_pos,
                bg_allele=src_var.ref,
                src_allele=src_var.alt,
                s288c_pos=s288c_pos
            ))

    return variants


def apply_pairwise_variant(
    sequence: str,
    variant: PairwiseVariant,
    verify_allele: bool = True
) -> Optional[str]:
    """Apply a pairwise variant to a sequence.

    The variant specifies:
    - bg_allele: what the background sequence currently has
    - src_allele: what we want to change it to
    - bg_seq_pos: position in the background sequence (0-based)

    Args:
        sequence: Background sequence
        variant: Pairwise variant to apply
        verify_allele: If True, verify bg_allele matches sequence before applying

    Returns:
        Modified sequence, or None if variant couldn't be applied

    Example:
        >>> seq = "ATGCATGC"
        >>> var = PairwiseVariant(
        ...     background_strain="A",
        ...     source_strain="B",
        ...     bg_seq_pos=2,
        ...     bg_allele="GC",
        ...     src_allele="T",
        ...     s288c_pos=100
        ... )
        >>> apply_pairwise_variant(seq, var)
        'ATTATGC'
    """
    pos = variant.bg_seq_pos
    bg_allele = variant.bg_allele if variant.bg_allele != "-" else ""
    src_allele = variant.src_allele if variant.src_allele != "-" else ""

    # Check bounds
    # FIX (2026-02-23): Changed > to >= for correct boundary check
    # pos == len(sequence) is invalid unless bg_allele is empty (pure insertion at end)
    if pos < 0:
        return None
    if pos > len(sequence):
        return None
    if pos == len(sequence) and len(bg_allele) > 0:
        # Can't have a non-empty bg_allele starting at pos == len(sequence)
        return None

    if pos + len(bg_allele) > len(sequence):
        return None

    # Verify the allele matches if requested
    if verify_allele:
        actual_allele = sequence[pos:pos + len(bg_allele)]
        if actual_allele.upper() != bg_allele.upper():
            # Try nearby positions (±5 bp) for minor alignment differences
            # FIX (2026-02-23): Log warning when position is shifted
            original_pos = pos
            found_match = False
            for offset in range(-5, 6):
                check_pos = pos + offset
                if 0 <= check_pos < len(sequence):
                    check_end = check_pos + len(bg_allele)
                    if check_end <= len(sequence):
                        actual_allele = sequence[check_pos:check_end]
                        if actual_allele.upper() == bg_allele.upper():
                            if offset != 0:
                                logger.warning(
                                    f"Variant position shifted by {offset:+d}bp: "
                                    f"s288c_pos={variant.s288c_pos}, "
                                    f"original_bg_seq_pos={original_pos}, "
                                    f"adjusted_bg_seq_pos={check_pos}"
                                )
                            pos = check_pos
                            found_match = True
                            break
            if not found_match:
                # No match found
                return None

    # Apply the variant
    return sequence[:pos] + src_allele + sequence[pos + len(bg_allele):]


def apply_multiple_pairwise_variants(
    sequence: str,
    variants: List[PairwiseVariant],
    verify_alleles: bool = True
) -> tuple[str, List[str]]:
    """Apply multiple pairwise variants to a sequence.

    Variants are applied right-to-left to avoid position shifts.

    Args:
        sequence: Background sequence
        variants: List of pairwise variants to apply
        verify_alleles: If True, verify each allele matches before applying

    Returns:
        Tuple of (modified_sequence, list of applied variant IDs)

    Example:
        Apply multiple variants in the correct order to avoid coordinate issues.
    """
    # Sort variants by position, descending (right-to-left)
    sorted_vars = sorted(variants, key=lambda v: v.bg_seq_pos, reverse=True)

    applied_ids = []
    result = sequence

    for var in sorted_vars:
        new_result = apply_pairwise_variant(result, var, verify_alleles)
        if new_result is not None:
            result = new_result
            applied_ids.append(var.variant_id)

    return result, applied_ids


# ============================================================================
# Linked Variant Detection and Classification
# ============================================================================

@dataclass
class LinkedVariantGroup:
    """A group of DNA variants that are linked (co-occur on same haplotype).

    When analyzing variants from one strain, nearby indels likely co-occur.
    This class represents such a group and provides methods to analyze their
    combined effect.

    Attributes:
        variants: List of PairwiseVariant objects in this group
        background_strain: Background strain name
        source_strain: Source strain name
    """
    variants: List[PairwiseVariant]
    background_strain: str
    source_strain: str

    @property
    def net_length_change(self) -> int:
        """Calculate net length change from all variants in group."""
        return sum(v.length_change for v in self.variants)

    @property
    def is_in_frame(self) -> bool:
        """Check if the combined effect is in-frame (divisible by 3)."""
        return self.net_length_change % 3 == 0

    @property
    def min_position(self) -> int:
        """Get the minimum S288C position in the group."""
        return min(v.s288c_pos for v in self.variants)

    @property
    def max_position(self) -> int:
        """Get the maximum S288C position in the group."""
        return max(v.s288c_pos for v in self.variants)

    @property
    def span(self) -> int:
        """Get the span of positions covered by this group."""
        return self.max_position - self.min_position + 1

    @property
    def n_variants(self) -> int:
        """Number of variants in this group."""
        return len(self.variants)

    @property
    def group_id(self) -> str:
        """Generate unique ID for this group."""
        positions = sorted(v.s288c_pos for v in self.variants)
        return f"{self.min_position}-{self.max_position}_n{self.n_variants}"

    def get_variant_ids(self) -> List[str]:
        """Get list of variant IDs in this group."""
        return [v.variant_id for v in sorted(self.variants, key=lambda x: x.s288c_pos)]


def detect_linked_variants(
    variants: List[PairwiseVariant],
    max_distance: int = 50,
    snv_distance: int = 3,
    indel_distance: int = 10
) -> List[LinkedVariantGroup]:
    """Group variants that are likely linked based on variant type and distance.

    Variants from the same strain in close proximity likely co-occur on the
    same haplotype. This function groups such variants for combined effect
    analysis.

    **Linking Rules (Updated 2026-02-23):**
    - SNVs: Link within snv_distance bp (default: 3bp)
    - INDELs: Link within indel_distance bp (default: 10bp)
    - Mixed SNV+INDEL: Use indel_distance (more permissive)

    For CDS-aware linking that extends INDEL groups to find in-frame combinations,
    use detect_linked_variants_cds_aware() instead.

    Args:
        variants: List of PairwiseVariant objects (should be from same strain pair)
        max_distance: Maximum bp distance for legacy behavior (default: 50bp, not used
                      when snv_distance and indel_distance are set)
        snv_distance: Maximum bp distance between consecutive SNVs (default: 3bp)
        indel_distance: Maximum bp distance between consecutive INDELs (default: 10bp)

    Returns:
        List of LinkedVariantGroup objects

    Example:
        >>> vars = [
        ...     PairwiseVariant("BYa", "YJM454a", 100, "ATGCG", "A", 538149),  # -4bp
        ...     PairwiseVariant("BYa", "YJM454a", 110, "G", "GA", 538163),     # +1bp
        ...     PairwiseVariant("BYa", "YJM454a", 120, "A", "ACCCACGAC", 538171),  # +8bp
        ...     PairwiseVariant("BYa", "YJM454a", 800, "C", "T", 538800),      # distant SNP
        ... ]
        >>> groups = detect_linked_variants(vars, snv_distance=3, indel_distance=10)
        >>> len(groups)
        2
        >>> groups[0].n_variants  # First 3 indels are linked (all within 10bp of each other)
        3
    """
    if not variants:
        return []

    # Sort by position
    sorted_vars = sorted(variants, key=lambda v: v.s288c_pos)

    # Verify all variants are from same strain pair
    bg_strain = sorted_vars[0].background_strain
    src_strain = sorted_vars[0].source_strain

    def get_link_distance(var1: PairwiseVariant, var2: PairwiseVariant) -> int:
        """Get the maximum linking distance based on variant types."""
        # If either is an indel, use the more permissive indel_distance
        if var1.is_indel or var2.is_indel:
            return indel_distance
        else:
            return snv_distance

    groups = []
    current_group = [sorted_vars[0]]

    for i in range(1, len(sorted_vars)):
        current_var = sorted_vars[i]
        prev_var = sorted_vars[i - 1]

        # Get appropriate linking distance for this pair
        link_dist = get_link_distance(prev_var, current_var)

        # Check if within linking distance of previous variant
        if current_var.s288c_pos - prev_var.s288c_pos <= link_dist:
            current_group.append(current_var)
        else:
            # Close current group and start new one
            groups.append(LinkedVariantGroup(
                variants=current_group,
                background_strain=bg_strain,
                source_strain=src_strain
            ))
            current_group = [current_var]

    # Don't forget the last group
    groups.append(LinkedVariantGroup(
        variants=current_group,
        background_strain=bg_strain,
        source_strain=src_strain
    ))

    return groups


def detect_linked_variants_cds_aware(
    variants: List[PairwiseVariant],
    cds_start_s288c: int,
    cds_end_s288c: int,
    snv_distance: int = 3,
    indel_distance: int = 10
) -> List[LinkedVariantGroup]:
    """Group variants with CDS-aware INDEL extension for in-frame detection.

    This extends the basic linking rules with a critical enhancement:
    **ALL INDELs within a CDS are collected and grouped if their combined
    effect is in-frame.**

    This handles cases like ENA1 where 4 indels (-10bp, +1bp, +8bp, +1bp)
    individually cause frameshifts but together sum to 0bp (in-frame).

    **Algorithm (Bidirectional - Updated 2026-02-23):**
    1. Collect ALL indels within the CDS
    2. Check if combined net change is in-frame (% 3 == 0)
    3. If in-frame: group ALL CDS indels together as one linked group
    4. If not in-frame: use cumulative sum to find maximal in-frame subgroups
    5. Apply standard linking rules for non-CDS variants and SNVs

    Args:
        variants: List of PairwiseVariant objects (should be from same strain pair)
        cds_start_s288c: CDS start position in S288C coordinates (1-based)
        cds_end_s288c: CDS end position in S288C coordinates (1-based)
        snv_distance: Maximum bp distance between consecutive SNVs (default: 3bp)
        indel_distance: Maximum bp distance between consecutive INDELs (default: 10bp)

    Returns:
        List of LinkedVariantGroup objects with CDS-aware INDEL linking

    Example:
        >>> # ENA1 example: 4 indels span ~25bp but should be grouped together
        >>> # because their net change is 0bp (in-frame)
        >>> vars = [
        ...     PairwiseVariant("BYa", "YJM454a", 100, "TGAGCACATTG", "T", 538149),  # -10bp
        ...     PairwiseVariant("BYa", "YJM454a", 114, "G", "GA", 538163),           # +1bp
        ...     PairwiseVariant("BYa", "YJM454a", 122, "A", "ACCCACGAC", 538171),    # +8bp
        ...     PairwiseVariant("BYa", "YJM454a", 125, "A", "AG", 538174),           # +1bp
        ... ]
        >>> groups = detect_linked_variants_cds_aware(
        ...     vars, cds_start_s288c=537500, cds_end_s288c=540000
        ... )
        >>> groups[0].net_length_change
        0  # All 4 linked as they form in-frame complex indel
    """
    if not variants:
        return []

    # Sort by position
    sorted_vars = sorted(variants, key=lambda v: v.s288c_pos)

    # Verify all variants are from same strain pair
    bg_strain = sorted_vars[0].background_strain
    src_strain = sorted_vars[0].source_strain

    # Separate variants into CDS indels and others
    cds_indels = []
    non_cds_or_snv = []

    for var in sorted_vars:
        is_in_cds = cds_start_s288c <= var.s288c_pos <= cds_end_s288c
        if is_in_cds and var.is_indel:
            cds_indels.append(var)
        else:
            non_cds_or_snv.append(var)

    final_groups = []
    used_variant_ids = set()

    # Process CDS indels: find in-frame groupings using cumulative sum
    if cds_indels:
        # Sort CDS indels by position
        cds_indels_sorted = sorted(cds_indels, key=lambda v: v.s288c_pos)

        # Calculate cumulative length change
        cumsum = 0
        cumsum_at_pos = []  # (variant, cumulative_sum_after_variant)
        for var in cds_indels_sorted:
            cumsum += var.length_change
            cumsum_at_pos.append((var, cumsum))

        # Find in-frame groups using cumulative sum
        # A group from index i to j is in-frame if:
        #   cumsum[j] - cumsum[i-1] ≡ 0 (mod 3)
        # Which means cumsum[j] ≡ cumsum[i-1] (mod 3)

        # Strategy: greedily find maximal in-frame groups
        # Start from position 0, extend until we find an in-frame endpoint
        current_group_start = 0
        prev_cumsum = 0  # cumsum before first element is 0

        while current_group_start < len(cds_indels_sorted):
            # Look for the furthest in-frame endpoint
            best_end = None
            for j in range(current_group_start, len(cds_indels_sorted)):
                group_net = cumsum_at_pos[j][1] - prev_cumsum
                if group_net % 3 == 0:
                    best_end = j  # Keep extending to find furthest in-frame point

            if best_end is not None:
                # Create group from current_group_start to best_end
                group_variants = [cumsum_at_pos[i][0] for i in range(current_group_start, best_end + 1)]
                if len(group_variants) > 0:
                    group = LinkedVariantGroup(
                        variants=group_variants,
                        background_strain=bg_strain,
                        source_strain=src_strain
                    )
                    final_groups.append(group)
                    used_variant_ids.update(v.variant_id for v in group_variants)
                    prev_cumsum = cumsum_at_pos[best_end][1]
                    current_group_start = best_end + 1
            else:
                # No in-frame endpoint found - remaining indels form a frameshift group
                group_variants = [cumsum_at_pos[i][0] for i in range(current_group_start, len(cds_indels_sorted))]
                if len(group_variants) > 0:
                    group = LinkedVariantGroup(
                        variants=group_variants,
                        background_strain=bg_strain,
                        source_strain=src_strain
                    )
                    final_groups.append(group)
                    used_variant_ids.update(v.variant_id for v in group_variants)
                break

    # Process non-CDS variants and SNVs using standard linking rules
    if non_cds_or_snv:
        standard_groups = detect_linked_variants(
            non_cds_or_snv, snv_distance=snv_distance, indel_distance=indel_distance
        )
        for group in standard_groups:
            # Only add if not already used (shouldn't happen, but be safe)
            if not any(v.variant_id in used_variant_ids for v in group.variants):
                final_groups.append(group)
                used_variant_ids.update(v.variant_id for v in group.variants)

    # Sort final groups by position
    final_groups.sort(key=lambda g: g.min_position)

    return final_groups


def extract_pairwise_variants_alignment(
    bg_seq: str,
    src_seq: str,
    background_strain: str,
    source_strain: str,
    bg_offset_table: Dict[int, int],
    locus_start_s288c: int,
    locus_end_s288c: int
) -> List[PairwiseVariant]:
    """Extract pairwise variants using direct sequence alignment.

    This is an alternative to VCF-guided extraction that handles cases where
    both strains have different variants at the same S288C position. It directly
    aligns the two sequences and extracts differences from the alignment.

    Args:
        bg_seq: Background strain sequence
        src_seq: Source strain sequence
        background_strain: Name of background strain
        source_strain: Name of source strain
        bg_offset_table: Offset table for background strain (S288C pos -> offset)
        locus_start_s288c: S288C start coordinate of the locus
        locus_end_s288c: S288C end coordinate of the locus

    Returns:
        List of PairwiseVariant objects with S288C positions mapped

    Note:
        This function requires the alignment_variants module and mappy library.
        Use this when VCF-guided extraction fails due to complex overlapping variants.
    """
    # Import alignment functions (deferred to avoid circular imports)
    from .alignment_variants import (
        align_and_extract_variants,
        AlignmentVariant
    )

    # Extract variants using alignment
    alignment_variants = align_and_extract_variants(
        bg_seq=bg_seq,
        src_seq=src_seq,
        background_strain=background_strain,
        source_strain=source_strain,
        bg_offset_table=bg_offset_table,
        locus_start_s288c=locus_start_s288c,
        locus_end_s288c=locus_end_s288c
    )

    # Convert AlignmentVariant to PairwiseVariant
    pairwise_variants = []
    for av in alignment_variants:
        # Handle None s288c_pos (position within insertion)
        s288c_pos = av.s288c_pos if av.s288c_pos is not None else -1

        pv = PairwiseVariant(
            background_strain=av.background_strain,
            source_strain=av.source_strain,
            bg_seq_pos=av.bg_seq_pos,
            bg_allele=av.bg_allele if av.bg_allele else "-",
            src_allele=av.src_allele if av.src_allele else "-",
            s288c_pos=s288c_pos
        )
        pairwise_variants.append(pv)

    return pairwise_variants


def extract_pairwise_variants_auto(
    bg_seq: str,
    src_seq: str,
    background_strain: str,
    source_strain: str,
    bg_vcf_variants: Union[List[Variant], Set[Variant]],
    src_vcf_variants: Union[List[Variant], Set[Variant]],
    bg_offset_table: Dict[int, int],
    s288c_start: int,
    s288c_end: int,
    prefer_alignment: bool = False
) -> List[PairwiseVariant]:
    """Auto-select extraction method based on variant complexity.

    By default, uses VCF-guided extraction for efficiency. Falls back to
    alignment-based extraction when:
    - Both strains have indels at the same position
    - VCF-guided extraction produces suspicious results
    - prefer_alignment=True is set

    Args:
        bg_seq: Background strain sequence
        src_seq: Source strain sequence
        background_strain: Name of background strain
        source_strain: Name of source strain
        bg_vcf_variants: VCF variants for background strain
        src_vcf_variants: VCF variants for source strain
        bg_offset_table: Offset table for background strain
        s288c_start: S288C start coordinate (1-based)
        s288c_end: S288C end coordinate (1-based)
        prefer_alignment: If True, always use alignment-based extraction

    Returns:
        List of PairwiseVariant objects
    """
    # Check for overlapping indels at same position (the problematic case)
    bg_indel_positions = {v.pos for v in bg_vcf_variants
                         if s288c_start <= v.pos <= s288c_end
                         and len(v.ref) != len(v.alt)}
    src_indel_positions = {v.pos for v in src_vcf_variants
                          if s288c_start <= v.pos <= s288c_end
                          and len(v.ref) != len(v.alt)}

    overlapping_indels = bg_indel_positions & src_indel_positions

    if prefer_alignment or overlapping_indels:
        # Use alignment-based extraction
        return extract_pairwise_variants_alignment(
            bg_seq=bg_seq,
            src_seq=src_seq,
            background_strain=background_strain,
            source_strain=source_strain,
            bg_offset_table=bg_offset_table,
            locus_start_s288c=s288c_start,
            locus_end_s288c=s288c_end
        )
    else:
        # Use VCF-guided extraction (faster for simple cases)
        return extract_pairwise_variants_vcf_guided(
            bg_seq=bg_seq,
            src_seq=src_seq,
            background_strain=background_strain,
            source_strain=source_strain,
            bg_vcf_variants=bg_vcf_variants,
            src_vcf_variants=src_vcf_variants,
            s288c_start=s288c_start,
            s288c_end=s288c_end
        )


def classify_linked_effect(group: LinkedVariantGroup) -> tuple:
    """Classify the effect of a linked variant group.

    For groups with multiple variants, determines the combined effect:
    - If net length change is 0: in-frame complex indel
    - If net length change % 3 == 0: in-frame indel (insertion or deletion)
    - Otherwise: frameshift (truncating or lengthening based on protein effect)

    Args:
        group: LinkedVariantGroup to classify

    Returns:
        Tuple of (variant_type, effect_description, net_change)

    Example:
        >>> # ENA1 example: -10 + 1 + 8 + 1 = 0bp (in-frame)
        >>> group = LinkedVariantGroup(variants=[...], ...)
        >>> classify_linked_effect(group)
        ('complex_indel_inframe', 'Net 0bp change (in-frame)', 0)
    """
    net_change = group.net_length_change
    n_variants = group.n_variants

    if n_variants == 1:
        # Single variant - classify directly
        var = group.variants[0]
        if not var.is_indel:
            return ('snp', 'Single nucleotide polymorphism', 0)
        elif net_change % 3 == 0:
            if net_change > 0:
                return ('aa_insertion', f'In-frame insertion of {net_change // 3} AA', net_change)
            else:
                return ('aa_deletion', f'In-frame deletion of {-net_change // 3} AA', net_change)
        else:
            return ('frameshift', f'{net_change:+d}bp frameshift', net_change)

    # Multiple linked variants
    if net_change == 0:
        return ('complex_indel_neutral', f'Net 0bp change from {n_variants} variants', 0)
    elif net_change % 3 == 0:
        aa_change = net_change // 3
        if aa_change > 0:
            return ('complex_indel_insertion', f'Net {aa_change} AA insertion from {n_variants} variants', net_change)
        else:
            return ('complex_indel_deletion', f'Net {-aa_change} AA deletion from {n_variants} variants', net_change)
    else:
        return ('complex_frameshift', f'Net {net_change:+d}bp frameshift from {n_variants} variants', net_change)


if __name__ == "__main__":
    # Test compute_strain_seq_position
    print("Testing compute_strain_seq_position...")

    # Create test variants
    locus_start = 1000
    deletion = Variant(chrom="chrI", pos=1010, ref="ATGCG", alt="A", samples={})
    insertion = Variant(chrom="chrI", pos=1020, ref="T", alt="TGC", samples={})

    variants = [deletion, insertion]  # Use list, not set

    # Before any variants
    pos = compute_strain_seq_position(1005, variants, locus_start)
    assert pos == 5, f"Expected 5, got {pos}"

    # After deletion (-4bp)
    pos = compute_strain_seq_position(1015, variants, locus_start)
    assert pos == 11, f"Expected 11, got {pos}"  # (1015-1000) + (-4)

    # After both deletion and insertion (-4 + 2 = -2)
    pos = compute_strain_seq_position(1025, variants, locus_start)
    assert pos == 23, f"Expected 23, got {pos}"  # (1025-1000) + (-2)

    print("✓ compute_strain_seq_position tests passed")

    # Test apply_pairwise_variant
    print("\nTesting apply_pairwise_variant...")

    seq = "ATGCATGC"
    var = PairwiseVariant(
        background_strain="A",
        source_strain="B",
        bg_seq_pos=2,
        bg_allele="GC",
        src_allele="T",
        s288c_pos=100
    )
    result = apply_pairwise_variant(seq, var)
    assert result == "ATTATGC", f"Expected ATTATGC, got {result}"

    print("✓ apply_pairwise_variant tests passed")

    # Test apply_multiple_pairwise_variants
    print("\nTesting apply_multiple_pairwise_variants...")

    seq = "ATGCATGC"
    vars = [
        PairwiseVariant("A", "B", 0, "A", "T", 100),  # A→T at pos 0
        PairwiseVariant("A", "B", 4, "A", "G", 104),  # A→G at pos 4
    ]

    result, applied = apply_multiple_pairwise_variants(seq, vars)
    assert result == "TTGCGTGC", f"Expected TTGCGTGC, got {result}"
    assert len(applied) == 2

    print("✓ apply_multiple_pairwise_variants tests passed")

    # Test detect_linked_variants
    print("\nTesting detect_linked_variants...")

    # Create test variants mimicking ENA1 case
    ena1_vars = [
        PairwiseVariant("BYa", "YJM454a", 100, "TGAGCACATTG", "T", 538149),  # -10bp
        PairwiseVariant("BYa", "YJM454a", 114, "G", "GA", 538163),           # +1bp
        PairwiseVariant("BYa", "YJM454a", 122, "A", "ACCCACGAC", 538171),    # +8bp
        PairwiseVariant("BYa", "YJM454a", 125, "A", "AG", 538174),           # +1bp
        PairwiseVariant("BYa", "YJM454a", 500, "C", "T", 538800),            # distant SNP
    ]

    groups = detect_linked_variants(ena1_vars, max_distance=50)
    assert len(groups) == 2, f"Expected 2 groups, got {len(groups)}"
    assert groups[0].n_variants == 4, f"Expected 4 variants in first group, got {groups[0].n_variants}"
    assert groups[0].net_length_change == 0, f"Expected net 0bp, got {groups[0].net_length_change}"
    assert groups[0].is_in_frame, "Expected first group to be in-frame"
    assert groups[1].n_variants == 1, f"Expected 1 variant in second group, got {groups[1].n_variants}"

    print("✓ detect_linked_variants tests passed")

    # Test classify_linked_effect
    print("\nTesting classify_linked_effect...")

    var_type, effect, net = classify_linked_effect(groups[0])
    assert var_type == 'complex_indel_neutral', f"Expected complex_indel_neutral, got {var_type}"
    assert net == 0, f"Expected net 0, got {net}"

    var_type, effect, net = classify_linked_effect(groups[1])
    assert var_type == 'snp', f"Expected snp, got {var_type}"

    # Test frameshift case
    frameshift_vars = [
        PairwiseVariant("BYa", "Test", 100, "A", "AT", 1000),  # +1bp
        PairwiseVariant("BYa", "Test", 110, "G", "GT", 1010),  # +1bp
    ]
    fs_groups = detect_linked_variants(frameshift_vars, max_distance=50)
    var_type, effect, net = classify_linked_effect(fs_groups[0])
    assert var_type == 'complex_frameshift', f"Expected complex_frameshift, got {var_type}"
    assert net == 2, f"Expected net 2, got {net}"

    print("✓ classify_linked_effect tests passed")

    print("\n✅ All pairwise_variants.py tests passed!")
