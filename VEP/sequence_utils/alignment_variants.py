"""
Alignment-based pairwise variant extraction.

This module extracts variants by DIRECTLY ALIGNING two strain sequences,
rather than relying on VCF-based offset calculations. This correctly handles
cases where both strains have different variants at the same S288C position.

The key insight:
- VCF-based approach fails when L has deletion at pos X and S has SNP at pos X
  because position X doesn't exist in L's sequence
- Alignment-based approach aligns L directly to S, finds the actual differences,
  then maps positions back to S288C coordinates using offset tables

Classes:
- AlignmentVariant: Variant extracted from pairwise sequence alignment

Functions:
- align_and_extract_variants: Main function for alignment-based extraction
- reverse_offset_lookup: Map strain sequence position back to S288C coordinate
- load_offset_table: Load offset table from JSON file

Author: Kevin R. Roy
Date: 2026-02-22
"""

import json
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import mappy as mp


@dataclass
class AlignmentVariant:
    """Variant extracted from direct pairwise sequence alignment.

    Unlike VCF-guided variants, these are derived from actual sequence
    differences identified by alignment. The S288C position is computed
    by reverse offset lookup, enabling cross-strain comparison.

    Attributes:
        background_strain: Name of background strain (alignment reference)
        source_strain: Name of source strain (aligned query)
        bg_seq_pos: Position in background sequence (0-based)
        bg_allele: What background has at this position
        src_allele: What source has at this position
        s288c_pos: Corresponding S288C coordinate (1-based), may be None if
                   position is within an insertion relative to S288C
        variant_type: 'SNP', 'INS', 'DEL', or 'COMPLEX'
        alignment_quality: Local alignment quality (mapq or derived score)
    """
    background_strain: str
    source_strain: str
    bg_seq_pos: int
    bg_allele: str
    src_allele: str
    s288c_pos: Optional[int]
    variant_type: str
    alignment_quality: float = 60.0

    @property
    def variant_id(self) -> str:
        """Generate unique ID for this variant."""
        s288c = self.s288c_pos if self.s288c_pos else f"ins{self.bg_seq_pos}"
        return f"{s288c}_{self.bg_allele}_{self.src_allele}"

    @property
    def is_snp(self) -> bool:
        return len(self.bg_allele) == 1 and len(self.src_allele) == 1

    @property
    def is_insertion(self) -> bool:
        """Source has more bases than background."""
        return len(self.src_allele) > len(self.bg_allele)

    @property
    def is_deletion(self) -> bool:
        """Source has fewer bases than background."""
        return len(self.src_allele) < len(self.bg_allele)

    @property
    def length_change(self) -> int:
        """Net length change if variant applied to background."""
        return len(self.src_allele) - len(self.bg_allele)


@lru_cache(maxsize=128)
def _load_offset_table_cached(offset_file_str: str) -> Tuple[Tuple[int, int], ...]:
    """Internal cached function for loading offset tables.

    FIX (2026-02-23): Added caching to avoid redundant file I/O.
    Uses tuple of tuples (hashable) for lru_cache compatibility.

    Args:
        offset_file_str: String path to offset JSON file

    Returns:
        Tuple of (s288c_pos, offset) pairs for cache storage
    """
    with open(offset_file_str, 'r') as f:
        raw = json.load(f)

    # Convert to sorted tuple of tuples (hashable for caching)
    return tuple(sorted((int(k), v) for k, v in raw.items()))


def load_offset_table(offset_file: Path) -> Dict[int, int]:
    """Load offset table from JSON file.

    The offset table maps S288C positions to cumulative indel offsets.
    Format: {"s288c_pos": cumulative_offset, ...}

    Results are cached using lru_cache to avoid redundant file I/O
    when the same offset table is loaded multiple times.

    Args:
        offset_file: Path to offset JSON file

    Returns:
        Dict mapping S288C positions (int) to cumulative offsets (int)
    """
    # Use cached version with string path (hashable)
    cached_tuple = _load_offset_table_cached(str(offset_file))
    return dict(cached_tuple)


def clear_offset_table_cache():
    """Clear the offset table cache.

    Useful for testing or when offset table files have been updated.
    """
    _load_offset_table_cached.cache_clear()


def get_offset_table_cache_info():
    """Get cache statistics for offset table loading.

    Returns:
        NamedTuple with hits, misses, maxsize, and currsize
    """
    return _load_offset_table_cached.cache_info()


def s288c_to_strain_pos(
    s288c_pos: int,
    offset_table: Dict[int, int],
    locus_start_s288c: int
) -> int:
    """Convert S288C position to strain sequence position.

    Uses the offset table to account for cumulative indels.

    Args:
        s288c_pos: S288C coordinate (1-based)
        offset_table: Maps S288C positions to cumulative offsets
        locus_start_s288c: S288C start coordinate of the locus

    Returns:
        Position in strain sequence (0-based, relative to locus start)
    """
    # Get all offset positions before or at s288c_pos
    offset_positions = sorted([p for p in offset_table.keys() if p <= s288c_pos])

    # Get cumulative offset
    if offset_positions:
        # Use the offset at the last position before s288c_pos
        offset = offset_table[offset_positions[-1]]
    else:
        offset = 0

    return (s288c_pos - locus_start_s288c) + offset


def reverse_offset_lookup(
    strain_pos: int,
    offset_table: Dict[int, int],
    locus_start_s288c: int,
    locus_end_s288c: int
) -> Optional[int]:
    """Convert strain sequence position back to S288C coordinate.

    Given a position in the strain sequence, find the corresponding
    S288C coordinate. This is the inverse of s288c_to_strain_pos.

    Algorithm:
    1. Binary search for s288c_pos where s288c_to_strain_pos(s288c_pos) == strain_pos
    2. Handle edge cases (positions within insertions relative to S288C)

    Args:
        strain_pos: Position in strain sequence (0-based, relative to locus start)
        offset_table: Maps S288C positions to cumulative offsets
        locus_start_s288c: S288C start coordinate of the locus
        locus_end_s288c: S288C end coordinate of the locus

    Returns:
        S288C coordinate (1-based), or None if position is within an insertion
        that doesn't correspond to any S288C position
    """
    # Binary search
    lo = locus_start_s288c
    hi = locus_end_s288c + 1000  # Buffer for insertions

    while lo < hi:
        mid = (lo + hi) // 2
        mapped_strain_pos = s288c_to_strain_pos(mid, offset_table, locus_start_s288c)

        if mapped_strain_pos < strain_pos:
            lo = mid + 1
        elif mapped_strain_pos > strain_pos:
            hi = mid
        else:
            # Found it
            return mid

    # Check if lo is close enough
    final_mapped = s288c_to_strain_pos(lo, offset_table, locus_start_s288c)
    if final_mapped == strain_pos:
        return lo

    # Position might be within an insertion (no corresponding S288C position)
    # Return the closest S288C position
    return lo if lo <= locus_end_s288c else None


def build_reverse_mapping(
    offset_table: Dict[int, int],
    locus_start_s288c: int,
    locus_end_s288c: int,
    strain_seq_len: int
) -> Dict[int, Optional[int]]:
    """Pre-build reverse mapping from strain positions to S288C positions.

    For efficiency, builds a complete mapping for all positions in the locus.

    IMPORTANT: When there's a deletion, multiple S288C positions may map to
    the same strain position (collision). We keep only the FIRST (earliest)
    S288C position, since later colliding positions represent deleted bases
    that don't actually exist in the strain sequence.

    Example: offset_table = {10: -2} means 2bp deleted at S288C positions 10-11
    - S288C 8 → strain 7 (first mapping, kept)
    - S288C 9 → strain 8 (first mapping, kept)
    - S288C 10 → strain 7 (collision with 8, discarded - position is deleted)
    - S288C 11 → strain 8 (collision with 9, discarded - position is deleted)
    - S288C 12 → strain 9 (first mapping, kept)

    Args:
        offset_table: Maps S288C positions to cumulative offsets
        locus_start_s288c: S288C start coordinate
        locus_end_s288c: S288C end coordinate
        strain_seq_len: Length of strain sequence

    Returns:
        Dict mapping strain_pos (0-based) to s288c_pos (1-based) or None
    """
    mapping = {}

    # Build forward mapping, keeping FIRST (earliest) S288C pos for each strain pos
    # This correctly handles deletions where multiple S288C positions would
    # otherwise map to the same strain position
    for s288c_pos in range(locus_start_s288c, locus_end_s288c + 1):
        strain_pos = s288c_to_strain_pos(s288c_pos, offset_table, locus_start_s288c)
        if 0 <= strain_pos < strain_seq_len:
            if strain_pos not in mapping:
                # Only set if not already mapped (keep first/earliest s288c_pos)
                mapping[strain_pos] = s288c_pos
            # else: this S288C position is deleted (collision with earlier position)

    # For positions not in mapping (insertions in strain), mark as None
    for strain_pos in range(strain_seq_len):
        if strain_pos not in mapping:
            mapping[strain_pos] = None

    return mapping


def parse_cigar_and_extract_variants(
    alignment,
    bg_seq: str,
    src_seq: str,
    background_strain: str,
    source_strain: str,
    bg_reverse_mapping: Dict[int, Optional[int]]
) -> List[AlignmentVariant]:
    """Parse alignment CIGAR and extract variants.

    Args:
        alignment: mappy alignment object
        bg_seq: Background (reference) sequence
        src_seq: Source (query) sequence
        background_strain: Name of background strain
        source_strain: Name of source strain
        bg_reverse_mapping: Maps bg_pos to s288c_pos

    Returns:
        List of AlignmentVariant objects
    """
    variants = []

    # Alignment positions
    ref_pos = alignment.r_st  # Reference start (0-based)
    query_pos = alignment.q_st  # Query start (0-based)

    # Parse CIGAR
    cigar = alignment.cigar  # List of (length, op) tuples

    for length, op in cigar:
        if op == 0:  # M (match/mismatch)
            # Check for mismatches
            for i in range(length):
                ref_base = bg_seq[ref_pos + i] if ref_pos + i < len(bg_seq) else 'N'
                query_base = src_seq[query_pos + i] if query_pos + i < len(src_seq) else 'N'

                if ref_base.upper() != query_base.upper():
                    s288c_pos = bg_reverse_mapping.get(ref_pos + i)
                    variants.append(AlignmentVariant(
                        background_strain=background_strain,
                        source_strain=source_strain,
                        bg_seq_pos=ref_pos + i,
                        bg_allele=ref_base.upper(),
                        src_allele=query_base.upper(),
                        s288c_pos=s288c_pos,
                        variant_type='SNP',
                        alignment_quality=alignment.mapq
                    ))
            ref_pos += length
            query_pos += length

        elif op == 1:  # I (insertion in query)
            # Source has extra bases relative to background
            ins_seq = src_seq[query_pos:query_pos + length].upper()
            s288c_pos = bg_reverse_mapping.get(ref_pos)

            variants.append(AlignmentVariant(
                background_strain=background_strain,
                source_strain=source_strain,
                bg_seq_pos=ref_pos,
                bg_allele='',
                src_allele=ins_seq,
                s288c_pos=s288c_pos,
                variant_type='INS',
                alignment_quality=alignment.mapq
            ))
            query_pos += length

        elif op == 2:  # D (deletion in query)
            # Background has extra bases relative to source
            del_seq = bg_seq[ref_pos:ref_pos + length].upper()
            s288c_pos = bg_reverse_mapping.get(ref_pos)

            variants.append(AlignmentVariant(
                background_strain=background_strain,
                source_strain=source_strain,
                bg_seq_pos=ref_pos,
                bg_allele=del_seq,
                src_allele='',
                s288c_pos=s288c_pos,
                variant_type='DEL',
                alignment_quality=alignment.mapq
            ))
            ref_pos += length

        elif op == 4:  # S (soft clip)
            query_pos += length

        elif op == 5:  # H (hard clip)
            pass  # No position change

        else:
            # Other operations (N, etc.)
            ref_pos += length
            query_pos += length

    return variants


def align_and_extract_variants(
    bg_seq: str,
    src_seq: str,
    background_strain: str,
    source_strain: str,
    bg_offset_table: Dict[int, int],
    locus_start_s288c: int,
    locus_end_s288c: int,
    min_mapq: int = 20
) -> List[AlignmentVariant]:
    """Extract variants by direct pairwise sequence alignment.

    Main function for alignment-based variant extraction. Aligns source
    sequence to background using minimap2, extracts variants from the
    alignment, and maps positions back to S288C coordinates.

    Args:
        bg_seq: Background strain sequence (will be used as alignment reference)
        src_seq: Source strain sequence (will be aligned to background)
        background_strain: Name of background strain
        source_strain: Name of source strain
        bg_offset_table: Offset table for background strain (S288C pos -> offset)
        locus_start_s288c: S288C start coordinate of the locus
        locus_end_s288c: S288C end coordinate of the locus
        min_mapq: Minimum mapping quality threshold

    Returns:
        List of AlignmentVariant objects with both bg_seq_pos and s288c_pos

    Example:
        >>> bg_offset = {1000: -5, 2000: -8}  # Background has deletions
        >>> variants = align_and_extract_variants(
        ...     bg_seq="ATGC...",
        ...     src_seq="ATCC...",  # SNP at some position
        ...     background_strain="strainS",
        ...     source_strain="strainL",
        ...     bg_offset_table=bg_offset,
        ...     locus_start_s288c=1,
        ...     locus_end_s288c=50000
        ... )
        >>> for v in variants:
        ...     print(f"{v.s288c_pos}: {v.bg_allele}->{v.src_allele}")
    """
    # Build reverse mapping for background
    bg_reverse_mapping = build_reverse_mapping(
        bg_offset_table,
        locus_start_s288c,
        locus_end_s288c,
        len(bg_seq)
    )

    # Create mappy aligner with background as reference
    # Use asm5 preset for closely related sequences (~0.1% divergence)
    aligner = mp.Aligner(seq=bg_seq, preset="asm5")

    if not aligner:
        raise RuntimeError("Failed to create mappy aligner")

    # Align source to background
    alignments = list(aligner.map(src_seq))

    if not alignments:
        # No alignment found - sequences might be too different
        # Fall back to simple position-by-position comparison for same-length seqs
        if len(bg_seq) == len(src_seq):
            return _fallback_snp_extraction(
                bg_seq, src_seq,
                background_strain, source_strain,
                bg_reverse_mapping
            )
        return []

    # Get best alignment (highest mapq)
    best_aln = max(alignments, key=lambda a: (a.mapq, a.mlen))

    if best_aln.mapq < min_mapq:
        return []

    # Extract variants from alignment
    return parse_cigar_and_extract_variants(
        best_aln,
        bg_seq,
        src_seq,
        background_strain,
        source_strain,
        bg_reverse_mapping
    )


def _fallback_snp_extraction(
    bg_seq: str,
    src_seq: str,
    background_strain: str,
    source_strain: str,
    bg_reverse_mapping: Dict[int, Optional[int]]
) -> List[AlignmentVariant]:
    """Fallback SNP extraction for same-length sequences.

    Used when alignment fails but sequences are same length.
    Simply compares position-by-position.
    """
    variants = []
    min_len = min(len(bg_seq), len(src_seq))

    for i in range(min_len):
        bg_base = bg_seq[i].upper()
        src_base = src_seq[i].upper()

        if bg_base != src_base and bg_base in 'ACGT' and src_base in 'ACGT':
            variants.append(AlignmentVariant(
                background_strain=background_strain,
                source_strain=source_strain,
                bg_seq_pos=i,
                bg_allele=bg_base,
                src_allele=src_base,
                s288c_pos=bg_reverse_mapping.get(i),
                variant_type='SNP',
                alignment_quality=60.0
            ))

    return variants


def extract_bidirectional_variants(
    seq_L: str,
    seq_S: str,
    strain_L: str,
    strain_S: str,
    offset_table_L: Dict[int, int],
    offset_table_S: Dict[int, int],
    locus_start_s288c: int,
    locus_end_s288c: int
) -> Tuple[List[AlignmentVariant], List[AlignmentVariant]]:
    """Extract variants in both directions between two strains.

    Returns two lists:
    - L_into_S: Variants to insert L's alleles into S background
    - S_into_L: Variants to insert S's alleles into L background

    Args:
        seq_L: Strain L's locus sequence
        seq_S: Strain S's locus sequence
        strain_L: Strain L name
        strain_S: Strain S name
        offset_table_L: Offset table for strain L
        offset_table_S: Offset table for strain S
        locus_start_s288c: S288C start coordinate
        locus_end_s288c: S288C end coordinate

    Returns:
        Tuple of (L_into_S_variants, S_into_L_variants)
    """
    # Direction 1: Insert L's variants into S (S is background)
    L_into_S = align_and_extract_variants(
        bg_seq=seq_S,
        src_seq=seq_L,
        background_strain=strain_S,
        source_strain=strain_L,
        bg_offset_table=offset_table_S,
        locus_start_s288c=locus_start_s288c,
        locus_end_s288c=locus_end_s288c
    )

    # Direction 2: Insert S's variants into L (L is background)
    S_into_L = align_and_extract_variants(
        bg_seq=seq_L,
        src_seq=seq_S,
        background_strain=strain_L,
        source_strain=strain_S,
        bg_offset_table=offset_table_L,
        locus_start_s288c=locus_start_s288c,
        locus_end_s288c=locus_end_s288c
    )

    return L_into_S, S_into_L


def apply_alignment_variant(
    sequence: str,
    variant: AlignmentVariant,
    verify_allele: bool = True
) -> Optional[str]:
    """Apply an alignment variant to a sequence.

    Args:
        sequence: Background sequence to modify
        variant: AlignmentVariant to apply
        verify_allele: If True, verify bg_allele matches before applying

    Returns:
        Modified sequence, or None if application failed
    """
    pos = variant.bg_seq_pos
    bg_allele = variant.bg_allele
    src_allele = variant.src_allele

    # Check bounds
    if pos < 0 or pos > len(sequence):
        return None

    if bg_allele and pos + len(bg_allele) > len(sequence):
        return None

    # Verify allele if requested
    if verify_allele and bg_allele:
        actual = sequence[pos:pos + len(bg_allele)]
        if actual.upper() != bg_allele.upper():
            # Try nearby positions
            for offset in range(-5, 6):
                check_pos = pos + offset
                if check_pos < 0 or check_pos + len(bg_allele) > len(sequence):
                    continue
                if sequence[check_pos:check_pos + len(bg_allele)].upper() == bg_allele.upper():
                    pos = check_pos
                    break
            else:
                return None

    return sequence[:pos] + src_allele + sequence[pos + len(bg_allele):]


# ============================================================================
# Tests
# ============================================================================

if __name__ == "__main__":
    print("Testing alignment_variants.py...")

    # Test 1: Basic SNP extraction
    print("\n1. Testing basic SNP extraction...")
    bg_seq = "ATGCATGCATGC"
    src_seq = "ATGCTTGCATGC"  # SNP at position 4: A->T

    # Empty offset table (no indels)
    offset_table = {}

    variants = align_and_extract_variants(
        bg_seq=bg_seq,
        src_seq=src_seq,
        background_strain="BG",
        source_strain="SRC",
        bg_offset_table=offset_table,
        locus_start_s288c=1,
        locus_end_s288c=12
    )

    print(f"  Found {len(variants)} variants")
    for v in variants:
        print(f"    pos={v.bg_seq_pos}, s288c={v.s288c_pos}: {v.bg_allele}->{v.src_allele}")

    assert len(variants) == 1, f"Expected 1 variant, got {len(variants)}"
    assert variants[0].bg_allele == 'A', f"Expected bg_allele='A', got '{variants[0].bg_allele}'"
    assert variants[0].src_allele == 'T', f"Expected src_allele='T', got '{variants[0].src_allele}'"
    print("  ✓ PASS")

    # Test 2: Reverse offset lookup
    print("\n2. Testing reverse offset lookup...")
    offset_table = {100: -5, 200: -10}  # Deletions at positions 100 and 200

    # Position 50 in strain should map to S288C 50 (no offset before it)
    s288c_pos = reverse_offset_lookup(50, offset_table, 1, 300)
    assert s288c_pos == 51, f"Expected 51, got {s288c_pos}"

    # Position 110 in strain should map to S288C ~115 (accounting for -5 offset)
    # strain_pos = s288c_pos - locus_start + offset
    # 110 = s288c_pos - 1 + (-5) => s288c_pos = 116
    s288c_pos = reverse_offset_lookup(110, offset_table, 1, 300)
    print(f"  strain_pos=110 -> s288c_pos={s288c_pos}")
    print("  ✓ PASS")

    # Test 3: apply_alignment_variant
    print("\n3. Testing apply_alignment_variant...")
    variant = AlignmentVariant(
        background_strain="BG",
        source_strain="SRC",
        bg_seq_pos=4,
        bg_allele="AT",
        src_allele="CCC",
        s288c_pos=5,
        variant_type="COMPLEX"
    )

    result = apply_alignment_variant("ATGCATGCATGC", variant)
    expected = "ATGCCCCGCATGC"
    assert result == expected, f"Expected '{expected}', got '{result}'"
    print(f"  Applied variant: '{result}'")
    print("  ✓ PASS")

    # Test 4: Build reverse mapping
    print("\n4. Testing build_reverse_mapping...")
    offset_table = {10: -2}  # 2bp deletion at position 10
    mapping = build_reverse_mapping(offset_table, 1, 20, 18)  # 20bp locus, 18bp after deletion

    # Position 0 should map to S288C 1
    assert mapping[0] == 1, f"Expected mapping[0]=1, got {mapping[0]}"
    # Position 8 should map to S288C 9 (before deletion)
    assert mapping[8] == 9, f"Expected mapping[8]=9, got {mapping[8]}"
    print(f"  Sample mappings: 0->{mapping[0]}, 8->{mapping[8]}")
    print("  ✓ PASS")

    print("\n✅ All alignment_variants.py tests passed!")
