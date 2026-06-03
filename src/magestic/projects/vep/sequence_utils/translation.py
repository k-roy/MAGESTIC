"""
DNA to protein translation utilities.

Functions:
- translate: Translate CDS to protein sequence
- translate_with_stops: Translate CDS including stop codons in output
- get_codon_at_position: Get codon and amino acid at a position
- find_aa_differences: Find AA changes between two proteins
- combine_nearby_aa_changes: Group AA changes within a distance

Classes:
- AAVariant: Single amino acid variant
- CombinedAAVariant: Group of nearby AA variants

Constants:
- CODON_TABLE: Standard genetic code

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple
from itertools import combinations

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(cds: str, remove_stop: bool = True) -> str:
    """
    Translate CDS to protein sequence.

    Args:
        cds: Coding DNA sequence (must be in-frame)
        remove_stop: If True, stop at first stop codon and exclude it from output.
                     If False, include '*' for stop codon and stop translation.

    Returns:
        Protein sequence. If remove_stop=True, excludes stop codon.
        If remove_stop=False, includes '*' for the stop codon.
    """
    protein = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            if remove_stop:
                break
            # Include stop codon when remove_stop=False
            protein.append(aa)
            continue
        protein.append(aa)
    return ''.join(protein)


def get_codon_at_position(cds: str, aa_position: int) -> tuple:
    """
    Get codon and amino acid at a position.

    Args:
        cds: Coding sequence
        aa_position: 0-based amino acid position

    Returns:
        (codon, amino_acid) tuple
    """
    start = aa_position * 3
    if start + 3 > len(cds):
        return None, None
    codon = cds[start:start+3].upper()
    aa = CODON_TABLE.get(codon, 'X')
    return codon, aa


def translate_with_stops(cds: str) -> str:
    """
    Translate CDS including stop codons in output.

    Unlike translate(), this function includes '*' for stop codons and
    continues translating after stop codons. Useful for analyzing frameshift
    variants that may introduce or remove stop codons.

    Args:
        cds: Coding DNA sequence (must be in-frame)

    Returns:
        Protein sequence with '*' for stop codons
    """
    protein = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)


# ============================================================================
# AA Variant Detection and Combining
# ============================================================================

@dataclass
class AAVariant:
    """Single amino acid variant.

    Represents a difference between two protein sequences at a single position.

    Attributes:
        position: 1-based AA position
        ref_aa: Reference (background) amino acid
        alt_aa: Alternate (source) amino acid
    """
    position: int  # 1-based
    ref_aa: str
    alt_aa: str

    def __str__(self) -> str:
        """Return string representation like 'V41I'."""
        return f"{self.ref_aa}{self.position}{self.alt_aa}"

    @property
    def is_missense(self) -> bool:
        """Check if this is a missense variant (AA substitution)."""
        return self.ref_aa != '*' and self.alt_aa != '*'

    @property
    def is_nonsense(self) -> bool:
        """Check if this introduces a stop codon."""
        return self.ref_aa != '*' and self.alt_aa == '*'

    @property
    def is_stop_loss(self) -> bool:
        """Check if this removes a stop codon."""
        return self.ref_aa == '*' and self.alt_aa != '*'


@dataclass
class CombinedAAVariant:
    """Group of nearby AA variants to be analyzed together.

    When multiple AA changes occur within a short distance (e.g., ≤3 positions),
    they may represent a single functional unit and should be tested together.

    Attributes:
        variants: List of AAVariant objects in this group
        positions: List of 1-based positions
    """
    variants: List[AAVariant]

    @property
    def positions(self) -> List[int]:
        """Get sorted list of positions."""
        return sorted(v.position for v in self.variants)

    @property
    def n_changes(self) -> int:
        """Number of AA changes in this group."""
        return len(self.variants)

    @property
    def ref_aas(self) -> str:
        """Reference AAs as a string (sorted by position)."""
        return ''.join(v.ref_aa for v in sorted(self.variants, key=lambda x: x.position))

    @property
    def alt_aas(self) -> str:
        """Alternate AAs as a string (sorted by position)."""
        return ''.join(v.alt_aa for v in sorted(self.variants, key=lambda x: x.position))

    def __str__(self) -> str:
        """Return string representation like 'V41I+A43T'."""
        sorted_vars = sorted(self.variants, key=lambda x: x.position)
        return '+'.join(str(v) for v in sorted_vars)

    @property
    def positions_str(self) -> str:
        """Positions as comma-separated string."""
        return ','.join(str(p) for p in self.positions)

    @property
    def variant_effect(self) -> str:
        """Classify the combined effect."""
        if any(v.is_nonsense for v in self.variants):
            return 'nonsense'
        elif any(v.is_stop_loss for v in self.variants):
            return 'stop_loss'
        else:
            return 'missense'


def find_aa_differences(
    ref_protein: str,
    alt_protein: str,
    include_length_diff: bool = False
) -> Tuple[List[AAVariant], Optional[int]]:
    """Find amino acid differences between two proteins.

    Compares two protein sequences position-by-position and returns all
    positions where they differ.

    Args:
        ref_protein: Reference (background) protein sequence
        alt_protein: Alternate (source) protein sequence
        include_length_diff: If True, also return length difference

    Returns:
        Tuple of (list of AAVariant objects, length_difference or None)
        Length difference is (alt_len - ref_len), positive if alt is longer.

    Example:
        >>> ref = "MAVIKEIVH"
        >>> alt = "MTVIKELVH"  # A2T, E7L
        >>> variants, _ = find_aa_differences(ref, alt)
        >>> len(variants)
        2
        >>> str(variants[0])
        'A2T'
    """
    variants = []
    min_len = min(len(ref_protein), len(alt_protein))

    for i in range(min_len):
        ref_aa = ref_protein[i]
        alt_aa = alt_protein[i]

        if ref_aa != alt_aa:
            variants.append(AAVariant(
                position=i + 1,  # 1-based
                ref_aa=ref_aa,
                alt_aa=alt_aa
            ))

    length_diff = None
    if include_length_diff:
        length_diff = len(alt_protein) - len(ref_protein)

    return variants, length_diff


def combine_nearby_aa_changes(
    variants: List[AAVariant],
    max_distance: int = 3,
    max_combo_size: int = 4
) -> List[CombinedAAVariant]:
    """Group AA changes that are within max_distance positions of each other.

    When multiple AA changes are close together, they may represent a single
    functional unit (e.g., two adjacent changes from a complex codon mutation).
    This function identifies such clusters and creates combined variants.

    Algorithm:
    1. Sort variants by position
    2. Find clusters where consecutive variants are within max_distance
    3. For each cluster, generate all valid combinations of size 2 to max_combo_size
    4. Validate that all consecutive pairs in each combination are within max_distance

    Args:
        variants: List of AAVariant objects
        max_distance: Maximum distance between consecutive variants to combine
                      (default: 3 positions)
        max_combo_size: Maximum number of variants to combine (default: 4)

    Returns:
        List of CombinedAAVariant objects (does NOT include single variants)

    Example:
        >>> variants = [
        ...     AAVariant(85, 'R', 'G'),
        ...     AAVariant(87, 'D', 'N'),  # 2 positions apart
        ...     AAVariant(100, 'A', 'V'),  # distant
        ... ]
        >>> combined = combine_nearby_aa_changes(variants, max_distance=3)
        >>> len(combined)
        1
        >>> str(combined[0])
        'R85G+D87N'
    """
    if len(variants) < 2:
        return []

    # Sort by position
    sorted_vars = sorted(variants, key=lambda v: v.position)

    combinations_list = []

    # Find clusters of nearby variants
    i = 0
    while i < len(sorted_vars):
        cluster = [sorted_vars[i]]
        j = i + 1

        # Extend cluster while variants are within max_distance
        while j < len(sorted_vars):
            if sorted_vars[j].position - cluster[-1].position <= max_distance:
                cluster.append(sorted_vars[j])
                j += 1
            else:
                break

        # Generate all valid combinations from cluster (size 2 to max_combo_size)
        if len(cluster) >= 2:
            max_size = min(len(cluster), max_combo_size)
            for size in range(2, max_size + 1):
                for combo in combinations(cluster, size):
                    combo_list = list(combo)
                    # Verify all consecutive pairs are within distance
                    valid = all(
                        combo_list[k + 1].position - combo_list[k].position <= max_distance
                        for k in range(len(combo_list) - 1)
                    )
                    if valid:
                        combinations_list.append(CombinedAAVariant(variants=combo_list))

        i = j if j > i + 1 else i + 1

    return combinations_list


def classify_frameshift_direction(
    wt_protein: str,
    mut_protein: str
) -> Tuple[str, str]:
    """Classify frameshift as truncating or lengthening based on protein length.

    Args:
        wt_protein: Wild-type protein sequence (may include '*')
        mut_protein: Mutant protein sequence (may include '*')

    Returns:
        Tuple of (variant_type, description)
        - 'truncating_frameshift' if mutant protein is shorter
        - 'lengthening_frameshift' if mutant protein is same length or longer
    """
    # Remove stop codons for length comparison
    wt_len = len(wt_protein.rstrip('*'))
    mut_len = len(mut_protein.rstrip('*'))

    if mut_len < wt_len:
        diff = wt_len - mut_len
        return "truncating_frameshift", f"Protein truncated by {diff} AA"
    else:
        diff = mut_len - wt_len
        if diff == 0:
            return "lengthening_frameshift", "Frameshift with same final length"
        else:
            return "lengthening_frameshift", f"Protein extended by {diff} AA"


def classify_start_loss(mut_cds: str) -> Tuple[str, str]:
    """Determine consequence of start codon loss.

    When the start codon (ATG → M) is mutated, translation would begin at
    the next ATG. This function determines if that next ATG is in-frame
    (position % 3 == 0) or out-of-frame.

    Args:
        mut_cds: Mutant CDS sequence (first codon should NOT be ATG)

    Returns:
        Tuple of (variant_type, description)
        - 'start_loss_inframe' if next ATG is in-frame
        - 'start_loss_outframe' if next ATG is out-of-frame
        - 'start_loss_no_atg' if no downstream ATG found
    """
    # Find next ATG in mutant CDS (skip first codon)
    for i in range(3, len(mut_cds) - 2, 1):
        if mut_cds[i:i + 3].upper() == 'ATG':
            if i % 3 == 0:
                aa_delay = i // 3
                return "start_loss_inframe", f"Translation delayed by {aa_delay} codons"
            else:
                return "start_loss_outframe", f"Next ATG at position {i} is out-of-frame"
    return "start_loss_no_atg", "No downstream ATG found"


if __name__ == "__main__":
    # Test translation
    assert translate("ATGGCGTAA") == "MA"  # ATG=M, GCG=A, TAA=stop
    assert translate("ATGTGA") == "M"      # TGA is stop
    assert translate("ATGGCG") == "MA"     # No stop codon
    assert translate("ATG") == "M"         # Single codon
    print("translate (remove_stop=True): PASS")

    # Test translate with remove_stop=False (includes stop codon as '*')
    assert translate("ATGGCGTAA", remove_stop=False) == "MA*"  # Should include '*'
    assert translate("ATGTGA", remove_stop=False) == "M*"       # TGA stop included
    assert translate("ATGGCG", remove_stop=False) == "MA"       # No stop, same result
    print("translate (remove_stop=False): PASS")

    # Test get_codon_at_position
    codon, aa = get_codon_at_position("ATGGCG", 0)
    assert codon == "ATG" and aa == "M"
    codon, aa = get_codon_at_position("ATGGCG", 1)
    assert codon == "GCG" and aa == "A"
    print("get_codon_at_position: PASS")

    # Test translate_with_stops
    result = translate_with_stops("ATGGCGTAAATG")  # M, A, *, M
    assert result == "MA*M", f"Expected 'MA*M', got '{result}'"
    print("translate_with_stops: PASS")

    # Test find_aa_differences
    ref = "MAVIKEIVH"
    alt = "MTVIKELVH"  # A2T at pos 2, I7L at pos 7 (1-based)
    variants, _ = find_aa_differences(ref, alt)
    assert len(variants) == 2, f"Expected 2 variants, got {len(variants)}"
    assert variants[0].position == 2, f"Expected position 2, got {variants[0].position}"
    assert str(variants[0]) == "A2T", f"Expected 'A2T', got '{variants[0]}'"
    assert variants[1].position == 7, f"Expected position 7, got {variants[1].position}"
    assert str(variants[1]) == "I7L", f"Expected 'I7L', got '{variants[1]}'"
    print("find_aa_differences: PASS")

    # Test with length difference
    ref = "MAVI"
    alt = "MAVIKEL"  # alt is 3 AA longer
    variants, length_diff = find_aa_differences(ref, alt, include_length_diff=True)
    assert len(variants) == 0, "Expected 0 variants (matching prefix)"
    assert length_diff == 3, f"Expected length_diff 3, got {length_diff}"
    print("find_aa_differences with length_diff: PASS")

    # Test combine_nearby_aa_changes
    test_variants = [
        AAVariant(85, 'R', 'G'),
        AAVariant(87, 'D', 'N'),  # 2 positions apart from 85
        AAVariant(100, 'A', 'V'),  # distant
    ]
    combined = combine_nearby_aa_changes(test_variants, max_distance=3)
    assert len(combined) == 1, f"Expected 1 combined variant, got {len(combined)}"
    assert str(combined[0]) == "R85G+D87N", f"Expected 'R85G+D87N', got '{combined[0]}'"
    assert combined[0].n_changes == 2
    assert combined[0].positions == [85, 87]
    print("combine_nearby_aa_changes: PASS")

    # Test with larger cluster
    cluster_variants = [
        AAVariant(10, 'A', 'V'),
        AAVariant(12, 'B', 'W'),  # 2 apart
        AAVariant(14, 'C', 'X'),  # 2 apart
        AAVariant(50, 'D', 'Y'),  # distant
    ]
    combined = combine_nearby_aa_changes(cluster_variants, max_distance=3)
    # Should get: (10,12), (12,14), (10,12,14) = 3 combinations
    assert len(combined) == 3, f"Expected 3 combinations, got {len(combined)}"
    print("combine_nearby_aa_changes with cluster: PASS")

    # Test classify_frameshift_direction
    vtype, desc = classify_frameshift_direction("MAVIKEIVHWVST", "MAVIK")  # truncated
    assert vtype == "truncating_frameshift", f"Expected truncating_frameshift, got {vtype}"

    vtype, desc = classify_frameshift_direction("MAVIK", "MAVIKEIVHWVST")  # extended
    assert vtype == "lengthening_frameshift", f"Expected lengthening_frameshift, got {vtype}"
    print("classify_frameshift_direction: PASS")

    # Test classify_start_loss
    # In-frame next ATG at position 9 (positions 0,1,2 are first codon, so pos 9 is codon 3)
    vtype, desc = classify_start_loss("TTTGCGATGAAA")  # ATG at position 6, which is 6%3=0 (in-frame)
    assert vtype == "start_loss_inframe", f"Expected start_loss_inframe, got {vtype}"

    # Out-of-frame next ATG
    vtype, desc = classify_start_loss("TTTATGAAA")  # ATG at position 3, 3%3=0 (in-frame!)
    # Wait, position 3 is in-frame. Let me fix this test.
    # Position 4 would be out of frame (4%3=1)
    vtype, desc = classify_start_loss("TTTTATGAAA")  # ATG at position 4, 4%3=1 (out-of-frame)
    assert vtype == "start_loss_outframe", f"Expected start_loss_outframe, got {vtype}"

    # No downstream ATG
    vtype, desc = classify_start_loss("TTTGCGTAA")
    assert vtype == "start_loss_no_atg", f"Expected start_loss_no_atg, got {vtype}"
    print("classify_start_loss: PASS")

    print("\n✅ All translation.py tests passed!")
