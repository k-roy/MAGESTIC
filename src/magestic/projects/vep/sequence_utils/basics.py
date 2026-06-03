"""
Basic sequence operations.

Functions:
- reverse_complement: Get reverse complement of DNA sequence
- pad_sequence: Pad sequence symmetrically to target length
- gc_content: Calculate GC content as fraction

Author: Kevin R. Roy
"""


def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of DNA sequence.

    Handles: A, T, G, C, N (returns N for N)

    Args:
        seq: DNA sequence

    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def pad_sequence(seq: str, target_length: int, pad_char: str = 'N') -> str:
    """
    Pad sequence symmetrically to target length.

    If seq is longer than target, returns seq unchanged.

    Args:
        seq: Input sequence
        target_length: Desired length
        pad_char: Character to use for padding (default 'N')

    Returns:
        Padded sequence
    """
    if len(seq) >= target_length:
        return seq
    pad_total = target_length - len(seq)
    pad_left = pad_total // 2
    pad_right = pad_total - pad_left
    return pad_char * pad_left + seq + pad_char * pad_right


def gc_content(seq: str) -> float:
    """
    Calculate GC content as fraction.

    Args:
        seq: DNA sequence

    Returns:
        GC content as float between 0 and 1
    """
    seq = seq.upper()
    gc = sum(1 for b in seq if b in 'GC')
    total = sum(1 for b in seq if b in 'ATGC')
    return gc / total if total > 0 else 0.0


if __name__ == "__main__":
    # Test reverse complement
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AANN") == "NNTT"
    assert reverse_complement("atgc") == "gcat"  # lowercase
    print("reverse_complement: PASS")

    # Test padding
    assert pad_sequence("ATCG", 8) == "NNATCGNN"
    assert pad_sequence("ATCG", 4) == "ATCG"  # no padding needed
    assert pad_sequence("ATCG", 3) == "ATCG"  # longer than target
    print("pad_sequence: PASS")

    # Test GC content
    assert gc_content("ATGC") == 0.5
    assert gc_content("GGGG") == 1.0
    assert gc_content("AAAA") == 0.0
    assert gc_content("ATGCN") == 0.5  # N ignored
    print("gc_content: PASS")

    print("\nAll basics.py tests passed!")
