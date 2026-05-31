"""
Sequence manipulation utilities.

Provides common functions for DNA sequence operations used across projects.

Author: Kevin R. Roy
"""

import re
from typing import List, Tuple


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


# Alias for backward compatibility
rev_comp = reverse_complement


def parse_cigar(cigar_str: str) -> List[Tuple[int, str]]:
    """Parse CIGAR string into list of (length, operation) tuples.

    CIGAR operations:
    - M: alignment match (can be match or mismatch)
    - I: insertion to reference
    - D: deletion from reference
    - N: skipped region from reference
    - S: soft clipping (sequence present but not aligned)
    - H: hard clipping (sequence not present)
    - P: padding
    - =: sequence match
    - X: sequence mismatch
    """
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    return [(int(length), op) for length, op in pattern.findall(cigar_str)]


def get_alignment_length(cigar_str: str) -> int:
    """Calculate the reference alignment length from a CIGAR string.

    Counts operations that consume reference: M, D, N, =, X
    """
    parsed = parse_cigar(cigar_str)
    ref_consuming = {'M', 'D', 'N', '=', 'X'}
    return sum(length for length, op in parsed if op in ref_consuming)


def hamming_distance(seq1: str, seq2: str) -> int:
    """Count mismatches between two equal-length sequences.

    Raises ValueError if sequences have different lengths.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must be equal length: {len(seq1)} vs {len(seq2)}")
    return sum(a != b for a, b in zip(seq1.upper(), seq2.upper()))


def gc_content(seq: str) -> float:
    """Calculate GC content of a sequence (0.0 to 1.0)."""
    seq = seq.upper()
    gc = sum(1 for base in seq if base in 'GC')
    total = sum(1 for base in seq if base in 'ACGT')
    return gc / total if total > 0 else 0.0


def find_kmer_positions(sequence: str, kmer: str, allow_overlap: bool = True) -> List[int]:
    """Find all positions of a k-mer in a sequence.

    Args:
        sequence: DNA sequence to search
        kmer: K-mer to find
        allow_overlap: If True, overlapping matches are returned

    Returns:
        List of 0-based start positions
    """
    positions = []
    seq_upper = sequence.upper()
    kmer_upper = kmer.upper()
    start = 0

    while True:
        pos = seq_upper.find(kmer_upper, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1 if allow_overlap else pos + len(kmer)

    return positions


def translate_codon(codon: str) -> str:
    """Translate a DNA codon to amino acid (single letter).

    Returns '*' for stop codons, 'X' for invalid codons.
    """
    codon_table = {
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
    return codon_table.get(codon.upper(), 'X')
