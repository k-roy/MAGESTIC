"""
DNA/RNA sequence utilities.

Refactored from the legacy bedgraph_computation.py module.
Removes Python 2 artifacts (imp.reload, circular imports).

Author: Kevin R. Roy
"""

from typing import Dict

# IUPAC complement map (includes ambiguous bases)
_COMPLEMENT: Dict[str, str] = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
    'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
    'D': 'H', 'H': 'D', 'V': 'B',
    'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
    'n': 'n', 'r': 'y', 'y': 'r', 's': 's',
    'w': 'w', 'k': 'm', 'm': 'k',
    ' ': ' ', '': '',
}

# IUPAC compatible base pairs for ambiguous base matching
_IUPAC_COMPATIBLE: Dict[str, str] = {
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}

# Standard genetic code (DNA codon → single-letter amino acid)
_CODON_TABLE: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def rev_comp(dna: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Handles IUPAC ambiguity codes and lowercase bases.

    Args:
        dna: DNA sequence string (IUPAC alphabet, upper or lower case)

    Returns:
        Reverse complement string

    Raises:
        KeyError: If an unrecognized base character is encountered
    """
    return ''.join(_COMPLEMENT[base] for base in reversed(dna))


def base_match(base1: str, base2: str, ambiguous_bases_allowed: bool = True) -> bool:
    """
    Check whether two IUPAC characters are compatible.

    Args:
        base1: First IUPAC base character
        base2: Second IUPAC base character
        ambiguous_bases_allowed: If True, IUPAC ambiguity codes are expanded
            and compared (e.g., R matches A or G)

    Returns:
        True if the bases are compatible, False otherwise
    """
    if base1 == base2:
        return True
    if ambiguous_bases_allowed:
        if base1 in _IUPAC_COMPATIBLE and base2 in _IUPAC_COMPATIBLE[base1]:
            return True
        if base2 in _IUPAC_COMPATIBLE and base1 in _IUPAC_COMPATIBLE[base2]:
            return True
    return False


def hamming_distance(
    string1: str,
    string2: str,
    ambiguous_bases_allowed: bool = True,
) -> int:
    """
    Count mismatches between two equal-length DNA strings.

    Args:
        string1: First DNA sequence
        string2: Second DNA sequence (must be same length as string1)
        ambiguous_bases_allowed: If True, IUPAC ambiguity codes are treated
            as compatible with any base they represent

    Returns:
        Number of mismatched positions
    """
    return sum(
        0 if base_match(b1, b2, ambiguous_bases_allowed) else 1
        for b1, b2 in zip(string1, string2)
    )


def edit_distance(v: str, w: str) -> tuple[int, int]:
    """
    Compute the best fitting alignment score between sequences v and w.

    Uses a dynamic programming algorithm that allows any 5' prefix of v
    to be unpenalized (fitting alignment — w must be fully consumed).
    Scoring: match = +1, mismatch penalty = 1, indel penalty = 1.

    Args:
        v: Longer reference sequence (length up to ~1000)
        w: Query sequence to fit within v (length up to ~100)

    Returns:
        Tuple of (best_score, start_index_in_v) where start_index_in_v
        is the 0-based position in v where the best alignment of w ends.
    """
    indel_penalty = 1
    num_rows = len(v) + 1
    num_columns = len(w) + 1
    import numpy as np
    path_array = np.zeros((num_rows, num_columns))
    backtracking_array = np.zeros((num_rows, num_columns))

    for i in range(1, num_rows):
        backtracking_array[i][0] = 1  # allow free 5' gaps in v

    for j in range(1, num_columns):
        path_array[0][j] = path_array[0][j - 1] - indel_penalty
        backtracking_array[0][j] = 2

    for i in range(1, num_rows):
        for j in range(1, num_columns):
            match_score = 1 if v[i - 1] == w[j - 1] else -1
            from_diag = path_array[i - 1][j - 1] + match_score
            from_up = path_array[i - 1][j] - indel_penalty
            from_left = path_array[i][j - 1] - indel_penalty
            best = max(from_diag, from_up, from_left)
            path_array[i][j] = best
            if best == from_diag:
                backtracking_array[i][j] = 0
            elif best == from_up:
                backtracking_array[i][j] = 1
            else:
                backtracking_array[i][j] = 2

    # Best score is the maximum value in the last column (w fully consumed)
    last_col = path_array[:, -1]
    best_row = int(np.argmax(last_col))
    best_score = int(last_col[best_row])
    return best_score, best_row


def translate_codon(codon: str) -> str:
    """
    Translate a DNA codon to a single-letter amino acid.

    Args:
        codon: Three-letter DNA codon (uppercase)

    Returns:
        Single-letter amino acid code, or '*' for stop codons

    Raises:
        KeyError: If the codon is not in the standard genetic code
    """
    return _CODON_TABLE[codon.upper()]


def translate_sequence(dna: str, codon_table: Dict[str, str] | None = None) -> str:
    """
    Translate a DNA sequence to a protein sequence.

    Args:
        dna: DNA sequence (length must be a multiple of 3)
        codon_table: Optional custom codon table. Defaults to standard code.

    Returns:
        Protein sequence string (stop codon represented as '*')
    """
    table = codon_table if codon_table is not None else _CODON_TABLE
    return ''.join(
        table.get(dna[i:i + 3].upper(), 'X')
        for i in range(0, len(dna) - 2, 3)
    )
