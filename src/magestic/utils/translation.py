"""
DNA to protein translation utilities.

Functions:
- translate: Translate CDS to protein sequence
- get_codon_at_position: Get codon and amino acid at a position

Constants:
- CODON_TABLE: Standard genetic code

Author: Kevin R. Roy
"""

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
        remove_stop: If True, stop at first stop codon (don't include it)

    Returns:
        Protein sequence (amino acids only, no stop codon)
    """
    protein = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            if remove_stop:
                break
            continue  # Skip stop codon
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


if __name__ == "__main__":
    # Test translation
    assert translate("ATGGCGTAA") == "MA"  # ATG=M, GCG=A, TAA=stop
    assert translate("ATGTGA") == "M"      # TGA is stop
    assert translate("ATGGCG") == "MA"     # No stop codon
    assert translate("ATG") == "M"         # Single codon
    print("translate: PASS")

    # Test get_codon_at_position
    codon, aa = get_codon_at_position("ATGGCG", 0)
    assert codon == "ATG" and aa == "M"
    codon, aa = get_codon_at_position("ATGGCG", 1)
    assert codon == "GCG" and aa == "A"
    print("get_codon_at_position: PASS")

    print("\nAll translation.py tests passed!")
