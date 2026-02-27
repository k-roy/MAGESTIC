# Sequence Utilities Module - Agent Notes

**Date:** 2026-02-10
**Purpose:** Common sequence operations for VEP pipeline

---

## Module Responsibility

Shared sequence utilities used by:
- Locus preparation (reverse complement, FASTA I/O)
- ESM1-v (translation)
- All models (one-hot encoding)

---

## Files to Create

```
sequence_utils/
├── __init__.py           # ✅ Created
├── AGENT_NOTES.md        # This file
├── fasta_io.py           # Read/write FASTA files
├── translation.py        # DNA → protein
├── one_hot.py            # One-hot encoding
└── basics.py             # reverse_complement, etc.
```

---

## Source Scripts to Extract From

### Primary Source: `20_run_trio_esm1v.py`

**Location:** `/path/to/projects/QTL/variant_effect_prediction/scripts/20260205_in_silico_QTN_dissection/20_run_trio_esm1v.py`

**Key Functions:**
- `reverse_complement()` → `basics.py`
- `translate()` → `translation.py`
- `CODON_TABLE` → `translation.py`

### Secondary Source: `19_run_trio_dissection_models.py`

**Key Functions:**
- `one_hot_encode()` → `one_hot.py`
- `pad_sequence()` → `basics.py`

---

## Function Specifications

### `basics.py`

```python
def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of DNA sequence.

    Handles: A, T, G, C, N (returns N for N)
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base.upper(), 'N') for base in reversed(seq))


def pad_sequence(seq: str, target_length: int, pad_char: str = 'N') -> str:
    """
    Pad sequence symmetrically to target length.

    If seq is longer than target, returns seq unchanged.
    """
    if len(seq) >= target_length:
        return seq
    pad_total = target_length - len(seq)
    pad_left = pad_total // 2
    pad_right = pad_total - pad_left
    return pad_char * pad_left + seq + pad_char * pad_right


def gc_content(seq: str) -> float:
    """Calculate GC content as fraction."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in 'GC')
    total = sum(1 for b in seq if b in 'ATGC')
    return gc / total if total > 0 else 0.0
```

### `translation.py`

```python
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
        (codon, amino_acid)
    """
    start = aa_position * 3
    codon = cds[start:start+3].upper()
    aa = CODON_TABLE.get(codon, 'X')
    return codon, aa
```

### `one_hot.py`

```python
import numpy as np


def one_hot_encode(seq: str) -> np.ndarray:
    """
    One-hot encode DNA sequence.

    Args:
        seq: DNA sequence (ATGCN)

    Returns:
        np.ndarray of shape (len(seq), 4) with order [A, C, G, T]
        N bases get 0.25 for all channels
    """
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    seq_upper = seq.upper()
    indices = np.array([mapping.get(base, 4) for base in seq_upper])

    one_hot = np.zeros((len(seq), 4), dtype=np.float32)
    valid_mask = indices < 4
    one_hot[valid_mask, indices[valid_mask]] = 1.0

    # N bases get 0.25 for all channels
    n_mask = indices == 4
    one_hot[n_mask] = 0.25

    return one_hot


def one_hot_decode(one_hot: np.ndarray) -> str:
    """
    Decode one-hot array back to sequence.

    Assumes order [A, C, G, T].
    """
    bases = ['A', 'C', 'G', 'T']
    indices = np.argmax(one_hot, axis=1)
    return ''.join(bases[i] for i in indices)
```

### `fasta_io.py`

```python
from pathlib import Path
from typing import Dict, Iterator, Tuple


def read_fasta(fasta_file: Path) -> Dict[str, str]:
    """
    Read FASTA file into dictionary.

    Returns:
        {header: sequence} for all entries
    """
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]  # First word after >
                current_seq = []
            else:
                current_seq.append(line)

        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def write_fasta(sequences: Dict[str, str], output_file: Path, line_width: int = 80):
    """
    Write sequences to FASTA file.

    Args:
        sequences: {header: sequence}
        output_file: Output path
        line_width: Characters per line (default 80)
    """
    with open(output_file, 'w') as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')


def iter_fasta(fasta_file: Path) -> Iterator[Tuple[str, str]]:
    """
    Iterate over FASTA entries (memory-efficient for large files).

    Yields:
        (header, sequence) tuples
    """
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    yield current_header, ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_header:
            yield current_header, ''.join(current_seq)
```

---

## Dependencies

- `numpy`: One-hot encoding

**Conda environment:** Any (pure Python + numpy)

---

## Testing

```python
if __name__ == "__main__":
    # Test reverse complement
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AANN") == "NNTT"
    print("reverse_complement: PASS")

    # Test translation
    assert translate("ATGGCGTAA") == "MA"  # ATG=M, GCG=A, TAA=stop
    assert translate("ATGTGA") == "M"      # TGA is stop
    print("translate: PASS")

    # Test one-hot
    oh = one_hot_encode("ACGT")
    assert oh.shape == (4, 4)
    assert oh[0, 0] == 1.0  # A
    assert oh[1, 1] == 1.0  # C
    print("one_hot_encode: PASS")

    # Test padding
    assert pad_sequence("ATCG", 8) == "NNATCGNN"
    print("pad_sequence: PASS")

    print("\nAll tests passed!")
```

---

*Last updated: 2026-02-10*
