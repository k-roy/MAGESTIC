"""
One-hot encoding utilities for neural models.

Functions:
- one_hot_encode: Encode DNA sequence as one-hot array
- one_hot_decode: Decode one-hot array back to sequence

Author: Kevin R. Roy
"""

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

    Args:
        one_hot: Array of shape (seq_len, 4) with order [A, C, G, T]

    Returns:
        DNA sequence string
    """
    bases = ['A', 'C', 'G', 'T']
    indices = np.argmax(one_hot, axis=1)
    return ''.join(bases[i] for i in indices)


if __name__ == "__main__":
    # Test one-hot encoding
    oh = one_hot_encode("ACGT")
    assert oh.shape == (4, 4)
    assert oh[0, 0] == 1.0  # A
    assert oh[1, 1] == 1.0  # C
    assert oh[2, 2] == 1.0  # G
    assert oh[3, 3] == 1.0  # T
    print("one_hot_encode basic: PASS")

    # Test N handling
    oh_n = one_hot_encode("AN")
    assert oh_n[0, 0] == 1.0  # A
    assert np.allclose(oh_n[1], [0.25, 0.25, 0.25, 0.25])  # N
    print("one_hot_encode N handling: PASS")

    # Test decode
    decoded = one_hot_decode(oh)
    assert decoded == "ACGT"
    print("one_hot_decode: PASS")

    # Round trip
    original = "ATGCATGC"
    decoded = one_hot_decode(one_hot_encode(original))
    assert decoded == original
    print("round trip: PASS")

    print("\nAll one_hot.py tests passed!")
