"""
Shorkie Model Module

Shorkie is an 8-fold ensemble model for predicting RNA expression effects
in yeast. It uses a Borzoi-like architecture with 16kb input sequences.

- 8-fold ensemble model
- Input: 16,384 bp DNA sequences
- Output: Per-position RNA-seq coverage predictions
- Package: baskerville-yeast

Author: Kevin R. Roy
"""

from .runner import (
    ShorkieModel,
    run_shorkie_batch,
    one_hot_encode,
    SEQ_LENGTH,
    NUM_FOLDS,
)

__all__ = [
    "ShorkieModel",
    "run_shorkie_batch",
    "one_hot_encode",
    "SEQ_LENGTH",
    "NUM_FOLDS",
]
