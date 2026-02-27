"""
Yorzoi Model Module

Yorzoi is a Borzoi-like model for yeast RNA expression prediction.
It uses a tiled prediction approach for sequences longer than 4992 bp.

- Uses flash-attention (requires Ampere+ GPU: A100/H100)
- HuggingFace model: tom-ellis-lab/yorzoi
- Input: 4992 bp DNA sequences (tiled for longer)
- Output: Per-position RNA coverage predictions

Author: Kevin R. Roy
"""

from .runner import (
    YorzoiModel,
    run_yorzoi_batch,
    one_hot_encode,
    MODEL_INPUT_LENGTH,
    MODEL_OUTPUT_LENGTH,
)

__all__ = [
    "YorzoiModel",
    "run_yorzoi_batch",
    "one_hot_encode",
    "MODEL_INPUT_LENGTH",
    "MODEL_OUTPUT_LENGTH",
]
