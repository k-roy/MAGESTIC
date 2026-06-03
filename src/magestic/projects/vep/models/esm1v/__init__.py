"""
ESM1-v Model Module

ESM1-v is a protein language model for scoring missense variant effects.
It uses masked language modeling to compute the effect of amino acid
substitutions on protein fitness.

- Model: esm1v_t33_650M_UR90S_1 (650M parameters)
- Input: Protein sequence + position + amino acid substitution
- Output: Delta log-likelihood (mutant - WT)
- Package: fair-esm

Author: Kevin R. Roy
"""

from .runner import (
    ESM1vModel,
    run_esm1v_batch,
    DELETERIOUS_THRESHOLD,
)

__all__ = [
    "ESM1vModel",
    "run_esm1v_batch",
    "DELETERIOUS_THRESHOLD",
]
