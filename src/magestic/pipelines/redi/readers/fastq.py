"""
magestic.pipelines.redi.readers.fastq — bc1 / REDI_bc extraction from collapsed reads.

Ports the bc1 / REDI_bc anchor-finding logic from `step_c`. Kept as pure
functions so step_c's pipeline orchestrator + unit tests can import them.

Inputs are the `*_R1_R2_counts.tsv` tables produced by step_a's
`collapsed_reads_to_counts_table.py`. Each row is (R1_seq, R2_seq, counts,
sample_name).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def extract_bc1_from_seq(
    seq: str | float | None,
    fwd_bc1_prefix: str = "CAGGTCATGCTC",
    bc1_length: int = 20,
    fallback_prefix_length: int = 16,
) -> str | None:
    """Extract a bc1 from a single trimmed merged read.

    Search for `fwd_bc1_prefix`; if found, take the next `bc1_length` bp.
    Otherwise fall back to a fixed slice (legacy step_c behavior for KR1965/
    KR1967 reads that present `TCGACAGGTCATGCTC` after 20 bp trimming).
    """
    if seq is None or (isinstance(seq, float) and np.isnan(seq)):
        return None
    idx = seq.find(fwd_bc1_prefix)
    if idx != -1:
        bc1_start = idx + len(fwd_bc1_prefix)
        return seq[bc1_start : bc1_start + bc1_length]
    return seq[fallback_prefix_length : fallback_prefix_length + bc1_length]


def extract_redi_bc_vectorized(
    seq_series: pd.Series,
    prefix_series: pd.Series,
    bc_length_series: pd.Series,
) -> pd.Series:
    """Vectorized REDI_bc extraction (step_c `extract_REDI_bc_vectorized`).

    `prefix_series` and `bc_length_series` come from the combined keyfile's
    `MAT_a_barcode_prefix` and `MAT_a_barcode_length` columns; they vary by
    inner primer set.
    """
    seq_arr = seq_series.fillna("").values
    prefix_arr = prefix_series.fillna("").values
    bc_length_arr = bc_length_series.fillna(0).astype(int).values
    prefix_length_arr = np.vectorize(len)(prefix_arr)

    idx_arr = np.array(
        [s.find(p) if p else -1 for s, p in zip(seq_arr, prefix_arr)]
    )
    bc_start_arr = idx_arr + prefix_length_arr
    bc_end_arr = bc_start_arr + bc_length_arr

    result = np.where(
        idx_arr != -1,
        [
            s[start:end] if start >= 0 and end <= len(s) else None
            for s, start, end in zip(seq_arr, bc_start_arr, bc_end_arr)
        ],
        None,
    )
    return pd.Series(result, index=seq_series.index)
