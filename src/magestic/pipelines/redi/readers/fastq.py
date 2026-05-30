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


_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _rc(seq: str) -> str:
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def extract_bc1_from_seq(
    seq: str | float | None,
    fwd_bc1_prefix: str = "CAGGTCATGCTC",
    bc1_length: int = 20,
    fallback_prefix_length: int = 16,
    fwd_bc1_suffix: str | None = None,
    length_tolerance: int = 1,
) -> str | None:
    """Extract a bc1 from a single trimmed merged read.

    Matches the `bc1_from_screen.snakemake.02_parse_bc1_reads` convention:
    forward orientation searches for `fwd_bc1_prefix`, then tries lengths in
    ``[bc1_length - length_tolerance, bc1_length + length_tolerance]``,
    preferring a length whose post-bc1 region equals `fwd_bc1_suffix` exactly;
    falls back to the expected length without suffix validation. Reverse
    orientation searches against ``rc(suffix)`` as the prefix-equivalent and
    returns the bc1 reverse-complemented back to forward orientation
    (skipped when `fwd_bc1_suffix` is None — there's no RC anchor to use).
    Legacy positional fallback when neither anchor is found (KR1965/KR1967
    reads with a 20 bp leading trim).
    """
    if seq is None or (isinstance(seq, float) and np.isnan(seq)):
        return None

    s = seq.upper()
    prefix = fwd_bc1_prefix.upper()
    suffix = fwd_bc1_suffix.upper() if fwd_bc1_suffix else None
    min_len = max(1, bc1_length - length_tolerance)
    max_len = bc1_length + length_tolerance

    # Forward orientation
    pidx = s.find(prefix)
    if pidx >= 0:
        start = pidx + len(prefix)
        if suffix:
            for try_len in range(min_len, max_len + 1):
                end = start + try_len
                if end + len(suffix) <= len(s) and s[end:end + len(suffix)] == suffix:
                    return s[start:end]
        end = start + bc1_length
        if end <= len(s):
            return s[start:end]

    # Reverse orientation (requires a suffix to derive the RC anchor)
    if suffix:
        rc_prefix = _rc(suffix)
        rc_suffix = _rc(prefix)
        pidx = s.find(rc_prefix)
        if pidx >= 0:
            start = pidx + len(rc_prefix)
            for try_len in range(min_len, max_len + 1):
                end = start + try_len
                if end + len(rc_suffix) <= len(s) and s[end:end + len(rc_suffix)] == rc_suffix:
                    return _rc(s[start:end])
            end = start + bc1_length
            if end <= len(s):
                return _rc(s[start:end])

    return s[fallback_prefix_length : fallback_prefix_length + bc1_length]


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
