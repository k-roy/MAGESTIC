"""Extract molecular barcodes with declared anchors and explicit ambiguity.

The configured cassette determines the bc1 flanks. A declared suffix is
mandatory evidence, not a preferred match followed by an unchecked slice.
Unanchored positional extraction is a separately opted-in legacy mode.
"""
from __future__ import annotations

from decimal import Decimal, InvalidOperation
from numbers import Integral
import re
import pandas as pd


def _rc(seq: str) -> str:
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def _integer(value, label, *, minimum=0):
    if isinstance(value, bool):
        raise ValueError(f'{label} must be an integer >= {minimum}')
    try:
        number = Decimal(str(value))
        if not number.is_finite() or number != number.to_integral_value() or number < minimum:
            raise ValueError
        return int(number)
    except (ValueError, TypeError, InvalidOperation):
        raise ValueError(f'{label} must be an integer >= {minimum}') from None


def _anchor(value, label):
    if not isinstance(value, str) or not re.fullmatch('[ACGT]+', value.upper()):
        raise ValueError(f'{label} must be a nonempty ACGT sequence')
    return value.upper()


def _starts(sequence, anchor):
    position = sequence.find(anchor)
    while position >= 0:
        yield position + len(anchor)
        position = sequence.find(anchor, position + 1)


def _candidates(sequence, prefix, suffix, sizes):
    found = set()
    for start in _starts(sequence, prefix):
        for size in sizes:
            end = start + size
            if end > len(sequence):
                continue
            if suffix is not None and sequence[end:end + len(suffix)] != suffix:
                continue
            barcode = sequence[start:end]
            if re.fullmatch('[ACGT]+', barcode):
                found.add(barcode)
    return found


def extract_bc1_from_seq(
    seq: str | float | None,
    fwd_bc1_prefix: str = 'CAGGTCATGCTC',
    bc1_length: int = 20,
    fallback_prefix_length: int | None = 16,
    fwd_bc1_suffix: str | None = None,
    length_tolerance: int = 1,
    allow_positional_fallback: bool = False,
) -> str | None:
    """Return a unique ACGT barcode, canonicalized to cassette orientation.

    With a suffix, both flanks must match exactly at a permitted length, in
    either orientation. Every anchor occurrence is considered; disagreeing
    candidates make identity uncallable. Without a suffix only exact expected
    length after the declared forward prefix is supported.

    Legacy positional mode requires an explicit True flag and is used only
    when neither molecular anchor occurs. Short or ambiguous slices never call.
    This legacy option is not appropriate evidence of cassette identity.
    """
    if not isinstance(bc1_length, Integral) or isinstance(bc1_length, bool) or bc1_length < 1:
        raise ValueError('bc1_length must be a positive integer')
    if not isinstance(length_tolerance, Integral) or isinstance(length_tolerance, bool) or length_tolerance < 0:
        raise ValueError('length_tolerance must be a nonnegative integer')
    if type(allow_positional_fallback) is not bool:
        raise ValueError('allow_positional_fallback must be boolean')
    prefix = _anchor(fwd_bc1_prefix, 'fwd_bc1_prefix')
    suffix = None if fwd_bc1_suffix is None else _anchor(fwd_bc1_suffix, 'fwd_bc1_suffix')
    if allow_positional_fallback:
        offset = _integer(fallback_prefix_length, 'fallback_prefix_length')
    if seq is None or (not isinstance(seq, str) and pd.isna(seq)):
        return None
    if not isinstance(seq, str):
        raise ValueError('read sequence must be a string or missing')
    sequence = seq.upper()
    sizes = range(max(1, bc1_length-length_tolerance), bc1_length+length_tolerance+1) if suffix else (bc1_length,)
    found = _candidates(sequence, prefix, suffix, sizes)
    if suffix:
        found.update(_rc(barcode) for barcode in _candidates(sequence, _rc(suffix), _rc(prefix), sizes))
    if len(found) == 1:
        return next(iter(found))
    if found:
        return None
    if allow_positional_fallback and prefix not in sequence and (suffix is None or _rc(suffix) not in sequence):
        barcode = sequence[offset:offset+bc1_length]
        if len(barcode) == bc1_length and re.fullmatch('[ACGT]+', barcode):
            return barcode
    return None


def extract_redi_bc_vectorized(seq_series: pd.Series, prefix_series: pd.Series,
                               bc_length_series: pd.Series) -> pd.Series:
    """One unique anchored REDI barcode per row with aligned metadata.

    Integral spreadsheet lengths are accepted, fractional lengths refused.
    Empty anchor plus zero length explicitly disables the assay. Multiple
    distinct candidates, unknown bases and incomplete reads remain uncalled.
    """
    if not seq_series.index.equals(prefix_series.index) or not seq_series.index.equals(bc_length_series.index):
        raise ValueError('sequence, prefix and length indices must align exactly')
    output = []
    for seq, prefix, length in zip(seq_series, prefix_series, bc_length_series):
        size = _integer(length, 'REDI barcode length')
        empty = pd.isna(prefix) or prefix == ''
        if empty and size == 0:
            output.append(None); continue
        if empty or size == 0:
            raise ValueError('enabled REDI extraction requires an anchor and positive length')
        anchor = _anchor(prefix, 'REDI barcode prefix')
        if seq is None or (not isinstance(seq, str) and pd.isna(seq)):
            output.append(None); continue
        if not isinstance(seq, str):
            raise ValueError('read sequence must be a string or missing')
        candidates = _candidates(seq.upper(), anchor, None, (size,))
        output.append(next(iter(candidates)) if len(candidates) == 1 else None)
    return pd.Series(output, index=seq_series.index, dtype=object)
