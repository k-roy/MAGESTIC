"""Strict identity keys and barcode validity; no registration inferred from calls."""
from __future__ import annotations
from collections.abc import Sequence
from decimal import Decimal, InvalidOperation
import pandas as pd

POSITION_COLUMNS = ('MAT_alpha_array_1536_barcode', 'REDI_plate_color_384', 'REDI_row_384', 'REDI_col_384')


def count_integer(value, label='counts') -> int:
    """Parse exact nonnegative counts without routing text through float."""
    try:
        if isinstance(value, bool) or (isinstance(value, float) and abs(value) >= 2**53):
            raise ValueError
        number=Decimal(str(value))
        if not number.is_finite() or number<0 or number!=number.to_integral_value():
            raise ValueError
        return int(number)
    except (ValueError,TypeError,InvalidOperation):
        raise ValueError(f'{label}: counts must be exact finite nonnegative integers') from None


def keyed_frame(frame: pd.DataFrame, keys: Sequence[str], *, label: str, unique: bool = True) -> pd.DataFrame:
    """Validate supplied keys; normalize integral spreadsheet scalars only.

    Callers must register plates into one dated physical namespace using transfer
    records first. This cannot select between alternative worklist generations.
    """
    keys = list(keys)
    if not keys or len(set(keys)) != len(keys):
        raise ValueError(f'{label}: position keys must be nonempty and unique')
    if frame.columns.duplicated().any():
        raise ValueError(f'{label}: duplicate column names')
    missing = set(keys)-set(frame.columns)
    if missing:
        raise ValueError(f'{label}: missing position/key columns {sorted(missing)}')
    out = frame.copy()
    for key in keys:
        values = out[key].astype('string').str.strip()
        if values.isna().any() or values.eq('').any():
            raise ValueError(f'{label}: missing physical position/key {key}')
        out[key] = values.str.replace(r'^(\d+)\.0+$', r'\1', regex=True)
    if unique and out.duplicated(keys).any():
        raise ValueError(f'{label}: duplicate physical position/key; unique registration required')
    return out


def valid_barcode(series: pd.Series, length: int = 20) -> pd.Series:
    if type(length) is not int or length <= 0:
        raise ValueError('barcode length must be a positive integer')
    return series.astype('string').str.fullmatch(f'[ACGT]{{{length}}}', na=False).astype(bool)
