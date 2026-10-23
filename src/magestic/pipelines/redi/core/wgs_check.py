"""Compare calls at independently registered physical positions.

WGS is another observation, not automatically ground truth. Missing and invalid
calls are unassessed; padding and barcode membership cannot establish identity.
"""
from __future__ import annotations
from collections.abc import Sequence
import pandas as pd
from .identity import POSITION_COLUMNS, keyed_frame, valid_barcode


def normalize_wgs_bc1(series: pd.Series, length: int = 20, pad_char: str = 'C') -> pd.Series:
    """Compatibility name: preserve native values, missing and off-length calls.

    Legacy length/pad arguments no longer edit sequence evidence. Validity is a
    separate measurement made by compare_to_wgs.
    """
    return series.copy()


def compare_to_wgs(redi_df: pd.DataFrame, wgs_bc1_call_df: pd.DataFrame,
                   wgs_plate_key_df: pd.DataFrame, array_a_name: str = 'SHA345',
                   array_b_name: str = 'yT182', bc1_length: int = 20, *,
                   position_keys: Sequence[str] = POSITION_COLUMNS) -> pd.DataFrame:
    """One row per registered WGS sample, with nullable comparisons.

    Multiple WGS samples may measure one position, but each sample must map once
    and each REDI position must have one call row. Extra shared metadata never
    becomes a join key. Callable flags provide each comparison's denominator.
    Native barcode strings and invalid-length observations remain in output.
    """
    if not array_a_name or not array_b_name or array_a_name == array_b_name:
        raise ValueError('two distinct array names are required')
    keys = list(position_keys)
    if 'bc1_from_wgs' in keys or any(k.endswith('_bc1') for k in keys):
        raise ValueError('position keys cannot be observed barcodes')
    wgs = keyed_frame(wgs_bc1_call_df, ['sample_name'], label='WGS calls')
    plate = keyed_frame(wgs_plate_key_df, ['sample_name'], label='WGS plate map')
    plate = keyed_frame(plate, keys, label='WGS plate map', unique=False)
    redi = keyed_frame(redi_df, keys, label='REDI calls')
    a, b = f'{array_a_name}_bc1', f'{array_b_name}_bc1'
    for frame, forbidden, label in (
        (redi, {'sample_name', 'bc1_from_wgs'}, 'REDI calls'),
        (plate, {a, b, 'bc1_from_wgs'}, 'physical plate map'),
        (wgs, {a, b, *keys}, 'WGS calls'),
    ):
        collision = forbidden & set(frame.columns)
        if collision:
            raise ValueError(f'{label}: reserved observation/key columns {sorted(collision)}')
    if 'bc1_from_wgs' not in wgs or a not in redi or b not in redi:
        raise ValueError('WGS/REDI barcode observation columns are required')
    unknown = set(wgs.sample_name)-set(plate.sample_name)
    if unknown:
        raise ValueError(f'WGS calls lack physical registration for samples {sorted(unknown)}')
    if 'bc1_from_wgs' in plate or 'sample_name' in keys:
        raise ValueError('physical plate map cannot include WGS observation or sample-name position key')
    merged = plate.merge(wgs, on='sample_name', how='left', validate='one_to_one', suffixes=('', '_wgs'))
    merged = merged.merge(redi, on=keys, how='left', validate='many_to_one', suffixes=('_wgs_key', ''))
    merged['wgs_has_call'] = merged['bc1_from_wgs'].notna() & merged['bc1_from_wgs'].ne('')
    merged['wgs_call_valid'] = valid_barcode(merged['bc1_from_wgs'], bc1_length)
    valid, agreements = {}, {}
    for name, column in ((array_a_name,a),(array_b_name,b)):
        valid[name] = valid_barcode(merged[column], bc1_length)
        merged[f'{name}_call_valid'] = valid[name]
        callable_mask = merged['wgs_call_valid'] & valid[name]
        merged[f'wgs_and_{name}_callable'] = callable_mask
        agreement = (merged['bc1_from_wgs'].astype('string') == merged[column].astype('string')).astype('boolean')
        agreements[name] = agreement.where(callable_mask, pd.NA)
        merged[f'wgs_and_{name}_bc1_agreement'] = agreements[name]
    both = merged['wgs_call_valid'] & valid[array_a_name] & valid[array_b_name]
    merged[f'wgs_and_{array_a_name}_{array_b_name}_callable'] = both
    merged[f'wgs_and_{array_a_name}_{array_b_name}_agreement'] = (agreements[array_a_name] & agreements[array_b_name]).where(both,pd.NA)
    redi_agree = merged[a].astype('string').eq(merged[b].astype('string')).fillna(False)
    conflict = valid[array_a_name] & valid[array_b_name] & ~redi_agree
    selected = merged[b].astype('string').where(valid[array_b_name], merged[a].astype('string'))
    merged['REDI_called_bc1'] = selected.where((valid[array_a_name] | valid[array_b_name]) & ~conflict,pd.NA)
    merged['REDI_call_conflict'] = conflict
    return merged
