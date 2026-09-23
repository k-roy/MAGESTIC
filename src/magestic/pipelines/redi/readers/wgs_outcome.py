"""
magestic.pipelines.redi.readers.wgs_outcome — wrapper around Shengdi's WGS
outcome table.

Canonical input:
    /oak/.../collaborator/projects/WGS_MAGESTIC_QTL_bc_calling/output3/summaries/
    wgs_outcome_full_table.v0.tsv

A newer `v20250701.tsv` exists alongside; configs choose which version.
Schema documented in MAGESTIC-SCORE repo.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from .keyfiles import _read_key


def load_wgs_outcome(
    path: str | Path,
    bc1_col: str = "bc1_from_wgs",
    keep_cols: list[str] | None = None,
) -> pd.DataFrame:
    """Read Shengdi's WGS outcome table.

    Parameters
    ----------
    path : path to `wgs_outcome_full_table.v?.tsv`.
    bc1_col : column name to rename → `bc1` (legacy is `bc1_from_wgs`).
    keep_cols : optional whitelist; if provided, only these columns are kept
        (besides the renamed `bc1`).
    """
    df = _read_key(path)
    if bc1_col not in df:
        raise ValueError(f'WGS outcome missing barcode observation column {bc1_col!r}')
    if bc1_col != 'bc1' and 'bc1' in df:
        raise ValueError('barcode observation rename collision with existing bc1')
    if bc1_col != "bc1":
        df = df.rename(columns={bc1_col: "bc1"})
    if keep_cols is not None:
        missing = set(keep_cols) - set(df)
        if missing:
            raise ValueError(f'WGS outcome missing requested columns: {sorted(missing)}')
        cols = list(dict.fromkeys(keep_cols))
        if "bc1" not in cols:
            cols.append("bc1")
        df = df[cols].drop_duplicates().reset_index(drop=True)
    return df


def load_wgs_plate_key(
    wgs_plate_key_path: str | Path,
    plate_384_to_96_format_path: str | Path | None = None,
) -> pd.DataFrame:
    """Load independent sample→physical-position registration.

    Prefer explicit sample_name and four position columns. Legacy plate-level
    tables require an explicit sample_number in the 384→96 format. Row order
    never defines identity.
    """
    from ..core.identity import POSITION_COLUMNS, keyed_frame
    wgs_plate_key=_read_key(wgs_plate_key_path)
    if 'sample_name' in wgs_plate_key:
        out=keyed_frame(wgs_plate_key,['sample_name'],label='WGS plate registration')
        return keyed_frame(out,POSITION_COLUMNS,label='WGS plate registration',unique=False)
    if plate_384_to_96_format_path is None:
        raise ValueError('legacy WGS registration requires format with explicit sample_number')
    fmt=_read_key(plate_384_to_96_format_path)
    # Require the real sequencing-index registration before choosing a geometry
    # schema. Neither a 96-well nor 384-well row order establishes this identity.
    fmt=keyed_frame(fmt,['sample_number'],label='384-to-96 format',unique=False)
    if 'plate_96' in wgs_plate_key and 'plate_96' in fmt:
        plate_column='plate_96'
    elif 'plate_384' in wgs_plate_key and 'plate_384' in fmt:
        plate_column='plate_384'
    else:
        raise ValueError('WGS key and geometry need the same explicit plate_96 or plate_384 registration')
    fmt=keyed_frame(fmt,[plate_column,'row_384','col_384'],label='384-to-96 format')
    fmt=keyed_frame(fmt,[plate_column,'sample_number'],label='format sequencing sample registration')
    required=[plate_column,'sequencing_date','Nextera outer primer set']
    wgs_plate_key=keyed_frame(wgs_plate_key,required,label='legacy WGS registration',unique=False)
    join=[plate_column]+[c for c in ('row_384','col_384') if c in wgs_plate_key]
    collisions=(set(wgs_plate_key)&set(fmt))-set(join)
    if collisions:
        raise ValueError(f'ambiguous legacy registration/format columns: {sorted(collisions)}')
    wgs_plate_key=wgs_plate_key.merge(fmt,on=join,how='left',indicator=True)
    if wgs_plate_key['_merge'].ne('both').any():
        raise ValueError('plate format missing for explicitly registered WGS join coordinates')
    wgs_plate_key=wgs_plate_key.drop(columns='_merge')
    wgs_plate_key=keyed_frame(wgs_plate_key,['sample_number'],label='expanded WGS registration',unique=False)
    wgs_plate_key["sample_name"] = (
        wgs_plate_key["sequencing_date"].astype(str)
        + "_"
        + wgs_plate_key["Nextera outer primer set"]
        + "_"
        + wgs_plate_key["sample_number"]
    )
    return keyed_frame(wgs_plate_key,['sample_name'],label='expanded WGS registration')
