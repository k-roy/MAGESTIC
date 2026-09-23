"""
magestic.pipelines.redi.readers.keyfiles — keyfile loaders for REDI screens.

Ports `step_b_combine_keyfiles.py` and the REDI_bc key consolidation logic
from `step_c`.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def _read_key(path: str | Path) -> pd.DataFrame:
    import csv
    with Path(path).open(newline='') as handle:
        reader = csv.reader(handle, delimiter='\t', strict=True)
        try:
            header = next(reader, [])
            if not header or any(not c for c in header) or len(header) != len(set(header)):
                raise ValueError(f'{path}: missing or duplicate key column names')
            # Pandas can silently promote extra leading fields to an index or
            # fill missing fields. Validate shape with constant auxiliary memory
            # before its typed parse so physical rows cannot shift or disappear.
            for row in reader:
                if not row:
                    raise ValueError(f'{path}: blank physical row {reader.line_num}')
                if len(row) != len(header):
                    raise ValueError(f'{path}: row {reader.line_num} has {len(row)} fields; expected {len(header)} columns')
        except csv.Error as error:
            raise ValueError(f'{path}: malformed TSV row {reader.line_num}: {error}') from error
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def _no_overlap(left, right, keys, label):
    overlap = (set(left) & set(right)) - set(keys)
    if overlap:
        raise ValueError(f'{label}: metadata column overlap/collision {sorted(overlap)}')


def combine_keyfiles(
    sample_key_path: str | Path,
    gdna_plate_key_path: str | Path,
    redi_plate_barcode_key_path: str | Path,
) -> pd.DataFrame:
    """Register samples to explicitly declared assays and physical arrays.

    Legacy sample plates may have several declared PCR assays. Per-well inner
    or outer primer columns restrict that expansion when present. The resulting
    outer-primer/sample identity must be unique. Only an explicit ``empty``
    sample is a control; missing or unregistered identifiers never join as NaN.
    """
    from ..core.identity import keyed_frame
    sample = _read_key(sample_key_path)
    gdna = _read_key(gdna_plate_key_path)
    plate = _read_key(redi_plate_barcode_key_path)
    join = ['gDNA_plate_name'] + [c for c in ('inner_primers', 'outer_primers') if c in sample]
    sample = keyed_frame(sample, ['gDNA_plate_name','sample_number'] + join[1:], label='sample key')
    sample = keyed_frame(sample, ['REDI_plate_barcode'], label='sample plate key', unique=False)
    gdna = keyed_frame(gdna, ['gDNA_plate_name','outer_primers','inner_primers'], label='gDNA assay key')
    plate = keyed_frame(plate, ['REDI_plate_barcode'], label='REDI plate key')
    physical = [c for c in ('MAT_alpha_array_1536_barcode','MAT_a_array_1536_barcode') if c in plate]
    if physical:
        plate = keyed_frame(plate, physical, label='registered physical arrays', unique=False)
    if plate.REDI_plate_barcode.str.casefold().eq('empty').any():
        raise ValueError('empty control cannot be a registered REDI plate key')
    _no_overlap(sample, gdna, join, 'sample/gDNA join')
    combined = sample.merge(gdna, on=join, how='left', indicator=True)
    if combined['_merge'].ne('both').any():
        raise ValueError('sample key lacks registered gDNA assay metadata')
    combined = combined.drop(columns='_merge')
    combined = keyed_frame(combined, ['outer_primers','sample_number'], label='output sample identity')
    _no_overlap(combined, plate, ['REDI_plate_barcode'], 'physical plate join')
    empty = combined.REDI_plate_barcode.str.casefold().eq('empty')
    unknown = set(combined.loc[~empty,'REDI_plate_barcode']) - set(plate.REDI_plate_barcode)
    if unknown:
        raise ValueError(f'sample key has unregistered REDI plates {sorted(unknown)}')
    combined = combined.merge(plate, on='REDI_plate_barcode', how='left', validate='many_to_one')
    expected_status = combined.REDI_plate_barcode.str.casefold().map(
        lambda value: 'declared_empty' if value == 'empty' else 'registered')
    if ('registration_status' in combined and
            not combined.registration_status.eq(expected_status).all()):
        raise ValueError('source registration_status conflicts with declared physical/control registration')
    combined['registration_status'] = expected_status
    return combined


def load_redi_bc_keys(
    v5_key_path: str | Path | None = None,
    yt182_yt183_key_path: str | Path | None = None,
) -> pd.DataFrame:
    """Load + consolidate the REDI_bc keyfile(s).

    Ports the V5 + yT182/yT183 keyfile consolidation from step_c.
    The output schema is a uniform 6-column table:

        REDI_library, REDI_bc, REDI_plate_color_384, REDI_row_384,
        REDI_col_384, validated_by_yT100_mating (optional),
        confirmed_contaminated_by_yT100_mating (optional),
        clone_removed_in_RBYG_signature_set (optional)
    """
    parts: list[pd.DataFrame] = []

    if v5_key_path is not None:
        v5 = _read_key(v5_key_path)
        v5["REDI_library"] = "SHA345"
        v5["REDI_bc"] = v5["REDI_bc_seq"].str.upper()
        v5_minimal = [
            "REDI_library",
            "REDI_bc",
            "REDI_plate_color_384",
            "REDI_row_384",
            "REDI_col_384",
        ]
        for opt in ["validated_by_yT100_mating", "confirmed_contaminated_by_yT100_mating"]:
            if opt in v5.columns:
                v5_minimal.append(opt)
        parts.append(v5[v5_minimal])

    if yt182_yt183_key_path is not None:
        y = _read_key(yt182_yt183_key_path)
        y["REDI_bc"] = y["REDI_bc1_seq"].str.upper()
        y_minimal = ["REDI_library", "REDI_bc", "plate_384", "row_384", "col_384"]
        if "clone_removed_in_RBYG_signature_set" in y.columns:
            y_minimal.append("clone_removed_in_RBYG_signature_set")
        y = y[y_minimal]
        y = y.rename(
            columns={
                "plate_384": "REDI_plate_color_384",
                "row_384": "REDI_row_384",
                "col_384": "REDI_col_384",
            }
        )
        parts.append(y)

    if not parts:
        return pd.DataFrame()
    from ..core.identity import keyed_frame
    combined = keyed_frame(pd.concat(parts, ignore_index=True), ['REDI_library','REDI_bc'],
                           label='REDI barcode registration')
    combined = keyed_frame(combined, ['REDI_plate_color_384','REDI_row_384','REDI_col_384'],
                           label='REDI barcode physical position', unique=False)
    if not combined.REDI_bc.str.fullmatch('[ACGT]+').all():
        raise ValueError('REDI barcode registration requires ACGT sequences')
    if not combined.REDI_row_384.isin(list('ABCDEFGHIJKLMNOP')).all():
        raise ValueError('REDI row must be A through P')
    if not combined.REDI_col_384.str.fullmatch(r'(?:[1-9]|1[0-9]|2[0-4])').all():
        raise ValueError('REDI column must be an integer from 1 through 24')
    if not combined.REDI_plate_color_384.isin(['RED_384','BLUE_384','YELLOW_384','GREEN_384']).all():
        raise ValueError('REDI quadrant must be a declared 384 color')
    return combined
