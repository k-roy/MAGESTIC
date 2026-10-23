#!/usr/bin/env python
"""
Snapshot the lab-shared agar plate barcode database (xlsx) to TSV + parquet.

Source (read-only — NEVER written by this script):
    /oak/.../lab_shared/SGTC_genome_editing_group/
    Rectangular_plate_barcode_database_for_PIXL_Spinnaker.xlsx

Outputs (default in this script's directory):
    plates.tsv                                — flat union of all 7 sheets
    plates/Source.parquet                     — sheet "Source 1000-1999"
    plates/Target.parquet                     — sheet "Target 2000-2999"
    plates/REDI_mating.parquet                — sheet "REDI_mating 3000-3999"
    plates/Misc.parquet                       — sheet "misc 4000-4999"
    plates/Rearray_Target.parquet             — sheet "Rearray_Target 5000-5999"
    plates/Consolidated_384.parquet           — sheet "Consolidated_384"
    plates/Misc_7000.parquet                 — optional sheet "misc_7000-7999"

Snapshot is read-only by design — answers OQ-27 / OQ-29 in MASTER_PLAN §12.
Multi-user write-contention on the live xlsx is avoided.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import tempfile

import pandas as pd
from magestic.pipelines.redi.core.io import provenance, publish_tables
from magestic.pipelines.redi.readers.plate_db import _plate_frame

DEFAULT_XLSX = (
    "/path/to/lab_shared/"
    "SGTC_genome_editing_group/"
    "Rectangular_plate_barcode_database_for_PIXL_Spinnaker.xlsx"
)

# Map snapshot stage name -> list of candidate sheet-name fragments.
# Excel sheet names sometimes have trailing whitespace / minor punctuation
# variants; match by case-insensitive substring against the SHEET name.
SHEET_FRAGMENTS: dict[str, list[str]] = {
    "Source": ["source"],
    "Target": ["target_plates_2000", "target_plates 2000", "target 2000", "target_2000"],
    "REDI_mating": ["redi_mating", "redi mating"],
    "Misc": ["misc 4000", "misc_4000", "misc4000"],
    "Rearray_Target": ["rearray_target", "rearray target", "rearray"],
    "Consolidated_384": ["consolidated_384", "consolidated 384"],
    "Misc_7000": ["misc 7000", "misc_7000", "misc7000"],
}
REQUIRED_STAGES = {"Source", "Target", "REDI_mating", "Misc", "Rearray_Target"}


def _match_sheet(sheet_names: list[str], fragments: list[str]) -> str | None:
    hits = [name for name in sheet_names
            if any(fragment.casefold() in name.casefold() for fragment in fragments)]
    if len(hits) > 1:
        raise ValueError(f'ambiguous plate snapshot sheet mapping: {hits}')
    return hits[0] if hits else None


def export(xlsx_path: Path, out_dir: Path) -> dict[str, Path]:
    """Export a complete, uniquely mapped snapshot with completion attestation.

    Required sheets and all primary identities are checked before publication.
    Optional older-era omissions are allowed, but cannot silently mix with old
    parquet files already present in the chosen destination.
    """
    xlsx_path, out_dir = Path(xlsx_path), Path(out_dir)
    flat_path, parquet_dir = out_dir/'plates.tsv', out_dir/'plates'
    prov = provenance('plate_db_export', flat_path)
    prov.add_input(xlsx_path, role='primary_plate_workbook')
    prov.add_input(Path(__file__).resolve(), role='exporter_implementation')
    frames = {}
    with pd.ExcelFile(xlsx_path, engine='openpyxl') as xl:
        mappings = {stage:_match_sheet(xl.sheet_names,fragments)
                    for stage,fragments in SHEET_FRAGMENTS.items()}
        missing = sorted(stage for stage in REQUIRED_STAGES if mappings[stage] is None)
        if missing:
            raise ValueError(f'missing required plate snapshot sheets: {missing}')
        selected = [sheet for sheet in mappings.values() if sheet is not None]
        if len(selected) != len(set(selected)):
            raise ValueError('ambiguous source sheet selected for multiple plate stages')
        for stage,sheet in mappings.items():
            if sheet is None:
                if (parquet_dir/f'{stage}.parquet').exists():
                    raise ValueError(f'existing snapshot contains omitted optional {stage}; use a new snapshot directory')
                continue
            frame = _plate_frame(stage,pd.read_excel(xl,sheet_name=sheet),'plate_id')
            if {'__stage','__source_sheet'} & set(frame):
                raise ValueError('reserved plate snapshot provenance columns in source')
            frame['__stage'], frame['__source_sheet'] = stage,sheet
            # Preserve heterogeneous descriptive cells as strings, while keeping
            # the validated primary plate_id nullable and exactly integral.
            for column in frame:
                if frame[column].dtype == object:
                    frame[column] = frame[column].astype('string')
            frames[stage] = frame
    union = pd.concat(frames.values(),ignore_index=True,sort=False)
    if union.loc[union.plate_id.notna(),'plate_id'].duplicated().any():
        raise ValueError('duplicate primary plate identifier across snapshot sheets')
    out_dir.mkdir(parents=True,exist_ok=True)
    written = {stage:parquet_dir/f'{stage}.parquet' for stage in frames}
    prov.set_parameters(sheet_mapping={stage:sheet for stage,sheet in mappings.items() if sheet is not None},
                        optional_sheets_omitted=[stage for stage,sheet in mappings.items() if sheet is None])
    prov.set_counts(sheets=len(frames),rows=len(union),primary_plate_ids=int(union.plate_id.notna().sum()),
                    **{stage+':rows':len(frame) for stage,frame in frames.items()})
    with tempfile.TemporaryDirectory(prefix='plate_snapshot_',dir=out_dir) as temporary:
        files=[]
        for stage,frame in frames.items():
            staged=Path(temporary)/f'{stage}.parquet'
            frame.to_parquet(staged,index=False)
            files.append((staged,written[stage],len(frame)))
        publish_tables([(flat_path,union)],prov,files=files)
    written['__union_tsv']=flat_path
    return written


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--xlsx",
        type=Path,
        default=Path(DEFAULT_XLSX),
        help=f"Source xlsx (default: {DEFAULT_XLSX})",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Output directory (default: this script's directory)",
    )
    args = p.parse_args(argv)
    if not args.xlsx.exists():
        raise SystemExit(f"Source xlsx not found: {args.xlsx}")
    export(args.xlsx, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
