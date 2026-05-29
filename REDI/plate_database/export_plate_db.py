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

Snapshot is read-only by design — answers OQ-27 / OQ-29 in MASTER_PLAN §12.
Multi-user write-contention on the live xlsx is avoided.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

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
    "Target": ["target 2000", "target_2000"],
    "REDI_mating": ["redi_mating", "redi mating"],
    "Misc": ["misc"],
    "Rearray_Target": ["rearray_target", "rearray target", "rearray"],
    "Consolidated_384": ["consolidated_384", "consolidated 384"],
}


def _match_sheet(sheet_names: list[str], fragments: list[str]) -> str | None:
    sheet_names_lc = [s.lower() for s in sheet_names]
    # Match exactly the highest-priority fragment first (e.g., 'rearray_target'
    # should match before plain 'rearray' which would also match the consolidated
    # 'Rearray_Target' if 'rearray' came first).
    for frag in fragments:
        for orig, lc in zip(sheet_names, sheet_names_lc):
            if frag in lc:
                return orig
    return None


def export(xlsx_path: Path, out_dir: Path) -> dict[str, Path]:
    """Read xlsx, write per-stage parquet + flat plates.tsv."""
    out_dir.mkdir(parents=True, exist_ok=True)
    parquet_dir = out_dir / "plates"
    parquet_dir.mkdir(exist_ok=True)

    xl = pd.ExcelFile(xlsx_path, engine="openpyxl")
    sheet_names = xl.sheet_names
    print(f"Loaded {xlsx_path} with sheets: {sheet_names}")

    written: dict[str, Path] = {}
    union_parts: list[pd.DataFrame] = []
    for stage, fragments in SHEET_FRAGMENTS.items():
        sheet = _match_sheet(sheet_names, fragments)
        if sheet is None:
            print(f"  WARN: no sheet matched for stage '{stage}'")
            continue
        df = pd.read_excel(xl, sheet_name=sheet)
        df["__stage"] = stage
        df["__source_sheet"] = sheet
        parquet_path = parquet_dir / f"{stage}.parquet"
        df.to_parquet(parquet_path, index=False)
        written[stage] = parquet_path
        union_parts.append(df)
        print(f"  {stage:18s}  ← '{sheet}'  rows={len(df)}")

    if union_parts:
        union = pd.concat(union_parts, ignore_index=True, sort=False)
        flat_path = out_dir / "plates.tsv"
        union.to_csv(flat_path, sep="\t", index=False)
        written["__union_tsv"] = flat_path
        print(f"Wrote union: {flat_path} ({len(union)} rows)")
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
