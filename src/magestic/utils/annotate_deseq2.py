"""
magestic.utils.annotate_deseq2 — Re-annotate DESeq2 results with the
current harmonized oligo table.

This is a thin CLI driver over `oligo_annotations.annotate_dataframe`.
It walks a directory of DESeq2 result TSVs, left-joins each on `oligo_name`
to the harmonized table at the canonical path (V574-V597-aware after the
2026-05-30 cAMP/Ras/PKA ingestion), and writes enriched copies.

Use case (HANDOFF-B yL410-yL433 cAMP/Ras/PKA):
    Existing Jun 2025 DESeq2 outputs already carry oligo_name but predate
    the V574-V597 V_aliases tagging. This CLI adds V_aliases + oligo_pool
    + all harmonized annotation columns without re-running the full
    bc0_donor_bc1 -> bc1_from_screen pipeline.

Author: Kevin R. Roy
Date: 2026-05-30
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

from .oligo_annotations import (
    HARMONIZED_OLIGO_TABLE,
    annotate_dataframe,
    load_oligo_annotations,
)

logger = logging.getLogger(__name__)


def annotate_one(
    in_tsv: Path,
    out_tsv: Path,
    annotation_df: pd.DataFrame,
    oligo_col: str = "oligo_name",
    overwrite: bool = False,
) -> dict:
    """Annotate one DESeq2 TSV; return a small summary dict."""
    if out_tsv.exists() and not overwrite:
        logger.info("SKIP (exists): %s", out_tsv.name)
        return {"file": in_tsv.name, "skipped": True}

    df = pd.read_csv(in_tsv, sep="\t")
    n_in = len(df)
    if oligo_col not in df.columns:
        logger.warning("  %s: no \"%s\" column; copying unchanged", in_tsv.name, oligo_col)
        df.to_csv(out_tsv, sep="\t", index=False)
        return {"file": in_tsv.name, "rows": n_in, "matched": 0, "no_oligo_col": True}

    # left-join — drop duplicate oligo_name column from annotation_df via merge on=oligo_col
    annot = annotation_df.copy()
    if oligo_col != "oligo_name" and "oligo_name" in annot.columns:
        annot = annot.rename(columns={"oligo_name": oligo_col})

    # Avoid colliding columns (DESeq2 already may have perfect_guide etc. from
    # the bc1_reference join). Add a suffix to the annotation side for collisions.
    overlap = (set(annot.columns) & set(df.columns)) - {oligo_col}
    if overlap:
        rename = {c: f"{c}__harmonized" for c in overlap}
        annot = annot.rename(columns=rename)

    merged = df.merge(annot, on=oligo_col, how="left")
    matched = merged["V_aliases"].notna().sum() if "V_aliases" in merged.columns else 0
    pct = 100.0 * matched / n_in if n_in else 0.0

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_tsv, sep="\t", index=False)
    logger.info(
        "  %s -> %s rows=%d matched=%d (%.1f%%)",
        in_tsv.name, out_tsv.name, n_in, matched, pct,
    )
    return {"file": in_tsv.name, "rows": n_in, "matched": int(matched), "pct": pct}


def _iter_inputs(input_dir: Path, glob: str) -> Iterable[Path]:
    for p in sorted(input_dir.glob(glob)):
        if p.is_file():
            yield p


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="magestic-annotate-deseq2",
        description=(
            "Re-annotate DESeq2 result TSVs by joining oligo_name to the "
            "current harmonized oligo table."
        ),
    )
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Directory containing DESeq2 TSVs")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Directory to write enriched TSVs")
    parser.add_argument("--glob", default="*DESeq2_results*.tsv",
                        help="Glob for input TSVs (default: %(default)s)")
    parser.add_argument("--harmonized", type=Path, default=None,
                        help="Override harmonized table path "
                             f"(default: {HARMONIZED_OLIGO_TABLE})")
    parser.add_argument("--oligo-col", default="oligo_name",
                        help="Column with oligo name (default: oligo_name)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing outputs")
    parser.add_argument("--summary-tsv", type=Path, default=None,
                        help="Optional per-file summary TSV path")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    src = args.input_dir
    dst = args.output_dir
    if not src.exists():
        logger.error("Input dir does not exist: %s", src)
        return 2

    logger.info("Loading harmonized oligo table (this may take ~30s)...")
    annotation_df = load_oligo_annotations(filepath=args.harmonized)
    logger.info("  %d oligos, %d columns", len(annotation_df), len(annotation_df.columns))

    inputs = list(_iter_inputs(src, args.glob))
    if not inputs:
        logger.error("No inputs matched %s in %s", args.glob, src)
        return 3
    logger.info("Found %d DESeq2 TSVs to re-annotate -> %s", len(inputs), dst)

    summary_rows = []
    for p in inputs:
        out = dst / p.name
        summary_rows.append(
            annotate_one(p, out, annotation_df,
                         oligo_col=args.oligo_col,
                         overwrite=args.overwrite)
        )

    if args.summary_tsv:
        args.summary_tsv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(summary_rows).to_csv(args.summary_tsv, sep="\t", index=False)
        logger.info("Wrote summary -> %s", args.summary_tsv)

    return 0


def _cli() -> None:
    sys.exit(main())


if __name__ == "__main__":
    _cli()
