"""
Step C — extract bc1 + REDI_bc from per-sample counts tables; produce the
`combined_REDI_bc_MAGESTIC_bc_merged_df.tsv` and
`MAGESTIC_bc_only_grouped_df.tsv` artifacts.

Ports legacy
`step_c_process_REDI_bc_and_MAGESTIC_bc.py` to the engine. The legacy script
glob'd `R1_R2_counts_*.tsv` from a hard-coded directory; this version uses
`config.paths.counts_tables`.

Read-pair attribution rules (which read holds bc1 vs REDI_bc) are config-driven
via `config.primer_rules` rather than hard-coded `KR1965`/`KR1967` strings.
"""

from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd

from ..core.config import REDIConfig
from ..readers.fastq import extract_bc1_from_seq, extract_redi_bc_vectorized
from ..readers.keyfiles import load_redi_bc_keys


def run_step_c(
    config: REDIConfig,
    combined_key_path: Path,
) -> tuple[Path, Path]:
    """Run step_c. Returns (path to combined_REDI_bc_MAGESTIC_bc_merged_df,
    path to MAGESTIC_bc_only_grouped_df)."""
    combined_data = config.paths.combined_data
    combined_data.mkdir(parents=True, exist_ok=True)
    config.paths.plots.mkdir(parents=True, exist_ok=True)

    combined_key = pd.read_csv(combined_key_path, sep="\t")
    combined_key["sequencing_date"] = config.sequencing_date
    combined_key["sample_name"] = (
        combined_key["sequencing_date"]
        + "_"
        + combined_key["outer_primers"]
        + "_"
        + combined_key["sample_number"]
    )

    filepaths = sorted(
        glob.glob(str(config.paths.counts_tables / "*.tsv"))
    )
    if not filepaths:
        raise FileNotFoundError(
            f"No *_R1_R2_counts.tsv files found in {config.paths.counts_tables}"
        )
    counts = pd.concat(
        (pd.read_csv(fp, sep="\t") for fp in filepaths),
        ignore_index=True,
    )

    bc = config.barcode
    counts["bc1_R1"] = counts["R1_seq"].apply(
        lambda s: extract_bc1_from_seq(
            s,
            fwd_bc1_prefix=bc.fwd_bc1_prefix,
            bc1_length=bc.bc1_length,
            fallback_prefix_length=bc.fallback_prefix_length,
            fwd_bc1_suffix=bc.fwd_bc1_suffix,
            length_tolerance=bc.bc1_length_tolerance,
        )
    )
    counts["bc1_R2"] = counts["R2_seq"].apply(
        lambda s: extract_bc1_from_seq(
            s,
            fwd_bc1_prefix=bc.fwd_bc1_prefix,
            bc1_length=bc.bc1_length,
            fallback_prefix_length=bc.fallback_prefix_length,
            fwd_bc1_suffix=bc.fwd_bc1_suffix,
            length_tolerance=bc.bc1_length_tolerance,
        )
    )

    counts = counts.merge(combined_key, on="sample_name", how="left")

    # Primer-rule dispatch: for each inner-primer match, pick bc1 from the
    # configured side.
    counts["bc1"] = None
    for primer, side in config.primer_rules.items():
        mask = counts["inner_primers"].str.contains(primer, na=False)
        col = f"bc1_{side}"
        counts.loc[mask, "bc1"] = counts.loc[mask, col]

    magestic_only = counts.query('purpose == "detect all bc1 on the plate"').copy()
    redi_df = counts.query('purpose != "detect all bc1 on the plate"').copy()

    # REDI_bc extraction via per-row prefix from the combined key
    redi_df["REDI_bc_R1"] = extract_redi_bc_vectorized(
        redi_df["R1_seq"],
        redi_df["MAT_a_barcode_prefix"],
        redi_df["MAT_a_barcode_length"],
    )
    redi_df["REDI_bc_R2"] = extract_redi_bc_vectorized(
        redi_df["R2_seq"],
        redi_df["MAT_a_barcode_prefix"],
        redi_df["MAT_a_barcode_length"],
    )

    # Per-primer side dispatch for REDI_bc (opposite of bc1 by design)
    redi_df["REDI_bc"] = None
    for primer, bc1_side in config.primer_rules.items():
        redi_side = "R2" if bc1_side == "R1" else "R1"
        mask = redi_df["inner_primers"].str.contains(primer, na=False)
        redi_df.loc[mask, "REDI_bc"] = redi_df.loc[mask, f"REDI_bc_{redi_side}"]

    redi_df = redi_df.drop(
        columns=[
            c
            for c in ["R1_seq", "R2_seq", "REDI_bc_R1", "REDI_bc_R2", "sample_number"]
            if c in redi_df.columns
        ]
    )

    redi_grouped = (
        redi_df.groupby(["bc1", "REDI_bc", "sample_name"])
        .agg({"counts": "sum"})
        .reset_index()
    )
    magestic_only_grouped = (
        magestic_only.groupby(["bc1", "sample_name"])
        .agg({"counts": "sum"})
        .reset_index()
    )

    magestic_only_path = combined_data / "MAGESTIC_bc_only_grouped_df.tsv"
    magestic_only_grouped.to_csv(magestic_only_path, sep="\t", index=False)

    # Load + consolidate REDI_bc keys (V5 + yT182/yT183 pattern from legacy step_c)
    redi_bc_key = _load_redi_bc_keys_from_config(config)

    annotated = redi_grouped.merge(combined_key, on="sample_name", how="left")
    if not redi_bc_key.empty:
        annotated = annotated.merge(redi_bc_key, on="REDI_bc", how="left")

    out_path = combined_data / "combined_REDI_bc_MAGESTIC_bc_merged_df.tsv"
    annotated.to_csv(out_path, sep="\t", index=False)
    return out_path, magestic_only_path


def _load_redi_bc_keys_from_config(config: REDIConfig) -> pd.DataFrame:
    """Resolve REDI_bc key paths from `config.redi_bc_keys` and consolidate."""
    if not config.redi_bc_keys:
        return pd.DataFrame()
    v5_path = None
    yt_path = None
    for keyfile in config.redi_bc_keys:
        p = Path(keyfile)
        if not p.is_absolute():
            p = config.paths.keyfiles / p
        name = p.name.lower()
        if "v5" in name or "redi_bc_v5" in name or "redi_bc_key" in name:
            v5_path = p
        elif "yt182" in name or "yt183" in name:
            yt_path = p
    return load_redi_bc_keys(v5_path, yt_path)
