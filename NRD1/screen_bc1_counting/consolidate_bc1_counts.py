#!/usr/bin/env python3
"""
Consolidate per-sample bc1 barcode counts into a single normalized matrix.

Reads per-sample barcode_counts.csv files from Umberto's bc_counts/ directory,
CPM-normalizes each sample using the total barcode count from the accompanying
pcr_duplicate_log.csv, and produces a wide-format CSV suitable as input for
fitness_slope_heatmap.py.

Output column naming: {treatment}_R{replicate}_Normalized_Count_{generation}
  e.g. YPD_R1_Normalized_Count_0, YPD_noIAA_R2_Normalized_Count_20

Annotations (amino_acid_position, ref_amino_acid, mutant_amino_acid, gene_name,
multiple_donor_bc0_fragments_passing_filter) are joined from the bc1 reference table.

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR = Path("/path/to")
PROJECT_DIR = BASE_DIR / "projects/NNS/20250813_NRD1_MAGESTIC_screen_reps_1-3"

DEFAULT_BC_COUNTS_DIR = "/path/to/20250813_AVITI_NRD1/bc_counts"
DEFAULT_SAMPLE_KEY    = "/path/to/20250813_AVITI_NRD1/bc_counts/sample_key.csv"
DEFAULT_REFERENCE_TABLE = str(
    BASE_DIR / "projects/NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1"
    / "processed_data/data_tables/20250820_NRD1_bc1_reference_table.tsv"
)
DEFAULT_OUTPUT = str(PROJECT_DIR / "processed_data/consolidated_normalized_counts.csv")

SCALING_FACTOR = 1_000_000  # CPM

# Rename treatments for clarity.
# All conditions except "YPD noIAA" were performed with IAA (NRD1-AID degron active).
# "YPD noIAA" (no depletion, NRD1 present) → "YPD"
# "YPD" (with IAA, NRD1 depleted) → "YPD+IAA"
# All stress conditions also had IAA → append "+IAA"
TREATMENT_RENAME = {
    "YPD noIAA":       "YPD",
    "YPD":             "YPD+IAA",
    "Actinomycin D":   "Actinomycin D+IAA",
    "Cobalt Chloride": "Cobalt Chloride+IAA",
    "Galactose":       "Galactose+IAA",
    "MMS":             "MMS+IAA",
    "Rapamycin":       "Rapamycin+IAA",
    "Sorbitol":        "Sorbitol+IAA",
}


# ============================================================================
# HELPERS
# ============================================================================

def sanitize_treatment(name: str) -> str:
    """Replace whitespace with underscores so column names are valid identifiers."""
    return re.sub(r"\s+", "_", name.strip())


def get_total_barcodes(log_path: Path) -> int:
    """Read 'Total barcodes: N' from the first line of a pcr_duplicate_log file."""
    with open(log_path) as fh:
        for line in fh:
            if line.startswith("Total barcodes:"):
                return int(line.split(":")[1].strip())
    raise ValueError(f"'Total barcodes:' not found in {log_path}")


def parse_filename(filename: str):
    """
    Extract (plate, sample_num) from e.g. '20250815_plate_3B_sample_1_barcode_counts.csv'.
    Returns (plate_str, sample_int) or (None, None) on failure.
    """
    m = re.search(r"plate_([^_]+)_sample_(\d+)_barcode_counts", filename)
    if m:
        return m.group(1), int(m.group(2))
    return None, None


# ============================================================================
# MAIN
# ============================================================================

def build_consolidated_counts(
    bc_counts_dir: str,
    sample_key_path: str,
    reference_table_path: str,
    output_path: str,
) -> None:
    bc_counts_dir = Path(bc_counts_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load sample key
    # ------------------------------------------------------------------
    sample_key = pd.read_csv(sample_key_path)
    sample_key["Plate"] = sample_key["Plate"].astype(str).str.strip()
    sample_key["Sample"] = sample_key["Sample"].astype(int)
    # Replicate number from Set column (R1 → 1, R2 → 2)
    sample_key["replicate"] = sample_key["Set"].str.extract(r"(\d+)$").astype(int)
    print(f"Sample key: {len(sample_key)} rows, "
          f"treatments: {sorted(sample_key['Treatment'].unique())}")

    # ------------------------------------------------------------------
    # 2. Load and CPM-normalize each per-sample count file
    # ------------------------------------------------------------------
    count_files = sorted(bc_counts_dir.glob("*_barcode_counts.csv"))
    # Exclude sample_key.csv if it happens to match the glob
    count_files = [f for f in count_files if "sample_key" not in f.name]
    print(f"Found {len(count_files)} barcode_counts files")

    # Collect one Series per sample, indexed by Barcode
    sample_series: list[pd.Series] = []
    skipped = 0

    for fpath in count_files:
        plate, sample_num = parse_filename(fpath.name)
        if plate is None:
            print(f"  SKIP (can't parse name): {fpath.name}")
            skipped += 1
            continue

        # Look up in sample key
        row = sample_key[
            (sample_key["Plate"] == plate) & (sample_key["Sample"] == sample_num)
        ]
        if row.empty:
            print(f"  SKIP (not in sample key): plate={plate} sample={sample_num}")
            skipped += 1
            continue

        treatment  = row["Treatment"].values[0]
        passage    = int(row["Passage"].values[0])
        replicate  = int(row["replicate"].values[0])

        # Get total barcodes from log file
        log_path = fpath.with_name(fpath.name.replace("_barcode_counts.csv", "_pcr_duplicate_log.csv"))
        if not log_path.exists():
            print(f"  SKIP (no log): {log_path.name}")
            skipped += 1
            continue

        try:
            total_barcodes = get_total_barcodes(log_path)
        except Exception as e:
            print(f"  SKIP (log error): {log_path.name}: {e}")
            skipped += 1
            continue

        if total_barcodes == 0:
            print(f"  SKIP (zero total barcodes): {fpath.name}")
            skipped += 1
            continue

        # Load counts and normalize
        counts = pd.read_csv(fpath)
        counts["cpm"] = (counts["Count"] / total_barcodes) * SCALING_FACTOR

        # Build column name: {treatment}_R{replicate}_Normalized_Count_{generation}
        treatment = TREATMENT_RENAME.get(treatment, treatment)
        col_name = f"{sanitize_treatment(treatment)}_R{replicate}_Normalized_Count_{passage}"

        s = counts.set_index("Barcode")["cpm"].rename(col_name)
        sample_series.append(s)

    print(f"Loaded {len(sample_series)} samples ({skipped} skipped)")

    # ------------------------------------------------------------------
    # 3. Build wide matrix: outer-join all samples, fill NA with 0
    # ------------------------------------------------------------------
    wide = pd.concat(sample_series, axis=1).fillna(0)
    wide.index.name = "bc1"
    wide = wide.reset_index()

    # Sort columns: bc1 first, then count columns sorted alphabetically
    norm_cols = sorted([c for c in wide.columns if c != "bc1"])
    wide = wide[["bc1"] + norm_cols]
    print(f"Wide matrix: {len(wide)} barcodes × {len(norm_cols)} samples")

    # ------------------------------------------------------------------
    # 4. Join with bc1 reference table for annotations
    # ------------------------------------------------------------------
    ref = pd.read_csv(reference_table_path, sep="\t", low_memory=False)

    # Select and rename annotation columns
    annotation_cols = {
        "bc1": "bc1",
        "common_ORF_name": "gene_name",
        "aa_nums_str": "amino_acid_position",
        "ref_aas_str": "ref_amino_acid",
        "alt_aas_str": "mutant_amino_acid",
        "multiple_donor_bc0_fragments_passing_filter": "multiple_donor_bc0_fragments_passing_filter",
    }
    ref_sub = ref[list(annotation_cols.keys())].rename(columns=annotation_cols)

    # amino_acid_position may be float (e.g. 94.0) — convert to nullable int
    ref_sub["amino_acid_position"] = pd.to_numeric(
        ref_sub["amino_acid_position"], errors="coerce"
    ).astype("Int64")

    n_before = len(wide)
    wide = wide.merge(ref_sub, on="bc1", how="left")
    print(f"After reference join: {len(wide)} rows (was {n_before}); "
          f"{wide['amino_acid_position'].notna().sum()} annotated with aa position")

    # Reorder: annotation columns first, then count columns
    annot_cols = ["bc1", "gene_name", "amino_acid_position", "ref_amino_acid",
                  "mutant_amino_acid", "multiple_donor_bc0_fragments_passing_filter"]
    # Only keep annotation columns that actually exist
    annot_cols = [c for c in annot_cols if c in wide.columns]
    wide = wide[annot_cols + norm_cols]

    # ------------------------------------------------------------------
    # 5. Save
    # ------------------------------------------------------------------
    wide.to_csv(output_path, index=False)
    print(f"\nSaved: {output_path}")
    print(f"  Dimensions: {len(wide)} rows × {len(wide.columns)} columns")
    print(f"  Treatments in output:")
    treatments_seen = set()
    for col in norm_cols:
        t = re.sub(r"_R\d+_Normalized_Count_\d+$", "", col)
        treatments_seen.add(t)
    for t in sorted(treatments_seen):
        print(f"    {t}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Consolidate per-sample bc1 counts into normalized wide-format CSV"
    )
    parser.add_argument(
        "--bc-counts-dir", default=DEFAULT_BC_COUNTS_DIR,
        help="Directory with per-sample *_barcode_counts.csv files"
    )
    parser.add_argument(
        "--sample-key", default=DEFAULT_SAMPLE_KEY,
        help="Path to sample_key.csv (Plate,Sample,Treatment,Passage,Set)"
    )
    parser.add_argument(
        "--reference-table", default=DEFAULT_REFERENCE_TABLE,
        help="Path to bc1 reference table TSV"
    )
    parser.add_argument(
        "--output", default=DEFAULT_OUTPUT,
        help="Output CSV path"
    )
    args = parser.parse_args()

    build_consolidated_counts(
        bc_counts_dir=args.bc_counts_dir,
        sample_key_path=args.sample_key,
        reference_table_path=args.reference_table,
        output_path=args.output,
    )
