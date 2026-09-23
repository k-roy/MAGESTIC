#!/usr/bin/env python3
"""
Assign bc1 barcodes to yeast sublibraries by integrating multiple amplicon data sources.

MULTI-SOURCE BC1 INTEGRATION
============================
This script integrates bc1-to-trafo assignments from multiple sources:

1. **bc1-donor-bc0 amplicons** (primary): Links bc1→donor→bc0→oligo with trafo_ID from sample metadata
2. **bc1-only amplicons** (supplementary): Per-trafo plate wash samples provide high-coverage bc1→trafo mapping
3. **Existing bc1_to_trafo_mapping** (pre-computed): Reuse previously computed mappings

The script combines these sources with confidence scoring based on:
- Read counts (higher counts = higher confidence)
- Source priority (configurable)
- Purity scores (fraction of counts in assigned trafo)

TERMINOLOGY CLARIFICATION
=========================
- yeast_sublibrary: The transformation batch level (trafo_1, trafo_2, trafo_3, trafo_4)
  Each yeast_sublibrary represents a unique combination of:
  - Strain/nuclease background (yKR1033=SpG, yKR1036=SpCas9, etc.)
  - V_library/oligo pool (V629, V631, V634)
  Replicate plates (trafo_1-1, trafo_1-2) belong to the same yeast_sublibrary.

- V_library: The plasmid oligo pool (V629, V631, V634)
  Determined by SPS (Sub-Pool Sequence) in the amplicon.
  NOTE: V634 is used by both trafo_3 (SpG) and trafo_4 (SpCas9), so SPS alone
  cannot distinguish these - sample metadata (trafo_ID) is required.

- nuclease_PAM: Descriptive name combining nuclease and PAM (SpG_NGNG, SpCas9_NGG, etc.)
  Derived from trafo_ID via trafo_key and strain_key.

WHY MULTI-SOURCE MATTERS
========================
Different amplicon types have different coverage:
- bc1-donor-bc0: Links bc1 to donors but may have lower bc1 coverage
- bc1-only plate wash: High bc1 coverage for trafo assignment but no donor info

By combining sources, we maximize both bc1 coverage AND donor linkage.

Usage:
    # Basic usage with bc1-donor-bc0 only (original behavior)
    python 05_assign_sublibrary.py \\
        --input <mapped_tsv> \\
        --sps-mapping <sps_library_mapping.tsv> \\
        --output-reference <output_tsv> \\
        --output-summary <summary_tsv>

    # With additional bc1-only trafo mapping
    python 05_assign_sublibrary.py \\
        --input <mapped_tsv> \\
        --sps-mapping <sps_library_mapping.tsv> \\
        --bc1-trafo-mapping <bc1_to_trafo_mapping.tsv> \\
        --trafo-key <trafo_key.tsv> \\
        --output-reference <output_tsv> \\
        --output-summary <summary_tsv>

    # With bc1-only parsed count files
    python 05_assign_sublibrary.py \\
        --input <mapped_tsv> \\
        --sps-mapping <sps_library_mapping.tsv> \\
        --bc1-counts-dir <parsed_counts_dir> \\
        --bc1-sample-key <sample_key.tsv> \\
        --trafo-condition "plate wash" \\
        --output-reference <output_tsv> \\
        --output-summary <summary_tsv>

Author: Kevin R. Roy
"""

import argparse
import pandas as pd
import yaml
from pathlib import Path
from collections import defaultdict
import re
from typing import Dict, List, Optional, Tuple
import logging
import json
import subprocess
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Mapping from V_library to PAM (used to construct nuclease_PAM)
V_LIBRARY_TO_PAM = {
    "V629": "NGNG",
    "V631": "NGNH",
    "V634": "NGG",
}


def load_sps_mapping(mapping_file):
    """
    Load SPS-to-V_library mapping from TSV or YAML file.

    TSV format: columns 'sps', 'library_id', optionally 'description'
    YAML format: dict with library_id as key, 'sps' in value

    Returns:
        dict: {sps_sequence (uppercase) -> V_library}
    """
    mapping_path = Path(mapping_file)

    if mapping_path.suffix in ['.yaml', '.yml']:
        # YAML format (legacy support)
        with open(mapping_path) as f:
            config = yaml.safe_load(f)
        libraries = config.get("libraries", config)
        return {v["sps"].upper(): k for k, v in libraries.items()}
    else:
        # TSV format (preferred)
        df = pd.read_csv(mapping_path, sep="\t")

        # Validate columns
        if "sps" not in df.columns or "library_id" not in df.columns:
            raise ValueError(
                f"SPS mapping file must have 'sps' and 'library_id' columns. "
                f"Found: {df.columns.tolist()}"
            )

        # Normalize SPS to uppercase
        return dict(zip(df["sps"].str.upper(), df["library_id"]))


def load_trafo_key(trafo_key_file):
    """
    Load trafo key mapping trafo_ID to strain (yKR), V_library, and other metadata.

    Returns:
        dict: {trafo_ID -> {'yKR': strain, 'V_library': V_library, 'yL': yL_number, ...}}
    """
    df = pd.read_csv(trafo_key_file, sep="\t")

    trafo_map = {}
    for _, row in df.iterrows():
        trafo_id = row.get("trafo_ID") or row.get("trafo_id")
        if not trafo_id or pd.isna(trafo_id):
            continue

        trafo_map[trafo_id] = {
            "yKR": row.get("yKR") or row.get("ykr"),
            "V_library": row.get("corresponding_step_1_V") or row.get("V_library"),
            "yL": row.get("yL") or row.get("yeast_library"),
            "nuclease": row.get("nuclease"),
            "nuclease_PAM": row.get("nuclease_PAM") or row.get("yeast_sublibrary"),
        }

    return trafo_map


def load_strain_key(strain_key_file):
    """Load strain key mapping yKR to nuclease name."""
    df = pd.read_csv(strain_key_file, sep="\t")
    return dict(zip(df["yKR"], df["nuclease"]))


def build_trafo_to_nuclease_pam(trafo_key, strain_key=None):
    """
    Build mapping from trafo_ID to nuclease_PAM.

    Returns:
        dict: {trafo_ID -> nuclease_PAM (e.g., 'SpG_NGNG')}
    """
    trafo_to_nuclease_pam = {}

    for trafo_id, info in trafo_key.items():
        # First check if nuclease_PAM is directly in the trafo_key
        if info.get("nuclease_PAM"):
            trafo_to_nuclease_pam[trafo_id] = info["nuclease_PAM"]
            continue

        # Otherwise, derive from nuclease + V_library
        yKR = info.get("yKR")
        V_library = info.get("V_library")
        nuclease = info.get("nuclease")

        # Get nuclease from strain_key if not in trafo_key
        if not nuclease and strain_key and yKR:
            nuclease = strain_key.get(yKR, "unknown")

        pam = V_LIBRARY_TO_PAM.get(V_library, "unknown") if V_library else "unknown"

        if nuclease:
            nuclease_clean = nuclease.replace(" ", "_").replace("WT_", "")
            trafo_to_nuclease_pam[trafo_id] = f"{nuclease_clean}_{pam}"
        else:
            trafo_to_nuclease_pam[trafo_id] = f"unknown_{pam}"

    return trafo_to_nuclease_pam


def load_bc1_trafo_mapping(mapping_file: Path) -> Dict[str, Dict]:
    """
    Load existing bc1→trafo mapping from TSV file.

    Expected columns: bc1, yeast_sublibrary (or trafo_ID), purity, total_counts

    Returns:
        dict: {bc1 -> {'trafo_ID': str, 'purity': float, 'counts': int, 'source': str}}
    """
    df = pd.read_csv(mapping_file, sep="\t")

    # Find trafo column
    trafo_col = None
    for col in ["yeast_sublibrary", "trafo_ID", "trafo_id", "sublibrary"]:
        if col in df.columns:
            trafo_col = col
            break

    if trafo_col is None:
        logger.warning(f"No trafo column found in {mapping_file}")
        return {}

    # Find purity and counts columns
    purity_col = "purity" if "purity" in df.columns else None
    counts_col = None
    for col in ["total_counts", "trafo_counts", "counts"]:
        if col in df.columns:
            counts_col = col
            break

    mapping = {}
    for _, row in df.iterrows():
        bc1 = row.get("bc1")
        if pd.isna(bc1):
            continue

        mapping[bc1] = {
            "trafo_ID": str(row[trafo_col]) if pd.notna(row[trafo_col]) else None,
            "purity": float(row[purity_col]) if purity_col and pd.notna(row.get(purity_col)) else 1.0,
            "counts": int(row[counts_col]) if counts_col and pd.notna(row.get(counts_col)) else 0,
            "source": "bc1_trafo_mapping",
        }

    logger.info(f"Loaded {len(mapping)} bc1→trafo mappings from {mapping_file.name}")
    return mapping


def build_bc1_trafo_from_parsed_counts(
    counts_dir: Path,
    sample_key: pd.DataFrame,
    trafo_condition: str = "plate wash"
) -> Dict[str, Dict]:
    """
    Build bc1→trafo mapping from per-sample parsed count files.

    This is used for bc1-only amplicons where each sample belongs to a specific trafo.

    Args:
        counts_dir: Directory containing parsed bc1 count files
        sample_key: DataFrame with sample_number, trafo_ID, condition columns
        trafo_condition: Filter samples to this condition (e.g., "plate wash")

    Returns:
        dict: {bc1 -> {'trafo_ID': str, 'purity': float, 'counts': int, 'source': str}}
    """
    # Find trafo column in sample key
    trafo_col = None
    for col in ["trafo_ID", "trafo_id", "yeast_sublibrary"]:
        if col in sample_key.columns:
            trafo_col = col
            break

    if trafo_col is None:
        logger.warning("No trafo column found in sample key")
        return {}

    # Filter to relevant condition
    condition_col = "condition" if "condition" in sample_key.columns else None
    if condition_col and trafo_condition:
        filtered_key = sample_key[
            sample_key[condition_col].str.lower().str.contains(trafo_condition.lower(), na=False)
        ]
        logger.info(f"Filtered to {len(filtered_key)} '{trafo_condition}' samples")
    else:
        filtered_key = sample_key
        logger.info(f"Using all {len(filtered_key)} samples (no condition filter)")

    # Build sample_number → trafo_ID mapping
    sample_to_trafo = {}
    for _, row in filtered_key.iterrows():
        sample_num = str(row.get("sample_number", "")).replace("sample_", "")
        trafo_id = row.get(trafo_col)
        if sample_num and pd.notna(trafo_id):
            # Extract base trafo (trafo_1-1 → trafo_1)
            trafo_base = str(trafo_id).split("-")[0]
            sample_to_trafo[sample_num] = trafo_base

    logger.info(f"Built sample→trafo mapping for {len(sample_to_trafo)} samples")

    # Aggregate bc1 counts by trafo
    bc1_trafo_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    # Find parsed count files
    parsed_files = list(counts_dir.glob("*_bc1_counts.tsv"))
    if not parsed_files:
        parsed_files = list(counts_dir.glob("*_bc1*.tsv"))

    logger.info(f"Found {len(parsed_files)} parsed count files in {counts_dir}")

    matched_files = 0
    for parsed_file in parsed_files:
        # Extract sample number from filename
        match = re.search(r'sample_(\d+)', parsed_file.name)
        if not match:
            continue

        sample_num = match.group(1)
        if sample_num not in sample_to_trafo:
            continue

        trafo_id = sample_to_trafo[sample_num]
        matched_files += 1

        # Load bc1 counts
        try:
            df = pd.read_csv(parsed_file, sep="\t")

            # Find bc1 and counts columns
            bc1_col = "bc1" if "bc1" in df.columns else df.columns[0]
            counts_col = "counts" if "counts" in df.columns else df.columns[-1]

            for _, row in df.iterrows():
                bc1 = str(row[bc1_col])
                count = int(row[counts_col])
                bc1_trafo_counts[bc1][trafo_id] += count

        except Exception as e:
            logger.debug(f"Failed to load {parsed_file.name}: {e}")

    logger.info(f"Processed {matched_files} files matching condition samples")
    logger.info(f"Found {len(bc1_trafo_counts)} unique bc1s")

    # Assign each bc1 to trafo with highest counts
    mapping = {}
    for bc1, trafo_counts in bc1_trafo_counts.items():
        total_counts = sum(trafo_counts.values())
        if total_counts == 0:
            continue

        best_trafo = max(trafo_counts.items(), key=lambda x: x[1])
        purity = best_trafo[1] / total_counts

        mapping[bc1] = {
            "trafo_ID": best_trafo[0],
            "purity": purity,
            "counts": total_counts,
            "source": "bc1_only_counts",
        }

    return mapping


def merge_bc1_trafo_sources(
    sources: List[Dict[str, Dict]],
    source_names: List[str],
    priority_weights: Optional[List[float]] = None
) -> pd.DataFrame:
    """
    Merge bc1→trafo mappings from multiple sources.

    Strategy:
    1. For each bc1, collect assignments from all sources
    2. If all sources agree, use that trafo
    3. If sources disagree, use weighted voting by (counts × purity × priority_weight)
    4. Record which source provided the final assignment

    Returns:
        DataFrame with columns: bc1, trafo_ID, purity, counts, source, confidence
    """
    if priority_weights is None:
        priority_weights = [1.0] * len(sources)

    # Collect all bc1s across sources
    all_bc1s = set()
    for source in sources:
        all_bc1s.update(source.keys())

    logger.info(f"Merging {len(all_bc1s)} unique bc1s from {len(sources)} sources")

    merged = []
    for bc1 in all_bc1s:
        assignments = []
        for i, (source, name, weight) in enumerate(zip(sources, source_names, priority_weights)):
            if bc1 in source:
                info = source[bc1]
                score = info["counts"] * info["purity"] * weight
                assignments.append({
                    "trafo_ID": info["trafo_ID"],
                    "purity": info["purity"],
                    "counts": info["counts"],
                    "source": name,
                    "score": score,
                })

        if not assignments:
            continue

        # Check if all sources agree
        trafos = set(a["trafo_ID"] for a in assignments)

        if len(trafos) == 1:
            # All sources agree - combine counts
            best = assignments[0]
            total_counts = sum(a["counts"] for a in assignments)
            if total_counts > 0:
                avg_purity = sum(a["purity"] * a["counts"] for a in assignments) / total_counts
            else:
                avg_purity = best["purity"]

            merged.append({
                "bc1": bc1,
                "trafo_ID": best["trafo_ID"],
                "purity": avg_purity,
                "counts": total_counts,
                "source": "+".join(sorted(set(a["source"] for a in assignments))),
                "confidence": "high",
            })
        else:
            # Sources disagree - use weighted voting
            best = max(assignments, key=lambda x: x["score"])
            merged.append({
                "bc1": bc1,
                "trafo_ID": best["trafo_ID"],
                "purity": best["purity"],
                "counts": best["counts"],
                "source": best["source"],
                "confidence": "medium" if best["purity"] > 0.8 else "low",
            })

    return pd.DataFrame(merged)


def assign_v_library_by_sps(row, sps_to_library):
    """Assign V_library based on SPS sequence."""
    sps = row.get("SPS") or row.get("sps")
    if pd.isna(sps):
        return None
    return sps_to_library.get(str(sps).upper())


def extract_yeast_sublibrary(trafo_id):
    """
    Extract yeast_sublibrary from trafo_ID.

    trafo_1-1, trafo_1-2 -> trafo_1 (same yeast_sublibrary)
    trafo_2 -> trafo_2
    """
    if pd.isna(trafo_id):
        return None

    # Handle replicate suffixes: trafo_1-1 -> trafo_1
    parts = str(trafo_id).split("-")
    return parts[0] if parts else None


def main():
    parser = argparse.ArgumentParser(
        description="Assign bc1 to yeast sublibraries using multiple data sources",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage (bc1-donor-bc0 only)
    python 05_assign_sublibrary.py \\
        --input mapped.tsv \\
        --sps-mapping sps_library_mapping.tsv \\
        --output-reference bc1_reference.tsv \\
        --output-summary summary.tsv

    # With pre-computed bc1→trafo mapping
    python 05_assign_sublibrary.py \\
        --input mapped.tsv \\
        --sps-mapping sps_library_mapping.tsv \\
        --bc1-trafo-mapping bc1_to_trafo_mapping.tsv \\
        --output-reference bc1_reference.tsv \\
        --output-summary summary.tsv

    # With bc1-only parsed counts
    python 05_assign_sublibrary.py \\
        --input mapped.tsv \\
        --sps-mapping sps_library_mapping.tsv \\
        --bc1-counts-dir processed_data/bc1_counts/parsed/ \\
        --bc1-sample-key keyfiles/bc1_from_screen/sample_key.tsv \\
        --trafo-condition "plate wash" \\
        --output-reference bc1_reference.tsv \\
        --output-summary summary.tsv
        """
    )

    # Required arguments
    parser.add_argument("--input", required=True, help="Input mapped bc1-donor-bc0 TSV file")
    parser.add_argument("--sps-mapping", required=True, help="SPS-to-V_library mapping (TSV or YAML)")
    parser.add_argument("--output-reference", required=True, help="Output bc1 reference table")
    parser.add_argument("--output-summary", required=True, help="Output purity summary")

    # Additional bc1 sources
    parser.add_argument("--bc1-trafo-mapping", help="Pre-computed bc1→trafo mapping TSV")
    parser.add_argument("--bc1-counts-dir", help="Directory with parsed bc1 count files")
    parser.add_argument("--bc1-sample-key", help="Sample key for bc1 counts (TSV)")
    parser.add_argument("--trafo-condition", default="plate wash",
                        help="Condition to filter for trafo assignment (default: 'plate wash')")

    # Metadata keys
    parser.add_argument("--trafo-key", help="Trafo key TSV file")
    parser.add_argument("--strain-key", help="Strain key TSV file")

    # Legacy argument
    parser.add_argument("--config", help="Deprecated: use --sps-mapping instead")

    args = parser.parse_args()

    # Handle legacy --config argument
    mapping_file = args.sps_mapping or args.config
    if not mapping_file:
        raise ValueError("--sps-mapping is required")

    # Load primary bc1-donor-bc0 data
    logger.info("Loading mapped bc1-donor-bc0 data...")
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    logger.info(f"Loaded {len(df)} rows with {df['bc1'].nunique()} unique bc1s")
    n_input_rows, n_input_bc1 = len(df), df["bc1"].nunique()
    link_key = [c for c in ("bc1", "donor", "bc0") if c in df.columns]
    if len(link_key) == 3 and df.duplicated(subset=link_key).any():
        raise SystemExit(
            f"ERROR: {int(df.duplicated(subset=link_key).sum())} duplicated (bc1, donor, bc0) links in "
            f"{args.input}. This is what per-run step-02 tables concatenated at step 04 look like "
            f"(review F2, 2026-09-10): the top link would come from ONE run. Run step 02 once over "
            f"every run with --manifest.")

    # Load SPS mapping
    logger.info(f"Loading SPS mapping from: {mapping_file}")
    sps_to_library = load_sps_mapping(mapping_file)
    logger.info(f"Loaded {len(sps_to_library)} SPS-to-V_library mappings")

    # Load trafo and strain keys if provided
    trafo_key = {}
    strain_key = {}
    trafo_to_nuclease_pam = {}

    if args.trafo_key:
        logger.info(f"Loading trafo key from: {args.trafo_key}")
        trafo_key = load_trafo_key(args.trafo_key)
        logger.info(f"Loaded {len(trafo_key)} trafo mappings")

        if args.strain_key:
            strain_key = load_strain_key(args.strain_key)

        trafo_to_nuclease_pam = build_trafo_to_nuclease_pam(trafo_key, strain_key)

    # Collect bc1→trafo mappings from all sources
    bc1_trafo_sources = []
    source_names = []

    # Source 1: bc1-donor-bc0 input (if trafo_ID column exists)
    if "trafo_ID" in df.columns or "trafo_id" in df.columns:
        trafo_col = "trafo_ID" if "trafo_ID" in df.columns else "trafo_id"
        count_col = "counts" if "counts" in df.columns else "total_counts"

        bdb_mapping = {}
        for _, row in df.iterrows():
            bc1 = row.get("bc1")
            trafo_id = row.get(trafo_col)
            if pd.notna(bc1) and pd.notna(trafo_id):
                trafo_base = extract_yeast_sublibrary(trafo_id)
                counts = int(row.get(count_col, 1))

                if bc1 not in bdb_mapping or counts > bdb_mapping[bc1]["counts"]:
                    bdb_mapping[bc1] = {
                        "trafo_ID": trafo_base,
                        "purity": 1.0,  # Assuming single source per bc1 in input
                        "counts": counts,
                        "source": "bc1_donor_bc0",
                    }

        if bdb_mapping:
            bc1_trafo_sources.append(bdb_mapping)
            source_names.append("bc1_donor_bc0")
            logger.info(f"Source 1: {len(bdb_mapping)} bc1s from bc1-donor-bc0 input")

    # Source 2: Pre-computed bc1→trafo mapping
    if args.bc1_trafo_mapping:
        mapping_path = Path(args.bc1_trafo_mapping)
        if mapping_path.exists():
            bc1_mapping = load_bc1_trafo_mapping(mapping_path)
            if bc1_mapping:
                bc1_trafo_sources.append(bc1_mapping)
                source_names.append("bc1_trafo_mapping")
                logger.info(f"Source 2: {len(bc1_mapping)} bc1s from pre-computed mapping")

    # Source 3: bc1-only parsed counts
    if args.bc1_counts_dir and args.bc1_sample_key:
        counts_dir = Path(args.bc1_counts_dir)
        sample_key_path = Path(args.bc1_sample_key)

        if counts_dir.exists() and sample_key_path.exists():
            sample_key_df = pd.read_csv(sample_key_path, sep="\t")
            bc1_counts_mapping = build_bc1_trafo_from_parsed_counts(
                counts_dir, sample_key_df, args.trafo_condition
            )
            if bc1_counts_mapping:
                bc1_trafo_sources.append(bc1_counts_mapping)
                source_names.append("bc1_only_counts")
                logger.info(f"Source 3: {len(bc1_counts_mapping)} bc1s from bc1-only counts")

    # Merge all bc1→trafo sources
    if bc1_trafo_sources:
        logger.info(f"\nMerging {len(bc1_trafo_sources)} bc1→trafo sources...")
        merged_trafo = merge_bc1_trafo_sources(bc1_trafo_sources, source_names)
        logger.info(f"Merged: {len(merged_trafo)} unique bc1s with trafo assignments")

        # Build bc1 → trafo lookup from merged data
        bc1_to_trafo = dict(zip(merged_trafo["bc1"], merged_trafo["trafo_ID"]))
        bc1_to_purity = dict(zip(merged_trafo["bc1"], merged_trafo["purity"]))
        bc1_to_source = dict(zip(merged_trafo["bc1"], merged_trafo["source"]))
        bc1_to_confidence = dict(zip(merged_trafo["bc1"], merged_trafo["confidence"]))

        # Report source distribution
        logger.info("\nTrafo assignment source distribution:")
        for src, count in merged_trafo["source"].value_counts().items():
            logger.info(f"  {src}: {count}")
    else:
        bc1_to_trafo = {}
        bc1_to_purity = {}
        bc1_to_source = {}
        bc1_to_confidence = {}

    # Assign V_library based on SPS
    if "SPS" in df.columns or "sps" in df.columns:
        logger.info("\nAssigning V_library based on SPS...")
        df["V_library"] = df.apply(lambda row: assign_v_library_by_sps(row, sps_to_library), axis=1)
        assigned = df["V_library"].notna().sum()
        logger.info(f"Assigned {assigned}/{len(df)} ({100*assigned/len(df):.1f}%) entries to V_libraries")
        df["library_ID"] = df["V_library"]  # Legacy column
    else:
        df["V_library"] = None
        df["library_ID"] = None

    # Assign yeast_sublibrary from merged trafo mapping
    logger.info("\nAssigning yeast_sublibrary from merged sources...")
    df["yeast_sublibrary"] = df["bc1"].map(bc1_to_trafo)
    df["sublibrary_purity"] = df["bc1"].map(bc1_to_purity)
    df["sublibrary_source"] = df["bc1"].map(bc1_to_source)
    df["sublibrary_confidence"] = df["bc1"].map(bc1_to_confidence)

    # Assign nuclease_PAM from trafo_key
    if trafo_to_nuclease_pam:
        df["nuclease_PAM"] = df["yeast_sublibrary"].map(trafo_to_nuclease_pam)
    else:
        df["nuclease_PAM"] = None

    assigned = df["yeast_sublibrary"].notna().sum()
    logger.info(f"Assigned {assigned}/{len(df)} ({100*assigned/len(df):.1f}%) entries to yeast_sublibraries")

    # Find counts column
    count_col = "counts" if "counts" in df.columns else "total_counts"

    # 2026-09-10: barcode<->donor UNCOUPLING, visible rather than silently first-wins
    # (the WGS lane's ask; QK-011 measured ~4.5x more uncoupling with Q5 than KAPA).
    # Distinct designed oligos over ALL of the bc1's links, computed BEFORE the
    # one-row-per-bc1 reduction below discards the minority links.
    if "oligo_name" in df.columns:
        df["maps_to_n_variants"] = df.groupby("bc1")["oligo_name"].transform("nunique").astype("int64")
        # maps_to_n_variants counts EVERY link, including 1-read error/chimera links, so on the 20210822 smoke it
        # was >= 2 for ~34% of bc1 -- too noisy to exclude on. Two better-behaved companions (29_):
        #   top_oligo_read_share: share of the bc1's reads on links naming the same oligo as its top link;
        #   maps_to_n_variants_confident: distinct oligos among links with >= 2 distinct PCR replicates.
        _cnt = "counts" if "counts" in df.columns else "total_counts"
        _top = df.sort_values([_cnt], ascending=False, kind="mergesort").drop_duplicates("bc1").set_index("bc1")["oligo_name"]
        _on_top = df["oligo_name"].notna() & (df["oligo_name"] == df["bc1"].map(_top))
        _tot = df.groupby("bc1")[_cnt].transform("sum")
        df["top_oligo_read_share"] = (df[_cnt].where(_on_top, 0).groupby(df["bc1"]).transform("sum")
                                      / _tot.where(_tot > 0)).round(4).where(df["bc1"].map(_top).notna())
        # (undefined -- NaN, not 0 -- when the bc1's top link maps to no designed oligo)
        if "support_tier" in df.columns:
            _conf = df["oligo_name"].where(df["support_tier"] == "confident")
            df["maps_to_n_variants_confident"] = _conf.groupby(df["bc1"]).transform("nunique").astype("int64")
    df["n_links"] = df.groupby("bc1")["bc1"].transform("size").astype("int64")

    # Keep top donor per bc1 for reference table. Deterministic on ties (the old
    # default quicksort was not stable, so a tie could pick a different link per run).
    tie_order = {count_col: False, "num_PCR_replicates": False, "donor": True, "bc0": True}
    tie_cols = [c for c in tie_order if c in df.columns]
    bc1_reference = df.sort_values(tie_cols, ascending=[tie_order[c] for c in tie_cols],
                                   kind="mergesort").drop_duplicates(subset=["bc1"], keep="first").copy()

    # Select columns for reference table
    ref_cols = ["bc1"]

    if "donor" in bc1_reference.columns:
        ref_cols.append("donor")
    elif "top_donor" in bc1_reference.columns:
        ref_cols.append("top_donor")

    if "bc0" in bc1_reference.columns:
        ref_cols.append("bc0")
    elif "top_bc0" in bc1_reference.columns:
        ref_cols.append("top_bc0")

    # Add assignment columns
    for col in ["SPS", "V_library", "yeast_sublibrary", "nuclease_PAM",
                "sublibrary_purity", "sublibrary_source", "sublibrary_confidence",
                "library_ID",  # Legacy
                "oligo_name", "guide", count_col,
                "num_PCR_replicates", "bc1_purity", "num_distinct_donors",
                # 2026-09-10 (planning 28_ Sec. 1e). The allow-list is why these must be
                # named: an unlisted column is dropped here silently.
                *TIER_COUNT_COLUMNS_05, "support_tier", "pcr_replicate_ids", "runs_observed",
                "yeast_libraries_observed",
                "total_bc1_counts", "maps_to_n_variants", "maps_to_n_variants_confident",
                "top_oligo_read_share", "n_links",
                "donor_suffix_ambiguous", "donor_lengths_observed",
                "guide_link_counts", "guide_link_n_guides", "bc1_orientation",
                # 2026-09-10 (planning 30_ arm_bdb-F1): the top link's match flags. oligo_name of an ambiguous
                # match is the sorted-first candidate; without these a wrong-edit row looks like a perfect one.
                "match_method", "is_unambiguous", "num_candidates", "num_candidate_variants", "variant",
                "candidate_oligos", "perfect_donor", "guide_match",
                # planning 218 item 10: what the read SHOWED about the designed edit. Without these two an
                # uncertified window candidate is indistinguishable from a certified match on this table.
                "edit_evidence", "num_edit_certified_variants",
                # planning 30_ step 6: guide_match=False conflated "no linked guide" with "the linked
                # guide disagrees"; guide_status separates them and num_guide_oligos sizes the conflict.
                "guide_status", "num_guide_oligos",
                # planning 30_ GDB-01: the linked guide was restored from gdb's 1-nt-early Cas12a extraction.
                "guide_repaired_gdb01", "guide_gdb_raw", "sps_scope_V", "sps_scoped"]:
        if col in bc1_reference.columns and col not in ref_cols:
            ref_cols.append(col)

    bc1_reference = bc1_reference[ref_cols]
    missing_contract = [c for c in ("support_tier", "bc1_orientation", *TIER_COUNT_COLUMNS_05)
                        if c not in bc1_reference.columns]
    if missing_contract:
        logger.warning(f"Reference table lacks contract columns {missing_contract} -- input predates "
                       f"the 2026-09-10 step 01/02 (tier counts / support_tier / orientation).")
    if all(c in bc1_reference.columns for c in TIER_COUNT_COLUMNS_05) and count_col in bc1_reference.columns:
        bad = int((bc1_reference[count_col] != bc1_reference[list(TIER_COUNT_COLUMNS_05)].sum(axis=1)).sum())
        if bad:
            raise SystemExit(f"ERROR: {bad} reference rows violate counts == sum(n_<tier>)")

    # Schema fidelity: the reference table has always called the top link's read count
    # ``total_counts`` (its 2026-05-27 build came via 02_simple_combine). The rebuilt
    # step 02 emits ``counts``; restore the established name so no consumer breaks.
    if count_col == "counts" and "total_counts" not in bc1_reference.columns:
        bc1_reference = bc1_reference.rename(columns={"counts": "total_counts"})

    # Save reference table
    bc1_reference.to_csv(args.output_reference, sep="\t", index=False)
    logger.info(f"\nSaved bc1 reference table ({len(bc1_reference)} bc1s) to: {args.output_reference}")

    # Generate purity summary
    group_col = "yeast_sublibrary" if df["yeast_sublibrary"].notna().any() else "V_library"

    summary_data = []
    for lib_id in df[group_col].dropna().unique():
        lib_df = df[df[group_col] == lib_id]
        lib_bc1s = lib_df.drop_duplicates(subset=["bc1"])

        summary_row = {
            group_col: lib_id,
            "num_bc1s": lib_bc1s["bc1"].nunique(),
            "total_counts": lib_df[count_col].sum(),
        }

        if "nuclease_PAM" in lib_df.columns:
            nuclease_pam = lib_df["nuclease_PAM"].dropna().unique()
            summary_row["nuclease_PAM"] = nuclease_pam[0] if len(nuclease_pam) > 0 else None

        if "V_library" in lib_df.columns:
            v_lib = lib_df["V_library"].dropna().unique()
            summary_row["V_library"] = v_lib[0] if len(v_lib) > 0 else None

        # Source distribution
        if "sublibrary_source" in lib_bc1s.columns:
            sources = lib_bc1s["sublibrary_source"].value_counts().to_dict()
            summary_row["source_distribution"] = str(sources)

        if "bc1_purity" in lib_bc1s.columns:
            summary_row["mean_purity"] = lib_bc1s["bc1_purity"].mean()
            summary_row["median_purity"] = lib_bc1s["bc1_purity"].median()

        summary_data.append(summary_row)

    # Add unassigned
    unassigned = df[df[group_col].isna()]
    if len(unassigned) > 0:
        summary_data.append({
            group_col: "unassigned",
            "num_bc1s": unassigned["bc1"].nunique(),
            "total_counts": unassigned[count_col].sum(),
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(args.output_summary, sep="\t", index=False)
    logger.info(f"Saved purity summary to: {args.output_summary}")
    # The sidecar is written LAST, after every output it checksums exists. (The 2026-09-11 smoke caught it
    # being written between the reference table and the summary: stat() on the summary raised.)
    _write_provenance(args, bc1_reference, n_input_rows, n_input_bc1)

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("BC1 REFERENCE TABLE SUMMARY")
    logger.info("=" * 60)
    print(summary_df.to_string(index=False))

    # Print multi-source statistics if applicable
    if len(bc1_trafo_sources) > 1:
        logger.info("\n" + "=" * 60)
        logger.info("MULTI-SOURCE INTEGRATION SUMMARY")
        logger.info("=" * 60)
        for name, source in zip(source_names, bc1_trafo_sources):
            logger.info(f"  {name}: {len(source)} bc1s")
        if bc1_trafo_sources:
            logger.info(f"  Merged total: {len(bc1_to_trafo)} bc1s")


# Imported unconditionally: a silent fallback to () would drop all seven tier columns
# through the allow-list and turn the invariant check into a false alarm (review F3).
from magestic.utils.bc1_extraction import TIER_COUNT_COLUMNS as TIER_COUNT_COLUMNS_05, BC1_ORIENTATION


def _git_state():
    """Commit sha, branch and dirty flag of the checkout that holds the running package.

    Delegates to magestic.utils.provenance.git_state, which runs git with cwd= rather than `git -C`:
    the cluster's system git predates -C, so the previous subprocess call exited 129 with an empty
    stdout and every 2026-09-10 reference-table sidecar shipped ``git_head: null``.
    """
    from magestic.utils.provenance import git_state
    g = git_state()
    return {"git_head": g.get("commit") or None, "git_branch": g.get("branch"), "git_repo": g.get("repo"),
            "git_tracked_dirty": bool(g.get("dirty")), "git_dirty_files": g.get("dirty_files") or []}


def _write_provenance(args, ref, n_input_rows, n_input_bc1):
    """<output_reference>.provenance.json -- what built this table, and its conventions."""
    from magestic.utils.provenance import environment, file_record, utc_now
    prov = {
        "step": "05_assign_sublibrary", "argv": sys.argv, **_git_state(),
        "environment": environment(), "finished": utc_now(),
        # inputs and outputs WITH md5 (PIPELINE_CORRECTNESS_POLICY.md Sec. 2): a path is not a version
        "inputs": [file_record(args.input, role="bc1_mapped_to_oligos"),
                   file_record(args.sps_mapping, role="sps_mapping")],
        "outputs": [file_record(args.output_reference, rows=int(len(ref))),
                    file_record(args.output_summary)],
        "input": args.input, "input_rows": int(n_input_rows), "input_bc1": int(n_input_bc1),
        "output_reference": args.output_reference, "rows": int(len(ref)),
        "unique_bc1": int(ref["bc1"].nunique()), "columns": list(ref.columns),
        "bc1_orientation": BC1_ORIENTATION,
        "bc1_orientation_meaning": ("bc1 string = the 20 nt read 5'->3' immediately downstream of "
                                    "CAGGTCATGCTC on whichever read strand carried it (the CATGCTC-[N20] "
                                    "read structure). Genome(+) holds the reverse complement. Join "
                                    "screen/WGS barcodes in THIS orientation or the join returns zero rows."),
        "support_tier_bc1": ref["support_tier"].value_counts().to_dict() if "support_tier" in ref.columns else None,
        "maps_to_n_variants_gt1": int((ref["maps_to_n_variants"] > 1).sum()) if "maps_to_n_variants" in ref.columns else None,
    }
    with open(f"{args.output_reference}.provenance.json", "w") as fh:
        json.dump(prov, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
