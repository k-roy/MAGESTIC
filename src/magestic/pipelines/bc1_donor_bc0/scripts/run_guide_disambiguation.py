#!/usr/bin/env python3
"""
Run guide-based disambiguation for multi-PAM oligo assignments.

This script runs the guide disambiguation on a QTL screen project's reference table.

Usage:
    python run_guide_disambiguation.py --project /path/to/project
    python run_guide_disambiguation.py --input ref_table.tsv --output ref_table_fixed.tsv

Author: Kevin Roy / Claude
Date: 2026-02-04
"""

import sys
import logging
from pathlib import Path
from datetime import datetime

BASE_DIR = Path("/path/to")

# Repointed 2026-06-02 off the common/scripts/bc1_donor_bc0_pipeline fork to the
# canonical magestic package (the fork is a Phase-4 delete candidate).
from magestic.pipelines.bc1_donor_bc0.core.guide_disambiguation import (
    run_guide_disambiguation,
    load_oligo_designs_for_disambiguation,
    reassign_oligos_by_guide,
)

# Default oligo design files
SPG_OLIGO_FILE = BASE_DIR / "scripts_and_keyfiles/guide_donor_design/SpCas9_oligo_pools/keyfiles/20240411_Twist_200mer_oligo_array_order.tsv"
SPCAS9_OLIGO_FILE = BASE_DIR / "common/oligo_designs/20210422_Twist_200mer_SpCas9_oligos_with_controls.tsv"

# Known screen configurations
SCREEN_CONFIGS = {
    "20250912_SpG_QTL_screen": {
        "data_dir": BASE_DIR / "processed_data/by_project/QTL/20250912_SpG_QTL_screen",
        "reference_table": "bc1_donor_bc0/data_tables/20250912_SpG_QTL_screen_bc1_reference_table_with_unmatched_annotations.tsv",
        "oligo_col": "bdb_oligo_name",
        "guide_col": "guide",
        "donor_col": "bdb_donor",
        "perfect_guide_col": "perfect_guide",
    },
    "20240806_yL406_SpG_QTL_screen": {
        "data_dir": BASE_DIR / "processed_data/by_project/QTL/20240806_yL406_SpG_QTL_screen",
        # Will need to find the actual reference table path
        "reference_table": "bc1_donor_bc0/data_tables/bc1_reference_table_with_annotations.tsv",
        "oligo_col": "bdb_oligo_name",
        "guide_col": "guide",
        "donor_col": "bdb_donor",
        "perfect_guide_col": "perfect_guide",
    },
}


def find_reference_table(data_dir: Path) -> Path:
    """Find the most recent reference table in a data directory."""
    patterns = [
        "bc1_donor_bc0/data_tables/*reference_table*with*annotations*.tsv",
        "bc1_donor_bc0/data_tables/*reference_table*.tsv",
        "bc1_donor_bc0/*reference_table*.tsv",
    ]

    for pattern in patterns:
        matches = list(data_dir.glob(pattern))
        if matches:
            # Return most recent
            return max(matches, key=lambda p: p.stat().st_mtime)

    raise FileNotFoundError(f"No reference table found in {data_dir}")


def run_for_screen(screen_name: str, force: bool = False) -> dict:
    """Run guide disambiguation for a known screen."""
    if screen_name not in SCREEN_CONFIGS:
        raise ValueError(f"Unknown screen: {screen_name}. Known screens: {list(SCREEN_CONFIGS.keys())}")

    config = SCREEN_CONFIGS[screen_name]
    data_dir = config["data_dir"]

    # Find reference table
    if "reference_table" in config and (data_dir / config["reference_table"]).exists():
        input_path = data_dir / config["reference_table"]
    else:
        input_path = find_reference_table(data_dir)

    logging.info(f"Found reference table: {input_path}")

    # Create output path
    timestamp = datetime.now().strftime("%Y%m%d")
    output_name = input_path.stem + "_guide_disambiguated.tsv"
    output_path = input_path.parent / output_name

    if output_path.exists() and not force:
        logging.warning(f"Output already exists: {output_path}")
        logging.warning("Use --force to overwrite")
        return {}

    # Run disambiguation
    stats = run_guide_disambiguation(
        reference_table_path=input_path,
        oligo_design_files=[SPG_OLIGO_FILE, SPCAS9_OLIGO_FILE],
        output_path=output_path,
        oligo_col=config.get("oligo_col", "bdb_oligo_name"),
        guide_col=config.get("guide_col", "guide"),
        donor_col=config.get("donor_col", "bdb_donor"),
        perfect_guide_col=config.get("perfect_guide_col", "perfect_guide"),
        max_guide_distance=2,
    )

    return stats


def main():
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    parser = argparse.ArgumentParser(
        description='Run guide-based multi-PAM disambiguation'
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--screen', '-s', choices=list(SCREEN_CONFIGS.keys()),
                       help='Run for a known screen')
    group.add_argument('--input', '-i', help='Input reference table path')

    parser.add_argument('--output', '-o', help='Output path (required with --input)')
    parser.add_argument('--force', '-f', action='store_true',
                        help='Overwrite existing output')
    parser.add_argument('--all-screens', action='store_true',
                        help='Run for all known screens')

    args = parser.parse_args()

    if args.all_screens:
        for screen in SCREEN_CONFIGS.keys():
            logging.info(f"\n{'='*60}")
            logging.info(f"Processing {screen}")
            logging.info('='*60)
            try:
                stats = run_for_screen(screen, force=args.force)
                if stats:
                    logging.info(f"Reassigned {stats.get('reassigned', 0):,} BC1s")
            except Exception as e:
                logging.error(f"Failed to process {screen}: {e}")

    elif args.screen:
        stats = run_for_screen(args.screen, force=args.force)
        if stats:
            print("\n=== Results ===")
            for k, v in stats.items():
                print(f"  {k}: {v}")

    elif args.input:
        if not args.output:
            parser.error("--output required with --input")

        stats = run_guide_disambiguation(
            reference_table_path=Path(args.input),
            oligo_design_files=[SPG_OLIGO_FILE, SPCAS9_OLIGO_FILE],
            output_path=Path(args.output),
            max_guide_distance=2,
        )

        print("\n=== Results ===")
        for k, v in stats.items():
            print(f"  {k}: {v}")


if __name__ == '__main__':
    main()
