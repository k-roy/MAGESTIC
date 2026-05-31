#!/usr/bin/env python3
"""
Add trafo_ID to existing bc1 reference tables via join.

This script updates bc1_reference_table.tsv files by joining trafo_ID from:
1. Existing sublibrary assignment tables (yL437/yL442)
2. Existing bc1_to_trafo_mapping + trafo_key (yL406)

No reprocessing required - just joins existing mapping data.

Usage:
    python 09_add_trafo_id_to_bc1_tables.py --project yL437
    python 09_add_trafo_id_to_bc1_tables.py --project yL406
    python 09_add_trafo_id_to_bc1_tables.py --all

Author: Kevin R. Roy
Date: 2026-03-23
"""

import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Tuple
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Base directories
BASE_DIR = Path("/path/to")
QTL_DIR = BASE_DIR / "projects/QTL"
MANIFEST_PATH = QTL_DIR / "QTL_PROJECT_MANIFEST.tsv"


def get_project_dir(project_name: str) -> Optional[Path]:
    """Find project directory by name pattern."""
    patterns = [
        f"*{project_name}*",
        f"*{project_name.replace('yL', '')}*",
    ]
    for pattern in patterns:
        matches = list(QTL_DIR.glob(pattern))
        if matches:
            return matches[0]
    return None


def process_yL437_yL442():
    """
    Process yL437/yL442 project using existing sublibrary assignment.

    Source: *_with_sublibrary_assignment.tsv already has trafo_ID
    """
    project_dir = get_project_dir("yL437")
    if project_dir is None:
        logger.error("yL437/yL442 project not found")
        return

    logger.info(f"Processing {project_dir.name}")
    logger.info("Source: *_with_sublibrary_assignment.tsv")

    # Load sublibrary assignment table
    sublibrary_file = project_dir / "processed_data/bc1_donor_bc0/data_tables/20250912_SpG_QTL_bc1_reference_table_min_2_reads_with_sublibrary_assignment.tsv"

    if not sublibrary_file.exists():
        logger.error(f"Sublibrary assignment file not found: {sublibrary_file}")
        return

    sublibrary_df = pd.read_csv(sublibrary_file, sep="\t", usecols=["bc1", "trafo_ID", "sublibrary"])
    logger.info(f"Loaded {len(sublibrary_df):,} bc1s with trafo_ID from sublibrary assignment")

    # Build lookup
    bc1_to_trafo = dict(zip(sublibrary_df["bc1"], sublibrary_df["trafo_ID"]))
    bc1_to_sublibrary = dict(zip(sublibrary_df["bc1"], sublibrary_df["sublibrary"]))

    # Load main bc1 reference table
    ref_file = project_dir / "processed_data/bc1_donor_bc0/bc1_reference/bc1_reference_table.tsv"
    ref_df = pd.read_csv(ref_file, sep="\t")
    logger.info(f"Main bc1 reference table: {len(ref_df):,} bc1s")

    # Add trafo_ID column
    ref_df["trafo_ID"] = ref_df["bc1"].map(bc1_to_trafo)
    ref_df["sublibrary"] = ref_df["bc1"].map(bc1_to_sublibrary)

    # Report coverage
    coverage = ref_df["trafo_ID"].notna().sum()
    logger.info(f"Coverage: {coverage:,} / {len(ref_df):,} ({100*coverage/len(ref_df):.1f}%)")

    if coverage > 0:
        logger.info(f"\ntrafo_ID distribution:")
        print(ref_df["trafo_ID"].value_counts())

    # Save updated table
    output_file = ref_file.parent / "bc1_reference_table_with_trafo_id.tsv"
    ref_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Saved: {output_file}")

    # Also check coverage for min_2_reads table
    min2_file = project_dir / "processed_data/bc1_donor_bc0/data_tables/20250912_SpG_QTL_bc1_reference_table_min_2_reads.tsv"
    if min2_file.exists():
        min2_df = pd.read_csv(min2_file, sep="\t")
        min2_df["trafo_ID"] = min2_df["bc1"].map(bc1_to_trafo)
        min2_coverage = min2_df["trafo_ID"].notna().sum()
        logger.info(f"\nCoverage for bc1s with ≥2 reads: {min2_coverage:,} / {len(min2_df):,} ({100*min2_coverage/len(min2_df):.1f}%)")

        # Save updated min2 table
        min2_output = min2_file.parent / f"{min2_file.stem}_with_trafo_id.tsv"
        min2_df["sublibrary"] = min2_df["bc1"].map(bc1_to_sublibrary)
        min2_df.to_csv(min2_output, sep="\t", index=False)
        logger.info(f"Saved: {min2_output}")


def process_yL406():
    """
    Process yL406 project using existing bc1_to_trafo_mapping + trafo_key.

    Source files:
    - bc1_to_trafo_mapping.tsv: bc1 → yeast_sublibrary (trafo_N)
    - trafo_key.tsv: trafo_N → yL, nuclease_PAM
    """
    project_dir = get_project_dir("yL406")
    if project_dir is None:
        logger.error("yL406 project not found")
        return

    logger.info(f"Processing {project_dir.name}")
    logger.info("Source: bc1_to_trafo_mapping.tsv + trafo_key.tsv")

    # Load bc1 to trafo mapping (from bc1_from_screen)
    bc1_trafo_file = project_dir / "processed_data/bc1_from_screen/bc1_to_trafo_mapping.tsv"
    if not bc1_trafo_file.exists():
        logger.error(f"bc1_to_trafo_mapping.tsv not found: {bc1_trafo_file}")
        return

    bc1_trafo_df = pd.read_csv(bc1_trafo_file, sep="\t")
    logger.info(f"Loaded {len(bc1_trafo_df):,} bc1s from bc1_to_trafo_mapping.tsv")

    # Build bc1 → trafo lookup (yeast_sublibrary column has trafo_N)
    bc1_to_trafo = dict(zip(bc1_trafo_df["bc1"], bc1_trafo_df["yeast_sublibrary"]))
    bc1_to_purity = dict(zip(bc1_trafo_df["bc1"], bc1_trafo_df["purity"]))

    logger.info(f"\ntrafo distribution in bc1_to_trafo_mapping:")
    print(bc1_trafo_df["yeast_sublibrary"].value_counts().head(15))

    # Load trafo_key for trafo → yL, nuclease_PAM mapping
    trafo_key_file = project_dir / "keyfiles/bc1_donor_bc0/trafo_key.tsv"
    if not trafo_key_file.exists():
        logger.error(f"trafo_key.tsv not found: {trafo_key_file}")
        return

    trafo_key_df = pd.read_csv(trafo_key_file, sep="\t")
    logger.info(f"\nLoaded trafo_key with {len(trafo_key_df)} entries")

    # Build trafo → metadata mapping (use first replica)
    trafo_to_yL = {}
    trafo_to_nuclease_PAM = {}
    trafo_to_nuclease = {}
    for _, row in trafo_key_df.iterrows():
        trafo = row["trafo_ID"]
        if trafo not in trafo_to_yL:  # First occurrence (rep_1)
            trafo_to_yL[trafo] = row["yL"]
            trafo_to_nuclease_PAM[trafo] = row["nuclease_PAM"]
            trafo_to_nuclease[trafo] = row["nuclease"]

    logger.info(f"\ntrafo → yL mapping:")
    for trafo, yL in sorted(trafo_to_yL.items()):
        logger.info(f"  {trafo} → {yL} ({trafo_to_nuclease_PAM.get(trafo, 'N/A')})")

    # Load main bc1 reference table
    ref_file = project_dir / "processed_data/bc1_donor_bc0/bc1_reference/bc1_reference_table.tsv"
    ref_df = pd.read_csv(ref_file, sep="\t")
    logger.info(f"\nMain bc1 reference table: {len(ref_df):,} bc1s")

    # Add trafo columns
    ref_df["trafo_ID"] = ref_df["bc1"].map(bc1_to_trafo)
    ref_df["trafo_purity"] = ref_df["bc1"].map(bc1_to_purity)
    ref_df["yL_from_trafo"] = ref_df["trafo_ID"].map(trafo_to_yL)
    ref_df["nuclease_from_trafo"] = ref_df["trafo_ID"].map(trafo_to_nuclease)

    # Report coverage
    coverage = ref_df["trafo_ID"].notna().sum()
    logger.info(f"\nCoverage: {coverage:,} / {len(ref_df):,} ({100*coverage/len(ref_df):.1f}%)")

    if coverage > 0:
        logger.info(f"\ntrafo_ID distribution in updated table:")
        print(ref_df["trafo_ID"].value_counts().head(15))

        logger.info(f"\nyL distribution:")
        print(ref_df["yL_from_trafo"].value_counts().head(15))

    # Save updated table
    output_file = ref_file.parent / "bc1_reference_table_with_trafo_id.tsv"
    ref_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"\nSaved: {output_file}")

    # Also update the main data_tables reference table
    main_ref_file = project_dir / "processed_data/bc1_donor_bc0/data_tables/MAGESTIC_SpG_QTL_library_yL406_bc1_reference_table.tsv"
    if main_ref_file.exists():
        main_df = pd.read_csv(main_ref_file, sep="\t")
        main_df["trafo_ID"] = main_df["bc1"].map(bc1_to_trafo)
        main_df["trafo_purity"] = main_df["bc1"].map(bc1_to_purity)
        main_df["yL_from_trafo"] = main_df["trafo_ID"].map(trafo_to_yL)
        main_df["nuclease_from_trafo"] = main_df["trafo_ID"].map(trafo_to_nuclease)

        main_output = main_ref_file.parent / f"{main_ref_file.stem}_with_trafo_id.tsv"
        main_df.to_csv(main_output, sep="\t", index=False)

        main_coverage = main_df["trafo_ID"].notna().sum()
        logger.info(f"\nMain reference table: {main_coverage:,} / {len(main_df):,} ({100*main_coverage/len(main_df):.1f}%)")
        logger.info(f"Saved: {main_output}")


def main():
    parser = argparse.ArgumentParser(
        description="Add trafo_ID to bc1 reference tables via join",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process yL437/yL442 (uses existing sublibrary assignment)
    python 09_add_trafo_id_to_bc1_tables.py --project yL437

    # Process yL406 (uses bc1_to_trafo_mapping + trafo_key)
    python 09_add_trafo_id_to_bc1_tables.py --project yL406

    # Process all projects
    python 09_add_trafo_id_to_bc1_tables.py --all

    # List available projects
    python 09_add_trafo_id_to_bc1_tables.py --list
        """
    )

    parser.add_argument("--project", help="Project to process (yL437, yL406, etc.)")
    parser.add_argument("--all", action="store_true", help="Process all QTL projects")
    parser.add_argument("--list", action="store_true", help="List available projects")

    args = parser.parse_args()

    if args.list:
        print("\nAvailable QTL projects:")
        for project_dir in sorted(QTL_DIR.glob("*")):
            if project_dir.is_dir():
                has_bc1 = bool(list(project_dir.glob("processed_data/bc1_donor_bc0/*")))
                status = "✓ has bc1_donor_bc0 data" if has_bc1 else "  (no bc1_donor_bc0 data)"
                print(f"  {project_dir.name} {status}")
        return

    if args.all:
        logger.info("Processing all projects...")

        logger.info("\n" + "=" * 60)
        process_yL437_yL442()

        logger.info("\n" + "=" * 60)
        process_yL406()

        logger.info("\n" + "=" * 60)
        logger.info("SUMMARY")
        logger.info("=" * 60)
        logger.info("""
yL437/yL442:
  - Source: *_with_sublibrary_assignment.tsv
  - trafo_ID: trafo_1, trafo_2, trafo_3, trafo_4
  - sublibrary: SpG_NGNG, SpG_NGNH, SpG_NGG, SpCas9_NGG

yL406:
  - Source: bc1_to_trafo_mapping.tsv + trafo_key.tsv
  - trafo_ID: trafo_1 through trafo_13
  - yL_from_trafo: yL386-yL405
  - nuclease_from_trafo: SpG, eSpCas9_SpG, Sniper2L_SpG, WT_SpCas9
        """)
        return

    if args.project:
        project = args.project.lower()

        if "437" in project or "442" in project:
            process_yL437_yL442()
        elif "406" in project:
            process_yL406()
        else:
            logger.error(f"Unknown project: {args.project}")
            logger.info("Supported projects: yL437/yL442, yL406")
        return

    parser.print_help()


if __name__ == "__main__":
    main()
