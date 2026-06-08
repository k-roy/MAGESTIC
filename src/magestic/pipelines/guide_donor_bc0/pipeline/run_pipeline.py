#!/usr/bin/env python3
"""
Guide-Donor-BC0 Pipeline Runner

Main entry point for running the guide-donor-bc0 plasmid library analysis pipeline.
Supports both local execution and SLURM jobs with Sherlock best practices.

Usage:
    # Basic usage
    magestic-gdb-pipeline --project-dir /path/to/project

    # With scratch (for SLURM jobs)
    magestic-gdb-pipeline --project-dir /path/to/project --use-scratch

    # Specific steps
    magestic-gdb-pipeline --project-dir /path/to/project --steps combine,match,purity

Author: Kevin R. Roy
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path
from typing import List, Optional

from ..config import (
    GuideDonorBC0Config,
    get_scratch_workdir,
)

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Pipeline Steps
# =============================================================================

AVAILABLE_STEPS = [
    "trim",       # Step 01: Trim, merge, collapse FASTQ
    "combine",    # Step 02: Combine collapsed FASTQ files
    "match",      # Step 03: Map to designed oligos
    "purity",     # Step 04: Calculate BC0 purity
    "represent",  # Step 05: Assess library representation
    "plots",      # Step 06: Plot library representation
]


def run_trim_merge_collapse(config: GuideDonorBC0Config) -> bool:
    """Run Step 01: Trim, merge, collapse FASTQ."""
    logger.info("Step 01: Trim, merge, collapse FASTQ...")
    logger.info("  (This step uses BBDuk/BBMerge - run via shell script)")
    logger.info("  Template: magestic/pipelines/guide_donor_bc0/snakemake/01_trim_merge_collapse.sh")
    logger.info("Step 01: Complete (run manually)")
    return True


def run_combine_fastq(config: GuideDonorBC0Config) -> bool:
    """Run Step 02: Parse collapsed FASTQ files into guide-donor-BC0 count tables."""
    logger.info("Step 02: Combining collapsed FASTQ files...")

    from ..core.fastq_processing import load_sample_keyfile, process_library

    # Resolve sample keyfile
    sample_keyfile = config.sample_keyfile
    if sample_keyfile is None:
        # Try default location
        sample_keyfile = (
            config.project_dir / "keyfiles" / "guide_donor_bc0"
            / "step_1_sample_key_unified.tsv"
        )
    if not sample_keyfile.exists():
        logger.error(f"  Sample keyfile not found: {sample_keyfile}")
        logger.error("  Set config.sample_keyfile or place keyfile at the expected path.")
        return False

    if not config.collapsed_fastq_dir or not config.collapsed_fastq_dir.exists():
        logger.error(f"  Collapsed FASTQ directory not configured or missing: "
                     f"{config.collapsed_fastq_dir}")
        return False

    counts_dir = config.get_working_dir("guide_donor_bc0")
    counts_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"  Keyfile: {sample_keyfile}")
    logger.info(f"  FASTQ dir: {config.collapsed_fastq_dir}")
    logger.info(f"  Counts dir: {counts_dir}")

    library_to_samples, all_dates = load_sample_keyfile(sample_keyfile)
    libraries = [lib for lib in library_to_samples if lib != "empty"]
    logger.info(f"  Libraries to process: {', '.join(libraries)}")

    for library_ID in libraries:
        stats = process_library(
            library_ID=library_ID,
            samples=library_to_samples[library_ID],
            collapsed_fastq_dir=config.collapsed_fastq_dir,
            counts_dir=counts_dir,
            all_dates=all_dates,
        )
        logger.info(
            f"  {library_ID}: {stats['unique_sequences']:,} unique sequences, "
            f"{stats['total_counts']:,} total counts, "
            f"{stats['files_processed']} files"
        )
        if stats["missing_files"]:
            logger.warning(f"    {len(stats['missing_files'])} FASTQ files not found")

    logger.info("Step 02: Complete")
    return True


def run_oligo_matching(config: GuideDonorBC0Config) -> bool:
    """Run Step 03: Map sequences to designed oligos.

    Two oligo-design sources:
      * raw 200mer order (default): guide/donor sliced at fixed SpCas9/SpG offsets
        (config.guide_slice / donor_slice). Preserved for bit-for-bit reproduction of
        legacy outputs; the fixed slice mis-extracts LbCas12a guide/donor.
      * harmonized table (config.use_harmonized_table): guide/donor already extracted
        per nuclease -- required for LbCas12a, recommended for the dense 2024 SpG design.
        Optional per-library scoping (config.library_scope) restricts each library's
        lookup to its nuclease/sublibrary, cutting cross-sublibrary ambiguity in the
        uniqueness-based rescue tiers.
    """
    logger.info("Step 03: Mapping to designed oligos...")

    from ..core.oligo_matching import (
        load_oligo_pool,
        load_oligo_pool_harmonized,
        build_lookup_dictionaries,
        process_counts_file,
    )

    working_dir = config.get_working_dir("guide_donor_bc0")
    working_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"  Counts directory: {working_dir}")

    def _build_lookups(v_libraries, nuclease_types):
        if config.use_harmonized_table:
            df = load_oligo_pool_harmonized(
                v_libraries=v_libraries,
                nuclease_types=nuclease_types,
            )
        else:
            df = load_oligo_pool(
                config.oligo_design_file,
                guide_slice=config.guide_slice,
                donor_slice=config.donor_slice,
            )
        lk = build_lookup_dictionaries(df)
        logger.info(
            f"    {len(df):,} oligos | exact guide+donor: {len(lk['guide_donor_to_oligo']):,}  "
            f"guides: {len(lk['guide_to_oligos']):,}  donors: {len(lk['perfect_donor_to_oligos']):,}  "
            f"partial: {len(lk['partial_donor_to_oligos']):,}"
        )
        return lk

    if config.use_harmonized_table:
        logger.info("  Oligo source: harmonized QTL table (per-nuclease guide/donor)")
    else:
        if not config.oligo_design_file or not config.oligo_design_file.exists():
            logger.error("  No oligo design file configured or file not found.")
            return False
        logger.info(f"  Oligo source: raw 200mer design {config.oligo_design_file.name}")

    # Raw path: one global lookup. Harmonized path: per-library (scoped) lookups.
    global_lookups = None
    if not config.use_harmonized_table:
        logger.info("  Building lookup dictionaries...")
        global_lookups = _build_lookups(None, None)

    counts_suffix = "_guide_donor_bc0_counts.tsv"
    # Exclude this step's own outputs: "*_matched_guide_donor_bc0_counts.tsv"
    # also ends in counts_suffix, so a naive glob re-ingests them on re-run
    # (creating "_matched_matched_" junk and doubling runtime).
    matched_suffix = "_matched_guide_donor_bc0_counts.tsv"
    input_files = sorted(
        f for f in working_dir.glob(f"*{counts_suffix}")
        if not f.name.endswith(matched_suffix)
    )
    if not input_files:
        logger.warning(
            f"  No files matching *{counts_suffix} found in {working_dir}. "
            "Run the 'combine' step first."
        )
        return False

    logger.info(f"  Found {len(input_files)} counts file(s)")
    all_ok = True
    for input_file in input_files:
        library_id = input_file.name.replace(counts_suffix, "")
        output_file = working_dir / f"{library_id}_matched_guide_donor_bc0_counts.tsv"
        logger.info(f"  Matching {library_id} -> {output_file.name}")
        try:
            if config.use_harmonized_table:
                scope = config.library_scope.get(library_id, {})
                lookups = _build_lookups(
                    scope.get("v_libraries", config.harmonized_v_libraries),
                    scope.get("nuclease_types", config.harmonized_nuclease_types),
                )
            else:
                lookups = global_lookups
            stats = process_counts_file(input_file, output_file, lookups)
            total = stats["total"]
            matched = stats["matched"]
            perfect = stats["perfect_match"]
            match_pct = 100 * matched / total if total else 0
            perfect_pct = 100 * perfect / total if total else 0
            logger.info(
                f"    {total:,} sequences, {matched:,} matched ({match_pct:.1f}%), "
                f"{perfect:,} perfect ({perfect_pct:.1f}%)"
            )
        except Exception:
            logger.exception(f"  Error matching {library_id}")
            all_ok = False

    logger.info("Step 03: Oligo matching complete")
    return all_ok


def run_bc0_purity(config: GuideDonorBC0Config) -> bool:
    """Run Step 04: Calculate BC0 purity for each matched-counts file."""
    logger.info("Step 04: Calculating BC0 purity...")

    from ..core.bc0_purity import (
        PurityConfig,
        find_input_files,
        calculate_bc0_purity,
        INPUT_SUFFIX,
        OUTPUT_SUFFIX,
    )

    counts_dir = config.get_working_dir("guide_donor_bc0")
    purity_dir = config.output_dir / "guide_donor_bc0_purity_filtered"
    purity_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"  Input dir:  {counts_dir}")
    logger.info(f"  Output dir: {purity_dir}")

    purity_config = PurityConfig()
    input_files = find_input_files(counts_dir, INPUT_SUFFIX)

    if not input_files:
        logger.warning(f"  No files matching *{INPUT_SUFFIX} found in {counts_dir}")
        logger.warning("  Run the 'match' step first.")
        return False

    logger.info(f"  Found {len(input_files)} input files")
    all_stats = []
    for input_file in input_files:
        library_id = input_file.stem.replace(INPUT_SUFFIX.replace(".tsv", ""), "")
        output_file = purity_dir / f"{library_id}{OUTPUT_SUFFIX}"
        stats = calculate_bc0_purity(input_file, output_file, purity_config)
        stats["library_id"] = library_id
        all_stats.append(stats)
        if stats["mean_purity"] is not None:
            logger.info(
                f"  {library_id}: {stats['unique_guide_bc0']:,} guide-BC0 pairs, "
                f"mean purity {stats['mean_purity']:.3f}"
            )

    # Write summary
    if all_stats:
        import pandas as pd
        pd.DataFrame(all_stats).to_csv(
            purity_dir / "purity_summary.tsv", sep="\t", index=False
        )

    logger.info("Step 04: Complete")
    return True


def _canonical_twist_designs():
    """Canonical Twist 200mer design files (with oligo_seq, hence SPS) for the
    represent/plots steps in harmonized mode. Resolved relative to the harmonized
    oligo table location so it tracks the same common/ tree the matcher uses."""
    try:
        from magestic.utils.oligo_annotations import HARMONIZED_OLIGO_TABLE
        common = Path(HARMONIZED_OLIGO_TABLE).parent.parent  # annotation_files -> common
    except Exception:
        return []
    candidates = [
        common / "oligo_designs" / "20210422_Twist_200mer.tsv",                    # 2021: SpCas9 + (imp)LbCas12a
        common / "oligo_designs" / "20240411_Twist_200mer_oligo_array_order.tsv",   # 2024: SpG
    ]
    return [c for c in candidates if c.exists()]


def run_library_representation(config: GuideDonorBC0Config) -> bool:
    """Run Step 05: Aggregate to oligo-name level and identify missing oligos."""
    logger.info("Step 05: Assessing library representation...")

    from ..core.library_representation import run_representation_analysis

    purity_dir = config.output_dir / "guide_donor_bc0_purity_filtered"
    summary_dir = purity_dir / "summary_plots"

    if not purity_dir.exists():
        logger.error(f"  Purity directory not found: {purity_dir}")
        logger.error("  Run the 'purity' step first.")
        return False

    # Resolve oligo design file(s)
    oligo_pool_files = []
    if config.oligo_design_file and config.oligo_design_file.exists():
        oligo_pool_files = [config.oligo_design_file]
    else:
        # Search keyfiles directory
        kf_dir = config.project_dir / "keyfiles" / "guide_donor_bc0"
        oligo_pool_files = list(kf_dir.glob("*Twist*.tsv"))

    if not oligo_pool_files and getattr(config, "use_harmonized_table", False):
        oligo_pool_files = _canonical_twist_designs()
        if oligo_pool_files:
            logger.info(
                "  Harmonized mode: using canonical Twist designs "
                f"{[f.name for f in oligo_pool_files]} "
                "(missing-oligo fill skipped unless library_sps_map is set)"
            )
    if not oligo_pool_files:
        logger.error("  No oligo design files found. Configure config.oligo_design_file.")
        return False

    library_sps_map = config.library_sps_map or None
    logger.info(f"  Purity dir: {purity_dir}")
    logger.info(f"  Summary dir: {summary_dir}")
    logger.info(f"  Oligo pools: {[f.name for f in oligo_pool_files]}")
    if library_sps_map:
        logger.info(f"  Library SPS map: {list(library_sps_map.keys())}")
    else:
        logger.info("  No library_sps_map set; missing-oligo filling skipped")

    oligo_df = run_representation_analysis(
        purity_dir=purity_dir,
        oligo_pool_files=oligo_pool_files,
        summary_dir=summary_dir,
        library_sps_map=library_sps_map,
    )
    logger.info(f"  Libraries: {list(oligo_df['library'].unique())}")
    logger.info(f"  Total oligo-library rows: {len(oligo_df):,}")
    logger.info("Step 05: Complete")
    return True


def run_plots(config: GuideDonorBC0Config) -> bool:
    """Run Step 06: Generate library QC plots."""
    logger.info("Step 06: Creating plots...")

    from ..core.library_representation import run_library_plots

    purity_dir = config.output_dir / "guide_donor_bc0_purity_filtered"
    plots_dir = purity_dir / "plots"

    if not purity_dir.exists():
        logger.error(f"  Purity directory not found: {purity_dir}")
        logger.error("  Run the 'purity' step first.")
        return False

    # Resolve oligo design file(s)
    oligo_pool_files = []
    if config.oligo_design_file and config.oligo_design_file.exists():
        oligo_pool_files = [config.oligo_design_file]
    else:
        kf_dir = config.project_dir / "keyfiles" / "guide_donor_bc0"
        oligo_pool_files = list(kf_dir.glob("*Twist*.tsv"))

    if not oligo_pool_files and getattr(config, "use_harmonized_table", False):
        oligo_pool_files = _canonical_twist_designs()
        if oligo_pool_files:
            logger.info(
                "  Harmonized mode: using canonical Twist designs "
                f"{[f.name for f in oligo_pool_files]} "
                "(missing-oligo fill skipped unless library_sps_map is set)"
            )
    if not oligo_pool_files:
        logger.error("  No oligo design files found. Configure config.oligo_design_file.")
        return False

    library_sps_map = config.library_sps_map or None
    library_order = sorted(library_sps_map.keys()) if library_sps_map else None

    logger.info(f"  Plots dir: {plots_dir}")
    run_library_plots(
        purity_dir=purity_dir,
        oligo_pool_files=oligo_pool_files,
        plots_dir=plots_dir,
        library_sps_map=library_sps_map,
        library_order=library_order,
    )
    logger.info("Step 06: Complete")
    return True


STEP_FUNCTIONS = {
    "trim": run_trim_merge_collapse,
    "combine": run_combine_fastq,
    "match": run_oligo_matching,
    "purity": run_bc0_purity,
    "represent": run_library_representation,
    "plots": run_plots,
}


# =============================================================================
# Data Staging (Sherlock Best Practices)
# =============================================================================

def stage_inputs_to_scratch(config: GuideDonorBC0Config) -> None:
    """Stage input files from Oak to Scratch."""
    if not config.use_scratch or config.scratch_dir is None:
        return

    logger.info("Staging inputs from Oak to Scratch...")

    scratch_input = config.scratch_dir / "input"
    scratch_input.mkdir(parents=True, exist_ok=True)

    # Stage collapsed FASTQ files
    if config.collapsed_fastq_dir and config.collapsed_fastq_dir.exists():
        dst = scratch_input / "collapsed_fastq"
        shutil.copytree(config.collapsed_fastq_dir, dst, dirs_exist_ok=True)
        config.collapsed_fastq_dir = dst
        logger.info(f"  Staged collapsed FASTQ to: {dst}")

    # Stage oligo design file
    if config.oligo_design_file and config.oligo_design_file.exists():
        dst = scratch_input / config.oligo_design_file.name
        shutil.copy2(config.oligo_design_file, dst)
        config.oligo_design_file = dst
        logger.info(f"  Staged oligo design: {dst.name}")

    logger.info(f"Inputs staged to: {scratch_input}")


def stage_outputs_to_oak(config: GuideDonorBC0Config, oak_output_dir: Path) -> None:
    """Stage output files from Scratch back to Oak."""
    if not config.use_scratch or config.scratch_dir is None:
        return

    logger.info("Staging outputs from Scratch to Oak...")

    # Copy results back to Oak
    dirs_to_copy = ["guide_donor_bc0", "plots"]

    for subdir in dirs_to_copy:
        src = config.scratch_dir / subdir
        if src.exists():
            dst = oak_output_dir / subdir
            dst.mkdir(parents=True, exist_ok=True)
            shutil.copytree(src, dst, dirs_exist_ok=True)
            logger.info(f"  Copied: {subdir}")

    logger.info(f"Outputs staged to: {oak_output_dir}")


# =============================================================================
# Main Pipeline
# =============================================================================

def run_pipeline(
    config: GuideDonorBC0Config,
    steps: Optional[List[str]] = None,
) -> bool:
    """
    Run the guide-donor-bc0 pipeline.

    Args:
        config: Pipeline configuration
        steps: List of steps to run (default: all steps)

    Returns:
        True if all steps completed successfully
    """
    if steps is None:
        steps = AVAILABLE_STEPS

    logger.info("=" * 70)
    logger.info("GUIDE-DONOR-BC0 PIPELINE")
    logger.info("=" * 70)
    logger.info(f"Project: {config.project_name}")
    logger.info(f"Use scratch: {config.use_scratch}")
    if config.use_scratch:
        logger.info(f"Scratch dir: {config.scratch_dir}")
    logger.info(f"Steps to run: {', '.join(steps)}")
    logger.info("=" * 70)

    # Validate configuration
    errors = config.validate()
    if errors:
        for error in errors:
            logger.error(f"Config error: {error}")
        return False

    # Save original output dir for staging back
    oak_output_dir = config.output_dir

    # Stage inputs if using scratch
    if config.use_scratch:
        stage_inputs_to_scratch(config)

    # Ensure output directories exist
    config.ensure_directories()

    # Run steps
    success = True
    for step in steps:
        if step not in STEP_FUNCTIONS:
            logger.warning(f"Unknown step: {step}")
            continue

        try:
            step_success = STEP_FUNCTIONS[step](config)
            if not step_success:
                logger.error(f"Step {step} failed")
                success = False
                break
        except Exception as e:
            logger.exception(f"Error in step {step}: {e}")
            success = False
            break

    # Stage outputs if using scratch
    if config.use_scratch and success:
        stage_outputs_to_oak(config, oak_output_dir)

    if success:
        logger.info("=" * 70)
        logger.info("PIPELINE COMPLETED SUCCESSFULLY")
        logger.info("=" * 70)
    else:
        logger.error("=" * 70)
        logger.error("PIPELINE FAILED")
        logger.error("=" * 70)

    return success


def main():
    parser = argparse.ArgumentParser(
        prog="magestic-gdb-pipeline",
        description="Guide-Donor-BC0 Pipeline Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline
  magestic-gdb-pipeline --project-dir /path/to/project

  # Run with scratch (for SLURM jobs)
  magestic-gdb-pipeline --project-dir /path/to/project --use-scratch

  # Run specific steps
  magestic-gdb-pipeline --project-dir /path/to/project --steps combine,match,purity
        """
    )

    parser.add_argument(
        "--project-dir",
        type=Path,
        required=True,
        help="Path to project directory"
    )
    parser.add_argument(
        "--use-scratch",
        action="store_true",
        help="Use $SCRATCH for job I/O (recommended for SLURM jobs)"
    )
    parser.add_argument(
        "--steps",
        type=str,
        help=f"Comma-separated list of steps. Available: {','.join(AVAILABLE_STEPS)}"
    )
    parser.add_argument(
        "--oligo-design-file",
        type=Path,
        default=None,
        help="Path to oligo design TSV (overrides auto-detection)"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode"
    )
    parser.add_argument(
        "--use-harmonized-table",
        action="store_true",
        help="Match against the harmonized QTL oligo table (per-nuclease guide/donor) "
             "instead of slicing the raw 200mer. Required for LbCas12a libraries."
    )
    parser.add_argument(
        "--nuclease",
        type=str,
        default=None,
        help="Comma-separated nuclease/PAM types to scope the harmonized lookup "
             "(e.g. SpG_NGNG,SpCas9_NGG,LbCas12a_TTTV). Requires --use-harmonized-table."
    )
    parser.add_argument(
        "--v-libraries",
        type=str,
        default=None,
        help="Comma-separated V libraries to scope the harmonized lookup. "
             "Requires --use-harmonized-table."
    )
    parser.add_argument(
        "--library-scope-json",
        type=Path,
        default=None,
        help="Optional JSON {library_ID: {nuclease_types:[...], v_libraries:[...]}} "
             "for per-library harmonized scoping (mixed-nuclease pools)."
    )

    args = parser.parse_args()

    # Create configuration
    config = GuideDonorBC0Config.from_project_dir(args.project_dir)

    # Override oligo design file if specified
    if args.oligo_design_file:
        config.oligo_design_file = args.oligo_design_file

    # Harmonized-table matching (gdb LbCas12a + dense-design fix)
    if args.use_harmonized_table:
        config.use_harmonized_table = True
    if args.nuclease:
        config.harmonized_nuclease_types = [x.strip() for x in args.nuclease.split(",") if x.strip()]
    if args.v_libraries:
        config.harmonized_v_libraries = [x.strip() for x in args.v_libraries.split(",") if x.strip()]
    if args.library_scope_json:
        import json
        with open(args.library_scope_json) as _fh:
            config.library_scope = json.load(_fh)

    # Apply scratch setting
    if args.use_scratch:
        config.use_scratch = True
        config.scratch_dir = get_scratch_workdir(config.project_name)

    # Parse steps
    steps = None
    if args.steps:
        steps = [s.strip() for s in args.steps.split(",")]

    # Enable debug logging
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Run pipeline
    success = run_pipeline(config, steps)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
