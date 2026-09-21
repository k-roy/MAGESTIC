#!/usr/bin/env python3
"""
Main pipeline runner for BC1-Donor-BC0 processing.

This script orchestrates the full pipeline:
1. Load data from available sources (2x150, 2x300, gdb)
2. Merge data sources at the earliest stage
3. Apply unified sliding window matching for oligo assignment
4. Output reference table with bc1 -> oligo mappings

Usage:
    python run_pipeline.py --project-dir /path/to/project [options]

The pipeline automatically detects available data sources and adjusts
processing accordingly. It supports three scenarios:
- Only 2x150 data
- Only 2x300 data
- Both 2x150 and 2x300 data (merged for maximum bc1 recovery)
"""

import argparse
import logging
import sys
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Optional
import pandas as pd
import numpy as np
from tqdm import tqdm

from ..config import PipelineConfig
from ..readers.data_loader import load_all_data_sources
from ..core.merge import merge_data_sources
from ..core.matching import (
    load_and_build_lookup,
    match_dataframe,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='BC1-Donor-BC0 Processing Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run with auto-detection of data sources
    python run_pipeline.py --project-dir /path/to/project

    # Run in debug mode with limited data
    python run_pipeline.py --project-dir /path/to/project --debug

    # Specify output directory
    python run_pipeline.py --project-dir /path/to/project --output-dir /path/to/output
        """
    )

    parser.add_argument(
        '--project-dir', '-p',
        type=Path,
        required=True,
        help='Root directory of the project'
    )

    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=None,
        help='Output directory (default: project_dir/processed_data/bc1_donor_bc0)'
    )

    parser.add_argument(
        '--common-dir', '-c',
        type=Path,
        default=Path('/path/to/common'),
        help='Directory containing common resources (oligo designs, etc.)'
    )

    parser.add_argument(
        '--min-reads',
        type=int,
        default=2,
        help='Minimum reads for a bc1-donor-bc0 to be included (default: 2)'
    )

    parser.add_argument(
        '--min-pcr-replicates',
        type=int,
        default=2,
        help='Minimum PCR replicates for confidence (default: 2)'
    )

    parser.add_argument(
        '--gdb-dir',
        type=Path,
        default=None,
        help='Guide-donor-bc0 table directory, EXPLICIT (overrides auto-detection). '
             'Auto-detection points at <project>/processed_data/guide_donor_bc0, whose '
             'top level holds no tables for some screens (0 files -> formerly a silent '
             'empty gdb) and whose subdirectory may hold SUPERSEDED tables (G3, 27_).'
    )

    parser.add_argument(
        '--no-gdb',
        action='store_true',
        help='Run explicitly WITHOUT guide-donor-bc0 data'
    )

    parser.add_argument(
        '--debug',
        action='store_true',
        help='Run in debug mode with limited data'
    )

    parser.add_argument(
        '--debug-n-samples',
        type=int,
        default=4,
        help='Number of samples to load in debug mode (default: 4)'
    )

    parser.add_argument(
        '--n-jobs',
        type=int,
        default=16,
        help='Number of parallel jobs (default: 16)'
    )

    parser.add_argument(
        '--annotate-unmatched',
        action='store_true',
        help='STEP 6: annotate unmatched donors via 06_annotate_unmatched_donors.py '
             '(uses the MAGESTIC_cs DESIGN reference, not the yKR1033 WGS assembly)'
    )

    parser.add_argument(
        '--integrate-annotations',
        action='store_true',
        help='After STEP 6, run 07 to integrate unmatched annotations into the final table'
    )

    parser.add_argument(
        '--spg-oligo-design', type=Path, default=None,
        help='SpG oligo design TSV for STEP 6 '
             '(default: common/oligo_designs/20240411_Twist_200mer_oligo_array_order.tsv)'
    )

    parser.add_argument(
        '--spcas9-oligo-design', type=Path, default=None,
        help='SpCas9 oligo design TSV for STEP 6 '
             '(default: common/oligo_designs/20210422_Twist_200mer.tsv)'
    )

    parser.add_argument(
        '--reference-genome', type=Path, default=None,
        help='MAGESTIC_cs DESIGN reference FASTA for STEP 6 '
             '(default: common/reference_genomes/MAGESTIC/MAGESTIC_background_strain.fasta)'
    )

    parser.add_argument(
        '--annotation-gff', type=Path, default=None,
        help='MAGESTIC_cs annotation GFF for STEP 6 '
             '(default: common/annotation_files/MAGESTIC_background_strain_annotations.gff)'
    )

    parser.add_argument(
        '--harmonized-oligo-key', type=Path, default=None,
        help='Harmonized oligo key for STEP 6 delta-variant computation '
             '(default: canonical harmonized table)'
    )

    # Explicit data-source overrides (added 2026-09-11). Their ABSENCE is why a 772-symlink shim was built
    # on Oak to fake the one directory layout the detector recognised: there was no way to point the
    # pipeline at a directory. Policy is to fix the tool, not to shim the filesystem.
    parser.add_argument(
        '--parsed-dir', type=Path, default=None,
        help='Directory of parsed *_bc1_donor_bc0.tsv files. Overrides auto-detection; the read layout is '
             'taken from the keyfile inner_product_size_bp rather than from directory names.'
    )
    parser.add_argument(
        '--data-2x150-dir', type=Path, default=None,
        help='Explicit 2x150 parsed directory (overrides auto-detection and --parsed-dir).'
    )
    parser.add_argument(
        '--data-2x300-dir', type=Path, default=None,
        help='Explicit 2x300 parsed directory (overrides auto-detection and --parsed-dir).'
    )

    return parser.parse_args()


def collapse_sequences_by_edit_distance(
    df: pd.DataFrame,
    group_cols: list,
    seq_col: str,
    count_cols: list,
    edit_distance_threshold: int = 3,
    show_progress: bool = True,
    set_cols: Optional[list] = None,
) -> pd.DataFrame:
    """
    Collapse similar sequences within groups by edit distance.

    Args:
        df: Input DataFrame
        group_cols: Columns to group by
        seq_col: Column containing sequences to collapse
        count_cols: Columns to SUM when collapsing (reads)
        edit_distance_threshold: Maximum edit distance for collapsing
        show_progress: Show progress bar
        set_cols: Columns holding *sets of identities* to UNION when collapsing
            (e.g. ``pcr_replicate_id``). These must never be summed --
            see the P2 note below.

    Returns:
        DataFrame with collapsed sequences

    Note (P2 fix, 2026-07-28):
        Replicate SUPPORT must be a **distinct count**, never a sum. If bc1
        variant A is seen in replicates {1,2} and its 1-nt sequencing-error twin
        B in replicate {2}, summing gives 2+1=3 when the true distinct count is
        2 -- inflating a single-tube chimera over a >=3-replicate threshold,
        which is precisely the artifact the recombination filter exists to
        remove. Carrying replicate-ID **sets** and taking the union makes the
        post-collapse count correct by construction.
    """
    from Levenshtein import distance
    set_cols = [c for c in (set_cols or []) if c in df.columns]

    def process_group(group):
        """Collapse sequences within a single group."""
        group = group.copy()
        result_rows = []

        while not group.empty:
            # Find row with highest count (ranking uses read counts only)
            top_idx = group[count_cols].sum(axis=1).idxmax()
            top_seq = group.loc[top_idx, seq_col]

            # Find all sequences within edit distance
            group['_edit_dist'] = group[seq_col].apply(
                lambda x: distance(x, top_seq) if pd.notna(x) else 999
            )
            collapse_mask = group['_edit_dist'] <= edit_distance_threshold

            # Aggregate counts (reads: sum)
            collapsed_counts = {
                col: group.loc[collapse_mask, col].sum()
                for col in count_cols
            }
            # Aggregate identity sets (replicates: UNION, never sum)
            collapsed_sets = {}
            for col in set_cols:
                acc = set()
                for v in group.loc[collapse_mask, col]:
                    if isinstance(v, (set, frozenset, list, tuple)):
                        acc |= set(v)
                    elif pd.notna(v):
                        acc.add(v)
                collapsed_sets[col] = acc

            # Create result row
            result_rows.append({
                **{col: group.iloc[0][col] for col in group_cols},
                seq_col: top_seq,
                **collapsed_counts,
                **collapsed_sets,
            })

            # Remove collapsed rows
            group = group.loc[~collapse_mask].drop(columns=['_edit_dist'])

        return result_rows

    # Process each group
    grouped = df.groupby(group_cols)
    all_results = []

    groups = list(grouped)
    iterator = tqdm(groups, desc=f"Collapsing {seq_col}") if show_progress else groups

    for _, group in iterator:
        all_results.extend(process_group(group))

    return pd.DataFrame(all_results)


def aggregate_and_filter(
    df: pd.DataFrame,
    min_reads: int = 2,
    min_pcr_replicates: int = 2,
    collapse_sequences: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """Aggregate observed link evidence without dropping its denominator scope."""
    from ..core.aggregation import aggregate_link_evidence
    return aggregate_link_evidence(df, min_reads, min_pcr_replicates,
                                   collapse_sequences, show_progress)


def _resolve_annot_paths(args):
    """Resolve STEP-6 reference/design paths, falling back to canonical MAGESTIC_cs
    DESIGN refs. gdb/bdb donors are compared against the design space they were built
    against -- the MAGESTIC_background_strain (MAGESTIC_cs) assembly -- NOT the more
    polished yKR1033-based assembly used for WGS/TRACE."""
    common = Path(args.common_dir)
    defaults = {
        "spg_oligo_design": common / "oligo_designs" / "20240411_Twist_200mer_oligo_array_order.tsv",
        "spcas9_oligo_design": common / "oligo_designs" / "20210422_Twist_200mer.tsv",
        "reference_genome": common / "reference_genomes" / "MAGESTIC" / "MAGESTIC_background_strain.fasta",
        "annotation_gff": common / "annotation_files" / "MAGESTIC_background_strain_annotations.gff",
        "harmonized_oligo_key": common / "annotation_files"
            / "20210422_and_20240411_Bloom_et_al_16_strains_QTL_harmonized_designed_variant_oligos.tsv",
    }
    resolved = {}
    for k, d in defaults.items():
        v = getattr(args, k, None)
        resolved[k] = Path(v) if v else d
    return resolved


def _run_unmatched_annotation(args, config, bc1_reference_table):
    """Flag-gated STEP 6: invoke 06_annotate_unmatched_donors.py (+ optional 07 integrate).

    Runs as a subprocess because 06 uses fork-based parallelism with its own main().
    06's raw output is not the deliverable; the integrated table (07) is, so 07 is
    chained when --integrate-annotations is set. Failures here never clobber the
    STEP-5 reference table, which is already on disk."""
    paths = _resolve_annot_paths(args)
    missing = [str(p) for p in paths.values() if not Path(p).exists()]
    if missing:
        logger.error("STEP 6 skipped -- required reference/design files not found:")
        for m in missing:
            logger.error(f"    {m}")
        logger.error("  (the STEP-5 bc1 reference table is unaffected)")
        return

    script = Path(__file__).resolve().parent.parent / "snakemake" / "06_annotate_unmatched_donors.py"
    out_dir = config.output_dir / "unmatched_annotation"
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(script),
        "--input", str(bc1_reference_table),
        "--output-dir", str(out_dir),
        "--spg-oligo-design", str(paths["spg_oligo_design"]),
        "--spcas9-oligo-design", str(paths["spcas9_oligo_design"]),
        "--reference-genome", str(paths["reference_genome"]),
        "--annotation-gff", str(paths["annotation_gff"]),
        "--harmonized-oligo-key", str(paths["harmonized_oligo_key"]),
        "--n-workers", str(args.n_jobs),
        "--n-chunks", str(args.n_jobs),
    ]
    logger.info("  Running 06_annotate_unmatched_donors.py:")
    logger.info("    " + " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        logger.error(f"STEP 6 (06_annotate_unmatched_donors) failed with exit {result.returncode}")
        return
    logger.info(f"  Unmatched-donor annotations written under: {out_dir}")

    if not getattr(args, "integrate_annotations", False):
        return
    annotations_file = out_dir / "unmatched_donor_annotations.tsv"
    if not annotations_file.exists():
        logger.error(f"STEP 6 integrate skipped -- {annotations_file} not found")
        return
    integrate_script = (Path(__file__).resolve().parent.parent / "snakemake"
                        / "07_integrate_unmatched_annotations.py")
    final_dir = config.output_dir / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    integrated = final_dir / "bc1_reference_table_with_unmatched_annotations.tsv"
    icmd = [
        sys.executable, str(integrate_script),
        "--bc1-reference", str(bc1_reference_table),
        "--unmatched-annotations", str(annotations_file),
        "--output", str(integrated),
    ]
    logger.info("  Running 07_integrate_unmatched_annotations.py:")
    logger.info("    " + " ".join(icmd))
    ires = subprocess.run(icmd)
    if ires.returncode != 0:
        logger.error(f"STEP 6 integrate (07) failed with exit {ires.returncode}")
    else:
        logger.info(f"  Integrated reference table: {integrated}")


def run_pipeline(args):
    """Main pipeline execution."""
    logger.info("=" * 60)
    logger.info("BC1-DONOR-BC0 PROCESSING PIPELINE")
    logger.info("=" * 60)

    # Initialize configuration
    logger.info(f"\nProject directory: {args.project_dir}")
    config = PipelineConfig.from_project_dir(
        args.project_dir, args.common_dir,
        parsed_dir=getattr(args, 'parsed_dir', None),
        data_2x150_dir=getattr(args, 'data_2x150_dir', None),
        data_2x300_dir=getattr(args, 'data_2x300_dir', None),
    )

    if args.output_dir:
        config.output_dir = args.output_dir

    config.min_reads = args.min_reads
    config.min_pcr_replicates = args.min_pcr_replicates

    if args.gdb_dir and args.no_gdb:
        logger.error("Give --gdb-dir or --no-gdb, not both.")
        return None
    if args.gdb_dir:
        config.gdb_dir = Path(args.gdb_dir)
        config.has_gdb_data = True
    elif args.no_gdb:
        config.gdb_dir = None
        config.has_gdb_data = False

    # Validate configuration
    errors = config.validate()
    if errors:
        for err in errors:
            logger.error(err)
        return None

    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Log detected data sources
    logger.info("\nDetected data sources:")
    logger.info(f"  2x150 data: {config.has_2x150_data}")
    logger.info(f"  2x300 data: {config.has_2x300_data}")
    logger.info(f"  gdb data: {config.has_gdb_data}")

    if not any([config.has_2x150_data, config.has_2x300_data]):
        logger.error("No data sources found. Please check project directory structure.")
        return None

    # Step 1: Load all available data sources
    logger.info("\n" + "=" * 60)
    logger.info("STEP 1: Loading data sources")
    logger.info("=" * 60)

    data_sources = load_all_data_sources(
        config,
        debug_mode=args.debug,
        debug_n_samples=args.debug_n_samples
    )
    # G3 (27_): a gdb directory that yields ZERO tables used to be a warning and an
    # empty frame -- the run then "succeeded" with no guide disambiguation at all.
    if config.has_gdb_data and 'gdb' not in data_sources:
        logger.error(f"gdb directory {config.gdb_dir} yielded no guide-donor-bc0 tables. Pass "
                     f"--gdb-dir <dir with the tables> or --no-gdb to run without them explicitly.")
        return None

    if not data_sources:
        logger.error("No data loaded from any source")
        return None

    for name, (df, info) in data_sources.items():
        logger.info(f"\n{name}: {info.num_rows} rows, {info.num_unique_bc1s} unique bc1s")

    # Step 2: Merge data sources
    logger.info("\n" + "=" * 60)
    logger.info("STEP 2: Merging data sources")
    logger.info("=" * 60)

    df_2x150 = data_sources.get('2x150', (None, None))[0]
    df_2x300 = data_sources.get('2x300', (None, None))[0]
    df_gdb = data_sources.get('gdb', (None, None))[0]

    merged_df = merge_data_sources(df_2x150, df_2x300, df_gdb)

    if merged_df.empty:
        logger.error("Merge resulted in empty DataFrame")
        return None

    logger.info(f"\nMerged dataset: {len(merged_df)} rows, {merged_df['bc1'].nunique()} unique bc1s")

    # Step 3: Aggregate and filter
    logger.info("\n" + "=" * 60)
    logger.info("STEP 3: Aggregating and filtering")
    logger.info("=" * 60)

    filtered_df = aggregate_and_filter(
        merged_df,
        min_reads=config.min_reads,
        min_pcr_replicates=config.min_pcr_replicates,
        collapse_sequences=not args.debug,  # Skip collapsing in debug mode for speed
        show_progress=True
    )

    logger.info(f"\nFiltered dataset: {len(filtered_df)} rows, {filtered_df['bc1'].nunique()} unique bc1s")

    # Step 4: Oligo matching
    logger.info("\n" + "=" * 60)
    logger.info("STEP 4: Matching to designed oligos")
    logger.info("=" * 60)

    # Build lookup tables from the harmonized oligo table (covers all four nuclease
    # classes -- SpG / SpCas9 / LbCas12a / impLbCas12a -- plus the cAMP/Ras/PKA
    # saturation oligos (V574-V597) added 2026-05-30). The previous two-file path
    # hard-coded `20210422_Twist_200mer_SpCas9_oligos_with_controls.tsv`, which has
    # never existed on /oak, so this step was a silent no-op (logged "Oligo design
    # files not found, skipping matching"). The harmonized table lives at the
    # canonical `common/annotation_files/` path resolved inside
    # `magestic.utils.oligo_annotations.HARMONIZED_OLIGO_TABLE`.
    logger.info("Building oligo lookup tables from harmonized table...")
    sublibrary_tables = load_and_build_lookup(use_harmonized_table=True)

    # Determine which columns to use for matching
    donor_col = 'donor' if 'donor' in filtered_df.columns else 'bdb_donor'
    guide_col = 'guide' if 'guide' in filtered_df.columns else 'gdb_guide'
    # Check for yeast_sublibrary column (or legacy "sublibrary" for backward compatibility)
    yeast_sublibrary_col = None
    if 'yeast_sublibrary' in filtered_df.columns:
        yeast_sublibrary_col = 'yeast_sublibrary'
    elif 'sublibrary' in filtered_df.columns:
        yeast_sublibrary_col = 'sublibrary'

    logger.info("Matching sequences to oligos...")
    matched_df = match_dataframe(
        filtered_df,
        sublibrary_tables,
        donor_col=donor_col,
        guide_col=guide_col if guide_col in filtered_df.columns else 'guide',
        prefix=''
    )

    # Step 5: Output results
    logger.info("\n" + "=" * 60)
    logger.info("STEP 5: Saving results")
    logger.info("=" * 60)

    timestamp = datetime.now().strftime("%Y%m%d")
    output_file = config.output_dir / f"{timestamp}_bc1_reference_table_min_{config.min_reads}_reads.tsv"

    matched_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved reference table to: {output_file}")

    # Summary statistics
    logger.info("\n" + "=" * 60)
    logger.info("PIPELINE SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total unique bc1s: {matched_df['bc1'].nunique()}")
    logger.info(f"Total rows: {len(matched_df)}")

    if 'oligo_name' in matched_df.columns:
        matched = matched_df['oligo_name'].notna().sum()
        logger.info(f"Matched to oligos: {matched} ({100*matched/len(matched_df):.1f}%)")

    if 'data_source' in matched_df.columns:
        logger.info("\nData source breakdown:")
        for source in matched_df['data_source'].unique():
            count = len(matched_df[matched_df['data_source'].str.contains(source, na=False)])
            logger.info(f"  {source}: {count}")

    # Step 6: Annotate unmatched donors (optional, flag-gated)
    if getattr(args, "annotate_unmatched", False):
        logger.info("\n" + "=" * 60)
        logger.info("STEP 6: Annotating unmatched donors")
        logger.info("=" * 60)
        _run_unmatched_annotation(args, config, output_file)

    return matched_df


def main():
    """Entry point."""
    args = parse_args()
    result = run_pipeline(args)

    if result is not None:
        logger.info("\nPipeline completed successfully!")
        return 0
    else:
        logger.error("\nPipeline failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
