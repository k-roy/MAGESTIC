#!/usr/bin/env python3
"""
Batch parser for bc1-donor-bc0 sequences.

Processes multiple samples in parallel using multiprocessing.
Designed to be called from SLURM array jobs where each task handles a batch of samples.

Usage:
    # With explicit sample list
    python 01_parse_batch.py --sample-list samples.txt --start 1 --end 32 --bc1-dir /path/to/fastqs --output-dir parsed/ --n-workers 8

    # With auto-discovery (finds all bc1_donor_bc0 FASTQs in project)
    python 01_parse_batch.py --project-dir /path/to/QTL/project --n-workers 8

    # Discovery-only mode (list samples without processing)
    python 01_parse_batch.py --discover --project-dir /path/to/QTL/project

Author: Kevin R. Roy
Date: 2026-03-22
"""

import argparse
import sys
import os
import subprocess
from pathlib import Path
from multiprocessing import Pool
from typing import List, Tuple, Optional
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Path to the parser script
SCRIPT_DIR = Path(__file__).parent
PARSER_SCRIPT = SCRIPT_DIR / "01_parse_bc1_donor_bc0.py"
PYTHON = sys.executable

# Base directories for auto-discovery
BASE_DIR = Path("/path/to")
QTL_PROJECTS_DIR = BASE_DIR / "projects/QTL"


def discover_bc1_donor_bc0_samples(project_dir: Optional[Path] = None) -> List[Tuple[Path, Path, Path]]:
    """
    Auto-discover bc1_donor_bc0 R1/R2 FASTQ pairs across QTL projects.

    Args:
        project_dir: Specific project directory to search, or None to search all QTL projects

    Returns:
        List of (R1_path, R2_path, output_path) tuples
    """
    samples = []

    if project_dir:
        search_dirs = [project_dir / "raw_data"]
    else:
        search_dirs = list(QTL_PROJECTS_DIR.glob("*/raw_data"))

    for raw_data_dir in search_dirs:
        if not raw_data_dir.exists():
            continue

        # Find bc1_donor_bc0 directories
        bc1_dirs = list(raw_data_dir.glob("*bc1_donor_bc0*"))
        for bc1_dir in bc1_dirs:
            # Find R1 files
            r1_files = list(bc1_dir.glob("**/*_R1.fastq.gz"))
            for r1_path in r1_files:
                r2_path = Path(str(r1_path).replace("_R1.fastq.gz", "_R2.fastq.gz"))
                if r2_path.exists():
                    # Determine output path
                    project_name = raw_data_dir.parent.name
                    sample_name = r1_path.stem.replace("_R1.fastq", "")
                    output_dir = raw_data_dir.parent / "processed_data" / "bc1_donor_bc0" / "parsed"
                    output_path = output_dir / f"{sample_name}_bc1_donor_bc0.tsv"
                    samples.append((r1_path, r2_path, output_path))

    return samples


def process_sample(args):
    """Process a single sample. Returns (sample_name, success, message)."""
    sample_name, r1_path, r2_path, output_path = args

    try:
        # Ensure output directory exists
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

        # Skip if output exists
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            return (sample_name, True, "skipped (exists)")

        # Call the parser script via subprocess
        cmd = [PYTHON, str(PARSER_SCRIPT), r1_path, r2_path, output_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            return (sample_name, False, result.stderr.strip() or "Unknown error")

        return (sample_name, True, "completed")
    except subprocess.TimeoutExpired:
        return (sample_name, False, "timeout (>10 min)")
    except Exception as e:
        return (sample_name, False, str(e))


def find_fastq(bc1_dir: Path, sample_name: str) -> tuple:
    """Find R1/R2 FASTQ files for a sample."""
    # Search in bc1_dir and subdirectories
    for r1_path in bc1_dir.rglob(f"{sample_name}_R1.fastq.gz"):
        r2_path = r1_path.parent / f"{sample_name}_R2.fastq.gz"
        if r2_path.exists():
            return str(r1_path), str(r2_path)
    return None, None


def main():
    parser = argparse.ArgumentParser(
        description="Batch process bc1-donor-bc0 samples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # With explicit sample list (for SLURM array jobs)
    python 01_parse_batch.py --sample-list samples.txt --start 1 --end 32 \\
        --bc1-dir /path/to/fastqs --output-dir parsed/ --n-workers 8

    # With auto-discovery (finds all bc1_donor_bc0 FASTQs in project)
    python 01_parse_batch.py --project-dir /path/to/QTL/project --n-workers 8

    # Discovery-only mode (list samples without processing)
    python 01_parse_batch.py --discover --project-dir /path/to/QTL/project
        """
    )

    # Mode selection
    mode_group = parser.add_argument_group("Mode selection (choose one)")
    mode_group.add_argument("--sample-list", help="File with sample names (one per line)")
    mode_group.add_argument("--project-dir", type=Path, help="Project directory for auto-discovery")
    mode_group.add_argument("--discover", action="store_true", help="Discovery-only mode (list samples, don't process)")

    # Sample list mode options
    list_group = parser.add_argument_group("Sample list options")
    list_group.add_argument("--start", type=int, help="Start index (1-based)")
    list_group.add_argument("--end", type=int, help="End index (1-based, inclusive)")
    list_group.add_argument("--bc1-dir", help="Directory containing bc1_donor_bc0 FASTQs")
    list_group.add_argument("--output-dir", help="Output directory for parsed TSVs")

    # Common options
    common_group = parser.add_argument_group("Common options")
    common_group.add_argument("--n-workers", type=int, default=8, help="Number of parallel workers")
    common_group.add_argument("--limit", type=int, help="Limit number of samples to process (for testing)")
    common_group.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")

    args = parser.parse_args()

    # Determine mode and validate arguments
    if args.discover or args.project_dir:
        # Auto-discovery mode
        logger.info("Using auto-discovery mode")
        samples = discover_bc1_donor_bc0_samples(args.project_dir)
        logger.info(f"Discovered {len(samples)} sample pairs")

        if args.discover:
            # Discovery-only: print summary by project and exit
            projects = {}
            for r1, r2, out in samples:
                # Extract project name from path
                parts = r1.parts
                try:
                    qtl_idx = parts.index("QTL")
                    project = parts[qtl_idx + 1] if qtl_idx + 1 < len(parts) else "unknown"
                except ValueError:
                    project = "unknown"

                if project not in projects:
                    projects[project] = []
                projects[project].append((r1, r2, out))

            print("\nDiscovered samples by project:")
            for project, project_samples in sorted(projects.items()):
                print(f"  {project}: {len(project_samples)} samples")
                # Show first few samples
                for r1, r2, out in project_samples[:3]:
                    print(f"    - {r1.name}")
                if len(project_samples) > 3:
                    print(f"    ... and {len(project_samples) - 3} more")
            return

        # Prepare tasks from discovered samples
        if args.limit:
            samples = samples[:args.limit]
            logger.info(f"Limited to {len(samples)} samples")

        tasks = []
        for r1_path, r2_path, output_path in samples:
            sample_name = output_path.stem.replace("_bc1_donor_bc0", "")
            if not args.overwrite and output_path.exists() and output_path.stat().st_size > 0:
                logger.info(f"{sample_name}: skipped (exists)")
                continue
            tasks.append((sample_name, str(r1_path), str(r2_path), str(output_path)))

    elif args.sample_list:
        # Sample list mode - validate required arguments
        if not all([args.start, args.end, args.bc1_dir, args.output_dir]):
            parser.error("--sample-list mode requires --start, --end, --bc1-dir, and --output-dir")

        # Read sample list
        with open(args.sample_list) as f:
            all_samples = [line.strip() for line in f if line.strip()]

        # Get samples for this batch (1-based indexing)
        start_idx = args.start - 1
        end_idx = min(args.end, len(all_samples))
        samples = all_samples[start_idx:end_idx]

        logger.info(f"Processing samples {args.start}-{end_idx} ({len(samples)} samples)")

        # Prepare tasks
        bc1_dir = Path(args.bc1_dir)
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        tasks = []
        for sample in samples:
            r1_path, r2_path = find_fastq(bc1_dir, sample)
            if r1_path is None:
                logger.warning(f"FASTQ not found for {sample}")
                continue
            output_path = str(output_dir / f"{sample}_bc1_donor_bc0.tsv")
            tasks.append((sample, r1_path, r2_path, output_path))
    else:
        parser.error("Must specify either --sample-list or --project-dir (or --discover)")

    logger.info(f"Found {len(tasks)} samples to process")
    logger.info(f"Using {args.n_workers} workers")

    if not tasks:
        logger.info("No samples to process")
        return

    # Process in parallel
    completed = 0
    failed = 0
    skipped = 0

    with Pool(processes=args.n_workers) as pool:
        for sample_name, success, message in pool.imap_unordered(process_sample, tasks):
            if success:
                if "skipped" in message:
                    skipped += 1
                else:
                    completed += 1
                logger.info(f"{sample_name}: {message}")
            else:
                failed += 1
                logger.error(f"{sample_name}: FAILED - {message}")

    logger.info(f"Batch complete: {completed} completed, {skipped} skipped, {failed} failed")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
