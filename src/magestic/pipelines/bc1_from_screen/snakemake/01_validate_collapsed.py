#!/usr/bin/env python3
"""
Step 01 Validation: Check if existing collapsed FASTQs are valid for reuse.

This script validates that existing collapsed FASTQ files were generated with
the same trimming parameters as the current configuration. If valid, step 01
can be safely skipped to save compute time.

The validation checks:
1. Read length distribution matches expected trimmed length
2. Parameter hash matches stored hash (if available)
3. File is not corrupted (gzip integrity)

Usage:
    # Validate a single sample
    python 01_validate_collapsed.py --sample sample_name --collapsed-dir /path/to/collapsed \
        --expected-length 125 --output validation_result.json

    # Validate with parameter hash check
    python 01_validate_collapsed.py --sample sample_name --collapsed-dir /path/to/collapsed \
        --expected-length 125 --param-hash abc123 --output validation_result.json

Author: Kevin R. Roy
"""

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple
import random


def compute_param_hash(
    expected_length: int,
    adapter: str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGA",
    min_counts: int = 1
) -> str:
    """
    Compute a hash of trimming parameters for validation.

    Args:
        expected_length: Expected read length after trimming
        adapter: Adapter sequence used for trimming
        min_counts: Minimum counts threshold

    Returns:
        MD5 hash of parameters
    """
    param_string = f"length={expected_length}|adapter={adapter}|min_counts={min_counts}"
    return hashlib.md5(param_string.encode()).hexdigest()[:12]


def sample_read_lengths(
    fastq_path: Path,
    n_samples: int = 1000,
    seed: int = 42
) -> Tuple[list, int]:
    """
    Sample read lengths from a gzipped FASTQ file.

    Args:
        fastq_path: Path to gzipped FASTQ file
        n_samples: Number of reads to sample
        seed: Random seed for reproducibility

    Returns:
        Tuple of (list of sampled lengths, total reads counted)
    """
    random.seed(seed)
    lengths = []
    total_reads = 0

    try:
        with gzip.open(fastq_path, 'rt') as f:
            # First pass: count total reads (sample every 4th line = sequence line)
            lines = []
            for i, line in enumerate(f):
                if i % 4 == 1:  # Sequence line
                    lines.append(len(line.strip()))
                    total_reads += 1
                if total_reads >= 100000:  # Cap at 100k for speed
                    break

            # Sample from collected lengths
            if total_reads <= n_samples:
                lengths = lines
            else:
                lengths = random.sample(lines, n_samples)

    except (gzip.BadGzipFile, EOFError) as e:
        raise ValueError(f"Corrupted gzip file: {fastq_path}") from e

    return lengths, total_reads


def validate_collapsed_fastq(
    collapsed_path: Path,
    expected_length: int,
    tolerance: float = 0.05,
    min_reads: int = 100
) -> Dict:
    """
    Validate a collapsed FASTQ file against expected parameters.

    Args:
        collapsed_path: Path to collapsed FASTQ file
        expected_length: Expected read length after trimming
        tolerance: Fraction of reads allowed to deviate from expected length
        min_reads: Minimum reads required for validation

    Returns:
        Dictionary with validation results
    """
    result = {
        "file": str(collapsed_path),
        "exists": collapsed_path.exists(),
        "valid": False,
        "reason": None,
        "stats": {}
    }

    if not collapsed_path.exists():
        result["reason"] = "File does not exist"
        return result

    if collapsed_path.stat().st_size == 0:
        result["reason"] = "File is empty"
        return result

    try:
        lengths, total_reads = sample_read_lengths(collapsed_path)
    except ValueError as e:
        result["reason"] = str(e)
        return result

    if total_reads < min_reads:
        result["reason"] = f"Too few reads: {total_reads} < {min_reads}"
        return result

    # Calculate length statistics
    if not lengths:
        result["reason"] = "No reads sampled"
        return result

    mean_length = sum(lengths) / len(lengths)
    matching_length = sum(1 for l in lengths if l == expected_length)
    match_fraction = matching_length / len(lengths)

    result["stats"] = {
        "total_reads_sampled": len(lengths),
        "total_reads_estimated": total_reads,
        "mean_length": round(mean_length, 2),
        "expected_length": expected_length,
        "matching_fraction": round(match_fraction, 4),
        "min_length": min(lengths),
        "max_length": max(lengths)
    }

    # Validation criteria:
    # 1. At least (1 - tolerance) fraction of reads should match expected length
    # 2. Mean length should be within 2bp of expected
    if match_fraction >= (1 - tolerance):
        result["valid"] = True
        result["reason"] = "Passed: read lengths match expected trimming"
    elif abs(mean_length - expected_length) <= 2:
        result["valid"] = True
        result["reason"] = "Passed: mean length within tolerance"
    else:
        result["valid"] = False
        result["reason"] = f"Failed: {match_fraction:.1%} reads match expected length {expected_length}, mean={mean_length:.1f}"

    return result


def check_param_hash_file(
    collapsed_dir: Path,
    current_hash: str
) -> Optional[bool]:
    """
    Check if a parameter hash file exists and matches current parameters.

    Args:
        collapsed_dir: Directory containing collapsed FASTQs
        current_hash: Hash of current trimming parameters

    Returns:
        True if hash matches, False if different, None if no hash file
    """
    hash_file = collapsed_dir / ".trimming_params_hash"

    if not hash_file.exists():
        return None

    stored_hash = hash_file.read_text().strip()
    return stored_hash == current_hash


def write_param_hash_file(
    collapsed_dir: Path,
    param_hash: str
) -> None:
    """Write parameter hash to file for future validation."""
    hash_file = collapsed_dir / ".trimming_params_hash"
    hash_file.write_text(param_hash)


def main():
    parser = argparse.ArgumentParser(
        description="Validate collapsed FASTQ files for step A skip"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="Sample name to validate"
    )
    parser.add_argument(
        "--collapsed-dir",
        type=Path,
        required=True,
        help="Directory containing collapsed FASTQ files"
    )
    parser.add_argument(
        "--expected-length",
        type=int,
        default=125,
        help="Expected read length after trimming (default: 125)"
    )
    parser.add_argument(
        "--adapter",
        type=str,
        default="AGATCGGAAGAGCGTCGTGTAGGGAAAGA",
        help="Adapter sequence used for trimming"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.05,
        help="Fraction of reads allowed to deviate from expected (default: 0.05)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file for validation results"
    )
    parser.add_argument(
        "--update-hash",
        action="store_true",
        help="Update the parameter hash file after validation"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only output pass/fail status"
    )

    args = parser.parse_args()

    # Compute current parameter hash
    current_hash = compute_param_hash(args.expected_length, args.adapter)

    # Check for existing parameter hash
    hash_check = check_param_hash_file(args.collapsed_dir, current_hash)

    # Validate the collapsed FASTQ
    collapsed_path = args.collapsed_dir / f"{args.sample}_collapsed.fastq.gz"
    result = validate_collapsed_fastq(
        collapsed_path,
        args.expected_length,
        args.tolerance
    )

    # Add hash check results
    result["param_hash"] = {
        "current": current_hash,
        "hash_file_exists": hash_check is not None,
        "hash_matches": hash_check if hash_check is not None else "N/A"
    }

    # If hash doesn't match, mark as invalid even if lengths look OK
    if hash_check is False:
        result["valid"] = False
        result["reason"] = "Parameter hash mismatch - trimming parameters changed"

    # Determine overall skip recommendation
    can_skip = result["valid"] and (hash_check is not False)
    result["can_skip_step_01"] = can_skip

    # Output results
    if args.quiet:
        print("PASS" if can_skip else "FAIL")
    else:
        print(json.dumps(result, indent=2))

    # Write output file
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)

    # Update hash file if requested and validation passed
    if args.update_hash and result["valid"]:
        write_param_hash_file(args.collapsed_dir, current_hash)
        if not args.quiet:
            print(f"\nUpdated parameter hash file: {current_hash}")

    # Exit with appropriate code
    sys.exit(0 if can_skip else 1)


if __name__ == "__main__":
    main()
