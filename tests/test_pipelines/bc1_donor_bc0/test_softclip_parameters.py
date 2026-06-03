#!/usr/bin/env python3
"""
Test soft-clipping parameters for the triple-aligner approach in step_g.

This script tests various alignment parameters to minimize soft-clipping
while maintaining alignment quality for donor sequence alignment.

The goal is to find optimal parameters for:
1. BWA-MEM: -L penalty for clipping
2. BBMap: local=f (global alignment mode)
3. Minimap2: --end-bonus parameter

Metrics evaluated:
- Soft-clip rate (lower is better)
- Perfect match rate (higher is better)
- Mapping rate (higher is better)
- Alignment score distribution

Usage:
    python test_softclip_parameters.py \\
        --query-fasta <donors.fa> \\
        --reference <genome.fa> \\
        --output-dir <output_dir> \\
        [--max-queries 10000]

Author: Kevin R. Roy
"""

import argparse
import subprocess
import tempfile
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
from collections import defaultdict

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False
    print("Error: pysam required for alignment parsing")
    exit(1)


@dataclass
class AlignmentStats:
    """Statistics for an alignment method."""
    method: str
    total_queries: int = 0
    mapped: int = 0
    unmapped: int = 0
    has_softclip: int = 0
    has_5prime_clip: int = 0
    has_3prime_clip: int = 0
    total_softclip_bp: int = 0
    max_softclip_bp: int = 0
    perfect_match: int = 0
    scores: List[int] = None

    def __post_init__(self):
        if self.scores is None:
            self.scores = []

    @property
    def map_rate(self) -> float:
        return 100 * self.mapped / self.total_queries if self.total_queries > 0 else 0

    @property
    def softclip_rate(self) -> float:
        return 100 * self.has_softclip / self.mapped if self.mapped > 0 else 0

    @property
    def perfect_rate(self) -> float:
        return 100 * self.perfect_match / self.mapped if self.mapped > 0 else 0

    @property
    def mean_softclip_bp(self) -> float:
        return self.total_softclip_bp / self.mapped if self.mapped > 0 else 0

    def to_dict(self) -> dict:
        return {
            "method": self.method,
            "total_queries": self.total_queries,
            "mapped": self.mapped,
            "unmapped": self.unmapped,
            "map_rate": f"{self.map_rate:.1f}%",
            "softclip_count": self.has_softclip,
            "softclip_rate": f"{self.softclip_rate:.1f}%",
            "has_5prime_clip": self.has_5prime_clip,
            "has_3prime_clip": self.has_3prime_clip,
            "mean_softclip_bp": f"{self.mean_softclip_bp:.1f}",
            "max_softclip_bp": self.max_softclip_bp,
            "perfect_match": self.perfect_match,
            "perfect_rate": f"{self.perfect_rate:.1f}%",
            "mean_score": np.mean(self.scores) if self.scores else 0,
            "median_score": np.median(self.scores) if self.scores else 0,
        }


def check_aligner(name: str) -> bool:
    """Check if aligner is available."""
    try:
        if name == "bwa":
            subprocess.run(["bwa"], capture_output=True)
            return True
        elif name == "bbmap":
            subprocess.run(["bbmap.sh", "--version"], capture_output=True)
            return True
        elif name == "minimap2":
            result = subprocess.run(["minimap2", "--version"], capture_output=True)
            return result.returncode == 0
    except FileNotFoundError:
        return False
    return False


def run_alignment(
    query_fasta: Path,
    reference: Path,
    output_sam: Path,
    method: str,
    threads: int = 4,
) -> bool:
    """Run alignment with specified method and parameters."""
    try:
        if method == "bwa_default":
            # Default BWA-MEM
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["bwa", "mem", "-t", str(threads), str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        elif method == "bwa_L50":
            # BWA-MEM with moderate clipping penalty
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["bwa", "mem", "-t", str(threads), "-L", "50,50", str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        elif method == "bwa_L100":
            # BWA-MEM with high clipping penalty
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["bwa", "mem", "-t", str(threads), "-L", "100,100", str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        elif method == "bwa_L200":
            # BWA-MEM with very high clipping penalty
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["bwa", "mem", "-t", str(threads), "-L", "200,200", str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        elif method == "bbmap_default":
            # BBMap default (local=t)
            subprocess.run([
                "bbmap.sh",
                f"in={query_fasta}",
                f"out={output_sam}",
                f"ref={reference}",
                "nmtag=t", "mdtag=t",
                f"threads={threads}",
            ], check=True, capture_output=True)

        elif method == "bbmap_global":
            # BBMap global alignment (local=f)
            subprocess.run([
                "bbmap.sh",
                f"in={query_fasta}",
                f"out={output_sam}",
                f"ref={reference}",
                "local=f",
                "nmtag=t", "mdtag=t",
                f"threads={threads}",
            ], check=True, capture_output=True)

        elif method == "bbmap_semiglobal":
            # BBMap semiglobal (local=f plus additional parameters)
            subprocess.run([
                "bbmap.sh",
                f"in={query_fasta}",
                f"out={output_sam}",
                f"ref={reference}",
                "local=f",
                "perfectmode=f",  # Allow mismatches
                "semiperfectmode=f",
                "nmtag=t", "mdtag=t",
                f"threads={threads}",
            ], check=True, capture_output=True)

        elif method == "minimap2_default":
            # Minimap2 default for short reads
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["minimap2", "-a", "-x", "sr", "-t", str(threads), str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        elif method == "minimap2_bonus25":
            # Minimap2 with end bonus 25
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["minimap2", "-a", "-x", "sr", "--end-bonus", "25", "-t", str(threads), str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        elif method == "minimap2_bonus50":
            # Minimap2 with end bonus 50
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["minimap2", "-a", "-x", "sr", "--end-bonus", "50", "-t", str(threads), str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        elif method == "minimap2_bonus100":
            # Minimap2 with end bonus 100
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["minimap2", "-a", "-x", "sr", "--end-bonus", "100", "-t", str(threads), str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        elif method == "minimap2_asm5":
            # Minimap2 with asm5 preset (for high-identity sequences)
            with open(output_sam, "w") as out:
                subprocess.run(
                    ["minimap2", "-a", "-x", "asm5", "--end-bonus", "50", "-t", str(threads), str(reference), str(query_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True
                )

        else:
            print(f"Unknown method: {method}")
            return False

        return True

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  Error running {method}: {e}")
        return False


def analyze_sam(sam_path: Path, method: str) -> AlignmentStats:
    """Analyze SAM file for alignment statistics."""
    stats = AlignmentStats(method=method)

    try:
        samfile = pysam.AlignmentFile(str(sam_path), "r")
    except Exception as e:
        print(f"  Error opening {sam_path}: {e}")
        return stats

    for read in samfile:
        stats.total_queries += 1

        if read.is_unmapped:
            stats.unmapped += 1
            continue

        stats.mapped += 1

        # Check for soft-clipping
        cigar = read.cigartuples or []
        has_5prime = False
        has_3prime = False
        softclip_bp = 0

        for i, (op, length) in enumerate(cigar):
            if op == 4:  # Soft-clip
                softclip_bp += length
                if i == 0:
                    has_5prime = True
                if i == len(cigar) - 1:
                    has_3prime = True

        if softclip_bp > 0:
            stats.has_softclip += 1
            stats.total_softclip_bp += softclip_bp
            stats.max_softclip_bp = max(stats.max_softclip_bp, softclip_bp)

        if has_5prime:
            stats.has_5prime_clip += 1
        if has_3prime:
            stats.has_3prime_clip += 1

        # Check for perfect match (no mismatches, no indels, no soft-clips)
        nm = read.get_tag('NM') if read.has_tag('NM') else -1
        if nm == 0 and softclip_bp == 0:
            stats.perfect_match += 1

        # Get alignment score
        if read.has_tag('AS'):
            stats.scores.append(read.get_tag('AS'))

    samfile.close()
    return stats


def main():
    parser = argparse.ArgumentParser(description="Test soft-clipping parameters for aligners")
    parser.add_argument("--query-fasta", "-q", required=True, type=Path,
                        help="Query FASTA file (donor sequences)")
    parser.add_argument("--reference", "-r", required=True, type=Path,
                        help="Reference genome FASTA")
    parser.add_argument("--output-dir", "-o", required=True, type=Path,
                        help="Output directory for results")
    parser.add_argument("--max-queries", type=int, default=10000,
                        help="Maximum queries to test (default: 10000)")
    parser.add_argument("--threads", "-t", type=int, default=4,
                        help="Threads per aligner (default: 4)")

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Check aligners
    print("Checking aligners...")
    aligners_available = {
        "bwa": check_aligner("bwa"),
        "bbmap": check_aligner("bbmap"),
        "minimap2": check_aligner("minimap2"),
    }
    for name, avail in aligners_available.items():
        print(f"  {name}: {'available' if avail else 'NOT FOUND'}")

    # Subset query FASTA if needed
    query_fasta = args.query_fasta
    if args.max_queries:
        print(f"\nSubsetting to {args.max_queries} queries...")
        subset_fasta = args.output_dir / "queries_subset.fa"
        count = 0
        with open(query_fasta) as inf, open(subset_fasta, "w") as outf:
            for line in inf:
                if line.startswith(">"):
                    count += 1
                    if count > args.max_queries:
                        break
                outf.write(line)
        query_fasta = subset_fasta
        print(f"  Wrote {count - 1} sequences")

    # Define methods to test
    methods = []

    if aligners_available["bwa"]:
        methods.extend([
            "bwa_default",
            "bwa_L50",
            "bwa_L100",
            "bwa_L200",
        ])

    if aligners_available["bbmap"]:
        methods.extend([
            "bbmap_default",
            "bbmap_global",
            "bbmap_semiglobal",
        ])

    if aligners_available["minimap2"]:
        methods.extend([
            "minimap2_default",
            "minimap2_bonus25",
            "minimap2_bonus50",
            "minimap2_bonus100",
            "minimap2_asm5",
        ])

    # Check BWA index
    if aligners_available["bwa"]:
        index_file = Path(str(args.reference) + ".bwt")
        if not index_file.exists():
            print("\nCreating BWA index...")
            subprocess.run(["bwa", "index", str(args.reference)], check=True, capture_output=True)

    # Run alignments
    print(f"\nRunning {len(methods)} alignment methods...")
    results = []

    for method in methods:
        print(f"\n  {method}...")
        sam_path = args.output_dir / f"{method}.sam"

        if run_alignment(query_fasta, args.reference, sam_path, method, args.threads):
            stats = analyze_sam(sam_path, method)
            results.append(stats)
            print(f"    Mapped: {stats.map_rate:.1f}%, Softclip: {stats.softclip_rate:.1f}%, Perfect: {stats.perfect_rate:.1f}%")

    # Create summary table
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    summary_data = [r.to_dict() for r in results]
    summary_df = pd.DataFrame(summary_data)

    # Sort by softclip rate (ascending), then perfect rate (descending)
    summary_df["softclip_rate_num"] = summary_df["softclip_rate"].str.rstrip('%').astype(float)
    summary_df["perfect_rate_num"] = summary_df["perfect_rate"].str.rstrip('%').astype(float)
    summary_df = summary_df.sort_values(["softclip_rate_num", "perfect_rate_num"], ascending=[True, False])
    summary_df = summary_df.drop(columns=["softclip_rate_num", "perfect_rate_num"])

    print("\nRanked by lowest soft-clip rate:")
    print(summary_df[["method", "map_rate", "softclip_rate", "perfect_rate", "mean_softclip_bp"]].to_string(index=False))

    # Save results
    output_file = args.output_dir / "softclip_parameter_comparison.tsv"
    summary_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nResults saved to: {output_file}")

    # Recommendations
    print("\n" + "=" * 80)
    print("RECOMMENDATIONS")
    print("=" * 80)

    best = summary_df.iloc[0]
    print(f"\nBest method: {best['method']}")
    print(f"  Soft-clip rate: {best['softclip_rate']}")
    print(f"  Perfect match rate: {best['perfect_rate']}")
    print(f"  Mapping rate: {best['map_rate']}")

    print("\nFor step_g triple-aligner approach:")
    print("  Consider using:")

    # Find best for each aligner
    for aligner in ["bwa", "bbmap", "minimap2"]:
        aligner_results = summary_df[summary_df["method"].str.startswith(aligner)]
        if len(aligner_results) > 0:
            best_for_aligner = aligner_results.iloc[0]
            print(f"  - {best_for_aligner['method']}: {best_for_aligner['softclip_rate']} softclip, {best_for_aligner['perfect_rate']} perfect")


if __name__ == "__main__":
    main()
