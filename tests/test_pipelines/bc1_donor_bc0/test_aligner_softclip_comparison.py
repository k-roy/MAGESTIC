#!/usr/bin/env python3
"""
Test script to compare alignment methods with different soft-clipping settings.

Compares 5 methods:
1. BBMap local=f (no soft-clipping)
2. BWA-MEM default (allows soft-clipping)
3. BWA-MEM -L 100,100 (discourages soft-clipping)
4. Minimap2 default (allows soft-clipping)
5. Minimap2 --end-bonus 50 (discourages soft-clipping)

For each method, we:
- Align observed donor sequences to designed donor sequences
- Parse CIGAR strings to identify soft-clipping
- Compare alignment rates and soft-clip frequencies

Author: Kevin R. Roy
"""

import argparse
import subprocess
import os
import pandas as pd
import pysam
from pathlib import Path
from collections import Counter
import re


def parse_cigar(cigar_string):
    """
    Parse CIGAR string and return statistics.

    Returns dict with:
    - total_len: total alignment length
    - matches: number of M/=/X operations
    - insertions: number of I operations
    - deletions: number of D operations
    - soft_clips: number of S operations (soft-clipping)
    - soft_clip_5prime: soft-clipping at 5' end
    - soft_clip_3prime: soft-clipping at 3' end
    """
    if cigar_string is None or cigar_string == "*":
        return None

    stats = {
        "total_len": 0,
        "matches": 0,
        "insertions": 0,
        "deletions": 0,
        "soft_clips": 0,
        "soft_clip_5prime": 0,
        "soft_clip_3prime": 0,
        "cigar": cigar_string
    }

    # Parse CIGAR operations
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)

    for i, (length, op) in enumerate(operations):
        length = int(length)
        stats["total_len"] += length

        if op in "M=X":
            stats["matches"] += length
        elif op == "I":
            stats["insertions"] += length
        elif op == "D":
            stats["deletions"] += length
        elif op == "S":
            stats["soft_clips"] += length
            if i == 0:  # First operation = 5' clip
                stats["soft_clip_5prime"] = length
            if i == len(operations) - 1:  # Last operation = 3' clip
                stats["soft_clip_3prime"] = length

    return stats


def run_bbmap(query_fasta, ref_fasta, output_sam, local=True):
    """Run BBMap alignment."""
    local_flag = "local=t" if local else "local=f"
    cmd = f"bbmap.sh in={query_fasta} out={output_sam} ref={ref_fasta} {local_flag} overwrite=t nmtag=t mdtag=t 2>&1"
    print(f"  Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False


def run_bwa_mem(query_fasta, ref_fasta, output_sam, clip_penalty=None):
    """Run BWA-MEM alignment."""
    # Index reference if needed
    if not os.path.exists(f"{ref_fasta}.bwt"):
        subprocess.run(f"bwa index {ref_fasta}", shell=True, capture_output=True)

    clip_flag = f"-L {clip_penalty},{clip_penalty}" if clip_penalty else ""
    cmd = f"bwa mem {clip_flag} {ref_fasta} {query_fasta} > {output_sam} 2>/dev/null"
    print(f"  Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False


def run_minimap2(query_fasta, ref_fasta, output_sam, end_bonus=None):
    """Run minimap2 alignment."""
    bonus_flag = f"--end-bonus {end_bonus}" if end_bonus else ""
    # Use -x sr for short reads (donors are ~129bp)
    cmd = f"minimap2 -ax sr {bonus_flag} {ref_fasta} {query_fasta} > {output_sam} 2>/dev/null"
    print(f"  Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False


def parse_sam_results(sam_file):
    """Parse SAM file and extract alignment statistics."""
    results = []

    try:
        sam = pysam.AlignmentFile(sam_file, "r")
        for read in sam:
            query_name = read.query_name
            ref_name = read.reference_name
            cigar = read.cigarstring
            mapq = read.mapping_quality
            is_unmapped = read.is_unmapped

            cigar_stats = parse_cigar(cigar) if cigar else None

            results.append({
                "query": query_name,
                "reference": ref_name,
                "mapq": mapq,
                "is_unmapped": is_unmapped,
                "cigar": cigar,
                **(cigar_stats if cigar_stats else {})
            })
        sam.close()
    except Exception as e:
        print(f"  Error parsing {sam_file}: {e}")

    return pd.DataFrame(results)


def summarize_results(df, method_name):
    """Generate summary statistics for alignment results."""
    total = len(df)
    mapped = (~df["is_unmapped"]).sum()
    unmapped = df["is_unmapped"].sum()

    mapped_df = df[~df["is_unmapped"]]

    summary = {
        "method": method_name,
        "total_queries": total,
        "mapped": mapped,
        "unmapped": unmapped,
        "map_rate": f"{100*mapped/total:.1f}%" if total > 0 else "N/A",
    }

    if len(mapped_df) > 0:
        # Soft-clipping stats
        has_softclip = (mapped_df["soft_clips"] > 0).sum()
        has_5prime_clip = (mapped_df["soft_clip_5prime"] > 0).sum()
        has_3prime_clip = (mapped_df["soft_clip_3prime"] > 0).sum()

        summary.update({
            "has_softclip": has_softclip,
            "softclip_rate": f"{100*has_softclip/len(mapped_df):.1f}%",
            "has_5prime_clip": has_5prime_clip,
            "has_3prime_clip": has_3prime_clip,
            "mean_softclip_bp": f"{mapped_df['soft_clips'].mean():.1f}",
            "max_softclip_bp": mapped_df["soft_clips"].max(),
            "perfect_match": (mapped_df["cigar"].str.match(r'^\d+M$')).sum() if "cigar" in mapped_df else 0,
        })

    return summary


def main():
    parser = argparse.ArgumentParser(description="Compare alignment methods for soft-clipping behavior")
    parser.add_argument("--query", required=True, help="Query FASTA (observed donors)")
    parser.add_argument("--reference", required=True, help="Reference FASTA (designed donors)")
    parser.add_argument("--output-dir", required=True, help="Output directory for SAM files and results")
    parser.add_argument("--max-queries", type=int, default=None, help="Limit number of queries for testing")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    query_fasta = args.query
    ref_fasta = args.reference

    # Optionally subset queries for faster testing
    if args.max_queries:
        subset_fasta = output_dir / "query_subset.fasta"
        print(f"Subsetting to {args.max_queries} queries...")
        with open(query_fasta) as f_in, open(subset_fasta, "w") as f_out:
            count = 0
            for line in f_in:
                f_out.write(line)
                if line.startswith(">"):
                    count += 1
                    if count > args.max_queries:
                        break
        query_fasta = str(subset_fasta)

    # Define alignment methods to test
    methods = [
        ("bbmap_no_softclip", lambda q, r, o: run_bbmap(q, r, o, local=False)),
        ("bwa_default", lambda q, r, o: run_bwa_mem(q, r, o, clip_penalty=None)),
        ("bwa_no_softclip", lambda q, r, o: run_bwa_mem(q, r, o, clip_penalty=100)),
        ("minimap2_default", lambda q, r, o: run_minimap2(q, r, o, end_bonus=None)),
        ("minimap2_no_softclip", lambda q, r, o: run_minimap2(q, r, o, end_bonus=50)),
    ]

    all_results = {}
    summaries = []

    print("\n" + "=" * 60)
    print("ALIGNER SOFT-CLIPPING COMPARISON TEST")
    print("=" * 60)
    print(f"Query: {query_fasta}")
    print(f"Reference: {ref_fasta}")
    print(f"Output: {output_dir}")
    print()

    for method_name, align_func in methods:
        print(f"\n--- {method_name} ---")
        sam_file = output_dir / f"{method_name}.sam"

        # Run alignment
        success = align_func(query_fasta, ref_fasta, str(sam_file))

        if success and sam_file.exists():
            # Parse results
            df = parse_sam_results(str(sam_file))
            all_results[method_name] = df

            # Summarize
            summary = summarize_results(df, method_name)
            summaries.append(summary)

            print(f"  Mapped: {summary['mapped']}/{summary['total_queries']} ({summary['map_rate']})")
            if "softclip_rate" in summary:
                print(f"  Soft-clipped: {summary['has_softclip']} ({summary['softclip_rate']})")
                print(f"  Perfect matches: {summary['perfect_match']}")
        else:
            print(f"  FAILED")

    # Save summary table
    summary_df = pd.DataFrame(summaries)
    summary_file = output_dir / "alignment_comparison_summary.tsv"
    summary_df.to_csv(summary_file, sep="\t", index=False)
    print(f"\n\nSummary saved to: {summary_file}")

    # Print comparison table
    print("\n" + "=" * 60)
    print("SUMMARY COMPARISON")
    print("=" * 60)
    print(summary_df.to_string(index=False))

    # Detailed soft-clip analysis
    print("\n" + "=" * 60)
    print("SOFT-CLIP DETAILS")
    print("=" * 60)

    for method_name, df in all_results.items():
        mapped_df = df[~df["is_unmapped"]]
        if len(mapped_df) > 0:
            softclip_df = mapped_df[mapped_df["soft_clips"] > 0]
            if len(softclip_df) > 0:
                print(f"\n{method_name}: {len(softclip_df)} alignments with soft-clipping")
                print("  Soft-clip length distribution:")
                print(softclip_df["soft_clips"].value_counts().head(10).to_string())

    # Save detailed results
    for method_name, df in all_results.items():
        detail_file = output_dir / f"{method_name}_details.tsv"
        df.to_csv(detail_file, sep="\t", index=False)

    print(f"\nDetailed results saved to {output_dir}/")


if __name__ == "__main__":
    main()
