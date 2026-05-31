#!/bin/bash
# ============================================================================
# Step 01: Trim, Merge, and Collapse BC1 from Screen FASTQ
# ============================================================================
#
# This script processes raw sequencing data for BC1 screen samples:
# 1. Trim reads for quality and TruSeq adapters (BBDuk)
# 2. Merge paired-end reads (BBMerge)
# 3. Collapse identical reads
#
# Can be run standalone or via Snakemake.
#
# Usage:
#   01_trim_merge_collapse.sh --r1 input_R1.fastq.gz --r2 input_R2.fastq.gz --output output_collapsed.fastq
#
#   # With scratch (recommended for SLURM jobs)
#   01_trim_merge_collapse.sh --r1 input_R1.fastq.gz --r2 input_R2.fastq.gz --output output_collapsed.fastq --use-scratch
#
# Author: Kevin R. Roy
# ============================================================================

set -e

# Defaults
USE_SCRATCH=false
THREADS=${SLURM_CPUS_PER_TASK:-8}
ADAPTER="AGATCGGAAGAGCGTCGTGTAGGGAAAGA"
MIN_COUNTS=1

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --r1)
            R1_INPUT="$2"
            shift 2
            ;;
        --r2)
            R2_INPUT="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --use-scratch)
            USE_SCRATCH=true
            shift
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --adapter)
            ADAPTER="$2"
            shift 2
            ;;
        --min-counts)
            MIN_COUNTS="$2"
            shift 2
            ;;
        --help)
            echo "Usage: 01_trim_merge_collapse.sh --r1 R1.fastq.gz --r2 R2.fastq.gz --output collapsed.fastq [--use-scratch] [--threads N]"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$R1_INPUT" ]] || [[ -z "$R2_INPUT" ]] || [[ -z "$OUTPUT" ]]; then
    echo "ERROR: --r1, --r2, and --output are required"
    echo "Usage: 01_trim_merge_collapse.sh --r1 R1.fastq.gz --r2 R2.fastq.gz --output collapsed.fastq"
    exit 1
fi

# Thread limiting
export OMP_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS

# Get sample name from R1
SAMPLE=$(basename "$R1_INPUT" | sed 's/_R1.fastq.gz//' | sed 's/_R1.fq.gz//')

echo "=============================================="
echo "Step 01: Trim, Merge, Collapse"
echo "Sample: $SAMPLE"
echo "Threads: $THREADS"
echo "Use Scratch: $USE_SCRATCH"
echo "Date: $(date)"
echo "=============================================="

# Set up working directory
if [[ "$USE_SCRATCH" == true ]] && [[ -n "$SCRATCH" ]]; then
    WORK_DIR="${SCRATCH}/step_01_$$"
    echo "Using scratch: $WORK_DIR"
else
    WORK_DIR=$(mktemp -d)
    echo "Using temp dir: $WORK_DIR"
fi

mkdir -p "$WORK_DIR"/{trimmed,merged,collapsed}

# Ensure output directory exists
mkdir -p "$(dirname "$OUTPUT")"

# Cleanup on exit
cleanup() {
    if [[ -d "$WORK_DIR" ]]; then
        rm -rf "$WORK_DIR"
    fi
}
trap cleanup EXIT

# ============================================================================
# Step 1: Trim reads for quality and adapters
# ============================================================================

echo ""
echo "Step 1: Trimming reads..."

bbduk.sh \
    in1="$R1_INPUT" \
    in2="$R2_INPUT" \
    out1="$WORK_DIR/trimmed/${SAMPLE}_trimmed_R1.fastq.gz" \
    out2="$WORK_DIR/trimmed/${SAMPLE}_trimmed_R2.fastq.gz" \
    literal="$ADAPTER" \
    overwrite=true \
    ftl=20 \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    tpe \
    tbo \
    threads=$THREADS \
    2>&1 | tail -5

# ============================================================================
# Step 2: Merge paired-end reads
# ============================================================================

echo ""
echo "Step 2: Merging reads..."

bbmerge.sh \
    in1="$WORK_DIR/trimmed/${SAMPLE}_trimmed_R1.fastq.gz" \
    in2="$WORK_DIR/trimmed/${SAMPLE}_trimmed_R2.fastq.gz" \
    out="$WORK_DIR/merged/${SAMPLE}_merged.fastq" \
    outu="$WORK_DIR/merged/${SAMPLE}_unmerged_R1.fastq" \
    outu2="$WORK_DIR/merged/${SAMPLE}_unmerged_R2.fastq" \
    overwrite=true \
    threads=$THREADS \
    2>&1 | tail -5

# ============================================================================
# Step 3: Collapse identical reads
# ============================================================================

echo ""
echo "Step 3: Collapsing identical reads..."

# Use Python to collapse reads
python3 - "$WORK_DIR/merged/${SAMPLE}_merged.fastq" "$WORK_DIR/collapsed/${SAMPLE}_collapsed.fastq" "$MIN_COUNTS" << 'PYTHON_SCRIPT'
#!/usr/bin/env python3
"""Collapse identical reads, keeping count in header."""
import sys
from collections import defaultdict

input_file = sys.argv[1]
output_file = sys.argv[2]
min_counts = int(sys.argv[3]) if len(sys.argv) > 3 else 1

# Count sequences
seq_counts = defaultdict(int)
seq_qual = {}  # Store best quality for each sequence

with open(input_file, 'r') as f:
    while True:
        header = f.readline().strip()
        if not header:
            break
        seq = f.readline().strip()
        plus = f.readline().strip()
        qual = f.readline().strip()

        seq_counts[seq] += 1
        # Keep quality string with highest average
        if seq not in seq_qual or len(qual) > len(seq_qual[seq]):
            seq_qual[seq] = qual

# Write collapsed sequences
with open(output_file, 'w') as f:
    seq_num = 0
    for seq, count in sorted(seq_counts.items(), key=lambda x: -x[1]):
        if count >= min_counts:
            seq_num += 1
            f.write(f"@seq_{seq_num}_count={count}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{seq_qual[seq]}\n")

print(f"  Collapsed {len(seq_counts)} unique sequences")
print(f"  Output {seq_num} sequences with count >= {min_counts}")
PYTHON_SCRIPT

# ============================================================================
# Copy output
# ============================================================================

echo ""
echo "Copying output..."

# Gzip if output expects gzipped file
if [[ "$OUTPUT" == *.gz ]]; then
    gzip -c "$WORK_DIR/collapsed/${SAMPLE}_collapsed.fastq" > "$OUTPUT"
else
    cp "$WORK_DIR/collapsed/${SAMPLE}_collapsed.fastq" "$OUTPUT"
fi

echo ""
echo "=============================================="
echo "Step A Complete"
echo "Output: $OUTPUT"
echo "Date: $(date)"
echo "=============================================="
