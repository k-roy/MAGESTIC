#!/bin/bash
#SBATCH --job-name=01_trim_merge_collapse
#SBATCH --partition=mylab
#SBATCH --account=mylab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/01_trim_merge_collapse_%j.out
#SBATCH --error=logs/01_trim_merge_collapse_%j.err

# ============================================================================
# Step 01: Trim, Merge, and Collapse Guide-Donor-BC0 FASTQ
# ============================================================================
#
# This script processes raw sequencing data:
# 1. Stage input FASTQs from Oak to Scratch (high-bandwidth I/O)
# 2. Trim reads for quality and TruSeq adapters (BBDuk)
# 3. Merge paired-end reads (BBMerge)
# 4. Collapse identical reads (custom Python script)
# 5. Copy final outputs back to Oak
#
# Input: Raw paired-end FASTQ files on Oak (*_R1.fastq.gz, *_R2.fastq.gz)
# Output: Collapsed FASTQ files on Oak with read counts in header
#
# Storage pattern per CLAUDE.md:
#   Oak (input) -> Scratch (processing) -> Oak (output)
#
# Usage:
#   sbatch 01_trim_merge_collapse.sh
#
# ============================================================================

set -e

# Thread limiting: Match to SLURM allocation
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
# CRITICAL: loky (joblib's backend) uses LOKY_MAX_CPU_COUNT, not JOBLIB_WORKERS
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK

echo "=============================================="
echo "Step 01: Trim, Merge, Collapse"
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}"
echo "Date: $(date)"
echo "=============================================="

# Initialize conda
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate base

# ============================================================================
# Configuration - BASE_DIR pattern per CLAUDE.md
# ============================================================================

BASE_DIR="/path/to"
PROJECT="QTL"
SCREEN_NAME="20250912_SpG_QTL_screen"
DATATYPE="guide_donor_bc0"

# Oak paths (long-term storage)
PROJECT_DIR="${BASE_DIR}/projects/${PROJECT}/${SCREEN_NAME}"
RAW_DATA_DIR="${PROJECT_DIR}/raw_data/${DATATYPE}"
PROCESSED_DATA_DIR="${PROJECT_DIR}/processed_data/${DATATYPE}"

# Scratch paths (high-bandwidth job I/O)
# Use SCRATCH for all intermediate processing per CLAUDE.md
SCRATCH_WORK_DIR="${SCRATCH}/${PROJECT}/${SCREEN_NAME}/${DATATYPE}/step_01_${SLURM_JOB_ID:-$$}"

# Common scripts
COMMON_SCRIPTS="${BASE_DIR}/common/scripts"
COLLAPSE_SCRIPT="${COMMON_SCRIPTS}/utils/collapse_identical_merged_reads_min_counts.py"

# Fallback to general_scripts if common script doesn't exist
if [[ ! -f "${COLLAPSE_SCRIPT}" ]]; then
    COLLAPSE_SCRIPT="${BASE_DIR}/scripts_and_keyfiles/general_scripts/collapse_identical_merged_reads_min_counts.py"
fi

# ============================================================================
# Setup scratch working directory
# ============================================================================

echo ""
echo "Setting up scratch working directory..."
echo "  Scratch: ${SCRATCH_WORK_DIR}"

mkdir -p "${SCRATCH_WORK_DIR}/input"
mkdir -p "${SCRATCH_WORK_DIR}/trimmed"
mkdir -p "${SCRATCH_WORK_DIR}/merged"
mkdir -p "${SCRATCH_WORK_DIR}/unmerged"
mkdir -p "${SCRATCH_WORK_DIR}/collapsed"

# Create output directory on Oak
mkdir -p "${PROCESSED_DATA_DIR}/collapsed_fastq"
mkdir -p "${PROJECT_DIR}/scripts/${DATATYPE}/logs"

# ============================================================================
# Stage input files from Oak to Scratch
# ============================================================================

FASTQ_DIR="${RAW_DATA_DIR}/fastq"

if [[ ! -d "${FASTQ_DIR}" ]]; then
    echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}"
    exit 1
fi

echo ""
echo "Staging input FASTQs from Oak to Scratch..."
echo "  From: ${FASTQ_DIR}"
echo "  To: ${SCRATCH_WORK_DIR}/input/"

# Copy all R1 and R2 files to scratch
cp "${FASTQ_DIR}"/*_R1.fastq.gz "${SCRATCH_WORK_DIR}/input/" 2>/dev/null || true
cp "${FASTQ_DIR}"/*_R2.fastq.gz "${SCRATCH_WORK_DIR}/input/" 2>/dev/null || true

# Count samples
n_samples=$(ls -1 "${SCRATCH_WORK_DIR}/input/"*_R1.fastq.gz 2>/dev/null | wc -l)
echo "  Staged ${n_samples} sample pairs"

if [[ ${n_samples} -eq 0 ]]; then
    echo "ERROR: No FASTQ files found"
    exit 1
fi

# ============================================================================
# Process samples on Scratch (high-bandwidth I/O)
# ============================================================================

echo ""
echo "Processing on Scratch..."
cd "${SCRATCH_WORK_DIR}/input" || exit 1

sample_num=0
for sample in *_R1.fastq.gz; do
    sample_num=$((sample_num + 1))
    prefix="${sample%_R1.fastq.gz}"

    echo ""
    echo "[${sample_num}/${n_samples}] Processing: ${prefix}"
    echo "  Started: $(date)"

    # Step 1: Trim reads for quality and TruSeq adapters
    echo "  Trimming reads..."
    bbduk.sh \
        in1="${SCRATCH_WORK_DIR}/input/${prefix}_R1.fastq.gz" \
        in2="${SCRATCH_WORK_DIR}/input/${prefix}_R2.fastq.gz" \
        out1="${SCRATCH_WORK_DIR}/trimmed/${prefix}_trimmed_R1.fastq.gz" \
        out2="${SCRATCH_WORK_DIR}/trimmed/${prefix}_trimmed_R2.fastq.gz" \
        literal=AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
        overwrite=true \
        ftl=20 \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        tpe \
        tbo \
        2>&1 | tail -5

    # Step 2: Merge paired-end reads
    echo "  Merging reads..."
    bbmerge.sh \
        in1="${SCRATCH_WORK_DIR}/trimmed/${prefix}_trimmed_R1.fastq.gz" \
        in2="${SCRATCH_WORK_DIR}/trimmed/${prefix}_trimmed_R2.fastq.gz" \
        out="${SCRATCH_WORK_DIR}/merged/${prefix}_trimmed_merged.fastq" \
        outu="${SCRATCH_WORK_DIR}/unmerged/${prefix}_trimmed_unmerged_R1.fastq" \
        outu2="${SCRATCH_WORK_DIR}/unmerged/${prefix}_trimmed_unmerged_R2.fastq" \
        overwrite=true \
        2>&1 | tail -5

    # Step 3: Collapse identical reads
    echo "  Collapsing identical reads..."
    python3 "${COLLAPSE_SCRIPT}" \
        "${SCRATCH_WORK_DIR}/merged/${prefix}_trimmed_merged.fastq" \
        "${SCRATCH_WORK_DIR}/collapsed/${prefix}_collapsed.fastq" \
        1

    echo "  Completed: $(date)"
done

# ============================================================================
# Copy final outputs from Scratch back to Oak
# ============================================================================

echo ""
echo "=============================================="
echo "Copying results from Scratch to Oak..."
echo "  From: ${SCRATCH_WORK_DIR}/collapsed/"
echo "  To: ${PROCESSED_DATA_DIR}/collapsed_fastq/"

cp "${SCRATCH_WORK_DIR}/collapsed/"*.fastq "${PROCESSED_DATA_DIR}/collapsed_fastq/"

n_output=$(ls -1 "${PROCESSED_DATA_DIR}/collapsed_fastq/"*.fastq 2>/dev/null | wc -l)
echo "  Copied ${n_output} collapsed FASTQ files"

# ============================================================================
# Cleanup scratch (optional - scratch auto-purges after 90 days)
# ============================================================================

echo ""
echo "Cleaning up scratch working directory..."
rm -rf "${SCRATCH_WORK_DIR}"

echo ""
echo "=============================================="
echo "Step 01 Complete"
echo "Processed ${sample_num} samples"
echo "Output: ${PROCESSED_DATA_DIR}/collapsed_fastq/"
echo "Date: $(date)"
echo "=============================================="
