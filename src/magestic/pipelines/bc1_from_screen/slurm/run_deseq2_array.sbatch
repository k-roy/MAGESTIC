#!/bin/bash
#SBATCH --job-name=deseq2_array
#SBATCH --partition=mylab
#SBATCH --account=mylab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=logs/deseq2_%A_%a.out
#SBATCH --error=logs/deseq2_%A_%a.err

# =============================================================================
# SLURM Job Array for DESeq2 Analysis
# =============================================================================
# Safe execution of DESeq2 on multiple comparisons using job arrays.
# Each array task processes ONE comparison with controlled resources.
#
# Usage:
#   # First, create comparison list
#   cut -f2 comparison_manifest.tsv | tail -n +2 > comparison_ids.txt
#
#   # Submit array job
#   sbatch --array=1-$(wc -l < comparison_ids.txt) run_deseq2_array.sbatch
#
# Required environment variables (set before sbatch or in this script):
#   PROJECT_DIR: Path to project directory
#   LIBRARY: Library name (e.g., yL406)
#
# Author: Kevin R. Roy
# =============================================================================

# CRITICAL: Limit threads BEFORE any Python imports
# This prevents joblib/numpy from spawning extra threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export LOKY_MAX_CPU_COUNT=1

# Configuration - override these via environment or edit here
PROJECT_DIR="${PROJECT_DIR:-/path/to/projects/QTL/20240806_yL406_SpG_QTL_screen}"
LIBRARY="${LIBRARY:-yL406}"

# Paths
CONDA_ENV="/path/to/anaconda3/envs/pydeseq2"
SCRIPT="/path/to/common/scripts/bc1_from_screen_pipeline/snakemake/08_deseq2_single.py"
COUNTS_DIR="${PROJECT_DIR}/processed_data/bc1_from_screen/counts_for_deseq/${LIBRARY}"
OUTPUT_DIR="${PROJECT_DIR}/processed_data/bc1_from_screen/DESeq2/${LIBRARY}"
MANIFEST="${COUNTS_DIR}/comparison_manifest.tsv"

# Create output and log directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "logs"

# Get comparison info for this array task
# SLURM_ARRAY_TASK_ID is 1-indexed
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))  # +1 to skip header
COMPARISON_INFO=$(sed -n "${LINE_NUM}p" "${MANIFEST}")

if [ -z "$COMPARISON_INFO" ]; then
    echo "ERROR: No comparison found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Parse manifest columns (tab-separated)
COMPARISON_ID=$(echo "$COMPARISON_INFO" | cut -f2)
REFERENCE_FILE=$(echo "$COMPARISON_INFO" | cut -f6)
TREATMENT_FILE=$(echo "$COMPARISON_INFO" | cut -f7)

echo "=========================================="
echo "DESeq2 Analysis - Array Task ${SLURM_ARRAY_TASK_ID}"
echo "=========================================="
echo "Comparison: ${COMPARISON_ID}"
echo "Reference: ${REFERENCE_FILE}"
echo "Treatment: ${TREATMENT_FILE}"
echo "Library: ${LIBRARY}"
echo "Output: ${OUTPUT_DIR}/${COMPARISON_ID}_deseq2_results.tsv"
echo ""

# Check if already completed
OUTPUT_FILE="${OUTPUT_DIR}/${COMPARISON_ID}_deseq2_results.tsv"
if [ -f "${OUTPUT_FILE}" ]; then
    echo "Output already exists, skipping"
    exit 0
fi

# Activate conda environment
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate pydeseq2

# Run DESeq2
${CONDA_ENV}/bin/python ${SCRIPT} \
    --reference "${COUNTS_DIR}/reference/${REFERENCE_FILE}" \
    --treatment "${COUNTS_DIR}/treatment/${TREATMENT_FILE}" \
    --output "${OUTPUT_FILE}" \
    --comparison "${COMPARISON_ID}" \
    --library "${LIBRARY}"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: ${COMPARISON_ID}"
else
    echo "FAILED: ${COMPARISON_ID} (exit code: ${EXIT_CODE})"
fi

exit $EXIT_CODE
