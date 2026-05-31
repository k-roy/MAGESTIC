#!/bin/bash
#SBATCH --job-name=unmatched_annot
#SBATCH --partition=mylab
#SBATCH --account=mylab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Unmatched Donor Annotation Pipeline Template
#
# This script runs the unmatched donor annotation pipeline (06) followed
# by integration (07).
#
# USAGE:
# 1. Copy this template to your project scripts directory
# 2. Edit the PROJECT CONFIGURATION section below
# 3. Submit with: sbatch run_unmatched_donor_annotation.sh
#
# Author: Kevin R. Roy

set -e

# =============================================================================
# PROJECT CONFIGURATION - EDIT THIS SECTION
# =============================================================================

# Base directories
BASE_DIR="/path/to"
PROJECT_DIR="${BASE_DIR}/projects/QTL/YOUR_PROJECT_NAME"  # <-- EDIT THIS

# Input bc1 reference table (output from 05)
BC1_REFERENCE="${PROJECT_DIR}/processed_data/bc1_donor_bc0/final/bc1_reference_table.tsv"  # <-- EDIT THIS

# Output directory for unmatched annotations
OUTPUT_DIR="${PROJECT_DIR}/processed_data/bc1_donor_bc0/unmatched_annotation"

# Oligo design files (adjust based on which libraries you're using)
SPG_OLIGO_DESIGN="${BASE_DIR}/common/oligo_designs/20240411_Twist_200mer_oligo_array_order.tsv"
SPCAS9_OLIGO_DESIGN="${BASE_DIR}/common/oligo_designs/20210422_Twist_200mer_oligo_array_order.tsv"

# Reference files
REFERENCE_GENOME="${BASE_DIR}/common/reference_genomes/MAGESTIC_background_strain.fasta"
ANNOTATION_GFF="${BASE_DIR}/common/annotation_files/MAGESTIC_background_strain_annotations.gff"

# Optional: Harmonized oligo key for delta variant computation
# Set to empty string if not available
HARMONIZED_OLIGO_KEY="${BASE_DIR}/common/annotation_files/20210422_and_20240411_Bloom_et_al_16_strains_QTL_harmonized_designed_variant_oligos.tsv"

# Pipeline scripts location
PIPELINE_DIR="${BASE_DIR}/common/scripts/bc1_donor_bc0_pipeline/snakemake"

# =============================================================================
# THREAD LIMITING - DO NOT EDIT
# =============================================================================

# Match thread limits to SLURM allocation
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK

# =============================================================================
# JOB INFO
# =============================================================================

echo "======================================"
echo "UNMATCHED DONOR ANNOTATION PIPELINE"
echo "======================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo ""
echo "Project: ${PROJECT_DIR}"
echo "Input: ${BC1_REFERENCE}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

# Activate conda environment
source ~/.bashrc
conda activate base

# Verify required packages
echo "Checking required packages..."
python -c "import pandas, numpy, pysam, Levenshtein, tqdm; print('All required packages available')" || {
    echo "Error: Missing required packages. Please install: pysam, Levenshtein, tqdm"
    exit 1
}

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# =============================================================================
# STEP G: ANNOTATE UNMATCHED DONORS
# =============================================================================

echo ""
echo "======================================"
echo "STEP G: Annotate Unmatched Donors"
echo "======================================"

# Build command with optional harmonized oligo key
CMD="python ${PIPELINE_DIR}/06_annotate_unmatched_donors.py \
    --input ${BC1_REFERENCE} \
    --output-dir ${OUTPUT_DIR} \
    --spg-oligo-design ${SPG_OLIGO_DESIGN} \
    --spcas9-oligo-design ${SPCAS9_OLIGO_DESIGN} \
    --reference-genome ${REFERENCE_GENOME} \
    --annotation-gff ${ANNOTATION_GFF} \
    --n-workers $SLURM_CPUS_PER_TASK \
    --n-chunks $SLURM_CPUS_PER_TASK"

if [[ -n "${HARMONIZED_OLIGO_KEY}" && -f "${HARMONIZED_OLIGO_KEY}" ]]; then
    CMD="$CMD --harmonized-oligo-key ${HARMONIZED_OLIGO_KEY}"
fi

echo "Running: $CMD"
$CMD

# =============================================================================
# STEP H: INTEGRATE ANNOTATIONS
# =============================================================================

echo ""
echo "======================================"
echo "STEP H: Integrate Unmatched Annotations"
echo "======================================"

ANNOTATIONS_FILE="${OUTPUT_DIR}/unmatched_donor_annotations.tsv"
OUTPUT_FILE="${OUTPUT_DIR}/../final/bc1_reference_table_with_unmatched_annotations.tsv"

if [[ -f "${ANNOTATIONS_FILE}" ]]; then
    python ${PIPELINE_DIR}/07_integrate_unmatched_annotations.py \
        --bc1-reference ${BC1_REFERENCE} \
        --unmatched-annotations ${ANNOTATIONS_FILE} \
        --output ${OUTPUT_FILE}

    echo ""
    echo "Final output: ${OUTPUT_FILE}"
else
    echo "Warning: Annotations file not found: ${ANNOTATIONS_FILE}"
    echo "Skipping integration step."
fi

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "======================================"
echo "PIPELINE COMPLETE"
echo "======================================"
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  Annotations: ${OUTPUT_DIR}/unmatched_donor_annotations.tsv"
echo "  Integrated: ${OUTPUT_FILE}"
