#!/bin/bash
#SBATCH --partition=mylab
#SBATCH --account=mylab
#SBATCH --job-name=bc1_pipeline
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

# BC1-From-Screen Full Pipeline Runner
# Author: Kevin R. Roy
# Date: 2026-02-08
#
# Usage:
#   # Run via SLURM
#   sbatch run_full_pipeline.sh --project 20250912 --library yL437 --t0-mode screen_date
#
#   # Run locally
#   bash run_full_pipeline.sh --project 20250912 --library yL437 --t0-mode screen_date
#
# Arguments:
#   --project: Project name (20240806 or 20250912)
#   --library: Library name (yL406, yL437, yL442)
#   --t0-mode: T0 consolidation mode (screen_date, pseudo_replicates, none)
#   --t0-pseudo-reps: Number of pseudo-replicates (for pseudo_replicates mode)
#   --start-step: Step to start from (default: 03b)
#   --end-step: Step to end at (default: 14)
#   --dry-run: Print commands without executing

set -e

# Default values
START_STEP="03b"
END_STEP="14"
DRY_RUN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --project) PROJECT="$2"; shift 2 ;;
        --library) LIBRARY="$2"; shift 2 ;;
        --t0-mode) T0_MODE="$2"; shift 2 ;;
        --t0-pseudo-reps) T0_PSEUDO_REPS="$2"; shift 2 ;;
        --start-step) START_STEP="$2"; shift 2 ;;
        --end-step) END_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "$PROJECT" || -z "$LIBRARY" || -z "$T0_MODE" ]]; then
    echo "Usage: $0 --project <project> --library <library> --t0-mode <mode> [options]"
    echo "  Required:"
    echo "    --project     Project name (20240806 or 20250912)"
    echo "    --library     Library name (yL406, yL437, yL442)"
    echo "    --t0-mode     T0 consolidation mode (screen_date, pseudo_replicates, none)"
    echo "  Optional:"
    echo "    --t0-pseudo-reps  Number of pseudo-replicates (for pseudo_replicates mode)"
    echo "    --start-step      Step to start from (default: 03b)"
    echo "    --end-step        Step to end at (default: 14)"
    echo "    --dry-run         Print commands without executing"
    exit 1
fi

# Set up paths
BASE_DIR="/path/to"
COMMON_SCRIPTS="$BASE_DIR/common/scripts/bc1_from_screen_pipeline/snakemake"
PYTHON="/path/to/anaconda3/envs/snakemake/bin/python"

if [[ "$PROJECT" == "20240806" ]]; then
    PROJECT_DIR="$BASE_DIR/projects/QTL/20240806_yL406_SpG_QTL_screen"
else
    PROJECT_DIR="$BASE_DIR/projects/QTL/20250912_SpG_QTL_screen"
fi

DATA_DIR="$PROJECT_DIR/processed_data/bc1_from_screen"
KEYFILE_DIR="$PROJECT_DIR/keyfiles/bc1_from_screen"
BC1_REF="$KEYFILE_DIR/bc1_reference_table.tsv"

# Set thread limits for SLURM compatibility
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export LOKY_MAX_CPU_COUNT=${SLURM_CPUS_PER_TASK:-16}

echo "=============================================="
echo "BC1-From-Screen Pipeline"
echo "=============================================="
echo "Project: $PROJECT"
echo "Library: $LIBRARY"
echo "T0 Mode: $T0_MODE"
echo "T0 Pseudo-Reps: ${T0_PSEUDO_REPS:-N/A}"
echo "Start Step: $START_STEP"
echo "End Step: $END_STEP"
echo "Dry Run: $DRY_RUN"
echo "=============================================="

# Function to run or print command
run_cmd() {
    echo ""
    echo ">>> $*"
    if [[ "$DRY_RUN" == "false" ]]; then
        "$@"
    fi
}

# Create output directories
mkdir -p "$DATA_DIR/bc1_counts"
mkdir -p "$DATA_DIR/error_detection"
mkdir -p "$DATA_DIR/count_matrix"
mkdir -p "$DATA_DIR/counts_for_deseq/$LIBRARY"
mkdir -p "$DATA_DIR/DESeq2/$LIBRARY"
mkdir -p logs

# Define file paths
LONG_COUNTS="$DATA_DIR/bc1_counts/${PROJECT}_*_bc1_counts_long_format.tsv"
LONG_COUNTS=$(ls $LONG_COUNTS 2>/dev/null | head -1)

if [[ -z "$LONG_COUNTS" ]]; then
    echo "ERROR: Long-format counts file not found"
    exit 1
fi
echo "Long counts file: $LONG_COUNTS"

# Step 03b: Error Detection
if [[ "$START_STEP" == "03b" || "$START_STEP" < "03b" ]]; then
    echo ""
    echo "=== Step 03b: Error Detection ==="
    run_cmd $PYTHON "$COMMON_SCRIPTS/03b_error_detection.py" \
        --counts "$LONG_COUNTS" \
        --output "$DATA_DIR/error_detection/${LIBRARY}_error_bc1_mapping.tsv" \
        --stats-output "$DATA_DIR/error_detection/${LIBRARY}_bc1_statistics.tsv" \
        --max-error-count 50 \
        --max-error-samples 3 \
        --min-count-ratio 100 \
        --max-hamming 2
fi

# Step 03c: Absorb Error Counts
if [[ "$START_STEP" <= "03c" && "$END_STEP" >= "03c" ]]; then
    echo ""
    echo "=== Step 03c: Absorb Error Counts ==="
    run_cmd $PYTHON "$COMMON_SCRIPTS/03c_absorb_error_counts.py" \
        --counts "$LONG_COUNTS" \
        --error-mappings "$DATA_DIR/error_detection/${LIBRARY}_error_bc1_mapping.tsv" \
        --output "$DATA_DIR/bc1_counts/${LIBRARY}_bc1_counts_long_format_cleaned.tsv" \
        --log-output "$DATA_DIR/error_detection/${LIBRARY}_absorption_log.tsv"
fi

# Step 04: Tier Classification
if [[ "$START_STEP" <= "04" && "$END_STEP" >= "04" ]]; then
    echo ""
    echo "=== Step 04: Tier Classification ==="
    CLEANED_COUNTS="$DATA_DIR/bc1_counts/${LIBRARY}_bc1_counts_long_format_cleaned.tsv"
    run_cmd $PYTHON "$COMMON_SCRIPTS/04_tier_classification.py" \
        --bc1-counts "$CLEANED_COUNTS" \
        --bc1-reference "$BC1_REF" \
        --library "$LIBRARY" \
        --output "$DATA_DIR/bc1_counts/${LIBRARY}_tier_mapping.tsv" \
        --stats-output "$DATA_DIR/bc1_counts/${LIBRARY}_tier_stats.tsv"
fi

# Step 05: Create Count Matrix
if [[ "$START_STEP" <= "05" && "$END_STEP" >= "05" ]]; then
    echo ""
    echo "=== Step 05: Create Count Matrix ==="
    CLEANED_COUNTS="$DATA_DIR/bc1_counts/${LIBRARY}_bc1_counts_long_format_cleaned.tsv"
    TIER_MAPPING="$DATA_DIR/bc1_counts/${LIBRARY}_tier_mapping.tsv"
    SAMPLE_KEY="$KEYFILE_DIR/combined_key.tsv"

    # Handle different file naming for yL406
    if [[ "$PROJECT" == "20240806" ]]; then
        SAMPLE_KEY="$KEYFILE_DIR/20240806_yL406_SpG_QTL_screen_combined_key.tsv"
    fi

    T0_ARGS="--t0-consolidation-mode $T0_MODE"
    if [[ "$T0_MODE" == "pseudo_replicates" && -n "$T0_PSEUDO_REPS" ]]; then
        T0_ARGS="$T0_ARGS --t0-pseudo-rep-count $T0_PSEUDO_REPS"
    fi

    run_cmd $PYTHON "$COMMON_SCRIPTS/05_create_count_matrix.py" \
        --bc1-counts "$CLEANED_COUNTS" \
        --tier-mapping "$TIER_MAPPING" \
        --sample-key "$SAMPLE_KEY" \
        --library "$LIBRARY" \
        --output "$DATA_DIR/count_matrix/${LIBRARY}_count_matrix.tsv" \
        --metadata-output "$DATA_DIR/count_matrix/${LIBRARY}_t0_consolidation_metadata.tsv" \
        $T0_ARGS
fi

# Step 06: Combine Technical Replicates
if [[ "$START_STEP" <= "06" && "$END_STEP" >= "06" ]]; then
    echo ""
    echo "=== Step 06: Combine Technical Replicates ==="
    run_cmd $PYTHON "$COMMON_SCRIPTS/06_combine_technical_replicates.py" \
        --count-matrix "$DATA_DIR/count_matrix/${LIBRARY}_count_matrix.tsv" \
        --sample-key "$SAMPLE_KEY" \
        --library "$LIBRARY" \
        --output "$DATA_DIR/count_matrix/${LIBRARY}_count_matrix_combined.tsv" \
        --qc-output "$DATA_DIR/qc/${LIBRARY}_technical_rep_correlation.tsv"
fi

# Step 07: Split Comparisons
if [[ "$START_STEP" <= "07" && "$END_STEP" >= "07" ]]; then
    echo ""
    echo "=== Step 07: Split Comparisons ==="
    run_cmd $PYTHON "$COMMON_SCRIPTS/07_split_comparisons.py" \
        --count-matrix "$DATA_DIR/count_matrix/${LIBRARY}_count_matrix_combined.tsv" \
        --sample-key "$SAMPLE_KEY" \
        --library "$LIBRARY" \
        --output-dir "$DATA_DIR/counts_for_deseq/$LIBRARY"
fi

# Step 08: DESeq2 (via SLURM job array)
if [[ "$START_STEP" <= "08" && "$END_STEP" >= "08" ]]; then
    echo ""
    echo "=== Step 08: DESeq2 ==="
    echo "Submit DESeq2 job array separately using:"
    echo "  sbatch $BASE_DIR/common/scripts/bc1_from_screen_pipeline/slurm/run_deseq2_array.sbatch \\"
    echo "    --library $LIBRARY \\"
    echo "    --project-dir $PROJECT_DIR"
    echo ""
    echo "Skipping automatic DESeq2 submission (requires SLURM job array)"
fi

# Steps 09-14 would continue here...
echo ""
echo "=== Pipeline steps 03b-07 complete ==="
echo ""
echo "Next steps:"
echo "1. Submit DESeq2 job array"
echo "2. Run steps 09-14 after DESeq2 completes"

echo ""
echo "Done!"
