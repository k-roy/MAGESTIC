#!/bin/bash
#SBATCH --job-name=calculate_bc0_purity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=mylab,owners
#SBATCH --account=mylab
#SBATCH --time=48:00:00

BASE_DIR="/path/to"
PROJECT_DIR="${BASE_DIR}/projects/NNS/20250628_repeat_step_1_library_cloning"
PYTHON="/path/to/anaconda3/bin/python"
SCRIPT="${BASE_DIR}/software/MAGESTIC/NRD1/guide_donor_bc0_linking/step_d_calculate_bc0_purity.py"

PROCESSED_DATA_DIR="${PROJECT_DIR}/processed_data"
GUIDE_DONOR_BC0_DIR="${PROCESSED_DATA_DIR}/guide_donor_bc0"
PURITY_FILTERED_DIR="${PROCESSED_DATA_DIR}/guide_donor_bc0_purity_filtered"
mkdir -p "$PURITY_FILTERED_DIR"

echo "Python: $($PYTHON --version)"

cd "$GUIDE_DONOR_BC0_DIR"
for sample in V*_matched_guide_donor_bc0_counts.tsv; do
    prefix=${sample%_matched_guide_donor_bc0_counts.tsv}
    echo "Processing: ${prefix}"
    $PYTHON "$SCRIPT" \
        "${GUIDE_DONOR_BC0_DIR}/${prefix}_matched_guide_donor_bc0_counts.tsv" \
        "${PURITY_FILTERED_DIR}/${prefix}_guide_donor_bc0_purity_filtered.tsv"
done
