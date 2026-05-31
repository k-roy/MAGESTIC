#!/bin/bash
# Activation script for bc1_from_screen pipeline
# Author: Kevin R. Roy
#
# Usage in SLURM scripts:
#   source /path/to/common/scripts/bc1_from_screen_pipeline/activate_env.sh
#
# Or for interactive use:
#   source activate_env.sh

# Conda location
CONDA_BASE="/path/to/anaconda3"
ENV_NAME="bc1_from_screen"

# Initialize conda (required for non-interactive shells)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Activate environment
conda activate "${ENV_NAME}" 2>/dev/null || {
    # Fallback to ihw_env if bc1_from_screen doesn't exist yet
    echo "Warning: ${ENV_NAME} environment not found, falling back to ihw_env"
    conda activate ihw_env
    ENV_NAME="ihw_env"
}

# Set R_HOME to use the environment's R (required for rpy2 + IHW)
export R_HOME="${CONDA_PREFIX}/lib/R"

# Ensure threading matches SLURM allocation (for SLURM jobs)
if [ -n "${SLURM_CPUS_PER_TASK}" ]; then
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
    export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
    export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
    export LOKY_MAX_CPU_COUNT="${SLURM_CPUS_PER_TASK}"
fi

echo "Activated ${ENV_NAME} environment"
echo "  R_HOME=${R_HOME}"
[ -n "${SLURM_CPUS_PER_TASK}" ] && echo "  Thread limit: ${SLURM_CPUS_PER_TASK}"
