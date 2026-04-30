#!/bin/bash
#SBATCH --job-name=map_guide_donor_bc0_to_designed_oligos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=mylab,owners
#SBATCH --account=mylab
#SBATCH --time=24:00:00

PYTHON="/path/to/anaconda3/bin/python"
SCRIPT_DIR="/path/to/software/MAGESTIC/NRD1/guide_donor_bc0_linking"

echo "Python: $($PYTHON --version)"
$PYTHON "${SCRIPT_DIR}/step_c_map_guide_donor_bc0_to_designed_oligos.py"
