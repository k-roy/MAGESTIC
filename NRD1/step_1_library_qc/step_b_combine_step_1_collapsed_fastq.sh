#!/bin/bash
#SBATCH --job-name=combine_step_1_collapsed_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=mylab,owners
#SBATCH --account=mylab
#SBATCH --time=24:00:00

PYTHON="/path/to/anaconda3/bin/python"
SCRIPT_DIR="/path/to/software/MAGESTIC/NRD1/step_1_library_qc"

echo "Python: $($PYTHON --version)"
$PYTHON "${SCRIPT_DIR}/step_b_combine_step_1_collapsed_fastq.py"
