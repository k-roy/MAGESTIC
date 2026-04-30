#!/bin/bash

# See `man sbatch` or https://slurm.schedmd.com/sbatch.html for descriptions
# of sbatch options.
#SBATCH --job-name=combine_step_1_collapsed_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=mylab
#SBATCH --account=mylab
#SBATCH --time=24:00:00

python combine_step_1_collapsed_fastq.py
python map_guide_donor_bc0_to_designed_oligos.py