#!/bin/bash

# See `man sbatch` or https://slurm.schedmd.com/sbatch.html for descriptions
# of sbatch options.

#SBATCH --job-name=step_b_process_guide_donor_bc0_collapsed_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=mylab
#SBATCH --account=mylab
#SBATCH --time=48:00:00

PARENT_DIR="/path/to/"
PROJECT_PATH="by_project/NNS/202250628_repeat_step_1_library_cloning/"
# scripts_and_keyfiles/by_project/NNS/202250628_repeat_step_1_library_cloning/scripts/step_d_calculate_bc0_purity.py
SCRIPT=$PARENT_DIR"scripts_and_keyfiles/"$PROJECT_PATH"scripts/step_d_calculate_bc0_purity.py"
PROCESSED_DATA_DIR=$PARENT_DIR"processed_data/"$PROJECT_PATH
PROCESSED_GUIDE_DONOR_BC0_DIR=$PROCESSED_DATA_DIR"guide_donor_bc0/"
GUIDE_DONOR_BC0_PURITY_FILTERED_DIR=$PROCESSED_DATA_DIR"guide_donor_bc0_purity_filtered/"
mkdir $GUIDE_DONOR_BC0_PURITY_FILTERED_DIR

cd $PROCESSED_GUIDE_DONOR_BC0_DIR
for sample in V*_matched_guide_donor_bc0_counts.tsv; 
do prefix=${sample%_matched_guide_donor_bc0_counts.tsv};
echo $prefix; 
python3 $SCRIPT $PROCESSED_GUIDE_DONOR_BC0_DIR$prefix"_matched_guide_donor_bc0_counts.tsv" $GUIDE_DONOR_BC0_PURITY_FILTERED_DIR$prefix"_guide_donor_bc0_purity_filtered.tsv" 
done