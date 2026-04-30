#!/bin/bash

# See `man sbatch` or https://slurm.schedmd.com/sbatch.html for descriptions
# of sbatch options.
#SBATCH --job-name=trim_merge_collapse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=mylab
#SBATCH --account=mylab
#SBATCH --time=2:00:00

PARENT_DIR="/path/to/"
mkdir $PARENT_DIR"intermediate_data/by_project/NNS/"
mkdir $PARENT_DIR"processed_data/by_project/NNS/"

PROJECT_FOLDER="20250813_NRD1_MAGESTIC_bc0_donor_bc1"

RAW_DATA_DIR=$PARENT_DIR"raw_data/by_project/NNS/$PROJECT_FOLDER/"
SCRIPTS_DIR=$PARENT_DIR"scripts_and_keyfiles/general_scripts/"

INTERMEDIATE_DATA_DIR=$PARENT_DIR"intermediate_data/by_project/NNS/$PROJECT_FOLDER/"
mkdir $INTERMEDIATE_DATA_DIR
TRIMMED_READS_DIR=$INTERMEDIATE_DATA_DIR"trimmed_fastq/"
mkdir $TRIMMED_READS_DIR
MERGED_READS_DIR=$INTERMEDIATE_DATA_DIR"merged_fastq/"
mkdir $MERGED_READS_DIR
UNMERGED_READS_DIR=$INTERMEDIATE_DATA_DIR"unmerged_fastq/"
mkdir $UNMERGED_READS_DIR

PROCESSED_DATA_DIR=$PARENT_DIR"processed_data/by_project/NNS/$PROJECT_FOLDER/"
mkdir $PROCESSED_DATA_DIR
COLLAPSED_READS_DIR=$PROCESSED_DATA_DIR"collapsed_fastq/"
mkdir $COLLAPSED_READS_DIR

## trim reads for quality and TruSeq adapters
FASTQ_DIR=$RAW_DATA_DIR"/fastq/"
cd $FASTQ_DIR
for sample in 20250815_plate_1C_*_R1.fastq.gz; 
do prefix=${sample%_R1.fastq.gz};
echo $prefix; 
OVERWRITE=false
bbduk.sh in1=$FASTQ_DIR$prefix"_R1.fastq.gz" in2=$FASTQ_DIR$prefix"_R2.fastq.gz" out1=$TRIMMED_READS_DIR$prefix"_trimmed_R1.fastq.gz" out2=$TRIMMED_READS_DIR$prefix"_trimmed_R2.fastq.gz" literal=AGATCGGAAGAGCGTCGTGTAGGGAAAGA overwrite=$OVERWRITE ftl=20 ktrim=r k=23 mink=11 hdist=1 tpe tbo; 
bbmerge.sh in1=$TRIMMED_READS_DIR$prefix"_trimmed_R1.fastq.gz" in2=$TRIMMED_READS_DIR$prefix"_trimmed_R2.fastq.gz" out=$MERGED_READS_DIR$prefix"_trimmed_merged.fastq" outu=$UNMERGED_READS_DIR$prefix"_trimmed_unmerged_R1.fastq" outu2=$UNMERGED_READS_DIR$prefix"_trimmed_unmerged_R2.fastq" overwrite=$OVERWRITE; 
## collapse identical reads and change read name to read counts
python3 $SCRIPTS_DIR"collapse_identical_merged_reads_min_counts.py" $MERGED_READS_DIR$prefix"_trimmed_merged.fastq" $COLLAPSED_READS_DIR$prefix"_collapsed.fastq" 1
done

