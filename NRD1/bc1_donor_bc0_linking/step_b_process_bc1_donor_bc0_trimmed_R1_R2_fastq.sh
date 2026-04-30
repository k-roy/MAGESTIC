#!/bin/bash

PARENT_DIR="/path/to/"
PROJECT_FOLDER="by_project/NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/"
SCRIPTS_DIR=$PARENT_DIR"scripts_and_keyfiles/"$PROJECT_FOLDER"scripts/"

INTERMEDIATE_DATA_DIR=$PARENT_DIR"intermediate_data/$PROJECT_FOLDER/"
TRIMMED_READS_DIR=$INTERMEDIATE_DATA_DIR"trimmed_fastq/"

PROCESSED_DATA_DIR=$PARENT_DIR"processed_data/"$PROJECT_FOLDER
PROCESSED_BC1_DONOR_BC0_DIR=$PROCESSED_DATA_DIR"bc1_donor_bc0_fragment_counts/"

mkdir -p $PROCESSED_BC1_DONOR_BC0_DIR

OVERWRITE=TRUE
cd $TRIMMED_READS_DIR

# Since the Python script does not use multiprocessing, 
# we can process them in parallel here with GNU parallel.
export SCRIPTS_DIR
export TRIMMED_READS_DIR
export PROCESSED_BC1_DONOR_BC0_DIR
export OVERWRITE

ls *_trimmed_R1.fastq.gz | sed 's/_trimmed_R1.fastq.gz$//' | \
    parallel --bar "
        prefix={}
        if [ -f \"\$PROCESSED_BC1_DONOR_BC0_DIR\${prefix}_processed_bc1_donor_bc0.tsv\" ] && [ \"\$OVERWRITE\" = \"FALSE\" ]; then
            echo \"Skipping \$prefix, already processed.\"
        else
            echo \"Processing \$prefix\"
            python3 \"\$SCRIPTS_DIR/step_b_process_bc1_donor_bc0_trimmed_R1_R2_fastq.py\" \
                \"\$TRIMMED_READS_DIR\${prefix}_trimmed_R1.fastq.gz\" \
                \"\$TRIMMED_READS_DIR\${prefix}_trimmed_R2.fastq.gz\" \
                \"\$PROCESSED_BC1_DONOR_BC0_DIR\${prefix}_bc1_donor_bc0_fragment_counts.tsv\" \
                \$OVERWRITE
        fi
    "
