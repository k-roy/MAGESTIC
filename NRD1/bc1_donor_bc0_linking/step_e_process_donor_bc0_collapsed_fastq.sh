#!/bin/bash

PARENT_DIR="/path/to/"
PROJECT_FOLDER="by_project/NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/"
SCRIPTS_DIR=$PARENT_DIR"scripts_and_keyfiles/"$PROJECT_FOLDER"scripts/"

COLLAPSED_READS_DIR=$PROCESSED_DATA_DIR"collapsed_fastq/"

PROCESSED_DATA_DIR=$PARENT_DIR"processed_data/"$PROJECT_FOLDER
PROCESSED_DONOR_BC0_DIR=$PROCESSED_DATA_DIR"complete_donor_bc0_counts/"

mkdir -p $PROCESSED_DONOR_BC0_DIR

OVERWRITE=TRUE
cd $COLLAPSED_READS_DIR

# Since the Python script does not use multiprocessing, 
# we can process them in parallel here with GNU parallel.
export SCRIPTS_DIR
export COLLAPSED_READS_DIR
export PROCESSED_DONOR_BC0_DIR
export OVERWRITE

ls *_collapsed.fastq | sed 's/_collapsed.fastq$//' | \
    parallel --bar "
        prefix={}
        if [ -f \"\$PROCESSED_DONOR_BC0_DIR\${prefix}_complete_donor_bc0_counts.tsv\" ] && [ \"\$OVERWRITE\" = \"FALSE\" ]; then
            echo \"Skipping \$prefix, already processed.\"
        else
            echo \"Processing \$prefix\"
            python3 \"\$SCRIPTS_DIR/step_b_process_donor_bc0_collapsed_fastq.py\" \
                \"\$COLLAPSED_READS_DIR\${prefix}_collapsed.fastq\" \
                \"\$PROCESSED_DONOR_BC0_DIR\${prefix}_complete_donor_bc0_counts.tsv\" \
                \$OVERWRITE
        fi
    "
