#!/bin/bash
#SBATCH --job-name=trim_merge_collapse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=mylab,owners
#SBATCH --account=mylab
#SBATCH --time=24:00:00

# Trim, merge, and collapse V625 guide-donor-bc0 FASTQ reads.
# Part of the guide → donor-bc0 → bc1 linkage workflow.
#
# Inputs:  raw_data/fastq/20250707_plate_5A_sample_{37-40}_R{1,2}.fastq.gz
# Outputs: processed_data/collapsed_fastq/*_collapsed.fastq
#          (intermediate trimmed/merged reads written to $SCRATCH, not Oak)

BASE_DIR="/path/to"
PROJECT_DIR="${BASE_DIR}/projects/NNS/20250628_repeat_step_1_library_cloning"
PYTHON="/path/to/anaconda3/bin/python"

FASTQ_DIR="${PROJECT_DIR}/raw_data/fastq"
COLLAPSED_READS_DIR="${PROJECT_DIR}/processed_data/collapsed_fastq"
mkdir -p "$COLLAPSED_READS_DIR"

# collapse utility — present in the MAGESTIC repo root
COLLAPSE_SCRIPT="${BASE_DIR}/software/MAGESTIC/collapse_identical_merged_reads_min_counts.py"

# Intermediate trimmed/merged reads go to $SCRATCH (high-bandwidth, auto-purged)
SCRATCH_DIR="${SCRATCH}/V625_trim_merge_$(date +%Y%m%d)"
TRIMMED_READS_DIR="${SCRATCH_DIR}/trimmed_fastq"
MERGED_READS_DIR="${SCRATCH_DIR}/merged_fastq"
UNMERGED_READS_DIR="${SCRATCH_DIR}/unmerged_fastq"
mkdir -p "$TRIMMED_READS_DIR" "$MERGED_READS_DIR" "$UNMERGED_READS_DIR"

echo "PROJECT_DIR:       ${PROJECT_DIR}"
echo "FASTQ_DIR:         ${FASTQ_DIR}"
echo "COLLAPSED_OUT:     ${COLLAPSED_READS_DIR}"
echo "SCRATCH_DIR:       ${SCRATCH_DIR}"
echo ""

cd "$FASTQ_DIR"
for sample in *_R1.fastq.gz; do
    prefix=${sample%_R1.fastq.gz}
    echo "Processing: ${prefix}"

    bbduk.sh \
        in1="${FASTQ_DIR}/${prefix}_R1.fastq.gz" \
        in2="${FASTQ_DIR}/${prefix}_R2.fastq.gz" \
        out1="${TRIMMED_READS_DIR}/${prefix}_trimmed_R1.fastq.gz" \
        out2="${TRIMMED_READS_DIR}/${prefix}_trimmed_R2.fastq.gz" \
        literal=AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
        overwrite=true ftl=20 ktrim=r k=23 mink=11 hdist=1 tpe tbo

    bbmerge.sh \
        in1="${TRIMMED_READS_DIR}/${prefix}_trimmed_R1.fastq.gz" \
        in2="${TRIMMED_READS_DIR}/${prefix}_trimmed_R2.fastq.gz" \
        out="${MERGED_READS_DIR}/${prefix}_trimmed_merged.fastq" \
        outu="${UNMERGED_READS_DIR}/${prefix}_trimmed_unmerged_R1.fastq" \
        outu2="${UNMERGED_READS_DIR}/${prefix}_trimmed_unmerged_R2.fastq" \
        overwrite=true

    $PYTHON "$COLLAPSE_SCRIPT" \
        "${MERGED_READS_DIR}/${prefix}_trimmed_merged.fastq" \
        "${COLLAPSED_READS_DIR}/${prefix}_collapsed.fastq" \
        1
done

echo ""
echo "Done. Collapsed FASTQs in: ${COLLAPSED_READS_DIR}"
echo "Scratch intermediates in:  ${SCRATCH_DIR}  (auto-purged after 90 days)"
