#!/bin/bash
#SBATCH --job-name=magestic-wgs-bc1
#SBATCH --account=mylab
#SBATCH --partition=mylab,owners
#SBATCH --array=1-288%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/wgs_bc1_%A_%a.out
#SBATCH --error=logs/wgs_bc1_%A_%a.err
#
# Stage 1 (bc1 calling) as an idempotent SLURM array — one SAMPLE per array task
# (ARCHITECTURE_POLICY §6). Submit with:
#
#   CONFIG=/path/to/<cohort>_wgs.yaml MANIFEST=/path/to/sample_manifest.tsv \
#     sbatch --array=1-$(( $(wc -l < "$MANIFEST") - 1 )) template.sbatch
#
# 🔴 This template deliberately carries NO CPU-generation node filter and no node
# exclusion, and none may be added. Two older repo files still pin an AVX-512 CPU
# generation list as a "SIGILL guard"; that advice was RETRACTED 2026-08-06. The lab
# partition is exactly two Milan nodes, so such a filter matches NO node in it — the
# partition hard-fails, and a `<lab>,owners` submission is SILENTLY EVICTED to the
# preemptible owners queue. A repo guard greps this file for those flags, so the flag
# names themselves are not written out even in this warning.

# Source ~/.bashrc BEFORE `set -u`: /etc/profile.d/256term.sh has an unbound COLORTERM
# that kills the job in ~3 s with no useful log.
source ~/.bashrc
set -euo pipefail

: "${CONFIG:?CONFIG (path to the cohort YAML) must be set}"
: "${MANIFEST:?MANIFEST (path to the sample manifest TSV) must be set}"

# Single-threaded per task; the fan-out is the array, not the BLAS pool. Exported BEFORE
# any Python import, or numpy/pandas grab every core on the node.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Row 1 of the manifest is the header, so task N reads line N+1.
LINE_NUM=$(( SLURM_ARRAY_TASK_ID + 1 ))
SAMPLE=$(awk -v n="$LINE_NUM" 'NR==n {print $1}' "$MANIFEST")
if [ -z "$SAMPLE" ]; then
    echo "no manifest row for array task ${SLURM_ARRAY_TASK_ID} (line ${LINE_NUM})" >&2
    exit 1
fi

# §6 skip-check: every artifact is written atomically, so a file that EXISTS is a file
# that is COMPLETE and re-running a requeued task is free. Checking the output rather
# than a sentinel is deliberate — a sentinel written after a loop is not a completion
# signal.
OUTPUT_DIR="${PROCESSED_DIR:-processed_data}/bc1_calls"
if [ -f "${OUTPUT_DIR}/${SAMPLE}.bc1_calls.tsv" ] && \
   [ -f "${OUTPUT_DIR}/${SAMPLE}.bc1_tracks.tsv" ]; then
    echo "${SAMPLE}: output already complete, skipping"
    exit 0
fi

# The SAME console script the Snakemake path calls, so there is one entry point per
# stage and dual-path parity stays structural (§4).
magestic-wgs-bc1 --config "$CONFIG" --sample "$SAMPLE"
