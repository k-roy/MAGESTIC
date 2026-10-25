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
#   export CONFIG=/path/to/<cohort>_wgs.yaml MANIFEST=/path/to/sample_manifest.tsv
#   sbatch --array=1-$(python -m magestic.pipelines.wgs.slurm.task_select \
#     --config "$CONFIG" --manifest "$MANIFEST" --count) template.sbatch
#
# Size the array with the selector, never `wc -l`: comment and blank manifest lines
# are not samples, so a line count creates tasks that have no sample.
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

# Task N is the Nth SAMPLE as the shared manifest reader enumerates it (the same one
# the direct and Snakemake paths use): columns by header name, comment and blank lines
# skipped, sample IDs kept whole. The selector refuses a task ID outside 1..n and a
# MANIFEST that is not the manifest CONFIG declares, and prints nothing when it refuses.
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set (submit this file as an array)}"
PYTHON="${PYTHON:-python3}"
SAMPLE=$("$PYTHON" -m magestic.pipelines.wgs.slurm.task_select \
    --config "$CONFIG" --manifest "$MANIFEST" --task "$SLURM_ARRAY_TASK_ID") || {
    echo "array task ${SLURM_ARRAY_TASK_ID} has no manifest sample; nothing was run" >&2
    exit 1
}

# §6 requeue safety lives in the stage, not here. A requeued task re-enters the stage
# CLI, whose cache reuses its outputs only when the recorded input, parameter and code
# identity still matches. An existence test on the output names would also "reuse" a
# stale or foreign file after the bridge, BAM or config changed, so there is none.

# The SAME console script the Snakemake path calls, so there is one entry point per
# stage and dual-path parity stays structural (§4).
magestic-wgs-bc1 --config "$CONFIG" --sample "$SAMPLE"
