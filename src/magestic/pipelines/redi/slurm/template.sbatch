#!/bin/bash
# REDI per-sample step_a SLURM template.
# Used by the Snakemake workflow OR can be submitted standalone via the
# `sherlock-sbatch` skill.
#
# Encodes the four documented Sherlock traps from CLAUDE.md:
#   - set -u AFTER source ~/.bashrc (COLORTERM is unbound in 256term.sh)
#   - AVX-512 SIGILL constraint for numpy/MKL-heavy Python
#   - LD_LIBRARY_PATH propagation when SSH-driving sbatch
#   - chunked into 20–60 min array tasks for owners-queue preemption safety
#
#SBATCH --job-name=redi_step_a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=mylab,owners
#SBATCH --account=mylab
#SBATCH --time=01:00:00
#SBATCH --constraint="CPU_GEN:SIE|CPU_GEN:GEN|CPU_GEN:BGM|CPU_GEN:SKX"
#SBATCH --output=logs/redi_step_a_%A_%a.out
#SBATCH --error=logs/redi_step_a_%A_%a.err

# Source ~/.bashrc BEFORE set -u (COLORTERM trap)
source ~/.bashrc
set -euo pipefail

# Activate conda env
conda activate magestic

# CONFIG and SAMPLE_PREFIX expected as env vars from submitter
: "${CONFIG:?CONFIG (path to YAML) must be set}"
: "${SAMPLE_PREFIX:?SAMPLE_PREFIX must be set}"

magestic-redi-trim --config "$CONFIG" --sample "$SAMPLE_PREFIX"
