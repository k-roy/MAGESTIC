#!/bin/bash
#
# download_references.sh
#
# Download S288C reference genome and annotations from SGD.
# These files are required for the VEP pipeline.
#
# Usage:
#   cd reference_data/
#   bash download_references.sh
#
# Author: Kevin R. Roy
# Date: 2026-02-26
#
# =============================================================================

set -euo pipefail

# Configuration
SGD_BASE_URL="https://downloads.yeastgenome.org/sequence/S288C_reference"
RELEASE="R64-5-1_20240529"
OUTPUT_DIR="${1:-.}"  # Default to current directory

echo "========================================"
echo "Downloading S288C Reference Files"
echo "========================================"
echo "Release: ${RELEASE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Create output directory
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# -----------------------------------------------------------------------------
# 1. Download reference genome FASTA
# -----------------------------------------------------------------------------
echo "[1/5] Downloading reference genome FASTA..."
FASTA_URL="${SGD_BASE_URL}/genome_releases/S288C_reference_genome_${RELEASE}.tgz"
FASTA_FILE="S288C_reference_sequence_${RELEASE}.fsa"

if [[ -f "${FASTA_FILE}" ]]; then
    echo "  Already exists: ${FASTA_FILE}"
else
    wget -q --show-progress "${FASTA_URL}" -O "S288C_reference_genome.tgz"
    tar -xzf "S288C_reference_genome.tgz"

    # Move files from extracted directory
    if [[ -d "S288C_reference_genome_${RELEASE}" ]]; then
        mv "S288C_reference_genome_${RELEASE}"/* .
        rmdir "S288C_reference_genome_${RELEASE}"
    fi
    rm -f "S288C_reference_genome.tgz"

    # Decompress FASTA if gzipped
    if [[ -f "${FASTA_FILE}.gz" ]] && [[ ! -f "${FASTA_FILE}" ]]; then
        gunzip -k "${FASTA_FILE}.gz"
    fi
    echo "  Downloaded: ${FASTA_FILE}"
fi

# -----------------------------------------------------------------------------
# 2. Download GFF annotations
# -----------------------------------------------------------------------------
echo "[2/5] Downloading GFF annotations..."
GFF_FILE="saccharomyces_cerevisiae_${RELEASE}.gff"

if [[ -f "${GFF_FILE}" ]]; then
    echo "  Already exists: ${GFF_FILE}"
else
    if [[ -f "${GFF_FILE}.gz" ]]; then
        gunzip -k "${GFF_FILE}.gz"
    else
        # GFF should be in the tarball, but download separately if needed
        GFF_URL="${SGD_BASE_URL}/genome_releases/S288C_reference_genome_${RELEASE}.tgz"
        echo "  GFF should be in the genome tarball"
    fi
    echo "  Available: ${GFF_FILE}"
fi

# -----------------------------------------------------------------------------
# 3. Create FASTA index
# -----------------------------------------------------------------------------
echo "[3/5] Creating FASTA index..."
FAI_FILE="${FASTA_FILE}.fai"

if [[ -f "${FAI_FILE}" ]]; then
    echo "  Already exists: ${FAI_FILE}"
else
    if command -v samtools &> /dev/null; then
        samtools faidx "${FASTA_FILE}"
        echo "  Created: ${FAI_FILE}"
    else
        echo "  WARNING: samtools not found, skipping index creation"
        echo "  Run 'samtools faidx ${FASTA_FILE}' manually"
    fi
fi

# -----------------------------------------------------------------------------
# 4. Download ORF coding sequences (optional, for protein analyses)
# -----------------------------------------------------------------------------
echo "[4/5] Downloading ORF coding sequences..."
ORF_FILE="orf_coding_all_${RELEASE}.fasta"

if [[ -f "${ORF_FILE}" ]] || [[ -f "${ORF_FILE}.gz" ]]; then
    echo "  Already exists: ${ORF_FILE}"
else
    # Should be in tarball
    if [[ -f "${ORF_FILE}.gz" ]]; then
        echo "  Available: ${ORF_FILE}.gz"
    fi
fi

# -----------------------------------------------------------------------------
# 5. Download protein sequences (optional, for ESM1-v)
# -----------------------------------------------------------------------------
echo "[5/5] Downloading protein sequences..."
PROTEIN_FILE="orf_trans_all_${RELEASE}.fasta"

if [[ -f "${PROTEIN_FILE}" ]] || [[ -f "${PROTEIN_FILE}.gz" ]]; then
    echo "  Already exists: ${PROTEIN_FILE}"
else
    if [[ -f "${PROTEIN_FILE}.gz" ]]; then
        echo "  Available: ${PROTEIN_FILE}.gz"
    fi
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Download Complete"
echo "========================================"
echo ""
echo "Required files:"
ls -lh "${FASTA_FILE}" 2>/dev/null || echo "  MISSING: ${FASTA_FILE}"
ls -lh "${GFF_FILE}" 2>/dev/null || echo "  MISSING: ${GFF_FILE}"
echo ""
echo "Optional files:"
ls -lh "${ORF_FILE}"* 2>/dev/null || echo "  Not downloaded: ${ORF_FILE}"
ls -lh "${PROTEIN_FILE}"* 2>/dev/null || echo "  Not downloaded: ${PROTEIN_FILE}"
echo ""
echo "To use these files with the VEP pipeline:"
echo "  1. Update your vep_config.yaml to point to this directory"
echo "  2. Or copy files to your project's reference_data/ folder"
echo ""
