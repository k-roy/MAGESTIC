#!/usr/bin/env bash
# Reproduce: Guide-Donor-BC0 matching + BC0 purity for the 2021 SpCas9/LbCas12a QTL library
#
# Usage (from repo root):
#   bash reproduce/20210422_SpCas9_LbCas12a_QTL_gdb/run_pipeline.sh
#
# Requires:
#   magestic-gdb-pipeline  (pip install -e .[dev] from repo root)
#
# Runtime: ~15-20 min on an interactive node (3M sequences).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REPRO_DIR="${REPO_ROOT}/reproduce/20210422_SpCas9_LbCas12a_QTL_gdb"
OUTPUT_DIR="${REPRO_DIR}/output"

BASE_DIR="/path/to"
PROJECT_DIR="${BASE_DIR}/projects/QTL/20210422_SpCas9_LbCas12a_QTL_guide_donor_libraries"

# Full Step 02 output (3M rows) produced by the legacy pipeline
COUNTS_SRC="${PROJECT_DIR}/processed_data/guide_donor_bc0/guide_donor_bc0/V455_guide_donor_bc0_counts.tsv"

OLIGO_DESIGN="${BASE_DIR}/common/oligo_designs/20210422_Twist_200mer.tsv"

echo "=== Reproduce: Guide-Donor-BC0 pipeline (2021 SpCas9/LbCas12a QTL) ==="
echo "Repo:          ${REPO_ROOT}"
echo "Source counts: ${COUNTS_SRC}"
echo "Oligo design:  ${OLIGO_DESIGN}"
echo "Output:        ${OUTPUT_DIR}"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Validate inputs
# ─────────────────────────────────────────────────────────────────────────────
if [[ ! -f "${COUNTS_SRC}" ]]; then
    echo "ERROR: Source counts file not found: ${COUNTS_SRC}"
    exit 1
fi

if [[ ! -f "${OLIGO_DESIGN}" ]]; then
    echo "ERROR: Oligo design file not found: ${OLIGO_DESIGN}"
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Stage counts file under the expected pipeline directory structure
#
# magestic-gdb-pipeline writes output to:
#   {project_dir}/processed_data/guide_donor_bc0/guide_donor_bc0/
#
# It searches that same directory for *_guide_donor_bc0_counts.tsv input files.
# ─────────────────────────────────────────────────────────────────────────────
COUNTS_DIR="${OUTPUT_DIR}/processed_data/guide_donor_bc0/guide_donor_bc0"
mkdir -p "${COUNTS_DIR}"
COUNTS_DEST="${COUNTS_DIR}/V455_guide_donor_bc0_counts.tsv"

echo "Staging counts file..."
cp "${COUNTS_SRC}" "${COUNTS_DEST}"
echo "  $(wc -l < "${COUNTS_DEST}") rows → ${COUNTS_DEST}"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Run Steps 03 (oligo matching) and 04 (BC0 purity)
# ─────────────────────────────────────────────────────────────────────────────
echo "Running Steps 03 (match) and 04 (purity)..."

magestic-gdb-pipeline \
    --project-dir "${OUTPUT_DIR}" \
    --oligo-design-file "${OLIGO_DESIGN}" \
    --steps match,purity

echo ""
echo "=== Pipeline complete ==="
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Verify checksums
# ─────────────────────────────────────────────────────────────────────────────
echo "Output files:"
ls -lh "${OUTPUT_DIR}/processed_data/guide_donor_bc0/guide_donor_bc0/" 2>/dev/null || true
ls -lh "${OUTPUT_DIR}/processed_data/guide_donor_bc0/guide_donor_bc0_purity_filtered/" 2>/dev/null || true
echo ""

echo "Verifying checksums..."
cd "${REPRO_DIR}"
md5sum -c expected_checksums.md5
