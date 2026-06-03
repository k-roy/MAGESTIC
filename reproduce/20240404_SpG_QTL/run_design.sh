#!/usr/bin/env bash
# MAGESTIC 3.0 — Reproduce SpG QTL oligo pool (2024-04-11 Twist order)
#
# Runtime: ~30 min on a 16-core interactive node (VCF processing + oligo design)
# Assumes you are in the MAGESTIC repository root and have MAGESTIC 3.0 installed.
#
# Usage:
#   bash reproduce/20240404_SpG_QTL/run_design.sh
#
# After completion, verify against the canonical Twist order:
#   md5sum -c reproduce/20240404_SpG_QTL/expected_checksums.md5

set -euo pipefail

PYTHON="/path/to/anaconda3/bin/python"
BASE_DIR="/path/to"
PROJ="$BASE_DIR/projects/QTL/20240404_SpG_QTL_guide_donor_libraries"
REF="$PROJ/final_reference_files"
COMMON="$BASE_DIR/common"
OUT_DIR="reproduce/20240404_SpG_QTL/output"

echo "=== MAGESTIC 3.0 — SpG QTL Oligo Design ==="
echo "Start: $(date)"
mkdir -p "$OUT_DIR"

# Step 1 — Process haplotype VCF
echo ""
echo "--- Step 1: VCF processing ---"
magestic-vcf-process \
    --haplotype-vcf "$REF/Bloom_et_al_16_strains_QTL_genes_by_haplotype_block_coord_adjusted.tsv" \
    --genome "$REF/MAGESTIC_background_strain.fasta" \
    --gene-list "$REF/Bloom_QTL_gene_level.tsv" \
    --max-fdr 0.10 \
    --exclude-background \
    --exclude-genes HO \
    --max-bp-between-snvs 3 \
    --max-bp-between-indels 10 \
    --output "$OUT_DIR/processed_variants.tsv"

echo "Processed variants written to $OUT_DIR/processed_variants.tsv"

# Step 2 — Guide selection (uses pre-computed off-target database)
# NOTE: re-run magestic-guide-selection only if the off-target database was regenerated.
echo ""
echo "--- Step 2: Using pre-computed allowed/disallowed guide files ---"
ALLOWED="$REF/MAGESTIC_base_strain_allowed_SpG_Cas9_20mer_guides_NGNN.txt"
DISALLOWED="$REF/MAGESTIC_base_strain_disallowed_SpG_Cas9_20mer_guides_NGNN.txt"

if [ ! -f "$ALLOWED" ]; then
    echo "ERROR: allowed guides file not found: $ALLOWED"
    echo "Run: magestic-pam-search --genome $REF/MAGESTIC_background_strain.fasta --nuclease SpG --output <off_target_file>"
    echo "Then: magestic-guide-selection --off-target-file <off_target_file> --allowed-out $ALLOWED --disallowed-out $DISALLOWED"
    exit 1
fi

# Step 3 — Oligo design
echo ""
echo "--- Step 3: Oligo design ---"
magestic-design-vcf \
    --vcf "$OUT_DIR/processed_variants.tsv" \
    --genome "$REF/MAGESTIC_background_strain.fasta" \
    --allowed-guides "$ALLOWED" \
    --disallowed-guides "$DISALLOWED" \
    --plus-bedgraph "$COMMON/reference_genomes/Pelechano_2013_PMID_23615609/plus_strand_coverage_ypd.bedgraph" \
    --minus-bedgraph "$COMMON/reference_genomes/Pelechano_2013_PMID_23615609/minus_strand_coverage_ypd.bedgraph" \
    --rev-primers "$REF/rev_subpool_priming_sequences.txt" \
    --nuclease SpG \
    --fasta-exclusion \
        "$REF/pS1380.fasta" \
        "$REF/pS1381.fasta" \
        "$REF/pF433.fasta" \
    --subpool-idx 10 \
    --oligo-length 200 \
    --minimum-homology 30 \
    --max-guides-per-donor 4 \
    --max-dist-to-pam 20 \
    --output-dir "$OUT_DIR" \
    --output-prefix "SpG_QTL_2024"

echo ""
echo "=== Design complete: $(date) ==="
echo "Outputs:"
echo "  $OUT_DIR/SpG_QTL_2024_oligo_pool.tsv"
echo "  $OUT_DIR/SpG_QTL_2024_info.tsv"
echo "  $OUT_DIR/SpG_QTL_2024_untargetable.tsv"
echo ""
echo "Verify against canonical output:"
echo "  md5sum -c reproduce/20240404_SpG_QTL/expected_checksums.md5"
