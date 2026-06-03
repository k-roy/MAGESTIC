#!/usr/bin/env bash
# MAGESTIC 3.0 — Reproduce SpCas9 + LbCas12a QTL oligo pool (2021-04-22 Twist order)
#
# Runtime: ~60 min on a 16-core interactive node
# Assumes you are in the MAGESTIC repository root and have MAGESTIC 3.0 installed.
#
# Usage:
#   bash reproduce/20210422_SpCas9_LbCas12a_QTL/run_design.sh
#
# After completion, verify against the canonical Twist order:
#   md5sum -c reproduce/20210422_SpCas9_LbCas12a_QTL/expected_checksums.md5
#
# NOTES on design complexity:
#   The 2021 oligo pool combined three nuclease sub-designs (SpCas9 + LbCas12a +
#   impLbCas12a) and saturation mutagenesis sub-pools for PDR5 and DHY214.
#   This script reproduces the QTL variant-editing component (the majority of
#   the pool). The saturation sub-pools (AA substitution controls) were generated
#   separately; see the legacy final_scripts/ for provenance.

set -euo pipefail

PYTHON="/path/to/anaconda3/bin/python"
BASE_DIR="/path/to"
PROJ="$BASE_DIR/projects/QTL/20210422_SpCas9_LbCas12a_QTL_guide_donor_libraries"
REF="$PROJ/final_reference_files"
COMMON="$BASE_DIR/common"
OUT_DIR="reproduce/20210422_SpCas9_LbCas12a_QTL/output"

echo "=== MAGESTIC 3.0 — SpCas9 + LbCas12a QTL Oligo Design ==="
echo "Start: $(date)"
mkdir -p "$OUT_DIR"

# Step 1 — Process haplotype VCF
echo ""
echo "--- Step 1: VCF processing ---"
magestic-vcf-process \
    --haplotype-vcf "$REF/Bloom_et_al_16_strains_QTL_genes_by_haplotype_block_coord_adjusted.tsv" \
    --genome "$REF/MAGESTIC_background_strain.fasta" \
    --exclude-background \
    --exclude-genes HO \
    --max-bp-between-snvs 5 \
    --max-bp-between-indels 5 \
    --output "$OUT_DIR/processed_variants.tsv"

echo "Processed variants written to $OUT_DIR/processed_variants.tsv"

# Step 2 — Check for pre-computed guide files
echo ""
echo "--- Step 2: Verifying pre-computed guide files ---"

CAS9_ALLOWED="$REF/MAGESTIC_base_strain_allowed_Cas9_20mer_guides_NGG.txt"
CAS9_DISALLOWED="$REF/MAGESTIC_base_strain_disallowed_Cas9_20mer_guides_NGG.txt"
CAS12A_ALLOWED="$REF/MAGESTIC_base_strain_allowed_Cas12a_TTTV_20mer_guides.txt"
CAS12A_DISALLOWED="$REF/MAGESTIC_base_strain_disallowed_Cas12a_TTTV_20mer_guides.txt"

for f in "$CAS9_ALLOWED" "$CAS9_DISALLOWED" "$CAS12A_ALLOWED" "$CAS12A_DISALLOWED"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: guide file not found: $f"
        exit 1
    fi
done
echo "All pre-computed guide files found."

# Step 3a — SpCas9 NGG oligo design (subpool 0)
echo ""
echo "--- Step 3a: SpCas9 NGG oligo design (subpool 0) ---"
magestic-design-vcf \
    --vcf "$OUT_DIR/processed_variants.tsv" \
    --genome "$REF/MAGESTIC_background_strain.fasta" \
    --allowed-guides "$CAS9_ALLOWED" \
    --disallowed-guides "$CAS9_DISALLOWED" \
    --plus-bedgraph "$COMMON/reference_genomes/Pelechano_2013_PMID_23615609/plus_strand_coverage_ypd.bedgraph" \
    --minus-bedgraph "$COMMON/reference_genomes/Pelechano_2013_PMID_23615609/minus_strand_coverage_ypd.bedgraph" \
    --rev-primers "$REF/rev_subpool_priming_sequences.txt" \
    --nuclease SpCas9 \
    --fasta-exclusion \
        "$REF/pS1380.fasta" \
        "$REF/pS1381.fasta" \
        "$REF/pF433.fasta" \
    --subpool-idx 0 \
    --oligo-length 200 \
    --minimum-homology 30 \
    --max-guides-per-donor 2 \
    --max-dist-to-pam 17 \
    --output-dir "$OUT_DIR" \
    --output-prefix "SpCas9_QTL_2021"

echo "SpCas9 oligos written to $OUT_DIR/SpCas9_QTL_2021_oligo_pool.tsv"

# Step 3b — LbCas12a TTTV oligo design (subpool 1)
echo ""
echo "--- Step 3b: LbCas12a TTTV oligo design (subpool 1) ---"
magestic-design-vcf \
    --vcf "$OUT_DIR/processed_variants.tsv" \
    --genome "$REF/MAGESTIC_background_strain.fasta" \
    --allowed-guides "$CAS12A_ALLOWED" \
    --disallowed-guides "$CAS12A_DISALLOWED" \
    --plus-bedgraph "$COMMON/reference_genomes/Pelechano_2013_PMID_23615609/plus_strand_coverage_ypd.bedgraph" \
    --minus-bedgraph "$COMMON/reference_genomes/Pelechano_2013_PMID_23615609/minus_strand_coverage_ypd.bedgraph" \
    --rev-primers "$REF/rev_subpool_priming_sequences.txt" \
    --nuclease LbCas12a \
    --fasta-exclusion \
        "$REF/pS1380.fasta" \
        "$REF/pS1381.fasta" \
        "$REF/pF433.fasta" \
    --subpool-idx 1 \
    --oligo-length 200 \
    --minimum-homology 30 \
    --max-guides-per-donor 2 \
    --max-dist-to-pam 20 \
    --output-dir "$OUT_DIR" \
    --output-prefix "LbCas12a_QTL_2021"

echo "LbCas12a oligos written to $OUT_DIR/LbCas12a_QTL_2021_oligo_pool.tsv"

# NOTE: impLbCas12a sub-pool (subpool 2) used an extended PAM set and was
# generated with a separate script invocation. See:
#   final_scripts/generate_guide_donor_sequences_from_vcf.py (NUCLEASE='impLbCas12a' branch)
# impLbCas12a is supported in MAGESTIC 3.0 via --nuclease impLbCas12a but is
# not invoked here because the canonical file also includes saturation sub-pools
# that are not part of the automated QTL variant design.

echo ""
echo "=== QTL variant design complete: $(date) ==="
echo ""
echo "NOTE: The canonical 2021 pool (48,012 oligos) also includes saturation"
echo "mutagenesis sub-pools for PDR5 and DHY214 that were generated separately."
echo "See reproduce/20210422_SpCas9_LbCas12a_QTL/config.yaml for full details."
echo ""
echo "Verify QTL variant components against reference counts:"
echo "  wc -l $OUT_DIR/SpCas9_QTL_2021_oligo_pool.tsv"
echo "  wc -l $OUT_DIR/LbCas12a_QTL_2021_oligo_pool.tsv"
echo ""
echo "Verify combined pool checksum:"
echo "  md5sum -c reproduce/20210422_SpCas9_LbCas12a_QTL/expected_checksums.md5"
