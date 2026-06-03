# Reproduce: SpG QTL Oligo Pool (2024-04-11 Twist Order)

## Overview

This directory contains everything needed to reproduce the SpG QTL natural-variant
oligo pool using the MAGESTIC 3.0 package.

```
bash reproduce/20240404_SpG_QTL/run_design.sh    # from repo root
md5sum -c reproduce/20240404_SpG_QTL/expected_checksums.md5
```

Runtime: ~10 min on a 16-core interactive node.

---

## Known Discrepancy vs Historical Twist Order

The MAGESTIC 3.0 output (`SpG_QTL_2024_oligo_pool.tsv`, 54,907 oligos) **differs from
the historical `final_subpools/...with_10_bp_linked...` file** (40,426 oligos) that was
submitted to Twist in April 2024.

### Why they differ

| Factor | Legacy design (Twist order) | MAGESTIC 3.0 |
|--------|----------------------------|--------------|
| Variant linking | 10 bp max for all variants | 3 bp (SNVs), 10 bp (INDELs) |
| FDR filter | localFDR < 0.10 | localFDR < 0.10 ✓ |
| Shared parent filter | Applied (removes cross-shared variants) | **Not implemented** |
| Background variants | Included (479 oligos) | Excluded |
| Reverse primer | Subpool-specific (rotated) | `--subpool-idx 10` |

The **shared parent filter** is the primary driver of the discrepancy.  The legacy
`remove_shared_variants_and_filter.py` removed variants present in a parent shared by
two adjacent crosses of a round-robin cross design.  This cross-specific logic has not
been ported to MAGESTIC 3.0 (planned for v3.1).

### Overlap analysis

Stripping the 3′ reverse primer (which differs by subpool assignment) and comparing
guide+donor sequences:

- **25,751 oligos** share identical guide+donor sequences with the historical design (63.7%)
- **14,676** in historical only (≈ shared-parent-filtered variants + background variants)
- **29,156** in MAGESTIC 3.0 only (≈ variants removed by shared parent filter in legacy)

The MAGESTIC 3.0 design is **larger** and **more complete** in terms of unique variant
coverage; the shared parent filter was a conservative exclusion step whose biological
rationale applies mainly to within-cross epistasis analysis.

### Impact on screen analysis

Oligos present in the MAGESTIC 3.0 design but absent from the historical Twist order
**cannot be assayed** with the 2024 library — they were never synthesized.  Use the
harmonized design files in `common/oligo_designs/` for comparisons against screen data
from the 2024 pool.

---

## Checksums

`expected_checksums.md5` reflects the MAGESTIC 3.0 pipeline output (not the historical
Twist order).  The checksum will change if the pipeline parameters, FDR threshold, or
allowed-guide database are updated.

---

## Reference Files

All reference files are in
`projects/QTL/20240404_SpG_QTL_guide_donor_libraries/final_reference_files/`:

| File | Description |
|------|-------------|
| `Bloom_et_al_16_strains_QTL_genes_by_haplotype_block_coord_adjusted.tsv` | Per-sample haplotype VCF (pipeline input) |
| `Bloom_QTL_gene_level.tsv` | Gene-level QTL table with localFDR (for FDR filtering) |
| `MAGESTIC_background_strain.fasta` | Reference genome |
| `MAGESTIC_base_strain_allowed_SpG_Cas9_20mer_guides_NGNN.txt` | Pre-computed allowed guides |
| `MAGESTIC_base_strain_disallowed_SpG_Cas9_20mer_guides_NGNN.txt` | Pre-computed disallowed guides |
| `plus_strand_coverage_ypd.bedgraph` | TIF-seq coverage (Pelechano 2013) |
| `minus_strand_coverage_ypd.bedgraph` | TIF-seq coverage (Pelechano 2013) |
| `rev_subpool_priming_sequences.txt` | Reverse subpool primer sequences |
| `pS1380.fasta`, `pS1381.fasta`, `pF433.fasta` | FASTA files to exclude reserved guides |
