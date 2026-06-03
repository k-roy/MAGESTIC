# MAGESTIC Reproducibility Wrappers

This directory contains configuration files and scripts to reproduce the two published
oligo pool designs generated with MAGESTIC.

## Requirements

- MAGESTIC ≥ 3.0.0 installed (`pip install magestic[design]`)
- Reference files from the project directories (see paths in each `config.yaml`)
- Run all commands from the **MAGESTIC repository root**

## Designs

### 2024 — SpG QTL oligo pool (180 k oligos)

**Directory:** `20240404_SpG_QTL/`

Published array order: `common/oligo_designs/20240411_Twist_200mer_oligo_array_order.tsv`

```bash
bash reproduce/20240404_SpG_QTL/run_design.sh
md5sum -c reproduce/20240404_SpG_QTL/expected_checksums.md5
```

**Design summary:**
- Nuclease: SpG (SpCas9 variant with relaxed NGNN PAM)
- Variants: Bloom et al. 16-strain QTL SNPs/INDELs
- Oligo length: 200 bp
- PAM preference: NGNG (subpool 10) > NGNH (subpool 11)
- ~180,000 oligos targeting ~1,000 QTL variants in 170 genes

### 2021 — SpCas9 + LbCas12a QTL oligo pool (48 k oligos)

**Directory:** `20210422_SpCas9_LbCas12a_QTL/`

Published array order: `projects/QTL/20210422_SpCas9_LbCas12a_QTL_guide_donor_libraries/final_oligo_pools/20210422_Twist_200mer_oligo_array_order.tsv`

```bash
bash reproduce/20210422_SpCas9_LbCas12a_QTL/run_design.sh
```

**Design summary:**
- Nucleases: SpCas9 (NGG, subpool 0) + LbCas12a (TTTV, subpool 1) + impLbCas12a (extended PAMs, subpool 2)
- Oligo length: 200 bp
- ~28,000 QTL variant oligos + ~20,000 saturation mutagenesis controls (PDR5, DHY214)

> **Note:** The canonical 2021 pool also contains saturation mutagenesis sub-pools
> for PDR5 and DHY214 that were generated with separate scripts. The automated
> `run_design.sh` reproduces only the QTL variant component. See `config.yaml`
> for full provenance documentation.

## Verification

After running a design, compare row counts and checksums:

```bash
# 2024 SpG design
wc -l reproduce/20240404_SpG_QTL/output/SpG_QTL_2024_oligo_pool.tsv
# Expected: 180001 (180000 oligos + header)
md5sum -c reproduce/20240404_SpG_QTL/expected_checksums.md5

# 2021 SpCas9+LbCas12a combined pool
wc -l projects/QTL/20210422_SpCas9_LbCas12a_QTL_guide_donor_libraries/final_oligo_pools/20210422_Twist_200mer_oligo_array_order.tsv
# Expected: 48013 (48012 oligos + header)
md5sum -c reproduce/20210422_SpCas9_LbCas12a_QTL/expected_checksums.md5
```

## Legacy Scripts

The original design scripts are preserved for provenance in each project's
`final_scripts/` directory:

- `projects/QTL/20240404_SpG_QTL_guide_donor_libraries/final_scripts/`
- `projects/QTL/20210422_SpCas9_LbCas12a_QTL_guide_donor_libraries/final_scripts/`

These scripts have been superseded by the MAGESTIC 3.0 package but are kept
for reference. See each project's `REPRODUCE.md` for the mapping between legacy
script parameters and MAGESTIC 3.0 CLI arguments.
