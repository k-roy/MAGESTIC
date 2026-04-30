# NRD1 Saturation Mutagenesis with MAGESTIC 3.0

Code and processed data for systematic NNS codon mutagenesis of the essential yeast transcription termination factor **NRD1** (YNL251C) using MAGESTIC 3.0 (Multiplexed Accurate Genome Editing with Short, Trackable, Integrated Cellular barcodes).

The library covers all 575 NRD1 residues with NNS codon substitutions, alanine-scan variants, synonymous-codon controls, and premature termination codon (PTC) controls. Designed mutations are tracked via unique genomic barcodes (bc1) integrated at the YBR209W locus and linked to their intended edits through amplicon sequencing.

---

## Repository Structure

```
NRD1/
├── README.md
│
├── # ── Library Design ──────────────────────────────────
├── generate_guide_donors_for_aa_substitutions_Cas9.py
├── guide_donor_functions.py
├── NRD1_alanine_scan_with_PTC_and_syn_controls.tsv       [23 KB]
├── NRD1_SpCas9_NGG_guide_donor_200bp_oligo_pool.tsv      [699 KB]
├── NRD1_SpCas9_NGG_guide_donor_200bp_info.tsv            [4 MB]
├── DNA_codon_to_AA.txt
├── DNA_codon_to_AA_suboptimal_codons_removed.txt
├── rev_subpool_priming_sequences.txt
│
├── # ── Sequence Maps ───────────────────────────────────
├── NRD1_MAGESTIC_3.0_sequence_maps/
│   ├── README.md
│   ├── YLS30.gb
│   ├── yT177.gb
│   ├── yT138.gb
│   ├── yT170.gb
│   ├── yKR1026_yT196.gb
│   ├── pS1381.gb
│   ├── pS1505.gb
│   ├── pS1439.gb
│   ├── pS1460.gb
│   ├── pS1461.gb
│   ├── pF833_cV437.gb
│   ├── pF834_pF835.gb
│   └── yT177_pF835_post_integration.gb
│
├── # ── Barcode Linking (bc1-donor-bc0) ─────────────────
├── bc1_donor_bc0_linking/
│   ├── step_a_trim_bc1_donor_bc0.sh
│   ├── step_b_process_bc1_donor_bc0_trimmed_R1_R2_fastq.py
│   ├── step_b_process_bc1_donor_bc0_trimmed_R1_R2_fastq.sh
│   ├── step_c_combine_all_bc1_donor_bc0_across_samples.py
│   ├── step_d_trim_merge_collapse_donor_bc0.sh
│   ├── step_e_process_donor_bc0_collapsed_fastq.py
│   ├── step_e_process_donor_bc0_collapsed_fastq.sh
│   ├── step_f_combine_donor_bc0_data_across_samples.py
│   ├── step_g_map_donor_bc0_to_designed_oligos.py
│   ├── step_h_link_complete_donor_bc0_and_guide_donor_bc0_to_bc1_donor_bc0_fragment.py
│   └── step_i_plot_bc1_reference_table.py
│
├── # ── Screen Processing ───────────────────────────────
├── screen_bc1_counting/
│   └── consolidate_bc1_counts.py
│
├── # ── Fitness Calling ─────────────────────────────────
├── fitness_slope_heatmap.py
├── fitness_slope_heatmap.R
│
├── # ── Figure Scripts ──────────────────────────────────
├── figure1_fitness_heatmap.py
├── figure_s1_screen_qc.py
│
└── # ── Processed Data (GitHub) ─────────────────────────
    data/
    ├── # Fitness results — input to figure scripts
    ├── slope_table_log2FC.csv                            [917 KB]
    ├── hits_summary.csv                                  [4 KB]
    ├── mean_log2FC_by_group.csv                          [2.2 MB]
    ├── depletion_log2FC.csv                              [446 KB]
    ├── barcode_support_flags.csv                         [644 KB]
    ├── hit_barcode_support_summary.csv                   [38 KB]
    └── # Lightweight counts — input to figure_s1_screen_qc.py replicate panels
        NRD1_bc1_counts_YPD_filtered.csv                 [9 MB]
        # (YPD + YPD+IAA conditions only; barcodes with mean YPD G0 ≥ 5 CPM; 14,967 barcodes)
        # Threshold on YPD G0 only — lethal/fitness-defective barcodes that deplete at later
        # timepoints are retained because they have normal counts at G0.
```

### Large files (Zenodo, DOI: TBD)

These files exceed GitHub's size limits and are deposited on Zenodo:

| File | Size | Description |
|------|------|-------------|
| `NRD1_bc1_reference_table.tsv` | 66 MB | Maps each bc1 barcode to its designed NRD1 mutation |
| `consolidated_normalized_counts.csv` | 373 MB | Full normalized barcode counts (all 8 conditions × 6 timepoints × 3 reps) |
| `MAGESTIC_background_strain.fasta` | ~12 MB | Genome sequence for the MAGESTIC background strain |
| `MAGESTIC_background_strain_annotations.gff` | ~5.6 MB | Gene annotations (GFF3) |
| `MAGESTIC_base_strain_allowed_Cas9_20mer_guides_NGG.txt` | ~34 MB | Pre-computed genome-wide SpCas9 allowed guides |
| `MAGESTIC_base_strain_disallowed_Cas9_20mer_guides_NGG.txt` | ~6.6 MB | Pre-computed SpCas9 disallowed guides |

---

## Pipeline Overview

The analysis proceeds in four stages:

```
1. Library Design          2. Barcode Linking         3. Screen Counting        4. Fitness Calling
   (oligo pool)     →      (bc1 ↔ mutation)    →     (bc1 counts/sample)  →   (slopes & hits)
                            amplicon-seq                amplicon-seq             computational
```

### Stage 1: Library Design

Designs SpCas9 guide RNA / donor DNA oligos for every desired NRD1 amino acid substitution. Each 200-nt oligo encodes a guide, a donor with the programmed mutation, and priming sites for subpool amplification.

**Input**: `NRD1_alanine_scan_with_PTC_and_syn_controls.tsv` (desired substitutions)
**Output**: `NRD1_SpCas9_NGG_guide_donor_200bp_oligo_pool.tsv` (oligo sequences for synthesis)

See [Library Design Details](#library-design-details) below.

### Stage 2: Barcode Linking (bc1-donor-bc0)

After MAGESTIC editing, each yeast clone carries a unique genomic barcode (bc1) at YBR209W that must be linked to the designed mutation it received. This is done by amplicon sequencing of a ~510 bp fragment spanning bc1, the donor region, and the donor barcode (bc0).

**Amplicons sequenced**:
- bc1-donor-bc0 (510 bp; primers KR1967 + KR1590)
- donor-bc0 (359 bp; primers KR1882 + KR1883)

**Processing pipeline** (`bc1_donor_bc0_linking/`):

| Step | Script | Description |
|------|--------|-------------|
| a | `step_a_trim_bc1_donor_bc0.sh` | Trim adapter sequences from paired-end reads |
| b | `step_b_process_bc1_donor_bc0_trimmed_R1_R2_fastq.py` | Extract bc1, donor, and bc0 from trimmed reads |
| c | `step_c_combine_all_bc1_donor_bc0_across_samples.py` | Merge results across multiplexed samples |
| d | `step_d_trim_merge_collapse_donor_bc0.sh` | Trim, merge, and collapse donor-bc0 amplicon reads |
| e | `step_e_process_donor_bc0_collapsed_fastq.py` | Extract donor sequence and bc0 from collapsed reads |
| f | `step_f_combine_donor_bc0_data_across_samples.py` | Merge donor-bc0 data across samples |
| g | `step_g_map_donor_bc0_to_designed_oligos.py` | Align observed donor-bc0 fragments to the designed oligo library |
| h | `step_h_link_complete_donor_bc0_and_guide_donor_bc0_to_bc1_donor_bc0_fragment.py` | Join bc1 ↔ donor-bc0 ↔ designed mutation into the final reference table |
| i | `step_i_plot_bc1_reference_table.py` | QC plots of the reference table (barcodes per variant, etc.) |

**Output**: `data/NRD1_bc1_reference_table.tsv` — maps every bc1 barcode to its designed NRD1 mutation.

### Stage 3: Screen bc1 Counting

Yeast libraries are grown competitively under 8 conditions (YPD ± IAA, Actinomycin D, CoCl2, Sorbitol, Galactose, MMS, Rapamycin) for 0–20 generations in 2 biological replicates. At each timepoint, genomic DNA is extracted and the bc1 region (232 bp; primers KR1966 + KR1967) is amplified and sequenced.

bc1 reads are counted per sample and joined to the reference table to produce a barcode × sample count matrix.

**Output**: `consolidated_normalized_counts_YPD.csv` — normalized counts per bc1 × treatment × generation.

*(Scripts TBD — see `screen_bc1_counting/`)*

### Stage 4: Fitness Calling

Fitness effects are estimated as log2 fold-change slopes over generations, median-centered within each treatment, with hit calling based on differential depletion between YPD+IAA (Nrd1 depleted) and YPD-IAA (Nrd1 present).

**Scripts**: `fitness_slope_heatmap.py` (Python) or `fitness_slope_heatmap.R` (R, original)

```bash
python fitness_slope_heatmap.py -i consolidated_normalized_counts_YPD.csv -o results/
```

**Pipeline steps**:
1. Load count matrix; classify mutations (missense / synonymous / PTC)
2. Convert to long format; compute per-barcode log2 fold-change vs generation 0
3. Collapse across barcodes (weighted mean by G0 count)
4. Trim depleted series: if k=3 consecutive post-G0 timepoints fall below FC < 0.01, truncate
5. Fit linear slope (log2FC ~ generation) on trimmed series
6. Median-center slopes within each treatment
7. Call hits: variants with (a) YPD+IAA adjusted slope < 2nd percentile, (b) slope_adj_YPD / 3 < slope_adj_YPD_noIAA, and (c) confirmed depletion (final FC < 1e-4)
8. Compute barcode-level support metrics for each hit

**Parameters** (defaults match the paper):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_total_count` | 10 | Minimum average G0 count to retain a barcode |
| `eps` | 1e-8 | Pseudocount to avoid log(0) |
| `use_weighted` | True | Weight collapse by G0 count |
| `fc_trim_cutoff` | 0.01 | Fold-change cutoff for series trimming |
| `sustained_k` | 3 | Consecutive points below cutoff to trigger trim |
| `hit_tail_prob` | 0.02 | Percentile threshold for YPD slope hit calling |
| `diff_factor_rule` | 3 | Required fold-difference between +IAA and -IAA slopes |
| `depletion_fc_cut` | 1e-4 | Final fold-change cutoff for depletion confirmation |

**Outputs** (6 CSV tables + 2 PDF plots):

| File | Description |
|------|-------------|
| `slope_table_log2FC.csv` | Per-variant, per-treatment fitness slopes with hit flags |
| `hits_summary.csv` | Variants called as fitness hits |
| `mean_log2FC_by_group.csv` | Mean log2FC trajectories by variant × treatment × generation |
| `depletion_log2FC.csv` | Depletion status (G0 vs final timepoint) |
| `barcode_support_flags.csv` | Per-barcode hit support flags |
| `hit_barcode_support_summary.csv` | Aggregated barcode support per variant |
| `heatmap_log2FC_slopes.pdf` | Heatmap of median-centered slopes by position and treatment |
| `dispersion_log2FC_slopes.pdf` | Scatter of YPD+IAA vs YPD-IAA adjusted slopes |

---

## Strain Background

All experiments use strain **yKR1026** (= yT196), which integrates four chromosomal modifications:

| Component | Locus | Description | Genetic Nomenclature |
|-----------|-------|-------------|----------------------|
| yT177 | YBR209W | Barcode landing pad | `natNT2-LYS2-CAN1-lox66-SceI_site-LEU2-SceI_site-kanR::ybr209w∆` |
| yT138 (WTC846) | YORWΔ17 | aTc-inducible system | `pRNR2-tetR-TUP1-pTDH3_tetO-tetR::YORW∆17` |
| yT170 | YNRCΔ9 | Gal/aTc-inducible SaCas9 guide-donor removal | `pGalL_SaCas9-pRPR1_tetO_HHR_SaCas9_X1_sgRNA::YNRC∆9` |
| YLS30 | YPRCΔ15 | Recoded NRD1 degron system | `pADH1_OsTIR1_FLAG_AID_NRD1rec_pRPB8B::YPRC∆15` |

In yKR1026/yT196 the LEU2 cassette at the yT177 landing pad is replaced with the full editing machinery:
`nat1-LYS2-CAN1-lox66-SceI_site-pTEF1-SpCas9-hphMX-pADH1-Ec86_RT-pDUT1-MCP-EcKlenow_DNA_polymerase-LexA-FHA-KlURA3-SceI_site-kanR::ybr209w∆`

Addition of IAA (indole-3-acetic acid) triggers OsTIR1-mediated degradation of the AID*-tagged recoded NRD1, making fitness of each mutant allele (at the endogenous locus) measurable by competitive growth.

YNRCΔ9, YORWΔ17, and YPRCΔ15 correspond to sites 16, 18, and 20, respectively, in Flagfeldt et al. 2009.

GenBank maps for all components are in `NRD1_MAGESTIC_3.0_sequence_maps/`.

---

## Screen Conditions

| Condition | Description | Biological Replicates |
|-----------|-------------|----------------------|
| YPD | Rich medium + IAA (Nrd1 depleted) | yL435, yL436 |
| YPD_noIAA | Rich medium, no IAA (Nrd1 present) | yL435, yL436 |
| Actinomycin_D | Transcription elongation stress | yL435, yL436 |
| CoCl2 | Hypoxia mimetic | yL435, yL436 |
| Sorbitol | Osmotic stress | yL435, yL436 |
| Galactose | Alternative carbon source | yL435, yL436 |
| MMS | DNA damage (alkylation) | yL435, yL436 |
| Rapamycin | TOR pathway inhibition | yL435, yL436 |

Timepoints: 0, 4, 8, 12, 16, 20 generations.

---

## Data Availability

- **Raw sequencing data**: NCBI Sequence Read Archive, BioProject [PRJNA1458821](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1458821)
- **Processed data**: Zenodo (DOI TBD, to be deposited during revision)

---

## Library Design Details

*(Preserved from the original README — full documentation of the oligo design pipeline.)*

### Contents (Library Design)

| File | Description |
|------|-------------|
| `generate_guide_donors_for_aa_substitutions_Cas9.py` | Main script — designs guide/donor oligos |
| `guide_donor_functions.py` | Utility functions (sequence manipulation, ORF parsing, donor assembly) |
| `NRD1_alanine_scan_with_PTC_and_syn_controls.tsv` | Input: list of desired amino acid substitutions |
| `NRD1_SpCas9_NGG_guide_donor_200bp_oligo_pool.tsv` | Output: pre-computed oligo library (names + sequences) |
| `NRD1_SpCas9_NGG_guide_donor_200bp_info.tsv` | Output: full per-oligo design details |
| `DNA_codon_to_AA.txt` | Standard genetic code (codon → amino acid) |
| `DNA_codon_to_AA_suboptimal_codons_removed.txt` | Genetic code with rare yeast codons removed |
| `rev_subpool_priming_sequences.txt` | Reverse priming sequences for each oligo subpool |

Guides that would target any of the integrated strain or vector sequences are
excluded from the library to prevent re-cutting during editing. The six relevant
maps used for guide exclusion are:

| File | Description | Genetic Nomenclature |
|------|-------------|----------------------|
| `yT138.gb` | aTc-inducible system @ YORWΔ17 | `pRNR2-tetR-TUP1-pTDH3_tetO-tetR::YORW∆17` |
| `yT170.gb` | Gal/aTc-inducible SaCas9 guide-donor removal @ YNRCΔ9 | `pGalL_SaCas9-pRPR1_tetO_HHR_SaCas9_X1_sgRNA::YNRC∆9` |
| `yKR1026_yT196.gb` | Editing machinery at barcode locus (YBR209W) | `nat1-LYS2-CAN1-lox66-SceI_site-pTEF1-SpCas9-hphMX-...::ybr209w∆` |
| `YLS30.gb` | Recoded NRD1 degron system @ YPRCΔ15 | `pADH1_OsTIR1_FLAG_AID_NRD1rec_pRPB8B::YPRC∆15` |
| `pS1439.gb` | MAGESTIC 3.0 guide-donor plasmid backbone | `CEN/ARS-AmpR-FCY1-pSNR52-sgRNA-BspQI-AscI-retron_msr/msd-CYC1t-X_site` |
| `pF834_pF835.gb` | Barcoded insert (pF834) and final insert for transformation (pF835) | `CYC1t-pGalL-SceI-kanR-HIS3-GAL7t-N20_bc1-lox66-X_site-msr/msd` |

### Large Reference Files (download from Zenodo)

Four large reference files are required to run the design script. They are not
bundled in this repository due to file size and will be available on Zenodo:

> **DOI: TBD**

| File | Description |
|------|-------------|
| `MAGESTIC_background_strain.fasta` | Genome sequence of the MAGESTIC background yeast strain (~12 MB) |
| `MAGESTIC_background_strain_annotations.gff` | Gene annotations in GFF3 format (~5.6 MB) |
| `MAGESTIC_base_strain_allowed_Cas9_20mer_guides_NGG.txt` | Pre-computed genome-wide SpCas9 NGG guide database — allowed guides (~34 MB) |
| `MAGESTIC_base_strain_disallowed_Cas9_20mer_guides_NGG.txt` | Pre-computed genome-wide SpCas9 NGG guide database — disallowed guides (~6.6 MB) |

The allowed/disallowed guide databases were generated by:
1. Enumerating all 20-nt NGG SpCas9 guide sequences in the genome
   (`process_NGG_guides_genomewide.py` in the repo root)
2. Filtering for guides with acceptable off-target profiles using Cas-OFFinder
   (`assign_allowable_guides.py` in the repo root)

This procedure is described in detail in [Roy et al. (2018) Nature Biotechnology](https://www.nature.com/articles/nbt.4137).

### Dependencies (Library Design)

Python 3.9+ (no third-party packages required beyond the standard library).

### Usage (Library Design)

After downloading the four large reference files into a directory (e.g., `./ref/`):

```bash
python generate_guide_donors_for_aa_substitutions_Cas9.py --ref-dir ./ref/
```

Output files are written to the script directory by default. Use `--output-dir`
to specify a different location:

```bash
python generate_guide_donors_for_aa_substitutions_Cas9.py --ref-dir ./ref/ --output-dir ./results/
```

### Input File Format

`NRD1_alanine_scan_with_PTC_and_syn_controls.tsv` is a tab-separated file
with a one-line header. Each data row specifies one amino acid substitution.
This library uses the 3-column format:

```
systematic_ORF_name  common_gene_name  aa_substitution
YNL251C              NRD1              Q2A
YNL251C              NRD1              M1*
```

`aa_substitution` uses one-letter amino acid codes: `<ref><position><alt>`
(e.g., `Q2A` = glutamine at position 2 → alanine; `M1*` = methionine → stop codon).

The script also supports 5-, 7-, 9-, and 11-column formats that allow explicit
specification of reference/alternate codons and additional metadata; see the
function-level docstring in `guide_donor_functions.py` for details.

### Output Files (Library Design)

**`*_oligo_pool.tsv`** — Two-column TSV (`oligo_name`, `oligo_seq`) ready for oligonucleotide synthesis.

**`*_info.tsv`** — Full per-oligo design details including guide sequence, PAM coordinates, donor sequence, homology arm lengths, synonymous changes introduced, and whether any synonymous codon substitutions or donor-arm shifts were needed to eliminate BspQI restriction sites at oligo junctions.

**`*_untargetable.tsv`** — Mutations for which no acceptable guide/donor design could be found, with the reason each mutation could not be targeted.

### Oligo Architecture

Each 200-nt oligo has the structure:

```
[FWD_PRIMING_SITE (20 nt)] [guide (20 nt)] [INTERNAL_CLONING_SITE (11 nt)] [donor (129 nt)] [REV_PRIMING_SITE (20 nt)]
```

The internal cloning site (`gtttgaagagc`) contains one BspQI (SapI) recognition
sequence used for ligation of the step-2 cloning insert, which consists of the
guide RNA scaffold, kanR cassette, and I-SceI sites that mediate linearization
of the guide-donor vector for MAGESTIC 3.0.

Synonymous codon changes are introduced in the guide seed / PAM-proximal region
of the donor to prevent Cas9 re-cutting after successful HDR, and to enable
editing outside of SpCas9 guide regions for full saturation mutagenesis.

---

## Dependencies (Full Pipeline)

### Library Design
- Python 3.9+ (standard library only)

### Barcode Linking
- Python 3.9+
- cutadapt
- bbmerge (BBTools)
- pandas, numpy, matplotlib

### Fitness Calling
- Python 3.9+
- pandas, numpy, scipy, matplotlib, seaborn

---

## Citation

If you use this code or data, please cite:

> Aiello U, Roy KR, Steinmetz LM. (2025). Mutational scanning by multiplexed genome editing of the essential transcription termination factor Nrd1. *bioRxiv* 2025.09.09.675026. https://doi.org/10.1101/2025.09.09.675026

For the MAGESTIC framework:

> Roy KR et al. (2018). Multiplexed precision genome editing with trackable genomic barcodes in yeast. *Nature Biotechnology* 36, 512–520. https://doi.org/10.1038/nbt.4137

> Roy KR, Smith JD, Li S, Vonesch SC, Nguyen M, Burnett WT, Orsley KM, Lee CS, Haber JE, St Onge RP, Steinmetz LM. (2024). Dissecting quantitative trait nucleotides by saturation genome editing. *bioRxiv* 2024.02.02.577784. https://doi.org/10.1101/2024.02.02.577784
