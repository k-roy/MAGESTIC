# MAGESTIC

**M**ultiplexed **A**ccurate **G**enome **E**diting with **S**hort, **T**rackable, **I**ntegrated **C**ellular barcodes

Tools and pipelines for CRISPR-based genome editing in yeast.

## Contents

### Guide-Donor Oligo Design

Scripts for designing MAGESTIC guide-donor oligonucleotides:

- `generate_guide_donor_oligos_for_natural_variants_from_vcf.py` - Design oligos for natural variants from VCF
- `generate_amino_acid_variant_guide_donor_oligos_for_saturation_editing.py` - Design oligos for saturation mutagenesis
- `guide_donor_oligo_functions.py` - Shared functions for oligo design

### Data Processing

- `process_batch_editing_data.py` - Process batch editing experiments
- `process_overlapping_reads.sh` - Handle overlapping paired-end reads
- `collapse_identical_merged_reads_min_counts.py` - Collapse duplicate reads

### Genome Utilities

- `GTF_GFF_manipulation.py` - GFF/GTF file parsing and manipulation
- `bedgraph_computation.py` - Generate bedgraph coverage files

### Variant Effect Prediction (VEP/)

A unified pipeline for running multiple variant effect prediction models on yeast sequences. See [VEP/README.md](VEP/README.md) for details.

**Supported models:**
- **Evo2 7B/40B** - DNA sequence log-probability scores
- **Yorzoi** - RNA expression predictions from DNA
- **Shorkie** - RNA expression (8-fold ensemble)
- **ESM1-v** - Protein variant effect scores

**Quick start:**
```bash
cd VEP/reference_data/
gunzip -k *.gz
samtools faidx S288C_reference_sequence_R64-5-1_20240529.fsa
```

## Requirements

- Python 3.8+
- Biopython
- pandas, numpy
- samtools (for FASTA indexing)
- Snakemake 8+ (for VEP pipeline)

## Citation

If you use MAGESTIC, please cite:

> Roy KR, et al. (2018). Multiplexed precision genome editing with trackable genomic barcodes in yeast. *Nature Biotechnology* 36:512-520.

## License

MIT License - see individual files for details.
