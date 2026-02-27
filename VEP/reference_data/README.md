# VEP Reference Data

This directory contains reference genome files needed for the VEP pipeline.

## Quick Start

### Option 1: Use included files (recommended)

The compressed reference files are included in this repo:

```bash
cd reference_data/

# Decompress the files you need
gunzip -k S288C_reference_sequence_R64-5-1_20240529.fsa.gz
gunzip -k saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz

# Create FASTA index (requires samtools)
samtools faidx S288C_reference_sequence_R64-5-1_20240529.fsa
```

### Option 2: Download fresh from SGD

```bash
bash download_references.sh
```

## Included Files

| File | Size | Description |
|------|------|-------------|
| `S288C_reference_sequence_R64-5-1_20240529.fsa.gz` | 3.7 MB | S288C genome FASTA (gzipped) |
| `saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz` | 5.0 MB | Gene annotations GFF3 (gzipped) |
| `orf_coding_all_R64-5-1_20240529.fasta.gz` | 3.7 MB | All ORF coding sequences |
| `orf_trans_all_R64-5-1_20240529.fasta.gz` | 2.6 MB | All protein sequences |
| `download_references.sh` | 4 KB | Download script |

## File Details

### S288C Reference Genome (Required)

**Release:** R64-5-1 (May 2024)
**Source:** [SGD Downloads](https://downloads.yeastgenome.org/sequence/S288C_reference/)

The reference genome contains 16 nuclear chromosomes plus the mitochondrial genome:
- chrI through chrXVI (nuclear)
- chrMito (mitochondrial)

### GFF3 Annotations (Required)

Contains gene, CDS, mRNA, and other feature annotations for all ~6,000 yeast genes.

Key feature types:
- `gene`: Gene loci
- `CDS`: Coding sequences (for protein-coding genes)
- `mRNA`: Transcripts
- `intron`: Intron positions (for ~300 intron-containing genes)

### ORF Sequences (Optional)

- `orf_coding_all`: DNA coding sequences for all verified/uncharacterized ORFs
- `orf_trans_all`: Protein sequences (useful for ESM1-v validation)

## Coordinate Systems

**IMPORTANT:** The VEP pipeline uses **1-based, inclusive** coordinates consistent with GFF3 format.

- Position 1 = first nucleotide of chromosome
- Start and end positions are both included in the feature
- This matches SGD, UCSC, and standard GFF3 conventions

## MAGESTIC Background Strain

For MAGESTIC-specific analyses, you also need the MAGESTIC background strain genome, which differs from S288C at ~few hundred positions. Contact the Steinmetz lab for access.

## Updating References

When SGD releases a new genome version:

1. Update the `RELEASE` variable in `download_references.sh`
2. Run the download script
3. Update any hardcoded paths in your config files

## Troubleshooting

### "File not found" errors

Make sure to decompress the `.gz` files:
```bash
gunzip -k *.gz
```

### FASTA index missing

Create with samtools:
```bash
module load samtools  # On Sherlock
samtools faidx S288C_reference_sequence_R64-5-1_20240529.fsa
```

### Permission denied

Files may be read-only. Copy to a writable location:
```bash
cp -r reference_data/ /path/to/your/project/
chmod -R u+w /path/to/your/project/reference_data/
```
