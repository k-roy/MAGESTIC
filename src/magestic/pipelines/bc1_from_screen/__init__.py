"""
magestic.pipelines.bc1_from_screen — BC1 screen analysis pipeline.

A unified pipeline for processing BC1 barcode counts from pooled CRISPR screens.
This pipeline:
- Processes collapsed FASTQ files to extract BC1 barcodes
- Classifies BC1s into tiers based on t0 read counts
- Creates DESeq2-ready count matrices
- Runs differential expression analysis (DESeq2)
- Performs tier-aware hit calling

Supports library-specific processing:
- Bottlenecked vs non-bottlenecked libraries
- Multi-timepoint vs single-timepoint screens

Migrated from common/scripts/bc1_from_screen_pipeline/ (v1.0.0).

Author: Kevin R. Roy
"""
