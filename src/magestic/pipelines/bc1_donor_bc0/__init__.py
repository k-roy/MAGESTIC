"""
magestic.pipelines.bc1_donor_bc0 — BC1-donor-BC0 integration mapping pipeline.

Links BC1 barcodes from screen sequencing to guide-donor identities via
shared BC0 barcodes. This pipeline:
- Parses BC1/donor/BC0 sequences from sequencing reads
- Links BC1 barcodes to guide-donor-BC0 plasmid libraries
- Maps to designed oligos
- Assigns sublibraries and annotates unmatched donors

Migrated from common/scripts/bc1_donor_bc0_pipeline/ (v1.1.0).

Author: Kevin R. Roy
"""
