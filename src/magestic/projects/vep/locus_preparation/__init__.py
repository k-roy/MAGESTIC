"""
Locus Preparation Module

Functions for building genomic loci and inserting variants:
- locus_builder: Build genomic loci for any strain
- variant_insertion: Insert variants into strain backgrounds
- coordinate_mapping: Handle offsets between strain assemblies
- variant_annotation: Annotate variants (codon changes, consequences)
- validation: GFF parsing, bcftools consensus checks

Author: Kevin R. Roy
"""

# Will be populated as functions are extracted
