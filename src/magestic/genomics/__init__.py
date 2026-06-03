"""
magestic.genomics — Core genomics utilities.

Provides low-level sequence manipulation, genome I/O, and GFF annotation
parsing. These utilities are used throughout the MAGESTIC design and analysis
pipelines.

Modules:
    sequence   — DNA/RNA sequence utilities (rev_comp, hamming_distance, etc.)
    genome_io  — Genome FASTA loading
    annotation — GFF/GTF annotation parsing (ORF coordinates, gene names)
"""

from .sequence import (
    rev_comp,
    base_match,
    hamming_distance,
    edit_distance,
    translate_codon,
)
from .genome_io import load_genome, load_custom_genome
from .annotation import (
    convert_chromosome_number_to_numeral,
    load_ORF_coordinates,
    get_ORF_seq_and_CDS_coords,
    get_aa_num_to_codon_coords_aa,
    load_systematic_to_common_gene_name_dict,
    load_common_to_systematic_gene_name_dict,
    get_annotation_fields_from_gff,
)

__all__ = [
    "rev_comp",
    "base_match",
    "hamming_distance",
    "edit_distance",
    "translate_codon",
    "load_genome",
    "load_custom_genome",
    "convert_chromosome_number_to_numeral",
    "load_ORF_coordinates",
    "get_ORF_seq_and_CDS_coords",
    "get_aa_num_to_codon_coords_aa",
    "load_systematic_to_common_gene_name_dict",
    "load_common_to_systematic_gene_name_dict",
    "get_annotation_fields_from_gff",
]
