"""
GFF/GTF annotation parsing for MAGESTIC design workflows.

Extracts the functions from the legacy GTF_GFF_manipulation.py that are
required for guide-donor oligo design: ORF coordinate loading, codon
coordinate mapping, and gene name lookups.

Legacy functions involving TIF-seq, bedgraph analysis, and RNA-seq annotation
types are NOT ported here — those were specific to a different analysis context.
The modern GFF loader with GeneFeature dataclass is in magestic.utils.gff_utils.

Author: Kevin R. Roy
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from .sequence import rev_comp

logger = logging.getLogger(__name__)

# Chromosome numeral lookup (Roman numeral chromosomes for S. cerevisiae)
_CHROMOSOME_NUMERALS = [
    'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII',
    'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI',
]
_NUMBER_TO_NUMERAL = {f'chr{i + 1}': f'chr{n}' for i, n in enumerate(_CHROMOSOME_NUMERALS)}


def convert_chromosome_number_to_numeral(chromosome_number: str) -> str:
    """
    Convert a numbered chromosome name to Roman numeral form.

    Args:
        chromosome_number: Chromosome name like 'chr1', 'chr16'

    Returns:
        Roman numeral chromosome name like 'chrI', 'chrXVI'

    Raises:
        KeyError: If the chromosome number is not in the S. cerevisiae range (chr1-chr16)
    """
    return _NUMBER_TO_NUMERAL[chromosome_number]


# ---------------------------------------------------------------------------
# GFF attribute field parsers
# ---------------------------------------------------------------------------

def get_annotation_fields_from_gff(annotation: str) -> Dict[str, str]:
    """
    Parse a GFF3 attribute column into a key-value dictionary.

    Handles URL-encoded values and semicolon-separated fields.

    Args:
        annotation: GFF3 attribute column string
            (e.g., 'ID=YAL001C;Name=TFC3;gene=TFC3;...')

    Returns:
        Dict mapping attribute names to their values
    """
    from urllib.parse import unquote
    fields: Dict[str, str] = {}
    for field in annotation.strip().split(';'):
        field = field.strip()
        if '=' in field:
            key, _, value = field.partition('=')
            fields[key.strip()] = unquote(value.strip())
    return fields


def _get_field(annotation: str, key: str) -> Optional[str]:
    """Extract a single GFF attribute field value by key."""
    fields = get_annotation_fields_from_gff(annotation)
    return fields.get(key)


def get_systematic_gene_name(annotation: str) -> str:
    """
    Extract the systematic gene name (ID=) from a GFF gene annotation string.

    Args:
        annotation: GFF3 attribute column string

    Returns:
        Systematic gene name (e.g., 'YAL001C')

    Raises:
        ValueError: If 'ID=' is not found in the annotation string
    """
    name = _get_field(annotation, 'ID')
    if name is None:
        raise ValueError(f"'ID=' not found in annotation: {annotation!r}")
    return name


def get_common_gene_name(annotation: str) -> str:
    """
    Extract the common gene name from a GFF annotation string.

    Prefers 'gene=' field, falls back to 'Name='.

    Args:
        annotation: GFF3 attribute column string

    Returns:
        Common gene name (e.g., 'TFC3'), or systematic name if common not available
    """
    fields = get_annotation_fields_from_gff(annotation)
    return fields.get('gene') or fields.get('Name') or fields.get('ID', '')


def get_gene_id_from_gff_annotation(annotation: str) -> Optional[str]:
    """Extract ID= field from GFF annotation. Returns None if not present."""
    return _get_field(annotation, 'ID')


def get_parent_from_gff_annotation(annotation: str) -> Optional[str]:
    """Extract Parent= field from GFF annotation. Returns None if not present."""
    return _get_field(annotation, 'Parent')


def get_gene_name_from_gff_annotation(annotation: str) -> Optional[str]:
    """Extract Name= field from GFF annotation. Returns None if not present."""
    return _get_field(annotation, 'Name')


# ---------------------------------------------------------------------------
# Gene name dictionary loaders
# ---------------------------------------------------------------------------

def load_systematic_to_common_gene_name_dict(
    gff_filename: Union[str, Path],
) -> Dict[str, str]:
    """
    Build a systematic → common gene name mapping from a GFF3 file.

    Args:
        gff_filename: Path to GFF3 annotation file

    Returns:
        Dict mapping systematic names (e.g., 'YAL001C') to common names
        (e.g., 'TFC3'). Genes without a common name map to their own
        systematic name.
    """
    result: Dict[str, str] = {}
    with open(gff_filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, _, feature_type = parts[0], parts[1], parts[2]
            if chrom in ('chrMito', 'chrmt'):
                continue
            if feature_type == 'gene':
                annotation = parts[8]
                systematic = get_systematic_gene_name(annotation)
                common = get_common_gene_name(annotation)
                result[systematic] = common
    return result


def load_common_to_systematic_gene_name_dict(
    gff_filename: Union[str, Path],
) -> Dict[str, str]:
    """
    Build a common → systematic gene name mapping from a GFF3 file.

    Args:
        gff_filename: Path to GFF3 annotation file

    Returns:
        Dict mapping common names (e.g., 'TFC3') to systematic names
        (e.g., 'YAL001C').
    """
    result: Dict[str, str] = {}
    with open(gff_filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, _, feature_type = parts[0], parts[1], parts[2]
            if chrom in ('chrMito', 'chrmt'):
                continue
            if feature_type == 'gene':
                annotation = parts[8]
                systematic = get_systematic_gene_name(annotation)
                common = get_common_gene_name(annotation)
                result[common] = systematic
    return result


# ---------------------------------------------------------------------------
# ORF coordinate loading
# ---------------------------------------------------------------------------

def load_ORF_coordinates(
    gff_filename: Union[str, Path],
) -> Dict[str, List[Tuple[str, int, int, str]]]:
    """
    Load CDS coordinates for all ORFs from a GFF3 file.

    For multi-exon genes, each exon appears as a separate tuple in the list.
    The '_CDS' suffix is stripped from feature names for consistency with
    the MAGESTIC design scripts.

    Args:
        gff_filename: Path to GFF3 annotation file

    Returns:
        Dict mapping ORF names (e.g., 'YAL001C_CDS' → 'YAL001C') to lists
        of (chrom, start_1based, end_1based, strand) tuples.
    """
    orf_dict: Dict[str, List[Tuple[str, int, int, str]]] = {}
    with open(gff_filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, _, feature_type, start, end, _, strand, _, annotation = parts
            if feature_type != 'CDS' or chrom == 'chrmt':
                continue
            start_i, end_i = int(start), int(end)
            gene_name_raw = get_gene_name_from_gff_annotation(annotation)
            if gene_name_raw is None:
                continue
            # Strip _CDS suffix if present (legacy format: YAL001C_CDS)
            gene_name = gene_name_raw[:-4] if gene_name_raw.endswith('_CDS') else gene_name_raw
            if gene_name not in orf_dict:
                orf_dict[gene_name] = []
            orf_dict[gene_name].append((chrom, start_i, end_i, strand))
    return orf_dict


# ---------------------------------------------------------------------------
# ORF sequence extraction and codon mapping
# ---------------------------------------------------------------------------

def get_ORF_seq_and_CDS_coords(
    orf: str,
    orf_dict: Dict[str, List[Tuple[str, int, int, str]]],
    genome_seq: Dict[str, str],
) -> Tuple[str, str, List[int], List[str]]:
    """
    Extract the codon sequence and genomic coordinates for an ORF.

    Handles multi-exon ORFs and correctly deals with reverse-strand genes.

    Args:
        orf: ORF name (key in orf_dict)
        orf_dict: ORF coordinate dictionary from :func:`load_ORF_coordinates`
        genome_seq: Genome sequence dict from :func:`magestic.genomics.genome_io.load_custom_genome`

    Returns:
        Tuple of (chrom, strand, exon_coords, codon_list) where:
          - chrom: Chromosome name
          - strand: '+' or '-'
          - exon_coords: List of 0-based genomic positions for each nucleotide
            in the ORF (length = 3 × number of codons)
          - codon_list: List of 3-letter codon strings

    Notes:
        - Uses 1-based coordinates from GFF internally, converts to 0-based
          for exon_coords output (consistent with Python slicing conventions).
        - Stop codons are included in the returned sequence.
    """
    exons = orf_dict[orf]
    orf_seq: List[str] = []
    exon_coords: List[int] = []
    strand = exons[0][3]
    frame_offset = 0
    stop_codons = {'TAA', 'TAG', 'TGA'}

    if strand == '+':
        for exon in exons:
            chrom, start, end, _ = exon
            # GFF is 1-based inclusive; Python slice is 0-based, exclusive end
            exon_seq = genome_seq[chrom][start - 1:end]
            if frame_offset > 0:
                prev = orf_seq[-1]
                rest = exon_seq[:frame_offset]
                orf_seq[-1] = prev + rest
                for idx in range(frame_offset):
                    exon_coords.append(start - 1 + idx)
            for idx in range(frame_offset, len(exon_seq), 3):
                frame_offset = 0
                codon = exon_seq[idx:idx + 3]
                orf_seq.append(codon)
                for exon_idx in range(len(codon)):
                    exon_coords.append(start - 1 + exon_idx + idx)
                if len(codon) < 3:
                    frame_offset = 3 - len(codon)
    else:
        for exon in reversed(exons):
            chrom, start, end, _ = exon
            exon_seq = rev_comp(genome_seq[chrom][start - 1:end])
            if frame_offset > 0:
                prev = orf_seq[-1]
                rest = exon_seq[:frame_offset]
                orf_seq[-1] = prev + rest
                for idx in range(frame_offset):
                    exon_coords.append(end - idx)
            for idx in range(frame_offset, len(exon_seq), 3):
                frame_offset = 0
                codon = exon_seq[idx:idx + 3]
                orf_seq.append(codon)
                for exon_idx in range(len(codon)):
                    exon_coords.append(end - exon_idx - idx)
                if len(codon) < 3:
                    frame_offset = 3 - len(codon)

    stop_count = sum(c in stop_codons for c in orf_seq)
    if stop_count != 1:
        logger.warning(
            "%s has %d stop codons (expected 1). Exons: %s", orf, stop_count, exons
        )

    return chrom, strand, exon_coords, orf_seq


def get_aa_num_to_codon_coords_aa(
    orf_chrom: str,
    orf_strand: str,
    orf_exon_coords: List[int],
    orf_seq: List[str],
    codon_to_aa: Dict[str, str],
) -> Tuple[Dict[int, Tuple[str, Tuple[int, ...], str]], Dict[int, List[int]]]:
    """
    Build per-amino-acid codon and coordinate dictionaries for an ORF.

    Args:
        orf_chrom: Chromosome name (unused but kept for API compatibility)
        orf_strand: Strand '+' or '-'
        orf_exon_coords: List of 0-based genomic coordinates (one per nucleotide)
            from :func:`get_ORF_seq_and_CDS_coords`
        orf_seq: List of codon strings from :func:`get_ORF_seq_and_CDS_coords`
        codon_to_aa: Codon to amino acid mapping dict

    Returns:
        Tuple of (orf_info, aa_num_to_exon) where:
          - orf_info: Dict mapping 1-based AA position → (codon, genomic_coords, aa)
          - aa_num_to_exon: Dict mapping 1-based AA position → list of exon numbers
            (1-based). Most AAs map to [1], split AAs map to [n, n+1].
    """
    current_exon_num = 1
    aa_num_to_exon: Dict[int, List[int]] = {}
    aa_length = len(orf_seq)
    previous_aa_end_coord: Optional[int] = None

    for aa_num in range(1, aa_length + 1):
        current_aa_coords = orf_exon_coords[(aa_num - 1) * 3: aa_num * 3]
        if previous_aa_end_coord is not None and abs(current_aa_coords[0] - previous_aa_end_coord) != 1:
            current_exon_num += 1
        aa_num_to_exon[aa_num] = [current_exon_num]
        if abs(current_aa_coords[2] - current_aa_coords[0]) != 2:
            aa_num_to_exon[aa_num] = [current_exon_num, current_exon_num + 1]
            current_exon_num += 1
        previous_aa_end_coord = current_aa_coords[2]

    orf_info: Dict[int, Tuple[str, Tuple[int, ...], str]] = {}
    for idx in range(len(orf_seq)):
        aa_num = idx + 1
        codon = orf_seq[idx].upper()
        aa = codon_to_aa[codon]
        aa_coords = tuple(orf_exon_coords[idx * 3: (idx + 1) * 3])
        orf_info[aa_num] = (codon, aa_coords, aa)

    return orf_info, aa_num_to_exon
