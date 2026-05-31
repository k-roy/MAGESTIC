"""
GFF file parsing utilities for yeast genome annotations.

This module provides functions for loading and parsing GFF3 format
annotation files, particularly for the MAGESTIC background strain.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import unquote


# Default GFF file path
DEFAULT_GFF_FILE = Path(
    "/path/to/common/annotation_files"
    "/MAGESTIC_background_strain_annotations.gff"
)


@dataclass
class GeneFeature:
    """
    Gene feature from GFF file.

    Attributes:
        chrom: Chromosome name (e.g., 'chrI', 'chrX')
        start: 1-based start position
        end: 1-based end position (inclusive)
        strand: '+' or '-'
        gene_id: Systematic gene ID (e.g., 'YGL003C')
        common_name: Common gene name if available (e.g., 'TOR1')
        feature_type: GFF feature type ('gene', 'CDS', etc.)
    """
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    common_name: Optional[str] = None
    feature_type: str = 'gene'

    @property
    def length(self) -> int:
        """Get the length of the feature in base pairs."""
        return self.end - self.start + 1

    def contains(self, position: int) -> bool:
        """Check if a position falls within this feature."""
        return self.start <= position <= self.end

    def overlaps(self, start: int, end: int) -> bool:
        """Check if a range overlaps with this feature."""
        return not (end < self.start or start > self.end)


def load_gff(
    gff_file: Path = DEFAULT_GFF_FILE,
    feature_types: Optional[List[str]] = None
) -> Dict[str, List[GeneFeature]]:
    """
    Load gene features from GFF file.

    Parses a GFF3 file and returns gene/CDS features organized by chromosome.
    Also builds a lookup table mapping systematic names to common names.

    Args:
        gff_file: Path to GFF3 file
        feature_types: List of feature types to load (default: ['gene', 'CDS'])

    Returns:
        Dict mapping chromosome names to lists of GeneFeature objects.
        Includes special key '_gene_name_lookup' with systematic->common name mapping.
    """
    if feature_types is None:
        feature_types = ['gene', 'CDS']

    features: Dict[str, List[GeneFeature]] = {}
    gene_name_lookup: Dict[str, str] = {}

    if not gff_file.exists():
        print(f"Warning: GFF file not found: {gff_file}")
        return features

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feature_type = parts[2]

            if feature_type not in feature_types:
                continue

            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            # Parse attributes
            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attrs[key] = unquote(value)

            gene_id = attrs.get('ID', attrs.get('Parent', '').replace('_mRNA', ''))
            common_name = attrs.get('gene', '')  # Only 'gene' features have this

            # For 'gene' features, build the lookup table
            if feature_type == 'gene' and common_name:
                systematic_name = attrs.get('ID', '')
                if systematic_name:
                    gene_name_lookup[systematic_name] = common_name
                    # Also map with _CDS suffix for CDS lookups
                    gene_name_lookup[f"{systematic_name}_CDS"] = common_name

            feature = GeneFeature(
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                gene_id=gene_id,
                common_name=common_name if common_name else None,
                feature_type=feature_type
            )

            if chrom not in features:
                features[chrom] = []
            features[chrom].append(feature)

    # Store lookup table in a special key
    features['_gene_name_lookup'] = gene_name_lookup  # type: ignore

    return features


def get_gene_by_name(
    features: Dict[str, List[GeneFeature]],
    gene_name: str
) -> Optional[GeneFeature]:
    """
    Find a gene by common name or systematic name.

    Args:
        features: Dict returned by load_gff()
        gene_name: Gene name to search for (common or systematic)

    Returns:
        GeneFeature if found, None otherwise
    """
    # Get the lookup table
    lookup = features.get('_gene_name_lookup', {})

    # Try to find the systematic name if given common name
    systematic_name = None
    for sys_name, common in lookup.items():
        if common == gene_name or sys_name == gene_name:
            systematic_name = sys_name
            break

    if systematic_name is None:
        systematic_name = gene_name

    # Search through all chromosomes
    for chrom, chrom_features in features.items():
        if chrom.startswith('_'):  # Skip special keys
            continue
        if not isinstance(chrom_features, list):
            continue
        for feature in chrom_features:
            if feature.feature_type == 'gene':
                if feature.gene_id == systematic_name or feature.common_name == gene_name:
                    return feature

    return None


def get_genes_in_region(
    features: Dict[str, List[GeneFeature]],
    chrom: str,
    start: int,
    end: int,
    feature_type: str = 'gene'
) -> List[GeneFeature]:
    """
    Get all genes/features overlapping a genomic region.

    Args:
        features: Dict returned by load_gff()
        chrom: Chromosome name
        start: Region start position
        end: Region end position
        feature_type: Feature type to return ('gene', 'CDS', or None for all)

    Returns:
        List of GeneFeature objects overlapping the region
    """
    if chrom not in features:
        return []

    result = []
    for feature in features[chrom]:
        if not isinstance(feature, GeneFeature):
            continue
        if feature_type and feature.feature_type != feature_type:
            continue
        if feature.overlaps(start, end):
            result.append(feature)

    return result


def get_systematic_to_common_name_map(
    features: Dict[str, List[GeneFeature]]
) -> Dict[str, str]:
    """
    Get the systematic name to common name mapping.

    Args:
        features: Dict returned by load_gff()

    Returns:
        Dict mapping systematic gene IDs to common names
    """
    return features.get('_gene_name_lookup', {})
