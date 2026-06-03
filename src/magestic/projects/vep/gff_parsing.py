"""
GFF3 Parsing Utilities for VEP Pipeline

Consolidated functions for parsing GFF3 files with proper URL decoding.
GFF3 spec requires URL-encoding for special characters (RFC 3986).

Functions:
- parse_gff_attributes: Parse GFF9 attributes with URL decoding
- load_gene_annotations: Load gene features from S288C GFF
- load_cds_info: Load CDS features with strain-specific info
- load_s288c_gene_structures: Load exon structure for multi-exon genes

Classes:
- GeneInfo: Gene annotation data
- CDSInfo: CDS annotation with strain coordinates
- ExonInfo: Single exon coordinates
- GeneStructure: Complete gene structure with exons

Author: Kevin R. Roy
Date: 2026-02-26
"""

import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union
from urllib.parse import unquote

# =============================================================================
# Constants
# =============================================================================

# Set VEP_BASE_DIR environment variable to override
import os
BASE_DIR = Path(os.environ.get("VEP_BASE_DIR", Path(__file__).parent))

# Default GFF paths
DEFAULT_S288C_GFF = (
    BASE_DIR / "common/reference_genomes/S288C_reference_genome_R64-5-1_20240529"
    / "saccharomyces_cerevisiae_R64-5-1_20240529.gff"
)

DEFAULT_STRAIN_GFF_DIR = (
    BASE_DIR / "projects/QTL/20260210_variant_effect_prediction"
    / "processed_data/20260205_in_silico_QTN_dissection/strain_genomes/gff"
)


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class GeneInfo:
    """Gene annotation data from GFF."""
    orf: str                    # Systematic name (e.g., YNL085W)
    common_name: str            # Common name (e.g., MKT1)
    chrom: str                  # Chromosome (e.g., chrXIV)
    start: int                  # Gene start (1-based)
    end: int                    # Gene end (1-based)
    strand: str                 # '+' or '-'
    cds_start: Optional[int] = None  # CDS start if known
    cds_end: Optional[int] = None    # CDS end if known


@dataclass
class CDSInfo:
    """CDS annotation with strain-specific coordinates."""
    orf: str                    # Systematic name
    strain: str                 # Strain name
    chrom: str                  # Chromosome
    strain_start: int           # Start in strain coordinates
    strain_end: int             # End in strain coordinates
    strand: str                 # '+' or '-'
    s288c_start: int            # Corresponding S288C start
    s288c_end: int              # Corresponding S288C end
    offset_start: int = 0       # Offset at start position
    offset_end: int = 0         # Offset at end position


@dataclass
class ExonInfo:
    """Single exon coordinates."""
    chrom: str
    start: int      # 1-based
    end: int        # 1-based
    strand: str
    phase: int = 0  # Reading frame (0, 1, or 2)


@dataclass
class GeneStructure:
    """Complete gene structure including exons."""
    orf: str
    common_name: str
    chrom: str
    strand: str
    gene_start: int
    gene_end: int
    exons: List[ExonInfo] = field(default_factory=list)

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    @property
    def is_spliced(self) -> bool:
        return len(self.exons) > 1

    @property
    def cds_length(self) -> int:
        """Total CDS length (sum of exon lengths)."""
        return sum(e.end - e.start + 1 for e in self.exons)


# =============================================================================
# Core Functions
# =============================================================================

def parse_gff_attributes(attr_string: str, url_decode: bool = True) -> Dict[str, str]:
    """Parse GFF9 attributes column into dictionary.

    GFF3 spec requires URL encoding for special characters:
    - %2C for comma, %3B for semicolon, %3D for equals, %26 for ampersand

    Args:
        attr_string: Attributes column from GFF (column 9)
        url_decode: If True, URL-decode values (recommended for GFF3)

    Returns:
        Dictionary of attribute key-value pairs
    """
    attrs = {}
    for item in attr_string.strip().split(';'):
        if not item.strip():
            continue
        if '=' in item:
            key, value = item.split('=', 1)
            key = key.strip()
            value = value.strip()
            if url_decode:
                # URL decode both key and value
                key = unquote(key)
                value = unquote(value)
            attrs[key] = value
    return attrs


def load_gene_annotations(
    gff_path: Optional[Path] = None,
    target_orfs: Optional[Union[List[str], Set[str]]] = None
) -> Dict[str, GeneInfo]:
    """Load gene annotations from GFF file.

    Args:
        gff_path: Path to GFF file. If None, uses DEFAULT_S288C_GFF.
        target_orfs: Optional list/set of ORFs to load. If None, loads all.

    Returns:
        Dictionary mapping ORF systematic name to GeneInfo.
    """
    if gff_path is None:
        gff_path = DEFAULT_S288C_GFF

    gff_path = Path(gff_path)
    target_set = set(target_orfs) if target_orfs else None

    genes = {}
    cds_coords = defaultdict(list)

    # Handle gzipped files
    is_gzipped = str(gff_path).endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'

    with opener(gff_path, mode, encoding='utf-8' if not is_gzipped else None) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attrs = parts
            start, end = int(start), int(end)
            attr_dict = parse_gff_attributes(attrs)

            orf_id = attr_dict.get('ID', attr_dict.get('Parent', ''))

            # Handle gene features
            if feature_type == 'gene':
                if target_set is None or orf_id in target_set:
                    common_name = attr_dict.get('gene', attr_dict.get('Name', orf_id))
                    genes[orf_id] = GeneInfo(
                        orf=orf_id,
                        common_name=common_name,
                        chrom=chrom,
                        start=start,
                        end=end,
                        strand=strand
                    )

            # Handle CDS features
            if feature_type == 'CDS':
                parent = attr_dict.get('Parent', '')
                # Parent might be "YNL085W_mRNA" - extract ORF
                for key in [parent, parent.replace('_mRNA', '')]:
                    if target_set is None or key in target_set:
                        cds_coords[key].append((start, end))

    # Add CDS coordinates to gene info
    for orf, coords in cds_coords.items():
        if orf in genes:
            genes[orf].cds_start = min(c[0] for c in coords)
            genes[orf].cds_end = max(c[1] for c in coords)

    return genes


def load_cds_info(
    strain: str,
    gff_dir: Optional[Path] = None,
    target_orfs: Optional[Set[str]] = None
) -> Dict[str, CDSInfo]:
    """Load CDS entries from strain-specific GFF file.

    Args:
        strain: Strain name (e.g., 'RMx', 'BYa')
        gff_dir: Directory containing strain GFF files.
                Default: DEFAULT_STRAIN_GFF_DIR
        target_orfs: Optional set of ORFs to load. If None, loads all.

    Returns:
        Dictionary mapping ORF name to CDSInfo.
    """
    if gff_dir is None:
        gff_dir = DEFAULT_STRAIN_GFF_DIR

    gff_dir = Path(gff_dir)
    gff_file = gff_dir / f"{strain}_CDS.gff"

    if not gff_file.exists():
        raise FileNotFoundError(f"Strain GFF not found: {gff_file}")

    cds_entries = {}

    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attributes = parts[:9]
            if feature != 'CDS':
                continue

            attrs = parse_gff_attributes(attributes)
            orf_id = attrs.get('ID', '').replace('_CDS', '')

            if target_orfs and orf_id not in target_orfs:
                continue

            # Parse S288C coordinates and offsets from attributes
            s288c_start = int(attrs.get('s288c_start', start))
            s288c_end = int(attrs.get('s288c_end', end))
            offset_start = int(attrs.get('offset_start', 0))
            offset_end = int(attrs.get('offset_end', 0))

            cds_entries[orf_id] = CDSInfo(
                orf=orf_id,
                strain=strain,
                chrom=chrom,
                strain_start=int(start),
                strain_end=int(end),
                strand=strand,
                s288c_start=s288c_start,
                s288c_end=s288c_end,
                offset_start=offset_start,
                offset_end=offset_end
            )

    return cds_entries


@lru_cache(maxsize=4)
def load_s288c_gene_structures(
    gff_path: Optional[Path] = None
) -> Dict[str, GeneStructure]:
    """Load complete gene structures including exon information.

    Useful for multi-exon (spliced) genes where CDS coordinates
    don't form a contiguous block.

    Args:
        gff_path: Path to S288C GFF file.

    Returns:
        Dictionary mapping ORF name to GeneStructure.
    """
    if gff_path is None:
        gff_path = DEFAULT_S288C_GFF

    gff_path = Path(gff_path)

    genes = {}
    exons = defaultdict(list)

    is_gzipped = str(gff_path).endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'

    with opener(gff_path, mode, encoding='utf-8' if not is_gzipped else None) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attrs = parts
            start, end = int(start), int(end)
            phase = int(phase) if phase != '.' else 0
            attr_dict = parse_gff_attributes(attrs)

            orf_id = attr_dict.get('ID', attr_dict.get('Parent', ''))

            if feature_type == 'gene':
                common_name = attr_dict.get('gene', attr_dict.get('Name', orf_id))
                genes[orf_id] = GeneStructure(
                    orf=orf_id,
                    common_name=common_name,
                    chrom=chrom,
                    strand=strand,
                    gene_start=start,
                    gene_end=end
                )

            elif feature_type == 'CDS':
                parent = attr_dict.get('Parent', '')
                orf = parent.replace('_mRNA', '')
                exons[orf].append(ExonInfo(
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    phase=phase
                ))

    # Attach exons to genes
    for orf, exon_list in exons.items():
        if orf in genes:
            # Sort exons by position
            genes[orf].exons = sorted(exon_list, key=lambda e: e.start)

    return genes


def get_spliced_genes(
    gff_path: Optional[Path] = None,
    min_exons: int = 2
) -> List[str]:
    """Get list of spliced genes (multi-exon).

    Args:
        gff_path: Path to GFF file.
        min_exons: Minimum number of exons to be considered spliced.

    Returns:
        List of ORF names for spliced genes.
    """
    structures = load_s288c_gene_structures(gff_path)
    return [orf for orf, gs in structures.items() if gs.n_exons >= min_exons]


def validate_gff_entry(parts: List[str]) -> Tuple[bool, Optional[str]]:
    """Validate a GFF entry has proper format.

    Args:
        parts: List of 9 tab-separated fields from GFF line.

    Returns:
        Tuple of (is_valid, error_message or None)
    """
    if len(parts) < 9:
        return False, f"Expected 9 fields, got {len(parts)}"

    chrom, source, feature_type, start, end, score, strand, phase, attrs = parts

    # Validate coordinates
    try:
        start_int = int(start)
        end_int = int(end)
        if start_int > end_int:
            return False, f"Start ({start}) > end ({end})"
        if start_int < 1:
            return False, f"Start ({start}) must be >= 1"
    except ValueError:
        return False, f"Invalid coordinates: start={start}, end={end}"

    # Validate strand
    if strand not in ['+', '-', '.']:
        return False, f"Invalid strand: {strand}"

    # Validate phase for CDS
    if feature_type == 'CDS' and phase not in ['0', '1', '2', '.']:
        return False, f"Invalid phase for CDS: {phase}"

    return True, None


# =============================================================================
# Test
# =============================================================================

if __name__ == "__main__":
    # Test parse_gff_attributes
    print("Testing parse_gff_attributes...")
    test_attrs = "ID=YNL085W;gene=MKT1;Note=Description%20with%20spaces"
    parsed = parse_gff_attributes(test_attrs)
    assert parsed['ID'] == 'YNL085W'
    assert parsed['gene'] == 'MKT1'
    assert parsed['Note'] == 'Description with spaces'  # URL decoded
    print("  parse_gff_attributes: PASS")

    # Test load_gene_annotations
    print("\nTesting load_gene_annotations...")
    if DEFAULT_S288C_GFF.exists():
        genes = load_gene_annotations(target_orfs=['YNL085W', 'YBR135W'])
        print(f"  Loaded {len(genes)} genes")
        if 'YNL085W' in genes:
            g = genes['YNL085W']
            print(f"  YNL085W: {g.common_name}, {g.chrom}:{g.start}-{g.end} ({g.strand})")
        print("  load_gene_annotations: PASS")
    else:
        print(f"  Skipping - GFF not found: {DEFAULT_S288C_GFF}")

    # Test load_s288c_gene_structures
    print("\nTesting load_s288c_gene_structures...")
    if DEFAULT_S288C_GFF.exists():
        structures = load_s288c_gene_structures()
        spliced = [orf for orf, gs in structures.items() if gs.is_spliced]
        print(f"  Loaded {len(structures)} genes, {len(spliced)} spliced")
        if spliced:
            example = structures[spliced[0]]
            print(f"  Example spliced gene: {example.orf} ({example.n_exons} exons)")
        print("  load_s288c_gene_structures: PASS")
    else:
        print(f"  Skipping - GFF not found: {DEFAULT_S288C_GFF}")

    print("\nAll gff_parsing tests completed!")
