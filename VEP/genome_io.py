"""
Genome I/O Utilities for VEP Pipeline

Consolidated functions for loading reference genomes and strain sequences.
This module centralizes genome loading to eliminate duplication across scripts.

Functions:
- load_s288c_genome: Load S288C reference genome
- load_magestic_genome: Load MAGESTIC background strain genome
- load_strain_chromosome: Load single chromosome from strain genome
- save_chromosome_fasta: Save chromosome sequence to FASTA
- load_strain_genome: Load complete strain genome

Constants:
- DEFAULT_S288C_PATH: Path to S288C reference FASTA
- DEFAULT_MAGESTIC_PATH: Path to MAGESTIC background strain FASTA
- CHROMOSOMES: List of yeast chromosome names
- CHROM_TO_NCBI: Mapping from chrI format to NC_* accession

Author: Kevin R. Roy
Date: 2026-02-26
"""

import gzip
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional, Union

# =============================================================================
# Constants
# =============================================================================

BASE_DIR = Path("/path/to")

# Default reference genome paths
DEFAULT_S288C_PATH = (
    BASE_DIR / "common/reference_genomes/S288C_reference_genome_R64-5-1_20240529"
    / "S288C_reference_sequence_R64-5-1_20240529.fsa"
)
DEFAULT_MAGESTIC_PATH = (
    BASE_DIR / "common/reference_genomes/MAGESTIC_background_strain.fasta"
)

# Default strain genomes directory (built from VCF)
DEFAULT_STRAIN_GENOMES_DIR = (
    BASE_DIR / "projects/QTL/20260210_variant_effect_prediction"
    / "processed_data/20260205_in_silico_QTN_dissection/strain_genomes"
)

# Yeast chromosome names (VCF format)
CHROMOSOMES = [
    'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII',
    'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI',
    'chrMito'
]

# Chromosome name mapping (S288C FASTA uses ref|NC_*| format)
CHROM_TO_NCBI = {
    'chrI': 'NC_001133', 'chrII': 'NC_001134', 'chrIII': 'NC_001135',
    'chrIV': 'NC_001136', 'chrV': 'NC_001137', 'chrVI': 'NC_001138',
    'chrVII': 'NC_001139', 'chrVIII': 'NC_001140', 'chrIX': 'NC_001141',
    'chrX': 'NC_001142', 'chrXI': 'NC_001143', 'chrXII': 'NC_001144',
    'chrXIII': 'NC_001145', 'chrXIV': 'NC_001146', 'chrXV': 'NC_001147',
    'chrXVI': 'NC_001148', 'chrMito': 'NC_001224'
}

# Reverse mapping
NCBI_TO_CHROM = {v: k for k, v in CHROM_TO_NCBI.items()}

# Bloom 16 strains
BLOOM_STRAINS = [
    'RMx', 'BYa', 'M22', 'CLIB219x', 'CBS2888a', 'YJM981x',
    '273614xa', 'PW5a', 'Y10x', 'I14a', 'YPS1009x',
    'YJM454a', 'YJM978x', 'CLIB413a', 'YJM145x', 'YPS163a'
]


# =============================================================================
# Core Functions
# =============================================================================

def load_s288c_genome(
    fasta_path: Optional[Path] = None,
    use_cache: bool = True
) -> Dict[str, str]:
    """Load S288C reference genome from FASTA.

    Args:
        fasta_path: Path to S288C FASTA file. If None, uses DEFAULT_S288C_PATH.
        use_cache: If True, cache the result for repeated calls with same path.

    Returns:
        Dictionary mapping chromosome names (chrI, chrII, ...) to sequences.
    """
    if fasta_path is None:
        fasta_path = DEFAULT_S288C_PATH

    fasta_path = Path(fasta_path)

    if use_cache:
        return _load_s288c_genome_cached(fasta_path)
    else:
        return _load_s288c_genome_impl(fasta_path)


@lru_cache(maxsize=4)
def _load_s288c_genome_cached(fasta_path: Path) -> Dict[str, str]:
    """Cached implementation of load_s288c_genome."""
    return _load_s288c_genome_impl(fasta_path)


def _load_s288c_genome_impl(fasta_path: Path) -> Dict[str, str]:
    """Implementation of load_s288c_genome.

    Handles both gzipped and plain FASTA files.
    Converts NC_* accessions to chrI format.
    """
    genome = {}
    current_chrom = None
    current_seq = []

    # Handle both gzipped and plain FASTA
    is_gzipped = str(fasta_path).endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'

    with opener(fasta_path, mode, encoding='utf-8' if not is_gzipped else None) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom is not None:
                    genome[current_chrom] = ''.join(current_seq).upper()
                    current_seq = []

                # Parse header: >ref|NC_001133| [org=...] [chromosome=I]
                header = line[1:]
                parts = header.split()
                ncbi_id = parts[0]

                # Extract NC_* ID from ref|NC_*| format
                if ncbi_id.startswith('ref|'):
                    ncbi_id = ncbi_id.split('|')[1]

                # Convert NC_* to chrI format
                if ncbi_id in NCBI_TO_CHROM:
                    current_chrom = NCBI_TO_CHROM[ncbi_id]
                else:
                    # Use as-is if not in mapping
                    current_chrom = ncbi_id
            else:
                current_seq.append(line)

    # Don't forget the last chromosome
    if current_chrom is not None:
        genome[current_chrom] = ''.join(current_seq).upper()

    return genome


def load_magestic_genome(
    fasta_path: Optional[Path] = None,
    use_cache: bool = True
) -> Dict[str, str]:
    """Load MAGESTIC background strain genome from FASTA.

    Args:
        fasta_path: Path to MAGESTIC FASTA file. If None, uses DEFAULT_MAGESTIC_PATH.
        use_cache: If True, cache the result for repeated calls.

    Returns:
        Dictionary mapping chromosome names (chrI, chrII, ...) to sequences.
    """
    if fasta_path is None:
        fasta_path = DEFAULT_MAGESTIC_PATH

    fasta_path = Path(fasta_path)

    if use_cache:
        return _load_magestic_genome_cached(fasta_path)
    else:
        return _load_simple_fasta(fasta_path)


@lru_cache(maxsize=4)
def _load_magestic_genome_cached(fasta_path: Path) -> Dict[str, str]:
    """Cached implementation of load_magestic_genome."""
    return _load_simple_fasta(fasta_path)


def _load_simple_fasta(fasta_path: Path) -> Dict[str, str]:
    """Load FASTA with simple >name headers."""
    genome = {}
    current_chrom = None
    current_seq = []

    is_gzipped = str(fasta_path).endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'

    with opener(fasta_path, mode, encoding='utf-8' if not is_gzipped else None) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom is not None:
                    genome[current_chrom] = ''.join(current_seq).upper()
                    current_seq = []
                current_chrom = line[1:].split()[0]
            else:
                current_seq.append(line)

    if current_chrom is not None:
        genome[current_chrom] = ''.join(current_seq).upper()

    return genome


def load_strain_chromosome(
    strain: str,
    chrom: str,
    strain_genomes_dir: Optional[Path] = None
) -> str:
    """Load single chromosome sequence from strain genome.

    Args:
        strain: Strain name (e.g., 'RMx', 'BYa')
        chrom: Chromosome name (e.g., 'chrVII')
        strain_genomes_dir: Directory containing strain genome files.
                           Default: DEFAULT_STRAIN_GENOMES_DIR

    Returns:
        Chromosome sequence as string.

    Raises:
        FileNotFoundError: If chromosome file doesn't exist.
    """
    if strain_genomes_dir is None:
        strain_genomes_dir = DEFAULT_STRAIN_GENOMES_DIR

    strain_genomes_dir = Path(strain_genomes_dir)
    chrom_file = strain_genomes_dir / "chromosomes" / strain / f"{chrom}.fasta"

    if not chrom_file.exists():
        raise FileNotFoundError(
            f"Strain chromosome file not found: {chrom_file}\n"
            f"Make sure strain genomes have been built with 28_build_strain_chromosomes.py"
        )

    # Load single chromosome from FASTA
    with open(chrom_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Skip header, join sequence lines
    seq_lines = [line.strip() for line in lines[1:] if line.strip()]
    return ''.join(seq_lines).upper()


def load_strain_genome(
    strain: str,
    strain_genomes_dir: Optional[Path] = None,
    chromosomes: Optional[list] = None
) -> Dict[str, str]:
    """Load complete strain genome (all chromosomes).

    Args:
        strain: Strain name (e.g., 'RMx', 'BYa')
        strain_genomes_dir: Directory containing strain genome files.
        chromosomes: List of chromosomes to load. Default: all 17 chromosomes.

    Returns:
        Dictionary mapping chromosome names to sequences.
    """
    if chromosomes is None:
        chromosomes = CHROMOSOMES

    genome = {}
    for chrom in chromosomes:
        try:
            genome[chrom] = load_strain_chromosome(strain, chrom, strain_genomes_dir)
        except FileNotFoundError:
            # Skip missing chromosomes (e.g., chrMito)
            pass

    return genome


def save_chromosome_fasta(
    sequence: str,
    output_path: Path,
    header: str,
    line_width: int = 80
) -> None:
    """Save chromosome sequence to FASTA file.

    Args:
        sequence: Chromosome sequence
        output_path: Path for output FASTA file
        header: Header line (without '>') for FASTA
        line_width: Characters per line (default 80)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(f">{header}\n")
        # Write sequence in chunks
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i + line_width] + "\n")


def save_genome_fasta(
    genome: Dict[str, str],
    output_path: Path,
    line_width: int = 80
) -> None:
    """Save complete genome to FASTA file.

    Args:
        genome: Dictionary mapping chromosome names to sequences
        output_path: Path for output FASTA file
        line_width: Characters per line (default 80)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as f:
        for chrom in CHROMOSOMES:
            if chrom in genome:
                f.write(f">{chrom}\n")
                seq = genome[chrom]
                for i in range(0, len(seq), line_width):
                    f.write(seq[i:i + line_width] + "\n")


# =============================================================================
# Test
# =============================================================================

if __name__ == "__main__":
    # Test S288C genome loading
    print("Testing load_s288c_genome...")
    if DEFAULT_S288C_PATH.exists():
        genome = load_s288c_genome()
        print(f"  Loaded {len(genome)} chromosomes")
        for chrom in ['chrI', 'chrVII', 'chrXVI']:
            if chrom in genome:
                print(f"  {chrom}: {len(genome[chrom]):,} bp")
        print("  load_s288c_genome: PASS")
    else:
        print(f"  Skipping - S288C file not found: {DEFAULT_S288C_PATH}")

    # Test MAGESTIC genome loading
    print("\nTesting load_magestic_genome...")
    if DEFAULT_MAGESTIC_PATH.exists():
        genome = load_magestic_genome()
        print(f"  Loaded {len(genome)} chromosomes")
        print("  load_magestic_genome: PASS")
    else:
        print(f"  Skipping - MAGESTIC file not found: {DEFAULT_MAGESTIC_PATH}")

    # Test strain chromosome loading
    print("\nTesting load_strain_chromosome...")
    strain_dir = DEFAULT_STRAIN_GENOMES_DIR / "chromosomes" / "RMx"
    if strain_dir.exists():
        seq = load_strain_chromosome("RMx", "chrVII")
        print(f"  RMx chrVII: {len(seq):,} bp")
        print("  load_strain_chromosome: PASS")
    else:
        print(f"  Skipping - strain genomes not built yet")

    print("\nAll genome_io tests completed!")
