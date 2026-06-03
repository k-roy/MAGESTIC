"""
Genome FASTA loading utilities.

Refactored from the legacy bedgraph_computation.py module.
Removes Python 2 artifacts (imp.reload, circular imports, print statements).

Author: Kevin R. Roy
"""

import logging
from pathlib import Path
from typing import Dict, Union

logger = logging.getLogger(__name__)


def load_genome(
    genome_file: Union[str, Path],
    include_mito: bool = False,
) -> Dict[str, str]:
    """
    Load a genome FASTA file into memory.

    Stops loading at the mitochondrial chromosome (any chromosome whose name
    contains 'm' or 'M') unless include_mito=True. This matches the legacy
    behavior for S. cerevisiae genomes where the mitochondrial chromosome
    typically comes last.

    For custom / non-yeast genomes, use :func:`load_custom_genome` instead,
    which loads all chromosomes unconditionally.

    Args:
        genome_file: Path to the genome FASTA file
        include_mito: If True, do not stop at mitochondrial chromosomes.
            Default is False (matches legacy behavior).

    Returns:
        Dict mapping chromosome names (FASTA headers without '>') to
        their uppercase sequence strings.
    """
    genome_file = Path(genome_file)
    genome_dict: Dict[str, str] = {}
    current_chromosome: str = ""

    with open(genome_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                chrom_name = line[1:]
                if not include_mito and ('m' in chrom_name or 'M' in chrom_name):
                    break
                current_chromosome = chrom_name
                genome_dict[current_chromosome] = ''
            else:
                genome_dict[current_chromosome] += line.upper()

    logger.info("Loaded genome from %s (%d chromosomes)", genome_file, len(genome_dict))
    return genome_dict


def load_custom_genome(
    genome_file: Union[str, Path],
) -> Dict[str, str]:
    """
    Load all chromosomes from a FASTA file unconditionally.

    Unlike :func:`load_genome`, this function does not skip mitochondrial
    or other chromosomes — it loads everything. Suitable for non-yeast
    genomes or any FASTA where all sequences should be retained.

    Args:
        genome_file: Path to the FASTA file

    Returns:
        Dict mapping chromosome names (FASTA headers without '>') to
        their uppercase sequence strings.
    """
    genome_file = Path(genome_file)
    genome_dict: Dict[str, str] = {}
    current_chromosome: str = ""

    with open(genome_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_chromosome = line[1:]
                genome_dict[current_chromosome] = ''
            else:
                genome_dict[current_chromosome] += line.upper()

    logger.info("Loaded genome from %s (%d sequences)", genome_file, len(genome_dict))
    return genome_dict
