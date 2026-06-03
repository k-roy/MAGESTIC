"""
magestic.projects.nrd1 — NRD1 gene MAGESTIC guide/donor design.

Migrated from NRD1/ directory.

Bundled reference data (accessible via get_nrd1_data_path):
    DNA_codon_to_AA.txt                    — Full codon table
    DNA_codon_to_AA_suboptimal_codons_removed.txt  — Codon table without suboptimal codons
    NRD1_alanine_scan_with_PTC_and_syn_controls.tsv — NRD1 alanine scan reference

Usage:
    from magestic.projects.nrd1 import get_nrd1_data_path
    codon_table = get_nrd1_data_path("DNA_codon_to_AA.txt")

Author: Kevin R. Roy
"""

from importlib.resources import files
from pathlib import Path


def get_nrd1_data_path(filename: str) -> Path:
    """
    Return the Path to a bundled NRD1 reference data file.

    Args:
        filename: Name of the data file (e.g., "DNA_codon_to_AA.txt")

    Returns:
        Path to the bundled file.

    Raises:
        FileNotFoundError: If the filename does not exist in the data directory.
    """
    resource = files("magestic.projects.nrd1.data").joinpath(filename)
    # Convert to a real Path for downstream code that expects os.PathLike
    path = Path(str(resource))
    if not path.exists():
        raise FileNotFoundError(
            f"NRD1 data file not found: {filename!r}. "
            f"Available files are in magestic.projects.nrd1.data"
        )
    return path


__all__ = ["get_nrd1_data_path"]
