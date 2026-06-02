"""
Common path utility functions for QTL screen pipelines.

Provides standardized path handling for the project structure:
- projects/QTL/{screen_name}/
  ├── processed_data/  (symlinked to legacy location for large data)
  ├── scripts/
  ├── keyfiles/
  ├── plots/
  └── reports/

Author: Kevin Roy Lab
"""

from pathlib import Path
from typing import Optional, Dict, List
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)

# Base directories
DEFAULT_BASE_DIR = Path("/path/to")
DEFAULT_COMMON_DIR = DEFAULT_BASE_DIR / "common"


@dataclass
class ProjectPaths:
    """
    Standardized paths for a QTL screen project.

    Uses the flattened directory structure (no nested scripts/scripts/).
    """
    project_name: str
    base_dir: Path = DEFAULT_BASE_DIR

    def __post_init__(self):
        """Initialize computed paths."""
        self.project_dir = self.base_dir / "projects" / "QTL" / self.project_name

    @property
    def processed_data_dir(self) -> Path:
        """Directory for processed data (symlinked to legacy location)."""
        return self.project_dir / "processed_data"

    @property
    def scripts_dir(self) -> Path:
        """Directory for analysis scripts."""
        return self.project_dir / "scripts"

    @property
    def keyfiles_dir(self) -> Path:
        """Consolidated keyfiles directory."""
        return self.project_dir / "keyfiles"

    @property
    def plots_dir(self) -> Path:
        """Directory for output plots."""
        return self.project_dir / "plots"

    @property
    def reports_dir(self) -> Path:
        """Directory for analysis reports."""
        return self.project_dir / "reports"

    @property
    def logs_dir(self) -> Path:
        """Directory for processing logs."""
        return self.project_dir / "logs"

    def analysis_dir(self, analysis_type: str) -> "AnalysisPaths":
        """
        Get paths for a specific analysis type.

        Args:
            analysis_type: One of 'bc1_from_screen', 'bc1_donor_bc0',
                          'guide_donor_bc0', 'bc1_from_outgrowth'

        Returns:
            AnalysisPaths object with analysis-specific directories
        """
        return AnalysisPaths(
            project_paths=self,
            analysis_type=analysis_type
        )

    def validate(self) -> List[str]:
        """
        Validate that expected directories exist.

        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []
        if not self.project_dir.exists():
            errors.append(f"Project directory does not exist: {self.project_dir}")
        return errors

    def ensure_directories(self) -> None:
        """Create all standard directories if they don't exist."""
        for dir_path in [
            self.processed_data_dir,
            self.scripts_dir,
            self.keyfiles_dir,
            self.plots_dir,
            self.reports_dir,
            self.logs_dir,
        ]:
            dir_path.mkdir(parents=True, exist_ok=True)


@dataclass
class AnalysisPaths:
    """
    Paths for a specific analysis type within a project.

    Analysis types:
    - bc1_from_screen: BC1 barcode counts from screen
    - bc1_donor_bc0: BC1-donor-bc0 associations
    - guide_donor_bc0: Guide-donor-bc0 from plasmid pool
    - bc1_from_outgrowth: BC1 counts from outgrowth cultures
    """
    project_paths: ProjectPaths
    analysis_type: str

    @property
    def processed_data_dir(self) -> Path:
        """Processed data directory for this analysis."""
        return self.project_paths.processed_data_dir / self.analysis_type

    @property
    def scripts_dir(self) -> Path:
        """Scripts directory for this analysis (flattened, no nested scripts/)."""
        return self.project_paths.scripts_dir / self.analysis_type

    @property
    def keyfiles_dir(self) -> Path:
        """Keyfiles directory for this analysis."""
        return self.project_paths.keyfiles_dir / self.analysis_type

    @property
    def plots_dir(self) -> Path:
        """Plots directory for this analysis."""
        # Analysis-specific plots can go in processed_data/analysis/plots/
        # or the project-level plots/ - depends on preference
        return self.processed_data_dir / "plots"

    def get_keyfile(self, filename: str) -> Path:
        """Get path to a specific keyfile."""
        return self.keyfiles_dir / filename

    def get_output_file(self, filename: str) -> Path:
        """Get path for an output file in processed data."""
        return self.processed_data_dir / filename

    def validate(self) -> List[str]:
        """Validate that expected directories exist."""
        errors = []
        if not self.scripts_dir.exists():
            errors.append(f"Scripts directory does not exist: {self.scripts_dir}")
        return errors


def get_common_dir() -> Path:
    """Get the common resources directory."""
    return DEFAULT_COMMON_DIR


def get_annotation_file(filename: str) -> Path:
    """
    Get path to a common annotation file.

    Args:
        filename: Name of the annotation file

    Returns:
        Full path to the annotation file
    """
    return DEFAULT_COMMON_DIR / "annotation_files" / filename


def get_oligo_design_file(filename: str) -> Path:
    """
    Get path to a common oligo design file.

    Args:
        filename: Name of the oligo design file

    Returns:
        Full path to the oligo design file
    """
    return DEFAULT_COMMON_DIR / "oligo_designs" / filename


def get_reference_genome(filename: str) -> Path:
    """
    Get path to a reference genome file.

    Args:
        filename: Name of the reference genome file

    Returns:
        Full path to the reference genome
    """
    return DEFAULT_COMMON_DIR / "reference_genomes" / filename


# Common annotation files
HARMONIZED_VARIANTS_FILE = get_annotation_file(
    "20210422_and_20240411_Bloom_et_al_16_strains_QTL_harmonized_designed_variant_oligos.tsv"
)

# Common oligo design files
SPG_OLIGO_DESIGN = get_oligo_design_file("20240411_Twist_200mer_oligo_array_order.tsv")
# 2021 SpCas9/LbCas12a pool: canonical AS-SUBMITTED (md5 79f34c86, 48000 rows).
# Repointed 2026-06-02 off 20210422_Twist_200mer.tsv (md5 048501cc, 48012 rows,
# +12 MIP1 over-rep pre-trim snapshot). See common/oligo_designs/OLIGO_POOLS_PROVENANCE.md.
SPCAS9_OLIGO_DESIGN = get_oligo_design_file("20210422_Twist_200mer_as_submitted.tsv")

# Reference genome
MAGESTIC_REFERENCE_GENOME = get_reference_genome("MAGESTIC_background_strain.fasta")
MAGESTIC_ANNOTATION_GFF = get_annotation_file("MAGESTIC_background_strain_annotations.gff")
