"""
Configuration management for the BC1-Donor-BC0 pipeline.

Provides dataclasses for pipeline configuration and automatic
detection of available data sources.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional
import logging

logger = logging.getLogger(__name__)


# Sequence constants
BC1_LENGTH = 20
BC0_LENGTH = 10
DONOR_BC0_FRAGMENT_LENGTH = 60  # ~30bp donor + ~30bp common sequence
TYPICAL_DONOR_LENGTH = 129
MIDDLE_DONOR_LENGTH = 60  # Window size for sliding match

# Oligo structure constants (200mer design)
GUIDE_START = 20
GUIDE_LENGTH = 20
SPCAS9_CLONING_SITE = "GTTTGAAGAGC"
DONOR_START = GUIDE_START + GUIDE_LENGTH + len(SPCAS9_CLONING_SITE)  # Position 51

# Common sequence in donor_bc0_fragment (appears after donor portion)
COMMON_SEQUENCE_START = "CGAT"  # Start of common region


@dataclass
class OligoDesignConfig:
    """Configuration for an oligo design file."""
    name: str
    file_path: Path
    sublibraries: List[str]  # Which sublibraries use this design


@dataclass
class PipelineConfig:
    """
    Configuration for the bc1-donor-bc0 pipeline.

    Attributes:
        project_name: Name of the project/screen
        project_dir: Root directory for the project

        # Data availability flags
        has_2x150_data: Whether 2x150 data is available
        has_2x300_data: Whether 2x300 data is available
        has_gdb_data: Whether guide-donor-bc0 data is available

        # Input paths (set based on availability)
        data_2x150_dir: Directory containing 2x150 processed data
        data_2x300_dir: Directory containing 2x300 processed data
        gdb_dir: Directory containing guide-donor-bc0 data

        # Oligo designs
        oligo_designs: List of oligo design configurations

        # Output
        output_dir: Directory for pipeline outputs

        # Processing parameters
        min_reads: Minimum reads for a bc1-donor-bc0 to be included
        min_pcr_replicates: Minimum PCR replicates for confidence
        bc1_purity_threshold: Minimum purity for bc1 assignment
    """
    project_name: str
    project_dir: Path

    # Data availability
    has_2x150_data: bool = False
    has_2x300_data: bool = False
    has_gdb_data: bool = False

    # Input paths
    data_2x150_dir: Optional[Path] = None
    data_2x300_dir: Optional[Path] = None
    gdb_dir: Optional[Path] = None

    # Oligo designs
    oligo_designs: List[OligoDesignConfig] = field(default_factory=list)

    # Output
    output_dir: Optional[Path] = None

    # Processing parameters
    min_reads: int = 2
    min_pcr_replicates: int = 1
    bc1_purity_threshold: float = 0.8

    def __post_init__(self):
        """Validate configuration after initialization."""
        if not any([self.has_2x150_data, self.has_2x300_data]):
            logger.warning("No data sources specified. Pipeline may not produce results.")

        if self.output_dir is None:
            self.output_dir = self.project_dir / "processed_data" / "bc1_donor_bc0"

    @classmethod
    def from_project_dir(cls, project_dir: Path,
                         common_dir: Optional[Path] = None) -> "PipelineConfig":
        """
        Create configuration by auto-detecting available data sources.

        Args:
            project_dir: Root directory of the project
            common_dir: Directory containing common resources (oligo designs, etc.)

        Returns:
            PipelineConfig with detected settings
        """
        project_dir = Path(project_dir)

        if common_dir is None:
            # Default to standard location
            common_dir = Path("/path/to/common")

        config = cls(
            project_name=project_dir.name,
            project_dir=project_dir,
        )

        # Detect 2x150 data
        possible_2x150_paths = [
            project_dir / "processed_data" / "bc1_donor_bc0" / "2x150",
            project_dir / "scripts" / "bc1_donor_bc0" / "scripts" / "20251125_2x150_AVITI_run",
        ]
        for path in possible_2x150_paths:
            if path.exists():
                config.has_2x150_data = True
                config.data_2x150_dir = path
                logger.info(f"Found 2x150 data at {path}")
                break

        # Detect 2x300 data
        possible_2x300_paths = [
            project_dir / "processed_data" / "bc1_donor_bc0" / "2x300",
            project_dir / "scripts" / "bc1_donor_bc0" / "scripts" / "20251212_2x300_AVITI_run",
        ]
        for path in possible_2x300_paths:
            if path.exists():
                config.has_2x300_data = True
                config.data_2x300_dir = path
                logger.info(f"Found 2x300 data at {path}")
                break

        # Detect gdb data
        possible_gdb_paths = [
            project_dir / "processed_data" / "guide_donor_bc0",
            project_dir / "scripts" / "guide_donor_bc0",
        ]
        for path in possible_gdb_paths:
            if path.exists():
                config.has_gdb_data = True
                config.gdb_dir = path
                logger.info(f"Found guide-donor-bc0 data at {path}")
                break

        # Set up default oligo designs
        spg_design = common_dir / "oligo_designs" / "20240411_Twist_200mer_oligo_array_order.tsv"
        spcas9_design = common_dir / "oligo_designs" / "20210422_Twist_200mer_SpCas9_oligos_with_controls.tsv"

        if spg_design.exists():
            config.oligo_designs.append(OligoDesignConfig(
                name="SpG_2024",
                file_path=spg_design,
                sublibraries=["SpG_NGNG", "SpG_NGNH"]
            ))

        if spcas9_design.exists():
            config.oligo_designs.append(OligoDesignConfig(
                name="SpCas9_2021",
                file_path=spcas9_design,
                sublibraries=["SpG_NGG", "SpCas9_NGG"]
            ))

        return config

    def validate(self) -> List[str]:
        """
        Validate the configuration.

        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []

        if not self.project_dir.exists():
            errors.append(f"Project directory does not exist: {self.project_dir}")

        if self.has_2x150_data and self.data_2x150_dir and not self.data_2x150_dir.exists():
            errors.append(f"2x150 data directory does not exist: {self.data_2x150_dir}")

        if self.has_2x300_data and self.data_2x300_dir and not self.data_2x300_dir.exists():
            errors.append(f"2x300 data directory does not exist: {self.data_2x300_dir}")

        if self.has_gdb_data and self.gdb_dir and not self.gdb_dir.exists():
            errors.append(f"Guide-donor-bc0 directory does not exist: {self.gdb_dir}")

        for design in self.oligo_designs:
            if not design.file_path.exists():
                errors.append(f"Oligo design file does not exist: {design.file_path}")

        return errors
