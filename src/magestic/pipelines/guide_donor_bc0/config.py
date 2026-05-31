"""
Configuration management for the Guide-Donor-BC0 pipeline.

Provides dataclasses for pipeline configuration and constants
for oligo matching.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


# Oligo structure constants (200mer design)
GUIDE_START = 20
GUIDE_END = 40
GUIDE_LENGTH = GUIDE_END - GUIDE_START  # 20bp

DONOR_START = 51
DONOR_END = 180
DONOR_LENGTH = DONOR_END - DONOR_START  # 129bp

# Middle donor region for partial matching
MIDDLE_DONOR_START = 30  # Relative to donor start
MIDDLE_DONOR_END = 90    # Relative to donor start
MIDDLE_DONOR_LENGTH = MIDDLE_DONOR_END - MIDDLE_DONOR_START  # 60bp

# Sequence constants
SPCAS9_CLONING_SITE = "GTTTGAAGAGC"
PERFECT_LEADER = "GCAGGCGCGCC"

# Sherlock scratch filesystem
SCRATCH = Path(os.environ.get("SCRATCH", "/tmp"))


def get_scratch_workdir(project_name: str) -> Path:
    """
    Return a per-project working directory under $SCRATCH.

    Args:
        project_name: Project identifier (used as subdirectory name)

    Returns:
        Path under $SCRATCH suitable for job I/O
    """
    return SCRATCH / "magestic_gdb" / project_name


@dataclass
class OligoDesignConfig:
    """Configuration for an oligo design file."""
    name: str
    file_path: Path
    sublibraries: List[str] = field(default_factory=list)

    def validate(self) -> List[str]:
        """Validate the configuration."""
        errors = []
        if not self.file_path.exists():
            errors.append(f"Oligo design file does not exist: {self.file_path}")
        return errors


@dataclass
class GuideDonorBC0Config:
    """
    Configuration for the guide-donor-bc0 pipeline.

    Attributes:
        project_name: Name of the project/screen
        project_dir: Root directory for the project

        # Input paths
        collapsed_fastq_dir: Directory containing collapsed FASTQ files
        oligo_design_file: Path to oligo design TSV

        # Output paths
        output_dir: Directory for pipeline outputs
        counts_dir: Directory for count tables
        plots_dir: Directory for output plots

        # Processing parameters
        min_reads: Minimum reads for a sequence to be counted
        bc0_length: Length of BC0 barcode (default 10)

        # Scratch I/O (for SLURM jobs)
        use_scratch: Stage I/O through $SCRATCH for better performance
        scratch_dir: Path to scratch working directory
    """
    project_name: str
    project_dir: Path

    # Input paths
    collapsed_fastq_dir: Optional[Path] = None
    oligo_design_file: Optional[Path] = None

    # Output paths
    output_dir: Optional[Path] = None
    counts_dir: Optional[Path] = None
    plots_dir: Optional[Path] = None

    # Processing parameters
    min_reads: int = 1
    bc0_length: int = 10

    # Oligo structure (can override defaults)
    guide_slice: tuple = (GUIDE_START, GUIDE_END)
    donor_slice: tuple = (DONOR_START, DONOR_END)
    middle_donor_start: int = MIDDLE_DONOR_START
    middle_donor_end: int = MIDDLE_DONOR_END

    # Step 02: Sample keyfile path
    sample_keyfile: Optional[Path] = None
    """Path to step_1_sample_key_unified.tsv (auto-detected if None)."""

    # Steps 05–06: Library-to-SPS mapping for missing-oligo filling
    library_sps_map: Dict[str, str] = field(default_factory=dict)
    """Maps library ID → expected SPS sequence (e.g. {"V629": "CTTGAC..."}).
    Leave empty to skip missing-oligo filling."""

    # Scratch I/O
    use_scratch: bool = False
    scratch_dir: Optional[Path] = None

    def __post_init__(self):
        """Set default paths after initialization."""
        if self.output_dir is None:
            self.output_dir = self.project_dir / "processed_data" / "guide_donor_bc0"

        if self.counts_dir is None:
            self.counts_dir = self.output_dir / "guide_donor_bc0"

        if self.plots_dir is None:
            self.plots_dir = self.output_dir / "plots"

    def get_working_dir(self, subdir: str) -> Path:
        """
        Return the working directory for a given output subdirectory.

        When use_scratch is True, returns a path under scratch_dir.
        Otherwise, returns a path under output_dir.

        Args:
            subdir: Subdirectory name (e.g., 'guide_donor_bc0', 'plots')

        Returns:
            Path to working directory
        """
        if self.use_scratch and self.scratch_dir is not None:
            return self.scratch_dir / subdir
        return self.output_dir / subdir

    @classmethod
    def from_project_dir(
        cls,
        project_dir: Path,
        common_dir: Optional[Path] = None,
    ) -> "GuideDonorBC0Config":
        """
        Create configuration by auto-detecting available data sources.

        Args:
            project_dir: Root directory of the project
            common_dir: Directory containing common resources

        Returns:
            GuideDonorBC0Config with detected settings
        """
        project_dir = Path(project_dir)

        if common_dir is None:
            common_dir = Path("/path/to/common")

        config = cls(
            project_name=project_dir.name,
            project_dir=project_dir,
        )

        # Detect collapsed FASTQ directory
        possible_fastq_paths = [
            project_dir / "processed_data" / "guide_donor_bc0" / "collapsed_fastq",
            project_dir / "scripts" / "guide_donor_bc0" / "collapsed_fastq",
        ]
        for path in possible_fastq_paths:
            if path.exists():
                config.collapsed_fastq_dir = path
                logger.info(f"Found collapsed FASTQ directory at {path}")
                break

        # Detect sample keyfile
        possible_keyfiles = [
            project_dir / "keyfiles" / "guide_donor_bc0" / "step_1_sample_key_unified.tsv",
            project_dir / "keyfiles" / "step_1_sample_key_unified.tsv",
        ]
        for path in possible_keyfiles:
            if path.exists():
                config.sample_keyfile = path
                logger.info(f"Found sample keyfile at {path}")
                break

        # Detect oligo design file
        possible_design_files = [
            project_dir / "keyfiles" / "guide_donor_bc0" / "20240411_Twist_200mer_oligo_array_order.tsv",
            common_dir / "oligo_designs" / "20240411_Twist_200mer_oligo_array_order.tsv",
        ]
        for path in possible_design_files:
            if path.exists():
                config.oligo_design_file = path
                logger.info(f"Found oligo design file at {path}")
                break

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

        if self.collapsed_fastq_dir and not self.collapsed_fastq_dir.exists():
            errors.append(f"Collapsed FASTQ directory does not exist: {self.collapsed_fastq_dir}")

        if self.oligo_design_file and not self.oligo_design_file.exists():
            errors.append(f"Oligo design file does not exist: {self.oligo_design_file}")

        return errors

    def ensure_directories(self) -> None:
        """Create output directories if they don't exist."""
        for dir_path in [self.output_dir, self.counts_dir, self.plots_dir]:
            if dir_path:
                dir_path.mkdir(parents=True, exist_ok=True)
