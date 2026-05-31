"""
Common configuration utilities for QTL screen pipelines.

Provides base configuration classes and utilities for:
- Library configuration
- Tier threshold configuration
- Screen date mappings

Author: Kevin Roy Lab
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any
import logging
import yaml
import json

logger = logging.getLogger(__name__)


@dataclass
class LibraryConfig:
    """
    Configuration for a single library within a screen.

    Libraries are typically strain-specific (e.g., yL437, yL442) and
    have different characteristics based on bottlenecking and timepoints.
    """
    name: str
    screen_dates: List[str]
    bottlenecked: bool = False
    timepoints: List[int] = field(default_factory=lambda: [0, 21])
    t0_consolidation_groups: int = 4  # How many consolidated t0 groups
    description: str = ""

    @property
    def is_multi_timepoint(self) -> bool:
        """Check if this library has multiple timepoints."""
        return len(self.timepoints) > 2

    @property
    def screen_date_group(self) -> str:
        """Get the combined screen date group identifier."""
        return "_".join(self.screen_dates)


@dataclass
class TierConfig:
    """
    Configuration for tiered BC1 classification.

    Tiers are based on t0 read counts to stratify BC1s by data quality.
    """
    tier_1_min_reads: int = 100  # Individual BC1 DESeq2
    tier_2_min_reads: int = 20   # Aggregated by oligo
    tier_3_min_reads: int = 1    # Flag only (low-confidence)

    # Hit calling thresholds per tier
    tier_1_z_threshold: float = 3.0
    tier_1_concordance_threshold: float = 0.75
    tier_2_z_threshold: float = 3.75
    tier_2_concordance_threshold: float = 0.85
    tier_2_min_bc1s_aggregated: int = 2

    def get_tier(self, t0_reads: int) -> str:
        """
        Classify a BC1 into a tier based on t0 read count.

        Args:
            t0_reads: Total reads in t0 samples

        Returns:
            Tier string: 'tier_1', 'tier_2', or 'tier_3'
        """
        if t0_reads >= self.tier_1_min_reads:
            return "tier_1"
        elif t0_reads >= self.tier_2_min_reads:
            return "tier_2"
        else:
            return "tier_3"


@dataclass
class ScreenConfig:
    """
    Configuration for a complete QTL screen.

    This is the top-level configuration that combines:
    - Project metadata
    - Library configurations
    - Tier thresholds
    - Path configuration
    """
    project_name: str
    screen_name: str
    libraries: List[LibraryConfig] = field(default_factory=list)
    tier_config: TierConfig = field(default_factory=TierConfig)

    # Sequencing parameters
    sequencing_dates: List[str] = field(default_factory=list)
    bc1_length: int = 20
    bc0_length: int = 10

    # Barcode detection sequences
    fwd_bc1_prefix: str = "CATGCTC"
    fwd_bc1_suffix: str = "CCCTAGG"

    # Paths (set after initialization)
    base_dir: Path = field(default_factory=lambda: Path("/path/to"))

    def __post_init__(self):
        """Initialize computed properties."""
        self._library_mapping: Dict[str, str] = {}
        for lib in self.libraries:
            for date in lib.screen_dates:
                self._library_mapping[date] = lib.name

    @property
    def library_mapping(self) -> Dict[str, str]:
        """Get mapping from screen date to library name."""
        return self._library_mapping

    @property
    def library_names(self) -> List[str]:
        """Get list of library names."""
        return [lib.name for lib in self.libraries]

    def get_library(self, name: str) -> Optional[LibraryConfig]:
        """Get library configuration by name."""
        for lib in self.libraries:
            if lib.name == name:
                return lib
        return None

    def get_library_for_date(self, screen_date: str) -> Optional[LibraryConfig]:
        """Get library configuration for a screen date."""
        lib_name = self._library_mapping.get(screen_date)
        if lib_name:
            return self.get_library(lib_name)
        return None

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> "ScreenConfig":
        """
        Load configuration from a YAML file.

        Args:
            yaml_path: Path to YAML configuration file

        Returns:
            ScreenConfig instance
        """
        with open(yaml_path) as f:
            data = yaml.safe_load(f)

        libraries = [
            LibraryConfig(**lib_data)
            for lib_data in data.get("libraries", [])
        ]

        tier_data = data.get("tier_config", {})
        tier_config = TierConfig(**tier_data)

        return cls(
            project_name=data["project_name"],
            screen_name=data["screen_name"],
            libraries=libraries,
            tier_config=tier_config,
            sequencing_dates=data.get("sequencing_dates", []),
            bc1_length=data.get("bc1_length", 20),
            bc0_length=data.get("bc0_length", 10),
            fwd_bc1_prefix=data.get("fwd_bc1_prefix", "CATGCTC"),
            fwd_bc1_suffix=data.get("fwd_bc1_suffix", "CCCTAGG"),
            base_dir=Path(data.get("base_dir", "/path/to")),
        )

    def to_yaml(self, yaml_path: Path) -> None:
        """
        Save configuration to a YAML file.

        Args:
            yaml_path: Path to save YAML configuration
        """
        data = {
            "project_name": self.project_name,
            "screen_name": self.screen_name,
            "libraries": [
                {
                    "name": lib.name,
                    "screen_dates": lib.screen_dates,
                    "bottlenecked": lib.bottlenecked,
                    "timepoints": lib.timepoints,
                    "t0_consolidation_groups": lib.t0_consolidation_groups,
                    "description": lib.description,
                }
                for lib in self.libraries
            ],
            "tier_config": {
                "tier_1_min_reads": self.tier_config.tier_1_min_reads,
                "tier_2_min_reads": self.tier_config.tier_2_min_reads,
                "tier_3_min_reads": self.tier_config.tier_3_min_reads,
                "tier_1_z_threshold": self.tier_config.tier_1_z_threshold,
                "tier_1_concordance_threshold": self.tier_config.tier_1_concordance_threshold,
                "tier_2_z_threshold": self.tier_config.tier_2_z_threshold,
                "tier_2_concordance_threshold": self.tier_config.tier_2_concordance_threshold,
                "tier_2_min_bc1s_aggregated": self.tier_config.tier_2_min_bc1s_aggregated,
            },
            "sequencing_dates": self.sequencing_dates,
            "bc1_length": self.bc1_length,
            "bc0_length": self.bc0_length,
            "fwd_bc1_prefix": self.fwd_bc1_prefix,
            "fwd_bc1_suffix": self.fwd_bc1_suffix,
            "base_dir": str(self.base_dir),
        }

        with open(yaml_path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)


# Pre-defined configuration for 20250912_SpG_QTL_screen
SPG_QTL_SCREEN_CONFIG = ScreenConfig(
    project_name="QTL",
    screen_name="20250912_SpG_QTL_screen",
    libraries=[
        LibraryConfig(
            name="yL437",
            screen_dates=["20251014", "20251018"],
            bottlenecked=True,
            timepoints=[0, 6, 12, 21],
            t0_consolidation_groups=4,
            description="Multi-timepoint bottlenecked library (~70k BC1s)",
        ),
        LibraryConfig(
            name="yL442",
            screen_dates=["20251212", "20251216"],
            bottlenecked=False,
            timepoints=[0, 21],
            t0_consolidation_groups=8,
            description="December non-bottlenecked library (~500k BC1s)",
        ),
    ],
    sequencing_dates=["20251125", "20260106", "20260107"],
)
