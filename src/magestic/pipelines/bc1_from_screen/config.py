"""
Configuration management for the BC1-From-Screen pipeline.

Provides dataclasses for pipeline configuration including:
- Library-specific settings
- Tier thresholds
- Hit calling parameters
- T0 consolidation settings
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
import os
import yaml

logger = logging.getLogger(__name__)

# Scratch filesystem path (uses $SCRATCH env variable if set, else /tmp)
SCRATCH = Path(os.environ.get("SCRATCH", "/tmp"))


def get_scratch_workdir(project_name: str) -> Path:
    """Return a per-project scratch working directory under $SCRATCH."""
    return SCRATCH / "magestic" / project_name


# Barcode constants
BC1_LENGTH = 20
FWD_BC1_PREFIX = "CATGCTC"
FWD_BC1_SUFFIX = "CCCTAGG"

# Primer pair patterns (forward and reverse complement orientations)
PRIMER_PAIRS = {
    "KR2470_KR2473": {
        "pattern": "TCATGCTCNNNNNNNNNNNNNNNNNNNNCCCTAggatAcg",
        "orientation": "forward",
    },
    "KR1966_KR1967": {
        "pattern": "TCGACAGGTCATGCTCNNNNNNNNNNNNNNNNNNNNCCCTAggatACGTGCATAAC",
        "orientation": "forward",
    },
}


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
    t0_consolidation_groups: int = 4
    description: str = ""

    @property
    def is_multi_timepoint(self) -> bool:
        """Check if this library has multiple timepoints."""
        return len(self.timepoints) > 2

    @property
    def screen_date_group(self) -> str:
        """Get the combined screen date group identifier."""
        return "_".join(self.screen_dates)

    @property
    def generations(self) -> List[int]:
        """Get list of generation timepoints."""
        return self.timepoints


@dataclass
class TierConfig:
    """
    Configuration for tiered BC1 classification based on read abundance.

    Tiers are based on t0 read counts to stratify BC1s by data quality:
    - abundance_high (≥100 reads): Individual BC1 analysis in DESeq2
    - abundance_medium (20-99 reads): Aggregate by oligo before DESeq2
    - abundance_low (<20 reads): Excluded from analysis (low confidence)

    NOTE: These "abundance tiers" are distinct from the "confidence_tier" column
    in bc1_donor_bc0_pipeline which classifies oligo attribution confidence.
    """
    # Read count thresholds (attribute names kept for YAML compatibility)
    tier_1_min_reads: int = 100  # abundance_high threshold
    tier_2_min_reads: int = 20   # abundance_medium threshold
    tier_3_min_reads: int = 1    # abundance_low threshold

    def get_tier(self, t0_reads: int) -> str:
        """
        Classify a BC1 into an abundance tier based on t0 read count.

        Args:
            t0_reads: Total reads in t0 samples

        Returns:
            Tier string: 'abundance_high', 'abundance_medium', or 'abundance_low'
        """
        if t0_reads >= self.tier_1_min_reads:
            return "abundance_high"
        elif t0_reads >= self.tier_2_min_reads:
            return "abundance_medium"
        else:
            return "abundance_low"


@dataclass
class HitCallingConfig:
    """
    Configuration for abundance-tier-aware hit calling.

    Different thresholds are applied to:
    - abundance_high (individual BC1s): More lenient thresholds
    - abundance_medium (oligo-aggregated pseudo-BC1s): Stricter thresholds

    NOTE: Attribute names kept as tier_1/tier_2 for YAML config compatibility.
    """
    # abundance_high thresholds (was tier_1)
    tier_1_z_threshold: float = 3.0
    tier_1_concordance_threshold: float = 0.75

    # abundance_medium thresholds (was tier_2) - stricter
    tier_2_z_threshold: float = 3.75
    tier_2_concordance_threshold: float = 0.85
    tier_2_min_bc1s_aggregated: int = 2

    # Adaptive settings
    use_adaptive_base_mean: bool = True
    base_mean_percentile: float = 5.0

    use_adaptive_min_bc1: bool = True
    strong_z_for_relaxed_bc1: float = 5.0
    low_cv_for_relaxed_bc1: float = 0.3
    high_concordance_for_relaxed_bc1: float = 0.9

    def get_thresholds(self, tier: str) -> Dict:
        """Get hit calling thresholds for a specific abundance tier."""
        # Accept both old and new tier names for backwards compatibility
        if tier in ("tier_1", "abundance_high"):
            return {
                "z_threshold": self.tier_1_z_threshold,
                "concordance_threshold": self.tier_1_concordance_threshold,
            }
        elif tier in ("tier_2", "abundance_medium"):
            return {
                "z_threshold": self.tier_2_z_threshold,
                "concordance_threshold": self.tier_2_concordance_threshold,
                "min_bc1s_aggregated": self.tier_2_min_bc1s_aggregated,
            }
        else:
            raise ValueError(f"Unknown tier: {tier}. Expected 'abundance_high', 'abundance_medium', 'tier_1', or 'tier_2'")


@dataclass
class BC1FromScreenConfig:
    """
    Configuration for the BC1-From-Screen pipeline.

    This is the top-level configuration that combines:
    - Project metadata
    - Library configurations
    - Tier thresholds
    - Hit calling parameters
    """
    project_name: str
    screen_name: str
    project_dir: Path

    # Library configurations
    libraries: List[LibraryConfig] = field(default_factory=list)

    # Processing configurations
    tier_config: TierConfig = field(default_factory=TierConfig)
    hit_calling_config: HitCallingConfig = field(default_factory=HitCallingConfig)

    # Sequencing parameters
    sequencing_dates: List[str] = field(default_factory=list)

    # Path configuration
    base_dir: Path = field(
        default_factory=lambda: Path("/path/to")
    )

    # Quality thresholds
    min_total_reads_per_sample: int = 50000
    min_replicates_per_condition: int = 2

    # Scratch filesystem support (for SLURM jobs)
    use_scratch: bool = False
    scratch_dir: Optional[Path] = None

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

    @property
    def processed_data_dir(self) -> Path:
        """Get processed data directory."""
        return self.project_dir / "processed_data" / "bc1_from_screen"

    @property
    def bc1_counts_dir(self) -> Path:
        """Get BC1 counts directory."""
        return self.processed_data_dir / "bc1_counts"

    @property
    def deseq2_dir(self) -> Path:
        """Get DESeq2 results directory."""
        return self.processed_data_dir / "DESeq2_consolidated_t0"

    @property
    def hit_calling_dir(self) -> Path:
        """Get hit calling results directory."""
        return self.processed_data_dir / "hit_calling_tiered"

    @property
    def keyfiles_dir(self) -> Path:
        """Get keyfiles directory."""
        return self.project_dir / "keyfiles" / "bc1_from_screen"

    def get_working_dir(self, subdir: str) -> Path:
        """
        Return the working directory for a pipeline sub-step.

        When use_scratch is True and scratch_dir is set, returns a path under
        scratch for fast I/O.  Otherwise returns the standard processed_data_dir.
        """
        if self.use_scratch and self.scratch_dir is not None:
            return self.scratch_dir / subdir
        return self.processed_data_dir / subdir

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> "BC1FromScreenConfig":
        """
        Load configuration from a YAML file.

        Args:
            yaml_path: Path to YAML configuration file

        Returns:
            BC1FromScreenConfig instance
        """
        with open(yaml_path) as f:
            data = yaml.safe_load(f)

        libraries = [
            LibraryConfig(**lib_data)
            for lib_data in data.get("libraries", [])
        ]

        tier_data = data.get("tier_config", {})
        tier_config = TierConfig(**tier_data)

        hit_calling_data = data.get("hit_calling_config", {})
        hit_calling_config = HitCallingConfig(**hit_calling_data)

        return cls(
            project_name=data["project_name"],
            screen_name=data["screen_name"],
            project_dir=Path(data["project_dir"]),
            libraries=libraries,
            tier_config=tier_config,
            hit_calling_config=hit_calling_config,
            sequencing_dates=data.get("sequencing_dates", []),
            base_dir=Path(data.get("base_dir", "/path/to")),
            min_total_reads_per_sample=data.get("min_total_reads_per_sample", 50000),
            min_replicates_per_condition=data.get("min_replicates_per_condition", 2),
        )

    def to_yaml(self, yaml_path: Path) -> None:
        """Save configuration to a YAML file."""
        data = {
            "project_name": self.project_name,
            "screen_name": self.screen_name,
            "project_dir": str(self.project_dir),
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
            },
            "hit_calling_config": {
                "tier_1_z_threshold": self.hit_calling_config.tier_1_z_threshold,
                "tier_1_concordance_threshold": self.hit_calling_config.tier_1_concordance_threshold,
                "tier_2_z_threshold": self.hit_calling_config.tier_2_z_threshold,
                "tier_2_concordance_threshold": self.hit_calling_config.tier_2_concordance_threshold,
                "tier_2_min_bc1s_aggregated": self.hit_calling_config.tier_2_min_bc1s_aggregated,
            },
            "sequencing_dates": self.sequencing_dates,
            "base_dir": str(self.base_dir),
            "min_total_reads_per_sample": self.min_total_reads_per_sample,
            "min_replicates_per_condition": self.min_replicates_per_condition,
        }

        with open(yaml_path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def validate(self) -> List[str]:
        """Validate the configuration."""
        errors = []

        if not self.project_dir.exists():
            errors.append(f"Project directory does not exist: {self.project_dir}")

        if not self.libraries:
            errors.append("No libraries configured")

        # Check for overlapping screen dates
        all_dates = []
        for lib in self.libraries:
            for date in lib.screen_dates:
                if date in all_dates:
                    errors.append(f"Screen date {date} appears in multiple libraries")
                all_dates.append(date)

        return errors


# Pre-defined configuration for 20250912_SpG_QTL_screen
def create_spg_qtl_screen_config(
    base_dir: Path = Path("/path/to")
) -> BC1FromScreenConfig:
    """Create configuration for the 20250912_SpG_QTL_screen project."""
    return BC1FromScreenConfig(
        project_name="QTL",
        screen_name="20250912_SpG_QTL_screen",
        project_dir=base_dir / "projects" / "QTL" / "20250912_SpG_QTL_screen",
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
        base_dir=base_dir,
    )
