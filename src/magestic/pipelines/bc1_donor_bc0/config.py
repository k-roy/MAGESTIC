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
    # P4 FIX (2026-07-28): was 1 here and 2 in run_pipeline.py's argparse, so the
    # effective threshold depended on the call path and the recombination filter
    # had probably never actually fired (>=1 passes everything). Aligned to the
    # documented CLI default of 2. NOTE this ACTIVATES a filter that was
    # previously a no-op -- see planning/QTL_analyses/BC1_LOADER_FIX_20260728.md.
    # Achievable ceiling (corrected 2026-07-28 by Kevin; my first read was wrong):
    # a bc1 can reach **4** independent PCR replicates, not 2. The bc1-donor-bc0
    # amplicons exist at TWO levels and the same bc1 appears in both:
    #   trafo-level amplicon (e.g. trafo_2-1)  x {plate_3F/Q5, plate_1A/KAPA} = 2
    #   pooled yL-level amplicon (e.g. yL387)  x {plate_3F/Q5, plate_1A/KAPA} = 2
    # For yL406: 45 distinct trafo-level + 48 pooled yL-level sample_numbers;
    # trafo_key maps yL386->trafo_1, yL387->trafo_2, ... Both are legitimate PCR
    # replicates for bc1->donor MATCHING (the bc1-donor linkage is fixed at
    # barcoding and does not change with outgrowth). They carry distinct
    # sample_numbers, so pcr_replicate_id already separates them correctly.
    # => a threshold of 3 IS achievable and meaningful. Measure the realized
    # distribution before choosing between 2 and 3.
    #
    # RATIFIED 2026-09-09 (Kevin, direct): N=2 -- **AS A LABEL, NOT A FILTER.**
    # "It is also useful to see what fraction of bc1 ... were only captured by a
    # single replicate. Those can be provisionally assigned."
    # Measured on yL406 (replicate_dist_20260729/03_postcollapse_minreads1.json):
    # 10,067,553 bc1-donor links; 9,208,614 (91.47%) have exactly 1 replicate;
    # 858,939 (8.53%) have >=2. Filtering at 2 therefore DISCARDS 9.2M links.
    # Required semantics: every link is RETAINED and labelled
    #   support_tier = 'confident'   if n_pcr_replicates >= 2
    #   support_tier = 'provisional' if n_pcr_replicates == 1
    # 🔴 `core/merge.py:630-632` DROPS rows when this value is >1. That function
    # currently has ZERO call sites (dead code, signature default 1), so nothing
    # is filtered today -- but wiring it up would silently delete 91.47% of the
    # data. When it IS wired, it must either gain mode="label" (writes
    # support_tier, returns every row; mode="filter" explicit and opt-in) or be
    # deleted. Do not leave a dead filter armed at 2.
    # Spec: planning/pipeline_provenance_rebuild/21_support_tier_requirement.md
    min_pcr_replicates: int = 2
    bc1_purity_threshold: float = 0.8

    def __post_init__(self):
        """Set derived defaults after initialization.

        FIX 2026-08-06: this used to emit
            "No data sources specified. Pipeline may not produce results."
        here. That is ALWAYS wrong inside ``from_project_dir``, which constructs
        the object FIRST and only detects data sources afterwards -- so the
        warning fired on every auto-detected run no matter what was on disk.
        A 2026-07-25 audit read that warning as proof the data-source detector
        had failed and attributed a partial yL406 load to it
        (AUDIT_20260725_SGTC_state.md:413-417); the real cause was per-file
        pattern resolution, fixed separately on 2026-07-28. Measured 2026-08-06
        with the current tree: yL406 detects 2x300 and loads 288/288 files while
        still printing this warning. The check now lives in
        ``warn_if_no_data_sources()`` and is called AFTER detection.
        """
        if self.output_dir is None:
            self.output_dir = self.project_dir / "processed_data" / "bc1_donor_bc0"

    def warn_if_no_data_sources(self) -> bool:
        """Warn when neither read layout was found. Returns True if none were.

        Call this AFTER data-source detection, never from ``__post_init__``.
        """
        if not any([self.has_2x150_data, self.has_2x300_data]):
            logger.warning(
                "No 2x150 or 2x300 data source found under "
                f"{self.project_dir / 'processed_data' / 'bc1_donor_bc0'} -- the pipeline "
                "will not produce results. Point --project-dir at a tree containing a "
                "2x150/ or 2x300/ directory, or pass the data directory explicitly."
            )
            return True
        return False

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

        # --- Oligo designs -----------------------------------------------------
        # FIX 2026-07-29: this looked for
        #   oligo_designs/20210422_Twist_200mer_SpCas9_oligos_with_controls.tsv
        # which DOES NOT EXIST on Oak -- the real file is `20210422_Twist_200mer.tsv`
        # (and run_pipeline.py's --spcas9-oligo-design help already documented that
        # correct name, so the two disagreed). The `if .exists()` had no `else`, so
        # the SpCas9 design was **silently dropped** and every SpCas9/NGG oligo went
        # unmatched. yL406 contains cV416 NGG sublibraries (yL391/yL397, trafo_6,
        # yKR1036 WT SpCas9), so this was actively wrong for that screen.
        # Same failure class as the loader defect: a path that quietly isn't there.
        # Now: try known filenames in order, and WARN LOUDLY if none resolve.
        design_specs = [
            ("SpG_2024", ["SpG_NGNG", "SpG_NGNH"], [
                "20240411_Twist_200mer_oligo_array_order.tsv",
            ]),
            ("SpCas9_2021", ["SpG_NGG", "SpCas9_NGG"], [
                "20210422_Twist_200mer.tsv",                          # actual file on Oak
                "20210422_Twist_200mer_SpCas9_oligos_with_controls.tsv",  # legacy name
                "20210422_Twist_200mer_as_submitted.tsv",             # as-submitted variant
            ]),
        ]
        for name, sublibs, candidates in design_specs:
            for fname in candidates:
                path = common_dir / "oligo_designs" / fname
                if path.exists():
                    config.oligo_designs.append(OligoDesignConfig(
                        name=name, file_path=path, sublibraries=sublibs))
                    logger.info(f"Found {name} oligo design: {path}")
                    break
            else:
                logger.warning(
                    f"NO oligo design found for {name} (sublibraries {sublibs}). Tried "
                    f"{candidates} under {common_dir / 'oligo_designs'}. Oligos from these "
                    f"sublibraries will NOT be matched -- pass --{name.split('_')[0].lower()}"
                    f"-oligo-design explicitly, or fix --common-dir."
                )

        # Now that detection has actually run, say something true about it
        # (see warn_if_no_data_sources / __post_init__ for why this moved here).
        config.warn_if_no_data_sources()

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
