"""
magestic.pipelines.redi.core.config — YAML config schema for REDI pipeline runs.

Defines the per-screen YAML schema (`REDI/configs/<screen>.yaml`) parsed into a
typed `REDIConfig` dataclass. Mirrors the constants hard-coded in the legacy
`step_a`–`step_i` scripts so a single config file replaces the
`DIR`/`SCREEN_NAME`/`PROJECT_SUBPATH` triplet that each script previously
duplicated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class REDIPaths:
    """Filesystem layout for one REDI screen."""

    project_root: Path  # /oak/.../projects/REDI/<screen>/
    raw_fastq: Path  # raw_data/fastq/
    intermediate: Path  # intermediate_data/
    processed: Path  # processed_data/
    keyfiles: Path  # keyfiles/

    @property
    def trimmed_fastq(self) -> Path:
        return self.intermediate / "trimmed_fastq"

    @property
    def merged_fastq(self) -> Path:
        return self.intermediate / "merged_fastq"

    @property
    def collapsed_fastq(self) -> Path:
        return self.processed / "collapsed_fastq"

    @property
    def counts_tables(self) -> Path:
        return self.processed / "read_pair_counts_tables"

    @property
    def combined_data(self) -> Path:
        return self.processed / "combined_data_tables"

    @property
    def plots(self) -> Path:
        return self.processed / "plots"


@dataclass
class REDIBarcodeSchema:
    """Sequence anchors used to locate bc1 and REDI_bc within reads.

    Defaults from `20250312_NNS_complex_satmut/scripts/step_c`. Per-screen
    overrides allowed (older primer sets KR1882/KR1883 differ).
    """

    bc1_length: int = 20
    fwd_bc1_prefix: str = "CAGGTCATGCTC"
    fwd_bc1_suffix: str = "CCCTAGGATACG"
    # When the prefix can't be found, fall back to a fixed offset that works
    # for KR1965/KR1967 primer sets (20-bp trim + 'TCGACAGGTCATGCTC').
    fallback_prefix_length: int = 16  # len('TCGACAGGTCATGCTC')
    # Synthesis/sequencing-error tolerance on bc1 length (matches bc1_from_screen).
    bc1_length_tolerance: int = 1


@dataclass
class REDIPurityConfig:
    """Purity calculation parameters.

    `top_n_for_purity`: number of bc1's per (sample, REDI_bc) summed for the
    purity denominator. Legacy step_d uses 4.
    `edit_distance_drop`: bc1 values whose Levenshtein distance to the top bc1
    is in this set are dropped (sequencing-error daughters). Legacy: {1, 2}.
    `min_counts`: optional read-count floor on `REDI_bc_bc1_counts`.
    """

    top_n_for_purity: int = 4
    edit_distance_drop: tuple[int, ...] = (1, 2)
    min_counts: int = 10


@dataclass
class REDIArrayPair:
    """One cross-array comparison (e.g., SHA345 vs yT182)."""

    name_a: str
    name_b: str
    purpose_a: str  # value of `purpose` column selecting array A rows
    purpose_b_not: str = ""  # if non-empty, B is everything where purpose != purpose_a


@dataclass
class REDIConfig:
    """Top-level config for one REDI screen run."""

    screen_name: str
    sequencing_date: str
    paths: REDIPaths
    barcode: REDIBarcodeSchema = field(default_factory=REDIBarcodeSchema)
    purity: REDIPurityConfig = field(default_factory=REDIPurityConfig)
    # REDI_bc key files; resolved relative to paths.keyfiles unless absolute
    redi_bc_keys: list[str] = field(default_factory=list)
    # Combined keyfile dependencies (legacy step_b inputs)
    sample_key: str = "sample_key.tsv"
    gdna_plate_key: str = "gDNA_plate_key.tsv"
    redi_plate_barcode_key: str = "REDI_plate_barcode_key.tsv"
    # Cross-array comparisons
    array_pairs: list[REDIArrayPair] = field(default_factory=list)
    # bc1-donor-bc0 reference table (output of bc1_donor_bc0 pipeline)
    bc1_reference_table: Path | None = None
    # WGS ground-truth table (Shengdi's output)
    wgs_outcome_table: Path | None = None
    # Inner-primer rules: which side of the read pair holds bc1 / REDI_bc
    primer_rules: dict[str, str] = field(default_factory=lambda: {
        "KR1965": "R2",  # KR1965 → bc1 on R2, REDI_bc on R1
        "KR1967": "R1",  # KR1967 → bc1 on R1, REDI_bc on R2
    })


def load_config(path: str | Path) -> REDIConfig:
    """Parse a YAML config into a `REDIConfig`."""
    path = Path(path)
    with path.open() as fh:
        data: dict[str, Any] = yaml.safe_load(fh)

    paths_d = data["paths"]
    paths = REDIPaths(
        project_root=Path(paths_d["project_root"]),
        raw_fastq=Path(paths_d["raw_fastq"]),
        intermediate=Path(paths_d["intermediate"]),
        processed=Path(paths_d["processed"]),
        keyfiles=Path(paths_d["keyfiles"]),
    )

    barcode = REDIBarcodeSchema(**data.get("barcode", {}))
    purity = REDIPurityConfig(**data.get("purity", {}))
    array_pairs = [REDIArrayPair(**ap) for ap in data.get("array_pairs", [])]

    return REDIConfig(
        screen_name=data["screen_name"],
        sequencing_date=str(data["sequencing_date"]),
        paths=paths,
        barcode=barcode,
        purity=purity,
        redi_bc_keys=data.get("redi_bc_keys", []),
        sample_key=data.get("sample_key", "sample_key.tsv"),
        gdna_plate_key=data.get("gdna_plate_key", "gDNA_plate_key.tsv"),
        redi_plate_barcode_key=data.get(
            "redi_plate_barcode_key", "REDI_plate_barcode_key.tsv"
        ),
        array_pairs=array_pairs,
        bc1_reference_table=Path(data["bc1_reference_table"])
        if data.get("bc1_reference_table")
        else None,
        wgs_outcome_table=Path(data["wgs_outcome_table"])
        if data.get("wgs_outcome_table")
        else None,
        primer_rules=data.get("primer_rules", {"KR1965": "R2", "KR1967": "R1"}),
    )
