"""
magestic.pipelines.redi.core.config — YAML config schema for REDI pipeline runs.

Defines the per-screen YAML schema (`REDI/configs/<screen>.yaml`) parsed into a
typed `REDIConfig` dataclass. Mirrors the constants hard-coded in the legacy
`step_a`–`step_i` scripts so a single config file replaces the
`DIR`/`SCREEN_NAME`/`PROJECT_SUBPATH` triplet that each script previously
duplicated.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
import math
from pathlib import Path
import re
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
    fwd_bc1_suffix: str | None = "CCCTAGGATACG"
    # Legacy positional extraction is opt-in and does not establish cassette identity.
    fallback_prefix_length: int | None = 16
    allow_positional_fallback: bool = False
    # Synthesis/sequencing-error tolerance on bc1 length (matches bc1_from_screen).
    bc1_length_tolerance: int = 1


@dataclass
class REDIPurityConfig:
    """Purity calculation parameters.

    Primary purity is the largest observed barcode read fraction, with all
    observed barcodes in its denominator; it is not an estimate of cell purity.
    `top_n_for_purity` and `edit_distance_drop` apply only to separately labelled
    legacy diagnostics. Distance alone does not establish a sequencing error.
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
    purpose_b_not: str = ""  # legacy field; now requires an explicit purpose_b selector
    purpose_b: str | list[str] | None = None
    require_validation_a: bool = True


@dataclass
class REDIPlotConfig:
    """Explicit plot filters; positions do not establish independent isolates."""

    min_counts: int = 20
    min_purity: float = 0.95
    # Positive values require an independently registered replicate-count field.
    min_replicates: int = 0


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
    library_by_purpose: dict[str, str] = field(default_factory=lambda: {
        'SHA345 REDI': 'SHA345', 'yT182 REDI': 'yT182', 'yT183 REDI': 'yT183',
    })
    allow_partial_samples: bool = False
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
    # Independent source_sample→sample_name namespace registration, when needed.
    wgs_sample_namespace: Path | None = None
    # Inner-primer rules: which side of the read pair holds bc1 / REDI_bc
    primer_rules: dict[str, str] = field(default_factory=lambda: {
        "KR1965": "R2",  # KR1965 → bc1 on R2, REDI_bc on R1
        "KR1967": "R1",  # KR1967 → bc1 on R1, REDI_bc on R2
    })
    plotting: REDIPlotConfig = field(default_factory=REDIPlotConfig)


class _UniqueSafeLoader(yaml.SafeLoader):
    def flatten_mapping(self, node):
        # SafeLoader recursively calls this for inline/aliased merge sources.
        # Validate each original mapping before flattening hides its duplicates.
        checked = getattr(self, '_checked_mappings', set())
        self._checked_mappings = checked
        if id(node) not in checked:
            declared = set()
            for key_node, _ in node.value:
                key = ('<<' if key_node.tag == 'tag:yaml.org,2002:merge'
                       else self.construct_object(key_node))
                if not isinstance(key, str):
                    raise ValueError("configuration mapping keys must be strings")
                if key in declared:
                    raise ValueError(f"duplicate configuration key: {key}")
                declared.add(key)
            checked.add(id(node))
        super().flatten_mapping(node)


def _unique_mapping(loader, node, deep=False):
    # Reject repeated keys written in the same mapping, while preserving YAML's
    # standard explicit override of a merged mapping (<<: *defaults).
    loader.flatten_mapping(node)
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str):
            raise ValueError("configuration mapping keys must be strings")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


_UniqueSafeLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping
)


def _schema(value: Any, cls: type, label: str) -> dict:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a mapping")
    unknown = set(value) - {f.name for f in fields(cls)}
    if unknown:
        raise ValueError(f"unknown {label} fields: {sorted(unknown)}")
    return value


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty string")
    return value


def _integer(value: Any, label: str, minimum: int = 0) -> None:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")


def _boolean(value: Any, label: str) -> None:
    if type(value) is not bool:
        raise ValueError(f"{label} must be a boolean")


def _relative_path(value: Any, base: Path, label: str) -> Path:
    result = Path(_text(value, label)).expanduser()
    return (result if result.is_absolute() else base / result).resolve()


def load_config(path: str | Path) -> REDIConfig:
    """Load an explicit schema, binding relative paths to the configuration.

    project_root is relative to the YAML directory. Other data paths are
    relative to project_root; key filenames remain relative to paths.keyfiles.
    Unknown and duplicate settings fail before any stage reads assay data.
    """
    path = Path(path).expanduser().resolve()
    with path.open() as fh:
        data = _schema(yaml.load(fh, Loader=_UniqueSafeLoader), REDIConfig, "config")
    for required in ("screen_name", "sequencing_date", "paths"):
        if required not in data:
            raise ValueError(f"missing configuration field: {required}")
    _text(data["screen_name"], "screen_name")
    date = data["sequencing_date"]
    if type(date) is int:
        date = str(date)
    _text(date, "sequencing_date")
    paths_d = _schema(data["paths"], REDIPaths, "paths")
    missing = {f.name for f in fields(REDIPaths)} - set(paths_d)
    if missing:
        raise ValueError(f"missing paths fields: {sorted(missing)}")
    project = _relative_path(paths_d["project_root"], path.parent, "paths.project_root")
    paths = REDIPaths(project_root=project, **{
        name: _relative_path(value, project, "paths." + name)
        for name, value in paths_d.items() if name != "project_root"
    })
    barcode = REDIBarcodeSchema(**_schema(data.get("barcode", {}), REDIBarcodeSchema, "barcode"))
    purity = REDIPurityConfig(**_schema(data.get("purity", {}), REDIPurityConfig, "purity"))
    plotting = REDIPlotConfig(**_schema(data.get("plotting", {}), REDIPlotConfig, "plotting"))
    _integer(barcode.bc1_length, "barcode.bc1_length", 1)
    _integer(barcode.bc1_length_tolerance, "barcode.bc1_length_tolerance")
    if barcode.bc1_length_tolerance >= barcode.bc1_length:
        raise ValueError("bc1 length tolerance must be smaller than bc1 length")
    _boolean(barcode.allow_positional_fallback, "barcode.allow_positional_fallback")
    if barcode.fallback_prefix_length is not None or barcode.allow_positional_fallback:
        _integer(barcode.fallback_prefix_length, "barcode.fallback_prefix_length")
    for name in ("fwd_bc1_prefix", "fwd_bc1_suffix"):
        if name == 'fwd_bc1_suffix' and barcode.fwd_bc1_suffix is None:
            continue  # Explicit prefix-only legacy assay; exact length only.
        value = _text(getattr(barcode, name), "barcode." + name).upper()
        if re.fullmatch("[ACGT]+", value) is None:
            raise ValueError(f"barcode.{name} must contain A/C/G/T only")
        setattr(barcode, name, value)
    _integer(purity.top_n_for_purity, "purity.top_n_for_purity", 1)
    _integer(purity.min_counts, "purity.min_counts", 1)
    if not isinstance(purity.edit_distance_drop, (list, tuple)):
        raise ValueError("purity.edit_distance_drop must be a list of integers")
    for distance in purity.edit_distance_drop:
        _integer(distance, "purity.edit_distance_drop", 1)
    purity.edit_distance_drop = tuple(purity.edit_distance_drop)
    _integer(plotting.min_counts, "plotting.min_counts")
    _integer(plotting.min_replicates, "plotting.min_replicates")
    if (type(plotting.min_purity) not in (int, float)
            or not math.isfinite(plotting.min_purity) or not 0 <= plotting.min_purity <= 1):
        raise ValueError("plotting.min_purity must be a finite number from 0 to 1")

    raw_pairs = data.get("array_pairs", [])
    if not isinstance(raw_pairs, list):
        raise ValueError("array_pairs must be a list")
    pairs = []
    pair_names = set()
    for raw_pair in raw_pairs:
        value = _schema(raw_pair, REDIArrayPair, "array_pairs")
        if not {"name_a", "name_b", "purpose_a"} <= set(value):
            raise ValueError("array_pairs requires name_a, name_b, purpose_a")
        pair = REDIArrayPair(**value)
        for name in (pair.name_a, pair.name_b):
            if not isinstance(name, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) is None:
                raise ValueError("array names must be safe nonempty identifiers")
        if pair.name_a == pair.name_b:
            raise ValueError("array names must be distinct")
        pair_key = pair.name_a + '_' + pair.name_b
        if pair_key in pair_names:
            raise ValueError("duplicate array comparison")
        pair_names.add(pair_key)
        _text(pair.purpose_a, "purpose_a")
        if not isinstance(pair.purpose_b_not, str):
            raise ValueError("purpose_b_not must be a string")
        if pair.purpose_b is not None:
            selectors = pair.purpose_b if isinstance(pair.purpose_b, list) else [pair.purpose_b]
            if not selectors:
                raise ValueError("purpose_b cannot be an empty list")
            for selector in selectors:
                _text(selector, "purpose_b")
            if len(selectors) != len(set(selectors)):
                raise ValueError("duplicate purpose_b selectors")
            if pair.purpose_a in selectors:
                raise ValueError("array purpose selectors must be disjoint")
        _boolean(pair.require_validation_a, "require_validation_a")
        pairs.append(pair)

    kwargs = dict(data)
    kwargs.update(paths=paths, sequencing_date=date, barcode=barcode, purity=purity,
                  plotting=plotting, array_pairs=pairs)
    for name in ("bc1_reference_table", "wgs_outcome_table", "wgs_sample_namespace"):
        value = kwargs.get(name)
        kwargs[name] = None if value is None else _relative_path(value, project, name)
    config = REDIConfig(**kwargs)
    _boolean(config.allow_partial_samples, "allow_partial_samples")
    for name in ("sample_key", "gdna_plate_key", "redi_plate_barcode_key"):
        _text(getattr(config, name), name)
    if not isinstance(config.redi_bc_keys, list):
        raise ValueError("redi_bc_keys must be a list")
    for key in config.redi_bc_keys:
        _text(key, "redi_bc_keys entry")
    for name in ("library_by_purpose", "primer_rules"):
        mapping = getattr(config, name)
        if not isinstance(mapping, dict) or not mapping:
            raise ValueError(f"{name} must be a nonempty mapping")
        for key, value in mapping.items():
            _text(key, name + " key")
            _text(value, name + " value")
            if name == "primer_rules" and value not in ("R1", "R2"):
                raise ValueError("primer_rules values must be R1 or R2")
    return config
