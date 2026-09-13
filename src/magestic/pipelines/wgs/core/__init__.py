"""magestic.pipelines.wgs.core — adapters between magestic config/paths and TRACE.

🔴 **This layer holds ADAPTERS, not algorithms.** ARCHITECTURE_POLICY §1 says `core/`
holds "ALL analysis logic" while §11 says the algorithms live in TRACE; the conflict is
named and resolved once here so no later reviewer "restores §1 compliance" by inlining
TRACE code.

**The test:** if a file in `wgs/` would contain a *threshold*, a *classification rule*, or
would *derive a new number from sequence or alignment data*, it belongs in TRACE. An
adapter that **re-represents** an existing quantity (96-sample VCF -> per-sample; `AD` in
FORMAT slot 4 -> slot 2; nested YAML -> flat sidecar) is orchestration. An adapter that
**derives** one (synthesizing reference depth from a BAM pileup) is an algorithm and goes
to TRACE.

Every module here takes primitive arguments — `Path`, `str`, `int`, mappings, open pysam
handles — and never imports `..config`. The config object is assembled and unpacked by
`pipeline/step_*.py`, exactly as `redi.core` does it, so these functions stay callable
from a test with no config at all.

TRACE itself is reached only through `._trace`, which defers the import to call time and
records which TRACE loaded.
"""

from . import _trace
from .emit import (
    count_data_rows,
    write_bc1_calls,
    write_plasmid_integration_summary,
    write_reconciliations,
    write_refused_table,
    write_sample_anno_table,
    write_table,
)
from .locus_config import (
    LOCUS_CONFIG_KEYS,
    LocusConfigError,
    build_tracks,
    format_locus,
    load_locus_sidecar,
    parse_window,
    validate_tracks,
    write_locus_sidecar,
)
from .manifest import (
    OPTIONAL_COLUMNS,
    REFUSED_COLUMNS,
    REQUIRED_COLUMNS,
    CrossReferenceError,
    ManifestError,
    ReferenceMismatchError,
    WGSSample,
    check_reference_consistency,
    group_by_reference,
    load_manifest,
    non_portable_contigs,
    refusal_reason,
    require_reference_consistency,
    require_single_reference,
)
from .reconcile_driver import (
    BC_CONFIDENCE_HIGH_TIERS,
    bc_confidence,
    load_oligo_windows,
    reconcile_sample,
    score_sites,
)
from .samplesheet import (
    AXIS_DISABLED,
    AXIS_ENABLED,
    AXIS_STAMP_KEY,
    SAMPLESHEET_COLUMNS,
    SYNTHETIC_DONOR_FROM_COLUMN,
    SYNTHETIC_DONOR_FROM_DESIGN,
    SamplesheetError,
    SamplesheetResult,
    axis_stamp,
    build_samplesheet,
    read_axis_stamp,
    stamp_outcome_table,
    underscore_contigs,
)
from .vcf_adapter import (
    PRECALL_GVCF_SUBDIR,
    PRECALL_GVCF_SUFFIX,
    PRECALL_VCF_SUBDIR,
    PRECALL_VCF_SUFFIX,
    TRACE_AD_SLOT,
    PrecallInputs,
    VcfAdapterError,
    VcfSplitResult,
    carries_alt,
    precall_inputs,
    reorder_format_field,
    split_plate_vcf,
)

__all__ = [
    # _trace — the single TRACE import site
    "_trace",
    # manifest — D1: reference is a per-sample attribute
    "WGSSample",
    "ManifestError",
    "ReferenceMismatchError",
    "CrossReferenceError",
    "REQUIRED_COLUMNS",
    "OPTIONAL_COLUMNS",
    "REFUSED_COLUMNS",
    "load_manifest",
    "group_by_reference",
    "check_reference_consistency",
    "require_reference_consistency",
    "require_single_reference",
    "non_portable_contigs",
    "refusal_reason",
    # locus_config — the flat sidecar TRACE can parse
    "LOCUS_CONFIG_KEYS",
    "LocusConfigError",
    "format_locus",
    "parse_window",
    "write_locus_sidecar",
    "load_locus_sidecar",
    "build_tracks",
    "validate_tracks",
    # samplesheet — TRACE's 7-column donor sheet, and GATE 1's stamp
    "SAMPLESHEET_COLUMNS",
    "SYNTHETIC_DONOR_FROM_DESIGN",
    "SYNTHETIC_DONOR_FROM_COLUMN",
    "AXIS_STAMP_KEY",
    "AXIS_DISABLED",
    "AXIS_ENABLED",
    "SamplesheetError",
    "SamplesheetResult",
    "build_samplesheet",
    "axis_stamp",
    "stamp_outcome_table",
    "read_axis_stamp",
    "underscore_contigs",
    # vcf_adapter — 96-sample plate VCF -> the per-sample layout precall hardcodes
    "PRECALL_VCF_SUBDIR",
    "PRECALL_VCF_SUFFIX",
    "PRECALL_GVCF_SUBDIR",
    "PRECALL_GVCF_SUFFIX",
    "TRACE_AD_SLOT",
    "VcfAdapterError",
    "VcfSplitResult",
    "PrecallInputs",
    "carries_alt",
    "reorder_format_field",
    "split_plate_vcf",
    "precall_inputs",
    # reconcile_driver — compose the reconcile primitives without argparse
    "BC_CONFIDENCE_HIGH_TIERS",
    "bc_confidence",
    "load_oligo_windows",
    "score_sites",
    "reconcile_sample",
    # emit — TSV writers whose schema comes from TRACE
    "write_table",
    "count_data_rows",
    "write_bc1_calls",
    "write_reconciliations",
    "write_plasmid_integration_summary",
    "write_sample_anno_table",
    "write_refused_table",
]
