"""magestic.pipelines.wgs.qc — provenance sidecars and the loud-failure guards.

`provenance` implements ARCHITECTURE_POLICY §8 ("every data output carries a
PROVENANCE.json") with an explicitly resolved code sha, because the shared helper
records a null sha when the code is staged rather than cloned.

`pipeline_qc` implements §7 ("no silent success"): each guard maps to a measured way
this module's inputs can produce a clean-looking but wrong answer — a header-only table
with exit code 0, a silently skipped sample, an unstamped gated axis.
"""

from .pipeline_qc import (
    QCFailure,
    require_axis_stamp,
    require_manifest_coverage,
    require_no_underscore_contigs,
    require_nonempty_output,
    require_samtools,
    summarize_bc1_table,
)
from .provenance import (
    RUN_SIDECAR_NAME,
    SIDECAR_SUFFIX,
    ProvenanceError,
    build_sha,
    code_state,
    write_run_provenance,
    write_step_provenance,
)

__all__ = [
    "RUN_SIDECAR_NAME",
    "SIDECAR_SUFFIX",
    "ProvenanceError",
    "QCFailure",
    "build_sha",
    "code_state",
    "require_axis_stamp",
    "require_manifest_coverage",
    "require_no_underscore_contigs",
    "require_nonempty_output",
    "require_samtools",
    "summarize_bc1_table",
    "write_run_provenance",
    "write_step_provenance",
]
