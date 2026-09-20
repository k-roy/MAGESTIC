"""
Step 2 — donor-window editing-outcome pre-calling. **GATED — spec §8.**

Orchestration over TRACE's `trace_crispr.wgs.donor_outcome` (a port of Shengdi
Li's `precalling_target_outcome.pl`). The stage is wired end-to-end and then
REFUSES TO RUN, because three measured gates stand between it and a correct
answer on this program's data.

🔴 **The gates are DETECTED, not declared.** There is no `enabled:` switch to
flip: each gate is a property of the inputs, checked against the filesystem every
run, so the day the inputs change the stage starts working with no code edit —
and until then it cannot be opened by editing a config.

GATE 2 (§8.2) — **no gVCF exists anywhere**, so the reference depth
    `precall_sample` needs to separate `designed_variant_unedited` from `unknown`
    has no source. Checked first because it is the cheapest check and, today,
    always true: `core.vcf_adapter.precall_inputs` reports `gvcf_missing`.
    `core/vcf_adapter.py` does the re-representation half of the adapter
    (96-sample bcftools merge -> per-sample plain-text `.fil.vcf`, AD remapped from
    FORMAT slot 4 to slot 2). Synthesizing per-position reference depth from a
    pileup would *derive* a quantity, which belongs in TRACE (spec D3), not in a
    magestic workaround. Until then every uncovered call collapses to `unknown`,
    and `unknown` must never be read as `unedited`.

GATE 3 (§8.3) — reproduced empirically: `precall_sample` re-parses its own
    `f"{chrom}_{pos}_{ref}_to_{alt}"` id with `split("_")` and indexes
    positionally, so any contig name containing an underscore raises `ValueError`
    mid-loop, after a partial output is already on disk. Every mandated composite
    reference carries `plasmid_{oligo}` / `integrated_plasmid_{oligo}` contigs, so
    the stage is incompatible with the reference standard this program mandates.
    Checked before the first sample, from the samplesheet.

GATE 1 (§8.1) — the SYNTHETIC_ERROR axis has no input: no per-oligo step-2 donor
    sequence exists for the 2021 pool and none is derivable from design data (it
    is an empirical per-oligo observation). This one does not block the run; it
    makes the output CONDITIONALLY MEANINGFUL, so `core/samplesheet.py` sets
    `synthetic_donor := design_donor` and stamps `synthetic_error_axis=disabled`,
    and this step re-reads the written file to confirm the stamp survived. A zero
    in the `synthetic_errors` column of an unstamped file is a wrong answer waiting
    to be quoted.

Author: Kevin R. Roy
"""

from __future__ import annotations

import logging
from pathlib import Path

from trace_crispr.wgs import donor_outcome

from ..config import DONOR_PRECALL_SUBDIR, WGSConfig
from ..core.manifest import WGSSample
from ..core.samplesheet import (
    build_samplesheet,
    read_axis_stamp,
    stamp_outcome_table,
    underscore_contigs,
)
from ..core.vcf_adapter import precall_inputs, split_plate_vcf
from ._common import (
    GateError,
    StepIncompleteError,
    cohort_table,
    collect_shards,
    manifest_samples,
    publish,
    record_for,
    require_row_count,
    shard_dir,
    write_provenance,
    write_tsv,
)

logger = logging.getLogger(__name__)

STEP = "step_2_donor_precall"
SHARD_SUFFIX = ".donor_outcome.tsv"
TABLE_NAME = "donor_outcome.tsv"

GATE_1 = (
    "GATE 1 (spec §8.1): no per-oligo step-2 synthetic donor exists for the 2021 pool, "
    "so the SYNTHETIC_ERROR axis ships DISABLED and every output is stamped."
)
GATE_2 = (
    "GATE 2 (spec §8.2): no gVCF exists for this cohort, so `designed_variant_unedited` "
    "cannot be separated from `unknown`. The missing piece is a reference-depth "
    "primitive in TRACE, not a magestic workaround."
)
GATE_3 = (
    "GATE 3 (spec §8.3): donor_outcome.precall_sample crashes on any contig name "
    "containing '_', and every mandated composite reference carries plasmid_{oligo} / "
    "integrated_plasmid_{oligo} contigs. The fix is TRACE-side (carry structured fields "
    "instead of re-parsing a delimited string)."
)


def shard_path(config: WGSConfig, sample: str) -> Path:
    """Per-sample donor-outcome shard."""
    return shard_dir(config, DONOR_PRECALL_SUBDIR) / f"{sample}{SHARD_SUFFIX}"


def table_path(config: WGSConfig) -> Path:
    """Cohort donor-outcome table."""
    return cohort_table(config, TABLE_NAME)


def samplesheet_path(config: WGSConfig, sample: str) -> Path:
    """TRACE's 7-column donor samplesheet for one sample."""
    return config.paths.samplesheets / f"{sample}.samplesheet.tsv"


def gate_reason() -> str:
    """The string `run_pipeline` records in place of an artifact path."""
    return f"GATED: {GATE_2} | {GATE_3} | {GATE_1}"


def run_step_2(config: WGSConfig, sample: str | None = None) -> Path:
    """Pre-call donor-window editing outcomes. Returns the path written.

    Raises `GateError` while any blocking gate is open — which is the state today
    for every cohort inventoried.
    """
    if sample is not None:
        return _run_one_sample(config, sample)
    return _combine(config)


def _run_one_sample(config: WGSConfig, sample: str) -> Path:
    """One sample. Idempotent: a finished shard short-circuits every gate check."""
    out = shard_path(config, sample)
    if out.exists():
        logger.info("step_2: %s already complete — skipping", sample)
        return out

    record = record_for(config, sample)
    workdir = config.paths.donor_precall_workdir

    _require_vcf_inputs(config, record, workdir)
    sheet, axis = _build_samplesheet(config, sample)

    raw = out.with_name(out.name + ".raw")
    results = donor_outcome.precall(
        sheet,
        workdir,
        Path(record.reference),
        raw,
        ne_threshold=config.donor_precall.ne_threshold,
        quantification_window_size=config.donor_precall.quantification_window_size,
        quantification_window_cov_perc=config.donor_precall.quantification_window_cov_perc,
        samples={sample},
    )
    if len(results) != 1:
        raise StepIncompleteError(
            f"step_2: precall returned {len(results)} outcomes for {sample}, expected 1. "
            "TRACE drops a sample whose vcf/gvcf is missing with a bare `continue` — no "
            "warning, no count, and an output file that looks complete."
        )
    stamp_outcome_table(raw, synthetic_error_axis=axis)
    publish(raw, out)
    logger.info("step_2: %s outcome=%s", sample, results[0].editing_outcome)
    return out


def _require_vcf_inputs(config: WGSConfig, record: WGSSample, workdir: Path) -> None:
    """Run the VCF adapter, then refuse on GATE 2 if the gVCF half is still missing."""
    for plate_vcf in config.donor_precall.plate_vcfs:
        split_plate_vcf(
            plate_vcf,
            workdir,
            samples=[record.sample],
            ad_index=config.donor_precall.ad_format_slot_in,
        )
    inputs = precall_inputs(workdir, [record.sample])
    if inputs.gvcf_missing:
        raise GateError(f"step_2: {GATE_2} (no gVCF for {', '.join(inputs.gvcf_missing)})")
    if inputs.vcf_missing:
        raise StepIncompleteError(
            f"step_2: no per-sample VCF for {', '.join(inputs.vcf_missing)} under {workdir}. "
            "precall would skip the sample silently and still write a file."
        )


def _build_samplesheet(config: WGSConfig, sample: str) -> tuple[Path, str]:
    """Build this sample's samplesheet and check GATE 3. Returns `(path, axis)`.

    `readers.oligo_tables` is imported inside the function, not at module scope:
    it is unreachable while GATE 2 is open, and deferring it keeps this module
    importable (and therefore testable) without it. See the integration note in
    the module report — the donor-row source is the one input `config` does not
    name today.
    """
    from ..readers.oligo_tables import donor_samplesheet_rows

    rows = list(donor_samplesheet_rows(config, [sample]))
    offenders = underscore_contigs(rows)
    if offenders:
        raise GateError(f"step_2: {GATE_3} Offending contig(s): {', '.join(offenders[:5])}")

    result = build_samplesheet(
        rows,
        samplesheet_path(config, sample),
        synthetic_donor_source=config.donor_precall.synthetic_donor_source,
    )
    if result.n_rows != 1:
        raise StepIncompleteError(
            f"step_2: samplesheet for {sample} has {result.n_rows} rows, expected 1"
        )
    return result.path, result.synthetic_error_axis


def _combine(config: WGSConfig) -> Path:
    """Fan-in: every sample, then the cohort table, then re-read the GATE-1 stamp."""
    records = manifest_samples(config)
    for record in records:
        _run_one_sample(config, record.sample)

    header, rows = collect_shards(
        records, lambda r: shard_path(config, r.sample), label="step_2 donor outcome"
    )
    out = write_tsv(table_path(config), header, rows)
    n = require_row_count(out, len(rows), len(records), label="step_2 donor outcome")
    _require_axis_stamp(out)

    unknown = sum(1 for row in rows if row.split("\t")[1].startswith("unknown"))
    logger.info(
        "step_2: %d samples; %d (%.1f%% of %d) outcome=unknown. unknown is NOT unedited "
        "(GATE 2 — there is no gVCF to prove coverage).",
        n,
        unknown,
        100.0 * unknown / n if n else 0.0,
        n,
    )

    write_provenance(
        STEP,
        out,
        inputs={"manifest": Path(config.paths.manifest)},
        parameters={
            "synthetic_error_axis": read_axis_stamp(out),
            "synthetic_donor_source": config.donor_precall.synthetic_donor_source,
            "ne_threshold": config.donor_precall.ne_threshold,
            "quantification_window_size": config.donor_precall.quantification_window_size,
            "quantification_window_cov_perc": config.donor_precall.quantification_window_cov_perc,
            "gates": [GATE_1, GATE_2, GATE_3],
        },
        counts={"samples": n, "outcome_unknown": unknown},
    )
    return out


def _require_axis_stamp(path: Path) -> None:
    """Re-read the written table and confirm the GATE-1 stamp survived the fan-in."""
    axis = read_axis_stamp(path)
    if axis is None:
        raise StepIncompleteError(
            f"step_2: {path} carries no synthetic_error_axis stamp. An unstamped Stage-2 "
            "output is indistinguishable from one where the axis genuinely ran, so a zero "
            "in `synthetic_errors` would read as 'no synthesis errors found'."
        )
