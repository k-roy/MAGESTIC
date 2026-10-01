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
from collections.abc import Sequence
from pathlib import Path

from trace_crispr.wgs import donor_outcome

from ..config import DONOR_PRECALL_SUBDIR, WGSConfig
from ..core.blacklist_adapter import reference_contig_lengths
from ..core.contig_alias import ContigAlias, alias_rows
from ..core.donor_precall import alias_fasta, alias_vcf_contigs, audit_vcf_window
from ..core.gvcf_adapter import expand_reference_blocks
from ..core.manifest import WGSSample
from ..core.outcome_reason import explain_outcome
from ..core.samplesheet import (
    build_samplesheet,
    read_axis_stamp,
    stamp_outcome_table,
    underscore_contigs,
)
from ..core.vcf_adapter import (
    CALLER_PLOIDY,
    PLOIDY_CAVEAT,
    PRECALL_GVCF_SUBDIR,
    PRECALL_GVCF_SUFFIX,
    PRECALL_VCF_SUBDIR,
    PRECALL_VCF_SUFFIX,
    caller_provenance,
    precall_inputs,
    split_plate_vcf,
)
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

    minor_alleles = _require_vcf_inputs(config, record, workdir)
    alias = _gate_3_alias(config, record)
    sheet, axis = _build_samplesheet(config, sample, alias=alias)

    genome = Path(record.reference)
    vcf = workdir / PRECALL_VCF_SUBDIR / f"{sample}{PRECALL_VCF_SUFFIX}"
    gvcf = workdir / PRECALL_GVCF_SUBDIR / f"{sample}{PRECALL_GVCF_SUFFIX}"
    if alias is not None:
        # PER-SAMPLE path, not `{stem}.aliased.fa`: step_2 fans out as a SLURM array,
        # so every sample on the same reference would otherwise write and index the
        # same FASTA concurrently — a torn file or a half-written `.fai`. Costs one
        # small FASTA per sample and removes the race entirely.
        genome = workdir / f"{sample}.{genome.stem}.aliased.fa"
        alias_fasta(record.reference, genome, alias)
        for path in (vcf, gvcf):
            staged = path.with_name(path.name + ".aliased")
            alias_vcf_contigs(path, staged, alias)
            staged.replace(path)
        _require_alias_consistency(alias, sheet, vcf, gvcf, genome)
        logger.info(
            "step_2: %s GATE 3 resolved by ContigAlias over %d aliased contig(s) "
            "(%s...); variant ids are DECODED before the shard is published",
            sample, len(alias.aliased), ", ".join(alias.aliased[:3]),
        )

    raw = out.with_name(out.name + ".raw")
    results = donor_outcome.precall(
        sheet,
        workdir,
        genome,
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
    if alias is not None:
        n_decoded = _decode_outcome_table(raw, alias)
        logger.info("step_2: %s decoded %d aliased variant id(s)", sample, n_decoded)
    _annotate_outcome_reason(
        config, record, workdir, raw, results[0], minor_alleles=minor_alleles
    )
    stamp_outcome_table(raw, synthetic_error_axis=axis)
    publish(raw, out)
    logger.info("step_2: %s outcome=%s", sample, results[0].editing_outcome)
    return out


#: Columns NN 251 appends to TRACE's 5-column Stage-2 table (defects 4 and 8).
#: Appended, never inserted, and `step_5` reads `editing_outcome` by NAME
#: (`step_5_join_outcomes.py:139`), so no consumer shifts.
OUTCOME_REASON_COLUMNS = (
    "callable",
    "outcome_reason",
    "qwin_covered",
    "qwin_size",
    "n_variants_discarded_outside_donor_window",
    "minor_allele_variants",
)


def _annotate_outcome_reason(
    config: WGSConfig, record: WGSSample, workdir: Path, table: Path, outcome: object,
    *, minor_alleles: Sequence[str] = (),
) -> None:
    """Append WHY this outcome is what it is to the shipped shard (NN 251 defect 4).

    🔴 **In the TABLE, not only in the log.** `unknown` is the only value in the
    Stage-2 table that means "I could not tell", and TRACE records nothing about
    which cause produced it — no gVCF, unexpanded reference blocks, depth below
    `ne_threshold`, or a genuinely unsequenced locus all collapse to the same five
    letters. A reason that reached only the run log would leave every downstream
    reader unable to tell "uncallable" from "unedited", which is the single most
    load-bearing distinction in this analysis.

    Also carries defect 8's count: how many VCF records fell outside the donor window
    and were therefore discarded in silence, so the VCF reconciles with the outcome.
    """
    sheet = samplesheet_path(config, record.sample)
    try:
        lines = sheet.read_text().splitlines()
        cols = lines[0].split("\t")
        row = lines[1].split("\t")
        chrom = row[cols.index("chrom")].strip()
        dstart = int(row[cols.index("donor_start_coord")])
        dend = int(row[cols.index("donor_end_coord")])
    except (OSError, IndexError, ValueError):
        # Even the failure is explicit in the table: a shard whose reason columns are
        # absent would be indistinguishable from one where the reason was computed.
        _append_columns(table, OUTCOME_REASON_COLUMNS,
                        ["", "samplesheet_unreadable", "", "", "",
                         ";".join(minor_alleles)])
        logger.warning("step_2: %s could not read its samplesheet to explain the "
                       "outcome; the flag stands but carries no reason", record.sample)
        return
    ex = explain_outcome(
        getattr(outcome, "editing_outcome", "unknown"),
        gvcf_path=workdir / PRECALL_GVCF_SUBDIR / f"{record.sample}{PRECALL_GVCF_SUFFIX}",
        chrom=chrom, donor_start=dstart, donor_end=dend, sample=record.sample,
        quantification_window_size=config.donor_precall.quantification_window_size,
        quantification_window_cov_perc=(
            config.donor_precall.quantification_window_cov_perc
        ),
        ne_threshold=config.donor_precall.ne_threshold,
    )
    audit = audit_vcf_window(
        workdir / PRECALL_VCF_SUBDIR / f"{record.sample}{PRECALL_VCF_SUFFIX}",
        chrom=chrom, donor_start=dstart, donor_end=dend,
    )
    _append_columns(
        table,
        OUTCOME_REASON_COLUMNS,
        [
            str(ex.callable_).lower(),
            ex.reason,
            str(ex.n_window_positions_covered),
            str(ex.n_window_positions),
            str(audit.n_discarded),
            ";".join(minor_alleles),
        ],
    )
    logger.info(
        "step_2: %s callable=%s reason=%s qwin=%d-%d covered=%d/%d",
        record.sample, ex.callable_, ex.reason, ex.qwin_start, ex.qwin_end,
        ex.n_window_positions_covered, ex.n_window_positions,
    )
    if audit.n_discarded:
        # NN 251 defect 8: distal unintended edits are invisible to Stage 2. Counted
        # here rather than discarded, so the VCF's records reconcile with the outcome.
        logger.info("step_2: %s %s", record.sample, audit.summary())


def _append_columns(table: Path, columns: Sequence[str], values: Sequence[str]) -> None:
    """Append fixed columns to TRACE's Stage-2 table, header and data rows alike.

    TRACE writes `#sample_name\\t...` as its header, so the header is the `#`-prefixed
    line and every other non-blank line is data. Appending (never inserting) keeps
    `step_5`'s `precall[sample]["editing_outcome"]` name-based lookup correct, and
    re-running is safe because a shard that already carries the columns is skipped.
    """
    lines = table.read_text().splitlines()
    if any(line.startswith("#") and columns[0] in line.split("\t") for line in lines):
        return  # idempotent: already annotated
    out: list[str] = []
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("#"):
            out.append(line + "\t" + "\t".join(columns))
        else:
            out.append(line + "\t" + "\t".join(values))
    with table.open("w", newline="\n") as fh:
        for line in out:
            fh.write(line + "\n")


def _require_vcf_inputs(
    config: WGSConfig, record: WGSSample, workdir: Path
) -> tuple[str, ...]:
    """Run the VCF adapter, then refuse on GATE 2 if the gVCF half is still missing.

    🔴 NN 251 defect 17 (found while fixing defect 1, not on the original inventory):
    this call previously read

        split_plate_vcf(..., ad_index=config.donor_precall.ad_format_slot_in)

    which has **two** bugs in one line. `split_plate_vcf` has no `ad_index`
    parameter — it is `ad_slot` — so the call was a guaranteed `TypeError`; and
    `ad_format_slot_in` is 4, the slot bcftools writes `AD` to, where the function
    wants the DESTINATION slot, which is `ad_format_slot_out` (2, where TRACE's
    positional gVCF parser reads reference depth). Passing 4 would have left `AD`
    where it already was and made TRACE read `PL`'s first comma element as the
    reference depth. Both halves are fixed here and asserted by
    `test_fixed_step_2_passes_the_destination_ad_slot`.

    🔴 NN 251 defect 14: `VcfSplitResult.n_multiallelic` was computed and never
    read. TRACE reads ALT whole, so a multi-allelic `"A,G"` reaches
    `bucket_variant` as a literal two-character allele and every Levenshtein
    distance from it is wrong. Refused here rather than discarded.
    """
    multiallelic = 0
    minor: tuple[str, ...] = ()
    for plate_vcf in config.donor_precall.plate_vcfs:
        _require_recoverable_caller(config, plate_vcf)
        split = split_plate_vcf(
            plate_vcf,
            workdir,
            samples=[record.sample],
            ad_slot=config.donor_precall.ad_format_slot_out,
            min_allele_fraction=config.donor_precall.min_allele_fraction,
            min_alt_reads=config.donor_precall.min_alt_reads,
        )
        if (
            split.n_rescued_by_af or split.n_no_ad or split.n_dropped_below_af
            or split.n_flagged_minor or split.n_dropped_below_min_alt_reads
            or split.n_no_ad_carrying_alt
        ):
            logger.info(
                "step_2: %s allele fraction — %d of %d record(s) had a usable AD; "
                "flagged MAG_MINOR %d (of which %d ALREADY CARRIED THE ALT IN GT, the "
                "NN 252 F11b over-call class); rescued %d, below fraction %d, below "
                "%d alt read(s) %d, no AD %d, no AD WITH an alt GT %d; "
                "minor-allele variants %s",
                record.sample, split.n_af_computed, split.n_records_in,
                split.n_flagged_minor, split.n_flagged_minor_carrying_alt,
                split.n_rescued_by_af, split.n_dropped_below_af,
                config.donor_precall.min_alt_reads,
                split.n_dropped_below_min_alt_reads, split.n_no_ad,
                split.n_no_ad_carrying_alt,
                split.minor_allele_variants.get(record.sample, ()),
            )
        multiallelic += split.n_multiallelic
        minor += split.minor_allele_variants.get(record.sample, ())
    if multiallelic:
        raise StepIncompleteError(
            f"step_2: {multiallelic} multi-allelic record(s) in the plate VCF(s) for "
            f"{record.sample}. TRACE reads the ALT column whole, so an ALT of 'A,G' "
            "reaches bucket_variant as a literal allele and every edit distance from "
            "it is wrong. Normalize upstream with `bcftools norm -m-` and re-run. "
            "(VcfSplitResult.n_multiallelic was computed and discarded before NN 251.)"
        )
    inputs = precall_inputs(workdir, [record.sample])
    if inputs.gvcf_missing:
        raise GateError(f"step_2: {GATE_2} (no gVCF for {', '.join(inputs.gvcf_missing)})")
    # 🔴 NN 251 defect 5, wired on the path Stage 2 actually takes. TRACE's
    # `flag_covered` compares only POS, so a NORMAL compressed gVCF — one record per
    # reference block carrying `INFO END=` — reports `unknown` for a well-covered
    # unedited cell. Expanding blocks to positions is pure re-representation, and it
    # must happen BEFORE `precall` reads the file. No gVCF exists in this cohort today
    # (GATE 2 above), so this is a landmine being defused, not a live correction.
    gvcf = workdir / PRECALL_GVCF_SUBDIR / f"{record.sample}{PRECALL_GVCF_SUFFIX}"
    if gvcf.exists():
        staged = gvcf.with_name(gvcf.name + ".expanded")
        res = expand_reference_blocks(gvcf, staged)
        if res.n_blocks_expanded:
            staged.replace(gvcf)
            logger.info(
                "step_2: %s expanded %d gVCF reference block(s) to %d position record(s) "
                "(widest %d bp); without this a compressed gVCF reports `unknown` for a "
                "covered unedited cell",
                record.sample, res.n_blocks_expanded, res.n_records_out,
                res.max_block_width,
            )
        else:
            staged.unlink(missing_ok=True)
        if res.n_malformed_end:
            raise StepIncompleteError(
                f"step_2: {record.sample}'s gVCF has {res.n_malformed_end} record(s) "
                "with an unparseable `END=`. Treating those as per-position would "
                "UNDER-report coverage and manufacture an `unknown`, so the run stops "
                "instead."
            )
    if inputs.vcf_missing:
        raise StepIncompleteError(
            f"step_2: no per-sample VCF for {', '.join(inputs.vcf_missing)} under {workdir}. "
            "precall would skip the sample silently and still write a file."
        )
    # 🔴 NN 251 defect 1, second half. A RESCUED minor allele reaches Stage 2 as a
    # 100% call, because TRACE reads only VCF columns 0/1/3/4 — never GT, FILTER,
    # INFO or FORMAT. Returning the ids so they reach the SHARD is what stops a
    # 10% allele from shipping as `designed_variant_installed`, which is the live
    # over-call class of 242 section 7a (`AD=29,1` read as installed).
    return tuple(dict.fromkeys(minor))


def _require_recoverable_caller(config: WGSConfig, plate_vcf: Path) -> None:
    """Read the plate VCF's OWN caller command line and refuse an unrecoverable combination.

    🔴 **Detected, not declared — NN 255.** Which ploidy these VCFs were called at is
    the premise the whole of defect 1 rests on; NN 251 asserted it in a docstring and
    NN 252 measured the opposite. `bcftools` writes `##bcftools_callCommand=` into every
    VCF it calls and `bcftools merge` preserves it, so the question is answerable from
    the file and is answered here, every run.

    Two measured facts drive the two refusals (NN 255 floor sweep, bcftools 1.21, flat
    40x, alt counts 0..40):

    * At `--ploidy 1` every minority allele is `GT 0`, so `carries_alt` drops it and the
      ALLELE-FRACTION RESCUE is the only thing that can see it. `--ploidy 1` with
      `min_allele_fraction=None` therefore DISCARDS every partial edit while looking
      like a clean haploid run — refused.
    * With `-v`, a site at which no sample was called non-reference is not emitted at
      all, so a mixture private to one well is unrecoverable at any threshold. Measured
      private-mixture emission floor at 40x: `-mv` 6/40 = 15%, `-mv --ploidy 1` 34/40 =
      85%, `-m --ploidy 1 -A` 1/40 = 2.5% (magnitudes are sequence-context dependent;
      the ORDERING is what is stable). **86.6% of wells are alone at their designed
      position within their Tn5 plate** — measured on the 20231213 plate-2348/2349
      cohort (n=768), the only cohort whose per-well table carries an oligo assignment;
      `spg_20240702_POST`'s keyfile names only a yL sublibrary, which pools many oligos,
      so its private fraction can only be HIGHER. Either way the PRIVATE floor is the
      one that applies — warned loudly, with the number, because it bounds the
      correct-edit rate and no code change can lift it.
    """
    prov = caller_provenance(plate_vcf)
    if prov.n_call_commands > 1:
        raise StepIncompleteError(
            f"step_2: {Path(plate_vcf).name} carries {prov.n_call_commands} different "
            "`##bcftools_callCommand=` lines, so its records were not all called the "
            "same way and no single detection floor describes the file. Re-call the "
            "plate in one pass."
        )
    if prov.ploidy == CALLER_PLOIDY and config.donor_precall.min_allele_fraction is None:
        raise StepIncompleteError(
            f"step_2: {Path(plate_vcf).name} was called at `--ploidy {prov.ploidy}` "
            "(read from its own header) but `donor_precall.min_allele_fraction` is "
            "None. At haploid ploidy EVERY minority allele is `GT 0`, so the "
            "allele-fraction rescue is the only thing that can see one, and with it "
            "disabled every partial edit is discarded while the run looks clean. Set "
            "`min_allele_fraction` (0.05 with the default "
            f"min_alt_reads={config.donor_precall.min_alt_reads} floor), or re-call "
            "the plate at the diploid default."
        )
    if prov.ploidy is None:
        logger.warning(
            "step_2: NOT VERIFIED: %s records no `##bcftools_callCommand=` (or an "
            "unparseable `--ploidy`), so the caller's ploidy is UNOBSERVED, not "
            "diploid. %s",
            Path(plate_vcf).name, PLOIDY_CAVEAT,
        )
    elif prov.variants_only:
        logger.warning(
            "step_2: %s was called with `-v` (`%s`), so a site at which NO sample was "
            "called non-reference carries no record and a minority allele PRIVATE to "
            "one well is unrecoverable at any threshold. Measured private-mixture "
            "emission floor at 40x: 6/40 (15%%) at the diploid default, 34/40 (85%%) at "
            "`--ploidy 1`; `-m --ploidy 1 -A` (no `-v`) reaches 1/40 (2.5%%). 86.6%% of "
            "wells are alone at their designed position within their Tn5 plate "
            "(measured on the 20231213 plate-2348/2349 cohort, n=768), so this is the "
            "floor that applies. Regenerate with `MPILEUP_CALL_WINDOWS_COMMAND` "
            "(all-sites RESTRICTED to the donor windows) to lift it; the unrestricted "
            "all-sites form is ~1.2e7 records per plate and is re-read once per sample.",
            Path(plate_vcf).name, prov.command,
        )
    else:
        logger.info(
            "step_2: %s caller provenance ploidy=%s variants_only=%s keep_alts=%s — "
            "private minor alleles ARE recoverable from this file",
            Path(plate_vcf).name, prov.ploidy, prov.variants_only, prov.keep_alts,
        )


def _build_samplesheet(
    config: WGSConfig, sample: str, *, alias: ContigAlias | None = None
) -> tuple[Path, str]:
    """Build this sample's samplesheet and check GATE 3. Returns `(path, axis)`.

    `alias` is the GATE-3 remedy from `_gate_3_alias`. When it is `None` the gate
    refuses on an underscore contig exactly as before; when it is supplied the rows'
    `chrom` is encoded through it and `_run_one_sample` decodes the output.

    `readers.oligo_tables` is imported inside the function, not at module scope:
    it is unreachable while GATE 2 is open, and deferring it keeps this module
    importable (and therefore testable) without it. See the integration note in
    the module report — the donor-row source is the one input `config` does not
    name today.
    """
    from ..readers.oligo_tables import donor_samplesheet_rows

    rows = list(donor_samplesheet_rows(config, [sample]))
    offenders = underscore_contigs(rows)
    if offenders and alias is None:
        raise GateError(f"step_2: {GATE_3} Offending contig(s): {', '.join(offenders[:5])}")
    if alias is not None:
        rows = alias_rows(rows, alias, key=config.columns.chrom_col)

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


def _gate_3_alias(config: WGSConfig, record: WGSSample) -> ContigAlias | None:
    """The contig alias GATE 3 needs for this reference, or `None` when unnecessary.

    🔴 NN 251 defect 3 + 14. `core/contig_alias.ContigAlias` is a bijective,
    self-asserting adapter that solves GATE 3, and before NN 251 it was imported by
    nothing — the gate simply refused. It is wired here.

    Built over **every** contig in the reference, not just this sample's donor
    contig: the alias must be a property of the reference so two samples on the same
    reference get the same mapping and their outputs are comparable.

    Returns `None` when the reference needs no aliasing, in which case Stage 2 runs
    byte-identically to before (`ContigAlias.is_noop`).
    """
    if not config.donor_precall.resolve_underscore_contigs:
        return None
    lengths = reference_contig_lengths(record.reference)
    if not lengths:
        raise StepIncompleteError(
            f"step_2: no contig names readable for {record.reference}; GATE 3 cannot be "
            "resolved without the reference's own contig list (an alias built from one "
            "sample's rows would not be comparable across samples)."
        )
    alias = ContigAlias.from_names(lengths)
    return None if alias.is_noop else alias


def _require_alias_consistency(
    alias: ContigAlias, sheet: Path, vcf: Path, gvcf: Path, fasta: Path
) -> None:
    """Refuse unless the alias was applied CONSISTENTLY to every input TRACE will read.

    🔴 This is the guard that keeps wiring the alias from being worse than the crash it
    replaces. `precall_sample` skips a record whose `v_chr != chrom` **silently**, so a
    half-applied alias produces no error, no variants, and a plausible
    `unedited`/`unknown` — a wrong answer where there used to be a loud `ValueError`.

    The samplesheet `chrom` must be an alias this reference produced, and must be in the
    aliased FASTA. For the two VCFs the test is **consistency, not presence**:

    * a RAW name the reference carries and the alias maps away is fatal (NN 255 F6 — set
      membership over CHROM values could not see this, because both forms can coexist);
    * a name the reference does not carry at all is counted and named (F6b);
    * **the donor contig being absent is NOT an error** (NN 255 F5). An unedited cell
      legitimately has variants elsewhere and none at its donor window, and the old
      check refused exactly that — the majority class, on the path NN 251 itself called
      "the riskiest thing in this change set".
    """
    import pysam

    header = sheet.read_text().splitlines()
    if len(header) < 2:
        raise StepIncompleteError(f"step_2: samplesheet {sheet} has no data row")
    cols = header[0].split("\t")
    chrom = header[1].split("\t")[cols.index("chrom")].strip()
    if chrom not in alias.decode_map:
        raise StepIncompleteError(
            f"step_2: samplesheet contig {chrom!r} is not an alias this reference "
            "produced — the alias was applied inconsistently, which would make "
            "precall_sample skip every variant in silence."
        )
    with pysam.FastaFile(str(fasta)) as fa:
        if chrom not in set(fa.references):
            raise StepIncompleteError(
                f"step_2: aliased contig {chrom!r} is absent from the aliased FASTA "
                f"{fasta}. precall_sample would raise KeyError, or worse skip silently."
            )
    aliasable = set(alias.aliased)
    for path, label in ((vcf, "vcf"), (gvcf, "gvcf")):
        seen = {
            line.split("\t", 1)[0]
            for line in path.read_text().splitlines()
            if line and not line.startswith("#")
        }
        # 🔴 NN 255 F5 + F6. The old check was `if seen and chrom not in seen: raise`,
        # which asks the WRONG QUESTION. "Is there a variant on the donor contig?" is
        # not the same question as "was the alias applied consistently": an UNEDITED
        # cell — the majority class, and the one `designed_variant_unedited` exists for
        # — legitimately carries variants elsewhere in the genome and NONE at its donor
        # window, and NN 252 measured that the guard hard-stopped exactly that input.
        # Its severity was also inverted by the `if seen` short-circuit: a header-only
        # VCF passed unchecked while one off-contig record was fatal.
        #
        # The right question is answerable from the alias itself. A name the reference
        # DOES carry and that needed an alias must never still be present in raw form:
        # that, and only that, is the half-applied alias whose records TRACE would skip
        # in silence (`v_chr != chrom`). It is also F6, which set membership over CHROM
        # values could not see, because the raw and aliased names can both be present.
        unaliased = sorted(seen & aliasable)
        if unaliased:
            raise StepIncompleteError(
                f"step_2: the aliased {label} {path.name} still carries RAW contig "
                f"name(s) {unaliased[:3]} that this reference's alias maps to "
                f"{[alias.encode(c) for c in unaliased[:3]]}. The alias was applied "
                "to some records and not others, so precall_sample would skip the "
                "un-aliased ones in silence (v_chr != chrom) and produce a plausible "
                "`unedited`/`unknown` instead of an error. Re-run alias_vcf_contigs "
                "over the whole file."
            )
        unknown = sorted(c for c in seen if c not in alias.decode_map)
        if unknown:
            # F6b: `alias.encode` is identity for a name absent from the reference's
            # `.fai`, so a VCF contig the composite reference does not carry keeps its
            # underscores. Harmless for TRACE's positional id parse only because
            # `precall_sample` filters `v_chr != chrom` first — a coincidence, not a
            # guard — so it is counted and named rather than passed over.
            logger.warning(
                "step_2: NOT VERIFIED: the aliased %s %s carries %d CHROM value(s) "
                "absent from the reference's own contig list, so the alias is identity "
                "for them and they keep any underscores: %s. TRACE skips them because "
                "they are not the donor contig, which is a coincidence rather than a "
                "guard; they belong in the reference or out of the VCF.",
                label, path.name, len(unknown), unknown[:5],
            )
        if chrom not in seen:
            # The UNEDITED class. Not an error, and the loudest thing it deserves is a
            # count: the alias has already been proven consistent above.
            logger.info(
                "step_2: %s carries no record on the donor contig %r (%d other CHROM "
                "value(s) present). That is the `designed_variant_unedited` class, not "
                "a mis-applied alias — the alias was verified consistent above.",
                path.name, chrom, len(seen),
            )


def _decode_outcome_table(path: Path, alias: ContigAlias) -> int:
    """Map every variant id in a Stage-2 outcome table back to real contig names.

    🔴 The half the peer who wrote `ContigAlias` deliberately stopped short of, and
    was right to: wiring the alias WITHOUT this would ship `ctg0003_131_A_to_G` into
    a shipped table — a silent wrong name in place of a loud crash.

    TRACE's schema is `#sample_name, editing_outcome, designed_variants,
    synthetic_errors, spontaneous_variants`; the last three are `;`-joined variant
    ids. Returns the number of ids rewritten, so the re-representation has a count.
    """
    lines = path.read_text().splitlines()
    n = 0
    out: list[str] = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            out.append(line)
            continue
        cols = line.split("\t")
        for i in (2, 3, 4):
            if i >= len(cols) or not cols[i].strip():
                continue
            ids = cols[i].split(";")
            decoded = [alias.decode_variant_id(v) for v in ids]
            n += sum(1 for a, b in zip(ids, decoded) if a != b)
            cols[i] = ";".join(decoded)
        out.append("\t".join(cols))
    with path.open("w", newline="\n") as fh:
        for line in out:
            fh.write(line + "\n")
    return n


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
