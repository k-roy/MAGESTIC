"""
magestic.pipelines.wgs.config — YAML config schema for one MAGESTIC WGS cohort run.

Defines the per-cohort YAML schema (`workflow/config/screens/<cohort>_wgs.yaml`) parsed
into a typed `WGSConfig` dataclass, and owns every pipeline-specific parameter that
otherwise squats as a default in the algorithm layer (ARCHITECTURE_POLICY §11; spec
`166_magestic_wgs_module.md` §2.6).

Both orchestration paths read this one file: the DIRECT path (`pipeline/run_pipeline.py`)
through `load_config()`, and the SNAKEMAKE path (`snakemake/Snakefile.rules`) by
`yaml.safe_load` on the same document — the Snakefile parses it only to compute file
names and the fan-out sample list, never to make an analysis decision (§4.1).

🔴 REFERENCE IS A PER-SAMPLE ATTRIBUTE (spec D1). There is deliberately no
`genome:` / `reference:` scalar on this config. 1137/1137 cohort BAMs were swept and one
directory spans two references, so a cohort-level reference is structurally wrong.
`load_config()` RAISES on either key rather than ignoring it. Per-sample reference comes
from the manifest (`paths.manifest`, resolved by `core/manifest.py`).

🔴 WHAT BELONGS HERE, AND WHAT DOES NOT. This config owns exactly four categories:

  (a) filesystem layout and input paths,
  (b) input-table column vocabulary (the schema of files *we* write or curate),
  (c) external binaries that cannot be expressed as Python dependencies,
  (d) the parameters spec §2.6 names explicitly — `bam_suffix`, `junction_distance`,
      `bp_to_extend`, `flexible_bases`, `ne_threshold`, `quantification_window_size`,
      `quantification_window_cov_perc`, `dup_mode`, the vcf/gvcf layout, and the
      `{oligo_name, chrom, donor_start_coord, donor_end_coord}` vocabulary.

Everything else keeps TRACE's default and is NOT re-declared here. Deliberately absent,
so a later reviewer does not "complete" the config by adding them: `classify_site`'s
`tau_full` / `tau_wt` / `min_cov`, `window_edit_vaf`'s `min_base_qual`,
`confident_barcodes`' `min_support` / `high_only`, `extract_track`'s `window_pad`, and
`BridgeIndex.from_tsv`'s `min_count`. Each is a threshold or a decision rule, and by
spec D3 a threshold lives in TRACE, not in the orchestration layer. The rule is
"named in §2.6 → here; not named → TRACE", applied with no exceptions.

Author: Kevin R. Roy
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

# --- Scratch filesystem ------------------------------------------------------
# Sherlock/SCG stage I/O through $SCRATCH; falls back to /tmp off-cluster.
SCRATCH = Path(os.environ.get("SCRATCH", "/tmp"))


def get_scratch_workdir(cohort_name: str) -> Path:
    """Return a per-cohort working directory under $SCRATCH."""
    return SCRATCH / "magestic_wgs" / cohort_name


# --- Output subdirectory names ----------------------------------------------
# 🔴 `snakemake/Snakefile.rules` re-derives these same directories from the raw YAML
# (it cannot import magestic). It MUST compose them from these literals. If the two
# derivations drift, paths (A) and (B) write to different directories and the §4
# byte-identity test fails looking like a *data* difference rather than a path bug.
BC1_SUBDIR = "bc1_calls"
DONOR_PRECALL_SUBDIR = "donor_precall"
SV_PLASMID_SUBDIR = "sv_plasmid"
RECONCILE_SUBDIR = "reconcile"
JOIN_SUBDIR = "joined_outcomes"
TABLES_SUBDIR = "tables"
LOCUS_CONFIG_SUBDIR = "locus_config"
SAMPLESHEET_SUBDIR = "samplesheets"
QC_SUBDIR = "qc"
LOGS_SUBDIR = "logs"

# --- The donor-precall workdir layout TRACE hardcodes ------------------------
# `donor_outcome.precall()` builds `{workdir}/vcf/{sample}.fil.vcf` and
# `{workdir}/gvcf/{sample}.g.vcf` from relative literals, and opens them with plain
# `open()` — so PLAIN TEXT, never bgzipped. These are not knobs: they mirror TRACE's
# hardcoded strings, and changing one requires a TRACE change, not a YAML edit.
# `core/vcf_adapter.py` writes into exactly this layout.
VCF_SUBDIR = "vcf"
GVCF_SUBDIR = "gvcf"
VCF_SUFFIX = ".fil.vcf"
GVCF_SUFFIX = ".g.vcf"

# --- Manifest schema (spec D1) ----------------------------------------------
# The manifest is the ONLY source of truth for which BAMs exist, which reference each
# was aligned against, and per-plate n (plate geometry differs by cohort: the v11 cohort
# is 96-complete, the rearray cohort is sparse — spec §8.6).
#
# Its COLUMN SCHEMA is deliberately NOT declared here: `core.manifest.REQUIRED_COLUMNS`
# and `core.manifest.OPTIONAL_COLUMNS` are its single home (ARCHITECTURE_POLICY §11 —
# exactly one home per definition). This config owns only the manifest's LOCATION,
# `paths.manifest`.

# --- Stage-2 SYNTHETIC_ERROR axis (spec §8.1, GATE 1) ------------------------
# No per-oligo step-2 / synthetic donor sequence exists for the 2021 pool, and none can
# be derived from design data (ruled out on all 40,191 rows). What ships is
# `synthetic_donor := design_donor`, which TRACE's own `test_precall_partial` does: it
# makes the SYNTHETIC_ERROR bucket unreachable and collapses cleanly to
# DESIGNED/SPONTANEOUS. The output must therefore be stamped as having the axis
# DISABLED — never as having found no synthesis errors.
SYNTHETIC_DONOR_FROM_DESIGN = "design_donor"
AXIS_DISABLED = "disabled"
AXIS_ENABLED = "enabled"

# --- Cluster resource keys ---------------------------------------------------
# One key per stage; `snakemake/Snakefile.rules` reads config["resources"][<key>].
STEP_RESOURCE_KEYS: tuple[str, ...] = (
    "step_1_bc1",
    "step_2_donor_precall",
    "step_3_sv_plasmid",
    "step_4_reconcile",
    "step_5_join",
)

#: Stage id (as used by `--steps` / `run_pipeline.STEPS`) -> its resource key. The ids
#: are the vocabulary `validate(steps=...)` scopes against, so they are declared here
#: beside the keys rather than re-derived from a string prefix.
STEP_RESOURCE_KEYS_BY_ID: dict[str, str] = {
    "1": "step_1_bc1",
    "2": "step_2_donor_precall",
    "3": "step_3_sv_plasmid",
    "4": "step_4_reconcile",
    "5": "step_5_join",
}

DEFAULT_RESOURCES: dict[str, dict[str, Any]] = {
    "step_1_bc1": {"mem_mb": 8000, "time": "01:00:00", "cpus": 1},
    "step_2_donor_precall": {"mem_mb": 16000, "time": "02:00:00", "cpus": 1},
    "step_3_sv_plasmid": {"mem_mb": 16000, "time": "02:00:00", "cpus": 2},
    "step_4_reconcile": {"mem_mb": 8000, "time": "01:00:00", "cpus": 1},
    "step_5_join": {"mem_mb": 8000, "time": "00:30:00", "cpus": 1},
}

# Top-level YAML keys this schema understands. An unrecognized key is a loud error
# (ARCHITECTURE_POLICY §7 — no silent success); a typo'd key must not silently no-op.
_KNOWN_TOP_LEVEL_KEYS: frozenset[str] = frozenset(
    {
        "cohort_name",
        "sequencing_date",
        "paths",
        "columns",
        "locus",
        "bc1",
        "donor_precall",
        "sv_plasmid",
        "reconcile",
        "resources",
        "python",
    }
)

# Keys that are forbidden outright, with the reason reported to the user.
_FORBIDDEN_TOP_LEVEL_KEYS: dict[str, str] = {
    "genome": (
        "reference is a PER-SAMPLE attribute carried in the manifest, never a "
        "cohort-level scalar (spec D1: one cohort directory was measured to span two "
        "references). Put a `reference` column in paths.manifest instead."
    ),
    "reference": (
        "reference is a PER-SAMPLE attribute carried in the manifest, never a "
        "cohort-level scalar (spec D1). Put a `reference` column in paths.manifest "
        "instead; use paths.reference_root only as a search root for that column."
    ),
    "bam_dir": (
        "BAMs are enumerated from paths.manifest, never from a directory glob — a glob "
        "whose suffix matches nothing writes a header-only TSV and exits 0 (spec F1)."
    ),
}


@dataclass
class WGSPaths:
    """Filesystem layout for one WGS cohort.

    Only `project_root` and `manifest` are required; every output location is derived
    (`@property`), never stored and never hardcoded at a call site.
    """

    project_root: Path  # cohort root, e.g. <projects>/QTL/WGS_barcoded_alignment/<cohort>/
    manifest: Path  # per-sample manifest — the ONLY source of BAM + reference (D1)
    intermediate: Path  # intermediate_data/ (adapters, sidecars, split VCFs)
    processed: Path  # processed_data/ (every shipped table)
    keyfiles: Path  # keyfiles/ (oligo tables, plate keys, blacklists)
    reference_root: Path | None = None
    """Search root for resolving RELATIVE `reference` values from the manifest.

    🔴 This is not "the" reference and must never be used as one. It exists so a manifest
    can carry `t2t/<strain>/<strain>_v1_T2T_full.fasta` instead of an absolute path;
    absolute manifest values ignore it entirely. `core/manifest.py` resolves per sample.
    """

    @property
    def bc1_calls(self) -> Path:
        """Stage 1 — 17-column bc1 emit tables."""
        return self.processed / BC1_SUBDIR

    @property
    def donor_precall(self) -> Path:
        """Stage 2 — 5-column donor editing-outcome tables (GATED, spec §8)."""
        return self.processed / DONOR_PRECALL_SUBDIR

    @property
    def sv_plasmid(self) -> Path:
        """Stage 3 — `.anno` files, cleaned BAMs, plasmid-integration summary."""
        return self.processed / SV_PLASMID_SUBDIR

    @property
    def reconcile(self) -> Path:
        """Stage 4 — 10-column per-cell reconciliation tables."""
        return self.processed / RECONCILE_SUBDIR

    @property
    def joined_outcomes(self) -> Path:
        """Stage 5 — the orchestration-owned Stage 2 x Stage 4 agreement matrix."""
        return self.processed / JOIN_SUBDIR

    @property
    def qc(self) -> Path:
        """QC artifacts and the pre-registered acceptance-test report."""
        return self.processed / QC_SUBDIR

    @property
    def logs(self) -> Path:
        """Per-rule logs (`2>&1 | tee {log}` in the Snakemake rules).

        🔴 Deliberately under `intermediate/`, NOT `processed/`. These logs are written
        by the SNAKEMAKE path (B) only and carry wall-clock text, so a copy inside
        `processed/` would make a literal `diff` of the two paths' processed trees
        differ on every run — reported as a §4 parity violation (fail condition F6)
        when it is nothing of the kind. The spec pre-registers exactly ONE parity
        carve-out (`*.provenance.json` / `PROVENANCE.json`); adding a second to
        accommodate a log location would be relaxing a pre-registered clause. Moving
        the logs out of the compared tree removes the need for a carve-out entirely.
        """
        return self.intermediate / LOGS_SUBDIR

    @property
    def locus_config(self) -> Path:
        """Flat key/value locus sidecars emitted by `core/locus_config.py`.

        TRACE's `load_locus_config` is NOT a YAML parser — it splits each line on the
        first `:` and ignores nesting entirely, so nested YAML is silently mis-read.
        The sidecar written here is flat by construction.
        """
        return self.intermediate / LOCUS_CONFIG_SUBDIR

    @property
    def samplesheets(self) -> Path:
        """TRACE's 7-column donor samplesheets, built by `core/samplesheet.py`."""
        return self.intermediate / SAMPLESHEET_SUBDIR

    @property
    def donor_precall_workdir(self) -> Path:
        """The `workdir` argument `donor_outcome.precall()` resolves vcf/gvcf under."""
        return self.intermediate / DONOR_PRECALL_SUBDIR

    @property
    def precall_vcf(self) -> Path:
        """`{workdir}/vcf/` — per-sample `<sample>.fil.vcf`, plain text."""
        return self.donor_precall_workdir / VCF_SUBDIR

    @property
    def precall_gvcf(self) -> Path:
        """`{workdir}/gvcf/` — per-sample `<sample>.g.vcf`, plain text."""
        return self.donor_precall_workdir / GVCF_SUBDIR


@dataclass
class WGSColumns:
    """The donor-window column vocabulary shared by `donor_outcome` and `reconcile`.

    Both TRACE modules default to these names, and `core/samplesheet.py` /
    `core/reconcile_driver.py` must agree with whatever the curated oligo table uses.
    Declared once here so the two stages cannot drift apart (spec §2.6).
    """

    oligo_col: str = "oligo_name"
    chrom_col: str = "chrom"
    start_col: str = "donor_start_coord"
    end_col: str = "donor_end_coord"


@dataclass
class WGSLocusConfig:
    """Barcode-locus geometry, flattened into the sidecar TRACE can actually parse.

    Loci are `contig:start-end`, 1-based inclusive, matching `TrackConfig.bc_start` /
    `bc_end`. `TrackConfig.validate_against_reference` asserts the locus is an N-run
    flanked by these junction sequences and fails loudly on reference drift — which is
    the whole point of carrying them, so never relax a mismatch into a warning.

    The plasmid fields are optional: a reference without a guide-donor plasmid contig
    yields a genome-only run. Note that TRACE's `tracks_from_config` indexes the plasmid
    keys unconditionally and `bc1_consensus.run` validates the plasmid track even when
    `use_plasmid=False`, so `pipeline/step_1_bc1_call.py` uses the LIBRARY entry
    (`call_sample(..., plasmid_track=None)`) rather than `run()`.
    """

    barcode_locus: str
    barcode_left_sequence: str
    barcode_right_sequence: str
    plasmid_bc_locus: str | None = None
    plasmid_bc_left_sequence: str | None = None
    plasmid_bc_right_sequence: str | None = None
    amplicon_left_window: str | None = None  # "start-end", 1-based inclusive
    amplicon_right_window: str | None = None  # "start-end", 1-based inclusive
    # The reference sequence each window covers. Required with its window: a window is a
    # raw coordinate pair, and `validate_tracks` compares this sequence with the FASTA
    # every run, as it does the junction flanks, so a re-pointed locus cannot leave the
    # windows silently at the old primer sites.
    amplicon_left_window_sequence: str | None = None
    amplicon_right_window_sequence: str | None = None

    @property
    def has_plasmid_track(self) -> bool:
        """True when all three plasmid keys are present and a plasmid track is buildable."""
        return all(
            v is not None
            for v in (
                self.plasmid_bc_locus,
                self.plasmid_bc_left_sequence,
                self.plasmid_bc_right_sequence,
            )
        )


@dataclass
class WGSBc1Config:
    """Stage 1 — de novo bc1 calling (`trace_crispr.wgs.bc1_consensus`)."""

    bam_suffix: str = ".bam"
    """Suffix stripped from a BAM basename to derive the sample name.

    🔴 TRACE defaults to `".sort.dedup.recal2.bam"`, which matches NONE of this
    program's BAMs (`plate_1B_clone_1.bam`). Stage 1 takes its BAM list from the
    manifest, not from a glob, so a wrong value here mis-names samples rather than
    silently emptying the run — but the header-only/short-row-count guard still lives in
    `qc/pipeline_qc.py` (spec F1), because exit code 0 proves nothing about row count.
    """

    use_plasmid: bool = False
    """Call the plasmid barcode track in addition to the genome track."""

    min_mapq: int | None = None
    """TRACE's default is conservative locus assignment (20); explicit 0 enables legacy sensitivity comparisons."""

    bridge_table: Path | None = None
    """Optional bc1 -> oligo bridge table. `None` runs BRIDGE-LESS, whose pre-registered
    consequence is `bridge_tier == "NONE"` on every row and a depressed `confidence`
    column — expected, not a defect. Running bridge-less isolates the caller from the
    bridge, which is what the §7 acceptance test measures."""

    bridge_bc1_col: str = "bc1"
    bridge_oligo_col: str = "oligo_name"
    bridge_count_col: str | None = None
    bridge_table_revcomp: bool = False
    """Set when the bridge table stores bc1 in the opposite orientation to the emit
    table. Orientation is the quantity acceptance criterion P3b pins down; do not guess
    it, record it."""


@dataclass
class WGSDonorPrecallConfig:
    """Stage 2 — donor-window editing outcome (`trace_crispr.wgs.donor_outcome`).

    🔴 GATED THREE WAYS (spec §8) and not part of the pre-registered acceptance test:
      GATE 1 — the SYNTHETIC_ERROR axis cannot run (no per-oligo step-2 donor exists);
      GATE 2 — no gVCF exists anywhere, so `designed_variant_unedited` cannot be
               separated from `unknown`, and `unknown` must NEVER be read as `unedited`;
      GATE 3 — `precall_sample` re-parses its own `chrom_pos_ref_to_alt` id with
               `split("_")`, so any contig name containing `_` raises ValueError. The
               mandated composite reference carries `plasmid_{oligo}` and
               `integrated_plasmid_{oligo}` contigs, so this is reachable here. The fix
               is TRACE-side; `qc/pipeline_qc.py` guards it loudly meanwhile.
    """

    donor_table: Path | None = None
    """Headered oligo donor definitions with explicit coordinates and sequences."""

    ne_threshold: int = 1
    """Minimum reference depth for a gVCF position to count as covered."""

    quantification_window_size: int = 10
    """Width of the coverage window centered on the donor midpoint, in bases."""

    quantification_window_cov_perc: float = 1.0
    """Fraction of the quantification window that must be covered before an unedited
    call is allowed; below it the outcome is `unknown`."""

    synthetic_donor_source: str = SYNTHETIC_DONOR_FROM_DESIGN
    """Where `core/samplesheet.py` gets the `synthetic_donor` column. Only
    `"design_donor"` is available today (GATE 1); any other value names a real per-oligo
    step-2 donor table and enables the axis."""

    ad_format_slot_in: int = 4
    """1-based FORMAT slot holding AD in the SOURCE VCFs. The existing per-plate
    bcftools merges are `GT:PL:DP:AD`, so AD is slot 4."""

    ad_format_slot_out: int = 2
    """1-based FORMAT slot holding AD that TRACE's POSITIONAL gVCF parser assumes
    (`sample_field.split(":")[1]`). `core/vcf_adapter.py` remaps 4 -> 2. This is pure
    re-representation, which is why it is orchestration and not an algorithm."""

    plate_vcfs: list[Path] = field(default_factory=list)
    """Per-plate multi-sample VCFs to split into the per-sample layout. Empty means
    Stage 2 has no input and must fail loudly rather than emit an empty table."""

    min_allele_fraction: float | None = None
    """Rescue a reference-GT record whose AD-derived allele fraction clears this.

    NN 251 defect 1. `None` keeps the GT-only behaviour, which is what every existing
    output was produced with. See `core/vcf_adapter.PLOIDY_CAVEAT` before setting it.
    planning/268 §2e (Kevin 2026-10-24): a rescued record is written with GT `0/1`
    (`vcf_adapter.RESCUED_GT`), with the caller's original kept in `MAG_RESCUED_GT`.
    The pinned TRACE refuses that as `unknown;ambiguous_genotype`, and step_5 reports a
    well whose ambiguous genotype is a `MAG_MINOR` allele as `mixture_signal`. A rescued
    minority allele therefore never reads as unedited, and is never asserted installed.

    ⚠️ NN 255: this governs the RESCUE only. Annotating a minority allele is now
    unconditional, so a 10% `GT 0/1` on a diploid-default plate VCF is flagged whatever
    this is set to — that was NN 252's F11b, a guard that only fired on a non-default
    branch. `None` is still required to reproduce a legacy run byte for byte, but a
    `--ploidy 1` plate VCF with `min_allele_fraction=None` is REFUSED by `step_2`,
    because at `--ploidy 1` every minority allele is `GT 0` and the rescue is the only
    thing that can see it."""

    min_alt_reads: int | None = None
    """Minimum ALT read count for a rescued record. NN 255; made depth-conditional by Kevin's
    2026-10-25 ruling (planning/284 §8 ruling 1, "even a single read is OK"): `None` (the default)
    takes TRACE's `reconcile.min_alt_reads(depth)` at the record's AD depth, which is one read up to
    51x; an integer pins a fixed floor, e.g. `2` to reproduce an NN 255-era run.

    The NN 255 reasoning for a fixed 2, kept for the record:

    A fraction threshold with no count floor selects for low coverage. Measured per-well
    depth, named by dataset: `spg_20240702_POST` median 20.1x / q25 9.1x (n=480);
    the 20231213 plate-2348/2349 cohort median 14.4x / q25 6.2x (n=768). At either,
    `min_allele_fraction=0.05` rounds to a ONE-read minimum whose per-position
    false-positive rate from sequencing error alone is ~1-2% (binomial, q=1e-3..5e-3).
    At two reads that falls below 1e-3 at 40x while sensitivity to a true 10% allele is
    still 0.92. Default matches `vcf_adapter.DEFAULT_MIN_ALT_READS`."""

    @property
    def synthetic_error_axis(self) -> str:
        """`"disabled"` while `synthetic_donor := design_donor`, else `"enabled"`.

        🔴 This one computed value is the single source for the stamp that
        `core/samplesheet.py` writes into the output header, that `qc/provenance.py`
        records in `PROVENANCE.json`, and that `qc/pipeline_qc.py` fails the run for
        missing. Read it; never hardcode the string in a step module.
        """
        if self.synthetic_donor_source == SYNTHETIC_DONOR_FROM_DESIGN:
            return AXIS_DISABLED
        return AXIS_ENABLED


@dataclass
class WGSSvPlasmidConfig:
    """Stage 3 — plasmid integration and donor-read cleaning (`trace_crispr.wgs.sv_plasmid`)."""

    donor_sequences: Path | None = None
    """Headered oligo_name/donor_sequence table for per-contig junction derivation."""

    min_mapq: int | None = None
    """Optional override; the standard-mode evidence default is owned by TRACE."""

    min_anchor: int | None = None
    """Minimum aligned bases per side of a junction; TRACE supplies the default."""

    blacklist: Path | None = None
    """`REDI_reads_mapped_regions.txt` — reference-specific REDI intervals.
    Standard mode validates intervals and reference identity; permissive historical
    parsing is available only in explicit parity mode."""

    blacklist_reference: Path | None = None
    """Reference on which blacklist coordinates were defined; required in standard mode."""

    samtools: str = "samtools"
    """🔴 UNDECLARED EXTERNAL BINARY. `sv_plasmid` shells out to `samtools` three times
    and this cannot be expressed as a Python dependency, so `pipeline/step_3_sv_plasmid.py`
    preflights it on PATH before any work."""

    flexible_bases: int = 6
    """Perl-parity slack when testing donor-window containment."""

    dup_mode: str = "standard"
    """`"standard"` decodes the SAM duplicate flag correctly; `"parity"` reproduces the
    original Perl bug. These genuinely differ (`decode_duplicate(1024)` is True vs
    False), so the choice is explicit rather than defaulted in the algorithm layer."""

    bp_to_extend: int = 1000
    """Flank length used when constructing the integrated-plasmid contig."""

    junction_distance: int | None = None
    """`None` DERIVES the junction per contig via `sv_plasmid.integrated_plasmid_junction`,
    which is what TRACE's own docstring says production should do. Setting an integer
    pins the `1129` parity literal instead and should only be done to reproduce a legacy
    run."""


@dataclass
class WGSReconcileConfig:
    """Stage 4 — per-cell edit reconciliation (`trace_crispr.wgs.reconcile`).

    Carries only the INPUT table; `classify_site`'s `tau_full` / `tau_wt` / `min_cov`
    and `window_edit_vaf`'s `min_base_qual` are deliberately absent — they are the
    algorithm's decision rule and live in TRACE (spec D3).

    `reconcile.run()` is the one TRACE entry point that takes an `argparse.Namespace`
    rather than kwargs, so `core/reconcile_driver.py` composes the primitives directly
    instead of synthesizing a Namespace.
    """

    require_frame_check: bool = True
    """Refuse a window whose `wt_seq` does not match the reference in either orientation.

    NN 251 defect 10 / NN 243 A11: a reference that already carries the design's ALT
    inverts every per-cell call **while the apparent aggregate rate lands on exactly
    truth** (per-cell accuracy 0.0%), so only a per-site check can detect it.

    Since planning/268 a window that declares `intended_variants` is scored by TRACE's
    allele-specific scorer, which validates each REF against the FASTA, so a wrong frame
    cannot invert it. A window that declares none is scored from read evidence against
    `fasta.fetch` (`core.edit_profile`), where it can. So `True` refuses a table with NO
    `wt_seq` column only when some window declares no intended allele, and otherwise
    records `frame_unverified_no_wt_seq` in the audit sidecar. When `wt_seq` IS present
    every window is verified with `core.window_table.verify_window_frame`; `True` refuses
    a mismatch or an empty `wt_seq`. `False` records each case as accepted-unverified,
    never a silent pass."""

    strict_window_table: bool = False
    """Also refuse a window table that repeats an oligo row.

    NN 251 defect 12. TRACE's `load_oligo_windows` now refuses a malformed row (a
    `3001.0`, an empty coordinate) and a conflicting duplicate itself, so Stage 4 cannot
    drop a window in silence at any setting; `core.window_table.load_windows_audited`
    runs first so the denominators are logged before such a refusal. What `True` adds is
    refusing IDENTICAL duplicate rows, which TRACE accepts. The real design table is
    only **50.2% filled** on `donor_start_coord`/`donor_end_coord`, so it must be
    curated before Stage 4 will accept it."""

    min_mapq: int | None = None
    """Optional intended-allele mapping-quality floor; TRACE owns default 20."""

    read_source: str = "auto"
    """Allele evidence: auto uses Stage 3 cleaning when its inputs are configured;
    raw preserves original alignments; donor_cleaned requires completed Half B.
    Barcode identity always uses the original BAM. Raw evidence is labelled and
    observed plasmid barcode evidence caps raw editing confidence at low.
    """

    oligo_windows: Path | None = None
    """Oligo-window table (`oligo_name`, `chrom`, `donor_start_coord`,
    `donor_end_coord`, optional `intended_variants` JSON `[{pos, ref, alt}]` in the BAM
    reference frame, optional `wt_seq` / `design_donor`). TRACE's `load_oligo_windows`
    refuses a malformed or conflicting row. A window with `intended_variants` is scored
    by TRACE per allele; one without is scored from its read evidence, gated at the
    designed base that `wt_seq` + `design_donor` give, and flagged
    `WINDOW_EVIDENCE_NO_INTENDED_ALLELES` (planning/268). A window with neither is
    REFUSED (`DesignedPositionsRequired`; Kevin 2026-10-25, planning/284 §8 ruling 4).

    The data that enables that (planning/268 §2g, agreed by Kevin 2026-10-24): the
    cohort's Stage-2 donor table (`donor_precall.donor_table`) already carries `chrom`,
    the donor coordinates, `wt_seq` and `design_donor` for every assigned oligo. Point
    this key at that table, or at a copy with `intended_variants` added. That turns on
    the frame check and the designed-base gate. The raw design table will not load: TRACE
    refuses its unfilled coordinate rows."""

    oligo_assignments: Path | None = None
    """Per-sample oligo assignment (`sample` -> one or more `oligo_name`), the input that
    turns an oligo-window table into this cohort's edit sites."""


@dataclass
class WGSConfig:
    """Top-level config for one MAGESTIC WGS cohort run.

    Attribute access paths the rest of the module depends on:

        config.paths.<derived dir>          filesystem layout (see `WGSPaths`)
        config.paths.manifest               per-sample BAM + reference (spec D1)
        config.columns.<col>_col            shared donor-window column vocabulary
        config.locus.<key>                  barcode-locus geometry
        config.bc1.<param>                  Stage 1
        config.donor_precall.<param>        Stage 2 (GATED)
        config.sv_plasmid.<param>           Stage 3
        config.reconcile.<param>            Stage 4
        config.resources_for("step_1_bc1")  cluster resources for one rule
        config.get_working_dir(BC1_SUBDIR)  scratch-aware write location
    """

    cohort_name: str
    sequencing_date: str
    paths: WGSPaths
    columns: WGSColumns = field(default_factory=WGSColumns)
    locus: WGSLocusConfig | None = None
    bc1: WGSBc1Config = field(default_factory=WGSBc1Config)
    donor_precall: WGSDonorPrecallConfig = field(default_factory=WGSDonorPrecallConfig)
    sv_plasmid: WGSSvPlasmidConfig = field(default_factory=WGSSvPlasmidConfig)
    reconcile: WGSReconcileConfig = field(default_factory=WGSReconcileConfig)
    resources: dict[str, dict[str, Any]] = field(default_factory=dict)
    python: str | None = None
    """Optional interpreter override for the Snakemake path; `None` uses `sys.executable`."""

    # Scratch I/O (the §9 standard `--use-scratch` flag)
    use_scratch: bool = False
    scratch_dir: Path | None = None

    def __post_init__(self) -> None:
        """Fill the scratch working directory when `--use-scratch` was requested."""
        if self.use_scratch and self.scratch_dir is None:
            self.scratch_dir = get_scratch_workdir(self.cohort_name)

    def get_working_dir(self, subdir: str) -> Path:
        """Return the write location for one stage's output subdirectory.

        With `use_scratch` this redirects under `scratch_dir`; otherwise it composes
        `paths.processed` with the same module-level constant the matching `WGSPaths`
        property uses, so the two agree by construction. `scratch_dir` alone changes
        nothing — it is inert unless `use_scratch` is also true.

        Args:
            subdir: One of the `*_SUBDIR` constants in this module.

        Returns:
            Path to the working directory (not created — see `ensure_directories`).
        """
        if self.use_scratch and self.scratch_dir is not None:
            return self.scratch_dir / subdir
        return self.paths.processed / subdir

    def resources_for(self, step_key: str) -> dict[str, Any]:
        """Return cluster resources for one stage, falling back to `DEFAULT_RESOURCES`.

        Args:
            step_key: One of `STEP_RESOURCE_KEYS`.

        Returns:
            Mapping with at least `mem_mb`, `time` and `cpus`.

        Raises:
            KeyError: If `step_key` is not a recognized stage.
        """
        if step_key not in STEP_RESOURCE_KEYS:
            raise KeyError(f"unknown step key {step_key!r}; expected one of {STEP_RESOURCE_KEYS}")
        merged = dict(DEFAULT_RESOURCES[step_key])
        merged.update(self.resources.get(step_key, {}))
        return merged

    def validate(self, steps: Sequence[str] | None = None) -> list[str]:
        """Validate the configuration against the filesystem.

        Parameters
        ----------
        steps : sequence of str, optional
            The stage ids actually about to run (`"1"`..`"5"`). Per-stage INPUT checks
            are applied only to these; global checks (project root, manifest, locus,
            resource keys) always run. `None` means "validate everything", which is the
            right behavior for a full run and for a config linter.

            🔴 **Scoping is required, not a convenience.** The pre-registered acceptance
            run is `magestic-wgs --config <yaml> --steps 1`, and this cohort has no
            Stage-3 blacklist and no Stage-4 oligo tables — those inputs are GATED and do
            not exist. Validating them unconditionally makes Stage 1 unrunnable on the
            only cohort it has, i.e. a config error is reported for a stage nobody asked
            to run. Note the failure mode this avoids is a FALSE failure: the unscoped
            check is loud, not silent, but it is loud about the wrong thing.

        Returns
        -------
        list of str
            Validation error messages; empty means valid. Callers MUST exit non-zero on
            a non-empty list (ARCHITECTURE_POLICY §7 — no silent success).
        """
        wanted = set(steps) if steps is not None else set(STEP_RESOURCE_KEYS_BY_ID)
        errors: list[str] = []

        if not self.paths.project_root.exists():
            errors.append(f"Project root does not exist: {self.paths.project_root}")

        # D1: without the manifest there is no per-sample reference, so nothing may run.
        if not self.paths.manifest.exists():
            errors.append(f"Sample manifest does not exist: {self.paths.manifest}")

        if self.paths.reference_root is not None and not self.paths.reference_root.exists():
            errors.append(f"Reference root does not exist: {self.paths.reference_root}")

        if self.locus is None:
            errors.append("No `locus:` block: Stage 1 cannot build a barcode track without one")

        if self.bc1.use_plasmid and (self.locus is None or not self.locus.has_plasmid_track):
            errors.append(
                "bc1.use_plasmid is true but locus.plasmid_bc_locus / "
                "plasmid_bc_left_sequence / plasmid_bc_right_sequence are not all set"
            )

        # Per-stage INPUT checks, applied only to the stages actually requested.
        # Stage 4 consumes the same oligo tables as Stage 3's Half B, so it is listed
        # against both; a path configured for neither stage is never checked.
        for label, path, owning_steps in (
            ("bc1.bridge_table", self.bc1.bridge_table, {"1", "4"}),
            ("sv_plasmid.blacklist", self.sv_plasmid.blacklist, {"3"}),
            ("sv_plasmid.blacklist_reference", self.sv_plasmid.blacklist_reference, {"3"}),
            ("sv_plasmid.donor_sequences", self.sv_plasmid.donor_sequences, {"3"}),
            ("donor_precall.donor_table", self.donor_precall.donor_table, {"2"}),
            ("reconcile.oligo_windows", self.reconcile.oligo_windows, {"3", "4"}),
            ("reconcile.oligo_assignments", self.reconcile.oligo_assignments, {"2", "3", "4"}),
        ):
            if not (owning_steps & wanted):
                continue
            if path is not None and not path.exists():
                errors.append(f"{label} does not exist: {path}")

        if "2" in wanted:
            for vcf in self.donor_precall.plate_vcfs:
                if not vcf.exists():
                    errors.append(f"donor_precall.plate_vcfs entry does not exist: {vcf}")

        for step_key in self.resources:
            if step_key not in STEP_RESOURCE_KEYS:
                errors.append(
                    f"resources key {step_key!r} is not a stage; expected one of "
                    f"{STEP_RESOURCE_KEYS}"
                )

        return errors

    def ensure_directories(self) -> None:
        """Create every output directory this run will write to."""
        for subdir in (
            BC1_SUBDIR,
            DONOR_PRECALL_SUBDIR,
            SV_PLASMID_SUBDIR,
            RECONCILE_SUBDIR,
            JOIN_SUBDIR,
            TABLES_SUBDIR,
            QC_SUBDIR,
            LOGS_SUBDIR,
        ):
            self.get_working_dir(subdir).mkdir(parents=True, exist_ok=True)
        for path in (
            self.paths.locus_config,
            self.paths.samplesheets,
            self.paths.precall_vcf,
            self.paths.precall_gvcf,
        ):
            path.mkdir(parents=True, exist_ok=True)


def _as_path(value: Any) -> Path | None:
    """Coerce an optional YAML scalar to a `Path`, preserving `None`/empty as `None`."""
    return Path(value).expanduser() if value else None


def _under_root_optional(value: Any, project_root: Path) -> Path | None:
    """Resolve an optional path against the project root when it is relative.

    A relative path left CWD-relative resolves DIFFERENTLY under the direct path (run
    from anywhere) and the Snakemake path (which sets its own working directory), which
    is a §4 parity break that surfaces as "file not found" rather than as a path bug.
    """
    path = _as_path(value)
    if path is None or path.is_absolute():
        return path
    return project_root / path


def load_config(
    path: str | Path,
    *,
    project_dir: str | Path | None = None,
    use_scratch: bool = False,
) -> WGSConfig:
    """Parse a per-cohort YAML config into a `WGSConfig`.

    Parameters
    ----------
    path
        Path to the cohort YAML (`workflow/config/screens/<cohort>_wgs.yaml`).
    project_dir
        Optional override for `paths.project_root`, wired to the §9 standard
        `--project-dir` flag. Relative sub-paths in the YAML are resolved against the
        effective project root, so the override moves the whole tree at once.
    use_scratch
        Wired to the §9 standard `--use-scratch` flag; stages I/O under `$SCRATCH`.

    Returns
    -------
    WGSConfig

    Raises
    ------
    ValueError
        If the document declares a cohort-level `genome:` / `reference:` / `bam_dir:`
        key (spec D1), or carries an unrecognized top-level key. Both raise rather than
        being ignored: a silently dropped key is a wrong answer waiting to be quoted.
    KeyError
        If a required key (`cohort_name`, `sequencing_date`, `paths.project_root`,
        `paths.manifest`) is missing.
    """
    path = Path(path).expanduser().resolve()
    with path.open() as fh:
        data: dict[str, Any] = yaml.safe_load(fh) or {}

    for key, reason in _FORBIDDEN_TOP_LEVEL_KEYS.items():
        if key in data:
            raise ValueError(f"{path}: `{key}:` is not a valid WGS config key — {reason}")

    unknown = sorted(set(data) - _KNOWN_TOP_LEVEL_KEYS)
    if unknown:
        raise ValueError(
            f"{path}: unrecognized top-level key(s) {unknown}; expected a subset of "
            f"{sorted(_KNOWN_TOP_LEVEL_KEYS)}"
        )

    paths_d: dict[str, Any] = data["paths"]
    if project_dir is not None:
        project_root = Path(project_dir).expanduser().resolve()
    else:
        project_root = Path(paths_d["project_root"]).expanduser()
        if not project_root.is_absolute():
            project_root = path.parent / project_root
        project_root = project_root.resolve()

    def _under_root(key: str, default: str) -> Path:
        """Resolve a path key against the effective project root when it is relative."""
        raw = Path(paths_d.get(key, default)).expanduser()
        return raw if raw.is_absolute() else project_root / raw

    paths = WGSPaths(
        project_root=project_root,
        manifest=_under_root("manifest", "keyfiles/sample_manifest.tsv"),
        intermediate=_under_root("intermediate", "intermediate_data"),
        processed=_under_root("processed", "processed_data"),
        keyfiles=_under_root("keyfiles", "keyfiles"),
        reference_root=_under_root_optional(paths_d.get("reference_root"), project_root),
    )

    columns = WGSColumns(**data.get("columns", {}))

    locus_d = data.get("locus")
    locus = WGSLocusConfig(**locus_d) if locus_d else None

    bc1_d = dict(data.get("bc1", {}))
    bc1_d["bridge_table"] = _under_root_optional(bc1_d.get("bridge_table"), project_root)
    value = bc1_d.get('min_mapq')
    if value is not None and (type(value) is not int or not 0 <= value <= 254):
        raise ValueError('bc1.min_mapq must be an integer from 0 to 254')
    bc1 = WGSBc1Config(**bc1_d)

    precall_d = dict(data.get("donor_precall", {}))
    precall_d["plate_vcfs"] = [_under_root_optional(p, project_root) for p in precall_d.get("plate_vcfs", [])]
    precall_d["donor_table"] = _under_root_optional(precall_d.get("donor_table"), project_root)
    donor_precall = WGSDonorPrecallConfig(**precall_d)

    sv_d = dict(data.get("sv_plasmid", {}))
    sv_d["blacklist"] = _under_root_optional(sv_d.get("blacklist"), project_root)
    sv_d["blacklist_reference"] = _under_root_optional(sv_d.get("blacklist_reference"), project_root)
    sv_d["donor_sequences"] = _under_root_optional(sv_d.get("donor_sequences"), project_root)
    if 'samtools' in sv_d and (not isinstance(sv_d['samtools'], str) or not sv_d['samtools'].strip()):
        raise ValueError('sv_plasmid.samtools must be a nonempty executable string')
    if sv_d.get('dup_mode', 'standard') not in ('standard', 'parity'):
        raise ValueError('sv_plasmid.dup_mode must be standard or parity')
    for key, lower in (('min_mapq', 0), ('min_anchor', 1)):
        value = sv_d.get(key)
        if value is not None and (type(value) is not int or value < lower):
            raise ValueError(f'sv_plasmid.{key} must be an integer >= {lower}')
    if sv_d.get('min_mapq') is not None and sv_d['min_mapq'] > 254:
        raise ValueError('sv_plasmid.min_mapq must be <= 254 (255 is unknown)')
    for key, lower in (('bp_to_extend', 1), ('flexible_bases', 0), ('junction_distance', 1)):
        value = sv_d.get(key)
        if value is not None and (type(value) is not int or value < lower):
            raise ValueError(f'sv_plasmid.{key} must be an integer >= {lower}')
    sv_plasmid = WGSSvPlasmidConfig(**sv_d)

    recon_d = dict(data.get("reconcile", {}))
    if recon_d.get('read_source', 'auto') not in ('auto', 'raw', 'donor_cleaned'):
        raise ValueError('reconcile.read_source must be auto, raw, or donor_cleaned')
    value = recon_d.get('min_mapq')
    if value is not None and (type(value) is not int or not 0 <= value <= 254):
        raise ValueError('reconcile.min_mapq must be an integer from 0 to 254')
    recon_d["oligo_windows"] = _under_root_optional(recon_d.get("oligo_windows"), project_root)
    recon_d["oligo_assignments"] = _under_root_optional(recon_d.get("oligo_assignments"), project_root)
    reconcile = WGSReconcileConfig(**recon_d)

    return WGSConfig(
        cohort_name=data["cohort_name"],
        sequencing_date=str(data["sequencing_date"]),
        paths=paths,
        columns=columns,
        locus=locus,
        bc1=bc1,
        donor_precall=donor_precall,
        sv_plasmid=sv_plasmid,
        reconcile=reconcile,
        resources=data.get("resources", {}),
        python=data.get("python"),
        use_scratch=use_scratch,
    )
