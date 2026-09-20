"""
magestic.pipelines.wgs.core.reconcile_driver — compose TRACE's reconcile primitives.

Stage 4 is the only TRACE entry point that is not library-callable with keyword
arguments: `reconcile.run(args: argparse.Namespace)` takes a Namespace and discovers its
BAMs with a bare `glob`. Synthesizing a Namespace to reach it would smuggle a CLI
boundary into the middle of the orchestration layer and would re-introduce the
silent-empty glob (a `bam_suffix` that matches nothing writes a header-only TSV and
returns 0). So this module composes the same primitives directly, from the manifest.

🔴 **It re-implements no algorithm.** `classify_site`, `window_edit_vaf`,
`build_edit_sites` and `reconcile_cell` are called, never copied, and every threshold
they own is left at TRACE's default unless the screen config explicitly overrides it —
hence the `None` defaults below, which mean "TRACE decides" rather than a number
restated in the wrong layer.

The one exception is documented in `bc_confidence()`, and it is a TRACE-side gap.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING, Any

from . import _trace

if TYPE_CHECKING:  # pragma: no cover - typing only
    import pysam
    from trace_crispr.wgs.bc1_consensus import BarcodeCall, BridgeIndex, TrackConfig
    from trace_crispr.wgs.reconcile import CellReconciliation, SiteCall

__all__ = [
    "BC_CONFIDENCE_HIGH_TIERS",
    "bc_confidence",
    "load_oligo_windows",
    "reconcile_sample",
    "score_sites",
]

#: Bridge/caller confidence tiers that make a cell's bc1 evidence count as "high".
#: Transcribed verbatim from `trace_crispr.wgs.reconcile.run` — see `bc_confidence`.
BC_CONFIDENCE_HIGH_TIERS: tuple[str, ...] = ("high", "medium")


def _only_set(**kwargs: Any) -> dict[str, Any]:
    """Drop `None` values, so TRACE's own defaults apply unless a caller overrode them."""
    return {k: v for k, v in kwargs.items() if v is not None}


def bc_confidence(genome_calls: Sequence[BarcodeCall]) -> str:
    """Collapse a cell's genome-track barcode calls to the `bc_confidence` label.

    🔴 **This is the one classification rule magestic owns that it should not.** It is
    transcribed verbatim from the body of `trace_crispr.wgs.reconcile.run`
    (`"high" if (gbcs and gbcs[0].confidence in ("high", "medium")) else "low"`), because
    TRACE exposes the rule only inside that argparse entry point and not as a function.
    It is isolated here, with the tier tuple as a module constant, so it is greppable and
    so replacing it is a one-line change once TRACE grows the primitive
    (`reconcile.bc_confidence_from(...)` is the upstream fix).

    It matters more than its size suggests: `reconcile_cell` defaults
    `bc_confidence="high"`, so a caller that omits it does not get a neutral default —
    it gets a wrong one that inflates confidence.
    """
    if genome_calls and genome_calls[0].confidence in BC_CONFIDENCE_HIGH_TIERS:
        return "high"
    return "low"


def load_oligo_windows(path: str, **column_names: str) -> dict[str, Any]:
    """Load the `oligo_name -> donor window` map with TRACE's own loader.

    Parameters
    ----------
    path : str
        Header'd TSV using the shared `{oligo_name, chrom, donor_start_coord,
        donor_end_coord}` column vocabulary. Pass `oligo_col=`, `chrom_col=`,
        `start_col=`, `end_col=` to override.

    Returns
    -------
    dict
        **Opaque.** Hand it straight to `score_sites` / `build_edit_sites` and never
        unpack a value: TRACE's element type is version-dependent, and this module's
        contract is that the two functions that produce and consume it always agree
        because they come from the same TRACE.

    Notes
    -----
    TRACE's loader silently skips rows that are too short or whose coordinates are not
    `.isdigit()`, and keeps the FIRST occurrence per oligo. Neither is reported, so a
    caller that needs a denominator must compare `len(...)` against the source row count.
    """
    return _trace.reconcile().load_oligo_windows(path, **column_names)


def score_sites(
    bam: pysam.AlignmentFile,
    fasta: pysam.FastaFile,
    oligos: Iterable[str],
    oligo_windows: dict[str, Any],
    *,
    tau_full: float | None = None,
    tau_wt: float | None = None,
    min_cov: int | None = None,
    min_base_qual: int | None = None,
) -> list[SiteCall]:
    """Measure and classify the edit VAF at each designed window this cell's oligos predict.

    Thin composition of `build_edit_sites` -> `window_edit_vaf` -> `classify_site`. The
    three threshold arguments default to `None`, meaning TRACE's defaults stand; pass a
    value only when the screen config deliberately overrides one.

    Coordinates are 1-based inclusive throughout (`EditSite.start`/`end` and
    `window_edit_vaf`), matching the samplesheet vocabulary — there is no conversion here.
    """
    tr = _trace.reconcile()
    site_kwargs = _only_set(tau_full=tau_full, tau_wt=tau_wt, min_cov=min_cov)
    vaf_kwargs = _only_set(min_base_qual=min_base_qual)
    calls: list[SiteCall] = []
    for site in tr.build_edit_sites(list(oligos), oligo_windows):
        vaf, cov = tr.window_edit_vaf(bam, fasta, site.chrom, site.start, site.end, **vaf_kwargs)
        calls.append(tr.SiteCall(site, vaf, cov, tr.classify_site(vaf, cov, **site_kwargs)))
    return calls


def reconcile_sample(
    bam: pysam.AlignmentFile,
    fasta: pysam.FastaFile,
    sample: str,
    *,
    genome_track: TrackConfig,
    oligo_windows: dict[str, Any],
    plasmid_track: TrackConfig | None = None,
    bridge: BridgeIndex | None = None,
    min_support: int | None = None,
    high_only: bool = True,
    tau_full: float | None = None,
    tau_wt: float | None = None,
    min_cov: int | None = None,
    min_base_qual: int | None = None,
) -> CellReconciliation:
    """Reconcile one cell: bc1 lineage set x per-site edit VAF -> a `CellReconciliation`.

    Parameters
    ----------
    bam, fasta : pysam handles
        Open, indexed. Both `extract_track` (via `bam.fetch`) and `window_edit_vaf` (via
        `bam.pileup`) require the index; TRACE documents neither requirement. The caller
        owns these handles, so a step can reuse one BAM across stages.
    sample : str
        The manifest sample name — NOT a BAM basename stripped of a suffix. Passing it in
        is what keeps the fan-out keyed to the manifest and the glob unreachable.
    genome_track, plasmid_track : TrackConfig
        From `core.locus_config.build_tracks`. `plasmid_track=None` runs genome-only and
        yields `n_plasmid_bc=0`.
    oligo_windows : dict
        From `load_oligo_windows`; opaque, passed straight through.
    bridge : BridgeIndex, optional
        bc1 -> oligo index. **Without it no oligo is attributed, so no site is scored and
        the cell classifies from zero sites** — a legitimate bridge-less mode, but never a
        "no edits found" result.
    min_support, high_only : int, bool
        Passed to `confident_barcodes`. `high_only=True` reproduces the choice
        `reconcile.run` makes at its own entry point; it is a pipeline parameter, so a
        screen config may override it.

    Returns
    -------
    CellReconciliation
        TRACE's own dataclass, unmodified. Serialize it with
        `core.emit.write_reconciliations`, which takes the column list from TRACE too.
    """
    bc1 = _trace.bc1_consensus()
    tr = _trace.reconcile()
    barcode_kwargs = _only_set(min_support=min_support)

    genome_calls = bc1.confident_barcodes(
        bc1.extract_track(bam, genome_track), bridge, high_only=high_only, **barcode_kwargs
    )
    plasmid_calls: list[BarcodeCall] = []
    if plasmid_track is not None:
        plasmid_calls = bc1.confident_barcodes(
            bc1.extract_track(bam, plasmid_track), bridge, high_only=high_only, **barcode_kwargs
        )

    oligos = [b.oligo for b in (list(genome_calls) + list(plasmid_calls)) if b.oligo]
    sites = score_sites(
        bam,
        fasta,
        oligos,
        oligo_windows,
        tau_full=tau_full,
        tau_wt=tau_wt,
        min_cov=min_cov,
        min_base_qual=min_base_qual,
    )
    return tr.reconcile_cell(
        sample,
        sites,
        n_genome_bc=len(genome_calls),
        n_plasmid_bc=len(plasmid_calls),
        bc_confidence=bc_confidence(genome_calls),
    )
