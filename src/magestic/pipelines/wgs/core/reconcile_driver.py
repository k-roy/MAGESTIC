"""Compose TRACE's allele-specific reconciliation without a glob/CLI boundary."""
from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING, Any
from . import _trace

if TYPE_CHECKING:
    import pysam
    from trace_crispr.wgs.bc1_consensus import BarcodeCall, BridgeIndex, TrackConfig
    from trace_crispr.wgs.reconcile import CellReconciliation, SiteCall

BC_CONFIDENCE_HIGH_TIERS = ("high",)


def _only_set(**kwargs: Any) -> dict[str, Any]:
    return {k: v for k, v in kwargs.items() if v is not None}


def bc_confidence(genome_calls: Sequence[BarcodeCall]) -> str:
    """Use barcode observation confidence, independently of bridge annotation."""
    return _trace.reconcile().bc_confidence_from(genome_calls)


def load_oligo_windows(path: str, **column_names: str) -> dict[str, Any]:
    """Load validated windows and optional intended_variants with TRACE's parser."""
    return _trace.reconcile().load_oligo_windows(path, **column_names)


def score_sites(bam: pysam.AlignmentFile, fasta: pysam.FastaFile,
                oligos: Iterable[str], oligo_windows: dict[str, Any], *,
                tau_full: float | None = None, tau_wt: float | None = None,
                min_cov: int | None = None, min_base_qual: int | None = None, min_mapq: int | None = None) -> list[SiteCall]:
    """Delegate allele counts and site classification entirely to TRACE.

    Missing intended variants yield unassessed window compatibility evidence.
    Missing candidate windows raise rather than silently changing the denominator.
    """
    tr = _trace.reconcile()
    kwargs = _only_set(tau_full=tau_full, tau_wt=tau_wt, min_cov=min_cov,
                       min_base_qual=min_base_qual, min_mapq=min_mapq)
    return [tr.score_edit_site(bam, fasta, site, **kwargs)
            for site in tr.build_edit_sites(list(oligos), oligo_windows, fasta=fasta)]


def reconcile_sample(bam: pysam.AlignmentFile, fasta: pysam.FastaFile, sample: str, *,
                     genome_track: TrackConfig, oligo_windows: dict[str, Any],
                     plasmid_track: TrackConfig | None = None,
                     bridge: BridgeIndex | None = None,
                     assigned_oligos: Sequence[str] | None = None,
                     allele_bam: pysam.AlignmentFile | None = None,
                     read_source_flags: Sequence[str] = (),
                     min_support: int | None = None, high_only: bool = True,
                     tau_full: float | None = None, tau_wt: float | None = None,
                     min_cov: int | None = None, min_base_qual: int | None = None, min_mapq: int | None = None) -> CellReconciliation:
    """Score the explicit assignment, or bridge candidates when none is supplied.

    Explicit assignments define the intended sites even with no bridge. Barcode
    observations still supply mixture evidence. Missing bridge identities cannot
    establish absence of additional edits and therefore lower confidence.
    """
    bc1 = _trace.bc1_consensus()
    tr = _trace.reconcile()
    kwargs = _only_set(min_support=min_support)
    genome_calls = bc1.confident_barcodes(bc1.extract_track(bam, genome_track), bridge,
                                         high_only=high_only, **kwargs)
    plasmid_calls = []
    if plasmid_track is not None:
        plasmid_calls = bc1.confident_barcodes(bc1.extract_track(bam, plasmid_track), bridge,
                                              high_only=high_only, **kwargs)
    calls = list(genome_calls) + list(plasmid_calls)
    bridged = {b.oligo for b in calls if b.oligo}
    flags = list(read_source_flags)
    if assigned_oligos is not None:
        if not assigned_oligos:
            raise ValueError(f"{sample}: explicit assignment is empty")
        oligos = list(dict.fromkeys(assigned_oligos))
        flags.append("EXPLICIT_ASSIGNMENT")
        if bridged - set(oligos):
            flags.append("BRIDGE_ASSIGNMENT_CONFLICT")
    else:
        oligos = sorted(bridged)
        if any(not b.oligo for b in calls):
            flags.append("UNRESOLVED_BC1")
    if not genome_calls:
        flags.append("NO_GENOME_BC1")
    sites = score_sites(allele_bam if allele_bam is not None else bam, fasta, oligos, oligo_windows, tau_full=tau_full,
                         tau_wt=tau_wt, min_cov=min_cov, min_base_qual=min_base_qual, min_mapq=min_mapq)
    confidence = bc_confidence(genome_calls)
    if "BRIDGE_ASSIGNMENT_CONFLICT" in flags:
        confidence = "low"
    if "RAW_DONOR_READS" in flags and plasmid_calls:
        flags.append("RAW_DONOR_ORIGIN_UNRESOLVED")
        confidence = "low"
    return tr.reconcile_cell(sample, sites, n_genome_bc=len(genome_calls),
                             n_plasmid_bc=len(plasmid_calls), bc_confidence=confidence,
                             flags=flags)
