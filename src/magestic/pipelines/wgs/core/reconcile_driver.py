"""Compose TRACE's allele-specific reconciliation without a glob/CLI boundary."""
from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any
from . import _trace, edit_profile

if TYPE_CHECKING:
    import pysam
    from trace_crispr.wgs.bc1_consensus import BarcodeCall, BridgeIndex, TrackConfig
    from trace_crispr.wgs.reconcile import CellReconciliation, SiteCall

BC_CONFIDENCE_HIGH_TIERS = ("high",)

#: Raised on every site scored from window evidence because its oligo declares no
#: intended allele. Provenance, not a penalty: it does not lower confidence.
FLAG_WINDOW_EVIDENCE = "WINDOW_EVIDENCE_NO_INTENDED_ALLELES"


class DesignedPositionsRequired(ValueError):
    """A window with no declared intended allele reached the read-evidence route without its
    designed base(s). Kevin 2026-10-25 (planning/284 §8 ruling 4): REFUSE, never score at the
    window maximum, where a background SNP becomes `SINGLE_EDIT/high` (the W07 regression)."""


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
                min_cov: int | None = None, min_base_qual: int | None = None, min_mapq: int | None = None,
                designed_positions: Mapping[str, Sequence[int]] | None = None,
                allele_evidence: str | None = None) -> list[SiteCall]:
    """Score each candidate window: by its declared alleles when it has them, else by
    its own read evidence.

    A window that declares `intended_variants` is scored by TRACE's allele-specific
    `score_edit_site`: fragment counts at each allele, REF validated against the FASTA.
    🔴 A window without them is NOT left `unassessed` (Kevin, 2026-10-24: "we should
    not require an intended edit in order to be able to assign confidence to variant
    calls from the WGS"; planning/268). It is classified from its read evidence by
    `core.edit_profile` at its designed base(s), which `designed_positions` supplies (from
    `wt_seq` vs `design_donor`, frame-verified by step_4). Without them it REFUSES
    (`DesignedPositionsRequired`; Kevin 2026-10-25, ruling 4), never reading the window
    maximum (NN 251 defects 6 and 11, the W07 regression). The site's `trace_status`
    records what TRACE would have said (`unassessed`), and `FLAG_WINDOW_EVIDENCE` names
    the route.

    Both routes classify from read COUNTS (ruling 1): presence from a single read up,
    absence (`wt`) only at `min_cov` reads; a `partial` resting on one alt read is flagged
    `PARTIAL_ON_A_SINGLE_ALT_READ`. Both pass the interim SV gate (ruling 3): a one-sided
    loss of signal in the window makes the site `nocov`, flagged
    `SV_SUSPECT_ONE_SIDED_COVERAGE_LOSS`.

    Missing candidate windows raise rather than silently changing the denominator.
    """
    tr = _trace.reconcile()
    kwargs = _only_set(tau_full=tau_full, tau_wt=tau_wt, min_cov=min_cov,
                       min_base_qual=min_base_qual, min_mapq=min_mapq)
    # planning/328: how TRACE reads an intended allele. None keeps TRACE's default (local realignment to the
    # site haplotypes since TRACE 6f5e439); "aligned" reproduces the aligned-base reader for parity.
    site_kwargs = dict(kwargs, **_only_set(evidence=allele_evidence))
    profile_kwargs = _only_set(min_base_qual=min_base_qual, min_mapq=min_mapq, min_cov=min_cov)
    calls = []
    for site in tr.build_edit_sites(list(oligos), oligo_windows, fasta=fasta):
        profile = edit_profile.window_edit_profile(
            bam, fasta, site.chrom, site.start, site.end, **profile_kwargs)
        if site.intended_alleles:
            call = tr.score_edit_site(bam, fasta, site, **site_kwargs)
            if call.status == tr.SITE_PARTIAL and any(e.alt_count == 1 for e in call.allele_evidence):
                call = dataclasses.replace(call, flags=tuple(call.flags) + (edit_profile.FLAG_SINGLE_ALT_READ,))
        else:
            dpos = tuple((designed_positions or {}).get(site.oligo or "", ()) or ())
            if not dpos:
                raise DesignedPositionsRequired(
                    f"{site.oligo or '?'} @ {site.chrom}:{site.start}-{site.end} declares no intended allele "
                    "and has no designed position (no frame-verified wt_seq + design_donor), so the window "
                    "cannot be scored. Declare `intended_variants`, or give the window table wt_seq and "
                    "design_donor; it is never scored at the window maximum (planning/268 §4).")
            wcall = edit_profile.classify_window(
                profile, designed_positions=dpos,
                **_only_set(tau_full=tau_full, tau_wt=tau_wt, min_cov=min_cov),
            )
            call = edit_profile.ScoredSite(
                site=site, vaf=profile.max_vaf, cov=profile.cov_at_max, status=wcall.status,
                reason=wcall.reason, variant_class=wcall.variant_class,
                min_depth=profile.min_depth, n_positions_uncovered=profile.n_positions_uncovered,
                depth_at_designed=wcall.depth_at_designed, trace_status=tr.SITE_UNASSESSED,
                flags=tuple(wcall.flags) + (FLAG_WINDOW_EVIDENCE,),
            )
        # 🔴 Interim SV gate on BOTH routes (planning/284 §8 ruling 3): a one-sided loss of signal
        # in the donor window excludes the site; it is never read as an edit or as NO_EDIT.
        side = edit_profile.one_sided_coverage_loss(
            bam, profile, **_only_set(min_mapq=min_mapq))
        if side is not None:
            extra = {"reason": f"sv_suspect_one_sided_loss_{side}"} if isinstance(call, edit_profile.ScoredSite) else {}
            call = dataclasses.replace(call, status=tr.SITE_NOCOV,
                                       flags=tuple(call.flags) + (edit_profile.FLAG_SV_ONE_SIDED_LOSS,), **extra)
        calls.append(call)
    return calls


def reconcile_sample(bam: pysam.AlignmentFile, fasta: pysam.FastaFile, sample: str, *,
                     genome_track: TrackConfig, oligo_windows: dict[str, Any],
                     plasmid_track: TrackConfig | None = None,
                     bridge: BridgeIndex | None = None,
                     assigned_oligos: Sequence[str] | None = None,
                     allele_bam: pysam.AlignmentFile | None = None,
                     read_source_flags: Sequence[str] = (),
                     min_support: int | None = None, high_only: bool = True,
                     tau_full: float | None = None, tau_wt: float | None = None,
                     min_cov: int | None = None, min_base_qual: int | None = None, min_mapq: int | None = None,
                     designed_positions: Mapping[str, Sequence[int]] | None = None,
                     allele_evidence: str | None = None) -> CellReconciliation:
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
                         tau_wt=tau_wt, min_cov=min_cov, min_base_qual=min_base_qual, min_mapq=min_mapq,
                         designed_positions=designed_positions, allele_evidence=allele_evidence)
    confidence = bc_confidence(genome_calls)
    if "BRIDGE_ASSIGNMENT_CONFLICT" in flags:
        confidence = "low"
    if "RAW_DONOR_READS" in flags and plasmid_calls:
        flags.append("RAW_DONOR_ORIGIN_UNRESOLVED")
        confidence = "low"
    return tr.reconcile_cell(sample, sites, n_genome_bc=len(genome_calls),
                             n_plasmid_bc=len(plasmid_calls), bc_confidence=confidence,
                             flags=flags)
