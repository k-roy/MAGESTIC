"""magestic.pipelines.wgs.core.edit_profile — per-position window measurement (defects 6, 11, 15).

WHAT TRACE MEASURES, AND THE TWO THINGS THAT COSTS
--------------------------------------------------
`trace_crispr.wgs.reconcile.window_edit_vaf` returns `(best_vaf, best_cov)`: the
**maximum** non-reference fraction anywhere in the donor window, and the depth at
that one position. `classify_site` then applies `min_cov` to `best_cov`.

Two defects follow, both proven on BAM fixtures in
`tests/test_pipelines/wgs/test_251_trace_fixes.py`:

**Defect 11 (A10/A12) — `min_cov` is gated on the window MAXIMUM depth.** A
soft-clipped deletion, or any dropout, can leave the designed base at depth 0 while
another position in the window sits at 42x. `best_cov` is 42, so `min_cov=4` passes
and the site classifies `SITE_WT`; the cell reports `NO_EDIT`, `n_nocov=0`, and —
because `reconcile_cell` only downgrades confidence when `n_nocov` is non-zero —
`confidence=high`. Measured: depth 0 across a 6 bp hole, reported
`NO_EDIT cov=14 n_nocov=0`. **A confident no-edit over a hole.**

**Defect 6 — one number cannot distinguish a substitution from an indel.** `is_del`,
`is_refskip`, a base mismatch and `indel > 0` all increment the same `nonref`
counter, so a 6 bp deletion and a 1 bp substitution return the same `(1.0, N)` and
both classify `SITE_FULL` -> `SINGLE_EDIT`. Nothing in `SiteCall` records which it
was, so **an SV can be reported as a substitution.**

WHAT THIS MODULE DOES
---------------------
One pileup pass over the window, retaining per-position counts instead of only the
argmax. From that, three things TRACE cannot currently say:

* `min_depth` / `n_positions_uncovered` / `depth_at(pos)` — so the gate can be the
  window MINIMUM, or (better, when the design is known) the designed base itself;
* `variant_class` — SUBSTITUTION / DELETION / INSERTION / REF_SKIP / COMPLEX / NONE,
  decided from which non-reference component dominates at the argmax position;
* `edit_evidence` — a depth-derived evidence tier, because TRACE's `confidence`
  column reflects only the bc1 call and never the edit (defect 15).

🔴 **PARITY IS ASSERTED, NOT ASSUMED.** `max_vaf` and `cov_at_max` are computed by
the *same* comparison TRACE uses, including its tie-break
(`vaf > best_vaf or (vaf == best_vaf and total > best_cov)`), and
`test_fixed_profile_is_vaf_parity_with_trace` asserts equality against
`window_edit_vaf` on three fixtures. The headline VAF provably does not move; this
module only adds what was discarded.

IS THIS AN ADAPTER?
-------------------
It re-reads the same pileup with the same parameters and keeps columns TRACE threw
away, so no new quantity is *defined* here — but the gating decision in
`classify_window` IS a policy choice, and it is magestic's deliberately: TRACE is
pinned, and shipping a half-measured rate is the failure mode this whole exercise is
about. THE TRACE-SIDE FIX: `window_edit_vaf` should return
`(max_vaf, cov_at_max, min_depth, counts_by_class)` and `classify_site` should take
`min_depth` and a `designed_positions` argument. Then `classify_window` collapses to
a call of `classify_site`.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any

from . import _trace

__all__ = [
    "EVIDENCE_MODERATE",
    "PositionOutsideWindowError",
    "EVIDENCE_NONE",
    "EVIDENCE_NO_SITES",
    "EVIDENCE_STRONG",
    "EVIDENCE_THIN",
    "FLAG_SINGLE_ALT_READ",
    "FLAG_SV_ONE_SIDED_LOSS",
    "FLAG_WINDOW_PARTIAL_COVERAGE",
    "PositionCount",
    "ScoredSite",
    "VARIANT_COMPLEX",
    "VARIANT_DELETION",
    "VARIANT_INSERTION",
    "VARIANT_NONE",
    "VARIANT_REF_SKIP",
    "VARIANT_SUBSTITUTION",
    "WindowCall",
    "WindowProfile",
    "classify_window",
    "edit_evidence",
    "one_sided_coverage_loss",
    "window_edit_profile",
]

VARIANT_NONE = "NONE"
VARIANT_SUBSTITUTION = "SUBSTITUTION"
VARIANT_DELETION = "DELETION"
VARIANT_INSERTION = "INSERTION"
VARIANT_REF_SKIP = "REF_SKIP"
VARIANT_COMPLEX = "COMPLEX"

EVIDENCE_NO_SITES = "no_sites_scored"
EVIDENCE_NONE = "none"
EVIDENCE_THIN = "thin"
EVIDENCE_MODERATE = "moderate"
EVIDENCE_STRONG = "strong"

#: Raised as a flag (never silently) when the designed base is callable but some
#: other window position is not. Without it, gating on the designed position would
#: *hide* the partial coverage that gating on the minimum would have surfaced.
FLAG_WINDOW_PARTIAL_COVERAGE = "WINDOW_PARTIAL_COVERAGE"
#: A `partial` call that rests on ONE alt read (the depth-conditional floor allows it, Kevin
#: 2026-10-25). Raised so its cost shows: at 20x a sequencing error does this ~2% of the time.
FLAG_SINGLE_ALT_READ = "PARTIAL_ON_A_SINGLE_ALT_READ"
#: The interim SV gate (Kevin 2026-10-25, planning/284 §8 ruling 3): until rearrangements are
#: detected, a one-sided loss of signal in the donor window EXCLUDES the site (`nocov`).
FLAG_SV_ONE_SIDED_LOSS = "SV_SUSPECT_ONE_SIDED_COVERAGE_LOSS"


class PositionOutsideWindowError(ValueError):
    """A designed position was read outside the profiled window (NN 255 F8).

    Its own class, not a bare `ValueError`, so a caller that genuinely wants to probe
    outside a window can catch exactly this and nothing else.
    """

#: TRACE's own defaults, re-stated so this module never invents a threshold.
#: `classify_site(vaf, cov, tau_full=0.85, tau_wt=0.15, min_cov=4)`.
TAU_FULL = 0.85
TAU_WT = 0.15
MIN_COV = 4
MIN_BASE_QUAL = 13
MIN_MAPQ = 20  # the pinned TRACE's `window_edit_vaf` default; MAPQ 255 is always rejected

#: Evidence tiers. Depth boundaries, nothing cleverer, and stated so a reader of the
#: shipped column knows exactly what it means.
EVIDENCE_STRONG_DEPTH = 100
EVIDENCE_MODERATE_DEPTH = 20


@dataclass(frozen=True)
class PositionCount:
    """Non-reference evidence at one window position.

    `depth` counts reads as the pinned TRACE's `window_edit_vaf` does (planning/268):
    primary, mapped, non-duplicate, QC-passing reads with 20 <= MAPQ < 255, deletions
    included, ref-skips and non-ACGT bases excluded. So `n_nonref / depth` reproduces
    TRACE's VAF exactly. A ref-skip is a MISSING observation, never an alternate
    allele: it is counted in `n_ref_skip` and nowhere else. A position the pileup never
    visited is present with `depth=0` rather than absent, which is what makes
    `min_depth` meaningful.
    """

    pos: int
    ref_base: str
    depth: int = 0
    n_substitution: int = 0
    n_deletion: int = 0
    n_insertion: int = 0
    n_ref_skip: int = 0

    @property
    def n_nonref(self) -> int:
        return self.n_substitution + self.n_deletion + self.n_insertion

    @property
    def vaf(self) -> float:
        return (self.n_nonref / self.depth) if self.depth else 0.0

    @property
    def variant_class(self) -> str:
        """Which non-reference component dominates here.

        A tie between two non-zero components is `COMPLEX`, never one of them: a
        deletion reported as a substitution is the exact failure defect 6 names.
        """
        parts = {
            VARIANT_SUBSTITUTION: self.n_substitution,
            VARIANT_DELETION: self.n_deletion,
            VARIANT_INSERTION: self.n_insertion,
        }
        top = max(parts.values())
        if top == 0:
            return VARIANT_NONE
        winners = [k for k, v in parts.items() if v == top]
        return winners[0] if len(winners) == 1 else VARIANT_COMPLEX


@dataclass(frozen=True)
class WindowProfile:
    """Every position of one donor window, measured once.

    Attributes
    ----------
    chrom, start, end : str, int, int
        1-based inclusive window, as `EditSite` states it.
    positions : tuple of PositionCount
        One entry per window base, in coordinate order. Length is always
        `end - start + 1`.
    max_vaf, cov_at_max, pos_at_max : float, int, int
        TRACE's `(best_vaf, best_cov)` plus WHERE it was found — the offset of the
        max from the designed base is itself a frame-error signal (NN 243 guard 3).
    """

    chrom: str
    start: int
    end: int
    positions: tuple[PositionCount, ...] = ()
    max_vaf: float = 0.0
    cov_at_max: int = 0
    pos_at_max: int = 0

    @property
    def min_depth(self) -> int:
        """Lowest depth anywhere in the window. 0 when any base is uncovered."""
        return min((p.depth for p in self.positions), default=0)

    @property
    def n_positions(self) -> int:
        return len(self.positions)

    @property
    def n_positions_uncovered(self) -> int:
        """Window bases below `MIN_COV`. Counted, so a hole is never invisible."""
        return sum(1 for p in self.positions if p.depth < MIN_COV)

    @property
    def variant_class(self) -> str:
        """The dominant non-reference class at the argmax position."""
        for p in self.positions:
            if p.pos == self.pos_at_max:
                return p.variant_class
        return VARIANT_NONE

    def contains(self, pos: int) -> bool:
        """Is `pos` inside the profiled window? 1-based inclusive."""
        return self.start <= int(pos) <= self.end

    def _require_contains(self, pos: int, caller: str) -> None:
        """🔴 NN 255 F8. A position outside the window is a COORDINATE BUG, not a zero.

        These three accessors used to `return 0` / `0.0` / `NONE` for an unprofiled
        position, silently. Measured consequence (NN 252): a fully-covered 40x window
        with one out-of-window designed position classified `nocov` with
        `reason='designed_position_depth_below_min_cov'` — **the reason named depth, not
        a coordinate problem** — and a wrong but IN-window position turned a real 100%
        edit into `wt / vaf=0.0 / reason='designed_position_wt_vaf'`, a confident
        no-edit.

        NN 252 could not build a live route to the in-window case and labelled it
        `SUSPECTED, NOT REPRODUCED:`; this is the one-assert-wide invariant it said was
        nonetheless real. The out-of-window case it DID reproduce now raises here.
        """
        if not self.contains(pos):
            raise PositionOutsideWindowError(
                f"{caller}({pos}) is outside the profiled window "
                f"{self.chrom}:{self.start}-{self.end}. That is a coordinate-frame "
                "error, not a depth of zero: returning 0 here made a fully-covered "
                "window read as NOCOV with a reason naming depth. Check that the "
                "designed positions and the window come from the same table row."
            )

    def depth_at(self, pos: int) -> int:
        self._require_contains(pos, "depth_at")
        for p in self.positions:
            if p.pos == pos:
                return p.depth
        return 0

    def vaf_at(self, pos: int) -> float:
        self._require_contains(pos, "vaf_at")
        for p in self.positions:
            if p.pos == pos:
                return p.vaf
        return 0.0

    def variant_class_at(self, pos: int) -> str:
        self._require_contains(pos, "variant_class_at")
        for p in self.positions:
            if p.pos == pos:
                return p.variant_class
        return VARIANT_NONE


@dataclass(frozen=True)
class WindowCall:
    """A classified window: TRACE's status, plus WHY, plus the variant type.

    `reason` is mandatory and machine-readable. No status in this module can be
    produced without one — that is the deliverable quality the audits asked for.
    """

    status: str
    reason: str
    variant_class: str
    vaf: float
    cov: int
    min_depth: int
    n_positions_uncovered: int
    depth_at_designed: int | None = None
    vaf_at_designed: float | None = None
    designed_positions: tuple[int, ...] = ()
    flags: tuple[str, ...] = field(default_factory=tuple)


@dataclass
class ScoredSite:
    """A `SiteCall` that also records WHY, and WHAT KIND of variant.

    🔴 **Duck-compatible with `trace_crispr.wgs.reconcile.SiteCall` on purpose.**
    `reconcile_cell` reads `.status`, `.flags` and `.allele_evidence`;
    `reconciliation_to_row` reads `.site`, `.vaf`, `.cov`, `.status` and
    `.allele_evidence`. Carrying the same attributes means this can be handed straight
    to TRACE's own aggregation and emit code, so the extra fields are additive and
    nothing in TRACE needs to change.

    `vaf` and `cov` are TRACE's `window_edit_vaf` values unchanged, so the headline
    numbers in an existing table do not move. `status` MAY differ from
    `classify_site`'s, and `reason` always says why.
    """

    site: Any
    vaf: float
    cov: int
    status: str
    reason: str = ""
    variant_class: str = VARIANT_NONE
    min_depth: int = 0
    n_positions_uncovered: int = 0
    depth_at_designed: int | None = None
    trace_status: str = ""
    flags: tuple[str, ...] = field(default_factory=tuple)
    #: Per-allele REF/ALT counts; empty here, because a window scored by this module has
    #: no declared intended allele. The pinned TRACE's `reconcile_cell` and
    #: `reconciliation_to_row` read it from every site (planning/268).
    allele_evidence: tuple[Any, ...] = field(default_factory=tuple)

    @property
    def status_changed_from_trace(self) -> bool:
        """True when the depth gate reclassified this site. Never silent."""
        return bool(self.trace_status) and self.status != self.trace_status


def window_edit_profile(
    bam: Any,
    fasta: Any,
    chrom: str,
    start: int,
    end: int,
    *,
    min_base_qual: int = MIN_BASE_QUAL,
    min_mapq: int = MIN_MAPQ,
    min_cov: int = MIN_COV,
) -> WindowProfile:
    """Measure every position of `[start, end]` (1-based inclusive) in one pass.

    The pileup call, the read filter, the `min_base_quality`, the non-reference
    predicate and the choice of the reported position are the pinned TRACE's
    `window_edit_vaf`, so `max_vaf`/`cov_at_max` are its return values (planning/268
    re-aligned them when the pin moved). Everything else is what TRACE discards.
    """
    ref = fasta.fetch(chrom, start - 1, end).upper()
    counts: dict[int, dict[str, int]] = {}
    best_vaf, best_cov, best_pos = 0.0, 0, start
    covered = False

    for col in bam.pileup(
        chrom, start - 1, end, truncate=True, stepper="samtools",
        min_base_quality=min_base_qual, flag_filter=0xF04,
    ):
        i = col.reference_pos - (start - 1)
        if i < 0 or i >= len(ref):
            continue
        refbase = ref[i]
        sub = dele = ins = skip = 0
        total = 0
        for pr in col.pileups:
            mapq = pr.alignment.mapping_quality
            if mapq == 255 or mapq < min_mapq:
                continue
            if pr.is_refskip:
                skip += 1  # a missing observation, never depth or an alternate allele
                continue
            if pr.is_del:
                total += 1
                dele += 1
                continue
            qpos = pr.query_position
            if qpos is None:
                continue
            base = pr.alignment.query_sequence[qpos].upper()
            if base not in "ACGT":
                continue
            total += 1
            if pr.indel > 0:
                # An insertion is scored at its anchor base, as TRACE does. When the
                # anchor ALSO mismatches, the read is counted once, as an insertion:
                # TRACE's `or` counts it once too, so the VAF is identical.
                ins += 1
            elif base != refbase:
                sub += 1
        counts[col.reference_pos + 1] = {
            "depth": total, "sub": sub, "del": dele, "ins": ins, "skip": skip,
        }
        if total == 0:
            continue
        vaf = (sub + dele + ins) / total
        # TRACE's exact choice: the highest (vaf, depth) among positions reaching
        # `min_cov`, else the deepest observed position, so low coverage stays explicit.
        if total >= min_cov:
            if not covered or (vaf, total) > (best_vaf, best_cov):
                best_vaf, best_cov, best_pos = vaf, total, col.reference_pos + 1
            covered = True
        elif not covered and (total, vaf) > (best_cov, best_vaf):
            best_vaf, best_cov, best_pos = vaf, total, col.reference_pos + 1

    positions = tuple(
        PositionCount(
            pos=p,
            ref_base=ref[p - start] if 0 <= p - start < len(ref) else "N",
            depth=counts.get(p, {}).get("depth", 0),
            n_substitution=counts.get(p, {}).get("sub", 0),
            n_deletion=counts.get(p, {}).get("del", 0),
            n_insertion=counts.get(p, {}).get("ins", 0),
            n_ref_skip=counts.get(p, {}).get("skip", 0),
        )
        for p in range(start, end + 1)
    )
    return WindowProfile(
        chrom=chrom, start=start, end=end, positions=positions,
        max_vaf=best_vaf, cov_at_max=best_cov, pos_at_max=best_pos,
    )


def classify_window(
    profile: WindowProfile,
    *,
    designed_positions: Sequence[int] | None = None,
    tau_full: float = TAU_FULL,
    tau_wt: float = TAU_WT,
    min_cov: int = MIN_COV,
) -> WindowCall:
    """Classify a window, gating depth where it actually matters.

    Gate, in order of preference:

    1. **At the designed positions**, when the caller knows them (they come from
       `window_table.designed_positions`, i.e. from `wt_seq` vs the design — defect
       10). This is the quantity that decides whether the edit was observable.
    2. **At the window minimum** otherwise. This is the drop-in that needs no new
       input and still closes defect 11.

    When (1) is used, a window position below `min_cov` elsewhere does not change the
    status but DOES raise `FLAG_WINDOW_PARTIAL_COVERAGE`, so the partial coverage the
    minimum would have surfaced is still reported rather than hidden by the narrower
    gate.

    The VAF is read at the designed position when known, and at the argmax otherwise.
    Reading it at the designed position is what stops a background polymorphism
    elsewhere in the window from being promoted to "the edit" (NN 243's A4b).
    """
    tr = _trace.reconcile()
    flags: list[str] = []
    dpos = tuple(int(p) for p in (designed_positions or ()))

    # 🔴 NN 255 F8 — the containment invariant, stated once, before anything is gated.
    # `designed_positions` and the profile window must come from the same window-table
    # row; `build_edit_sites` takes `(chrom, start, end)` straight from
    # `oligo_windows`, and `designed_positions_by_oligo` derives its positions from the
    # same `start`, so a disagreement here means the two sources have diverged. That is
    # worth stopping for: NN 243 measured that a wrong frame flips every per-cell call
    # while the apparent aggregate rate lands on exactly truth, so no rate can detect it.
    outside = tuple(p for p in dpos if not profile.contains(p))
    if outside:
        raise PositionOutsideWindowError(
            f"designed position(s) {outside} fall outside the profiled window "
            f"{profile.chrom}:{profile.start}-{profile.end} (designed: {dpos}). "
            "Gating there would read depth 0 and VAF 0.0 and classify a fully-covered "
            "window as NOCOV, with a reason naming depth rather than the coordinate "
            "frame. The designed positions and the window must come from the same "
            "window-table row."
        )

    if dpos:
        # 2026-10-25 rulings (planning/268 §4): each designed base is classified from its COUNTS by
        # TRACE's `classify_counts` -- presence from a single read up when the non-reference count
        # beats sequencing error at that depth (FLOOR), absence only at `min_cov` reads (UNCALLABLE).
        depth_at_designed = min(profile.depth_at(p) for p in dpos)
        vaf_at_designed = max(profile.vaf_at(p) for p in dpos)
        vclass = next(
            (
                profile.variant_class_at(p)
                for p in dpos
                if profile.variant_class_at(p) != VARIANT_NONE
            ),
            VARIANT_NONE,
        )
        if profile.min_depth < min_cov:
            flags.append(FLAG_WINDOW_PARTIAL_COVERAGE)
        counts = {p: next(c for c in profile.positions if c.pos == p) for p in dpos}
        statuses = {p: tr.classify_counts(c.n_nonref, c.depth, tau_full=tau_full, min_cov=min_cov)
                    for p, c in counts.items()}
        values = set(statuses.values())
        if tr.SITE_PARTIAL in values:
            status, reason = tr.SITE_PARTIAL, "designed_position_minority_allele"
            if any(counts[p].n_nonref == 1 for p, st in statuses.items() if st == tr.SITE_PARTIAL):
                flags.append(FLAG_SINGLE_ALT_READ)
        elif tr.SITE_NOCOV in values:
            status, reason = tr.SITE_NOCOV, "designed_position_depth_below_min_cov"
        elif values == {tr.SITE_FULL}:
            status, reason = tr.SITE_FULL, "designed_position_full_vaf"
        elif values == {tr.SITE_WT}:
            status, reason = tr.SITE_WT, "designed_position_wt_vaf"
        else:
            status, reason = tr.SITE_INCOMPLETE, "designed_positions_partly_installed"
        return WindowCall(
            status=status, reason=reason, variant_class=vclass, vaf=vaf_at_designed,
            cov=depth_at_designed, min_depth=profile.min_depth,
            n_positions_uncovered=profile.n_positions_uncovered,
            depth_at_designed=depth_at_designed, vaf_at_designed=vaf_at_designed,
            designed_positions=dpos, flags=tuple(flags),
        )

    if profile.min_depth < min_cov:
        return WindowCall(
            status=tr.SITE_NOCOV, reason="window_min_depth_below_min_cov",
            variant_class=profile.variant_class, vaf=profile.max_vaf,
            cov=profile.cov_at_max, min_depth=profile.min_depth,
            n_positions_uncovered=profile.n_positions_uncovered, flags=tuple(flags),
        )
    if profile.max_vaf >= tau_full:
        status, reason = tr.SITE_FULL, "window_max_full_vaf"
    elif profile.max_vaf <= tau_wt:
        status, reason = tr.SITE_WT, "window_max_wt_vaf"
    else:
        status, reason = tr.SITE_PARTIAL, "window_max_intermediate_vaf"
    return WindowCall(
        status=status, reason=reason, variant_class=profile.variant_class,
        vaf=profile.max_vaf, cov=profile.cov_at_max, min_depth=profile.min_depth,
        n_positions_uncovered=profile.n_positions_uncovered, flags=tuple(flags),
    )


def one_sided_coverage_loss(
    bam: Any,
    profile: WindowProfile,
    *,
    min_run: int = 5,
    min_flank_depth: int = 3,
    cluster_bp: int = 2,
    min_mapq: int = MIN_MAPQ,
) -> str | None:
    """`"left"` / `"right"` when the donor window loses all signal on ONE side, else `None`.

    The interim SV gate (Kevin 2026-10-25, planning/284 §8 ruling 3): until inversions and
    translocations are detected, a site whose reads stop dead on one side is excluded rather
    than read as `NO_EDIT`. All three must hold, so ordinary coverage noise is not flagged:
    a zero-depth run of at least `min_run` bases touching that window edge (and not the other);
    at least `min_flank_depth` reads at the last covered base; and at least half of those reads,
    and at least 3, ENDING (or, for a left-side loss, starting) within `cluster_bp` of that base,
    which is the pile-up a breakpoint leaves and random read ends do not. A designed deletion
    inside the window does not touch an edge, so it is never flagged.
    """
    depths = [p.depth for p in profile.positions]
    if not depths:
        return None
    lead = next((i for i, d in enumerate(depths) if d), len(depths))
    trail = next((i for i, d in enumerate(reversed(depths)) if d), len(depths))
    if lead == len(depths) or (lead >= min_run) == (trail >= min_run):
        return None
    side = "left" if lead >= min_run else "right"
    edge_index = lead if side == "left" else len(depths) - 1 - trail
    edge = profile.start + edge_index  # the last (right) or first (left) covered base, 1-based
    flank = depths[edge_index]
    if flank < min_flank_depth:
        return None
    clustered = 0
    for read in bam.fetch(profile.chrom, max(edge - 1 - cluster_bp, 0), edge + cluster_bp):
        if (read.flag & 0xF04 or read.mapping_quality == 255 or read.mapping_quality < min_mapq
                or read.reference_end is None):
            continue
        tip = read.reference_end if side == "right" else read.reference_start + 1
        clustered += abs(tip - edge) <= cluster_bp
    return side if clustered >= max(3, flank / 2) else None


def edit_evidence(sites: Iterable[Any]) -> str:
    """A depth-derived evidence tier for a cell's scored sites (defect 15).

    `reconcile_cell`'s `confidence` column is `"high" if bc_confidence == "high"`,
    downgraded only for `LOW_CONF` or an uncovered `NO_EDIT`. It therefore reflects
    the **barcode** call and never the edit: a 4-read site and a 400-read site are
    both `high`. NN 243 measured exactly this — every `confidence=high` in its
    harness was inherited from `bc_confidence`, and in a bridge-less cohort the
    column is uniformly `low` and carries no information either way.

    This is a separate, explicitly-named quantity. It does not overwrite
    `confidence`; `emit.RECONCILE_COLUMNS` ships both, with TRACE's renamed
    `bc_confidence` so the ambiguous bare name cannot mislead.
    """
    sites = list(sites)
    if not sites:
        return EVIDENCE_NO_SITES
    tr = _trace.reconcile()
    # An `unassessed` site (no declared intended alleles, pinned TRACE) has a depth but
    # measured no edit, so it is evidence of nothing here, exactly like `nocov`.
    unmeasured = {tr.SITE_NOCOV, getattr(tr, "SITE_UNASSESSED", tr.SITE_NOCOV)}
    covered = [s for s in sites if getattr(s, "status", None) not in unmeasured]
    if not covered:
        return EVIDENCE_NONE
    best = max(int(getattr(s, "cov", 0) or 0) for s in covered)
    if best >= EVIDENCE_STRONG_DEPTH:
        return EVIDENCE_STRONG
    if best >= EVIDENCE_MODERATE_DEPTH:
        return EVIDENCE_MODERATE
    return EVIDENCE_THIN
