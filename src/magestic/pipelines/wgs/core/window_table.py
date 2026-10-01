"""magestic.pipelines.wgs.core.window_table — denominators and frame for Stage 4 (defects 10, 12).

DEFECT 10 — THE ALLELE INVERSION, AND WHY IT IS THE HIGHEST-VALUE FIX HERE
---------------------------------------------------------------------------
Stage 4 scores a window by comparing reads to `fasta.fetch(chrom, start-1, end)`.
It is **never given the design's `wt_seq`**, and its only reference check,
`locus_config.validate_tracks`, validates the bc1 locus alone. So when the reference
frame is wrong — the reference already carries the design's ALT at the designed base
— every per-cell call FLIPS, and NN 243 measured the consequence:

* whole-frame error: **0.0% / 100.0%** apparent rate, no warning (A1);
* planted 50% truth: **50.0% apparent — exactly truth — with per-cell accuracy 0.0%**
  (A11);
* mixed frame damage: 55.0% apparent vs 50.0% truth, 45.0% accuracy (A11b), i.e.
  **the rate looks BETTER the worse the damage is.**

🔴 **The aggregate edit rate is therefore not a detector for reference-frame error.**
Any validation that compares rates passes while every call is wrong. The only check
that works is a per-site one, and `trace_crispr.wgs.donor_outcome.orient_to_genome`
already *is* that check — NN 243 confirmed it catches frame error in Stage 2. It is
simply never reached from Stage 4, because the window table carries no `wt_seq`.

`verify_window_frame` reuses that same comparison (including its minus-strand
fallback) per site. This converts an undetectable failure into a loud one.

⚠️ **Scope, preserved from NN 243 and not to be restated as a live error:** the live
22.6% is a **Stage-2** number, and Stage 4 runs on no shipped config today
(`reconcile.oligo_windows` is a `/path/to/` placeholder). These are pre-deployment
landmines in the scoring core shared with TRACE, not wrong numbers in tonight's rates.

DEFECT 12 — A DROPPED WINDOW IS A CONFIDENT UNDERCOUNT
-------------------------------------------------------
`reconcile.load_oligo_windows` silently skips a row whose coordinates fail
`.isdigit()` and keeps the FIRST occurrence per oligo, reporting neither. An oligo
with no window is then simply **absent** from `build_edit_sites`, so `MULTI_EDIT`
becomes `SINGLE_EDIT` with no flag and without being counted as `nocov`.

This is live-adjacent, not hypothetical: `readers/oligo_tables` measures
`donor_start_coord`/`donor_end_coord` as only **50.2% filled (20,177/40,191)** over
the 2021-pool rows. `3001.0`, `3,001`, a bare CR and a duplicate oligo row all
trigger it.

`load_windows_audited` returns TRACE's exact window dict — byte-for-byte the same
acceptance rule, so Stage 4's numbers cannot move — **plus every denominator**, and
`strict=True` refuses. `unwindowed_oligos` names the oligos that must be counted as
`SITE_NOCOV` rather than omitted.

THE TRACE-SIDE FIX
------------------
`load_oligo_windows` should return `(windows, LoadReport)` with `n_rows_read`,
`n_dropped`, `n_duplicate_oligos` and the dropped/duplicate oligo names, and should
accept an optional `wt_seq_col`. `build_edit_sites` should return the oligos it could
not place alongside the sites it could. `reconcile_driver`'s docstring already
specifies the comparison; nothing performed it.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import _trace

__all__ = [
    "DesignOrientation",
    "DesignOrientationError",
    "FRAME_UNVERIFIED",
    "FRAME_VERIFIED_MINUS",
    "FRAME_VERIFIED_PLUS",
    "FrameCheck",
    "FrameMismatchError",
    "WindowAudit",
    "WindowTableError",
    "designed_positions",
    "designed_positions_by_oligo",
    "load_windows_audited",
    "unwindowed_oligos",
    "verify_design_orientation",
    "verify_window_frame",
]

FRAME_UNVERIFIED = "frame_unverified_no_wt_seq"
FRAME_VERIFIED_PLUS = "frame_verified_plus_strand"
FRAME_VERIFIED_MINUS = "frame_verified_minus_strand"


class WindowTableError(ValueError):
    """The window table dropped rows, or collapsed duplicate oligos, under `strict`."""


class FrameMismatchError(ValueError):
    """The reference window does not match the design's `wt_seq` in either orientation."""


@dataclass(frozen=True)
class WindowAudit:
    """TRACE's window dict plus every denominator it discards.

    Attributes
    ----------
    windows : dict
        `{oligo: (chrom, start, end)}` — **exactly** what `load_oligo_windows`
        returns for the same file, including its first-wins and `.isdigit()` rules.
        Parity is the point: this module reports, it does not reinterpret.
    n_rows_read : int
        Data rows in the file. The denominator for everything below.
    n_windows : int
        `len(windows)`.
    n_dropped_bad_coords : int
        Rows whose start or end is not a bare integer (`3001.0`, `3,001`, empty).
    n_dropped_short_row : int
        Rows with fewer columns than the header declares.
    n_duplicate_oligos : int
        Rows whose oligo was already placed. **The first wins, silently, in TRACE.**
    frame_status : str
        `FRAME_UNVERIFIED` when the table carries no `wt_seq` column — so a reader
        can never mistake "the guard did not run" for "the guard passed".
    """

    path: Path
    windows: dict[str, tuple[str, int, int]]
    wt_seqs: dict[str, str]
    design_seqs: dict[str, str] = field(default_factory=dict)
    n_rows_read: int = 0
    n_dropped_bad_coords: int = 0
    n_dropped_short_row: int = 0
    n_dropped_no_oligo: int = 0
    n_duplicate_oligos: int = 0
    dropped_oligos: tuple[str, ...] = ()
    duplicate_oligos: tuple[str, ...] = ()
    has_wt_seq: bool = False

    @property
    def n_windows(self) -> int:
        return len(self.windows)

    @property
    def n_dropped(self) -> int:
        return (
            self.n_dropped_bad_coords + self.n_dropped_short_row + self.n_dropped_no_oligo
        )

    @property
    def frame_status(self) -> str:
        return FRAME_VERIFIED_PLUS if self.has_wt_seq else FRAME_UNVERIFIED

    def summary(self) -> str:
        return (
            f"{self.path.name}: {self.n_windows} window(s) from {self.n_rows_read} row(s) "
            f"— dropped {self.n_dropped} (bad coords {self.n_dropped_bad_coords}, "
            f"short rows {self.n_dropped_short_row}, no oligo {self.n_dropped_no_oligo}), "
            f"duplicate oligos {self.n_duplicate_oligos}, wt_seq column "
            f"{'present' if self.has_wt_seq else 'ABSENT'}"
        )


def load_windows_audited(
    path: str | Path,
    *,
    oligo_col: str = "oligo_name",
    chrom_col: str = "chrom",
    start_col: str = "donor_start_coord",
    end_col: str = "donor_end_coord",
    wt_seq_col: str = "wt_seq",
    design_seq_col: str = "design_donor",
    strict: bool = False,
) -> WindowAudit:
    """`load_oligo_windows` with every denominator returned, and an optional refusal.

    The acceptance rule is TRACE's, deliberately unchanged: a coordinate must be a
    bare `.isdigit()` integer and the first row per oligo wins. Loosening it would
    silently change which sites Stage 4 scores, which is a different (and worse)
    defect than the one being fixed.

    Fields are `.strip()`ed, which also removes a stray `\\r` from a CRLF/CR file —
    an invisible line ending must not look like a bad coordinate.
    """
    path = Path(path)
    windows: dict[str, tuple[str, int, int]] = {}
    wt_seqs: dict[str, str] = {}
    design_seqs: dict[str, str] = {}
    n_rows = n_bad = n_short = n_no_oligo = n_dupe = 0
    dropped: list[str] = []
    dupes: list[str] = []

    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = [h.strip() for h in next(reader)]
        except StopIteration as exc:
            raise WindowTableError(f"{path} is empty — it has not even a header") from exc
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in (oligo_col, chrom_col, start_col, end_col) if c not in idx]
        if missing:
            raise WindowTableError(
                f"{path}: window table is missing required column(s) {missing}; "
                f"header is {header}"
            )
        oi, ci, si, ei = idx[oligo_col], idx[chrom_col], idx[start_col], idx[end_col]
        wi = idx.get(wt_seq_col)
        di = idx.get(design_seq_col)
        has_wt = wi is not None

        for row in reader:
            if not any(f.strip() for f in row):
                continue
            n_rows += 1
            if len(row) <= max(oi, ci, si, ei):
                n_short += 1
                dropped.append(row[oi].strip() if len(row) > oi else "<short row>")
                continue
            oligo = row[oi].strip()
            start_s, end_s = row[si].strip(), row[ei].strip()
            if not oligo:
                n_no_oligo += 1
                continue
            if not start_s.isdigit() or not end_s.isdigit():
                n_bad += 1
                dropped.append(oligo)
                continue
            if oligo in windows:
                n_dupe += 1
                dupes.append(oligo)
                continue
            windows[oligo] = (row[ci].strip(), int(start_s), int(end_s))
            if has_wt and len(row) > wi:
                wt_seqs[oligo] = row[wi].strip().upper()
            if di is not None and len(row) > di:
                design_seqs[oligo] = row[di].strip().upper()

    audit = WindowAudit(
        path=path, windows=windows, wt_seqs=wt_seqs, design_seqs=design_seqs,
        n_rows_read=n_rows,
        n_dropped_bad_coords=n_bad, n_dropped_short_row=n_short,
        n_dropped_no_oligo=n_no_oligo, n_duplicate_oligos=n_dupe,
        dropped_oligos=tuple(dropped), duplicate_oligos=tuple(dupes),
        has_wt_seq=has_wt,
    )
    if strict and (audit.n_dropped or audit.n_duplicate_oligos):
        raise WindowTableError(
            f"refusing a window table with silent losses — {audit.summary()}. "
            f"{audit.n_dropped} of {audit.n_rows_read} row(s) produced no window and "
            f"{audit.n_duplicate_oligos} duplicate oligo row(s) were collapsed. Each "
            "lost window makes an edit at that site invisible, which reads downstream "
            "as a confident undercount (MULTI_EDIT -> SINGLE_EDIT) rather than as "
            "missing data. Fix the table, or pass strict=False and record the counts."
        )
    return audit


def unwindowed_oligos(
    oligos: Iterable[str], audit: WindowAudit
) -> tuple[str, ...]:
    """Requested oligos with no window — the set that must count as `SITE_NOCOV`.

    `build_edit_sites` simply omits them, which shrinks the denominator. Naming them
    is what lets `step_4` report them as uncovered instead of as absent.
    """
    return tuple(dict.fromkeys(o for o in oligos if o and o not in audit.windows))


@dataclass(frozen=True)
class FrameCheck:
    """The result of comparing a reference window to the design's `wt_seq`."""

    chrom: str
    start: int
    end: int
    oriented: str  # "plus" / "minus"
    n_mismatch: int = 0
    status: str = FRAME_VERIFIED_PLUS


def verify_window_frame(
    fasta: Any,
    chrom: str,
    start: int,
    end: int,
    wt_seq: str,
) -> FrameCheck:
    """Assert the reference window IS the design's `wt_seq`, in either orientation.

    This is `donor_outcome.orient_to_genome`'s comparison, applied per site at Stage
    4 where it is currently absent. A minus-strand donor is accepted (and reported),
    because the design tables legitimately store either orientation; anything else
    raises.

    Raises
    ------
    FrameMismatchError
        Neither orientation matches. The message carries the mismatch count and the
        first differing offset, because a 1-base mismatch at the designed position
        is the allele-inversion signature and a wholesale mismatch is a wrong contig
        or a wrong coordinate frame — different fixes.
    """
    do = _trace.donor_outcome()
    wt = str(wt_seq).upper()
    window = fasta.fetch(chrom, start - 1, end).upper()
    if window == wt:
        return FrameCheck(chrom, start, end, "plus", 0, FRAME_VERIFIED_PLUS)
    rc = do.reverse_complement(wt)
    if window == rc:
        return FrameCheck(chrom, start, end, "minus", 0, FRAME_VERIFIED_MINUS)

    def _diff(a: str, b: str) -> tuple[int, int]:
        n = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
        first = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), -1)
        return n, first

    n_plus, first_plus = _diff(window, wt)
    n_minus, _ = _diff(window, rc)
    n, first = (n_plus, first_plus) if n_plus <= n_minus else (n_minus, -1)
    raise FrameMismatchError(
        f"reference frame does not match the design at {chrom}:{start}-{end}: "
        f"{n} mismatch(es) against wt_seq"
        + (f", first at window offset {first} (genomic {start + first})" if first >= 0 else "")
        + f" [plus {n_plus}, minus {n_minus}, len ref {len(window)} vs wt {len(wt)}]. "
        "🔴 A reference that already carries the design's ALT inverts every per-cell "
        "call while leaving the aggregate rate on top of truth, so the rate cannot "
        "detect this. Rebuild the composite reference for this strain background "
        "(WGS_REFERENCE_MAPPING_POLICY) before scoring."
    )


def designed_positions_by_oligo(
    audit: WindowAudit,
    orientations: Mapping[str, str] | None = None,
) -> dict[str, tuple[int, ...]]:
    """`{oligo: designed genomic positions}`, for oligos whose STRAND is KNOWN.

    This is what lets `reconcile_driver.score_sites` gate depth and VAF at the base
    the edit actually touches (defect 11's better half) instead of at the window
    minimum.

    🔴 **`orientations` is required in practice, and an oligo without one is OMITTED
    rather than assumed plus-strand.** The design tables legitimately store either
    orientation — `verify_window_frame` accepts both and reports which — and
    `designed_positions` maps string index `i` to genomic `window_start + i`, which is
    true only on the plus strand. On the minus strand `wt[i]` is genomic `end - i`, so
    assuming plus would gate depth and read VAF at the **mirror base**: a hole at the
    real edit would pass as covered and a real edit would read as WT. That is the A11
    class of silent wrong answer, and guessing the strand would reintroduce it while
    fixing defect 11.

    Omission is the safe failure: `score_sites` then gates on the window MINIMUM for
    that oligo, which is conservative, never on TRACE's window maximum.

    Parameters
    ----------
    audit : WindowAudit
        From `load_windows_audited`, with both `wt_seq` and `design_donor` columns.
    orientations : mapping, optional
        `{oligo: "plus" | "minus"}`, from `verify_window_frame(...).oriented`. An
        oligo absent from it — or the mapping being `None` — omits that oligo.
    """
    out: dict[str, tuple[int, ...]] = {}
    for oligo, (_chrom, start, end) in audit.windows.items():
        wt = audit.wt_seqs.get(oligo)
        design = audit.design_seqs.get(oligo)
        oriented = (orientations or {}).get(oligo)
        if not wt or not design or oriented not in ("plus", "minus"):
            continue
        # 🔴 NN 255 F9 — the strand of `design_donor` was NEVER checked.
        # `verify_window_frame` compares only `wt_seq` to the reference, and the
        # orientation it returns is then applied by `designed_positions` to BOTH
        # columns. With `wt_seq` stored minus and `design_donor` stored plus — two
        # different columns of a design table, so a real possibility — NN 252 measured
        # `designed_positions` returning **10 positions, the whole window, where truth
        # is 1**, which degrades the gate back to the window minimum and the VAF back to
        # TRACE's window maximum (re-admitting NN 243's A4b, a background polymorphism
        # promoted to "the edit") with no flag distinguishing it from a legitimately
        # narrow gate.
        verify_design_orientation(oligo, wt, design)
        pos = designed_positions(start, wt, design, oriented=oriented)
        if pos:
            out[oligo] = pos
    return out


class DesignOrientationError(ValueError):
    """`design_donor` and `wt_seq` are not written on the same strand (NN 255 F9)."""


@dataclass(frozen=True)
class DesignOrientation:
    """How `design_donor` compares to `wt_seq` as stored and reverse-complemented."""

    oligo: str
    n_mismatch_as_stored: int
    n_mismatch_reverse_complement: int
    n_compared: int

    @property
    def same_strand(self) -> bool:
        """As-stored must be at least as close as the reverse complement.

        Ties pass as same-strand: that is the conservative direction, because a false
        alarm on a palindromic window would refuse a legitimate design.
        """
        return self.n_mismatch_as_stored <= self.n_mismatch_reverse_complement


def verify_design_orientation(
    oligo: str, wt_seq: str, design_seq: str
) -> DesignOrientation:
    """Refuse a `(wt_seq, design_donor)` pair whose two strands disagree.

    The test does not need the reference: a design donor differs from its own `wt_seq`
    at the designed positions and nowhere else, so **as stored it must be at least as
    similar to `wt_seq` as its reverse complement is.** When the reverse complement is
    STRICTLY closer, the two columns are written on opposite strands.

    A second, independent refusal catches the measured symptom directly: a design that
    differs from `wt_seq` at EVERY compared position is not a design, it is a different
    sequence — and that is exactly what the mixed-strand case produces.

    Raises
    ------
    DesignOrientationError
        Either condition. Loud, because the alternative is a gate silently widened to
        the whole window with no flag — see the comment in
        `designed_positions_by_oligo`.
    """
    rc = _trace.donor_outcome().reverse_complement
    wt = str(wt_seq).upper()
    des = str(design_seq).upper()

    def _diff(a: str, b: str) -> int:
        return sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))

    n_same = _diff(wt, des)
    n_rc = _diff(wt, rc(des))
    n_cmp = min(len(wt), len(des))
    result = DesignOrientation(oligo, n_same, n_rc, n_cmp)
    if not result.same_strand:
        raise DesignOrientationError(
            f"oligo {oligo}: `design_donor` is closer to the REVERSE COMPLEMENT of "
            f"`wt_seq` ({n_rc} mismatch(es)) than to `wt_seq` as stored ({n_same}), so "
            "the two design columns are written on opposite strands. Nothing in "
            "`verify_window_frame` compares `design_donor` to anything, so the "
            "orientation verified from `wt_seq` would be applied to both and the "
            "designed positions would span the whole window — degrading the depth gate "
            "to the window minimum and the VAF to TRACE's window MAXIMUM, which "
            "re-admits a background polymorphism as 'the edit'. Reverse-complement one "
            "column in the window table."
        )
    if n_cmp >= 3 and n_same >= n_cmp:
        raise DesignOrientationError(
            f"oligo {oligo}: `design_donor` differs from `wt_seq` at ALL {n_cmp} "
            "compared position(s), so every base of the window would be treated as a "
            "designed position and the gate would collapse to the window minimum. That "
            "is a wrong pair of columns, not a design."
        )
    return result


def designed_positions(
    window_start: int, wt_seq: str, design_seq: str, *, oriented: str = "plus"
) -> tuple[int, ...]:
    """1-based genomic positions where the design differs from `wt_seq`.

    These are the only positions at which the edit is observable, so they are what
    `edit_profile.classify_window` gates depth and VAF on (defect 11's better half).

    Parameters
    ----------
    oriented : {"plus", "minus"}
        Which strand `wt_seq`/`design_seq` are written on, as reported by
        `verify_window_frame(...).oriented`. 🔴 **A minus-oriented pair is
        reverse-complemented before diffing**, because otherwise string index `i`
        would map to `window_start + i` when the actual genomic position is
        `window_end - i` — gating on the mirror base. Both strings are flipped, so
        the diff is taken in genomic orientation and the arithmetic below is correct
        for either input.

    Notes
    -----
    Equal-length sequences are compared base by base. An indel design (unequal
    lengths) returns the common-prefix divergence point onward up to the shorter
    length, which is deliberately conservative: it never claims a position outside
    the window, and a caller that needs exact indel coordinates should take them
    from the design table rather than from a diff.
    """
    if oriented not in ("plus", "minus"):
        raise ValueError(
            f"oriented must be 'plus' or 'minus', got {oriented!r}; the strand cannot "
            "be guessed — a wrong guess gates the mirror base"
        )
    wt = str(wt_seq).upper()
    des = str(design_seq).upper()
    if oriented == "minus":
        rc = _trace.donor_outcome().reverse_complement
        wt, des = rc(wt), rc(des)
    n = min(len(wt), len(des))
    out = [window_start + i for i in range(n) if wt[i] != des[i]]
    if len(wt) != len(des) and not out:
        out = [window_start + n - 1] if n else []
    return tuple(out)
