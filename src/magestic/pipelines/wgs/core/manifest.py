"""
magestic.pipelines.wgs.core.manifest — the per-sample manifest, and D1's home.

**Reference is a per-SAMPLE attribute, never a cohort-level `config.genome`.** That is
the dominant structural constraint on this module and it overturns the obvious design:
a BAM header sweep over the program's 1137 clonal BAMs (100%, reference read from
`@PG ID:bwa`) found one BAM *directory* spanning two references (315 + 342). A
`genome:` scalar in the screen config would be structurally wrong, so every stage
resolves its reference from this manifest and runs within a reference group.

This module is pure: it reads a TSV, validates it, and compares the declared reference
against an observed one supplied by the caller. Reading `@PG`/`@RG` out of a BAM belongs
to `readers.bam_inventory`, which hands the observation here as a plain mapping — so the
cross-check is unit-testable without a BAM.

It also owns the coordinate-portability refusal (§8.4): chrVIII is CUP1-contracted and
its length differs across the three references in use (568,236 / 544,282 / 560,249),
with no single offset. Cross-reference coordinate comparison REFUSES rather than lifts,
following the existing 435-row `designed_donor_intervals_PREcomposite_REFUSED_chrVIII`
precedent.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "CrossReferenceError",
    "ManifestError",
    "OPTIONAL_COLUMNS",
    "REFUSED_COLUMNS",
    "REQUIRED_COLUMNS",
    "ReferenceMismatchError",
    "WGSSample",
    "check_reference_consistency",
    "group_by_reference",
    "load_manifest",
    "non_portable_contigs",
    "refusal_reason",
    "require_reference_consistency",
    "require_single_reference",
]

#: Columns every manifest must carry. `reference` is per-sample by design (D1).
REQUIRED_COLUMNS: tuple[str, ...] = ("sample", "bam", "reference")

#: Columns a manifest may carry. Absent ones default to the empty string — they are
#: bookkeeping, not analysis inputs, so a missing `plate` must not fail a run.
OPTIONAL_COLUMNS: tuple[str, ...] = ("cohort", "plate", "background")

#: Schema of the refused-comparison table (§8.4), matching the existing
#: `*_REFUSED_chrVIII.tsv` precedent on Oak: chrom, start, end, reason.
REFUSED_COLUMNS: tuple[str, ...] = ("chrom", "start", "end", "reason")


class ManifestError(ValueError):
    """The manifest is malformed, incomplete, or self-contradictory."""


class ReferenceMismatchError(ManifestError):
    """A sample's declared reference disagrees with the reference its BAM was aligned to."""


class ReferenceUnobservedError(ManifestError):
    """A sample's BAM was never inventoried, so its reference could not be observed.

    Its own class, distinct from `ReferenceMismatchError`, because the two have
    different fixes and the distinction is load-bearing: a BAM aligned with a tool
    that writes no `@PG ID:bwa` line, or one merged afterwards so the line is lost,
    is UNOBSERVED, not mismatched. Collapsing them would make a legitimate pipeline
    look like a wrong-reference alignment, and a reviewer who has been told "the
    reference check raises on this cohort" would draw the wrong conclusion.
    """


class CrossReferenceError(ManifestError):
    """An operation that is only meaningful within one reference spans several."""


@dataclass(frozen=True)
class WGSSample:
    """One clonal WGS sample: the fan-out unit for every stage.

    Attributes
    ----------
    sample : str
        Sample name. Must match the BAM basename minus the configured `bam_suffix`,
        because that is the key TRACE emits (`Bc1Call.sample`) and the key every
        downstream join uses.
    bam : Path
        Coordinate-sorted BAM. Stages 1, 3 and 4 all `fetch`, so it must be indexed;
        `readers.bam_inventory` asserts the `.bai`.
    reference : Path
        The composite reference this BAM was aligned to. Per-sample by design (D1).
    cohort, plate, background : str
        Bookkeeping. `plate` is NOT a geometry promise: the v11 cohort is 96-complete
        per plate while the rearray cohort is sparse (group n as low as 1), so anything
        assuming 96 wells is right for one cohort and wrong for the other. The manifest
        is the only source of truth for per-plate n.
    """

    sample: str
    bam: Path
    reference: Path
    cohort: str = ""
    plate: str = ""
    background: str = ""


def load_manifest(path: str | Path, *, base_dir: str | Path | None = None,
                  reference_root: str | Path | None = None) -> list[WGSSample]:
    """Read a header'd TSV manifest into `WGSSample`s, failing loudly on every defect.

    Parameters
    ----------
    path : str or Path
        Manifest TSV. Header required; blank lines and lines whose first field starts
        with `#` are skipped.
    base_dir : str or Path, optional
        Directory relative `bam` / `reference` paths resolve against. Defaults to the
        manifest's own parent, so a manifest can be relocated with its data.

    Returns
    -------
    list of WGSSample
        In file order. Order is preserved because it becomes the fan-out order and
        therefore the row order of every emitted table (dual-path parity, §4).

    Raises
    ------
    ManifestError
        Missing required column, missing header, blank required value, short row, or a
        duplicate sample name. A duplicate is fatal rather than deduplicated: it would
        silently halve a cohort's denominator.
    """
    path = Path(path)
    root = Path(base_dir) if base_dir is not None else path.parent
    ref_root = Path(reference_root) if reference_root is not None else root
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ManifestError(f"{path}: empty manifest (no header row)") from None
        header = [c.strip() for c in header]
        if any(not c for c in header) or len(set(header)) != len(header):
            raise ManifestError(f"{path}: blank or duplicate header column")
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise ManifestError(
                f"{path}: manifest is missing required column(s) {missing}; "
                f"found {header}"
            )
        idx = {name: header.index(name) for name in header}
        samples: list[WGSSample] = []
        seen: dict[str, int] = {}
        for lineno, row in enumerate(reader, start=2):
            if not row or not any(v.strip() for v in row) or row[0].lstrip().startswith("#"):
                continue
            if len(row) != len(header):
                raise ManifestError(
                    f"{path}:{lineno}: row has {len(row)} fields, header has {len(header)}"
                )
            values = {name: row[i].strip() for name, i in idx.items()}
            for column in REQUIRED_COLUMNS:
                if not values[column]:
                    raise ManifestError(f"{path}:{lineno}: required column {column!r} is blank")
            name = values["sample"]
            if name in (".", "..") or "/" in name or "\\" in name:
                raise ManifestError(f"{path}:{lineno}: sample is not a safe file identity: {name!r}")
            if name in seen:
                raise ManifestError(
                    f"{path}:{lineno}: duplicate sample {name!r} (first seen on line {seen[name]})"
                )
            seen[name] = lineno
            samples.append(
                WGSSample(
                    sample=name,
                    bam=_resolve(values["bam"], root),
                    reference=_resolve(values["reference"], ref_root),
                    cohort=values.get("cohort", ""),
                    plate=values.get("plate", ""),
                    background=values.get("background", ""),
                )
            )
    if not samples:
        raise ManifestError(f"{path}: manifest has a header but no sample rows")
    return samples


def _resolve(value: str, root: Path) -> Path:
    """Resolve a manifest path against `root`, leaving absolute paths untouched."""
    p = Path(value)
    return p if p.is_absolute() else root / p


def group_by_reference(samples: Iterable[WGSSample]) -> dict[Path, list[WGSSample]]:
    """Group samples by reference path, preserving first-seen order within each group.

    Every stage runs within a group. Stage 5 may join *across* groups on `sample`, but
    must never compare coordinates across them (see `non_portable_contigs`).
    """
    groups: dict[Path, list[WGSSample]] = {}
    for s in samples:
        groups.setdefault(s.reference, []).append(s)
    return groups


def check_reference_consistency(
    samples: Iterable[WGSSample],
    observed: Mapping[str, str],
    *,
    compare_basename: bool = True,
) -> list[str]:
    """Compare each sample's declared reference to the one its BAM was aligned to.

    Parameters
    ----------
    samples : iterable of WGSSample
        The manifest.
    observed : mapping
        `{sample: reference}` as read from the BAM's own `@PG ID:bwa` line by
        `readers.bam_inventory`. A sample absent from this mapping is reported — an
        unobserved BAM is the silent-skip failure mode, not a neutral one.
    compare_basename : bool, default True
        Compare basenames rather than full paths. The `@PG` line records the path used
        at alignment time, which legitimately differs in prefix from the manifest's
        (different mount, staged scratch copy). Comparing full paths would raise false
        alarms; comparing basenames still catches a genuinely different reference,
        which is the failure this guards.

    Returns
    -------
    list of str
        One human-readable message per problem, empty when every sample agrees. Returning
        messages rather than raising lets a caller report *all* mismatches at once.
    """
    problems: list[str] = []
    for s in samples:
        seen = observed.get(s.sample)
        if seen is None:
            problems.append(f"{s.sample}: no observed reference (BAM not inventoried)")
            continue
        declared_key = Path(s.reference).name if compare_basename else str(s.reference)
        observed_key = Path(seen).name if compare_basename else str(seen)
        if declared_key != observed_key:
            problems.append(
                f"{s.sample}: manifest declares reference {declared_key!r} but the BAM "
                f"was aligned to {observed_key!r}"
            )
    return problems


def require_reference_consistency(
    samples: Iterable[WGSSample],
    observed: Mapping[str, str],
    *,
    compare_basename: bool = True,
) -> None:
    """`check_reference_consistency`, raising on any problem.

    Raises
    ------
    ReferenceUnobservedError
        At least one sample has no observed reference. Checked FIRST and raised as
        its own class: `observed` is built from the BAM's `@PG ID:bwa` line, which a
        BAM aligned by another tool — or merged after alignment — legitimately lacks.
        That is "could not check", not "checked and wrong", and conflating the two
        would make a clean cohort look mis-aligned (NN 251, defect 14).
    ReferenceMismatchError
        Every sample was observed and at least one disagrees with the manifest.
    """
    problems = check_reference_consistency(samples, observed, compare_basename=compare_basename)
    if not problems:
        return
    unobserved = [p for p in problems if "no observed reference" in p]
    mismatched = [p for p in problems if p not in set(unobserved)]
    # 🔴 A REAL MISMATCH IS REPORTED FIRST. In a mixed cohort — some BAMs with a `@PG`
    # line, some without — raising Unobserved first would hide a genuine wrong-reference
    # alignment behind a "could not check" message, which is the louder-but-wronger
    # failure. Mismatch is the one that means data is wrong; unobserved only means the
    # check could not run.
    if mismatched:
        raise ReferenceMismatchError(
            "manifest/BAM reference mismatch for "
            f"{len(mismatched)} sample(s):\n  " + "\n  ".join(mismatched)
            + (
                f"\n  (plus {len(unobserved)} sample(s) whose reference could not be "
                "observed at all; fix the mismatch first, then re-run to see those.)"
                if unobserved else ""
            )
        )
    if unobserved:
        raise ReferenceUnobservedError(
            f"{len(unobserved)} sample(s) have no observed reference, so the "
            "manifest/BAM reference check could not run:\n  "
            + "\n  ".join(unobserved)
            + "\n  (A BAM with no `@PG ID:bwa` line is UNOBSERVED, not mismatched — "
            "inventory it with `readers.bam_inventory`, or compare `@SQ` lengths with "
            "`incomparable_contigs` instead.)"
        )
    raise ReferenceMismatchError(  # pragma: no cover - unreachable; kept as a backstop
        "manifest/BAM reference mismatch for "
        f"{len(problems)} sample(s):\n  " + "\n  ".join(problems)
    )


def incomparable_contigs(
    lengths_a: Mapping[str, int], lengths_b: Mapping[str, int]
) -> list[str]:
    """Contigs present in BOTH references whose LENGTHS disagree.

    🔴 **NN 251 correction to NN 243.** NN 243 listed this function among the
    "guards that are written, documented and called by nothing". It was in fact
    **never written** — `grep` over the package finds only prose references to it.
    So this is a new guard, not a wiring job, and the audit's "4 free wins" are
    really 3 wirings plus 1 implementation. Recorded here rather than silently
    fixed, because a reviewer told the guard already existed would not re-check it.

    The comparison is deliberately restricted to the INTERSECTION of names: a
    composite reference legitimately adds contigs (`native_2micron`,
    `integrated_plasmid_<oligo>`) that a plain reference lacks, and treating those
    as a problem would make every correctly-built composite fail. A shared name
    whose LENGTH differs, though, is never benign — it is the chrVIII/CUP1 case,
    where the same name denotes a contracted and an uncontracted contig and **no
    single coordinate offset exists between them**, so every coordinate read in the
    other frame is silently wrong.

    Returns
    -------
    list of str
        One message per disagreeing contig, empty when every shared name agrees.
        Messages rather than a raise, matching `check_reference_consistency`, so a
        caller can report them all at once.
    """
    shared = set(lengths_a) & set(lengths_b)
    # 🔴 A FULLY DISJOINT PAIR IS A PROBLEM, NOT A PASS. Restricting the comparison to
    # the intersection is right (a composite reference adds contigs), but with an EMPTY
    # intersection the loop below finds nothing and the guard would return "all clear"
    # for the worst case of all — a `chr15`-named BAM against `chrXV`-named references,
    # which is the very naming mismatch defect 2 was about. Reported explicitly.
    if not shared and (lengths_a or lengths_b):
        return [
            "no contig name is shared between the two references "
            f"({len(lengths_a)} vs {len(lengths_b)} contigs; first names "
            f"{sorted(lengths_a)[:3]} vs {sorted(lengths_b)[:3]}) — they are not the "
            "same naming convention (chr15 vs chrXV), so no coordinate is comparable"
        ]
    problems: list[str] = []
    for name in sorted(shared):
        a, b = int(lengths_a[name]), int(lengths_b[name])
        if a != b:
            problems.append(f"{name}: length {a} vs {b}")
    return problems


def require_single_reference(samples: Sequence[WGSSample], *, operation: str) -> Path:
    """Return the one reference these samples share, or refuse.

    Any operation that compares or merges coordinates is only meaningful within a single
    reference (§8.4). Call this at the top of such an operation rather than discovering
    the problem as a plausible-looking number.
    """
    references = {str(s.reference) for s in samples}
    if len(references) != 1:
        raise CrossReferenceError(
            f"{operation}: refusing — samples span {len(references)} references "
            f"{sorted(references)}. Coordinates are not portable across references "
            "(chrVIII is CUP1-contracted; no single offset exists). Run per reference "
            "group, or record the refusal."
        )
    return Path(next(iter(references)))


def non_portable_contigs(
    lengths_a: Mapping[str, int],
    lengths_b: Mapping[str, int],
) -> list[str]:
    """Contigs whose coordinates cannot be compared between two references.

    Parameters
    ----------
    lengths_a, lengths_b : mapping
        `{contig: length}` for each reference, as read from its `.fai` by
        `readers.bam_inventory`.

    Returns
    -------
    list of str
        Sorted contigs that differ in length or are absent from one reference. A contig
        with equal length in both is *not* proof of portability — equal length with
        internal rearrangement would pass — so this is a necessary, not sufficient, test
        and callers should refuse rather than lift on any hit.

    Notes
    -----
    Measured on the three references in use: chrVIII is 568,236 / 544,282 / 560,249
    (the CUP1 contraction) while chrIII is identical at 317,783 across all three. Contig
    *names* drift too — one reference carries `unplaced_contig_34` where another carries
    `34a` + `34b` — which is why absence counts as non-portable here.
    """
    out = set()
    for contig, length in lengths_a.items():
        if contig not in lengths_b or lengths_b[contig] != length:
            out.add(contig)
    out.update(contig for contig in lengths_b if contig not in lengths_a)
    return sorted(out)


def refusal_reason(contig: str, reference_a: str, reference_b: str) -> str:
    """The `reason` cell for a refused cross-reference row (`REFUSED_COLUMNS`)."""
    return (
        f"{contig} differs between {Path(reference_a).name} and {Path(reference_b).name}; "
        "no single offset exists"
    )
