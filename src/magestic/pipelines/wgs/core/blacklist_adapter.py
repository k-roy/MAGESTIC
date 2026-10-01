"""magestic.pipelines.wgs.core.blacklist_adapter — refuse an inert blacklist (defect 2).

THE DEFECT
----------
`step_3_sv_plasmid._run_half_b` checks the REDI blacklist with `Path(blacklist).exists()`
and nothing else, then hands it to `trace_crispr.wgs.sv_plasmid.load_blacklist`, whose
parsing is **deliberately permissive** (it reproduces Shengdi's Perl `split("\\t")`,
including coercing a non-numeric coordinate to `0`, which makes that row inert).

Measured on the real inputs: **one copy of `REDI_reads_mapped_regions.txt` is
`chr15`-named while the composite references are `chrXV`-named, so 0 of its 9 contig
names exist in the reference and the blacklist is 100% inert.** A second copy is in
PRE coordinates. One row is space-separated rather than tab-separated and so parses
to the interval `(0, 0)`, which can never match.

Nothing raised. `clean_donor_reads` ran, reported frag counts, and the output looked
filtered. 🔴 **A blacklist that matches nothing is worse than no blacklist**, because
the run carries the authority of having been filtered.

WHAT THIS MODULE ADDS, AND WHAT IT CANNOT CATCH
-----------------------------------------------
`audit_blacklist` parses STRICTLY and counts every reject against a stated
denominator; `require_blacklist_applicable` raises when the blacklist cannot apply to
the reference at all.

⚠️ **Stated plainly: a same-name blacklist in the wrong COORDINATE FRAME, whose
intervals still fall inside the contig, remains undetectable here.** The name check
catches the `chr15`/`chrXV` copy; the bounds check catches a PRE-coordinate interval
only when it runs off the end of the contig. A 200 bp frame shift inside a 1 Mb
contig passes every check in this module. The only real defence is to version the
blacklist with the reference it was built against — recorded as text, per the
no-symlinks policy, in the screen's `PROJECT.md`.

This is an adapter: it reads the same file TRACE reads and reports on it. TRACE's
permissive parse is faithful to the Perl and should NOT be changed — the parity it
buys is real. THE TRACE-SIDE FIX is additive: `load_blacklist` should return a
`(Blacklist, LoadReport)` pair rather than a bare dict, so a caller can see the
rejects without re-parsing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "BlacklistAudit",
    "BlacklistInapplicableError",
    "audit_blacklist",
    "reference_contig_lengths",
    "require_blacklist_applicable",
]


class BlacklistInapplicableError(ValueError):
    """The blacklist cannot apply to this reference, so filtering would be a fiction."""


@dataclass(frozen=True)
class BlacklistAudit:
    """Every count, against a stated denominator.

    Attributes
    ----------
    path : Path
    n_rows : int
        Non-blank rows read. The denominator for everything below.
    n_well_formed : int
        Rows with >= 3 tab-separated fields whose fields 2 and 3 are integers.
    n_malformed : int
        `n_rows - n_well_formed`. These are INERT in TRACE's parse (coordinates
        coerce to 0), not errors — counted so "inert" is never invisible.
    n_space_separated : int
        Rows with no tab but with whitespace: the known `yT177_...` row's shape.
    n_present_contigs, n_absent_contigs : int
        Distinct contig names in the blacklist that are / are not in the reference.
    n_applicable_rows : int
        Well-formed rows on a present contig whose interval is within its length.
        **This is the number that decides whether filtering happened.**
    n_out_of_bounds : int
        Rows on a present contig whose interval exceeds the contig length — the
        wrong-coordinate-frame signal.
    fraction_applicable : float
        `n_applicable_rows / n_rows`, or 0.0 for an empty file.
    """

    path: Path
    n_rows: int = 0
    n_well_formed: int = 0
    n_malformed: int = 0
    n_space_separated: int = 0
    n_applicable_rows: int = 0
    n_out_of_bounds: int = 0
    present_contigs: tuple[str, ...] = ()
    absent_contigs: tuple[str, ...] = ()

    @property
    def n_present_contigs(self) -> int:
        return len(self.present_contigs)

    @property
    def n_absent_contigs(self) -> int:
        return len(self.absent_contigs)

    @property
    def fraction_applicable(self) -> float:
        return (self.n_applicable_rows / self.n_rows) if self.n_rows else 0.0

    def summary(self) -> str:
        return (
            f"{self.path.name}: {self.n_applicable_rows} of {self.n_rows} row(s) can "
            f"apply ({self.fraction_applicable:.1%}); contigs present "
            f"{self.n_present_contigs}, absent {self.n_absent_contigs} "
            f"(first absent: {', '.join(self.absent_contigs[:5]) or 'none'}); "
            f"malformed {self.n_malformed} (space-separated {self.n_space_separated}); "
            f"out-of-bounds {self.n_out_of_bounds}"
        )


def reference_contig_lengths(fai_or_fasta: str | Path) -> dict[str, int]:
    """`{contig: length}` from a `.fai` (or a FASTA whose `.fai` sits beside it).

    Reads the index rather than the FASTA, so this is cheap on a composite reference
    and needs no pysam handle.
    """
    p = Path(fai_or_fasta)
    fai = p if p.suffix == ".fai" else p.with_name(p.name + ".fai")
    if not fai.exists():
        raise BlacklistInapplicableError(
            f"no FASTA index at {fai}; a blacklist cannot be checked against a "
            "reference whose contig names and lengths are unknown"
        )
    out: dict[str, int] = {}
    for line in fai.read_text().splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) >= 2 and cols[1].strip().isdigit():
            out[cols[0].strip()] = int(cols[1])
    return out


def audit_blacklist(
    path: str | Path, reference_contigs: dict[str, int]
) -> BlacklistAudit:
    """Parse the blacklist strictly and report every reject.

    `reference_contigs` is `{contig: length}` — from `reference_contig_lengths`, or
    from `pysam.AlignmentFile(...).header` / a `.fai` the caller already has.
    """
    path = Path(path)
    n_rows = n_well = n_space = n_applicable = n_oob = 0
    seen: dict[str, bool] = {}
    for raw in path.read_text().splitlines():
        line = raw.rstrip("\r")
        if not line.strip():
            continue
        n_rows += 1
        cols = line.split("\t")
        if len(cols) < 2 and " " in line.strip():
            n_space += 1
        chrom = cols[0].strip()
        seen.setdefault(chrom, chrom in reference_contigs)
        well = (
            len(cols) >= 3
            and cols[1].strip().lstrip("-").isdigit()
            and cols[2].strip().lstrip("-").isdigit()
        )
        if not well:
            continue
        n_well += 1
        start, end = int(cols[1]), int(cols[2])
        length = reference_contigs.get(chrom)
        if length is None:
            continue
        if start < 0 or end > length:
            n_oob += 1
            continue
        n_applicable += 1
    return BlacklistAudit(
        path=path,
        n_rows=n_rows,
        n_well_formed=n_well,
        n_malformed=n_rows - n_well,
        n_space_separated=n_space,
        n_applicable_rows=n_applicable,
        n_out_of_bounds=n_oob,
        present_contigs=tuple(sorted(c for c, ok in seen.items() if ok)),
        absent_contigs=tuple(sorted(c for c, ok in seen.items() if not ok)),
    )


def require_blacklist_applicable(
    path: str | Path,
    reference_contigs: dict[str, int],
    *,
    min_fraction_applicable: float = 0.0,
    strict_rows: bool = False,
) -> BlacklistAudit:
    """`audit_blacklist`, raising when the blacklist cannot do its job.

    Parameters
    ----------
    min_fraction_applicable : float, default 0.0
        The audit must EXCEED this fraction. The default refuses only a blacklist
        that is **entirely** inert (the measured `chr15` case), which is the
        bright-line failure; a screen may raise the bar.
    strict_rows : bool, default False
        Also refuse any malformed row. Default False **on purpose**: the real
        `REDI_reads_mapped_regions.txt` contains one known space-separated row that
        the Perl treats as inert, and refusing it would block every existing cohort
        over a row that has never mattered. It is counted either way.

    Raises
    ------
    BlacklistInapplicableError
    """
    audit = audit_blacklist(path, reference_contigs)
    if audit.n_rows == 0:
        raise BlacklistInapplicableError(
            f"{audit.path} is empty — 0 of 0 rows can apply. An empty blacklist makes "
            "the run look filtered while filtering nothing."
        )
    if audit.fraction_applicable <= min_fraction_applicable:
        raise BlacklistInapplicableError(
            f"blacklist is inert against this reference: {audit.summary()}. "
            f"{audit.n_applicable_rows} of {audit.n_rows} rows can apply. A blacklist "
            "matching nothing is worse than no blacklist, because the run looks "
            "filtered. Check the contig-naming convention (chr15 vs chrXV) and the "
            "coordinate frame (PRE vs POST) before re-running."
        )
    if strict_rows and audit.n_malformed:
        raise BlacklistInapplicableError(
            f"{audit.n_malformed} of {audit.n_rows} blacklist rows are malformed and "
            f"therefore inert: {audit.summary()}"
        )
    return audit
