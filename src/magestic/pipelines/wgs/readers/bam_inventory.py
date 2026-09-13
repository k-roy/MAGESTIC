"""
magestic.pipelines.wgs.readers.bam_inventory — BAM discovery and the D1 cross-check.

BAMs are enumerated from the MANIFEST, never from a directory glob. That is not a style
preference: TRACE's batch entry point globs `bam_dir + "/*" + bam_suffix` with a default
suffix (`.sort.dedup.recal2.bam`) that matches none of this program's BAMs, and an empty
glob writes a header-only TSV and returns 0. Taking the BAM list from the manifest makes
that failure unreachable rather than merely guarded.

This module reads the observations; `core.manifest` owns the comparison. Keeping the two
apart is what makes the D1 reference cross-check unit-testable without a BAM: the reader
hands `core.manifest.check_reference_consistency` a plain `{sample: reference}` mapping.

Author: Kevin R. Roy
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from pathlib import Path

from ..core.manifest import WGSSample

__all__ = [
    "BamInventoryError",
    "index_path",
    "observed_references",
    "reference_contig_lengths",
    "require_indexed",
]

#: `@PG` lines record the command line; the reference is the last `.fa`/`.fasta` token.
_FASTA_TOKEN = re.compile(r"\S+\.(?:fa|fasta|fna)(?:\.gz)?\b")


class BamInventoryError(RuntimeError):
    """A BAM is missing, unindexed, or carries no usable alignment provenance."""


def index_path(bam: str | Path) -> Path | None:
    """The htslib index beside `bam`, or None.

    htslib resolves `<bam>.bai` / `<bam>.csi` — NOT the picard-style `<stem>.bai`, so
    accepting the latter would pass a precondition check and then fail inside `fetch`.
    """
    bam = Path(bam)
    for ext in (".bai", ".csi"):
        candidate = bam.with_name(bam.name + ext)
        if candidate.exists():
            return candidate
    return None


def require_indexed(samples: Iterable[WGSSample]) -> list[Path]:
    """Assert every sample's BAM exists and is indexed; return the BAM paths.

    Stages 1, 3 and 4 all reach the BAM through `fetch`/`pileup`, which require an
    index. TRACE documents neither requirement, so without this a missing `.bai`
    surfaces as an opaque pysam error deep inside a per-sample loop.
    """
    bams: list[Path] = []
    problems: list[str] = []
    for s in samples:
        bam = Path(s.bam)
        if not bam.exists():
            problems.append(f"{s.sample}: BAM not found: {bam}")
            continue
        if index_path(bam) is None:
            problems.append(f"{s.sample}: no {bam.name}.bai or .csi beside {bam}")
            continue
        bams.append(bam)
    if problems:
        raise BamInventoryError(
            f"{len(problems)} BAM precondition failure(s):\n  " + "\n  ".join(problems)
        )
    return bams


def observed_references(samples: Iterable[WGSSample]) -> dict[str, str]:
    """`{sample: reference}` read from each BAM's own `@PG ID:bwa` line.

    This is the observation `core.manifest.check_reference_consistency` compares against
    the manifest's declared reference (decision D1). A BAM with no parseable `@PG`
    reference is simply ABSENT from the returned mapping — and `check_reference_consistency`
    reports an absent sample as a problem, so an unreadable header is never silently
    treated as agreement.
    """
    import pysam

    out: dict[str, str] = {}
    for s in samples:
        bam = Path(s.bam)
        if not bam.exists():
            continue
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            ref = _reference_from_header(handle.header.to_dict())
        if ref:
            out[s.sample] = ref
    return out


def _reference_from_header(header: dict) -> str | None:
    """Pull the reference FASTA out of a BAM header's `@PG` records.

    Prefers `ID:bwa` (the aligner that recorded the reference for this program's BAMs)
    and falls back to any `@PG` whose command line names a FASTA.
    """
    programs: Sequence[dict] = header.get("PG", []) or []
    ordered = [pg for pg in programs if str(pg.get("ID", "")).startswith("bwa")] + list(programs)
    for pg in ordered:
        match = _FASTA_TOKEN.search(str(pg.get("CL", "")))
        if match:
            return match.group(0)
    return None


def reference_contig_lengths(reference: str | Path) -> dict[str, int]:
    """`{contig: length}` from a reference's `.fai`.

    Feeds `core.manifest.non_portable_contigs`, which refuses cross-reference coordinate
    comparison — chrVIII is CUP1-contracted and differs across the references in use
    (568,236 / 544,282 / 560,249) with no single offset.
    """
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        raise BamInventoryError(
            f"{reference}: no .fai beside the reference. Index it with `samtools faidx`; "
            "contig lengths are needed for the coordinate-portability refusal (§8.4)."
        )
    lengths: dict[str, int] = {}
    for line in fai.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) >= 2:
            lengths[fields[0]] = int(fields[1])
    return lengths
