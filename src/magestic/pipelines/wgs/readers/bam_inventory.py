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
import shlex
import hashlib
from functools import lru_cache
from collections.abc import Iterable, Sequence
from pathlib import Path

from ..core.manifest import WGSSample, ReferenceMismatchError, require_reference_consistency

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

    Includes the stem.bai/stem.csi alternatives recognized by htslib. Existence
    is only discovery; require_indexed also opens the index through htslib.
    """
    bam = Path(bam)
    for ext in (".bai", ".csi"):
        for candidate in (bam.with_name(bam.name + ext), bam.with_suffix(ext)):
            if candidate.is_file():
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
        try:
            import pysam
            with pysam.AlignmentFile(str(bam), "rb") as handle:
                handle.check_index()
        except (ValueError, OSError) as exc:
            problems.append(f"{s.sample}: BAM/index cannot be opened: {exc}")
            continue
        bams.append(bam)
    if problems:
        raise BamInventoryError(
            f"{len(problems)} BAM precondition failure(s):\n  " + "\n  ".join(problems)
        )
    return bams


@lru_cache(maxsize=8)
def _reference_md5(path: str, identity: tuple) -> dict[str, str]:
    import pysam
    with pysam.FastaFile(path) as fa:
        return {name: hashlib.md5(fa.fetch(name).upper().encode()).hexdigest()
                for name in fa.references}


def validate_bam_reference(sample: WGSSample) -> None:
    """Refuse contradictory or absent alignment reference evidence.

    Contig dictionaries must match. SQ M5 values, when present, must match the
    actual FASTA. Complete M5 evidence also permits renamed/staged references.
    Otherwise a parsed PG reference basename must agree. Basename agreement is
    a provenance check, not a sequence-identity proof when M5 is absent.
    """
    import pysam
    reference = Path(sample.reference)
    if not reference.is_file() or not Path(str(reference)+".fai").is_file():
        raise BamInventoryError(f"{sample.sample}: reference FASTA/.fai missing: {reference}")
    with pysam.AlignmentFile(str(sample.bam), "rb") as bam:
        header = bam.header.to_dict()
        observed = dict(zip(bam.references, bam.lengths))
    with pysam.FastaFile(str(reference)) as fa:
        declared = dict(zip(fa.references, fa.lengths))
    if observed != declared:
        raise ReferenceMismatchError(f"{sample.sample}: BAM SQ dictionary differs from {reference}")
    m5 = {sq['SN']: str(sq['M5']).lower() for sq in header.get('SQ', []) if sq.get('M5')}
    if m5:
        stat = reference.stat()
        hashes = _reference_md5(str(reference.resolve()),
                                (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns))
        wrong = [name for name, value in m5.items() if hashes.get(name) != value]
        if wrong:
            raise ReferenceMismatchError(f"{sample.sample}: reference sequence M5 mismatch: {wrong}")
        if len(m5) == len(declared):
            return
    ref = _reference_from_header(header)
    require_reference_consistency([sample], {sample.sample: ref} if ref else {})


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
        try:
            tokens = shlex.split(str(pg.get("CL", "")))
        except ValueError:
            continue
        for token in tokens:
            candidate = token.split("=", 1)[-1]
            if re.search(r"\.(?:fa|fasta|fna)(?:\.gz)?$", candidate):
                return candidate
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
