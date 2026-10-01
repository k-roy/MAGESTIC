"""magestic.pipelines.wgs.core.gvcf_adapter — expand gVCF reference BLOCKS to positions.

WHY THIS IS AN ADAPTER AND NOT AN ALGORITHM CHANGE (defect 5)
--------------------------------------------------------------
`trace_crispr.wgs.donor_outcome.precall_sample`'s coverage test is

    if v_chr == chrom and ref_depth >= ne_threshold and qwin_start <= v_pos <= qwin_end:
        flag_covered += 1

— it compares only the record's **POS**. A *compressed* gVCF (what GATK
`-ERC BP_RESOLUTION`'s opposite, `-ERC GVCF`, and `bcftools mpileup -g` both emit)
represents a run of equally-covered reference bases as ONE record carrying
`INFO END=<last position>`. So a single reference block that covers the entire
quantification window contributes `flag_covered == 1`, the test
`flag_covered >= quantification_window_size * quantification_window_cov_perc`
fails, and a perfectly well-covered unedited cell is reported `unknown`.

Measured on a synthetic fixture (`tests/.../test_251_trace_fixes.py`): donor window
101-160, quantification window [126,135], one block `POS=101 END=160 AD=30,0` ->
`unknown`. After expansion -> `designed_variant_unedited`.

**Expanding a block to its positions derives no new quantity.** The `END` key is
already in the record and the per-position depth is already stated by that record's
own `AD`/`DP`; writing one line per covered position is a pure *re-representation*,
the same class of work `vcf_adapter.reorder_format_field` already does. It is
therefore magestic's to own while TRACE is pinned.

THE TRACE-SIDE FIX (do this when the pin moves)
-----------------------------------------------
In `precall_sample`'s gVCF pass, read `END` from column 7 and count the overlap
instead of the single position::

    end = v_pos
    for kv in cols[7].split(";"):
        if kv.startswith("END="):
            end = int(kv[4:])
    if v_chr == chrom and ref_depth >= ne_threshold:
        flag_covered += max(0, min(end, qwin_end) - max(v_pos, qwin_start) + 1)

That is 5 lines and removes the need for this module. Until then, every gVCF handed
to Stage 2 must pass through `expand_reference_blocks`.

NO SILENT DROP
--------------
`ExpandResult` states every denominator: records in, records out, blocks expanded,
positions added, and the longest block seen. A block longer than
`max_expand_positions` RAISES rather than being truncated, because a truncated block
would under-report coverage and so manufacture an `unknown`, which is the very
failure this module exists to remove.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "DEFAULT_MAX_EXPAND_POSITIONS",
    "ExpandResult",
    "GvcfAdapterError",
    "block_end",
    "expand_reference_blocks",
    "window_reference_coverage",
]

#: A single block wider than this raises. A whole-chromosome block (which some
#: callers emit for an unmapped contig) would otherwise expand to millions of lines.
DEFAULT_MAX_EXPAND_POSITIONS = 1_000_000


class GvcfAdapterError(ValueError):
    """The gVCF is malformed, or a reference block is implausibly wide."""


@dataclass
class ExpandResult:
    """What `expand_reference_blocks` did, with every denominator stated.

    Attributes
    ----------
    src, dst : Path
        Input and output.
    n_records_in : int
        Data records read (the denominator for everything below).
    n_records_out : int
        Data records written. Equals `n_records_in` when no record carried `END`.
    n_blocks_expanded : int
        Records whose `END` exceeded their `POS`. **Zero means this adapter was a
        no-op and the gVCF was already per-position** — not that it failed.
    n_positions_added : int
        `n_records_out - n_records_in`.
    max_block_width : int
        Widest block seen, in positions. 1 when nothing was a block.
    n_malformed_end : int
        Records with an unparseable `END=` value. Counted, and the record is written
        through unchanged at its own POS rather than discarded.
    """

    src: Path
    dst: Path
    n_records_in: int = 0
    n_records_out: int = 0
    n_blocks_expanded: int = 0
    n_positions_added: int = 0
    max_block_width: int = 1
    n_malformed_end: int = 0


def _open_text(path: Path):
    """Open a plain or bgzipped/gzipped gVCF as text (bgzip is gzip-compatible)."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def block_end(info_field: str, pos: int) -> tuple[int, bool]:
    """`(end, malformed)` for one INFO column.

    Returns `pos` when there is no `END` key, which is the per-position case. A
    present-but-unparseable `END` returns `(pos, True)` so the caller can count it:
    silently treating it as per-position would under-report coverage.
    """
    for kv in str(info_field).split(";"):
        if kv.startswith("END="):
            try:
                return int(kv[4:]), False
            except ValueError:
                return pos, True
    return pos, False


def expand_reference_blocks(
    src: str | Path,
    dst: str | Path,
    *,
    restrict: tuple[str, int, int] | None = None,
    max_expand_positions: int = DEFAULT_MAX_EXPAND_POSITIONS,
) -> ExpandResult:
    """Rewrite `src` so every reference block becomes one record per covered position.

    Parameters
    ----------
    src : str or Path
        The gVCF. Plain, gzipped or bgzipped.
    dst : str or Path
        Output. **Always plain text**, because TRACE opens it with a bare `open()`.
    restrict : (chrom, start, end), optional
        Expand only within this 1-based inclusive span, and pass every other record
        through at its own POS. Use it to bound the output to the donor/quantification
        window, which is all `precall_sample` reads. Records outside are NOT dropped.
    max_expand_positions : int
        A single block wider than this raises `GvcfAdapterError`.

    Returns
    -------
    ExpandResult

    Notes
    -----
    A record with no `END` is written **byte-identically**, so a per-position gVCF
    round-trips unchanged and Stage 2's numbers cannot move on an already-expanded
    input. The `END` key is stripped from an expanded record's INFO (replaced by
    `.` if it was the only key) because a per-position record carrying `END=` from
    its parent block would be self-contradictory.
    """
    src = Path(src)
    dst = Path(dst)
    res = ExpandResult(src=src, dst=dst)
    dst.parent.mkdir(parents=True, exist_ok=True)

    with _open_text(src) as fh, dst.open("w", newline="\n") as out:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                out.write(line + "\n")
                continue
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                raise GvcfAdapterError(
                    f"{src}: record {line[:60]!r} has {len(cols)} columns, "
                    "a gVCF data record needs at least 8"
                )
            res.n_records_in += 1
            try:
                pos = int(cols[1])
            except ValueError as exc:
                raise GvcfAdapterError(f"{src}: non-integer POS in {line[:60]!r}") from exc

            end, malformed = block_end(cols[7], pos)
            if malformed:
                res.n_malformed_end += 1
            if end <= pos:
                out.write(line + "\n")  # per-position: byte-identical passthrough
                res.n_records_out += 1
                continue

            width = end - pos + 1
            res.max_block_width = max(res.max_block_width, width)
            if width > max_expand_positions:
                raise GvcfAdapterError(
                    f"{src}: reference block {cols[0]}:{pos}-{end} spans {width} "
                    f"positions, above max_expand_positions={max_expand_positions}. "
                    "Truncating it would under-report coverage and manufacture an "
                    "`unknown`; pass `restrict=` to bound the region instead."
                )

            lo, hi = pos, end
            if restrict is not None:
                r_chrom, r_start, r_end = restrict
                if cols[0] != r_chrom:
                    out.write(line + "\n")
                    res.n_records_out += 1
                    continue
                lo, hi = max(pos, r_start), min(end, r_end)
                if lo > hi:
                    out.write(line + "\n")
                    res.n_records_out += 1
                    continue

            res.n_blocks_expanded += 1
            info = ";".join(
                kv for kv in cols[7].split(";") if kv and not kv.startswith("END=")
            ) or "."
            head = cols[:7]
            tail = cols[8:]
            for p in range(lo, hi + 1):
                row = head[:1] + [str(p)] + head[2:] + [info] + tail
                out.write("\t".join(row) + "\n")
                res.n_records_out += 1

    res.n_positions_added = res.n_records_out - res.n_records_in
    return res


def window_reference_coverage(
    gvcf_path: str | Path,
    *,
    chrom: str,
    window_start: int,
    window_end: int,
    ne_threshold: int = 1,
    ad_slot: int = 2,
) -> int:
    """Count DISTINCT positions in `[window_start, window_end]` covered at depth
    `>= ne_threshold`, honouring `END`.

    This is the quantity `precall_sample`'s `flag_covered` is trying to be, computed
    from the same fields (`AD`'s first value at FORMAT slot `ad_slot`, and `END`).
    It is used by `outcome_reason.explain_outcome` to say *why* a sample is
    uncallable, and the test suite asserts it never contradicts TRACE's own flag.

    Counting DISTINCT positions rather than records also fixes a second-order
    hazard: an already-expanded gVCF double-counted a position would inflate
    `flag_covered` above the window size.
    """
    covered: set[int] = set()
    path = Path(gvcf_path)
    if not path.exists():
        return 0
    with _open_text(path) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 10 or cols[0] != chrom:
                continue
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            info = cols[9].split(":")
            index = ad_slot - 1
            if len(info) <= index:
                continue
            try:
                ref_depth = int(info[index].split(",")[0])
            except ValueError:
                continue
            if ref_depth < ne_threshold:
                continue
            end, _ = block_end(cols[7], pos)
            lo, hi = max(pos, window_start), min(max(end, pos), window_end)
            if lo <= hi:
                covered.update(range(lo, hi + 1))
    return len(covered)
