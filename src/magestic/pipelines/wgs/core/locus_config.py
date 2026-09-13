"""
magestic.pipelines.wgs.core.locus_config — the flat locus sidecar TRACE can actually parse.

TRACE's `bc1_consensus.load_locus_config` is documented as reading "the magestic
`input.yaml` subset", but it is **not a YAML parser**: it splits each line on the first
`:`, strips quotes, and ignores nesting entirely. Handed a real nested YAML it does not
error — it silently produces wrong track configs, which is the worst available outcome
for a barcode caller. So magestic emits a **flat key/value sidecar** derived from the
screen config, and TRACE parses that. Writing it is orchestration; parsing it stays in
TRACE so there is exactly one parser.

Two TRACE behaviors this module works around, both measured:

* `tracks_from_config` unconditionally indexes the three `plasmid_bc_*` keys, so a
  genome-only sidecar raises `KeyError` inside TRACE. `write_locus_sidecar` therefore
  requires all six keys and names the missing ones.
* `bc1_consensus.run()` validates the plasmid track against the reference even when
  `use_plasmid=False`, so a reference lacking the plasmid contig fails a genome-only
  run. `validate_tracks` takes the tracks to validate, so the caller chooses.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from . import _trace

if TYPE_CHECKING:  # pragma: no cover - typing only
    import pysam
    from trace_crispr.wgs.bc1_consensus import TrackConfig

__all__ = [
    "LOCUS_CONFIG_KEYS",
    "LocusConfigError",
    "build_tracks",
    "format_locus",
    "load_locus_sidecar",
    "parse_window",
    "validate_tracks",
    "write_locus_sidecar",
]

#: Every key `tracks_from_config` indexes. All six are required — see module docstring.
LOCUS_CONFIG_KEYS: tuple[str, ...] = (
    "barcode_locus",
    "barcode_left_sequence",
    "barcode_right_sequence",
    "plasmid_bc_locus",
    "plasmid_bc_left_sequence",
    "plasmid_bc_right_sequence",
)

_SIDECAR_BANNER = (
    "# magestic.pipelines.wgs locus sidecar — flat key/value, one pair per line.\n"
    "# Consumed by trace_crispr.wgs.bc1_consensus.load_locus_config, which is NOT a YAML\n"
    "# parser: it splits on the first ':' and ignores nesting. Do not add structure here.\n"
)


class LocusConfigError(ValueError):
    """The locus sidecar is missing keys or is not in the flat form TRACE can parse."""


def format_locus(contig: str, start: int, end: int) -> str:
    """Render a `contig:start-end` locus string (1-based inclusive, as TRACE parses it)."""
    return f"{contig}:{start}-{end}"


def parse_window(value: str | None) -> tuple[int, int] | None:
    """Parse an amplicon-telltale window `"10224-10235"` into a 1-based inclusive tuple.

    Mirrors TRACE's `_parse_window`, which is private and therefore not API. Kept to three
    lines rather than imported so a TRACE refactor of a private helper cannot break us.
    """
    if not value:
        return None
    a, b = str(value).split("-")
    return int(a), int(b)


def write_locus_sidecar(path: str | Path, values: Mapping[str, str]) -> Path:
    """Write the flat locus sidecar, requiring all of `LOCUS_CONFIG_KEYS`.

    Parameters
    ----------
    path : str or Path
        Destination. Overwritten if present (steps are idempotent, §6).
    values : mapping
        Must supply all six `LOCUS_CONFIG_KEYS`. Extra keys are written through, since
        `load_locus_config` ignores what it does not look up.

    Returns
    -------
    Path
        The path written.

    Raises
    ------
    LocusConfigError
        A required key is missing or blank, or a value contains a newline (which would
        split into two bogus config lines). Failing here is the whole point: the
        alternative is TRACE silently building a track from a mis-parsed value.

    Notes
    -----
    Output is deterministic — fixed key order, no timestamp — because both orchestration
    paths write it and §4 requires byte-identical outputs.

    The locus coordinates come from the screen config, NOT from the reference `.fai`: a
    `.fai` carries contig names and lengths only, so the position of the barcode N-run
    cannot be derived from it. The reference's role here is *validation*
    (`validate_tracks`), not derivation.
    """
    missing = [k for k in LOCUS_CONFIG_KEYS if not str(values.get(k, "")).strip()]
    if missing:
        raise LocusConfigError(
            f"locus sidecar is missing required key(s) {missing}. All six are required even "
            "for a genome-only run: trace_crispr.wgs.bc1_consensus.tracks_from_config "
            "indexes the plasmid_bc_* keys unconditionally."
        )
    ordered = list(LOCUS_CONFIG_KEYS) + [k for k in values if k not in LOCUS_CONFIG_KEYS]
    lines = []
    for key in ordered:
        value = str(values[key]).strip()
        if "\n" in value or "\r" in value:
            raise LocusConfigError(f"locus sidecar value for {key!r} contains a newline")
        lines.append(f"{key}: {value}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="\n") as fh:
        fh.write(_SIDECAR_BANNER)
        fh.write("\n".join(lines) + "\n")
    return path


def load_locus_sidecar(path: str | Path) -> dict[str, str]:
    """Read the sidecar back through TRACE's own parser, then require the six keys.

    Round-tripping through `load_locus_config` rather than a local parser means the
    file is validated against the code that will actually consume it.
    """
    cfg = dict(_trace.bc1_consensus().load_locus_config(str(path)))
    missing = [k for k in LOCUS_CONFIG_KEYS if not cfg.get(k)]
    if missing:
        raise LocusConfigError(f"{path}: sidecar parsed but missing key(s) {missing}")
    return cfg


def build_tracks(
    sidecar: str | Path | Mapping[str, str],
    *,
    amplicon_left_window: str | tuple[int, int] | None = None,
    amplicon_right_window: str | tuple[int, int] | None = None,
) -> tuple[TrackConfig, TrackConfig]:
    """Build the `(genome, plasmid)` TRACE tracks from a sidecar path or parsed mapping.

    Parameters
    ----------
    sidecar : str, Path or mapping
        A sidecar file, or an already-parsed flat mapping.
    amplicon_left_window, amplicon_right_window : str or tuple, optional
        The genome-track bc1-amplicon contaminant telltale window, 1-based inclusive,
        as `"10224-10235"` or `(10224, 10235)`. Omit to disable the filter. This is a
        per-screen parameter and belongs in the screen config, not in a default here.

    Returns
    -------
    tuple
        `(genome_track, plasmid_track)` exactly as TRACE builds them — genome
        `revcomp=True`, plasmid `revcomp=False`. Orientation is TRACE's call, not ours.
    """
    cfg = dict(sidecar) if isinstance(sidecar, Mapping) else load_locus_sidecar(sidecar)
    left = amplicon_left_window if isinstance(amplicon_left_window, tuple) else parse_window(
        amplicon_left_window
    )
    right = amplicon_right_window if isinstance(amplicon_right_window, tuple) else parse_window(
        amplicon_right_window
    )
    return _trace.bc1_consensus().tracks_from_config(
        cfg, amplicon_left_window=left, amplicon_right_window=right
    )


def validate_tracks(
    tracks: Sequence[TrackConfig],
    fasta: str | Path | pysam.FastaFile,
) -> None:
    """Assert each track's locus is an N-run flanked by its configured junctions.

    Delegates to `TrackConfig.validate_against_reference`, which fails loudly on
    reference drift — the bug class that silently broke hard-coded coordinates upstream.

    Parameters
    ----------
    tracks : sequence of TrackConfig
        The tracks to validate. Pass `[genome]` for a genome-only run; TRACE's own
        `run()` validates the plasmid track unconditionally, which fails on a reference
        that legitimately has no plasmid contig.
    fasta : str, Path or pysam.FastaFile
        An indexed reference FASTA, or an already-open handle. A path is opened and
        closed here; an open handle is left to its owner.
    """
    if isinstance(fasta, (str, Path)):
        import pysam as _pysam  # local: keeps pysam off the import path of this module

        with _pysam.FastaFile(str(fasta)) as handle:
            _validate_all(tracks, handle)
    else:
        _validate_all(tracks, fasta)


def _validate_all(tracks: Sequence[TrackConfig], handle: Any) -> None:
    for track in tracks:
        track.validate_against_reference(handle)
