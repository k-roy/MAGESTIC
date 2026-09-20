"""
magestic.pipelines.wgs.readers.wgs_calls — THE REDI SEAM (spec §6).

One loader that normalizes **either** a TRACE `bc1_consensus` emit table **or** the
legacy wide WGS table into a single schema, so that `redi/wgs_check` CAN consume the
shared module instead of re-implementing it (ARCHITECTURE_POLICY §11: exactly one home
per algorithm).

🔴 **NOT YET WIRED, as of 2026-09-13. This file states a capability, not a fact.**
`redi/wgs_check` still defines its own `normalize_wgs_bc1` and `compare_to_wgs`, and
nothing outside this package imports `load_wgs_bc1_calls`. **§11's "exactly one home"
is therefore NOT met on the REDI seam** — two homes exist today. An earlier revision of
this docstring asserted the consumption in the present tense; that was false and is
corrected here. Wiring the consumer is the remaining work, and until it lands this
module is the intended seam rather than the live one.

The returned frame, one row per sample:

    sample_name  join key, already in the consumer's namespace
    bc1          canonical orientation, NATIVE LENGTH — no truncation, no padding
    has_call     bool — distinguishes "no call" from "disagrees"
    confidence   TRACE's high/medium/low/none
    oligo        from the bridge, empty when bridge-less

🔴 **`has_call` is the column that fixes REDI defect D1.** REDI's `compare_to_wgs`
compares with `==` and no na-handling, so a NaN barcode yields `False` — a sample with
NO WGS call is recorded identically to one whose WGS call CONTRADICTS REDI. Measured on
the legacy artifact: 301 of 576 rows (52.3%) had no WGS bc1 at all and all 301 were
recorded as disagreement; headline agreement read 37.2%/39.4% where the call-bearing
subset was 77.8%/82.5%. A >2x understatement produced entirely by conflating absence
with disagreement. Carrying `has_call` explicitly is what lets the consumer state its
denominator.

🔴 **No truncate-and-pad.** REDI's `normalize_wgs_bc1` slices to 20 and right-pads with
`C`, which FABRICATES BASES: in the legacy 874-row table, 8 barcodes were 19 nt and 4
were 21 nt (12/874 = 1.4% altered), and a 19-mer's invented 20th base is correct only
~1 time in 4. TRACE guarantees length STRUCTURALLY instead — `bc_len = bc_end - bc_start
+ 1`, with `validate_against_reference` asserting the locus is an N-run flanked by the
configured junctions — so a call that comes from this module never needs normalizing and
an off-length one is a REAL defect that must stay visible. `normalize_wgs_bc1` should be
deleted from REDI, not moved here.

**Orientation.** `bc1_consensus` emits `TrackConfig.to_canonical(...)` output: the genome
track is built `revcomp=True` and the plasmid `revcomp=False` (the plasmid cassette is
inverted relative to the genome integrant — verified on the v11 reference, where the
plasmid's left/right flanks are the exact reverse complements of the genome track's). So
an emit-table barcode is ALREADY canonical and this module re-orients nothing.

Author: Kevin R. Roy
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

__all__ = [
    "CALL_COLUMNS",
    "LEGACY_BC1_COLUMN",
    "WgsCallsError",
    "load_wgs_bc1_calls",
]

#: The one schema both input shapes are normalized to.
CALL_COLUMNS: tuple[str, ...] = ("sample_name", "bc1", "has_call", "confidence", "oligo")

#: The legacy wide table's barcode column (REDI's `load_wgs_outcome` renames it to `bc1`).
LEGACY_BC1_COLUMN = "bc1_from_wgs"


class WgsCallsError(RuntimeError):
    """The source table is not a shape this loader recognizes, or lacks its key column."""


def load_wgs_bc1_calls(
    source: str | Path,
    *,
    sample_namespace: Callable[[str], str] | Mapping[str, str] | None = None,
) -> Any:
    """Load bc1 calls from a TRACE emit table or the legacy wide table.

    Parameters
    ----------
    source : str or Path
        Either a `bc1_consensus` emit TSV (detected by its `sample` + `barcode` columns)
        or the legacy wide WGS table (detected by `bc1_from_wgs`).
    sample_namespace : callable or mapping, optional
        Maps the source's sample key into the CONSUMER's namespace. **This is the hard
        part, and this module owns it** (spec §6.3): REDI's key is
        `sequencing_date + "_Nextera outer primer set_" + sample_number` where
        `sample_number` is POSITIONALLY SYNTHESIZED from row order in a rearray file,
        while TRACE's `Bc1Call.sample` is a BAM basename. They are different namespaces
        with no mapping anywhere — the already-lossy 523/576 overlap measured between
        them is live breakage, not a hypothetical. A mapping that drops a sample raises;
        silently dropping it would shrink the denominator invisibly.

    Returns
    -------
    pandas.DataFrame
        Columns exactly `CALL_COLUMNS`, one row per sample, source order preserved.

    Raises
    ------
    WgsCallsError
        Unrecognized shape, a duplicated sample key, or a `sample_namespace` mapping
        that has no entry for a sample present in the source.
    """
    import pandas as pd

    path = Path(source)
    if not path.exists():
        raise WgsCallsError(f"WGS bc1 call source not found: {path}")
    frame = pd.read_csv(path, sep="\t", comment=None, dtype=str, keep_default_na=False)
    frame.columns = [c.lstrip("#") for c in frame.columns]

    if "barcode" in frame.columns and "sample" in frame.columns:
        out = _from_trace_emit(frame)
    elif LEGACY_BC1_COLUMN in frame.columns:
        out = _from_legacy_wide(frame)
    else:
        raise WgsCallsError(
            f"{path}: not a recognized WGS bc1 source. Expected either a TRACE emit "
            f"table (columns 'sample' + 'barcode') or the legacy wide table (column "
            f"{LEGACY_BC1_COLUMN!r}). Found: {list(frame.columns)[:12]}"
        )

    if sample_namespace is not None:
        out["sample_name"] = [_map_sample(s, sample_namespace) for s in out["sample_name"]]

    duplicated = out["sample_name"][out["sample_name"].duplicated()].tolist()
    if duplicated:
        raise WgsCallsError(
            f"{path}: duplicate sample key(s) after namespace mapping: "
            f"{sorted(set(duplicated))[:5]}. A duplicate silently multiplies rows in the "
            "consumer's merge and changes its denominator."
        )
    return out[list(CALL_COLUMNS)]


def _from_trace_emit(frame: Any) -> Any:
    """Normalize a `bc1_consensus` emit table (17 columns)."""
    barcode = frame["barcode"].fillna("").astype(str).str.strip()
    out = frame[["sample"]].copy()
    out.columns = ["sample_name"]
    out["bc1"] = barcode
    out["has_call"] = barcode.str.len() > 0
    out["confidence"] = frame.get("confidence", "").astype(str) if "confidence" in frame else ""
    out["oligo"] = frame.get("oligo", "").astype(str) if "oligo" in frame else ""
    return out


def _from_legacy_wide(frame: Any) -> Any:
    """Normalize the legacy wide WGS table.

    `sample_name` is required. The legacy `wgs_outcome_full_table.v20250701.tsv` carries
    `sample_name_short` instead and NO barcode column at all, which is why REDI's step_g
    raises `KeyError` against the only such table that exists on Oak — so a missing key
    column is named explicitly here rather than surfacing three frames down.
    """
    if "sample_name" not in frame.columns:
        alt = "sample_name_short" if "sample_name_short" in frame.columns else None
        raise WgsCallsError(
            "legacy WGS table has no 'sample_name' column"
            + (
                f"; it carries {alt!r}, which is the same namespace under a different "
                "name — pass it through sample_namespace explicitly rather than "
                "renaming silently"
                if alt
                else ""
            )
        )
    barcode = frame[LEGACY_BC1_COLUMN].fillna("").astype(str).str.strip()
    out = frame[["sample_name"]].copy()
    out["bc1"] = barcode
    out["has_call"] = barcode.str.len() > 0
    out["confidence"] = ""
    out["oligo"] = ""
    return out


def _map_sample(sample: str, namespace: Callable[[str], str] | Mapping[str, str]) -> str:
    """Apply the namespace mapping, raising on an unmapped sample."""
    if callable(namespace):
        return namespace(sample)
    if sample not in namespace:
        raise WgsCallsError(
            f"sample {sample!r} has no entry in the sample_namespace mapping. Dropping "
            "it would shrink the consumer's denominator invisibly (spec §6.3)."
        )
    return namespace[sample]
