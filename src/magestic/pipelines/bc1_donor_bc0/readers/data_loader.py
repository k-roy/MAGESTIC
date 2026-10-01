"""
Data loader module for BC1-Donor-BC0 pipeline.

Handles loading processed bc1-donor-bc0 data from various sources:
- 2x150 bp sequencing (partial donor in bc1-donor-bc0, full donor via separate donor-bc0)
- 2x300 bp sequencing (full donor in merged reads)
- Guide-donor-bc0 (gdb) data

All loaders return DataFrames with standardized columns for downstream processing.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re
import pandas as pd
import logging
from tqdm import tqdm

from magestic.pipelines.bc1_donor_bc0.config import (
    FragmentKeyNotBuildable,
    build_donor_bc0_fragment,
)

logger = logging.getLogger(__name__)


class BridgeRecoveredNothing(RuntimeError):
    """A bridge join between two NON-EMPTY tables recovered zero rows.

    🔴 NN 250 / C17. This has to be an exception, not a log line. The gdb bridge used to log
    ``Recovered full donor for 0/N (0.0%)`` at INFO and return normally when its fragment key
    disagreed with stage 01's, so a total pipeline failure was indistinguishable from a measured
    biological zero -- the defect shape that produced wrong conclusions repeatedly in the
    2026-10-01 audits. An EMPTY lookup is a different condition and stays quiet: NN 240's
    non-vacuity control deliberately runs with no bridge.
    """


def ensure_donor_bc0_fragment(df, *, source: str, require: bool = True):
    """Add ``donor_bc0_fragment`` to ``df`` if absent, using the ONE canonical builder.

    Replaces two hand-written fallbacks that had drifted to a 40-nt SPS-less key (NN 250 / C17).
    Every row that cannot produce a joinable key is counted and logged -- never silently
    dropped or silently shortened.

    Args:
        df: a bc1-donor-bc0 or guide-donor-bc0 frame.
        source: 'bdb' / 'gdb' / a filename, for the log line.
        require: if True (default), raise when no row can produce a key at all. A frame that
            cannot be keyed bridges nothing, which must not look like an absence of overlap.
    """
    if "donor_bc0_fragment" in df.columns or df.empty:
        return df
    if "donor" not in df.columns or "bc0" not in df.columns:
        return df
    if "SPS" not in df.columns:
        msg = (
            f"{source}: cannot rebuild donor_bc0_fragment -- the table has donor and bc0 but no "
            f"SPS column, and stage 01's key is donor[-30:] + SPS + bc0 (60 nt). The SPS-less "
            f"key is 40 nt and joins NOTHING while reading as a clean zero (NN 250 / C17). "
            f"Re-run the step that produces donor_bc0_fragment, or supply the SPS."
        )
        if require:
            raise FragmentKeyNotBuildable(msg)
        logger.warning(msg)
        return df

    df = df.copy()
    n_bad = 0

    def _build(r):
        nonlocal n_bad
        try:
            return build_donor_bc0_fragment(r["donor"], r["SPS"], r["bc0"])
        except FragmentKeyNotBuildable:
            n_bad += 1
            return None

    df["donor_bc0_fragment"] = df.apply(_build, axis=1)
    n_ok = int(df["donor_bc0_fragment"].notna().sum())
    logger.info(
        f"{source}: rebuilt donor_bc0_fragment for {n_ok:,}/{len(df):,} rows "
        f"({n_bad:,} lacked a buildable donor/SPS/bc0 and are NaN, not dropped)"
    )
    if n_ok == 0 and require:
        raise FragmentKeyNotBuildable(
            f"{source}: not one of {len(df):,} rows could produce a donor_bc0_fragment; this "
            f"bridges nothing. Refusing to continue and report it as a zero (NN 250 / C17)."
        )
    return df


# Fields composing the PCR-replicate identity, in order.
#
# EXPLICIT-OVER-IMPLICIT (2026-07-28, Kevin): the identity previously relied on
# two *implicit* redundancies that happened to hold for yL406 but are not
# guaranteed by anything:
#
#   * `plate_polymerase` was implied by `outer_primers` (plate_3F=Q5,
#     plate_1A=KAPA). True today; diverges the moment one plate mixes chemistries.
#   * `condition` and `amplicon_level` were implied by `sample_number` (each
#     sample_number happens to sit in exactly one condition and one level).
#     True today; diverges the moment a sample_number is reused across
#     conditions -- and the two would then silently MERGE into one replicate.
#
# Naming them costs nothing (each is a function of a field already present, so
# the distinct-id count is unchanged -- asserted by the loader) and makes the
# identity self-describing plus robust to those keyfile changes.
#
# `sequencing_date` is intentionally NOT here: re-sequencing the same PCR
# product is a technical sequencing replicate, not an independent amplification.
PCR_REPLICATE_ID_FIELDS = (
    'gDNA_plate_name',      # physical gDNA prep / amplification event
    'inner_primers',        # inner PCR primer pair
    'outer_primers',        # outer PCR primer plate (plate_3F / plate_1A)
    'plate_polymerase',     # EXPLICIT: from plate_key, not inferred from the plate
    'condition',            # EXPLICIT: plate wash vs 5FC outgrowth
    'amplicon_level',       # EXPLICIT: trafo-level vs pooled yL-level amplicon
    'sample_number',        # the well
)

# Values for the derived `amplicon_level` column.
AMPLICON_LEVEL_TRAFO = 'trafo'       # e.g. trafo_ID == "trafo_2-1"
AMPLICON_LEVEL_POOLED = 'pooled_yL'  # e.g. trafo_ID == "yL387" (trafos pooled BEFORE barcoding)
AMPLICON_LEVEL_UNKNOWN = 'unknown'


def derive_amplicon_level(trafo_id) -> str:
    """Classify a bc1-donor-bc0 amplicon as trafo-level or pooled yL-level.

    Both are legitimate, INDEPENDENT PCR replicates for bc1->donor **matching**:
    the bc1<->donor linkage is fixed at barcoding and does not change during
    outgrowth, so a per-trafo amplicon and a pooled yL amplicon of the same
    material are independent amplifications of the same linkage. (They are NOT
    interchangeable for abundance/LFC work -- different biological material.)

    yL406 example: 45 distinct trafo-level sample_numbers (``trafo_1-1`` ...)
    and 48 pooled yL-level ones (``yL386``..``yL405``); the trafo_key maps
    yL386->trafo_1, yL387->trafo_2, ... so the same bc1 appears in both. That is
    what lets a bc1 reach 4 replicates (2 levels x {Q5, KAPA}) rather than 2.

    Args:
        trafo_id: The ``trafo_ID`` cell (``"trafo_2-1"``, ``"yL387"``, or empty).

    Returns:
        One of ``trafo`` / ``pooled_yL`` / ``unknown``.
    """
    s = str(trafo_id).strip() if pd.notna(trafo_id) else ''
    if re.fullmatch(r'yL\d+', s):
        return AMPLICON_LEVEL_POOLED
    if s.lower().startswith('trafo'):
        return AMPLICON_LEVEL_TRAFO
    return AMPLICON_LEVEL_UNKNOWN


def build_pcr_replicate_id(row) -> str:
    """PCR-replicate identity = one independent amplification of the gDNA, per well.

    Identity = :data:`PCR_REPLICATE_ID_FIELDS` — the keyfile's own
    ``PCR_plate_name`` (gDNA_plate + inner + outer) at well resolution, plus the
    three attributes made **explicit** on 2026-07-28 (``plate_polymerase``,
    ``condition``, ``amplicon_level``). Each of those three is a function of a
    field already in the identity for the yL406 layout, so naming them does not
    change the distinct-id count (the loader asserts this); they exist so the id
    is self-describing and stays correct if a future keyfile breaks the
    coincidence. Rationale for the original four:

    - ``gDNA_plate_name`` identifies the physical prep/amplification event, so an
      independent re-prep or re-amplification that *coincidentally* reuses the
      same outer primers is still distinguished (different gDNA plate → a new
      replicate) instead of being silently collapsed.
    - ``inner_primers`` / ``outer_primers`` distinguish independent amplifications
      of one prep (the same sample split across primer plates = real technical
      replicates).
    - ``sample_number`` is the well.
    - ``sequencing_date`` is deliberately EXCLUDED: re-sequencing the *same*
      library (same gDNA_plate + primers + well, new flowcell/date) is a
      sequencing replicate, not an independent PCR replicate — a PCR/chimera
      artifact reappears in every re-sequencing, so it is not independent
      evidence and must collapse to a single replicate.

    Absent/NaN fields render as ``"NA"`` so the id degrades gracefully (never
    worse than the prior sample_number-only behaviour). When ``gDNA_plate_name``
    is unavailable the re-prep guard is inactive; :func:`check_replicate_metadata`
    warns so this is never silent.

    Args:
        row: A keyfile row (pandas Series or mapping).

    Returns:
        The fields of :data:`PCR_REPLICATE_ID_FIELDS` joined by ``_``.
    """
    def _f(key: str) -> str:
        val = row.get(key) if hasattr(row, 'get') else (row[key] if key in row else None)
        return str(val) if pd.notna(val) else 'NA'

    return "_".join(_f(k) for k in PCR_REPLICATE_ID_FIELDS)


def check_replicate_metadata(df: pd.DataFrame, source_label: str = "data") -> None:
    """Sanity-check PCR-replicate metadata after loading; warn (never silently).

    Two conditions are surfaced:

    1. ``gDNA_plate_name`` absent / all-NA → the re-prep guard is INACTIVE; the
       replicate id degrades to ``(inner_primers, outer_primers, sample_number)``,
       so an independent re-prep that reused the same outer primers would be
       under-counted as one replicate. Reported so the user can supply
       ``gDNA_plate_name`` (it lives in the plate_key).
    2. The same ``(sample_number, outer_primers)`` maps to >1 ``gDNA_plate_name``
       → genuine independent preps/amplifications reused the same outer primers.
       These are now (correctly) counted as separate PCR replicates; reported so
       the count change is visible rather than silent.

    Args:
        df: A loaded source frame (post-concat, pre-merge).
        source_label: Tag for the log line (e.g. "2x300", "2x150").
    """
    if df is None or len(df) == 0:
        return

    cols = df.columns
    have_gdna = (
        'gDNA_plate_name' in cols
        and df['gDNA_plate_name'].notna().any()
        and (df['gDNA_plate_name'].astype(str) != 'NA').any()
    )
    if not have_gdna:
        logger.warning(
            f"[{source_label}] gDNA_plate_name unavailable -> PCR-replicate re-prep guard "
            f"INACTIVE; replicate id degrades to (inner_primers, outer_primers, sample_number). "
            f"An independent re-prep that reused the same outer primers would be under-counted "
            f"as one replicate. Supply gDNA_plate_name (plate_key) to activate the guard."
        )
        return

    if 'sample_number' in cols and 'outer_primers' in cols:
        combos = df[['sample_number', 'outer_primers', 'gDNA_plate_name']].drop_duplicates()
        per = combos.groupby(['sample_number', 'outer_primers'])['gDNA_plate_name'].nunique()
        multi = per[per > 1]
        if len(multi) > 0:
            logger.warning(
                f"[{source_label}] {len(multi)} (sample_number, outer_primers) combo(s) span "
                f">1 gDNA_plate_name -> independent re-preps reused the same outer primers; each "
                f"gDNA plate is now counted as a separate PCR replicate. Examples: "
                f"{list(multi.index[:5])}"
            )
        else:
            logger.info(
                f"[{source_label}] replicate-metadata OK: no outer-primer reuse across gDNA plates."
            )


def _read_keyfiles_concat(
    keyfile_dir: Path,
    explicit_name: str,
    glob_pattern: str,
    dedup_subset: Optional[List[str]] = None,
) -> Optional[pd.DataFrame]:
    """Load and concatenate ALL keyfiles matching ``explicit_name`` or ``glob_pattern``.

    Recent screens ship **per-sequencing-date** plate_keys / sample_keys (e.g.
    ``20251125_*_plate_key.tsv`` AND ``20251212_*_plate_key.tsv``). The prior
    code took only the first glob match (``matches[0]``), so the other date's
    rows — and their ``gDNA_plate_name`` / ``inner_primers`` — were dropped,
    silently degrading the PCR-replicate re-prep guard (which needs every date's
    plate_key to tell re-sequencing from an independent re-prep). This loads
    every match and de-duplicates, so single-keyfile screens are unchanged while
    multi-keyfile screens get the full set.

    Args:
        keyfile_dir: Directory to search.
        explicit_name: Exact filename to prefer if present (e.g. "plate_key.tsv").
        glob_pattern: Glob for the per-date variants (e.g. "*_plate_key.tsv").
        dedup_subset: Columns to de-duplicate on after concat (None = full row).

    Returns:
        Concatenated/de-duplicated DataFrame, or None if nothing matched.
    """
    paths: List[Path] = []
    explicit = keyfile_dir / explicit_name
    if explicit.exists():
        paths.append(explicit)
    for m in sorted(keyfile_dir.glob(glob_pattern)):
        if m not in paths:
            paths.append(m)
    if not paths:
        return None

    frames = [pd.read_csv(p, sep='\t') for p in paths]
    out = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0].copy()
    if len(frames) > 1:
        subset = [c for c in dedup_subset if c in out.columns] if dedup_subset else None
        out = out.drop_duplicates(subset=subset).reset_index(drop=True)
        logger.info(
            f"Concatenated {len(paths)} keyfiles matching '{glob_pattern}' "
            f"({[p.name for p in paths]}) -> {len(out)} rows after de-dup"
        )
    return out


def _resolve_unique_file(
    data_dir: Path,
    patterns: List[str],
    context: str,
    allow_ambiguous: bool = False,
) -> Optional[Path]:
    """Resolve ONE data file from ordered candidate patterns, refusing to guess.

    Shared by every per-sample loader. Extracted 2026-07-28 because the same
    ``matches[0]`` construct existed in three places and one of them silently
    discarded two entire re-sequencing runs for months
    (``planning/QTL_analyses/BC1_LOADER_DEFECT_20260727.md``).

    Patterns are tried in order, most specific first. Exact (glob-free) patterns
    are resolved by existence. A glob that matches more than one file is an
    ERROR, not a coin flip: in this layout a loose glob matches one file per
    sequencing run, so picking the first silently selects the oldest run.

    Args:
        data_dir: Directory to resolve within.
        patterns: Ordered candidate filename patterns (globs allowed).
        context: Human-readable identifier for error messages (e.g. the row key).
        allow_ambiguous: If True, warn and take the first match instead of raising.

    Returns:
        The resolved Path, or None if nothing matched.

    Raises:
        ValueError: if a glob matched >1 file and allow_ambiguous is False.
    """
    for pattern in patterns:
        if '*' in pattern:
            matches = sorted(data_dir.glob(pattern))
            if not matches:
                continue
            if len(matches) > 1:
                names = [m.name for m in matches]
                msg = (f"Ambiguous file match for pattern '{pattern}' ({context}): "
                       f"{len(matches)} candidates {names}. This usually means the exact "
                       f"pattern failed and a loose fallback matched every sequencing run. "
                       f"Refusing to guess -- silently taking one is how an entire "
                       f"re-sequencing run gets dropped.")
                if allow_ambiguous:
                    logger.warning(msg + " [allow_ambiguous: taking first]")
                    return matches[0]
                raise ValueError(msg)
            return matches[0]
        candidate = data_dir / pattern
        if candidate.exists():
            return candidate
    return None


def _assert_full_coverage(
    data_dir: Path,
    loaded_files: List[Path],
    label: str,
    glob_pattern: str = "*.tsv",
    require_full_coverage: bool = True,
    debug_mode: bool = False,
) -> None:
    """Fail loudly when files present in ``data_dir`` were never loaded.

    Any one of the three 2x300 loader defects would have been caught instantly by
    this check. Runs even when zero files loaded -- "loaded 0 of 288" is the worst
    case, not an exempt one.
    """
    available = sorted(data_dir.glob(glob_pattern))
    n_avail, n_loaded = len(available), len(set(loaded_files))
    if n_loaded != len(loaded_files):
        logger.warning(f"  {label}: {len(loaded_files) - n_loaded} file(s) loaded more than once")
    if debug_mode:
        logger.info(f"  {label} coverage: {n_loaded}/{n_avail} files (debug_mode -- assertion skipped)")
        return
    if n_loaded == n_avail:
        logger.info(f"  {label} coverage OK: {n_loaded}/{n_avail} files loaded")
        return
    missed = sorted({p.name for p in available} - {p.name for p in loaded_files})
    msg = (f"{label} COVERAGE FAILURE: loaded {n_loaded} of {n_avail} files matching "
           f"'{glob_pattern}' in {data_dir}. Missing {len(missed)}: "
           f"{missed[:12]}{' ...' if len(missed) > 12 else ''}. Every downstream number "
           f"would be computed on partial data. Fix the keyfile/filename mapping, or pass "
           f"require_full_coverage=False with an explicit reason.")
    if require_full_coverage:
        raise ValueError(msg)
    logger.warning(msg + " [require_full_coverage=False]")


@dataclass
class DataSourceInfo:
    """Information about a loaded data source."""
    source_type: str  # "2x150", "2x300", "gdb"
    num_rows: int
    num_unique_bc1s: int
    num_samples: int
    file_paths: List[Path]


def load_keyfile(keyfile_path: Path, filters: Optional[Dict] = None) -> pd.DataFrame:
    """
    Load and optionally filter a keyfile.

    Args:
        keyfile_path: Path to the keyfile TSV
        filters: Dict of column -> value(s) to filter on

    Returns:
        Filtered keyfile DataFrame
    """
    df = pd.read_csv(keyfile_path, sep='\t')

    if filters:
        for col, val in filters.items():
            if col in df.columns:
                if isinstance(val, (list, tuple)):
                    df = df[df[col].isin(val)]
                else:
                    df = df[df[col] == val]

    return df


def load_2x300_data(
    data_dir: Path,
    keyfile_dir: Path,
    plate_key_name: str = "plate_key.tsv",
    sample_key_name: str = "sample_key.tsv",
    debug_mode: bool = False,
    debug_n_samples: int = 4,
    allow_ambiguous: bool = False,
    require_full_coverage: bool = True,
) -> Tuple[pd.DataFrame, DataSourceInfo]:
    """
    Load 2x300 bp bc1-donor-bc0 data.

    The 2x300 bp reads merge to capture the full donor sequence directly.

    Args:
        data_dir: Directory containing processed TSV files
        keyfile_dir: Directory containing keyfiles
        plate_key_name: Name of plate keyfile
        sample_key_name: Name of sample keyfile
        debug_mode: If True, only load a few samples
        debug_n_samples: Number of samples to load in debug mode
        allow_ambiguous: If True, a glob matching >1 file logs a warning and
            takes the first instead of raising. Default False -- silently
            taking one match is exactly how an entire re-sequencing run was
            dropped for months (BC1_LOADER_DEFECT_20260727.md).
        require_full_coverage: If True (default), raise when the number of
            files loaded does not equal the number of candidate .tsv files in
            ``data_dir``. Suppressed automatically in debug_mode.

    Returns:
        Tuple of (DataFrame, DataSourceInfo)

    Raises:
        ValueError: on ambiguous filename matches (unless allow_ambiguous) or
            on incomplete coverage (unless require_full_coverage is False).
    """
    logger.info("Loading 2x300 bp data...")

    # Load ALL matching plate keyfiles (recent screens ship one per sequencing
    # date; each carries that date's gDNA_plate_name needed by the replicate
    # re-prep guard). Concatenate rather than taking only the first match.
    plate_key = _read_keyfiles_concat(keyfile_dir, plate_key_name, "*_plate_key.tsv")

    if plate_key is None:
        logger.warning(f"No plate keyfile found in {keyfile_dir}")
        return pd.DataFrame(), DataSourceInfo("2x300", 0, 0, 0, [])

    # Filter to bc1-donor-bc0 amplicons
    if 'description' in plate_key.columns:
        plate_key = plate_key[plate_key['description'] == 'bc1-donor-bc0']

    # --- DEFECT 1 FIX (2026-07-28) -------------------------------------------
    # Real plate_keys spell this column ``sequencing_dates`` (PLURAL). Every
    # downstream guard tested for the singular, so on real data the multi-date
    # explode below AND the exact filename patterns were both dead code: every
    # row fell through to the loose ``*_{n}_*.tsv`` fallback, which matched all
    # three sequencing runs and silently took ``matches[0]`` -- the 2024 run --
    # for all 288 rows. Normalise to the singular so the rest of the function
    # works as written. See planning/QTL_analyses/BC1_LOADER_DEFECT_20260727.md.
    if 'sequencing_date' not in plate_key.columns and 'sequencing_dates' in plate_key.columns:
        plate_key = plate_key.rename(columns={'sequencing_dates': 'sequencing_date'})
        logger.info("  Normalised plate_key column 'sequencing_dates' -> 'sequencing_date'")
    # -------------------------------------------------------------------------

    # Expand rows whose sequencing_date encodes multiple re-sequencing dates as
    # "d1,d2" into one row per date, so every date's file loads and shares the
    # same gDNA_plate_name -> re-sequencing collapses to one PCR replicate (the
    # filename is per-date, so a comma-joined date would otherwise match no file).
    if 'sequencing_date' in plate_key.columns:
        plate_key = plate_key.copy()
        plate_key['sequencing_date'] = (
            plate_key['sequencing_date'].astype(str).str.split(r'\s*,\s*')
        )
        plate_key = plate_key.explode('sequencing_date', ignore_index=True)
        plate_key['sequencing_date'] = plate_key['sequencing_date'].str.strip()

    # Load ALL matching sample keyfiles (de-dup on sample_number so duplicate
    # per-date sample keys don't double-load the same file via the cross join).
    sample_key = _read_keyfiles_concat(
        keyfile_dir, sample_key_name, "*_sample_key.tsv", dedup_subset=['sample_number'])

    # --- EXPLICIT COLUMN DISAMBIGUATION (2026-07-28) --------------------------
    # `polymerase` and `plate_category` exist in BOTH keyfiles. A cross-merge
    # silently suffixes them to `_x`/`_y`, so a later `row['polymerase']` does
    # not resolve at all -- the exact class of implicit breakage that hid the
    # loader defect. Name both sides explicitly BEFORE merging:
    #   plate_key.polymerase  -> plate_polymerase          (AUTHORITATIVE: the
    #                            chemistry of the physical PCR plate; yL406 has
    #                            plate_3F=Q5 n=96, plate_1A=KAPA n=18)
    #   sample_key.polymerase -> sample_listed_polymerase  (which chemistries a
    #                            sample was *listed* under; enumerated, so a
    #                            sample can appear under both)
    plate_key = plate_key.rename(columns={
        'polymerase': 'plate_polymerase',
        'plate_category': 'plate_category_from_plate_key',
    })
    if sample_key is not None:
        sample_key = sample_key.rename(columns={
            'polymerase': 'sample_listed_polymerase',
            'plate_category': 'plate_category_from_sample_key',
        })
        overlap = set(plate_key.columns) & set(sample_key.columns)
        if overlap:
            logger.warning(
                f"  plate_key/sample_key still share column(s) {sorted(overlap)}; "
                f"a cross-merge will suffix them _x/_y and they will not resolve by name."
            )
        combined_key = plate_key.merge(sample_key, how='cross')
    else:
        combined_key = plate_key

    # Derive the amplicon level explicitly rather than leaving it implied by
    # sample_number (trafo-level vs trafos-pooled-before-barcoding).
    if 'trafo_ID' in combined_key.columns:
        combined_key['amplicon_level'] = combined_key['trafo_ID'].apply(derive_amplicon_level)
        lvl = combined_key['amplicon_level'].value_counts().to_dict()
        logger.info(f"  Amplicon levels (from trafo_ID): {lvl}")
    else:
        combined_key['amplicon_level'] = AMPLICON_LEVEL_UNKNOWN
        logger.warning(
            "  No trafo_ID column -> amplicon_level='unknown' for every row. "
            "Trafo-level and pooled-yL amplicons cannot be distinguished, so a "
            "bc1 present in both will still be counted correctly ONLY because "
            "they carry different sample_numbers."
        )

    # One physical file == one (sequencing_date, outer_primers, sample_number);
    # the filename does not encode gDNA_plate_name, so two plate-key rows that
    # differ only in gDNA but share (date, outer) point at the SAME file and must
    # not be loaded twice. Genuine independent re-preps appear as a DIFFERENT
    # sequencing_date (hence a different file) and are preserved.
    file_id_cols = [c for c in ('sequencing_date', 'outer_primers', 'sample_number')
                    if c in combined_key.columns]
    if file_id_cols:
        before = len(combined_key)
        # Rows sharing a file can DIFFER in the descriptive columns -- yL406's
        # sample_key lists each sample_number twice (once per polymerase) and the
        # second row is often blank (`trafo_ID='empty'`). drop_duplicates keeps
        # the FIRST row, so an arbitrary ordering decided whether the surviving
        # row carried a real trafo_ID/condition or 'empty'. Sort the informative
        # rows first so the choice is deterministic and never loses metadata.
        _info_cols = [c for c in ('trafo_ID', 'condition', 'sample_listed_polymerase',
                                  'plate_replicate', 'generations') if c in combined_key.columns]
        if _info_cols:
            _blank = combined_key[_info_cols].isna() | combined_key[_info_cols].astype(str).isin(
                ['', 'empty', 'NA', 'nan', 'None'])
            combined_key = combined_key.assign(_n_blank=_blank.sum(axis=1)) \
                                       .sort_values('_n_blank', kind='mergesort') \
                                       .drop(columns='_n_blank')
        combined_key = combined_key.drop_duplicates(subset=file_id_cols).reset_index(drop=True)
        if len(combined_key) < before:
            logger.info(f"  De-duplicated combined_key on {file_id_cols}: "
                        f"{before} -> {len(combined_key)} file rows")

    if debug_mode:
        combined_key = combined_key.head(debug_n_samples)

    dfs = []
    loaded_files = []

    for _, row in tqdm(combined_key.iterrows(), total=len(combined_key),
                       desc="Loading 2x300 samples", disable=len(combined_key) < 10):
        # Try various file naming patterns, most specific first.
        # --- DEFECT 2 FIX (2026-07-28) ---------------------------------------
        # Real filenames carry a literal ``sample_`` token:
        #   {date}_{outer_primers}_sample_{n}_bc1_donor_bc0.tsv
        # e.g. 20240821_plate_3F_sample_12_bc1_donor_bc0.tsv
        # The previous exact patterns omitted it and so never matched anything,
        # even once the column-name defect above was fixed. Legacy spellings are
        # retained after the correct ones for other screens.
        patterns = []
        if 'sequencing_date' in row and 'outer_primers' in row and 'sample_number' in row:
            d, o, n = row['sequencing_date'], row['outer_primers'], row['sample_number']
            patterns.append(f"{d}_{o}_sample_{n}_bc1_donor_bc0.tsv")
            patterns.append(f"{d}_{o}_sample_{n}_processed_bc1_donor_bc0.tsv")
            patterns.append(f"{d}_{o}_{n}_processed_bc1_donor_bc0.tsv")  # legacy
            patterns.append(f"{d}_{o}_{n}_bc1_donor_bc0.tsv")            # legacy
        # Loose fallback ONLY when no exact pattern could be constructed (i.e.
        # this screen's keyfile lacks sequencing_date/outer_primers). If exact
        # patterns exist and none matched, the file is genuinely MISSING -- do
        # not then glob across sequencing runs, or a missing file is reported as
        # an ambiguous match and the real cause is masked. The end-of-load
        # coverage assertion reports the missing file properly.
        if not patterns and 'sample_number' in row:
            # Anchored on sample_ and .tsv so it cannot match a different
            # sample_number by substring (e.g. _1_ vs _12_).
            patterns.append(f"*_sample_{row['sample_number']}_*.tsv")
            patterns.append(f"*_{row['sample_number']}_*.tsv")  # legacy, unanchored

        filepath = _resolve_unique_file(
            data_dir, patterns,
            context=f"sample_number={row.get('sample_number')}, "
                    f"date={row.get('sequencing_date')}, outer={row.get('outer_primers')}",
            allow_ambiguous=allow_ambiguous,
        )

        if filepath and filepath.exists():
            df = pd.read_csv(filepath, sep='\t')

            # Standardize columns
            df['sample_number'] = row.get('sample_number', 'unknown')
            # Carry the replicate-identity fields so the merge step can count
            # PCR replicates on the full PCR_plate_name identity
            # (gDNA_plate_name, inner_primers, outer_primers, sample_number)
            # instead of sample_number alone (R2 fix); see build_pcr_replicate_id.
            df['outer_primers'] = row.get('outer_primers') if pd.notna(row.get('outer_primers')) else 'NA'
            df['inner_primers'] = row.get('inner_primers') if pd.notna(row.get('inner_primers')) else 'NA'
            df['gDNA_plate_name'] = row.get('gDNA_plate_name') if pd.notna(row.get('gDNA_plate_name')) else 'NA'
            df['sequencing_date'] = row.get('sequencing_date') if pd.notna(row.get('sequencing_date')) else 'NA'
            # Explicit replicate-identity attributes (2026-07-28). Each is a
            # function of a field already carried, so naming them does not change
            # the distinct-id count -- asserted after the load.
            for col in ('plate_polymerase', 'condition', 'amplicon_level',
                        'trafo_ID', 'sample_listed_polymerase'):
                df[col] = row.get(col) if pd.notna(row.get(col)) else 'NA'
            df['pcr_replicate_id'] = build_pcr_replicate_id(row)
            df['sequencing_run'] = '2x300'
            df['data_source'] = '2x300'

            # Ensure donor_bc0_fragment exists
            # ONE canonical builder (NN 250 / C17). The old else-branch here built an SPS-less
            # 40-nt key that cannot join stage 01's 60-nt key.
            df = ensure_donor_bc0_fragment(df, source=str(filepath.name))

            # Select standard columns
            std_cols = ['counts', 'bc1', 'donor', 'bc0', 'donor_bc0_fragment',
                       'sample_number', 'outer_primers', 'inner_primers',
                       'gDNA_plate_name', 'sequencing_date', 'pcr_replicate_id',
                       # explicit replicate-identity attributes (2026-07-28)
                       'plate_polymerase', 'condition', 'amplicon_level',
                       'trafo_ID', 'sample_listed_polymerase',
                       'sequencing_run', 'data_source']
            available_cols = [c for c in std_cols if c in df.columns]
            dfs.append(df[available_cols])
            loaded_files.append(filepath)

    # DEFECT 4 FIX (2026-07-28): coverage assertion. Runs BEFORE the empty-result
    # early return on purpose -- 'loaded 0 of 288' is the worst case, not an exempt one.
    _assert_full_coverage(data_dir, loaded_files, "2x300",
                          require_full_coverage=require_full_coverage,
                          debug_mode=debug_mode)

    if not dfs:
        logger.warning("No 2x300 data files loaded")
        return pd.DataFrame(), DataSourceInfo("2x300", 0, 0, 0, [])

    combined = pd.concat(dfs, ignore_index=True)
    check_replicate_metadata(combined, "2x300")

    # Naming the previously-implicit attributes must not CHANGE the identity --
    # each added field is a function of one already present, so the distinct-id
    # count has to match the legacy 4-field identity. If it does not, an added
    # field is varying independently (e.g. a sample_number reused across
    # conditions) and the replicate counts would silently inflate.
    if 'pcr_replicate_id' in combined.columns:
        legacy = ('gDNA_plate_name', 'inner_primers', 'outer_primers', 'sample_number')
        if all(c in combined.columns for c in legacy):
            n_legacy = combined[list(legacy)].drop_duplicates().shape[0]
            n_now = combined['pcr_replicate_id'].nunique()
            if n_now != n_legacy:
                logger.warning(
                    f"  [2x300] EXPLICIT-FIELD INFLATION: pcr_replicate_id yields {n_now} distinct "
                    f"ids vs {n_legacy} under the legacy 4-field identity. One of "
                    f"plate_polymerase/condition/amplicon_level varies independently of "
                    f"(gDNA_plate_name, inner_primers, outer_primers, sample_number). This is "
                    f"NOT automatically wrong -- it may be a genuinely finer replicate axis -- "
                    f"but it CHANGES replicate counts, so confirm it is intended."
                )
            else:
                logger.info(
                    f"  [2x300] explicit replicate identity: {n_now} distinct ids, "
                    f"unchanged vs the legacy 4-field identity (no inflation)"
                )

    info = DataSourceInfo(
        source_type="2x300",
        num_rows=len(combined),
        num_unique_bc1s=combined['bc1'].nunique() if 'bc1' in combined.columns else 0,
        num_samples=combined['sample_number'].nunique() if 'sample_number' in combined.columns else 0,
        file_paths=loaded_files
    )

    logger.info(f"  Loaded {info.num_rows} rows, {info.num_unique_bc1s} unique bc1s from {info.num_samples} samples")

    return combined, info


def load_2x150_bc1_donor_bc0_data(
    data_dir: Path,
    keyfile_dir: Path,
    combined_key_name: str = "bc1_donor_bc0_combined_key.tsv",
    debug_mode: bool = False,
    debug_n_samples: int = 4,
    allow_ambiguous: bool = False,
) -> Tuple[pd.DataFrame, DataSourceInfo]:
    """
    Load 2x150 bp bc1-donor-bc0 fragment data.

    The 2x150 bp reads capture bc1 + partial donor (donor_bc0_fragment).
    Full donor must be recovered via separate donor-bc0 amplicon.

    Args:
        data_dir: Directory containing processed TSV files
        keyfile_dir: Directory containing keyfiles
        combined_key_name: Name of combined keyfile
        debug_mode: If True, only load a few samples
        debug_n_samples: Number of samples to load in debug mode

    Returns:
        Tuple of (DataFrame, DataSourceInfo)
    """
    logger.info("Loading 2x150 bp bc1-donor-bc0 data...")

    keyfile_path = keyfile_dir / combined_key_name
    if not keyfile_path.exists():
        # Try to find any keyfile
        matches = list(keyfile_dir.glob("*combined_key*.tsv")) + list(keyfile_dir.glob("*bc1_donor_bc0*.tsv"))
        if matches:
            keyfile_path = matches[0]
        else:
            logger.warning(f"No keyfile found in {keyfile_dir}")
            return pd.DataFrame(), DataSourceInfo("2x150", 0, 0, 0, [])

    combined_key = pd.read_csv(keyfile_path, sep='\t')

    # Filter to bc1-donor-bc0 amplicons from 2x150 run
    if 'description' in combined_key.columns:
        combined_key = combined_key[combined_key['description'] == 'bc1-donor-bc0']
    if 'sequencing_date' in combined_key.columns:
        # Try to identify 2x150 runs (typically earlier dates or specific patterns)
        pass  # Keep all if no clear filter

    if debug_mode:
        combined_key = combined_key.head(debug_n_samples)

    dfs = []
    loaded_files = []

    for _, row in tqdm(combined_key.iterrows(), total=len(combined_key),
                       desc="Loading 2x150 bc1-donor-bc0", disable=len(combined_key) < 10):
        patterns = []
        if 'sequencing_date' in row and 'outer_primers' in row and 'sample_number' in row:
            patterns.append(f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}_bc1_donor_bc0_fragment_counts.tsv")
        if 'sample_number' in row:
            patterns.append(f"*_{row['sample_number']}_*fragment*.tsv")

        # Shared resolver: refuses to guess between multiple glob matches.
        # (2026-07-28 -- same defect class as the 2x300 loader; see
        # planning/QTL_analyses/BC1_LOADER_DEFECT_20260727.md.)
        filepath = _resolve_unique_file(
            data_dir, patterns,
            context=f"sample_number={row.get('sample_number')}",
            allow_ambiguous=allow_ambiguous,
        )

        if filepath and filepath.exists():
            df = pd.read_csv(filepath, sep='\t')

            df['sample_number'] = row.get('sample_number', 'unknown')
            # Carry replicate-identity fields for the full PCR_plate_name
            # replicate id (gDNA_plate_name, inner_primers, outer_primers,
            # sample_number); sequencing_date for provenance only (R2 fix).
            df['outer_primers'] = row.get('outer_primers') if pd.notna(row.get('outer_primers')) else 'NA'
            df['inner_primers'] = row.get('inner_primers') if pd.notna(row.get('inner_primers')) else 'NA'
            df['gDNA_plate_name'] = row.get('gDNA_plate_name') if pd.notna(row.get('gDNA_plate_name')) else 'NA'
            df['sequencing_date'] = row.get('sequencing_date') if pd.notna(row.get('sequencing_date')) else 'NA'
            df['pcr_replicate_id'] = build_pcr_replicate_id(row)
            df['sequencing_run'] = '2x150'
            df['data_source'] = '2x150'

            # Extract bc0 from donor_bc0_fragment if not present
            if 'bc0' not in df.columns and 'donor_bc0_fragment' in df.columns:
                df['bc0'] = df['donor_bc0_fragment'].str[-10:]

            std_cols = ['counts', 'bc1', 'donor_bc0_fragment', 'bc0',
                       'sample_number', 'outer_primers', 'inner_primers',
                       'gDNA_plate_name', 'sequencing_date', 'pcr_replicate_id',
                       'sequencing_run', 'data_source']
            available_cols = [c for c in std_cols if c in df.columns]
            dfs.append(df[available_cols])
            loaded_files.append(filepath)

    if not dfs:
        logger.warning("No 2x150 bc1-donor-bc0 data files loaded")
        return pd.DataFrame(), DataSourceInfo("2x150", 0, 0, 0, [])

    combined = pd.concat(dfs, ignore_index=True)
    check_replicate_metadata(combined, "2x150")

    info = DataSourceInfo(
        source_type="2x150",
        num_rows=len(combined),
        num_unique_bc1s=combined['bc1'].nunique() if 'bc1' in combined.columns else 0,
        num_samples=combined['sample_number'].nunique() if 'sample_number' in combined.columns else 0,
        file_paths=loaded_files
    )

    logger.info(f"  Loaded {info.num_rows} rows, {info.num_unique_bc1s} unique bc1s from {info.num_samples} samples")

    return combined, info


def load_donor_bc0_lookup(
    data_dir: Path,
    keyfile_dir: Path,
    combined_key_name: str = "bc1_donor_bc0_combined_key.tsv",
    debug_mode: bool = False,
    debug_n_samples: int = 4,
    allow_ambiguous: bool = False,
) -> pd.DataFrame:
    """
    Load donor-bc0 data to create a lookup table for donor recovery.

    The donor-bc0 amplicon captures the full donor sequence. We use this
    to recover full donors for 2x150 bc1-donor-bc0 data by linking on
    donor_bc0_fragment.

    Args:
        data_dir: Directory containing donor-bc0 TSV files
        keyfile_dir: Directory containing keyfiles
        combined_key_name: Name of combined keyfile
        debug_mode: If True, only load a few samples
        debug_n_samples: Number of samples to load in debug mode

    Returns:
        DataFrame with columns: donor_bc0_fragment, donor, bc0, counts
    """
    logger.info("Loading donor-bc0 data for donor recovery...")

    keyfile_path = keyfile_dir / combined_key_name
    if not keyfile_path.exists():
        matches = list(keyfile_dir.glob("*combined_key*.tsv"))
        if matches:
            keyfile_path = matches[0]
        else:
            logger.warning(f"No keyfile found in {keyfile_dir}")
            return pd.DataFrame()

    combined_key = pd.read_csv(keyfile_path, sep='\t')

    # Filter to donor-bc0 amplicons
    if 'description' in combined_key.columns:
        combined_key = combined_key[combined_key['description'] == 'donor-bc0']

    if debug_mode:
        combined_key = combined_key.head(debug_n_samples)

    dfs = []

    for _, row in tqdm(combined_key.iterrows(), total=len(combined_key),
                       desc="Loading donor-bc0", disable=len(combined_key) < 10):
        patterns = []
        if 'sequencing_date' in row and 'outer_primers' in row and 'sample_number' in row:
            patterns.append(f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}_complete_donor_bc0_counts.tsv")
            patterns.append(f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}_donor_bc0_counts.tsv")
        if 'sample_number' in row:
            patterns.append(f"*_{row['sample_number']}_*donor_bc0*.tsv")

        # Shared resolver: refuses to guess between multiple glob matches.
        # (2026-07-28 -- same defect class as the 2x300 loader; see
        # planning/QTL_analyses/BC1_LOADER_DEFECT_20260727.md.)
        filepath = _resolve_unique_file(
            data_dir, patterns,
            context=f"sample_number={row.get('sample_number')}",
            allow_ambiguous=allow_ambiguous,
        )

        if filepath and filepath.exists():
            df = pd.read_csv(filepath, sep='\t')
            df['sample_number'] = row.get('sample_number', 'unknown')

            cols = ['counts', 'donor', 'donor_bc0_fragment', 'bc0', 'sample_number']
            available_cols = [c for c in cols if c in df.columns]
            dfs.append(df[available_cols])

    if not dfs:
        logger.warning("No donor-bc0 data files loaded")
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)

    # Aggregate to create lookup table (NN 250 / D5: support-based and counted, see below)
    logger.info("Aggregating donor-bc0 data into lookup table...")
    agg_df, _stats = collapse_donor_bc0_lookup(combined)

    logger.info(f"  Created lookup table with {len(agg_df)} unique donor_bc0_fragments")

    return agg_df


#: Columns searched, in order, for per-row read support when collapsing the donor-bc0 lookup.
#: Same real-names-first rule as _GDB_SUPPORT_COLUMNS: the donor-bc0 tables the loader selects
#: carry `counts`, the gdb-shaped ones carry `total_counts` / `grouped_guide_bc0_counts`.
_DONOR_BC0_SUPPORT_COLUMNS = (
    "counts",
    "total_counts",
    "grouped_guide_bc0_counts",
    "total_fwd_counts",
)


def collapse_donor_bc0_lookup(
    combined: pd.DataFrame,
    *,
    key: str = "donor_bc0_fragment",
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Collapse the donor-bc0 table to one row per ``key``, BY SUPPORT, counting the ambiguity.

    🔴 NN 250 / D5. This replaces
    ``groupby(key).agg({'donor':'first','bc0':'first','counts':'sum', ...})``, which had three
    measured failure modes:

    1. **Support-blind.** ``'first'`` takes whichever row came first in the file. Reproduced: a
       fragment whose ``donor_WRONG`` had 1 read and whose ``donor_RIGHT`` had 180 returned
       ``donor_WRONG``. ``merge.py`` already picked its best gdb match by ``max(counts)``, so a
       support-blind collapse here DEMOTED a correct tie-break to an arbitrary one -- the same
       defect class already fixed in :func:`collapse_gdb_bridge`.
    2. **Uncounted.** Nothing in the output said the fragment had been ambiguous, so the rate
       could not be audited. A fragment with two plausible donors is exactly the condition that
       shows up downstream as a bc1 assigned to the wrong edit.
    3. **Chimeric pair.** ``groupby.agg('first')`` skips NaN *independently per column*, so
       ``donor`` could come from one row and ``bc0`` from another. Reproduced exactly: rows
       ``(NaN, BC0_A)`` and ``(donor_B, NaN)`` emitted ``('donor_B','BC0_A')`` -- a pair present
       in no read. Taking the WHOLE winning row makes that impossible.

    Support is summed per ``(key, donor, bc0)`` across samples first, so a donor supported by
    many small samples beats one supported by a single large read pile-up only if it really has
    more reads. Ties resolve on ``(donor, bc0)`` ascending as terminal keys, so the result is
    independent of input order (tested by permutation invariance, not by inspecting the sort).

    Returns ``(collapsed, stats)``. Added columns: ``donor_support`` (reads behind the winner),
    ``n_distinct_donors`` (how many donors the fragment was seen with),
    ``num_samples_donor_bc0`` (unchanged meaning).
    """
    stats: Dict[str, object] = {
        "rows_in": int(len(combined)),
        "rows_out": 0,
        "dropped_null_key": 0,
        "unique_fragments": 0,
        "ambiguous_fragments": 0,
        "support_column": None,
    }
    if combined.empty or key not in combined.columns:
        return combined, stats

    has_key = combined[key].notna() & (combined[key].astype(str).str.len() > 0)
    stats["dropped_null_key"] = int((~has_key).sum())
    if stats["dropped_null_key"]:
        logger.warning(
            f"  donor-bc0 collapse: {stats['dropped_null_key']:,} of {len(combined):,} rows have "
            f"no {key} and cannot be bridged (counted, not silently dropped)"
        )
    work = combined.loc[has_key].copy()
    if work.empty:
        return work.reset_index(drop=True), stats

    stats["unique_fragments"] = int(work[key].nunique())

    support = next((c for c in _DONOR_BC0_SUPPORT_COLUMNS if c in work.columns), None)
    stats["support_column"] = support
    if support is None:
        logger.warning(
            "  donor-bc0 collapse: no support column found among "
            f"{_DONOR_BC0_SUPPORT_COLUMNS}; falling back to a deterministic (donor, bc0) order. "
            "The winner is reproducible but NOT support-based -- record this in the run's notes."
        )
        work["_support"] = 0
    else:
        work["_support"] = pd.to_numeric(work[support], errors="coerce").fillna(0)

    group_cols = [key] + [c for c in ("donor", "bc0") if c in work.columns]
    # Support per distinct observation, summed across samples.
    obs = work.groupby(group_cols, dropna=False, observed=True)["_support"].sum().reset_index()
    obs = obs.rename(columns={"_support": "donor_support"})

    # Rank observations: most-supported first, then the terminal keys so the order is TOTAL.
    terminal = [c for c in ("donor", "bc0") if c in obs.columns]
    obs = obs.sort_values(
        ["donor_support"] + terminal,
        ascending=[False] + [True] * len(terminal),
        kind="mergesort",
    )
    winners = obs.drop_duplicates(subset=[key], keep="first").set_index(key)

    # Per-fragment ambiguity and totals, from the full table rather than the winner's row.
    per_frag = work.groupby(key, dropna=False, observed=True)
    n_donors = (
        per_frag["donor"].nunique(dropna=True)
        if "donor" in work.columns
        else pd.Series(dtype="int64")
    )
    total_support = per_frag["_support"].sum()
    n_samples = (
        per_frag["sample_number"].nunique(dropna=True)
        if "sample_number" in work.columns
        else pd.Series(1, index=total_support.index)
    )

    out = winners.copy()
    out["n_distinct_donors"] = n_donors.reindex(out.index).fillna(0).astype("int64")
    out["num_samples_donor_bc0"] = n_samples.reindex(out.index).fillna(0).astype("int64")
    if support is not None:
        # Keep the ORIGINAL support column name carrying the fragment's TOTAL, as before.
        out[support] = total_support.reindex(out.index).fillna(0)
    out = out.reset_index()

    stats["rows_out"] = int(len(out))
    stats["ambiguous_fragments"] = int((out["n_distinct_donors"] > 1).sum())
    if stats["ambiguous_fragments"]:
        logger.warning(
            f"  donor-bc0 collapse: {stats['ambiguous_fragments']:,} of {stats['rows_out']:,} "
            f"fragments map to >1 distinct donor; the best-supported one wins and the count is "
            f"carried in n_distinct_donors / donor_support (NN 250 / D5)"
        )
    return out, stats


def recover_full_donors(
    bc1_df: pd.DataFrame,
    donor_lookup: pd.DataFrame
) -> pd.DataFrame:
    """
    Recover full donor sequences for 2x150 bc1-donor-bc0 data.

    Links bc1-donor-bc0 records to full donors via donor_bc0_fragment.

    Args:
        bc1_df: DataFrame with bc1 and donor_bc0_fragment columns
        donor_lookup: Lookup table from load_donor_bc0_lookup()

    Returns:
        DataFrame with donor column added
    """
    if donor_lookup.empty or bc1_df.empty:
        return bc1_df

    logger.info("Recovering full donors for 2x150 data...")

    merged = bc1_df.merge(
        donor_lookup[['donor_bc0_fragment', 'donor']],
        on='donor_bc0_fragment',
        how='left'
    )

    recovered = merged['donor'].notna().sum()
    total = len(merged)
    logger.info(f"  Recovered full donor for {recovered}/{total} ({100*recovered/total:.1f}%) entries")

    # 🔴 NN 250 / C17: zero recovery from two NON-EMPTY tables is a FAILURE, never a result.
    # Both tables having rows while the join produces nothing means the two sides' fragment keys
    # do not agree -- which is exactly how the drifted 40-nt key presented: as "0.0%" at INFO.
    # An empty lookup returns early above, so NN 240's no-bridge non-vacuity control is unaffected.
    if recovered == 0 and total > 0:
        raise BridgeRecoveredNothing(
            f"the gdb bridge recovered 0 of {total:,} donors while the lookup held "
            f"{len(donor_lookup):,} fragments. Both sides are non-empty, so their "
            f"donor_bc0_fragment keys disagree -- check that both were built by "
            f"config.build_donor_bc0_fragment (donor[-30:] + SPS + bc0, 60 nt). Refusing to "
            f"report this as a zero recovery rate (NN 250 / C17)."
        )

    return merged


#: Columns searched, in order, for per-row read support when collapsing the gdb bridge.
#:
#: 🔴 The first three are the names the real gdb tables actually use. An earlier version of
#: this tuple guessed at generic names ("count", "n", "reads", …), none of which exist in a
#: gdb table, so no sort ran and the collapse fell back to first-in-order. That was a
#: regression: ``merge.py`` already picked the best gdb match by ``max(counts)`` (using
#: ``grouped_guide_bc0_counts`` when present), so a support-blind collapse DEMOTED a
#: correct tie-break to an arbitrary one. Caught by the NN 239 pipeline audit, 2026-10-01.
#: Keep the real column names first; the generic ones are a harmless fallback.
_GDB_SUPPORT_COLUMNS = (
    "grouped_guide_bc0_counts",
    "total_counts",
    "counts",
    "count",
    "n_reads",
    "num_reads",
    "read_count",
    "reads",
    "n",
)


def collapse_gdb_bridge(
    combined: pd.DataFrame,
    *,
    key: str = "donor_bc0_fragment",
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Collapse a read-level gdb table to one row per ``key`` for use as a bridge.

    The gdb table is read-level, so a library contributes on the order of tens of
    millions of rows while the thing the bridge actually needs is the mapping
    ``donor_bc0_fragment -> (oligo_name, guide)``, which is at most design size. Merging
    the read-level frame is what made the full-scale bridge run exhaust 128 GB
    (HANDOFF_20260608_gdb_bdb_pipeline.md §3).

    Collapsing is not merely a size optimization, because a fragment can map to more than
    one ``oligo_name`` across reads. The most-supported row wins, and **the ambiguity is
    counted and returned rather than silently resolved** — a fragment with two plausible
    oligos is exactly the condition that shows up downstream as a bc1 assigned to the
    wrong edit, so the rate belongs in the run's statistics.

    Returns ``(collapsed, stats)``.
    """
    stats: Dict[str, int] = {
        "rows_in": int(len(combined)),
        "rows_out": 0,
        "dropped_null_key": 0,
        "unique_keys": 0,
        "ambiguous_keys": 0,
    }
    if combined.empty or key not in combined.columns:
        return combined, stats

    has_key = combined[key].notna() & (combined[key].astype(str).str.len() > 0)
    stats["dropped_null_key"] = int((~has_key).sum())
    work = combined.loc[has_key]
    if work.empty:
        return work.reset_index(drop=True), stats

    stats["unique_keys"] = int(work[key].nunique())

    if "oligo_name" in work.columns:
        n_oligos = work.groupby(key, observed=True)["oligo_name"].nunique(dropna=True)
        stats["ambiguous_keys"] = int((n_oligos > 1).sum())

    support = next((c for c in _GDB_SUPPORT_COLUMNS if c in work.columns), None)
    if support is not None:
        # mergesort is stable, so ties collapse deterministically run to run.
        work = work.sort_values(support, ascending=False, kind="mergesort")

    collapsed = work.drop_duplicates(subset=[key], keep="first").reset_index(drop=True)
    stats["rows_out"] = int(len(collapsed))
    return collapsed, stats


def load_gdb_data(
    gdb_dir: Path,
    debug_mode: bool = False,
    debug_n_rows: int = 100000,
    *,
    collapse_bridge: bool = True,
) -> Tuple[pd.DataFrame, DataSourceInfo]:
    """
    Load guide-donor-bc0 (gdb) data.

    The gdb data provides guide sequences which can help disambiguate
    oligo assignments.

    Args:
        gdb_dir: Directory containing gdb TSV files
        debug_mode: If True, only load subset of rows
        debug_n_rows: Number of rows to load in debug mode
        collapse_bridge: Collapse to one row per ``donor_bc0_fragment`` before
            returning. Default True; the read-level frame is what OOMs the full-scale
            bridge. Pass False only to inspect the raw read-level table.

    Returns:
        Tuple of (DataFrame, DataSourceInfo)
    """
    logger.info("Loading guide-donor-bc0 data...")

    # Find gdb files
    gdb_files = list(gdb_dir.glob("*guide_donor_bc0*.tsv")) + list(gdb_dir.glob("*gdb*.tsv"))

    if not gdb_files:
        logger.warning(f"No gdb files found in {gdb_dir}")
        return pd.DataFrame(), DataSourceInfo("gdb", 0, 0, 0, [])

    dfs = []
    for filepath in tqdm(gdb_files, desc="Loading gdb files", disable=len(gdb_files) < 5):
        if debug_mode:
            df = pd.read_csv(filepath, sep='\t', nrows=debug_n_rows)
        else:
            df = pd.read_csv(filepath, sep='\t')

        df['data_source'] = 'gdb'
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)

    # Ensure donor_bc0_fragment exists, via the ONE canonical builder (NN 250 / C17). This site
    # had no SPS branch at all and built a 40-nt key, so a gdb table lacking the column bridged
    # zero bc1s and the run read as a biological null.
    combined = ensure_donor_bc0_fragment(combined, source='gdb')

    n_samples = (
        combined['sample_number'].nunique() if 'sample_number' in combined.columns else 0
    )
    n_read_level = len(combined)

    if collapse_bridge:
        combined, stats = collapse_gdb_bridge(combined)
        logger.info(
            "  gdb bridge collapsed: %s read-level rows -> %s unique donor_bc0_fragment",
            f"{stats['rows_in']:,}", f"{stats['rows_out']:,}",
        )
        if stats["dropped_null_key"]:
            logger.warning(
                "  %s gdb rows had no donor_bc0_fragment and cannot bridge to a bc1",
                f"{stats['dropped_null_key']:,}",
            )
        if stats["ambiguous_keys"]:
            pct = 100.0 * stats["ambiguous_keys"] / max(stats["unique_keys"], 1)
            logger.warning(
                "  %s of %s donor_bc0_fragments (%.2f%%) map to MORE THAN ONE oligo_name; "
                "the most-supported was kept. This is an upper bound on bc1->oligo "
                "mis-assignment from the bridge.",
                f"{stats['ambiguous_keys']:,}", f"{stats['unique_keys']:,}", pct,
            )

    info = DataSourceInfo(
        source_type="gdb",
        num_rows=len(combined),
        num_unique_bc1s=0,  # gdb doesn't have bc1
        num_samples=n_samples,
        file_paths=gdb_files
    )

    logger.info(
        f"  Loaded {n_read_level:,} rows from {len(gdb_files)} files"
        + (f"; returning {info.num_rows:,} bridge rows" if collapse_bridge else "")
    )

    return combined, info


def load_all_data_sources(
    config,
    debug_mode: bool = False,
    debug_n_samples: int = 4,
) -> Dict[str, Tuple[pd.DataFrame, DataSourceInfo]]:
    """
    Load all available data sources based on configuration.

    Args:
        config: PipelineConfig object with data source paths
        debug_mode: If True, only load subset of data
        debug_n_samples: Number of samples to load in debug mode

    Returns:
        Dict mapping source name to (DataFrame, DataSourceInfo) tuple
    """
    results = {}

    # Load 2x300 data if available
    if config.has_2x300_data and config.data_2x300_dir:
        keyfile_dir = config.project_dir / "keyfiles"
        if not keyfile_dir.exists():
            keyfile_dir = config.data_2x300_dir.parent / "keyfiles"

        df, info = load_2x300_data(
            config.data_2x300_dir,
            keyfile_dir,
            debug_mode=debug_mode,
            debug_n_samples=debug_n_samples
        )
        if not df.empty:
            results['2x300'] = (df, info)

    # Load 2x150 data if available
    if config.has_2x150_data and config.data_2x150_dir:
        keyfile_dir = config.project_dir / "keyfiles"
        if not keyfile_dir.exists():
            keyfile_dir = config.data_2x150_dir.parent / "keyfiles"

        # Load bc1-donor-bc0 fragment data
        df_bc1, info_bc1 = load_2x150_bc1_donor_bc0_data(
            config.data_2x150_dir,
            keyfile_dir,
            debug_mode=debug_mode,
            debug_n_samples=debug_n_samples
        )

        # Load donor-bc0 for donor recovery
        donor_bc0_dir = config.data_2x150_dir.parent / "complete_donor_bc0_counts"
        if not donor_bc0_dir.exists():
            donor_bc0_dir = config.data_2x150_dir

        donor_lookup = load_donor_bc0_lookup(
            donor_bc0_dir,
            keyfile_dir,
            debug_mode=debug_mode,
            debug_n_samples=debug_n_samples
        )

        # Recover full donors
        if not df_bc1.empty and not donor_lookup.empty:
            df_bc1 = recover_full_donors(df_bc1, donor_lookup)

        if not df_bc1.empty:
            results['2x150'] = (df_bc1, info_bc1)

    # Load gdb data if available
    if config.has_gdb_data and config.gdb_dir:
        df, info = load_gdb_data(
            config.gdb_dir,
            debug_mode=debug_mode,
            debug_n_rows=100000 if debug_mode else None
        )
        if not df.empty:
            results['gdb'] = (df, info)

    return results
