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
import pandas as pd
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


# Fields composing the PCR-replicate identity, in order. This is the keyfile's
# own PCR_plate_name (gDNA_plate_name + inner_primers + outer_primers) at well
# (sample_number) resolution. sequencing_date is intentionally NOT here.
PCR_REPLICATE_ID_FIELDS = ('gDNA_plate_name', 'inner_primers', 'outer_primers', 'sample_number')


def build_pcr_replicate_id(row) -> str:
    """PCR-replicate identity = one independent amplification of the gDNA, per well.

    Identity = ``(gDNA_plate_name, inner_primers, outer_primers, sample_number)``
    — the keyfile's own ``PCR_plate_name`` (gDNA_plate + inner + outer) at well
    resolution. Rationale:

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
        ``"{gDNA_plate_name}_{inner_primers}_{outer_primers}_{sample_number}"``.
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

    Returns:
        Tuple of (DataFrame, DataSourceInfo)
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

    if sample_key is not None:
        combined_key = plate_key.merge(sample_key, how='cross')
    else:
        combined_key = plate_key

    # One physical file == one (sequencing_date, outer_primers, sample_number);
    # the filename does not encode gDNA_plate_name, so two plate-key rows that
    # differ only in gDNA but share (date, outer) point at the SAME file and must
    # not be loaded twice. Genuine independent re-preps appear as a DIFFERENT
    # sequencing_date (hence a different file) and are preserved.
    file_id_cols = [c for c in ('sequencing_date', 'outer_primers', 'sample_number')
                    if c in combined_key.columns]
    if file_id_cols:
        before = len(combined_key)
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
        # Try various file naming patterns
        patterns = []
        if 'sequencing_date' in row and 'outer_primers' in row and 'sample_number' in row:
            patterns.append(f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}_processed_bc1_donor_bc0.tsv")
            patterns.append(f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}_bc1_donor_bc0.tsv")
        if 'sample_number' in row:
            patterns.append(f"*_{row['sample_number']}_*.tsv")

        filepath = None
        for pattern in patterns:
            if '*' in pattern:
                matches = list(data_dir.glob(pattern))
                if matches:
                    filepath = matches[0]
                    break
            else:
                candidate = data_dir / pattern
                if candidate.exists():
                    filepath = candidate
                    break

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
            df['pcr_replicate_id'] = build_pcr_replicate_id(row)
            df['sequencing_run'] = '2x300'
            df['data_source'] = '2x300'

            # Ensure donor_bc0_fragment exists
            if 'donor_bc0_fragment' not in df.columns and 'donor' in df.columns and 'bc0' in df.columns:
                sps = df.get('SPS', '')
                if 'SPS' in df.columns:
                    df['donor_bc0_fragment'] = df.apply(
                        lambda r: r['donor'][-30:] + r['SPS'] + r['bc0']
                        if pd.notna(r['donor']) and len(str(r['donor'])) >= 30 else None,
                        axis=1
                    )
                else:
                    df['donor_bc0_fragment'] = df.apply(
                        lambda r: r['donor'][-30:] + r['bc0']
                        if pd.notna(r['donor']) and len(str(r['donor'])) >= 30 else None,
                        axis=1
                    )

            # Select standard columns
            std_cols = ['counts', 'bc1', 'donor', 'bc0', 'donor_bc0_fragment',
                       'sample_number', 'outer_primers', 'inner_primers',
                       'gDNA_plate_name', 'sequencing_date', 'pcr_replicate_id',
                       'sequencing_run', 'data_source']
            available_cols = [c for c in std_cols if c in df.columns]
            dfs.append(df[available_cols])
            loaded_files.append(filepath)

    if not dfs:
        logger.warning("No 2x300 data files loaded")
        return pd.DataFrame(), DataSourceInfo("2x300", 0, 0, 0, [])

    combined = pd.concat(dfs, ignore_index=True)
    check_replicate_metadata(combined, "2x300")

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

        filepath = None
        for pattern in patterns:
            if '*' in pattern:
                matches = list(data_dir.glob(pattern))
                if matches:
                    filepath = matches[0]
                    break
            else:
                candidate = data_dir / pattern
                if candidate.exists():
                    filepath = candidate
                    break

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

        filepath = None
        for pattern in patterns:
            if '*' in pattern:
                matches = list(data_dir.glob(pattern))
                if matches:
                    filepath = matches[0]
                    break
            else:
                candidate = data_dir / pattern
                if candidate.exists():
                    filepath = candidate
                    break

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

    # Aggregate to create lookup table
    logger.info("Aggregating donor-bc0 data into lookup table...")
    agg_df = combined.groupby('donor_bc0_fragment').agg({
        'donor': 'first',
        'bc0': 'first',
        'counts': 'sum',
        'sample_number': lambda x: len(set(x))
    }).reset_index()
    agg_df.rename(columns={'sample_number': 'num_samples_donor_bc0'}, inplace=True)

    logger.info(f"  Created lookup table with {len(agg_df)} unique donor_bc0_fragments")

    return agg_df


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

    return merged


def load_gdb_data(
    gdb_dir: Path,
    debug_mode: bool = False,
    debug_n_rows: int = 100000,
) -> Tuple[pd.DataFrame, DataSourceInfo]:
    """
    Load guide-donor-bc0 (gdb) data.

    The gdb data provides guide sequences which can help disambiguate
    oligo assignments.

    Args:
        gdb_dir: Directory containing gdb TSV files
        debug_mode: If True, only load subset of rows
        debug_n_rows: Number of rows to load in debug mode

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

    # Ensure donor_bc0_fragment exists
    if 'donor_bc0_fragment' not in combined.columns and 'donor' in combined.columns and 'bc0' in combined.columns:
        combined['donor_bc0_fragment'] = combined.apply(
            lambda r: r['donor'][-30:] + r['bc0']
            if pd.notna(r['donor']) and len(str(r['donor'])) >= 30 else None,
            axis=1
        )

    info = DataSourceInfo(
        source_type="gdb",
        num_rows=len(combined),
        num_unique_bc1s=0,  # gdb doesn't have bc1
        num_samples=combined['sample_number'].nunique() if 'sample_number' in combined.columns else 0,
        file_paths=gdb_files
    )

    logger.info(f"  Loaded {info.num_rows} rows from {len(gdb_files)} files")

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
