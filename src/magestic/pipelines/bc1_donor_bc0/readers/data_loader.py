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

    # Try to find keyfiles with various naming conventions
    plate_key_paths = [
        keyfile_dir / plate_key_name,
        keyfile_dir / f"*_plate_key.tsv",
    ]

    plate_key = None
    for pattern in plate_key_paths:
        if '*' in str(pattern):
            matches = list(keyfile_dir.glob(pattern.name))
            if matches:
                plate_key = pd.read_csv(matches[0], sep='\t')
                break
        elif pattern.exists():
            plate_key = pd.read_csv(pattern, sep='\t')
            break

    if plate_key is None:
        logger.warning(f"No plate keyfile found in {keyfile_dir}")
        return pd.DataFrame(), DataSourceInfo("2x300", 0, 0, 0, [])

    # Filter to bc1-donor-bc0 amplicons
    if 'description' in plate_key.columns:
        plate_key = plate_key[plate_key['description'] == 'bc1-donor-bc0']

    # Load sample key
    sample_key_path = keyfile_dir / sample_key_name
    if not sample_key_path.exists():
        matches = list(keyfile_dir.glob("*_sample_key.tsv"))
        if matches:
            sample_key_path = matches[0]

    if sample_key_path.exists():
        sample_key = pd.read_csv(sample_key_path, sep='\t')
        combined_key = plate_key.merge(sample_key, how='cross')
    else:
        combined_key = plate_key

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
                       'sample_number', 'sequencing_run', 'data_source']
            available_cols = [c for c in std_cols if c in df.columns]
            dfs.append(df[available_cols])
            loaded_files.append(filepath)

    if not dfs:
        logger.warning("No 2x300 data files loaded")
        return pd.DataFrame(), DataSourceInfo("2x300", 0, 0, 0, [])

    combined = pd.concat(dfs, ignore_index=True)

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
            df['sequencing_run'] = '2x150'
            df['data_source'] = '2x150'

            # Extract bc0 from donor_bc0_fragment if not present
            if 'bc0' not in df.columns and 'donor_bc0_fragment' in df.columns:
                df['bc0'] = df['donor_bc0_fragment'].str[-10:]

            std_cols = ['counts', 'bc1', 'donor_bc0_fragment', 'bc0',
                       'sample_number', 'sequencing_run', 'data_source']
            available_cols = [c for c in std_cols if c in df.columns]
            dfs.append(df[available_cols])
            loaded_files.append(filepath)

    if not dfs:
        logger.warning("No 2x150 bc1-donor-bc0 data files loaded")
        return pd.DataFrame(), DataSourceInfo("2x150", 0, 0, 0, [])

    combined = pd.concat(dfs, ignore_index=True)

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
