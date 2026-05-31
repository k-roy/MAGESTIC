"""
Early merge module for combining 2x150 and 2x300 bc1-donor-bc0 data.

This module implements the early merge strategy that:
1. Combines data from both sequencing depths at the bc1-donor-bc0 level
2. Deduplicates records to avoid redundant processing
3. Tracks data source(s) for each record
4. Handles conflicts when donors differ between sources
5. Counts PCR replicates across all sources for increased confidence

Key insight: The 30bp donor portion + 10bp bc0 = 40bp unique identifier
that can link 2x150 fragments to 2x300 complete sequences.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict
import pandas as pd
import logging

from ..config import COMMON_SEQUENCE_START, BC0_LENGTH

logger = logging.getLogger(__name__)


@dataclass
class BC1Record:
    """
    Unified record for a bc1-donor-bc0 combination.

    This dataclass represents the merged view of data from multiple sources.
    A single bc1 may have data from 2x150, 2x300, and/or gdb sources.
    """
    bc1: str
    bc0: str

    # Donor information - may come from different sources
    donor_full: Optional[str] = None           # Full donor from 2x300
    donor_bc0_fragment: Optional[str] = None   # 60bp fragment from 2x150
    donor_from_gdb: Optional[str] = None       # Donor from guide-donor-bc0

    # Guide information (from gdb)
    guide: Optional[str] = None

    # Data source tracking
    has_2x150: bool = False
    has_2x300: bool = False
    has_gdb: bool = False

    # Read counts from each source
    counts_2x150: int = 0
    counts_2x300: int = 0
    counts_gdb: int = 0

    # PCR replicate counts from each source
    pcr_replicates_2x150: int = 0
    pcr_replicates_2x300: int = 0

    # Conflict tracking
    has_donor_conflict: bool = False
    conflict_details: Optional[str] = None

    # Additional metadata
    library_id: Optional[str] = None
    scaffold: Optional[str] = None

    @property
    def total_counts(self) -> int:
        """Total read counts across all sources."""
        return self.counts_2x150 + self.counts_2x300

    @property
    def total_pcr_replicates(self) -> int:
        """Total PCR replicates across all sources."""
        return self.pcr_replicates_2x150 + self.pcr_replicates_2x300

    @property
    def best_donor(self) -> Optional[str]:
        """
        Get the best available donor sequence.

        Priority:
        1. Full donor from 2x300 (most complete)
        2. Donor from gdb (if available and consistent)
        3. None if only fragment available (need to infer from matching)
        """
        if self.donor_full:
            return self.donor_full
        if self.donor_from_gdb:
            return self.donor_from_gdb
        return None

    @property
    def donor_bc0_key(self) -> Optional[str]:
        """
        Get the donor_bc0 key for linking records.

        This is the 30bp donor portion + 10bp bc0 that uniquely identifies
        a donor-bc0 combination and can link 2x150 fragments to 2x300 data.
        """
        if self.donor_bc0_fragment:
            # Extract donor portion (before common sequence)
            if COMMON_SEQUENCE_START in self.donor_bc0_fragment:
                donor_portion = self.donor_bc0_fragment[:self.donor_bc0_fragment.index(COMMON_SEQUENCE_START)]
            else:
                # Fallback: assume last 10bp is bc0, rest is donor portion
                donor_portion = self.donor_bc0_fragment[:-BC0_LENGTH]
            return donor_portion + self.bc0

        if self.donor_full:
            # Extract last 30bp of donor + bc0
            donor_portion = self.donor_full[-30:] if len(self.donor_full) >= 30 else self.donor_full
            return donor_portion + self.bc0

        return None

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame creation."""
        return {
            'bc1': self.bc1,
            'bc0': self.bc0,
            'donor': self.best_donor,
            'donor_full': self.donor_full,
            'donor_bc0_fragment': self.donor_bc0_fragment,
            'donor_from_gdb': self.donor_from_gdb,
            'guide': self.guide,
            'has_2x150': self.has_2x150,
            'has_2x300': self.has_2x300,
            'has_gdb': self.has_gdb,
            'counts_2x150': self.counts_2x150,
            'counts_2x300': self.counts_2x300,
            'counts_gdb': self.counts_gdb,
            'total_counts': self.total_counts,
            'pcr_replicates_2x150': self.pcr_replicates_2x150,
            'pcr_replicates_2x300': self.pcr_replicates_2x300,
            'total_pcr_replicates': self.total_pcr_replicates,
            'has_donor_conflict': self.has_donor_conflict,
            'conflict_details': self.conflict_details,
            'library_id': self.library_id,
            'scaffold': self.scaffold,
        }


def extract_donor_portion_from_fragment(fragment: str) -> str:
    """
    Extract the donor portion from a donor_bc0_fragment.

    The fragment structure is: ~30bp donor + ~20bp common + 10bp bc0
    Common sequence starts with 'CGAT'.

    Args:
        fragment: The 60bp donor_bc0_fragment

    Returns:
        The donor portion (before common sequence)
    """
    if COMMON_SEQUENCE_START in fragment:
        return fragment[:fragment.index(COMMON_SEQUENCE_START)]
    # Fallback: assume structure is 30bp donor + 20bp common + 10bp bc0
    return fragment[:30]


def build_donor_bc0_key(donor: str, bc0: str, use_last_n_bp: int = 30) -> str:
    """
    Build a donor_bc0 key for linking records.

    Args:
        donor: Full or partial donor sequence
        bc0: The bc0 barcode
        use_last_n_bp: How many bp from end of donor to use

    Returns:
        Concatenated key: donor_portion + bc0
    """
    if len(donor) >= use_last_n_bp:
        donor_portion = donor[-use_last_n_bp:]
    else:
        donor_portion = donor
    return donor_portion + bc0


def merge_2x150_data(df_2x150: pd.DataFrame, show_progress: bool = True) -> Dict[str, BC1Record]:
    """
    Process 2x150 data into BC1Record objects.

    PERFORMANCE: Uses itertuples instead of iterrows for ~10x faster iteration.

    Expected columns:
    - bc1: The bc1 barcode
    - donor_bc0_fragment: The 60bp fragment (or 'donor' if that's the column name)
    - bc0: The bc0 barcode
    - total_counts or counts: Read counts
    - num_PCR_replicates: PCR replicate count

    Args:
        df_2x150: DataFrame with 2x150 bc1-donor_bc0_fragment data
        show_progress: Whether to show a progress bar

    Returns:
        Dictionary mapping bc1 -> BC1Record
    """
    from tqdm import tqdm

    records: Dict[str, BC1Record] = {}

    # Handle column name variations
    fragment_col = 'donor_bc0_fragment' if 'donor_bc0_fragment' in df_2x150.columns else 'donor'
    counts_col = 'total_bc1_donor_bc0_counts' if 'total_bc1_donor_bc0_counts' in df_2x150.columns else 'counts'
    pcr_col = 'num_PCR_replicates' if 'num_PCR_replicates' in df_2x150.columns else None
    # 2026-05-31 fix: when the parsed TSV does not carry a pre-aggregated
    # num_PCR_replicates column (the current parser does not), derive the
    # n-distinct-samples count from sample_number (tagged per-row by
    # data_loader.load_2x150_bc1_donor_bc0_data ~line 273).
    sample_col = 'sample_number' if pcr_col is None and 'sample_number' in df_2x150.columns else None
    has_library_id = 'library_ID' in df_2x150.columns

    # Prepare subset of columns for efficient iteration
    cols_needed = ['bc1', 'bc0', fragment_col, counts_col]
    if pcr_col:
        cols_needed.append(pcr_col)
    elif sample_col:
        cols_needed.append(sample_col)
    if has_library_id:
        cols_needed.append('library_ID')

    df_subset = df_2x150[cols_needed]

    # samples-per-bc1 accumulator: drives pcr_replicates_2x150 when the input
    # lacks a pre-aggregated num_PCR_replicates column.
    samples_per_bc1: Dict[str, set] = defaultdict(set) if sample_col else {}

    # Use itertuples for faster iteration
    iterator = df_subset.itertuples(index=False)
    if show_progress:
        iterator = tqdm(iterator, total=len(df_subset), desc="Processing 2x150")

    for row in iterator:
        bc1 = row.bc1
        bc0 = row.bc0
        # HANDOFF_20260528 fix (Agent B, yL437/yL442): removed dead getattr line that threw
        # AttributeError when pandas-mangled `donor_bc0_fragment` → `donorbc0fragment` did
        # not match the actual itertuples attribute name. Position-based access below is the
        # correct path; the getattr was always overridden by line 231.
        fragment = row[2]  # Third column is fragment_col
        counts = row[3] if pd.notna(row[3]) else 0

        if bc1 not in records:
            records[bc1] = BC1Record(bc1=bc1, bc0=bc0)

        record = records[bc1]
        record.has_2x150 = True
        record.donor_bc0_fragment = fragment
        # 2026-05-31 fix: SUM counts across all rows seen for this bc1 (one
        # row per (bc1, sample) in the concatenated parsed TSVs). Prior code
        # overwrote on each iteration, so final value was just the last row's.
        record.counts_2x150 += int(counts)

        if pcr_col:
            pcr_val = row[4]
            # 2026-05-31 fix: SUM pre-aggregated PCR replicate counts across rows.
            record.pcr_replicates_2x150 += int(pcr_val) if pd.notna(pcr_val) else 0
        elif sample_col:
            samples_per_bc1[bc1].add(row[4])

        if has_library_id:
            lib_idx = 5 if (pcr_col or sample_col) else 4
            lib_id = row[lib_idx] if len(row) > lib_idx else None
            if pd.notna(lib_id):
                record.library_id = lib_id

    # 2026-05-31 fix: when num_PCR_replicates was absent, derive
    # pcr_replicates_2x150 from distinct sample_number count per bc1.
    if sample_col:
        for bc1, samples in samples_per_bc1.items():
            records[bc1].pcr_replicates_2x150 = len(samples)

    logger.info(f"Processed {len(records)} bc1 records from 2x150 data")
    return records


def merge_2x300_data(df_2x300: pd.DataFrame,
                     existing_records: Optional[Dict[str, BC1Record]] = None,
                     show_progress: bool = True) -> Dict[str, BC1Record]:
    """
    Process 2x300 data into BC1Record objects, optionally merging with existing records.

    PERFORMANCE: Uses vectorized operations where possible and itertuples for iteration.

    Expected columns:
    - bc1: The bc1 barcode
    - donor: The full donor sequence
    - bc0: The bc0 barcode
    - total_bc1_donor_bc0_counts or counts: Read counts
    - num_PCR_replicates: PCR replicate count

    Args:
        df_2x300: DataFrame with 2x300 bc1-donor-bc0 data
        existing_records: Optional dict of existing records to merge into
        show_progress: Whether to show a progress bar

    Returns:
        Dictionary mapping bc1 -> BC1Record
    """
    from tqdm import tqdm

    records = existing_records if existing_records is not None else {}

    # Handle column name variations
    donor_col = 'donor' if 'donor' in df_2x300.columns else 'bdb_donor'
    counts_col = 'total_bc1_donor_bc0_counts' if 'total_bc1_donor_bc0_counts' in df_2x300.columns else 'counts'
    pcr_col = 'num_PCR_replicates' if 'num_PCR_replicates' in df_2x300.columns else None
    # 2026-05-31 fix: derive pcr_replicates_2x300 from distinct sample_number
    # values per bc1 when the parsed TSV lacks a pre-aggregated
    # num_PCR_replicates column (the current 2x300 parser does not produce one;
    # data_loader.load_2x300_data tags each row with sample_number ~line 153).
    sample_col = 'sample_number' if pcr_col is None and 'sample_number' in df_2x300.columns else None
    has_library_id = 'library_ID' in df_2x300.columns

    new_records = 0
    merged_records = 0

    # Get numpy arrays for faster access
    bc1_vals = df_2x300['bc1'].values
    bc0_vals = df_2x300['bc0'].values
    donor_vals = df_2x300[donor_col].values
    counts_vals = df_2x300[counts_col].values
    pcr_vals = df_2x300[pcr_col].values if pcr_col else None
    sample_vals = df_2x300[sample_col].values if sample_col else None
    lib_vals = df_2x300['library_ID'].values if has_library_id else None

    # samples-per-bc1 accumulator for the no-num_PCR_replicates path.
    samples_per_bc1: Dict[str, set] = defaultdict(set) if sample_col else {}

    iterator = range(len(df_2x300))
    if show_progress:
        iterator = tqdm(iterator, desc="Processing 2x300")

    for i in iterator:
        bc1 = bc1_vals[i]
        bc0 = bc0_vals[i]
        donor = donor_vals[i]

        if bc1 not in records:
            records[bc1] = BC1Record(bc1=bc1, bc0=bc0)
            new_records += 1
        else:
            merged_records += 1

        record = records[bc1]
        record.has_2x300 = True
        record.donor_full = donor
        # 2026-05-31 fix: SUM counts across all rows for this bc1 (one row per
        # (bc1, sample) in the concatenated parsed TSVs). Prior code overwrote
        # so final value was just the last row's count.
        record.counts_2x300 += int(counts_vals[i]) if pd.notna(counts_vals[i]) else 0

        if pcr_vals is not None:
            # 2026-05-31 fix: SUM pre-aggregated PCR replicate counts across rows.
            record.pcr_replicates_2x300 += int(pcr_vals[i]) if pd.notna(pcr_vals[i]) else 0
        elif sample_vals is not None:
            samples_per_bc1[bc1].add(sample_vals[i])

        # Check for bc0 consistency
        if record.bc0 != bc0:
            logger.warning(f"bc0 mismatch for bc1 {bc1}: 2x150={record.bc0}, 2x300={bc0}")

        # Extract library_id if available
        if lib_vals is not None and pd.notna(lib_vals[i]):
            record.library_id = lib_vals[i]

    # 2026-05-31 fix: when num_PCR_replicates was absent, derive
    # pcr_replicates_2x300 from distinct sample_number count per bc1.
    if sample_vals is not None:
        for bc1, samples in samples_per_bc1.items():
            records[bc1].pcr_replicates_2x300 = len(samples)

    logger.info(f"Processed 2x300 data: {new_records} new records, {merged_records} merged with existing")
    return records


def merge_gdb_data(df_gdb: pd.DataFrame,
                   existing_records: Dict[str, BC1Record]) -> Dict[str, BC1Record]:
    """
    Merge guide-donor-bc0 data into existing records.

    The merge is done via the donor_bc0 key (30bp donor portion + 10bp bc0).

    Expected columns:
    - guide: The guide sequence
    - donor: The donor sequence
    - bc0: The bc0 barcode
    - donor_bc0_fragment: The 60bp fragment
    - oligo_name: The matched oligo (if available)
    - counts or grouped_guide_bc0_counts: Read counts

    Args:
        df_gdb: DataFrame with guide-donor-bc0 data
        existing_records: Dict of existing bc1 records to merge into

    Returns:
        Updated dictionary of BC1Records
    """
    # Build lookup from donor_bc0_key -> gdb rows
    gdb_lookup: Dict[str, List[pd.Series]] = defaultdict(list)

    for _, row in df_gdb.iterrows():
        donor = row['donor']
        bc0 = row['bc0']
        key = build_donor_bc0_key(donor, bc0)
        gdb_lookup[key].append(row)

    # Also build lookup from fragment -> gdb rows
    fragment_lookup: Dict[str, List[pd.Series]] = defaultdict(list)
    if 'donor_bc0_fragment' in df_gdb.columns:
        for _, row in df_gdb.iterrows():
            fragment = row['donor_bc0_fragment']
            if pd.notna(fragment):
                fragment_lookup[fragment].append(row)

    # Merge gdb data into existing records
    matched_count = 0
    for bc1, record in existing_records.items():
        # Try to match via donor_bc0_key
        key = record.donor_bc0_key
        gdb_matches = []

        if key and key in gdb_lookup:
            gdb_matches = gdb_lookup[key]
        elif record.donor_bc0_fragment and record.donor_bc0_fragment in fragment_lookup:
            gdb_matches = fragment_lookup[record.donor_bc0_fragment]

        if gdb_matches:
            # Take the best match (highest counts)
            counts_col = 'grouped_guide_bc0_counts' if 'grouped_guide_bc0_counts' in gdb_matches[0].index else 'counts'
            best_match = max(gdb_matches, key=lambda x: x.get(counts_col, 0))

            record.has_gdb = True
            record.guide = best_match.get('guide')
            record.donor_from_gdb = best_match.get('donor')
            record.counts_gdb = int(best_match.get(counts_col, 0))

            if 'library_ID' in best_match.index and pd.notna(best_match.get('library_ID')):
                record.library_id = best_match['library_ID']
            if 'scaffold' in best_match.index and pd.notna(best_match.get('scaffold')):
                record.scaffold = best_match['scaffold']

            matched_count += 1

    logger.info(f"Merged gdb data: {matched_count}/{len(existing_records)} records matched")
    return existing_records


def resolve_donor_conflicts(records: Dict[str, BC1Record]) -> Dict[str, BC1Record]:
    """
    Resolve conflicts when donors differ between data sources.

    Strategy:
    1. If 2x300 and gdb donors match: high confidence, no conflict
    2. If 2x300 and gdb donors differ: flag as conflict, keep both
    3. If only 2x150 fragment available: check if it matches 2x300/gdb donor end

    Args:
        records: Dictionary of BC1Records to check

    Returns:
        Updated dictionary with conflict flags set
    """
    conflict_count = 0

    for bc1, record in records.items():
        if not record.has_donor_conflict:
            # Check for conflicts between sources
            donors_to_compare = []

            if record.donor_full:
                donors_to_compare.append(('2x300', record.donor_full))
            if record.donor_from_gdb:
                donors_to_compare.append(('gdb', record.donor_from_gdb))

            if len(donors_to_compare) >= 2:
                # Compare donors
                d1_source, d1 = donors_to_compare[0]
                d2_source, d2 = donors_to_compare[1]

                if d1 != d2:
                    record.has_donor_conflict = True
                    record.conflict_details = f"{d1_source} ({len(d1)}bp) != {d2_source} ({len(d2)}bp)"
                    conflict_count += 1

            # Check fragment consistency if available
            if record.donor_bc0_fragment and record.donor_full:
                donor_portion = extract_donor_portion_from_fragment(record.donor_bc0_fragment)
                if not record.donor_full.endswith(donor_portion):
                    if not record.has_donor_conflict:
                        record.has_donor_conflict = True
                        record.conflict_details = "Fragment doesn't match 2x300 donor end"
                        conflict_count += 1

    logger.info(f"Found {conflict_count} records with donor conflicts")
    return records


def merge_data_sources(
    df_2x150: Optional[pd.DataFrame] = None,
    df_2x300: Optional[pd.DataFrame] = None,
    df_gdb: Optional[pd.DataFrame] = None,
    resolve_conflicts: bool = True
) -> pd.DataFrame:
    """
    Main entry point: merge all available data sources into a unified DataFrame.

    This function:
    1. Processes 2x150 data (if available)
    2. Merges 2x300 data (if available)
    3. Merges gdb data (if available)
    4. Resolves donor conflicts
    5. Returns a unified DataFrame

    Args:
        df_2x150: DataFrame with 2x150 bc1-donor_bc0_fragment data
        df_2x300: DataFrame with 2x300 bc1-donor-bc0 data
        df_gdb: DataFrame with guide-donor-bc0 data
        resolve_conflicts: Whether to check and flag donor conflicts

    Returns:
        Unified DataFrame with all merged data
    """
    records: Dict[str, BC1Record] = {}

    # Process each data source
    if df_2x150 is not None and len(df_2x150) > 0:
        logger.info(f"Processing {len(df_2x150)} rows of 2x150 data")
        records = merge_2x150_data(df_2x150)

    if df_2x300 is not None and len(df_2x300) > 0:
        logger.info(f"Processing {len(df_2x300)} rows of 2x300 data")
        records = merge_2x300_data(df_2x300, existing_records=records)

    if df_gdb is not None and len(df_gdb) > 0:
        logger.info(f"Processing {len(df_gdb)} rows of gdb data")
        records = merge_gdb_data(df_gdb, existing_records=records)

    # Resolve conflicts
    if resolve_conflicts:
        records = resolve_donor_conflicts(records)

    # Convert to DataFrame
    if not records:
        logger.warning("No records after merging - returning empty DataFrame")
        return pd.DataFrame()

    df = pd.DataFrame([r.to_dict() for r in records.values()])

    # Sort by total counts descending
    df = df.sort_values('total_counts', ascending=False).reset_index(drop=True)

    # Log summary statistics
    logger.info(f"Merged data summary:")
    logger.info(f"  Total records: {len(df)}")
    logger.info(f"  With 2x150 data: {df['has_2x150'].sum()}")
    logger.info(f"  With 2x300 data: {df['has_2x300'].sum()}")
    logger.info(f"  With gdb data: {df['has_gdb'].sum()}")
    logger.info(f"  With both 2x150 and 2x300: {((df['has_2x150']) & (df['has_2x300'])).sum()}")
    logger.info(f"  With donor conflicts: {df['has_donor_conflict'].sum()}")

    return df


def filter_by_quality(
    df: pd.DataFrame,
    min_reads: int = 2,
    min_pcr_replicates: int = 1,
    min_purity: Optional[float] = None,
    purity_column: Optional[str] = None
) -> pd.DataFrame:
    """
    Filter merged data by quality thresholds.

    Args:
        df: Merged DataFrame from merge_data_sources
        min_reads: Minimum total read counts
        min_pcr_replicates: Minimum total PCR replicates
        min_purity: Minimum purity value (if purity column exists)
        purity_column: Name of the purity column

    Returns:
        Filtered DataFrame
    """
    initial_count = len(df)

    # Filter by read counts
    df = df[df['total_counts'] >= min_reads]
    logger.info(f"After min_reads={min_reads} filter: {len(df)} records (removed {initial_count - len(df)})")

    # Filter by PCR replicates
    if min_pcr_replicates > 1:
        before = len(df)
        df = df[df['total_pcr_replicates'] >= min_pcr_replicates]
        logger.info(f"After min_pcr_replicates={min_pcr_replicates} filter: {len(df)} records (removed {before - len(df)})")

    # Filter by purity
    if min_purity is not None and purity_column and purity_column in df.columns:
        before = len(df)
        df = df[df[purity_column] >= min_purity]
        logger.info(f"After min_purity={min_purity} filter: {len(df)} records (removed {before - len(df)})")

    return df
