"""Merge plasmid-assay evidence at the observed, scoped bc1/donor/bc0 link.

Read support is summed; PCR identities are unioned. Competing links and
unresolved fragments remain separate observations. GDB reads are annotation
evidence only and never enter BDB read or PCR totals.
"""
from collections import defaultdict
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
import json
import logging
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
import pandas as pd
from ..config import COMMON_SEQUENCE_START, BC0_LENGTH

logger = logging.getLogger(__name__)
# These namespaces stay separate. In particular, SPS is not transformation ID.
SCOPE_COLUMNS = ('yeast_library', 'yeast_sublibrary', 'trafo_ID', 'library_ID', 'V_library')


def _text(value):
    return str(value) if value is not None and not pd.isna(value) and str(value) else None


def _count(value, label='read count'):
    try:
        parsed = Decimal(str(value))
        if not parsed.is_finite() or parsed < 0 or parsed != parsed.to_integral_value():
            raise ValueError
        return int(parsed)
    except (InvalidOperation, ValueError, TypeError):
        raise ValueError(f'Invalid {label}: {value!r}; expected a finite nonnegative integer') from None


def _ids(value, *, serialized=False):
    if isinstance(value, (set, frozenset, tuple, list)):
        return {str(v) for v in value if _text(v) is not None}
    if _text(value) is None:
        return set()
    if serialized and isinstance(value, str) and value.lstrip().startswith('['):
        decoded = json.loads(value)
        if not isinstance(decoded, list):
            raise ValueError('PCR identities must be a JSON list')
        return _ids(decoded)
    if serialized and isinstance(value, str):
        return {part for part in value.split(';') if part}
    return {str(value)}


def _fragment_suffix(fragment, sps, bc0):
    """Decode only an explicitly supplied SPS/bc0 boundary; never first-CGAT."""
    if not all(isinstance(v, str) and v for v in (fragment, sps, bc0)):
        return None
    tail = sps + bc0
    if not fragment.endswith(tail) or len(fragment) <= len(tail):
        return None
    return fragment[:-len(tail)]


@dataclass
class BC1Record:
    """One observed link, possibly with a uniquely attached short observation."""
    bc1: Optional[str]
    bc0: Optional[str]
    donor_full: Optional[str] = None
    donor_bc0_fragment: Optional[str] = None
    donor_from_gdb: Optional[str] = None
    guide: Optional[str] = None
    has_2x150: bool = False
    has_2x300: bool = False
    has_gdb: bool = False
    counts_2x150: int = 0
    counts_2x300: int = 0
    counts_gdb: int = 0
    pcr_replicates_2x150: int = 0
    pcr_replicates_2x300: int = 0
    has_donor_conflict: bool = False
    conflict_details: Optional[str] = None
    library_id: Optional[str] = None
    scaffold: Optional[str] = None
    donor_recovered: Optional[str] = None
    SPS: Optional[str] = None
    scope: Tuple = field(default_factory=lambda: (None,) * len(SCOPE_COLUMNS))
    pcr_ids_2x150: Set[str] = field(default_factory=set)
    pcr_ids_2x300: Set[str] = field(default_factory=set)
    anonymous_pcr_2x150: int = 0
    anonymous_pcr_2x300: int = 0
    support_unknown: bool = False
    fragments: Set[str] = field(default_factory=set)
    short_link_candidates: int = 0
    gdb_candidates: List[Dict] = field(default_factory=list)

    @property
    def total_counts(self):
        return self.counts_2x150 + self.counts_2x300

    @property
    def total_pcr_replicates(self):
        # Without identities the largest supplied scalar is only a lower bound.
        return max(len(self.pcr_ids_2x150 | self.pcr_ids_2x300),
                   self.anonymous_pcr_2x150, self.anonymous_pcr_2x300,
                   self.pcr_replicates_2x150, self.pcr_replicates_2x300)

    @property
    def best_donor(self):
        return self.donor_full or self.donor_recovered or self.donor_from_gdb

    @property
    def donor_bc0_key(self):
        if not self.bc0:
            return None
        if self.best_donor:
            return build_donor_bc0_key(self.best_donor, self.bc0)
        suffix = _fragment_suffix(self.donor_bc0_fragment, self.SPS, self.bc0)
        return suffix[-30:] + self.bc0 if suffix else None

    def to_dict(self):
        ids = self.pcr_ids_2x150 | self.pcr_ids_2x300
        result = dict(bc1=self.bc1, bc0=self.bc0, donor=self.best_donor,
            donor_full=self.donor_full, donor_recovered=self.donor_recovered,
            donor_bc0_fragment=self.donor_bc0_fragment, donor_from_gdb=self.donor_from_gdb,
            guide=self.guide, has_2x150=self.has_2x150, has_2x300=self.has_2x300, has_gdb=self.has_gdb,
            counts_2x150=self.counts_2x150, counts_2x300=self.counts_2x300, counts_gdb=self.counts_gdb,
            total_counts=self.total_counts, pcr_replicates_2x150=self.pcr_replicates_2x150,
            pcr_replicates_2x300=self.pcr_replicates_2x300, total_pcr_replicates=self.total_pcr_replicates,
            pcr_replicate_ids=';'.join(sorted(ids)),
            pcr_ids_2x150=json.dumps(sorted(self.pcr_ids_2x150)),
            pcr_ids_2x300=json.dumps(sorted(self.pcr_ids_2x300)),
            pcr_support_exact=not self.support_unknown,
            has_donor_conflict=self.has_donor_conflict, conflict_details=self.conflict_details,
            library_id=self.library_id, scaffold=self.scaffold, SPS=self.SPS,
            input_library_scope=json.dumps(dict(zip(SCOPE_COLUMNS, self.scope)), sort_keys=True),
            observed_fragments=json.dumps(sorted(self.fragments)),
            short_link_candidates=self.short_link_candidates,
            gdb_candidates=json.dumps(self.gdb_candidates, sort_keys=True),
            gdb_ambiguous=len(self.gdb_candidates) > 1,
            link_resolved=bool(self.bc1 and self.bc0 and self.best_donor))
        result.update(dict(zip(SCOPE_COLUMNS, self.scope)))
        return result


def extract_donor_portion_from_fragment(fragment: str) -> str:
    """Legacy convenience heuristic; not used to establish a merge assignment."""
    if COMMON_SEQUENCE_START in fragment:
        return fragment[:fragment.index(COMMON_SEQUENCE_START)]
    return fragment[:30]


def build_donor_bc0_key(donor: str, bc0: str, use_last_n_bp: int = 30) -> str:
    """Candidate lookup key only: a suffix match is not a unique donor call."""
    if not isinstance(donor, str) or not isinstance(bc0, str):
        return ''
    return donor[-use_last_n_bp:] + bc0


def _key(record):
    known = record.donor_full or record.donor_recovered
    return (record.scope, record.bc1, record.bc0, known,
            None if known else record.donor_bc0_fragment, record.SPS)


def _refresh_support(record):
    for source in ('2x150', '2x300'):
        setattr(record, 'pcr_replicates_' + source,
                max(len(getattr(record, 'pcr_ids_' + source)),
                    getattr(record, 'anonymous_pcr_' + source)))


def _combine(target, other):
    """Transfer one observation exactly once, with SUM reads / UNION IDs."""
    for source in ('2x150', '2x300'):
        setattr(target, 'counts_' + source, getattr(target, 'counts_' + source) + getattr(other, 'counts_' + source))
        setattr(target, 'has_' + source, getattr(target, 'has_' + source) or getattr(other, 'has_' + source))
        getattr(target, 'pcr_ids_' + source).update(getattr(other, 'pcr_ids_' + source))
        setattr(target, 'anonymous_pcr_' + source,
                max(getattr(target, 'anonymous_pcr_' + source), getattr(other, 'anonymous_pcr_' + source)))
    target.support_unknown |= other.support_unknown
    target.donor_full = target.donor_full or other.donor_full
    target.donor_recovered = target.donor_recovered or other.donor_recovered
    target.SPS = target.SPS or other.SPS
    target.fragments.update(other.fragments)
    target.donor_bc0_fragment = min(target.fragments) if target.fragments else None
    _refresh_support(target)


def _merge_source(df, source, existing=None, show_progress=True):
    """Return a composite-key dictionary; never key evidence by bc1 alone."""
    records = {_key(record): record for record in (existing or {}).values()}
    count_cols = [c for c in ('total_bc1_donor_bc0_counts', 'counts', 'total_counts') if c in df]
    if not count_cols:
        raise ValueError('Missing BDB read count column')
    for required in ('bc1', 'bc0'):
        if required not in df:
            raise ValueError(f'Missing BDB column: {required}')
    iterator = df.itertuples(index=False, name=None)
    if show_progress:
        from tqdm import tqdm
        iterator = tqdm(iterator, total=len(df), desc='Merging observed ' + source + ' links')
    columns = list(df.columns)
    for values in iterator:
        row = dict(zip(columns, values))
        counts = [_count(row[c]) for c in count_cols]
        if len(set(counts)) != 1:
            raise ValueError(f'Conflicting read count columns: {count_cols}')
        count = counts[0]
        bc1, bc0 = _text(row['bc1']), _text(row['bc0'])
        scope = tuple(_text(row.get(c)) for c in SCOPE_COLUMNS)
        fragment = _text(row.get('donor_bc0_fragment'))
        if source == '2x150' and 'donor_bc0_fragment' not in df:
            fragment = _text(row.get('donor'))
        donor = _text(row.get('donor', row.get('bdb_donor')))
        record = BC1Record(bc1, bc0, scope=scope, SPS=_text(row.get('SPS')))
        if source == '2x300':
            record.donor_full = donor
        elif 'donor_bc0_fragment' in df:
            record.donor_recovered = donor
        record.donor_bc0_fragment = fragment
        if fragment:
            record.fragments.add(fragment)
        setattr(record, 'counts_' + source, count)
        setattr(record, 'has_' + source, True)
        if count:
            identities = set()
            for field_name in ('pcr_replicate_ids', 'pcr_replicate_id', 'sample_number'):
                if field_name in row:
                    identities = _ids(row[field_name], serialized=field_name == 'pcr_replicate_ids')
                    if identities:
                        break
            if any(';' in identity or identity.lstrip().startswith('[') for identity in identities):
                raise ValueError('PCR identity contains reserved semicolon or JSON-list prefix')
            getattr(record, 'pcr_ids_' + source).update(identities)
            scalar = row.get('num_PCR_replicates')
            if not identities:
                record.support_unknown = True
                if _text(scalar) is not None:
                    setattr(record, 'anonymous_pcr_' + source, _count(scalar, 'PCR support count'))
            elif _text(scalar) is not None and _count(scalar, 'PCR support count') > len(identities):
                record.support_unknown = True
                setattr(record, 'anonymous_pcr_' + source, _count(scalar, 'PCR support count'))
        _refresh_support(record)
        key = _key(record)
        if key in records:
            _combine(records[key], record)
        else:
            records[key] = record
    return records


def _attach_unique_short_links(records):
    candidates = defaultdict(list)
    for record in records.values():
        if record.best_donor and record.bc1 and record.bc0:
            candidates[(record.scope, record.bc1, record.bc0)].append(record)
    remove = []
    for key, record in records.items():
        if record.best_donor or not record.bc1 or not record.bc0:
            continue
        suffix = _fragment_suffix(record.donor_bc0_fragment, record.SPS, record.bc0)
        if not suffix:
            continue
        matches = [parent for parent in candidates[(record.scope, record.bc1, record.bc0)]
                   if parent.best_donor.endswith(suffix)
                   and parent.SPS == record.SPS]
        record.short_link_candidates = len(matches)
        if len(matches) == 1:
            _combine(matches[0], record)
            remove.append(key)
    for key in remove:
        del records[key]
    return {_key(record): record for record in records.values()}


def merge_2x150_data(df_2x150, show_progress=True):
    """Merge exact observed short links; dictionary keys are scoped link tuples."""
    return _merge_source(df_2x150, '2x150', show_progress=show_progress)


def merge_2x300_data(df_2x300, existing_records=None, show_progress=True):
    """Merge exact full links, then attach only uniquely compatible short links."""
    records = _merge_source(df_2x300, '2x300', existing_records, show_progress)
    return _attach_unique_short_links(records)


def merge_gdb_data(df_gdb, existing_records):
    """Annotate compatible links without multiplying or adding BDB evidence."""
    lookup = defaultdict(list)
    fragments = defaultdict(list)
    for values in df_gdb.itertuples(index=False, name=None):
        row = dict(zip(df_gdb.columns, values))
        donor, bc0 = _text(row.get('donor')), _text(row.get('bc0'))
        if not donor or not bc0:
            continue
        count_col = 'grouped_guide_bc0_counts' if 'grouped_guide_bc0_counts' in row else 'counts'
        candidate = {c: _text(row.get(c)) for c in ('donor', 'bc0', 'guide', 'library_ID', 'scaffold', 'SPS', *SCOPE_COLUMNS)}
        candidate['counts'] = _count(row.get(count_col, 0), 'GDB annotation count')
        lookup[build_donor_bc0_key(donor, bc0)].append(candidate)
        fragment = _text(row.get('donor_bc0_fragment'))
        if fragment:
            fragments[fragment].append(candidate)
    for record in existing_records.values():
        possible = lookup.get(record.donor_bc0_key, [])
        if not possible and record.donor_bc0_fragment:
            possible = fragments.get(record.donor_bc0_fragment, [])
        matches = []
        for candidate in possible:
            if candidate['bc0'] != record.bc0:
                continue
            if record.best_donor and candidate['donor'] != record.best_donor:
                continue
            if record.SPS and candidate.get('SPS') and record.SPS != candidate['SPS']:
                continue
            if any(value and candidate.get(col) and value != candidate[col]
                   for col, value in zip(SCOPE_COLUMNS, record.scope)):
                continue
            suffix = _fragment_suffix(record.donor_bc0_fragment, record.SPS, record.bc0)
            if not record.best_donor and suffix and not candidate['donor'].endswith(suffix):
                continue
            matches.append(candidate)
        # Stable candidate evidence; no counts from a short read are duplicated.
        record.gdb_candidates = sorted(matches, key=lambda row: (-row['counts'], json.dumps(row, sort_keys=True)))
        if (not matches or len({row['donor'] for row in matches}) != 1
                or record.short_link_candidates > 1):
            continue
        best = record.gdb_candidates[0]
        record.has_gdb = True
        record.donor_from_gdb = best['donor']
        guides = {row['guide'] for row in matches if row['guide']}
        # Metadata is one real observation, not independently coalesced fields.
        record.guide = best['guide'] if len(guides) == 1 else None
        record.counts_gdb = best['counts']
        record.library_id, record.scaffold = best['library_ID'], best['scaffold']
    return existing_records


def resolve_donor_conflicts(records):
    """Flag conflicts only where both sequence identities/boundaries are known."""
    for record in records.values():
        donor = record.donor_full or record.donor_recovered
        suffix = _fragment_suffix(record.donor_bc0_fragment, record.SPS, record.bc0)
        if donor and record.donor_from_gdb and donor != record.donor_from_gdb:
            record.has_donor_conflict = True
            record.conflict_details = 'Observed/recovered donor differs from GDB donor'
        if donor and suffix and not donor.endswith(suffix):
            record.has_donor_conflict = True
            record.conflict_details = 'Explicit fragment boundary disagrees with donor suffix'
    return records


def merge_data_sources(df_2x150=None, df_2x300=None, df_gdb=None, resolve_conflicts=True):
    """Return all observed links, including ambiguous/unresolved count evidence.

    This is an evidence table, not a unique bc1-to-construct assignment table.
    The caller must retain scope, PCR identities and unresolved denominator mass.
    """
    records = {}
    if df_2x150 is not None and len(df_2x150):
        records = merge_2x150_data(df_2x150)
    if df_2x300 is not None and len(df_2x300):
        records = merge_2x300_data(df_2x300, records)
    if df_gdb is not None and len(df_gdb):
        records = merge_gdb_data(df_gdb, records)
    if resolve_conflicts:
        records = resolve_donor_conflicts(records)
    columns = list(BC1Record(None, None).to_dict())
    return pd.DataFrame([record.to_dict() for record in records.values()], columns=columns).sort_values(
        'total_counts', ascending=False, kind='stable').reset_index(drop=True)


def filter_by_quality(
    df: pd.DataFrame,
    min_reads: int = 2,
    min_pcr_replicates: int = 1,
    min_purity: Optional[float] = None,
    purity_column: Optional[str] = None,
    mode: str = "label",
) -> pd.DataFrame:
    """
    Filter merged data by quality thresholds.

    2026-09-10 (21_ guard): ``mode="label"`` (DEFAULT) applies the PCR-replicate
    threshold as a ``support_tier`` label and drops nothing on replicates; the
    read and purity thresholds still filter. ``mode="filter"`` is the old
    row-dropping behaviour and must be asked for by name. This function had zero
    callers while armed at 2 by config -- wiring it in would have silently deleted
    the singleton tier.

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

    # PCR replicates: label (default) or filter (explicit)
    if mode not in ("label", "filter"):
        raise ValueError(f"mode must be 'label' or 'filter', got {mode!r}")
    if mode == "label":
        df = df.copy()
        df['support_tier'] = np.where(df['total_pcr_replicates'] >= max(min_pcr_replicates, 1),
                                      'confident', 'provisional')
    elif min_pcr_replicates > 1:
        before = len(df)
        df = df[df['total_pcr_replicates'] >= min_pcr_replicates]
        logger.info(f"After min_pcr_replicates={min_pcr_replicates} filter: {len(df)} records (removed {before - len(df)})")

    # Filter by purity
    if min_purity is not None and purity_column and purity_column in df.columns:
        before = len(df)
        df = df[df[purity_column] >= min_purity]
        logger.info(f"After min_purity={min_purity} filter: {len(df)} records (removed {before - len(df)})")

    return df
