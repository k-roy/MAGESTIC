"""Aggregate observed links while retaining scoped denominator and PCR evidence."""
from collections import defaultdict
import json

import numpy as np
import pandas as pd
from rapidfuzz.distance.Levenshtein import distance

from .merge import (POPULATION_COLUMNS, SCOPE_COLUMNS, _count, _ids, _sample_population, _text,
                    decode_population_counts, encode_population_counts)
from ..readers.data_loader import count_distinct_pcrs

# POPULATION_COLUMNS identify the sampled biological population, trafo_ID being the
# transformation (planning 218 patch 45): a barcode's purity denominator is taken within
# it. V_library/library_ID/SPS are retained as link evidence, never used to partition a
# purity denominator. Each link keeps its reads per sample label in population_counts.


def _flag(value, default=False):
    if value is None or pd.isna(value):
        return default
    if value in (True, 1, 'True', 'true', '1'):
        return True
    if value in (False, 0, 'False', 'false', '0'):
        return False
    raise ValueError(f'Invalid evidence flag: {value!r}')


def _scope(row, columns):
    return tuple(_text(row.get(column)) for column in columns)


def _union(values):
    result = set()
    for value in values:
        result.update(value)
    return result


def _gdb_records(value):
    if not isinstance(value, list):
        if _text(value) is None:
            return []
        value = json.loads(value)
    if not isinstance(value, list) or any(not isinstance(item, dict) for item in value):
        raise ValueError('GDB candidates must be a JSON list of evidence records')
    return value


def _combine(group, count_columns, representative=None):
    """Retain a real representative and combine only declared evidence fields."""
    records = group.to_dict('records')
    row = dict(representative if representative is not None else
               max(records, key=lambda item: item['_representative_reads']))
    for column in count_columns:
        row[column] = sum(item[column] for item in records)
    populations = {}
    for item in records:
        for population, count in item['_populations'].items():
            populations[population] = populations.get(population, 0) + count
    row['_populations'] = populations
    row['population_counts'] = encode_population_counts(populations)
    for column in ('_pcr_ids', '_original_bc1s', '_guide_candidates', '_pcr_ids_2x150', '_pcr_ids_2x300', '_gdb_candidates'):
        row[column] = _union(item[column] for item in records)
    for column in ('_pcr_lower_bound', '_pcr_lower_bound_2x150', '_pcr_lower_bound_2x300'):
        row[column] = max(item[column] for item in records)
    row['pcr_support_exact'] = all(item['pcr_support_exact'] for item in records if item['counts'] > 0)
    for column in ('has_2x150', 'has_2x300', 'has_gdb', 'has_donor_conflict', 'gdb_ambiguous'):
        if column in group:
            row[column] = any(_flag(item[column]) for item in records)
    if 'gdb_candidates' in group:
        row['gdb_candidates'] = '[' + ', '.join(sorted(row['_gdb_candidates'])) + ']'
        row['gdb_ambiguous'] = _flag(row.get('gdb_ambiguous')) or len(row['_gdb_candidates']) > 1
    if 'guide' in group:
        guides = row['_guide_candidates']
        row['guide_candidates'] = json.dumps(sorted(guides))
        conflict = len(guides) > 1 or any(_flag(item.get('guide_evidence_conflict')) for item in records)
        row['guide_evidence_conflict'] = conflict
        if conflict:
            row['guide'] = None
        elif len(guides) == 1 and _text(row.get('guide')) is None:
            # Carry a whole annotation bundle only for the same observed donor.
            carriers = [item for item in records if _text(item.get('guide'))
                        and _text(item.get('donor')) == _text(row.get('donor'))]
            if carriers:
                carrier = max(carriers, key=lambda item: item['_representative_reads'])
                # Structured evidence belongs to the whole group. Only scalar
                # display metadata comes from this coherent carrier record.
                for column in ('guide', 'donor_from_gdb', 'library_id', 'scaffold', 'counts_gdb'):
                    if column in carrier:
                        row[column] = carrier[column]
    return row


def _collapse(frame, group_columns, sequence_column, count_columns, threshold):
    """Greedy existing-distance policy, preserving each row's evidence bundle."""
    rows = []
    for _, group in frame.groupby(group_columns, dropna=False, sort=False):
        group = group.reset_index(drop=True)
        while len(group):
            top_index = group['counts'].idxmax()
            representative = group.loc[top_index].to_dict()
            sequence = _text(representative[sequence_column])
            if sequence is None:
                # Unknown sequences cannot be error-collapsed; still consume the row.
                mask = group.index == top_index
            else:
                mask = group[sequence_column].map(
                    lambda value: _text(value) is not None and distance(str(value), sequence) <= threshold)
            rows.append(_combine(group.loc[mask], count_columns, representative))
            group = group.loc[~mask].reset_index(drop=True)
    return pd.DataFrame(rows, columns=frame.columns.union(pd.Index(['guide_candidates', 'guide_evidence_conflict']), sort=False))


def aggregate_link_evidence(df, min_reads=2, min_pcr_replicates=2, collapse_sequences=True,
                            show_progress=True):
    """Keep all count evidence until scoped barcode denominators are calculated.

    A raw barcode sent to multiple representatives by separate link collapses has
    ambiguous identity: its affected purity is undefined rather than invented.
    Min-PCR support is a label. Only the explicit final min-read rule drops rows.
    Inputs must be complete observed-link evidence, not a filtered prior output.
    """
    frame = df.copy()
    if any(column in frame for column in ('original_bc1s', 'total_bc1_counts', 'bc1_purity')):
        raise ValueError('Cannot reaggregate a prior aggregation output; use complete observed-link evidence')
    for column in ('bc1', 'donor', 'bc0'):
        if column not in frame:
            if len(frame):
                raise ValueError(f'Missing link identity column: {column}')
            frame[column] = pd.Series(dtype=object)
    if 'counts' not in frame:
        if 'total_counts' in frame:
            frame['counts'] = frame['total_counts']
        else:
            source_columns = [column for column in ('counts_2x150', 'counts_2x300') if column in frame]
            if not source_columns and len(frame):
                raise ValueError('Missing BDB read count column')
            frame['counts'] = frame[source_columns].sum(axis=1) if source_columns else pd.Series(dtype='int64')
    count_columns = [column for column in ('counts', 'total_counts', 'counts_2x150', 'counts_2x300') if column in frame]
    for column in count_columns:
        frame[column] = frame[column].map(_count)
    for column in ('bc1', 'donor', 'bc0'):
        frame[column] = frame[column].map(_text)
    scopes = [column for column in POPULATION_COLUMNS if column in frame]
    # SPS identifies construct evidence, not an independently transformed pool.
    # It can constrain link collapsing but cannot split a barcode denominator.
    link_scopes = [column for column in (*SCOPE_COLUMNS, 'input_library_scope', 'SPS') if column in frame]
    if frame.empty:
        for column, dtype in [('total_bc1_counts', float), ('bc1_purity', float),
                              ('unresolved_bc1_counts', float), ('num_distinct_donors', int),
                              ('num_PCR_replicates', int), ('support_tier', object),
                              ('pcr_replicate_ids', object), ('barcode_identity_ambiguous', bool),
                              ('population_counts', object)]:
            frame[column] = pd.Series(dtype=dtype)
        return frame

    original_totals = defaultdict(int)
    original_unresolved = defaultdict(int)
    ids, lower, exact, origins, guides, gdb_evidence, populations = [], [], [], [], [], [], []
    source_ids = {source: [] for source in ('2x150', '2x300')}
    source_lower = {source: [] for source in source_ids}
    source_fields = {source for source in source_ids if
                     any(column in frame for column in ('pcr_ids_' + source, 'pcr_replicates_' + source))}
    for row in frame.to_dict('records'):
        key = (_scope(row, scopes), row['bc1'])
        original_totals[key] += row['counts']
        if row['donor'] is None or row['bc0'] is None:
            original_unresolved[key] += row['counts']
        identities = set()
        if row['counts']:
            for column in ('pcr_replicate_ids', 'pcr_replicate_id', 'sample_number'):
                if column in row:
                    identities = _ids(row[column], serialized=column == 'pcr_replicate_ids')
                    if identities:
                        break
        source_bounds = []
        for source in source_ids:
            supported = row['counts'] > 0 and row.get('counts_' + source, row['counts']) > 0
            values = _ids(row.get('pcr_ids_' + source), serialized=True) if supported else set()
            source_ids[source].append(values)
            bound_value = row.get('pcr_replicates_' + source)
            bound = _count(bound_value, 'source PCR support count') if supported and _text(bound_value) is not None else 0
            source_lower[source].append(bound)
            source_bounds.append(bound)
            identities.update(values)
        if any(';' in identity or identity.lstrip().startswith('[') for identity in identities):
            raise ValueError('PCR identity contains reserved semicolon or JSON-list prefix')
        scalar = max([_count(row[column], 'PCR support count') for column in
                      ('num_PCR_replicates', 'total_pcr_replicates')
                      if column in row and _text(row[column]) is not None] + source_bounds + [0]) if row['counts'] else 0
        ids.append(identities)
        lower.append(scalar)
        exact.append(bool((identities and scalar <= len(identities) and _flag(row.get('pcr_support_exact'), True))
                          or row['counts'] == 0))
        origins.append({row['bc1']} if row['bc1'] is not None else set())
        guide_values = {_text(row.get('guide'))} - {None}
        guide_values.update(_ids(row.get('guide_candidates'), serialized=True))
        candidates = _gdb_records(row.get('gdb_candidates'))
        guide_values.update({_text(item.get('guide')) for item in candidates} - {None})
        gdb_evidence.append({json.dumps(item, sort_keys=True) for item in candidates})
        guides.append(guide_values)
        observed = decode_population_counts(row.get('population_counts'))
        if observed and sum(observed.values()) != row['counts']:
            raise ValueError('population_counts does not sum to the row read count')
        scope = tuple(_text(row.get(column)) for column in POPULATION_COLUMNS)
        populations.append(observed or {_sample_population(scope, _text(row.get('trafo_label'))): row['counts']})
    frame['_populations'] = populations
    frame['_pcr_ids'], frame['_pcr_lower_bound'] = ids, lower
    frame['pcr_support_exact'], frame['_original_bc1s'] = exact, origins
    frame['_guide_candidates'] = guides
    frame['_gdb_candidates'] = gdb_evidence
    for source in source_ids:
        frame['_pcr_ids_' + source] = source_ids[source]
        frame['_pcr_lower_bound_' + source] = source_lower[source]
    frame['_representative_reads'] = frame['counts']
    frame['_unresolved_fragment'] = [row.get('donor_bc0_fragment') if row['donor'] is None else None
                                      for row in frame.to_dict('records')]
    groups = link_scopes + ['bc1', 'donor', 'bc0', '_unresolved_fragment']
    result = pd.DataFrame([_combine(group, count_columns) for _, group in frame.groupby(groups, dropna=False, sort=False)])
    if collapse_sequences:
        # Unknown link components cannot establish that two observations belong
        # to the same construct. Keep their exact rows and read mass intact.
        complete = result[['bc1', 'donor', 'bc0']].notna().all(axis=1)
        resolved = result.loc[complete]
        unresolved_links = result.loc[~complete]
        resolved = _collapse(resolved, link_scopes + ['bc1', 'bc0'], 'donor', count_columns, 5)
        resolved = _collapse(resolved, link_scopes + ['bc0', 'donor'], 'bc1', count_columns, 1)
        result = pd.concat([resolved, unresolved_links], ignore_index=True)

    # Map original identity to final identity before applying the read filter.
    destinations = defaultdict(set)
    destination_origins = defaultdict(set)
    for row in result.to_dict('records'):
        scope = _scope(row, scopes)
        for original in row['_original_bc1s']:
            destinations[(scope, original)].add(row['bc1'])
            destination_origins[(scope, row['bc1'])].add(original)
    totals, unresolved, ambiguities = [], [], []
    for row in result.to_dict('records'):
        scope = _scope(row, scopes)
        originals = destination_origins[(scope, row['bc1'])]
        ambiguous = row['bc1'] is None or any(len(destinations[(scope, original)]) != 1 for original in originals)
        ambiguities.append(ambiguous)
        totals.append(np.nan if ambiguous else sum(original_totals[(scope, original)] for original in originals))
        unresolved.append(np.nan if ambiguous else sum(original_unresolved[(scope, original)] for original in originals))
    result['total_bc1_counts'] = totals
    result['unresolved_bc1_counts'] = unresolved
    result['barcode_identity_ambiguous'] = ambiguities
    result['bc1_purity'] = result['counts'] / result['total_bc1_counts'].replace(0, np.nan)
    result['num_distinct_donors'] = result.groupby(scopes + ['bc1'], dropna=False)['donor'].transform('nunique')
    # Identity strings from routes that declare different fields can name one PCR, so
    # support is the fewest PCRs consistent with them; an inferred merge is inexact.
    support = [count_distinct_pcrs(value) for value in result['_pcr_ids']]
    per_source = {source: [max(count_distinct_pcrs(value)[0], bound) for value, bound in
                           zip(result['_pcr_ids_' + source], result['_pcr_lower_bound_' + source])]
                  for source in source_fields}
    # An exact minimum cover never shrinks when identities are added, but the search's
    # fallback lower bound can, so the total is floored at every per-source count, as
    # BC1Record._pcr_support does in merge.
    result['num_PCR_replicates'] = [max([count, bound] + [per_source[source][i] for source in per_source])
                                    for i, ((count, _), bound) in enumerate(zip(support, result['_pcr_lower_bound']))]
    result['pcr_support_exact'] = [bool(flag) and identities_exact for flag, (_, identities_exact)
                                   in zip(result['pcr_support_exact'], support)]
    result['total_pcr_replicates'] = result['num_PCR_replicates']
    result['pcr_replicate_ids'] = [';'.join(sorted(value)) for value in result['_pcr_ids']]
    if 'pcr_replicate_id' in result:
        result['pcr_replicate_id'] = [next(iter(value)) if len(value) == 1 else None for value in result['_pcr_ids']]
    for source in source_fields:
        result['pcr_replicates_' + source] = per_source[source]
        result['pcr_ids_' + source] = [json.dumps(sorted(value)) for value in result['_pcr_ids_' + source]]
    result['original_bc1s'] = [json.dumps(sorted(value)) for value in result['_original_bc1s']]
    result['support_tier'] = np.where(result['num_PCR_replicates'] >= max(min_pcr_replicates, 1), 'confident', 'provisional')
    if sum(result['counts']) != sum(frame['counts']):
        raise ValueError('BDB aggregation did not conserve input read counts')
    removed = result[result['counts'] < min_reads]
    kept = result[result['counts'] >= min_reads].copy()
    kept.attrs['count_accounting'] = dict(reads_in=int(sum(frame['counts'])),
        reads_out=int(sum(kept['counts'])), rows_below_min_reads=len(removed),
        reads_below_min_reads=int(sum(removed['counts'])), rows_missing_bc1=int(frame['bc1'].isna().sum()))
    return kept.drop(columns=[column for column in kept if column.startswith('_')])
