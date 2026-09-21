"""Candidate-preserving recovery from a separate plasmid donor/bc0 amplicon."""
from collections import defaultdict
import json

import pandas as pd

from .merge import SCOPE_COLUMNS, _count, _text


def aggregate_donor_candidates(frame):
    """Sum only identical observed candidates; never pick the first donor."""
    required = ('donor_bc0_fragment', 'donor', 'bc0', 'counts')
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f'Missing donor-recovery evidence columns: {sorted(missing)}')
    frame = frame.copy()
    frame['counts'] = frame['counts'].map(_count)
    keys = [column for column in (*SCOPE_COLUMNS, 'SPS') if column in frame]
    keys += ['donor_bc0_fragment', 'donor', 'bc0']
    aggregations = {'counts': lambda values: sum(int(value) for value in values)}
    if 'sample_number' in frame:
        aggregations['sample_number'] = 'nunique'
    if 'source_file' in frame:
        aggregations['source_file'] = lambda values: json.dumps(sorted(set(values)))
    result = frame.groupby(keys, dropna=False, sort=False).agg(aggregations).reset_index()
    return result.rename(columns={'sample_number': 'num_samples_donor_bc0',
                                  'source_file': 'donor_recovery_source_files'})


def recover_unique_donors(observations, candidates):
    """Annotate once per input row; unknown/conflicting evidence stays visible.

    Candidate counts are separate reference-amplicon read support. They never
    enter an observation's counts or PCR support. Known scopes match exactly;
    absent scope is not a wildcard for another library or transformation.
    """
    output = observations.copy()
    if 'donor' not in output:
        output['donor'] = None
    for column in ('donor_bc0_fragment', 'bc0'):
        if len(output) and column not in output:
            raise ValueError(f'Missing BDB recovery key: {column}')
    if len(candidates):
        missing = {'donor_bc0_fragment', 'bc0', 'donor', 'counts'} - set(candidates.columns)
        if missing:
            raise ValueError(f'Missing donor-recovery evidence columns: {sorted(missing)}')
    scopes = [column for column in (*SCOPE_COLUMNS, 'SPS')
              if column in output or column in candidates]

    def key(row):
        fragment, bc0 = _text(row.get('donor_bc0_fragment')), _text(row.get('bc0'))
        if fragment is None or bc0 is None:
            return None
        return (fragment, bc0, *(_text(row.get(column)) for column in scopes))

    grouped = defaultdict(lambda: dict(donors=set(), unresolved=0, reads=0))
    for row in candidates.to_dict('records'):
        count = _count(row['counts'])
        identity = key(row)
        if identity is None or count == 0:
            continue
        group = grouped[identity]
        group['reads'] += count
        donor = _text(row['donor'])
        if donor is None:
            group['unresolved'] += count
        else:
            group['donors'].add(donor)

    donors, statuses, encoded, totals, unresolved = [], [], [], [], []
    for row in output.to_dict('records'):
        group = grouped.get(key(row), dict(donors=set(), unresolved=0, reads=0))
        alternatives = sorted(group['donors'])
        observed = _text(row.get('donor'))
        donor = observed
        if observed is not None:
            status = ('conflict_with_observed' if any(value != observed for value in alternatives)
                      else 'observed_consistent' if alternatives else 'observed_only')
        elif group['unresolved']:
            status = 'unresolved_reference'
        elif len(alternatives) > 1:
            status = 'ambiguous'
        elif len(alternatives) == 1:
            donor, status = alternatives[0], 'recovered_unique'
        else:
            status = 'no_positive_candidate'
        donors.append(donor); statuses.append(status)
        encoded.append(json.dumps(alternatives))
        totals.append(group['reads']); unresolved.append(group['unresolved'])
    output['donor'] = donors
    output['donor_recovery_status'] = statuses
    output['donor_recovery_candidates'] = encoded
    output['donor_recovery_support_reads'] = totals
    output['donor_recovery_unresolved_reads'] = unresolved
    return output
