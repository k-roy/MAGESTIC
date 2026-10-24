"""Candidate-preserving recovery from a separate plasmid donor/bc0 amplicon, ranked and tiered by rule R*."""
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


#: Rule R* (planning/222_triumph_embl_library_design/267_bc1_oligo_assignment_audit.md section 20.4; Kevin
#: 2026-10-24, planning/259 section 7). The support floor and the lead ratio are TIER BOUNDARIES, not filters:
#: a lead below them is still assigned, as `provisional`. R* was flat across floors {5, 10, 20} x ratios
#: {2, 3, 5, 10}, so these are the pre-registered point, not a tuned optimum.
RSTAR_FLOOR_READS = 10
RSTAR_LEAD_RATIO = 3


def _tier(lead, runner_up, floor, ratio):
    return 'confident' if lead >= floor and lead >= ratio * runner_up else 'provisional'


def recover_ranked_donors(observations, candidates, *, floor=RSTAR_FLOOR_READS, ratio=RSTAR_LEAD_RATIO,
                          is_perfect_donor=None):
    """Annotate once per input row with rule R*: every candidate kept, the support leader assigned, tiered.

    Assignment (R* S1+S2+S6, unconditional): the donor with strictly the most reference-amplicon reads for
    the row's (fragment, bc0, scope) key. Reads whose donor is unknown compete as one anonymous candidate:
    they can block an unclear lead but can never be assigned. An exact tie at the top is not guessed.

    Tier (R* S3+S4): `confident` when the lead holds at least `floor` reads and `ratio` times the runner-up.
    With `is_perfect_donor` (a predicate on a candidate donor; the BDB route passes the oligo lookup's
    perfect-donor test), the tier is computed on perfect-donor reads only, and is `confident` only when the
    perfect-donor leader is also the assigned donor. Without it, the tier is computed on all reads.

    Candidate counts are separate reference-amplicon read support. They never enter an observation's counts
    or PCR support, and an observation is never duplicated. Known scopes match exactly; absent scope is not a
    wildcard for another library or transformation, and a blank scope value matches only a blank one.
    """
    output = observations.copy()
    if 'donor' not in output:
        output['donor'] = None
    if len(output) and len(candidates):
        for column in ('donor_bc0_fragment', 'bc0'):
            if column not in output:
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

    grouped = defaultdict(lambda: dict(support=defaultdict(int), unresolved=0, reads=0))
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
            group['support'][donor] += count

    perfect_cache = {}

    def perfect(donor):
        if donor not in perfect_cache:
            perfect_cache[donor] = bool(is_perfect_donor(donor))
        return perfect_cache[donor]

    empty = dict(support={}, unresolved=0, reads=0)
    columns = defaultdict(list)
    for row in output.to_dict('records'):
        group = grouped.get(key(row), empty)
        # Total order: most-supported first, then the donor string, so the result is input-order free.
        ranked = sorted(group['support'].items(), key=lambda item: (-item[1], item[0]))
        lead = ranked[0][1] if ranked else 0
        runner_up = max(ranked[1][1] if len(ranked) > 1 else 0, group['unresolved'])
        observed = _text(row.get('donor'))
        donor, tier, perfect_reads = observed, 'unassigned', None
        if observed is not None:
            status = ('conflict_with_observed' if any(value != observed for value, _ in ranked)
                      else 'observed_consistent' if ranked else 'observed_only')
            tier = 'observed'
        elif not ranked:
            status = 'unresolved_reference' if group['unresolved'] else 'no_positive_candidate'
        elif lead > runner_up:
            donor = ranked[0][0]
            status = 'recovered_unique' if len(ranked) == 1 and not group['unresolved'] else 'recovered_by_support'
            if is_perfect_donor is None:
                tier = _tier(lead, runner_up, floor, ratio)
            else:
                exact = sorted(((d, n) for d, n in ranked if perfect(d)), key=lambda item: (-item[1], item[0]))
                perfect_reads = dict(exact).get(donor, 0)
                second = exact[1][1] if len(exact) > 1 else 0
                tier = ('confident' if exact and exact[0][0] == donor and exact[0][1] > second
                        and _tier(exact[0][1], second, floor, ratio) == 'confident' else 'provisional')
        elif len(ranked) > 1 and ranked[1][1] == lead:
            status = 'ambiguous'
        else:
            status = 'unresolved_reference'
        columns['donor'].append(donor)
        columns['donor_recovery_status'].append(status)
        columns['donor_recovery_tier'].append(tier)
        columns['donor_recovery_candidates'].append(json.dumps(sorted(group['support'])))
        columns['donor_recovery_candidate_support'].append(json.dumps([[d, n] for d, n in ranked]))
        columns['donor_recovery_lead_reads'].append(lead)
        columns['donor_recovery_runner_up_reads'].append(runner_up)
        columns['donor_recovery_support_reads'].append(group['reads'])
        columns['donor_recovery_unresolved_reads'].append(group['unresolved'])
        columns['donor_recovery_perfect_reads'].append(perfect_reads)
    for name in ('donor', 'donor_recovery_status', 'donor_recovery_tier', 'donor_recovery_candidates',
                 'donor_recovery_candidate_support', 'donor_recovery_lead_reads',
                 'donor_recovery_runner_up_reads', 'donor_recovery_support_reads',
                 'donor_recovery_unresolved_reads'):
        output[name] = columns[name]
    if is_perfect_donor is not None:
        # Emitted only when it was measured: a column that is constant on every row names nothing.
        output['donor_recovery_perfect_reads'] = columns['donor_recovery_perfect_reads']
    return output
