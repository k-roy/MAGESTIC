"""
Data loader module for BC1-Donor-BC0 pipeline.

Handles loading processed bc1-donor-bc0 data from various sources:
- 2x150 bp sequencing (partial donor in bc1-donor-bc0, full donor via separate donor-bc0)
- 2x300 bp sequencing (full donor in merged reads)
- Guide-donor-bc0 (gdb) data

All loaders return DataFrames with standardized columns for downstream processing.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import json
import re
from urllib.parse import quote, unquote
import pandas as pd
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


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


# Cell spellings that mean "no value": pandas' default ``na_values``. The loaders read
# keyfiles with default NA parsing and numbered step 02 reads literal strings
# (keep_default_na=False); normalizing both through this list makes them agree.
MISSING_VALUE_SPELLINGS = frozenset((
    '', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND',
    '1.#QNAN', '<NA>', 'N/A', 'NA', 'NULL', 'NaN', 'None', 'n/a', 'nan', 'null'))


def declared_value(value):
    """Return a metadata cell as stripped text, or None when it holds no value.

    An integer-valued float (``1.0``, which pandas produces for an integer column
    containing a blank) is rendered as ``1``.
    """
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    text = str(value).strip()
    return None if text in MISSING_VALUE_SPELLINGS else text


def build_pcr_replicate_id(row) -> str:
    """PCR-replicate identity = one independent amplification of the gDNA, per well.

    BDB reads the guide-donor plasmid carried in that gDNA prep, so this identifies
    a PCR event, not a genomic integration.

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

    Every producer (both loaders and numbered step 02) builds the identity here,
    and the cells are normalized here so that producers agree. Surrounding
    whitespace is removed. An absent field, or a cell spelled as in
    :data:`MISSING_VALUE_SPELLINGS` (including a literal ``"NA"``), is missing and
    encodes as JSON null. ``sample_number`` drops a ``sample_`` prefix and leading
    zeros, so ``sample_01`` and ``1`` are the same well. A missing
    ``amplicon_level`` is derived from ``trafo_ID`` with
    :func:`derive_amplicon_level` (``unknown`` when no trafo_ID is declared),
    exactly as the loaders derive it.

    The ordered tuple is serialized as JSON behind a ``pcr:v1:`` prefix and
    percent-escaped, so field boundaries survive underscores and the ``;`` list
    delimiter used by merge and aggregation. When ``gDNA_plate_name`` is
    unavailable the re-prep guard is inactive; :func:`check_replicate_metadata`
    warns so this is never silent.

    Args:
        row: A keyfile row (pandas Series or mapping).

    Returns:
        A ``pcr:v1:`` identifier. Rebuild all identities together: a legacy
        underscore-joined id and a ``pcr:v1:`` id for the same physical PCR are
        different strings and would be counted as two replicates.
    """
    def _get(key):
        if hasattr(row, 'get'):
            return row.get(key)
        return row[key] if key in row else None

    values = {key: declared_value(_get(key)) for key in PCR_REPLICATE_ID_FIELDS}
    well = values['sample_number']
    if well is not None:
        well = re.sub(r'^sample_', '', well)
        values['sample_number'] = str(int(well)) if well.isdigit() else well
    if values['amplicon_level'] is None:
        values['amplicon_level'] = derive_amplicon_level(declared_value(_get('trafo_ID')))
    encoded = json.dumps([values[key] for key in PCR_REPLICATE_ID_FIELDS],
                         separators=(',', ':'))
    return 'pcr:v1:' + quote(encoded, safe='')


_PCR_V1_PREFIX = 'pcr:v1:'
_AMPLICON_LEVEL_POSITION = PCR_REPLICATE_ID_FIELDS.index('amplicon_level')
# Branch-and-bound steps allowed for one group of compatible identities. Past this the
# count falls back to a proven lower bound and is reported as inexact.
PCR_COVER_SEARCH_LIMIT = 100_000
# Components larger than this skip the exact search (it also bounds recursion depth).
_PCR_COVER_EXACT_MAX = 40


@lru_cache(maxsize=65536)
def pcr_identity_fields(identity: str) -> Optional[Tuple[Optional[str], ...]]:
    """Invert :func:`build_pcr_replicate_id`.

    Returns the :data:`PCR_REPLICATE_ID_FIELDS` values of a ``pcr:v1:`` identity, with
    None for a field the producer did not declare. Any other encoding (a legacy
    underscore-joined id, a bare sample_number, a test label) returns None and is
    treated by callers as an opaque identity, equal only to itself.
    """
    if not isinstance(identity, str) or not identity.startswith(_PCR_V1_PREFIX):
        return None
    try:
        values = json.loads(unquote(identity[len(_PCR_V1_PREFIX):]))
    except ValueError:
        return None
    if (not isinstance(values, list) or len(values) != len(PCR_REPLICATE_ID_FIELDS)
            or any(value is not None and not isinstance(value, str) for value in values)):
        return None
    return tuple(values)


def _declared_pcr_fields(identity):
    """The fields an identity declares; None where nothing was declared.

    A derived ``amplicon_level`` of ``unknown`` means no trafo_ID was declared, so it
    counts as undeclared too.
    """
    fields = pcr_identity_fields(identity)
    if fields is not None and fields[_AMPLICON_LEVEL_POSITION] == AMPLICON_LEVEL_UNKNOWN:
        fields = (fields[:_AMPLICON_LEVEL_POSITION] + (None,)
                  + fields[_AMPLICON_LEVEL_POSITION + 1:])
    return fields


def _can_be_one_pcr(a, b) -> bool:
    """Two declarations can describe one PCR when every field both declare agrees."""
    return all(x is None or y is None or x == y for x, y in zip(a, b))


def _fields_order(fields):
    return tuple((value is None, value or '') for value in fields)


def _signature(fields):
    return tuple(value is None for value in fields)


def _proven_lower_bound(vertices) -> int:
    """Size of a set of pairwise incompatible identities: each needs its own PCR."""
    degree = {v: sum(_can_be_one_pcr(v, w) for w in vertices if w != v) for v in vertices}
    chosen = []
    for fields in sorted(vertices, key=lambda v: (degree[v], _fields_order(v))):
        if not any(_can_be_one_pcr(fields, other) for other in chosen):
            chosen.append(fields)
    # Identities declaring the same fields are pairwise incompatible, so the largest such
    # class is also a set of provably distinct PCRs.
    return max(len(chosen), max(Counter(map(_signature, vertices)).values()))


def _fewest_pcrs(vertices) -> int:
    """Minimum number of groups of mutually compatible identities (a minimum clique cover).

    Exact branch and bound: each identity joins a group it is compatible with or opens a
    new one, and a branch stops once it cannot beat the best cover found. Being
    compatible with a group's merged declaration is the same as being compatible with
    every member, so no valid grouping is missed. When the search budget runs out, a
    proven lower bound is returned instead of the best cover found so far, because a
    cover found by an incomplete search can exceed the true minimum.
    """
    order = sorted(vertices, key=lambda v: (-sum(x is not None for x in v), _fields_order(v)))
    best, steps = len(order), 0

    def place(index, groups):
        nonlocal best, steps
        steps += 1
        if steps > PCR_COVER_SEARCH_LIMIT or len(groups) >= best:
            return
        if index == len(order):
            best = len(groups)
            return
        fields = order[index]
        for position, merged in enumerate(groups):
            if _can_be_one_pcr(fields, merged):
                groups[position] = tuple(a if a is not None else b for a, b in zip(merged, fields))
                place(index + 1, groups)
                groups[position] = merged
        groups.append(fields)
        place(index + 1, groups)
        groups.pop()

    if len(order) <= _PCR_COVER_EXACT_MAX:
        place(0, [])
        if steps <= PCR_COVER_SEARCH_LIMIT:
            return best
    return _proven_lower_bound(order)


def count_distinct_pcrs(identities: Iterable[str]) -> Tuple[int, bool]:
    """Fewest physical PCRs consistent with a set of PCR identities, and whether that is exact.

    Two ``pcr:v1:`` identities can be one physical PCR when they agree on every field
    both declare: a field one producer did not declare is unknown, not different. This
    matters when one PCR is sequenced by two routes whose keyfiles declare different
    fields (the yL437/yL442 2x150 combined key has no gDNA plate, inner primers or
    polymerase; its 2x300 plate keys do). Counting the identity strings would count that
    one PCR twice.

    The count is the minimum number of PCRs that explains every identity, so it is a
    lower bound and never raises support above what the declarations prove. It is exact
    when no two distinct identities could be the same PCR. When a merge had to be
    inferred, the under-declared PCR may still be a different one, so the second value
    is False. Identities in any other encoding are opaque and count once each.

    Returns:
        ``(count, exact)``.
    """
    atoms, declared = 0, set()
    for identity in set(identities):
        fields = _declared_pcr_fields(identity)
        if fields is None:
            atoms += 1
        else:
            declared.add(fields)
    by_signature = defaultdict(list)
    for fields in declared:
        by_signature[_signature(fields)].append(fields)
    # Within one signature, compatible means equal, so only pairs across signatures
    # need checking.
    if len(by_signature) < 2:
        return atoms + len(declared), True
    vertices = sorted(declared, key=_fields_order)
    index = {fields: i for i, fields in enumerate(vertices)}
    parent = list(range(len(vertices)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    groups = sorted(by_signature.values(), key=lambda group: _fields_order(min(group, key=_fields_order)))
    for g, first in enumerate(groups):
        for second in groups[g + 1:]:
            for a in first:
                for b in second:
                    if _can_be_one_pcr(a, b):
                        parent[root(index[a])] = root(index[b])
    components = defaultdict(list)
    for fields in vertices:
        components[root(index[fields])].append(fields)
    count, exact = atoms, True
    for members in components.values():
        if len(members) == 1:
            count += 1
        else:
            # A component of two or more contains a compatible pair, so a merge is inferred.
            exact = False
            count += _fewest_pcrs(members)
    return count, exact


# Declarations that identify the sample, population or PCR event. When both the plate
# key and the sample key carry one, they must agree for every file that is loaded.
_KEY_IDENTITY_COLUMNS = frozenset((*PCR_REPLICATE_ID_FIELDS, 'sequencing_date',
                                   'yeast_library', 'yeast_sublibrary', 'trafo_ID',
                                   'library_ID', 'V_library', 'SPS'))
_KEY_CONFLICTS = '_key_conflicts'
# Columns that name the physical file ({date}_{outer}_sample_{n}_...). When both keys
# carry one, a plate row and a sample row that disagree on it describe different files.
_KEY_FILE_COLUMNS = ('sequencing_date', 'outer_primers', 'sample_number')
_PLAIN_NUMBER = re.compile(r'^-?\d+(\.\d+)?$')


def _same_declaration(a: str, b: str) -> bool:
    """Equal declarations; plain decimal numerals compare by value (1 and 1.0 agree).

    Compared as decimals, not floats: two integers past float precision (for example
    9007199254740992 and 9007199254740993) are different declarations.
    """
    if a == b:
        return True
    if _PLAIN_NUMBER.match(a) and _PLAIN_NUMBER.match(b):
        try:
            return Decimal(a) == Decimal(b)
        except InvalidOperation:
            return False
    return False


def _cross_merge_keys(plate_key: pd.DataFrame, sample_key: pd.DataFrame) -> pd.DataFrame:
    """Pair plate-key rows with sample-key rows, resolving columns that both declare.

    A plain cross-merge renames a shared column to ``_x``/``_y``, so an equal
    declaration in both keys vanishes and a contradictory one passes unnoticed.

    * A shared file-naming column (:data:`_KEY_FILE_COLUMNS`) is a join key: a plate
      row and a sample row that disagree on it describe different files and are not
      paired. Rows left with no partner are logged; a file on disk that no pair
      names is then reported by the loader's coverage check.
    * Every other shared column becomes one column holding the value both keys agree
      on, or one key's value when the other is blank.
    * A disagreement in an identity column (:data:`_KEY_IDENTITY_COLUMNS`) is listed
      in ``_key_conflicts``. The loader combines the lists of all rows that name one
      file and refuses that file when it loads. Other disagreeing columns keep both
      values as ``<column>__plate_key`` / ``<column>__sample_key`` and leave the
      shared column empty.

    Resolved shared values are stripped text (:func:`declared_value`), so a blank
    never turns an integer date column into floats.
    """
    shared = sorted(set(plate_key.columns) & set(sample_key.columns))
    combined = plate_key.rename(columns={c: f'{c}__plate_key' for c in shared}).assign(
        _plate_row=range(len(plate_key))).merge(
        sample_key.rename(columns={c: f'{c}__sample_key' for c in shared}).assign(
            _sample_row=range(len(sample_key))), how='cross')
    paired = [True] * len(combined)
    conflicts = [[] for _ in range(len(combined))]
    unresolved = set()
    for column in shared:
        resolved = []
        pairs = zip(combined[f'{column}__plate_key'], combined[f'{column}__sample_key'])
        for i, (plate_value, sample_value) in enumerate(pairs):
            a, b = declared_value(plate_value), declared_value(sample_value)
            if a is None or b is None or _same_declaration(a, b):
                resolved.append(b if b is not None else a)
                continue
            resolved.append(None)
            if column in _KEY_FILE_COLUMNS:
                paired[i] = False
            elif column in _KEY_IDENTITY_COLUMNS:
                conflicts[i].append(column)
            else:
                unresolved.add(column)
        combined[column] = pd.Series(resolved, index=combined.index, dtype=object)
    combined[_KEY_CONFLICTS] = [';'.join(c) for c in conflicts]
    kept = combined[paired]
    if len(kept) < len(combined):
        join_columns = [c for c in shared if c in _KEY_FILE_COLUMNS]
        logger.info(f"  {len(combined) - len(kept):,} plate/sample row combinations disagree on "
                    f"{join_columns} and are not paired")
        for side, other, column in (('sample', 'plate', '_sample_row'), ('plate', 'sample', '_plate_row')):
            lost = sorted(set(combined[column]) - set(kept[column]))
            if lost:
                logger.warning(f"  {len(lost)} {side}_key row(s) pair with no {other}_key row "
                               f"(0-based rows {lost[:10]})")
    if shared:
        logger.info(f"  plate_key and sample_key both declare {shared}; resolved to one column each")
    if unresolved:
        logger.warning(
            f"  plate_key and sample_key disagree on non-identity column(s) {sorted(unresolved)}; "
            f"both values kept as <column>__plate_key / <column>__sample_key")
    return kept.drop(columns=['_plate_row', '_sample_row']).reset_index(drop=True)


def _declared_or_derived_amplicon_level(key: pd.DataFrame) -> pd.Series:
    """The key's declared amplicon_level, else the level derived from trafo_ID."""
    declared = (key['amplicon_level'].map(declared_value) if 'amplicon_level' in key.columns
                else pd.Series(None, index=key.index, dtype=object))
    derived = (key['trafo_ID'].map(derive_amplicon_level) if 'trafo_ID' in key.columns
               else pd.Series(AMPLICON_LEVEL_UNKNOWN, index=key.index, dtype=object))
    return declared.where(declared.notna(), derived)


# Transformation scope (planning 218 patch 45; Kevin, 2026-09-23). The bc1 PCR fragment is joined to
# the I-SceI-cut guide-donor plasmid in yeast at each transformation (utils/harmonize_oligos.py
# construction chain), so a bc1-donor-bc0 linkage is one assembly event in one transformation. Its
# PCR support counts separate PCR wells WITHIN that transformation; another transformation is a
# question of attribution, not support. A key row's trafo_ID is a SAMPLE label: a plate of a
# transformation ('trafo_2-1'), a pooled outgrowth of one transformation ('yL387'), a pool of several
# ('yL406', 'yL442'), or 'empty'. The loaders keep that label in `trafo_label` and set trafo_ID to the
# transformation it belongs to, which is what merge, aggregation and donor recovery key on.
TRAFO_LABEL_COLUMN = 'trafo_label'
# 'trafo_N-M' is plate M of transformation N: the yL406 notebook (20240702_SpG_QTL_trafo, entry of
# 20240705) washes each transformation's three plates separately and pools the washes for outgrowth.
_TRAFO_PLATE = re.compile(r'^(trafo_\d+)-\d+$')


def load_transformation_map(keyfile_dir: Path) -> Dict[str, str]:
    """yL -> the one transformation it pools, from every ``*trafo_key.tsv`` in ``keyfile_dir``.

    Each file needs ``yL`` and ``trafo_ID`` columns; identity rows (yL equal to trafo_ID) add nothing.
    A yL mapped to two different transformations is refused rather than resolved by file order.
    """
    mapping, sources = {}, {}
    for path in sorted(Path(keyfile_dir).glob('*trafo_key.tsv')):
        table = pd.read_csv(path, sep='\t', dtype=str)
        if not {'yL', 'trafo_ID'} <= set(table.columns):
            continue
        for yl, trafo in zip(table['yL'], table['trafo_ID']):
            yl, trafo = declared_value(yl), declared_value(trafo)
            if yl is None or trafo is None or yl == trafo:
                continue
            if mapping.get(yl, trafo) != trafo:
                raise ValueError(f'Trafo keys map {yl} to both {mapping[yl]} ({sources[yl]}) '
                                 f'and {trafo} ({path.name})')
            mapping[yl], sources[yl] = trafo, path.name
    return mapping


def transformation_of(label, trafo_map: Dict[str, str]) -> Optional[str]:
    """The transformation a sample label belongs to (see TRAFO_LABEL_COLUMN).

    A yL the trafo key maps to one transformation takes it; 'trafo_N-M' becomes 'trafo_N'; any other
    label (a transformation, a pool of several, 'empty') is its own scope.
    """
    value = declared_value(label)
    if value is None:
        return None
    if value in trafo_map:
        return trafo_map[value]
    plate = _TRAFO_PLATE.fullmatch(value)
    return plate.group(1) if plate else value


def _stamp_transformation(frame, row, trafo_map):
    """Keep the row's sample label and scope the frame by its transformation.

    Called after the PCR identity is built from the raw row, so a pooled outgrowth keeps its
    'pooled_yL' amplicon level.
    """
    label = declared_value(row.get('trafo_ID')) if 'trafo_ID' in row else None
    frame[TRAFO_LABEL_COLUMN] = label
    frame['trafo_ID'] = transformation_of(label, trafo_map)
    return frame


def _attach_sample_scope(frame, row, source):
    """Keep sequence evidence; validate physical scope against sample metadata."""
    from ..core.merge import SCOPE_COLUMNS, _text
    frame = frame.copy()
    population = {'yeast_library', 'yeast_sublibrary', 'trafo_ID'}
    for column in (*SCOPE_COLUMNS, 'SPS'):
        declared = _text(row.get(column))
        observed = {_text(value) for value in frame[column]} - {None} if column in frame else set()
        if declared is None:
            if column in population and observed:
                raise ValueError(f'Sample metadata lacks {column} supplied by amplicon table: {source}')
        else:
            if observed - {declared}:
                raise ValueError(f'Sample metadata conflicts with table {column}: {source}')
            frame[column] = declared
    return frame


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
        missing = set(filters) - set(df.columns)
        if missing:
            raise ValueError(f'Unknown keyfile filter column(s): {sorted(missing)}')
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
    trafo_map = load_transformation_map(keyfile_dir)

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
        combined_key = _cross_merge_keys(plate_key, sample_key)
    else:
        combined_key = plate_key

    # Derive the amplicon level explicitly rather than leaving it implied by
    # sample_number (trafo-level vs trafos-pooled-before-barcoding). A level the key
    # declares is kept; only missing levels are derived, as build_pcr_replicate_id does.
    combined_key['amplicon_level'] = _declared_or_derived_amplicon_level(combined_key)
    if 'trafo_ID' in combined_key.columns:
        lvl = combined_key['amplicon_level'].value_counts().to_dict()
        logger.info(f"  Amplicon levels (declared, else from trafo_ID): {lvl}")
    elif (combined_key['amplicon_level'] == AMPLICON_LEVEL_UNKNOWN).any():
        logger.warning(
            "  No trafo_ID column -> amplicon_level='unknown' where the key declares none. "
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
        if _KEY_CONFLICTS in combined_key.columns:
            # A conflicting row must not be de-duplicated away behind an agreeing
            # sibling that names the same file: every row for the file carries the
            # conflicts of all of them, and the file is refused when it loads.
            combined_key[_KEY_CONFLICTS] = combined_key.groupby(
                file_id_cols, dropna=False, sort=False)[_KEY_CONFLICTS].transform(
                lambda values: ';'.join(sorted({c for value in values for c in value.split(';') if c})))
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
            conflicts = row.get(_KEY_CONFLICTS)
            if conflicts:
                raise ValueError(
                    f"plate_key and sample_key declare conflicting {conflicts.split(';')} "
                    f"for {filepath.name}; correct the keyfiles before loading")
            df = pd.read_csv(filepath, sep='\t')

            df = _attach_sample_scope(df, row, filepath)

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
                        'sample_listed_polymerase'):
                df[col] = row.get(col) if pd.notna(row.get(col)) else 'NA'
            df['pcr_replicate_id'] = build_pcr_replicate_id(row)
            df = _stamp_transformation(df, row, trafo_map)
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

            # Preserve input sequence/annotation evidence and validated sample
            # scope; a whitelist previously erased donor/guide/SPS/population.
            dfs.append(df)
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


def _key_row_patterns(row, exact_suffixes, loose_pattern):
    """Exact ``{date}_{outer}_{sample}{suffix}`` names, or the loose glob only when none can be built.

    A loose ``*_{sample}_*`` glob accepts a single match from any other date or plate, so a
    row whose own file is absent would silently load another row's file under its metadata.
    It is used only when the key cannot name the file (no date, outer primers or sample).
    """
    fields = [row.get(column) if column in row else None
              for column in ('sequencing_date', 'outer_primers', 'sample_number')]
    values = [declared_value(value) for value in fields]
    if all(value is not None for value in values):
        # declared_value renders an integer-valued float date (20251125.0) as 20251125.
        stem = '_'.join(values)
        return [stem + suffix for suffix in exact_suffixes]
    if 'sample_number' in row and declared_value(row.get('sample_number')) is not None:
        return [loose_pattern.format(sample=row['sample_number'])]
    return []


def _claim_physical_file(claims, filepath, row, description):
    """Register ``filepath`` for one key row; True when it was already loaded for the same row.

    One physical file (the same device and inode, however it is named) loads once. A second
    key row naming it with different metadata is refused, because loading it again would
    count its reads twice under two sets of metadata.
    """
    status = filepath.stat()
    identity = (status.st_dev, status.st_ino)
    metadata = tuple(sorted((str(column), declared_value(value)) for column, value in row.items()))
    if identity in claims:
        if claims[identity] != metadata:
            raise ValueError(f'The combined key names one physical file under conflicting '
                             f'{description} rows: {filepath}')
        return True
    claims[identity] = metadata
    return False


def load_2x150_bc1_donor_bc0_data(
    data_dir: Path,
    keyfile_dir: Path,
    combined_key_name: str = "bc1_donor_bc0_combined_key.tsv",
    debug_mode: bool = False,
    debug_n_samples: int = 4,
    allow_ambiguous: bool = False,
    require_full_coverage: bool = True,
) -> Tuple[pd.DataFrame, DataSourceInfo]:
    """
    Load 2x150 bp bc1-donor-bc0 fragment data.

    The 2x150 bp reads capture bc1 + partial donor (donor_bc0_fragment).
    Full donor must be recovered via separate donor-bc0 amplicon.

    Each key row loads its own ``{date}_{outer}_{sample}_bc1_donor_bc0_fragment_counts.tsv`` or
    nothing; a row naming a file that is not in ``data_dir`` (for example another sequencing
    run) is logged and skipped. Every fragment table in ``data_dir`` must be named by a row.

    Args:
        data_dir: Directory containing processed TSV files
        keyfile_dir: Directory containing keyfiles
        combined_key_name: Name of combined keyfile
        debug_mode: If True, only load a few samples
        debug_n_samples: Number of samples to load in debug mode
        require_full_coverage: If True (default), raise when a fragment table in
            ``data_dir`` is named by no key row. Suppressed in debug_mode.

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
    trafo_map = load_transformation_map(keyfile_dir)

    if 'plate_polymerase' not in combined_key and 'polymerase' in combined_key:
        combined_key = combined_key.rename(columns={'polymerase': 'plate_polymerase'})
    combined_key['amplicon_level'] = _declared_or_derived_amplicon_level(combined_key)

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
    claims = {}
    unnamed_rows = []

    for _, row in tqdm(combined_key.iterrows(), total=len(combined_key),
                       desc="Loading 2x150 bc1-donor-bc0", disable=len(combined_key) < 10):
        patterns = _key_row_patterns(row, ['_bc1_donor_bc0_fragment_counts.tsv'],
                                     '*_{sample}_*fragment*.tsv')

        # Shared resolver: refuses to guess between multiple glob matches.
        # (2026-07-28 -- same defect class as the 2x300 loader; see
        # planning/QTL_analyses/BC1_LOADER_DEFECT_20260727.md.)
        filepath = _resolve_unique_file(
            data_dir, patterns,
            context=f"sample_number={row.get('sample_number')}",
            allow_ambiguous=allow_ambiguous,
        )

        if not (filepath and filepath.exists()):
            unnamed_rows.append(patterns[0] if patterns else f"sample_number={row.get('sample_number')}")
        elif not _claim_physical_file(claims, filepath, row, 'bc1-donor-bc0'):
            df = pd.read_csv(filepath, sep='\t')

            df = _attach_sample_scope(df, row, filepath)
            df['sample_number'] = row.get('sample_number', 'unknown')
            # Carry replicate-identity fields for the full PCR_plate_name
            # replicate id (gDNA_plate_name, inner_primers, outer_primers,
            # sample_number); sequencing_date for provenance only (R2 fix).
            df['outer_primers'] = row.get('outer_primers') if pd.notna(row.get('outer_primers')) else 'NA'
            df['inner_primers'] = row.get('inner_primers') if pd.notna(row.get('inner_primers')) else 'NA'
            df['gDNA_plate_name'] = row.get('gDNA_plate_name') if pd.notna(row.get('gDNA_plate_name')) else 'NA'
            df['sequencing_date'] = row.get('sequencing_date') if pd.notna(row.get('sequencing_date')) else 'NA'
            for col in ('plate_polymerase', 'condition', 'amplicon_level'):
                df[col] = row.get(col) if pd.notna(row.get(col)) else None
            df['pcr_replicate_id'] = build_pcr_replicate_id(row)
            df = _stamp_transformation(df, row, trafo_map)
            df['sequencing_run'] = '2x150'
            df['data_source'] = '2x150'

            # Extract bc0 from donor_bc0_fragment if not present
            if 'bc0' not in df.columns and 'donor_bc0_fragment' in df.columns:
                df['bc0'] = df['donor_bc0_fragment'].str[-10:]

            dfs.append(df)
            loaded_files.append(filepath)

    if unnamed_rows:
        logger.info(f"  {len(unnamed_rows)} combined-key row(s) name no file in {data_dir} "
                    f"(for example {unnamed_rows[:3]}); they are not loaded")
    _assert_full_coverage(data_dir, loaded_files, "2x150", glob_pattern="*fragment_counts.tsv",
                          require_full_coverage=require_full_coverage, debug_mode=debug_mode)

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
    from ..core.donor_recovery import aggregate_donor_candidates
    from ..core.merge import SCOPE_COLUMNS, _text
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
    trafo_map = load_transformation_map(keyfile_dir)

    # Filter to donor-bc0 amplicons
    if 'description' in combined_key.columns:
        combined_key = combined_key[combined_key['description'] == 'donor-bc0']

    if debug_mode:
        combined_key = combined_key.head(debug_n_samples)

    dfs = []
    registrations = {}

    for _, row in tqdm(combined_key.iterrows(), total=len(combined_key),
                       desc="Loading donor-bc0", disable=len(combined_key) < 10):
        # A row whose own table is absent loads nothing rather than another date's table.
        patterns = _key_row_patterns(row, ['_complete_donor_bc0_counts.tsv', '_donor_bc0_counts.tsv'],
                                     '*_{sample}_*donor_bc0*.tsv')

        # Shared resolver: refuses to guess between multiple glob matches.
        # (2026-07-28 -- same defect class as the 2x300 loader; see
        # planning/QTL_analyses/BC1_LOADER_DEFECT_20260727.md.)
        filepath = _resolve_unique_file(
            data_dir, patterns,
            context=f"sample_number={row.get('sample_number')}",
            allow_ambiguous=allow_ambiguous,
        )

        if filepath and filepath.exists():
            # One physical file however it is named: a symlink or a hardlink under a
            # second keyfile row shares the device and inode and is read once.
            status = filepath.stat()
            identity = (status.st_dev, status.st_ino)
            metadata = tuple(_text(row.get(column)) for column in ('sample_number', *SCOPE_COLUMNS, 'SPS'))
            if identity in registrations:
                if registrations[identity] != metadata:
                    raise ValueError(f'Conflicting donor-recovery registration for one physical file: {filepath}')
                continue
            registrations[identity] = metadata
            df = pd.read_csv(filepath, sep='\t')
            df['sample_number'] = row.get('sample_number', 'unknown')
            for column in (*SCOPE_COLUMNS, 'SPS'):
                value = _text(row.get(column))
                if value is not None:
                    if column in df and any(_text(item) not in (None, value) for item in df[column]):
                        raise ValueError(f'Donor-recovery metadata conflicts for {column}: {filepath}')
                    df[column] = value
            # Donor candidates pool within a transformation: a pooled outgrowth may recover
            # from its own transformation's donor-bc0 amplicon (planning 218 patch 45).
            df = _stamp_transformation(df, row, trafo_map)
            df['source_file'] = str(filepath.resolve())
            cols = ['counts', 'donor', 'donor_bc0_fragment', 'bc0', 'sample_number',
                    *SCOPE_COLUMNS, 'SPS', 'source_file']
            available_cols = [c for c in cols if c in df.columns]
            dfs.append(df[available_cols])

    if not dfs:
        logger.warning("No donor-bc0 data files loaded")
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)

    # Aggregate to create lookup table
    logger.info("Aggregating donor-bc0 data into lookup table...")
    agg_df = aggregate_donor_candidates(combined)

    logger.info(f"  Created lookup table with {len(agg_df)} scoped donor candidates; "
                "num_samples_donor_bc0 counts distinct sample_number labels, not biological replicates")

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
    from ..core.donor_recovery import recover_unique_donors
    result = recover_unique_donors(bc1_df, donor_lookup)
    recovered = int((result['donor_recovery_status'] == 'recovered_unique').sum())
    logger.info(f"  Uniquely recovered full donor for {recovered}/{len(result)} input BDB rows")
    return result


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
    gdb_files = sorted(set(gdb_dir.glob("*guide_donor_bc0*.tsv")) |
                       set(gdb_dir.glob("*gdb*.tsv")))

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
