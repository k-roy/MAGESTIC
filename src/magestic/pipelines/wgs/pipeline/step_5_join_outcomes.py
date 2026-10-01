"""
Step 5 — join Stage 2's editing outcome to Stage 4's cell reconciliation.

**Orchestration-owned, and it exists because nothing in TRACE joins these two.**
Stage 2 calls the edit from variant records inside the donor window; Stage 4 calls
it from a pileup VAF over the same window. Two independent estimates of one
quantity, from different evidence, and their disagreement list is the QC artifact
worth putting in front of a human. No TRACE call: a join, a cross-tab and a count.

Cohort-level, so there is no `--sample`. Joining ACROSS reference groups on
`sample` is fine; comparing COORDINATES across them is not (chrVIII is
CUP1-contracted — 568,236 / 544,282 / 560,249 — and no single offset exists), and
that refusal lives in `core/manifest.py`. This step joins on the sample name only.

The agreement rule is declared once, at module scope, and is deliberately
THREE-valued. `unknown` (Stage 2) and `LOW_CONF` (Stage 4) are absence of
evidence, not evidence of no edit — GATE 2 means every in-window-uncovered sample
collapses to `unknown` today, so folding `unknown` into "disagrees" would
manufacture a disagreement rate out of missing gVCFs, and folding it into "agrees"
would do the reverse. Both are wrong; `undetermined` is counted separately
everywhere, and the agreement percentage states its denominator.

Author: Kevin R. Roy
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..config import WGSConfig
from ..core.emit import write_table
from ..core.sv_status import require_no_sv_rate, require_no_sv_rate_in_columns
from ._common import (
    GateError,
    StepIncompleteError,
    atomic,
    cohort_table,
    read_table,
    write_provenance,
)

logger = logging.getLogger(__name__)

STEP = "step_5_join_outcomes"
TABLE_NAME = "outcome_join.tsv"
MATRIX_NAME = "outcome_agreement_matrix.tsv"
DISAGREEMENT_NAME = "outcome_disagreements.tsv"

# The three comparable edit states, plus the explicit absence-of-evidence state.
EDIT = "edit"
NO_EDIT = "no_edit"
MIXED = "mixed"
UNDETERMINED = "undetermined"
STATES = (EDIT, NO_EDIT, MIXED, UNDETERMINED)

#: Stage 2 `editing_outcome` (a `;spontaneous_mutations` / `;synthetic_error_installed`
#: suffix is stripped first) -> comparable state. Anything unlisted — `unknown`
#: included — is UNDETERMINED.
OUTCOME_TO_STATE = {
    "designed_variant_installed": EDIT,
    "designed_variant_partially_installed": MIXED,
    "designed_variant_unedited": NO_EDIT,
}

#: Stage 4 `CellReconciliation.label` -> comparable state.
LABEL_TO_STATE = {
    "SINGLE_EDIT": EDIT,
    "MULTI_EDIT": EDIT,
    "NO_EDIT": NO_EDIT,
    "CONTAMINATION": MIXED,
    "LOW_CONF": UNDETERMINED,
}

JOIN_COLUMNS = (
    "sample",
    "editing_outcome",
    "precall_state",
    "label",
    "reconcile_confidence",
    "reconcile_state",
    "agreement",
)
MATRIX_COLUMNS = ("precall_state", "reconcile_state", "n")
AGREEMENT_INDEX = JOIN_COLUMNS.index("agreement")

# --------------------------------------------------------------------------- #
# NN 255 F7 — the mixture signal has to ARRIVE                                  #
# --------------------------------------------------------------------------- #
#
# NN 251 appended six columns to the Stage-2 table (defects 4 and 8) and named
# `minor_allele_variants` as *"the only signal that it was a mixture"*. NN 252 measured
# that **none of the six appears anywhere in this module**: the joined row was
# `[sample, outcome, p_state, label, confidence, r_state, agreement]`, so a reader of
# the cross-stage table — the table the TRIUMPH rates are read off — saw
# `designed_variant_installed` for a 10% allele with nothing beside it.
#
# The passthrough is derived FROM `step_2.OUTCOME_REASON_COLUMNS` rather than restated,
# so a seventh column cannot be added upstream and silently dropped here; and
# `test_f7_every_step_2_column_reaches_step_5` fails if one ever is.
#: Reason codes for `precall_state`, so the override below is never silent.
PRECALL_STATE_FROM_OUTCOME = "from_editing_outcome"
PRECALL_STATE_MINOR_ALLELE = "downgraded_to_mixed_minor_allele_at_designed_position"


def precall_passthrough() -> tuple[str, ...]:
    """The step-2 columns this join carries through, taken from step_2 itself."""
    from . import step_2_donor_precall

    return tuple(step_2_donor_precall.OUTCOME_REASON_COLUMNS)


def join_columns() -> tuple[str, ...]:
    """The full joined-table header, including the step-2 passthrough.

    A function rather than a constant because the step-2 import is deliberately lazy
    (this module makes no TRACE call and must stay importable without the `wgs` extra).
    """
    return JOIN_COLUMNS + ("precall_state_reason",) + precall_passthrough()

# --------------------------------------------------------------------------- #
# NN 255 F1 — a renamed upstream column must not be a KeyError                  #
# --------------------------------------------------------------------------- #
#
# NN 251 renamed `confidence` -> `bc_confidence` in `emit.RECONCILE_COLUMNS`
# (defect 15) and recorded that it had grepped `step_5` and found no consumer.
# There was exactly one, `record["confidence"]` on the row build below, and it made
# **Stage 5 unrunnable** against the table `step_4` actually writes — the stage that
# produces the Stage-2 <-> Stage-4 agreement matrix the TRIUMPH rates are read off.
#
# The remedy is deliberately a TRANSITION, not a flag-day rename: this reader accepts
# BOTH names, in preference order, and raises naming every candidate when none is
# present. A `.get(name, "")` default would have turned the KeyError into a column of
# empty strings, which is the silent-wrong-answer class this whole exercise is about.
#: `logical field -> the header names that may carry it`, newest first.
FIELD_ALIASES: dict[str, tuple[str, ...]] = {
    "reconcile_confidence": ("bc_confidence", "confidence"),
    "sample": ("sample",),
    "label": ("label",),
    "editing_outcome": ("editing_outcome",),
}


def read_field(record: dict[str, str], field: str) -> str:
    """One field of an upstream row, tolerating a rename, refusing an absence.

    Raises `StepIncompleteError` — never returns a default — when no alias is
    present: an empty string here would be published as a real value.
    """
    for name in FIELD_ALIASES.get(field, (field,)):
        if name in record:
            return record[name]
    raise StepIncompleteError(
        f"step_5: no column for {field!r} in the upstream row. Tried "
        f"{list(FIELD_ALIASES.get(field, (field,)))}; the row carries "
        f"{sorted(record)}. A renamed upstream column must be added to "
        "`step_5_join_outcomes.FIELD_ALIASES`, not defaulted to an empty value."
    )


def table_path(config: WGSConfig) -> Path:
    """The joined per-sample table (a Snakemake rule target)."""
    return cohort_table(config, TABLE_NAME)


def matrix_path(config: WGSConfig) -> Path:
    """The agreement cross-tab."""
    return cohort_table(config, MATRIX_NAME)


def disagreement_path(config: WGSConfig) -> Path:
    """The disagreements-only subset — the review artifact (a rule target)."""
    return cohort_table(config, DISAGREEMENT_NAME)


def run_step_5(config: WGSConfig) -> list[Path]:
    """Join Stages 2 and 4. Returns every path written.

    Raises `GateError` when Stage 2 produced no table, because Stage 2 is GATED by
    default (§8) — a join with one side missing is not a result.
    """
    # Imported here, not at module scope: this step makes no TRACE call, and a
    # module-scope import of its two neighbours would make `import step_5` require
    # the `magestic[wgs]` extra for nothing but two file names.
    from . import step_2_donor_precall, step_4_reconcile

    precall_table = step_2_donor_precall.table_path(config)
    reconcile_table = step_4_reconcile.table_path(config)
    if not precall_table.exists():
        raise GateError(
            f"step_5 needs Stage 2's donor-outcome table, which does not exist "
            f"({precall_table}). {step_2_donor_precall.gate_reason()}"
        )
    if not reconcile_table.exists():
        raise StepIncompleteError(
            f"step_5 needs Stage 4's reconciliation table ({reconcile_table}); run step 4."
        )

    precall = {r["sample_name"]: r for r in read_table(precall_table)}
    reconcile_rows = read_table(reconcile_table)
    shared = [r for r in reconcile_rows if r["sample"] in precall]
    if not shared:
        raise StepIncompleteError(
            f"step_5: Stage 2 ({len(precall)} samples) and Stage 4 "
            f"({len(reconcile_rows)} samples) share 0 sample names. An empty join is "
            "not 100% agreement."
        )

    passthrough = precall_passthrough()
    rows: list[list[str]] = []
    matrix: dict[tuple[str, str], int] = {}
    n_minor_downgraded = 0
    n_minor_elsewhere = 0
    for record in shared:  # Stage 4 order, which is manifest order
        sample = read_field(record, "sample")
        p_row = precall[sample]
        outcome = read_field(p_row, "editing_outcome")
        label = read_field(record, "label")
        p_state = OUTCOME_TO_STATE.get(outcome.split(";", 1)[0], UNDETERMINED)
        r_state = LABEL_TO_STATE.get(label, UNDETERMINED)
        p_reason = PRECALL_STATE_FROM_OUTCOME
        minor_at_design, minor_elsewhere = _minor_allele_split(p_row)
        if minor_elsewhere:
            n_minor_elsewhere += 1
        if minor_at_design and p_state == EDIT:
            # 🔴 A 10% allele reaches TRACE as `(chrom, pos, ref, alt)` and is bucketed
            # as a full designed variant, because `_iter_vcf_records` reads neither INFO
            # nor FORMAT. The allele FRACTION is therefore knowable only from
            # `minor_allele_variants`, and a cross-stage table that drops it publishes a
            # mixture as a clean edit. The downgrade is recorded in its own column, never
            # applied silently, and `editing_outcome` is passed through verbatim.
            #
            # Restricted to the INTERSECTION with this row's `designed_variants`: a
            # clean 100% designed edit in a cell that also carries a 30% SPONTANEOUS
            # variant elsewhere in the genome must not be downgraded, and
            # `split_plate_vcf` scans the whole plate VCF, so the column is genome-wide.
            p_state, p_reason = MIXED, PRECALL_STATE_MINOR_ALLELE
            n_minor_downgraded += 1
        rows.append(
            [
                sample,
                outcome,
                p_state,
                label,
                read_field(record, "reconcile_confidence"),
                r_state,
                _agreement(p_state, r_state),
                p_reason,
                *[p_row.get(c, "") for c in passthrough],
            ]
        )
        matrix[p_state, r_state] = matrix.get((p_state, r_state), 0) + 1
    if n_minor_downgraded or n_minor_elsewhere:
        logger.info(
            "step_5: %d sample(s) moved EDIT -> MIXED because a MAG_MINOR allele sits "
            "at a DESIGNED position (reason=%s); a further %d carry a minority allele "
            "only ELSEWHERE in the genome and were deliberately NOT downgraded",
            n_minor_downgraded, PRECALL_STATE_MINOR_ALLELE, n_minor_elsewhere,
        )

    columns = join_columns()
    # 🔴 NN 255, F2's missing half: `require_no_sv_rate`'s docstring claimed a step-5
    # call site and `grep` found none (NN 252 B4). It is wired here, over BOTH the
    # written HEADER and the provenance counts, because a column name reaches a file
    # exactly as readily as a count does.
    require_no_sv_rate_in_columns(columns, source=f"{STEP} join table header")
    require_no_sv_rate_in_columns(MATRIX_COLUMNS, source=f"{STEP} agreement matrix")
    out = table_path(config)
    with atomic(out) as tmp:
        write_table(tmp, columns, rows)

    matrix_out = matrix_path(config)
    with atomic(matrix_out) as tmp:
        write_table(
            tmp,
            MATRIX_COLUMNS,
            [[p, r, matrix.get((p, r), 0)] for p in STATES for r in STATES],
        )

    # Indexed by COLUMN POSITION, not by `r[-1]`: the step-2 passthrough appends
    # columns after `agreement`, and a negative index would silently start reading
    # `minor_allele_variants` as the agreement verdict.
    disagreements = [r for r in rows if r[AGREEMENT_INDEX] == "disagree"]
    disagreement_out = disagreement_path(config)
    with atomic(disagreement_out) as tmp:
        write_table(tmp, columns, disagreements)

    agree = sum(1 for r in rows if r[AGREEMENT_INDEX] == "agree")
    undetermined = sum(1 for r in rows if r[AGREEMENT_INDEX] == UNDETERMINED)
    comparable = len(rows) - undetermined
    logger.info(
        "step_5: %d samples joined (Stage 2 had %d, Stage 4 had %d); %d comparable, "
        "%d agree (%s), %d disagree, %d undetermined",
        len(rows),
        len(precall),
        len(reconcile_rows),
        comparable,
        agree,
        f"{100.0 * agree / comparable:.1f}% of the {comparable} comparable"
        if comparable
        else "n/a — nothing comparable",
        len(disagreements),
        undetermined,
    )

    write_provenance(
        STEP,
        out,
        inputs={"donor_outcome": precall_table, "reconcile": reconcile_table},
        parameters={
            "agreement_rule": "three-valued; `unknown` (Stage 2) and `LOW_CONF` "
            "(Stage 4) are undetermined, never counted as agreement or disagreement",
        },
        counts=_counts(
            joined=len(rows),
            precall_rows=len(precall),
            reconcile_rows=len(reconcile_rows),
            comparable=comparable,
            agree=agree,
            disagree=len(disagreements),
            undetermined=undetermined,
            minor_allele_downgraded_to_mixed=n_minor_downgraded,
            minor_allele_elsewhere_not_downgraded=n_minor_elsewhere,
        ),
    )
    return [out, disagreement_out, matrix_out]


def _counts(**kwargs: int) -> dict[str, int]:
    """Provenance counts, refused if any key would publish an SV rate (F2 / NN 252 B4)."""
    require_no_sv_rate(kwargs)
    return kwargs


def _minor_allele_split(precall_row: dict[str, str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Split `minor_allele_variants` into (at a DESIGNED position, elsewhere).

    `split_plate_vcf` scans the whole plate VCF, so the column is genome-wide: a cell
    with a clean 100% designed edit and a 30% SPONTANEOUS variant on another chromosome
    carries an id in this column, and downgrading it would be wrong. The split is by
    intersection with this row's own `designed_variants`.

    Both id strings are `{chrom}_{pos}_{ref}_to_{alt}` over RAW contig names:
    `split_plate_vcf` runs BEFORE `alias_vcf_contigs` in `step_2._run_one_sample`, and
    `_decode_outcome_table` maps TRACE's aliased ids back, so the two meet in the same
    namespace. `test_f7_the_minor_allele_id_matches_the_designed_variant_id` pins that.
    """
    minor = tuple(
        v.strip() for v in (precall_row.get("minor_allele_variants") or "").split(";")
        if v.strip()
    )
    designed = {
        v.strip() for v in (precall_row.get("designed_variants") or "").split(";")
        if v.strip()
    }
    return (
        tuple(v for v in minor if v in designed),
        tuple(v for v in minor if v not in designed),
    )


def _agreement(precall_state: str, reconcile_state: str) -> str:
    """agree / disagree / undetermined for one sample."""
    if UNDETERMINED in (precall_state, reconcile_state):
        return UNDETERMINED
    return "agree" if precall_state == reconcile_state else "disagree"


def outcome_states(outcomes: Sequence[str]) -> list[str]:
    """Map raw `editing_outcome` strings to comparable states (exposed for tests)."""
    return [OUTCOME_TO_STATE.get(o.split(";", 1)[0], UNDETERMINED) for o in outcomes]
