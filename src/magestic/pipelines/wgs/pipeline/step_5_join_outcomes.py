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

    rows: list[list[str]] = []
    matrix: dict[tuple[str, str], int] = {}
    for record in shared:  # Stage 4 order, which is manifest order
        sample = record["sample"]
        outcome = precall[sample]["editing_outcome"]
        label = record["label"]
        p_state = OUTCOME_TO_STATE.get(outcome.split(";", 1)[0], UNDETERMINED)
        r_state = LABEL_TO_STATE.get(label, UNDETERMINED)
        rows.append(
            [
                sample,
                outcome,
                p_state,
                label,
                record["confidence"],
                r_state,
                _agreement(p_state, r_state),
            ]
        )
        matrix[p_state, r_state] = matrix.get((p_state, r_state), 0) + 1

    out = table_path(config)
    with atomic(out) as tmp:
        write_table(tmp, JOIN_COLUMNS, rows)

    matrix_out = matrix_path(config)
    with atomic(matrix_out) as tmp:
        write_table(
            tmp,
            MATRIX_COLUMNS,
            [[p, r, matrix.get((p, r), 0)] for p in STATES for r in STATES],
        )

    disagreements = [r for r in rows if r[-1] == "disagree"]
    disagreement_out = disagreement_path(config)
    with atomic(disagreement_out) as tmp:
        write_table(tmp, JOIN_COLUMNS, disagreements)

    agree = sum(1 for r in rows if r[-1] == "agree")
    undetermined = sum(1 for r in rows if r[-1] == UNDETERMINED)
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
        counts={
            "joined": len(rows),
            "precall_rows": len(precall),
            "reconcile_rows": len(reconcile_rows),
            "comparable": comparable,
            "agree": agree,
            "disagree": len(disagreements),
            "undetermined": undetermined,
        },
    )
    return [out, disagreement_out, matrix_out]


def _agreement(precall_state: str, reconcile_state: str) -> str:
    """agree / disagree / undetermined for one sample."""
    if UNDETERMINED in (precall_state, reconcile_state):
        return UNDETERMINED
    return "agree" if precall_state == reconcile_state else "disagree"


def outcome_states(outcomes: Sequence[str]) -> list[str]:
    """Map raw `editing_outcome` strings to comparable states (exposed for tests)."""
    return [OUTCOME_TO_STATE.get(o.split(";", 1)[0], UNDETERMINED) for o in outcomes]
