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
    manifest_samples,
    read_table,
    write_provenance,
)

logger = logging.getLogger(__name__)

STEP = "step_5_join_outcomes"
TABLE_NAME = "outcome_join.tsv"
MATRIX_NAME = "outcome_agreement_matrix.tsv"
DISAGREEMENT_NAME = "outcome_disagreements.tsv"

# The comparable edit states, plus the explicit absence-of-evidence state.
#
# PARTIAL and MIXTURE were one `mixed` state until 2026-09-23. Collapsing them made
# Stage 2's `designed_variant_partially_installed` and Stage 4's `CONTAMINATION` score
# as AGREEMENT, but they are different claims about different things: one clone whose
# donor installed incompletely is not the same observation as signal from more than one
# genotype in a well. Stage 4 says as much itself by flagging `MIXTURE_SIGNAL_NOT_PROOF`.
# They are kept apart, and the pair gets its own verdict rather than being scored as
# agreement or as contradiction (Kevin 2026-09-23: "your call for now, likely we need to
# revisit when the data comes in").
EDIT = "edit"
NO_EDIT = "no_edit"
PARTIAL = "partial_install"
MIXTURE = "mixture_signal"
UNDETERMINED = "undetermined"
STATES = (EDIT, NO_EDIT, PARTIAL, MIXTURE, UNDETERMINED)

#: The one cross-state pair that is neither agreement nor contradiction.
MIXED_STATE_MISMATCH = "mixed_state_mismatch"

#: Stage 2 `editing_outcome` (a `;spontaneous_mutations` / `;synthetic_error_installed`
#: suffix is stripped first) -> comparable state. Anything unlisted — `unknown`
#: included — is UNDETERMINED.
OUTCOME_TO_STATE = {
    "designed_variant_installed": EDIT,
    "designed_variant_partially_installed": PARTIAL,
    "designed_variant_unedited": NO_EDIT,
}

#: Stage 4 `CellReconciliation.label` -> comparable state.
LABEL_TO_STATE = {
    "SINGLE_EDIT": EDIT,
    "MULTI_EDIT": EDIT,
    "NO_EDIT": NO_EDIT,
    "CONTAMINATION": MIXTURE,
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

    manifest = manifest_samples(config)
    expected = {r.sample for r in manifest}
    precall_rows = read_table(precall_table)
    reconcile_rows = read_table(reconcile_table)
    def by_sample(rows, column, label):
        mapping = {}
        for row in rows:
            name = row.get(column, '')
            if not name or name in mapping:
                raise StepIncompleteError(f'step_5: {label} has blank or duplicate sample {name!r}')
            mapping[name] = row
        observed = set(mapping)
        if observed != expected:
            raise StepIncompleteError(f'step_5: {label} sample set differs from manifest; '
                                      f'missing={sorted(expected-observed)}, extra={sorted(observed-expected)}')
        return mapping
    precall = by_sample(precall_rows, 'sample_name', 'Stage 2')
    reconciliation = by_sample(reconcile_rows, 'sample', 'Stage 4')
    shared = [reconciliation[r.sample] for r in manifest]

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

    # The review artifact carries both true contradictions and the mixed-state pair; the
    # `agreement` column keeps them distinguishable.
    disagreements = [r for r in rows if r[-1] in ("disagree", MIXED_STATE_MISMATCH)]
    disagreement_out = disagreement_path(config)
    with atomic(disagreement_out) as tmp:
        write_table(tmp, JOIN_COLUMNS, disagreements)

    agree = sum(1 for r in rows if r[-1] == "agree")
    undetermined = sum(1 for r in rows if r[-1] == UNDETERMINED)
    mixed_mismatch = sum(1 for r in rows if r[-1] == MIXED_STATE_MISMATCH)
    contradictions = sum(1 for r in rows if r[-1] == "disagree")
    comparable = len(rows) - undetermined
    # A `low` confidence Stage 4 call is still compared, because dropping it would
    # silently shrink the denominator, but an agreement resting on one is weaker evidence
    # than the headline count suggests, so it is reported on its own line.
    agree_low_confidence = sum(1 for r in rows if r[-1] == "agree" and r[4] == "low")
    logger.info(
        "step_5: %d samples joined (Stage 2 had %d, Stage 4 had %d); %d comparable, "
        "%d agree (%s; %d resting on a low-confidence Stage 4 call), "
        "%d disagree, %d mixed-state mismatch, %d undetermined",
        len(rows),
        len(precall),
        len(reconcile_rows),
        comparable,
        agree,
        f"{100.0 * agree / comparable:.1f}% of the {comparable} comparable"
        if comparable
        else "n/a — nothing comparable",
        agree_low_confidence,
        contradictions,
        mixed_mismatch,
        undetermined,
    )

    # Every published table gets its own receipt. The join, its disagreement subset
    # and the cross-tab are separate artifacts; a reader of any one of them must be
    # able to prove which Stage 2/4 tables and which manifest (the sample scope and
    # row order) produced it.
    inputs = {"donor_outcome": precall_table, "reconcile": reconcile_table,
              "manifest": config.paths.manifest}
    counts = {
        "joined": len(rows),
        "precall_rows": len(precall),
        "reconcile_rows": len(reconcile_rows),
        "comparable": comparable,
        "agree": agree,
        "agree_low_confidence": agree_low_confidence,
        "disagree": contradictions,
        "mixed_state_mismatch": mixed_mismatch,
        "review_rows": len(disagreements),
        "undetermined": undetermined,
    }
    for artifact, artifact_rows in ((out, len(rows)), (disagreement_out, len(disagreements)),
                                    (matrix_out, len(STATES) ** 2)):
        write_provenance(
            STEP,
            artifact,
            inputs=inputs,
            parameters={
                "artifact": artifact.name,
                "agreement_rule": (
                "`unknown` (Stage 2) and `LOW_CONF` (Stage 4) are undetermined, never counted as "
                "agreement or disagreement. Partial installation against mixture signal is "
                "`mixed_state_mismatch`, counted separately and listed for review. A low-confidence "
                "Stage 4 call stays comparable; agreements resting on one are counted separately in "
                "`agree_low_confidence`."
            ),
            },
            counts={**counts, "artifact_rows": artifact_rows},
        )
    return [out, disagreement_out, matrix_out]


def _agreement(precall_state: str, reconcile_state: str) -> str:
    """agree / disagree / mixed_state_mismatch / undetermined for one sample.

    `mixed_state_mismatch` is partial installation on one side against mixture signal on
    the other. Calling that agreement would assert the two describe one phenomenon;
    calling it disagreement would assert they contradict. Neither is supported, so it is
    counted on its own and listed for review.
    """
    if UNDETERMINED in (precall_state, reconcile_state):
        return UNDETERMINED
    if precall_state == reconcile_state:
        return "agree"
    if {precall_state, reconcile_state} == {PARTIAL, MIXTURE}:
        return MIXED_STATE_MISMATCH
    return "disagree"


def outcome_states(outcomes: Sequence[str]) -> list[str]:
    """Map raw `editing_outcome` strings to comparable states (exposed for tests)."""
    return [OUTCOME_TO_STATE.get(o.split(";", 1)[0], UNDETERMINED) for o in outcomes]
