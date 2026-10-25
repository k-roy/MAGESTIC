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
    manifest_samples,
    preflight_provenance,
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
#: planning/268 §2e (Kevin 2026-10-24). The pinned TRACE refuses a heterozygous genotype
#: in a clonal donor as `unknown;ambiguous_genotype`, which a rescued minority allele now
#: always is (`vcf_adapter.RESCUED_GT`). When that same donor-window record carries
#: `MAG_MINOR`, the well holds more than one genotype: a mixture, reported as one.
PRECALL_STATE_MINOR_ALLELE_AMBIGUOUS = "mixture_ambiguous_genotype_with_minor_allele_in_donor_window"


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
    preflight_provenance(STEP)
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

    passthrough = precall_passthrough()
    rows: list[list[str]] = []
    matrix: dict[tuple[str, str], int] = {}
    n_minor_downgraded = 0
    n_minor_ambiguous = 0
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
        ambiguous_mixture = p_state == UNDETERMINED and _minor_allele_at_ambiguous_genotype(p_row)
        if minor_elsewhere and not ambiguous_mixture:
            n_minor_elsewhere += 1
        if minor_at_design and p_state == EDIT:
            # 🔴 A 10% allele called as the clone's genotype (GT `1`) reaches TRACE as a
            # full designed variant: TRACE selects alleles by GT and never reads AD or
            # INFO for the call. The allele FRACTION is therefore knowable only from
            # `minor_allele_variants`, and a cross-stage table that drops it publishes a
            # mixture as a clean edit. The downgrade is recorded in its own column, never
            # applied silently, and `editing_outcome` is passed through verbatim.
            #
            # Restricted to the INTERSECTION with this row's `designed_variants`: a
            # clean 100% designed edit in a cell that also carries a 30% SPONTANEOUS
            # variant elsewhere in the genome must not be downgraded, and
            # `split_plate_vcf` scans the whole plate VCF, so the column is genome-wide.
            #
            # Merge note (planning/268): master downgraded to the collapsed `MIXED`, which the
            # 205 lane split into `partial_install` (Stage 2: the donor installed incompletely)
            # and `mixture_signal` (more than one genotype in the well). A minority allele
            # FRACTION at a designed base is the second, so it maps to MIXTURE.
            p_state, p_reason = MIXTURE, PRECALL_STATE_MINOR_ALLELE
            n_minor_downgraded += 1
        elif ambiguous_mixture:
            p_state, p_reason = MIXTURE, PRECALL_STATE_MINOR_ALLELE_AMBIGUOUS
            n_minor_ambiguous += 1
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
    if n_minor_downgraded or n_minor_ambiguous or n_minor_elsewhere:
        logger.info(
            "step_5: %d sample(s) moved EDIT -> mixture_signal because a MAG_MINOR allele sits "
            "at a DESIGNED position (reason=%s); %d moved undetermined -> mixture_signal "
            "because the ambiguous donor-window genotype is itself a MAG_MINOR allele "
            "(reason=%s); a further %d carry a minority allele only ELSEWHERE in the genome "
            "and were deliberately NOT downgraded",
            n_minor_downgraded, PRECALL_STATE_MINOR_ALLELE, n_minor_ambiguous,
            PRECALL_STATE_MINOR_ALLELE_AMBIGUOUS, n_minor_elsewhere,
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
    # `minor_allele_variants` as the agreement verdict. The review artifact carries both
    # true contradictions and the mixed-state pair; the `agreement` column keeps them
    # distinguishable.
    disagreements = [r for r in rows if r[AGREEMENT_INDEX] in ("disagree", MIXED_STATE_MISMATCH)]
    disagreement_out = disagreement_path(config)
    with atomic(disagreement_out) as tmp:
        write_table(tmp, columns, disagreements)

    agree = sum(1 for r in rows if r[AGREEMENT_INDEX] == "agree")
    undetermined = sum(1 for r in rows if r[AGREEMENT_INDEX] == UNDETERMINED)
    mixed_mismatch = sum(1 for r in rows if r[AGREEMENT_INDEX] == MIXED_STATE_MISMATCH)
    contradictions = sum(1 for r in rows if r[AGREEMENT_INDEX] == "disagree")
    comparable = len(rows) - undetermined
    # A `low` confidence Stage 4 call is still compared, because dropping it would
    # silently shrink the denominator, but an agreement resting on one is weaker evidence
    # than the headline count suggests, so it is reported on its own line.
    # Positional, as above: the step-2 passthrough sits after `agreement`, so `r[-1]` is
    # not the verdict on a merged row (a 205-lane line that git merged without a conflict).
    confidence_index = JOIN_COLUMNS.index("reconcile_confidence")
    agree_low_confidence = sum(
        1 for r in rows if r[AGREEMENT_INDEX] == "agree" and r[confidence_index] == "low"
    )
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
    counts = _counts(
        joined=len(rows),
        precall_rows=len(precall),
        reconcile_rows=len(reconcile_rows),
        comparable=comparable,
        agree=agree,
        agree_low_confidence=agree_low_confidence,
        disagree=contradictions,
        mixed_state_mismatch=mixed_mismatch,
        review_rows=len(disagreements),
        undetermined=undetermined,
        minor_allele_downgraded_to_mixed=n_minor_downgraded,
        minor_allele_ambiguous_to_mixed=n_minor_ambiguous,
        minor_allele_elsewhere_not_downgraded=n_minor_elsewhere,
    )
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

    Both id strings are `{chrom}_{pos}_{ref}_to_{alt}` over RAW contig names: step_2
    no longer aliases contigs (planning/268; the pinned TRACE builds its ids from
    structured fields), and reads the minor ids back from the split VCF's `MAG_MINOR`
    tags. `test_f7_the_minor_allele_id_matches_the_designed_variant_id` pins that.
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


def _minor_allele_at_ambiguous_genotype(precall_row: dict[str, str]) -> bool:
    """True when Stage 2 withheld the call for an ambiguous genotype AND a `MAG_MINOR`
    allele sits at one of those ambiguous donor-window positions (planning/268 §2e).

    Matching is by `(chrom, pos)` against TRACE's `ambiguous_genotypes` column, so a
    minority allele elsewhere in the genome never turns an ambiguous call into a
    mixture. An `unresolved_allele` entry (a symbolic or missing ALT) is not a genotype
    mixture and does not count.
    """
    import json

    outcome = precall_row.get("editing_outcome") or ""
    if not outcome.startswith("unknown") or "ambiguous_genotype" not in outcome.split(";"):
        return False
    try:
        entries = json.loads(precall_row.get("ambiguous_genotypes") or "[]")
    except ValueError:
        return False
    ambiguous = {(str(e.get("chrom")), int(e.get("pos")))
                 for e in entries if isinstance(e, dict) and e.get("reason") != "unresolved_allele"
                 and e.get("pos") is not None}
    for v in (precall_row.get("minor_allele_variants") or "").split(";"):
        v = v.strip()
        if "_to_" not in v:
            continue
        left = v.rsplit("_to_", 1)[0]
        parts = left.rsplit("_", 2)
        if len(parts) == 3 and parts[1].isdigit() and (parts[0], int(parts[1])) in ambiguous:
            return True
    return False


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
