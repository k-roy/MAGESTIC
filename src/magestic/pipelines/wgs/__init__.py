"""
magestic.pipelines.wgs — clonal-WGS orchestration over the pinned TRACE algorithms.

This module is the ORCHESTRATION/INTEGRATION layer for barcoded clonal whole-genome
sequencing (ARCHITECTURE_POLICY §11). The algorithms live in TRACE
(``trace_crispr.wgs``), which stays standalone and never imports magestic; the
dependency is one-directional and pinned by the ``magestic[wgs]`` extra.

Four analysis stages plus an orchestration-owned join. They are a DIAMOND, not a
chain: Stage 1 feeds Stage 4, while Stages 2, 3 and 4 are mutually independent.

    Stage 1  bc1 calling          de novo lineage barcode per clone
                                  TRACE: trace_crispr.wgs.bc1_consensus
    Stage 2  donor pre-calling    on-target editing outcome in the donor window
                                  TRACE: trace_crispr.wgs.donor_outcome
                                  *** GATED — see below; does not run today ***
    Stage 3  SV / plasmid         plasmid integration + donor-read cleaning
                                  TRACE: trace_crispr.wgs.sv_plasmid
    Stage 4  reconciliation       bc1 lineage set x edit VAF -> per-cell call
                                  TRACE: trace_crispr.wgs.reconcile
    Stage 5  join                 Stage 2 x Stage 4 agreement matrix (no TRACE call)

Stage 2 is GATED three ways and MUST NOT be reported as a result until all three
close (planning/166 §8):

    GATE 1  no per-oligo step-2 (synthetic) donor sequence exists for the 2021 pool,
            so the SYNTHETIC_ERROR axis is unreachable. Outputs carry an explicit
            ``synthetic_error_axis=disabled`` stamp; a zero in ``synthetic_errors``
            on an unstamped file is a wrong answer waiting to be quoted.
    GATE 2  no gVCF exists anywhere in the cohort, so ``designed_variant_unedited``
            cannot be separated from ``unknown``. ``unknown`` is never ``unedited``.
    GATE 3  TRACE's ``precall_sample`` re-parses its own variant id with ``split("_")``
            and crashes on an underscore-bearing contig — which every mandated
            composite reference carries (``plasmid_*`` / ``integrated_plasmid_*``).

Reference is a PER-SAMPLE attribute carried in the manifest, never a cohort-level
``config.genome``: one BAM directory in this program spans two references.

Author: Kevin R. Roy
"""

__all__ = ["config", "core", "readers", "pipeline", "qc"]
