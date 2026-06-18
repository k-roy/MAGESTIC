#!/usr/bin/env python3
"""One-shot generator for Phase-1 documentation stubs.

Creates a real (non-empty) placeholder for every nav entry that does not yet
have hand-written content, so `mkdocs build` renders the full navigation. Each
stub carries a one-line scope summary and a "being expanded" admonition; Phase 2
replaces these with full content. Safe to re-run — it never overwrites a page
that has already grown past the stub marker.
"""
from __future__ import annotations

from pathlib import Path

DOCS = Path(__file__).resolve().parent
STUB_MARKER = "<!-- phase1-stub -->"

# path -> (H1 title, scope sentence)
PAGES: dict[str, tuple[str, str]] = {
    # Quick start
    "quickstart/design.md": (
        "Quick Start: Library Design",
        "End-to-end walkthrough from a VCF (or a target ORF) to a synthesizable "
        "guide-donor oligo pool, with example commands and a sample config.",
    ),
    "quickstart/screen.md": (
        "Quick Start: Screen Analysis",
        "End-to-end walkthrough from barcode reads to per-variant fitness: BC0 QC, "
        "BC1 linking, and DESeq2 screen analysis.",
    ),
    # User guide
    "user_guide/overview.md": (
        "User Guide Overview",
        "Orientation to the MAGESTIC command-line tools and how the design and "
        "analysis commands chain together.",
    ),
    "user_guide/commands/pam_search.md": (
        "magestic-pam-search",
        "Enumerate all PAM sites genome-wide and write a 0/1-mismatch off-target "
        "analysis file (run once per genome/nuclease).",
    ),
    "user_guide/commands/design_vcf.md": (
        "magestic-design-vcf",
        "Design guide-donor oligos for natural variants supplied as a VCF.",
    ),
    "user_guide/commands/design_saturation.md": (
        "magestic-design-saturation",
        "Design saturation-mutagenesis guide-donor oligos covering every amino-acid "
        "substitution in a target ORF.",
    ),
    "user_guide/commands/gdb_pipeline.md": (
        "magestic-gdb-pipeline",
        "Plasmid library BC0 QC: map reads to designed oligos, compute BC0 purity, "
        "and assess library representation.",
    ),
    "user_guide/commands/bc1_link.md": (
        "magestic-bc1-link",
        "Link BC1 barcodes to guide-donor identities via shared BC0 barcodes and "
        "integrate multiple sequencing runs.",
    ),
    "user_guide/commands/screen_analysis.md": (
        "magestic-screen-analysis",
        "Pooled-screen analysis: BC1 counting, error correction, tier classification, "
        "DESeq2, and tier-aware hit calling.",
    ),
    "user_guide/commands/redi.md": (
        "magestic-redi",
        "REDI clone-array verification: colony purity, cross-array agreement, and WGS "
        "edit correspondence. Umbrella plus per-step CLIs.",
    ),
    "user_guide/commands/harmonize_oligos.md": (
        "magestic-harmonize-oligos",
        "Merge per-nuclease oligo designs into a single annotated oligo table.",
    ),
    "user_guide/commands/annotate_deseq2.md": (
        "magestic-annotate-deseq2",
        "Re-annotate DESeq2 screen results with the harmonized oligo table.",
    ),
    "user_guide/hpc_slurm.md": (
        "HPC / SLURM",
        "Running the pipelines as Snakemake workflows and SLURM job arrays on a "
        "cluster; idempotent step design for preemptible queues.",
    ),
    # Pipelines (conceptual)
    "pipelines/overview.md": (
        "Pipelines Overview",
        "The shared shape of the four analysis pipelines: config dataclass, numbered "
        "idempotent steps, and optional Snakemake workflow.",
    ),
    "pipelines/guide_donor_bc0.md": (
        "Guide-Donor-BC0 (plasmid QC)",
        "How the plasmid library is QC'd: read trimming/merging/collapse, oligo "
        "mapping over the middle-donor window, BC0 purity, and representation.",
    ),
    "pipelines/bc1_donor_bc0.md": (
        "BC1-Donor-BC0 (linking)",
        "How BC1 is linked to the guide-donor identity via the shared BC0, including "
        "donor→oligo matching, guide disambiguation, and multi-run integration.",
    ),
    "pipelines/bc1_from_screen.md": (
        "BC1 Screen Analysis",
        "How pooled-screen BC1 counts become per-variant fitness: error correction, "
        "tiering, DESeq2, slope fitness, and hit calling.",
    ),
    "pipelines/redi.md": (
        "REDI (clone arrays)",
        "How arrayed clones are verified: per-barcode counts, colony purity, "
        "cross-array agreement, and WGS correspondence.",
    ),
    # Algorithms
    "algorithms/overview.md": (
        "Algorithms Overview",
        "Index of the core algorithms behind MAGESTIC design and analysis, with "
        "pointers to each detailed page.",
    ),
    "algorithms/pam_search.md": (
        "PAM Search & Off-Target Analysis",
        "Genome-wide PAM enumeration and the 0/1-mismatch seed-region off-target "
        "model used to classify guides as unique or multi-targeting.",
    ),
    "algorithms/guide_selection.md": (
        "Guide Selection",
        "Partitioning guides into allowed (uniquely targeting) and disallowed sets, "
        "with seed-region mismatch tolerance.",
    ),
    "algorithms/donor_assembly.md": (
        "Donor Assembly",
        "Constructing the donor: homology arms, the edit, guide-disruption changes, "
        "restriction-site and poly-T avoidance, and oligo packing.",
    ),
    "algorithms/saturation_design.md": (
        "Saturation Design",
        "Generating guide-donor oligos for every amino-acid substitution across a "
        "target ORF, including subpool windowing.",
    ),
    "algorithms/barcode_collapse.md": (
        "Barcode Collapse & Error Correction",
        "Edit-distance collapsing of barcode/donor sequences and Hamming-based "
        "absorption of error barcodes into their parents.",
    ),
    "algorithms/oligo_matching.md": (
        "Donor-to-Oligo Matching",
        "Matching observed donors back to designed oligos: perfect-donor, "
        "sliding-window, and guide-aware disambiguation.",
    ),
    "algorithms/hit_calling.md": (
        "Tiered Hit Calling (DESeq2)",
        "Abundance-tiered statistics, DESeq2 normalization, multi-timepoint slope "
        "fitness, and tier-aware thresholds for calling hits.",
    ),
    "algorithms/redi_purity.md": (
        "REDI Colony Purity",
        "Scoring colony purity (top barcode / top-4 sum) and cross-array agreement "
        "for arrayed clone validation.",
    ),
    # Projects
    "projects/nrd1.md": (
        "NRD1 Alanine Scan",
        "The NRD1 alanine-scan project: bundled reference data and the saturation "
        "design specialized for NRD1.",
    ),
    "projects/vep.md": (
        "Variant Effect Prediction",
        "The VEP framework: Evo2 (DNA LM), ESM1v (protein LM), Shorkie and Yorzoi "
        "(RNA expression) runners, locus preparation, and result handling.",
    ),
    # Reference
    "coordinate_system.md": (
        "Coordinate Systems",
        "Coordinate conventions used across MAGESTIC, including MAGESTIC ↔ S288C "
        "genome coordinate conversion.",
    ),
    "troubleshooting.md": (
        "Troubleshooting",
        "Common failure modes and fixes across design and analysis.",
    ),
}


STUB_TEMPLATE = """\
# {title}

{marker}

{summary}

!!! info "Documentation in progress"
    This page is being expanded to the full MAGESTIC documentation standard
    (worked examples, parameters, and schematics). The summary above describes
    its intended scope. For the authoritative behavior today, see the
    [API Reference](../api/index.md) and the module source.
"""


def write_stub(relpath: str, title: str, summary: str) -> str:
    path = DOCS / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and STUB_MARKER not in path.read_text():
        return f"skip (has content): {relpath}"
    # adjust the relative API link depth
    api_link = "../" * relpath.count("/") + "api/index.md"
    body = STUB_TEMPLATE.format(title=title, marker=STUB_MARKER, summary=summary)
    body = body.replace("../api/index.md", api_link)
    path.write_text(body)
    return f"wrote: {relpath}"


def main() -> None:
    for relpath, (title, summary) in PAGES.items():
        print(write_stub(relpath, title, summary))


if __name__ == "__main__":
    main()
