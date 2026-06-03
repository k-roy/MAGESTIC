# Changelog

All notable changes to the MAGESTIC package are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Version scheme follows [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

### Added
- Unified `magestic` Python package combining all MAGESTIC 3.0 pipelines
- `magestic.genomics` — core genomics utilities refactored from legacy `bedgraph_computation.py` and `GTF_GFF_manipulation.py` (removed Python 2 `imp.reload` artifacts, removed hardcoded paths)
- `magestic.utils` — 16 cross-pipeline shared utility modules (migrated from `common/scripts/utils/`)
- `magestic.design` — guide-donor oligo design for VCF variants and saturation mutagenesis
- `magestic.pipelines.guide_donor_bc0` — plasmid library BC0 QC pipeline (from `guide_donor_bc0_pipeline` v1.0.0)
- `magestic.pipelines.bc1_donor_bc0` — BC1-donor-BC0 triplet linking pipeline (from `bc1_donor_bc0_pipeline` v1.1.0)
- `magestic.pipelines.bc1_from_screen` — pooled CRISPR screen analysis pipeline (from `bc1_from_screen_pipeline` v1.0.0)
- `magestic.plotting` — locus visualization: gene tracks, VEP score panels (migrated from `locus_plotting`)
- `magestic.projects.nrd1` — NRD1 alanine-scan project with bundled reference data
- `magestic.projects.vep` — Variant Effect Prediction (Evo2, Yorzoi, Shorkie, ESM1v)
- CLI entry points: `magestic-design-vcf`, `magestic-design-saturation`, `magestic-gdb-pipeline`, `magestic-bc1-link`, `magestic-screen-analysis`
- `pyproject.toml` with hatchling build system, optional dependency groups (`[align]`, `[screen]`, `[ihw]`, `[full]`)

### Changed
- `guide_donor_bc0_pipeline` → `magestic.pipelines.guide_donor_bc0` (algorithm version preserved as `__component_version__ = "1.0.0"`)
- `bc1_donor_bc0_pipeline` → `magestic.pipelines.bc1_donor_bc0` (algorithm version preserved as `__component_version__ = "1.1.0"`)
- `bc1_from_screen_pipeline` → `magestic.pipelines.bc1_from_screen` (algorithm version preserved as `__component_version__ = "1.0.0"`)
- All `sys.path.insert` hacks for shared utilities replaced with proper `magestic.utils` imports

### Fixed
- Removed Python 2 `import imp; imp.reload()` calls from genomics utilities
- Removed hardcoded absolute paths from all modules

### Known Issues
- **Overlapping gene coordinate assignment bug:** ~522 variants in overlapping gene regions are silently dropped by the legacy haplotype-block extraction step due to a flat dict overwrite when two gene ±2 kb windows overlap. Affects both published libraries (20210422 SpCas9+LbCas12a, 20240411 SpG). Fix deferred to v3.1 to preserve exact reproducibility of the published arrays. See `docs/overlapping_gene_coordinate_bug.md` for full analysis and proposed fix.

---

## [Component History — Pre-Packaging]

### guide_donor_bc0_pipeline v1.0.0
- Plasmid library BC0 QC: match guide-donor sequences to designed oligos, calculate BC0 purity
- Key constants: GUIDE_START=20, GUIDE_END=40, DONOR_START=51, DONOR_END=180
- Auto-inference of SpCas9/SpG vs LbCas12a nuclease type

### bc1_donor_bc0_pipeline v1.1.0
- BC1-donor-BC0 triplet linking: maps BC1 barcodes to designed oligos via donor sequences
- Sliding window matching (60bp), de novo structure inference, sublibrary-aware matching
- Supports 2×150 and 2×300 read data

### bc1_from_screen_pipeline v1.0.0
- Pooled CRISPR screen phenotype analysis: 18-step pipeline
- Tier classification (Tier 1: ≥100 reads → individual DESeq2, Tier 2: 20–99 reads → aggregated, Tier 3: <20 → excluded)
- Multi-timepoint fitness via linear regression, variance calibration using synonymous controls
- Optional IHW multiple testing correction (falls back to BH if R/rpy2 not available)
