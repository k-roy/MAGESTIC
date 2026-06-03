# Reproduce: Guide-Donor-BC0 Pipeline — 2021 SpCas9/LbCas12a QTL Library

## Overview

Reproduces Steps 03 (oligo matching) and 04 (BC0 purity) of the guide-donor-bc0
pipeline for library V455 (2021 SpCas9/LbCas12a QTL Twist order).

```bash
bash reproduce/20210422_SpCas9_LbCas12a_QTL_gdb/run_pipeline.sh  # from repo root
```

Runtime: ~15–20 min on a 16-core interactive node.

---

## Inputs

| File | Location | Role |
|------|----------|------|
| `V455_guide_donor_bc0_counts.tsv` | `projects/QTL/20210422_.../processed_data/guide_donor_bc0/guide_donor_bc0/` | Step 02 output (3,007,770 rows) |
| `20210422_Twist_200mer.tsv` | `common/oligo_designs/` | 48,012 designed oligos |

## Outputs

| File | Description |
|------|-------------|
| `V455_matched_guide_donor_bc0_counts.tsv` | Per-sequence oligo assignments |
| `V455_guide_donor_bc0_purity_filtered.tsv` | BC0 purity per guide-BC0 pair |
| `purity_summary.tsv` | Library-level purity statistics |

---

## Key Statistics

| Metric | Value |
|--------|-------|
| Total sequences | 3,007,770 |
| Matched | 53,786 (1.8%) |
| Perfect matches (exact guide+donor) | 18,946 (0.6%) |
| Unique guide-BC0 pairs | 2,628,120 |
| Mean BC0 purity | 0.973 |

The 1.8% match rate reflects that most collapsed FASTQ sequences are off-target
amplicons, primer dimers, or chimeras; only sequences with a recognizable
guide+SpCas9-cloning-site+donor structure are matched to the designed oligo pool.

---

## Known Discrepancy vs Legacy Pipeline Output

The MAGESTIC 3.0 `V455_matched_guide_donor_bc0_counts.tsv` differs slightly from
the legacy output in
`projects/QTL/20210422_.../processed_data/guide_donor_bc0/guide_donor_bc0/V455_matched_guide_donor_bc0_counts.tsv`.

### Concordance analysis

| Category | Rows |
|----------|------|
| Matched in both (same bc0 + guide + oligo_name) | 53,724 |
| In legacy only | 817 |
| In MAGESTIC 3.0 only | 62 |
| Agreement rate | 98.6% of legacy matched rows |

### Why they differ

| Factor | Legacy pipeline | MAGESTIC 3.0 |
|--------|----------------|--------------|
| Match priority | Unknown (legacy script not fully documented) | Exact → guide∩donor → unambiguous guide → unambiguous donor → partial donor |
| Middle donor window | Unknown | positions 30–90 within donor (60 bp) |
| Output schema | 15 columns (no `match_type`) | 16 columns (adds `match_type`) |

The 817 legacy-only and 62 MAGESTIC-3.0-only rows represent marginal cases where
the matching heuristics produce different results. The core matched set (53,724 rows)
is identical between the two implementations.

---

## Checksums

`expected_checksums.md5` reflects the MAGESTIC 3.0 pipeline output (not the legacy).
