# REDI pipeline reference (step_a → step_i)

This document mirrors the per-step prose from the legacy scripts at
`/oak/.../projects/REDI/20250312_NNS_complex_saturation_mutagenesis/scripts/`,
formalized for the `magestic.pipelines.redi` engine.

## Step A — Trim, merge, collapse, counts table

Per-sample:

1. **bbduk**: TruSeq adapter trim + quality trim. Args: `literal=AGATCGGAAGAGCGTCGTGTAGGGAAAGA ftl=20 ktrim=r k=23 mink=11 hdist=1 tpe tbo`.
2. **bbmerge**: merge overlapping paired reads.
3. **collapse_identical_merged_reads_min_counts.py**: dedupe + rename header to count.
4. **collapsed_reads_to_counts_table.py**: emit `<sample>_R1_R2_counts.tsv` with columns `[R1_seq, R2_seq, counts, sample_name]`.

## Step B — Combine keyfiles

Merges `sample_key.tsv` ← `gDNA_plate_key.tsv` ← `REDI_plate_barcode_key.tsv` into one `combined_key.tsv`. Critical detail: `REDI_plate_barcode` is coerced via `pd.to_numeric(errors='coerce')` on both sides to avoid join failures from non-integer placeholder strings.

## Step C — bc1 / REDI_bc extraction

For each row of the concatenated counts tables, finds bc1 anchored by the 12-bp prefix `CAGGTCATGCTC`; falls back to a fixed 16-bp offset slice when the prefix isn't found (KR1965/KR1967 fixed-position case). Per-row dispatch based on `inner_primers`:

- `KR1965` → bc1 from **R2**, REDI_bc from **R1**
- `KR1967` → bc1 from **R1**, REDI_bc from **R2**

The REDI_bc prefix and length come from the per-row `MAT_a_barcode_prefix` / `MAT_a_barcode_length` columns of the combined key. Output:

- `combined_REDI_bc_MAGESTIC_bc_merged_df.tsv` — grouped `(bc1, REDI_bc, sample_name) → counts` plus all metadata columns from `combined_key`.
- `MAGESTIC_bc_only_grouped_df.tsv` — `(bc1, sample_name) → counts` for samples whose `purpose == "detect all bc1 on the plate"`.

## Step D — Colony purity

For each `(sample_name, REDI_bc)` group:

1. Identify the dominant bc1 (`top_bc1`) by sorting `counts` descending.
2. Compute Levenshtein distance from every observed bc1 to `top_bc1`. Drop rows where distance is 1 or 2 (sequencing-error daughters of the dominant bc1).
3. Take the top 4 bc1's; sum their counts → `total_counts`.
4. `colony_purity = top_bc1.counts / total_counts`.

Output: `top_REDI_bc_MAGESTIC_bc1_df.tsv` (one row per `(sample_name, REDI_bc)`).

## Step E — Annotate with bc1 reference

Left-merge against the bc1-donor-bc0 reference table (output of `magestic.pipelines.bc1_donor_bc0`). Overlapping column names other than `bc1`/`counts` are prefixed `REDI_bc_bc1_` on the REDI side; `counts` → `REDI_bc_bc1_counts`. `amino_acid_position` is coerced to int with NaN → -1 sentinel.

Output: `annotated_top_REDI_bc_MAGESTIC_bc1_df.tsv`.

## Step F — Cross-array agreement

For each configured array pair (e.g., `SHA345` vs `yT182`):

1. Filter to `REDI_bc_bc1_counts >= min_counts` (default 10).
2. Build `plate_and_position` key from `MAT_alpha_array_1536_barcode_<plate>_<row>_<col>`.
3. Group by plate-and-position, sum read counts, deduplicate by position keeping the top-counts row.
4. Outer-merge A and B on plate-and-position; emit `<A>_bc1_and_<B>_bc1_both_present` and `bc1_agreement` flags.

Optional join with `1536_PIXL_rearray_format.tsv` adds 384/96-well coordinates.

Output: `<A>_<B>_merged_REDI_df.tsv`.

## Step G — WGS comparison (type-(i) ground truth)

Loads Shengdi's WGS outcome table, truncates/right-pads `bc1_from_wgs` to 20 bp with `C` fill (legacy convention), merges with step_f output on plate-and-position. Output:

- `wgs_and_<A>_bc1_agreement` — per-clone WGS-vs-array-A agreement
- `wgs_and_<B>_bc1_agreement`
- `wgs_and_<A>_<B>_agreement` — agreement when both arrays + WGS all match
- `REDI_called_bc1` — array-B bc1 fallback to array-A
- Output file: `WGS_bc1_call_outer_merge_with_REDI_df_and_plate_key.tsv`

## Step H — Plots + filtering

High-quality clone filter (legacy): `(A_counts >= 20 AND A_purity >= 0.95) OR (B_counts >= 20 AND B_purity >= 0.95)`, plus `num_independent_replicates > 1`. Plots: `bc1_clones_distribution.png`, `sublibrary_ID_counts_barchart.png`.

## Step I — Clone selection (generic primitives only)

Engine provides `select_top_n_per_group` (generic). Per-screen selection logic (NRD1/NAB3/SEN1/auxin/etc) lives in `REDI/screens/<screen>/driver.py`.
