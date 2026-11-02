# NRD1 2025 screen bc1 reference table (library yL435)

Rebuilds the table that maps every genomic barcode (bc1) of the NRD1 screen library to its designed
guide-donor, from the sequencing reads, with one package command:

```bash
magestic-bc1-reference --config config.json --jobs 16          # stages b,e,c,f,g,h
```

Inputs (fill `config.template.json`): the 96 trimmed 2x150 R1/R2 pairs of the bc1-donor-bc0 amplicon
(KR1967-KR1590, plate 5B; `*_trimmed_R1.fastq.gz`), the 96 merged and collapsed donor-bc0 reads (KR1882-KR1883,
plate 1C; `*_collapsed.fastq`, read name = count), the sample key, the GDB plasmid map of library V625, the
2024-04-11 oligo order and the NRD1 oligo annotation. Reads: BioProject PRJNA1458821.

Trimming and merging (stage a/d, bbduk + bbmerge + collapse) are upstream of this command and unchanged from
`NRD1/bc1_donor_bc0_linking/step_a_*.sh` and `step_d_*.sh`.

## What reproduces, measured (Kevin Roy lab notes, planning/378)

* Stages b, e, f and g reproduce the surviving 2025 intermediates byte for byte (96 + 96 per-sample tables,
  4,101,755 and 2,419,499 rows).
* With `tie_order: stable` (the default) the output is identical on AMD and Intel CPUs:
  `md5(NRD1_screen_bc1_reference_table_2025.tsv) = 56e0dace5d9bb1ebec0d29505be674ca`
  (pandas 2.2.2, numpy 1.26.4).
* The 2025-08-20 table that the screen analysis used (`5c02f460ebb77a2f17b1865849e59456`) cannot be rebuilt bit
  for bit by any code that survives. The original scripts choose among tied rows by the order of an unstable
  sort, which differs between CPUs with and without AVX-512 (1,758 barcodes differ at stage c between the two),
  and the surviving stage c output matches neither. Measured against it, the stable rebuild has the same
  34,500 barcodes and the same 72 columns; 32,455 barcodes (94.1%) are identical in every column, 2,045 differ,
  and 725 of those are assigned a different designed oligo. In the cases inspected the difference is which of
  several equally supported fragments, often one or two reads each, represents the barcode.

`tie_order: legacy` keeps the original sorts and equals the original scripts run on the same CPU class.
