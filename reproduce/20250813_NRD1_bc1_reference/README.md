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

## What reproduces, measured (Kevin Roy lab notes, planning/378 and 379)

* Stages b, e, f and g reproduce the surviving 2025 intermediates byte for byte (96 + 96 per-sample tables,
  4,101,755 and 2,419,499 rows).
* With `tie_order: designed` (the default) the output is identical on AMD and Intel CPUs:
  `md5(NRD1_screen_bc1_reference_table_2025.tsv) = 341b11ce2bc6c53acb0efccc6020394c`
  (pandas 2.2.2, numpy 1.26.4). Among fragments tied on read support it prefers the one the donor-bc0 or
  guide-donor-bc0 map resolves to a designed donor; untied choices are unchanged.
* The 2025-08-20 table that the screen analysis first used (`5c02f460ebb77a2f17b1865849e59456`) cannot be rebuilt
  bit for bit by any code that survives. The original scripts choose among tied rows by the order of an unstable
  sort, which differs between CPUs with and without AVX-512 (1,758 barcodes differ at stage c between the two),
  and the surviving stage c output matches neither.
* Against that table, the rebuild holds 34,714 barcodes (34,500 before). 30,491 keep the same designed oligo,
  3,914 that had none now have one, none that stays loses one, and 77 change between designs.
* Supplementary Table 1 recomputed from the rebuild (same consolidation and fitness script) has the same 23 hits,
  the same 24 edits in Figure 1D and the same depletion calls; 19 edits gain a measurement (none a hit), the
  fitness index changes by a median of 0.0004 (at most 0.12), and the barcodes used rise from 5,945 to 6,484.

`tie_order: stable` drops the designed preference; `tie_order: legacy` keeps the original sorts and equals the
original scripts run on the same CPU class.
