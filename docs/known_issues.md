# Known Issues

Defects and caveats that affect interpretation or reproducibility. Where a fix
would change the output of an already-published library, MAGESTIC deliberately
**preserves the original behavior** (and documents it here) rather than silently
altering historical designs.

## Overlapping-gene coordinate assignment bug

A flat-dictionary overwrite in the legacy haplotype-block extraction step
silently drops ~522 variants that fall in overlapping gene regions (where two
gene ±2 kb windows overlap). It affects both published arrays —
`20210422_Twist_200mer` (SpCas9 + LbCas12a) and `20240411_Twist_200mer`
(SpG). The fix is deferred to preserve exact reproducibility of the published
oligo pools.

**Full analysis, affected files, and the proposed fix:**
[Overlapping-Gene Coordinate Bug](overlapping_gene_coordinate_bug.md).

---

*As additional issues are characterized they will be listed here with status,
affected outputs, and mitigation.*
