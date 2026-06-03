# Overlapping Gene Coordinate Assignment Bug

**Status:** Known historical defect — present in all libraries generated with pre-v3.0 scripts  
**Affects:** `20210422_Twist_200mer.tsv` (SpCas9 + LbCas12a) and `20240411_Twist_200mer_oligo_array_order.tsv` (SpG)  
**Estimated impact:** ~522 variants silently dropped  
**Fix status:** Deferred to v3.1 (fixing would break exact reproducibility of published arrays)

---

## Summary

The first step of the legacy VCF-to-oligo pipeline builds a flat
`coordinate → (ORF, gene_name)` dictionary to tag each variant with the gene it
falls in.  When two QTL genes have ±2 kb flanking windows that overlap in genomic
coordinates, the dictionary only retains the *last-processed* gene for any
coordinate in the overlap zone.  Variants in that zone are subsequently assigned
to the wrong gene and, if the winning gene does not pass the FDR threshold, are
**silently dropped** from the final oligo design with no error or warning.

---

## Pipeline Context

The bug spans two legacy scripts in the 2021 and 2024 QTL `final_scripts/`
directories.  The full five-step pipeline is:

```
modify_ref_seq_with_background_variants.py
  ↓
get_haplotype_blocks_from_vcf_with_coord_offset.py   ← BUG INTRODUCED HERE
  ↓
combine_linked_variants.py
  ↓
remove_shared_variants_and_filter.py                  ← BUG TAKES EFFECT HERE
  ↓
generate_guide_donor_sequences_from_vcf.py
```

---

## Root Cause: Step 1 — Flat Dictionary Overwrites

In `get_haplotype_blocks_from_vcf_with_coord_offset.py`, a coordinate-to-gene
lookup is constructed by iterating over the QTL gene table row by row and writing
each gene's ±2 kb flanking window into a flat dict:

```python
# Lines 178–185 (2024 final_scripts version)
chrom_to_coord_to_gene_region[chrom] = chrom_to_coord_to_gene_region.get(chrom, {})
for coord in range(adjusted_start_flanking, adjusted_end_flanking):
    chrom_to_coord_to_gene_region[chrom][coord] = ORF, NAME  # last writer wins
```

For any two genes whose ±2 kb windows overlap, the gene processed **later in
the DataFrame** overwrites the earlier gene's coordinates.  After this loop,
every coordinate in the overlap zone maps to exactly one gene — whichever was
iterated last — regardless of FDR status or biological relevance.

Later in the same script, every variant is assigned to whichever gene currently
"owns" its coordinate:

```python
# Line 277
elif POS in chrom_to_coord_to_gene_region[chrom]:
    ORF, NAME = chrom_to_coord_to_gene_region[chrom][POS]
```

A variant that physically belongs to gene A (FDR-passing) but falls in the
overlap zone may be stamped with gene B's `ORF` and `NAME` instead.

---

## How Variants Are Dropped: Step 4 — FDR Filter

In `remove_shared_variants_and_filter.py`, the FDR filter (line 86) is:

```python
if NAME == 'HO' or ORF not in ORFs_passing_filter:
    skip_variant = True
```

`ORFs_passing_filter` contains only genes with `localFDR < 0.10`.  A variant
that was mis-stamped with a non-passing gene's ORF in step 1 will be skipped
here, even if it physically resides in an FDR-passing gene.  No warning or count
is emitted — the variant simply disappears from the output.

---

## Complete Failure Chain

```
Gene A (FDR-passing)  ──┐
                         ├─ overlap zone in genomic coordinates
Gene B (not FDR-passing)─┘

1. Both genes have ±2 kb flanking windows that overlap.
2. Gene B is iterated after Gene A in the QTL gene table.
3. Gene B's coordinates overwrite Gene A's in chrom_to_coord_to_gene_region.
4. Variants in the overlap zone are stamped with Gene B's ORF and NAME.
5. FDR filter: Gene B ∉ ORFs_passing_filter → skip_variant = True.
6. Those variants are absent from the VCF fed to the oligo generator.
   No warning is raised. No count of dropped variants is recorded.
```

---

## Known Examples

### MKT1 / END3 (2024 library)

MKT1 (`YNL085W`) is explicitly acknowledged in `remove_shared_variants_and_filter.py`
(line 92) as requiring "manual design" for the `MKT1 G30D` variant, but the
root cause — overlap with a neighbouring gene's window — is not named.

### Other Affected Gene Pairs

Approximately 19 gene pairs were identified as having overlapping ±2 kb windows
in the QTL gene table, totalling roughly 522 variants across the 2021 and 2024
libraries.  The exact list depends on which gene is processed last per pair (i.e.
row order in `Bloom_QTL_gene_level.tsv`), which is not deterministic across runs
if the input DataFrame is not consistently sorted.

---

## Why MAGESTIC 3.0 Does Not Fix This

`magestic.design.vcf_processing` reads the **haplotype TSV already produced by
step 1** (`load_haplotype_vcf`).  The ORF and NAME columns in that TSV are
whatever the legacy script wrote.  The MAGESTIC package does not currently
re-implement the coordinate-to-gene assignment logic, so it inherits the wrong
gene labels without knowing it.

Fixing the assignment logic in `vcf_processing.py` would cause
`magestic-design-vcf` to include ~522 additional oligos that are absent from the
historical Twist orders.  The `expected_checksums.md5` files in `reproduce/`
would then fail to verify.  Because the stated goal of v3.0 is **exact
reproduction** of those published arrays, the fix is deferred.

---

## Proposed Fix (for v3.1+)

### Step 1: Store all genes per coordinate, not just the last

Replace the overwriting assignment with a list append:

```python
# In the haplotype block extraction step (or its MAGESTIC 3.1 equivalent):
chrom_to_coord_to_gene_region[chrom].setdefault(coord, []).append((ORF, NAME))
```

### Step 2: Resolve conflicts at variant-assignment time, preferring FDR-passing genes

```python
candidates = chrom_to_coord_to_gene_region[chrom].get(POS, [])
passing = [(o, n) for o, n in candidates if o in ORFs_passing_filter]
if passing:
    ORF, NAME = passing[0]   # FDR-passing gene takes priority
elif candidates:
    ORF, NAME = candidates[0]
else:
    continue  # variant outside all gene windows
```

This ensures that a variant in an overlap zone is attributed to an FDR-passing
gene whenever one is available, preventing silent drops.

### Step 3: Emit a warning for every conflict

```python
if len(candidates) > 1:
    logger.warning(
        "Coordinate %s:%d claimed by multiple genes: %s — assigned to %s",
        chrom, POS,
        [n for _, n in candidates],
        NAME,
    )
```

### Step 4: Update expected checksums

Because v3.1 will by design produce a different (more complete) output than the
historical libraries, the `reproduce/` checksums must be updated and the release
notes must document that v3.1 designs supersede the v3.0-equivalent historical
arrays.

---

## Impact Assessment for Published Libraries

| Library | Estimated variants dropped | Approximate oligos missing |
|---|---|---|
| 20210422 SpCas9 + LbCas12a | ~300 | ~900–1,500 (3–5 oligos/variant × two nucleases) |
| 20240411 SpG | ~220 | ~660–1,100 (3–5 oligos/variant) |

These estimates assume an average of 3–5 guide-donor pairs per variant.  Variants
in high-guide-density regions may contribute more missing oligos; those with no
valid guides near the overlap zone may contribute fewer.

The missing variants do not invalidate the published screens.  The affected
variants are a small fraction of each library, and their absence is systematic
(not random), so any gene whose variants were preferentially dropped would simply
appear as a gene with fewer data points rather than a confounded data point.

---

*First documented: 2026-05-09*  
*Author: Kevin R. Roy*
