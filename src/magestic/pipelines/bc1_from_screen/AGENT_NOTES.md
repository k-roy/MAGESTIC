# Notes for bc1_from_screen_pipeline Agent

These notes were left by another agent working on the annotation pipeline.

## Relevant Work Completed

### 1. Annotation Module Created
I created an annotation module at `core/annotation.py` that:
- Parses oligo names to extract variant information (gene, AA change, codon, etc.)
- Classifies variant types (missense, nonsense, synonymous, control)
- Provides `annotate_deseq2_results()` to merge DESeq2 results with variant annotations
- Provides `create_gene_summary()` for gene-level aggregation

**Important**: The annotation module was added to `core/__init__.py` exports.

### 2. run_annotation Function Implemented
The `run_annotation()` function in `pipeline/run_pipeline.py` was implemented to:
- Load oligo design files from `common/oligo_designs/`
- Create annotation lookup from parsed oligo names
- Process DESeq2 result files and add variant annotations
- Create gene-level summaries

### 3. Oligo Name Format
The oligo names encode variant information in a specific format:

**SpG format**: `SpG_Cas9_{GENE_SYS}_{GENE_COMMON}_{CODON_POS}_{AA_REF}_{AA_ALT}_{CODON_REF}_{CODON_ALT}_{US_CHANGES}_US_codon_changes_{DS_CHANGES}_DS_codon_changes_{STRAND}_{CHR}_{GUIDE_STRAND}_{GUIDE_POS}_{PAM}_{OFFSET}_{EDIT_START}_{EDIT_END}`

**SpCas9 format**: `SpCas9_{GENE_SYS}_{GENE_COMMON}_{CODON_POS}_{AA_REF}_{AA_ALT}_{CODON_REF}_{CODON_ALT}_{US_CHANGES}_US_syn_changes_{DS_CHANGES}_DS_syn_changes_{STRAND}_{CHR}_{GUIDE_STRAND}_{GUIDE_POS}_{PAM}_{OFFSET}_{EDIT_START}_{EDIT_END}`

## Related Pipeline: variant_annotation_harmonization

There's a more comprehensive annotation pipeline at:
`/path/to/projects/QTL/20260116_variant_annotation_harmonization/`

Key script: `scripts/combine_and_annotate_designed_variants.py`

This pipeline does more advanced annotation including:
- SCORE-SV prediction assignment and imputation
- Indel length calculation and frameshift detection
- CDS overlap checking
- Sub-variant relationship identification
- Incomplete frameshift flagging

The harmonized output file contains all designed variants with comprehensive annotations.

## Sherlock Best Practices Applied

The pipelines have been updated to use:
- Environment variables (`$OAK`, `$SCRATCH`) instead of hardcoded paths
- Data staging workflow: Oak → Scratch → process → Scratch → Oak
- Conda initialization in sbatch scripts (not ~/.bashrc)

## Files Modified

1. `core/annotation.py` - NEW: Annotation module
2. `core/__init__.py` - Added annotation exports
3. `pipeline/run_pipeline.py` - Implemented `run_annotation()`

## Testing Status

- Module imports: ✓ Tested
- Oligo name parsing: ✓ Basic testing done
- Full pipeline run: NOT TESTED (needs DESeq2 results to exist)

## Recommendations

1. The annotation module provides basic functionality. For more comprehensive variant annotation, consider using the harmonized variant table from the `variant_annotation_harmonization` pipeline.

2. The current annotation parsing in `core/annotation.py` may need adjustment for edge cases (control oligos, unusual naming patterns).

3. Consider adding the harmonized variant table as a standard annotation source for DESeq2 results.

---
*Left by: Agent working on bc1_donor_bc0, guide_donor_bc0, and variant_annotation_harmonization pipelines*
*Date: 2026-02-02*
