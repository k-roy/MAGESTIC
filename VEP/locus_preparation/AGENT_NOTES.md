# Locus Preparation Module - Agent Notes

**Date:** 2026-02-10
**Updated:** 2026-02-22 (alignment-based variant extraction)
**Purpose:** Build genomic loci and insert variants for VEP

---

## Module Responsibility

This module handles:
1. **Building genomic loci** from reference + VCF
2. **Inserting variants** into strain backgrounds
3. **Coordinate mapping** between strain assemblies
4. **Variant annotation** (codon changes, consequences)
5. **Validation** via GFF parsing and bcftools

---

## CRITICAL: Variant Extraction Methods (Updated 2026-02-22)

### Problem with VCF-Guided Approach

The original VCF-guided variant extraction (`extract_pairwise_variants_vcf_guided`) **fails when both strains have different variants at the same S288C position**:

```
Example failure:
- S288C position 1000: Reference = "ATGCG"
- Strain L: deletion "ATGCG" → "A" (offset -4)
- Strain S: SNP "ATGCG" → "TTGCG"
- Problem: Position 1000 doesn't exist in L's sequence!
```

### Solution: Use `sequence_utils` Module

The `sequence_utils` module (in `common/scripts/variant_effect_prediction/sequence_utils/`) provides three extraction methods:

| Method | Use Case | Function |
|--------|----------|----------|
| **VCF-guided** | Simple SNPs, no overlapping indels | `extract_pairwise_variants_vcf_guided()` |
| **Alignment-based** | Complex regions, overlapping indels | `extract_pairwise_variants_alignment()` |
| **Auto** | Recommended - auto-selects | `extract_pairwise_variants_auto()` |

### When to Use Alignment-Based Extraction

Use alignment-based extraction when:
- Both strains have different variants at the same S288C position
- One strain has a deletion where another has a SNP
- Complex indel regions (e.g., ENA1 with 4 linked indels)
- Uncertain about variant complexity

### Example Usage

```python
from sequence_utils import (
    extract_pairwise_variants_auto,
    extract_pairwise_variants_alignment,
    detect_linked_variants,
    classify_linked_effect,
    load_offset_table,
)

# Load offset table for background strain
bg_offsets = load_offset_table(Path("strain_genomes/offset_tables/BYa_chrIV_offsets.json"))

# Auto-select extraction method (recommended)
variants = extract_pairwise_variants_auto(
    bg_seq=strain_L_seq,
    src_seq=strain_S_seq,
    background_strain="StrainL",
    source_strain="StrainS",
    bg_vcf_variants=L_variants,
    src_vcf_variants=S_variants,
    bg_offset_table=bg_offsets,
    s288c_start=locus_start,
    s288c_end=locus_end
)

# Detect and classify linked variants
groups = detect_linked_variants(variants, max_distance=50)
for group in groups:
    var_type, effect, net_change = classify_linked_effect(group)
    print(f"{group.n_variants} variants, net {net_change:+d}bp: {var_type}")
```

### How Alignment-Based Extraction Works

1. **Direct alignment** using minimap2 (via `mappy` Python bindings)
2. **CIGAR parsing** to extract SNPs, insertions, deletions
3. **Reverse offset lookup** to map strain positions back to S288C coordinates
4. **Collision handling** for deleted positions (keeps earliest S288C position)

### Files in `sequence_utils/`

- `alignment_variants.py` - Alignment-based extraction (NEW 2026-02-22)
- `pairwise_variants.py` - VCF-guided extraction + linked detection
- `vcf_parsing.py` - VCF loading and Variant class
- `variant_ops.py` - Apply variants to sequences

---

## Files to Create

```
locus_preparation/
├── __init__.py           # ✅ Created
├── AGENT_NOTES.md        # This file
├── locus_builder.py      # Build genomic loci
├── variant_insertion.py  # Insert variants into backgrounds
├── coordinate_mapping.py # Handle offsets
├── variant_annotation.py # Annotate variants
└── validation.py         # GFF/bcftools validation
```

---

## Source Scripts to Extract From

### Primary Source: `12_build_all_bloom_loci_s288c.py`

**Location:** `/path/to/projects/QTL/variant_effect_prediction/scripts/20260205_in_silico_QTN_dissection/12_build_all_bloom_loci_s288c.py`

**Key Functions:**
- `build_strain_sequence()` → `locus_builder.py`
- `apply_variants_to_sequence()` → `variant_insertion.py`
- `calculate_offset()` → `coordinate_mapping.py`

### Secondary Source: `18_generate_trio_dissection_tests.py`

**Location:** Same directory

**Key Functions:**
- `get_unique_variants()` → `variant_insertion.py`
- `insert_variant()` → `variant_insertion.py`
- Round-robin strain logic → Keep in project scripts

---

## Key Data Structures

### Locus Metadata (from `all_loci_metadata_reclassified.tsv`)

| Column | Description |
|--------|-------------|
| orf | Gene name (e.g., YOR128C) |
| strain | Strain name (e.g., RMx, BYa) |
| chrom | Chromosome |
| start | S288C start coordinate |
| end | S288C end coordinate |
| strand | + or - |
| offset | Coordinate offset for this strain |
| classification | VALID_REFERENCE, VALID_EXTENDED, etc. |

### Variant Format (from VCF)

```python
variant = {
    'chrom': 'chrXIV',
    'pos': 461485,      # 1-based
    'ref': 'A',
    'alt': 'G',
    'strain': 'CBS2888a'
}
```

---

## Function Specifications

### `locus_builder.py`

```python
def build_locus(
    reference_fasta: Path,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    vcf_file: Path = None,
    strain: str = None,
    flank: int = 0
) -> str:
    """
    Build a genomic locus sequence.

    Args:
        reference_fasta: Path to S288C reference FASTA
        chrom: Chromosome name
        start: Start coordinate (1-based)
        end: End coordinate (1-based)
        strand: + or -
        vcf_file: Optional VCF with strain variants
        strain: Strain name for VCF extraction
        flank: Additional flanking sequence on each side

    Returns:
        DNA sequence string (uppercase)
    """
    pass


def build_loci_batch(
    metadata_df: pd.DataFrame,
    reference_fasta: Path,
    vcf_file: Path,
    output_dir: Path,
    seq_length: int = 15000
) -> pd.DataFrame:
    """Build multiple loci from metadata table."""
    pass
```

### `variant_insertion.py`

```python
def insert_variant(
    sequence: str,
    position: int,  # 0-based within sequence
    ref: str,
    alt: str
) -> str:
    """
    Insert a single variant into a sequence.

    Validates ref matches sequence at position.
    Handles SNPs and indels.
    """
    pass


def get_unique_variants(
    strain_a: str,
    strain_b: str,
    vcf_file: Path,
    region: str = None  # e.g., "chrXIV:460000-465000"
) -> List[dict]:
    """
    Get variants in strain_a that are NOT in strain_b.

    Used for trio-based dissection.
    """
    pass
```

### `coordinate_mapping.py`

```python
def calculate_offset(
    strain: str,
    chrom: str,
    vcf_file: Path
) -> int:
    """
    Calculate coordinate offset for a strain relative to S288C.

    Offset accumulates from indels upstream of the region.
    """
    pass


def s288c_to_strain_coord(
    s288c_pos: int,
    strain: str,
    chrom: str,
    offset_map: dict
) -> int:
    """Convert S288C coordinate to strain-specific coordinate."""
    pass
```

---

## Validation Requirements

### bcftools Consensus Check

```bash
# Build sequence using bcftools (gold standard)
bcftools consensus -f S288C.fa -s STRAIN -r chrN:start-end bloom.vcf.gz > expected.fa

# Compare to our rebuilt sequence
diff expected.fa our_rebuilt.fa
```

### GFF Parsing for CDS Boundaries

```python
def parse_gff_cds(gff_file: Path, orf: str) -> dict:
    """
    Extract CDS coordinates from GFF.

    Returns:
        {'chrom': str, 'start': int, 'end': int, 'strand': str}
    """
    pass
```

---

## Dependencies

- `pysam`: FASTA/VCF reading
- `pandas`: Metadata handling
- `subprocess`: bcftools calls (validation only)

**Conda environment:** `pysam` or `snakemake`

---

## Testing

```python
if __name__ == "__main__":
    from pathlib import Path

    BASE = Path("/path/to")
    REF = BASE / "common/reference_genomes/S288C_reference.fasta"
    VCF = BASE / "common/published_data/bloom_et_al_2019/bloom.vcf.gz"

    # Test building one locus
    seq = build_locus(REF, "chrXIV", 460000, 475000, "+", VCF, "RMx")
    print(f"Built locus: {len(seq)} bp")

    # Test variant insertion
    test_seq = "ATCGATCGATCG"
    modified = insert_variant(test_seq, 3, "G", "A")
    assert modified == "ATCAATCGATCG"
    print("Variant insertion test passed")
```

---

*Last updated: 2026-02-22 (added alignment-based variant extraction)*
