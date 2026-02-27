# Variant Effect Prediction (VEP) Pipeline

A unified pipeline for running multiple variant effect prediction models on yeast sequences.

## Supported Models

| Model | Input | Output | GPU Required |
|-------|-------|--------|--------------|
| **Evo2 7B** | 50kb DNA sequences | Log-probability scores | 1x H100 (80GB) |
| **Evo2 40B** | 50kb DNA sequences | Log-probability scores | 4x H100 (80GB each) |
| **Yorzoi** | 50kb DNA sequences | RNA expression predictions | 1x A100/H100 |
| **Shorkie** | 16kb DNA sequences | RNA expression (8-fold ensemble) | 1x GPU |
| **ESM1-v** | Protein sequences + variants | ΔLL scores | 1x GPU |

## Quick Start

### 1. Copy the project template

```bash
cp /path/to/common/scripts/variant_effect_prediction/templates/project_config.yaml \
   /path/to/your/project/vep_config.yaml
```

### 2. Edit the configuration

Update paths in `vep_config.yaml`:
- `project_dir`: Your project directory
- `fasta_50kb_dir`: Directory with 50kb FASTA files
- `fasta_16kb_dir`: Directory with 16kb FASTA files (for Shorkie)
- `protein_variants_file`: TSV with protein variants (for ESM1-v)

### 3. Run with Snakemake

```bash
# Dry run
snakemake -n -p --configfile vep_config.yaml

# Run with SLURM
snakemake \
    --configfile vep_config.yaml \
    --executor slurm \
    --default-resources slurm_partition=gpu slurm_account=mylab \
    --jobs 10
```

## Directory Structure

```
variant_effect_prediction/
├── __init__.py              # Package init
├── config.py                # VEPConfig dataclass
├── exceptions.py            # Custom exception hierarchy (2026-02-23)
├── pytest.ini               # Pytest configuration (2026-02-23)
├── README.md                # This file
│
├── core/                    # Shared utilities
│   └── __init__.py          # FASTA I/O, metadata, results, SLURM
│
├── sequence_utils/          # Sequence manipulation (v0.4.0)
│   ├── __init__.py          # Package exports
│   ├── basics.py            # reverse_complement, pad_sequence, gc_content
│   ├── translation.py       # CODON_TABLE, translate, codon operations
│   ├── one_hot.py           # One-hot encoding for models
│   ├── fasta_io.py          # read_fasta, write_fasta, iter_fasta
│   ├── vcf_parsing.py       # Variant class, VCF loading
│   ├── variant_ops.py       # apply_variants_to_sequence, filters
│   ├── pairwise_variants.py # VCF-guided extraction + linked detection
│   ├── alignment_variants.py # Alignment-based extraction + cached offset loading
│   └── variant_inserter.py  # Variant application with copy semantics
│
├── models/                  # Model implementations
│   ├── base.py              # VEPModel abstract base class
│   ├── evo2/
│   │   ├── __init__.py
│   │   └── runner.py        # Evo2Model class
│   ├── yorzoi/
│   │   ├── __init__.py
│   │   └── runner.py        # YorzoiModel + GPU memory management
│   ├── shorkie/
│   │   ├── __init__.py
│   │   └── runner.py        # ShorkieModel + streaming mode
│   └── esm1v/
│       ├── __init__.py
│       └── runner.py        # ESM1vModel + batch processing
│
├── snakemake/               # Snakemake rules and scripts
│   ├── Snakefile.rules      # Include-able rules
│   ├── 02_run_evo2.py
│   ├── 03_run_yorzoi.py
│   ├── 04_run_shorkie.py
│   ├── 05_run_esm1v.py
│   └── 06_merge_results.py
│
├── tests/                   # Test suite
│   ├── unit/                # Pytest unit tests (2026-02-23)
│   │   ├── __init__.py
│   │   ├── conftest.py      # Shared fixtures
│   │   ├── test_pairwise_variants.py
│   │   └── test_translation.py
│   └── ...                  # Snakemake integration test data
│
└── templates/
    └── project_config.yaml  # Example configuration
```

## Sequence Utilities (sequence_utils)

The `sequence_utils` module provides reusable functions for sequence manipulation and variant extraction.

### Variant Extraction Methods

Two approaches are available for extracting variants between strains:

```python
from sequence_utils import (
    # VCF-guided (fast, simple cases)
    extract_pairwise_variants_vcf_guided,

    # Alignment-based (handles complex overlapping indels)
    extract_pairwise_variants_alignment,
    align_and_extract_variants,

    # Auto-select based on complexity (recommended)
    extract_pairwise_variants_auto,
)
```

| Method | When to Use | Requires |
|--------|-------------|----------|
| `extract_pairwise_variants_vcf_guided` | SNP-only regions, simple cases | VCF variants for both strains |
| `extract_pairwise_variants_alignment` | Complex regions, overlapping indels | Offset table, mappy installed |
| `extract_pairwise_variants_auto` | **Recommended** - auto-selects | Both VCF variants and offset table |

**When to use alignment-based extraction:**
- Both strains have different variants at the same S288C position
- One strain has a deletion where the other has a SNP
- Complex indel regions (like ENA1 with 4 linked indels)

### Linked Variant Detection

```python
from sequence_utils import (
    detect_linked_variants,   # Group nearby variants (within 50bp)
    classify_linked_effect,   # Determine net effect (in-frame, frameshift)
    LinkedVariantGroup,       # Data class for groups
)

# Detect linked variants
groups = detect_linked_variants(variants, max_distance=50)

# Classify each group
for group in groups:
    var_type, effect, net_change = classify_linked_effect(group)
    print(f"{group.n_variants} variants, net {net_change:+d}bp: {var_type}")
```

### S288C Position Mapping

```python
from sequence_utils import (
    s288c_to_strain_pos,     # S288C → strain position
    reverse_offset_lookup,    # Strain → S288C position
    build_reverse_mapping,    # Pre-compute full mapping
    load_offset_table,        # Load JSON offset table
)
```

## Programmatic Usage

```python
from variant_effect_prediction import VEPConfig
from variant_effect_prediction.models.yorzoi import YorzoiModel, run_yorzoi_batch

# Load config
config = VEPConfig.from_yaml("vep_config.yaml")

# Run single model
with YorzoiModel() as model:
    results = model.predict_batch(sequences)

# Or use the batch function with WT caching
results = run_yorzoi_batch(
    sequences=sequences,
    metadata=metadata_df,
    wt_cache_path=Path("wt_cache.pkl")
)
```

## SLURM Resources

| Model | Partition | GPUs | Memory | Time |
|-------|-----------|------|--------|------|
| Evo2 7B | gpu | 1x H100 | 64 GB | 12 hr |
| Evo2 40B | gpu | 4x H100 | 256 GB | 48 hr |
| Yorzoi | gpu | 1x A100/H100 | 64 GB | 12 hr |
| Shorkie | gpu | 1x any | 64 GB | 12 hr |
| ESM1-v | gpu | 1x any | 32 GB | 2 hr |

## Conda Initialization

SLURM scripts use proper conda initialization for batch jobs:

```bash
eval "$(/path/to/anaconda3/bin/conda shell.bash hook)"
conda activate <env_name>
```

**Do not use** `source ~/.bashrc && conda activate` - this doesn't properly update `$PATH` in batch jobs.

## CDS Validation and Extraction

### Overview

When extracting CDS sequences from strain loci, several coordinate and biological issues can arise. The `validate_cds_sequences_v3.py` script handles:

1. **Strain coordinate calculation** using offset tables (GFF attributes are unreliable)
2. **Multi-exon gene splicing** with correct minus-strand handling
3. **Frameshift detection and fixing** by extending to next stop codon
4. **Start codon loss** detection

### Key Principle: Use Offset Tables

**Critical:** Strain GFF files may have incorrect offset values. Always compute strain coordinates from S288C coordinates + offset table:

```python
from sequence_utils import get_offset_at_position, load_offset_table

# Load offset table
offset_table = load_offset_table(Path(f"offset_tables/{strain}_{chrom}_offsets.json"))

# Compute strain position from S288C + offset
offset = get_offset_at_position(offset_table, s288c_pos)
strain_pos = s288c_pos + offset
```

### Multi-Exon Genes (Minus Strand)

For minus-strand spliced genes, reverse complement each exon individually before joining:

```python
# CORRECT for minus strand:
rc_exon_seqs = [reverse_complement(seq) for seq in exon_seqs]
cds_seq = ''.join(reversed(rc_exon_seqs))

# WRONG - don't RC the entire joined sequence:
# cds_seq = reverse_complement(''.join(exon_seqs))
```

### Frameshift Handling

When CDS length is not divisible by 3 (frameshift variant):

1. Extend downstream to find next in-frame stop codon
2. Use large extension window (2000 bp recommended)
3. Classify as lengthening (protein extended) or truncating (premature stop)

```python
def fix_frameshift(cds_seq: str, downstream_seq: str) -> str:
    """Extend CDS to next stop codon after frameshift."""
    extended = cds_seq + downstream_seq[:2000]
    for i in range(0, len(extended) - 2, 3):
        if extended[i:i+3] in {'TAA', 'TAG', 'TGA'}:
            return extended[:i+3]  # Proper CDS ending at stop
    return cds_seq  # No stop found
```

### Validation Results (Bloom et al. 127 QTL Genes × 16 Strains)

| Metric | Result |
|--------|--------|
| Total loci | 2,032 |
| Valid (original) | 94.1% |
| Valid (after fixing) | 100.0% |
| Unfixable | 1 (0.05%) |

The single unfixable case is YIR020W-A in CLIB219x, where a start codon mutation (ATG→ATT) creates a CDS with zero ATG codons anywhere - the gene is genuinely non-functional.

### Script Location

```
projects/QTL/20260210_variant_effect_prediction/scripts/validate_cds_sequences_v3.py
```

Output: `loci_50kb/metadata/cds_validation_results_v3.tsv`

## Complete Workflow: Locus Reconstruction → Variant Injection → VEP

This section documents the full data flow from the Bloom et al. VCF through variant effect prediction, with S288C coordinate back-references.

### Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│ SOURCE OF TRUTH: Bloom VCF (parents_w_svar_sorted.vcf)              │
│ - 16 strain samples, ~320,000 variants                              │
│ - ALL positions in S288C coordinates (1-based)                      │
│ - Created by aligning each strain to S288C reference                │
└─────────────────────────────────────────────────────────────────────┘
                                 │
                                 ↓
┌─────────────────────────────────────────────────────────────────────┐
│ STEP 1: Reconstruct Strain Chromosomes                              │
│ Script: build_strain_chromosomes.py                                 │
│                                                                     │
│ For each strain:                                                    │
│   1. Start with S288C reference sequence                            │
│   2. Apply VCF variants (insertions/deletions/SNPs) in order        │
│   3. Track cumulative offset at each S288C position                 │
│                                                                     │
│ OUTPUT:                                                             │
│   strain_genomes/chromosomes/{strain}_{chrom}.fasta                 │
│   strain_genomes/offset_tables/{strain}_{chrom}_offsets.json        │
└─────────────────────────────────────────────────────────────────────┘
                                 │
                                 ↓
┌─────────────────────────────────────────────────────────────────────┐
│ STEP 2: Extract 50kb/16kb Loci Windows                              │
│ Script: extract_loci.py                                             │
│                                                                     │
│ For each QTL gene (e.g., YNL085W/MKT1):                             │
│   1. Define window in S288C coordinates (gene center ± 25kb)        │
│   2. For each strain:                                               │
│      a. Look up offset at window start/end                          │
│      b. Extract corresponding region from strain chromosome         │
│   3. Store metadata with both S288C and strain coordinates          │
│                                                                     │
│ OUTPUT:                                                             │
│   loci_50kb/fasta/{ORF}_strain_loci_50kb.fasta                      │
│   loci_50kb/metadata/loci_metadata_50kb.tsv                         │
└─────────────────────────────────────────────────────────────────────┘
                                 │
                                 ↓
┌─────────────────────────────────────────────────────────────────────┐
│ STEP 3: Pairwise Variant Extraction                                 │
│ Script: pairwise_extraction.py                                      │
│                                                                     │
│ For each pair of strains (e.g., RM11 ↔ M22):                        │
│                                                                     │
│ METHOD A: VCF-Guided (fast, simple cases)                           │
│   - Load VCF variants for both strains in window                    │
│   - Compute differences: variants in L not in S, vice versa         │
│   - Map S288C positions to strain sequence positions                │
│                                                                     │
│ METHOD B: Alignment-Based (complex cases)                           │
│   - Align strain sequences directly using minimap2                  │
│   - Parse CIGAR to extract SNPs/indels                              │
│   - Use reverse offset lookup for S288C back-reference              │
│                                                                     │
│ AUTO-SELECT: extract_pairwise_variants_auto() chooses based on      │
│ whether overlapping indels exist at the same S288C position         │
│                                                                     │
│ OUTPUT: PairwiseVariant objects with:                               │
│   - bg_seq_pos (0-based position in background strain sequence)     │
│   - s288c_pos (1-based S288C coordinate for cross-strain comparison)│
│   - bg_allele, src_allele (what each strain has)                    │
│   - variant_type (SNP, INS, DEL, COMPLEX)                           │
└─────────────────────────────────────────────────────────────────────┘
                                 │
                                 ↓
┌─────────────────────────────────────────────────────────────────────┐
│ STEP 4: Generate Test Sequences (Variant Injection)                 │
│ Script: generate_test_sequences.py                                  │
│                                                                     │
│ For each PairwiseVariant:                                           │
│   1. Take background strain locus sequence                          │
│   2. Apply single variant at bg_seq_pos                             │
│   3. Retain s288c_pos for back-reference to original VCF            │
│                                                                     │
│ Test sequence naming convention:                                    │
│   {ORF}_{bg_strain}_from_{src_strain}_v{variant_index}              │
│   Example: YNL085W_RMx_from_M22_v15                                 │
│                                                                     │
│ CRITICAL: Each test sequence stores:                                │
│   - test_id (unique identifier)                                     │
│   - s288c_pos (links back to original Bloom VCF entry)              │
│   - bg_allele, src_allele (for VCF cross-reference)                 │
│                                                                     │
│ OUTPUT:                                                             │
│   pairwise_tests_50kb/fasta_50kb/{ORF}_tests_50kb.fasta             │
│   pairwise_tests_50kb/metadata/test_metadata.tsv                    │
└─────────────────────────────────────────────────────────────────────┘
                                 │
                                 ↓
┌─────────────────────────────────────────────────────────────────────┐
│ STEP 5: Run VEP Models                                              │
│                                                                     │
│ Evo2/Yorzoi: 50kb sequences → log-prob scores / expression          │
│ Shorkie: 16kb sequences → RNA expression (8-fold ensemble)          │
│ ESM1-v: Protein variants → ΔLL scores                               │
│                                                                     │
│ OUTPUT:                                                             │
│   pairwise_tests_50kb/evo2_results/scores.tsv                       │
│   pairwise_tests_50kb/yorzoi_results/predictions.tsv                │
│   pairwise_tests_50kb/shorkie_results/predictions.tsv               │
└─────────────────────────────────────────────────────────────────────┘
                                 │
                                 ↓
┌─────────────────────────────────────────────────────────────────────┐
│ STEP 6: Merge Results with VCF Back-Reference                       │
│                                                                     │
│ Final output includes:                                              │
│   - test_id, VEP scores (all models)                                │
│   - s288c_pos, chrom → original VCF lookup key                      │
│   - bg_allele, src_allele → VCF REF/ALT verification                │
│   - variant_type, codon_variant_type → annotation                   │
│                                                                     │
│ VCF BACK-REFERENCE:                                                 │
│   bcftools query parents_w_svar_sorted.vcf \                        │
│       -r chrVII:234567 -s RMx,M22 -f '%CHROM\t%POS\t%REF\t%ALT\n'   │
│                                                                     │
│ This confirms the variant in our test sequence matches the          │
│ original VCF entry, enabling integration with other Bloom analyses. │
└─────────────────────────────────────────────────────────────────────┘
```

### Offset Tables: The Key Data Structure

Offset tables track cumulative indel effects for each strain relative to S288C:

```python
# Example offset table (strain has net -3bp at position 10000)
{
    "1": 0,        # Start of chromosome: no offset
    "5432": -2,    # 2bp deletion at S288C 5432
    "8901": -3,    # 1bp deletion at S288C 8901 (cumulative: -3)
    "15000": -1,   # 2bp insertion at S288C 15000 (cumulative: -1)
}
```

**Forward mapping (S288C → strain):**
```python
def s288c_to_strain_pos(s288c_pos: int, offset_table: dict) -> int:
    """Convert S288C position to strain sequence position."""
    offset = get_offset_at_position(offset_table, s288c_pos)
    return s288c_pos + offset
```

**Reverse mapping (strain → S288C):**
```python
def reverse_offset_lookup(strain_pos: int, offset_table: dict,
                          locus_start_s288c: int, locus_end_s288c: int) -> int:
    """Convert strain position back to S288C coordinate."""
    # Build reverse mapping for the locus region
    reverse_map = {}
    for s288c_pos in range(locus_start_s288c, locus_end_s288c + 1):
        strain_p = s288c_to_strain_pos(s288c_pos, offset_table)
        if strain_p not in reverse_map:  # Keep first (handles deletions)
            reverse_map[strain_p] = s288c_pos
    return reverse_map.get(strain_pos)
```

### Why S288C Back-References Matter

1. **Cross-strain comparison**: Same S288C position = same biological location across all strains
2. **VCF integration**: Look up any variant in the original Bloom VCF by (chrom, s288c_pos)
3. **Annotation consistency**: Gene coordinates, functional annotations use S288C
4. **QTL mapping**: Bloom QTL peak positions are in S288C coordinates
5. **Literature cross-reference**: Published yeast genomics uses S288C as standard

### Handling Complex Cases

**Case 1: Overlapping indels at same S288C position**

When strain L has a deletion and strain S has a different variant at the same S288C position, VCF-guided extraction fails. Solution: use alignment-based extraction.

```python
# Auto-select handles this automatically
variants = extract_pairwise_variants_auto(
    bg_seq, src_seq,
    bg_strain, src_strain,
    bg_vcf_variants, src_vcf_variants,
    bg_offset_table,
    locus_start_s288c, locus_end_s288c
)
```

**Case 2: Multi-variant linkage (ENA1 example)**

Four indels at S288C positions 538149-538174 sum to in-frame (−10 + 1 + 8 + 1 = 0):

```python
# Detect linked variants within 50bp
groups = detect_linked_variants(variants, max_distance=50)

for group in groups:
    var_type, effect, net_change = classify_linked_effect(group)
    # group.n_variants=4, net_change=0 → "in_frame_complex"
```

**Case 3: Insertion-only positions**

When alignment reveals an insertion in the source strain, there's no corresponding S288C position (it doesn't exist in the reference). The s288c_pos is set to None, and the variant is handled specially during analysis.

### Data Provenance Trail

Every VEP prediction can be traced back to the original VCF:

```
VEP result for YNL085W_RMx_from_M22_v15
    ↓
test_metadata.tsv: s288c_pos=217234, bg_allele=G, src_allele=A
    ↓
bcftools query -r chrXIV:217234 -s M22 parents_w_svar_sorted.vcf
    ↓
Confirms: M22 has A at this position (matches src_allele)
    ↓
Annotation: D30G missense in MKT1 (known causal QTN)
```

This complete provenance enables:
- Verification of variant extraction accuracy
- Integration with experimental fitness data (MAGESTIC screens)
- Cross-referencing with published QTL studies
- Debugging of unexpected VEP predictions

## Control Sequence Generation

Control sequences serve as validation for VEP predictions. They include designed MAGESTIC edits with known effects.

### Control Categories

| Category | Description | Count | Target Genes |
|----------|-------------|-------|--------------|
| **PTC_del_syn_controls** | Premature stop codons, CDS deletions, synonymous | 1,206 | 77 QTL genes |
| **PDR5_controls** | Drug resistance saturation mutagenesis | 549 | YOR153W (PDR5) |
| **DHY214_controls** | MKT1 locus saturation mutagenesis | 58 | YNL085W + 3 others |

### MAGESTIC vs S288C Reference

**Critical:** MAGESTIC has engineered "fixes" that differ from S288C. These affect coordinate mapping:

| Locus | S288C Allele | MAGESTIC Allele | Notes |
|-------|--------------|-----------------|-------|
| MKT1 position 30 | D (GAT) | **G (GGT)** | MAGESTIC has RM "wild-type" |
| HAP1 | Frameshift (Ty1) | Full-length | Restored functional |
| ... | ... | ... | ~80 engineered changes |

This means **oligo donors designed for MAGESTIC may not match S288C coordinates**.

### DHY214 Donor-Based Alignment

For DHY214 controls, we use **donor sequence alignment** rather than coordinate-based matching:

```python
def align_donor_to_magestic(donor: str, magestic_chrom: str) -> Tuple[int, str]:
    """
    Align donor to MAGESTIC using uppercase flanks.

    Donor structure:
      - UPPERCASE: Matches MAGESTIC WT exactly (homology arms)
      - lowercase: Contains mutation + synonymous PAM disruption

    Returns:
        (donor_start_position, strand)
    """
    # Find first uppercase region
    upper_parts = re.finditer(r'[A-Z]+', donor)
    first_upper = next(upper_parts).group()[:30]

    # Search MAGESTIC for match
    if first_upper in magestic_chrom:
        pos = magestic_chrom.index(first_upper)
        # Verify all uppercase regions match
        ...
```

**Why this approach:**
1. `leftmost_codon_coord` in oligo design file may refer to S288C or design reference
2. MAGESTIC has different sequence due to engineered fixes
3. Donor uppercase regions ALWAYS match target genome (that's how HDR works)

### Applying Donor Mutations

For DHY214 controls, the entire donor sequence replaces the corresponding WT region:

```python
def apply_dhy214_edit(locus_seq: str, edit: ControlEdit) -> str:
    """Apply DHY214 edit by replacing WT with donor."""
    donor = edit.donor_seq.upper()  # Convert lowercase to uppercase

    # Position in locus window
    pos = edit.donor_start_genomic - locus.full_start_s288c

    # Replace WT region with donor (includes mutation + syn changes)
    return locus_seq[:pos] + donor + locus_seq[pos + len(donor):]
```

The donor includes:
- **Main mutation** (e.g., G30D = GGT→GAC)
- **Synonymous changes** for PAM disruption (prevents re-cutting)

### Script Location

```
projects/QTL/20260210_variant_effect_prediction/scripts/
    20260205_in_silico_QTN_dissection/bloom_trio_variants/
        33_extract_control_sequences.py   # Main extraction script
        34_run_vep_controls.py            # VEP job submission
```

### Output Structure

```
control_tests/
├── fasta_50kb/
│   ├── PTC_del_syn_controls_tests_50kb.fasta
│   ├── PDR5_controls_tests_50kb.fasta
│   ├── DHY214_controls_tests_50kb.fasta
│   └── WT_tests_50kb.fasta
├── fasta_16kb/
│   └── (same structure for Shorkie)
├── metadata/
│   └── control_test_metadata.tsv
├── evo2_results/
├── yorzoi_results/
└── shorkie_results/
```

### Validation: Expected Control Effects

| Variant Type | Expected VEP Signal |
|--------------|-------------------|
| **synonymous** | Minimal/zero effect |
| **stop_gain** (PTC) | Strong negative (loss of function) |
| **missense** | Variable (depends on position, conservation) |
| **complete_CDS_deletion** | Strong negative |

DHY214 G30D specifically is a known causal QTN affecting high-temperature growth.

## Changelog: 2026-02-23 Audit Fixes

A comprehensive audit was performed on 2026-02-23, identifying and fixing 6 critical bugs, 3 efficiency issues, and adding infrastructure improvements.

### Critical Bug Fixes

| Bug | Location | Fix |
|-----|----------|-----|
| Off-by-one boundary check | `pairwise_variants.py:249` | Position at `len(sequence)` now correctly rejected for non-empty bg_allele |
| Silent ±5bp position shift | `pairwise_variants.py:256-271` | Added `logger.warning()` when variant position is adjusted |
| Missing bounds validation | `pairwise_variants.py:73-120` | Added `validate_bounds` parameter to `compute_strain_seq_position()` |
| Variant object mutation | `variant_inserter.py:684` | Creates copies instead of mutating original PairwiseVariant objects |
| VCF filtering misses large indels | `31_pairwise_50kb_extraction.py:426` | Now includes variants that span into the region from before start |
| Multi-exon genes: only first exon | `21_validate_cds_gff.py:377` | Two-pass GFF parsing to collect all exons per ORF |

### Efficiency Improvements

| Model | Issue | Fix | Impact |
|-------|-------|-----|--------|
| **Shorkie** | Loaded all 8 models (~800MB) simultaneously | Added streaming mode in `_predict_streaming()` | ~50% memory reduction |
| **Yorzoi** | GPU memory fragmentation during tiled prediction | Added `torch.cuda.empty_cache()` every 10 tiles | Prevents OOM on long jobs |
| **ESM1-v** | Processed variants one at a time | Added `predict_variants_batch()` with batch_size=16 | ~4x throughput |

### Infrastructure Additions

#### Pytest Test Suite

```bash
# Run tests
cd common/scripts/variant_effect_prediction
pytest tests/unit/ -v

# Run with coverage
pytest tests/unit/ --cov=. --cov-report=html
```

**Test files:**
- `tests/unit/conftest.py` - Shared fixtures
- `tests/unit/test_pairwise_variants.py` - Variant application and boundary tests
- `tests/unit/test_translation.py` - DNA-to-protein translation tests

#### Custom Exception Hierarchy

New `exceptions.py` provides specific exception classes:

```python
from exceptions import (
    VEPError,              # Base exception
    SequenceError,         # Invalid/missing sequences
    VariantError,          # Invalid variants, position errors
    CoordinateError,       # Offset table, bounds errors
    TranslationError,      # Frameshift, start/stop codon errors
    ModelError,            # Model not loaded, inference errors
    FileError,             # FASTA/GFF/VCF parsing errors
)

# Usage example
try:
    result = apply_variant(sequence, variant)
except VariantPositionError as e:
    print(f"Position error: {e.message} [position={e.context['position']}]")
```

#### Offset Table Caching

Offset tables are now cached using `lru_cache` to avoid redundant file I/O:

```python
from sequence_utils import (
    load_offset_table,              # Automatically cached
    clear_offset_table_cache,       # Clear cache (for testing)
    get_offset_table_cache_info,    # Check cache stats
)

# Cache is transparent - just use load_offset_table()
offset_table = load_offset_table(Path("offsets/RM11_chrI_offsets.json"))

# Check cache performance
info = get_offset_table_cache_info()
print(f"Cache hits: {info.hits}, misses: {info.misses}")
```

### Updated Directory Structure

```
variant_effect_prediction/
├── ...
├── exceptions.py            # NEW: Custom exception hierarchy
├── pytest.ini               # NEW: Pytest configuration
└── tests/
    ├── unit/                # NEW: Unit test suite
    │   ├── __init__.py
    │   ├── conftest.py      # Shared fixtures
    │   ├── test_pairwise_variants.py
    │   └── test_translation.py
    └── ...                  # Existing Snakemake test data
```

### Migration Notes

**For existing code using `compute_strain_seq_position()`:**
```python
# Old usage (still works)
pos = compute_strain_seq_position(s288c_pos, variants, locus_start)

# New usage with validation (recommended)
pos = compute_strain_seq_position(
    s288c_pos, variants, locus_start,
    locus_end=locus_end,
    seq_length=len(sequence),
    validate_bounds=True  # Raises ValueError on invalid positions
)
```

**For existing code using `run_esm1v_batch()`:**
```python
# New batch_size parameter available (default: 16)
results = run_esm1v_batch(
    variants_df,
    output_file,
    batch_size=32,  # Increase for better throughput on large GPUs
)
```

## Author

Kevin R. Roy
Date: 2026-02-23 (comprehensive audit fixes, testing infrastructure, efficiency improvements)
