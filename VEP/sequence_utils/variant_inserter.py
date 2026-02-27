"""
Unified variant insertion and protein analysis module.

This module provides a standardized interface for inserting variants from one
strain into another and analyzing the resulting protein effects. It handles:
- Complete variants (full linked haplotype units)
- Sub-variants (individual variants from linked groups)
- Combined AA variants (nearby AA changes tested together)
- Cross-reference position mapping via protein alignment

Classes:
- VariantInserter: Main class for variant insertion and analysis
- VariantAnalysisResult: Container for analysis results

Author: Kevin R. Roy
Date: 2026-02-21
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set, Union
import re

# Handle both module import and direct execution
try:
    from .pairwise_variants import (
        PairwiseVariant,
        LinkedVariantGroup,
        detect_linked_variants,
        classify_linked_effect,
        extract_pairwise_variants_vcf_guided,
        apply_pairwise_variant,
        apply_multiple_pairwise_variants,
    )
    from .translation import (
        translate,
        translate_with_stops,
        AAVariant,
        CombinedAAVariant,
        find_aa_differences,
        combine_nearby_aa_changes,
        classify_frameshift_direction,
        classify_start_loss,
    )
    from .vcf_parsing import (
        Variant,
        load_variants_for_region,
        get_strain_variants,
    )
    from .fasta_io import read_fasta, write_fasta
except ImportError:
    # Direct execution
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from pairwise_variants import (
        PairwiseVariant,
        LinkedVariantGroup,
        detect_linked_variants,
        classify_linked_effect,
        extract_pairwise_variants_vcf_guided,
        apply_pairwise_variant,
        apply_multiple_pairwise_variants,
    )
    from translation import (
        translate,
        translate_with_stops,
        AAVariant,
        CombinedAAVariant,
        find_aa_differences,
        combine_nearby_aa_changes,
        classify_frameshift_direction,
        classify_start_loss,
    )
    from vcf_parsing import (
        Variant,
        load_variants_for_region,
        get_strain_variants,
    )
    from fasta_io import read_fasta, write_fasta


# ============================================================================
# Position Mapping for Cross-Strain Comparison
# ============================================================================

@dataclass
class PositionMapping:
    """Mapping between strain and reference protein positions.

    Attributes:
        strain_pos: Position in strain protein (1-based)
        ref_pos: Corresponding position in reference protein (1-based)
        strain_aa: Amino acid at strain position
        ref_aa: Amino acid at reference position (may differ due to variant)
        alignment_quality: Local alignment quality score (0-1)
    """
    strain_pos: int
    ref_pos: int
    strain_aa: str
    ref_aa: str
    alignment_quality: float = 1.0


def build_position_mapping(
    strain_protein: str,
    ref_protein: str,
) -> Dict[int, PositionMapping]:
    """Build position mapping from strain protein to reference protein.

    Uses simple gapless alignment for proteins of similar length, or
    pairwise alignment for proteins with indels.

    Args:
        strain_protein: Strain protein sequence
        ref_protein: Reference (S288C) protein sequence

    Returns:
        Dict mapping strain positions (1-based) to PositionMapping objects
    """
    mapping = {}

    # For proteins of same length, use direct mapping
    if len(strain_protein) == len(ref_protein):
        for i in range(len(strain_protein)):
            mapping[i + 1] = PositionMapping(
                strain_pos=i + 1,
                ref_pos=i + 1,
                strain_aa=strain_protein[i],
                ref_aa=ref_protein[i],
                alignment_quality=1.0
            )
        return mapping

    # For proteins of different length, use Needleman-Wunsch-like alignment
    # Simple implementation: find longest common subsequence regions
    return _align_proteins_simple(strain_protein, ref_protein)


def _align_proteins_simple(
    strain_protein: str,
    ref_protein: str,
) -> Dict[int, PositionMapping]:
    """Simple protein alignment for position mapping.

    Uses a gap-aware alignment approach. For each position in the strain
    protein, finds the corresponding reference position accounting for
    insertions/deletions.

    This is a simplified implementation that assumes:
    - Proteins are mostly similar (>80% identity)
    - Indels are relatively small (<10 AA)
    """
    mapping = {}

    # Use dynamic programming alignment
    # Match score = 2, Mismatch = -1, Gap = -2
    m, n = len(strain_protein), len(ref_protein)

    # Score matrix
    score = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[None] * (n + 1) for _ in range(m + 1)]

    # Initialize
    for i in range(1, m + 1):
        score[i][0] = -2 * i
        traceback[i][0] = 'U'  # Up
    for j in range(1, n + 1):
        score[0][j] = -2 * j
        traceback[0][j] = 'L'  # Left

    # Fill matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = 2 if strain_protein[i-1] == ref_protein[j-1] else -1
            diag = score[i-1][j-1] + match
            up = score[i-1][j] - 2
            left = score[i][j-1] - 2

            if diag >= up and diag >= left:
                score[i][j] = diag
                traceback[i][j] = 'D'
            elif up >= left:
                score[i][j] = up
                traceback[i][j] = 'U'
            else:
                score[i][j] = left
                traceback[i][j] = 'L'

    # Traceback to build mapping
    i, j = m, n
    alignment_pairs = []

    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == 'D':
            alignment_pairs.append((i, j))
            i -= 1
            j -= 1
        elif i > 0 and traceback[i][j] == 'U':
            alignment_pairs.append((i, None))  # Gap in reference
            i -= 1
        else:
            alignment_pairs.append((None, j))  # Gap in strain
            j -= 1

    alignment_pairs.reverse()

    # Build mapping from alignment
    for strain_pos, ref_pos in alignment_pairs:
        if strain_pos is not None:
            if ref_pos is not None:
                mapping[strain_pos] = PositionMapping(
                    strain_pos=strain_pos,
                    ref_pos=ref_pos,
                    strain_aa=strain_protein[strain_pos - 1],
                    ref_aa=ref_protein[ref_pos - 1],
                    alignment_quality=1.0 if strain_protein[strain_pos - 1] == ref_protein[ref_pos - 1] else 0.5
                )
            else:
                # Strain has insertion - no corresponding ref position
                # Map to nearest ref position
                nearest_ref = _find_nearest_ref_pos(alignment_pairs, strain_pos)
                mapping[strain_pos] = PositionMapping(
                    strain_pos=strain_pos,
                    ref_pos=nearest_ref,
                    strain_aa=strain_protein[strain_pos - 1],
                    ref_aa='-',
                    alignment_quality=0.0
                )

    return mapping


def _find_nearest_ref_pos(alignment_pairs: List[Tuple], strain_pos: int) -> int:
    """Find nearest reference position for a strain position with no direct mapping."""
    # Look backward for nearest mapped reference position
    for sp, rp in reversed(alignment_pairs):
        if sp is not None and sp < strain_pos and rp is not None:
            # Return the next reference position
            return rp + (strain_pos - sp)
    # Fallback: return strain position
    return strain_pos


def map_strain_position_to_reference(
    strain_protein: str,
    ref_protein: str,
    strain_position: int
) -> Tuple[int, float]:
    """Map a single strain AA position to reference position.

    Convenience function that builds full alignment and returns the
    mapping for a specific position.

    Args:
        strain_protein: Strain protein sequence
        ref_protein: Reference (S288C) protein sequence
        strain_position: Position in strain protein (1-based)

    Returns:
        Tuple of (reference_position, alignment_quality)
        - reference_position: Corresponding position in reference (1-based)
        - alignment_quality: Score 0-1 indicating mapping confidence
    """
    mapping = build_position_mapping(strain_protein, ref_protein)

    if strain_position in mapping:
        pm = mapping[strain_position]
        return pm.ref_pos, pm.alignment_quality
    else:
        # Position outside protein - return with low quality
        return strain_position, 0.0


# ============================================================================
# Variant Test Record Classes
# ============================================================================

@dataclass
class VariantTestRecord:
    """Record for a single variant test.

    Attributes:
        test_id: Unique identifier for this test
        gene_name: Gene/ORF name
        background_strain: Strain receiving the variant
        source_strain: Strain providing the variant
        variant_type: single, complete, or sub_variant
        dna_variants: List of DNA-level PairwiseVariants
        linked_group_id: ID if part of a linked group (None for single)
        is_complete_variant: True if testing all variants in linked group
        background_cds: Background CDS sequence
        mutant_cds: Mutant CDS sequence after variant insertion
        background_protein: Background protein sequence
        mutant_protein: Mutant protein sequence
        aa_variants: List of AA-level changes
        combined_aa_variants: List of combined nearby AA changes
        variant_classification: missense, nonsense, frameshift, etc.
        frameshift_direction: truncating or lengthening (if applicable)
        start_loss_consequence: inframe or outframe (if applicable)
    """
    test_id: str
    gene_name: str
    background_strain: str
    source_strain: str
    variant_type: str  # 'single', 'complete', 'sub_variant'
    dna_variants: List[PairwiseVariant]
    linked_group_id: Optional[str] = None
    is_complete_variant: bool = False
    background_cds: str = ""
    mutant_cds: str = ""
    background_protein: str = ""
    mutant_protein: str = ""
    aa_variants: List[AAVariant] = field(default_factory=list)
    combined_aa_variants: List[CombinedAAVariant] = field(default_factory=list)
    variant_classification: str = ""
    frameshift_direction: Optional[str] = None
    start_loss_consequence: Optional[str] = None
    ref_position_mapping: Dict[int, int] = field(default_factory=dict)  # strain pos -> ref pos

    @property
    def net_dna_length_change(self) -> int:
        """Net DNA length change from all variants."""
        return sum(v.length_change for v in self.dna_variants)

    @property
    def protein_length_change(self) -> int:
        """Protein length change."""
        return len(self.mutant_protein) - len(self.background_protein)

    @property
    def s288c_positions(self) -> List[int]:
        """List of S288C DNA positions affected."""
        return sorted(v.s288c_pos for v in self.dna_variants)

    @property
    def aa_positions(self) -> List[int]:
        """List of AA positions changed (1-based)."""
        return sorted(v.position for v in self.aa_variants)


@dataclass
class VariantAnalysisResult:
    """Container for complete analysis results.

    Attributes:
        gene_name: Gene/ORF being analyzed
        background_strain: Strain receiving variants
        source_strain: Strain providing variants
        single_variant_tests: Individual variant tests
        complete_variant_tests: Tests with all linked variants together
        sub_variant_tests: Tests with subsets of linked variants
        linked_groups: Linked variant groups detected
        position_mapping: Protein position mapping to reference
    """
    gene_name: str
    background_strain: str
    source_strain: str
    single_variant_tests: List[VariantTestRecord] = field(default_factory=list)
    complete_variant_tests: List[VariantTestRecord] = field(default_factory=list)
    sub_variant_tests: List[VariantTestRecord] = field(default_factory=list)
    linked_groups: List[LinkedVariantGroup] = field(default_factory=list)
    position_mapping: Dict[int, PositionMapping] = field(default_factory=dict)

    @property
    def all_tests(self) -> List[VariantTestRecord]:
        """All test records combined."""
        return self.single_variant_tests + self.complete_variant_tests + self.sub_variant_tests

    @property
    def n_total_tests(self) -> int:
        """Total number of tests."""
        return len(self.all_tests)

    @property
    def n_missense(self) -> int:
        """Count of missense variants."""
        return sum(1 for t in self.all_tests if t.variant_classification == 'missense')

    @property
    def n_synonymous(self) -> int:
        """Count of synonymous variants."""
        return sum(1 for t in self.all_tests if t.variant_classification == 'synonymous')

    @property
    def n_frameshift(self) -> int:
        """Count of frameshift variants."""
        return sum(1 for t in self.all_tests if 'frameshift' in t.variant_classification)


# ============================================================================
# Main VariantInserter Class
# ============================================================================

class VariantInserter:
    """Unified interface for inserting variants and analyzing protein effects.

    This class provides a standardized procedure for:
    - Extracting variants between two strains (from VCF or sequence comparison)
    - Identifying linked/complete variants (nearby indels that co-occur)
    - Inserting variants into background sequences
    - Translating to protein and finding AA differences
    - Combining nearby AA changes for ESM1-v scoring
    - Mapping positions to a reference (S288C) for cross-strain comparison

    Usage:
        inserter = VariantInserter(
            background_genome=Path("magestic.fasta"),
            source_strain_name="YJM454a",
            vcf_file=Path("parents.vcf.gz"),
        )

        result = inserter.analyze_cds(
            cds_name="ENA1",
            cds_start=538000,
            cds_end=541500,
            ref_protein=s288c_ena1_protein,  # For position mapping
        )
    """

    def __init__(
        self,
        background_genome: Optional[Path] = None,
        background_strain_name: str = "background",
        source_strain_name: str = "source",
        vcf_file: Optional[Path] = None,
        max_linked_distance_bp: int = 50,
        max_combined_distance_aa: int = 3,
    ):
        """Initialize VariantInserter.

        Args:
            background_genome: Path to background strain genome FASTA
            background_strain_name: Name of background strain (for VCF lookup)
            source_strain_name: Name of source strain (for VCF lookup)
            vcf_file: Path to VCF file with variant calls
            max_linked_distance_bp: Max bp distance for DNA-level variant linking
            max_combined_distance_aa: Max AA distance for protein-level combining
        """
        self.background_genome_path = background_genome
        self.background_strain_name = background_strain_name
        self.source_strain_name = source_strain_name
        self.vcf_file = vcf_file
        self.max_linked_distance_bp = max_linked_distance_bp
        self.max_combined_distance_aa = max_combined_distance_aa

        # Load genome if provided
        self.background_genome: Dict[str, str] = {}
        if background_genome and background_genome.exists():
            self.background_genome = read_fasta(background_genome)

    def extract_variants_from_vcf(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> Tuple[List[Variant], List[Variant]]:
        """Extract variants for background and source strains from VCF.

        Args:
            chrom: Chromosome name
            start: Start position (1-based)
            end: End position (1-based)

        Returns:
            Tuple of (background_variants, source_variants)
        """
        if self.vcf_file is None:
            raise ValueError("VCF file not provided")

        # Load all variants in region
        all_variants = load_variants_for_region(
            self.vcf_file,
            chrom=chrom,
            start=start,
            end=end,
            samples_to_load={self.background_strain_name, self.source_strain_name}
        )

        # Get strain-specific variants
        bg_variants = get_strain_variants(all_variants, self.background_strain_name)
        src_variants = get_strain_variants(all_variants, self.source_strain_name)

        return bg_variants, src_variants

    def extract_pairwise_variants(
        self,
        bg_seq: str,
        src_seq: str,
        bg_vcf_variants: List[Variant],
        src_vcf_variants: List[Variant],
        s288c_start: int,
        s288c_end: int,
    ) -> List[PairwiseVariant]:
        """Extract pairwise variants between background and source.

        Uses VCF-guided coordinate mapping to properly handle indels.

        Args:
            bg_seq: Background strain sequence
            src_seq: Source strain sequence
            bg_vcf_variants: VCF variants for background strain
            src_vcf_variants: VCF variants for source strain
            s288c_start: S288C start coordinate
            s288c_end: S288C end coordinate

        Returns:
            List of PairwiseVariant objects
        """
        return extract_pairwise_variants_vcf_guided(
            bg_seq=bg_seq,
            src_seq=src_seq,
            background_strain=self.background_strain_name,
            source_strain=self.source_strain_name,
            bg_vcf_variants=bg_vcf_variants,
            src_vcf_variants=src_vcf_variants,
            s288c_start=s288c_start,
            s288c_end=s288c_end,
        )

    def detect_linked_variants(
        self,
        variants: List[PairwiseVariant],
    ) -> List[LinkedVariantGroup]:
        """Detect linked variant groups (nearby variants that co-occur).

        Args:
            variants: List of pairwise variants

        Returns:
            List of LinkedVariantGroup objects
        """
        return detect_linked_variants(variants, max_distance=self.max_linked_distance_bp)

    def insert_variant(
        self,
        sequence: str,
        variant: PairwiseVariant,
    ) -> Optional[str]:
        """Insert a single variant into a sequence.

        Args:
            sequence: Background sequence
            variant: Variant to insert

        Returns:
            Modified sequence, or None if insertion failed
        """
        return apply_pairwise_variant(sequence, variant)

    def insert_variants(
        self,
        sequence: str,
        variants: List[PairwiseVariant],
    ) -> Tuple[str, List[str]]:
        """Insert multiple variants into a sequence.

        Variants are applied right-to-left to maintain coordinate consistency.

        Args:
            sequence: Background sequence
            variants: List of variants to insert

        Returns:
            Tuple of (modified_sequence, list_of_applied_variant_ids)
        """
        return apply_multiple_pairwise_variants(sequence, variants)

    def analyze_protein_effect(
        self,
        background_cds: str,
        mutant_cds: str,
        ref_protein: Optional[str] = None,
    ) -> Tuple[str, List[AAVariant], List[CombinedAAVariant], Dict[int, int]]:
        """Analyze protein-level effect of CDS variants.

        Args:
            background_cds: Background CDS sequence
            mutant_cds: Mutant CDS sequence
            ref_protein: Reference protein for position mapping (e.g., S288C)

        Returns:
            Tuple of:
            - variant_classification: missense, nonsense, frameshift, etc.
            - aa_variants: List of AA changes
            - combined_variants: List of combined nearby AA changes
            - position_mapping: Dict mapping background AA pos to ref AA pos
        """
        # Translate both
        bg_protein = translate(background_cds)
        mut_protein = translate(mutant_cds)

        # Check for frameshift
        net_change = len(mutant_cds) - len(background_cds)

        if net_change % 3 != 0:
            # Frameshift
            direction, desc = classify_frameshift_direction(bg_protein, mut_protein)
            return direction, [], [], {}

        # Check for start loss
        if not mutant_cds.upper().startswith('ATG'):
            consequence, desc = classify_start_loss(mutant_cds)
            return consequence, [], [], {}

        # Find AA differences
        aa_variants, length_diff = find_aa_differences(bg_protein, mut_protein, include_length_diff=True)

        # Classify based on AA changes
        if not aa_variants and length_diff == 0:
            classification = 'synonymous'
        elif any(v.is_nonsense for v in aa_variants):
            classification = 'nonsense'
        elif any(v.is_stop_loss for v in aa_variants):
            classification = 'stop_loss'
        elif length_diff != 0:
            if length_diff > 0:
                classification = 'aa_insertion'
            else:
                classification = 'aa_deletion'
        else:
            classification = 'missense'

        # Combine nearby AA changes
        combined = combine_nearby_aa_changes(aa_variants, max_distance=self.max_combined_distance_aa)

        # Build position mapping if reference provided
        pos_mapping = {}
        if ref_protein:
            full_mapping = build_position_mapping(bg_protein, ref_protein)
            for aa_var in aa_variants:
                if aa_var.position in full_mapping:
                    pos_mapping[aa_var.position] = full_mapping[aa_var.position].ref_pos

        return classification, aa_variants, combined, pos_mapping

    def analyze_cds(
        self,
        gene_name: str,
        background_cds: str,
        source_cds: str,
        pairwise_variants: List[PairwiseVariant],
        ref_protein: Optional[str] = None,
        cds_start_in_locus: int = 0,
    ) -> VariantAnalysisResult:
        """Full analysis pipeline for a CDS.

        1. Detect linked variants
        2. For each single variant and linked group:
           a. Insert variant(s) into background CDS
           b. Translate to protein
           c. Find AA differences
           d. Combine nearby AA changes
           e. Classify variant effect
        3. Build position mapping to reference

        Args:
            gene_name: Gene/ORF name
            background_cds: Background strain CDS
            source_cds: Source strain CDS (for validation)
            pairwise_variants: List of variants between strains
            ref_protein: Reference protein for position mapping
            cds_start_in_locus: Start position of CDS within the locus sequence

        Returns:
            VariantAnalysisResult with all test records
        """
        result = VariantAnalysisResult(
            gene_name=gene_name,
            background_strain=self.background_strain_name,
            source_strain=self.source_strain_name,
        )

        # Filter variants to those within CDS
        # FIX (2026-02-23): Create copies instead of mutating original objects
        # This prevents coordinate corruption when variants are reused
        cds_variants = []
        for v in pairwise_variants:
            if cds_start_in_locus <= v.bg_seq_pos < cds_start_in_locus + len(background_cds):
                # Create a copy with adjusted position (don't mutate original)
                adjusted_variant = PairwiseVariant(
                    background_strain=v.background_strain,
                    source_strain=v.source_strain,
                    bg_seq_pos=v.bg_seq_pos - cds_start_in_locus,  # Adjusted to CDS-relative
                    bg_allele=v.bg_allele,
                    src_allele=v.src_allele,
                    s288c_pos=v.s288c_pos
                )
                cds_variants.append(adjusted_variant)

        if not cds_variants:
            return result

        # Detect linked variants
        linked_groups = self.detect_linked_variants(cds_variants)
        result.linked_groups = linked_groups

        # Process each variant/group
        for group in linked_groups:
            if group.n_variants == 1:
                # Single variant test
                test = self._create_test_record(
                    gene_name=gene_name,
                    variants=[group.variants[0]],
                    background_cds=background_cds,
                    variant_type='single',
                    ref_protein=ref_protein,
                )
                result.single_variant_tests.append(test)
            else:
                # Linked variant group
                # First, create complete variant test (all variants together)
                complete_test = self._create_test_record(
                    gene_name=gene_name,
                    variants=group.variants,
                    background_cds=background_cds,
                    variant_type='complete',
                    linked_group_id=group.group_id,
                    is_complete=True,
                    ref_protein=ref_protein,
                )
                result.complete_variant_tests.append(complete_test)

                # Then, create sub-variant tests (each variant individually)
                for var in group.variants:
                    sub_test = self._create_test_record(
                        gene_name=gene_name,
                        variants=[var],
                        background_cds=background_cds,
                        variant_type='sub_variant',
                        linked_group_id=group.group_id,
                        is_complete=False,
                        ref_protein=ref_protein,
                    )
                    result.sub_variant_tests.append(sub_test)

        # Build position mapping for background protein
        if ref_protein:
            bg_protein = translate(background_cds)
            result.position_mapping = build_position_mapping(bg_protein, ref_protein)

        return result

    def _create_test_record(
        self,
        gene_name: str,
        variants: List[PairwiseVariant],
        background_cds: str,
        variant_type: str,
        linked_group_id: Optional[str] = None,
        is_complete: bool = False,
        ref_protein: Optional[str] = None,
    ) -> VariantTestRecord:
        """Create a test record for a variant or variant group.

        Args:
            gene_name: Gene/ORF name
            variants: List of variants to test
            background_cds: Background CDS sequence
            variant_type: 'single', 'complete', or 'sub_variant'
            linked_group_id: ID of linked group (if applicable)
            is_complete: Whether this tests all variants in a linked group
            ref_protein: Reference protein for position mapping

        Returns:
            VariantTestRecord with full analysis
        """
        # Generate test ID
        var_ids = '_'.join(str(v.s288c_pos) for v in sorted(variants, key=lambda x: x.s288c_pos))
        test_id = f"{gene_name}_{self.background_strain_name}_from_{self.source_strain_name}_v{var_ids}"

        # Insert variants into CDS
        mutant_cds, applied_ids = self.insert_variants(background_cds, variants)

        # Analyze protein effect
        classification, aa_variants, combined, pos_mapping = self.analyze_protein_effect(
            background_cds, mutant_cds, ref_protein
        )

        # Get additional classification details
        frameshift_dir = None
        start_loss_cons = None

        if 'frameshift' in classification:
            bg_protein = translate(background_cds)
            mut_protein = translate(mutant_cds)
            frameshift_dir, _ = classify_frameshift_direction(bg_protein, mut_protein)
        elif 'start_loss' in classification:
            start_loss_cons, _ = classify_start_loss(mutant_cds)

        return VariantTestRecord(
            test_id=test_id,
            gene_name=gene_name,
            background_strain=self.background_strain_name,
            source_strain=self.source_strain_name,
            variant_type=variant_type,
            dna_variants=variants,
            linked_group_id=linked_group_id,
            is_complete_variant=is_complete,
            background_cds=background_cds,
            mutant_cds=mutant_cds,
            background_protein=translate(background_cds),
            mutant_protein=translate(mutant_cds),
            aa_variants=aa_variants,
            combined_aa_variants=combined,
            variant_classification=classification,
            frameshift_direction=frameshift_dir,
            start_loss_consequence=start_loss_cons,
            ref_position_mapping=pos_mapping,
        )


# ============================================================================
# Convenience Functions
# ============================================================================

def analyze_strain_pair(
    background_cds: str,
    source_cds: str,
    background_strain: str,
    source_strain: str,
    gene_name: str,
    bg_vcf_variants: List[Variant],
    src_vcf_variants: List[Variant],
    s288c_start: int,
    s288c_end: int,
    ref_protein: Optional[str] = None,
    max_linked_distance_bp: int = 50,
    max_combined_distance_aa: int = 3,
) -> VariantAnalysisResult:
    """Convenience function to analyze a strain pair without creating VariantInserter.

    Args:
        background_cds: Background CDS sequence
        source_cds: Source CDS sequence
        background_strain: Background strain name
        source_strain: Source strain name
        gene_name: Gene/ORF name
        bg_vcf_variants: VCF variants for background strain
        src_vcf_variants: VCF variants for source strain
        s288c_start: S288C start coordinate
        s288c_end: S288C end coordinate
        ref_protein: Reference protein for position mapping
        max_linked_distance_bp: Max bp distance for variant linking
        max_combined_distance_aa: Max AA distance for combining

    Returns:
        VariantAnalysisResult with all tests
    """
    # Create inserter
    inserter = VariantInserter(
        background_strain_name=background_strain,
        source_strain_name=source_strain,
        max_linked_distance_bp=max_linked_distance_bp,
        max_combined_distance_aa=max_combined_distance_aa,
    )

    # Need locus sequences (not just CDS) for VCF-guided extraction
    # For this convenience function, assume CDS = locus
    pairwise_variants = inserter.extract_pairwise_variants(
        bg_seq=background_cds,
        src_seq=source_cds,
        bg_vcf_variants=bg_vcf_variants,
        src_vcf_variants=src_vcf_variants,
        s288c_start=s288c_start,
        s288c_end=s288c_end,
    )

    return inserter.analyze_cds(
        gene_name=gene_name,
        background_cds=background_cds,
        source_cds=source_cds,
        pairwise_variants=pairwise_variants,
        ref_protein=ref_protein,
        cds_start_in_locus=0,
    )


# ============================================================================
# Tests
# ============================================================================

if __name__ == "__main__":
    print("Testing VariantInserter module...")

    # Test position mapping for same-length proteins
    print("\n1. Testing position mapping (same length)...")
    ref = "MAVIKEIVH"
    strain = "MTVIKELVH"  # A2T, E7L
    mapping = build_position_mapping(strain, ref)
    assert len(mapping) == len(strain), f"Expected {len(strain)} mappings, got {len(mapping)}"
    assert mapping[2].ref_pos == 2, f"Expected pos 2 to map to ref 2"
    assert mapping[2].strain_aa == 'T', f"Expected strain AA at pos 2 to be T"
    assert mapping[2].ref_aa == 'A', f"Expected ref AA at pos 2 to be A"
    print("✓ Same-length position mapping: PASS")

    # Test position mapping for different-length proteins
    print("\n2. Testing position mapping (different length)...")
    ref = "MAVIKEIVH"
    strain = "MAVIKE"  # Shorter
    mapping = build_position_mapping(strain, ref)
    assert len(mapping) == len(strain), f"Expected {len(strain)} mappings, got {len(mapping)}"
    assert mapping[1].ref_pos == 1
    print("✓ Different-length position mapping: PASS")

    # Test map_strain_position_to_reference
    print("\n3. Testing map_strain_position_to_reference...")
    ref_pos, quality = map_strain_position_to_reference("MAVIKEIVH", "MAVIKEIVH", 5)
    assert ref_pos == 5, f"Expected ref_pos 5, got {ref_pos}"
    assert quality == 1.0, f"Expected quality 1.0, got {quality}"
    print("✓ map_strain_position_to_reference: PASS")

    # Test VariantInserter basic functionality
    print("\n4. Testing VariantInserter...")
    inserter = VariantInserter(
        background_strain_name="BYa",
        source_strain_name="YJM454a",
        max_linked_distance_bp=50,
        max_combined_distance_aa=3,
    )

    # Create test variants (mimicking ENA1)
    test_variants = [
        PairwiseVariant("BYa", "YJM454a", 100, "TGAGCACATTG", "T", 538149),  # -10bp
        PairwiseVariant("BYa", "YJM454a", 101, "G", "GA", 538163),           # +1bp
        PairwiseVariant("BYa", "YJM454a", 102, "A", "ACCCACGAC", 538171),    # +8bp
        PairwiseVariant("BYa", "YJM454a", 103, "A", "AG", 538174),           # +1bp
    ]

    # Detect linked variants
    groups = inserter.detect_linked_variants(test_variants)
    assert len(groups) == 1, f"Expected 1 linked group, got {len(groups)}"
    assert groups[0].net_length_change == 0, f"Expected net 0bp, got {groups[0].net_length_change}"
    assert groups[0].is_in_frame, "Expected in-frame"
    print("✓ Linked variant detection: PASS")

    # Test classify_linked_effect via the group
    var_type, effect, net = classify_linked_effect(groups[0])
    assert var_type == 'complex_indel_neutral', f"Expected complex_indel_neutral, got {var_type}"
    print("✓ Linked effect classification: PASS")

    # Test VariantTestRecord
    print("\n5. Testing VariantTestRecord...")
    test_record = VariantTestRecord(
        test_id="test_1",
        gene_name="ENA1",
        background_strain="BYa",
        source_strain="YJM454a",
        variant_type="single",
        dna_variants=test_variants[:1],
        background_cds="ATGCATGC",
        mutant_cds="ATGC",
        background_protein="MH",
        mutant_protein="M",
    )
    assert test_record.net_dna_length_change == -10
    assert test_record.protein_length_change == -1
    print("✓ VariantTestRecord: PASS")

    # Test VariantAnalysisResult
    print("\n6. Testing VariantAnalysisResult...")
    result = VariantAnalysisResult(
        gene_name="ENA1",
        background_strain="BYa",
        source_strain="YJM454a",
    )
    result.single_variant_tests.append(test_record)
    assert result.n_total_tests == 1
    print("✓ VariantAnalysisResult: PASS")

    # Test analyze_protein_effect
    print("\n7. Testing analyze_protein_effect...")
    bg_cds = "ATGGCGGTGTAA"  # MAV*
    mut_cds = "ATGTCGGTGTAA"  # MSV* (A2S missense)

    classification, aa_variants, combined, pos_mapping = inserter.analyze_protein_effect(
        bg_cds, mut_cds
    )
    assert classification == 'missense', f"Expected missense, got {classification}"
    assert len(aa_variants) == 1, f"Expected 1 AA variant, got {len(aa_variants)}"
    assert str(aa_variants[0]) == "A2S", f"Expected A2S, got {aa_variants[0]}"
    print("✓ analyze_protein_effect (missense): PASS")

    # Test frameshift detection
    print("\n8. Testing frameshift detection...")
    bg_cds = "ATGGCGGTGTAA"  # MAV*
    mut_cds = "ATGGCGTGTAA"   # Frameshift (removed one base)

    classification, _, _, _ = inserter.analyze_protein_effect(bg_cds, mut_cds)
    assert 'frameshift' in classification, f"Expected frameshift, got {classification}"
    print("✓ analyze_protein_effect (frameshift): PASS")

    # Test synonymous detection
    print("\n9. Testing synonymous detection...")
    bg_cds = "ATGGCGTAA"  # MA*
    mut_cds = "ATGGCATAA"  # MA* (GCA also codes for A)

    classification, aa_variants, _, _ = inserter.analyze_protein_effect(bg_cds, mut_cds)
    assert classification == 'synonymous', f"Expected synonymous, got {classification}"
    assert len(aa_variants) == 0, f"Expected 0 AA variants for synonymous"
    print("✓ analyze_protein_effect (synonymous): PASS")

    print("\n✅ All variant_inserter.py tests passed!")
