"""
VCF file parsing utilities.

Classes:
- Variant: Dataclass representing a VCF variant

Functions:
- parse_vcf_genotype: Parse VCF genotype string
- load_variants_for_region: Load variants from VCF within a genomic region
- get_strain_variants: Filter variants to get strain-specific alleles

Author: Kevin R. Roy
"""

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union


@dataclass
class Variant:
    """Variant from VCF.

    Attributes:
        chrom: Chromosome name
        pos: 1-based genomic position
        ref: Reference allele
        alt: Alternate allele
        samples: Dict mapping sample names to their alleles
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    samples: Dict[str, str] = field(default_factory=dict)

    @property
    def variant_type(self) -> str:
        """Get variant type: SNP, MNP, or INDEL."""
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "SNP"
        elif len(self.ref) == len(self.alt):
            return "MNP"
        else:
            return "INDEL"

    @property
    def variant_id(self) -> str:
        """Get unique variant identifier."""
        return f"{self.pos}_{self.ref}_{self.alt}_{self.variant_type}"

    @property
    def length_change(self) -> int:
        """Get length change from ref to alt (positive for insertion)."""
        return len(self.alt) - len(self.ref)


def parse_vcf_genotype(gt_str: str, ref: str, alt_list: List[str]) -> str:
    """Parse VCF genotype string and return the corresponding allele.

    Args:
        gt_str: Genotype string (e.g., "0/1", "1|1", "./.")
        ref: Reference allele
        alt_list: List of alternate alleles

    Returns:
        The allele for this genotype (ref or one of alt_list)
    """
    if gt_str in ['.', './.', '.|.', './1', '1/.']:
        return ref

    gt_str = gt_str.replace('|', '/')

    if '/' in gt_str:
        alleles = gt_str.split('/')
    else:
        alleles = [gt_str]

    allele_indices = []
    for a in alleles:
        if a.isdigit():
            allele_indices.append(int(a))

    if not allele_indices:
        return ref

    max_idx = max(allele_indices)

    if max_idx == 0:
        return ref
    elif max_idx <= len(alt_list):
        return alt_list[max_idx - 1]
    else:
        return ref


def load_variants_for_region(
    vcf_path: Union[str, Path],
    chrom: str,
    start: int,
    end: int,
    samples_to_load: Optional[Set[str]] = None
) -> List[Variant]:
    """Load variants from VCF within a genomic region.

    Handles both gzipped and uncompressed VCF files.
    VCF positions are 1-based.

    Args:
        vcf_path: Path to VCF file
        chrom: Chromosome name (must match VCF format)
        start: Start position (1-based, inclusive)
        end: End position (1-based, inclusive)
        samples_to_load: Optional set of sample names to load. If None, loads all.

    Returns:
        List of Variant objects
    """
    vcf_path = Path(vcf_path)
    variants = []
    sample_names = []

    # Determine if gzipped
    is_gzipped = str(vcf_path).endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'

    with open_func(vcf_path, mode) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                sample_names = parts[9:]
                continue

            parts = line.strip().split('\t')
            if parts[0] != chrom:
                continue

            vcf_pos = int(parts[1])
            if vcf_pos < start or vcf_pos > end:
                continue

            ref = parts[3]
            alt_str = parts[4]
            alt_list = alt_str.split(',')

            # Parse format field
            format_fields = parts[8].split(':')
            gt_idx = format_fields.index('GT') if 'GT' in format_fields else 0

            # Parse samples
            samples = {}
            for i, sample in enumerate(sample_names):
                if samples_to_load and sample not in samples_to_load:
                    continue

                sample_data = parts[9 + i].split(':')
                gt = sample_data[gt_idx] if len(sample_data) > gt_idx else './.'
                allele = parse_vcf_genotype(gt, ref, alt_list)
                samples[sample] = allele

            # Skip if no samples have alternate allele
            alleles_used = set(samples.values())
            alleles_used.discard(ref)

            if not alleles_used:
                continue

            # Create variant for each alternate allele used
            for alt in alleles_used:
                # Skip symbolic/special alleles
                if alt in ('DEL', 'INS', 'DUP', 'INV', 'CNV', '*'):
                    continue
                if alt.startswith('<'):
                    continue

                variant = Variant(
                    chrom=chrom,
                    pos=vcf_pos,
                    ref=ref,
                    alt=alt,
                    samples=samples
                )
                variants.append(variant)

    return variants


def get_strain_variants(
    variants: List[Variant],
    strain: str,
    skip_ref_matches: bool = True
) -> List[Variant]:
    """Filter variants to get strain-specific alleles.

    Args:
        variants: List of all variants
        strain: Strain/sample name
        skip_ref_matches: If True, skip variants where strain matches reference

    Returns:
        List of variants where strain has alternate allele
    """
    strain_variants = []

    for var in variants:
        if strain not in var.samples:
            continue

        strain_allele = var.samples[strain]

        # Skip symbolic alleles
        if strain_allele in ('DEL', 'INS', 'DUP', 'INV', 'CNV', '*'):
            continue
        if strain_allele.startswith('<'):
            continue

        # Skip if strain matches reference
        if skip_ref_matches and strain_allele == var.ref:
            continue

        # Create new variant with strain's specific allele
        strain_var = Variant(
            chrom=var.chrom,
            pos=var.pos,
            ref=var.ref,
            alt=strain_allele,
            samples=var.samples
        )
        strain_variants.append(strain_var)

    return strain_variants


if __name__ == "__main__":
    # Test parse_vcf_genotype
    assert parse_vcf_genotype("0/0", "A", ["T"]) == "A"
    assert parse_vcf_genotype("0/1", "A", ["T"]) == "T"
    assert parse_vcf_genotype("1/1", "A", ["T"]) == "T"
    assert parse_vcf_genotype("./.", "A", ["T"]) == "A"
    assert parse_vcf_genotype("1|1", "A", ["T"]) == "T"  # phased
    assert parse_vcf_genotype("0/2", "A", ["T", "G"]) == "G"  # multi-allelic
    print("parse_vcf_genotype: PASS")

    # Test Variant dataclass
    v = Variant(chrom="chrI", pos=100, ref="A", alt="T")
    assert v.variant_type == "SNP"
    assert v.length_change == 0
    print("Variant SNP: PASS")

    v_indel = Variant(chrom="chrI", pos=100, ref="AT", alt="A")
    assert v_indel.variant_type == "INDEL"
    assert v_indel.length_change == -1
    print("Variant INDEL: PASS")

    v_ins = Variant(chrom="chrI", pos=100, ref="A", alt="ATG")
    assert v_ins.variant_type == "INDEL"
    assert v_ins.length_change == 2
    print("Variant insertion: PASS")

    print("\nAll vcf_parsing.py tests passed!")
