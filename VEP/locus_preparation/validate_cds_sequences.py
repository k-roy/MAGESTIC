#!/usr/bin/env python3
"""
Validate CDS sequences extracted from strain loci - Version 3.

Handles:
1. Multi-exon genes (spliced genes with introns)
2. Frameshift variants (extend to next stop codon)
3. Start codon variants (find next ATG)

Key improvements in v3:
- Frameshifts are FIXED by extending to next stop codon
- Invalid start codons are FIXED by finding next ATG
- Outputs both original and corrected CDS info

Author: Kevin R. Roy
Date: 2026-02-22
"""

import json
import re
from pathlib import Path
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass, field
import sys

# Add sequence_utils to path
BASE_DIR = Path("/path/to")
sys.path.insert(0, str(BASE_DIR / "common/scripts/variant_effect_prediction"))

from sequence_utils import reverse_complement, translate

# Project paths
PROJECT_DIR = BASE_DIR / "projects/QTL/20260210_variant_effect_prediction"
DATA_DIR = PROJECT_DIR / "processed_data/20260205_in_silico_QTN_dissection"
LOCI_DIR = DATA_DIR / "loci_50kb/fasta_50kb"
GFF_DIR = DATA_DIR / "strain_gff"
OFFSET_DIR = DATA_DIR / "strain_genomes/offset_tables"
METADATA_FILE = DATA_DIR / "loci_50kb/metadata/loci_metadata_50kb.tsv"
S288C_GFF = BASE_DIR / "common/reference_genomes/S288C_reference_genome_R64-5-1_20240529/saccharomyces_cerevisiae_R64-5-1_20240529.gff"

# Valid codons
START_CODONS = {'ATG'}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


@dataclass
class ExonInfo:
    """Single exon from GFF."""
    start: int  # 1-based, S288C coordinates
    end: int    # 1-based, inclusive


@dataclass
class GeneStructure:
    """Complete gene structure from S288C GFF."""
    orf: str
    chrom: str
    strand: str
    exons: List[ExonInfo] = field(default_factory=list)
    introns: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def is_spliced(self) -> bool:
        return len(self.exons) > 1

    @property
    def total_cds_length(self) -> int:
        return sum(e.end - e.start + 1 for e in self.exons)

    @property
    def cds_start(self) -> int:
        """5' end of CDS in S288C coordinates."""
        if self.strand == '+':
            return min(e.start for e in self.exons)
        else:
            return max(e.end for e in self.exons)

    @property
    def cds_end(self) -> int:
        """3' end of CDS in S288C coordinates."""
        if self.strand == '+':
            return max(e.end for e in self.exons)
        else:
            return min(e.start for e in self.exons)


def load_s288c_gene_structures(orfs: set) -> Dict[str, GeneStructure]:
    """Load gene structures from S288C GFF for specified ORFs.

    Uses a two-pass approach to avoid matching ORF names in Note fields.
    """
    structures = {}

    # First pass: Find gene features
    with open(S288C_GFF) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            attrs = dict(item.split('=', 1) for item in parts[8].split(';') if '=' in item)
            gene_id = attrs.get('ID', '')
            if gene_id in orfs:
                structures[gene_id] = GeneStructure(
                    orf=gene_id, chrom=parts[0], strand=parts[6]
                )

    # Second pass: Find CDS and intron features
    with open(S288C_GFF) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] not in ('CDS', 'intron'):
                continue

            start, end = int(parts[3]), int(parts[4])
            attrs = dict(item.split('=', 1) for item in parts[8].split(';') if '=' in item)
            parent = attrs.get('Parent', '')

            for orf in orfs:
                if orf in structures:
                    # Match Parent patterns like YOR140W_mRNA or YOR140W_id001
                    if parent.startswith(orf + '_') or (',' + orf + '_') in parent:
                        if parts[2] == 'CDS':
                            structures[orf].exons.append(ExonInfo(start=start, end=end))
                        elif parts[2] == 'intron':
                            structures[orf].introns.append((start, end))
                        break

    # Sort exons by genomic position
    for struct in structures.values():
        struct.exons.sort(key=lambda e: e.start)

    return structures


def load_offset_table(strain: str, chrom: str) -> Dict[int, int]:
    """Load offset table for a strain/chromosome pair."""
    offset_file = OFFSET_DIR / f"{strain}_{chrom}_offsets.json"
    if not offset_file.exists():
        return {}
    with open(offset_file) as f:
        raw_offsets = json.load(f)
    return {int(k): v for k, v in raw_offsets.items()}


def get_offset_at_position(offset_table: Dict[int, int], s288c_pos: int) -> int:
    """Get cumulative offset at an S288C position."""
    if not offset_table:
        return 0
    keys = sorted(offset_table.keys())
    offset = 0
    for k in keys:
        if k <= s288c_pos:
            offset = offset_table[k]
        else:
            break
    return offset


def read_fasta_sequences(fasta_file: Path) -> Dict[str, str]:
    """Read all sequences from a FASTA file."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def load_loci_metadata():
    """Load loci metadata."""
    import pandas as pd
    df = pd.read_csv(METADATA_FILE, sep='\t')
    metadata = {}
    for _, row in df.iterrows():
        key = (row['orf'], row['strain'])
        metadata[key] = row
    return metadata


def s288c_to_strain_pos(s288c_pos: int, offset_table: Dict[int, int]) -> int:
    """Convert S288C position to strain position."""
    return s288c_pos + get_offset_at_position(offset_table, s288c_pos)


def extract_extended_cds(
    locus_seq: str,
    locus_s288c_start: int,
    gene_struct: GeneStructure,
    offset_table: Dict[int, int],
    extension_bp: int = 2000
) -> Tuple[str, str, int, int]:
    """Extract CDS with extension capability for frameshift handling.

    Returns:
        (cds_seq, downstream_seq, cds_start_in_locus, cds_end_in_locus)

    The cds_seq is the initial extraction based on S288C coordinates.
    The downstream_seq provides extra sequence for extending frameshifts.
    """
    offset_at_locus = get_offset_at_position(offset_table, locus_s288c_start)
    locus_strain_start = locus_s288c_start + offset_at_locus

    # For single-exon genes (non-spliced)
    if not gene_struct.is_spliced:
        exon = gene_struct.exons[0]

        strain_cds_start = s288c_to_strain_pos(exon.start, offset_table)
        strain_cds_end = s288c_to_strain_pos(exon.end, offset_table)

        cds_start_in_locus = strain_cds_start - locus_strain_start
        cds_end_in_locus = strain_cds_end - locus_strain_start

        if cds_start_in_locus < 0 or cds_end_in_locus >= len(locus_seq):
            return "", "", -1, -1

        cds_seq = locus_seq[cds_start_in_locus:cds_end_in_locus + 1]

        # Get downstream sequence for extension
        if gene_struct.strand == '+':
            downstream_end = min(cds_end_in_locus + extension_bp, len(locus_seq) - 1)
            downstream_seq = locus_seq[cds_end_in_locus + 1:downstream_end + 1]
        else:
            # For minus strand, "downstream" is actually upstream in genomic coords
            downstream_start = max(cds_start_in_locus - extension_bp, 0)
            downstream_seq = locus_seq[downstream_start:cds_start_in_locus]

        # Handle strand
        if gene_struct.strand == '-':
            cds_seq = reverse_complement(cds_seq)
            downstream_seq = reverse_complement(downstream_seq)

        return cds_seq, downstream_seq, cds_start_in_locus, cds_end_in_locus

    # For multi-exon (spliced) genes
    else:
        sorted_exons = sorted(gene_struct.exons, key=lambda e: e.start)
        exon_seqs = []

        for exon in sorted_exons:
            strain_exon_start = s288c_to_strain_pos(exon.start, offset_table)
            strain_exon_end = s288c_to_strain_pos(exon.end, offset_table)

            exon_start_in_locus = strain_exon_start - locus_strain_start
            exon_end_in_locus = strain_exon_end - locus_strain_start

            if exon_start_in_locus < 0 or exon_end_in_locus >= len(locus_seq):
                return "", "", -1, -1

            exon_seqs.append(locus_seq[exon_start_in_locus:exon_end_in_locus + 1])

        # Get downstream sequence (after last exon)
        last_exon = sorted_exons[-1]
        strain_last_end = s288c_to_strain_pos(last_exon.end, offset_table)
        last_end_in_locus = strain_last_end - locus_strain_start

        if gene_struct.strand == '+':
            downstream_end = min(last_end_in_locus + extension_bp, len(locus_seq) - 1)
            downstream_seq = locus_seq[last_end_in_locus + 1:downstream_end + 1]
        else:
            first_exon = sorted_exons[0]
            strain_first_start = s288c_to_strain_pos(first_exon.start, offset_table)
            first_start_in_locus = strain_first_start - locus_strain_start
            downstream_start = max(first_start_in_locus - extension_bp, 0)
            downstream_seq = locus_seq[downstream_start:first_start_in_locus]

        # Handle strand-specific joining
        if gene_struct.strand == '-':
            rc_exon_seqs = [reverse_complement(seq) for seq in exon_seqs]
            cds_seq = ''.join(reversed(rc_exon_seqs))
            downstream_seq = reverse_complement(downstream_seq)
        else:
            cds_seq = ''.join(exon_seqs)

        # Return positions of first exon
        first_start = s288c_to_strain_pos(sorted_exons[0].start, offset_table) - locus_strain_start
        last_end = s288c_to_strain_pos(sorted_exons[-1].end, offset_table) - locus_strain_start

        return cds_seq, downstream_seq, first_start, last_end


def find_stop_codon_in_frame(seq: str) -> Optional[int]:
    """Find first stop codon in reading frame (starting from position 0)."""
    for i in range(0, len(seq) - 2, 3):
        if seq[i:i+3] in STOP_CODONS:
            return i
    return None


def find_next_atg(seq: str, start_pos: int = 0) -> Optional[int]:
    """Find next ATG start codon."""
    for i in range(start_pos, len(seq) - 2):
        if seq[i:i+3] == 'ATG':
            return i
    return None


def fix_cds(cds_seq: str, downstream_seq: str) -> Tuple[str, List[str], List[str]]:
    """Fix CDS issues (frameshift, invalid start/stop).

    Returns:
        (fixed_cds, issues_found, fixes_applied)
    """
    issues = []
    fixes = []
    fixed_cds = cds_seq

    # Check start codon
    if len(cds_seq) >= 3 and cds_seq[:3] not in START_CODONS:
        issues.append(f'Invalid start codon: {cds_seq[:3]}')
        next_atg = find_next_atg(cds_seq, 0)
        if next_atg is not None and next_atg < len(cds_seq) // 2:
            frame = "in-frame" if next_atg % 3 == 0 else "out-of-frame"
            fixes.append(f'Start at next ATG +{next_atg} ({frame})')
            fixed_cds = cds_seq[next_atg:]
        else:
            fixes.append('No suitable ATG found')

    # Check if length is divisible by 3
    if len(fixed_cds) % 3 != 0:
        issues.append(f'Frameshift: length {len(fixed_cds)} not divisible by 3')

        # Try to extend to next stop codon
        extended = fixed_cds + downstream_seq
        stop_pos = find_stop_codon_in_frame(extended)

        if stop_pos is not None:
            proper_cds = extended[:stop_pos + 3]
            if len(proper_cds) % 3 == 0:
                original_len = len(fixed_cds)
                fixed_cds = proper_cds
                extension = len(fixed_cds) - original_len
                if extension > 0:
                    fixes.append(f'Extended by {extension}bp to stop codon')
                else:
                    fixes.append(f'Truncated to stop codon at {stop_pos}')
            else:
                fixes.append(f'Stop codon found but still not in frame')
        else:
            fixes.append('No stop codon found in extended region')

    # Check stop codon
    if len(fixed_cds) >= 3:
        if fixed_cds[-3:] not in STOP_CODONS:
            if 'Extended' not in str(fixes) and 'Truncated' not in str(fixes):
                issues.append(f'Invalid stop codon: {fixed_cds[-3:]}')

    return fixed_cds, issues, fixes


def validate_cds(cds_seq: str) -> Tuple[bool, str]:
    """Check if CDS is valid (ATG start, stop codon, divisible by 3)."""
    if len(cds_seq) < 6:
        return False, "Too short"

    problems = []

    if cds_seq[:3] not in START_CODONS:
        problems.append(f"start={cds_seq[:3]}")

    if cds_seq[-3:] not in STOP_CODONS:
        problems.append(f"stop={cds_seq[-3:]}")

    if len(cds_seq) % 3 != 0:
        problems.append(f"len={len(cds_seq)}")

    if problems:
        return False, "; ".join(problems)

    return True, "Valid"


def main():
    """Run CDS validation with frameshift fixing."""
    import pandas as pd

    print("Loading loci metadata...")
    metadata = load_loci_metadata()
    print(f"Found {len(metadata)} loci entries")

    orfs = sorted(set(k[0] for k in metadata.keys()))
    strains = sorted(set(k[1] for k in metadata.keys()))
    print(f"ORFs: {len(orfs)}, Strains: {len(strains)}")

    # Load gene structures from S288C
    print("\nLoading S288C gene structures...")
    gene_structures = load_s288c_gene_structures(set(orfs))
    print(f"Loaded {len(gene_structures)} gene structures")

    # Count spliced genes
    spliced = [orf for orf, gs in gene_structures.items() if gs.is_spliced]
    print(f"Spliced genes (multi-exon): {len(spliced)}")
    for orf in spliced:
        gs = gene_structures[orf]
        print(f"  {orf}: {len(gs.exons)} exons, strand {gs.strand}")

    results = []
    stats = {
        'total': 0,
        'valid_original': 0,
        'valid_after_fix': 0,
        'unfixable': 0
    }

    # Process each locus
    for orf in orfs:
        fasta_file = LOCI_DIR / f"{orf}_strain_loci_50kb.fasta"
        if not fasta_file.exists():
            continue

        gene_struct = gene_structures.get(orf)
        if not gene_struct:
            print(f"Warning: No gene structure for {orf}")
            continue

        sequences = read_fasta_sequences(fasta_file)

        for seq_id, locus_seq in sequences.items():
            parts = seq_id.split('_')
            strain = '_'.join(parts[1:])

            if (orf, strain) not in metadata:
                continue

            meta = metadata[(orf, strain)]
            chrom = meta['chrom']
            locus_s288c_start = meta['full_start_s288c']

            offset_table = load_offset_table(strain, chrom)

            # Extract CDS with downstream sequence for extension
            cds_seq, downstream_seq, _, _ = extract_extended_cds(
                locus_seq, locus_s288c_start, gene_struct, offset_table
            )

            stats['total'] += 1

            if not cds_seq:
                results.append({
                    'orf': orf,
                    'strain': strain,
                    'original_valid': False,
                    'fixed_valid': False,
                    'original_length': 0,
                    'fixed_length': 0,
                    'issues': 'Extraction failed',
                    'fixes': '',
                    'protein_length': 0,
                    'is_spliced': gene_struct.is_spliced
                })
                stats['unfixable'] += 1
                continue

            # Check original validity
            original_valid, original_status = validate_cds(cds_seq)
            if original_valid:
                stats['valid_original'] += 1

            # Try to fix
            fixed_cds, issues, fixes = fix_cds(cds_seq, downstream_seq)
            fixed_valid, _ = validate_cds(fixed_cds)

            if fixed_valid:
                stats['valid_after_fix'] += 1
            else:
                stats['unfixable'] += 1

            # Calculate protein length
            protein_length = 0
            if fixed_valid:
                try:
                    protein = translate(fixed_cds)
                    protein_length = len(protein.rstrip('*'))
                except:
                    pass

            results.append({
                'orf': orf,
                'strain': strain,
                'original_valid': original_valid,
                'fixed_valid': fixed_valid,
                'original_length': len(cds_seq),
                'fixed_length': len(fixed_cds),
                'original_start': cds_seq[:3] if len(cds_seq) >= 3 else '',
                'original_stop': cds_seq[-3:] if len(cds_seq) >= 3 else '',
                'fixed_start': fixed_cds[:3] if len(fixed_cds) >= 3 else '',
                'fixed_stop': fixed_cds[-3:] if len(fixed_cds) >= 3 else '',
                'issues': '; '.join(issues),
                'fixes': '; '.join(fixes),
                'protein_length': protein_length,
                'is_spliced': gene_struct.is_spliced
            })

    # Summary
    print(f"\n{'='*60}")
    print("CDS Validation Summary (v3 - with fixing)")
    print(f"{'='*60}")
    print(f"Total loci: {stats['total']}")
    print(f"Valid (original): {stats['valid_original']} ({100*stats['valid_original']/stats['total']:.1f}%)")
    print(f"Valid (after fix): {stats['valid_after_fix']} ({100*stats['valid_after_fix']/stats['total']:.1f}%)")
    print(f"Unfixable: {stats['unfixable']} ({100*stats['unfixable']/stats['total']:.1f}%)")

    # Save results
    results_df = pd.DataFrame(results)
    output_file = DATA_DIR / "loci_50kb/metadata/cds_validation_results_v3.tsv"
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults saved to: {output_file}")

    # Show examples of fixes applied
    fixed_cases = results_df[results_df['original_valid'] == False]
    if len(fixed_cases) > 0:
        print(f"\n{'='*60}")
        print("Cases with issues (showing first 10)")
        print(f"{'='*60}")
        for _, row in fixed_cases.head(10).iterrows():
            print(f"\n{row['orf']} in {row['strain']}:")
            print(f"  Issues: {row['issues']}")
            print(f"  Fixes: {row['fixes']}")
            print(f"  Length: {row['original_length']} -> {row['fixed_length']}")
            print(f"  Fixed valid: {row['fixed_valid']}")

    # Analyze unfixable cases
    unfixable = results_df[(results_df['fixed_valid'] == False)]
    if len(unfixable) > 0:
        print(f"\n{'='*60}")
        print(f"Unfixable cases: {len(unfixable)}")
        print(f"{'='*60}")
        for _, row in unfixable.head(10).iterrows():
            print(f"  {row['orf']} in {row['strain']}: {row['issues']}")

    return results_df


if __name__ == '__main__':
    main()
