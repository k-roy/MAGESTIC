"""
Alignment parsing and analysis utilities.

Provides functions for parsing SAM files and comparing alignments
from multiple aligners.

Author: Kevin R. Roy
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .sequence_utils import parse_cigar


@dataclass
class AlignmentScore:
    """Unified alignment score for comparing aligners."""
    aligner: str
    ref_start: int = 0
    ref_end: int = 0
    mismatches: int = 0
    soft_clips: int = 0
    insertions: int = 0
    deletions: int = 0
    is_aligned: bool = False
    cigar: str = ""
    aligned_seq: str = ""
    md_tag: str = ""
    mate_unmapped: bool = False  # For paired-end: is the mate unmapped?

    @property
    def penalty_score(self) -> int:
        """Score for ranking alignments. Lower is better.
        Soft clips are penalized equally to mismatches."""
        if not self.is_aligned:
            return float('inf')
        return self.mismatches + self.soft_clips

    @property
    def priority(self) -> int:
        """Priority for tie-breaking: bwa=0 (highest), bbmap=1, minimap2=2"""
        priority_map = {'bwa': 0, 'bbmap': 1, 'minimap2': 2}
        return priority_map.get(self.aligner, 99)


@dataclass
class AlignmentInfo:
    """Complete information about a read alignment."""
    read_name: str
    flag: int
    ref_name: str
    ref_start: int
    mapq: int
    cigar: str
    sequence: str
    quality: str
    is_aligned: bool
    is_reverse: bool
    is_read1: bool = True
    is_read2: bool = False
    mate_unmapped: bool = False
    md_tag: str = ""


def extract_alignment_info(cigar: str, seq: str, ref_pos: int, md_tag: str = "") -> AlignmentScore:
    """
    Extract alignment details from CIGAR and MD tag.

    Returns AlignmentScore with positions, mismatches, soft clips, etc.
    """
    result = AlignmentScore(aligner="", ref_start=ref_pos, is_aligned=True, cigar=cigar, md_tag=md_tag)

    parsed = parse_cigar(cigar)
    ref_consumed = 0
    soft_clips = 0
    mismatches = 0
    insertions = 0
    deletions = 0

    for length, op in parsed:
        if op == 'M' or op == '=':
            ref_consumed += length
        elif op == 'X':  # Mismatch in extended CIGAR
            ref_consumed += length
            mismatches += length
        elif op == 'I':
            insertions += length
        elif op == 'D' or op == 'N':
            ref_consumed += length
            deletions += length
        elif op == 'S':
            soft_clips += length
        elif op == 'H':
            pass  # Hard clip doesn't consume anything

    result.ref_end = ref_pos + ref_consumed
    result.soft_clips = soft_clips
    result.insertions = insertions
    result.deletions = deletions

    # Parse MD tag for mismatches if available
    if md_tag and mismatches == 0:
        # MD tag format: numbers for matches, letters for reference bases at mismatches
        # e.g., "50A20" means 50 matches, then A (ref was A, read has mismatch), then 20 matches
        md_mismatches = len(re.findall(r'[ACGT]', md_tag.upper()))
        mismatches = md_mismatches

    result.mismatches = mismatches

    return result


def parse_sam_file(sam_path: Path, aligner: str) -> Dict[str, AlignmentScore]:
    """
    Parse SAM file and extract alignment scores for each read.

    For paired-end SAM files, uses FLAG bits to distinguish R1 from R2:
    - FLAG 0x40 (64) = first in pair (R1)
    - FLAG 0x80 (128) = second in pair (R2)

    Returns dict of read_id -> AlignmentScore, where read_id includes _R1 or _R2 suffix
    for paired-end reads.
    """
    results = {}

    with open(sam_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            # Normalize read ID (some aligners include full Illumina header with spaces)
            base_read_id = fields[0].split()[0]
            flag = int(fields[1])
            ref_name = fields[2]
            ref_pos = int(fields[3]) - 1  # Convert to 0-based
            cigar = fields[5]
            seq = fields[9]

            # For paired-end reads, add suffix based on FLAG bits
            # FLAG 0x40 = first in pair (R1), FLAG 0x80 = second in pair (R2)
            if flag & 0x1:  # paired-end read
                if flag & 0x40:
                    read_id = f"{base_read_id}_R1"
                elif flag & 0x80:
                    read_id = f"{base_read_id}_R2"
                else:
                    read_id = base_read_id  # shouldn't happen in proper paired-end
            else:
                read_id = base_read_id  # single-end read

            # Skip unmapped reads
            if flag & 4 or ref_name == '*' or cigar == '*':
                results[read_id] = AlignmentScore(aligner=aligner, is_aligned=False)
                continue

            # Check if mate is mapped (for paired-end quality filtering)
            mate_unmapped = bool(flag & 0x8) if (flag & 0x1) else False

            # Extract MD tag if present
            md_tag = ""
            for field in fields[11:]:
                if field.startswith('MD:Z:'):
                    md_tag = field[5:]
                    break

            # Parse alignment
            score = extract_alignment_info(cigar, seq, ref_pos, md_tag)
            score.aligner = aligner
            score.aligned_seq = seq
            score.mate_unmapped = mate_unmapped  # Track mate status

            results[read_id] = score

    return results


def select_best_alignment(
    bwa_results: Dict[str, AlignmentScore],
    bbmap_results: Dict[str, AlignmentScore],
    minimap2_results: Dict[str, AlignmentScore]
) -> Dict[str, AlignmentScore]:
    """
    Select the best alignment for each read across all three aligners.

    Selection criteria:
    - Score = mismatches + soft_clips (lower is better)
    - Ties go to bwa (priority 0), then bbmap (priority 1), then minimap2 (priority 2)
    """
    all_read_ids = set(bwa_results.keys()) | set(bbmap_results.keys()) | set(minimap2_results.keys())

    best_results = {}

    for read_id in all_read_ids:
        candidates = []

        if read_id in bwa_results:
            candidates.append(bwa_results[read_id])
        if read_id in bbmap_results:
            candidates.append(bbmap_results[read_id])
        if read_id in minimap2_results:
            candidates.append(minimap2_results[read_id])

        # Filter to aligned candidates
        aligned_candidates = [c for c in candidates if c.is_aligned]

        if not aligned_candidates:
            # No aligner could align this read
            best_results[read_id] = AlignmentScore(aligner='none', is_aligned=False)
            continue

        # Sort by (penalty_score, priority) - lower is better for both
        aligned_candidates.sort(key=lambda x: (x.penalty_score, x.priority))

        best_results[read_id] = aligned_candidates[0]

    return best_results


def select_best_alignment_paired(
    bwa_r1: AlignmentScore, bwa_r2: AlignmentScore,
    bbmap_r1: AlignmentScore, bbmap_r2: AlignmentScore,
    minimap2_r1: AlignmentScore, minimap2_r2: AlignmentScore
) -> Tuple[str, AlignmentScore, AlignmentScore]:
    """
    Select best aligner based on combined R1+R2 penalty scores.

    Returns (aligner_name, best_r1, best_r2)
    """
    # Priority order for tie-breaking: bwa (0) > bbmap (1) > minimap2 (2)
    ALIGNER_PRIORITY = {'bwa': 0, 'bbmap': 1, 'minimap2': 2}

    def combined_score(r1: AlignmentScore, r2: AlignmentScore) -> int:
        """Calculate combined penalty score for a read pair."""
        if not r1.is_aligned or not r2.is_aligned:
            return float('inf')
        return r1.penalty_score + r2.penalty_score

    scores = [
        ('bwa', combined_score(bwa_r1, bwa_r2), bwa_r1, bwa_r2),
        ('bbmap', combined_score(bbmap_r1, bbmap_r2), bbmap_r1, bbmap_r2),
        ('minimap2', combined_score(minimap2_r1, minimap2_r2), minimap2_r1, minimap2_r2),
    ]

    # Sort by score (lower is better), then by priority for tie-breaking
    scores.sort(key=lambda x: (x[1], ALIGNER_PRIORITY[x[0]]))

    best_aligner, _, best_r1, best_r2 = scores[0]
    return best_aligner, best_r1, best_r2
