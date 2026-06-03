#!/usr/bin/env python3
"""
FASTA Deduplication Utility for VEP Pipeline

Identifies and removes duplicate sequences before VEP job submission.
Sequences are considered duplicates if their actual sequence content is identical,
regardless of header/ID differences.

Usage:
    python deduplicate_fasta.py <input.fasta> <output.fasta>
    python deduplicate_fasta.py --check <input.fasta>  # Just report stats
    python deduplicate_fasta.py --dir <fasta_dir>       # Process all FASTAs

Author: Kevin R. Roy
"""

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def hash_sequence(seq: str) -> str:
    """Create a hash of a sequence for comparison."""
    return hashlib.md5(seq.upper().encode()).hexdigest()


def read_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    """Read FASTA file, return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences


def deduplicate_sequences(sequences: List[Tuple[str, str]]) -> Tuple[List[Tuple[str, str]], Dict[str, List[str]]]:
    """
    Deduplicate sequences, keeping first occurrence.

    Returns:
        - List of unique (header, sequence) tuples
        - Dict mapping each unique sequence hash to list of duplicate headers
    """
    seen_hashes = {}  # hash -> (header, sequence)
    duplicates = {}   # hash -> [list of duplicate headers]
    unique = []

    for header, seq in sequences:
        seq_hash = hash_sequence(seq)

        if seq_hash not in seen_hashes:
            seen_hashes[seq_hash] = (header, seq)
            unique.append((header, seq))
            duplicates[seq_hash] = [header]
        else:
            duplicates[seq_hash].append(header)

    return unique, duplicates


def write_fasta(sequences: List[Tuple[str, str]], output_path: Path, line_width: int = 80):
    """Write sequences to FASTA file."""
    with open(output_path, 'w') as f:
        for header, seq in sequences:
            f.write(f">{header}\n")
            # Write sequence in chunks
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')


def check_duplicates(fasta_path: Path) -> dict:
    """Check for duplicates and return statistics."""
    sequences = read_fasta(fasta_path)
    unique, duplicates = deduplicate_sequences(sequences)

    total = len(sequences)
    unique_count = len(unique)
    dup_count = total - unique_count

    # Find most duplicated sequences
    most_duplicated = sorted(
        [(h, len(headers)) for h, headers in duplicates.items()],
        key=lambda x: -x[1]
    )[:10]

    return {
        'path': str(fasta_path),
        'total_sequences': total,
        'unique_sequences': unique_count,
        'duplicate_sequences': dup_count,
        'duplication_ratio': total / unique_count if unique_count > 0 else 0,
        'most_duplicated': [(duplicates[h][0], count) for h, count in most_duplicated if count > 1]
    }


def process_directory(dir_path: Path, output_dir: Path = None, check_only: bool = False):
    """Process all FASTA files in a directory."""
    fasta_files = list(dir_path.glob("*.fasta")) + list(dir_path.glob("*.fa"))

    if not fasta_files:
        print(f"No FASTA files found in {dir_path}")
        return

    total_input = 0
    total_unique = 0

    for fasta_file in sorted(fasta_files):
        stats = check_duplicates(fasta_file)
        total_input += stats['total_sequences']
        total_unique += stats['unique_sequences']

        if stats['duplicate_sequences'] > 0:
            print(f"{fasta_file.name}: {stats['total_sequences']} -> {stats['unique_sequences']} "
                  f"({stats['duplication_ratio']:.1f}x duplication)")

            if not check_only and output_dir:
                sequences = read_fasta(fasta_file)
                unique, _ = deduplicate_sequences(sequences)
                output_file = output_dir / fasta_file.name
                write_fasta(unique, output_file)

    print(f"\nTotal: {total_input} sequences -> {total_unique} unique "
          f"({total_input/total_unique:.1f}x overall duplication)")


def main():
    parser = argparse.ArgumentParser(description="FASTA deduplication utility")
    parser.add_argument("input", nargs="?", help="Input FASTA file")
    parser.add_argument("output", nargs="?", help="Output FASTA file (deduplicated)")
    parser.add_argument("--check", action="store_true", help="Just check for duplicates, don't write")
    parser.add_argument("--dir", type=str, help="Process all FASTAs in directory")
    parser.add_argument("--output-dir", type=str, help="Output directory for deduplicated files")
    args = parser.parse_args()

    if args.dir:
        dir_path = Path(args.dir)
        output_dir = Path(args.output_dir) if args.output_dir else None
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
        process_directory(dir_path, output_dir, check_only=args.check or not output_dir)

    elif args.input:
        input_path = Path(args.input)

        if not input_path.exists():
            print(f"ERROR: File not found: {input_path}")
            sys.exit(1)

        stats = check_duplicates(input_path)

        print(f"Input:     {stats['total_sequences']} sequences")
        print(f"Unique:    {stats['unique_sequences']} sequences")
        print(f"Duplicate: {stats['duplicate_sequences']} sequences")
        print(f"Ratio:     {stats['duplication_ratio']:.2f}x")

        if stats['most_duplicated']:
            print(f"\nMost duplicated sequences:")
            for header, count in stats['most_duplicated'][:5]:
                print(f"  {count}x: {header[:60]}...")

        if not args.check and args.output:
            output_path = Path(args.output)
            sequences = read_fasta(input_path)
            unique, _ = deduplicate_sequences(sequences)
            write_fasta(unique, output_path)
            print(f"\nWrote {len(unique)} unique sequences to {output_path}")

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
