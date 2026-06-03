#!/usr/bin/env python3
"""
FASTA Deduplication with ID Mapping

Creates:
1. A deduplicated FASTA file with unique sequences
2. A mapping file (TSV) that tracks: seq_hash -> [list of original IDs]

After prediction on deduplicated sequences, use expand_results.py to
map scores back to all original variant IDs.

Usage:
    # Deduplicate
    python deduplicate_with_mapping.py input.fasta output_dedup.fasta mapping.tsv

    # Check status
    python deduplicate_with_mapping.py --check input.fasta

Author: Kevin R. Roy
"""

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


def hash_sequence(seq: str) -> str:
    """Create MD5 hash of sequence content (case-insensitive)."""
    return hashlib.md5(seq.upper().strip().encode()).hexdigest()


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


def deduplicate_with_mapping(
    sequences: List[Tuple[str, str]]
) -> Tuple[List[Tuple[str, str, str]], Dict[str, List[str]]]:
    """
    Deduplicate sequences and create ID mapping.

    Returns:
        - List of (canonical_id, sequence, seq_hash) for unique sequences
        - Dict mapping seq_hash -> [list of all IDs with that sequence]
    """
    hash_to_ids = defaultdict(list)  # hash -> [id1, id2, ...]
    hash_to_seq = {}                  # hash -> sequence
    hash_to_canonical = {}            # hash -> first ID (canonical)

    for header, seq in sequences:
        # Extract ID (first field before any whitespace)
        seq_id = header.split()[0]
        seq_hash = hash_sequence(seq)

        hash_to_ids[seq_hash].append(seq_id)

        if seq_hash not in hash_to_seq:
            hash_to_seq[seq_hash] = seq
            hash_to_canonical[seq_hash] = seq_id

    # Build unique sequences list
    unique_seqs = []
    for seq_hash, seq in hash_to_seq.items():
        canonical_id = hash_to_canonical[seq_hash]
        unique_seqs.append((canonical_id, seq, seq_hash))

    return unique_seqs, dict(hash_to_ids)


def write_fasta(sequences: List[Tuple[str, str, str]], output_path: Path, line_width: int = 80):
    """Write deduplicated sequences to FASTA, with hash in header."""
    with open(output_path, 'w') as f:
        for seq_id, seq, seq_hash in sequences:
            # Include hash in header for easy lookup
            f.write(f">{seq_id} hash={seq_hash}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')


def write_mapping(hash_to_ids: Dict[str, List[str]], output_path: Path):
    """
    Write mapping file (TSV format):
    seq_hash<TAB>canonical_id<TAB>all_ids_json
    """
    with open(output_path, 'w') as f:
        f.write("seq_hash\tcanonical_id\tn_duplicates\tall_ids\n")
        for seq_hash, ids in sorted(hash_to_ids.items()):
            canonical = ids[0]
            all_ids_json = json.dumps(ids)
            f.write(f"{seq_hash}\t{canonical}\t{len(ids)}\t{all_ids_json}\n")


def load_mapping(mapping_path: Path) -> Dict[str, List[str]]:
    """Load mapping file back into dict."""
    hash_to_ids = {}
    with open(mapping_path) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            seq_hash = parts[0]
            all_ids = json.loads(parts[3])
            hash_to_ids[seq_hash] = all_ids
    return hash_to_ids


def check_duplicates(fasta_path: Path):
    """Check and report duplication statistics."""
    sequences = read_fasta(fasta_path)
    unique_seqs, hash_to_ids = deduplicate_with_mapping(sequences)

    total = len(sequences)
    unique = len(unique_seqs)

    print(f"Input file:  {fasta_path}")
    print(f"Total IDs:   {total}")
    print(f"Unique seqs: {unique}")
    print(f"Duplicates:  {total - unique}")
    print(f"Ratio:       {total/unique:.2f}x")
    print()

    # Show most duplicated
    by_count = sorted(hash_to_ids.items(), key=lambda x: -len(x[1]))[:10]
    duplicated = [(h, ids) for h, ids in by_count if len(ids) > 1]

    if duplicated:
        print("Most duplicated sequences:")
        for seq_hash, ids in duplicated[:5]:
            print(f"  {len(ids)}x: {ids[0][:50]}... (+ {len(ids)-1} others)")
            if len(ids) <= 5:
                for other_id in ids[1:]:
                    print(f"       -> {other_id[:60]}")


def main():
    parser = argparse.ArgumentParser(
        description="Deduplicate FASTA with ID mapping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Deduplicate and create mapping
    python deduplicate_with_mapping.py input.fasta output.fasta mapping.tsv

    # Check duplication level only
    python deduplicate_with_mapping.py --check input.fasta

    # Process directory
    python deduplicate_with_mapping.py --dir fasta_dir/ --output-dir dedup_dir/
"""
    )
    parser.add_argument("input", nargs="?", help="Input FASTA file")
    parser.add_argument("output", nargs="?", help="Output deduplicated FASTA")
    parser.add_argument("mapping", nargs="?", help="Output mapping TSV file")
    parser.add_argument("--check", action="store_true", help="Check duplicates only")
    parser.add_argument("--dir", type=str, help="Process all FASTAs in directory")
    parser.add_argument("--output-dir", type=str, help="Output directory")
    args = parser.parse_args()

    if args.check and args.input:
        check_duplicates(Path(args.input))
        return

    if args.dir:
        input_dir = Path(args.dir)
        output_dir = Path(args.output_dir) if args.output_dir else input_dir / "deduplicated"
        output_dir.mkdir(parents=True, exist_ok=True)

        fasta_files = list(input_dir.glob("*.fasta")) + list(input_dir.glob("*.fa"))

        total_input = 0
        total_unique = 0

        for fasta_file in sorted(fasta_files):
            sequences = read_fasta(fasta_file)
            unique_seqs, hash_to_ids = deduplicate_with_mapping(sequences)

            total_input += len(sequences)
            total_unique += len(unique_seqs)

            # Write outputs
            out_fasta = output_dir / fasta_file.name
            out_mapping = output_dir / (fasta_file.stem + "_mapping.tsv")

            write_fasta(unique_seqs, out_fasta)
            write_mapping(hash_to_ids, out_mapping)

            if len(sequences) > len(unique_seqs):
                print(f"{fasta_file.name}: {len(sequences)} -> {len(unique_seqs)} "
                      f"({len(sequences)/len(unique_seqs):.1f}x)")

        print(f"\nTotal: {total_input} -> {total_unique} ({total_input/total_unique:.1f}x)")
        print(f"Output: {output_dir}")

    elif args.input and args.output and args.mapping:
        input_path = Path(args.input)
        output_path = Path(args.output)
        mapping_path = Path(args.mapping)

        print(f"Reading {input_path}...")
        sequences = read_fasta(input_path)

        print(f"Deduplicating {len(sequences)} sequences...")
        unique_seqs, hash_to_ids = deduplicate_with_mapping(sequences)

        print(f"Writing {len(unique_seqs)} unique sequences to {output_path}...")
        write_fasta(unique_seqs, output_path)

        print(f"Writing mapping to {mapping_path}...")
        write_mapping(hash_to_ids, mapping_path)

        print(f"\nDone! Reduced {len(sequences)} -> {len(unique_seqs)} "
              f"({len(sequences)/len(unique_seqs):.1f}x)")

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
