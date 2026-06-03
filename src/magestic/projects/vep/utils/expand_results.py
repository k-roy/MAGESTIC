#!/usr/bin/env python3
"""
Expand Deduplicated Results Back to All Original IDs

After running VEP predictions on deduplicated sequences, use this script
to expand results back to all original variant IDs using the mapping file.

Usage:
    python expand_results.py results.tsv mapping.tsv expanded_results.tsv

The results file must have a column that matches the canonical_id from deduplication.
Common column names: test_id, seq_id, sequence_id, id

Author: Kevin R. Roy
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd


def load_mapping(mapping_path: Path) -> Dict[str, List[str]]:
    """Load mapping file: canonical_id -> [all original IDs]."""
    canonical_to_all = {}

    with open(mapping_path) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            canonical_id = parts[1]
            all_ids = json.loads(parts[3])
            canonical_to_all[canonical_id] = all_ids

    return canonical_to_all


def expand_results(
    results_df: pd.DataFrame,
    mapping: Dict[str, List[str]],
    id_column: str = "test_id"
) -> pd.DataFrame:
    """
    Expand results to all original IDs.

    For each row with a canonical_id, create duplicate rows for all
    other IDs that had the same sequence.
    """
    expanded_rows = []

    for _, row in results_df.iterrows():
        canonical_id = row[id_column]

        if canonical_id in mapping:
            all_ids = mapping[canonical_id]
            for original_id in all_ids:
                new_row = row.copy()
                new_row[id_column] = original_id
                expanded_rows.append(new_row)
        else:
            # ID not in mapping, keep as-is
            expanded_rows.append(row)

    return pd.DataFrame(expanded_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Expand deduplicated results back to all original IDs"
    )
    parser.add_argument("results", help="Results TSV file from VEP predictions")
    parser.add_argument("mapping", help="Mapping TSV file from deduplication")
    parser.add_argument("output", help="Output expanded results TSV")
    parser.add_argument(
        "--id-column", default="test_id",
        help="Column name containing sequence IDs (default: test_id)"
    )
    args = parser.parse_args()

    results_path = Path(args.results)
    mapping_path = Path(args.mapping)
    output_path = Path(args.output)

    print(f"Loading mapping from {mapping_path}...")
    mapping = load_mapping(mapping_path)
    print(f"  {len(mapping)} canonical IDs")

    print(f"Loading results from {results_path}...")
    results_df = pd.read_csv(results_path, sep='\t')
    print(f"  {len(results_df)} rows")

    if args.id_column not in results_df.columns:
        print(f"ERROR: Column '{args.id_column}' not found in results.")
        print(f"Available columns: {list(results_df.columns)}")
        sys.exit(1)

    print(f"Expanding results using '{args.id_column}' column...")
    expanded_df = expand_results(results_df, mapping, args.id_column)
    print(f"  {len(expanded_df)} rows after expansion")

    print(f"Writing to {output_path}...")
    expanded_df.to_csv(output_path, sep='\t', index=False)

    expansion_ratio = len(expanded_df) / len(results_df)
    print(f"\nDone! Expanded {len(results_df)} -> {len(expanded_df)} ({expansion_ratio:.1f}x)")


if __name__ == "__main__":
    main()
