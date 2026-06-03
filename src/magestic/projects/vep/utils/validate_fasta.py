#!/usr/bin/env python3
"""
FASTA Validation Script for VEP Pipeline

Validates FASTA files before submitting VEP jobs.
Run this BEFORE submitting to GPU queue to catch issues early.

Usage:
    python validate_fasta.py <fasta_file> [fasta_file2 ...]

Author: Kevin R. Roy
"""

import sys
from pathlib import Path

def validate_fasta(fasta_path: str) -> bool:
    """Validate a FASTA file. Returns True if valid."""
    path = Path(fasta_path)
    
    if not path.exists():
        print(f"ERROR: File does not exist: {fasta_path}")
        return False
    
    if path.stat().st_size == 0:
        print(f"ERROR: File is empty: {fasta_path}")
        return False
    
    seq_count = 0
    errors = []
    
    with open(path) as f:
        first_line = f.readline().strip()
        if not first_line.startswith('>'):
            print(f"ERROR: First line must start with '>' (header)")
            print(f"  Found: {first_line[:50]}...")
            return False
        
        seq_count = 1
        f.seek(0)
        
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line.startswith('>'):
                seq_count += 1
    
    seq_count = seq_count // 2 + 1  # Rough count
    
    # Recount properly
    with open(path) as f:
        seq_count = sum(1 for line in f if line.startswith('>'))
    
    print(f"VALID: {path.name} ({seq_count} sequences)")
    return True

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: validate_fasta.py <fasta_file> [...]")
        sys.exit(1)
    
    all_valid = all(validate_fasta(f) for f in sys.argv[1:])
    sys.exit(0 if all_valid else 1)
