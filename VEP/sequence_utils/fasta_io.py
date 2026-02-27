"""
FASTA file I/O utilities.

Functions:
- read_fasta: Read FASTA file into dictionary
- write_fasta: Write sequences to FASTA file
- iter_fasta: Memory-efficient iterator over FASTA entries

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Dict, Iterator, Tuple, Union


def read_fasta(fasta_file: Union[str, Path]) -> Dict[str, str]:
    """
    Read FASTA file into dictionary.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        {header: sequence} for all entries
    """
    fasta_file = Path(fasta_file)
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]  # First word after >
                current_seq = []
            else:
                current_seq.append(line)

        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def write_fasta(sequences: Dict[str, str], output_file: Union[str, Path],
                line_width: int = 80) -> None:
    """
    Write sequences to FASTA file.

    Args:
        sequences: {header: sequence}
        output_file: Output path
        line_width: Characters per line (default 80)
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')


def iter_fasta(fasta_file: Union[str, Path]) -> Iterator[Tuple[str, str]]:
    """
    Iterate over FASTA entries (memory-efficient for large files).

    Args:
        fasta_file: Path to FASTA file

    Yields:
        (header, sequence) tuples
    """
    fasta_file = Path(fasta_file)
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    yield current_header, ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_header:
            yield current_header, ''.join(current_seq)


if __name__ == "__main__":
    import tempfile
    import os

    # Create test FASTA content
    test_content = """>seq1
ATGCATGC
ATGCATGC
>seq2 extra description
GGGGCCCC
"""

    # Write to temp file and test reading
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(test_content)
        temp_path = f.name

    try:
        # Test read_fasta
        seqs = read_fasta(temp_path)
        assert len(seqs) == 2
        assert seqs['seq1'] == "ATGCATGCATGCATGC"
        assert seqs['seq2'] == "GGGGCCCC"
        print("read_fasta: PASS")

        # Test iter_fasta
        seqs_iter = list(iter_fasta(temp_path))
        assert len(seqs_iter) == 2
        assert seqs_iter[0] == ('seq1', 'ATGCATGCATGCATGC')
        print("iter_fasta: PASS")

        # Test write_fasta
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            out_path = f.name

        write_fasta({'test': 'ATGC' * 30}, out_path, line_width=80)
        written = read_fasta(out_path)
        assert written['test'] == 'ATGC' * 30
        print("write_fasta: PASS")

        os.unlink(out_path)

    finally:
        os.unlink(temp_path)

    print("\nAll fasta_io.py tests passed!")
