"""
Genome / gene-model layer for domain-insertion design.

Design consideration 5 in one sentence: **the donor must match the GENOME (it is a homologous
repair template), but codons must be read from the SPLICED CDS.** Every coordinate helper here
exists to keep those two frames straight. 55 of the 1,065 library genes have introns, including
NBL1 (YHR199C-A), one of the nine genes added in artifact 02.

Coordinates are 1-based inclusive throughout, matching GFF and GenBank.
"""

from dataclasses import dataclass
from functools import cached_property
import gzip
import re

_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv",
                      "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII",
         "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]


@dataclass
class GeneModel:
    """A protein-coding gene: its chromosome, strand, and CDS exons."""

    orf: str
    chrom: str
    strand: str
    exons: tuple          # in TRANSCRIPT order; each (start, end) with start <= end genomic

    @cached_property
    def cds_length(self) -> int:
        return sum(e - s + 1 for s, e in self.exons)

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    @property
    def has_intron(self) -> bool:
        return len(self.exons) > 1

    def cds_to_genomic(self, pos: int) -> int:
        """1-based CDS nucleotide -> genomic coordinate."""
        if not 1 <= pos <= self.cds_length:
            raise IndexError(f"CDS position {pos} outside 1..{self.cds_length} for {self.orf}")
        remaining = pos - 1
        for s, e in self.exons:
            span = e - s + 1
            if remaining < span:
                return s + remaining if self.strand == "+" else e - remaining
            remaining -= span
        raise AssertionError("unreachable")

    def genomic_to_cds(self, coord: int):
        """Genomic coordinate -> 1-based CDS nucleotide, or None if intronic/outside."""
        offset = 0
        for s, e in self.exons:
            if s <= coord <= e:
                return offset + ((coord - s) if self.strand == "+" else (e - coord)) + 1
            offset += e - s + 1
        return None

    def codon_cds_range(self, residue: int) -> tuple:
        """1-based residue -> (first, last) CDS nucleotide positions."""
        if residue < 1:
            raise IndexError(f"residue {residue} < 1")
        return (residue - 1) * 3 + 1, residue * 3

    def codon_genomic_coords(self, residue: int) -> list:
        """The three genomic coordinates of a codon, in transcript order.

        For a codon spanning an intron these are NOT contiguous -- that is the point.
        """
        a, b = self.codon_cds_range(residue)
        return [self.cds_to_genomic(p) for p in range(a, b + 1)]

    def codon_is_split(self, residue: int) -> bool:
        """True if the codon straddles an intron."""
        c = self.codon_genomic_coords(residue)
        step = 1 if self.strand == "+" else -1
        return not (c[1] - c[0] == step and c[2] - c[1] == step)

    def spliced_cds(self, chrom_seq: str) -> str:
        """Concatenated exon sequence in transcript orientation."""
        parts = []
        for s, e in self.exons:
            seg = chrom_seq[s - 1:e]
            parts.append(seg if self.strand == "+" else rc(seg))
        return "".join(parts).upper()

    def genomic_window(self, chrom_seq: str, center: int, flank: int) -> tuple:
        """
        Unspliced genomic sequence around `center`, returned in TRANSCRIPT orientation.

        Returns (sequence, offset_of_center). Used for PAM search and donor-arm construction,
        both of which operate on genomic (intron-containing) sequence.
        """
        lo, hi = center - flank, center + flank
        lo, hi = max(1, lo), min(len(chrom_seq), hi)
        seg = chrom_seq[lo - 1:hi].upper()
        if self.strand == "+":
            return seg, center - lo
        return rc(seg), hi - center

    def window_with_coords(self, chrom_seq: str, center: int, flank: int) -> tuple:
        """
        Like `genomic_window`, but also returns the genomic coordinate of every window index
        and its CDS position (None where intronic or outside the CDS).

        Returns (sequence, center_index, genomic_coords, cds_positions) -- all in transcript
        orientation, so index arithmetic runs left-to-right in gene sense regardless of strand.
        This is what lets the designer recode CODONS while writing GENOMIC donor arms.
        """
        lo, hi = max(1, center - flank), min(len(chrom_seq), center + flank)
        seg = chrom_seq[lo - 1:hi].upper()
        if self.strand == "+":
            W = seg
            coords = list(range(lo, hi + 1))
            center_idx = center - lo
        else:
            W = rc(seg)
            coords = list(range(hi, lo - 1, -1))
            center_idx = hi - center
        cds = [self.genomic_to_cds(c) for c in coords]
        return W, center_idx, coords, cds


_PARENT_RE = re.compile(r"Parent=([^;]+)")


def _orf_from_parent(attr: str):
    m = _PARENT_RE.search(attr)
    if not m:
        return None
    name = m.group(1).split(",")[0]
    name = re.sub(r"_mRNA$", "", name)
    name = re.sub(r"_id\d+$", "", name)
    return name or None


def load_gene_models(gff_path: str) -> dict:
    """Parse CDS features from an SGD GFF into {ORF: GeneModel}."""
    raw = {}
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            orf = _orf_from_parent(f[8])
            if not orf:
                continue
            rec = raw.setdefault(orf, {"chrom": f[0], "strand": f[6], "exons": []})
            rec["exons"].append((int(f[3]), int(f[4])))

    models = {}
    for orf, rec in raw.items():
        exons = sorted(set(rec["exons"]))
        if rec["strand"] == "-":
            exons = exons[::-1]                      # transcript order
        models[orf] = GeneModel(orf=orf, chrom=rec["chrom"], strand=rec["strand"],
                                exons=tuple(exons))
    return models


def load_chromosomes(fasta_path: str) -> dict:
    """
    Load chromosome sequences keyed by SGD-style name ('chrI'...'chrXVI', 'chrmt').

    SGD's reference FASTA uses RefSeq ids with a '[chromosome=IV]' tag in the defline;
    the GFF uses 'chrIV'. This reconciles the two.
    """
    seqs, name, buf = {}, None, []
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                m = re.search(r"\[chromosome=([^\]]+)\]", line)
                if m:
                    name = "chrmt" if m.group(1) == "mito" else "chr" + m.group(1)
                else:
                    m2 = re.search(r"\[location=mitochondrion\]", line)
                    name = "chrmt" if m2 else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf)
    return {k: v.upper() for k, v in seqs.items()}
