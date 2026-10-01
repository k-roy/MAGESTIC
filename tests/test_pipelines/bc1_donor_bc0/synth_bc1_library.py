"""Synthetic bc1-donor-bc0 library with PLANTED TRUTH, for a known-answer round trip.

Why a generator and not a code read: on 2026-10-01 a defect asserted from reading
`domain_insertion/validate.py` was falsified by a test. The evidence for a pipeline's
correctness is recovery against planted truth.

THE TEMPLATE IS ANCHORED IN REAL DATA, not inferred from the parser. Measured
2026-10-01 against
`projects/QTL/20240806_yL406_SpG_QTL_screen/.../mapped_to_oligos/bc1_mapped_to_oligos.tsv`
and the harmonized oligo table:

  * a real parsed `donor` field is 238 nt = [98 nt constant pad][GTTTGAAGAGC][129 nt design donor]
  * `trim_to_variable_region(real parsed donor) == harmonized donor`  -> True
  * `donor_bc0_fragment == donor[-30:] + SPS(20) + bc0(10)`, 60 nt -> matches byte for byte
  * the real 98-nt pad is PAD98 below (plasmid scaffold between the bc1 block and the donor)
  * the bc1-donor-bc0 amplicon does NOT carry the guide. The guide reaches the bc1 only
    through the gdb (plasmid) bridge on `donor_bc0_fragment`. That is why guide-redundant
    oligos are the interesting case.

Amplicon, 400 nt, cut into R1 = amp[:300] and R2 = revcomp(amp)[:300]:

    [PS1 20][CAGGTCATGCTC][bc1 20][tail 10][CGTGCATAAC][MID 28]
    [PAD98][GTTTGAAGAGC][donor 129][SPS 20][bc0 10][AGGAAAACAGAC][flank 20]

With those lengths `revcomp(amp)[:300]` is exactly revcomp(amp[100:]), so stage 01 parses
a 238-nt donor and the anchor lands at position 98 -- `DEFAULT_ANCHOR_POS`, as in real data.
"""

from __future__ import annotations

import gzip
import json
import random
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# --- constants measured from real data / the canonical extractor -------------------
UPSTREAM_ANCHOR = "CAGGTCATGCTC"
UNIVERSAL_DOWNSTREAM = "CGTGCATAAC"
#: layout C (MAGESTIC 3.0 / yL406 era) places the universal anchor at offset +10
TAIL_C = "CCCTAGGATA"          # 10 nt -> offset +10
TAIL_A = "CCCTAGGATACCGTAGGCT"[:18]   # 18 nt -> offset +18 (layout A, 1.0/2.5)
TAIL_B = "C"                   # 1 nt  -> offset +1  (layout B/D, no spacer)
PS1_PAD = "GCTAGCTTCAGGCATTGACC"      # 20 nt
MID = "ACGTTGCATCGGATCCAGTTACGCATGA"  # 28 nt
#: the real 98-nt constant pad upstream of the donor anchor (yL406 2x300, verified)
PAD98 = ("ATCCCCTATTTCATAGGGCACTATCCATCCTTACCCACTCCTCCGTGGGTTTCTTCCACTC"
         "AGGTGCGCACCCTTACGTCAGAAGAAACGGGTTTCCT")
DONOR_ANCHOR = "GTTTGAAGAGC"
MSD = "AGGAAAACAGAC"
FLANK3 = "GGATCCTTAAGGCTCGAATT"   # 20 nt
DONOR_LEN = 129
#: the real V545 subpool priming site
SPS_V545 = "CTTGACTCGACAGTGAGAGC"
SPS_V429 = "CGATAGACGACTGGACAGCA"

_COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _rand_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


# ---------------------------------------------------------------------------------
# The design
# ---------------------------------------------------------------------------------

@dataclass
class Oligo:
    oligo_name: str
    guide: str
    donor: str                 # the 129-nt variable region == the harmonized `donor` column
    bc0: str
    sps: str
    individual_variants: str
    library: str
    locus: str
    nuclease_pam: str = "SpG_NGNG"
    codon_variant_type: str = "missense"

    @property
    def fragment(self) -> str:
        """donor[-30:] + SPS + bc0, 60 nt -- the gdb<->bdb bridge key."""
        return self.donor[-30:] + self.sps + self.bc0

    @property
    def oligo_sequence(self) -> str:
        return "ctgcgattggcaggcgcgcc" + self.guide + "gtttgaagagc" + self.donor + self.sps.lower()


@dataclass
class Design:
    oligos: List[Oligo]
    #: bc1 -> list of oligo_names it was planted against (len>1 only for bc1 collisions)
    links: Dict[str, List[str]] = field(default_factory=dict)
    notes: Dict[str, str] = field(default_factory=dict)

    def by_name(self, name: str) -> Optional[Oligo]:
        for o in self.oligos:
            if o.oligo_name == name:
                return o
        return None

    def variant_of(self, name: str) -> Optional[str]:
        o = self.by_name(name)
        return o.individual_variants if o else None


def make_design(
    seed: int = 20261001,
    n_loci: int = 8,
    variants_per_locus: int = 3,
    guides_per_variant: int = 2,
    library: str = "V545",
    sps: str = SPS_V545,
) -> Design:
    """A realistic design: loci with several SNP variants sharing flanks, each variant
    reachable by several guides (so guide-redundant oligos exist, which is the only place
    the gdb bridge can change an answer), one unique bc0 per oligo."""
    rng = random.Random(seed)
    oligos: List[Oligo] = []
    used_bc0: set = set()
    used_guide: set = set()
    for li in range(n_loci):
        locus = f"chr_syn_{li}_YSY{li:03d}W"      # underscores on purpose
        flank_l = _rand_seq(rng, 64)
        flank_r = _rand_seq(rng, 64)
        base = flank_l + "A" + flank_r            # 129 nt, edit at the centre (index 64)
        assert len(base) == DONOR_LEN
        for vi in range(variants_per_locus):
            alt = "CGT"[vi % 3]
            donor = flank_l + alt + flank_r
            pos = 100000 + li * 1000 + vi
            variant = f"{pos}_A_{alt}_SNP"
            for gi in range(guides_per_variant):
                while True:
                    g = _rand_seq(rng, 20)
                    if g not in used_guide:
                        used_guide.add(g)
                        break
                while True:
                    b = _rand_seq(rng, 10)
                    if b not in used_bc0:
                        used_bc0.add(b)
                        break
                oligos.append(Oligo(
                    oligo_name=f"SpG_Cas9_{locus}_{pos}_A_{alt}_g{gi}",
                    guide=g, donor=donor, bc0=b, sps=sps,
                    individual_variants=variant, library=library, locus=locus,
                ))
    return Design(oligos=oligos)


def plant_bc1s(design: Design, seed: int = 7, min_hd: int = 4) -> Dict[str, List[str]]:
    """One bc1 per oligo, pairwise Hamming >= min_hd so 02's --bc1-edit-distance 1 cannot
    legitimately merge any planted pair. Deliberate HD-1 / HD-2 pairs are added separately."""
    rng = random.Random(seed)
    chosen: List[str] = []
    links: Dict[str, List[str]] = {}
    for o in design.oligos:
        for _ in range(20000):
            b = _rand_seq(rng, 20)
            if all(_hamming(b, c) >= min_hd for c in chosen):
                chosen.append(b)
                links[b] = [o.oligo_name]
                break
        else:
            raise RuntimeError("could not draw a well-separated bc1")
    design.links = links
    return links


# ---------------------------------------------------------------------------------
# Molecules -> reads
# ---------------------------------------------------------------------------------

def build_amplicon(bc1: str, oligo: Oligo, tail: str = TAIL_C,
                   donor_override: Optional[str] = None) -> str:
    donor = donor_override if donor_override is not None else oligo.donor
    return (PS1_PAD + UPSTREAM_ANCHOR + bc1 + tail + UNIVERSAL_DOWNSTREAM + MID
            + PAD98 + DONOR_ANCHOR + donor + oligo.sps + oligo.bc0 + MSD + FLANK3)


def cut_reads(amp: str, l1: int = 300, l2: int = 300) -> Tuple[str, str]:
    """R1 from the 5' end, R2 from the 3' end (reverse complement) -- as a sequencer does."""
    return amp[:l1], rc(amp)[:l2]


def _mutate(seq: str, n: int, rng: random.Random, lo: int = 0, hi: Optional[int] = None) -> str:
    hi = len(seq) if hi is None else hi
    s = list(seq)
    for i in rng.sample(range(lo, hi), n):
        s[i] = rng.choice([c for c in "ACGT" if c != s[i]])
    return "".join(s)


@dataclass
class Molecule:
    bc1: str
    oligo: Oligo
    depth: int = 5
    tail: str = TAIL_C
    donor_override: Optional[str] = None
    #: substitutions applied to the AMPLICON (synthesis error -- both reads see it)
    synth_mismatches: int = 0
    #: nt deleted from the donor region of the amplicon (synthesis indel/truncation)
    synth_deletion: int = 0
    #: substitutions applied to R1 only, inside the bc1 (sequencing error)
    bc1_seq_errors: int = 0
    bc1_len_override: Optional[int] = None


def emit_reads(mols: List[Molecule], seed: int = 11) -> List[Tuple[str, str, str]]:
    """-> list of (read_name, r1, r2). Order is shuffled: `detect_orientation` and PASS-1
    sample only the head of the file, so a sorted file would hide a layout problem."""
    rng = random.Random(seed)
    out: List[Tuple[str, str, str]] = []
    i = 0
    for m in mols:
        bc1 = m.bc1
        if m.bc1_len_override is not None and m.bc1_len_override != len(bc1):
            bc1 = (bc1[:m.bc1_len_override] if m.bc1_len_override < len(bc1)
                   else bc1 + _rand_seq(rng, m.bc1_len_override - len(bc1)))
        donor = m.donor_override if m.donor_override is not None else m.oligo.donor
        if m.synth_deletion:
            donor = donor[:-m.synth_deletion]
        amp = build_amplicon(bc1, m.oligo, tail=m.tail, donor_override=donor)
        if m.synth_mismatches:
            # confine synthesis errors to the donor variable region
            d0 = amp.find(DONOR_ANCHOR) + len(DONOR_ANCHOR)
            amp = _mutate(amp, m.synth_mismatches, rng, d0, d0 + len(donor))
        for _ in range(m.depth):
            r1, r2 = cut_reads(amp)
            if m.bc1_seq_errors:
                b0 = len(PS1_PAD) + len(UPSTREAM_ANCHOR)
                r1 = _mutate(r1, m.bc1_seq_errors, rng, b0, b0 + len(bc1))
            i += 1
            out.append((f"synth_{i}", r1, r2))
    rng.shuffle(out)
    return out


def write_fastq_pair(reads: List[Tuple[str, str, str]], r1_path: Path, r2_path: Path,
                     gzipped: bool = False, header_style: str = "plain") -> None:
    """header_style: 'plain' (@name) | 'int' (@12) | 'count' (@seq_1_count=N) -- the two
    collapsed-fastq producer formats that the bc1 caller is documented to confuse."""
    op = (lambda p: gzip.open(p, "wt")) if gzipped else (lambda p: open(p, "w"))
    with op(r1_path) as f1, op(r2_path) as f2:
        for n, (name, r1, r2) in enumerate(reads, start=1):
            if header_style == "int":
                h = str(n)
            elif header_style == "count":
                h = f"seq_{n}_count=1"
            else:
                h = name
            f1.write(f"@{h}\n{r1}\n+\n{'I' * len(r1)}\n")
            f2.write(f"@{h}\n{r2}\n+\n{'I' * len(r2)}\n")


# ---------------------------------------------------------------------------------
# Reference tables the stages need
# ---------------------------------------------------------------------------------

HARMONIZED_NAME = ("20210422_and_20240411_Bloom_et_al_16_strains_QTL_harmonized_"
                   "designed_variant_oligos.tsv")
HARMONIZED_COLS = [
    "oligo_name", "oligo_sequence", "guide", "donor", "guide_nuclease_PAM", "SPS",
    "guide_donor_bc0_plasmid_library", "codon_variant_type",
    "systematic_gene_name_locus", "common_gene_name_locus", "individual_variants",
]


def write_common_dir(design: Design, common: Path) -> Path:
    """Plant the design as the harmonized oligo table under a synthetic MAGESTIC_COMMON_DIR.
    `path_utils` binds HARMONIZED_VARIANTS_FILE at import, so stages must run as
    subprocesses with the env var set."""
    (common / "annotation_files").mkdir(parents=True, exist_ok=True)
    p = common / "annotation_files" / HARMONIZED_NAME
    with open(p, "w") as f:
        f.write("\t".join(HARMONIZED_COLS) + "\n")
        for o in design.oligos:
            f.write("\t".join([
                o.oligo_name, o.oligo_sequence, o.guide, o.donor, o.nuclease_pam, o.sps,
                o.library, o.codon_variant_type, o.locus, o.locus, o.individual_variants,
            ]) + "\n")
    return p


#: the REAL gdb counts schema, measured 2026-10-01 from
#: 20210822_yL274-yL284_QTL_screen/.../V437_guide_donor_bc0_counts.tsv
GDB_REAL_COLS = ["library_ID", "total_counts", "total_fwd_counts", "total_rev_counts",
                 "leader", "guide", "donor", "SPS", "bc0", "msd", "donor_bc0_fragment", "seq"]


def write_gdb_dir(design: Design, gdb_dir: Path, *,
                  reads_per_fragment: int = 3,
                  count_col: str = "total_counts",
                  ambiguous: Optional[Dict[str, str]] = None,
                  read_level: bool = True,
                  omit_fragment_column: bool = False) -> Path:
    """A read-level gdb table in the REAL column schema.

    `ambiguous` maps fragment -> a wrong oligo_name whose row is written FIRST in file
    order with LOWER support, so "most-supported wins" and "first in file wins" give
    different answers. That is the only way the collapse on/off comparison is not vacuous.
    """
    gdb_dir.mkdir(parents=True, exist_ok=True)
    ambiguous = ambiguous or {}
    p = gdb_dir / f"{design.oligos[0].library}_guide_donor_bc0_counts.tsv"
    cols = list(GDB_REAL_COLS)
    if count_col not in cols:
        cols.insert(1, count_col)
    cols.append("oligo_name")
    if omit_fragment_column:
        cols.remove("donor_bc0_fragment")

    def row(o: Oligo, n: int, oligo_name: str) -> List[str]:
        d = {"library_ID": o.library, "total_counts": str(n), "total_fwd_counts": str(n),
             "total_rev_counts": "0", "leader": "NA", "guide": o.guide, "donor": o.donor,
             "SPS": o.sps, "bc0": o.bc0, "msd": MSD,
             "donor_bc0_fragment": o.fragment, "seq": "NA", "oligo_name": oligo_name,
             count_col: str(n)}
        return [d.get(c, "NA") for c in cols]

    with open(p, "w") as f:
        f.write("\t".join(cols) + "\n")
        for o in design.oligos:
            wrong = ambiguous.get(o.fragment)
            if wrong:
                # written FIRST, with LOWER support: distinguishes the two tie-break rules
                f.write("\t".join(row(o, 1, wrong)) + "\n")
            n_rows = reads_per_fragment if read_level else 1
            for k in range(n_rows):
                f.write("\t".join(row(o, 10 + k, o.oligo_name)) + "\n")
    return p


def write_sps_mapping(design: Design, path: Path) -> Path:
    """Satisfies BOTH loaders: core.matching.load_sps_mapping wants columns
    (sps, library_id); 05_assign_sublibrary has its own reader."""
    seen = {}
    for o in design.oligos:
        seen[o.sps] = o.library
    with open(path, "w") as f:
        f.write("sps\tlibrary_id\toligo_pool\tdescription\n")
        for s, lib in seen.items():
            f.write(f"{s}\t{lib}\tsynth_pool\tsynthetic round-trip subpool\n")
    return path


MANIFEST02_COLS = ["sample_id", "include", "exclude_reason", "keyfile_note", "raw_r1",
                   "raw_r2", "r1_bytes", "parsed_path", "sequencing_date", "run_id",
                   "sample_number", "outer_primers", "gDNA_plate_name", "inner_primers",
                   "plate_polymerase", "condition", "amplicon_level", "yeast_library",
                   "source_rule"]


def write_manifest02(samples: List[dict], path: Path) -> Path:
    """The real 02 manifest schema, measured from
    20210822_.../bc1_donor_bc0_REBUILD_20260910/manifest_02.tsv"""
    with open(path, "w") as f:
        f.write("\t".join(MANIFEST02_COLS) + "\n")
        for s in samples:
            f.write("\t".join(str(s.get(c, "")) for c in MANIFEST02_COLS) + "\n")
    return path


def write_truth(design: Design, path: Path) -> Path:
    truth = {
        "links": design.links,
        "oligos": {o.oligo_name: asdict(o) for o in design.oligos},
        "fragment_to_oligo": {o.fragment: o.oligo_name for o in design.oligos},
        "notes": design.notes,
    }
    path.write_text(json.dumps(truth, indent=1))
    return path
