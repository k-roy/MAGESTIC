"""
Coding-consequence annotation for guide-donor design.

Guide-donor design depends on knowing what an edit does to a protein: the choice
of a synonymous recode, the PTC arm, the missense arm and the decision of which
variants deserve extra guides all key on the consequence. That annotation
therefore belongs beside the custom-genome build and liftover, not in per-project
scripts.

WHY THIS MODULE EXISTS
----------------------
The ad-hoc annotators used before it only ever emitted a call for *single-base
substitutions*. Everything else was left EMPTY, coding or not. Measured on a
shipped 24,000-oligo library (15,641 natural designs):

    variant class   total    blank consequence
    substitution   13,868      882  (6.4 %, essentially all multi-base)
    deletion          914      914  (100 %)
    insertion         859      859  (100 %)

Of those 2,655 blanks, **365 sit inside a CDS** — 257 substitutions (120 of them
4 nt spanning two codons), 75 frameshift indels and 33 in-frame indels — where a
consequence is fully determined by the sequence and was simply never computed.
The other 2,290 are outside any CDS and should have read ``non_coding`` rather
than nothing. A multi-base MNP in a coding region produces a definite protein
change, and an indel is predictably in-frame or a frameshift, so "the annotator
only handled SNVs" is a defect and not a limitation.

CONTRACT
--------
``annotate_variant`` NEVER returns an empty consequence. Every call lands on one
of ``CONSEQUENCES``. An edit outside the named gene's CDS is ``non_coding``,
which is a positive statement, not a blank.

Consequences are judged by translating the whole spliced CDS before and after the
edit, so multi-codon MNPs, in-frame indels and frameshifts are all handled by the
same path rather than by special cases. For a change that alters several residues
at once, ``aa_changes`` lists all of them and ``aa_change`` carries the single
most severe one.

A position can sit inside two CDSs at once (14 designs in the library above do).
``annotate_overlapping`` returns one result per containing gene; a caller that
needs a scalar must choose, and must not assume there is only one.

Author: Kevin R. Roy
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

logger = logging.getLogger(__name__)

__all__ = [
    "CONSEQUENCES",
    "VariantConsequence",
    "translate",
    "build_cds",
    "cds_index",
    "load_cds_table",
    "load_cds_from_gff",
    "build_gene_index",
    "genes_at",
    "annotate_variant",
    "annotate_overlapping",
]

# ---------------------------------------------------------------- codon table

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE: Dict[str, str] = {
    a + b + c: _AAS[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES
    )
}

#: Every value ``annotate_variant`` can return. There is deliberately no empty
#: string here: a blank consequence is the defect this module was written to fix.
CONSEQUENCES = frozenset(
    {
        "non_coding",
        "synonymous",
        "missense",
        "nonsense",
        "stop_lost",
        "start_lost",
        "frameshift",
        "inframe_insertion",
        "inframe_deletion",
    }
)

# Severity order used to pick the single `aa_change` when a change alters several
# residues. Documented rather than incidental, so callers can rely on it.
_SEVERITY = ("nonsense", "start_lost", "stop_lost", "missense", "synonymous")

_CDS_NAME_RE = re.compile(r"Name=([^;_]+(?:-[A-Z0-9]+)?)_CDS")
_GFF_NAME_RE = re.compile(r"(?:^|;)(?:Name|Parent|ID)=([^;]+)")

#: ``{gene: (chrom, strand, [(start0, end0), ...] in genome order)}``
CdsModel = Dict[str, Tuple[str, str, List[Tuple[int, int]]]]


@dataclass(frozen=True)
class VariantConsequence:
    """What an edit does to one gene's protein.

    ``consequence`` is always set and always a member of :data:`CONSEQUENCES`.
    ``aa_change`` is the single most severe residue change in HGVS-ish short form
    (``W976*``), or ``""`` where no single residue describes the event
    (``synonymous``, ``non_coding``, and frameshifts, which change everything
    downstream). ``aa_changes`` lists every altered residue, which is what a
    multi-codon MNP needs.
    """

    consequence: str
    aa_change: str = ""
    aa_changes: Tuple[str, ...] = ()
    gene: str = ""
    detail: str = ""

    def __post_init__(self) -> None:
        if self.consequence not in CONSEQUENCES:
            raise ValueError(
                f"{self.consequence!r} is not a known consequence; "
                f"expected one of {sorted(CONSEQUENCES)}"
            )


# ---------------------------------------------------------------- primitives


def _rev_comp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(nt: str) -> str:
    """Translate a coding sequence. Unknown or partial codons become ``X``."""
    return "".join(
        CODON_TABLE.get(nt[i : i + 3].upper(), "X") for i in range(0, len(nt) - 2, 3)
    )


def build_cds(seq: str, exons: Sequence[Tuple[int, int]], strand: str) -> str:
    """Splice the exons out of ``seq`` and orient them coding-strand 5'->3'."""
    spliced = "".join(seq[a:b] for a, b in exons).upper()
    return _rev_comp(spliced) if strand == "-" else spliced


def cds_index(
    exons: Sequence[Tuple[int, int]], strand: str, pos0: int
) -> Optional[int]:
    """Offset of genome position ``pos0`` within the spliced CDS, else ``None``.

    ``pos0`` is 0-based. Returns ``None`` when the position is intronic or
    outside the gene, which is what makes ``non_coding`` a measured answer
    rather than a guess.
    """
    run = 0
    total = sum(b - a for a, b in exons)
    for a, b in exons:
        if a <= pos0 < b:
            off = run + (pos0 - a)
            return (total - 1 - off) if strand == "-" else off
        run += b - a
    return None


# ---------------------------------------------------------------- CDS loading


def load_cds_table(path: Union[str, Path]) -> CdsModel:
    """Load a HEADERLESS CDS table: chrom, start(1-based), end, strand, phase, attrs.

    This is the layout emitted by the design-frame CDS exports (for example
    ``310_r64_5_1_cds.tsv`` and ``467ad_designframe_cds.tsv``). It has NO header
    row; reading one with a header-consuming parser silently yields zero
    intervals, which makes every downstream membership test answer "outside a
    CDS". That failure is invisible unless the number of intervals is checked, so
    this function raises rather than returning an empty model.
    """
    out: Dict[str, List] = defaultdict(lambda: [None, None, []])
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            match = _CDS_NAME_RE.search(parts[5])
            if not match:
                continue
            gene = match.group(1)
            out[gene][0] = parts[0]
            out[gene][1] = parts[3]
            out[gene][2].append((int(parts[1]) - 1, int(parts[2])))
    for gene in out:
        out[gene][2].sort()
    model = {g: (v[0], v[1], v[2]) for g, v in out.items()}
    if not model:
        raise ValueError(
            f"{path} yielded no CDS features. This table is HEADERLESS; a parser "
            "that consumes the first row as a header produces exactly this, and "
            "an empty model makes every position read as non_coding."
        )
    logger.info("loaded %d genes from CDS table %s", len(model), path)
    return model


def load_cds_from_gff(path: Union[str, Path]) -> CdsModel:
    """Load CDS features from a GFF3 (column 3 == ``CDS``).

    Genes are keyed by the ``Name=<gene>_CDS`` attribute where present, else by
    ``Parent``/``ID`` with any ``_CDS``/``_mRNA`` suffix stripped.
    """
    out: Dict[str, List] = defaultdict(lambda: [None, None, []])
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            match = _CDS_NAME_RE.search(parts[8]) or _GFF_NAME_RE.search(parts[8])
            if not match:
                continue
            gene = match.group(1)
            for suffix in ("_CDS", "_mRNA"):
                if gene.endswith(suffix):
                    gene = gene[: -len(suffix)]
            out[gene][0] = parts[0]
            out[gene][1] = parts[6]
            out[gene][2].append((int(parts[3]) - 1, int(parts[4])))
    for gene in out:
        out[gene][2].sort()
    model = {g: (v[0], v[1], v[2]) for g, v in out.items()}
    if not model:
        raise ValueError(f"{path} yielded no CDS features")
    logger.info("loaded %d genes from GFF %s", len(model), path)
    return model


def build_gene_index(model: CdsModel) -> Dict[str, List[Tuple[int, int, str]]]:
    """Index of ``{chrom: [(start0, end0, gene), ...]}`` for containment lookups."""
    index: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for gene, (chrom, _strand, exons) in model.items():
        for a, b in exons:
            index[chrom].append((a, b, gene))
    for chrom in index:
        index[chrom].sort()
    return dict(index)


def genes_at(
    index: Dict[str, List[Tuple[int, int, str]]], chrom: str, pos0: int
) -> List[str]:
    """Every gene whose CDS contains ``pos0``. Can legitimately return 2+."""
    return sorted({g for a, b, g in index.get(chrom, ()) if a <= pos0 < b})


# ---------------------------------------------------------------- the classifier


def _residue_changes(prot_ref: str, prot_alt: str) -> List[str]:
    """``A123V``-style descriptions of every residue that differs, 1-based."""
    changes = []
    for i in range(min(len(prot_ref), len(prot_alt))):
        if prot_ref[i] != prot_alt[i]:
            changes.append(f"{prot_ref[i]}{i + 1}{prot_alt[i]}")
    return changes


def _first_stop(prot: str) -> int:
    idx = prot.find("*")
    return idx if idx >= 0 else len(prot)


def _pick(changes: Sequence[str], kind: str) -> str:
    """The representative ``aa_change`` for a consequence of this ``kind``."""
    if not changes:
        return ""
    if kind == "nonsense":
        return next((c for c in changes if c.endswith("*")), changes[0])
    if kind == "stop_lost":
        return next((c for c in changes if c.startswith("*")), changes[0])
    if kind == "start_lost":
        return next((c for c in changes if c[1:-1] == "1"), changes[0])
    return changes[0]


def _classify_substitution(
    prot_ref: str, prot_alt: str
) -> Tuple[str, str, Tuple[str, ...]]:
    """A length-preserving change of any width, including a multi-codon MNP."""
    changes = _residue_changes(prot_ref, prot_alt)
    if not changes and len(prot_ref) == len(prot_alt):
        return "synonymous", "", ()
    if not changes:
        # Same nucleotide length but a different protein length should not happen;
        # report it as missense rather than silently calling it synonymous.
        return "missense", "", ()

    gained_stop = any(c.endswith("*") for c in changes)
    lost_stop = any(c.startswith("*") for c in changes)
    lost_start = any(c[1:-1] == "1" and c.startswith("M") for c in changes)

    # Documented precedence: a gained stop dominates, then a destroyed start,
    # then a destroyed stop, then an ordinary residue change.
    if gained_stop:
        kind = "nonsense"
    elif lost_start:
        kind = "start_lost"
    elif lost_stop:
        kind = "stop_lost"
    else:
        kind = "missense"
    return kind, _pick(changes, kind), tuple(changes)


def _classify_indel(
    prot_ref: str, prot_alt: str, ref: str, alt: str, codon: int
) -> Tuple[str, str, Tuple[str, ...], str]:
    """A length-changing change: frameshift, or an in-frame insertion/deletion."""
    net = len(alt) - len(ref)
    if net % 3:
        return (
            "frameshift",
            "",
            (),
            f"net {net:+d} nt from residue {codon + 1}; every downstream residue changes",
        )

    kind = "inframe_insertion" if net > 0 else "inframe_deletion"
    detail = f"net {net:+d} nt ({abs(net) // 3} residue(s)) at residue {codon + 1}"

    # An in-frame indel can still introduce a premature stop, which matters far
    # more than its length does. The comparison MUST allow for the length the
    # indel itself adds or removes: a clean 3 nt deletion moves the terminal stop
    # one residue earlier, and testing against the unadjusted reference position
    # reports every in-frame deletion as nonsense.
    if _first_stop(prot_alt) < _first_stop(prot_ref) + net // 3:
        return (
            "nonsense",
            f"*{_first_stop(prot_alt) + 1}",
            (),
            f"premature stop introduced by an {kind} ({detail})",
        )
    changes = _residue_changes(prot_ref[:codon], prot_alt[:codon])
    return kind, "", tuple(changes), detail


def shift_exons_for_indel(
    exons: Sequence[Tuple[int, int]], pos0: int, net: int
) -> Optional[List[Tuple[int, int]]]:
    """The exon intervals as they are AFTER an indel of ``net`` nt at ``pos0``.

    🔴 **WITHOUT THIS, AN INSERTION NEVER LENGTHENS THE PROTEIN.** Splicing the edited genome with
    the ORIGINAL exon coordinates keeps the CDS window the same width, so inserted bases displace
    the same number of downstream bases and the terminal stop can fall out of the window entirely.
    Measured on 10 real designs with net +3 to +27 nt: ``len(prot_alt) == len(prot_ref)`` on every
    one, 4 of them lost their stop codon, and 3 landed a "premature" stop at exactly the reference
    stop position, which made the premature-stop test true by construction.

    So every consequence computed for a length-changing allele this way was unsound, and the
    refinement built on top of it reported truncations that are artifacts of the window.

    An interval entirely before the edit is unchanged; the interval containing the edit grows or
    shrinks by ``net``; an interval entirely after it shifts by ``net``.

    Returns ``None`` when the edit SPANS an exon boundary, where this arithmetic does not apply and
    guessing would be worse than declining: the caller must treat that as unresolvable rather than
    receive a plausible wrong answer.
    """
    out: List[Tuple[int, int]] = []
    containing = 0
    for a, b in exons:
        if b <= pos0:
            out.append((a, b))
        elif a <= pos0 < b:
            containing += 1
            if b + net < a:
                return None          # the deletion removes more than the exon holds
            out.append((a, b + net))
        else:
            out.append((a + net, b + net))
    if containing != 1:
        return None                  # at a boundary, or outside every exon
    return out


def classify_translations(
    prot_ref: str,
    prot_alt: str,
    ref: str,
    alt: str,
    codon: Optional[int] = None,
) -> VariantConsequence:
    """The consequence of ``ref``->``alt`` read from the two TRANSLATIONS alone.

    :func:`annotate_variant` is the right entry point when you hold the genome and the CDS model.
    This one exists for the callers that have already spliced and translated both proteins and only
    need the judgement, so that **the judgement has one implementation**.

    🔴 **WHY IT IS PUBLIC.** A project copy of this logic shipped a 24,000-oligo pool in which 7
    designs that destroy a gene's start methionine were labelled ``missense``. Its classifier took
    the FIRST differing residue and had no start-codon case, so ``M1V`` and ``M1I`` read as ordinary
    substitutions. Measured against this module on the same 14,748 coding rows, the two agreed on
    14,741 and the 7 they disagreed on were all codon 1. Destroying the start is not a missense
    change, and the project had no way to say so. One implementation, imported, is the fix.

    ``codon`` is the 0-based residue index of the edit. It is optional because a caller holding only
    the translations does not always know it; it refines an indel's ``detail`` and ``aa_changes``
    and changes no consequence.

    Precedence is the same as :func:`annotate_variant`: a gained stop dominates, then a destroyed
    start, then a destroyed stop, then an ordinary residue change.
    """
    ref, alt = (ref or "").upper(), (alt or "").upper()
    if len(ref) == len(alt):
        kind, aa, changes = _classify_substitution(prot_ref, prot_alt)
        return VariantConsequence(kind, aa, changes)

    net = len(alt) - len(ref)
    # 🔴 THE GUARD THAT CATCHES THE CALLER'S WINDOW BUG, because a caller supplying translations has
    # already done the splicing and can have done it wrong. A correctly spliced in-frame indel MUST
    # change the protein's length, by net//3 residues; equal lengths mean the edited CDS was cut
    # with the UNSHIFTED exon coordinates, so the inserted bases displaced downstream ones and the
    # terminal stop may have fallen out of the window. Measured on 10 real designs at net +3 to +27:
    # equal lengths on all 10, 4 with no stop left at all, and 3 whose "premature" stop sat exactly
    # at the reference stop, which makes the premature-stop test true by construction.
    # Refusing is the only safe answer: a plausible wrong consequence is worse than none, and
    # `shift_exons_for_indel` plus `annotate_variant` is the path that gets it right.
    if net % 3 == 0 and net != 0 and len(prot_alt) == len(prot_ref):
        return VariantConsequence(
            "non_coding",
            detail=(f"unresolvable: net {net:+d} nt but both translations are {len(prot_ref)} "
                    f"residues, so the edited CDS was spliced with the UNSHIFTED exon coordinates. "
                    f"Use annotate_variant, or shift_exons_for_indel before translating."),
        )
    kind, aa, changes, detail = _classify_indel(
        prot_ref, prot_alt, ref, alt, 0 if codon is None else codon
    )
    if codon is None:
        # the residue number in an indel's detail is the one thing `codon` was needed for
        detail = detail.split(" at residue ")[0].split("; every downstream")[0]
    return VariantConsequence(kind, aa, changes, detail=detail)


def annotate_variant(
    genome: Dict[str, str],
    model: CdsModel,
    gene: str,
    chrom: str,
    pos0: int,
    ref: str,
    alt: str,
    *,
    verify_ref: bool = True,
) -> VariantConsequence:
    """The consequence of ``ref``->``alt`` at ``pos0`` in ``gene``'s frame.

    ``pos0`` is 0-based. ``ref``/``alt`` are the genomic alleles on the plus
    strand; either may be empty for a pure insertion or deletion.

    Returns ``non_coding`` when the edit is not inside this gene's CDS. Never
    returns an empty consequence.

    With ``verify_ref`` (the default), the reference allele is asserted against
    the genome. A coordinate that is right in one frame and wrong in another is
    the most common way this annotation goes silently wrong, so the sequence is
    checked rather than trusted.
    """
    ref, alt = (ref or "").upper(), (alt or "").upper()
    seq = genome.get(chrom)
    if seq is None:
        raise KeyError(f"{chrom} is not in the supplied genome")
    if verify_ref and ref and seq[pos0 : pos0 + len(ref)].upper() != ref:
        raise ValueError(
            f"reference allele {ref!r} does not match the genome at {chrom}:{pos0} "
            f"(found {seq[pos0:pos0 + len(ref)].upper()!r}). Assert against sequence, "
            "never coordinates alone."
        )

    entry = model.get(gene)
    if entry is None:
        return VariantConsequence("non_coding", gene=gene, detail="gene not in CDS model")
    g_chrom, strand, exons = entry
    if g_chrom != chrom:
        return VariantConsequence(
            "non_coding", gene=gene, detail=f"gene is on {g_chrom}, edit on {chrom}"
        )
    offset = cds_index(exons, strand, pos0)
    if offset is None:
        return VariantConsequence(
            "non_coding", gene=gene, detail="edit is not inside this gene's CDS"
        )

    prot_ref = translate(build_cds(seq, exons, strand))
    edited = seq[:pos0] + alt + seq[pos0 + len(ref) :]
    codon = offset // 3

    if len(ref) == len(alt):
        prot_alt = translate(build_cds(edited, exons, strand))
        kind, change, changes = _classify_substitution(prot_ref, prot_alt)
        detail = "" if len(changes) <= 1 else f"{len(changes)} residues altered"
        return VariantConsequence(kind, change, changes, gene, detail)

    # 🔴 A LENGTH-CHANGING ALLELE MOVES EVERY DOWNSTREAM EXON BOUNDARY. Splicing the edited genome
    # with the original coordinates gives a CDS of the SAME width, so an insertion displaces
    # downstream bases instead of extending the protein and can push the terminal stop out of the
    # window. Measured on 10 real designs at net +3 to +27: len(prot_alt) == len(prot_ref) on all
    # of them. The shift is what makes the edited protein real.
    shifted = shift_exons_for_indel(exons, pos0, len(alt) - len(ref))
    if shifted is None:
        return VariantConsequence(
            "non_coding", gene=gene,
            detail="a length-changing allele spanning an exon boundary; the edited reading frame "
                   "is not resolvable from this model, so no consequence is claimed",
        )
    prot_alt = translate(build_cds(edited, shifted, strand))
    kind, change, changes, detail = _classify_indel(prot_ref, prot_alt, ref, alt, codon)
    return VariantConsequence(kind, change, changes, gene, detail)


def annotate_overlapping(
    genome: Dict[str, str],
    model: CdsModel,
    index: Dict[str, List[Tuple[int, int, str]]],
    chrom: str,
    pos0: int,
    ref: str,
    alt: str,
    **kwargs,
) -> List[VariantConsequence]:
    """One result per gene whose CDS contains the edit; ``[]`` when none does.

    A position can sit inside two CDSs at once, and a single scalar consequence
    cannot express that. Callers that need one value must choose explicitly.
    """
    return [
        annotate_variant(genome, model, g, chrom, pos0, ref, alt, **kwargs)
        for g in genes_at(index, chrom, pos0)
    ]


# ---------------------------------------------------------------- self-test

# Runnable without pytest, as `python src/magestic/genomics/consequence.py`,
# mirroring the `--selftest` convention used elsewhere in the package. Every
# check below is a case the previous annotators left BLANK, plus two planted
# failures that must raise.


def _selftest() -> int:  # pragma: no cover - exercised by tests and by hand
    failures = 0

    def chk(name, got, want):
        nonlocal failures
        ok = got == want
        print(f"  {'PASS' if ok else 'FAIL'} {name}: got {got!r} want {want!r}")
        if not ok:
            failures += 1

    # A synthetic gene: ATG AAA CCC GGG TTT TGG TAA on the plus strand.
    plus = "ATGAAACCCGGGTTTTGGTAA"
    genome = {"chrT": "GGGG" + plus + "GGGG"}
    start = 4
    model: CdsModel = {"PLUS": ("chrT", "+", [(start, start + len(plus))])}
    idx = build_gene_index(model)

    chk("ref protein", translate(build_cds(genome["chrT"], model["PLUS"][2], "+")),
        "MKPGFW*")

    # --- single-base substitution, the only case the old annotators handled
    chk("SNV synonymous (AAA->AAG, K)",
        annotate_variant(genome, model, "PLUS", "chrT", start + 5, "A", "G").consequence,
        "synonymous")
    chk("SNV missense (AAA->GAA, K2E)",
        annotate_variant(genome, model, "PLUS", "chrT", start + 3, "A", "G").aa_change,
        "K2E")

    # --- multi-base MNP spanning TWO codons: left blank before, and the single
    #     most common coding blank class (120 of 365 were 4 nt over 2 codons).
    #     AAC->TTG at plus[4:7] straddles codon 2 (AAA->ATT) and codon 3 (CCC->GCC).
    mnp = annotate_variant(genome, model, "PLUS", "chrT", start + 4, "AAC", "TTG")
    chk("MNP over two codons is classified", mnp.consequence, "missense")
    chk("MNP reports BOTH changed residues", mnp.aa_changes, ("K2I", "P3A"))

    # --- MNP that creates a stop
    chk("MNP gaining a stop is nonsense",
        annotate_variant(genome, model, "PLUS", "chrT", start + 3, "AAA", "TAA").consequence,
        "nonsense")

    # --- start lost, stop lost
    chk("start codon destroyed",
        annotate_variant(genome, model, "PLUS", "chrT", start, "A", "G").consequence,
        "start_lost")
    chk("terminal stop destroyed",
        annotate_variant(genome, model, "PLUS", "chrT", start + 18, "T", "C").consequence,
        "stop_lost")

    # --- indels: 100 % blank before this module
    chk("3 nt deletion is an in-frame deletion",
        annotate_variant(genome, model, "PLUS", "chrT", start + 3, "AAA", "").consequence,
        "inframe_deletion")
    chk("1 nt deletion is a frameshift",
        annotate_variant(genome, model, "PLUS", "chrT", start + 3, "A", "").consequence,
        "frameshift")
    chk("3 nt insertion is an in-frame insertion",
        annotate_variant(genome, model, "PLUS", "chrT", start + 3, "", "GGG").consequence,
        "inframe_insertion")
    chk("2 nt insertion is a frameshift",
        annotate_variant(genome, model, "PLUS", "chrT", start + 3, "", "GG").consequence,
        "frameshift")

    # --- minus strand must be handled, not merely accepted
    rc_seq = _rev_comp(plus)
    genome_m = {"chrM": "GGGG" + rc_seq + "GGGG"}
    model_m: CdsModel = {"MINUS": ("chrM", "-", [(4, 4 + len(rc_seq))])}
    chk("minus-strand ref protein",
        translate(build_cds(genome_m["chrM"], model_m["MINUS"][2], "-")), "MKPGFW*")
    chk("minus-strand missense is K2E",
        annotate_variant(genome_m, model_m, "MINUS", "chrM",
                         4 + len(rc_seq) - 1 - 3, "T", "C").aa_change, "K2E")

    # --- non_coding is a POSITIVE answer, never a blank
    nc = annotate_variant(genome, model, "PLUS", "chrT", 0, "G", "A")
    chk("an edit outside the CDS is non_coding", nc.consequence, "non_coding")
    chk("non_coding carries no residue change", nc.aa_change, "")
    chk("a gene absent from the model is non_coding",
        annotate_variant(genome, model, "ABSENT", "chrT", start + 3, "A", "G").consequence,
        "non_coding")

    # --- overlapping genes: one position, two frames
    model2 = dict(model)
    model2["OTHER"] = ("chrT", "-", [(start, start + len(plus))])
    idx2 = build_gene_index(model2)
    chk("both overlapping genes are found", genes_at(idx2, "chrT", start + 3),
        ["OTHER", "PLUS"])
    chk("one result per containing gene",
        len(annotate_overlapping(genome, model2, idx2, "chrT", start + 3, "A", "G")), 2)

    # --- PLANTED FAILURE 1: a wrong reference allele must raise, not annotate.
    try:
        annotate_variant(genome, model, "PLUS", "chrT", start + 3, "C", "G")
        print("  FAIL planted wrong ref allele was accepted")
        failures += 1
    except ValueError:
        print("  PASS planted wrong ref allele raises")

    # --- PLANTED FAILURE 2: an unknown consequence string must be rejected.
    try:
        VariantConsequence("")
        print("  FAIL a blank consequence was accepted")
        failures += 1
    except ValueError:
        print("  PASS a blank consequence is rejected by the dataclass")

    # --- PLANTED FAILURE 3: a headerless table read as if it had a header.
    import tempfile, os as _os
    fd, tmp = tempfile.mkstemp(suffix=".tsv")
    with _os.fdopen(fd, "w") as fh:
        fh.write("chrom\tstart\tend\tstrand\tphase\tattrs\n")  # a header, no data
    try:
        load_cds_table(tmp)
        print("  FAIL an empty CDS model was returned instead of raising")
        failures += 1
    except ValueError:
        print("  PASS an empty CDS table raises instead of making everything non_coding")
    finally:
        _os.unlink(tmp)

    total = 22
    print(
        f"consequence selftest: {total - failures}/{total} checks passed, "
        "including 3 planted failures"
        if failures
        else f"consequence selftest: all {total} checks passed, including 3 planted failures"
    )
    return 1 if failures else 0


if __name__ == "__main__":  # pragma: no cover
    import sys as _sys

    logging.basicConfig(level=logging.INFO)
    _sys.exit(_selftest())
