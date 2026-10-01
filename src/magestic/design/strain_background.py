"""
magestic.design.strain_background — de novo background-strain reference preparation.

This module packages the step that was previously project-specific and undistributed:
``get_haplotype_blocks_from_vcf_with_coord_offset.py``. Its absence is the root cause of
two independent classes of design error seen in collaborator reimplementations
(EMBL/TRIUMPH, 2026-09), because every reimplementation has to rediscover the same
four traps.

What it does, in the order the design pipeline needs it:

1. :func:`load_strain_variants` — read a multi-sample VCF and extract one strain's
   variants, **consulting FILTER, QUAL and genotype** rather than blindly taking
   ``ALT[0]``.
2. :func:`build_strain_reference` — substitute those variants into the reference FASTA
   to produce the strain's own design reference, applying variants right-to-left and
   **retaining every contig** (mitochondrion, 2-micron, plasmids).
3. :class:`CoordinateMap` — the bidirectional reference↔strain coordinate mapping that
   indels make necessary. This is the "coord_adjusted" in
   ``Bloom_et_al_16_strains_QTL_genes_by_haplotype_block_coord_adjusted.tsv``.
4. :func:`remap_gff` — shift an annotation onto the new coordinates, flagging any
   feature whose span contains an indel.
5. :func:`project_variants` — express another strain's variants (called against the
   shared reference) in *this* strain's coordinates. This is what makes an all-by-all
   variant swap across N strains tractable: N references, and for each one the other
   N-1 strains' alleles as edits.
6. :func:`verify_designs_against_assembly` — exact-substring check of a finished guide
   and its donor homology arms against a de novo assembly, so a locally wrong
   substituted reference is caught before synthesis rather than after.

The four traps, each of which was observed as a real defect rather than imagined:

- **Dropping contigs.** A FASTA reader that stops at the mitochondrial defline loses
  chrM, the 2-micron and any plasmid. Project policy requires the A364A 2-micron variant
  in *every* reference, even for [cir0], because a generic 2-micron reads 0x against
  DHY214 and the mis-call is silent.
- **Position-keyed variant dicts.** ``variants[chrom][pos] = record`` silently discards
  one of two records at the same position (multi-allelic sites, or an SNV and an indel
  anchored together).
- **Left-to-right application.** Applying an indel shifts every downstream coordinate,
  so the next variant lands in the wrong place. Right-to-left is the only order that
  needs no bookkeeping.
- **Ignoring genotype.** Applying ``ALT[0]`` with no GT check writes heterozygous and
  reference-genotype calls into the reference as if they were homozygous alternates.

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

import gzip
import logging
import re
from bisect import bisect_right
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

logger = logging.getLogger(__name__)

__all__ = [
    "StrainVariant",
    "CoordinateMap",
    "VariantRejection",
    "load_strain_variants",
    "build_strain_reference",
    "remap_gff",
    "project_variants",
    "annotate_variants_with_features",
    "verify_designs_against_assembly",
    "read_fasta",
    "write_fasta",
]


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StrainVariant:
    """One variant to substitute into a reference.

    Coordinates are **1-based inclusive** on ``chrom``, matching VCF and GFF.
    ``ref`` is the reference allele as written in the VCF, so ``pos`` through
    ``pos + len(ref) - 1`` is the span this variant replaces.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    sample: str = ""
    genotype: str = ""
    orf: str = ""
    name: str = ""

    @property
    def end(self) -> int:
        """Last reference position (1-based, inclusive) that this variant replaces."""
        return self.pos + len(self.ref) - 1

    @property
    def length_change(self) -> int:
        """Net length change in the strain relative to the reference."""
        return len(self.alt) - len(self.ref)

    @property
    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt)


@dataclass(frozen=True)
class VariantRejection:
    """A variant that was not applied, and why. Rejections are reported, never silent."""

    variant: StrainVariant
    reason: str


@dataclass
class CoordinateMap:
    """Bidirectional coordinate mapping between a reference and a strain sequence.

    Built from the applied variants of one chromosome. Internally it stores, for each
    applied variant, the reference start and the cumulative length change *after* that
    variant, so a lookup is a binary search rather than a per-base table. That matters:
    a per-base dict over 12 Mb of yeast genome is ~1 GB in Python and is the usual
    reason a naive implementation cannot run on a laptop.
    """

    #: chrom -> sorted list of reference start positions of applied variants
    _starts: Dict[str, List[int]] = field(default_factory=dict)
    #: chrom -> cumulative offset in effect *after* the variant at the same index
    _cum_offsets: Dict[str, List[int]] = field(default_factory=dict)
    #: chrom -> (ref_start, ref_end, length_change) for each applied variant
    _spans: Dict[str, List[Tuple[int, int, int]]] = field(default_factory=dict)

    def add_chromosome(self, chrom: str, applied: Sequence[StrainVariant]) -> None:
        """Register the applied variants for one chromosome, in reference order."""
        ordered = sorted(applied, key=lambda v: v.pos)
        starts: List[int] = []
        cum: List[int] = []
        spans: List[Tuple[int, int, int]] = []
        running = 0
        for v in ordered:
            running += v.length_change
            starts.append(v.pos)
            cum.append(running)
            spans.append((v.pos, v.end, v.length_change))
        self._starts[chrom] = starts
        self._cum_offsets[chrom] = cum
        self._spans[chrom] = spans

    def total_length_change(self, chrom: str) -> int:
        cum = self._cum_offsets.get(chrom)
        return cum[-1] if cum else 0

    def ref_to_strain(self, chrom: str, ref_pos: int) -> Optional[int]:
        """Map a 1-based reference position to the strain coordinate.

        Returns ``None`` when the reference base was deleted in the strain, which is a
        real answer and not an error: a design anchored on a deleted base has no target.
        """
        starts = self._starts.get(chrom)
        if not starts:
            return ref_pos

        idx = bisect_right(starts, ref_pos) - 1
        if idx < 0:
            return ref_pos

        span_start, span_end, delta = self._spans[chrom][idx]
        offset_before = self._cum_offsets[chrom][idx - 1] if idx > 0 else 0

        if ref_pos > span_end:
            # Strictly downstream of this variant: the full cumulative offset applies.
            return ref_pos + self._cum_offsets[chrom][idx]

        # Inside the replaced span.
        alt_len = (span_end - span_start + 1) + delta
        offset_in_span = ref_pos - span_start
        if offset_in_span < alt_len:
            return span_start + offset_before + offset_in_span
        # Reference base falls in the deleted tail of the span.
        return None

    def strain_to_ref(self, chrom: str, strain_pos: int) -> Optional[int]:
        """Map a 1-based strain position back to the reference coordinate."""
        starts = self._starts.get(chrom)
        if not starts:
            return strain_pos

        # Walk the variants in order, tracking the strain coordinate of each span.
        offset = 0
        for (span_start, span_end, delta) in self._spans[chrom]:
            strain_span_start = span_start + offset
            alt_len = (span_end - span_start + 1) + delta
            if strain_pos < strain_span_start:
                return strain_pos - offset
            if strain_pos < strain_span_start + alt_len:
                # Inside an inserted/substituted span; anchor to the span start.
                return span_start + min(strain_pos - strain_span_start,
                                        span_end - span_start)
            offset += delta
        return strain_pos - offset


# ---------------------------------------------------------------------------
# FASTA I/O — retains every contig, which is the point
# ---------------------------------------------------------------------------


def _open_maybe_gzip(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path: str | Path) -> Dict[str, str]:
    """Read a FASTA into ``{contig_name: sequence}``, keeping **every** contig.

    The contig name is the first whitespace-delimited token of the defline. No defline
    is treated as a terminator: a reader that stops at the mitochondrial or 2-micron
    record silently produces a reference that cannot call those loci.
    """
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    logger.info("read_fasta: %d contigs from %s", len(seqs), path)
    return seqs


def write_fasta(seqs: Dict[str, str], path: str | Path, width: int = 60) -> None:
    """Write ``{contig: sequence}`` to FASTA with fixed line width."""
    with open(path, "w") as out:
        for name, seq in seqs.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


# ---------------------------------------------------------------------------
# VCF loading — genotype-aware
# ---------------------------------------------------------------------------

_HET_GT = re.compile(r"^([0-9])[/|]([0-9])$")


def _genotype_allele_index(gt: str, require_homozygous: bool) -> Optional[int]:
    """Return the ALT index to apply for a genotype string, or ``None`` to skip.

    ``0`` means the reference allele, so a ``0/0`` call contributes no variant. A
    heterozygous call is skipped when *require_homozygous* is set, because writing one
    allele of a heterozygous site into a haploid design reference asserts something the
    data does not support.
    """
    gt = gt.split(":")[0].strip()
    if gt in {".", "./.", ".|.", ""}:
        return None
    m = _HET_GT.match(gt)
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        if a != b:
            return None if require_homozygous else max(a, b) or None
        return a if a > 0 else None
    if gt.isdigit():
        idx = int(gt)
        return idx if idx > 0 else None
    return None


def load_strain_variants(
    vcf_path: str | Path,
    sample: Optional[str] = None,
    *,
    require_pass: bool = True,
    min_qual: Optional[float] = None,
    require_homozygous: bool = True,
    chroms: Optional[Set[str]] = None,
) -> Tuple[List[StrainVariant], List[VariantRejection]]:
    """Load one strain's variants from a VCF, consulting FILTER, QUAL and genotype.

    Returns ``(variants, rejections)``. Rejections carry a reason so the caller can
    report how many calls were dropped and why — a count that a position-keyed dict
    cannot produce, because its losses are invisible.

    Parameters
    ----------
    sample
        Sample column to read. ``None`` uses the first sample column; for a
        single-sample VCF with no genotype columns, every record is taken.
    require_pass
        Drop records whose FILTER is neither ``PASS`` nor ``.``.
    min_qual
        Drop records with QUAL below this value. ``None`` disables the check.
    require_homozygous
        Drop heterozygous calls rather than guessing which allele to write.
    """
    variants: List[StrainVariant] = []
    rejections: List[VariantRejection] = []
    sample_idx: Optional[int] = None
    sample_name = sample or ""

    with _open_maybe_gzip(vcf_path) as fh:
        for raw in fh:
            if raw.startswith("##"):
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if raw.startswith("#CHROM"):
                if len(fields) > 9:
                    samples = fields[9:]
                    if sample is None:
                        sample_idx, sample_name = 0, samples[0]
                    elif sample in samples:
                        sample_idx, sample_name = samples.index(sample), sample
                    else:
                        raise KeyError(
                            f"sample {sample!r} not in VCF; available: {samples}"
                        )
                continue
            if len(fields) < 5:
                continue

            chrom, pos_s, _id, ref, alt_s = fields[0], fields[1], fields[2], fields[3], fields[4]
            if chroms is not None and chrom not in chroms:
                continue
            flt = fields[6] if len(fields) > 6 else "."
            qual_s = fields[5] if len(fields) > 5 else "."

            alts = alt_s.split(",")
            alt_index = 1
            genotype = ""
            if sample_idx is not None and len(fields) > 9 + sample_idx:
                genotype = fields[9 + sample_idx]
                idx = _genotype_allele_index(genotype, require_homozygous)
                if idx is None:
                    rejections.append(VariantRejection(
                        StrainVariant(chrom, int(pos_s), ref, alts[0], sample_name, genotype),
                        "genotype is reference, missing, or heterozygous",
                    ))
                    continue
                alt_index = idx

            if alt_index > len(alts):
                rejections.append(VariantRejection(
                    StrainVariant(chrom, int(pos_s), ref, alts[0], sample_name, genotype),
                    f"genotype allele index {alt_index} exceeds ALT count {len(alts)}",
                ))
                continue
            alt = alts[alt_index - 1]

            cand = StrainVariant(chrom, int(pos_s), ref.upper(), alt.upper(),
                                 sample_name, genotype)

            if require_pass and flt not in {"PASS", ".", ""}:
                rejections.append(VariantRejection(cand, f"FILTER={flt}"))
                continue
            if min_qual is not None and qual_s not in {".", ""}:
                try:
                    if float(qual_s) < min_qual:
                        rejections.append(VariantRejection(cand, f"QUAL={qual_s} < {min_qual}"))
                        continue
                except ValueError:
                    pass
            if alt in {".", "*", ""} or not re.fullmatch(r"[ACGTN]+", alt):
                rejections.append(VariantRejection(cand, f"non-nucleotide ALT {alt!r}"))
                continue

            variants.append(cand)

    logger.info(
        "load_strain_variants: sample=%s kept=%d rejected=%d",
        sample_name or "<single>", len(variants), len(rejections),
    )
    return variants, rejections


# ---------------------------------------------------------------------------
# Reference construction
# ---------------------------------------------------------------------------


def build_strain_reference(
    reference: Dict[str, str],
    variants: Iterable[StrainVariant],
    *,
    check_ref_allele: bool = True,
) -> Tuple[Dict[str, str], CoordinateMap, List[VariantRejection]]:
    """Substitute *variants* into *reference* and return the strain reference.

    Returns ``(strain_sequences, coordinate_map, rejections)``.

    Every contig of *reference* appears in the output, whether or not it carries a
    variant. Variants are applied **right-to-left** within each contig so that an indel
    never displaces a variant that has not yet been applied. Overlapping variants are
    rejected rather than applied, because the second of two overlapping edits would be
    written against sequence the first one already replaced.

    ``check_ref_allele`` verifies each record's REF against the reference sequence. A
    mismatch almost always means the VCF was called against a different assembly or a
    different coordinate convention, and it is the cheapest available detector for that.
    """
    by_chrom: Dict[str, List[StrainVariant]] = {}
    rejections: List[VariantRejection] = []

    for v in variants:
        if v.chrom not in reference:
            rejections.append(VariantRejection(v, f"contig {v.chrom!r} absent from reference"))
            continue
        by_chrom.setdefault(v.chrom, []).append(v)

    out: Dict[str, str] = dict(reference)
    cmap = CoordinateMap()

    for chrom, seq in reference.items():
        chrom_vars = sorted(by_chrom.get(chrom, []), key=lambda v: (v.pos, v.end))
        if not chrom_vars:
            cmap.add_chromosome(chrom, [])
            continue

        accepted: List[StrainVariant] = []
        last_end = 0
        for v in chrom_vars:
            if v.pos < 1 or v.end > len(seq):
                rejections.append(VariantRejection(v, "span outside contig bounds"))
                continue
            if check_ref_allele:
                observed = seq[v.pos - 1:v.end].upper()
                if observed != v.ref:
                    rejections.append(VariantRejection(
                        v, f"REF mismatch: VCF {v.ref!r} vs reference {observed!r}"
                    ))
                    continue
            if v.pos <= last_end:
                rejections.append(VariantRejection(
                    v, f"overlaps a preceding variant ending at {last_end}"
                ))
                continue
            accepted.append(v)
            last_end = v.end

        # Right-to-left application: no coordinate bookkeeping required.
        chars = list(seq)
        for v in sorted(accepted, key=lambda v: v.pos, reverse=True):
            chars[v.pos - 1:v.end] = list(v.alt)
        out[chrom] = "".join(chars)

        cmap.add_chromosome(chrom, accepted)
        logger.info(
            "build_strain_reference: %s applied=%d rejected=%d len %d -> %d",
            chrom, len(accepted), len(chrom_vars) - len(accepted), len(seq), len(out[chrom]),
        )

    return out, cmap, rejections


# ---------------------------------------------------------------------------
# GFF remapping
# ---------------------------------------------------------------------------


def remap_gff(
    gff_path: str | Path,
    cmap: CoordinateMap,
    out_path: str | Path,
) -> Dict[str, int]:
    """Rewrite a GFF3 onto strain coordinates. Returns a counts summary.

    GFF is 1-based inclusive on both ends, the same convention as
    :class:`CoordinateMap`, so no off-by-one conversion is applied here. Mixing this up
    with 0-based BED or Python slicing is the single most common source of a design that
    is uniformly one base out.

    A feature whose span contains an indel is still remapped, but counted under
    ``spans_indel`` so it can be reviewed: its length changes, and for a CDS that may
    mean a frameshift relative to the annotation.
    """
    counts = {"written": 0, "spans_indel": 0, "dropped": 0, "comment": 0}

    with _open_maybe_gzip(gff_path) as fh, open(out_path, "w") as out:
        for raw in fh:
            if raw.startswith("#"):
                out.write(raw)
                counts["comment"] += 1
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 9:
                out.write(raw)
                continue

            chrom, start_s, end_s = fields[0], fields[3], fields[4]
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                out.write(raw)
                continue

            new_start = cmap.ref_to_strain(chrom, start)
            new_end = cmap.ref_to_strain(chrom, end)
            if new_start is None or new_end is None or new_end < new_start:
                counts["dropped"] += 1
                continue

            if (new_end - new_start) != (end - start):
                counts["spans_indel"] += 1

            fields[3], fields[4] = str(new_start), str(new_end)
            out.write("\t".join(fields) + "\n")
            counts["written"] += 1

    logger.info("remap_gff: %s", counts)
    return counts


# ---------------------------------------------------------------------------
# Cross-strain variant projection — the all-by-all enabler
# ---------------------------------------------------------------------------


def project_variants(
    variants: Iterable[StrainVariant],
    cmap: CoordinateMap,
    recipient_reference: Dict[str, str],
    *,
    drop_if_identical: bool = True,
) -> Tuple[List[StrainVariant], List[VariantRejection]]:
    """Express donor-strain *variants* in the recipient strain's coordinates.

    For an all-by-all swap across N strains you build N references and, for each
    recipient, project the other N-1 strains' alleles onto it. Two things have to happen
    and both are easy to get wrong:

    - the coordinate moves, because the recipient's own indels shifted it;
    - **a variant the recipient already carries is not an edit at all.** Those are
      dropped when *drop_if_identical* is set. This is the same concern as
      ``filter_variants(exclude_background=True)``, applied pairwise instead of against
      a single background, and leaving them in produces oligos that edit nothing.
    """
    projected: List[StrainVariant] = []
    rejections: List[VariantRejection] = []

    for v in variants:
        if v.chrom not in recipient_reference:
            rejections.append(VariantRejection(v, f"contig {v.chrom!r} absent from recipient"))
            continue
        new_pos = cmap.ref_to_strain(v.chrom, v.pos)
        if new_pos is None:
            rejections.append(VariantRejection(v, "position deleted in recipient"))
            continue

        seq = recipient_reference[v.chrom]
        new_end = new_pos + len(v.ref) - 1
        if new_pos < 1 or new_end > len(seq):
            rejections.append(VariantRejection(v, "projected span outside recipient contig"))
            continue

        recipient_allele = seq[new_pos - 1:new_end].upper()
        if drop_if_identical and recipient_allele == v.alt.upper():
            rejections.append(VariantRejection(v, "recipient already carries the ALT allele"))
            continue

        projected.append(StrainVariant(
            chrom=v.chrom, pos=new_pos, ref=recipient_allele, alt=v.alt.upper(),
            sample=v.sample, genotype=v.genotype, orf=v.orf, name=v.name,
        ))

    logger.info("project_variants: projected=%d rejected=%d", len(projected), len(rejections))
    return projected, rejections


# ---------------------------------------------------------------------------
# Feature annotation — the fix for silently dropped non-coding variants
# ---------------------------------------------------------------------------


def annotate_variants_with_features(
    variants: Iterable[StrainVariant],
    gff_path: str | Path,
    *,
    feature_types: Sequence[str] = ("gene", "CDS", "mRNA"),
    flank: int = 1000,
) -> List[StrainVariant]:
    """Assign ``orf``/``name`` to each variant, labeling non-coding ones explicitly.

    Downstream, ``filter_variants`` keeps a variant only when its ``ORF`` or ``NAME`` is
    in the supplied gene list, so a variant left with an empty ``ORF`` is dropped with no
    message. That is how a natural-variant library ends up with no non-coding members
    while every individual step looks like it worked.

    This function therefore never leaves ``orf`` empty. A variant inside a feature gets
    that feature; one within *flank* bp gets ``upstream:<gene>`` or ``downstream:<gene>``;
    anything else gets ``intergenic``. Those labels are real identifiers, so a caller can
    deliberately include or exclude non-coding classes instead of losing them by default.
    """
    feats: Dict[str, List[Tuple[int, int, str, str, str]]] = {}
    wanted = set(feature_types)

    with _open_maybe_gzip(gff_path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            f = raw.rstrip("\n\r").split("\t")
            if len(f) < 9 or f[2] not in wanted:
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            attrs = f[8]
            orf = _gff_attr(attrs, "ID") or _gff_attr(attrs, "Name") or ""
            name = _gff_attr(attrs, "gene") or _gff_attr(attrs, "Name") or orf
            feats.setdefault(f[0], []).append((start, end, orf, name, f[6]))

    for chrom in feats:
        feats[chrom].sort()

    out: List[StrainVariant] = []
    for v in variants:
        chrom_feats = feats.get(v.chrom, [])
        orf, name = "intergenic", "intergenic"
        best_dist: Optional[int] = None

        for (start, end, f_orf, f_name, strand) in chrom_feats:
            if start - flank > v.end:
                break
            if end + flank < v.pos:
                continue
            if v.pos <= end and v.end >= start:
                orf, name, best_dist = f_orf, f_name, 0
                break
            if v.end < start:
                dist = start - v.end
                label = "upstream" if strand == "+" else "downstream"
            else:
                dist = v.pos - end
                label = "downstream" if strand == "+" else "upstream"
            if best_dist is None or dist < best_dist:
                best_dist, orf, name = dist, f"{label}:{f_orf}", f_name

        out.append(StrainVariant(
            chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt,
            sample=v.sample, genotype=v.genotype, orf=orf, name=name,
        ))
    return out


def _gff_attr(attrs: str, key: str) -> str:
    for kv in attrs.split(";"):
        kv = kv.strip()
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return ""


# ---------------------------------------------------------------------------
# Design verification against a de novo assembly
# ---------------------------------------------------------------------------


def verify_designs_against_assembly(
    designs: Iterable[Dict[str, str]],
    assembly: Dict[str, str],
    *,
    guide_key: str = "guide",
    left_arm_key: str = "left_homology_arm",
    right_arm_key: str = "right_homology_arm",
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    """Split *designs* into (passed, failed) by exact presence in a de novo assembly.

    This is the recommended guard when a design reference was built by substituting
    variants into a distant reference: the substituted reference can be locally wrong
    wherever the variant calls were incomplete, and a guide or homology arm built on
    wrong sequence fails in a way that is indistinguishable from a biological negative.

    A design passes when its guide and both homology arms each occur as an exact
    substring of the assembly, on either strand. Each failed design carries a
    ``failure_reason``. Nothing here is a substitute for alignment; it is a cheap
    sufficient check that catches the case that matters.
    """
    haystack = "\n".join(assembly.values()).upper()
    haystack_rc = _revcomp(haystack)

    passed: List[Dict[str, str]] = []
    failed: List[Dict[str, str]] = []

    for design in designs:
        missing: List[str] = []
        for key in (guide_key, left_arm_key, right_arm_key):
            seq = (design.get(key) or "").upper().replace("-", "")
            if not seq:
                continue
            if seq not in haystack and seq not in haystack_rc:
                missing.append(key)
        if missing:
            rec = dict(design)
            rec["failure_reason"] = "absent from assembly: " + ",".join(missing)
            failed.append(rec)
        else:
            passed.append(dict(design))

    logger.info(
        "verify_designs_against_assembly: passed=%d failed=%d", len(passed), len(failed)
    )
    return passed, failed


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]
