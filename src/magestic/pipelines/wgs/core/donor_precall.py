"""magestic.pipelines.wgs.core.donor_precall — alias contigs IN, decode them OUT (defects 3, 8).

DEFECT 3 — WIRING `ContigAlias`, WITH THE DECODE THAT MAKES IT SAFE
--------------------------------------------------------------------
`trace_crispr.wgs.donor_outcome.precall_sample` formats its own variant id as
``f"{chrom}_{pos}_{ref}_to_{alt}"`` and then re-parses it with `split("_")`,
indexing positionally (`sp[1]` pos, `sp[2]` ref, `sp[4]` alt) at
`donor_outcome.py:325` and `:334`. Any underscore in `chrom` shifts every index, so
`int(sp[1])` raises `ValueError` **mid-loop, after a partial output is on disk**.

Measured: **4 of the 5 contigs the WGS reference-mapping policy mandates crash it** —
`unplaced_contig_34`, `blank_MAGESTIC3_pF835_guide_donor_plasmid`, `native_2micron`,
`integrated_plasmid_o1`. Only `chrXV` survives. The stage is incompatible with the
reference standard, not with one odd dataset.

`core/contig_alias.py` (a peer's work, 17 tests) solves the representation half and
was **deliberately not wired in**, with a correct reason: aliasing without decoding
would ship `ctg0003_131_A_to_G` into a shipped table, replacing a loud crash with a
silent wrong name. This module is the other half.

🔴 **Every input TRACE sees is encoded, and every output is decoded.** The three
names `precall_sample` reads are the spec's `chrom`, the `.fil.vcf` CHROM column and
the `.g.vcf` CHROM column; `precall`'s driver additionally reads FASTA headers, which
`alias_fasta` handles. Miss one and `v_chr != chrom` skips every variant with **no
crash and a wrong outcome** — a strictly worse failure than today's. The test suite
asserts (a) variants are actually bucketed on an underscore contig and (b) no string
beginning `ctg` reaches any shipped field.

THE TRACE-SIDE FIX (precise, for when the pin moves)
-----------------------------------------------------
`precall_sample` should carry structured fields instead of re-parsing a string it
just formatted. Failing that, parse RIGHT-anchored — unambiguous because `ref`/`alt`
are nucleotide strings and `to` is a literal::

    sp = v_id.split("_")
    alt   = sp[-1]
    ref   = sp[-3]
    pos   = int(sp[-4])
    chrom = "_".join(sp[:-4])

That is `donor_outcome.py:325-327` and `:334-336`, four lines each, and it removes
the need for this module entirely.

DEFECT 8 — DISTAL VARIANTS ARE DISCARDED INVISIBLY
---------------------------------------------------
`precall_sample`'s VCF filter is `donor_outcome.py:265`::

    if v_chr != chrom or pos0 < 0 or (pos0 + len(v_ref)) >= wt_len: continue

Every variant outside the donor window — a distal unintended edit, an off-target, a
variant on another contig — is dropped with no count and no log line, so **distal
unintended edits are invisible** and the VCF's record count never reconciles with the
outcome's variant count. `audit_vcf_window` reproduces TRACE's arithmetic exactly and
reports what it discards, with the ids; `flank` makes the window configurable so a
screen can widen it deliberately rather than by patching TRACE.
"""

from __future__ import annotations

import dataclasses
import gzip
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import _trace
from .contig_alias import ContigAlias

__all__ = [
    "AliasedPrecall",
    "VcfWindowAudit",
    "alias_fasta",
    "alias_vcf_contigs",
    "audit_vcf_window",
    "decode_outcome",
    "precall_sample_aliased",
]


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


# --------------------------------------------------------------------------- #
# Defect 3 — encode on the way in, decode on the way out                       #
# --------------------------------------------------------------------------- #


def alias_vcf_contigs(src: str | Path, dst: str | Path, alias: ContigAlias) -> int:
    """Rewrite a VCF/gVCF's CHROM column (and `##contig=<ID=...>`) to aliases.

    Returns the number of DATA records rewritten. Header `##contig` lines are mapped
    too, so the file stays self-consistent for any tool that reads them; TRACE reads
    only the data columns, so this is belt and braces rather than load-bearing.

    Output is always plain text, because TRACE opens these with a bare `open()`.
    """
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with _open_text(src) as fh, dst.open("w", newline="\n") as out:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("##contig=<ID="):
                head, _, rest = line.partition("##contig=<ID=")
                name, sep, tail = rest.partition(",")
                if not sep:
                    name, sep, tail = rest.partition(">")
                    sep = ">"
                out.write(f"##contig=<ID={alias.encode(name)}{sep}{tail}\n")
                continue
            if line.startswith("#") or not line.strip():
                out.write(line + "\n")
                continue
            cols = line.split("\t")
            cols[0] = alias.encode(cols[0])
            out.write("\t".join(cols) + "\n")
            n += 1
    return n


def alias_fasta(src: str | Path, dst: str | Path, alias: ContigAlias) -> int:
    """Rewrite a FASTA's headers to aliases and index the result.

    Needed only on the `donor_outcome.precall` DRIVER path, which loads a FASTA by
    `display_id` (the first whitespace token of the header). `precall_sample` takes
    the sequence string directly and so needs no FASTA aliasing.
    """
    import pysam

    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with src.open() as fh, dst.open("w", newline="\n") as out:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0] if line[1:].strip() else ""
                out.write(f">{alias.encode(name)}\n")
                n += 1
            else:
                out.write(line)
    pysam.faidx(str(dst))
    return n


def decode_outcome(outcome: Any, alias: ContigAlias) -> Any:
    """Return a copy of a TRACE `SampleOutcome` with real contig names restored.

    Every variant-id list is mapped through `ContigAlias.decode_variant_id`, which
    parses right-anchored and understands the `(LC)` less-confident prefix. An id
    that does not have the expected shape is returned unchanged rather than mangled.
    """
    return dataclasses.replace(
        outcome,
        designed_variants=[alias.decode_variant_id(v) for v in outcome.designed_variants],
        synthetic_errors=[alias.decode_variant_id(v) for v in outcome.synthetic_errors],
        spontaneous_variants=[
            alias.decode_variant_id(v) for v in outcome.spontaneous_variants
        ],
    )


@dataclass
class AliasedPrecall:
    """A precall run through the alias adapter, with the mapping kept for provenance.

    Attributes
    ----------
    outcome : SampleOutcome
        TRACE's own dataclass, with **decoded** (real) contig names.
    alias : ContigAlias
        The bijection used. `alias.is_noop` is True when the reference needed no
        aliasing, in which case the outputs are byte-identical to an unaliased run.
    n_vcf_records, n_gvcf_records : int
        Records rewritten on the way in, so the re-representation has a denominator.
    """

    outcome: Any
    alias: ContigAlias
    n_vcf_records: int = 0
    n_gvcf_records: int = 0
    aliased_contigs: tuple[str, ...] = ()


def precall_sample_aliased(
    spec: Any,
    vcf_path: str | Path,
    gvcf_path: str | Path,
    chrom_seq: str,
    *,
    contigs: Iterable[str],
    workdir: str | Path | None = None,
    **precall_kwargs: Any,
) -> AliasedPrecall:
    """`precall_sample`, made safe for underscore-bearing contigs.

    Parameters
    ----------
    spec : DonorSpec
        TRACE's spec. Its `chrom` is encoded before the call and never mutated.
    vcf_path, gvcf_path : path
        Rewritten to alias-named copies beside the originals (or under `workdir`).
        The originals are left untouched.
    chrom_seq : str
        The donor contig's sequence, passed through unchanged — it carries no name.
    contigs : iterable of str
        **Every** contig name in the reference, not just this sample's. The alias
        must be built over the whole reference or two runs could assign the same
        alias to different contigs and the mapping would not be comparable.

    Returns
    -------
    AliasedPrecall
    """
    do = _trace.donor_outcome()
    names = list(dict.fromkeys(list(contigs) + [spec.chrom]))
    alias = ContigAlias.from_names(names)

    vcf_path, gvcf_path = Path(vcf_path), Path(gvcf_path)
    base = Path(workdir) if workdir is not None else vcf_path.parent
    if alias.is_noop:
        a_vcf, a_gvcf, n_v, n_g = vcf_path, gvcf_path, 0, 0
    else:
        a_vcf = base / (vcf_path.name + ".aliased.vcf")
        a_gvcf = base / (gvcf_path.name + ".aliased.vcf")
        n_v = alias_vcf_contigs(vcf_path, a_vcf, alias)
        n_g = alias_vcf_contigs(gvcf_path, a_gvcf, alias)

    aliased_spec = dataclasses.replace(spec, chrom=alias.encode(spec.chrom))
    outcome = do.precall_sample(aliased_spec, a_vcf, a_gvcf, chrom_seq, **precall_kwargs)
    return AliasedPrecall(
        outcome=decode_outcome(outcome, alias),
        alias=alias,
        n_vcf_records=n_v,
        n_gvcf_records=n_g,
        aliased_contigs=alias.aliased,
    )


# --------------------------------------------------------------------------- #
# Defect 8 — count what the donor-window filter throws away                    #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class VcfWindowAudit:
    """What `precall_sample`'s donor-window filter keeps and discards.

    Attributes
    ----------
    n_records : int
        Data records in the VCF. The denominator.
    n_kept : int
        Records inside the window AS CONFIGURED (donor window widened by `flank`).
        With `flank=0` this is exactly what TRACE buckets.
    n_kept_donor_window : int
        Records TRACE would bucket, whatever `flank` is — the parity number.
    n_discarded : int
        `n_records - n_kept`: records outside the configured window. With `flank=0`,
        **every one of these is invisible today.**
    n_distal_same_contig : int
        Records on the donor contig outside the DONOR window, whatever `flank` is —
        candidate distal unintended edits.
    n_other_contig : int
        Records on another contig (never in any flank of this window).
    n_in_flank : int
        Records outside the donor window but inside `flank`, i.e. the ones the wider
        window recovers. `n_kept == n_kept_donor_window + n_in_flank`.
    discarded_ids : tuple of str
        `{chrom}_{pos}_{ref}_to_{alt}` for each, so they can be inspected rather
        than merely counted. Capped by `max_ids` to bound memory on a whole-genome
        VCF; `n_discarded` is always exact.
    """

    n_records: int = 0
    n_kept: int = 0
    n_kept_donor_window: int = 0
    n_discarded: int = 0
    n_distal_same_contig: int = 0
    n_other_contig: int = 0
    n_in_flank: int = 0
    n_ids_truncated: int = 0
    flank: int = 0
    discarded_ids: tuple[str, ...] = field(default_factory=tuple)

    def summary(self) -> str:
        return (
            f"{self.n_kept_donor_window} of {self.n_records} VCF record(s) fall in the "
            f"donor window ({self.n_kept} within flank={self.flank}); "
            f"{self.n_discarded} discarded — {self.n_distal_same_contig} distal on the "
            f"donor contig, {self.n_other_contig} on another contig, {self.n_in_flank} "
            "recovered by the configured flank"
        )


def audit_vcf_window(
    vcf_path: str | Path,
    *,
    chrom: str,
    donor_start: int,
    donor_end: int,
    flank: int = 0,
    max_ids: int = 1000,
) -> VcfWindowAudit:
    """Count what the donor-window filter discards, reproducing TRACE's arithmetic.

    TRACE keeps a record when, with `offset = donor_start - 1`,
    `wt_len = donor_end - donor_start + 1` and `pos0 = pos - offset - 1`::

        v_chr == chrom and pos0 >= 0 and (pos0 + len(ref)) < wt_len

    (The `>=` on the right edge makes the final window base uncallable — NN 243's A7,
    a separate TRACE-side one-character fix `>= -> >`. It is reproduced here, not
    corrected, because this function's job is to report TRACE's behaviour.)

    `flank=0` reproduces TRACE exactly. A non-zero `flank` WIDENS the window that
    counts as kept and reports the recovered records separately as `n_in_flank`, so a
    screen can both see how many variants a wider window recovers and choose to widen
    it — without patching TRACE. `n_kept_donor_window` always stays the parity number.
    """
    offset = donor_start - 1
    wt_len = donor_end - donor_start + 1
    lo, hi = donor_start - flank, donor_end + flank
    n = kept_donor = distal = other = in_flank = truncated = 0
    ids: list[str] = []
    path = Path(vcf_path)
    with _open_text(path) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            n += 1
            v_chr, v_ref, v_alt = cols[0], cols[3], cols[4]
            try:
                v_pos = int(cols[1])
            except ValueError:
                continue
            pos0 = v_pos - offset - 1
            if v_chr == chrom and pos0 >= 0 and (pos0 + len(v_ref)) < wt_len:
                kept_donor += 1
                continue
            recovered = False
            if v_chr == chrom:
                distal += 1
                if flank and lo <= v_pos <= hi:
                    in_flank += 1
                    recovered = True
            else:
                other += 1
            if recovered:
                continue
            if len(ids) < max_ids:
                ids.append(f"{v_chr}_{v_pos}_{v_ref}_to_{v_alt}")
            else:
                truncated += 1
    kept = kept_donor + in_flank
    return VcfWindowAudit(
        n_records=n, n_kept=kept, n_kept_donor_window=kept_donor,
        n_discarded=n - kept, n_distal_same_contig=distal, n_other_contig=other,
        n_in_flank=in_flank, n_ids_truncated=truncated, flank=flank,
        discarded_ids=tuple(ids),
    )
