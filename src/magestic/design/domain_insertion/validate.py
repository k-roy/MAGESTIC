"""
Independent validator for domain-insertion guide-donor oligos.

This deliberately does NOT import the designer. It takes a finished oligo, parses it back into
its parts, simulates the step-3 PaqCI digest and LOV ligation, and asserts that the resulting
allele is what was intended. That independence is the point: a designer bug and a validator bug
would have to agree to pass, which is far less likely than either alone being wrong.

Checks (see `validate_oligo`):
  C1  exactly one BspQI site, in the expected cloning position
  C2  exactly one inward-facing PaqCI pair, and no other PaqCI site
  C3  no AscI/NotI internal to the donor
  C4  the donor starts on a codon boundary (design consideration 4)
  C5  both arms match the reference, allowing only synonymous designed changes
  C6  recoding is translation-neutral (design consideration 2)
  C7  simulated PaqCI digest + ligation gives the native protein with the domain inserted,
      in frame, at exactly the intended residue
  C8  the guide does not re-cut the finished donor
"""

from dataclasses import dataclass, field

from .spec import (
    ASLOV2_AA,
    ASLOV2_CDS,
    ASCI_SITE,
    BSPQI_SITE,
    LEFT_OVERHANG,
    NOTI_SITE,
    PAQCI_CUT_BOTTOM,
    PAQCI_CUT_TOP,
    PAQCI_SITE,
    RIGHT_OVERHANG,
    OligoSpec,
)

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTT"
        "NNKKSSRRVVVVAAAADDEEGGGG")
CODON_TABLE = {a + b + c: _AAS[i]
               for i, (a, b, c) in enumerate((x, y, z) for x in _BASES
                                             for y in _BASES for z in _BASES)}


def translate(seq: str, strict: bool = True) -> str:
    seq = seq.upper()
    if strict and len(seq) % 3:
        raise ValueError(f"sequence length {len(seq)} is not a multiple of 3")
    out = []
    for i in range(0, len(seq) - 2, 3):
        out.append(CODON_TABLE.get(seq[i:i + 3], "X"))
    return "".join(out)


def find_all(hay: str, needle: str):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


@dataclass
class PaqCISites:
    """0-based indices into the donor, and the resulting cut/overhang coordinates."""

    bottom_site: int          # index of the rc(PaqCI) match (left, cuts leftward)
    top_site: int             # index of the PaqCI match (right, cuts rightward)
    left_cut_top: int         # top strand cut; donor[:left_cut_top] is retained
    right_cut_top: int        # top strand resumes at donor[right_cut_top:]
    left_overhang: str
    right_overhang: str


def locate_paqci(donor: str) -> PaqCISites:
    """
    Find the inward-facing PaqCI pair and compute cut coordinates.

    PaqCI is CACCTGC(4/8). For the top-strand site the top cut falls 4 nt past the end of the
    recognition sequence. For the bottom-strand site (rc = GCAGGTG) the cutting runs toward
    lower coordinates: the top-strand cut falls 8 nt before the start of the match.
    """
    top_hits = find_all(donor, PAQCI_SITE)
    bot_hits = find_all(donor, rc(PAQCI_SITE))
    if len(top_hits) != 1 or len(bot_hits) != 1:
        raise ValueError(
            f"expected exactly one PaqCI site per strand, got "
            f"{len(top_hits)} top and {len(bot_hits)} bottom")

    bot, top = bot_hits[0], top_hits[0]
    if bot >= top:
        raise ValueError("PaqCI sites are not in the expected inward-facing arrangement")

    left_cut_top = bot - PAQCI_CUT_BOTTOM
    right_cut_top = top + len(PAQCI_SITE) + PAQCI_CUT_TOP
    if left_cut_top < 0:
        raise ValueError("left PaqCI cut falls outside the donor")

    return PaqCISites(
        bottom_site=bot,
        top_site=top,
        left_cut_top=left_cut_top,
        right_cut_top=right_cut_top,
        left_overhang=donor[left_cut_top:left_cut_top + 4],
        right_overhang=donor[right_cut_top:right_cut_top + 4],
    )


@dataclass
class ValidationResult:
    ok: bool
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    info: dict = field(default_factory=dict)

    def fail(self, code: str, msg: str):
        self.ok = False
        self.errors.append(f"{code}: {msg}")

    def warn(self, code: str, msg: str):
        self.warnings.append(f"{code}: {msg}")


def extract_donor(oligo: str, spec: OligoSpec) -> str:
    """Donor = everything between the internal cloning site and the 3' priming sequence."""
    i = oligo.upper().find(spec.cloning_site.upper())
    if i == -1:
        raise ValueError("internal cloning site not found in oligo")
    start = i + len(spec.cloning_site)
    if spec.rev_primer:
        j = oligo.upper().rfind(rc(spec.rev_primer).upper())
        if j == -1:
            j = oligo.upper().rfind(spec.rev_primer.upper())
        if j != -1 and j > start:
            return oligo[start:j]
    return oligo[start:]


def extract_guide(oligo: str, spec: OligoSpec) -> str:
    i = oligo.upper().find(spec.cloning_site.upper())
    if i == -1:
        raise ValueError("internal cloning site not found in oligo")
    if i < spec.guide_length:
        raise ValueError("oligo too short upstream of the cloning site to hold a guide")
    return oligo[i - spec.guide_length:i]


def _splice(arm: str, cds_positions=None) -> str:
    """
    Drop intronic nucleotides from a genomic arm before translating it.

    Design consideration 5: the donor arm must MATCH THE GENOME (it is a repair template, so it
    carries introns), but its CODONS live in the spliced CDS. Translating a genomic arm that
    spans an intron produces nonsense -- which is exactly how the two NBL1 (YHR199C-A) failures
    surfaced during scoping. `cds_positions` is the CDS coordinate of each arm nucleotide, or
    None where intronic / outside the CDS.
    """
    if cds_positions is None:
        return arm
    return "".join(b for b, cp in zip(arm, cds_positions) if cp is not None)


def validate_oligo(
    oligo: str,
    spec: OligoSpec,
    native_protein: str = None,
    insertion_after_residue: int = None,
    lov_cds: str = ASLOV2_CDS,
    lov_aa: str = ASLOV2_AA,
    left_cds_positions=None,
    right_cds_positions=None,
) -> ValidationResult:
    """
    Validate one designed oligo.

    `native_protein` is the full unmodified protein of the target gene; `insertion_after_residue`
    is the 1-based residue after which the domain is inserted. When both are supplied, C7 checks
    the reconstructed allele against the real gene rather than merely against itself.
    """
    r = ValidationResult(ok=True)
    oligo = oligo.upper()
    spec.validate()

    # ---- C1: the BspQI sites the spec calls for -------------------------
    # Count comes from `spec.expected_bspqi_sites` (derived from the cloning site) rather than being
    # hard-coded to 1, so the two-BspQI variant Alex uses is representable. Kevin's current default
    # is ONE site (his decision, 2026-07-26); the two-site cassette is under evaluation.
    _want_bspqi = spec.expected_bspqi_sites
    bspqi = find_all(oligo, BSPQI_SITE) + find_all(oligo, rc(BSPQI_SITE))
    if len(bspqi) != _want_bspqi:
        r.fail("C1", f"expected exactly {_want_bspqi} BspQI site(s), "
                     f"found {len(bspqi)} at {sorted(bspqi)}")

    donor = extract_donor(oligo, spec)
    guide = extract_guide(oligo, spec)
    r.info["guide"] = guide
    r.info["donor_length"] = len(donor)

    # ---- C2: exactly one inward PaqCI pair -----------------------------
    try:
        sites = locate_paqci(donor)
        r.info["left_overhang"] = sites.left_overhang
        r.info["right_overhang"] = sites.right_overhang
        if sites.left_overhang != spec.left_overhang:
            r.fail("C2", f"left overhang {sites.left_overhang!r} != {spec.left_overhang!r}; "
                         "a single pooled LOV insert will not ligate")
        if sites.right_overhang != spec.right_overhang:
            r.fail("C2", f"right overhang {sites.right_overhang!r} != {spec.right_overhang!r}; "
                         "a single pooled LOV insert will not ligate")
    except ValueError as e:
        r.fail("C2", str(e))
        return r

    # ---- C3: no AscI / NotI inside the donor ---------------------------
    for name, site in (("AscI", ASCI_SITE), ("NotI", NOTI_SITE)):
        if find_all(donor, site) or find_all(donor, rc(site)):
            r.fail("C3", f"{name} site present inside the donor")

    # ---- C4: donor starts on a codon boundary --------------------------
    left_arm_plus_linker = donor[:sites.left_cut_top] + sites.left_overhang
    if left_cds_positions is None and len(left_arm_plus_linker) % 3:
        r.fail("C4", f"donor start is not codon-phased: {len(left_arm_plus_linker)} nt "
                     "from donor start to the end of the left linker is not a multiple of 3")

    # ---- C7: simulate step 3 and rebuild the final allele --------------
    insert_top = spec.insert_top_strand(lov_cds)
    if not insert_top.startswith(spec.left_overhang):
        r.fail("C7", "LOV insert does not begin with the left overhang")
    final_donor = donor[:sites.left_cut_top] + insert_top + donor[sites.right_cut_top:]
    r.info["final_donor_length"] = len(final_donor)

    # The finished donor is GENOMIC: it may carry introns and UTR. Translate the CODING portion.
    # Without CDS annotation we fall back to translating the donor directly, which is correct
    # only for an intronless arm that stays inside the CDS.
    if left_cds_positions is not None or right_cds_positions is not None:
        _left_coding = _splice(donor[:sites.left_cut_top - (len(spec.left_linker) - 4)],
                               left_cds_positions)
        _right_coding = _splice(donor[sites.right_cut_top + len(spec.right_linker):],
                                right_cds_positions)
        # The arm may begin mid-codon: once the donor is allowed to start in 5' UTR or inside an
        # intron, the first codon it touches is usually incomplete (its remaining bases lie
        # outside the donor). Drop that leading partial codon before translating -- the coding
        # portion still ENDS on a codon boundary, because the insertion point is one.
        _lead = len(_left_coding) % 3
        _left_coding = _left_coding[_lead:]
        coding = (_left_coding + spec.left_linker + spec.amplicon_linker_left + lov_cds
                  + spec.amplicon_linker_right + spec.right_linker + _right_coding)
    else:
        coding = final_donor

    trimmed = coding[:len(coding) - (len(coding) % 3)]
    if len(coding) % 3:
        r.warn("C7", f"coding length {len(coding)} is not a whole number of codons; "
                     "translating the in-frame prefix only")
    edited_aa = translate(trimmed)
    r.info["edited_protein"] = edited_aa

    # Stop-codon handling is deferred until we know how much of the right arm lies past the
    # end of the CDS: for a C-terminal insertion the arm legitimately spans the native stop and
    # continues into the 3' UTR, so a '*' there is correct, not premature. See below.
    _stop_check_pending = True

    if lov_aa not in edited_aa:
        r.fail("C7", "AsLOV2 protein sequence is not present in frame in the edited allele")
    else:
        r.info["lov_offset_in_donor_aa"] = edited_aa.index(lov_aa) + 1

    # ---- C5/C6: arms against the real protein --------------------------
    #
    # The left arm runs from the donor start to the beginning of the left linker. The linker's
    # last 4 nt are the PaqCI overhang, so the linker begins at left_cut_top - (len(linker) - 4).
    if native_protein and insertion_after_residue:
        # left_cut_top falls INSIDE the left linker: the overhang is the linker's last 4 nt,
        # so the arm ends (len(linker) - 4) nt before the cut.
        left_arm_end = sites.left_cut_top - (len(spec.left_linker) - 4)
        # right_cut_top is where the right linker STARTS (its first 4 nt are the overhang),
        # so the right arm begins a full linker later -- not (len - 4) later.
        right_arm_start = sites.right_cut_top + len(spec.right_linker)
        left_arm = donor[:left_arm_end]
        right_arm = donor[right_arm_start:]
        r.info["left_arm_nt"] = len(left_arm)
        r.info["right_arm_nt"] = len(right_arm)

        # An arm may legitimately run off the end of the CDS: into the 5' UTR for insertions
        # near the N-terminus, or through the stop codon near the C-terminus. Kevin's
        # consideration 4 says so explicitly ("only applies when we are already ~20+ codons
        # into the CDS"). So compare only the portion that overlaps the coding sequence.
        left_coding = _splice(left_arm, left_cds_positions)
        right_coding = _splice(right_arm, right_cds_positions)
        if left_cds_positions is not None:
            r.info["left_arm_intronic_nt"] = len(left_arm) - len(left_coding)
            r.info["right_arm_intronic_nt"] = len(right_arm) - len(right_coding)

        left_aa = ""
        lead = len(left_coding) % 3          # partial leading codon; see note above
        if lead:
            r.info["left_arm_partial_leading_codon_nt"] = lead
            left_coding = left_coding[lead:]
        if len(left_coding) % 3:
            r.fail("C4", f"left arm coding portion is {len(left_coding)} nt, "
                         "not a whole number of codons")
        else:
            left_aa = translate(left_coding)
            avail = min(len(left_aa), insertion_after_residue)
            expected_left = native_protein[insertion_after_residue - avail:insertion_after_residue]
            if left_aa[len(left_aa) - avail:] != expected_left:
                r.fail("C5", f"left arm translates to {left_aa!r} but the native protein has "
                             f"{expected_left!r} - recoding is not translation-neutral, or the "
                             "arm is placed at the wrong residue")
            if avail < len(left_aa):
                r.warn("C5", f"left arm extends {len(left_aa) - avail} codons upstream of the "
                             "start codon (5' UTR) - expected for N-terminal insertions")
            r.info["left_arm_utr_codons"] = len(left_aa) - avail

        right_aa = translate(right_coding[:len(right_coding) - (len(right_coding) % 3)])
        avail_r = min(len(right_aa), len(native_protein) - insertion_after_residue)
        expected_right = native_protein[insertion_after_residue:
                                        insertion_after_residue + avail_r]
        if right_aa[:avail_r] != expected_right:
            r.fail("C5", f"right arm translates to {right_aa!r} but the native protein has "
                         f"{expected_right!r}")
        if avail_r < len(right_aa):
            r.warn("C5", f"right arm extends {len(right_aa) - avail_r} codons past the last "
                         "residue (stop codon / 3' UTR) - expected for C-terminal insertions")
        r.info["right_arm_utr_codons"] = len(right_aa) - avail_r

        # C7 (strict form): the finished allele must read
        # [left arm][GS][ampL][AsLOV2][ampR][SG][right arm] contiguously.
        # Built from the arms themselves (already anchored to the native protein by C5 above),
        # so this specifically tests frame and junction assembly rather than re-testing C5.
        want = left_aa \
            + translate(spec.left_linker) \
            + translate(spec.amplicon_linker_left) \
            + lov_aa \
            + translate(spec.amplicon_linker_right) \
            + translate(spec.right_linker) \
            + right_aa
        if want and want != edited_aa:
            r.fail("C7", "edited allele does not match the intended "
                         "[arm][GS][linker][AsLOV2][linker][SG][arm] product")

        # The junction assembly is the part C5 cannot see: AsLOV2 must be immediately flanked
        # by the linker residues, in frame, with nothing inserted or lost at either seam.
        want_left_flank = translate(spec.left_linker) + translate(spec.amplicon_linker_left)
        want_right_flank = translate(spec.amplicon_linker_right) + translate(spec.right_linker)
        i = edited_aa.find(lov_aa)
        if i < 0:
            r.fail("C7", "AsLOV2 is not present in frame in the edited allele")
        else:
            if edited_aa[i - len(want_left_flank):i] != want_left_flank:
                r.fail("C7", f"AsLOV2 is not preceded by {want_left_flank!r} - "
                             "the left PaqCI junction is malformed")
            if edited_aa[i + len(lov_aa):i + len(lov_aa) + len(want_right_flank)] \
                    != want_right_flank:
                r.fail("C7", f"AsLOV2 is not followed by {want_right_flank!r} - "
                             "the right PaqCI junction is malformed")

        r.info["insertion_after_residue"] = insertion_after_residue
        r.info["native_context"] = (native_protein[insertion_after_residue - 8:insertion_after_residue]
                                    + "^" + native_protein[insertion_after_residue:
                                                           insertion_after_residue + 8])

    # ---- C7 (deferred): internal stop codons ---------------------------
    # Only the coding portion may not contain a stop. Codons past the end of the CDS are the
    # native stop plus 3' UTR, carried because the right arm is genomic.
    if _stop_check_pending:
        utr = r.info.get("right_arm_utr_codons", 0)
        coding_part = edited_aa[:len(edited_aa) - utr] if utr else edited_aa
        if "*" in coding_part:
            r.fail("C7", f"premature stop codon in the coding portion of the edited allele "
                         f"at residue {coding_part.index('*') + 1}")
        elif utr:
            r.warn("C7", f"right arm extends {utr} codons past the last residue and spans the "
                         "native stop codon - expected for a C-terminal insertion")

    # ---- C8: guide must not re-cut the finished donor ------------------
    if guide:
        for hay, label in ((final_donor, "final donor"), (donor, "oligo donor")):
            for strand_seq, strand in ((hay, "+"), (rc(hay), "-")):
                pos = strand_seq.find(guide)
                while pos != -1:
                    pam = strand_seq[pos + len(guide):pos + len(guide) + 3]
                    if len(pam) == 3 and pam[1:] == "GG":
                        r.fail("C8", f"guide + NGG PAM still present in the {label} "
                                     f"({strand} strand) - the donor would be re-cut")
                    pos = strand_seq.find(guide, pos + 1)

    r.info["residues_added"] = spec.residues_added(lov_aa)
    return r
