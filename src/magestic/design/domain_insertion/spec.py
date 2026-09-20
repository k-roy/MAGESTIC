"""
Oligo architecture for MAGESTIC domain-insertion libraries (AsLOV2 into essential genes).

Every constant here was decoded from, and is regression-tested against, the hand-designed
reference oligo `pol3-example-oligo-250-nt.gb` (Wilson lab / K. Roy, 2026-07).

Layout of the 250-nt designed oligo
-----------------------------------
    [5' const ..AscI] [guide 20] [GTTTGAAGAGC] [left arm] [GS] [PaqCI stuffer] [SG] [right arm] [rev primer]
             20            20          11         ~66       6         38         6      ~63          20

Fixed overhead is 121 nt, leaving 129 nt of homology arms at a 250-nt oligo.

Why the GS / SG linkers matter
------------------------------
The two linker codons are NOT native gene sequence -- verified at the protein level against
SGD R64 (POL3/YDL102W: the native left and right arms are contiguous at residue 773, and the
`...KGS|SGHP...` junction does not occur in POL3). They are added by the design.

That is what makes a *pooled* third cloning step possible: PaqCI cuts inside the linker codons,
so the 4-nt overhangs are CONSTANT across the whole library --

    left overhang  TTCT   (last 4 nt of GGT TCT)
    right overhang TCCG   (first 4 nt of TCC GGA)

Had the overhangs been gene-derived, one LOV insert could not ligate into ~580k members.

Any ADDITIONAL linker length is carried on the AsLOV2 PCR amplicon, added by the amplification
primers -- so it costs zero oligo budget (Kevin, 2026-07-24). Final protein product:

    ...left arm | G S | [linkerL] AsLOV2(136 aa) [linkerR] | S G | right arm...
                 └─ oligo ─┘      └──── LOV amplicon ────┘   └ oligo ┘
"""

from dataclasses import dataclass

# --- type IIS enzymes -------------------------------------------------------
# BspQI: step-1 guide-donor cloning. PaqCI: step-3 domain insertion.
BSPQI_SITE = "GCTCTTC"
BSPQI_CUT_TOP, BSPQI_CUT_BOTTOM = 1, 4      # GCTCTTC(1/4), 3-nt 5' overhang

PAQCI_SITE = "CACCTGC"
PAQCI_CUT_TOP, PAQCI_CUT_BOTTOM = 4, 8      # CACCTGC(4/8), 4-nt 5' overhang

# Sites that must NOT occur inside a donor arm.
#
# NOTE: the legacy saturation-editing code screens only
#   RESTRICTION_SITES_DISALLOWED = ['GGCGCGCC', 'GCGGCCGC', 'GCTCTTC']
# which OMITS PaqCI. A donor arm containing CACCTGC would be cut during step 3 and the
# member lost. Adding PaqCI here is the single most important change for this library.
ASCI_SITE = "GGCGCGCC"
NOTI_SITE = "GCGGCCGC"
DISALLOWED_IN_DONOR = (BSPQI_SITE, PAQCI_SITE, ASCI_SITE, NOTI_SITE)

# --- fixed oligo elements ---------------------------------------------------
INTERNAL_CLONING_SITE = "GTTTGAAGAGC"       # scaffold GTTT + BspQI (rev-comp orientation)

# --- OPTIONAL two-BspQI variant (Kevin, 2026-07-27; Alex's rationale, 2026-07-27) --------------
# Alex: "we added the second BSPQI site for the reason you mentioned of the insert potentially being
# inserted head to tail, and it has always been our least efficient Gibson so we wanted to add the
# site to see if that would help there also."
#
# Transcribed from Alex's Benchling map (Kevin, 2026-07-27). The segment between the guide and the
# upstream homology reads, 5'->3':
#
#     GTTTGAAGAGC  GTAGTCAGATGC  GCTCTTC  GCCTGCA
#     `----------' `----------'  `-----'  `-----'
#      scaffold +   spacer (12)   BspQI    cut(1) + 3-nt overhang, + 3 trailing
#      rc-BspQI                   site 2
#
# 37 nt total, i.e. **+26 nt** over the single-site cassette -- NOT the +9 nt of my first
# provisional guess, which used a 2-nt spacer. That guess was wrong and unrealistic: two type-IIS
# sites 2 nt apart cannot both be bound and cut efficiently, and Alex's 12-nt spacer is the real
# constraint. The corrected cost changes the conclusion (see planning/28).
#
# Kevin, 2026-07-27: "I don't think the terminal gca is needed as it's not involved in the overlaps.
# Perhaps that was homology that he mislabeled as stuffer." The cut geometry agrees -- BspQI
# GCTCTTC(1/4) needs only 4 nt after the site (1 nt spacer + the 3-nt overhang) -- so GCCTGCA can be
# trimmed to GCCT, giving 34 nt.
#
# 34 leaves arm_budget 106, which is NOT divisible by 3 (see `validate`), so one inert base is added
# to the SPACER (the safe place to pad) to reach 35 -> arm_budget 105 -> arms 54/51.
ALEX_BSPQI_CLONING_SITE_AS_DRAWN = "GTTTGAAGAGC" + "GTAGTCAGATGC" + BSPQI_SITE + "GCCTGCA"   # 37
# ⛔ BLOCKER FOUND AND FIXED (adversary panel + verification, 2026-07-27).
# The trailing 4 nt after the second BspQI site DEFINE THE VECTOR'S DONOR-SIDE OVERHANG, because
# BspQI is GCTCTTC(1/4): the 3-nt 5' protrusion is bases +2..+4 after the site.
#   `GCCT`  -> protrusion CCT.  cV451's insert right end presents 5'-AAC-3' (bottom strand), which
#              pairs with GTT, NOT CCT.  A fully digested two-site vector COULD NOT RECEIVE cV451.
#   `GGTT`  -> protrusion GTT.  Matches, and is the same overhang the validated single-site
#              cassette produces (verified: single-site vector end = 5'-GTT-3' top / AAC bottom,
#              identical to cV451's own ends -- which is also why cV451 can concatenate).
# Length is unchanged at 35 nt, so the arm-budget arithmetic in `validate()` is untouched.
# STILL UNCONFIRMED against Alex's Benchling map, which has never been seen -- the map's trailing
# `GCCTGCA` was transcribed from a screenshot and its last 3 nt were already suspected to be
# homology mislabelled as stuffer (Kevin, 2026-07-27).
DOUBLE_BSPQI_CLONING_SITE = "GTTTGAAGAGC" + "GTAGTCAGATGCA" + BSPQI_SITE + "GGTT"            # 35

# DECISION (Kevin, 2026-07-25): the linkers stay as they are. `GGTTCT` begins with GG, so a guide
# whose protospacer ends exactly at the arm/linker boundary is handed a functional NGG PAM by the
# linker (PAM = arm[-1] + linker[0] + linker[1]). No codon choice avoids this -- every glycine codon
# is GG* -- and swapping to Ser-Gly would change the left overhang TTCT -> TGGT and force a redesign
# of the AsLOV2 amplicon primers. Measured across 1,200 positions, the swap changed designability,
# validator cleanliness and self-disrupting fraction by nothing at all (artifact 14).
#
# Instead the hazard is blocked where it arises: such guides are NOT treated as self-disrupting
# (splitting the PAM is not disruption -- the protospacer survives intact), so they fall through to
# normal SPACER SYNONYMOUS RECODING and must reach the disruption threshold like any other guide.
# Validator check C8 confirms per design that the finished donor carries no protospacer + NGG.
# Measured: 9.8% of designs take this path, median 8 codons recoded, median disruption score 6.
LEFT_LINKER = "GGTTCT"                      # Gly-Ser, immediately 5' of the stuffer
RIGHT_LINKER = "TCCGGA"                     # Ser-Gly, immediately 3' of the stuffer

LEFT_OVERHANG = "TTCT"                      # LEFT_LINKER[-4:]
RIGHT_OVERHANG = "TCCG"                     # RIGHT_LINKER[:4]

# The 38-nt stuffer excised in step 3. Carries an inward-facing PaqCI pair:
# bottom-strand site at stuffer offsets 5..11, top-strand site at offsets 28..34.
STUFFER = "gagtgcaggtggtagtgctaacaggatcacctgcgcta".upper()

GUIDE_LENGTH = 20
PAM = "NGG"                                 # SpCas9, per project decision (robust WT enzyme)

# --- payload ----------------------------------------------------------------
# AsLOV2 (A. sativa phototropin 1, residues ~404-546), 408 nt / 136 aa.
# Byte-identical to CDS 3922..4329 of pKR497 (ADE2-disrupt donor).
# Verified free of BspQI, PaqCI, AscI and NotI sites.
ASLOV2_CDS = (
    "TTGGAACGTATCGAAAAGAATTTTGTCATCACGGATCCGCGTCTTCCCGACAATCCGATTATCTTCGCGTCAGAC"
    "TCTTTCTTACAACTGACTGAGTATAGTAGAGAGGAGATATTGGGGCGTAACTGTAGATTTCTTCAGGGGCCAGAA"
    "ACTGATCGGGCTACCGTTCGCAAGATACGTGACGCAATAGACAACCAGACCGAGGTGACGGTGCAGCTGATTAAC"
    "TACACAAAGTCTGGGAAGAAGTTCTGGAACCTGTTTCATTTACAACCTATGAGAGACCAAAAAGGTGACGTTCAA"
    "TATTTCATCGGGGTTCAGTTAGATGGGACTGAGCACGTGAGAGATGCAGCAGAAAGAGAGGGTGTAATGCTTATT"
    "AAAAAAACAGCCGAGAATATCGACGAAGCCGCT"
)
ASLOV2_AA = (
    "LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRATVRKIRDAIDNQTEVTVQLINY"
    "TKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLDGTEHVRDAAEREGVMLIKKTAENIDEAA"
)


@dataclass(frozen=True)
class OligoSpec:
    """Tunable design parameters. Defaults reproduce the POL3 reference oligo."""

    oligo_length: int = 250
    guide_length: int = GUIDE_LENGTH
    fwd_const: str = ""                     # 5' constant ending in AscI; set per library
    rev_primer: str = ""                    # 3' subpool priming sequence
    left_linker: str = LEFT_LINKER
    right_linker: str = RIGHT_LINKER
    stuffer: str = STUFFER
    cloning_site: str = INTERNAL_CLONING_SITE
    min_arm: int = 20                       # MINIMUM_HOMOLOGY in the legacy code
    left_arm_must_be_codon_aligned: bool = True   # donor STARTS on a codon boundary (readability)
    # Additional linker carried on the LOV amplicon (must be a codon multiple).
    amplicon_linker_left: str = ""
    amplicon_linker_right: str = ""

    @property
    def left_overhang(self) -> str:
        """PaqCI leaves the linker's last 4 nt as the left 5' overhang."""
        return self.left_linker[-4:].upper()

    @property
    def right_overhang(self) -> str:
        """...and the right linker's first 4 nt as the right overhang."""
        return self.right_linker[:4].upper()

    @property
    def creates_junction_pam(self) -> bool:
        """
        Does the LEFT linker hand a functional NGG PAM to a guide whose protospacer ends exactly
        at the arm/linker boundary?

        The reconstituted PAM spans the junction as [last arm base][linker[0]][linker[1]], so the
        condition is simply `left_linker[:2] == "GG"`. Every glycine codon is GG*, so a
        Gly-FIRST linker always triggers it -- no codon choice avoids it, only a different
        residue order. See artifact 14.
        """
        return self.left_linker[:2].upper() == "GG"

    @property
    def fixed_overhead(self) -> int:
        return (len(self.fwd_const) + self.guide_length + len(self.cloning_site)
                + len(self.left_linker) + len(self.stuffer) + len(self.right_linker)
                + len(self.rev_primer))

    @property
    def arm_budget(self) -> int:
        """Total nt available for the two homology arms."""
        return self.oligo_length - self.fixed_overhead

    def split_arms(self) -> tuple:
        """
        (left, right) arm lengths, as close to even as the frame allows.

        The LEFT arm length must be a whole number of codons: the insertion point sits on a
        codon boundary, so the donor only starts on a codon boundary (design consideration 4)
        if the left arm spans complete codons. The right arm takes whatever is left -- it
        starts on a boundary regardless, so its length is unconstrained.

        With the reference 250-nt oligo this yields (66, 63), matching the hand design exactly.
        """
        b = self.arm_budget
        left = 3 * round(b / 6)
        left = min(left, (b // 3) * 3)
        return left, b - left

    def insert_top_strand(self, lov_cds: str = ASLOV2_CDS) -> str:
        """
        Top strand contributed by the ligated LOV amplicon, from the left overhang up to
        (but not including) the right overhang -- i.e. exactly what replaces the excised
        stuffer region between the two PaqCI cuts.
        """
        return (self.left_overhang + self.amplicon_linker_left + lov_cds
                + self.amplicon_linker_right)

    def residues_added(self, lov_aa: str = ASLOV2_AA) -> int:
        """Number of residues the finished edit inserts into the target protein."""
        return (len(self.left_linker) // 3 + len(self.amplicon_linker_left) // 3
                + len(lov_aa)
                + len(self.amplicon_linker_right) // 3 + len(self.right_linker) // 3)

    @property
    def expected_bspqi_sites(self) -> int:
        """How many BspQI sites the assembled oligo should carry -- read from the cloning site.

        Previously hard-coded as 1 in `_oligo_sites_are_canonical` and validator C1, which made the
        two-site variant unrepresentable.
        """
        from .genome import rc as _rc
        cs = self.cloning_site.upper()
        return cs.count(BSPQI_SITE) + cs.count(_rc(BSPQI_SITE))

    def validate(self) -> None:
        if len(self.left_linker) % 3 or len(self.right_linker) % 3:
            raise ValueError("linkers must be whole codons")
        if len(self.amplicon_linker_left) % 3 or len(self.amplicon_linker_right) % 3:
            raise ValueError("amplicon linkers must be whole codons")
        if len(self.left_linker) < 4 or len(self.right_linker) < 4:
            raise ValueError("each linker must be at least 4 nt to carry a PaqCI overhang")
        if self.arm_budget < 2 * self.min_arm:
            raise ValueError(f"arm budget {self.arm_budget} below 2 x min_arm")
        # A3 (adversarial panel, 2026-07-26) used to RAISE here unless `arm_budget % 3 == 0`, on the
        # grounds that the right arm would otherwise end mid-codon "in every design -- a silent
        # missense". **That rationale was wrong and the guard is REMOVED (Kevin, 2026-07-28):**
        #
        #   "I don't think the padding is necessary. We don't need the right arms to end exactly at
        #    the end of a codon. HDR machinery doesn't care. I only said it's helpful for us
        #    eyeballing each donor if they START at a codon when possible. No need to have them also
        #    end at one, and it's moot for UTR and introns anyways."
        #
        # He is right, and the error is worth naming because it cost a design decision. The right
        # arm is a COPY OF GENOMIC SEQUENCE. Ending mid-codon does not mutate anything -- the
        # genome simply supplies the rest of that codon after the crossover. There is no missense
        # to hide from C5/C7 because no codon is altered. The guard confused "the arm boundary
        # falls inside a codon" with "a codon is changed".
        #
        # Consequence: the guard rejected oligo_length=300 outright (arm_budget 179, remainder 2)
        # for BOTH cloning sites, which is why artifact `37` never evaluated the 300-nt single-BspQI
        # configuration and recommended 250 nt. See `62`.
        #
        # What DOES still hold: `split_arms()` keeps the LEFT arm a whole number of codons, so the
        # donor starts on a codon boundary -- which is the part Kevin actually asked for, and it is
        # a readability convenience, not a correctness requirement.
        if self.arm_budget % 3 and self.left_arm_must_be_codon_aligned:
            left, _ = self.split_arms()
            if left % 3:
                raise ValueError(
                    f"left arm {left} is not a whole number of codons; the donor would not start "
                    f"on a codon boundary")
