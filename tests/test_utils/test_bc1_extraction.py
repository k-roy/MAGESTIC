"""Tests for the canonical bc1 extractor. Run: python3 test_bc1_extraction.py"""
import sys, random, itertools
pass
from magestic.utils.bc1_extraction import (extract_bc1, reverse_complement, _find_upstream, _seed_positions,
                            UPSTREAM_ANCHOR, UNIVERSAL_DOWNSTREAM, DOWNSTREAM_FALLBACK, BC1_LENGTH)

P=F=0
def chk(name, cond, extra=""):
    global P,F
    if cond: P+=1; print(f"  PASS  {name}")
    else:    F+=1; print(f"  FAIL  {name} {extra}")

random.seed(1337)
BC = "ACGTACGTACGTACGTACGT"
TAILS = {"A_1.0/2.5":("ATAACAGGGTAATGGAAA",18), "B_no_spacer":("",1), "C_3.0":("CCCTAGGATA",10)}
def build(layout, bc=BC, lead="GAGCTCCTACGACGTCGA"):
    spacer,off = TAILS[layout]
    # universal anchor sits at `off`; for B the tail starts with 'A' then the anchor
    pre = "A" if layout=="B_no_spacer" else spacer
    return lead + UPSTREAM_ANCHOR + bc + pre + UNIVERSAL_DOWNSTREAM + "CGTATAATG", off

print("\n== layout calling from the observed offset ==")
for lay,(_,off) in TAILS.items():
    read,exp = build(lay)
    r = extract_bc1(read)
    chk(f"{lay}: bc1 correct",   r.bc1==BC, f"got {r.bc1}")
    chk(f"{lay}: tier confirmed", r.anchor_tier=="confirmed", f"got {r.anchor_tier}")
    chk(f"{lay}: offset {exp}",   r.anchor_offset==exp, f"got {r.anchor_offset}")

print("\n== reverse strand ==")
read,_ = build("A_1.0/2.5")
r = extract_bc1(reverse_complement(read))
chk("rev: bc1 recovered", r.bc1==BC, f"got {r.bc1}")
chk("rev: orientation flagged", r.orientation=="rev", f"got {r.orientation}")
chk("rev not searched when disabled",
    extract_bc1(reverse_complement(read), search_both_strands=False).anchor_tier=="none")

print("\n== bounded, COUNTED error tolerance on the upstream anchor ==")
def mut(s,k,protect_last=1):
    idx=random.sample(range(len(s)-protect_last),k); s=list(s)
    for i in idx: s[i]=random.choice([c for c in "ACGT" if c!=s[i]])
    return "".join(s)
for k_data in (0,1,2):
    ok=0
    for _ in range(300):
        a=mut(UPSTREAM_ANCHOR,k_data)
        read="GAGCTCCTACGACGTCGA"+a+BC+TAILS["A_1.0/2.5"][0]+UNIVERSAL_DOWNSTREAM+"CG"
        if extract_bc1(read,k=1).bc1==BC: ok+=1
    if k_data<=1: chk(f"k=1 accepts {k_data} mismatch(es): {ok}/300", ok==300)
    else:         chk(f"k=1 REJECTS {k_data} mismatches: {ok}/300 accepted", ok==0)
chk("k=0 rejects 1 mismatch",
    extract_bc1("GAGCTCCTACGACGTCGA"+mut(UPSTREAM_ANCHOR,1)+BC+TAILS["A_1.0/2.5"][0]+UNIVERSAL_DOWNSTREAM,k=0).bc1 is None)

print("\n== protected 3' base (a mismatch there would shift the barcode start) ==")
bad = UPSTREAM_ANCHOR[:-1] + ("A" if UPSTREAM_ANCHOR[-1]!="A" else "T")
r = extract_bc1("GAGCTCCTACGACGTCGA"+bad+BC+TAILS["A_1.0/2.5"][0]+UNIVERSAL_DOWNSTREAM, k=1, protect=1)
chk("mismatch at the barcode-adjacent base is refused", r.bc1 is None, f"got tier={r.anchor_tier}")

print("\n== pigeonhole seeding is COMPLETE (never misses a <=k match) ==")
miss=0
for _ in range(2000):
    a=mut(UPSTREAM_ANCHOR,1)
    seq="".join(random.choice("ACGT") for _ in range(40))+a+"".join(random.choice("ACGT") for _ in range(40))
    if _find_upstream(seq,UPSTREAM_ANCHOR,1,0) is None: miss+=1
chk(f"2000 random k=1 placements, misses={miss}", miss==0)

print("\n== nothing is returned unlabelled ==")
noanchor="".join(random.choice("ACGT") for _ in range(150))
chk("no anchor -> tier 'none', bc1 None", extract_bc1(noanchor).anchor_tier=="none")
prefix_only="GAGCTC"+UPSTREAM_ANCHOR+BC+"".join(random.choice("TG") for _ in range(40))
r=extract_bc1(prefix_only)
chk("prefix but no downstream -> tier 'prefix_only', bc1 KEPT", r.anchor_tier=="prefix_only" and r.bc1==BC,
    f"got {r.anchor_tier}/{r.bc1}")
short="GAGCTC"+UPSTREAM_ANCHOR+BC[:12]
chk("truncated read -> length_mismatch, bc1 None", extract_bc1(short).anchor_tier=="length_mismatch")

print("\n== downstream confirmation is OFFSET-CONSTRAINED, never a scan ==")
LEAD="GAGCTCCTACGACGTCGA"
def mk(tail): return LEAD+UPSTREAM_ANCHOR+BC+tail
# truncated mid-anchor: longest fitting prefix must match exactly, >= 6 nt
chk("truncated read -> confirmed_truncated at the established offset",
    extract_bc1(mk("A"+UNIVERSAL_DOWNSTREAM[:6]), established_offset=1).anchor_tier=="confirmed_truncated")
chk("truncated below the 6-nt floor is REFUSED",
    extract_bc1(mk("A"+UNIVERSAL_DOWNSTREAM[:5]), established_offset=1).anchor_tier=="prefix_only")
chk("truncated 6-mer prefix with a mismatch is REFUSED (below the mm floor)",
    extract_bc1(mk("A"+"T"+UNIVERSAL_DOWNSTREAM[1:6]), established_offset=1).anchor_tier=="prefix_only")
_t8 = "T"+UNIVERSAL_DOWNSTREAM[1:8] if UNIVERSAL_DOWNSTREAM[0]!="T" else "A"+UNIVERSAL_DOWNSTREAM[1:8]
chk("truncated 8-mer prefix with ONE mismatch -> confirmed_truncated_mismatch",
    extract_bc1(mk("A"+_t8), established_offset=1).anchor_tier=="confirmed_truncated_mismatch")
_t8b = "TT"+UNIVERSAL_DOWNSTREAM[2:8]
chk("truncated 8-mer prefix with TWO mismatches is REFUSED",
    extract_bc1(mk("A"+_t8b), established_offset=1).anchor_tier=="prefix_only")
# single substitution inside a full-length anchor
mm = UNIVERSAL_DOWNSTREAM[:4]+("A" if UNIVERSAL_DOWNSTREAM[4]!="A" else "T")+UNIVERSAL_DOWNSTREAM[5:]
chk("1 mismatch at the established offset -> confirmed_mismatch",
    extract_bc1(mk("A"+mm+"TTTT"), established_offset=1).anchor_tier=="confirmed_mismatch")
mm2 = "A"+mm[1:]
mm2 = ("A" if UNIVERSAL_DOWNSTREAM[0]!="A" else "T")+mm[1:]
chk("2 mismatches are REFUSED",
    extract_bc1(mk("A"+mm2+"TTTT"), established_offset=1).anchor_tier=="prefix_only")
chk("no established offset -> no loosening at all",
    extract_bc1(mk("A"+mm+"TTTT")).anchor_tier=="prefix_only")
chk("loosening does NOT fire at a different offset",
    extract_bc1(mk("A"+mm+"TTTT"), established_offset=10).anchor_tier=="prefix_only")
# widened window still finds a genuinely shifted anchor, exactly
chk("shifted anchor within the widened window is found EXACTLY",
    extract_bc1(mk("T"*24+UNIVERSAL_DOWNSTREAM+"TT")).anchor_tier=="confirmed")
chk("anchor beyond the window is NOT found",
    extract_bc1(mk("T"*60+UNIVERSAL_DOWNSTREAM)).anchor_tier=="prefix_only")

print("\n== bc1 is never padded or normalized ==")
chk("no N-padding anywhere", all("N" not in (extract_bc1(x).bc1 or "") for x in (short, prefix_only, noanchor)))

print(f"\n{P} passed, {F} failed")
sys.exit(1 if F else 0)
