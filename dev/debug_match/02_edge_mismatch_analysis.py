#!/usr/bin/env python3
"""
B. Where are step-2 perfect_middle_region-but-not-perfect_donor mismatches?

Run the matcher on yL406 sample 10. For rows where perfect_middle_region=True
but perfect_donor=False, retrieve the candidate designed-donor sequence and
align it positionally against the trimmed observed donor. Bucket each
mismatch position into (a) within first 5bp from 5' end, (b) within last 5bp
from 3' end, (c) middle. If most diffs concentrate at edges => boundary
problem. If scattered => real HDR end-resection / partial conversion noise.
"""
import sys, os
from pathlib import Path
from collections import Counter
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src" / "MAGESTIC_src"))
import pandas as pd, numpy as np
from magestic.pipelines.bc1_donor_bc0.core.matching import (
    build_lookup_tables, match_dataframe, trim_to_variable_region,
)

DATA = REPO / "data"
HARMONIZED = DATA / "20210422_and_20240411_Bloom_et_al_16_strains_QTL_harmonized_designed_variant_oligos.tsv"

print("--- load designs ---")
oligos = pd.read_csv(HARMONIZED, sep="\t",
                     usecols=["oligo_name","guide","donor","individual_variants","donor_changes"],
                     low_memory=False)
oligos = oligos.dropna(subset=["guide","donor"]).reset_index(drop=True)
oligo_name_to_donor = dict(zip(oligos["oligo_name"], oligos["donor"]))
print(f"  {len(oligos):,} oligos")

print("--- load + match sample 10 ---")
tsv = sorted(DATA.glob("*sample_10_bc1_donor_bc0.tsv"))[0]
obs = pd.read_csv(tsv, sep="\t", low_memory=False)
tables = build_lookup_tables(oligos, window_size=60)
matched = match_dataframe(obs, tables, donor_col="donor", guide_col="guide",
                          prefix="", n_jobs=1, chunk_size=5000, pre_trim=True)
print(f"  matched {len(matched):,} rows")

# Subset: perfect_middle_region=True AND perfect_donor=False AND has a single
# unambiguous candidate (so we know which design donor to compare against).
mask = (matched["perfect_middle_region"] == True) & (matched["perfect_donor"] == False)
sub = matched[mask].copy()
print(f"\nperfect_middle_region=True & perfect_donor=False: {len(sub):,} rows")
print(f"  (read-weighted: {sub['counts'].sum():,} reads)")

# A row's matched_oligo set is in 'matched_oligos' column (list[str]) or 'oligo_name'.
# Use the first oligo's design donor as the reference for position-diff.
def first_oligo(row):
    for c in ("oligo_name", "matched_oligos"):
        v = row.get(c)
        if isinstance(v, str) and v: return v.split(",")[0] if "," in v else v
        if isinstance(v, list) and v: return v[0]
    return None

print("\ncols on matcher output:")
print([c for c in matched.columns if "oligo" in c.lower() or "donor" in c.lower() or "perfect" in c.lower()])

# Compare trimmed observed donor vs design donor positionally
edge_5p = 0; edge_3p = 0; middle = 0; length_diff = 0; multi_diff = 0
diff_positions_5p_relative = []   # position from 5' end of design donor
diff_positions_3p_relative = []   # position from 3' end of design donor
length_deltas = []                # observed_len - design_len
examples = []

EDGE = 5

sample = sub.sample(min(5000, len(sub)), random_state=0)
for _, row in sample.iterrows():
    obs_raw = row["donor"]
    if not isinstance(obs_raw, str): continue
    obs_trimmed = trim_to_variable_region(obs_raw)
    if not obs_trimmed: continue
    name = first_oligo(row)
    if not name or name not in oligo_name_to_donor: continue
    design = oligo_name_to_donor[name]
    if len(obs_trimmed) != len(design):
        length_diff += 1
        length_deltas.append(len(obs_trimmed) - len(design))
        if len(examples) < 5:
            examples.append(("LENGTH_DIFF", obs_trimmed[:80], design[:80], len(obs_trimmed), len(design)))
        continue
    diff_idxs = [i for i in range(len(design)) if obs_trimmed[i] != design[i]]
    if not diff_idxs: continue
    if len(diff_idxs) > 1:
        multi_diff += 1
    for idx in diff_idxs:
        if idx < EDGE:
            edge_5p += 1
            diff_positions_5p_relative.append(idx)
        elif idx >= len(design) - EDGE:
            edge_3p += 1
            diff_positions_3p_relative.append(len(design) - 1 - idx)
        else:
            middle += 1
    if len(examples) < 5 and diff_idxs:
        examples.append((f"SUB({len(diff_idxs)})", obs_trimmed, design, len(obs_trimmed), diff_idxs[:5]))

total_mismatches = edge_5p + edge_3p + middle
print(f"\n--- positional mismatch analysis (sampled {len(sample):,} rows; perfect_middle+!perfect_donor) ---")
print(f"  length-mismatch rows:      {length_diff:>6,}  ({length_diff/len(sample):.1%})")
print(f"  multi-diff rows:           {multi_diff:>6,}  ({multi_diff/len(sample):.1%})")
print(f"  total mismatch positions:  {total_mismatches:>6,}")
if total_mismatches:
    print(f"    5' edge (first {EDGE}bp):  {edge_5p:>6,}  ({edge_5p/total_mismatches:.1%})")
    print(f"    3' edge (last  {EDGE}bp):  {edge_3p:>6,}  ({edge_3p/total_mismatches:.1%})")
    print(f"    middle:                    {middle:>6,}  ({middle/total_mismatches:.1%})")

if diff_positions_5p_relative:
    print(f"\n  5' diff position histogram (0=first bp): {Counter(diff_positions_5p_relative).most_common()}")
if diff_positions_3p_relative:
    print(f"  3' diff position histogram (0=last bp):  {Counter(diff_positions_3p_relative).most_common()}")
if length_deltas:
    c = Counter(length_deltas)
    print(f"  length-delta distribution (obs - design): {c.most_common(10)}")

print(f"\nexamples:")
for ex in examples:
    print(f"  {ex[0]}: obs len={ex[3]} design len={ex[-2] if isinstance(ex[-2],int) else ex[3]}")
    if isinstance(ex[1], str):
        print(f"    obs:    {ex[1][:120]}")
        print(f"    design: {ex[2][:120]}")
    print(f"    diff_idxs/lens: {ex[-1]}")
