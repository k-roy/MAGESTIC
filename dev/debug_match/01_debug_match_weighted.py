#!/usr/bin/env python3
"""Re-run the matcher with read weighting to compare against step-1 V545 numbers."""
import sys, argparse
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src" / "MAGESTIC_src"))
import pandas as pd
from magestic.pipelines.bc1_donor_bc0.core.matching import build_lookup_tables, match_dataframe

DATA = REPO / "data"
HARMONIZED = DATA / "20210422_and_20240411_Bloom_et_al_16_strains_QTL_harmonized_designed_variant_oligos.tsv"

p = argparse.ArgumentParser()
p.add_argument("--n-samples", type=int, default=1)
p.add_argument("--top-n", type=int, default=200_000)
p.add_argument("--no-cap", action="store_true")
args = p.parse_args()

print("--- load designs + tsv ---")
oligos = pd.read_csv(HARMONIZED, sep="\t", usecols=["oligo_name","guide","donor","individual_variants","donor_changes"], low_memory=False)
oligos = oligos.dropna(subset=["guide","donor"]).reset_index(drop=True)
print(f"oligos: {len(oligos):,}")

tsvs = sorted(DATA.glob("*_bc1_donor_bc0.tsv"))[:args.n_samples]
dfs = [pd.read_csv(p, sep="\t", low_memory=False) for p in tsvs]
obs = pd.concat(dfs, ignore_index=True)
if (not args.no_cap) and len(obs) > args.top_n:
    obs = obs.head(args.top_n)
print(f"obs rows: {len(obs):,}    total reads (counts): {obs['counts'].sum():,}")

tables = build_lookup_tables(oligos, window_size=60)
matched = match_dataframe(obs, tables, donor_col="donor", guide_col="guide", prefix="", n_jobs=1, chunk_size=5000, pre_trim=True)
print(f"\nmatcher output rows: {len(matched):,}")

# parse-failure denom
donor_str = matched["donor"].astype("string")
is_parse_fail = donor_str.isna() | (donor_str.str.len()==0)
n_total_rows = len(matched); n_total_reads = matched["counts"].sum()
m = matched[~is_parse_fail]
print(f"parse failures: {is_parse_fail.sum():,} rows, {matched.loc[is_parse_fail,'counts'].sum():,} reads")
n_rows = len(m); n_reads = m["counts"].sum()
print(f"matcheable denom: {n_rows:,} rows, {n_reads:,} reads\n")

cols = ["perfect_donor","perfect_middle_region","is_unambiguous"]
for col in cols:
    rt = int(m[col].sum()); rd = int(m.loc[m[col],"counts"].sum())
    print(f"  {col:30s}  rows {rt:>7,} ({rt/n_rows:6.1%})   reads {rd:>12,} ({rd/n_reads:6.1%})")

print("\nmatch_method (rows | reads):")
for k, sub in matched.groupby("match_method"):
    r = len(sub); rd = sub["counts"].sum()
    print(f"  {k:30s}  rows {r:>7,} ({r/n_total_rows:6.1%})   reads {rd:>12,} ({rd/n_total_reads:6.1%})")
