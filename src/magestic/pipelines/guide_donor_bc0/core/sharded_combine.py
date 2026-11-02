"""
Step 02, sharded: parse each collapsed FASTQ as its own cluster task, then merge per library.

Why this exists. ``process_library`` parses all of a library's FASTQs in one process and holds
the whole library in one dict; for the 2021 libraries that is up to 33 GB of input and 36 M
unique sequences per library, which is a multi-hour monolithic job that preemptible queues
requeue from zero. This module splits step 02 into idempotent pieces, each with its own rc
sentinel and provenance sidecar:

  parse  one collapsed FASTQ -> one shard TSV (one row per unique read, oriented + parsed)
  merge  all shards of one library -> ``{library}_guide_donor_bc0_counts.tsv`` in exactly the
         schema ``fastq_processing.write_counts_file`` produces (so step 03 is unchanged).

The parse step is stateless per read (a collapsed FASTQ holds each read once), so its memory
is O(1). The merge keys on the oriented read (``seq``) -- every other parsed field is a pure
function of ``seq`` -- and does two passes: pass 1 sums counts per ``seq`` into one array per
key; pass 2 re-streams the shards and writes each key the first time it is seen, popping it,
so peak memory is one dict of ``seq -> counts`` (~15 GB for the 36 M-key V454), not the
per-key 8-tuples ``process_library`` keeps.

Both steps write ``<output>.provenance.json`` (magestic.utils.provenance).

CLI (no entry point needed; run as a module)::

    python -m magestic.pipelines.guide_donor_bc0.core.sharded_combine manifest \
        --keyfile K --collapsed-fastq-dir D --out MANIFEST.tsv [--libraries V429,V454]
    python -m ... parse --manifest MANIFEST.tsv --shard-id N --shard-dir S [--md5 M]
    python -m ... merge --manifest MANIFEST.tsv --library V454 --shard-dir S --out-dir O

Author: Kevin R. Roy
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from array import array
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import pysam

from magestic.utils.collapsed_fastq import parse_collapsed_header_count
from magestic.utils.provenance import StepProvenance, md5_file, write_json_atomic
from .fastq_processing import (
    SEQ_CONST,
    load_sample_keyfile,
    new_parse_stats,
    parse_parameters,
    parse_sequence,
    rev_comp,
    _REVERSE_PRIMER_PAIRS,
)

SHARD_FIELDS = ("seq", "count", "leader", "guide", "donor", "SPS", "bc0", "msd", "donor_bc0_fragment")
MANIFEST_FIELDS = ("shard_id", "library_ID", "sample_number", "inner_primers", "outer_primers",
                   "sequencing_date", "direction", "counts_col", "fastq", "bytes")
COUNTS_SEQ_COLS = ("leader", "guide", "donor", "SPS", "bc0", "msd", "donor_bc0_fragment", "seq")


# ============================================================================
# Manifest
# ============================================================================

def build_manifest(keyfile: Path, collapsed_fastq_dir: Path, libraries: Optional[Sequence[str]] = None,
                   skip_libraries: Sequence[str] = ("empty",)) -> List[Dict[str, object]]:
    """One row per (sample, sequencing date) = one collapsed FASTQ = one parse shard.

    Raises if any referenced FASTQ is missing: a silently skipped input is exactly the kind
    of accounting hole the provenance policy forbids.
    """
    library_to_samples, _all_dates = load_sample_keyfile(keyfile)
    rows: List[Dict[str, object]] = []
    missing: List[str] = []
    sid = 0
    for library_ID, samples in library_to_samples.items():
        if library_ID in skip_libraries:
            continue
        if libraries and library_ID not in libraries:
            continue
        for s in samples:
            for date in s.sequencing_dates:
                fq = collapsed_fastq_dir / f"{date}_{s.outer_primers}_{s.sample_number}_collapsed.fastq"
                if not fq.exists():
                    missing.append(str(fq))
                    continue
                direction = "rev" if s.inner_primers in _REVERSE_PRIMER_PAIRS else "fwd"
                rows.append({
                    "shard_id": sid, "library_ID": library_ID, "sample_number": s.sample_number,
                    "inner_primers": s.inner_primers, "outer_primers": s.outer_primers,
                    "sequencing_date": date, "direction": direction,
                    "counts_col": f"{direction}_{date[-4:]}", "fastq": str(fq), "bytes": fq.stat().st_size,
                })
                sid += 1
    if missing:
        raise FileNotFoundError(f"{len(missing)} keyfile FASTQ(s) missing, e.g. {missing[:3]}")
    return rows


def write_manifest(rows: Iterable[Dict[str, object]], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def read_manifest(path: Path) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def all_counts_cols(manifest: Sequence[Dict[str, str]]) -> List[str]:
    """Every ``{fwd,rev}_{MMDD}`` column across the whole manifest, so all libraries share one
    schema (``process_library`` does the same via ``all_dates``).

    Order is ``sorted()``, matching ``fastq_processing.write_counts_file``'s
    ``sorted(per_run_cols)`` exactly, so the sharded table is column-identical to the
    in-process one (fwd_0906, fwd_1012, ..., rev_0906, ...)."""
    return sorted({f"{d}_{r['sequencing_date'][-4:]}" for r in manifest for d in ("fwd", "rev")})


# ============================================================================
# Parse: one FASTQ -> one shard
# ============================================================================

def shard_paths(shard_dir: Path, shard_id: int | str) -> Dict[str, Path]:
    base = shard_dir / f"shard_{int(shard_id):05d}"
    return {"tsv": base.with_suffix(".tsv"), "stats": base.with_suffix(".stats.json"),
            "rc": base.with_suffix(".rc")}


def orient_read(seq: str, inner_primers: str) -> str:
    """Apply the step-02 read-orientation rule: reverse-primer reads are reverse-complemented,
    and KR1882-KR1953 reads (which end 3 nt short) get 'GAC' appended."""
    if inner_primers in _REVERSE_PRIMER_PAIRS:
        seq = rev_comp(seq)
        if inner_primers == "KR1882-KR1953":
            seq += "GAC"
    return seq


def input_fingerprint(path: str | os.PathLike, edge: int = 1 << 20) -> Dict[str, object]:
    """Identity of an input file for resume decisions: size, mtime_ns, inode, and the md5 of its first and
    last MiB. The content part matters on Oak, whose mtimes have one-second resolution, so a same-size
    rewrite inside one second leaves every stat field unchanged."""
    st = os.stat(path)
    h = hashlib.md5()
    with open(path, "rb") as fh:
        h.update(fh.read(edge))
        if st.st_size > edge:
            fh.seek(max(st.st_size - edge, edge))
            h.update(fh.read())
    return {"size": st.st_size, "mtime_ns": st.st_mtime_ns, "ino": st.st_ino, "edge_md5": h.hexdigest()}


def parse_shard(row: Dict[str, str], shard_dir: Path, fastq_md5: Optional[str] = None) -> Dict[str, object]:
    """Parse one collapsed FASTQ into ``shard_<id>.tsv`` + stats + provenance sidecar.

    Idempotent ON INPUT IDENTITY: a shard is reused only if its rc sentinel reads ``0``, its TSV exists,
    and the FASTQ it was parsed from still has the fingerprint recorded at parse time (and the md5, when
    one is supplied). A changed FASTQ at the same manifest path is re-parsed; the old rule (sentinel +
    output exist) kept a stale shard silently (GDB-RESUME-R1, planning/288).
    """
    paths = shard_paths(shard_dir, row["shard_id"])
    if paths["rc"].exists() and paths["rc"].read_text().strip() == "0" and paths["tsv"].exists() \
            and paths["stats"].exists():
        try:
            st = json.load(open(paths["stats"]))
        except (OSError, ValueError):
            st = {}
        same = st.get("fastq_fingerprint") == input_fingerprint(row["fastq"])
        if fastq_md5 is not None:
            same = same and st.get("fastq_md5") == fastq_md5
        if same:
            return {"skipped": True, "shard_id": row["shard_id"], "tsv": str(paths["tsv"])}
        paths["rc"].unlink()

    shard_dir.mkdir(parents=True, exist_ok=True)
    fq = Path(row["fastq"])
    prov = StepProvenance(pipeline="guide_donor_bc0", step="02_combine.parse_shard", output=str(paths["tsv"]),
                          description="one collapsed FASTQ parsed into oriented reads with guide/donor/SPS/bc0 "
                                      "fields; one row per unique read; counts kept per row")
    prov.add_input(fq, role="collapsed_fastq", md5=fastq_md5, library_ID=row["library_ID"],
                   sample_number=row["sample_number"], sequencing_date=row["sequencing_date"],
                   inner_primers=row["inner_primers"], counts_col=row["counts_col"])
    prov.set_parameters(**parse_parameters(), shard_id=int(row["shard_id"]), counts_col=row["counts_col"],
                        direction=row["direction"])

    stats = new_parse_stats()
    n_rows = 0
    reads = 0
    tmp = paths["tsv"].with_suffix(".tsv.tmp")
    with pysam.FastxFile(str(fq)) as fin, open(tmp, "w") as out:
        out.write("\t".join(SHARD_FIELDS) + "\n")
        for entry in fin:
            count, _ = parse_collapsed_header_count(entry.name)
            reads += count
            seq = orient_read(entry.sequence, row["inner_primers"])
            leader, guide, donor, SPS, bc0, msd, frag = parse_sequence(seq, stats=stats)
            out.write(f"{seq}\t{count}\t{leader}\t{guide}\t{donor}\t{SPS}\t{bc0}\t{msd}\t{frag}\n")
            n_rows += 1
    os.replace(tmp, paths["tsv"])

    counts = {"unique_reads_in": n_rows, "reads_in": reads, "rows_out": n_rows, "parse": stats}
    if stats["reads"] != n_rows or stats["reads"] != stats["no_msd"] + stats["no_re_site"] + stats["bad_frame"] + stats["parsed"]:
        raise RuntimeError(f"parse accounting broken for shard {row['shard_id']}: {stats} vs rows {n_rows}")
    prov.set_counts(**counts)
    prov.add_output(paths["tsv"], rows=n_rows)
    prov.write()
    write_json_atomic({"shard_id": int(row["shard_id"]), "library_ID": row["library_ID"], "fastq": str(fq),
                       "counts_col": row["counts_col"], **counts, "tsv_md5": prov.outputs[-1]["md5"],
                       "fastq_fingerprint": input_fingerprint(fq), "fastq_md5": fastq_md5}, paths["stats"])
    paths["rc"].write_text("0\n")
    return {"skipped": False, "shard_id": row["shard_id"], "tsv": str(paths["tsv"]), **counts}


# ============================================================================
# Merge: all shards of one library -> counts table
# ============================================================================

def _iter_shard_rows(tsv: Path):
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if tuple(header) != SHARD_FIELDS:
            raise ValueError(f"{tsv}: unexpected shard header {header}")
        for line in fh:
            yield line.rstrip("\n").split("\t", 2)   # seq, count, rest-of-fields


def merge_library(library_ID: str, manifest: Sequence[Dict[str, str]], shard_dir: Path, out_dir: Path,
                  counts_cols: Optional[List[str]] = None) -> Dict[str, object]:
    """Sum the library's shards per oriented read and write ``{library}_guide_donor_bc0_counts.tsv``.

    Output schema (identical to ``fastq_processing.write_counts_file``)::

        library_ID  total_counts  total_fwd_counts  total_rev_counts  <fwd_MMDD/rev_MMDD ...>
        leader  guide  donor  SPS  bc0  msd  donor_bc0_fragment  seq

    Row order = first appearance across shards in manifest order (collapsed FASTQs are
    abundance-sorted, so this is roughly abundance order). Two passes; see module docstring.
    """
    rows = [r for r in manifest if r["library_ID"] == library_ID]
    if not rows:
        raise ValueError(f"{library_ID}: no manifest rows")
    if counts_cols is None:
        counts_cols = all_counts_cols(manifest)
    col_index = {c: i for i, c in enumerate(counts_cols)}
    ncol = len(counts_cols)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{library_ID}_guide_donor_bc0_counts.tsv"
    rc_path = out_dir / f".{library_ID}_guide_donor_bc0_counts.rc"
    sidecar = Path(str(out_path) + ".provenance.json")
    if rc_path.exists() and rc_path.read_text().strip() == "0" and out_path.exists() and sidecar.exists():
        # Reuse only if every shard is the one this merge consumed (a re-parsed shard has a new md5).
        try:
            used = {int(i["shard_id"]): i.get("md5") for i in json.load(open(sidecar)).get("inputs", [])
                    if i.get("role") == "shard"}
            now = {int(r["shard_id"]): json.load(open(shard_paths(shard_dir, r["shard_id"])["stats"])).get("tsv_md5")
                   for r in rows}
        except (OSError, ValueError, KeyError):
            used, now = None, {}
        if used == now:
            return {"skipped": True, "library_ID": library_ID, "output": str(out_path)}

    prov = StepProvenance(pipeline="guide_donor_bc0", step="02_combine.merge_library", output=str(out_path),
                          description="per-library guide_donor_bc0 counts table: shard rows summed per oriented "
                                      "read into per-run columns; schema identical to fastq_processing.write_counts_file")
    prov.set_parameters(**parse_parameters(), counts_cols=counts_cols, row_order="first appearance across shards "
                        "in manifest order", key="seq (every other parsed field is a pure function of seq)")

    # Inputs: the shards (with their recorded md5, verified against the bytes on disk) and,
    # through them, the FASTQs.
    shard_stats: List[Dict[str, object]] = []
    for r in rows:
        p = shard_paths(shard_dir, r["shard_id"])
        if not (p["rc"].exists() and p["rc"].read_text().strip() == "0"):
            raise RuntimeError(f"{library_ID}: shard {r['shard_id']} has no rc=0 sentinel ({p['rc']})")
        st = json.load(open(p["stats"]))
        actual = md5_file(p["tsv"])
        if actual != st["tsv_md5"]:
            raise RuntimeError(f"{library_ID}: shard {r['shard_id']} md5 {actual} != recorded {st['tsv_md5']}")
        prov.add_input(p["tsv"], role="shard", md5=actual, shard_id=int(r["shard_id"]), fastq=r["fastq"],
                       counts_col=r["counts_col"])
        shard_stats.append(st)

    # Pass 1: counts per seq.
    counts: Dict[str, array] = {}
    reads_in = 0
    rows_in = 0
    for r in rows:
        ci = col_index[r["counts_col"]]
        for seq, count, _rest in _iter_shard_rows(shard_paths(shard_dir, r["shard_id"])["tsv"]):
            c = int(count)
            reads_in += c
            rows_in += 1
            a = counts.get(seq)
            if a is None:
                a = array("q", [0]) * ncol
                counts[seq] = a
            a[ci] += c
    unique = len(counts)

    # Pass 2: write each key on first sight, popping it.
    fwd_idx = [i for i, c in enumerate(counts_cols) if c.startswith("fwd_")]
    rev_idx = [i for i, c in enumerate(counts_cols) if c.startswith("rev_")]
    header = ["library_ID", "total_counts", "total_fwd_counts", "total_rev_counts"] + counts_cols + list(COUNTS_SEQ_COLS)
    tmp = out_path.with_suffix(".tsv.tmp")
    rows_out = 0
    reads_out = 0
    with open(tmp, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            for seq, _count, rest in _iter_shard_rows(shard_paths(shard_dir, r["shard_id"])["tsv"]):
                a = counts.pop(seq, None)
                if a is None:
                    continue
                tf = sum(a[i] for i in fwd_idx)
                tr = sum(a[i] for i in rev_idx)
                reads_out += tf + tr
                rows_out += 1
                out.write(f"{library_ID}\t{tf + tr}\t{tf}\t{tr}\t" + "\t".join(str(x) for x in a)
                          + f"\t{rest}\t{seq}\n")
    if counts:
        raise RuntimeError(f"{library_ID}: {len(counts)} keys never written in pass 2")
    if rows_out != unique or reads_out != reads_in:
        raise RuntimeError(f"{library_ID}: accounting broken rows {rows_out}/{unique} reads {reads_out}/{reads_in}")
    os.replace(tmp, out_path)

    parse_tot = new_parse_stats()
    for st in shard_stats:
        for k, v in st["parse"].items():
            parse_tot[k] = parse_tot.get(k, 0) + v
    c = {"shards": len(rows), "shard_rows_in": rows_in, "reads_in": reads_in, "unique_seq": unique,
         "rows_out": rows_out, "reads_out": reads_out, "parse_reads_unique": parse_tot}
    prov.set_counts(**c)
    prov.add_output(out_path, rows=rows_out)
    prov.write()
    rc_path.write_text("0\n")
    return {"skipped": False, "library_ID": library_ID, "output": str(out_path), **c}


# ============================================================================
# CLI
# ============================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="sharded_combine", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("manifest", help="build the shard manifest from a keyfile")
    m.add_argument("--keyfile", required=True, type=Path)
    m.add_argument("--collapsed-fastq-dir", required=True, type=Path)
    m.add_argument("--out", required=True, type=Path)
    m.add_argument("--libraries", default=None, help="comma-separated subset of library_IDs")
    p = sub.add_parser("parse", help="parse one manifest shard")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--shard-id", required=True, type=int)
    p.add_argument("--shard-dir", required=True, type=Path)
    p.add_argument("--md5", default=None, help="raw md5 of the FASTQ if already known (else computed)")
    g = sub.add_parser("merge", help="merge one library's shards")
    g.add_argument("--manifest", required=True, type=Path)
    g.add_argument("--library", required=True)
    g.add_argument("--shard-dir", required=True, type=Path)
    g.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args(argv)

    if a.cmd == "manifest":
        libs = [x for x in a.libraries.split(",") if x] if a.libraries else None
        rows = build_manifest(a.keyfile, a.collapsed_fastq_dir, libs)
        write_manifest(rows, a.out)
        print(f"manifest: {len(rows)} shards, {sum(int(r['bytes']) for r in rows):,} bytes -> {a.out}")
        return 0
    manifest = read_manifest(a.manifest)
    if a.cmd == "parse":
        row = next((r for r in manifest if int(r["shard_id"]) == a.shard_id), None)
        if row is None:
            print(f"shard {a.shard_id} not in manifest", file=sys.stderr)
            return 2
        res = parse_shard(row, a.shard_dir, fastq_md5=a.md5)
        print(json.dumps(res, default=str))
        return 0
    if a.cmd == "merge":
        res = merge_library(a.library, manifest, a.shard_dir, a.out_dir)
        print(json.dumps(res, default=str))
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
