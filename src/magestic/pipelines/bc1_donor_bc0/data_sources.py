"""
Locate the parsed bc1-donor-bc0 files of a screen, whatever layout it uses on disk.

WHY THIS MODULE EXISTS. ``PipelineConfig.from_project_dir`` detected exactly two directory names --
``processed_data/bc1_donor_bc0/2x150`` and ``.../2x300`` -- plus two hard-coded dated paths. Only ONE screen on
Oak is laid out that way (20240806_yL406). Four others keep their parsed TSVs in
``processed_data/bc1_donor_bc0/parsed/``, so the detector found nothing for them, and the response was to build a
tree of SYMLINKS faking the expected layout:

    processed_data/bc1_donor_bc0_PACKAGE_RERUN/<yL>/processed_data/bc1_donor_bc0/2x300/*.tsv -> ../../parsed/*.tsv

That shim is a standing violation of the no-links-on-Oak policy, it was removed from yL406 and reappeared 2.7x
larger (772 links) in 20250912, and -- measured 2026-09-11 -- **it is also silently incomplete: 192 of the 768
real parsed files have no link pointing at them, all of them the 20251125 run.** A consumer reading through the
shim sees 576 files and no error. That is the whole argument for fixing the tool instead of shimming the tree:
the shim is a second, hand-maintained copy of the truth, and it had already drifted.

THE READ LAYOUT IS RECORDED IN TEXT, so it never needs to be encoded in directory names. Every screen's plate
keyfile carries ``sequencing_date`` and ``inner_product_size_bp``; the layout follows from whether the amplicon
fits the read. Measured on 20250912: the 20251125 run is ``inner_product_size_bp=152`` (2x150) and the
20251212/20260213/20260302 runs are 451-460 (2x300). Parsed filenames begin with their sequencing date, so the
keyfile classifies them without any filesystem aliasing.

Author: Kevin Roy Lab
"""

from __future__ import annotations

import csv
import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

logger = logging.getLogger(__name__)

PARSED_SUFFIX = "_bc1_donor_bc0.tsv"
#: A parsed filename starts with its sequencing date: ``20251125_plate_2A_sample_25_bc1_donor_bc0.tsv``.
DATE_PREFIX_RE = re.compile(r"^(\d{8})_")
#: 2x150 reads merge across an insert up to ~300 nt; the real screens sit far from this boundary
#: (152 vs 451-460), so the exact cut only has to be defensible, not tuned.
MAX_INSERT_FOR_2X150 = 300
LAYOUT_2X150 = "2x150"
LAYOUT_2X300 = "2x300"


@dataclass
class DataSources:
    """Where a screen's parsed bc1-donor-bc0 files actually are, and how they were found."""

    files_by_layout: Dict[str, List[Path]] = field(default_factory=dict)
    dirs_by_layout: Dict[str, Path] = field(default_factory=dict)
    unclassified: List[Path] = field(default_factory=list)
    searched: List[str] = field(default_factory=list)
    how: str = ""

    @property
    def total_files(self) -> int:
        return sum(len(v) for v in self.files_by_layout.values()) + len(self.unclassified)

    def describe(self) -> str:
        parts = [f"{k}: {len(v)} files" for k, v in sorted(self.files_by_layout.items())]
        if self.unclassified:
            parts.append(f"unclassified: {len(self.unclassified)}")
        return f"{self.how} -> " + (", ".join(parts) if parts else "NOTHING FOUND")


def read_layout_map(keyfile: str | os.PathLike) -> Dict[str, str]:
    """``sequencing_date -> '2x150'|'2x300'`` from a plate/combined keyfile.

    Uses ``inner_product_size_bp``: an insert a 2x150 pair cannot span needs 2x300. A ``sequencing_date``
    cell may hold several comma-separated dates (a library sequenced more than once); each is mapped.
    Conflicting sizes for one date resolve to the LARGER requirement, because a run that must cover the
    longest amplicon was necessarily the longer read.
    """
    out: Dict[str, str] = {}
    with open(keyfile, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return out
    cols = rows[0].keys()
    date_col = next((c for c in cols if c and c.startswith("sequencing_date")), None)
    size_col = next((c for c in cols if c == "inner_product_size_bp"), None)
    if not date_col or not size_col:
        return out
    for row in rows:
        raw_size = (row.get(size_col) or "").strip()
        if not raw_size:
            continue
        try:
            size = int(float(raw_size))
        except ValueError:
            continue
        layout = LAYOUT_2X150 if size <= MAX_INSERT_FOR_2X150 else LAYOUT_2X300
        for date in (d.strip() for d in (row.get(date_col) or "").split(",")):
            if not date:
                continue
            if out.get(date) == LAYOUT_2X300:
                continue          # the larger requirement already won
            out[date] = layout
    return out


def find_keyfiles(project_dir: Path) -> List[Path]:
    """Plate/combined keyfiles of a screen, most specific first."""
    kf = project_dir / "keyfiles"
    if not kf.is_dir():
        return []
    found = [p for p in sorted(kf.rglob("*.tsv"))
             if "plate_key" in p.name or "combined_key" in p.name]
    return found


def _parsed_files(directory: Path) -> List[Path]:
    return sorted(p for p in directory.glob(f"*{PARSED_SUFFIX}") if p.is_file())


def discover(project_dir: str | os.PathLike,
             parsed_dir: Optional[str | os.PathLike] = None,
             data_2x150_dir: Optional[str | os.PathLike] = None,
             data_2x300_dir: Optional[str | os.PathLike] = None,
             keyfiles: Optional[Sequence[str | os.PathLike]] = None) -> DataSources:
    """Find a screen's parsed bc1-donor-bc0 files. Explicit arguments always win.

    Order: explicit per-layout directories, then explicit ``parsed_dir``, then the legacy ``2x150``/``2x300``
    directories (back-compat for yL406), then ``parsed/`` classified by the keyfile. **Symlinks are followed
    but never required**, and nothing here creates one.

    A result whose ``total_files`` is 0 carries the full ``searched`` list, so "no data" is always
    reportable as "looked here, here and here" rather than as a bare warning.
    """
    project_dir = Path(project_dir)
    ds = DataSources()
    base = project_dir / "processed_data" / "bc1_donor_bc0"

    explicit = {LAYOUT_2X150: data_2x150_dir, LAYOUT_2X300: data_2x300_dir}
    if any(explicit.values()):
        for layout, d in explicit.items():
            if not d:
                continue
            d = Path(d)
            ds.searched.append(str(d))
            files = _parsed_files(d)
            if files:
                ds.files_by_layout[layout] = files
                ds.dirs_by_layout[layout] = d
        ds.how = "explicit --data-2x150-dir/--data-2x300-dir"
        if ds.total_files:
            return ds

    if parsed_dir:
        pd_ = Path(parsed_dir)
        ds.searched.append(str(pd_))
        files = _parsed_files(pd_)
        if files:
            ds.dirs_by_layout["parsed"] = pd_
            _classify(ds, files, project_dir, keyfiles)
            ds.how = f"explicit --parsed-dir {pd_}"
            return ds

    for layout in (LAYOUT_2X150, LAYOUT_2X300):
        d = base / layout
        ds.searched.append(str(d))
        if d.is_dir():
            files = _parsed_files(d)
            if files:
                ds.files_by_layout[layout] = files
                ds.dirs_by_layout[layout] = d
    if ds.total_files:
        ds.how = f"legacy {base}/{{2x150,2x300}} layout"
        return ds

    d = base / "parsed"
    ds.searched.append(str(d))
    if d.is_dir():
        files = _parsed_files(d)
        if files:
            ds.dirs_by_layout["parsed"] = d
            _classify(ds, files, project_dir, keyfiles)
            ds.how = f"{d} classified by keyfile inner_product_size_bp"
            return ds

    ds.how = "no data source found"
    return ds


def _classify(ds: DataSources, files: List[Path], project_dir: Path,
              keyfiles: Optional[Sequence[str | os.PathLike]]) -> None:
    """Split a flat parsed directory into read layouts using the keyfile record."""
    kfs = [Path(k) for k in keyfiles] if keyfiles else find_keyfiles(project_dir)
    layout_map: Dict[str, str] = {}
    for k in kfs:
        try:
            layout_map.update({d: l for d, l in read_layout_map(k).items() if d not in layout_map})
        except Exception as exc:                       # a malformed keyfile must not kill discovery
            logger.warning("could not read layout map from %s: %s", k, exc)
    for f in files:
        m = DATE_PREFIX_RE.match(f.name)
        layout = layout_map.get(m.group(1)) if m else None
        if layout:
            ds.files_by_layout.setdefault(layout, []).append(f)
        else:
            ds.unclassified.append(f)
    for v in ds.files_by_layout.values():
        v.sort()
    ds.unclassified.sort()
