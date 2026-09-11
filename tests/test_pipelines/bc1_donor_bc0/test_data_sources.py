"""Layout-agnostic discovery of a screen's parsed bc1-donor-bc0 files.

The detector used to recognise only `2x150/` and `2x300/`, which is ONE screen's layout out of five. The
response on Oak was a tree of 772 symlinks faking that layout -- a standing policy violation that had also
silently drifted: 192 of 768 real parsed files had no link pointing at them. These tests pin the behaviour
that makes the shim unnecessary.
"""
import pytest

from magestic.pipelines.bc1_donor_bc0 import data_sources as dsrc

KEY_HEADER = "sequencing_date\tgDNA_plate_name\tinner_primers\tinner_product_size_bp\tcycles\n"


def _keyfile(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(KEY_HEADER + "".join(
        f"{d}\tplate_{i}\tKR1-KR2\t{bp}\t9+9\n" for i, (d, bp) in enumerate(rows)))
    return path


def _parsed(directory, names):
    directory.mkdir(parents=True, exist_ok=True)
    for n in names:
        (directory / n).write_text("bc1\tdonor\tbc0\tcounts\n")
    return directory


def _screen(tmp_path):
    """The real 20250912 shape: a flat parsed/ dir mixing a 2x150 run and three 2x300 runs."""
    p = tmp_path / "20250912_screen"
    _keyfile(p / "keyfiles" / "bc1_donor_bc0" / "plate_key.tsv",
             [("20251125", 152), ("20251212", 455), ("20260213", 457), ("20260302", 460)])
    _parsed(p / "processed_data" / "bc1_donor_bc0" / "parsed", [
        "20251125_plate_2A_sample_1_bc1_donor_bc0.tsv",
        "20251125_plate_2A_sample_2_bc1_donor_bc0.tsv",
        "20251212_plate_4C_sample_1_bc1_donor_bc0.tsv",
        "20260213_plate_4C_sample_2_bc1_donor_bc0.tsv",
        "20260302_plate_4C_sample_3_bc1_donor_bc0.tsv",
    ])
    return p


def test_flat_parsed_dir_is_found_and_split_by_the_keyfile_not_by_directory_names(tmp_path):
    ds = dsrc.discover(_screen(tmp_path))
    assert ds.total_files == 5 and not ds.unclassified
    assert len(ds.files_by_layout["2x150"]) == 2      # the 152 bp run
    assert len(ds.files_by_layout["2x300"]) == 3      # the 451-460 bp runs
    assert "parsed" in ds.how and "inner_product_size_bp" in ds.how


def test_legacy_2x150_2x300_layout_still_works(tmp_path):
    """yL406 is laid out this way and must not regress."""
    p = tmp_path / "yL406"
    base = p / "processed_data" / "bc1_donor_bc0"
    _parsed(base / "2x300", ["20240806_plate_1A_sample_1_bc1_donor_bc0.tsv"])
    ds = dsrc.discover(p)
    assert ds.files_by_layout["2x300"] and "legacy" in ds.how
    assert "2x150" not in ds.files_by_layout


def test_explicit_directories_win_over_detection(tmp_path):
    p = _screen(tmp_path)
    other = _parsed(tmp_path / "elsewhere", ["20991231_plate_9Z_sample_9_bc1_donor_bc0.tsv"])
    ds = dsrc.discover(p, data_2x300_dir=other)
    assert ds.total_files == 1 and "explicit" in ds.how
    ds2 = dsrc.discover(p, parsed_dir=other)
    assert ds2.total_files == 1 and "explicit --parsed-dir" in ds2.how


def test_an_unmapped_sequencing_date_is_reported_not_silently_dropped(tmp_path):
    """The failure that made the shim dangerous: files that no keyfile explains must be VISIBLE."""
    p = _screen(tmp_path)
    (p / "processed_data" / "bc1_donor_bc0" / "parsed"
     / "20770101_plate_9Z_sample_1_bc1_donor_bc0.tsv").write_text("x\n")
    ds = dsrc.discover(p)
    assert len(ds.unclassified) == 1
    assert ds.total_files == 6
    assert "unclassified: 1" in ds.describe()


def test_nothing_found_reports_everywhere_it_looked(tmp_path):
    p = tmp_path / "empty_screen"
    (p / "processed_data" / "bc1_donor_bc0").mkdir(parents=True)
    ds = dsrc.discover(p)
    assert ds.total_files == 0 and ds.how == "no data source found"
    assert any(s.endswith("2x150") for s in ds.searched)
    assert any(s.endswith("2x300") for s in ds.searched)
    assert any(s.endswith("parsed") for s in ds.searched)
    assert "NOTHING FOUND" in ds.describe()


def test_read_layout_map_handles_multi_date_cells_and_conflicts(tmp_path):
    k = _keyfile(tmp_path / "k.tsv", [("20251212,20260213", 457), ("20251125", 152)])
    m = dsrc.read_layout_map(k)
    assert m == {"20251212": "2x300", "20260213": "2x300", "20251125": "2x150"}
    # a date carrying both a short and a long amplicon takes the LONGER requirement
    k2 = _keyfile(tmp_path / "k2.tsv", [("20260106", 135), ("20260106", 455)])
    assert dsrc.read_layout_map(k2) == {"20260106": "2x300"}


def test_keyfile_without_the_needed_columns_is_not_fatal(tmp_path):
    k = tmp_path / "bad.tsv"
    k.write_text("a\tb\n1\t2\n")
    assert dsrc.read_layout_map(k) == {}
    k2 = tmp_path / "empty.tsv"
    k2.write_text("")
    assert dsrc.read_layout_map(k2) == {}


def test_discovery_follows_symlinks_but_never_needs_them(tmp_path):
    """The shim can be deleted without breaking anything: real files are found directly."""
    p = _screen(tmp_path)
    shim = p / "processed_data" / "bc1_donor_bc0_PACKAGE_RERUN" / "yL442" / "processed_data" / "bc1_donor_bc0" / "2x300"
    shim.mkdir(parents=True)
    real = p / "processed_data" / "bc1_donor_bc0" / "parsed" / "20251212_plate_4C_sample_1_bc1_donor_bc0.tsv"
    (shim / real.name).symlink_to(real)
    ds = dsrc.discover(p)
    assert ds.total_files == 5
    assert not any("PACKAGE_RERUN" in str(f) for v in ds.files_by_layout.values() for f in v)
