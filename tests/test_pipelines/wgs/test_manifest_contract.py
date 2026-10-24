"""One manifest parser for both execution paths (audit WF-WGS-HASH-PATH).

At 9d654b8 the Snakefile parsed the manifest itself with `pd.read_csv(..., comment="#")`,
which truncates any field at its first `#`, while the direct path kept the field whole. The
workflow now imports `core.manifest.load_manifest`, so the two paths cannot disagree.
"""
from pathlib import Path

import magestic.pipelines.wgs.snakemake as workflow
from magestic.pipelines.wgs.core.manifest import load_manifest


def test_hash_inside_a_field_is_data_and_only_a_leading_hash_comments_a_row(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample\tbam\treference\n"
                        "# a comment row\tx.bam\tx.fa\n"
                        "s1\tbams/plate#1/s1.bam\tref.fa\n")
    samples = load_manifest(manifest)
    assert [s.sample for s in samples] == ["s1"]
    assert samples[0].bam == tmp_path / "bams/plate#1/s1.bam"


def test_workflow_reads_the_manifest_with_the_shared_loader():
    text = (Path(workflow.__file__).parent / "Snakefile").read_text()
    assert "from magestic.pipelines.wgs.core.manifest import load_manifest" in text
    assert "read_csv" not in text and "comment=" not in text
