"""A cached reference backbone is only used while it still matches the run.

scripts/CheckReferenceBackbone.py gates params.ref_set_aligned. The shipped HCV
and RABV caches must pass for their profiles, and each way a cache goes stale must fail.
"""

import shutil
from pathlib import Path

import pytest

import CheckReferenceBackbone as CRB

REPO_ROOT = Path(__file__).resolve().parents[2]
HCV_CACHE = REPO_ROOT / "generic" / "hcv" / "ref_set_aligned"
HCV_REF_LIST = REPO_ROOT / "generic" / "hcv" / "ref_list_subtype_genotype.txt"

RABV_CACHE = REPO_ROOT / "generic" / "rabv" / "ref_set_aligned"
RABV_REF_LIST = REPO_ROOT / "generic" / "rabv" / "ref_list_clades.txt"

requires_cache = pytest.mark.skipif(not (HCV_CACHE / CRB.MANIFEST).exists(), reason="HCV backbone cache not present")
requires_rabv_cache = pytest.mark.skipif(not (RABV_CACHE / CRB.MANIFEST).exists(), reason="RABV backbone cache not present")


def _tiny_cache(tmp_path, accessions=("NC_1", "R1", "R2")):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "refset_1_aln.fasta").write_text("".join(f">{a}\nACGT\n" for a in accessions))
    ref_list = tmp_path / "ref_list.txt"
    ref_list.write_text("NC_1\tmaster\t1\n" + "".join(f"{a}\treference\t1\n" for a in accessions[1:]))
    CRB.write_manifest(str(cache), str(ref_list), 2, "unit test")
    return cache, ref_list


@requires_cache
def test_shipped_hcv_cache_matches_its_profiles():
    failures, warnings, notes = CRB.check(str(HCV_CACHE), str(HCV_REF_LIST), 2)
    assert failures == []
    assert warnings == [], "BuildReferenceAlignment.py changed since the HCV cache was made: rebuild it"
    assert "238 references" in notes[0]


@requires_rabv_cache
def test_shipped_rabv_cache_matches_its_profile():
    failures, warnings, notes = CRB.check(str(RABV_CACHE), str(RABV_REF_LIST), 2)
    assert failures == []
    assert warnings == [], "BuildReferenceAlignment.py changed since the RABV cache was made: rebuild it"
    assert "497 references" in notes[0]


@requires_cache
def test_cache_carries_everything_downstream_reads():
    names = {p.name for p in HCV_CACHE.iterdir()}
    assert {"refset_1_aln.fasta", "dropped_insertions.tsv", "build_report.tsv", CRB.MANIFEST} <= names


@requires_rabv_cache
def test_rabv_cache_carries_everything_downstream_reads():
    names = {p.name for p in RABV_CACHE.iterdir()}
    assert {"refset_1_aln.fasta", "dropped_insertions.tsv", "build_report.tsv", CRB.MANIFEST} <= names


def test_fresh_cache_passes(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path)
    assert CRB.check(str(cache), str(ref_list), 2)[0] == []


def test_reference_added_to_list_fails(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path)
    ref_list.write_text(ref_list.read_text() + "R3\treference\t1\n")
    failures = CRB.check(str(cache), str(ref_list), 2)[0]
    assert failures and "R3" in failures[0]


def test_duplicate_rows_in_the_list_are_not_a_change(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path)
    ref_list.write_text(ref_list.read_text() + "R1\treference\t1\n")
    assert CRB.check(str(cache), str(ref_list), 2)[0] == []


def test_exclusion_list_rows_are_not_expected_in_the_backbone(tmp_path):
    """The flu list names influenza B/C/D as exclusion_list; the builder leaves them out."""
    cache, ref_list = _tiny_cache(tmp_path)
    ref_list.write_text(ref_list.read_text() + "B_FLU\texclusion_list\t1\n")
    assert CRB.check(str(cache), str(ref_list), 2)[0] == []
    CRB.write_manifest(str(cache), str(ref_list), 2, "with exclusions")


def test_different_insertion_support_fails(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path)
    failures = CRB.check(str(cache), str(ref_list), 1)[0]
    assert any("min_insertion_support" in f for f in failures)


def test_changed_builder_warns(tmp_path, monkeypatch):
    cache, ref_list = _tiny_cache(tmp_path)
    fake = tmp_path / "BuildReferenceAlignment.py"
    fake.write_text("# edited\n")
    monkeypatch.setattr(CRB, "BUILDER", str(fake))
    failures, warnings, _ = CRB.check(str(cache), str(ref_list), 2)
    assert failures == [] and warnings


def test_curated_backbone_holding_every_reference_passes(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path)
    (cache / CRB.MANIFEST).unlink()
    failures, _, notes = CRB.check(str(cache), str(ref_list), 2)
    assert failures == [] and "coverage checked" in notes[0]


def test_curated_backbone_missing_a_reference_fails(tmp_path):
    """No manifest is no exemption: a run must never align against a backbone lacking a reference."""
    cache, ref_list = _tiny_cache(tmp_path)
    (cache / CRB.MANIFEST).unlink()
    ref_list.write_text(ref_list.read_text() + "R3\treference\t1\n")
    failures = CRB.check(str(cache), str(ref_list), 2)[0]
    assert failures and "R3" in failures[0]


def test_curated_backbone_extra_rows_only_warn(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path, accessions=("NC_1", "R1", "R2", "UNLISTED"))
    ref_list.write_text("NC_1\tmaster\t1\nR1\treference\t1\nR2\treference\t1\n")
    (cache / CRB.MANIFEST).unlink()
    failures, warnings, _ = CRB.check(str(cache), str(ref_list), 2)
    assert failures == [] and "UNLISTED" in warnings[0]


def test_cache_with_rows_no_longer_listed_fails(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path)
    ref_list.write_text("NC_1\tmaster\t1\nR1\treference\t1\n")
    failures = CRB.check(str(cache), str(ref_list), 2)[0]
    assert failures and "R2" in failures[0]


def _segmented(tmp_path, seg2_rows):
    cache = tmp_path / "seg"
    cache.mkdir()
    (cache / "refset_1_aln.fasta").write_text(">M1\nACGT\n>R1\nACGT\n")
    (cache / "refset_2_aln.fasta").write_text("".join(f">{a}\nACGT\n" for a in seg2_rows))
    ref_list = tmp_path / "seg_list.txt"
    ref_list.write_text("M1\tmaster\t1\nM2\tmaster\t2\nR1\treference\t1\nR2\treference\t2\n")
    return cache, ref_list


def test_segmented_backbone_with_each_reference_in_its_segment_passes(tmp_path):
    cache, ref_list = _segmented(tmp_path, ["M2", "R2"])
    assert CRB.check(str(cache), str(ref_list), 2, "Y")[0] == []


def test_segmented_reference_in_the_wrong_segment_file_fails(tmp_path):
    cache, ref_list = _segmented(tmp_path, ["M2"])
    (cache / "refset_1_aln.fasta").write_text(">M1\nACGT\n>R1\nACGT\n>R2\nACGT\n")
    failures = CRB.check(str(cache), str(ref_list), 2, "Y")[0]
    assert failures and "R2 (segment 2" in failures[0]


def test_segment_with_no_backbone_file_fails(tmp_path):
    cache, ref_list = _segmented(tmp_path, ["M2", "R2"])
    (cache / "refset_2_aln.fasta").unlink()
    failures = CRB.check(str(cache), str(ref_list), 2, "Y")[0]
    assert failures and "no backbone file" in failures[0]


INFLUENZA_BACKBONE = REPO_ROOT / "generic" / "influenza" / "ref_set_aligned"
INFLUENZA_REF_LIST = REPO_ROOT / "generic" / "influenza" / "ref_list_refmast.txt"


@pytest.mark.skipif(not INFLUENZA_BACKBONE.exists(), reason="influenza backbone not present")
def test_shipped_influenza_backbone_holds_every_reference_in_its_segment():
    failures, warnings, _ = CRB.check(str(INFLUENZA_BACKBONE), str(INFLUENZA_REF_LIST), 2, "Y")
    assert failures == []


def test_manifest_refuses_a_mismatched_cache(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path)
    ref_list.write_text("NC_1\tmaster\t1\n")
    with pytest.raises(ValueError):
        CRB.write_manifest(str(cache), str(ref_list), 2, "x")


def test_cli_exit_codes(tmp_path):
    cache, ref_list = _tiny_cache(tmp_path)
    args = ["--backbone_dir", str(cache), "--ref_list", str(ref_list)]
    assert CRB.main(args + ["--min_insertion_support", "2"]) == 0
    assert CRB.main(args + ["--min_insertion_support", "3"]) == 1


def test_workflow_checks_a_supplied_backbone_before_using_it():
    text = (REPO_ROOT / "vgtk-init.nf").read_text()
    assert "process CHECK_REFERENCE_BACKBONE" in text
    assert text.count("ref_backbone_ch = CHECK_REFERENCE_BACKBONE.out.backbone") == 2, \
        "both a built and a supplied backbone must pass the check"
    assert "CHECK_REFERENCE_BACKBONE(BUILD_REFERENCE_ALIGNMENT.out.ref_set_aligned" in text
    assert "--is_segmented !{params.is_segmented}" in text
    assert "in ['', 'null', 'false']" in text
    config = (REPO_ROOT / "nextflow.config").read_text()
    assert config.count('ref_set_aligned   = "${projectDir}/generic/hcv/ref_set_aligned"') == 3
    # setup_rabv_full only: the `test` profile keeps building so CI exercises the builder.
    assert config.count('ref_set_aligned   = "${projectDir}/generic/rabv/ref_set_aligned"') == 1
    bash = (REPO_ROOT / "vgtk-rabv.sh").read_text()
    assert bash.count("CheckReferenceBackbone.py") == 2
    assert 'ValidateRefListAgainstDb.py" --ref_list "$REF_LIST" --db "$UPDATE_DB"' in bash
