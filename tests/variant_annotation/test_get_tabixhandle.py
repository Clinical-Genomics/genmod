import pytest
from genmod.annotate_variants.read_tabix_files import TabixError, get_tabixhandle


def test_get_tabixhandle(thousand_g_path):
    assert get_tabixhandle(thousand_g_path) is not None


def test_get_tabixhandle_not_gzipped():
    with pytest.raises(TabixError):
        get_tabixhandle("tests/fixtures/annotate_variant/small_1000G_maxAF.vcf")


def test_get_tabixhandle_missing_index(tmp_path):
    path = tmp_path / "no_index.vcf.gz"
    path.write_bytes(b"")
    with pytest.raises(TabixError):
        get_tabixhandle(str(path))


def test_get_cadd_scores():
    from genmod.annotate_variants import get_cadd_scores

    handle = get_tabixhandle("tests/fixtures/annotate_variant/small_CADD.tsv.gz")
    assert get_cadd_scores(handle, "1", "10001", "C") == {
        "cadd_raw": "-0.642756",
        "cadd_phred": "1.165",
    }
    assert get_cadd_scores(handle, "1", "10001", "T") == {}
