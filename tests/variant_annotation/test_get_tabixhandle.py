from pathlib import Path

import pytest
from genmod.annotate_variants import get_cadd_scores
from genmod.annotate_variants.read_tabix_files import TabixError, get_tabixhandle

CADD_PATH = "tests/fixtures/annotate_variant/small_CADD.tsv.gz"
UNZIPPED_VCF_PATH = "tests/fixtures/annotate_variant/small_1000G_maxAF.vcf"


def test_get_tabixhandle(thousand_g_path):
    # GIVEN a bgzipped file that exists and has a .tbi index
    assert Path(thousand_g_path).exists()
    assert Path(thousand_g_path + ".tbi").exists()

    # WHEN opening it
    handle = get_tabixhandle(thousand_g_path)

    # THEN a tabix handle is returned
    assert handle is not None


def test_get_tabixhandle_not_gzipped():
    # GIVEN an existing file that is not gzipped
    assert Path(UNZIPPED_VCF_PATH).exists()

    # WHEN opening it
    # THEN a TabixError is raised
    with pytest.raises(TabixError):
        get_tabixhandle(UNZIPPED_VCF_PATH)


def test_get_tabixhandle_missing_index(tmp_path):
    # GIVEN an existing .vcf.gz file without a .tbi index
    path = tmp_path / "no_index.vcf.gz"
    path.write_bytes(b"")
    assert path.exists()
    assert not Path(str(path) + ".tbi").exists()

    # WHEN opening it
    # THEN a TabixError is raised
    with pytest.raises(TabixError):
        get_tabixhandle(str(path))


def test_get_cadd_scores():
    # GIVEN an indexed CADD file with a C and no T alternative at 1:10001
    assert Path(CADD_PATH).exists()
    handle = get_tabixhandle(CADD_PATH)

    # WHEN fetching the scores for the C and the T alternative
    scores_c = get_cadd_scores(handle, "1", "10001", "C")
    scores_t = get_cadd_scores(handle, "1", "10001", "T")

    # THEN the C scores are returned and the T alternative has none
    assert scores_c == {"cadd_raw": "-0.642756", "cadd_phred": "1.165"}
    assert scores_t == {}
