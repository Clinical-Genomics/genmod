from pathlib import Path
from typing import Any

import pysam
import pytest
from genmod.annotate_variants import get_cosmic, get_spidex_score
from genmod.annotate_variants.read_tabix_files import get_tabixhandle


def make_tabix_file(path: Path, lines: list[str], **index_args: Any) -> pysam.TabixFile:
    """Write rows to a bgzipped, tabix indexed file and return a handle to it"""
    path.write_text("".join(f"{line}\n" for line in lines))
    pysam.tabix_index(str(path), force=True, **index_args)
    return get_tabixhandle(f"{path}.gz")


@pytest.fixture
def spidex_handle(tmp_path: Path) -> pysam.TabixFile:
    """Return a handle to a SPIDEX tsv: chrom, pos, ref, alt, dpsi_max_tissue, dpsi_zscore"""
    lines = [
        "#Chromosome\tPosition\tRef\tAlt\tdpsi_max_tissue\tdpsi_zscore",
        "1\t100\tA\tC\t-1.5\t-0.75",
        "1\t100\tA\tG\t2.5\t1.25",
    ]
    return make_tabix_file(tmp_path / "spidex.tsv", lines, seq_col=0, start_col=1, end_col=1)


@pytest.fixture
def cosmic_handle(tmp_path: Path) -> pysam.TabixFile:
    """Return a handle to a COSMIC vcf with one record at 1:100 (A>C)"""
    lines = [
        "##fileformat=VCFv4.1",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "1\t100\tCOSM1\tA\tC\t.\t.\t.",
    ]
    return make_tabix_file(tmp_path / "cosmic.vcf", lines, preset="vcf")


def test_get_spidex_score(spidex_handle: pysam.TabixFile) -> None:
    # GIVEN a SPIDEX file with scores for the C and G alternatives at 1:100
    # WHEN fetching the score for the G alternative
    score = get_spidex_score(spidex_handle, "1", "100", "G")

    # THEN the z score of the G row is returned as a float
    assert score == 1.25


def test_get_spidex_score_unknown_alt(spidex_handle: pysam.TabixFile) -> None:
    # GIVEN a SPIDEX file without a T alternative at 1:100
    # WHEN fetching the score for the T alternative
    score = get_spidex_score(spidex_handle, "1", "100", "T")

    # THEN there is no score
    assert score is None


def test_get_spidex_score_unknown_position(spidex_handle: pysam.TabixFile) -> None:
    # GIVEN a SPIDEX file without a record at 1:200
    # WHEN fetching the score at that position
    score = get_spidex_score(spidex_handle, "1", "200", "C")

    # THEN there is no score
    assert score is None


def test_get_cosmic(cosmic_handle: pysam.TabixFile) -> None:
    # GIVEN a COSMIC file with a record at 1:100 for the C alternative
    # WHEN checking the C alternative
    in_cosmic = get_cosmic(cosmic_handle, "1", "100", "C")

    # THEN the variant is in COSMIC
    assert in_cosmic is True


def test_get_cosmic_unknown_alt(cosmic_handle: pysam.TabixFile) -> None:
    # GIVEN a COSMIC file with a record at 1:100 for the C alternative only
    # WHEN checking the G alternative
    in_cosmic = get_cosmic(cosmic_handle, "1", "100", "G")

    # THEN the variant is not in COSMIC
    assert in_cosmic is False


def test_get_cosmic_unknown_position(cosmic_handle: pysam.TabixFile) -> None:
    # GIVEN a COSMIC file without a record at 1:200
    # WHEN checking that position
    in_cosmic = get_cosmic(cosmic_handle, "1", "200", "C")

    # THEN the variant is not in COSMIC
    assert in_cosmic is False
