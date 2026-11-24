from genmod.annotate_regions.get_features import get_region
from genmod.annotate_regions.parse_annotations import build_region_trees
from genmod.annotate_regions.region_index import RegionIndex


def test_build_region_trees(bed_lines: list[str]) -> None:
    # GIVEN some lines in the bed format
    # WHEN building region trees
    region_trees = build_region_trees(bed_lines, padding=0)
    # THEN assert there are region indexes in the result
    assert len(region_trees) > 0
    for chrom in region_trees:
        assert isinstance(region_trees[chrom], RegionIndex)


def test_query_region_trees() -> None:
    # GIVEN some bed formated lines
    lines = [
        "13\t1\t100\tHMGA1P6\tENSG00000233440\n",
        "13\t50\t200\tRNY3P4\tENSG00000207157\n",
        "13\t300\t1000\tLINC00362\tENSG00000229483\n",
    ]
    # WHEN building region trees
    region_trees = build_region_trees(lines, padding=0)

    # THEN make sure that we get the correct info when querying
    index = region_trees["13"]

    # There should be two results here
    assert {r.symbol for r in index.at(50)} == {"HMGA1P6", "RNY3P4"}

    result = index.at(100)
    # Intervals are half opened
    assert [r.symbol for r in result] == ["RNY3P4"]

    result = index.query(100, 100)
    # Intervals without lenght does not match anythin
    assert result == []

    result = index.query(250, 600)
    assert [r.symbol for r in result] == ["LINC00362"]

    # Test none overlapping
    assert index.query(10000, 11000) == []


def test_get_region_returns_all_overlapping_symbols() -> None:
    # GIVEN overlapping genes on chromosome 13 and one on X
    lines = [
        "13\t1\t100\tOUTER\tENSG1\n",
        "13\t50\t80\tINNER\tENSG2\n",
        "X\t10\t20\tXGENE\tENSG3\n",
    ]
    region_trees = build_region_trees(lines, padding=0)

    # WHEN looking up a position in both overlapping genes
    # THEN a variant in two genes gets both symbols
    assert get_region("13", 60, 61, region_trees) == {"OUTER", "INNER"}
    # AND a position outside is empty
    assert get_region("13", 500, 501, region_trees) == set()
    # AND an unknown chromosome is empty
    assert get_region("7", 60, 61, region_trees) == set()


def test_padding_extends_regions_upstream_but_not_below_zero() -> None:
    # GIVEN a gene close to the chromosome start
    lines = ["1\t1000\t2000\tGENE\tENSG1\n"]

    # WHEN building with padding larger than the start
    region_trees = build_region_trees(lines, padding=4000)

    # THEN the region starts at 0 and covers upstream positions
    assert get_region("1", 0, 1, region_trees) == {"GENE"}
    assert get_region("1", 1999, 2000, region_trees) == {"GENE"}
    assert get_region("1", 2000, 2001, region_trees) == set()


def test_chr_prefix_and_mitochondrial_are_normalised() -> None:
    # GIVEN lines with a chr prefix and the mitochondrial chromosome named MT
    lines = ["chr1\t10\t20\tG1\tx\n", "MT\t10\t20\tMTG\tx\n"]
    region_trees = build_region_trees(lines, padding=0)

    # THEN chromosomes are stored as '1' and 'M'
    assert set(region_trees) == {"1", "M"}


# def test_build_all(ensembl_file):
#     with open(ensembl_file, 'r') as ensembl_handle:
#         region_trees = build_region_trees(ensembl_handle, padding=4000)
