"""Tests for RegionIndex: overlap queries on half-open [start, end) regions that may overlap."""

import random

import pytest
from genmod.annotate_regions.region_index import Region, RegionIndex


def symbols(regions: list[Region]) -> set[str]:
    return {region.symbol for region in regions}


def brute_force(regions: list[Region], start: int, end: int) -> set[str]:
    """Reference implementation: check every region. A range without length matches nothing."""
    if end <= start:
        return set()
    return {r.symbol for r in regions if r.start < end and r.end > start}


# Four genes: B is nested in A, C partly overlaps A, D is separate
GENES = [
    Region(100, 500, "A"),
    Region(200, 300, "B"),
    Region(400, 800, "C"),
    Region(1000, 1200, "D"),
]


def test_empty_index_has_no_hits():
    # GIVEN an index without regions
    index = RegionIndex([])

    # WHEN querying
    # THEN nothing is found and there is no error
    assert index.query(0, 1000) == []
    assert len(index) == 0


def test_single_region_hit():
    # GIVEN an index with one region
    index = RegionIndex([Region(100, 200, "A")])

    # WHEN querying inside it
    # THEN the region is returned
    assert symbols(index.query(150, 151)) == {"A"}


def test_region_start_is_inclusive_and_end_is_exclusive():
    # GIVEN an index with the region [100, 200)
    index = RegionIndex([Region(100, 200, "A")])

    # THEN the first base is in the region and the end coordinate is not
    assert symbols(index.at(100)) == {"A"}
    assert symbols(index.at(199)) == {"A"}
    assert index.at(99) == []
    assert index.at(200) == []


def test_adjacent_regions_do_not_overlap():
    # GIVEN two regions that touch at 200
    index = RegionIndex([Region(100, 200, "A"), Region(200, 300, "B")])

    # THEN position 200 is only in the second region
    assert symbols(index.at(200)) == {"B"}
    assert symbols(index.at(199)) == {"A"}
    # AND a query ending at 200 does not reach the second region
    assert symbols(index.query(150, 200)) == {"A"}


def test_nested_regions_return_both():
    # GIVEN the GENES, where B is nested in A
    index = RegionIndex(GENES)

    # WHEN querying a position in B
    # THEN both the inner and the outer gene are returned
    assert symbols(index.at(250)) == {"A", "B"}


def test_partially_overlapping_regions_return_both():
    # GIVEN the GENES, where C partly overlaps A (400-500)
    index = RegionIndex(GENES)

    # THEN a position in the overlap is in both genes
    assert symbols(index.at(450)) == {"A", "C"}
    # AND positions outside the overlap are only in one
    assert symbols(index.at(350)) == {"A"}
    assert symbols(index.at(600)) == {"C"}


def test_variant_walk_back_stops_without_missing_earlier_long_region():
    # GIVEN a long region followed by short ones that end before the query
    regions = [Region(0, 10_000, "LONG"), Region(100, 110, "S1"), Region(200, 210, "S2")]
    index = RegionIndex(regions)

    # WHEN querying after the short regions
    # THEN the long region is still found
    assert symbols(index.at(5000)) == {"LONG"}


def test_query_between_and_outside_regions_is_empty():
    # GIVEN the GENES
    index = RegionIndex(GENES)

    # THEN gaps, before the first and after the last region are empty
    assert index.at(900) == []
    assert index.at(50) == []
    assert index.at(5000) == []
    assert index.query(10_000, 11_000) == []


def test_range_query_returns_all_overlapping_regions():
    # GIVEN the GENES
    index = RegionIndex(GENES)

    # WHEN querying a range spanning several genes
    # THEN all genes touching the range are returned
    assert symbols(index.query(250, 1100)) == {"A", "B", "C", "D"}
    assert symbols(index.query(300, 400)) == {"A"}


@pytest.mark.parametrize("start, end", [(100, 100), (500, 400), (-5, -5)])
def test_empty_or_inverted_query_matches_nothing(start, end):
    # GIVEN the GENES
    index = RegionIndex(GENES)

    # WHEN querying without length
    # THEN nothing matches, as for the previous interval tree
    assert index.query(start, end) == []


def test_identical_regions_are_both_returned():
    # GIVEN two regions with the same coordinates but different symbols
    index = RegionIndex([Region(100, 200, "A"), Region(100, 200, "B")])

    # THEN both are found
    assert symbols(index.at(150)) == {"A", "B"}


def test_duplicate_symbol_on_different_regions():
    # GIVEN one symbol on two separate regions (e.g. two exons)
    index = RegionIndex([Region(100, 200, "A"), Region(300, 400, "A")])

    # THEN each region is found on its own
    assert len(index.at(150)) == 1
    assert len(index.at(350)) == 1
    assert len(index) == 2


def test_input_order_does_not_matter():
    # GIVEN the same regions in reverse order
    index = RegionIndex(list(reversed(GENES)))

    # THEN the result is the same
    assert symbols(index.at(450)) == {"A", "C"}
    assert symbols(index.at(250)) == {"A", "B"}


def test_negative_start_coordinates_are_supported():
    # GIVEN a region starting below zero
    index = RegionIndex([Region(-100, 50, "A")])

    # THEN it is found
    assert symbols(index.at(-1)) == {"A"}
    assert symbols(index.at(0)) == {"A"}


@pytest.mark.parametrize("start, end", [(100, 100), (200, 100)])
def test_region_without_length_is_rejected(start, end):
    # GIVEN a region with end <= start, as the previous interval tree rejected
    # WHEN building the index
    # THEN a clear error is raised
    with pytest.raises(ValueError, match="start < end"):
        RegionIndex([Region(start, end, "BAD")])


def test_matches_brute_force_on_random_overlapping_regions():
    # GIVEN many random regions, with plenty of overlap and some very long ones
    rng = random.Random(1234)
    regions = []
    for i in range(300):
        start = rng.randint(0, 2000)
        length = rng.choice([1, 5, 20, 100, 1000, 5000])
        regions.append(Region(start, start + length, f"G{i}"))
    index = RegionIndex(regions)

    # WHEN querying many random positions and ranges
    for _ in range(2000):
        start = rng.randint(-50, 3000)
        end = start + rng.choice([0, 1, 1, 2, 10, 500])
        # THEN the result is the same as checking every region
        assert symbols(index.query(start, end)) == brute_force(regions, start, end)
