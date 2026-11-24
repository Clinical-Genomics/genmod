"""A small index for finding which (possibly overlapping) regions overlap a position or range."""

from __future__ import annotations

from bisect import bisect_left
from collections.abc import Iterable
from typing import NamedTuple


class Region(NamedTuple):
    """A half-open region ``[start, end)`` with a name, e.g. a gene from a BED file."""

    start: int
    end: int
    symbol: str


class RegionIndex:
    """Static index of regions on one chromosome that answers overlap queries.

    Regions may overlap each other (nested or partly overlapping genes) and every
    region overlapping a query is returned. The index is built once and not changed.

    Regions are sorted on ``start``. For each position in that order we also keep the
    largest ``end`` seen so far (``max_end``). A query for ``[start, end)`` then:

    1. uses ``bisect`` to find the regions that begin before ``end`` (only those can overlap),
    2. walks backwards through them and stops as soon as ``max_end`` shows that no
       earlier region reaches ``start``.

    Lookups are fast for ordinary gene annotations. A very long region keeps ``max_end``
    high and makes queries after it walk further back, which is still correct.
    """

    def __init__(self, regions: Iterable[Region]) -> None:
        """Build the index.

        Args:
            regions: The regions to index, in any order.

        Raises:
            ValueError: If a region does not have ``start < end``.
        """
        sorted_regions = sorted(regions, key=lambda region: region.start)
        for region in sorted_regions:
            if not region.start < region.end:
                raise ValueError(
                    f"Region {region.symbol} must have start < end, got {region.start}-{region.end}"
                )
        self._regions: list[Region] = sorted_regions
        self._starts: list[int] = [region.start for region in sorted_regions]
        # _max_end[i] is the largest end among regions 0..i
        self._max_end: list[int] = []
        largest = sorted_regions[0].end if sorted_regions else 0
        for region in sorted_regions:
            largest = max(largest, region.end)
            self._max_end.append(largest)

    def __len__(self) -> int:
        return len(self._regions)

    def query(self, start: int, end: int) -> list[Region]:
        """Return all regions overlapping the half-open range ``[start, end)``.

        A range without length (``end <= start``) overlaps nothing.

        Args:
            start: First position of the range (inclusive).
            end: End of the range (exclusive).

        Returns:
            The overlapping regions, in order of region start.
        """
        if end <= start:
            return []
        # Regions starting at or after `end` cannot overlap
        upper = bisect_left(self._starts, end)
        hits: list[Region] = []
        for i in range(upper - 1, -1, -1):
            if self._max_end[i] <= start:
                # No region at or before i reaches `start`
                break
            if self._regions[i].end > start:
                hits.append(self._regions[i])
        hits.reverse()
        return hits

    def at(self, position: int) -> list[Region]:
        """Return all regions containing a single position.

        Args:
            position: The position to look up.

        Returns:
            The regions with ``start <= position < end``.
        """
        return self.query(position, position + 1)


# One index per chromosome, e.g. {"1": RegionIndex, "X": RegionIndex}
RegionTrees = dict[str, RegionIndex]
