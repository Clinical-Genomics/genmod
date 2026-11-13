#!/usr/bin/env python
# encoding: utf-8
"""
annotation_parser.py

This script will parse a file with intervals in .bed format and build
one region index for each chromosome.
These intervals typically represents genes

So self.chromosomes will look like:

{'1':RegionIndex, '2':RegionIndex, ..., 'X':RegionIndex}

The intervals represent features that are annotated in the infile.

Created by Måns Magnusson on 2016-12-23.
Copyright (c) 2016 __MonsoInc__. All rights reserved.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Iterator
from typing import NamedTuple

from genmod.annotate_regions.region_index import Region, RegionIndex, RegionTrees

logger = logging.getLogger(__name__)


class BedRegion(NamedTuple):
    """A region from a bed file, with the chromosome it is on."""

    chrom: str
    start: int
    end: int
    symbol: str


def build_region_trees(bed_lines: Iterable[str], padding: int) -> RegionTrees:
    """Build a region index for each chromosome

    Build a dictionary with chromosomes as keys and region indexes as
    values. Regions may overlap each other.

    Args:
        bed_lines: An iterable with bed formated lines
        padding: Defines what should be considered upstream and downstream variants

    Returns:
        The region index of each chromosome, keyed on chromosome
    """
    regions_by_chrom: dict[str, list[Region]] = {}
    for bed_region in bed_parser(bed_lines, padding):
        regions_by_chrom.setdefault(bed_region.chrom, []).append(
            Region(bed_region.start, bed_region.end, bed_region.symbol)
        )

    return {chrom: RegionIndex(regions) for chrom, regions in regions_by_chrom.items()}


def bed_parser(bed_lines: Iterable[str], padding: int = 4000) -> Iterator[BedRegion]:
    """Parse a file in the bed format.

    Comment lines (starting with ``#``) and empty lines are skipped. The chromosome is
    stored without a ``chr`` prefix and the mitochondrial chromosome (``MT``, ``M``) as
    ``M``. Regions are extended ``padding`` bases upstream, but not on ``M``, and never
    below position 0. If a line has no fourth column the line index is used as symbol.

    Args:
        bed_lines: An iterable with bed formated lines
        padding: Defines what should be considered upstream and downstream variants

    Yields:
        The parsed regions
    """
    for index, line in enumerate(bed_lines):
        if line.startswith("#") or len(line) <= 1:
            continue
        columns = line.rstrip().split()
        # Get the coordinates for the region:
        chrom = columns[0].lstrip("chr")
        feature_stop = int(columns[2])
        if chrom in ["MT", "M"]:
            # Only represent the mitochondrial chromosome as M
            chrom = "M"
            feature_start = int(columns[1])
        else:
            feature_start = max(int(columns[1]) - padding, 0)

        # Get the feature id
        feature_id = columns[3] if len(columns) > 3 else str(index)

        yield BedRegion(chrom, feature_start, feature_stop, feature_id)
