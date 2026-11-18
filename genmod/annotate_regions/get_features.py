#!/usr/bin/env python
# encoding: utf-8
"""
get_features.py

Check exon trees and gene trees which features that are overlapped by a variant.

Created by Måns Magnusson on 2015-09-16.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping

from genmod.annotate_regions.region_index import RegionIndex

logger = logging.getLogger(__name__)


def get_region(
    chrom: str, start: int, end: int, region_trees: Mapping[str, RegionIndex]
) -> set[str]:
    """Check if a position overlaps any regions

    Args:
        chrom: The chromosome
        start: The start position for the feature (inclusive)
        end: The end position for the feature (exclusive)
        region_trees: The region index of each chromosome, keyed on chromosome

    Returns:
        The symbols of the regions that the variant overlaps, empty if the chromosome
        has no regions
    """
    if chrom not in region_trees:
        return set()
    return {region.symbol for region in region_trees[chrom].query(start, end)}
