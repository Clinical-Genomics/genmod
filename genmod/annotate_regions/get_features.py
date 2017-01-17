#!/usr/bin/env python
# encoding: utf-8
"""
get_features.py

Check exon trees and gene trees which features that are overlapped by a variant.

Created by Måns Magnusson on 2015-09-16.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

import logging

logger = logging.getLogger(__name__)

def get_region(chrom, start, end, region_trees):
    """Check if a position overlapps any regions
    
    Arguments:
       chrom (str): The chromosome
       start (int): The start position for the feature
       end (int): The stop position for the feature
       region_trees (dict): A dictionary with chromosomes as keys and interval trees as values
    
    Returns:
        regions (set): The regions that the variant ovelapps
    """
    regions = set()
    if chrom in region_trees:
        tree = region_trees[chrom]
        result = tree[start:end]
        for interval in result:
            regions.add(interval.data)
    return regions
