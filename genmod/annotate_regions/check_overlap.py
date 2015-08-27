#!/usr/bin/env python
# encoding: utf-8
"""
check_overlap.py

Function to check what intervals a that are ovelapped by a region
 
Created by Måns Magnusson on 2015-08-25.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

import logging

def check_overlap(interval, interval_tree):
    """
    Function to test what features that are overlapped by a region in a
    interval tree
    
    Arguments:
        interval (list): A list with integers that defines start and stop
        interval_tree (IntervalTree): An interval tree with regions
    
    Returns:
        features_overlapped (set): A set with id:s of the features that are
                                    overlapped
    """
    logger = logging.getLogger(__name__)
    logger.debug("Finding overlap for interval {0}".format(interval))
    
    return set(interval_tree.find_range(interval))
    