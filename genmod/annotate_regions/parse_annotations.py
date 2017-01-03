#!/usr/bin/env python
# encoding: utf-8
"""
annotation_parser.py

This script will parse a file with intervals in .bed format and build 
one interval tree for each chromosome.
These intervals typically represents genes

So self.chromosomes will look like:

{'1':intervalTree, '2':intervalTree, ..., 'X':intervalTree}

The intervals represent features that are annotated in the infile.

Created by Måns Magnusson on 2016-12-23.
Copyright (c) 2016 __MonsoInc__. All rights reserved.
"""
import logging

try:
    import cPickle as pickle
except:
    import pickle

from intervaltree import (Interval, IntervalTree)

logger = logging.getLogger(__name__)

def get_interval(start, stop, value):
    """Create an interval instance
    
        Args:
            start(int)
            stop(int)
            value
        
        Returns:
            interval(intervaltree.Interval)
    
    """
    interval = Interval(start, stop, value)
    return interval

def build_region_trees(bed_lines, padding):
    """Build region trees for each chromosome
    
        Build a dictionary with chromosomes as keys and interval trees as 
        values.
    
        Args:
            bed_lines(iterable): An iterable with bed formated lines
            padding (int): Defines what should be considered upstream 
                           and downstream variants
    """
    region_trees = {}
    for region in bed_parser(bed_lines, padding):
        chrom = region['chrom']
        start = region['start']
        stop = region['stop']
        symbol = region['symbol']

        if chrom not in region_trees:
            region_trees[chrom] = IntervalTree()

        region_trees[chrom].add(Interval(start, stop, symbol))

    return region_trees


def bed_parser(bed_lines, padding=4000):
    """
    Parse a file in the bed format.
    
    Arguments:
        bed_lines(iterable): An iterable with bed formated lines
        padding (int): Defines what should be considered upstream 
                       and downstream variants
    
    Yields:
        region(dict): 
                    {
                        'chrom': str,
                        'start': int,
                        'stop': int,
                        'symbol': str
                    }
    """
    genes = {}
    for index, line in enumerate(bed_lines):
        if not line.startswith('#') and len(line) > 1:
            line = line.rstrip().split()
            feature_id = str(index)
            # Get the coordinates for the region:
            chrom = line[0].lstrip('chr')
            if chrom == 'MT':
                feature_start = int(line[1])
                feature_stop = int(line[2])
            else:
                feature_start = max(int(line[1]) - padding, 0)
                feature_stop = int(line[2])
                
            # Get the feature id
            if len(line) > 3:
                feature_id = line [3]
            
            region = {
                'chrom': chrom,
                'start': feature_start,
                'stop': feature_stop,
                'symbol': feature_id
            }
            
            yield region
