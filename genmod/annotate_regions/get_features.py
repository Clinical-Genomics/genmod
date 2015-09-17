#!/usr/bin/env python
# encoding: utf-8
"""
get_features.py

Check exon trees and gene trees which features that are overlapped by a variant.

Created by Måns Magnusson on 2015-09-16.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import os
import logging

from genmod.annotate_regions import check_overlap



def check_exonic(chrom, start, stop, exon_trees):
    """
    Check if variant overlapps any exon
    
    Arguments:
       chrom (str): The chromosome
       start (int): The start position for the feature
       stop (int): The stop position for the feature
       exon_trees (dict): A dictionary with chromosomes as keys and interval trees as values
    
    Returns:
        bool : If the variant overlapps any exons
    """
    logger = logging.getLogger(__name__)
    
    variant_interval = [start, stop]
    
    exon_features = None
    
    try:
        exon_tree = exon_trees[chrom]
        exon_features = check_overlap(variant_interval, exon_tree)
    except KeyError:
        logger.warning("Chromosome {0} is not in annotation file".format(chrom))
    
    if exon_features:
        return True
    
    return False


def get_genes(chrom, start, stop, gene_trees):
    """Check if variant overlapps genes
    
    Arguments:
       chrom (str): The chromosome
       start (int): The start position for the feature
       stop (int): The stop position for the feature
       gene_trees (dict): A dictionary with chromosomes as keys and interval trees as values
    
    Returns:
        gene_features (set): The genes that the variant ovelapps
    """
    logger = logging.getLogger(__name__)
    
    variant_interval = [start, stop]
    gene_features = set()
    
    try:
        gene_tree = gene_trees[chrom]
        gene_features = check_overlap(variant_interval, gene_tree)
    except KeyError:
        logger.warning("Chromosome {0} is not in annotation file".format(chrom))
    
    return gene_features
    
    
    
