#!/usr/bin/env python
# encoding: utf-8
"""
load_annotation.py

Load annotations built by genmod.
Return gene_trees and exon_trees which are interval trees.

Created by Måns Magnusson on 2015-01-12.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import os
import logging

try:
    import cPickle as pickle
except:
    import pickle

from codecs import open


def load_annotations(annotation_dir):
    """
    Load the annotations found in the indata path.
    These are pickled interval trees that are returned as dictionaries.
    
    Arguments:
        annotation_dir (str): Path to folder where the annotation trees are
    
    Returns:
        gene_trees (dict): A dictionary with {<chr>: IntervalTree}
        exon_trees (dict): A dictionary with {<chr>: IntervalTree}
    """
    logger = logging.getLogger(__name__)
    gene_trees = {}
    exon_trees = {}

    logger.info('Reading annotations...')

    gene_db = os.path.join(annotation_dir, 'genes.db')
    logger.debug("Reading gene trees from {0}".format(gene_db))

    exon_db = os.path.join(annotation_dir, 'exons.db')
    logger.debug("Reading exon trees from {0}".format(exon_db))

    try:
        with open(gene_db, 'rb') as f:
            gene_trees = pickle.load(f)
            logger.info("Gene trees loaded")
        with open(exon_db, 'rb') as g:
            exon_trees = pickle.load(g)
            logger.info("Exon trees loaded")
    except IOError as e:
        logger.warning("No annotations found. You need to build annotations!"\
                        " See documentation.")
        raise e
 
    return gene_trees, exon_trees
