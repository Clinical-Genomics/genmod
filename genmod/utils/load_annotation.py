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

import sys
import os
import logging

try:
    import cPickle as pickle
except:
    import pickle

from codecs import open

from genmod.errors import warning


def load_annotations(annotation_dir, verbose=False):
    """
    Load the annotations found in the indata path.
    These are pickled interval trees that are returned as dictionaries.
    """
    logger = logging.getLogger(__name__)
    gene_trees = {}
    exon_trees = {}
    
    logger.info('Reading annotations...')
    
    gene_db = os.path.join(annotation_dir, 'genes.db')
    exon_db = os.path.join(annotation_dir, 'exons.db')
    
    try:
        with open(gene_db, 'rb') as f:
            gene_trees = pickle.load(f)
        with open(exon_db, 'rb') as g:
            exon_trees = pickle.load(g)
    except IOError as e:
        if verbose:
            logger.warning('No annotations found.')
            logger.warning('You need to build annotations! See documentation.')
            # It is possible to continue the analysis without annotation files
        pass
    
    logger.info('Annotations used found in: {0}, {1}'.format(gene_db, exon_db))
     
    return gene_trees, exon_trees
