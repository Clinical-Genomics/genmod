#!/usr/bin/env python
# encoding: utf-8
"""
load_annotation.py

Load annotations built by genmod.
Return gene_trees and exon_trees which are interval trees.

Created by Måns Magnusson on 2015-01-12.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os

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
    
    gene_trees = {}
    exon_trees = {}
    
    if verbose:
        print('Reading annotations...\n', file=sys.stderr)
    
    gene_db = os.path.join(annotation_dir, 'genes.db')
    exon_db = os.path.join(annotation_dir, 'exons.db')
    
    try:
        with open(gene_db, 'rb') as f:
            gene_trees = pickle.load(f)
        with open(exon_db, 'rb') as g:
            exon_trees = pickle.load(g)
    except IOError as e:
        if verbose:
            warning('No annotations found.')
            warning('You need to build annotations! See documentation.')
            # It is possible to continue the analysis without annotation files
        pass
    
    if verbose:
        print('Annotations used found in: %s, %s\n' % (gene_db, exon_db), file=sys.stderr)
     
    return gene_trees, exon_trees