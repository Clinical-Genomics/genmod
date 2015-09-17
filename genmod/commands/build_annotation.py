#!/usr/bin/env python
# encoding: utf-8
"""
build_annotation.py

Command line tool for building new annotation databases for genmod.

The annotation database is a data structure for fast feature annotation with 
genmod. This gives the user the opportunity to use their own set of genes, or
other features, to annotate their variants with.


Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import sys
import os

import logging
import click
import pkg_resources
import logging

try:
    import cPickle as pickle
except:
    import pickle

from genmod.annotate_regions import parse_annotations


###        This is for building new annotations     ###

@click.command()
@click.argument('annotation_file', 
                nargs=1, 
                type=click.Path(exists=True),
)
@click.option('-o', '--outdir', 
                    type=click.Path(exists=True),
                    help="Specify the path to a folder where the annotation"\
                         " files should be stored."
)
@click.option('-t' ,'--annotation_type',
                type=click.Choice(['bed', 'ccds', 'gtf', 'gene_pred', 'gff']), 
                default='gene_pred',
                help='Specify the format of the annotation file.'
)
@click.option('--splice_padding',
                    type=int, nargs=1, default=2,
                    help='Specify the the number of bases that the exons should be padded with. Default is 2 bases.'
)
@click.option('-v', '--verbose', 
                count=True,
                help='Increase output verbosity.'
)
def build(annotation_file, annotation_type, outdir, splice_padding, verbose):
    """Build a new annotation database
    
        Build an annotation database from an annotation file.
    
    """
    logger = logging.getLogger(__name__)
    
    if not outdir:
        raise IOError("Please give a outdir with flag -o/--outdir")
    
    logger.info("Building new annotation databases from {0} into {1}.".format(
        annotation_file, outdir))
    
    logger.info("Initializing annotation parser.")
    
    gene_trees, exon_trees = parse_annotations(
                            annotation_file, 
                            annotation_type, 
                            splice_padding = splice_padding, 
                        )
    
    logger.info("Creating files for dumping the data structures")
    gene_db = os.path.join(outdir, 'genes.db')
    exon_db = os.path.join(outdir, 'exons.db')
    logger.info("Gene annotation files will be dumped on {0}".format(gene_db))
    logger.info("Exon annotation files will be dumped on {0}".format(exon_db))
    
    
    with open(gene_db, 'wb') as f:
        logger.info("Dumping gene database to {0}.".format(gene_db))
        pickle.dump(gene_trees, f)
        logger.debug("Dumping successful.")
    
    with open(exon_db, 'wb') as g:
        logger.info("Dumping exon database to {0}.".format(exon_db))
        pickle.dump(exon_trees, g)
        logger.debug("Dumping successful.")
    

if __name__ == '__main__':
    build()
