#!/usr/bin/env python
# encoding: utf-8
"""
build_annotation.py

Command line tool for building new annotation databases for genmod.

Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import sys
import os

import logging
import click
import pkg_resources

try:
    import cPickle as pickle
except:
    import pickle

from genmod.annotation import AnnotationParser


###        This is for building new annotations     ###

@click.command()
@click.argument('annotation_file', 
                nargs=1, 
                type=click.Path(exists=True),
)
@click.option('-t' ,'--annotation_type',
                type=click.Choice(['bed', 'ccds', 'gtf', 'gene_pred', 'gff']), 
                default='gene_pred',
                help='Specify the format of the annotation file.'
)
@click.option('-o', '--outdir', 
                    type=click.Path(exists=True),
                    default=pkg_resources.resource_filename('genmod', 'annotations'),
                    help=("""Specify the path to a folder where the annotation files should be stored. 
                            Default is the annotations dir of the ditribution.""")
)
@click.option('--splice_padding',
                    type=int, nargs=1, default=2,
                    help='Specify the the number of bases that the exons should be padded with. Default is 2 bases.'
)
@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def build_annotation(annotation_file, annotation_type, outdir, splice_padding, verbose):
    """Build a new annotation database."""
    logger = logging.getLogger(__name__)
    logger.info("Building new annotation databases from {0} into {1}.".format(
        annotation_file, outdir))
    
    logger.info("Initializing annotation parser.")
    
    anno_parser = AnnotationParser(
                            annotation_file, 
                            annotation_type, 
                            splice_padding = splice_padding, 
                            verbosity=verbose
                        )
    logger.info("Creating files for dumping the data structures")
    gene_db = os.path.join(outdir, 'genes.db')
    exon_db = os.path.join(outdir, 'exons.db')
    logger.info("Gene annotation files will be dumped on {0}".format(gene_db))
    logger.info("Exon annotation files will be dumped on {0}".format(exon_db))
    
    
    with open(gene_db, 'wb') as f:
        logger.info("Dumping gene database to {0}.".format(gene_db))
        pickle.dump(anno_parser.gene_trees, f)
        logger.debug("Dumping successful.")
    
    with open(exon_db, 'wb') as g:
        logger.info("Dumping exon database to {0}.".format(exon_db))
        pickle.dump(anno_parser.exon_trees, g)
        logger.debug("Dumping successful.")
    

if __name__ == '__main__':
    build_annotation()
