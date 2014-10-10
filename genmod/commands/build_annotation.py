#!/usr/bin/env python
# encoding: utf-8
"""
build_annotation.py

Command line tool for building new annotation databases for genmod.

Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import click
import pkg_resources

try:
    import cPickle as pickle
except:
    import pickle

import genmod


###        This is for building new annotations     ###

@click.command()
@click.argument('annotation_file', 
                nargs=1, 
                type=click.Path(exists=True),
)
@click.option('-t' ,'--type',
                type=click.Choice(['bed', 'ccds', 'gtf', 'gene_pred']), 
                default='gene_pred',
                help='Specify the format of the annotation file.'
)
@click.option('-o', '--outdir', 
                    type=click.Path(exists=True),
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
def build_annotation(annotation_file, type, outdir, splice_padding, verbose):
    """Build a new annotation database."""
    if verbose:
        click.echo('Building new annotation databases from %s into %s.' % (annotation_file, outdir))
    
    anno_parser = genmod.annotation_parser.AnnotationParser(annotation_file, type, 
                            splice_padding = splice_padding, verbosity=verbosity)
    
    gene_db = pkg_resources.resource_filename('genmod', 'annotations/genes.db')
    exon_db = pkg_resources.resource_filename('genmod', 'annotations/exons.db')
    
    if outdir:
        gene_db = os.path.join(outdir, 'genes.db')
        exon_db = os.path.join(outdir, 'exons.db')
    
    with open(gene_db, 'wb') as f:
        pickle.dump(anno_parser.gene_trees, f)
    
    with open(exon_db, 'wb') as g:
        pickle.dump(anno_parser.exon_trees, g)
    

if __name__ == '__main__':
    build_annotation()