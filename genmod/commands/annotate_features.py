#!/usr/bin/env python
# encoding: utf-8
"""
annotate_features.py

Command line tool for annotating variants with what features they belong to.
The features are searced in a data structure constructed by 'genmod build_annotations'.
There is a default annotation set that is the refSeq genes, these will be used if no
other is specified as input

Created by Måns Magnusson on 2015-08-25.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import (print_function)

import logging
import pkg_resources

import click

from genmod import __version__
from genmod.annotate_regions import load_annotations, check_overlap
from genmod.vcf_tools import (HeaderParser, add_vcf_info, add_metadata, print_headers,
                              print_variant)

@click.command()
@click.argument('variant_file', 
                    nargs=1, 
                    type=click.File('rb'),
                    metavar='<vcf_file> or -'
)
@click.option('-a' ,'--annotation_dir', 
                    type=click.Path(exists=True),
                    default=pkg_resources.resource_filename('genmod', 'annotations'),
                    help="""Specify the path to the directory where the annotation 
                    databases are. 
                    Default is the gene pred files that comes with the distribution."""
)
@click.option('-o', '--outfile', 
                    type=click.File('w'),
                    help='Specify the path to a file where results should be stored.'
)
@click.option('-k', '--keyword', 
                    default="Annotation",
                    help="What keyword should be used for annotation in vcf."\
                         "Default is Annotation"
)
def annotate_features(variant_file, annotation_dir, outfile, keyword):
    """
    Annotate vcf variants with what features they overlap and if they are exonic.
    """

    logger = logging.getLogger(__name__)
    #For testing
    logger = logging.getLogger("genmod.commands.annotate_features")
    logger.info("Running genmod annotate_features version {0}".format(__version__))

    logger.info("Initializing a Header Parser")
    head = HeaderParser()
    
    headers_done = False
    
    logger.info("Loading annotations")
    
    gene_trees, exon_trees = load_annotations(annotation_dir)

    for line in variant_file:
        line = line.rstrip()
        
        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            if not headers_done:
                add_metadata(
                    head,
                    'info',
                    'Annotation',
                    annotation_number='.',
                    entry_type='String',
                    description='Annotates what feature(s) this variant belongs to.'
                )
                add_metadata(
                    head,
                    'info',
                    'Exonic',
                    annotation_number='0',
                    entry_type='Flag',
                    description='Indicates if the variant is exonic.'
                )
                logger.info("Printing vcf header")
                print_headers(head, outfile)
                
                headers_done = True
            
            gene_features = None
            exon_features = None

            variant_line = line.rstrip().split('\t')

            chrom = variant_line[0].strip('chr')
            position = int(variant_line[1])
            longest_alt = max([
                len(alt) for alt in variant_line[4].split(',')])
            variant_interval = [position, (position + longest_alt-1)]

            try:
                gene_tree = gene_trees[chrom]
                gene_features = check_overlap(variant_interval, gene_tree)
            except KeyError:
                logger.warning("Chromosome {0} is not in annotation file".format(chrom))

            try:
                exon_tree = exon_trees[chrom]
                exon_features = check_overlap(variant_interval, exon_tree)
            except KeyError:
                logger.warning("Chromosome {0} is not in annotation file".format(chrom))

            if gene_features:
                line = add_vcf_info(
                    variant_line = line,
                    keyword = keyword,
                    annotation = ','.join(gene_features)
                )
            if exon_features:
                line = add_vcf_info(
                    variant_line = line,
                    keyword = "Exonic"
                )
            
            print_variant(
                variant_line = line,
                outfile = outfile
            )

if __name__ == '__main__':
    from genmod.log import init_log
    from genmod import logger
    init_log(logger, loglevel="INFO")
    annotate_features()
