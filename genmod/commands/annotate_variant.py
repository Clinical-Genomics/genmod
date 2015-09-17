#!/usr/bin/env python
# encoding: utf-8
"""
annotate_variants.py

Command line tool for annotating vcf variants.
How and what that should be annotated is specified on the command line or with a config file.


Created by Måns Magnusson on 2015-08-25.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import (print_function)

import sys
import logging
import pkg_resources

import click
import tabix

from codecs import open

from genmod import __version__
# from genmod.annotate_regions import load_annotations, check_overlap
from genmod.vcf_tools import (HeaderParser, add_vcf_info, add_metadata, print_headers,
                              print_variant)

from genmod.annotate_variants import annotate_thousand_g
from genmod.annotate_regions import get_genes, check_exonic

@click.command()
@click.argument('variant_file',
                    nargs=1,
                    type=click.File('r'),
                    metavar='<vcf_file> or -'
)
@click.option('-r', '--annotate_regions', 
                is_flag=True,
                help='Increase output verbosity.'
)
@click.option('--cadd_file', 
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped cadd file (with index) with variant scores."""
)
@click.option('--cadd_1000g',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped cadd file (with index) with variant scores 
                            for all 1000g variants."""
)
@click.option('--cadd_exac',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped cadd file (with index) with variant scores 
                            for all ExAC variants."""
)
@click.option('--cadd_esp',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped cadd file (with index) with variant scores 
                            for all ESP6500 variants."""
)
@click.option('--cadd_indels',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped cadd file (with index) with variant scores 
                            for all CADD InDel variants."""
)
@click.option('--thousand_g',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped vcf file (with index) with 1000g variants"""
)
@click.option('--exac',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped vcf file (with index) with exac variants."""
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
@click.option('-s', '--silent',
                is_flag=True,
                help='Do not print the variants.'
)
def annotate(variant_file, annotate_regions, cadd_file, cadd_1000g, 
cadd_exac, cadd_esp, cadd_indels, thousand_g, exac, annotation_dir, outfile, silent):
    """
    Annotate vcf variants.
    
    Annotate variants with a number of different sources.
    Please use --help for more info.
    """

    logger = logging.getLogger(__name__)
    #For testing
    logger = logging.getLogger("genmod.commands.annotate_variants")
    
    logger.info("Running genmod annotate_variant version {0}".format(__version__))

    logger.info("Initializing a Header Parser")
    head = HeaderParser()
    
    for line in variant_file:
        line = line.rstrip()

        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            break
    
    if annotate_regions:
        logger.info("Loading annotations")
        gene_trees, exon_trees = load_annotations(annotation_dir)
        
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
        
    variant_file.seek(0)
    
    any_cadd_file = False

    if exac:
        logger.debug("Opening EXaC frequency file with tabix open")
        exac_handle = tabix.open(exac)
        logger.debug("EXaC frequency file opened")
    
    if thousand_g:
        logger.debug("Opening 1000G frequency file with tabix open")
        thousand_g_handle = tabix.open(thousand_g)
        logger.debug("1000G frequency file opened")

    if cadd_file:
        logger.debug("Opening CADD file with tabix open")
        cadd_handle = tabix.open(cadd_file)
        any_cadd_file = True
        logger.debug("CADD file opened")

    if cadd_1000g:
        logger.debug("Opening CADD 1000G file with tabix open")
        cadd_thousand_g_handle = tabix.open(cadd_1000g)
        any_cadd_file = True
        logger.debug("CADD 1000G file opened")

    if cadd_exac:
        logger.debug("Opening CADD EXaC file with tabix open")
        cadd_exac_handle = tabix.open(cadd_exac)
        any_cadd_file = True
        logger.debug("CADD EXaC file opened")

    if cadd_esp:
        logger.debug("Opening CADD ESP file with tabix open")
        cadd_esp_handle = tabix.open(cadd_esp)
        any_cadd_file = True
        logger.debug("CADD ESP file opened")

    if cadd_indels:
        logger.debug("Opening CADD INDELS file with tabix open")
        cadd_indels_handle = tabix.open(cadd_indels)
        any_cadd_file = True
        logger.debug("CADD INDELS file opened")
    
    for line in variant_file:
        line = line.rstrip()
        
        if not line.startswith('#'):
            splitted_line = line.split()
            chrom = variant_line[0].strip('chr')
            position = int(variant_line[1])
            ref = splitted_line[3]
            alternatives = splitted_line[4]
            
            longest_alt = max([
                len(alt) for alt in alternatives.split(',')])
            
            if annotate_regions:
                if check_exonic(
                    chrom = chrom, 
                    start = position, 
                    stop = (position+longest_alt)-1, 
                    exon_trees = exon_trees):
                    line = add_vcf_info(
                        keyword = 'Exonic', 
                        variant_line=line, 
                        annotation=None
                    )
                
                genes = get_genes(
                    chrom = chrom, 
                    start = position, 
                    stop = (position+longest_alt)-1, 
                    gene_trees = gene_trees
                )
                if genes:
                    line = add_vcf_info(
                        keyword = "Annotation",
                        variant_line = line,
                        annotation = ','.join(gene_features)
                    )
                    
                
            # annotated_line = annotate_thousand_g(
            #     variant_line = line,
            #     thousand_g = thousand_g_handle
            # )
            print(line)
    #         if not headers_done:
    #             add_annotation_header(head)
    #             add_metadata(
    #                 head,
    #                 'info',
    #                 'Exonic',
    #                 annotation_number='0',
    #                 entry_type='Flag',
    #                 description='Indicates if the variant is exonic.'
    #             )
    #             logger.info("Printing vcf header")
    #             print_headers(head, outfile)
    #
    #             headers_done = True
    #
    #         gene_features = None
    #         exon_features = None
    #
    #         variant_line = line.rstrip().split('\t')
    #
    #         chrom = variant_line[0].strip('chr')
    #         position = int(variant_line[1])
    #         longest_alt = max([
    #             len(alt) for alt in variant_line[4].split(',')])
    #         variant_interval = [position, (position + longest_alt-1)]
    #
    #         try:
    #             gene_tree = gene_trees[chrom]
    #             gene_features = check_overlap(variant_interval, gene_tree)
    #         except KeyError:
    #             logger.warning("Chromosome {0} is not in annotation file".format(chrom))
    #
    #         try:
    #             exon_tree = exon_trees[chrom]
    #             exon_features = check_overlap(variant_interval, exon_tree)
    #         except KeyError:
    #             logger.warning("Chromosome {0} is not in annotation file".format(chrom))
    #
    #         if gene_features:
    #             line = add_vcf_info(
    #                 variant_line = line,
    #                 keyword = keyword,
    #                 annotation = ','.join(gene_features)
    #             )
    #         if exon_features:
    #             line = add_vcf_info(
    #                 variant_line = line,
    #                 keyword = "Exonic"
    #             )
    #
    #         print_variant(
    #             variant_line = line,
    #             outfile = outfile
    #         )

if __name__ == '__main__':
    from genmod.log import init_log
    from genmod import logger
    init_log(logger, loglevel="INFO")
    annotate_variants()
