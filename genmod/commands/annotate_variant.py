#!/usr/bin/env python
# encoding: utf-8
"""
annotate_variants.py

Command line tool for annotating vcf variants.
How and what that should be annotated is specified on the command line or with a config file.


Created by Måns Magnusson on 2015-08-25.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

import sys
import logging
import pkg_resources
import itertools

import click

from codecs import open
from datetime import datetime

from tabix import TabixError

from genmod import __version__

from genmod.vcf_tools import (HeaderParser, print_headers, print_variant)

from genmod.annotations import (ensembl_path_37, ensembl_path_38)

from genmod.annotate_regions.parse_annotations import (build_region_trees)
from genmod.annotate_variants.add_annotations import (add_regions, add_exac, 
    add_exac_max, add_thousandg, add_thousandg_max, add_spidex, add_cadd, 
    add_cadd_raw, add_cosmic)

from genmod.annotate_variants.read_tabix_files import (get_tabixhandle)
from genmod.annotate_variants.annotate import annotate_variant

from genmod.commands.utils import (outfile, silent, temp_dir, 
                                   variant_file, get_file_handle)

logger = logging.getLogger(__name__)

@click.command('annotate', short_help="Annotate vcf variants")
@variant_file
@click.option('-r', '--regions', '--annotate_regions','--annotate-regions', 
                is_flag=True,
                help='Annotate what regions a variant belongs to (eg. genes).'
)
@click.option('--region-file','--region_file', 
                type=click.Path(exists=True),
                show_default=True,
                help='Choose a bed file with regions that should be used.'
)
@click.option('--genome-build','-b',
                type=click.Choice(['37','38']),
                default='37',
                show_default=True,
                help='Choose what genome build to use.'
)
@click.option('-c', '--cadd-file', '--cadd_file', 
                    multiple = True,
                    type=click.Path(exists=True), 
                    help="Specify the path to a bgzipped cadd file"\
                    " (with index) with variant scores. This command can be"\
                    " used multiple times if multiple cadd files."
)
@click.option('--thousand-g', '--thousand_g',
                    type=click.Path(exists=True), 
                    help="Specify the path to a bgzipped vcf file"\
                            " (with index) with 1000g variants"
)
@click.option('--exac',
                    type=click.Path(exists=True),
                    help="Specify the path to a bgzipped vcf file"\
                            " (with index) with exac variants."
)
@click.option('--cosmic',
                    type=click.Path(exists=True), 
                    help="Specify the path to a bgzipped vcf file"\
                            " (with index) with COSMIC variants."
)
@click.option('--max-af', '--max_af',
                    is_flag=True,
                    help="If the MAX AF should be annotated"
)
@click.option('--spidex',
                    type=click.Path(exists=True), 
                    help="Specify the path to a bgzipped tsv file"\
                            " (with index) with spidex information."
)
@click.option('--cadd-raw', '--cadd_raw',
                    is_flag=True,
                    help="""If the raw cadd scores should be annotated."""
)
@outfile
@silent
@temp_dir
@click.pass_context
def annotate(context, variant_file, annotate_regions, region_file, cadd_file, 
             thousand_g, exac, spidex, outfile, silent, cadd_raw, cosmic, 
             max_af, temp_dir, genome_build):
    """
    Annotate vcf variants.
    
    Annotate variants with a number of different sources.
    Please use --help for more info.
    """
    regions = annotate_regions
    logger.info("Running genmod annotate_variant version {0}".format(__version__))
    
    if not region_file:
        if genome_build == '37':
            region_file = ensembl_path_37
        elif genome_build == '38':
            region_file = ensembl_path_38
    
    start_time_analysis = datetime.now()
    annotation_arguments = {}
    
    variants = get_file_handle(variant_file)
    
    logger.info("Initializing a Header Parser")
    head = HeaderParser()
    
    line = None
    for line in variants:
        line = line.rstrip()

        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            break
    
    #Add the first variant back to the iterator
    # If the vcf has no variants the last line will be a header
    if not line.startswith('#'):
        variants = itertools.chain([line], variants)
    else:
        print_headers(head, outfile, silent)
        sys.exit(0)
    
    header_line = head.header
    annotation_arguments['header_line'] = header_line
    
    try:
        if regions:
            logger.info("Loading annotations")
            logger.info("Use annotations file: {0}".format(region_file))
            add_regions(head)
            regions_handle = get_file_handle(region_file)
            logger.debug("Adding region trees to arguments")
            annotation_arguments['region_trees'] = build_region_trees(regions_handle, padding=4000)
        
        if exac:
            logger.info("Annotating ExAC frequencies")
            logger.debug("Using ExAC file: {0}".format(exac))
            annotation_arguments['exac'] = get_tabixhandle(exac)
            add_exac(head)
        
        if thousand_g:
            logger.info("Annotating 1000G frequencies")
            logger.debug("Using 1000G file: {0}".format(thousand_g))
            annotation_arguments['thousand_g'] = get_tabixhandle(thousand_g)
            add_thousandg(head)
        
        if spidex:
            logger.info("Annotating Spidex z scores")
            logger.debug("Using Spidex file: {0}".format(spidex))
            annotation_arguments['spidex'] = get_tabixhandle(spidex)
            add_spidex(head)
        
        if cadd_file:
            logger.info("Annotating CADD scores")
            logger.debug("Using CADD file(s): {0}".format(', '.join(cadd_file)))
            annotation_arguments['cadd_files'] = [get_tabixhandle(cadd) for cadd in cadd_file]
            
            add_cadd(head)
        
            if cadd_raw:
                annotation_arguments['cadd_raw'] = cadd_raw
                add_cadd_raw(head)
        
        if max_af:
            annotation_arguments['max_af'] = max_af
            if thousand_g:
                add_thousandg_max(head)
            if exac:
                add_exac_max(head)
        
        if cosmic:
            logger.info("Annotating if variant is in COSMIC")
            logger.debug("Using COSMIC file: {0}".format(cosmic))
            annotation_arguments['cosmic'] = get_tabixhandle(cosmic)
            add_cosmic(head)
    except TabixError as err:
        logger.warning(err)
        context.abort()
    
    print_headers(head, outfile, silent)
    
    for variant in variants:
        print_variant(
            variant_line = annotate_variant(variant, annotation_arguments),
            outfile = outfile,
            silent = silent
        )
        
