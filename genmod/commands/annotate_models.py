#!/usr/bin/env python
# encoding: utf-8
"""
annotate_models.py

Command line tool for annotating what genetic models that are followed for variants in the vcf format.
The variants needs to be annotated with what genetic regions they belong to and preferably if they are exonic or not.
This can be done with 'genmod annotate_regions', vep or other tools. Please specify on the command line.

Created by Måns Magnusson on 2015-08-25.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import (print_function)

import sys
import os
import click
import inspect
import logging

from multiprocessing import JoinableQueue, Manager, cpu_count
from codecs import open
from datetime import datetime
from tempfile import mkdtemp, TemporaryFile, NamedTemporaryFile

import shutil
import pkg_resources

from ped_parser import FamilyParser

from genmod import (__version__, VariantPrinter)

from genmod.log import init_log
from genmod.variant_annotation import (VariantAnnotator)
from genmod.vcf_tools import (HeaderParser, add_vcf_info, add_version_header, 
add_genetic_models_header, add_model_score_header, print_headers, print_variant,
add_compounds_header,)
from genmod.utils import check_individuals

@click.command()
@click.argument('variant_file', 
                    nargs=1, 
                    type=click.File('rb'),
                    metavar='<vcf_file> or -'
)
@click.option('-f', '--family_file',
                    nargs=1, 
                    type=click.File('r'),
                    metavar='<ped_file>'
)
@click.option('-t' ,'--family_type', 
                type=click.Choice(['ped', 'alt', 'cmms', 'mip']), 
                default='ped',
                help='If the analysis use one of the known setups, please specify which one.'
)
@click.option('-a', '--annotation',
                    default="Annotation",
                    help="Specify what field to search for annotations in."\
                         "Default is Annotation"
)
@click.option('--vep', 
                    is_flag=True,
                    help='If variants are annotated with the Variant Effect Predictor.'
)
@click.option('--phased', 
                    is_flag=True,
                    help='If data is phased use this flag.'
)
@click.option('-s' ,'--strict', 
                    is_flag=True,
                    help='If strict model annotations should be used(see documentation).'
)
@click.option('-split' ,'--split_variants', 
                    is_flag=True,
                    help='If the variants should be splitted.'
)
@click.option('--silent', 
                    is_flag=True,
                    help='Do not print the variants.'
)
@click.option('-g' ,'--whole_gene', 
                    is_flag=True,
                    help="""If compounds should be checked in the whole gene regions. 
                    Not only exonic/splice sites."""
)
@click.option('-o', '--outfile', 
                    type=click.File('w'),
                    help='Specify the path to a file where results should be stored.'
)
def annotate_models(variant_file, family_file, family_type, annotation, vep,
                    split_variants, phased, strict, silent, whole_gene, outfile):
    """
    Annotate variants with what genetic models that are followed in a VCF file.
    The analysis is family based so each family that are specified in the family
    file and exists in the variant file will get it's own annotation.
    """
    logger = logging.getLogger(__name__)
    # For testing only:
    logger = logging.getLogger("genmod.commands.annotate")
    init_log(logger, loglevel="INFO")
    
    ######### This is for logging the command line string #########
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [
        i+'='+str(values[i]) for i in values if values[i] and 
        i not in ['frame']
    ]
    annotator_arguments = {}
    annotator_arguments['phased'] = phased
    annotator_arguments['strict'] = strict
    ##############################################################
    
    logger.info("Running GENMOD annotate version {0}".format(__version__))
    logger.debug("Arguments: {0}".format(', '.join(argument_list)))
    
    if not family_file:
        logger.warning("Please provide a family file with -f/--family_file")
        sys.exit(0)
    
    logger.info("Setting up a family parser")
    family_parser = FamilyParser(family_file, family_type)
    logger.debug("Family parser done")
    logger.info("Families used in analysis: {0}".format(
                    ','.join(list(family_parser.families.keys()))))
    logger.info("Individuals included in analysis: {0}".format(
                    ','.join(list(family_parser.individuals.keys()))))
    
    
    if vep:
        logger.info("Using VEP annotation")
    
    # The individuals in the ped file must be present in the variant file:
    families = family_parser.families
    
    start_time_analysis = datetime.now()
    
    head = HeaderParser()
    
    headers_done = False
    
    for line in variant_file:
        line = line.rstrip()
        
        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            if not headers_done:
                add_version_header(
                    head, 
                    command_line_string=' '.join(argument_list)
                )
                add_genetic_models_header(head)
                add_model_score_header(head)
                add_compounds_header(head)
                
                vcf_individuals = head.individuals
                logger.info("Individuals found in vcf file: {}".format(
                    ', '.join(vcf_individuals)))
                print_headers(head, outfile)
                headers_done = True
                check_individuals(family_parser, vcf_individuals)
            else:
                break
                
                
    #
    # ###################################################################
    # ### The task queue is where all jobs(in this case batches that  ###
    # ### represents variants in a region) is put. The consumers will ###
    # ### then pick their jobs from this queue.                       ###
    # ###################################################################
    #
    # logger.debug("Setting up a JoinableQueue for storing variant batches")
    # variant_queue = JoinableQueue(maxsize=1000)
    # logger.debug("Setting up a Queue for storing results from workers")
    # results = Manager().Queue()
    #
    #
    # num_model_checkers = processes
    # #Adapt the number of processes to the machine that run the analysis
    # logger.info('Number of CPU:s {}'.format(cpu_count()))
    # logger.info('Number of model checkers: {}'.format(num_model_checkers))
    #
    # # We use a temp file to store the processed variants
    # logger.debug("Build a tempfile for printing the variants")
    # temp_file = NamedTemporaryFile(delete=False)
    # temp_file.close()
    # # Open the temp file with codecs
    # temporary_variant_file = open(
    #                             temp_file.name,
    #                             mode='w',
    #                             encoding='utf-8',
    #                             errors='replace'
    #                             )
    #
    #
    # # These are the workers that do the heavy part of the analysis
    # logger.info('Seting up the workers')
    # model_checkers = [
    #     VariantAnnotator(
    #         variant_queue,
    #         results,
    #         **annotator_arguments
    #     )
    #     for i in range(num_model_checkers)
    # ]
    

if __name__ == '__main__':
    annotate_models()