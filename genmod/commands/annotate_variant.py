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
import os
import logging
import pkg_resources
import itertools

import click
import tabix

from multiprocessing import JoinableQueue, Manager, cpu_count
from codecs import open
from tempfile import NamedTemporaryFile
from datetime import datetime

from genmod import __version__
# from genmod.annotate_regions import load_annotations, check_overlap
from genmod.vcf_tools import (HeaderParser, add_vcf_info, add_metadata, 
print_headers, print_variant, sort_variants)

from genmod.annotate_variants import VariantAnnotator
from genmod.annotate_regions import (get_genes, check_exonic, load_annotations)
from genmod.utils import VariantPrinter

from .utils import (outfile, silent, processes, temp_dir, variant_file)

logger = logging.getLogger(__name__)

@click.command()
@variant_file
@click.option('-r', '--annotate_regions', 
                is_flag=True,
                help='Increase output verbosity.'
)
@click.option('-c', '--cadd_file', 
                    multiple = True,
                    type=click.Path(exists=True), 
                    help="Specify the path to a bgzipped cadd file"\
                    " (with index) with variant scores. This command can be"\
                    " used multiple times if multiple cadd files."
)
@click.option('--thousand_g',
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
@click.option('--max_af',
                    is_flag=True,
                    help="If the MAX AF should be annotated"
)
@click.option('--spidex',
                    type=click.Path(exists=True), 
                    help="Specify the path to a bgzipped tsv file"\
                            " (with index) with spidex information."
)
@click.option('-a' ,'--annotation_dir',
                    type=click.Path(exists=True),
                    default=pkg_resources.resource_filename('genmod', 'annotations'),
                    help="""Specify the path to the directory where the annotation 
                    databases are. 
                    Default is the gene pred files that comes with the distribution."""
)
@click.option('--cadd_raw', 
                    is_flag=True,
                    help="""If the raw cadd scores should be annotated."""
)
@processes
@outfile
@silent
@temp_dir
def annotate(variant_file, annotate_regions, cadd_file, thousand_g, exac, 
spidex,annotation_dir, outfile, silent, cadd_raw, cosmic, max_af, processes,
temp_dir):
    """
    Annotate vcf variants.
    
    Annotate variants with a number of different sources.
    Please use --help for more info.
    """

    logger.info("Running genmod annotate_variant version {0}".format(__version__))
    
    start_time_analysis = datetime.now()
    annotator_arguments = {}
    
    logger.info("Initializing a Header Parser")
    head = HeaderParser()
    
    line = None
    for line in variant_file:
        line = line.rstrip()

        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            break
    
    #Add the first variant to the iterator
    if line:
        variant_file = itertools.chain([line], variant_file)
    
    header_line = head.header
    annotator_arguments['header_line'] = header_line
    
    if annotate_regions:
        logger.info("Loading annotations")
        gene_trees, exon_trees = load_annotations(annotation_dir)
        annotator_arguments['gene_trees'] = gene_trees
        annotator_arguments['exon_trees'] = exon_trees
        
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

    if exac:
        logger.info("Annotating ExAC frequencies")
        logger.debug("Using ExAC file: {0}".format(exac))
        annotator_arguments['exac'] = exac
        add_metadata(
            head,
            'info',
            'ExACAF',
            annotation_number='1',
            entry_type='Float',
            description="Frequency in the ExAC database."
        )
        
    if thousand_g:
        logger.info("Annotating 1000G frequencies")
        logger.debug("Using 1000G file: {0}".format(thousand_g))
        annotator_arguments['thousand_g'] = thousand_g
        add_metadata(
            head,
            'info',
            '1000GAF',
            annotation_number='1',
            entry_type='Float',
            description="Frequency in the 1000G database."
        )

    if spidex:
        logger.info("Annotating Spidex z scores")
        logger.debug("Using Spidex file: {0}".format(spidex))
        annotator_arguments['spidex'] = spidex
        add_metadata(
            head,
            'info',
            'SPIDEX',
            annotation_number='1',
            entry_type='Float',
            description="Z score from the spidex database."
        )
    
    if cadd_file:
        logger.info("Annotating CADD scores")
        logger.debug("Using CADD file(s): {0}".format(', '.join(cadd_file)))
        annotator_arguments['cadd_files'] = cadd_file
        any_cadd_file = True

        add_metadata(
            head,
            'info',
            'CADD',
            annotation_number='1',
            entry_type='Integer',
            description="The CADD relative score for this alternative."
        )
        if cadd_raw:
            annotator_arguments['cadd_raw'] = cadd_raw
            logger.debug("Adding vcf metadata for CADD raw score")
            add_metadata(
                head,
                'info',
                'CADD_raw',
                annotation_number='1',
                entry_type='Float',
                description="The CADD raw score(s) for this alternative(s)."
            )

    if max_af:
        annotator_arguments['max_af'] = max_af
        if thousand_g:
            add_metadata(
                head,
                'info',
                '1000G_MAX_AF',
                annotation_number='1',
                entry_type='Float',
                description="The max af for thousand genomes populations."
            )
        if exac:
            add_metadata(
                head,
                'info',
                'ExAC_MAX_AF',
                annotation_number='1',
                entry_type='Float',
                description="The max af for ExAC populations."
            )

    if cosmic:
        logger.info("Annotating if variant is in COSMIC")
        logger.debug("Using COSMOC file: {0}".format(cosmic))
        annotator_arguments['cosmic'] = cosmic
        add_metadata(
            head,
            'info',
            'COSMIC',
            annotation_number='0',
            entry_type='Flag',
            description="If variant is in COSMIC database."
        )
    
    ###################################################################
    ### The task queue is where all jobs(in this case batches that  ###
    ### represents variants in a region) is put. The consumers will ###
    ### then pick their jobs from this queue.                       ###
    ###################################################################

    logger.debug("Setting up a JoinableQueue for storing variant batches")
    variant_queue = JoinableQueue(maxsize=1000)
    logger.debug("Setting up a Queue for storing results from workers")
    results = Manager().Queue()

    num_annotators = processes
    #Adapt the number of processes to the machine that run the analysis
    if cadd_file or spidex:
        # We need more power when annotating cadd scores:
        # But if flag is used that overrides
        if num_annotators == min(4, cpu_count()):
            num_annotators = min(8, cpu_count())

    logger.info('Number of CPU:s {}'.format(cpu_count()))
    logger.info('Number of model checkers: {}'.format(num_annotators))


    # These are the workers that do the heavy part of the analysis
    logger.info('Setting up the workers')
    annotators = [
        VariantAnnotator(
            variant_queue, 
            results, 
            **annotator_arguments
        )
        for i in range(num_annotators)
    ]

    logger.info('Starting the workers')
    for worker in annotators:
        logger.debug('Starting worker {0}'.format(worker))
        worker.start()

    # This process prints the variants to temporary files
    # If there is only one annotation process we can print the results as soon
    # as they are done
    logger.info('Setting up the variant printer')
    if len(annotators) == 1:
        print_headers(head, outfile, silent)
        var_printer = VariantPrinter(
                        task_queue = results, 
                        head = head, 
                        mode='normal', 
                        outfile = outfile
                        )
    else:
        # We use a temp file to store the processed variants
        logger.debug("Build a tempfile for printing the variants")
        if temp_dir:
            temp_file = NamedTemporaryFile(delete=False, dir=temp_dir)
        else:
            temp_file = NamedTemporaryFile(delete=False)
            
        temp_file.close()
        
        var_printer = VariantPrinter(
                        task_queue = results, 
                        head = head, 
                        mode='chromosome', 
                        outfile = temp_file.name
                        )
    
    logger.info('Starting the variant printer process')
    var_printer.start()

    start_time_variant_parsing = datetime.now()
    start_time_twenty = datetime.now()
    nr_of_lines = 0
    # This process parses the original vcf and create batches to put in the variant queue:
    logger.info('Start parsing the variants')
    
    for line in variant_file:
        line = line.rstrip()
        
        if not line.startswith('#'):
            variant_queue.put(line)
            
            nr_of_lines += 1
            
            if nr_of_lines % 20000 == 0:
                logger.info('{0} variants parsed'.format(nr_of_lines))
                logger.info('Last 20000 took {0} to parse'.format(
                    datetime.now()-start_time_twenty))
                start_time_twenty = datetime.now()
    
    logger.info('Put stop signs in the variant queue')
    
    for i in range(num_annotators):
        variant_queue.put(None)

    variant_queue.join()
    results.put(None)
    var_printer.join()

    if len(annotators) > 1:
        logger.info("Start sorting the variants")
        sort_variants(temp_file.name, mode='chromosome')

        logger.info("Print the headers")
        print_headers(head, outfile, silent)

        with open(temp_file.name, 'r', encoding='utf-8') as f:
            for line in f:
                print_variant(
                    variant_line=line,
                    outfile=outfile,
                    mode='modified',
                    silent=silent
                )

        logger.info("Removing temp file")
        os.remove(temp_file.name)
        logger.debug("Temp file removed")

    logger.info('Time for whole analyis: {0}'.format(
        str(datetime.now() - start_time_analysis)))
    

if __name__ == '__main__':
    from genmod.log import init_log
    from genmod import logger
    init_log(logger, loglevel="INFO")
    annotate()
