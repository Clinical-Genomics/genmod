#!/usr/bin/env python
# encoding: utf-8
"""
score_compounds.py

Script for scoring compound variants in VCF files.

Created by Måns Magnusson on 2015-09-03.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import sys
import os
import click
import logging
import itertools

from multiprocessing import JoinableQueue, Manager, cpu_count, util
from codecs import open
from datetime import datetime
from tempfile import NamedTemporaryFile

from genmod.vcf_tools import (HeaderParser, add_metadata, print_headers, 
sort_variants, print_variant)
from genmod.utils import (get_batches, VariantPrinter)
from genmod.score_variants import CompoundScorer

from genmod import __version__

from .utils import (variant_file, silent, outfile, processes, temp_dir, 
                    get_file_handle)

logger = logging.getLogger(__name__)
util.abstract_sockets_supported = False

@click.command('compound', short_help="Score compounds")
@variant_file
@silent
@outfile
@processes
@temp_dir
@click.option('--vep', 
                    is_flag=True,
                    help='If variants are annotated with the Variant Effect Predictor.'
)
@click.option('--threshold', type=int, help="Threshold for model-dependent penalty if no compounds with passing score", default=9)
@click.option('--penalty', type=int, help="Penalty applied together with --threshold", default=6)
@click.pass_context
def compound(context, variant_file, silent, outfile, vep, threshold: int, penalty: int, processes, temp_dir):
    """
    Score compound variants in a vcf file based on their rank score.
    """
    logger.info('Running GENMOD score_compounds, version: {0}'.format(__version__))
    
    variant_file = get_file_handle(variant_file)
    
    start_time_analysis = datetime.now()
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

    logger.info("Headers parsed")
    
    if not line.startswith('#'):
        variant_file = itertools.chain([line], variant_file)
    else:
        print_headers(head=head, outfile=outfile, silent=silent)
        sys.exit(0)

    header_line = head.header
    individuals = head.individuals

    add_metadata(head,
                 'info',
                 'CompoundsNormalized',
                 annotation_number='.',
                 entry_type='String',
                 description='Rank score as provided by compound analysis, based on RankScoreNormalized. family_id:rank_score')

    ###################################################################
    ### The task queue is where all jobs(in this case batches that  ###
    ### represents variants in a region) is put. The consumers will ###
    ### then pick their jobs from this queue.                       ###
    ###################################################################

    logger.debug("Setting up a JoinableQueue for storing variant batches")
    variant_queue = JoinableQueue(maxsize=1000)
    logger.debug("Setting up a Queue for storing results from workers")
    results = Manager().Queue()

    num_scorers = processes
    #Adapt the number of processes to the machine that run the analysis
    logger.info('Number of CPU:s {}'.format(cpu_count()))
    logger.info('Number of model checkers: {}'.format(num_scorers))

    # These are the workers that do the heavy part of the analysis
    logger.info('Seting up the workers')
    compound_scorers = [
        CompoundScorer(
            task_queue=variant_queue,
            results_queue=results,
            individuals=individuals,
            threshold=threshold,
            penalty=penalty,
        )
        for i in range(num_scorers)
    ]
    
    try:
        logger.info('Starting the workers')
        for worker in compound_scorers:
            logger.debug('Starting worker {0}'.format(worker))
            worker.start()
        
        # This process prints the variants to temporary files
        logger.info('Seting up the variant printer')
        
        # We use a temp file to store the processed variants
        logger.debug("Build a tempfile for printing the variants")
        if temp_dir:
            temp_file = NamedTemporaryFile(delete=False, dir=temp_dir)
        else:
            temp_file = NamedTemporaryFile(delete=False)
        temp_file.close()
        
        variant_printer = VariantPrinter(
            task_queue=results,
            head=head,
            mode='chromosome',
            outfile = temp_file.name
        )
        
        logger.info('Starting the variant printer process')
        variant_printer.start()
        
        start_time_variant_parsing = datetime.now()
        
        # This process parses the original vcf and create batches to put in the variant queue:
        chromosome_list = get_batches(
                                    variants = variant_file,
                                    batch_queue = variant_queue,
                                    header = head,
                                    vep = vep,
                                    results_queue=results
                                )
        
        logger.debug("Put stop signs in the variant queue")
        for i in range(num_scorers):
            variant_queue.put(None)
        
        variant_queue.join()
        results.put(None)
        variant_printer.join()
        
        sort_variants(infile=temp_file.name, mode='chromosome')
        
        print_headers(head=head, outfile=outfile, silent=silent)
        
        with open(temp_file.name, 'r', encoding='utf-8') as f:
            for line in f:
                print_variant(
                    variant_line=line,
                    outfile=outfile,
                    mode='modified',
                    silent=silent
                )
    except Exception as e:
        logger.warning(e)
        for worker in compound_scorers:
            worker.terminate()
        variant_printer.terminate()
        context.abort()
    finally:
        logger.info("Removing temp file")
        os.remove(temp_file.name)
        logger.debug("Temp file removed")
    

    logger.info('Time for whole analyis: {0}'.format(
        str(datetime.now() - start_time_analysis)))
