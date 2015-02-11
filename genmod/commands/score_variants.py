#!/usr/bin/env python
# encoding: utf-8
"""
score_variants.py

Script for scoring genetic variants in VCF files.

Created by Henrik Stranneheim and Måns Magnusson on 2015-01-08.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os
import argparse
import inspect
import pkg_resources
import click
import genmod

from multiprocessing import JoinableQueue, Manager, cpu_count
from codecs import open, getwriter
from tempfile import NamedTemporaryFile
from datetime import datetime

from pprint import pprint as pp

from genmod import (VariantScorer, VariantPrinter, get_batches, collectKeys, 
                    load_annotations, add_metadata, warning, print_headers,
                    sort_variants, print_variants)

from ped_parser import FamilyParser
from vcf_parser import VCFParser

# Import third party library
# https://github.com/mitsuhiko/logbook
from logbook import Logger, StderrHandler
log = Logger('Logbook')
log_handler = StderrHandler()

if sys.version_info < (3,0):
    sys.stdout = getwriter('UTF-8')(sys.stdout)

VERSION = pkg_resources.require("genmod")[0].version


def check_plugin(config_file, variant_parser, verbose=False):
    """
    Collect keys from config_file and check vcf compatibility.
    
    Args:
        config_file   (file) : plugin file
        variant_parser (object) : vcf_parser object
        verbose       (boolean, optional) : Enable informative print
    
    Returns:
        dict:   Dictionary of alternatives
    """
    ## Collect supplied plugin
    alt_dict, score_dict, value_dict, operation_dict = collectKeys(
                                                            config_file,
                                                            variant_parser,
                                                            verbose
                                                        )
    if verbose == 2:
        log.info("Plugin file: " + config_file)
        log.info("alt_dict:" + str(alt_dict))
        log.info("score_dict: " + str(score_dict))
        log.info("value_dict: " + str(value_dict))
        log.info("operation_dict" + str(operation_dict))
    
    return alt_dict, score_dict, value_dict, operation_dict

def get_genetic_models(family_file, family_type):
    """
    Return the genetic models found for the family(families).
    
    Args:
        family_file (file): A file with family information 
                            in ped or ped like format.
    
    Returns:
        inheritance_models  : A set with the expected inheritance models
        family_id   : A string that represents the family id
    """
    inheritance_models = set([])
    my_family_parser = FamilyParser(family_file, family_type)
    family_id = None
    for family in my_family_parser.families:
        family_id = family
        for model in my_family_parser.families[family].models_of_inheritance:
            if model not in ['NA', 'na', 'Na']:
                inheritance_models.add(model)
    # Stupid thing but for now when we only look at one family
    return inheritance_models, family_id


@click.command()
@click.argument('variant_file',
                nargs=1,
                type=click.Path(),
                metavar='<vcf_file> or -'
)
@click.option('-f', '--family_file',
                nargs=1, 
                type=click.Path(exists=True),
                metavar='<ped_file>'
)
@click.option('-t' ,'--family_type', 
                type=click.Choice(['ped', 'alt', 'cmms', 'mip']), 
                default='ped',
                help='If the analysis use one of the known setups, please specify which one.'
)
@click.option('-a' ,'--annotation_dir', 
                    type=click.Path(exists=True),
                    default=pkg_resources.resource_filename('genmod', 'annotations'),
                    help="""Specify the path to the directory where the annotation 
                    databases are. 
                    Default is the gene pred files that comes with the distribution."""
)
@click.option('--vep', 
                    is_flag=True,
                    help='If variants are annotated with the Variant Effect Predictor.'
)
@click.option('-s', '--silent',
                is_flag=True,
                help='Do not print the variants.'
)
@click.option('-o', '--outfile',
                type=click.Path(exists=False),
                help='Specify the path to a file where results should be stored.'
)
@click.option('-pi', '--plugin_file',
              type=click.Path(exists=True),
              default=pkg_resources.resource_filename('genmod', 'configs/rank_model_test_v1.4.ini'),
              help="The plug-in config file(.ini)"
)
@click.option('-p', '--processes', 
                default=min(4, cpu_count()),
                help='Define how many processes that should be use for annotation.'
)
@click.option('-v', '--verbose',
              count=True,
              help='Increase output verbosity. If -vv all scores will be printed'
)

def score(family_file, variant_file, family_type, annotation_dir, vep,
                       plugin_file, processes, silent, outfile, verbose):
    """
    Score variants in a vcf file using Weighted Sum Model.
    The specific scores should be defined in a config file, see examples in 
    genmod/configs
    """
    
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [i+'='+str(values[i]) for i in values if values[i] and 
                            i != 'args' and i != 'frame' and i != 'parser']
    
    start_time_analysis = datetime.now()
    
    
    if verbose:
        log.info('Running GENMOD score, version: %s \n' % VERSION)
    
    ## Start by parsing the pedigree file:
    prefered_models = []
    family_id = None
    
    if family_file:
        prefered_models, family_id = get_genetic_models(
                                                    family_file, 
                                                    family_type
                                                    )
    
    if verbose:
        log.info('Prefered model found in family file: %s \n' % 
                prefered_models)
    
    ######### Read to the annotation data structures #########
    
    gene_trees = {}
    exon_trees = {}
    
    # If the variants are already annotated we do not need to redo the annotation
    if not vep:
        gene_trees, exon_trees = load_annotations(annotation_dir, verbose)
    else:
        if verbose:
            log.info('Using VEP annotation')
    
    ## Check the variants:
    
    if variant_file == '-':
        variant_parser = VCFParser(fsock=sys.stdin)
    else:
        variant_parser = VCFParser(infile=variant_file)
    
    head = variant_parser.metadata
    
    add_metadata(
        head,
        'version',    
        'genmod_score', 
        version=VERSION, 
        command_line_string=' '.join(argument_list)
    )
    
    add_metadata(
        head,
        'info',
        'RankScore',
        annotation_number='1', 
        entry_type='Integer', 
        description="Combined rank score for the variant in this family."'GeneticModels'
    )
    
    alt_dict, score_dict, value_dict, operation_dict = check_plugin(
                                                            plugin_file, 
                                                            variant_parser, 
                                                            verbose
                                                            )
    
    ####################################################################
    ### The variant queue is where all jobs(in this case batches that###
    ### represents variants in a region) is put. The consumers will  ###
    ### then pick their jobs from this queue.                        ###
    ####################################################################
    
    variant_queue = JoinableQueue(maxsize=1000)
    # The consumers will put their results in the results queue
    results = Manager().Queue()
    
    
    num_model_scorers = processes
    
    if verbose:
        log.info('Number of CPU:s %s' % cpu_count(), file=sys.stderr)
        log.info('Number of model scorers: %s' % num_model_scorers, file=sys.stderr)
    
    
    temp_file = NamedTemporaryFile(delete=False)
    temp_file.close()
    
    # We open a variant file to print the variants before sorting:
    temporary_variant_file = open(
                                temp_file.name, 
                                mode='w', 
                                encoding='utf-8', 
                                errors='replace'
                                )
    
    model_scorers = [
                        VariantScorer(
                            variant_queue,
                            results,
                            variant_parser.header,
                            prefered_models,
                            family_id,
                            alt_dict, 
                            score_dict, 
                            value_dict,
                            operation_dict, 
                            verbose
                        )
                        for i in range(num_model_scorers)
                    ]
    
    
    for proc in model_scorers:
        proc.start()
    
    # This process prints the variants to temporary files
    var_printer = VariantPrinter(
                        results,
                        temporary_variant_file,
                        head,
                        mode='score',
                        verbosity=verbose
                    )
    
    var_printer.start()
    
    start_time_variant_parsing = datetime.now()
    
    if verbose:
        log.info('Start parsing the variants ... \n')
    
    # get_batches put the variants in the queue and returns all chromosomes
    # found among the variants
    chromosome_list = get_batches(
                        variant_parser, 
                        variant_queue,
                        individuals = [],
                        gene_trees = gene_trees, 
                        exon_trees = exon_trees, 
                        phased = False, 
                        vep = vep, 
                        whole_genes = True, 
                        verbosity = verbose
                    )
    
    
    # Put stop signs in the variant queue
    for i in range(num_model_scorers):
        variant_queue.put(None)
    
    variant_queue.join()
    
    results.put(None)
    var_printer.join()
    
    temporary_variant_file.close()
    
    if verbose:
        log.info('Cromosomes found in variant file: %s \n' % ','.join(chromosome_list))
        log.info('Variants scored!\n')
    
    
    sort_variants(
            infile=temp_file.name,
            mode='rank',
            verbose=verbose
    )
    
    print_headers(head, outfile, silent)
    
    print_variants(temp_file.name, outfile, mode='modified',  silent=silent)
    
    os.remove(temp_file.name)
    
    if verbose:
        log.info('Time for whole analyis: %s' % str(datetime.now() - start_time_analysis))
    
    
if __name__ == '__main__':
    score()
