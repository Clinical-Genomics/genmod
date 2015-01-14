#!/usr/bin/env python
# encoding: utf-8
"""
score.py

Script for scoring genetic variants in VCF files.

Created by Henrik Stranneheim and Måns Magnusson on 2015-01-08.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import argparse
import inspect
import pkg_resources
import click
import genmod

from multiprocessing import JoinableQueue
from codecs import open, getwriter
from tempfile import NamedTemporaryFile
from datetime import datetime

from pprint import pprint as pp

from genmod import (VariantScorer, FileSort, get_batches)
from genmod.utils import (collectKeys, load_annotations)

from ped_parser import parser as ped_parser
from vcf_parser import parser as vcf_parser

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
    if verbose:
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
        family_file (file): A file with family information in ped or ped like format
    
    Returns:
        inheritance_models  (set): A list with the expected inheritance models
    """
    inheritance_models = set([])
    my_family_parser = ped_parser.FamilyParser(family_file, family_type)
    for family in my_family_parser.families:
        for model in my_family_parser.families[family].models_of_inheritance:
            if model not in ['NA', 'na', 'Na']:
                inheritance_models.add(model)
    # Stupid thing but for now when we only look at one family
    return inheritance_models


def add_metadata(head, command_line_string):
    """Add metadata for the information added by this script."""

    head.add_info('RankScore', '1', 'Integer',
                  "Combined rank score for the variant in this family.")
    head.add_version_tracking('score_mip_variants', VERSION,
                              str(datetime.now()), command_line_string)
    
    return


def print_headers(head, outfile=None, silent=False):
    """Print the headers to a results file."""
    if outfile:
        with open(outfile, 'w', encoding='utf-8') as f:
            for head_count in head.print_header():
                f.write(head_count+'\n')
    else:
        if not silent:
            for line in head.print_header():
                print(line)
    return


def print_version(ctx, param, value):
    """Callback function for printing version and exiting.

    Args:
        ctx   (object) : Current context
        param (object) : Click parameter(s)
        value (boolean) : Click parameter was supplied or not

    Returns:
        None:
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo('score_mip_variants version: ' + VERSION)
    ctx.exit()


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
              default=pkg_resources.resource_filename('genmod', 'configs/rank_model_test_v1.0.ini'),
              help="The plug-in config file(.ini)"
)
@click.option('-v', '--verbose',
              is_flag=True,
              help='Increase output verbosity.'
)
def score(family_file, variant_file, family_type, annotation_dir, vep,
                       plugin_file, verbose, silent, outfile):
    """
    Score variants in a vcf file using Weighted Sum Model.
    The specific scores should be defined in a config file, see examples in 
    genmod/configs
    """
    
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [i+'='+str(values[i]) for i in values if values[i] and i != 'args' and i != 'frame' and i != 'parser']
    
    start_time_analysis = datetime.now()
    
    if verbose:
        print('\nRunning GENMOD score, version: %s \n' % VERSION, file=sys.stderr)
    
    ## Start by parsing the pedigree file:
    prefered_models = []
    if family_file:
        prefered_models = get_genetic_models(family_file, family_type)
    
    if verbose:
        print('Prefered model found in family file: %s \n' % prefered_models, file=sys.stderr)
    
    ######### Read to the annotation data structures #########
    
    gene_trees = {}
    exon_trees = {}
    
    # If the variants are already annotated we do not need to redo the annotation
    if not vep:
        gene_trees, exon_trees = load_annotations(annotation_dir, verbose)
    else:
        if verbose:
            print('Using VEP annotation', file=sys.stderr)
    
    ## Check the variants:
    
    if variant_file == '-':
        variant_parser = vcf_parser.VCFParser(fsock=sys.stdin)
    else:
        variant_parser = vcf_parser.VCFParser(infile=variant_file)
    
    head = variant_parser.metadata
    
    alt_dict, score_dict, value_dict, operation_dict = check_plugin(plugin_file, variant_parser, verbose)
    
    
    variant_queue = JoinableQueue(maxsize=1000)
    temp_file = NamedTemporaryFile(delete=False)
    temp_file.close()
    
    try:
        temporary_variant_file = open(temp_file.name, mode='w', encoding='utf-8', errors='replace')
        
        scorer = VariantScorer(
                                variant_queue,
                                temporary_variant_file,
                                variant_parser.header,
                                prefered_models,
                                alt_dict, 
                                score_dict, 
                                value_dict,
                                operation_dict, 
                                verbose
                            )
        
        scorer.start()
        
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
        
        variant_queue.put(None)
        
        variant_queue.join()
        
        temporary_variant_file.seek(0)
        
        ## Add the new metadata to the headers:
        add_metadata(head, ','.join(argument_list))
        print_headers(head, outfile, silent)
        
        var_sorter = FileSort(
                            infile=temporary_variant_file,
                            mode='rank',
                            outfile=outfile, 
                            silent=silent
                        )
        var_sorter.sort()
        
    finally:
        
        temporary_variant_file.close()
        os.remove(temp_file.name)
    
if __name__ == '__main__':
    score()
