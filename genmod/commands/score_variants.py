#!/usr/bin/env python
# encoding: utf-8
"""
score_variants.py

Script for scoring genetic variants in VCF files.

Created by Måns Magnusson on 2015-09-03.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import sys
import os
import click
import logging
<<<<<<< HEAD
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


from genmod import __version__ as VERSION


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
    logger = logging.getLogger(__name__)
    ## Collect supplied plugin
    alt_dict, score_dict, value_dict, operation_dict = collectKeys(
                                                            config_file,
                                                            variant_parser,
                                                            verbose
                                                        )
    logging.debug("Plugin file: " + config_file)
    logging.debug("alt_dict:" + str(alt_dict))
    logging.debug("score_dict: " + str(score_dict))
    logging.debug("value_dict: " + str(value_dict))
    logging.debug("operation_dict" + str(operation_dict))
    
    return alt_dict, score_dict, value_dict, operation_dict
=======


from codecs import open
from datetime import datetime
from validate import ValidateError
>>>>>>> feature/fix_compounds_single

from genmod.vcf_tools import (add_metadata, print_variant_dict, add_vcf_info,
print_headers, HeaderParser, get_variant_dict, get_info_dict)
from genmod.score_variants import (ConfigParser, score_variant, check_plugins)

from genmod import __version__

@click.command()
@click.argument('variant_file',
                nargs=1,
                type=click.File('r'),
                metavar='<vcf_file> or -'
)
@click.option('-f', '--family_id',
                default='1', 
)
@click.option('-s', '--silent',
                is_flag=True,
                help='Do not print the variants.'
)
@click.option('-o', '--outfile',
                type=click.File('w'),
                help='Specify the path to a file where results should be stored.'
)
@click.option('-c', '--score_config',
              type=click.Path(exists=True),
              help="The plug-in config file(.ini)"
)
@click.option('-v', '--verbose',
              count=True,
              help='Increase output verbosity. If -vv all scores will be printed'
)
def score(variant_file, family_id, score_config, silent, outfile, verbose):
    """
    Score variants in a vcf file using Weighted Sum Model.
<<<<<<< HEAD
    The specific scores should be defined in a config file, see examples in 
    genmod/configs
    """
    from genmod.log import init_log, LEVELS
    from genmod import logger as root_logger
    loglevel = LEVELS.get(min(verbose,2), "WARNING")
    init_log(root_logger, loglevel=loglevel)
    
    logger = logging.getLogger(__name__)
    
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [i+'='+str(values[i]) for i in values if values[i] and 
                            i != 'args' and i != 'frame' and i != 'parser']
    
    start_time_analysis = datetime.now()
    
    
    logger.info('Running GENMOD score, version: {0}'.format(VERSION))
    
    ## Start by parsing the pedigree file:
    prefered_models = []
    family_id = None
    
    if family_file:
        prefered_models, family_id = get_genetic_models(
                                                    family_file, 
                                                    family_type
                                                    )
    else:
        logger.critical("Please provide a family file")
        sys.exit()
    
    logger.info('Prefered model found in family file:{0}'.format(prefered_models))
    
    if not plugin_file:
        logger.critical("Please provide a plugin file")
        sys.exit()
=======
    
    The specific scores should be defined in a config file, see examples on 
    github.
    """
    # from genmod import logger as root_logger
    # from genmod.log import init_log, LEVELS
    # loglevel = LEVELS.get(min(verbose,2), "WARNING")
    # init_log(root_logger, loglevel=loglevel)
    
    logger = logging.getLogger(__name__)
    # logger = logging.getLogger("genmod.commands.score")
    
    logger.info('Running GENMOD score, version: {0}'.format(__version__))
>>>>>>> feature/fix_compounds_single
    
    ## Check the score config:
    if not score_config:
        logger.warning("Please provide a score config file.")
        sys.exit(1)
    
    logger.debug("Parsing config file")
    
<<<<<<< HEAD
    # If the variants are already annotated we do not need to redo the annotation
    if not vep:
        gene_trees, exon_trees = load_annotations(annotation_dir, verbose)
    else:
        logger.info('Using VEP annotation')
=======
    try:
        config_parser = ConfigParser(score_config)
    except ValidateError as e:
        logger.error("Something wrong in plugin file, please see log")
        logger.info("Exiting")
        sys.exit(1)

    logger.debug("Config parsed succesfully")

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
>>>>>>> feature/fix_compounds_single
    
    variant_file.seek(0)
    header_line = head.header
    
<<<<<<< HEAD
    if variant_file == '-':
        variant_parser = VCFParser(
            fsock = sys.stdin, 
            )
    else:
        variant_parser = VCFParser(
            infile = variant_file, 
            )
=======
    if "RankScore" in head.info_dict:
        logger.warning("Variants already scored according to VCF header")
        logger.info("Please check VCF file")
        logger.info("Exiting...")
        sys.exit(1)
>>>>>>> feature/fix_compounds_single
    
    
    add_metadata(
        head,
        'info',
        'IndividualRankScore',
        annotation_number='.', 
        entry_type='String', 
        description="Individual rank score for the variant in this family. "\
        "This score is NOT corrected for compounds"
    )
    
    add_metadata(
        head,
        'info',
        'RankScore',
<<<<<<< HEAD
        annotation_number='.', 
        entry_type='String', 
        description="Combined rank score for the variant in this family. "\
        "This score is corrected for compounds"
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
    
    logger.info('Number of CPU:s {0}'.format(cpu_count()))
    logger.info('Number of model scorers: {0}'.format(num_model_scorers))
    
    
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
    
    logger.info('Start parsing the variants ... ')
    
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
    
    logger.info('Cromosomes found in variant file: {0}'.format(','.join(chromosome_list)))
    logger.info('Variants scored!')
    
    
    sort_variants(
            infile=temp_file.name,
            mode='rank',
            verbose=verbose
    )
    
    print_headers(head, outfile, silent)
    
    print_variants(temp_file.name, outfile, mode='modified',  silent=silent)
    
    os.remove(temp_file.name)
    
    logger.info('Time for whole analyis: {0}'.format(str(datetime.now() - start_time_analysis)))
    
    
=======
        annotation_number='1',
        entry_type='Integer',
        description="The rank score for this variant."
    )
    print_headers(
        head=head,
        outfile=outfile,
        silent=silent
    )
    start_scoring = datetime.now()
    last_twenty = datetime.now()
    nr_of_variants = 1

    for line in variant_file:
        if not line.startswith('#'):
            variant = get_variant_dict(line, header_line)
            variant['info_dict'] = get_info_dict(variant['INFO'])

            rank_score = score_variant(variant, config_parser)
            
            variant = add_vcf_info(
                keyword = 'RankScore',
                variant_dict=variant,
                annotation="{0}:{1}".format(family_id, rank_score))

            print_variant_dict(
                variant=variant,
                header_line=header_line,
                outfile=outfile,
                silent=silent)

            nr_of_variants += 1

            if nr_of_variants % 20000 == 0:
                logger.info("{0} variants scored.".format(nr_of_variants))
                logger.info("Last 20000 took {0} to score.".format(datetime.now()-last_twenty))
                last_twenty = datetime.now()

    logger.info("Variants scored. Number of variants: {0}".format(nr_of_variants))
    logger.info("Time to score variants: {0}".format(datetime.now()-start_scoring))
        

>>>>>>> feature/fix_compounds_single
if __name__ == '__main__':
    score()
