#!/usr/bin/env python
# encoding: utf-8
"""
annotate.py

Command line tool for annotating variants with genmod.

Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import (print_function, unicode_literals)

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
from vcf_parser import VCFParser

from genmod import (__version__, VariantPrinter)

from genmod.variant_annotation import (VariantAnnotator)

# , annotation_parser,
#                     , get_batches, sort_variants,
#                     load_annotations, print_headers, print_variants,
#                     add_metadata, warning, __version__)

from genmod.utils import (get_batches, load_annotations, add_metadata, 
                          sort_variants, print_variants, print_headers)

from genmod.log import init_log


def check_tabix_index(compressed_file, file_type='cadd', verbose=False):
    """
    Check if a compressed file have a tabix index, if not build one.
    
    Args:
        compressed_file : Path to a file that is assumed to be compressed.
        file_type   : The type of the file.
        verobe      : Increase output verbosity
    
    Returns:
        0 if everything went ok.
    
    """
    if file_type == 'cadd':
        try:
            tabix_index(compressed_file, seq_col=0, start_col=1, end_col=1, meta_char='#')
        except IOError as e:
            pass
    elif file_type == 'vcf':
        try:
            tabix_index(compressed_file, preset='vcf')
        except IOError as e:
            pass
    return 0



###           This is for annotating the variants       ###


@click.command()
@click.argument('variant_file', 
                    nargs=1, 
                    type=click.Path(),
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
@click.option('--silent', 
                    is_flag=True,
                    help='Do not print the variants.'
)
@click.option('-split' ,'--split_variants', 
                    is_flag=True,
                    help='If the variants should be splitted.'
)
@click.option('-g' ,'--whole_gene', 
                    is_flag=True,
                    help="""If compounds should be checked in the whole gene regions. 
                    Not only exonic/splice sites."""
)
@click.option('-r' ,'--cadd_raw', 
                    is_flag=True,
                    help="""If the raw cadd scores should be annotated."""
)
@click.option('-a' ,'--annotation_dir', 
                    type=click.Path(exists=True),
                    default=pkg_resources.resource_filename('genmod', 'annotations'),
                    help="""Specify the path to the directory where the annotation 
                    databases are. 
                    Default is the gene pred files that comes with the distribution."""
)
@click.option('-o', '--outfile', 
                    type=click.Path(exists=False),
                    help='Specify the path to a file where results should be stored.'
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
@click.option('--dbnfsp',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped dbNSFP file (with index). THIS IS CURRENTLY NOT DOING ANYTHING!"""
)
@click.option('-p', '--processes', 
                default=min(4, cpu_count()),
                help='Define how many processes that should be use for annotation.'
)

@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def annotate(family_file, variant_file, family_type, vep, silent, phased, strict, cadd_raw, whole_gene, 
                annotation_dir, cadd_file, cadd_1000g, cadd_exac, cadd_esp, cadd_indels, thousand_g, exac, outfile,
                split_variants, processes, dbnfsp, verbose):
    """
    Annotate variants in a VCF file.\n
    
    The main function with genmod is to annotate genetic inheritance 
    patterns for variants in families. 
    Use flag --family together with a .ped file to describe which individuals
     in the vcf you wish to check inheritance patterns for.
    Individuals that are not present in the ped file will not be considered 
    in the analysis.\n
    It is also possible to use genmod without a family file. In this case the 
    variants will be annotated with a variety of options seen below.
    Please see docuentation on github.com/moonso/genmod or 
    genmod/examples/readme.md for more information.
    """    
    logger = logging.getLogger(__name__)
    # logger = logging.getLogger("genmod.commands.annotate")
    ######### This is for logging the command line string #########
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [
        i+'='+str(values[i]) for i in values if values[i] and 
        i not in ['config','frame','logger']
    ]
    annotator_arguments = {}
    annotator_arguments['phased'] = phased
    annotator_arguments['strict'] = strict
    annotator_arguments['cadd_raw'] = cadd_raw
    annotator_arguments['whole_gene'] = whole_gene
    ##############################################################
    
    logger.info("Running GENMOD annotate version {0}".format(__version__))
    logger.debug("Arguments: {0}".format(', '.join(argument_list)))
    
    start_time_analysis = datetime.now()
    
    ######### Setup a variant parser #########
    init_log(logging.getLogger("vcf_parser"), loglevel="INFO")
    
    logger.debug("Setting up a variant parser")
    if variant_file == '-':
        variant_parser = VCFParser(
            fsock = sys.stdin, 
            split_variants=split_variants, 
            skip_info_check=True
            )
    else:
        variant_parser = VCFParser(
            infile = variant_file, 
            split_variants=split_variants,
            skip_info_check=True
            )
    logger.debug("Variant parser setup")
    
    vcf_individuals = variant_parser.individuals
    logger.info("Individuals found in vcf file: {}".format(', '.join(vcf_individuals)))
    
    head = variant_parser.metadata
    
    # Update version logging
    add_metadata(
        head,
        'version',    
        'genmod', 
        version=__version__, 
        command_line_string=' '.join(argument_list)
    )
    
    ######### Parse the ped file (if there is one) #########
    
    families = {}
    # These are the individuals that should be used in the analysis
    individuals = []
    if family_file:
        logger.info("Annotating Genetic Models")
        family_parser = FamilyParser(family_file, family_type)
        # The individuals in the ped file must be present in the variant file:
        families = family_parser.families
        
        for individual in family_parser.individuals:
            if individual not in vcf_individuals:
                logger.warning("All individuals in ped file must be in vcf "\
                                "file!")
                logger.warning("Individuals in PED file: {0}".format(
                                ', '.join(list(family_parser.individuals.keys()))))
                logger.warning("Individuals in VCF file: {0}".format(', '.join(vcf_individuals)))
                raise IOError() # Raise proper exception here
            else:
                individuals.append(individual)
        
        annotator_arguments['individuals'] = individuals
        annotator_arguments['families'] = families
        
        logger.info("Families used in analysis: {0}".format(
                    ','.join(list(families.keys()))))
        logger.info("Individuals included in analysis: {0}".format(
                    ','.join(list(family_parser.individuals.keys()))))
        
        
        logger.debug("Adding vcf metadata for GeneticModels")
        add_metadata(
            head,
            'info',
            'GeneticModels',
            annotation_number='.',
            entry_type='String',
            description="':'-separated list of genetic models for this variant."
        )

        logger.debug("Adding vcf metadata for ModelScore")
        add_metadata(
            head,
            'info',
            'ModelScore',
            annotation_number='1',
            entry_type='Integer',
            description="PHRED score for genotype models."
        )

        logger.debug("Adding vcf metadata for compounds")
        add_metadata(
            head,
            'info',
            'Compounds',
            annotation_number='.',
            entry_type='String',
            description=("List of compound pairs for this variant."
            "The list is splitted on ',' family id is separated with compounds"
            "with ':'. Compounds are separated with '|'.")
        )
    else:
        logger.info("No families included in analysis")

    ######### Read to the annotation data structures #########

    gene_trees = {}
    exon_trees = {}

    # If the variants are already annotated we do not need to redo the annotation
    if not vep:

        gene_trees, exon_trees = load_annotations(annotation_dir)
        annotator_arguments['exon_trees'] = exon_trees
        logger.debug("Adding vcf metadata for annotation")
        add_metadata(
            head,
            'info',
            'Annotation',
            annotation_number='.',
            entry_type='String',
            description='Annotates what feature(s) this variant belongs to.'
        )
    else:
        logger.info("Using VEP annotation")

    ######### Check which other annotations files that should be used in the analysis #########

    cadd_annotation = False

    if cadd_file:
        logger.info('Cadd file! {0}'.format(cadd_file))
        cadd_annotation = True
        annotator_arguments['cadd_file'] = cadd_file
    if cadd_1000g:
        logger.info('Cadd 1000G file! {0}'.format(cadd_1000g))
        cadd_annotation = True
        annotator_arguments['cadd_1000g'] = cadd_1000g
    if cadd_esp:
        logger.info('Cadd ESP6500 file! {0}'.format(cadd_esp))
        cadd_annotation = True
        annotator_arguments['cadd_ESP'] = cadd_esp
    if cadd_indels:
        logger.info('Cadd InDel file! {0}'.format(cadd_indels))
        cadd_annotation = True
        annotator_arguments['cadd_InDels'] = cadd_indels
    if cadd_exac:
        logger.info('Cadd ExAC file! {0}'.format(cadd_exac))
        cadd_annotation = True
        annotator_arguments['cadd_exac'] = cadd_exac


    if cadd_annotation:
        logger.debug("Adding vcf metadata for CADD score")
        add_metadata(
            head,
            'info',
            'CADD',
            annotation_number='A',
            entry_type='Float',
            description="The CADD relative score for this alternative."
        )
        if cadd_raw:
            logger.debug("Adding vcf metadata for CADD raw score")
            add_metadata(
                head,
                'info',
                'CADD_raw',
                annotation_number='A',
                entry_type='Float',
                description="The CADD raw score(s) for this alternative(s)."
            )

    if thousand_g:
        logger.info('1000G frequency file! {0}'.format(thousand_g))
        logger.debug("Adding vcf metadata for 1000G_freq")
        add_metadata(
            head,
            'info',
            '1000G_freq',
            annotation_number='A',
            entry_type='Float',
            description="Frequency in the 1000G database."
        )
        annotator_arguments['thousand_g'] = thousand_g

    if exac:
        logger.info('ExAC frequency file! {0}'.format(exac))
        add_metadata(
            head,
            'info',
            'ExAC_freq',
            annotation_number='A',
            entry_type='Float',
            description="Frequency in the ExAC database."
        )
        annotator_arguments['exac'] = exac


    ###################################################################
    ### The task queue is where all jobs(in this case batches that  ###
    ### represents variants in a region) is put. The consumers will ###
    ### then pick their jobs from this queue.                       ###
    ###################################################################

    logger.debug("Setting up a JoinableQueue for storing variant batches")
    variant_queue = JoinableQueue(maxsize=1000)
    logger.debug("Setting up a Queue for storing results from workers")
    results = Manager().Queue()


    num_model_checkers = processes
    #Adapt the number of processes to the machine that run the analysis
    if cadd_annotation:
        # We need more power when annotating cadd scores:
        # But if flag is used that overrides
        if num_model_checkers == min(4, cpu_count()):
            num_model_checkers = min(8, cpu_count())

    logger.info('Number of CPU:s {}'.format(cpu_count()))
    logger.info('Number of model checkers: {}'.format(num_model_checkers))

    # We use a temp file to store the processed variants
    logger.debug("Build a tempfile for printing the variants")
    temp_file = NamedTemporaryFile(delete=False)
    temp_file.close()
    # Open the temp file with codecs
    temporary_variant_file = open(
                                temp_file.name,
                                mode='w',
                                encoding='utf-8',
                                errors='replace'
                                )


    # These are the workers that do the heavy part of the analysis
    logger.info('Seting up the workers')
    model_checkers = [
        VariantAnnotator(
            variant_queue, 
            results, 
            **annotator_arguments
        )
        for i in range(num_model_checkers)
    ]

    logger.info('Starting the workers')
    for worker in model_checkers:
        logger.debug('Starting worker {0}'.format(worker))
        worker.start()

    # This process prints the variants to temporary files
    logger.info('Seting up the variant printer')
    var_printer = VariantPrinter(
                            results,
                            temporary_variant_file,
                            head,
                            mode='chromosome',
                            verbosity=verbose
                        )
    logger.info('Starting the variant printer process')
    var_printer.start()

    start_time_variant_parsing = datetime.now()

    # This process parses the original vcf and create batches to put in the variant queue:
    logger.info('Start parsing the variants')
    chromosome_list = get_batches(
                                variant_parser,
                                variant_queue,
                                gene_trees,
                            )

    logger.info('Put stop signs in the variant queue')
    for i in range(num_model_checkers):
        variant_queue.put(None)

    variant_queue.join()
    results.put(None)
    var_printer.join()

    temporary_variant_file.close()

    logger.info('Cromosomes found in variant file: {0}'.format(
        ','.join(chromosome_list)))
    logger.info('Models checked')

    logger.info("Start sorting the variants")
    sort_variants(temp_file.name, mode='chromosome')

    logger.info("Print the headers")
    print_headers(head, outfile, silent)

    logger.info("Start printing the variants")
    print_variants(temp_file.name, outfile, mode='modified',  silent=silent)

    logger.info("Removing temp file")
    os.remove(temp_file.name)
    logger.debug("Temp file removed")

    logger.info('Time for whole analyis: {0}'.format(
        str(datetime.now() - start_time_analysis)))



if __name__ == '__main__':
    from genmod import logger
    init_log(logger, loglevel="INFO")
    annotate()
