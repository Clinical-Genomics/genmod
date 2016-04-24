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
import shutil
import itertools

from multiprocessing import JoinableQueue, Manager, cpu_count
from codecs import open
from datetime import datetime
from tempfile import NamedTemporaryFile

from ped_parser import FamilyParser

from genmod import (__version__)

from genmod.utils import (get_batches, VariantPrinter, check_individuals)
from genmod.annotate_models import (VariantAnnotator)
from genmod.vcf_tools import (add_metadata, print_headers, sort_variants, 
print_variant, HeaderParser)

from .utils import (temp_dir, silent, outfile, processes, variant_file, 
                    family_file, family_type)

@click.command()
@variant_file
@family_file
@family_type
@click.option('-r', '--reduced_penetrance',
                    nargs=1, 
                    type=click.File('r'),
                    metavar='<tsv_file>',
                    help='File with gene ids that have reduced penetrance.'
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
@processes
@silent
@click.option('-w', '--whole_gene',
                    is_flag=True,
                    help='If compounds should be checked over the whole gene.'
)
@click.option('-k' ,'--keyword', 
                    default="Annotation",
                    help="""What annotation keyword that should be used when 
                    searching for features."""
)
@outfile
@temp_dir
def models(variant_file, family_file, family_type, reduced_penetrance, vep,
keyword, phased, strict, silent, processes, whole_gene, outfile, temp_dir):
    """
    Annotate genetic models for vcf variants. 
    
    Checks what patterns of inheritance that are followed in a VCF file.
    The analysis is family based so each family that are specified in the family
    file and exists in the variant file will get it's own annotation.
    """
    logger = logging.getLogger(__name__)
    
    ######### This is for logging the command line string #########
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [
        i+'='+str(values[i]) for i in values if values[i] and 
        i not in ['frame']
    ]
    
    ###########################################################################
    
    logger.info("Running GENMOD annotate version {0}".format(__version__))
    logger.debug("Arguments: {0}".format(', '.join(argument_list)))
    
    reduced_penetrance_genes = set()
    nr_reduced_penetrance_genes = 0
    if reduced_penetrance:
        logger.info("Found file with genes that have reduced penetrance")
        for line in reduced_penetrance:
            if not line.startswith('#'):
                nr_reduced_penetrance_genes += 1
                gene_id = line.rstrip().split()[0]
                logger.debug("Adding gene {0} to reduced penetrance genes".format(
                    gene_id
                ))
                reduced_penetrance_genes.add(
                    gene_id
                )
    
        logger.info("Found {0} genes with reduced penetrance".format(
            nr_reduced_penetrance_genes))
    
    
    if not family_file:
        logger.warning("Please provide a family file with -f/--family_file")
        logger.info("Exiting")
        sys.exit(1)
    
    logger.info("Setting up a family parser")
    family_parser = FamilyParser(family_file, family_type)
    logger.debug("Family parser done")
    
    families = {}
    logger.info("Check if the familys have any affected")
    for family_id in family_parser.families:
        found_affected = False
        family_obj = family_parser.families[family_id]
        for ind_id in family_obj.individuals:
            ind_obj = family_obj.individuals[ind_id]
            if ind_obj.affected:
                found_affected = True
        
        if found_affected:
            families[family_id] = family_obj
        else:
            logger.warning("No affected individuals found for family {0}."\
                           " Skipping family.".format(family_id))
    
    if not families:
        logger.warning("Please provide at least one family with affected individuals")
        sys.exit(0)
    # The individuals in the ped file must be present in the variant file:
    logger.info("Families used in analysis: {0}".format(
                    ','.join(list(families.keys()))))
    logger.info("Individuals included in analysis: {0}".format(
                    ','.join(list(family_parser.individuals.keys()))))
    
    
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
    
    #Add the first variant to the iterator
    variant_file = itertools.chain([line], variant_file)
    
    if vep:
        if not "CSQ" in head.info_dict:
            logger.warning("vep flag is used but there is no CSQ field specified in header")
            logger.info("Please check VCF file")
            logger.info("Exiting...")
            sys.exit(1)
        else:
            logger.info("Using VEP annotation")
    else:
        if not keyword in head.info_dict:
            logger.warning("Annotation key {0} could not be found in VCF header".format(keyword))
            logger.info("Please check VCF file")
            logger.info("Exiting...")
            sys.exit(1)
        else:
            logger.info("Using {0} annotation".format(keyword))
        
    
    if "GeneticModels" in head.info_dict:
        logger.warning("Genetic models are already annotated according to vcf"\
        " header.")
        logger.info("Exiting...")
        sys.exit(1)
    
    logger.info("Adding genmod version to vcf header")
    head.add_version_tracking(
                    info_id='genmod',
                    version=__version__,
                    date=datetime.now().strftime("%Y-%m-%d %H:%M"),
                    command_line=' '.join(argument_list)
                )
    
    logger.debug("Version added")
    logger.info("Adding genetic models to vcf header")
    add_metadata(
        head,
        'info',
        'GeneticModels',
        annotation_number='.',
        entry_type='String',
        description="':'-separated list of genetic models for this variant."
    )
    
    logger.debug("Genetic models added")
    logger.info("Adding model score to vcf header")
    add_metadata(
        head,
        'info',
        'ModelScore',
        annotation_number='.',
        entry_type='String',
        description="PHRED score for genotype models."
    )
    logger.debug("Model score added")
    
    logger.info("Adding Compounds to vcf header")
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
    logger.debug("Compounds added")
    
    vcf_individuals = head.individuals
    logger.debug("Individuals found in vcf file: {}".format(', '.join(vcf_individuals)))
    

    start_time_analysis = datetime.now()
    
    try:
        check_individuals(family_parser.individuals, vcf_individuals)
    except IOError as e:
        logger.error(e)
        logger.info("Individuals in PED file: {0}".format(
                        ', '.join(family_parser.individuals)))
        logger.info("Individuals in VCF file: {0}".format(', '.join(vcf_individuals)))
        logger.info("Exiting...")
        sys.exit(1)

    analysis_individuals = list(family_parser.individuals.keys())
    
    logger.info("Individuals used in analysis: {0}".format(
        ', '.join(analysis_individuals)))
    
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
    logger.info('Number of CPU:s {}'.format(cpu_count()))
    logger.info('Number of model checkers: {}'.format(num_model_checkers))


    # These are the workers that do the heavy part of the analysis
    logger.info('Seting up the workers')
    model_checkers = [
        VariantAnnotator(
            task_queue=variant_queue,
            results_queue=results,
            families=families,
            individuals=analysis_individuals,
            phased=phased,
            strict=strict,
            whole_gene=whole_gene,
            vep=vep,
            reduced_penetrance_genes = reduced_penetrance_genes
        )
        for i in range(num_model_checkers)
    ]
    logger.info('Starting the workers')
    for worker in model_checkers:
        logger.debug('Starting worker {0}'.format(worker))
        worker.start()

    # This process prints the variants to temporary files
    logger.info('Seting up the variant printer')
    if len(model_checkers) == 1:
        print_headers(head=head, outfile=outfile, silent=silent)
        variant_printer = VariantPrinter(
                task_queue=results,
                head=head,
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
    logger.info('Start parsing the variants')
    chromosome_list = get_batches(
                                variants = variant_file,
                                batch_queue = variant_queue,
                                header = head,
                                vep = vep,
                                annotation_keyword = keyword
                            )
    
    logger.debug("Put stop signs in the variant queue")
    for i in range(num_model_checkers):
        variant_queue.put(None)
    
    variant_queue.join()
    results.put(None)
    variant_printer.join()
    
    if len(model_checkers) > 1:
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
    
        logger.debug("Removing temp file")
        os.remove(temp_file.name)
        logger.debug("Temp file removed")

    logger.info('Time for whole analyis: {0}'.format(
        str(datetime.now() - start_time_analysis)))
    

if __name__ == '__main__':
    from genmod import logger as root_logger
    from genmod.log import init_log
    init_log(root_logger, loglevel="INFO")
    
    models()