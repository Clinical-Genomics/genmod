#!/usr/bin/env python
# encoding: utf-8
"""
annotate.py

Command line tool for annotating variants with genmod.

Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os
import click
import inspect


from multiprocessing import JoinableQueue, Manager, cpu_count
from codecs import open
from datetime import datetime
from tempfile import mkdtemp, TemporaryFile, NamedTemporaryFile
from pprint import pprint as pp

import shutil
import pkg_resources

from ped_parser import FamilyParser
from vcf_parser import VCFParser

from genmod import (VariantConsumer, annotation_parser, 
                    VariantPrinter, get_batches, sort_variants,
                    load_annotations, print_headers, print_variants,
                    add_metadata, warning)

VERSION = pkg_resources.require("genmod")[0].version


def check_tabix_index(compressed_file, file_type='cadd', verbose=False):
    """
    Check if a compressed file have a tabix index, if not build one.
    
    Args:
        compressed_file : Path to a file that is assumed to be compressed.
        file_type   : The type of the file.
        verobe      : Increase output verbosity
    
    Returns:
        0 if everythong went ok.
    
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
                    type=click.Path(exists=True),
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
    """Annotate variants in a VCF file.\n
        The main function with genmod is to annotate genetic inheritance patterns for variants in families. 
        Use flag --family together with a .ped file to describe which individuals in the vcf you wish to check inheritance for in the analysis.
        Individuals that are not present in the ped file will not be considered in the analysis.\n
        It is also possible to use genmod without a family file. In this case the variants will be annotated with a variety of options seen below.
        Please see docuentation on github.com/moonso/genmod or genmod/examples/readme.md for more information.
    """    
    
    ######### This is for logging the command line string #########
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [i+'='+str(values[i]) for i in values if values[i] and i != 'config' and i != 'frame']
    
    if verbose:
        print('\nRunning GENMOD annotate version %s \n' % VERSION ,file=sys.stderr)
    
    start_time_analysis = datetime.now()
    
    ######### Setup a variant parser #########
    
    if variant_file == '-':
        variant_parser = VCFParser(fsock = sys.stdin, split_variants=split_variants)
    else:
        variant_parser = VCFParser(infile = variant_file, split_variants=split_variants)
    
    # These are the individuals in from the vcf file
    individuals = variant_parser.individuals
    
    head = variant_parser.metadata
    
    # Update version logging
    add_metadata(
        head,
        'version',    
        'genmod', 
        version=VERSION, 
        command_line_string=' '.join(argument_list)
    )
    
    ######### Parse the ped file (if there is one) #########
    
    families = {}
    
    if family_file:
        family_parser = FamilyParser(family_file, family_type)
        # The individuals in the ped file must be present in the variant file:
        families = family_parser.families
        
        for individual in family_parser.individuals:
            if individual not in individuals:
                warning('All individuals in ped file must be in vcf file! Aborting...')
                warning('Individuals in PED file: %s' % ' '.join(list(family_parser.individuals.keys())))
                warning('Individuals in VCF file: %s' % ' '.join(individuals))
                print('Exiting...', file=sys.stderr)
                sys.exit()
        
        add_metadata(
            head,
            'info',
            'GeneticModels', 
            annotation_number='.', 
            entry_type='String', 
            description="':'-separated list of genetic models for this variant."
        )
        add_metadata(
            head,
            'info',
            'ModelScore', 
            annotation_number='1', 
            entry_type='Integer', 
            description="PHRED score for genotype models."
        )
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
        
    
    if verbose:
        if family_file:
            print('Starting analysis of families: %s' % 
                    ','.join(list(families.keys())), file=sys.stderr)
            print('Individuals included in analysis: %s\n' % 
                    ','.join(list(family_parser.individuals.keys())), file=sys.stderr)
    ######### Read to the annotation data structures #########
    
    gene_trees = {}
    exon_trees = {}
    
    # If the variants are already annotated we do not need to redo the annotation
    if not vep:
        
        gene_trees, exon_trees = load_annotations(annotation_dir, verbose)
        
        add_metadata(
            head,
            'info',
            'Annotation', 
            annotation_number='.', 
            entry_type='String', 
            description='Annotates what feature(s) this variant belongs to.'
        )
    else:
        if verbose:
            print('Using VEP annotation', file=sys.stderr)
    
    
    ######### Check which other annotations files that should be used in the analysis #########
    
    cadd_annotation = False
    
    if cadd_file:
        if verbose:
            print('Cadd file! %s' % cadd_file, file=sys.stderr)
        cadd_annotation = True
    if cadd_1000g:
        if verbose:
            print('Cadd 1000G file! %s' % cadd_1000g, file=sys.stderr)
        cadd_annotation = True
    if cadd_esp:
        if verbose:
            print('Cadd ESP6500 file! %s' % cadd_esp, file=sys.stderr)
        cadd_annotation = True
    if cadd_indels:
        if verbose:
            print('Cadd InDel file! %s' % cadd_indels, file=sys.stderr)
        cadd_annotation = True
    if cadd_exac:
        if verbose:
            print('Cadd ExAC file! %s' % cadd_exac, file=sys.stderr)
        cadd_annotation = True
    
    
    if cadd_annotation:
        add_metadata(
            head,
            'info',
            'CADD', 
            annotation_number='A',
            entry_type='Float', 
            description="The CADD relative score for this alternative."
        )
        if cadd_raw:
            add_metadata(
                head,
                'info',
                'CADD_raw', 
                annotation_number='A',
                entry_type='Float', 
                description="The CADD raw score(s) for this alternative(s)."
            )
        
    if thousand_g:
        if verbose:
            print('1000G frequency file! %s' % thousand_g, file=sys.stderr)
        add_metadata(
            head,
            'info',
            '1000G_freq', 
            annotation_number='A', 
            entry_type='Float', 
            description="Frequency in the 1000G database."
        )
        
    if exac:
        if verbose:
            print('ExAC frequency file! %s' % exac, file=sys.stderr)
        add_metadata(
            head,
            'info',
            'ExAC_freq', 
            annotation_number='A', 
            entry_type='Float', 
            description="Frequency in the ExAC database."
        )
        
    if dbnfsp:
        if verbose:
            print('dbNFSP file! %s' % dbnfsp, file=sys.stderr)
    
    
    ###################################################################
    ### The task queue is where all jobs(in this case batches that  ###
    ### represents variants in a region) is put. The consumers will ###
    ### then pick their jobs from this queue.                       ###
    ###################################################################
    
    variant_queue = JoinableQueue(maxsize=1000)
    # The consumers will put their results in the results queue
    results = Manager().Queue()
    
    num_model_checkers = processes
    #Adapt the number of processes to the machine that run the analysis
    if cadd_annotation:
        # We need more power when annotating cadd scores:
        # But if flag is used that overrides
        if num_model_checkers == min(4, cpu_count()):
            num_model_checkers = min(8, cpu_count())
    
    if verbose:
        print('Number of CPU:s %s' % cpu_count(), file=sys.stderr)
        print('Number of model checkers: %s' % num_model_checkers, file=sys.stderr)
    
    # We use a temp file to store the processed variants
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
    model_checkers = [
                    VariantConsumer(
                                variant_queue, 
                                results,
                                families,
                                phased,
                                vep,
                                cadd_raw,
                                cadd_file,
                                cadd_1000g,
                                cadd_exac,
                                cadd_esp,
                                cadd_indels,
                                thousand_g,
                                exac,
                                dbnfsp,
                                strict,
                                verbose
                            )
                        for i in range(num_model_checkers)
                        ]
    
    for w in model_checkers:
        w.start()
    
    # This process prints the variants to temporary files
    var_printer = VariantPrinter(
                            results, 
                            temporary_variant_file, 
                            head,
                            mode='chromosome',
                            verbosity=verbose
                        )
    var_printer.start()
    
    start_time_variant_parsing = datetime.now()
    
    if verbose:
        print('Start parsing the variants ... \n', file=sys.stderr)
    
    # This process parses the original vcf and create batches to put in the variant queue:
    
    chromosome_list = get_batches(
                                variant_parser, 
                                variant_queue,
                                individuals,
                                gene_trees, 
                                exon_trees, 
                                phased, 
                                vep, 
                                whole_gene, 
                                verbose
                            )
    
    # Put stop signs in the variant queue
    for i in range(num_model_checkers):
        variant_queue.put(None)
    
    variant_queue.join()
    results.put(None)
    var_printer.join()
    
    temporary_variant_file.close()
        
    if verbose:
        print('Cromosomes found in variant file: %s \n' % ','.join(chromosome_list), file=sys.stderr)
        print('Models checked!\n', file=sys.stderr)
    
    sort_variants(temp_file.name, mode='chromosome', verbose=verbose)
    
    print_headers(head, outfile, silent)
    
    print_variants(temp_file.name, outfile, mode='modified',  silent=silent)
    
    # Remove all temp files:
    os.remove(temp_file.name)
    
    if verbose:
        print('Time for whole analyis: %s' % str(datetime.now() - start_time_analysis), file=sys.stderr)
    


if __name__ == '__main__':
    annotate()