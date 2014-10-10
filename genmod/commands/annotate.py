#!/usr/bin/env python
# encoding: utf-8
"""
annotate.py

Command line tool for annotating variants with genmod.

Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import click
import inspect

try:
    import cPickle as pickle
except:
    import pickle

from multiprocessing import JoinableQueue, Manager, cpu_count
from codecs import open
from datetime import datetime
from tempfile import mkdtemp, TemporaryFile, NamedTemporaryFile
from pprint import pprint as pp

# import vcf

import shutil
import pkg_resources

from pysam import tabix_index, tabix_compress

from ped_parser import parser as ped_parser
from vcf_parser import parser as vcf_parser

import genmod
from genmod import variant_consumer, variant_sorter, annotation_parser, variant_printer, variant_annotator, warning

version = pkg_resources.require("genmod")[0].version


def get_family(family_file, family_type):
    """Return the family"""
    
    my_family_parser = ped_parser.FamilyParser(family_file, family_type)
    # Stupid thing but for now when we only look at one family
    
    try:
        return my_family_parser.families.popitem()[1]
    except KeyError as e:
        raise SyntaxError("Something wrong with the pedigree file? No families found")
    
    return

def add_metadata(head, annotate_models=False, vep=False, cadd_annotation=False, cadd_raw=False, thousand_g=None, command_line_string=''):
    """Add metadata for the information added by this script."""
    # Update INFO headers
    if not vep:
        head.add_info('Annotation', '.', 'String', 'Annotates what feature(s) this variant belongs to.')
    if annotate_models:
        head.add_info('Compounds', '.', 'String', "':'-separated list of compound pairs for this variant.")
        head.add_info('GeneticModels', '.', 'String', "':'-separated list of genetic models for this variant.")
        head.add_info('ModelScore', '1', 'Integer', "PHRED score for genotype models.")
    if cadd_annotation:
        head.add_info('CADD', 'A', 'Float', "The CADD relative score for this alternative.")
        if cadd_raw:
            head.add_info('CADD_raw', 'A', 'Float', "The CADD raw score(s) for this alternative(s).")
    
    if thousand_g:
        head.add_info('1000GMAF', 'A', 'Float', "Frequency in the 1000G database.")
    
    # Update version logging
    head.add_version_tracking('genmod', version, str(datetime.now()), command_line_string)
    return

def print_headers(head, outfile, silent=False):
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


def print_variants(sorted_variants, outfile, silent=False):
    """Print the variants to a results file or stdout."""
    
    with open(sorted_variants, mode='r', encoding='utf-8') as f:
        if outfile:
            with open(outfile, 'a', encoding='utf-8') as g:
                for variant in f:
                    g.write(variant)
        else:
            if not silent:
                for line in f:
                    print(line.rstrip())
    return


def check_tabix_index(compressed_file, file_type='cadd', verbose=False):
    """Check if a compressed file have a tabix index. If not build one."""
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
    return



###           This is for annotating the variants       ###


@click.command()
@click.argument('variant_file', 
                    nargs=1, 
                    type=click.Path(),
                    metavar='<vcf_file> or -'
)
@click.option('--family_file', '-fam',
                    nargs=1, 
                    type=click.Path(exists=True),
                    metavar='<ped_file>'
)
@click.option('-f' ,'--family_type', 
                type=click.Choice(['ped', 'alt', 'cmms', 'mip']), 
                default='ped',
                help='If the analysis use one of the known setups, please specify which one.'
)
@click.option('--vep', 
                    is_flag=True,
                    help='If variants are annotated with the Variant Effect Predictor.'
)
@click.option('--chr_prefix', 
                    is_flag=True,
                    help='If chr prefix is used in vcf file.'
)
@click.option('-p' ,'--phased', 
                    is_flag=True,
                    help='If data is phased use this flag.'
)
@click.option('-strict' ,'--strict', 
                    is_flag=True,
                    help='If strict model annotations should be used(see documentation).'
)
@click.option('-s' ,'--silent', 
                    is_flag=True,
                    help='Do not print the variants.'
)
@click.option('-g' ,'--whole_gene', 
                    is_flag=True,
                    help="""If compounds should be checked in the whole gene regions. 
                    Not only exonic/splice sites."""
)
@click.option('-craw' ,'--cadd_raw', 
                    is_flag=True,
                    help="""If the raw cadd scores should be annotated."""
)
@click.option('-a' ,'--annotation_dir', 
                    type=click.Path(exists=True), 
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
                    help="""Specify the path to a bgzipped cadd file with variant scores.
                            If no index is present it will be created."""
)
@click.option('--cadd_1000g',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped cadd file with variant scores 
                            for all 1000g variants. If no index is present a new index
                            will be created."""
)
@click.option('--cadd_esp',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped cadd file with variant scores 
                            for all ESP6500 variants. If no index is present a new index
                            will be created."""
)
@click.option('--cadd_indels',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped cadd file with variant scores 
                            for all CADD InDel variants. If no index is present a new index
                            will be created."""
)
@click.option('--thousand_g',
                    type=click.Path(exists=True), 
                    help="""Specify the path to a bgzipped vcf file frequency info of all 
                            1000g variants. If no index is present a new index
                            will be created."""
)
@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def annotate(family_file, variant_file, family_type, vep, silent, phased, strict, cadd_raw, whole_gene, 
                annotation_dir, cadd_file, cadd_1000g, cadd_esp, cadd_indels, thousand_g, outfile,
                chr_prefix, verbose):
    """Annotate variants in a VCF file.
        It is possible to annotate from sources shown as options only
        If a ped file is provided then the genetic inheritance patterns for all individuals are followed.
        Individuals that are not present in ped file will not be considered in the analysis.
    """    
    verbosity = verbose
    
    # This is for logging the command line string:
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [i+'='+str(values[i]) for i in values if values[i] and i != 'config' and i != 'frame']
    gene_db = pkg_resources.resource_filename('genmod', 'annotations/genes.db')
    exon_db = pkg_resources.resource_filename('genmod', 'annotations/exons.db')
    
    cadd_annotation = False
    
    # Connect to the annotation directorues that should be used:
    if annotation_dir:
        gene_db = os.path.join(annotation_dir, 'genes.db')
        exon_db = os.path.join(annotation_dir, 'exons.db')
    
    # Setup a variant parser
    if variant_file == '-':
        variant_parser = vcf_parser.VCFParser(fsock = sys.stdin)
    else:
        variant_parser = vcf_parser.VCFParser(infile = variant_file)
    
    head = variant_parser.metadata
    
    if verbosity:
        start_time_analysis = datetime.now()
    try:
        with open(gene_db, 'rb') as f:
            gene_trees = pickle.load(f)
        with open(exon_db, 'rb') as g:
            exon_trees = pickle.load(g)
    except IOError as e:
        if verbosity:
            warning('You need to build annotations! See documentation.')
            # It is possible to continue the analysis without annotation files
        pass
        
        
    if family_file:
        family = get_family(family_file, family_type)
        # There has to be same individuals in ped file and variant file:
        if set(family.individuals.keys()) != set(variant_parser.individuals):
            warning.warning('There must be same individuals in ped file and vcf file! Aborting...')
            warning.warning('Individuals in PED file: %s' % '\t'.join(list(family.individuals.keys())))
            warning.warning('Individuals in VCF file: %s' % '\t'.join(list(variant_parser.individuals)))
            sys.exit()
    else:
        family = False
    
    
    if cadd_file:
        if verbosity:
            click.echo('Cadd file! %s' % cadd_file)
        check_tabix_index(cadd_file, 'cadd', verbosity)
        cadd_annotation = True
    if cadd_1000g:
        if verbosity:
            click.echo('Cadd 1000G file! %s' % cadd_1000g)
        check_tabix_index(cadd_1000g, 'cadd', verbosity)
        cadd_annotation = True
    if cadd_esp:
        if verbosity:
            click.echo('Cadd ESP6500 file! %s' % cadd_esp)
        check_tabix_index(cadd_esp, 'cadd', verbosity)
        cadd_annotation = True
    if cadd_indels:
        if verbosity:
            click.echo('Cadd InDel file! %s' % cadd_indels)
        check_tabix_index(cadd_indels, 'cadd', verbosity)
        cadd_annotation = True
    if thousand_g:
        if verbosity:
            click.echo('1000G frequency file! %s' % thousand_g)
        check_tabix_index(thousand_g, 'vcf', verbosity)
    
    
    ###################################################################
    ### The task queue is where all jobs(in this case batches that ###
    ### represents variants in a region) is put. The consumers will ###
    ### then pick their jobs from this queue. ###
    ###################################################################
    
    variant_queue = JoinableQueue(maxsize=1000)
    # The consumers will put their results in the results queue
    results = Manager().Queue()
    
    # Create a directory to keep track of temp files
    temp_dir = mkdtemp()
    #Adapt the number of processes to the machine that run the analysis
    if cadd_annotation:
        num_model_checkers = cpu_count()
    else:
        num_model_checkers = min(4, cpu_count())
    
    if verbosity:
        print('Number of CPU:s %s' % cpu_count())
    
    # These are the workers that do the analysis
    model_checkers = [variant_consumer.VariantConsumer(variant_queue, results, family,
                        phased, vep, cadd_raw, cadd_file, cadd_1000g, cadd_esp, cadd_indels,
                        thousand_g, chr_prefix, strict, verbosity) for i in range(num_model_checkers)]
    
    for w in model_checkers:
        w.start()
    
    # This process prints the variants to temporary files
    var_printer = variant_printer.VariantPrinter(results, temp_dir, head, verbosity)
    var_printer.start()
    
    if verbosity:
        print('Start parsing the variants ...')
        print('')
        start_time_variant_parsing = datetime.now()
    
    # For parsing the vcf:
    var_annotator = variant_annotator.VariantAnnotator(variant_parser, variant_queue,
                        gene_trees, exon_trees, phased, vep, whole_gene, verbosity)
    
    var_annotator.annotate()
    
    for i in range(num_model_checkers):
        variant_queue.put(None)
    
    variant_queue.join()
    results.put(None)
    var_printer.join()
    
    chromosome_list = var_annotator.chromosomes
    
    if verbosity:
        print('Cromosomes found in variant file: %s' % ','.join(chromosome_list))
        print('Models checked!')
        print('Start sorting the variants:')
        print('')
        start_time_variant_sorting = datetime.now()
    
    # Add the new metadata to the headers:
    add_metadata(head, family, vep, cadd_annotation, cadd_raw, thousand_g, ' '.join(argument_list))
    
    print_headers(head, outfile, silent)
    
    sorted_temp = NamedTemporaryFile(delete=False)
    sorted_temp.close()
    sorted_file = sorted_temp.name
    for chromosome in chromosome_list:
        for temp_file in os.listdir(temp_dir):
            if temp_file.split('_')[0] == chromosome:
                var_sorter = variant_sorter.FileSort(os.path.join(temp_dir, temp_file), outfile=sorted_file, silent=silent)
                var_sorter.sort()
    
    if verbosity:
        print('Sorting done!')
        print('Time for sorting: %s' % str(datetime.now()-start_time_variant_sorting))
        print('')
        print('Time for whole analyis: %s' % str(datetime.now() - start_time_analysis))
    
    print_variants(sorted_file, outfile, silent)
    
    # Remove all temp files:
    shutil.rmtree(temp_dir)
    


if __name__ == '__main__':
    annotate()