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

import shutil
import pkg_resources

from ped_parser import parser as ped_parser
from vcf_parser import parser as vcf_parser

from genmod import variant_consumer, variant_sorter, annotation_parser, variant_printer, variant_annotator, warning

version = pkg_resources.require("genmod")[0].version


def add_metadata(head, annotate_models=False, vep=False, cadd_annotation=False, cadd_raw=False, thousand_g=None,
                 exac=None, command_line_string=''):
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
        head.add_info('1000G_freq', 'A', 'Float', "Frequency in the 1000G database.")
    if exac:
        head.add_info('ExAC_freq', 'A', 'Float', "Frequency in the ExAC database.")
    
    # Update version logging
    head.add_version_tracking('genmod', version, datetime.now().strftime("%Y-%m-%d %H:%M"), command_line_string)
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
                chr_prefix, split_variants, processes, dbnfsp, verbose):
    """Annotate variants in a VCF file.\n
        The main function with genmod is to annotate genetic inheritance patterns for variants in families. 
        Use flag --family together with a .ped file to describe which individuals in the vcf you wish to check inheritance for in the analysis.
        Individuals that are not present in the ped file will not be considered in the analysis.\n
        It is also possible to use genmod without a family file. In this case the variants will be annotated with a variety of options seen below.
        Please see docuentation on github.com/moonso/genmod or genmod/examples/readme.md for more information.
    """    
    verbosity = verbose
    
    ######### This is for logging the command line string #########
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [i+'='+str(values[i]) for i in values if values[i] and i != 'config' and i != 'frame']
    
    
    if verbosity:
        start_time_analysis = datetime.now()
    
    
    ######### Setup a variant parser #########
    
    if variant_file == '-':
        variant_parser = vcf_parser.VCFParser(fsock = sys.stdin, split_variants=split_variants)
    else:
        variant_parser = vcf_parser.VCFParser(infile = variant_file, split_variants=split_variants)
    
    # These are the individuals in from the vcf file
    individuals = variant_parser.individuals
    
    head = variant_parser.metadata
    
    ######### Start to parse the ped file (if there is one) #########
    
    families = {}
    
    if family_file:
        family_parser = ped_parser.FamilyParser(family_file, family_type)
        # The individuals in the ped file must be present in the variant file:
        families = family_parser.families
        
        for individual in family_parser.individuals:
            if individual not in individuals:
                warning.warning('All individuals in ped file must be in vcf file! Aborting...')
                warning.warning('Individuals in PED file: %s' % ' '.join(list(family_parser.individuals.keys())))
                warning.warning('Individuals in VCF file: %s' % ' '.join(individuals))
                sys.exit()
    
    if verbosity:
        if family_file:
            print('Starting analysis of families: %s' % ','.join(list(families.keys())), file=sys.stderr)
            print('Individuals included in analysis: %s\n' % ','.join(list(family_parser.individuals.keys())), file=sys.stderr)
    ######### Connect to the annotations #########
    
    gene_trees = {}
    exon_trees = {}

    if not vep:
        
        if verbosity:
            print('Reading annotations...\n', file=sys.stderr)
        
        gene_db = pkg_resources.resource_filename('genmod', 'annotations/genes.db')
        exon_db = pkg_resources.resource_filename('genmod', 'annotations/exons.db')
        
        # If the user have provided an annotation use this
        if annotation_dir:
            gene_db = os.path.join(annotation_dir, 'genes.db')
            exon_db = os.path.join(annotation_dir, 'exons.db')
        
        try:
            with open(gene_db, 'rb') as f:
                gene_trees = pickle.load(f)
            with open(exon_db, 'rb') as g:
                exon_trees = pickle.load(g)
        except IOError as e:
            if verbosity:
                warning.warning('You need to build annotations! See documentation.')
                # It is possible to continue the analysis without annotation files
            pass
    
        if verbosity:
            print('Annotations used found in: %s, %s\n' % (gene_db, exon_db), file=sys.stderr)
    else:
        if verbosity:
            print('Using VEP annotation', file=sys.stderr)
        
    # sys.exit()
    
    
    ######### Check which other annotations files that should be used in the analysis #########
    
    cadd_annotation = False
    
    if cadd_file:
        if verbosity:
            print('Cadd file! %s' % cadd_file, file=sys.stderr)
        cadd_annotation = True
    if cadd_1000g:
        if verbosity:
            print('Cadd 1000G file! %s' % cadd_1000g, file=sys.stderr)
        cadd_annotation = True
    if cadd_esp:
        if verbosity:
            print('Cadd ESP6500 file! %s' % cadd_esp, file=sys.stderr)
        cadd_annotation = True
    if cadd_indels:
        if verbosity:
            print('Cadd InDel file! %s' % cadd_indels, file=sys.stderr)
        cadd_annotation = True
    if cadd_exac:
        if verbosity:
            print('Cadd ExAC file! %s' % cadd_exac, file=sys.stderr)
        cadd_annotation = True
    if thousand_g:
        if verbosity:
            print('1000G frequency file! %s' % thousand_g, file=sys.stderr)
    if exac:
        if verbosity:
            print('ExAC frequency file! %s' % exac, file=sys.stderr)
    if dbnfsp:
        if verbosity:
            print('dbNFSP file! %s' % dbnfsp, file=sys.stderr)
    
    
    ###################################################################
    ### The task queue is where all jobs(in this case batches that ###
    ### represents variants in a region) is put. The consumers will ###
    ### then pick their jobs from this queue. ###
    ###################################################################
    
    variant_queue = JoinableQueue(maxsize=1000)
    # The consumers will put their results in the results queue
    results = Manager().Queue()
    
    # Create a directory to keep track of temp files. The annotated variants will be printed to these files
    # They are then sorted and printed to a proper vcf
    temp_dir = mkdtemp()
    
    num_model_checkers = processes
    #Adapt the number of processes to the machine that run the analysis
    if cadd_annotation:
        # We need more power when annotating cadd scores:
        # But if flag is used that overrides
        if num_model_checkers == min(4, cpu_count()):
            num_model_checkers = min(8, cpu_count())
    
    if verbosity:
        print('Number of CPU:s %s' % cpu_count(), file=sys.stderr)
        print('Number of model checkers: %s' % num_model_checkers, file=sys.stderr)
    
    # These are the workers that do the heavy part of the analysis
    model_checkers = [variant_consumer.VariantConsumer(variant_queue, 
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
                                                        chr_prefix, 
                                                        strict, 
                                                        verbosity) 
                                                        for i in range(num_model_checkers)
                        ]
    
    for w in model_checkers:
        w.start()
    
    # This process prints the variants to temporary files
    var_printer = variant_printer.VariantPrinter(results, 
                                                    temp_dir, 
                                                    head, 
                                                    chr_prefix, 
                                                    verbosity
                                                )
    var_printer.start()
    
    
    # This process parses the original vcf and annotate the variants with regions:
    var_annotator = variant_annotator.VariantAnnotator(variant_parser, 
                                                        variant_queue,
                                                        individuals,
                                                        gene_trees, 
                                                        exon_trees, 
                                                        phased, 
                                                        vep, 
                                                        whole_gene, 
                                                        chr_prefix, 
                                                        verbosity
                                                    )
    
    
    start_time_variant_parsing = datetime.now()
    if verbosity:
        print('Start parsing the variants ... \n', file=sys.stderr)
    
    var_annotator.annotate()
    
    # Put stop signs in the variant queue
    for i in range(num_model_checkers):
        variant_queue.put(None)
    
    variant_queue.join()
    results.put(None)
    var_printer.join()
    
    
    chromosome_list = var_annotator.chromosomes
    
    if verbosity:
        print('Cromosomes found in variant file: %s \n' % ','.join(chromosome_list), file=sys.stderr)
        print('Models checked!\n', file=sys.stderr)
        print('Start sorting the variants:\n', file=sys.stderr)
        start_time_variant_sorting = datetime.now()
    
    # Add the new metadata to the headers:
    add_metadata(head, 
                 families, 
                 vep, 
                 cadd_annotation, 
                 cadd_raw, 
                 thousand_g, 
                 exac, 
                 ' '.join(argument_list)
                 )
    
    print_headers(head, outfile, silent)
    
    # We need to open and close files to read and write proper with codecs:
    sorted_temp = NamedTemporaryFile(delete=False)
    sorted_temp.close()
    sorted_file = sorted_temp.name
    for chromosome in chromosome_list:
        for temp_file in os.listdir(temp_dir):
            if temp_file.split('_')[0] == chromosome:
                var_sorter = variant_sorter.FileSort(os.path.join(temp_dir, temp_file), 
                                                        outfile=sorted_file, 
                                                        silent=silent
                                                    )
                var_sorter.sort()
    
    if verbosity:
        print('Sorting done!', file=sys.stderr)
        print('Time for sorting: %s \n' % str(datetime.now()-start_time_variant_sorting), file=sys.stderr)
        print('Time for whole analyis: %s' % str(datetime.now() - start_time_analysis), file=sys.stderr)
    
    print_variants(sorted_file, outfile, silent)
    
    # Remove all temp files:
    shutil.rmtree(temp_dir)
    


if __name__ == '__main__':
    annotate()