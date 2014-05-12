#!/usr/bin/env python
# encoding: utf-8
"""
run_genmod.py

Script for annotating genetic models in variant files.

Created by Måns Magnusson on 2014-01-21.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
from multiprocessing import JoinableQueue, Manager, cpu_count
from codecs import open
from datetime import datetime
from tempfile import mkdtemp
import shutil
import pkg_resources
import genmod

try:
    import cPickle as pickle
except:
    import pickle


from pysam import tabix_index, tabix_compress

from ped_parser import parser
from vcf_parser import vcf_parser

from genmod.utils import variant_consumer, variant_sorter, annotation_parser, variant_printer, variant_annotator


def get_annotation(args):
    """Initialize the annotation parser and return the gene and exon trees."""
    
    gene_trees = {}
    exon_trees = {}
    
    if args.annotation_file:
        anno_file = args.annotation_file[0]
        
        if args.verbose:
            print('Parsing annotation ...')
            print('')
            start_time_annotation = datetime.now()
        
        anno_file_name, anno_file_extension = os.path.splitext(anno_file)
        anno_zipped = False
        annotation = 'gene_pred'
        
        if anno_file_extension == '.gz':
            anno_zipped = True
            anno_file_name, anno_file_extension = os.path.splitext(anno_file_name)
        
        if args.annotation_type:
            annotation = args.annotation_type[0]
        
        else:
            
            if  anno_file_extension[1:] == 'bed':
                annotation = 'bed'
            if  file_extension[1:] == 'ccds':
                annotation = 'ccds'
            if anno_file_extension[1:] == 'gtf':
                annotation = 'gtf'
            if anno_file_extension[1:] in ['ref_gene', 'gene_pred']:
                annotation = 'gene_pred'
        
        anno_parser = annotation_parser.AnnotationParser(anno_file, annotation,
                     zipped=anno_zipped, splice_padding = args.splice_padding[0], verbosity=args.verbose)
        
        gene_trees = anno_parser.gene_trees
        exon_trees = anno_parser.exon_trees
        
        if args.verbose:
            print('Annotation Parsed!')
            print('Cromosomes found in annotation file: %s' % ','.join(list(gene_trees.keys())))
            print('Time to parse annotation: %s' % (datetime.now() - start_time_annotation))
            print('')
            
    
    else:
        annopath = os.path.join(os.path.split(os.path.dirname(genmod.__file__))[0], 'annotations/')
        gene_db = os.path.join(annopath, 'genes.db')
        exon_db = os.path.join(annopath, 'exons.db')
        
        try:
            with open(gene_db, 'rb') as f:
                gene_trees = pickle.load(f)
            with open(exon_db, 'rb') as g:
                exon_trees = pickle.load(g)
        except FileNotFoundError:
            print('You need to build annotations! See documentation.')
            pass
                
    return gene_trees, exon_trees 

def check_file_existence(infile):
    """Check is the file exists. Quit if something is wrong."""
    if not os.path.isfile(infile):
        print('The file %s does not exist!!!\nPlease check what is wrong and rerun.\nExiting...' % infile)
        sys.exit()

def get_family(args):
    """Return the family"""
    family_type = args.family_type[0]
    family_file = args.family_file[0]
    
    my_family_parser = parser.FamilyParser(family_file, family_type)
    # Stupid thing but for now when we only look at one family
    return my_family_parser.families.popitem()[1]

def add_metadata(head, args):
    """Add metadata for the information added by this script."""
    head.add_info('ANN', '.', 'String', 'Annotates what feature(s) this variant belongs to.')
    head.add_info('Comp', '.', 'String', "':'-separated list of compound pairs for this variant.")
    head.add_info('GM', '.', 'String', "':'-separated list of genetic models for this variant.")
    head.add_info('MS', '1', 'Integer', "PHRED score for genotype models.")
    if args.cadd_file[0] or args.cadd_1000g[0]:
        head.add_info('CADD', '1', 'Float', "The CADD relative score for this alternative.")
    if args.thousand_g[0]:
        head.add_info('1000G_freq', '1', 'Float', "Frequency in the 1000G database.")
    return

def print_headers(head,args):
    """Print the headers to a results file."""
    if args.outfile[0]:
        with open(args.outfile[0], 'w', encoding='utf-8') as f: 
            for head_count in head.print_header():
                f.write(head_count+'\n')
    else:
        if not args.silent:
            for line in head.print_header():
                print(line)
    return

def main():
    
    info_string = """Individuals that are not present in ped file will not be considered in the analysis."""
    
    parser = argparse.ArgumentParser(description="Annotate genetic models in variant files..")
    
    parser.add_argument('family_file', 
        type=str, nargs=1, 
        help='A pedigree file in .ped format.'
    )
    parser.add_argument('variant_file', 
        type=str, nargs=1, 
        help='A variant file. Default is vcf format.'
    )

    parser.add_argument('-an', '--annotation_file', 
        type=str, nargs=1, 
        help='An annotation file. Default is gene_pred format.'
    )
    
    parser.add_argument('-at', '--annotation_type',  
        type=str, nargs=1, 
        choices=['bed', 'ccds', 'gtf', 'gene_pred'], 
        default=['gene_pred'],
        help='Specify the format of the annotation file. gene_pred is default (this is the format of the refgene files.)'
    )    
    
    parser.add_argument('-vep', '--vep', 
        action="store_true", 
        help='If variants are annotated with the Variant Effect Predictor.'
    )
    
    parser.add_argument('--version', 
        action="version", 
        version=pkg_resources.require("genmod")[0].version
    )
    
    parser.add_argument('-v', '--verbose', 
        action="store_true", 
        help='Increase output verbosity.'
    )
    
    parser.add_argument('-chr', '--chr_prefix', 
        action="store_true", 
        help='If chr prefix is used in vcf.'
    )
    
    parser.add_argument('-s', '--silent', 
        action="store_true", 
        help='Do not print the variants.'
    )
    
    parser.add_argument('-phased', '--phased', 
        action="store_true", 
        help='If data is phased use this flag.'
    )
    
    parser.add_argument('-gene', '--whole_gene', 
        action="store_true", 
        help='If compounds should be checked in the whole gene regions. Not only exonic/splice sites.'
    )
    
    parser.add_argument('-o', '--outfile', 
        type=str, nargs=1, default=[None],
        help='Specify the path to a file where results should be stored.'
    )
    
    parser.add_argument('-splice', '--splice_padding', 
        type=int, nargs=1, default=[2],
        help='Specify the the number of bases that the exons should be padded with.'
    )
    
    parser.add_argument('-cadd', '--cadd_file', 
        type=str, nargs=1, default=[None],
        help='Specify the path to a bgzipped cadd file with variant scores.\
            If no index is present it will be created.'
    )    
    
    parser.add_argument('-c1kg', '--cadd_1000g', 
        type=str, nargs=1, default=[None],
        help='Specify the path to a bgzipped cadd file with variant scores for all 1000g variants.\
            If no index is present it will be created.'
    )    

    parser.add_argument('-kg', '--thousand_g', 
        type=str, nargs=1, default=[None],
        help='Specify the path to a bgzipped vcf file frequency info of all 1000g variants.\
            If no index is present it will be created.'
    )
    parser.add_argument('-family', '--family_type', 
        type=str, nargs=1, default=['ped'], 
        choices=['ped', 'alt', 'cmms', 'mip'],
        help='If the analysis use one of the known setups, please specify which.'
    )
    
    args = parser.parse_args()
    var_file = args.variant_file[0]
    file_name, file_extension = os.path.splitext(var_file)
    
    start_time_analysis = datetime.now()
    
    check_file_existence(var_file)
    check_file_existence(args.family_file[0])
    # Start by parsing at the pedigree file:
    
    my_family = get_family(args)
        
    # Parse the annotation file and make annotation trees:
    
    gene_trees = {}
    exon_trees = {}
    
    #If files are annotated with vep we do not need to annotate again:
    if not args.vep:
        # If annotation file is provided, use this annotation and make new .db files
        gene_trees, exon_trees = get_annotation(args)
        
    # Check if the cadd file is compressed and indexed:
    
    if args.cadd_file[0]:
        check_file_existence(args.cadd_file[0])
        if args.verbose:
            print('Cadd file! %s' % args.cadd_file[0])
        try:
            tabix_index(args.cadd_file[0], seq_col=0, start_col=1, end_col=1, meta_char='#')
        except IOError as e:
            if args.verbose:
                print(e)
        
    if args.cadd_1000g[0]:
        check_file_existence(args.cadd_1000g[0])
        if args.verbose:
            print('Cadd 1000g file! %s' % args.cadd_1000g[0])
        try:
            tabix_index(args.cadd_1000g[0], seq_col=0, start_col=1, end_col=1, meta_char='#')
        except IOError as e:
            if args.verbose:
                print(e)
    
    if args.thousand_g[0]:
        check_file_existence(args.thousand_g[0])
        if args.verbose:
            print('1000g frequency file! %s' % args.thousand_g[0])
        try:
            tabix_index(args.thousand_g[0], preset='vcf')
        except IOError as e:
            if args.verbose:
                print(e)
    
    # # Check the variants:
    
    my_vcf_parser = vcf_parser.VCFParser(var_file)
    head = my_vcf_parser.metadata
    
    if set(my_family.individuals.keys()) != set(my_vcf_parser.individuals):
        
        print('There must be same individuals in ped file and vcf file! Aborting...')
        print('Individuals in PED file: %s' % '\t'.join(list(my_family.individuals.keys())))
        print('Individuals in VCF file: %s' % '\t'.join(list(my_vcf_parser.individuals)))
        sys.exit()
        
    # The task queue is where all jobs(in this case batches that represents variants in a region) is put
    # the consumers will then pick their jobs from this queue.
    variant_queue = JoinableQueue(maxsize=1000)
    # The consumers will put their results in the results queue
    results = Manager().Queue()
    
    # Create a directory to keep track of temp files
    temp_dir = mkdtemp()
        
    num_model_checkers = (cpu_count()*2-1)
    
    if args.verbose:
        print('Number of CPU:s %s' % cpu_count())
    
    # These are the workers that do the analysis
    model_checkers = [variant_consumer.VariantConsumer(variant_queue, results, my_family, args) for i in range(num_model_checkers)]
    
    for w in model_checkers:
        w.start()
    
    # This process prints the variants to temporary files
    var_printer = variant_printer.VariantPrinter(results, temp_dir, head, args.verbose)
    var_printer.start()
    
    if args.verbose:
        print('Start parsing the variants ...')
        print('')
        start_time_variant_parsing = datetime.now()    
    
    # For parsing the vcf:
    var_annotator = variant_annotator.VariantAnnotator(my_vcf_parser, variant_queue, args, gene_trees, exon_trees)
    var_annotator.annotate()
    
    for i in range(num_model_checkers):
        variant_queue.put(None)
    
    variant_queue.join()
    results.put(None)
    var_printer.join()
    
    chromosome_list = var_annotator.chromosomes
        
    if args.verbose:
        print('Cromosomes found in variant file: %s' % ','.join(chromosome_list))
        print('Models checked!')
        print('Start sorting the variants:')
        print('')
        start_time_variant_sorting = datetime.now()
    
    # Add the new metadata to the headers:
    add_metadata(head, args)
    print_headers(head, args)
    
    for chromosome in chromosome_list:
        for temp_file in os.listdir(temp_dir):
            if temp_file.split('_')[0] == chromosome:
                var_sorter = variant_sorter.FileSort(os.path.join(temp_dir, temp_file), 
                                                outFile=args.outfile[0], silent=args.silent)
                var_sorter.sort()
    
    if args.verbose:
        print('Sorting done!')
        print('Time for sorting: %s' % str(datetime.now()-start_time_variant_sorting))
        print('')
        print('Time for whole analyis: %s' % str(datetime.now() - start_time_analysis))
    
    # Remove all temp files:
    shutil.rmtree(temp_dir)
    

if __name__ == '__main__':
    main()

