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
from multiprocessing import JoinableQueue, Manager, cpu_count, Lock
from datetime import datetime
from tempfile import mkdtemp
import shutil
import pkg_resources
if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict


from pprint import pprint as pp

# pp(sys.path)

from ped_parser import parser

from genmod.utils import is_number, variant_consumer, variant_sorter, annotation_parser, variant_printer
from genmod.vcf import vcf_header, vcf_parser

def get_family(args):
    """Return the family"""
    family_type = 'ped'
    family_file = args.family_file[0]
    
    my_family_parser = parser.FamilyParser(family_file, family_type)
    # Stupid thing but for now when we only look at one family
    return my_family_parser.families.popitem()[1]

def get_header(variant_file):
    """Return a fixed header parser"""
    head = vcf_header.VCFParser(variant_file)
    head.parse()
    return head

def add_metadata(head):
    """Add metadata for the information added by this script."""
    head.metadataparser.add_info('ANN', '.', 'String', 'Annotates what feature(s) this variant belongs to.')
    head.metadataparser.add_info('Comp', '.', 'String', "':'-separated list of compound pairs for this variant.")
    head.metadataparser.add_info('GM', '.', 'String', "':'-separated list of genetic models for this variant.")
    return

def print_headers(args, header_object):
    """Print the headers to a results file."""
    if args.outfile[0]:
        with open(args.outfile[0], 'w') as f: 
            for head_count in header_object.print_header():
                f.write(head_count+'\n')
    else:
        if not args.silent:
            for line in header_object.print_header():
                print line
    return

def main():
    parser = argparse.ArgumentParser(description="Annotate genetic models in variant files..")
    
    parser.add_argument('family_file', 
        type=str, nargs=1, 
        help='A pedigree file in .ped format.'
    )
    parser.add_argument('variant_file', 
        type=str, nargs=1, 
        help='A variant file. Default is vcf format.'
    )

    parser.add_argument('annotation_file', 
        type=str, nargs=1, 
        help='A annotations file. Default is ref_gene format.'
    )
    
    parser.add_argument('-at', '--annotation_type',  
        type=str, nargs=1, choices=['bed', 'ccds', 'gtf', 'ref_gene'],
        default=['ref_gene'], help='Specify the format of the annotation file.'
    )    
    
    parser.add_argument('--version', 
        action="version", 
        version=pkg_resources.require("genmod")[0].version)
    
    
    parser.add_argument('-v', '--verbose', 
        action="store_true", 
        help='Increase output verbosity.'
    )
    
    parser.add_argument('-s', '--silent', 
        action="store_true", 
        help='Do not print the variants.'
    )
    
    parser.add_argument('-o', '--outfile', 
        type=str, nargs=1, default=[None],
        help='Specify the path to a file where results should be stored.'
    )
    
    
    args = parser.parse_args()
    var_file = args.variant_file[0]
    file_name, file_extension = os.path.splitext(var_file)
    anno_file = args.annotation_file[0]
        
    start_time_analysis = datetime.now()
    
            
    # Start by parsing at the pedigree file:

    my_family = get_family(args)
    
    # Parse the header of the vcf:
    
    head = get_header(var_file)
    add_metadata(head)
    # Parse the annotation file and make annotation trees:

    if args.verbose:
        print 'Parsing annotation ...'
        print ''
        start_time_annotation = datetime.now()
    
    annotation_trees = annotation_parser.AnnotationParser(anno_file, args.annotation_type[0])
            
    # # Check the variants:
    
    if args.verbose:
        print 'Annotation Parsed!'
        print 'Time to parse annotation:', datetime.now() - start_time_annotation
        print ''
        print 'Parsing variants ...'
        print ''
    
    start_time_variant_parsing = datetime.now()
        
    # The task queue is where all jobs(in this case batches that represents variants in a region) is put
    # the consumers will then pick their jobs from this queue.
    variant_queue = JoinableQueue()
    # The consumers will put their results in the results queue
    results = Manager().Queue()
    
    temp_dir = mkdtemp()
    # Create a temporary file for the variants:
        
    num_model_checkers = (cpu_count()*2-1)
    
    model_checkers = [variant_consumer.VariantConsumer(variant_queue, results, my_family, 
                     args.verbose) for i in xrange(num_model_checkers)]
    
    for w in model_checkers:
        w.start()
    
    var_printer = variant_printer.VariantPrinter(results, temp_dir, args.verbose)
    var_printer.start()
    
    if args.verbose:
        print 'Start parsing the variants ...'
        print ''
        start_time_variant_parsing = datetime.now()    
    
    var_parser = vcf_parser.VariantFileParser(var_file, variant_queue, head, annotation_trees, args.verbose)
    var_parser.parse()
    
    for i in xrange(num_model_checkers):
        variant_queue.put(None)
    
    variant_queue.join()
    results.put(None)
    var_printer.join()
    
    chromosome_list = var_parser.chromosomes
        
    if args.verbose:
        print 'Cromosomes in variant file:', ','.join(chromosome_list)
        print 'Models checked!'
        print 'Start sorting the variants:'
        print ''
        start_time_variant_sorting = datetime.now()
    
    print_headers(args, head)
    
    for chromosome in chromosome_list:
        for temp_file in os.listdir(temp_dir):
            if temp_file.split('_')[0] == chromosome:
                var_sorter = variant_sorter.FileSort(os.path.join(temp_dir, temp_file), outFile=args.outfile[0], silent=args.silent)
                var_sorter.sort()
    
    if args.verbose:
        print 'Sorting done!'
        print 'Time for sorting:', datetime.now()-start_time_variant_sorting
        print ''
        print 'Time for analyis:', datetime.now() - start_time_analysis
    
    shutil.rmtree(temp_dir)
    

if __name__ == '__main__':
    main()

