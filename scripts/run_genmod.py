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
from datetime import datetime
from tempfile import NamedTemporaryFile
import pkg_resources

from pprint import pprint as pp

# pp(sys.path)

from ped_parser import parser

from genmod.utils import is_number, variant_consumer, variant_sorter
from genmod.variants import variant_parser
from genmod.vcf import vcf_header

def get_family(args):
    """Return the family"""
    family_type = 'ped'
    family_file = args.family_file[0]
    
    my_family_parser = parser.FamilyParser(family_file, family_type)
    # Stupid thing but for now when we only look at one family
    return my_family_parser.families.popitem()[1]

def get_header(variant_file):
    """Return a fixed header parser"""
    head = vcf_header.VcfParser(variant_file)
    return head

def print_headers(args, header_object):
    """Print the headers to a results file."""
    header_object.header.append('Inheritance_model')
    header_object.header.append('Compounds')
    header_object.header.append('Rank_score')
    if args.outfile[0]:
        with open(args.outfile[0], 'w') as f: 
            for head_count in header_object.metadata:
                f.write(header_object.metadata[head_count]+'\n')
            f.write('#' + '\t'.join(header_object.header) + '\n')
    else:
        if not args.silent:
            for head_count in header_object.metadata:
                print header_object.metadata[head_count]
            print '#' + '\t'.join(header_object.header)
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
    
    parser.add_argument('--version', 
        action="version", 
        version=pkg_resources.require("genmod")[0].version)
    
    
    parser.add_argument('-v', '--verbose', 
        action="store_true", 
        help='Increase output verbosity.'
    )
    
    
    parser.add_argument('-o', '--output', 
        type=str, nargs=1, 
        help='Specify the path to a file where results should be stored.'
    )
    
    
    args = parser.parse_args()
    
    start_time_analysis = datetime.now()
    
            
    # Start by parsing at the pedigree file:

    my_family = get_family(args)
    
    head = get_header(var_file)        
    # # Check the variants:
    
    if args.verbose:
        print 'Parsing variants ...'
        print ''
    
    
    start_time_variant_parsing = datetime.now()
        
    
    var_file = args.variant_file[0]
    file_name, file_extension = os.path.splitext(var_file)
    
    individuals = head.individuals

    # # The task queue is where all jobs(in this case batches that represents variants in a region) is put
    # # the consumers will then pick their jobs from this queue.
    variant_queue = JoinableQueue()
    # # The consumers will put their results in the results queue
    # results = Manager().Queue()
    # # We will need a lock so that the consumers can print their results to screen
    # lock = Lock()
    
    # num_consumers = cpu_count() * 2
    # consumers = [variant_consumer.VariantConsumer(lock, tasks, results, my_family, args.verbose, file_handle) for i in xrange(num_consumers)]
    # 
    # for w in consumers:
    #     w.start()
    
    var_parser = variant_parser.VariantParser(var_file, variant_queue, head, args.verbose)
    
    # for i in xrange(num_consumers):
    #     tasks.put(None)
    
    # tasks.join()
    # file_handle.close()
    # temp_file.seek(0)
    
    # var_sorter = variant_sorter.FileSort(temp_file)
    # var_sorter.sort()
    
    # for line in temp_file:
    #     print line
        
    # os.remove(temp_file)    
    # # 
    # # while num_jobs:
    # #     result = results.get()
    # #     print 'Result: ', result
    # #     num_jobs -= 1
    # 
    # 
    if args.verbose:
        print 'Variants done!. Time to parse variants: ', (datetime.now() - start_time_variant_parsing)
        print ''
    # # 
    # # # Add info about variant file:
    # # new_headers = my_variant_parser.header_lines 
    # # 
    # # # Add new headers:
    # # 
    # # new_headers.append('Inheritance_model')
    # # new_headers.append('Compounds')
    # # new_headers.append('Rank_score')
    # # 
    # # 
    # # if args.verbose:
    # #     print 'Checking genetic models...'
    # #     print ''
    # # 
    # # for data in my_variant_parser.metadata:
    # #     print data
    # # 
    # # print '#'+'\t'.join(new_headers)
    # # 
    # # if not args.position:
    # #     all_variants = {}
    # # 
    # # # Check the genetic models
    # # 
    # # jobs=[]
    # # 
    # # for chrom in my_variant_parser.chrom_shelves:
    # #     
    # #     shelve_directory = os.path.split(my_variant_parser.chrom_shelves[chrom])[0]
    # #     current_shelve = my_variant_parser.chrom_shelves[chrom]
    # #     
    # #     p = multiprocessing.Process(target=check_variants, args=(current_shelve, my_family, gene_annotation, args, preferred_models))
    # #     jobs.append(p)
    # #     p.start()
    # # 
    # # for job in jobs:
    # #     job.join()
    # # 
    # # # Print all variants:
    # # 
    # # for chrom in my_variant_parser.chrom_shelves:
    # #     
    # #     variants = []        
    # #     variant_db = shelve.open(my_variant_parser.chrom_shelves[chrom])
    # #     
    # #     for var_id in variant_db:
    # #         variants.append(variant_db[var_id])
    # #         
    # #     for variant in sorted(variants, key=lambda genetic_variant:genetic_variant.start):
    # #         pass
    # #         # print '\t'.join(variant.get_cmms_variant())
    # # 
    # # 
    # #     os.remove(my_variant_parser.chrom_shelves[chrom])
    # # os.removedirs(shelve_directory)
    # # 
    # # # Else print by rank score:
    # # # if not args.position:
    # # #     for variant in sorted(all_variants.iteritems(), key=lambda (k,v): int(operator.itemgetter(-1)(v)), reverse=True):
    # # #         if args.treshold:
    # # #             rank_score = int(variant[-1][-1])
    # # #             if rank_score >= args.treshold[0]:
    # # #                 print '\t'.join(variant[1])
    # # #         else:
    # #             # print '\t'.join(variant[1])
    # # if args.verbose:
    # #     print 'Finished analysis!'
    # #     print 'Time for analysis', (datetime.now() - start_time_variant_parsing)
    pass


if __name__ == '__main__':
    main()

