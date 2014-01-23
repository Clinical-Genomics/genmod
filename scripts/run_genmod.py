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
from multiprocessing import JoinableQueue, Queue, Lock, cpu_count
from datetime import datetime

from pprint import pprint as pp

# pp(sys.path)

from genmod.utils import is_number, variant_consumer
from genmod.family import family_parser
from genmod.variants import variant_parser



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
    
    parser.add_argument('-var_type', '--variant_type', 
        type=str, choices=['CMMS', 'VCF'], nargs=1, default=['CMMS'],
        help='Specify the format on the variant file.'
    )
    
    parser.add_argument('-v', '--verbose', 
        action="store_true", 
        help='Increase output verbosity.'
    )
    
    parser.add_argument('-ga', '--gene_annotation', 
        type=str, choices=['Ensembl', 'HGNC'], nargs=1, default=['HGNC'],
        help='What gene annotation should be used, HGNC or Ensembl.'
    )
    
    parser.add_argument('-o', '--output', 
        type=str, nargs=1, 
        help='Specify the path to a file where results should be stored.'
    )
    
    # parser.add_argument('-pos', '--position', 
    #     action="store_true", 
    #     help='If output should be sorted by position. Default is sorted on rank score'
    # )
    
    parser.add_argument('-tres', '--treshold', 
        type=int, nargs=1,  
        help='Specify the lowest rank score to be outputted.'
    )
    
    args = parser.parse_args()
    
    
    new_headers = []    
        
    # Start by parsing at the pedigree file:
    family_type = 'CMMS'    
    my_family_parser = family_parser.FamilyParser(args.family_file[0], family_type)
    
    print my_family_parser
    # # Stupid thing but for now when we only look at one family
    my_family = my_family_parser.families.popitem()[1]
    
    preferred_models = my_family.models_of_inheritance
        
    # # Check the variants:
    
    if args.verbose:
        print 'Parsing variants ...'
        print ''
    
    
    start_time_variant_parsing = datetime.now()
    
    var_file = args.variant_file[0]
    file_name, file_extension = os.path.splitext(var_file)
    
    individuals = [ind.individual_id for ind in my_family.individuals]

    var_type = 'cmms'        
    # header_line = []
    # metadata = []
    # 
    # # The task queue is where all jobs(in this case batches that represents variants in a region) is put
    # # the consumers will then pick their jobs from this queue.
    tasks = JoinableQueue()
    # # The consumers will put their results in the results queue
    results = Queue()
    # # We will need a lock so that the consumers can print their results to screen
    lock = Lock()
    
    num_consumers = cpu_count() * 2
    consumers = [variant_consumer.VariantConsumer(lock, tasks, results, my_family, args.verbose) for i in xrange(num_consumers)]
    
    for w in consumers:
        w.start()
    
    var_parser = variant_parser.VariantParser(var_file, tasks, individuals)
    
    for i in xrange(num_consumers):
        tasks.put(None)
    
    tasks.join()
    # # 
    # # while num_jobs:
    # #     result = results.get()
    # #     print 'Result: ', result
    # #     num_jobs -= 1
    # 
    # 
    # if args.verbose:
    #     print 'Variants done!. Time to parse variants: ', (datetime.now() - start_time_variant_parsing)
    #     print ''
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

