#!/usr/bin/env python
# encoding: utf-8
"""
variant_consumer.py

Class that takes a list of objects and return all unordered pairs as a generator.

If only one object? Raise Exception
 
Created by Måns Magnusson on 2013-03-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import multiprocessing
from pprint import pprint as pp

from genmod.models import genetic_models

class VariantConsumer(multiprocessing.Process):
    """Yeilds all unordered pairs from a list of objects as tuples, like (obj_1, obj_2)"""
    
    def __init__(self, task_queue, results_queue, family, verbosity = False):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.family = family
        self.results_queue = results_queue
        self.verbosity = verbosity
    
    def run(self):
        """Run the consuming"""
        proc_name = self.name
        if self.verbosity:
            print proc_name, 'Starting!'
        while True:
            # A batch is a dictionary on the form {gene:{variant_id:variant_dict}}
            next_batch = self.task_queue.get()
            if self.verbosity:
                # if self.task_queue.empty():
                #     print 'No variants to parse!', proc_name
                if self.results_queue.full():
                    print 'Batch results queue Full!', proc_name
                if self.task_queue.full():
                    print 'Variant queue full!', proc_name
            if next_batch is None:
                self.task_queue.task_done()
                if self.verbosity:
                    print '%s: Exiting' % proc_name
                break
            # print '%s: %s' % (proc_name, next_batch)
            variant_batch = genetic_models.check_genetic_models(next_batch, self.family, self.verbosity, proc_name)
            fixed_variants = {}
            # Make shore we only have one copy of each variant:
            for gene in variant_batch:
                #Make one dictionary for each gene:
                variant_dict = dict((variant_id, variant_info) for variant_id, variant_info in variant_batch[gene].items())
                for variant_id in variant_dict:
                    #Remove the 'Genotypes' post since we will not need them for now
                    variant_dict[variant_id].pop('Genotypes', 0)
                    if variant_id in fixed_variants:
                        if len(variant_dict[variant_id]['Compounds']) > 0:
                            for compound in variant_dict[variant_id]['Compounds']:
                                fixed_variants[variant_id]['Compounds'].append(compound)
                            fixed_variants[variant_id]['Inheritance_model']['AR_compound'] = True
                    else:
                        fixed_variants[variant_id] = variant_dict[variant_id]
            
            # fixed_variants = score_variants.score_variant(fixed_variants, self.family.models_of_inheritance)
            # # with self.lock:
            # #     for variant_id, variant in fixed_variants.items():
            # #         print variant_id
            # #         print variant
            # for variant_id in fixed_variants:
            #     model_list = []
            #     if len(fixed_variants[variant_id]['Compounds']) > 0:
            #         #We do not want reference to itself as a compound:
            #         fixed_variants[variant_id]['Compounds'].pop(variant_id, 0)
            #         compounds_list = []
            #         for compound_id in fixed_variants[variant_id]['Compounds']:
            #             compound_score = int(fixed_variants[variant_id]['Rank_score']) + int(fixed_variants[compound_id]['Rank_score'])
            #             fixed_variants[variant_id]['Compounds'][compound_id] = compound_score
            #             compounds_list.append(compound_id + '=' + str(compound_score))
            #         fixed_variants[variant_id]['Compounds'] = ':'.join(compounds_list)
            #     else:
            #         fixed_variants[variant_id]['Compounds'] = '-'
            #     for model in fixed_variants[variant_id]['Inheritance_model']:
            #         if fixed_variants[variant_id]['Inheritance_model'][model]:
            #             model_list.append(model)
            #     if len(model_list) == 0:
            #         fixed_variants[variant_id]['Inheritance_model'] = 'NA'
            #     else:
            #         fixed_variants[variant_id]['Inheritance_model'] = ':'.join(model_list)
            #     fixed_variants[variant_id]['Rank_score'] = str(fixed_variants[variant_id]['Rank_score'])
            #             
            #               
            # # print next_batch
            # self.results_queue.put(fixed_variants)
            self.task_queue.task_done()
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()