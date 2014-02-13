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
            for feature in variant_batch:
                #Make one dictionary for each feature:
                variant_dict = dict((variant_id, variant_info) for variant_id, variant_info in variant_batch[feature].items())
                for variant_id in variant_dict:
                    #Remove the 'Genotypes' post since we will not need them for now
                    variant_dict[variant_id].pop('Genotypes', 0)
                    if variant_id in fixed_variants:
                        if len(variant_dict[variant_id]['Compounds']) > 0:
                            fixed_variants[variant_id]['Compounds'] = dict(variant_dict[variant_id]['Compounds'].items() +
                                                                    fixed_variants[variant_id]['Compounds'].items())
                            fixed_variants[variant_id]['Inheritance_model']['AR_compound'] = True
                    else:
                        fixed_variants[variant_id] = variant_dict[variant_id]
                        
            # Now we want to make versions of the variants that are ready for printing.
            for variant_id in fixed_variants:
                model_list = []
                compounds_list = []
                feature_list = fixed_variants[variant_id]['Annotation']
                if len(fixed_variants[variant_id]['Compounds']) > 0:
                    #We do not want reference to itself as a compound:
                    fixed_variants[variant_id]['Compounds'].pop(variant_id, 0)
                    compounds_list = fixed_variants[variant_id]['Compounds'].keys()
                else:
                    compounds_list = ['-']
                for model in fixed_variants[variant_id]['Inheritance_model']:
                    if fixed_variants[variant_id]['Inheritance_model'][model]:
                        model_list.append(model)
                if len(model_list) == 0:
                    model_list = ['NA']
                fixed_variants[variant_id].pop('Compounds',0)
                fixed_variants[variant_id].pop('Inheritance_model',0)
                fixed_variants[variant_id].pop('Annotation',0)
                vcf_info = fixed_variants[variant_id]['INFO'].split(';')
                # if we should include the annotation:
                vcf_info.append('ANN=' + ':'.join(feature_list))
                # if we should include compounds:
                vcf_info.append('Comp=' + ':'.join(compounds_list))
                # if we should include genetic models:
                vcf_info.append('GM=' + ':'.join(model_list))
                fixed_variants[variant_id]['INFO'] = ';'.join(vcf_info)
                # pp(fixed_variants[variant_id])
            self.results_queue.put(fixed_variants)
            self.task_queue.task_done()
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()