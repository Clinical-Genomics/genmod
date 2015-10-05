#!/usr/bin/env python
# encoding: utf-8
"""
variant_consumer.py

Consumes batches of variants and annotates them. Each batch is a dictionary 
with variant_id:s as keys and dictionaries with variant information.
The variants will get different annotations depending on input
 
Created by Måns Magnusson on 2013-03-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import (division, print_function, absolute_import)

import sys
import os
import operator
import logging

from pprint import pprint as pp
from multiprocessing import Process
from math import log10

from genmod.utils import check_exonic
from genmod.vcf_tools import get_genotypes
from . import (get_haploblocks, check_genetic_models, get_model_score, 
make_print_version)
                

class VariantAnnotator(Process):
    """
    Annotates variant in batches from the task queue and puts the result in 
    the results queue.
    """
    
    def __init__(self, task_queue, results_queue, families, individuals, 
                phased=False, strict=False, whole_gene=False, vep=False,
                reduced_penetrance_genes = set()):
        """
        Initialize the VariantAnnotator
        
        Consume variant batches from the task queue, annotate them with the 
        genetic inheritance patterns that they follow and put them in the 
        results queue.
        
        Arguments:
            task_queue (Queue)
            results_queue (Queue)
            families (dict)
            individuals (list)
            phased (bool)
            strict (bool)
            whole_gene (bool)
            vep (bool)
            reduced_penetrance_genes (set): Set of reduced penetrance genes
        """
        Process.__init__(self)
        self.logger = logging.getLogger(__name__)
        
        self.proc_name = self.name
        
        self.logger.info("Setting up variant_annotator: {0}".format(
            self.proc_name))
        
        self.logger.debug("Setting up task queue")
        self.task_queue = task_queue
        
        self.logger.debug("Setting up results queue")
        self.results_queue = results_queue
        
        # The families that should be annotated
        self.families = families
        self.logger.debug("Families found: {0}".format(self.families))
        self.individuals = individuals
        self.logger.debug("Individuals found: {0}".format(self.individuals))
        
        # Settings for the annotation
        self.phased = phased
        self.logger.debug("Setting phased to {0}".format(self.phased))
        self.strict = strict
        self.logger.debug("Setting strict to {0}".format(self.strict))
        self.vep = vep
        self.logger.debug("Setting vep to {0}".format(self.vep))
        self.whole_gene = whole_gene
        self.logger.debug("Setting whole_gene to {0}".format(self.whole_gene))
        self.reduced_penetrance = reduced_penetrance_genes
        

    def run(self):
        """Run the consuming"""
        self.logger.info('%s: Starting!' % self.proc_name)
        # Check if there are any batches in the queue
        while True:
            # A batch is a dictionary with varints on the form {variant_id:variant_dict}
            self.logger.debug("Getting task from task_queue")
            variant_batch = self.task_queue.get()
            
            if variant_batch is None:
                self.logger.info('No more batches')
                self.task_queue.task_done()
                self.logger.info('{0}: Exiting'.format(self.proc_name))
                break
                
            # We are now going to check the genetic models for the variants in
            # the batch
            
            for variant_id in variant_batch:
                variant = variant_batch[variant_id]
                variant['genotypes'] = get_genotypes(variant, self.individuals)
                
                # Check if the variant is in a gene with reduced penetrance
                if variant.get('annotation', set()).intersection(self.reduced_penetrance):
                    self.logger.debug("Setting reduced_penetrance to True for"\
                    " variant: {0}".format(variant_id))
                    
                    variant['reduced_penetrance'] = True
            
            if len(variant_batch) > 1:
                #If the variant are phased we want to find out which 
                #haploblocks they belong to for compounds
                if self.phased:
                    self.logger.debug("Get haploblocks for variant batch")
                    variant_batch['haploblocks'] = get_haploblocks(
                        variant_batch, self.individuals
                    )
                
                # We only need to check compound candidates if there is 
                # more than one variant in the batch
                for variant_id in variant_batch:
                    self.logger.debug("Check compound candidates")
                    variant = variant_batch[variant_id]
                    
                    variant['compound_candidate'] = False
                    
                    if self.whole_gene:
                        if variant['annotation']:
                            variant['compound_candidate'] = True
                            self.logger.debug("Set compound_candidate to True")
                    elif check_exonic(variant, self.vep):
                            variant['compound_candidate'] = True
                            self.logger.debug("Set compound_candidate to True")
                    

            # Check the genetic models for all variants in the batch
            check_genetic_models(
                variant_batch = variant_batch,
                families = self.families,
                phased = self.phased,
                strict = self.strict,
                )

            # # Now we want to make versions of the variants that are ready for printing.
            for variant_id in variant_batch:
                variant = make_print_version(
                    variant=variant_batch[variant_id],
                    families=self.families
                )
                self.logger.debug("Putting variant in results_queue")
                self.results_queue.put(variant)
                
            self.task_queue.task_done()
        
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()