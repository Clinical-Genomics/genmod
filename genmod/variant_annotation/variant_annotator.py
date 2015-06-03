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

from __future__ import (division, print_function, unicode_literals, absolute_import)

import sys
import os
import operator
import tabix
import logging

from pprint import pprint as pp
from multiprocessing import Process
from math import log10

from genmod.utils import check_exonic
from . import (get_model_score, annotate_cadd_score, 
                get_haploblocks, annotate_frequency, check_genetic_models,
                make_print_version)

class VariantAnnotator(Process):
    """
    Annotates variant in batches from the task queue and puts the result in 
    the results queue.
    """
    
    def __init__(self, task_queue, results_queue, **kwargs):
        
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
        self.families = kwargs.get('families', {})
        self.logger.debug("Families found: {0}".format(self.families))
        self.individuals = kwargs.get('individuals', [])
        self.logger.debug("Individuals found: {0}".format(self.individuals))
        # Settings for the annotation
        self.phased = kwargs.get('phased', False)
        self.logger.debug("Setting phased to {0}".format(self.phased))
        self.strict = kwargs.get('strict', False)
        self.logger.debug("Setting strict to {0}".format(self.strict))
        self.cadd_raw = kwargs.get('cadd_raw', None)
        self.logger.debug("Setting cadd raw to {0}".format(self.cadd_raw))
        
        self.exon_trees = kwargs.get('exon_trees', {})
        self.whole_gene = kwargs.get('whole_gene', False)
        ######### Annotation files #########
        # Cadd files #
        self.cadd_file = kwargs.get('cadd_file', None)
        self.logger.debug("Cadd file {0}".format(self.cadd_file))
        self.cadd_1000g = kwargs.get('cadd_1000g', None)
        self.logger.debug("Cadd 1000G file {0}".format(self.cadd_1000g))
        self.cadd_exac = kwargs.get('cadd_exac', None)
        self.logger.debug("Cadd exac file {0}".format(self.cadd_exac))
        self.cadd_ESP = kwargs.get('cadd_ESP', None)
        self.logger.debug("Cadd ESP file {0}".format(self.cadd_ESP))
        self.cadd_InDels = kwargs.get('cadd_InDels', None)
        self.logger.debug("Cadd InDels file {0}".format(self.cadd_InDels))
        # Frequency files #
        self.thousand_g = kwargs.get('thousand_g', None)
        self.logger.debug("1000G frequency file {0}".format(self.thousand_g))
        self.exac = kwargs.get('exac', None)
        self.logger.debug("Exac frequency file {0}".format(self.thousand_g))
        self.dbNSFP = kwargs.get('dbNSFP', None)
        self.any_cadd_info = False
        # Setup file handles to the annotation files
        if self.cadd_file:
            self.logger.debug("Opening cadd file with tabix open")
            self.cadd_file = tabix.open(self.cadd_file)
            self.logger.debug("Cadd file opened")
            self.any_cadd_info = True
        if self.cadd_1000g:
            self.logger.debug("Opening cadd 1000G file with tabix open")
            self.cadd_1000g = tabix.open(self.cadd_1000g)
            self.logger.debug("Cadd 1000G file opened")
            self.any_cadd_info = True
        if self.cadd_exac:
            self.logger.debug("Opening cadd exac file with tabix open")
            self.cadd_exac = tabix.open(self.cadd_exac)
            self.logger.debug("Cadd exac file opened")
            self.any_cadd_info = True
        if self.cadd_ESP:
            self.logger.debug("Opening cadd ESP file with tabix open")
            self.cadd_ESP = tabix.open(self.cadd_ESP)
            self.logger.debug("Cadd ESP file opened")
            self.any_cadd_info = True
        if self.cadd_InDels:
            self.logger.debug("Opening cadd InDels file with tabix open")
            self.cadd_InDels = tabix.open(self.cadd_InDels)
            self.logger.debug("Cadd InDels file opened")
            self.any_cadd_info = True
        if self.thousand_g:
            self.logger.debug("Opening 1000G frequency file with tabix open")
            self.thousand_g = tabix.open(self.thousand_g)
            self.logger.debug("1000G frequency file opened")
        if self.exac:
            self.logger.debug("Opening ExAC frequency file with tabix open")
            self.exac = tabix.open(self.exac)
            self.logger.debug("ExAC frequency file opened")
        if self.dbNSFP:
            self.dbNSFP = tabix.open(self.exac)
        self.logger.debug("Setting any cadd info to {0}".format(self.any_cadd_info))

    def run(self):
        """Run the consuming"""
        self.logger.info('%s: Starting!' % self.proc_name)
        # Check if there are any batches in the queue
        while True:
            # A batch is a dictionary on the form {variant_id:variant_dict}
            variant_batch = self.task_queue.get()
            
            if variant_batch is None:
                self.logger.info('No more batches')
                self.task_queue.task_done()
                self.logger.info('{0}: Exiting'.format(self.proc_name))
                break
                
            # If there are families we will annotate the genetic inheritance 
            # patterns that are followed
            if self.families:
                if len(variant_batch) > 1:
                    self.logger.debug("Get haploblocks for variant batch")
                    variant_batch['haploblocks'] = get_haploblocks(variant_batch, self.individuals)
                
                for variant_id in variant_batch:
                    variant = variant_batch[variant_id]
                    variant['compound_candidate'] = False
                    if variant['annotation']:
                        if self.whole_gene:
                            variant['compound_candidate'] = True
                        else:
                            variant['compound_candidate'] = check_exonic(variant)
                    variant_batch[variant_id] = variant
                
                check_genetic_models(
                                variant_batch = variant_batch, 
                                families = self.families, 
                                phased = self.phased, 
                                strict = self.strict, 
                            )

            # These are family independent annotations which will be done annyway:
            for variant_id in variant_batch:
                variant = variant_batch[variant_id]
                
                if self.any_cadd_info:
                    variant = annotate_cadd_score(
                        variant,
                        self.cadd_raw,
                        self.cadd_file, 
                        self.cadd_1000g, 
                        self.cadd_exac, 
                        self.cadd_ESP, 
                        self.cadd_InDels
                        )
                
                if self.thousand_g or self.exac:
                    variant = annotate_frequency(
                        variant,
                        self.thousand_g,
                        self.exac
                        )
                
                # Now we want to make versions of the variants that are ready for printing.
                variant = make_print_version(
                    variant=variant,
                    families=self.families
                )
                
                variant_batch[variant_id] = variant
            
            
            self.results_queue.put(variant_batch)
            self.task_queue.task_done()
        
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()