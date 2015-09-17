#!/usr/bin/env python
# encoding: utf-8
"""
variant_annotator.py

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
import tabix
import logging

from multiprocessing import Process
from math import log10

from genmod.vcf_tools import get_variant_dict, add_vcf_info
from genmod.annotate_regions import (get_genes, check_exonic)
from genmod.annotate_variants import annotate_frequency, get_cadd_scores

# from . import (annotate_cadd_score, annotate_frequency, check_genetic_models)

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
        
        self.header_line = kwargs.get('header_line')
        if not self.header_line:
            raise IOError("VariantAnnotator needs a vcf header line")

        self.cadd_raw = kwargs.get('cadd_raw', None)
        self.logger.debug("Setting cadd raw to {0}".format(self.cadd_raw))

        self.exon_trees = kwargs.get('exon_trees', {})
        self.gene_trees = kwargs.get('gene_trees', {})

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
        # self.dbNSFP = kwargs.get('dbNSFP', None)
        self.any_cadd_info = False

        # Setup file handles to the annotation files
        self.cadd_handles = []
        if self.cadd_file:
            self.logger.debug("Opening cadd file with tabix open")
            self.cadd_handles.append(tabix.open(self.cadd_file))
            self.logger.debug("Cadd file opened")
            self.any_cadd_info = True
        if self.cadd_1000g:
            self.logger.debug("Opening cadd 1000G file with tabix open")
            self.cadd_handles.append(tabix.open(self.cadd_1000g))
            self.logger.debug("Cadd 1000G file opened")
            self.any_cadd_info = True
        if self.cadd_exac:
            self.logger.debug("Opening cadd exac file with tabix open")
            self.cadd_handles.append(tabix.open(self.cadd_exac))
            self.logger.debug("Cadd exac file opened")
            self.any_cadd_info = True
        if self.cadd_ESP:
            self.logger.debug("Opening cadd ESP file with tabix open")
            self.cadd_handles.append(tabix.open(self.cadd_ESP))
            self.logger.debug("Cadd ESP file opened")
            self.any_cadd_info = True
        if self.cadd_InDels:
            self.logger.debug("Opening cadd InDels file with tabix open")
            self.cadd_handles.append(tabix.open(self.cadd_InDels))
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
        # if self.dbNSFP:
        #     self.dbNSFP = tabix.open(self.exac)
        self.logger.debug("Setting any cadd info to {0}".format(self.any_cadd_info))

    def run(self):
        """Run the consuming"""
        self.logger.info('%s: Starting!' % self.proc_name)
        # Check if there are any batches in the queue
        while True:
            # A batch is a dictionary on the form {variant_id:variant_dict}
            variant_line = self.task_queue.get()
            
            if variant_line is None:
                self.logger.info('No more variants')
                self.task_queue.task_done()
                self.logger.info('{0}: Exiting'.format(self.proc_name))
                break

            variant_dict = get_variant_dict(
                variant_line = variant_line, 
                header_line = self.header_line
            )
            
            chrom = variant_dict['CHROM']
            position = int(variant_dict['POS'])
            ref = variant_dict['REF']
            alternatives = variant_dict['ALT']

            longest_alt = max([
                len(alt) for alt in alternatives.split(',')])

            if self.exon_trees:
                if check_exonic( chrom = chrom, start = position,
                    stop = (position+longest_alt)-1, exon_trees = self.exon_trees):

                    variant_dict = add_vcf_info(
                            keyword = 'Exonic',
                            variant_dict=variant_dict,
                            annotation=None
                            )
            if self.gene_trees:
                genes = get_genes(
                    chrom = chrom,
                    start = position,
                    stop = (position+longest_alt)-1,
                    gene_trees = self.gene_trees
                )
                if genes:
                    variant_dict = add_vcf_info(
                        keyword = "Annotation",
                        variant_dict = variant_dict,
                        annotation = ','.join(genes)
                    )
            if self.thousand_g:
                
                frequencies = annotate_frequency(
                    chrom=chrom, 
                    position=position, 
                    alternatives=alternatives,
                    tabix_handle=self.thousand_g, 
                    
                )
                for frequency in frequencies:
                    if frequency != '0':    
                        variant_dict = add_vcf_info(
                            keyword = "1000GAF",
                            variant_dict = variant_dict,
                            annotation = ','.join(frequencies)
                        )

            if self.exac:
                
                frequencies = annotate_frequency(
                    chrom=chrom, 
                    position=position, 
                    alternatives=alternatives,
                    tabix_handle=self.exac
                    
                )
                for frequency in frequencies:
                    if frequency != '0':    
                        variant_dict = add_vcf_info(
                            keyword = "ExACAF",
                            variant_dict = variant_dict,
                            annotation = ','.join(frequencies)
                        )
            if self.any_cadd_info:
                cadd_phred = []
                cadd_raw = []
                for alternative in alternatives:
                    cadd_found = False
                    for cadd_handle in self.cadd_handles:
                        cadd_scores = get_cadd_scores(
                            tabix_reader=cadd_handle, 
                            chrom=chrom, 
                            start=position,
                            alt=alternative 
                        )                        
                        if cadd_scores['cadd_phred']:
                            cadd_found = True
                            cadd_phred.append(cadd_scores['cadd_phred'])
                            cadd_raw.append(cadd_scores['cadd_raw'])
                            break
                    if not cadd_found:
                        cadd_phred.append('0')
                        cadd_raw.append('0')
                for cadd_score in cadd_phred:
                    if cadd_score != '0':
                        variant_dict = add_vcf_info(
                            keyword = "CADD",
                            variant_dict = variant_dict,
                            annotation = ','.join(cadd_phred)
                        )
                
                if self.cadd_raw:
                    for cadd_score in cadd_raw:
                        if cadd_score != '0':
                            variant_dict = add_vcf_info(
                                keyword = "CADD_raw",
                                variant_dict = variant_dict,
                                annotation = ','.join(cadd_raw)
                            )
                    
            #     annotated_line = annotate_thousand_g(
            #         variant_line = line,
            #         thousand_g = thousand_g_handle
            #     )
            # if self.any_cadd_info:
            #     variant = annotate_cadd_score(
            #         variant=variant,
            #         cadd_raw=self.cadd_raw,
            #         cadd_file=self.cadd_file,
            #         cadd_1000g=self.cadd_1000g,
            #         cadd_exac=self.cadd_exac,
            #         cadd_ESP=self.cadd_ESP,
            #         cadd_InDels=self.cadd_InDels
            #         )
            #
            # if self.thousand_g or self.exac:
            #     variant = annotate_frequency(
            #         variant=variant,
            #         thousand_g=self.thousand_g,
            #         exac=self.exac
            #         )

            self.results_queue.put(variant_dict)
            self.task_queue.task_done()

        return



def main():
    pass

if __name__ == '__main__':
    main()