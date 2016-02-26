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
from genmod.annotate_variants import (get_frequencies, get_cadd_scores, 
get_spidex_score, get_cosmic)

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

        self.max_af = kwargs.get('max_af', None)
        self.logger.debug("Setting max_af to {0}".format(self.max_af))

        self.spidex = kwargs.get('spidex', None)
        self.logger.debug("Setting spidex to {0}".format(self.spidex))

        self.cosmic = kwargs.get('cosmic', None)
        self.logger.debug("Setting cosmic to {0}".format(self.spidex))

        self.exon_trees = kwargs.get('exon_trees', {})
        self.gene_trees = kwargs.get('gene_trees', {})

        ######### Annotation files #########
        # Cadd files #
        self.cadd_handles = []
        for cadd_path in kwargs.get('cadd_files', []):
            self.logger.debug("Opening cadd file {0} with tabix open".format(
                cadd_path
            ))
            self.cadd_handles.append(tabix.open(cadd_path))
            self.logger.debug("Cadd file opened")
            
        # Frequency files #
        self.thousand_g = kwargs.get('thousand_g', None)
        self.logger.debug("1000G frequency file {0}".format(self.thousand_g))
        self.exac = kwargs.get('exac', None)
        self.logger.debug("Exac frequency file {0}".format(self.thousand_g))
        # self.dbNSFP = kwargs.get('dbNSFP', None)
        
        if self.thousand_g:
            self.logger.debug("Opening 1000G frequency file with tabix open")
            self.thousand_g = tabix.open(self.thousand_g)
            self.logger.debug("1000G frequency file opened")
        if self.exac:
            self.logger.debug("Opening ExAC frequency file with tabix open")
            self.exac = tabix.open(self.exac)
            self.logger.debug("ExAC frequency file opened")
        if self.spidex:
            self.logger.debug("Opening Spidex database with tabix open")
            self.spidex = tabix.open(self.spidex)
            self.logger.debug("Spidex file opened")
        if self.cosmic:
            self.logger.debug("Opening COSMIC database with tabix open")
            self.cosmic = tabix.open(self.cosmic)
            self.logger.debug("COSMIC file opened")

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
            
            chrom = variant_dict['CHROM'].lstrip('chr')
            position = int(variant_dict['POS'])
            ref = variant_dict['REF']
            alternative = variant_dict['ALT']

            longest_alt = max([
                len(alt) for alt in alternative.split(',')])

            if self.exon_trees:
                exonic = check_exonic(
                    chrom = chrom, 
                    start = position,
                    stop = (position+longest_alt)-1, 
                    exon_trees = self.exon_trees
                )
                if exonic:
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
                
                frequencies = get_frequencies(
                    tabix_reader = self.thousand_g,
                    chrom=chrom, 
                    start=position, 
                    alt=alternative
                )
                if frequencies.get('AF'):    
                    variant_dict = add_vcf_info(
                        keyword = "1000GAF",
                        variant_dict = variant_dict,
                        annotation = frequencies['AF']
                    )
                if self.max_af:
                    if frequencies.get('MAX_AF'):
                        variant_dict = add_vcf_info(
                            keyword = "1000G_MAX_AF",
                            variant_dict = variant_dict,
                            annotation = frequencies['MAX_AF']
                            )

            if self.exac:
                
                frequencies = get_frequencies(
                    tabix_reader = self.exac,
                    chrom=chrom, 
                    start=position, 
                    alt=alternative
                )
                if frequencies.get('AF'):    
                    variant_dict = add_vcf_info(
                        keyword = "ExACAF",
                        variant_dict = variant_dict,
                        annotation = frequency
                    )
                if self.max_af:
                    if frequencies.get('MAX_AF'):
                        variant_dict = add_vcf_info(
                            keyword = "ExAC_MAX_AF",
                            variant_dict = variant_dict,
                            annotation = frequency
                        )

            if self.spidex:
                spidex_score = get_spidex_score(
                        tabix_reader=self.spidex, 
                        chrom=chrom, 
                        start=position, 
                        alt=alternative
                    )
                if spidex_score:
                    variant_dict = add_vcf_info(
                        keyword = "SPIDEX",
                        variant_dict = variant_dict,
                        annotation = spidex_score
                    )

            if self.cosmic:
                in_cosmic = get_cosmic(
                        tabix_reader=self.cosmic, 
                        chrom=chrom, 
                        start=position, 
                        alt=alternative
                    )
                if in_cosmic:
                    variant_dict = add_vcf_info(
                        keyword = "COSMIC",
                        variant_dict = variant_dict,
                    )
                
            if self.cadd_handles:
                
                cadd_phred = None
                cadd_raw = None
                for cadd_handle in self.cadd_handles:
                    
                    cadd_scores = get_cadd_scores(
                        tabix_reader=cadd_handle, 
                        chrom=chrom, 
                        start=position,
                        alt=alternative 
                    )                        
                    if cadd_scores['cadd_phred']:
                        cadd_phred = cadd_scores['cadd_phred']
                        cadd_raw = cadd_scores['cadd_raw']

                        variant_dict = add_vcf_info(
                            keyword = "CADD",
                            variant_dict = variant_dict,
                            annotation = cadd_phred
                        )

                        if self.cadd_raw and cadd_raw:
                            variant_dict = add_vcf_info(
                                keyword = "CADD_raw",
                                variant_dict = variant_dict,
                                annotation = cadd_raw
                            )
                        break
                    

            self.results_queue.put(variant_dict)
            self.task_queue.task_done()

        return



def main():
    pass

if __name__ == '__main__':
    main()