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

from __future__ import (division, print_function, unicode_literals, relative_import)

import sys
import os
import operator
import tabix
import logging

from pprint import pprint as pp
from multiprocessing import Process
from math import log10

from . import check_genetic_models
from genmod.errors import warning

from . import (get_model_score, get_frequency, get_cadd_scores)

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
        # Settings for the annotation
        self.phased = kwargs.get('phased', False)
        self.logger.debug("Setting phased to {0}".format(self.phased))
        self.strict = kwargs.get('strict', False)
        self.logger.debug("Setting strict to {0}".format(self.strict))
        self.cadd_raw = kwargs.get('cadd_raw', None)
        self.logger.debug("Setting cadd raw to {0}".format(self.cadd_raw))
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
            self.exac = tabix.open(self.exac)
        self.logger.debug("Setting any cadd info to {0}".format(self.any_cadd_info))
    
    def add_cadd_score(self, variant):
        """Add the CADD relative score to this variant."""
        cadd_score = None
        cadd_relative = None
        cadd_absolute = None
        cadd_relative_scores = []
        cadd_absolute_scores = []
        #Check CADD file(s):
        for alt in variant['ALT'].split(','):
            if self.cadd_file:
                cadd_scores = get_cadd_scores(
                                        self.cadd_file, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
                cadd_relative = cadd_score['cadd_phred']
                cadd_absolute = cadd_score['cadd_raw']
            # If variant not found in big CADD file check the 1000G file:
            if not (cadd_relative and cadd_absolute) and self.cadd_1000g:
                cadd_scores = get_cadd_scores(
                                        self.cadd_1000g, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
                cadd_relative = cadd_score['cadd_phred']
                cadd_absolute = cadd_score['cadd_raw']
            
            if not (cadd_relative and cadd_absolute) and self.cadd_exac:
                cadd_scores = get_cadd_scores(
                                        self.cadd_exac, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
                cadd_relative = cadd_score['cadd_phred']
                cadd_absolute = cadd_score['cadd_raw']
            
            if not (cadd_relative and cadd_absolute) and self.cadd_ESP:
                cadd_scores = get_cadd_scores(
                                        self.cadd_ESP, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
                cadd_relative = cadd_score['cadd_phred']
                cadd_absolute = cadd_score['cadd_raw']
            
            if not (cadd_relative and cadd_absolute) and self.cadd_InDels:
                cadd_scores = get_cadd_scores(
                                        self.cadd_ESP, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
                cadd_relative = cadd_score['cadd_phred']
                cadd_absolute = cadd_score['cadd_raw']
            
            if cadd_relative:
                cadd_relative_scores.append(str(cadd_relative))
            if cadd_absolute:
                cadd_absolute_scores.append(str(cadd_absolute))
        
        if len(cadd_relative_scores) > 0:
            variant['CADD'] = ','.join(cadd_relative_scores)
            if self.cadd_raw:
                variant['CADD_raw'] = ','.join(cadd_absolute_scores)
        return

    
    def add_frequency(self, variant):
        """Add the thousand genome frequency if present."""
        #Check 1000G frequency:
        thousand_g_freq = None
        exac_freq = None
        if self.thousand_g:
            thousand_g_freq = get_frequency(
                                            self.thousand_g, 
                                            variant['CHROM'], 
                                            variant['POS'],
                                            variant['ALT'].split(',')[0]
                                            )
            if thousand_g_freq:
                variant['1000G_freq'] = thousand_g_freq
        if self.exac:
            exac_freq = get_frequency(
                                        self.exac, 
                                        variant['CHROM'], 
                                        variant['POS'],
                                        variant['ALT'].split(',')[0]
                                    )
            if exac_freq:
                variant['ExAC_freq'] = exac_freq
        return
    
    
    def make_print_version(self, variant_dict):
        """Get the variants ready for printing"""
        
        for variant_id in variant_dict:
            variant = variant_dict[variant_id]
                        
            vcf_info = variant_dict[variant_id]['INFO'].split(';')
            
            feature_list = variant_dict[variant_id].get('annotation', set())
            
            # variant[compounds] is a dictionary with family id as keys and a set of compounds as values
            compounds = variant.get('compounds', dict())
            # Here we store the compound strings that should be added to the variant:
            family_compound_strings = []
            
            for family_id in compounds:
                compound_string = ''
                compound_set = compounds[family_id]
                #We do not want reference to itself as a compound:
                compound_set.discard(variant_id)
                # If there are any compounds for the family:
                if compounds[family_id]:
                    compound_string = '|'.join(compound_set)
                    family_compound_strings.append(':'.join([family_id, compound_string]))

            # We need to check if compounds have already been annotated.
            if 'Compounds' not in variant['info_dict']:
                if len(family_compound_strings) > 0:
                    vcf_info.append('Compounds=' + ','.join(family_compound_strings))
            
            # Check if any genetic models are followed
            if 'GeneticModels' not in variant['info_dict'] and self.families:
                # Here we store the compound strings that should be added to the variant:
                family_model_strings = []
                model_scores = {}
                genetic_models = variant.get('inheritance_models', {})
                for family_id in genetic_models:
                    model_string = ''
                    model_list = []
                    for model in genetic_models[family_id]:
                        if genetic_models[family_id][model]:
                            model_list.append(model)
                        model_string = '|'.join(model_list)
                    if len(model_list) > 0:
                        family_model_strings.append(':'.join(
                                                    [
                                                        family_id, 
                                                        model_string
                                                    ]
                                                    )
                                                )
                        model_scores[family_id] = str(get_model_score(
                                                    self.families[
                                                        family_id
                                                    ].individuals, variant))
                
                if len(family_model_strings) > 0:
                    vcf_info.append(
                                'GeneticModels=' + 
                                ','.join(family_model_strings)
                                )
                    model_score_list = []
                    for family_id in model_scores:
                        if model_scores[family_id]:
                            if float(model_scores[family_id]) > 0:
                                model_score_list.append(
                                    ':'.join(
                                                [
                                                    family_id, 
                                                    model_scores[family_id]
                                                ]
                                            )
                                        )
                    if len(model_score_list) > 0:
                        vcf_info.append(
                                    'ModelScore=' +  
                                    ','.join(model_score_list)
                                    )
            
            
            # We only want to include annotations where we have a value
            if 'Annotation' not in variant['info_dict']:
                if len(feature_list) != 0 and feature_list != ['-']:
                    vcf_info.append(
                                'Annotation=' + 
                                ','.join(feature_list)
                                )
            
            if variant.get('CADD', None):
                if 'CADD' not in variant['info_dict']:
                    vcf_info.append(
                                'CADD=%s' % str(variant.pop('CADD', '.'))
                                )
            
            if self.cadd_raw:
                if 'CADD_raw' not in variant['info_dict']:
                    if variant.get('CADD_raw', None):
                        vcf_info.append(
                                    'CADD_raw=%s' % 
                                    str(variant.pop('CADD_raw', '.'))
                                    )
            
            if variant.get('1000G_freq', None):
                if '1000G_freq' not in variant['info_dict']:
                    vcf_info.append(
                                '1000G_freq=%s' % 
                                str(variant.pop('1000G_freq', '.'))
                                )

            if variant.get('ExAC_freq', None):
                if 'ExAC_freq' not in variant['info_dict']:
                    vcf_info.append(
                                'ExAC_freq=%s' % 
                                str(variant.pop('ExAC_freq', '.'))
                                )
            
            variant_dict[variant_id]['INFO'] = ';'.join(vcf_info)
            
        return
    
    def run(self):
        """Run the consuming"""
        if self.verbosity:
            log.info('%s: Starting!' % self.proc_name)
        # Check if there are any batches in the queue
        while True:
            # A batch is a dictionary on the form {variant_id:variant_dict}
            variant_batch = self.task_queue.get()
            
            if variant_batch is None:
                self.logger.info('No more batches')
                self.task_queue.task_done()
                self.logger.info('{0}: Exiting'.format(proc_name))
                break

            
            if self.families:
                check_genetic_models(
                                variant_batch, 
                                self.families, 
                                self.verbosity,
                                self.phased, 
                                self.strict, 
                                proc_name
                            )
            
            # We can now free som space by removing the haploblocks
            variant_batch.pop('haploblocks', None)
            
            # These are family independent annotations which will be done annyway:
            for variant_id in variant_batch:
                if self.any_cadd_info:
                    self.add_cadd_score(variant_batch[variant_id])
                if self.thousand_g or self.exac:
                    self.add_frequency(variant_batch[variant_id])
            
            # Now we want to make versions of the variants that are ready for printing.
            self.make_print_version(variant_batch)
            self.results_queue.put(variant_batch)
            self.task_queue.task_done()
        
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()