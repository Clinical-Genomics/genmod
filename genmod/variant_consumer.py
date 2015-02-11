#!/usr/bin/env python
# encoding: utf-8
"""
variant_consumer.py

Class that takes a list of objects and return all unordered pairs as a generator.

If only one object? Raise Exception
 
Created by Måns Magnusson on 2013-03-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import division, print_function, unicode_literals

import sys
import os
import operator
import tabix

from pprint import pprint as pp
from functools import reduce
from multiprocessing import Process
from math import log10

from genmod import check_genetic_models
from genmod.errors import warning

# Import third party library
# https://github.com/mitsuhiko/logbook
from logbook import Logger, StderrHandler
log = Logger('Logbook')
log_handler = StderrHandler()

class VariantConsumer(Process):
    """
    Yeilds all unordered pairs from a list of objects as tuples, 
    like (obj_1, obj_2)
    """
    
    def __init__(self, task_queue, results_queue, families={}, phased=False, 
                vep=False, cadd_raw=False, cadd_file=None, cadd_1000g=None, 
                cadd_exac=None, cadd_ESP=None, cadd_InDels=None, 
                thousand_g=None, exac=None, dbNSFP=None, strict=False, 
                verbosity=False):
        Process.__init__(self)
        self.task_queue = task_queue
        self.families = families
        self.results_queue = results_queue
        self.verbosity = verbosity
        self.phased = phased
        self.vep = vep
        self.cadd_raw = cadd_raw
        self.cadd_file = cadd_file
        self.cadd_1000g = cadd_1000g
        self.cadd_exac = cadd_exac
        self.cadd_ESP = cadd_ESP
        self.cadd_InDels = cadd_InDels
        self.thousand_g = thousand_g
        self.exac = exac
        self.dbNSFP = dbNSFP
        self.strict = strict
        self.any_cadd_info = False
        if self.cadd_file:
            self.cadd_file = tabix.open(self.cadd_file)
            self.any_cadd_info = True
        if self.cadd_1000g:
            self.cadd_1000g = tabix.open(self.cadd_1000g)
            self.any_cadd_info = True
        if self.cadd_exac:
            self.cadd_exac = tabix.open(self.cadd_exac)
            self.any_cadd_info = True
        if self.cadd_ESP:
            self.cadd_ESP = tabix.open(self.cadd_ESP)
            self.any_cadd_info = True
        if self.cadd_InDels:
            self.cadd_InDels = tabix.open(self.cadd_InDels)
            self.any_cadd_info = True
        if self.thousand_g:
            self.thousand_g = tabix.open(self.thousand_g)
        if self.exac:
            self.exac = tabix.open(self.exac)
        if self.dbNSFP:
            self.exac = tabix.open(self.exac)
    
    def get_model_score(self, individuals, variant):
        """Return the model score for this variant."""
        model_score = None
        genotype_scores = []
        
        for individual in individuals:
            gt_call = variant['genotypes'].get(individual, None)
            if gt_call:
                if gt_call.genotype_quality > 0:
                    genotype_scores.append(10**-(float(gt_call.genotype_quality)/10))
        if len(genotype_scores) > 0:
            model_score = (str(
                            round(-10*log10(
                                    1-reduce(
                                        operator.mul, 
                                            [
                                            1-score for score in genotype_scores]
                                            )
                                        )
                                    )
                                )
                            )
        
        return model_score
    
    def get_cadd_score(self, tabix_reader, chrom, start, alt=None):
        """Return the record from a cadd file."""
        score = None
        # CADD values are only for snps:
        cadd_key = int(start)
        try:
            for record in tabix_reader.query(chrom, cadd_key-1, cadd_key):
                if record[3] == alt:
                    #We need to send both cadd values
                    cadd_raw = record[-2]
                    cadd_phred = record[-1]
                    return (cadd_raw, cadd_phred)
        except TypeError:
            for record in tabix_reader.query(str(chrom), cadd_key-1, cadd_key):
                record = record.split('\t')
                if record[3] == alt:
                    #We need to send both cadd values
                    cadd_raw = record[-2]
                    cadd_phred = record[-1]
                    return (cadd_raw, cadd_phred)
        except:
            pass
        
        return score
    
    def add_cadd_score(self, variant):
        """Add the CADD relative score to this variant."""
        cadd_score = None
        cadd_relative_scores = []
        cadd_absolute_scores = []
        #Check CADD file(s):
        for alt in variant['ALT'].split(','):
            if self.cadd_file:
                cadd_score = self.get_cadd_score(
                                        self.cadd_file, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
            # If variant not found in big CADD file check the 1000G file:
            if not cadd_score and self.cadd_1000g:
                cadd_score = self.get_cadd_score(
                                        self.cadd_1000g, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
            
            if not cadd_score and self.cadd_exac:
                cadd_score = self.get_cadd_score(
                                        self.cadd_exac, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
            
            if not cadd_score and self.cadd_ESP:
                cadd_score = self.get_cadd_score(
                                        self.cadd_ESP, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
            
            if not cadd_score and self.cadd_ESP:
                cadd_score = self.get_cadd_score(
                                        self.cadd_ESP, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
            
            if not cadd_score and self.cadd_ESP:
                cadd_score = self.get_cadd_score(
                                        self.cadd_ESP, 
                                        variant['CHROM'], 
                                        variant['POS'], 
                                        alt
                                    )
            
            if cadd_score:
                cadd_relative_scores.append(str(cadd_score[1]))
                cadd_absolute_scores.append(str(cadd_score[0]))
        
        if len(cadd_relative_scores) > 0:
            variant['CADD'] = ','.join(cadd_relative_scores)
            if self.cadd_raw:
                variant['CADD_raw'] = ','.join(cadd_absolute_scores)
        return

    def get_freq(self, tabix_reader, chrom, start, alt):
        """Return the record from a cadd file."""
        freq = None
        # CADD values are only for snps:
        cadd_key = int(start)
        try:
            for record in tabix_reader.query(chrom, cadd_key-1, cadd_key):
                i = 0
                #We can get multiple rows so need to check each one
                #We also need to check each one of the alternatives per row
                for alternative in record[4].split(','):
                    if alternative == alt:
                        for info in record[7].split(';'):
                            info = info.split('=')
                            if info[0] == 'AF':
                                frequencies = info[-1].split(',')
                                return frequencies[i]
                    i += 1
        except TypeError:            
            for record in tabix_reader.query(str(chrom), cadd_key-1, cadd_key):
                i = 0
                #We can get multiple rows so need to check each one
                #We also need to check each one of the alternatives per row
                for alternative in record[4].split(','):
                    if alternative == alt:
                        for info in record[7].split(';'):
                            info = info.split('=')
                            if info[0] == 'AF':
                                frequencies = info[-1].split(',')
                                return frequencies[i]
                    i += 1
        except:
            pass
        
        return freq
    
    def add_frequency(self, variant):
        """Add the thousand genome frequency if present."""
        #Check 1000G frequency:
        thousand_g_freq = None
        exac_freq = None
        if self.thousand_g:
            thousand_g_freq = self.get_freq(self.thousand_g, variant['CHROM'], variant['POS'],
                                                              variant['ALT'].split(',')[0])
            if thousand_g_freq:
                variant['1000G_freq'] = thousand_g_freq
        if self.exac:
            exac_freq = self.get_freq(self.exac, variant['CHROM'], variant['POS'],
                                                              variant['ALT'].split(',')[0])
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
                        model_scores[family_id] = self.get_model_score(
                                                    self.families[
                                                        family_id
                                                    ].individuals, variant)
                
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
            
            if not self.vep:
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
        proc_name = self.name
        if self.verbosity:
            log.info('%s: Starting!' % proc_name)
        while True:
            # A batch is a dictionary on the form {variant_id:variant_dict}
            variant_batch = self.task_queue.get()
            
            if variant_batch is None:
                self.task_queue.task_done()
                if self.verbosity:
                    log.info('%s: Exiting' % proc_name)
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