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

from __future__ import (division, print_function)

import sys
import os
import logging

from multiprocessing import Process

from genmod.vcf_tools import (replace_vcf_info, add_vcf_info)

logger = logging.getLogger(__name__)

class CompoundScorer(Process):
    """
    Annotates variant in batches from the task queue and puts the result in 
    the results queue.
    """
    
    def __init__(self, task_queue, results_queue, individuals):
        """
        Initialize the VariantAnnotator
        
        Consume variant batches from the task queue, annotate them with the 
        genetic inheritance patterns that they follow and put them in the 
        results queue.
        
        Arguments:
            task_queue (Queue)
            results_queue (Queue)
            individuals (list)
        """
        Process.__init__(self)
        
        self.proc_name = self.name
        
        logger.info("Setting up variant_annotator: {0}".format(
            self.proc_name))
        
        logger.debug("Setting up task queue")
        self.task_queue = task_queue
        
        logger.debug("Setting up results queue")
        self.results_queue = results_queue

        logger.debug("Setting up individuals")
        self.individuals = individuals
        
        if len(self.individuals) == 1:
            self.models = ['AR_comp', 'AR_comp_dn', 'AD', 'AD_dn']
        else:
            self.models = ['AR_comp', 'AR_comp_dn']
            

    def run(self):
        """Run the consuming"""
        logger.info('%s: Starting!' % self.proc_name)
        # Check if there are any batches in the queue
        while True:
            # A batch is a dictionary with varints on the form {variant_id:variant_dict}
            logger.debug("Getting task from task_queue")
            variant_batch = self.task_queue.get()
            
            if variant_batch is None:
                logger.info('No more batches')
                self.task_queue.task_done()
                logger.info('{0}: Exiting'.format(self.proc_name))
                break
            
            # We need to save the compound scores in a dict and group them by family
            # This is a dictionary on the form {'variant_id: rank_score}
            rank_scores = {}
            # We are now going to score the compounds
            for variant_id in variant_batch:
                # First we store the scores for each variant in a dictionary
                variant = variant_batch[variant_id]
                rank_score_entry = variant['info_dict'].get('RankScore', '')
                
                # We need to loop through the families
                # This entry looks like <family_id>:<rank_score>, <family_id>:<rank_score>
                for family_rank_score in rank_score_entry.split(','):
                    family_rank_score = family_rank_score.split(':')
                    
                    #TODO check if correct family id
                    # Right now we assume that there is only one family in the vcf
                    family_id = family_rank_score[0]
                    rank_score = int(family_rank_score[-1])
                    
                    rank_scores[variant_id] = rank_score
            
            #We now have a dictionary with variant ids and rank scores
            for variant_id in variant_batch:
                # If the variants only follow AR_comp (and AD for single individual families)
                # we want to pennalise the score if the compounds have low scores
                variant = variant_batch[variant_id]
                raw_compounds = variant['info_dict'].get('Compounds', None)
                
                if raw_compounds:
                    logger.debug("Scoring compound for variant {0}".format(variant_id))
                    variant = variant_batch[variant_id]
                    #Variable to see if we should correct the rank score
                    correct_score = True
                    # First we check if the rank score should be corrected:
                    for family in variant['info_dict'].get('GeneticModels', '').split(','):
                        for model in family.split(':')[-1].split('|'):
                            # If the variant follows any model more than the specified it should
                            # not be corrected
                            if model not in self.models:
                                correct_score = False
                    
                    logger.debug("Setting correct_score to {0}".format(correct_score))
                    
                    current_rank_score = rank_scores[variant_id]
                    
                    logger.debug("Current rank score is {0}".format(current_rank_score))
                    
                    scored_compound_list = []
                    only_low = True
                    
                    # One entry per family, splitted on ','
                    # family_id and compounds splitted with ':'
                    # list of compounds splitted on '|'
                    
                    #TODO Only checks first family now
                    family_compound_entry = raw_compounds.split(',')[0]
                    splitted_entry = family_compound_entry.split(':')
                    compound_family_id = splitted_entry[0]
                    compound_list = splitted_entry[-1].split('|')
                    
                    logger.debug("Checking compounds for family {0}".format(
                        compound_family_id))
                    
                    #Loop through compounds to check if they are only low scored
                    for compound_id in compound_list:
                        compound_rank_score = rank_scores[compound_id]
                        if compound_rank_score > 9:
                            only_low = False
                        
                    logger.debug("Setting only_low to {0}".format(only_low))
                    
                    if (correct_score and only_low):
                        logger.debug("correcting rank score for {0}".format(
                            variant_id))
                        current_rank_score -= 6
                        
                    for compound_id in compound_list:
                        logger.debug("Checking compound {0}".format(compound_id))
                        compound_rank_score = rank_scores[compound_id]
                        # This is the combined score for current variant and 
                        # its compound:
                        compound_score = current_rank_score + compound_rank_score
                        # This is the new compound
                        new_compound = "{0}>{1}".format(compound_id, compound_score)
                        scored_compound_list.append(new_compound)

                    new_compound_string = "{0}:{1}".format(
                        compound_family_id, '|'.join(scored_compound_list))

                    new_rank_score_string = "{0}:{1}".format(compound_family_id, current_rank_score)
                    
                    # variant['info_dict']['IndividualRankScore'] = current_rank_score_string
                    variant['info_dict']['RankScore'] = new_rank_score_string
                    variant['info_dict']['Compounds'] = new_compound_string
                    
                    variant = replace_vcf_info(
                        keyword = 'RankScore',
                        annotation = new_rank_score_string, 
                        variant_dict=variant
                    )
                    
                    variant = replace_vcf_info(
                        keyword = 'Compounds',
                        annotation = new_compound_string, 
                        variant_dict=variant
                    )
                logger.debug("Putting variant in results_queue")
                self.results_queue.put(variant)
            
            self.task_queue.task_done()
        
        return
