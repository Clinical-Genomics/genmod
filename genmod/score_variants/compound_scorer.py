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
        self.logger = logging.getLogger(__name__)
        
        self.proc_name = self.name
        
        self.logger.info("Setting up variant_annotator: {0}".format(
            self.proc_name))
        
        self.logger.debug("Setting up task queue")
        self.task_queue = task_queue
        
        self.logger.debug("Setting up results queue")
        self.results_queue = results_queue

        self.logger.debug("Setting up individuals")
        self.individuals = individuals
        
        if len(self.individuals) == 1:
            self.models = ['AR_comp', 'AR_comp_dn', 'AD', 'AD_dn']
        else:
            self.models = ['AR_comp', 'AR_comp_dn']
            

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
            
            # We need to save the compound scores in a dict and group them by family
            # This is a dictionary on the form {'variant_id: rank_score}
            compound_scores = {}
            # We are now going to score the compounds
            for variant_id in variant_batch:
                # If the variants only follow AR_comp (and AD for single individual families)
                # We want to pennalise the score if the compounds have low scores
                variant = variant_batch[variant_id]
                correct_score = True
                for family in variant['info_dict'].get('GeneticModels', '').split(','):
                    for model in family.split(':')[-1].split('|'):
                        if model not in self.models:
                            correct_score = False
                
                
                current_rank_score_entry = variant['info_dict'].get('RankScore', None)
                
                if current_rank_score_entry:
                    # We need to loop through the families
                    for current_family_rank_score in current_rank_score_entry.split(','):
                        current_family_rank_score = current_family_rank_score.split(':')
                        current_family_id = current_family_rank_score[0]
                        current_rank_score = int(current_family_rank_score[-1])
                
                raw_compounds = variant['info_dict'].get('Compounds', None)
                
                new_compound_list = []
                only_low = True
                
                # We know that there are compounds since they are here
                for family_compounds in raw_compounds.split(','):
                    splitted_entry = family_compounds.split(':')
                    # This is the id of the family compounds we are looking at
                    campound_family_id = splitted_entry[0]
                    # We loop over the compounds of this family:
                    for compound in splitted_entry[-1].split('|'):
                        if compound in compound_scores:
                            compound_rank_score = compound_scores[compound]
                        else:
                            try:
                                compound_rank_score_entry = variant_batch[compound]['info_dict'].get('RankScore', None)
                            except KeyError as e:
                                logger.error("The compounds have been annotated with wrong annotation set")
                                raise e
                            
                            if compound_rank_score_entry:
                                for family_score in compound_rank_score_entry.split(','):
                                    family_score = family_score.split(':')
                                    compound_family_id = family_score[0]
                                    compound_rank_score = int(family_score[-1])
                                    compound_scores[compound] = compound_rank_score
                                    if compound_rank_score > 10:
                                        only_low = False
                                
                        compound_score = current_rank_score + compound_rank_score
                        new_compound = "{0}>{1}".format(compound, compound_score)
                        new_compound_list.append(new_compound)
                        
                new_compound_string = "{0}:{1}".format(current_family_id, '|'.join(new_compound_list))
                
                # This is the rank score string before it is corrected:
                # current_rank_score_string = "{0}:{1}".format(current_family_id, current_rank_score)
                
                if correct_score:
                    if only_low:
                        current_rank_score -= 6
                
                new_rank_score_string = "{0}:{1}".format(current_family_id, current_rank_score)

                # variant['info_dict']['IndividualRankScore'] = current_rank_score_string
                variant['info_dict']['RankScore'] = new_rank_score_string
                variant['info_dict']['Compounds'] = new_compound_string
                
                # variant = add_vcf_info(
                #     keyword = 'IndividualRankScore',
                #     variant_dict=variant,
                #     annotation=current_rank_score_string)
                
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
                self.logger.debug("Putting variant in results_queue")
                self.results_queue.put(variant)
            
            self.task_queue.task_done()
        
        return
        
    

def main():
    pass

if __name__ == '__main__':
    main()