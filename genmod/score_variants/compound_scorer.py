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
from typing import Dict, Tuple, Union, List
from multiprocessing import Process

from genmod.vcf_tools import (replace_vcf_info, add_vcf_info)

from genmod.score_variants.score_variant import as_normalized_max_min, MIN_SCORE_NORMALIZED, MAX_SCORE_NORMALIZED
from genmod.score_variants.rank_score_variant_definitions import RANK_SCORE_TYPE_NAMES
from genmod.score_variants.cap_rank_score_to_min_bound import cap_rank_score_to_min_bound

logger = logging.getLogger(__name__)


def get_rank_score(rank_score_type: str,
                   threshold: Union[int, float],
                   min_rank_score_value: float,
                   max_rank_score_value: float) -> Union[int, float]:
    """
    Return raw rank score or normalized rank score.

    Args:
        rank_score_type: Type according to RANK_SCORE_TYPE_NAMES
        threshold: A rank score like value
        min_rank_score_value: input to as_normalized_max_min
        max_rank_score_value: input to as_normalized_max_min
    Returns:
        A rank score like value, possibly normalized
    """
    if rank_score_type == 'RankScore':
        return threshold
    elif rank_score_type == 'RankScoreNormalized':
        # Normalize raw rank score
        return as_normalized_max_min(score=float(threshold),
                                     min_score_value=min_rank_score_value,
                                     max_score_value=max_rank_score_value)
    raise ValueError('Unknown RANK_SCORE_TYPE_NAMES config', rank_score_type)


def get_rank_score_as_magnitude(rank_score_type: str,
                                rank_score: Union[int, float],
                                min_rank_score_value: float,
                                max_rank_score_value: float) -> float:
    """
    Returns rank score as a magnitude (delta), to make the rank score
    suitable for addition/subtraction operations.

    In case of raw rank score, return the input value as the magnitude.

    In case of normalized rank score, return the input value as the
    percentage of the rank score dynamic range.

    Args:
        rank_score_type: Type according to RANK_SCORE_TYPE_NAMES
        rank_score: A rank score like value, to be used in an addition/subtraction operation
        min_rank_score_value: input to as_normalized_max_min
        max_rank_score_value: input to as_normalized_max_min
    Returns:
        A value, magnitude, compatible with raw or normalized rank score values

    """
    if rank_score_type == 'RankScore':
        return rank_score
    elif rank_score_type == 'RankScoreNormalized':
        normalized_rank_score: float = rank_score / (max_rank_score_value - min_rank_score_value)
        if not (MIN_SCORE_NORMALIZED <= normalized_rank_score <= MAX_SCORE_NORMALIZED):
            raise ValueError(f'Failed to normalize to within expected bounds {normalized_rank_score}')
        return normalized_rank_score
    raise ValueError(f'Unknown rank score type {rank_score_type}')

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

    @staticmethod
    def _get_rankscore_normalization_bounds(variant_batch: Dict[str, Dict]) -> Dict[str, Tuple]:
        """
        For all variants in a variant batch, find the rank score normalization
        min-max bounds.

        Returns:
            A dict containing tuple min-max keyed on variant id
        """
        variant_rankscore_normalization_bounds: dict = {}
        for variant_id in variant_batch:
            entry_minmax: List[str] = variant_batch[variant_id]['info_dict']['RankScoreMinMax'].split(':')
            rankscore_normalization_min_max: tuple = (float(entry_minmax[1]), float(entry_minmax[2]))
            if not rankscore_normalization_min_max[0] <= rankscore_normalization_min_max[1]:
                raise ValueError(f'Invalid min-max normalization value expected MIN-MAX \
                {rankscore_normalization_min_max}')
            if variant_id in variant_rankscore_normalization_bounds.keys():
                raise KeyError(f'Cannot add variant ID to normalization data dict since it\'s already present \
                               {variant_id}, {variant_rankscore_normalization_bounds}')
            variant_rankscore_normalization_bounds.update({variant_id: rankscore_normalization_min_max})
        return variant_rankscore_normalization_bounds

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
            for rank_score_type in RANK_SCORE_TYPE_NAMES:
                # Prepare rank_scores dict to contain raw RankScore and RankScoreNormalized compound scores
                rank_scores.update({rank_score_type: dict()})

                for variant_id in variant_batch:
                    # First we store the scores for each variant in a dictionary
                    variant = variant_batch[variant_id]
                    rank_score_entry = variant['info_dict'].get(f'{rank_score_type}', '')

                    # We need to loop through the families
                    # This entry looks like <family_id>:<rank_score>, <family_id>:<rank_score>
                    for family_rank_score in rank_score_entry.split(','):
                        family_rank_score = family_rank_score.split(':')

                        #TODO check if correct family id
                        # Right now we assume that there is only one family in the vcf
                        family_id = family_rank_score[0]
                        rank_score = float(family_rank_score[-1])

                        rank_scores[rank_score_type][variant_id] = rank_score

            # Per variant, find rank score max min values used for normalization
            variant_rankscore_normalization_bounds: Dict[str, Tuple] = \
                self._get_rankscore_normalization_bounds(variant_batch)
            
            #We now have a dictionary with variant ids and rank scores, per rank_score_type
            for variant_id in variant_batch:
                # If the variants only follow AR_comp (and AD for single individual families)
                # we want to pennalise the score if the compounds have low scores
                variant = variant_batch[variant_id]
                raw_compounds = variant['info_dict'].get('Compounds', None)
                for rank_score_type in RANK_SCORE_TYPE_NAMES:
                    if raw_compounds:
                        logger.debug("Scoring compound for variant %s" % variant_id)
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

                        current_rank_score = rank_scores[rank_score_type][variant_id]

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
                            compound_rank_score = rank_scores[rank_score_type][compound_id]
                            if compound_rank_score > get_rank_score(rank_score_type=rank_score_type,
                                                                    threshold=9,
                                                                    min_rank_score_value=variant_rankscore_normalization_bounds[variant_id][0],
                                                                    max_rank_score_value=variant_rankscore_normalization_bounds[variant_id][1]
                                                                    ):
                                only_low = False
                        logger.debug("Setting only_low to {0}".format(only_low))

                        if (correct_score and only_low):
                            logger.debug("correcting rank score for {0}".format(
                                variant_id))
                            current_rank_score -= get_rank_score_as_magnitude(rank_score_type=rank_score_type,
                                                                              rank_score=6,
                                                                              min_rank_score_value=variant_rankscore_normalization_bounds[variant_id][0],
                                                                              max_rank_score_value=variant_rankscore_normalization_bounds[variant_id][1]
                                                                              )
                            # In case the current_rank_score falls outside normalization bounds after modification,
                            # cap it to within the MIN normalization bound.
                            current_rank_score = cap_rank_score_to_min_bound(rank_score_type=rank_score_type,
                                                                             rank_score=current_rank_score,
                                                                             min_rank_score_value=variant_rankscore_normalization_bounds[variant_id][0]
                                                                             )

                        for compound_id in compound_list:
                            logger.debug("Checking compound {0}".format(compound_id))
                            compound_rank_score = rank_scores[rank_score_type][compound_id]
                            # This is the combined score for current variant and
                            # its compound:
                            compound_score = current_rank_score + compound_rank_score
                            # This is the new compound
                            new_compound = "{0}>{1}".format(compound_id, compound_score)
                            scored_compound_list.append(new_compound)

                        new_compound_string = "{0}:{1}".format(
                            compound_family_id, '|'.join(scored_compound_list))

                        current_rank_score = float(current_rank_score)  # Export rank score as float type
                        new_rank_score_string = "{0}:{1}".format(compound_family_id, current_rank_score)

                        # variant['info_dict']['IndividualRankScore'] = current_rank_score_string
                        variant['info_dict'][f'{rank_score_type}'] = new_rank_score_string
                        variant['info_dict'][f'Compounds{rank_score_type.strip("RankScore")}'] = new_compound_string

                        variant = replace_vcf_info(
                            keyword=f'{rank_score_type}',
                            annotation = new_rank_score_string,
                            variant_dict=variant
                        )

                        # CompoundsNormalized is not previously added to VCF.
                        # For this case, perform an VCF INFO ADD operation, rather than a REPLACE
                        keyword_compounds = f'Compounds{rank_score_type.strip("RankScore")}'
                        fn_add_replace_vcf_info = replace_vcf_info
                        if not (keyword_compounds in variant['INFO'] and
                                keyword_compounds in variant['info_dict']):
                            # In case INFO subfield is not previously added to VCF,
                            # there's a need to do so now.
                            fn_add_replace_vcf_info = add_vcf_info
                        variant = fn_add_replace_vcf_info(
                            keyword=keyword_compounds,
                            annotation=new_compound_string,
                            variant_dict=variant
                        )
                logger.debug("Putting variant in results_queue")

                self.results_queue.put(variant)
            
            self.task_queue.task_done()
        
        return
