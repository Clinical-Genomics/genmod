import logging
import operator
import functools
import math

def get_model_score(individuals, variant):
    """
    Return the model score for this variant.
    
    The score is a estimation on how likely the genetic inheritance pattern
    is for these individuals.
    Scores are based on the genotype calls from variant calling.
    
    Arguments:
        individuals (list): A list with individual id:s
        variant (dict): A variant dictionary
    
    Returns:
        model_score (float): The model score
    
    """
    logger = logging.getLogger(__name__)
    
    model_score = 0
    genotype_scores = []
    
    for individual in individuals:
        logger.debug("Checking gt call for individual {0}".format(
            individual
        ))
        
        gt_call = variant.get('genotypes', {}).get(individual, None)
        if gt_call:
            if gt_call.genotype_quality > 0:
                genotype_scores.append(10**-(float(gt_call.genotype_quality)/10))
    
    if len(genotype_scores) > 0:
        model_score = (
            round(-10*math.log10(1-functools.reduce(
                operator.mul, [1-score for score in genotype_scores]))))
    
    return model_score