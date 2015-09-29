"""
Score a variant line based on the score functions in config_parser

"""

import logging

def score_variant(variant, config_parser):
    """Score a variant line
    
        Args:
            variant (dict): A variant dictionary
            config_parser (ConfigParser): A config parser object with score functions
        Returns:
            rank_score (float): The rank score for this variant
    """
    logger = logging.getLogger(__name__)
    
    rank_score = 0
    
    for category in config_parser.categories:
        logger.debug("Checking scores for category {0}".format(category))
        #We save all the scores for a category
        category_aggregate = config_parser.categories[category]['category_aggregation']
        category_scores = []
        category_score = None
        
        for plugin_name in config_parser.categories[category]['plugins']:
            logger.debug("Checking scores for plugin {0}".format(plugin_name))
            plugin = config_parser.plugins[plugin_name]
            score_function = config_parser.score_functions[plugin_name]
            value = plugin.get_value(variant_dict=variant)
            logger.debug("Found value {0} for plugin {1}".format(value, plugin_name))
            score = score_function.get_score(value)
            logger.debug("Score is {0} for plugin {1}".format(score, plugin_name))
            category_scores.append(score)

        if category_scores:
        
            if category_aggregate == 'max' and category_scores:
                logger.debug("Take the max score for category {0}".format(
                    category))
                category_score = max(category_scores)
                logger.debug("Max value is {0}".format(
                    category_score))
            elif category_aggregate == 'min' and category_scores:
                logger.debug("Take the min score for category {0}".format(
                    category))
                category_score = min(category_scores)
                logger.debug("Min value is {0}".format(
                    category_score))
            elif category_aggregate == 'sum' and category_scores:
                logger.debug("Take the sum of scores score for category {0}".format(
                    category))
                category_score = sum(category_scores)
                logger.debug("Sum of scores is {0}".format(
                    category_score))
            
            logger.debug("Adding category score {0} to rank_score".format(category_score))
            rank_score += category_score
            logger.debug("Updating rank score to {0}".format(rank_score))
        else:
            logger.debug("No scores found for category {0}".format(category))
            
        
    # logger.debug("Rank score was found to be {0} for variant {1}".format(
    #     rank_score, variant.get('variant_id', 'unknown variant id')))
            
    return rank_score