"""
Score a variant line based on the score functions in config_parser

"""

import logging

logger = logging.getLogger(__name__)

def get_plugin_score(variant, plugin_name, config_parser):
    """Return the score found for a plugin
    
        Args:
            variant (dict): A variant dictionary
            plugin_name (str): A plugin name
            config_parser (ConfigParser): A config parser object with score functions
        Returns:
            category_score (float): The rank score for this variant
    
    """
    logger.debug("Checking scores for plugin {0}".format(plugin_name))
    # This is a plugin to collect the correct value from the vcf
    plugin = config_parser.plugins[plugin_name]
    # This is the score function for this plugin
    score_function = config_parser.score_functions[plugin_name]
    
    value = plugin.get_value(variant_dict=variant)
    logger.debug("Found value {0} for plugin {1}".format(value, plugin_name))
    
    # Score is allways a number
    score = score_function.get_score(value)
    logger.debug("Score is {0} for plugin {1}".format(score, plugin_name))
    
    return score

def get_category_score(variant, category, config_parser):
    """Return the score for a given category
    
        Args:
            variant (dict): A variant dictionary
            category (str): A score category
            config_parser (ConfigParser): A config parser object with score functions
        Returns:
            category_score (float): The rank score for this variant
    """
    
    category_score = 0
    
    logger.debug("Checking scores for category {0}".format(category))
    #We save all the scores for a category
    category_aggregate = config_parser.categories[category]['category_aggregation']
    category_scores = []
    
    for plugin_name in config_parser.categories[category]['plugins']:
        category_scores.append(get_plugin_score(
            variant = variant, 
            plugin_name = plugin_name, 
            config_parser = config_parser
            )
        )

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
        
    else:
        logger.debug("No scores found for category {0}".format(category))
    
    return category_score
    

def score_variant(variant, config_parser):
    """Score a variant line
    
        Args:
            variant (dict): A variant dictionary
            config_parser (ConfigParser): A config parser object with score functions
        Returns:
            rank_score (float): The rank score for this variant
    """
    
    rank_score = 0
    
    for category in config_parser.categories:
        logger.debug("Checking scores for category {0}".format(category))
        category_score = get_category_score(
            variant = variant,
            category = category,
            config_parser = config_parser
        )

        logger.debug("Adding category score {0} to rank_score".format(category_score))
        rank_score += category_score
        logger.debug("Updating rank score to {0}".format(rank_score))

    return rank_score
