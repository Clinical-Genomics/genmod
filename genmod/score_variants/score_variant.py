"""
Score a variant line based on the score functions in config_parser

"""

import logging
from typing import Any, Tuple, List

logger = logging.getLogger(__name__)

MIN_SCORE_NORMALIZED: float = 0.0
MAX_SCORE_NORMALIZED: float = 1.0


def get_plugin_score(variant, plugin_name, config_parser, csq_format=None) -> Tuple[Any, float, float]:
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
    
    value = plugin.get_value(variant_dict=variant, csq_format=csq_format)
    logger.debug("Found value {0} for plugin {1}".format(value, plugin_name))
    
    # Score is allways a number
    score = score_function.get_score(value)
    logger.debug("Score is {0} for plugin {1}".format(score, plugin_name))

    plugin_score_min: float = score_function.score_min
    plugin_score_max: float = score_function.score_max
    logger.debug(f'Minmax scores for plugin {plugin_name} is ({plugin_score_min},{plugin_score_max})')
    
    return score, plugin_score_min, plugin_score_max


def get_category_score(variant, category, config_parser, csq_format=None) -> Tuple[int, float, float]:
    """Return the score for a given category.

       A category (such as 'allele_frequency') groups multiple plugin scores.
       This method selects final score based on category_aggregation [Categories]
       rule (see genmod_example.ini)

        Args:
            variant (dict): A variant dictionary
            category (str): A score category
            config_parser (ConfigParser): A config parser object with score functions
        Returns:
            category_score, sum of min and max scores for this category
    """

    category_score = 0

    logger.debug("Checking scores for category {0}".format(category))
    #We save all the scores for a category
    category_aggregate = config_parser.categories[category]['category_aggregation']
    category_scores = []
    category_score_mins: List[float] = []
    category_score_maxs: List[float] = []

    for plugin_name in config_parser.categories[category]['plugins']:
        plugin_score, plugin_score_min, plugin_score_max = get_plugin_score(
            variant = variant,
            plugin_name = plugin_name,
            config_parser = config_parser,
            csq_format = csq_format
            )
        category_scores.append(plugin_score)
        # Add the maximal and minimal score points that can be provided for this category
        category_score_mins.append(plugin_score_min)
        category_score_maxs.append(plugin_score_max)

    # Provide default plugin bounds in case there were no plugins in category
    category_score_min: float = 0.0
    category_score_max: float = 0.0

    if category_scores:

        if category_aggregate == 'max' and category_scores:
            logger.debug("Take the max score for category {0}".format(
                category))
            category_score = max(category_scores)
            category_score_min = max(category_score_mins)
            category_score_max = max(category_score_maxs)
            logger.debug("Max value is {0}".format(
                category_score))
        elif category_aggregate == 'min' and category_scores:
            logger.debug("Take the min score for category {0}".format(
                category))
            category_score = min(category_scores)
            category_score_min = min(category_score_mins)
            category_score_max = min(category_score_maxs)
            logger.debug("Min value is {0}".format(
                category_score))
        elif category_aggregate == 'sum' and category_scores:
            logger.debug("Take the sum of scores score for category {0}".format(
                category))
            category_score = sum(category_scores)
            category_score_min = sum(category_score_mins)
            category_score_max = sum(category_score_maxs)
            logger.debug("Sum of scores is {0}".format(
                category_score))

    else:
        logger.debug("No scores found for category {0}".format(category))

    if not (category_score_min <= category_score <= category_score_max):
        raise ValueError('Category score outside expected category score range',
                          category_score_min, category_score, category_score_max)
    return category_score, category_score_min, category_score_max


def as_normalized_max_min(score: float,
                          min_score_value: float,
                          max_score_value: float) -> float:
    """
    Performs max-min normalization of score to range (0, 1).
    Args:
        score: Rank score
        min_score_value: Lower bound on score
        max_score_value: Upper bound on score
    Returns:
        Score in range (0, 1)
    """
    for key, value in {'score': score,
                       'min_score_value': min_score_value,
                       'max_score_value': max_score_value}.items():
        if not isinstance(value, float):
            raise TypeError(f'Bad type {key} {type(value)} {value}')
    if not max_score_value >= min_score_value:
        raise ValueError('Inverted minmax values for normalization', max_score_value, min_score_value)
    score_normalized: float = (score - min_score_value) / (max_score_value - min_score_value)
    if not (MIN_SCORE_NORMALIZED <= score_normalized <= MAX_SCORE_NORMALIZED):
        raise ValueError('Failed to normalize to within expected bounds',
                                        min_score_value, max_score_value, score, score_normalized)
    return score_normalized
