"""Vendored from extract_vcf (https://github.com/moonso/extract_vcf, MIT License, Måns Magnusson)."""

import logging

logger = logging.getLogger(__name__)


def split_strings(string, separators):
    """
    Split a string with arbitrary number of separators.
    Return a list with the splitted values

    Arguments:
        string (str): ex. "a:1|2,b:2"
        separators (list): ex. [',',':','|']

    Returns:
         results (list) : ex. ['a','1','2','b','2']
    """
    logger = logging.getLogger("genmod.plugin.split_strings")
    logger.debug("splitting string '{0}' with separators {1}".format(string, separators))
    results = []

    def recursion(recursive_string, separators, i=1):
        """
        Split a string with arbitrary number of separators.
        Add the elements of the string to global list result.

        Arguments:
            string : ex. "a:1|2,b:2"
            separators (list): ex. [',',':','|']

        Returns:
             Adds splitted string to results. ex. ['a','1','2','b','2']
        """
        if i == len(separators):
            for value in recursive_string.split(separators[i - 1]):
                logger.debug("Adding {0} to results".format(value))
                results.append(value)
        else:
            for value in recursive_string.split(separators[i - 1]):
                recursion(value, separators, i + 1)

    if len(separators) > 0:
        recursion(string, separators)
    else:
        results = [string]

    return results
