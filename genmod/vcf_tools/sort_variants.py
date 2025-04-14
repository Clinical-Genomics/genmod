#!/usr/bin/env python
# encoding: utf-8
"""
sort_test.py

Try to sort with unix sort.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import logging
from datetime import datetime
from subprocess import call

logger = logging.getLogger(__name__)


def sort_variants(infile, mode="chromosome"):
    """
    Sort variants based on rank score or chromosome.

    Uses unix sort to sort the variants and overwrites the infile.

    Args:
        infile : A string that is the path to a file
        mode : 'chromosome' or 'rank'
        outfile : The path to an outfile where the variants should be printed

    Returns:
        0 if sorting was performed
        1 if variants where not sorted
    """
    command = [
        "sort",
    ]
    if mode == "chromosome":
        command.append("-k1,1V")  # Version sorting to deal with e.g. Un_* contigs
        command.append("-k3,3n")  # Sorting positions numerically
        command.append("-s")

    elif mode == "rank":
        command.append("-rn")
        command.append("-k1")

    command = command + [infile, "-o", infile]

    logger.info("Start sorting variants...")
    logger.info("Sort command: {0}".format(" ".join(command)))
    sort_start = datetime.now()

    try:
        call(command)
    except OSError as e:
        logger.warning("unix program 'sort' does not seem to exist on your system...")
        logger.warning("genmod needs unix sort to provide a sorted output.")
        logger.warning("Output VCF will not be sorted since genmod can not findunix sort")
        raise e

    logger.info("Sorting done. Time to sort: {0}".format(datetime.now() - sort_start))

    return
