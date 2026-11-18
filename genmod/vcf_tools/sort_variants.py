#!/usr/bin/env python
# encoding: utf-8
"""
sort_variants.py

Sort a file of variants with unix sort.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import annotations

import logging
from datetime import datetime
from subprocess import call
from typing import Literal

logger = logging.getLogger(__name__)

SORT_OPTIONS = {
    # Version sorting (-V) deals with e.g. Un_* contigs
    "chromosome": ["-k1,1", "-k3,3", "-s", "-V"],
    "vcf": ["-k1,1", "-k2,2", "-s", "-V"],
    "rank": ["-rn", "-k1"],
}


def sort_variants(infile: str, mode: Literal["chromosome", "vcf", "rank"] = "chromosome") -> None:
    """Sort variants based on rank score or chromosome.

    Uses unix sort to sort the variants and overwrites the infile.

    Args:
        infile: The path to the file to sort, which is overwritten
        mode: ``chromosome`` sorts on the first and third column (a temporary file with a
            sort column first), ``vcf`` on chromosome and position of a plain vcf, and
            ``rank`` on the first column, highest number first

    Raises:
        ValueError: If the sort mode is not known
        OSError: If unix ``sort`` is not available
    """
    if mode not in SORT_OPTIONS:
        raise ValueError(f"Unknown sort mode '{mode}', use one of {list(SORT_OPTIONS)}")

    command = ["sort", *SORT_OPTIONS[mode], infile, "-o", infile]

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
