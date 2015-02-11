#!/usr/bin/env python
# encoding: utf-8
"""
sort_test.py

Try to sort with unix sort.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os
import click


from subprocess import call
from datetime import datetime

# Import third party library
# https://github.com/mitsuhiko/logbook
from logbook import Logger, StderrHandler
log = Logger('Logbook')
log_handler = StderrHandler()


def sort_variants(infile, mode='chromosome', verbose=False):
    """
    Sort variants based on rank score or chromosome.
    
    Uses unix sort to sort the variants and overwrites the infile.
    
    Args:
        infile : A string that is the path to a file
        mode : 'chromosome' or 'rank'
        outfile : The path to an outfile where the variants should be printed
        verbose : Increase output verbosity
    
    Returns:
        0 if sorting was performed
        1 if variants where not sorted
    """
    command = [
            'sort',
            ]
    if mode == 'chromosome':
        command.append('-n')
        command.append('-k1')
        command.append('-k3')

    elif mode == 'rank':
        command.append('-rn')
        command.append('-k1')

    command = command + [infile, '-o', infile]

    if verbose:
        log.info("Start sorting variants...")
        log.info("Sort command: %s" % ' '.join(command))
        sort_start = datetime.now()
    
    try:
        call(command)
    except OSError:
        if verbose:
            log.warn("unix command sort does not seem to exist on your system...")
            log.warn("genmod needs unix sort to provide a sorted output.")
        log.warn("""Output VCF will not be sorted since genmod can not find
                unix sort""")
        return 1

    if verbose:
        log.info("Sorting done!")
        log.info("Time to sort %s" % (datetime.now()-sort_start))
    
    return 0

@click.command()
@click.argument('csv_file',
                nargs=1,
                type=click.Path(exists=True)
)
@click.option('-m' ,'--mode', 
                type=click.Choice(['chromosome', 'score']), 
                default='chromosome',
                help='Specify the mode on how we should sort the variants.'
)
@click.option('-v', '--verbose',
              is_flag=True,
              help='Increase output verbosity.'
)
def cli(csv_file,mode,verbose):
    sort_variants(csv_file, mode, verbose)

if __name__ == '__main__':
    cli()

