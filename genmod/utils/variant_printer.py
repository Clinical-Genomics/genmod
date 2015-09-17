#!/usr/bin/env python
# encoding: utf-8
"""
variant_printer.py

Print the variants of a results queue to a file.

Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import (print_function, absolute_import)

from multiprocessing import Process
from codecs import open

from . import get_chromosome_priority

import logging

class VariantPrinter(Process):
    """
    Print variants to a temporary file.
    
    There are two modes that will be used when sorting the file.
    'chromosome': In this case the priority order of the chromosome is printed
                    to the first position in the results file
    
    'score': Here the Individual score of the variant is printed to the first
            position of the results file
    
    Args:
        task_queue : A variants queue with batches of variants
        head : The header line to specify what from the variant object to print
        mode : 'chromosome' or 'score'. See above.
        chr_map : If mode='chromosome' we need a map to specify the sort order
                    of the chromosomes. (Soon to come?)
        outfile : File that all variants should be printed to
    
    """
    def __init__(self, task_queue, head, mode='chromosome', outfile = None):
        Process.__init__(self)
        self.logger = logging.getLogger(__name__)
        self.task_queue = task_queue
        self.outfile = outfile
        self.header = head.header
        self.mode = mode
    
    def run(self):
        """Starts the printing"""
        # Print the results to a temporary file:
        number_of_finished = 0
        proc_name = self.name
        self.logger.info(('{0}: starting'.format(proc_name)))
        if self.outfile:
            f = open(self.outfile, 'w', encoding="utf-8")
        
        while True:
            
            
            # A task is a variant dictionary
            self.logger.debug(('{0} fetching next variant'.format(proc_name)))
            variant = self.task_queue.get()
            
            if self.task_queue.full():
                self.logger.warning('Variant queue full')
            
            if variant is None:
                self.logger.info('All variants printed.')
                if self.outfile:
                    f.close()
                break
            
            self.logger.debug("Printing variant {0}".format(variant.get('variant_id', 'unknown')))
                
            priority = get_chromosome_priority(variant['CHROM'])
            
            print_line = [priority] + [variant.get(entry, '.') 
                            for entry in self.header]
                
            if self.outfile:
                f.write('\t'.join(print_line) + '\n')
            else:
                print('\t'.join(print_line))
        
        return

