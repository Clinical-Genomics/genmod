#!/usr/bin/env python
# encoding: utf-8
"""
variant_printer.py

Print the variants of a results queue to a file.

Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

from multiprocessing import Process

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
        temp_file : The temporary file that all variants should be printed to
        head : The header line to specify what from the variant object to print
        mode : 'chromosome' or 'score'. See above.
        chr_map : If mode='chromosome' we need a map to specify the sort order
                    of the chromosomes
        verbosity : Increase output verbosity
    
    """
    def __init__(self, task_queue, temp_file, head, mode='chromosome', 
                verbosity=False):
        Process.__init__(self)
        self.logger = logging.getLogger(__name__)
        self.task_queue = task_queue
        self.verbosity = verbosity
        self.temp_file = temp_file
        self.header = head.header
        self.mode = mode
    
    def run(self):
        """Starts the printing"""
        # Print the results to a temporary file:
        number_of_finished = 0
        proc_name = self.name
        self.logger.info(('{0}: starting'.format(proc_name)))
        
        while True:
            self.logger.debug(('{0} fetching next variant batch'.format(proc_name)))
            # A task is a dictionary with variants where variant_id is key and
            # a variant dict is value
            next_result = self.task_queue.get()
            
            if self.task_queue.full():
                self.logger.warning('Variant queue full')
            
            if next_result is None:
                self.logger.info('All variants printed.')
                self.temp_file.close()
                break
                
            else:
                
                for variant_id in next_result:
                    variant = next_result[variant_id]
                    
                    if self.mode == 'score':
                        try:
                            priority = variant['Individual_rank_score']
                        except KeyError:
                            priority = '0'
                    elif self.mode == 'chromosome':
                        chrom = variant['CHROM']
                        if chrom.startswith('chr'):
                            chrom = chrom[3:]
                        try:
                            priority = int(chrom)
                        except ValueError:
                            if chrom == 'X':
                                priority = 23
                            elif chrom == 'Y':
                                priority = 24
                            elif chrom == 'MT':
                                priority = 25
                            else:
                                priority = 26
                    else:
                        raise SyntaxError("""Need to specify priority mode for 
                                            printing the variants""")
                    
                    print_line = [str(priority)] + [variant.get(entry, '-')
                                                 for entry in self.header]
                    
                    self.temp_file.write('\t'.join(print_line) + '\n')
    
def main():
    pass

if __name__ == '__main__':
    main()
