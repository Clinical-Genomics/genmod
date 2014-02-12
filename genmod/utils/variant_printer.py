#!/usr/bin/env python
# encoding: utf-8
"""
variant_printer.py


Print the variants of a results queue to a file.


Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import multiprocessing
from pprint import pprint as pp


class VariantPrinter(multiprocessing.Process):
    """docstring for VariantPrinter"""
    def __init__(self, task_queue, outfile, chromosomes, verbosity=False):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.outfile = outfile
        self.verbosity = verbosity
        # self.lock = lock
        self.chromosomes = chromosomes
    
    def run(self):
        """Starts the printing"""
        # Print the results to a temporary file:
        number_of_finished = 0
        proc_name = self.name
        if self.verbosity:
            print proc_name ,'starting!'
        while True:
            next_result = self.task_queue.get()
            if self.verbosity:
                if self.task_queue.full():
                    print 'Printing queue full'
            if next_result is None:
                if self.verbosity:
                    print 'All variants printed!'
                self.outfile.close()
                break
            else:
                for variant_id in next_result:
                    # print '\t'.join(next_result[variant_id].values())
                    # self.chromosomes[next_result[variant_id['CHROM']]].write(
                    self.outfile.write('\t'.join(next_result[variant_id].values())+'\n')
        return

def main():
    pass

if __name__ == '__main__':
    main()
