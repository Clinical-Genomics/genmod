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
from tempfile import NamedTemporaryFile
from pprint import pprint as pp


class VariantPrinter(multiprocessing.Process):
    """docstring for VariantPrinter"""
    def __init__(self, task_queue, temp_dir,verbosity=False):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.verbosity = verbosity
        self.file_handles = {}
        self.temp_dir = temp_dir
    
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
                for chromosome in self.file_handles:
                    self.file_handles[chromosome].close()
                break
                
            else:
                for variant_id in next_result:
                    variant_chrom = next_result[variant_id]['CHROM']
                    if variant_chrom in self.file_handles:
                        self.file_handles[variant_chrom].write('\t'.join(next_result[variant_id].values()) + '\n')
                    else:
                        self.file_handles[variant_chrom] = NamedTemporaryFile(prefix=variant_chrom+'_', dir=self.temp_dir, delete=False)
        return
    
    def get_temp_files(self):
        """Return the dictionary with temp files"""
        print self.file_handles
        return self.file_handles

def main():
    pass

if __name__ == '__main__':
    main()
