#!/usr/bin/env python
# encoding: utf-8
"""
print_headers.py

Print the lines in a header object to a file or stdout.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

from codecs import open

def print_headers(head, outfile=None, silent=False):
    """
    Print the vcf headers.
    
    If a result file is provided headers will be printed here, otherwise
    they are printed to stdout.
    
    Args:
        head    : A vcf header object
        outfile : Path to outfile or None
        silent  : Bool. If nothing should be printed.
        
    """
    if outfile:
        with open(outfile, 'w', encoding='utf-8') as f:
            for head_count in head.print_header():
                f.write(head_count+'\n')
    else:
        if not silent:
            for line in head.print_header():
                print(line)
    return

