#!/usr/bin/env python
# encoding: utf-8
"""
vcf_header.py


Parse the header of a vcf file.

Create a variant objects and a dictionary with individuals that have a dictionary with genotypes for each variant.

Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

### TODO make a proper vcf parser ###


class HeaderParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile):
        super(HeaderParser, self).__init__()
        self.metadata=OrderedDict()
        self.header=[]
        self.metadata_counter = 1
        self.line_counter = 0
        start_symbol = '#'
        with open(infile, 'rb') as f:
            for line in f:
                self.line_counter += 1
                start_symbol = line[0]
                if line[0] == '#':
                    if line [1] == '#':
                        line = line[2:].rstrip()
                        self.metadata[self.metadata_counter] = line
                        self.metadata_counter += 1
                    else:
                        self.header = line[1:].rstrip().split('\t')
                else:
                    break
                        


def main():
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    args = parser.parse_args()
    infile = args.variant_file[0]
    my_parser = HeaderParser(infile)
    for line in my_parser.metadata:
        print line, my_parser.metadata[line]
    print '\t'.join(my_parser.header)
    print my_parser.line_counter
    # for individual in my_parser.individuals:
    #     for genotype in my_parser.individuals[individual]:
    #         print individual, genotype, my_parser.individuals[individual][genotype]
    


if __name__ == '__main__':
    main()
