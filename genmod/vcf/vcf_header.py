#!/usr/bin/env python
# encoding: utf-8
"""
vcf_parser.py


Parse a vcf file.

Includes a header class for storing information about the headers.
Create variant objects and a dictionary with individuals that have a dictionary with genotypes for each variant.


Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import gzip
import re
import argparse
if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from pprint import pprint as pp
### TODO make a proper vcf parser ###


class HeaderParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self):
        super(HeaderParser, self).__init__()
        self.info_lines=[]
        self.filter_lines=[]
        self.contig_lines=[]
        self.format_lines=[]
        self.alt_lines=[]
        self.other_lines=[]
        self.header=[]
        self.header_keys={'info' : ['ID', 'Number', 'Type', 'Description'], 
                            'form' : ['ID', 'Number', 'Type', 'Description'], 
                            'filt' : ['ID', 'Description'],
                            'alt' : ['ID', 'Description'],
                            'contig' : ['ID', 'length']}
        self.fileformat = ''
        self.metadata_counter = 1
        self.line_counter = 0
        self.individuals = []
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.contig_pattern = re.compile(r'''\#\#contig=<
            ID=(?P<id>[^,]+),
            .*
            length=(?P<length>-?\d+)
            .*
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.alt_pattern = re.compile(r'''\#\#ALT=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')
    
    def parse_meta_data(self, line):
        """Parse a vcf metadataline"""
        line = line.rstrip()
        line_info = line[2:].split('=')
        match = False
        
        if line_info[0] == 'fileformat':
            self.fileformat = line_info[1]
        elif line_info[0] == 'INFO':
            match = self.info_pattern.match(line)
            if not match:
                raise SyntaxError("One of the INFO lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('number'), match.group('type'), match.group('desc')]
            self.info_lines.append(dict(zip(self.header_keys['info'],matches)))
        elif line_info[0] == 'FILTER':
            match = self.filter_pattern.match(line)
            if not match:
                raise SyntaxError("One of the FILTER lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('desc')]
            self.filter_lines.append(dict(zip(self.header_keys['filt'],matches)))
        elif line_info[0] == 'contig':
            match = self.contig_pattern.match(line)
            if not match:
                raise SyntaxError("One of the contig lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('length')]
            self.contig_lines.append(dict(zip(self.header_keys['contig'],matches)))
        elif line_info[0] == 'FORMAT':
            match = self.format_pattern.match(line)
            if not match:
                raise SyntaxError("One of the FORMAT lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('number'), match.group('type'), match.group('desc')]
            self.format_lines.append(dict(zip(self.header_keys['form'],matches)))
        elif line_info[0] == 'ALT':
            match = self.alt_pattern.match(line)
            if not match:
                raise SyntaxError("One of the ALT lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('desc')]
            self.alt_lines.append(dict(zip(self.header_keys['alt'],matches)))
        else:
            match = self.meta_pattern.match(line)
            if not match:
                raise SyntaxError("One of the meta data lines is malformed: %s" % line)
            self.other_lines.append({match.group('key'): match.group('val')})
        # self.metadata[self.metadata_counter] = line
        # self.metadata_counter += 1
    
    def parse_header_line(self, line):
        """docstring for parse_header_line"""
        self.header = line[1:].rstrip().split('\t')
        self.individuals = self.header[9:]
        
    def print_header(self):
        """Returns a list with the header lines if proper format"""
        pass


class VCFParser(object):
    """docstring for VCFParser"""
    def __init__(self, infile):
        super(VCFParser, self).__init__()
        self.infile = infile
        self.metadataparser = HeaderParser()
        self.individuals = []
        file_name, file_extension = os.path.splitext(self.infile)
        if file_extension == '.gz':
            self.vcf = gzip.open(self.infile)
        elif file_extension == '.vcf':
            self.vcf = open(self.infile, 'rb')
        else:
            raise SyntaxError('File is not in a supported format.')
        print file_name, file_extension
            
    def parse(self):
        """Start parsing the vcf"""
        for line in self.vcf:
            if line.startswith('##'):
                self.metadataparser.parse_meta_data(line)
            elif line.startswith('#'):
                self.metadataparser.parse_header_line(line)
            else:
                # Here comes the variant lines
                pass
        self.individuals = self.metadataparser.individuals
        self.header = self.metadataparser.header
        # pp(self.metadataparser.info_lines)
        # print ''
        # pp(self.metadataparser.filter_lines)
        # print ''
        # pp(self.metadataparser.format_lines)
        # print ''
        # pp(self.metadataparser.contig_lines)
        # print ''
        # pp(self.metadataparser.alt_lines)
        # print ''
        # pp(self.metadataparser.other_lines)
        # print ''
        # pp(self.metadataparser.header)
        # print ''
        # pp(self.metadataparser.individuals)
        
        

def main():
    parser = argparse.ArgumentParser(description="Parse vcf headers.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    args = parser.parse_args()
    infile = args.variant_file[0]
    my_parser = VCFParser(infile)
    my_parser.parse()
    # for line in my_parser.metadata:
    #     print line, my_parser.metadata[line]
    # print '\t'.join(my_parser.header)
    # print my_parser.line_counter
    # print my_parser.individuals
    # for individual in my_parser.individuals:
    #     for genotype in my_parser.individuals[individual]:
    #         print individual, genotype, my_parser.individuals[individual][genotype]
    


if __name__ == '__main__':
    main()
