#!/usr/bin/env python
# encoding: utf-8
"""
annotation_parser.py

This script will parse a file with intervals in .ccds or .bed format and build one interval tree for each chromosome.
So self.chromosomes will look like:

{'1':intervalTree, '2':intervalTree, ..., 'X':intervalTree}

Create a family object and its family members from different types of input file
Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
import gzip

from pprint import pprint as pp

from genmod.utils import interval_tree, is_number

class AnnotationParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile, annotation_type, zipped = False):
        super(AnnotationParser, self).__init__()
        self.annotation_type = annotation_type
                
        self.interval_trees = {}# A dictionary with {<chr>:<intervalTree>}
        
        chromosomes = {} # A dictionary with {<chr>: [feature_1, feature_2, ...]} 
        chromosome_stops = {}# A dictionary with information about the last positions on each chromosome:
        
        if zipped:
            f = gzip.open(infile, 'r')
        else: 
            f = open(infile, 'r')
        line_count = 0
        for line in f:
            if line[0] != '#' and len(line) > 1:
                line = line.rstrip()
                info = {'chrom':'Na', 'start':0, 'stop':0, 'transcript_id':'0', 'gene_id':'0', 'feature_id':str(line_count)}
                
                if self.annotation_type == 'ccds':
                    info = self.ccds_parser(line, info, line_count)
                elif self.annotation_type == 'bed':
                    info = self.bed_parser(line, info, line_count)
                elif self.annotation_type == 'gtf':
                    info = self.gtf_parser(line, info, line_count)
                elif self.annotation_type == 'ref_gene':
                    info = self.ref_gene_parser(line, info, line_count)
                #Feature is a list with [start, stop, id]
                feature = [info['start'], info['stop'], info['feature_id']]
                
                if info['chrom'] in chromosomes:
                    chromosomes[info['chrom']].append(feature)
                else:
                    chromosomes[info['chrom']] = [feature]
                # Update the last end position, is it is bigger than before:
                if info['stop'] > chromosome_stops.get(info['chrom'], 0):
                    chromosome_stops[info['chrom']] = info['stop'] + 1
            line_count += 1
        
        number_of_intervals = 0
        # for chrom, stop in self.chromosome_stops.items():
        #     print chrom, stop, len(self.chromosomes[chrom])
        #     number_of_intervals += len(self.chromosomes[chrom])
        
        #Build one interval tree for each chromosome:
        for chrom in chromosomes:
            self.interval_trees[chrom] = interval_tree.intervalTree(chromosomes[chrom], 
                                        0, 1, 1, chromosome_stops[chrom])
                    
    def bed_parser(self, line, info, line_count):
        """Parse a .bed."""
        line = line.split()
        if 'hr' in line[0]:
            info['chrom'] = line[0][3:]
        else:
            info['chrom'] = line[0]
        info['start'] = int(line[1])
        info['stop'] = int(line[2])
        if len(line) > 3:
            info['gene_id'] = line [3]
            info['feature_id'] = line [3]
        return info
    
    def ccds_parser(self, line, info, line_count):
        """Parse a ccds line"""
        line = line.split('\t')
        if 'hr' in line[0]:
            info['chrom'] = line[0][3:]
        else:
            info['chrom'] = line[0]
        info['transcript_id'] = line[1]
        info['gene_id'] = line[2]
        info['feature_id'] = info['gene_id']
        if is_number.is_number(line[7]) and is_number.is_number(line[8]):
            info['start'] = int(line[7])
            info['stop'] = int(line[8])
        return info
    
    def gtf_parser(self, line, info, line_count):
        """Parse a gtf line"""
        line = line.split('\t')
        if 'hr' in line[0]:
            info['chrom'] = line[0][3:]
        else:
            info['chrom'] = line[0]
        if is_number.is_number(line[3]) and is_number.is_number(line[4]):
            info['start'] = int(line[3])
            info['stop'] = int(line[4])
        info_field = line[8].split(';')[:-1]
        for information in info_field:
            entry = information.split()
            if entry[0] == 'transcript_id':
                info['transcript_id'] = entry[1][1:-1]
            if entry[0] == 'gene_id':
                info['gene_id'] = entry[1][1:-1]
        info['feature_id'] = info['gene_id']
        return info
    
    def ref_gene_parser(self, line, info, line_count):
        """Parse a file in the refGene format"""
        line = line.split('\t')
        if 'hr' in line[2]:
            info['chrom'] = line[2][3:]
        else:
            info['chrom'] = line[2]
        if is_number.is_number(line[4]) and is_number.is_number(line[5]):
            info['start'] = int(line[4])
            info['stop'] = int(line[5])
        info['transcript_id'] = line[1] #???
        info['gene_id'] = line[12]
        info['feature_id'] = info['gene_id']
        return info
        
        

def main():
    parser = argparse.ArgumentParser(description="Parse different kind of annotation files.")
    parser.add_argument('annotation_file', type=str, nargs=1 , help='A file with anotations.')
    parser.add_argument('-bed', '--bed', action="store_true", help='Annotation file is in bed format.')
    parser.add_argument('-ccds', '--ccds', action="store_true", help='Annotation file is in ccds format.')
    parser.add_argument('-gtf', '--gtf', action="store_true", help='Annotation file is in gtf format.')
    parser.add_argument('-ref_gene', '--ref_gene', action="store_true", help='Annotation file is in gtf format.')
    args = parser.parse_args()
    infile = args.annotation_file[0]
    # TODO write a check for zipped files
    file_type = 'ccds'
    if args.bed:
        file_type = 'bed'
    if args.ccds:
        file_type = 'ccds'
    if args.gtf:
        file_type = 'gtf'
    if args.ref_gene:
        file_type = 'ref_gene'
        
    my_parser = AnnotationParser(infile, file_type)
    pp(my_parser.interval_trees)


if __name__ == '__main__':
    main()
