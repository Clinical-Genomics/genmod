#!/usr/bin/env python
# encoding: utf-8
"""
annotation_parser.py

This script will parse a file with intervals in .ccds or .bed format and build one interval tree for each chromosome.
So self.chromosomes will look like:

{'1':intervalTree, '2':intervalTree, ..., 'X':intervalTree}

The intervals represent features that are annotated in the infile.


###UPDATE 2014-04-15###

This parser has been a bit naive since files can have one transcript per line and some one gene per line.
It is hard to write a tool to parse everything. I will now focus on the refGene format that is one transcript per line and each line is like:

from http://genome.ucsc.edu/FAQ/FAQformat.html#format9:

table refFlat
"A gene prediction with additional geneName field."
    (
    string  geneName;           "Name of gene as it appears in Genome Browser."
    string  name;               "Name of gene"
    string  chrom;              "Chromosome name"
    char[1] strand;             "+ or - for strand"
    uint    txStart;            "Transcription start position"
    uint    txEnd;              "Transcription end position"
    uint    cdsStart;           "Coding region start"
    uint    cdsEnd;             "Coding region end"
    uint    exonCount;          "Number of exons"
    uint[exonCount] exonStarts; "Exon start positions"
    uint[exonCount] exonEnds;   "Exon end positions"
    )



Create a family object and its family members from different types of input file
Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import unicode_literals
import sys
import os
import argparse
import gzip
from codecs import open, getreader

from pprint import pprint as pp

from genmod.utils import is_number

from interval_tree import interval_tree

class AnnotationParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile, annotation_type, zipped = False):
        super(AnnotationParser, self).__init__()
        self.annotation_type = annotation_type
                
        self.gene_trees = {}# A dictionary with {<chr>:<intervalTree>} the trees are intervals with genes
        
        self.exon_trees = {}# A dictionary with {<chr>:<intervalTree>} the trees are intervals with exons
        
        file_name, file_extension = os.path.splitext(infile)
        nr_of_genes = 0
        if file_extension == '.gz':
            zipped = True
                
        if zipped:
            f = getreader('utf-8')(gzip.open(infile), errors='replace')
        else: 
            f = open(infile, mode='r', encoding='utf-8')
        line_count = 0
        
        if self.annotation_type == 'ref_gene':
            chromosomes,exons, chromosome_stops = self.ref_gene_parser(f)

        elif self.annotation_type == 'gtf':
            chromosomes,exons, chromosome_stops = self.gtf_parser(f)
        
        else:
            chromosomes = {} # A dictionary with {<chr>: [feature_1, feature_2, ...]} 
            chromosome_stops = {}# A dictionary with information about the last positions on each chromosome:
            
            for line in f:
                if not line.startswith('#') and len(line) > 1:
                    line = line.rstrip()
                    # Info allways contain the important information that is needed to create a feature
                    info = {'chrom':'Na', 'start':0, 'stop':0, 'transcript_id':'0', 'gene_id':'0', 'feature_id':str(line_count)}
                    # print(line)
                    if self.annotation_type == 'ccds':
                        info = self.ccds_parser(line, info, line_count)
                    elif self.annotation_type == 'bed':
                        info = self.bed_parser(line, info, line_count)
                    elif self.annotation_type == 'ref_gene':
                        info = self.ref_gene_parser(line, info, line_count)
                    #Feature is a list with [start, stop, id]
                    feature = [info['start'], info['stop'], info['feature_id']]
                    
                    if info['chrom'] in chromosomes:
                        chromosomes[info['chrom']].append(feature)
                    else:
                        chromosomes[info['chrom']] = [feature]
                    # Update the last end position, if it is bigger than before:
                    
                    if info['stop'] > chromosome_stops.get(info['chrom'], 0):
                        chromosome_stops[info['chrom']] = info['stop'] + 1
                line_count += 1
        
        number_of_intervals = 0

        #Build one interval tree for each chromosome:
        for chrom in chromosomes:
            self.gene_trees[chrom] = interval_tree.IntervalTree(chromosomes[chrom], 1, chromosome_stops[chrom])
        # for chrom in chromosomes:
            # self.exon_trees[chrom] = interval_tree.IntervalTree(exons[chrom], 1, chromosome_stops[chrom])
                    
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
    
    def gtf_parser(self, gtf_file_handle):
        """Parse a gtf file"""
        genes = {}
        exons = {}
        chromosomes = {}
        chromosome_stops = {}
        line_count = 0
        
        for line in gtf_file_handle:
            if not line.startswith('#') and len(line) > 1:
                line_count += 1
                line = line.split('\t')
                chrom = line[0]
                if 'hr' in line[0]:
                    chrom = line[0][3:]
                if is_number.is_number(line[3]) and is_number.is_number(line[4]):
                    feature_start = int(line[3])
                    feature_stop = int(line[4])
                info_field = line[8].split(';')[:-1]
                for information in info_field:
                    entry = information.split()
                    if entry[0] == 'transcript_id':
                        transcript_id = entry[1][1:-1]
                    if entry[0] == 'gene_id':
                        gene_id = entry[1][1:-1]
                
                if line[2] == 'gene':
                    if chrom not in genes:
                        genes[chrom] = {}
                        # raw_exons[chrom] = {}
                
                    if gene_id in genes[chrom]:
                        # If this transcript starts before update the start of the gene:
                        if feature_start < genes[chrom][gene_id]['gene_start']:
                            genes[chrom][gene_id]['gene_start'] = feature_start
                        # If this transcript ends after update the stop of the gene:
                        if feature_stop > genes[chrom][gene_id]['gene_stop']:
                            genes[chrom][gene_id]['gene_stop'] = feature_stop
                    else:
                        genes[chrom][gene_id] = {'gene_start':feature_start, 'gene_stop':feature_stop}
                        
        for chrom in genes:
            # prepare the intervals for the tree:
            if chrom not in chromosomes:
                chromosomes[chrom] = []
                exons[chrom] = []
            for gene_id in genes[chrom]:
                feature = [genes[chrom][gene_id]['gene_start'],
                        genes[chrom][gene_id]['gene_stop'], gene_id]
                chromosomes[chrom].append(feature)
                
                # Update the end position of the interval
                if genes[chrom][gene_id]['gene_stop'] > chromosome_stops.get(chrom, 0):
                    chromosome_stops[chrom] = genes[chrom][gene_id]['gene_stop'] + 1
        
        return chromosomes,exons,chromosome_stops
    
    def ref_gene_parser(self, ref_file_handle):
        """Parse a file in the refGene format, we should add the information about gene or transcript here"""
        genes = {} # A dictionary with {<chr>: [feature_1, feature_2, ...]} 
        raw_exons = {} # A dictionary with {exon_start:{exon_stop:[<feature_id_1>...]}}
        exons = {} # A dictionary with {<chr>: [feature_1, feature_2, ...]} 
        chromosomes = {}# A dictionary with information about the last positions on each chromosome:
        chromosome_stops = {}# A dictionary with information about the last positions on each chromosome:
        
        for line in ref_file_handle:
            line = line.split('\t')
            transcript = line[1]
            chrom = line[2]
            if 'hr' in chrom:
                chrom = chrom[3:]
            direction = line[3]
            transc_start = int(line[4])
            transc_stop = int(line[5])
            coding_start = int(line[6])
            coding_stop = int(line[7])
            nr_of_exons = int(line[8])
            exon_starts = [int(boundary) for boundary in line[9].split(',')[:-1]]
            exon_stops = [int(boundary) for boundary in line[10].split(',')[:-1]]
            gene_id = line[12]
            transcript_id = ':'.join([transcript, gene_id])
            
            if chrom not in genes:
                genes[chrom] = {}
                raw_exons[chrom] = {}
            
            if gene_id in genes[chrom]:
                # If this transcript starts before update the start of the gene:
                if transc_start < genes[chrom][gene_id]['gene_start']:
                    genes[chrom][gene_id]['gene_start'] = transc_start
                # If this transcript ends after update the stop of the gene:
                if transc_stop > genes[chrom][gene_id]['gene_stop']:
                    genes[chrom][gene_id]['gene_stop'] = transc_stop
            else:
                genes[chrom][gene_id] = {'gene_start':transc_start, 'gene_stop':transc_stop}
            
            
            for i in range(len(exon_starts)):
                start = exon_starts[i]
                stop = exon_stops[i]
                # Check if the exon is already seen
                if start in raw_exons[chrom]:
                    if stop in raw_exons[chrom][start]:
                        raw_exons[chrom][start][stop].append(transcript_id)
                    else:
                        raw_exons[chrom][start][stop] = [transcript_id]
                else:
                    raw_exons[chrom][start] = {stop:[transcript_id]}
        
        for chrom in genes:
            # prepare the intervals for the tree:
            if chrom not in chromosomes:
                chromosomes[chrom] = []
                exons[chrom] = []
            for gene_id in genes[chrom]:
                feature = [genes[chrom][gene_id]['gene_start'],
                        genes[chrom][gene_id]['gene_stop'], gene_id]
                chromosomes[chrom].append(feature)
                
                # Update the end position of the interval
                if genes[chrom][gene_id]['gene_stop'] > chromosome_stops.get(chrom, 0):
                    chromosome_stops[chrom] = genes[chrom][gene_id]['gene_stop'] + 1
            
            for start in raw_exons[chrom]:
                for stop in raw_exons[chrom][start]:
                    feature = [start, stop, ';'.join(raw_exons[chrom][start][stop])]
                    exons[chrom].append(feature)
        
        return chromosomes,exons,chromosome_stops
        
        

def main():
    parser = argparse.ArgumentParser(description="Parse different kind of annotation files.")
    parser.add_argument('annotation_file', type=str, nargs=1 , help='A file with anotations.')
    parser.add_argument('-bed', '--bed', action="store_true", help='Annotation file is in bed format.')
    parser.add_argument('-ccds', '--ccds', action="store_true", help='Annotation file is in ccds format.')
    parser.add_argument('-gtf', '--gtf', action="store_true", help='Annotation file is in gtf format.')
    parser.add_argument('-gene_pred', '--gene_pred', action="store_true", help='Annotation file is in gene pred format. This is the format for the ref_seq and ensembl gene files.')
    args = parser.parse_args()
    infile = args.annotation_file[0]
    file_name, file_extension = os.path.splitext(infile)
    if file_extension == '.gz':
        file_name, file_extension = os.path.splitext(infile)
        
    # TODO write a check for zipped files
    file_type = 'gene_pred'
    if args.bed or file_extension[1:] == 'bed':
        file_type = 'bed'
    if args.ccds or file_extension[1:] == 'ccds':
        file_type = 'ccds'
    if args.gtf or file_extension[1:] == 'gtf':
        file_type = 'gtf'
    if args.gene_pred or file_extension[1:] in ['ref_gene', 'gene_pred']:
        file_type = 'gene_pred'
        
    my_parser = AnnotationParser(infile, file_type)
    pp(my_parser.gene_trees)
    pp(my_parser.gene_trees['1'])
    pp(my_parser.gene_trees['1'].find_range([11900, 11901]))
    pp(my_parser.gene_trees['1'].find_range([721290, 721291]))
    # pp(my_parser.exon_trees['1'].find_range([721190, 821290]))
    pp(my_parser.gene_trees['1'].find_range([721190, 821290]))


if __name__ == '__main__':
    main()
