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

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import argparse
import gzip
import genmod
try:
    import cPickle as pickle
except:
    import pickle

from datetime import datetime
from codecs import open, getreader

from pprint import pprint as pp

from genmod import is_number

from interval_tree import interval_tree

class AnnotationParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile, annotation_type, zipped = False, splice_padding = 2, verbosity=False):
        super(AnnotationParser, self).__init__()
        self.verbosity = verbosity
        self.annotation_type = annotation_type
        
        self.gene_trees = {}# A dictionary with {<chr>:<intervalTree>} the trees are intervals with genes
        self.exon_trees = {}# A dictionary with {<chr>:<intervalTree>} the trees are intervals with exons
        chromosome_stops = {}## A dictionary with information about the last positions on each chromosome:
        infile_name, infile_extension = os.path.splitext(infile)
        
        genes = {}
        exons = {}
        nr_of_genes = 0
        zipped = False
    
        if infile_extension == '.gz':
            zipped = True        
        
        if zipped:
            f = getreader('utf-8')(gzip.open(infile), errors='replace')
        else: 
            f = open(infile, mode='r', encoding='utf-8', errors='replace')
        
        
        if self.annotation_type == 'gene_pred':
            genes,exons = self.gene_pred_parser(f, splice_padding)
        
        elif self.annotation_type == 'gtf':
            genes,exons = self.gtf_parser(f, splice_padding)
        
        elif self.annotation_type == 'ccds':
            genes,exons = self.ccds_parser(f, splice_padding)
        
        elif self.annotation_type == 'bed':
            genes,exons = self.bed_parser(f, splice_padding)
        
        for chrom in dict(genes):
            # prepare the intervals for the tree:
            gene_intervals = []
            exon_intervals = []
            for gene_id in dict(genes[chrom]):
                feature = [genes[chrom][gene_id]['gene_start'],
                        genes[chrom][gene_id]['gene_stop'], gene_id]
                gene_intervals.append(feature)
                
                # Update the end position of the interval
                if genes[chrom][gene_id]['gene_stop'] > chromosome_stops.get(chrom, 0):
                    chromosome_stops[chrom] = genes[chrom][gene_id]['gene_stop'] + 1
                
            genes[chrom] = gene_intervals
            
            for exon in dict(exons.get(chrom,{})):
                exon_intervals.append([exon[0], exon[1], ':'.join(exons[chrom][exon])])
                
            exons[chrom] = exon_intervals
        
        #Build one interval tree for each chromosome:
                
        nr_of_genes = 0
        
        for chrom in genes:
            nr_of_genes += len(genes[chrom])
            self.gene_trees[chrom] = interval_tree.IntervalTree(genes[chrom], 1, chromosome_stops[chrom])
        for chrom in exons:
            if len(exons[chrom]) > 0:
                self.exon_trees[chrom] = interval_tree.IntervalTree(exons[chrom], 1, chromosome_stops[chrom])
            # If no exons found
            else:
                self.exon_trees[chrom] = interval_tree.IntervalTree([[0,0,None]], 0, 0)
        
        if self.verbosity:
            print('Number of genes in annotation file: %s' % nr_of_genes)
                
        return
    
    def add_gene(self, genes, chrom, start, stop, gene_id):
        """Add a gene to a gene dict"""
        if chrom not in genes:
            genes[chrom] = {}
        
        if gene_id in genes[chrom]:
            # If this transcript starts before the previous, update the start of the gene:
            if start < genes[chrom][gene_id]['gene_start']:
                genes[chrom][gene_id]['gene_start'] = start
            # If this transcript ends after the previous update the stop of the gene:
            if stop > genes[chrom][gene_id]['gene_stop']:
                genes[chrom][gene_id]['gene_stop'] = stop
        else:
            genes[chrom][gene_id] = {'gene_start':start, 'gene_stop':stop}
        return
    
    def add_exon(self, exons, chrom, start, stop, feature_id, splice_padding):
        """Add a exon to a exon dict"""
        if chrom not in exons:
            exons[chrom] = {}
        exon = (start - splice_padding, stop + splice_padding)
        if exon in exons[chrom]:
            exons[chrom][exon].append(feature_id)
        else:
            exons[chrom][exon] = [feature_id]
        return
    
    def bed_parser(self, bed_file_handle, splice_padding):
        """Parse a .bed."""
        genes = {}
        exons = {}        
        for line in bed_file_handle:
            if not line.startswith('#') and len(line) > 1:
                line = line.split()
                feature_id = 'None'
                chrom = line[0].lstrip('chr')
                feature_start = int(line[1])
                feature_stop = int(line[2])
                if len(line) > 3:
                    feature_id = line [3]
            
                self.add_gene(genes, chrom, feature_start, feature_stop, feature_id)
                self.add_exon(exons, chrom, feature_start, feature_stop, feature_id, splice_padding)
            
        return genes, exons
    
    def ccds_parser(self, ccds_file_handle, splice_padding):
        """Parse a ccds line"""
        genes = {}
        exons = {}
        
        for line in ccds_file_handle:
            if not line.startswith('#') and len(line) > 1:
                line = line.split('\t')                
                chrom = line[0].lstrip('chr')
                transcript_id = line[1]
                gene_id = line[2]
                if is_number.is_number(line[7]) and is_number.is_number(line[8]):
                    feature_start = int(line[7])
                    feature_stop = int(line[8])
                    #TODO raise exception?
                
                    self.add_gene(genes, chrom, feature_start, feature_stop, gene_id)
                
                    for interval in line[9][1:-1].split(','):
                        boundaries = (interval.split('-'))
                        exon_start = int(boundaries[0].lstrip())
                        exon_stop = int(boundaries[1].lstrip())
        
        return genes, exons
    
    def gtf_parser(self, gtf_file_handle, splice_padding):
        """Parse a gtf file"""
        genes = {}
        exons = {}
        
        for line in gtf_file_handle:
            if not line.startswith('#') and len(line) > 1:
                transcript_id = ''
                gene_id = ''
                gene_name = ''
                
                line = line.split('\t')
                if len(line) < 5:
                    line = line.split()
                chrom = line[0].lstrip('chr')
                if is_number.is_number(line[3]) and is_number.is_number(line[4]):
                    feature_start = int(line[3])
                    feature_stop = int(line[4])
                #TODO Raise exception?
                info_field = line[8].split(';')[:-1]
                for information in info_field:
                    entry = information.split()
                    if entry[0] == 'transcript_id':
                        transcript_id = entry[1][1:-1]
                    if entry[0] == 'gene_id':
                        gene_id = entry[1][1:-1]
                    if entry[0] == 'gene_name':
                        gene_name = entry[1][1:-1]
                
                if line[2] == 'gene':
                    self.add_gene(genes, chrom, feature_start, feature_stop, gene_id)
                 
                elif line[2] == 'exon':
                    self.add_exon(exons, chrom, feature_start, feature_stop, transcript_id, splice_padding)
                
        
        return genes,exons
    
    def gene_pred_parser(self, ref_file_handle, splice_padding):
        """Parse a file in the refGene format, we should add the information about gene or transcript here"""
        genes = {} # A dictionary with {<chr>: {gene_id:{'gene_start':0,'gene_stop':0}}
        exons = {} # A dictionary with {<chr>: {(exon_start, exon_stop, transcript_id_1):[gene_id_1, ..]
        ### Features look like [start,stop, feature_id] ###
        nr_of_lines = 0
        for line in ref_file_handle:
            if not line.startswith('#') and len(line) > 1:
                line = line.split('\t')
                nr_of_lines += 1
                transcript = line[1]
                chrom = line[2].lstrip('chr')
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
                                
                self.add_gene(genes, chrom, transc_start, transc_stop, gene_id)
                for i in range(len(exon_starts)):
                    self.add_exon(exons, chrom, exon_starts[i], exon_stops[i], transcript_id, splice_padding)
                
        return genes,exons
        
        

def main():
    parser = argparse.ArgumentParser(description="Parse different kind of annotation files. If no annotation file given, use genes.db and exons.db")
    parser.add_argument('-an','--annotation_file', type=str, nargs=1, default=[None], help='A file with anotations.')
    parser.add_argument('-bed', '--bed', action="store_true", help='Annotation file is in bed format.')
    parser.add_argument('-ccds', '--ccds', action="store_true", help='Annotation file is in ccds format.')
    parser.add_argument('-gtf', '--gtf', action="store_true", help='Annotation file is in gtf format.')
    parser.add_argument('-gene_pred', '--gene_pred', action="store_true", help='Annotation file is in gene pred format. This is the format for the ref_seq and ensembl gene files.')
    parser.add_argument('-v', '--verbose', action="store_true", help='Increase output verbosity.')
    args = parser.parse_args()
    anno_file = args.annotation_file[0]
    gene_trees = {}
    exon_trees = {}
    
    outpath = os.path.join(os.path.split(os.path.dirname(genmod.__file__))[0], 'annotations/')
    gene_db = os.path.join(outpath, 'genes.db')
    exon_db = os.path.join(outpath, 'exons.db')
    
    if args.verbose:
        print(gene_db)
        print(exon_db)
        print('Start parsing annotation!')
        start = datetime.now()
    
    if anno_file:
        file_name, file_extension = os.path.splitext(anno_file)
        zipped = False
        if file_extension == '.gz':
            zipped = True
            file_name, file_extension = os.path.splitext(file_name)
        
        file_type = 'gene_pred'
        if args.bed or file_extension[1:] == 'bed':
            file_type = 'bed'
        if args.ccds or file_extension[1:] == 'ccds':
            file_type = 'ccds'
        if args.gtf or file_extension[1:] == 'gtf':
            file_type = 'gtf'
        if args.gene_pred or file_extension[1:] in ['ref_gene', 'gene_pred']:
            file_type = 'gene_pred'
        
        my_parser = AnnotationParser(anno_file, file_type, zipped=zipped, verbosity=args.verbose)
        
        gene_trees = my_parser.gene_trees
        exon_trees = my_parser.exon_trees
    
    else:
        try:
            with open(gene_db, 'rb') as f:
                gene_trees = pickle.load(f)
            with open(exon_db, 'rb') as g:
                exon_trees = pickle.load(g)
        except FileNotFoundError:
            print('You need to build annotations! See documentation.')
            pass
            # raise FileNotFoundError("You need to make annotations for genmod with an annotation file!")
            #TODO write help line
    
    if args.verbose:
        print('Gene chromosomes:%s' % sorted(gene_trees.keys()))
        print('Exon chromosomes:%s' % sorted(exon_trees.keys()))
        pp(gene_trees.get('1',interval_tree.IntervalTree([[0,0,None]], 1, 1)))
        pp(gene_trees.get('1',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([11900, 11901]))
        pp(gene_trees.get('8',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([6913070, 6913071]))
        pp(exon_trees.get('8',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([6913070, 6913071]))
        pp(gene_trees.get('1',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([721190, 821290]))
        pp(exon_trees.get('1',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([721190, 821290]))
        print('Time to parse file: %s' % str(datetime.now()-start))


if __name__ == '__main__':
    main()
