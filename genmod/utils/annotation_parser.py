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
try:
    import cPickle as pickle
except:
    import pickle

from codecs import open, getreader

from pprint import pprint as pp

from genmod.utils import is_number

from interval_tree import interval_tree

class AnnotationParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile, annotation_type, splice_padding = 2, zipped = False):
        super(AnnotationParser, self).__init__()
        self.annotation_type = annotation_type
        
        self.gene_trees = {}# A dictionary with {<chr>:<intervalTree>} the trees are intervals with genes
        self.exon_trees = {}# A dictionary with {<chr>:<intervalTree>} the trees are intervals with exons
        chromosome_stops = {}## A dictionary with information about the last positions on each chromosome:
        
        
        nr_of_genes = 0
                        
        if zipped:
            f = getreader('utf-8')(gzip.open(infile), errors='replace')
        else: 
            f = open(infile, mode='r', encoding='utf-8', errors='replace')
        
        line_count = 0
        
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
        
        for chrom in genes:
            self.gene_trees[chrom] = interval_tree.IntervalTree(genes[chrom], 1, chromosome_stops[chrom])
        for chrom in exons:
            if len(exons[chrom]) > 0:
                self.exon_trees[chrom] = interval_tree.IntervalTree(exons[chrom], 1, chromosome_stops[chrom])
            # If no exons found
            else:
                self.exon_trees[chrom] = interval_tree.IntervalTree([[0,0,None]], 0, 0)
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
                if 'hr' in line[0]:
                    chrom = line[0][3:]
                else:
                    chrom = line[0]
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
            
                if 'hr' in line[0]:
                    chrom = line[0][3:]
                else:
                    chrom = line[0]
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
                chrom = line[0]
                if 'hr' in line[0]:
                    chrom = line[0][3:]
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
        
        for line in ref_file_handle:
            if not line.startswith('#') and len(line) > 1:
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
    args = parser.parse_args()
    anno_file = args.annotation_file[0]
    gene_trees = {}
    exon_trees = {}
    
    gene_db = 'genes.db'
    exon_db = 'exons.db'
    
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
        
        my_parser = AnnotationParser(anno_file, file_type, zipped)
        
        gene_trees = my_parser.gene_trees
        exon_trees = my_parser.exon_trees
        
        with open(gene_db, 'wb') as f:
            pickle.dump(gene_trees, f)
        
        with open(exon_db, 'wb') as g:
            pickle.dump(exon_trees, g)
        
    
    else:
        
        with open(gene_db, 'rb') as f:
            gene_trees = pickle.load(f)
        with open(exon_db, 'rb') as g:
            exon_trees = pickle.load(g)
        
    print('Gene_trees:')
    pp(gene_trees)
    print('Exon_trees:')
    pp(exon_trees)
    pp(gene_trees.get('1',interval_tree.IntervalTree([[0,0,None]], 1, 1)))
    pp(gene_trees.get('1',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([11900, 11901]))
    pp(gene_trees.get('8',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([6913070, 6913071]))
    pp(exon_trees.get('8',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([6913070, 6913071]))
    pp(gene_trees.get('1',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([721190, 821290]))
    pp(exon_trees.get('1',interval_tree.IntervalTree([[0,0,None]], 1, 1)).find_range([721190, 821290]))


if __name__ == '__main__':
    main()
