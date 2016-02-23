#!/usr/bin/env python
# encoding: utf-8
"""
annotation_parser.py

This script will parse a file with intervals in .ccds or .bed format and build 
one interval tree for each chromosome.
So self.chromosomes will look like:

{'1':intervalTree, '2':intervalTree, ..., 'X':intervalTree}

The intervals represent features that are annotated in the infile.


###UPDATE 2014-04-15###

This parser has been a bit naive since files can have one transcript per line 
and some one gene per line.
It is hard to write a tool to parse everything. It will now focus on the refGene 
format that is one transcript per line and each line is like:

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

import sys
import os
import gzip
import logging

# For testing:
import click

try:
    import cPickle as pickle
except:
    import pickle

from datetime import datetime
from codecs import open, getreader

from pprint import pprint as pp

from genmod.utils import is_number

from interval_tree import interval_tree

logger = logging.getLogger(__name__)

def parse_annotations(infile, annotation_type, zipped = False, 
                    splice_padding = 2):
    """
    Parses a file with genetic region annotations.
    
    Go through the file and creates one interval tree with genes and one with
    exons.
    
    Arguments:
        infile (str): Path to the annotation file
        annotation_type (str): Specifies the file format
        zipped (bool): If the annotation file is zipped
        splice_padding (int): How many bases around each exon should be added
                              (These are the canonical splice regions)
    """
    
    # A dictionary with {<chr>:<intervalTree>} 
    # the trees are intervals with genes
    logger.info("Initializing gene trees")
    gene_trees = {}
    # A dictionary with {<chr>:<intervalTree>} 
    # the trees are intervals with exons
    logger.info("Initializing exon trees")
    exon_trees = {}
    # A dictionary with information about the last 
    # positions on each chromosome:
    logger.info("Initializing chromosome stops")
    chromosome_stops = {}
    
    # Chromosome stops are necessary to build interval trees
    infile_name, infile_extension = os.path.splitext(infile)
    
    genes = {}
    exons = {}
    nr_of_genes = 0
    zipped = False
    
    if infile_extension == '.gz':
        zipped = True        
    
    if zipped:
        annotation_handle = getreader('utf-8')(gzip.open(infile), errors='replace')
    else: 
        annotation_handle = open(infile, mode='r', encoding='utf-8', errors='replace')
    
    if annotation_type == 'gene_pred':
        genes,exons = gene_pred_parser(annotation_handle, splice_padding)
    
    # elif annotation_type == 'gtf':
    #     genes,exons = gtf_parser(annotation_handle, splice_padding)
    
    elif annotation_type == 'ccds':
        genes,exons = ccds_parser(annotation_handle, splice_padding)
        
    
    elif annotation_type == 'bed':
        genes,exons = bed_parser(annotation_handle, splice_padding)
    
    elif annotation_type == 'gff':
        genes,exons = gff_parser(annotation_handle, splice_padding)

    
    for chrom in genes:
        # prepare the intervals for the tree:
        gene_intervals = []
        exon_intervals = []
        
        if chrom not in chromosome_stops:
            chromosome_stops[chrom] = 0
        
        for gene_id in genes[chrom]:
            feature = genes[chrom][gene_id]
            feature_stop = feature[1]
            
            gene_intervals.append(feature)
            
            # Update the end position of the interval
            if feature_stop > chromosome_stops.get(chrom, 0):
                chromosome_stops[chrom] = feature_stop + 1

        genes[chrom] = gene_intervals

        if chrom in exons:
            for exon_id in exons[chrom]:
                exon_feature = exons[chrom][exon_id]
                exon_intervals.append(exon_feature)
            
        exons[chrom] = exon_intervals
    
    #Build one interval tree for each chromosome:
            
    nr_of_genes = 0
    
    for chrom in genes:
        nr_of_genes += len(genes[chrom])
        gene_trees[chrom] = interval_tree.IntervalTree(
                                                    genes[chrom], 
                                                    1, 
                                                    chromosome_stops[chrom]
                                                    )
    
    logger.info("Gene interval trees created")
    
    for chrom in exons:
        if len(exons[chrom]) > 0:
            exon_trees[chrom] = interval_tree.IntervalTree(
                                                    exons[chrom], 
                                                    1, 
                                                    chromosome_stops[chrom]
                                                    )
        # If no exons found
        else:
            exon_trees[chrom] = interval_tree.IntervalTree(
                                                    [[0,0,None]], 
                                                    0, 
                                                    0
                                                    )

    logger.info("Exon interval trees created")
    
    logger.info("Number of genes in annotation file: {0}".format(nr_of_genes))
    
    return gene_trees, exon_trees
        

def get_coordinates(genes, start, stop, gene_id):
    """
    Return the new coordinates for a gene.
    Since the start and stop can describe different transcripts of the same
    gene we need to check if the gene is already added.
    
    Arguments:
        genes (dict): A dictionary with gene id:s as keys and
                      a list with [start, stop, id] as value
        start (int): An integer with the start of a gene
        stop (int): The stop position of a gene
        gene_id (str): The gene id
    
    Returns:
        gene_info: A list on the form [start, stop, id]
    """
    
    gene_info = genes.get(gene_id)
    
    if gene_info:
        # If this transcript starts after the previous, 
        # we keep the previous start:
        if start > gene_info[0]:
            start = gene_info[0]
        # If this transcript ends befort the previous update 
        # we keep the old one:
        if stop < gene_info[1]:
            stop = gene_info[1]
    
    return [start, stop, gene_id]


def ccds_parser(ccds_file_handle, splice_padding=2):
    """
    Parser a file in the ccds format.
    
    Arguments:
        ccds_file_handle (file_handle): An opened file in ccds format
        splice_padding (int): An integer that describes how we should expand 
                              the exons
    Returns:
        genes (dict): A dictionary chromosome ids as keys and gene dictionarys
                      as values.
        exons (dict): A dictionary chromosome ids as keys and exon dictionarys
                      as values.
    """
    genes = {}
    exons = {}
    
    for line in ccds_file_handle:
        
        if not line.startswith('#') and len(line) > 1:
            
            line = line.split('\t')
            chrom = line[0].lstrip('chr')
            transcript_id = line[1]
            gene_id = line[2]
            
            if chrom not in genes:
                genes[chrom] = {}
            
            if chrom not in exons:
                exons[chrom] = {}
            
            if is_number(line[7]) and is_number(line[8]):
                feature_start = int(line[7])
                feature_stop = int(line[8])
                #TODO raise exception?
                
                genes[chrom][gene_id] = get_coordinates(
                    genes[chrom],
                    feature_start,
                    feature_stop,
                    gene_id
                )
                
                for interval in line[9][1:-1].split(','):
                    boundaries = (interval.split('-'))
                    exon_start = int(boundaries[0].lstrip())
                    exon_stop = int(boundaries[1].lstrip())
                    exon_id = str(exon_start)+str(exon_stop)
                    
                    if not exon_id in exons[chrom]:
                        exons[chrom][exon_id] = [
                            exon_start - splice_padding, 
                            exon_stop + splice_padding, 
                            exon_id
                        ]
    return genes, exons


def gene_pred_parser(ref_file_handle, splice_padding=2):
    """
    Parse a file in the refGene format.
    
    Arguments:
        ref_file_handle (file_handle): An opened file in ref_gene format
        splice_padding (int): An integer that describes how we should expand 
                              the exons
    Returns:
        genes (dict): A dictionary chromosome ids as keys and gene dictionarys
                      as values.
        exons (dict): A dictionary chromosome ids as keys and exon dictionarys
                      as values.
    """
    # A dictionary with {<chr>: {gene_id:[start, stop, gene_id]}
    genes = {}
    # A dictionary with {<chr>: {exon_id:[exon_start, exon_stop, exon_id]}}
    exons = {} 
    ### Features look like [start,stop, feature_id] ###
    
    for line in ref_file_handle:
        if not line.startswith('#') and len(line) > 1:
            line = line.split('\t')
            transcript = line[1]
            chrom = line[2].lstrip('chr')
            direction = line[3]
            transc_start = int(line[4])
            transc_stop = int(line[5])
            coding_start = int(line[6])
            coding_stop = int(line[7])
            nr_of_exons = int(line[8])
            exon_starts = [
                    int(boundary) for boundary in line[9].split(',')[:-1]
                    ]
            exon_stops = [
                    int(boundary) for boundary in line[10].split(',')[:-1]
                    ]
            gene_id = line[12]
            transcript_id = ':'.join([transcript, gene_id])
            
            if chrom not in genes:
                genes[chrom] = {}
            
            if chrom not in exons:
                exons[chrom] = {}
            
            genes[chrom][gene_id] = get_coordinates(
                        genes[chrom], 
                        transc_start, 
                        transc_stop, 
                        gene_id
                        )
            
            for i in range(len(exon_starts)):
                exon_id = str(exon_starts[i]) + str(exon_stops[i])
                exons[chrom][exon_id] = [
                    exon_starts[i] - splice_padding,
                    exon_stops[i] + splice_padding,
                    exon_id
                ]

    return genes,exons

def gff_parser(gff_file_handle, splice_padding=2, id_key='ID'):
    """
    Parse a file in the gff format.
    
    Arguments:
        gff_file_handle (file_handle): An opened file in gff file format
        splice_padding (int): An integer that describes how we should expand 
                              the exons
        id_key (str): The key that defines the gene id in the info string
    
    Returns:
        genes (dict): A dictionary chromosome ids as keys and gene dictionarys
                      as values.
        exons (dict): A dictionary chromosome ids as keys and exon dictionarys
                      as values.
    """
    genes = {}
    exons = {}

    for line in gff_file_handle:
        line = line.rstrip()
        if not line.startswith('#') and len(line) > 1:
            if line.startswith('>'):
                break
            transcript_id = ''
            gene_id = ''
            gene_name = ''
            feature_start = 0
            feature_stop = 0
            
            line = line.split('\t')
            if len(line) < 5:
                line = line.split()
            chrom = line[0].lstrip('chr')
            
            feature_type = line[2]

            if is_number(line[3]) and is_number(line[4]):
                feature_start = int(line[3])
                feature_stop = int(line[4])
            #TODO Raise exception?

            for information in line[8].split(';'):
                entry = information.split('=')
                if entry[0] == 'transcript_id':
                    transcript_id = entry[1][1:-1]
                if entry[0] == 'gene_id':
                    gene_id = entry[1][1:-1]
                if entry[0] == 'gene_name':
                    gene_name = entry[1][1:-1]
                if entry[0] == id_key:
                    gene_id = entry[1]
            
            if feature_type in ['gene', 'CDS', 'exon']:
                if chrom not in genes:
                    genes[chrom] = {}
                if chrom not in exons:
                    exons[chrom] = {}

            if feature_type in ['gene', 'CDS']:
                genes[chrom][gene_id] = get_coordinates(
                            genes[chrom],
                            feature_start,
                            feature_stop,
                            gene_id
                            )

            elif line[2] == 'exon':
                exon_id = str(feature_start) + str(feature_stop)
                exon_feature = [
                    feature_start - splice_padding,
                    feature_stop + splice_padding,
                    exon_id
                ]
                print(exon_feature)
                exons[chrom][exon_id] = exon_feature


    return genes,exons


def bed_parser(bed_file_handle, splice_padding=2):
    """
    Parse a file in the bed format.
    
    Arguments:
        bed_file_handle (file_handle): An opened file in gff file format
        splice_padding (int): An integer that describes how we should expand 
                              the exons
    
    Returns:
        genes (dict): A dictionary chromosome ids as keys and gene dictionarys
                      as values.
        exons (dict): A dictionary chromosome ids as keys and exon dictionarys
                      as values.
    
    """
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
            
            if chrom not in genes:
                genes[chrom] = {}
            if chrom not in exons:
                exons[chrom] = {}
            
            genes[chrom][gene_id] = get_coordinates(
                        genes[chrom],
                        feature_start,
                        feature_stop,
                        feature_id
                        )
            
            exon_id = str(feature_start) + str(feature_stop)
            
            exons[chrom][exon_id] = [
                feature_start - splice_padding,
                feature_stop + splice_padding,
                exon_id
                ]
    
    return genes, exons

def gtf_parser(gtf_file_handle, splice_padding=2):
   """
   Parse a file in the gtf format.
   
   Arguments:
       gtf_file_handle (file_handle): An opened file in gtf file format
       splice_padding (int): An integer that describes how we should expand 
                             the exons
   
   Returns:
       genes (dict): A dictionary chromosome ids as keys and gene dictionarys
                     as values.
       exons (dict): A dictionary chromosome ids as keys and exon dictionarys
                     as values.
   
   """
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
           feature_type = line[2]
           
           if is_number(line[3]) and is_number(line[4]):
               feature_start = int(line[3])
               feature_stop = int(line[4])

           info_field = line[8].split(';')[:-1]

           for information in info_field:
               entry = information.split()

               if entry[0] == 'transcript_id':
                   transcript_id = entry[1].strip('"')
               if entry[0] == 'gene_id':
                   gene_id = entry[1].strip('"')
               if entry[0] == 'gene_name':
                   gene_name = entry[1].strip('"')

           if feature_type == 'gene':
               
               if chrom not in genes:
                   genes[chrom] = {}
               
               genes[chrom][gene_id] = get_coordinates(
                           genes[chrom],
                           feature_start,
                           feature_stop,
                           gene_id
                           )

           elif feature_type == 'exon':
               
               if chrom not in exons:
                   exons[chrom] = {}
               
               exon_id = str(feature_start) + str(feature_stop)
               exons[chrom][exon_id] = [
                   feature_start - splice_padding,
                   feature_stop + splice_padding,
                   exon_id
                   ]
               

   return genes,exons

        
        
@click.command()
@click.argument('annotation_file',
                nargs=1,
                type=click.Path(exists=True),
)
@click.option('-t' ,'--annotation_type',
                type=click.Choice(['bed', 'ccds', 'gtf', 'gene_pred', 'gff']),
                default='gene_pred',
                help='Specify the format of the annotation file.'
)
@click.option('-v', '--verbose',
                is_flag=True,
                help='Increase output verbosity.'
)
def main(annotation_file, annotation_type, verbose):
    from genmod import logger
    from genmod.log import init_log
    init_log(logger, loglevel="DEBUG")
    print(annotation_file)
    gene_trees, exon_trees = parse_annotations(annotation_file, annotation_type)
    print(gene_trees.keys())
    print(exon_trees.keys())


if __name__ == '__main__':
    main()
