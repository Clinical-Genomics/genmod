#!/usr/bin/env python
# encoding: utf-8
"""
vcf_parser.py

Parse a file with variant info in vcf format.

Creates batches and put them into a queue object.
The batches are dictionary objects with overlapping features where the feature id:s are keys and the values are dictionarys with variants.


Batch = {feature_1_id:{variant_1_id:variant_1_info, variant_2_id: variant_2_info}, feature_2_id:... }   

Created by Måns Magnusson on 2014-03-17.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
from datetime import datetime
from codecs import open

from pprint import pprint as pp

from genmod.variants import genotype

from interval_tree import interval_tree



class VariantFileParser(object):
    """Creates parser objects for parsing variant files"""
    def __init__(self, variant_file, batch_queue, head, gene_trees, args):
        super(VariantFileParser, self).__init__()
        self.variant_file = variant_file
        self.batch_queue = batch_queue
        self.verbosity = args.verbose
        self.phased = args.phased
        self.individuals = head.individuals
        self.header_line = head.header
        self.gene_trees  = gene_trees
        self.chromosomes = []
    
    def parse(self):
        """Start the parsing"""        
        beginning = True
        batch = {}
        new_chrom = None
        current_chrom = None
        current_features = []
        haploblock_id = 1
        # Haploblocks is a dictionary with list of lists like {ind_id:[[start, stop, id],[start, stop,id],...], ...}
        haploblocks = {ind_id:[] for ind_id in self.individuals}
        nr_of_batches = 0
        # Parse the vcf file:
        with open(self.variant_file, mode='r', encoding='utf-8') as f:
            
            if self.verbosity:
                start_parsing_time = datetime.now()
                start_chrom_time = start_parsing_time
                start_twenty_time = start_parsing_time
                nr_of_variants = 0
                if self.batch_queue.full():
                    print('Queue full!!')
            
            for line in f:
                # Only metadata lines start with '#'
                if not line.startswith('#'):
                    splitted_line = line.rstrip().split('\t')
                    if len(splitted_line) < 8:
                        splitted_line = line.rstrip().split()
                    variant = self.vcf_variant(splitted_line)
                    new_chrom = variant['CHROM']
                    new_features = variant['Annotation']
                    
                    if self.verbosity:
                        nr_of_variants += 1
                        if nr_of_variants % 20000 == 0:
                            print('%s variants parsed!' % nr_of_variants)
                            print('Last 20.000 took %s to parse.\n' % str(datetime.now() - start_twenty_time))
                            start_twenty_time = datetime.now()
                    
                    # If we look at the first variant, setup boundary conditions:
                    if beginning:
                        current_features = new_features
                        beginning = False
                        # Add the variant to each of its features in a batch
                        batch = self.add_variant(batch, variant, new_features)
                        current_chrom = new_chrom
                        if self.phased:
                            haploblock_starts = {ind_id:int(variant['POS']) for ind_id in self.individuals}
                            batch['haploblocks'] = {}
                    else:
                        # If we should put the batch in the queue:
                        send = True
                        
                        if self.phased:
                            for ind_id in self.individuals:
                                #A new haploblock is indicated by '/' if the data is phased
                                if '/' in variant.get(ind_id, './.'):
                                #If call is not passed we consider it to be on same haploblock(GATK recommendations)
                                    if variant.get('FILTER', '.') == 'PASS':
                                        haploblocks[ind_id].append([haploblock_starts[ind_id], int(variant['POS']) - 1,
                                                                     str(haploblock_id)])
                                        haploblock_id += 1
                                        haploblock_starts[ind_id] = int(variant['POS'])
                        
                    # Check if we are in a space between features:
                        if len(new_features) == 0:
                            if len(current_features) == 0:
                                send = False
                    #If not check if we are in a region with overlapping genes
                        elif len(set.intersection(set(new_features),set(current_features))) > 0:
                            send = False
                        
                        if send:
                            if self.phased:
                            # Create an interval tree for each individual with the phaing intervals 
                                for ind_id in self.individuals:
                                    #Check if we have just finished an interval
                                    if haploblock_starts[ind_id] != int(variant['POS']):                                        
                                        haploblocks[ind_id].append([haploblock_starts[ind_id], int(variant['POS']) - 1, 
                                                                    str(haploblock_id)])
                                        haploblock_id += 1
                                    batch['haploblocks'][ind_id] = interval_tree.IntervalTree(haploblocks[ind_id], 
                                                        0, 1, haploblocks[ind_id][0][0], haploblocks[ind_id][-1][1])
                                haploblocks = {ind_id:[] for ind_id in self.individuals}
                            # Put the job in the queue
                            self.batch_queue.put(batch)
                            nr_of_batches += 1
                            #Reset the variables
                            current_features = new_features
                            batch = self.add_variant({}, variant, new_features)
                            batch['intervals'] = {}
                        else:
                            current_features = list(set(current_features) | set(new_features))
                            batch = self.add_variant(batch, variant, new_features) # Add variant batch

                    if new_chrom != current_chrom:
                        self.chromosomes.append(current_chrom)

                        if self.verbosity:
                            print('Chromosome %s parsed!' % current_chrom)
                            print('Time to parse chromosome %s' % str(datetime.now()-start_chrom_time))
                            current_chrom = new_chrom
                            start_chrom_time = datetime.now()
        
        self.chromosomes.append(current_chrom)
        
        if self.verbosity:
            print('Chromosome %s parsed!' % current_chrom)
            print('Time to parse chromosome %s \n' % str(datetime.now()-start_chrom_time))
            print('Variants parsed!')
            print('Time to parse variants:%s' % str(datetime.now() - start_parsing_time))
        
        if self.phased:
        # Create an interval tree for each individual with the phaing intervals
            for ind_id in self.individuals:
                #check if we have just finished an interval
                if haploblock_starts[ind_id] != int(variant['POS']):
                    intervals[ind_id].append([haploblock_starts[ind_id], int(variant['POS']) - 1, str(haploblock_id)])
                    haploblock_id += 1
                try:
                    batch['haploblocks'][ind_id] = interval_tree.IntervalTree(intervals[ind_id], 
                                        0, 1, intervals[ind_id][0][0], intervals[ind_id][-1][1])
                except IndexError:
                    pass
        
        self.batch_queue.put(batch)
        nr_of_batches += 1
        return nr_of_batches
    
    def add_variant(self, batch, variant, features):
        """Adds the variant to the proper gene(s) in the batch."""
        variant_id = [variant['CHROM'], variant['POS'], variant['REF'], variant['ALT']]
        variant_id = '_'.join(variant_id)
        # If we are in a region between features:
        if len(features) == 0:
            if len(batch) == 0:
                batch['-'] = {variant_id:variant}
            else:
                batch['-'][variant_id] = variant
        for feature in features:
            if feature in batch:
                batch[feature][variant_id] = variant
            else:
                batch[feature] = {variant_id:variant}
        return batch
    
    def vcf_variant(self, splitted_variant_line):
        """Returns a variant object in the cmms format."""
        my_variant = dict(zip(self.header_line, splitted_variant_line))
        variant_chrom = my_variant['CHROM']
        alternatives = my_variant['ALT'].split(',')
        
        # When checking what features that are overlapped we use the longest alternative
        longest_alt = max([len(alternative) for alternative in alternatives])
        
        # Internally we never use 'chr' in the chromosome names:
        variant_chrom = variant_chrom.lstrip('chr')
        
        variant_interval = [int(my_variant['POS']), (int(my_variant['POS']) + 
                            longest_alt - 1)]
        
        # if self.exon_annotation:
        #     try:
        #         my_variant['Annotation'] = self.gene_trees.exon_trees[variant_chrom].find_range(variant_interval)
        #     except KeyError:
        #         if self.verbosity:
        #             print('Chromosome', variant_chrom, 'is not in annotation file!')
        #         my_variant['Annotation'] = []
        
        # The feature files does not have to include all chromosomes that are in the vcf:
        try:
            my_variant['Annotation'] = self.gene_trees.gene_trees[variant_chrom].find_range(variant_interval)
        except KeyError:
            if self.verbosity:
                print('Chromosome', variant_chrom, 'is not in annotation file!')
            my_variant['Annotation'] = []
        
        return my_variant

def main():
    from multiprocessing import JoinableQueue
    from genmod.vcf import vcf_header
    from genmod.utils import annotation_parser
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    parser.add_argument('annotation_file', type=str, nargs=1 , help='A file with feature annotations.')
    parser.add_argument('-phased', '--phased', action="store_true", help='If variant file is phased.')    
    parser.add_argument('-v', '--verbose', action="store_true", help='Increase output verbosity.')
    
    args = parser.parse_args()
    infile = args.variant_file[0]
    if args.verbose:
        print('Parsing annotationfile...')
        start_time_annotation = datetime.now()
    my_anno_parser = annotation_parser.AnnotationParser(args.annotation_file[0], 'ref_gene')
    
    if args.verbose:
        print('annotation parsed. Time to parse annotation: %s\n' % str(datetime.now() - start_time_annotation))
    
    my_head_parser = vcf_header.VCFParser(infile)
    my_head_parser.parse()
    print(my_head_parser.__dict__)
    variant_queue = JoinableQueue()
    start_time = datetime.now()        
    
    my_parser = VariantFileParser(infile, variant_queue, my_head_parser, my_anno_parser, args)
    nr_of_batches = my_parser.parse()
    print(nr_of_batches)
    for i in range(nr_of_batches):
        variant_queue.get()
        variant_queue.task_done()
    
    variant_queue.join()
    print('Time to parse variants: %s ' % str(datetime.now()-start_time))

if __name__ == '__main__':
    main()
