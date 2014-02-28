#!/usr/bin/env python
# encoding: utf-8
"""
variant_parser.py


Parse a file with variant info, this can be a .vcf file, an annotated annovar file, 
a annotated .txt cmms file, a annotated .txt cmms_ranked .

Create a variant objects and a dictionary with individuals that have a dictionary with genotypes for each variant.

Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
from datetime import datetime
from pysam import Tabixfile, asTuple
import sqlite3

from pprint import pprint as pp

from genmod.variants import genotype


class VariantFileParser(object):
    """docstring for VariantParser"""
    def __init__(self, variant_file, batch_queue, head, interval_trees, cadd_db = None, cadd_file = None, verbosity = False):
        super(VariantFileParser, self).__init__()
        self.variant_file = variant_file
        self.batch_queue = batch_queue
        self.verbosity = verbosity
        self.individuals = head.individuals
        self.header_line = head.header
        self.interval_trees = interval_trees
        self.chromosomes = []
        self.cadd_db = None
        if cadd_db:
            self.conn = sqlite3.connect(cadd_db)
            self.cadd_db = self.conn.cursor() 
    
    def parse(self):
        """Start the parsing"""        
        start_parsing = datetime.now()
        start_chrom = start_parsing
        start_twenty = start_parsing
        beginning = True
        batch = {}
        new_chrom = None
        current_chrom = None
        current_features = []
        nr_of_variants = 0
        nr_of_batches = 0
        with open(self.variant_file, 'r') as f:
            for line in f:
                # Variant lines do not start with '#'
                if not line.startswith('#'):
                    variant, new_features = self.vcf_variant(line.rstrip().split('\t'))
                    new_chrom = variant['CHROM']
                    if self.verbosity:
                        nr_of_variants += 1
                        if nr_of_variants % 20000 == 0:
                            print('%s variants parsed!' % nr_of_variants)
                            print('Last 20.000 took %s to parse.\n' % str(datetime.now() - start_twenty))
                            start_twenty = datetime.now()
                    # If we look at the first variant, setup boundary conditions:
                    if beginning:
                        current_features = new_features
                        beginning = False
                        # Add the variant to each of its features in a batch
                        batch = self.add_variant(batch, variant, new_features)
                        current_chrom = new_chrom
                    else:
                        send = True
                    # Check if we are in a space between features:
                        if len(new_features) == 0:
                            if len(current_features) == 0:
                                send = False
                    #If not check if we are in a region with overlapping genes
                        elif len(set.intersection(set(new_features),set(current_features))) > 0:
                            send = False
                        
                        if send:
                            nr_of_batches += 1
                            self.batch_queue.put(batch)
                            current_features = new_features
                            batch = self.add_variant({}, variant, new_features)
                        else:
                            current_features = list(set(current_features) | set(new_features))
                            batch = self.add_variant(batch, variant, new_features) # Add variant batch

                    if new_chrom != current_chrom:
                        self.chromosomes.append(current_chrom)

                        if self.verbosity:
                            print('Chromosome %s parsed!' % current_chrom)
                            print('Time to parse chromosome %s' % str(datetime.now()-start_chrom))
                            current_chrom = new_chrom
                            start_chrom = datetime.now()
        
        self.chromosomes.append(current_chrom)
        if self.verbosity:
            print('Chromosome %s parsed!' % current_chrom)
            print('Time to parse chromosome %s \n' % str(datetime.now()-start_chrom))
            print('Variants parsed!')
            print('Time to parse variants:%s' % str(datetime.now() - start_parsing))
        nr_of_batches += 1
        self.batch_queue.put(batch)
        return nr_of_batches
    
    def add_variant(self, batch, variant, features):
        """Adds the variant to the proper gene(s) in the batch."""
        variant_id = [variant['CHROM'], variant['POS'], variant['REF'], variant['ALT']]
        variant_id = '_'.join(variant_id)
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
        longest_alt = max([len(alternative) for alternative in alternatives])
        
        if variant_chrom.startswith('chr') or variant_chrom.startswith('Chr'):
            variant_chrom = variant_chrom[3:]
        
        variant_interval = [int(my_variant['POS']), (int(my_variant['POS']) + 
                            longest_alt - 1)]
        features_overlapped = []
        
        try:
            features_overlapped = self.interval_trees.interval_trees[variant_chrom].findRange(variant_interval)
        except KeyError:
            if self.verbosity:
                print('Chromosome', variant_chrom, 'is not in annotation file!')

        my_variant['Annotation'] = features_overlapped

        
        cadd_score = 0
        if self.cadd_db:
            #Only check for snps:
            if longest_alt == 1 and len(my_variant['REF']) == 1:
                cadd_key = '_'.join([variant_chrom, my_variant['REF']])
                try:
                    cad_score = self.cadd_db.execute(reference=variant_chrom, start=int(my_variant['POS'])-1, end=int(my_variant['POS'])):
                        if my_variant['ALT'].split(',')[0] == tpl[3]:
        #                     cad_score = float(tpl[5])
        #         except KeyError:
        #             pass
        # my_variant['CADD'] = cad_score
        # if self.verbosity:
        #     if cad_score > 30:
        #         print('High sore! Cadd score %s' % str(cad_score))
        #         print('%s \n' % '\t'.join([key+':'+my_variant[key] for key in self.header_line]))
                
        return my_variant, features_overlapped

def main():
    from multiprocessing import JoinableQueue
    from genmod.vcf import vcf_header
    from genmod.utils import annotation_parser
    from pysam import tabix_index
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    parser.add_argument('annotation_file', type=str, nargs=1 , help='A file with feature annotations.')
    parser.add_argument('-cadd_db', '--cadd_db', type=str, nargs=1 , default=[None], help='A sqlite db with cadd values.')
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
    
    if args.cadd_file[0]:
        try:
            tabix_index(args.cadd_file[0], seq_col=0, start_col=1, end_col=1, meta_char='#')
        except IOError as e:
            if args.verbose:
                print e
            pass
        
    
    my_parser = VariantFileParser(infile, variant_queue, my_head_parser, my_anno_parser, args.cadd_file[0], args.verbose)
    nr_of_batches = my_parser.parse()
    for i in xrange(nr_of_batches):
        variant_queue.get()
        variant_queue.task_done()
    
    variant_queue.join()
    print('Time to parse variants: %s ' % str(datetime.now()-start_time))

if __name__ == '__main__':
    main()
