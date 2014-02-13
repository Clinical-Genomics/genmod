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
from tempfile import gettempdir


if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict
    
from pprint import pprint as pp

from genmod.variants import genotype


class VariantFileParser(object):
    """docstring for VariantParser"""
    def __init__(self, variant_file, batch_queue, head, interval_tree, verbosity = False):
        super(VariantFileParser, self).__init__()
        self.variant_file = variant_file
        self.batch_queue = batch_queue
        self.verbosity = verbosity
        self.individuals = head.individuals
        self.header_line = head.header
        self.interval_tree = interval_tree
        self.chromosomes = []
    
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
        temp_dir = gettempdir()
        with open(self.variant_file, 'rb') as f:
            for line in f:
                
                if not line.startswith('#'):
                    variant, new_features = self.vcf_variant(line.rstrip().split('\t'))
                    new_chrom = variant['CHROM']
                    if self.verbosity:
                        nr_of_variants += 1
                        if nr_of_variants % 20000 == 0:
                            print nr_of_variants, 'variants parsed!'
                            print 'Last 20.000 took', datetime.now() - start_twenty, 'to parse.'
                            print ''
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
                        # print current_features, new_features
                        if len(new_features) == 0:
                            if len(current_features) == 0:
                                send = False
                    #If not check if we are in a consecutive region
                        elif len(set.intersection(set(new_features),set(current_features))) > 0:
                            send = False
                        
                        if send:
                            # If there is an intergenetic region we do not look at the compounds.
                            # The tasks are tuples like (variant_list, bool(if compounds))
                            self.batch_queue.put(batch)
                            current_features = new_features
                            batch = self.add_variant({}, variant, new_features)
                        else:
                            current_features = list(set(current_features) | set(new_features))
                            batch = self.add_variant(batch, variant, new_features) # Add variant batch

                    if new_chrom != current_chrom:
                        self.chromosomes.append(current_chrom)

                        if self.verbosity:
                            print 'Chromosome', current_chrom, 'parsed!'
                            print 'Time to parse chromosome', datetime.now()-start_chrom
                            current_chrom = new_chrom
                            start_chrom = datetime.now()
        
        self.chromosomes.append(current_chrom)
        if self.verbosity:
            print 'Chromosome', current_chrom, 'parsed!'
            print 'Time to parse chromosome', datetime.now()-start_chrom
            print ''
            print 'Variants parsed!'
            print 'Time to parse variants:', datetime.now() - start_parsing
        self.batch_queue.put(batch)
        return
    
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
        my_variant = OrderedDict(zip(self.header_line, splitted_variant_line))
        variant_chrom = my_variant['CHROM']
        variant_interval = [int(my_variant['POS']), int(my_variant['POS'])]
        features_overlapped = []
        
        try:
            features_overlapped = self.interval_tree.interval_trees[variant_chrom].findRange(variant_interval)
        except KeyError:
            if self.verbosity:
                print 'Chromosome', variant_chrom, 'is not in annotation file!'
        my_variant['Annotation'] = features_overlapped
        return my_variant, features_overlapped

def main():
    from multiprocessing import JoinableQueue
    from genmod.vcf import vcf_header
    from genmod.utils import annotation_parser
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    parser.add_argument('annotation_file', type=str, nargs=1 , help='A file with feature annotations.')
    parser.add_argument('-v', '--verbose', action="store_true", help='Increase output verbosity.')
    
    args = parser.parse_args()
    infile = args.variant_file[0]
    anno_file = args.annotation_file[0]
    my_anno_parser = annotation_parser.AnnotationParser(anno_file, 'ref_gene')
    for tree in my_anno_parser.interval_trees:
        print tree
    my_head_parser = vcf_header.VCFParser(infile)
    my_head_parser.parse()
    print my_head_parser.__dict__
    variant_queue = JoinableQueue()
    start_time = datetime.now()
    my_parser = VariantFileParser(infile, variant_queue, my_head_parser,my_anno_parser, args.verbose)
    my_parser.parse()

    print 'Time to parse variants:', datetime.now()-start_time

if __name__ == '__main__':
    main()
