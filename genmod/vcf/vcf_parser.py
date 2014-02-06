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
import shelve
from datetime import datetime
if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict
    
from pprint import pprint as pp

from genmod.variants import genotype
from genmod.utils import get_genes



class VariantFileParser(object):
    """docstring for VariantParser"""
    def __init__(self, variant_file, batch_queue, head, verbosity = False):
        super(VariantFileParser, self).__init__()
        self.variant_file = variant_file
        self.batch_queue = batch_queue
        self.verbosity = verbosity
        self.individuals = head.individuals
        self.header_line = head.header
    
    def parse(self):
        """Start the parsing"""        
        start_parsing = datetime.now()
        start_chrom = start_parsing
        start_twenty = start_parsing
        beginning = True
        batch = {}
        new_chrom = None
        current_chrom = None
        current_genes = []
        nr_of_variants = 0
        with open(self.variant_file, 'rb') as f:
            for line in f:
                
                if not line.startswith('#'):
                    variant, new_genes = self.vcf_variant(line.rstrip().split('\t'))
                    if self.verbosity:
                        nr_of_variants += 1
                        new_chrom = variant['CHROM']
                        if nr_of_variants % 20000 == 0:
                            print nr_of_variants, 'variants parsed!'
                            print 'Last 20.000 took', datetime.now() - start_twenty, 'to parse.'
                            print ''
                            start_twenty = datetime.now()
                    # If we look at the first variant, setup boundary conditions:
        #             if beginning:
        #                 current_genes = new_genes
        #                 beginning = False
        #                 # Add the variant to each of its genes in a batch
        #                 batch = self.add_variant(batch, variant, new_genes)
        #                 if self.verbosity:
        #                     current_chrom = new_chrom
        #             else:
        #                 send = True
        #             
        #             # Check if we are in a space between genes:
        #                 # print current_genes, new_genes
        #                 if len(new_genes) == 0:
        #                     if len(current_genes) == 0:
        #                         send = False
        #             #If not check if we are in a consecutive region
        #                 elif len(set.intersection(set(new_genes),set(current_genes))) > 0:
        #                     send = False
        #                 
        #                 if send:
        #                     # If there is an intergenetic region we do not look at the compounds.
        #                     # The tasks are tuples like (variant_list, bool(if compounds))
        #                     self.batch_queue.put(batch)
        #                     current_genes = new_genes
        #                     batch = self.add_variant({}, variant, new_genes)
        #                 else:
        #                     current_genes = list(set(current_genes) | set(new_genes))
        #                     batch = self.add_variant(batch, variant, new_genes) # Add variant batch
        #             
        #             if self.verbosity:
        #                 if new_chrom != current_chrom:
        #                     print 'Chromosome', current_chrom, 'parsed!'
        #                     print 'Time to parse chromosome', datetime.now()-start_chrom
        #                     current_chrom = new_chrom
        #                     start_chrom = datetime.now()
        #                 
        # if self.verbosity:
        #     print 'Chromosome', current_chrom, 'parsed!'
        #     print 'Time to parse chromosome', datetime.now()-start_chrom
        # self.batch_queue.put(batch)
        return
    
    def add_variant(self, batch, variant, genes):
        """Adds the variant to the proper gene(s) in the batch."""
        variant_id = [variant['Chromosome'], variant['Variant_start'], variant['Reference_allele'], variant['Alternative_allele']]
        variant_id = '_'.join(variant_id)
        if len(genes) == 0:
            if len(batch) == 0:
                batch['-'] = {variant_id:variant}
            else:
                batch['-'][variant_id] = variant
        for gene in genes:
            if gene in batch:
                batch[gene][variant_id] = variant
            else:
                batch[gene] = {variant_id:variant}
        return batch
    
    
    def vcf_variant(self, splitted_variant_line):
        """Returns a variant object in the cmms format."""
    
        # These must be parsed separately
        hgnc_genes = []
        # ensemble_genes = get_genes.get_genes(ensemble_entry, 'Ensemble')
        
        return OrderedDict(zip(self.header_line, splitted_variant_line)), hgnc_genes
    
        # for entry in range(len(splitted_variant_line)):
        #     
        #     if 'IDN' in self.header_line[entry]:
        #         # Looks like IDN:11-1-2A
        #         gt_info = splitted_variant_line[entry].split(':')
        #         individual_genotypes[gt_info[0]] = genotype.Genotype(GT=gt_info[1].split('=')[1])  
        #         individual = gt_info[0]
        #         if individual not in self.individuals:
        #             raise SyntaxError('One of the individuals in the variant file \
        #                         is not in the ped file: %s' % individual)
        #     variant_info[self.header_line[entry]] = splitted_variant_line[entry]
    
        # chrom = variant_info['Chromosome']
        # start = variant_info['Variant_start']
        # stop = variant_info['Variant_stop']
        # alternative = variant_info['Alternative_allele']
        # reference = variant_info['Reference_allele']
        # identity = variant_info['Dbsnp_rs_nr']
        # my_variant = genetic_variant.Variant(chrom , start, stop, reference, 
        #                 alternative, identity, genes=hgnc_genes, all_info=variant_info)
        #     
        # Add the genotypes to variant:
            
        # for individual in self.individuals:
        #     genotype_arguments = {} # args for genotype class
        #     key = 'IDN:' + individual
        #     # gt_info looks like 11-1-2A:GT=0/1:PL=32,3,2:...
        #     for gt_info in variant_info[key].split(':')[1:]:
        #         value_pair = gt_info.split('=')
        #         genotype_arguments[value_pair[0]] = value_pair[-1]
        #     my_genotype = genotype.Genotype(GT=genotype_arguments.get('GT','./.'), 
        #                                     AD=genotype_arguments.get('AD','.,.'), 
        #                                     DP=genotype_arguments.get('DP','0'), 
        #                                     GQ=genotype_arguments.get('GQ','0'))
        #     my_variant.genotypes[individual] = my_genotype
        # return variant_info, individual_genotypes, hgnc_genes
        # return my_variant


def main():
    from multiprocessing import JoinableQueue
    from genmod.vcf import vcf_header
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    parser.add_argument('-v', '--verbose', action="store_true", help='Increase output verbosity.')
    
    args = parser.parse_args()
    infile = args.variant_file[0]
    head = vcf_header.VCFParser(infile)
    variant_queue = JoinableQueue()
    start_time = datetime.now()
    my_parser = VariantFileParser(infile, variant_queue, head, args.verbose)

    print 'Time to parse variants:', datetime.now()-start_time

if __name__ == '__main__':
    main()
