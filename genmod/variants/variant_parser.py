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
import tempfile
from genmod.variants import genetic_variant, genotype
from genmod.utils import get_genes
from collections import OrderedDict



class VariantParser(object):
    """docstring for VariantParser"""
    def __init__(self, variant_file, tasks_queue, individuals):
        super(VariantParser, self).__init__()
        self.variant_file = variant_file
        self.tasks_queue = tasks_queue
        self.individuals = individuals
        self.header_line = []
        self.metadata = []
            
        variants = []
    
        with open(variant_file, 'r') as f:
            
            beginning = True
            batch = {} # This is a dictionary to store the variant lines of a batch
            current_genes = []  # These are lists to keep track of the regions that we look at
            new_region = []
            for line in f:
                line = line.rstrip()
                if line[:2] == '##':
                #This is the metadata information
                    self.metadata.append(line)
                elif line[:1] == '#':
                    self.header_line = line[1:].split('\t')
                else:
                    #These are variant lines                    
                    variant = self.cmms_variant(line)
                    new_genes = variant.genes
                    # If we look at the first variant, setup boundary conditions:
                    if beginning:
                        current_genes = new_genes
                        beginning = False
                        batch = self.add_variant(batch, variant) # Add variant batch

                    else:
                        send = True
                    
                    #Check if we are in a space between genes:
                        if len(new_genes) == 0:
                            if len(current_genes) == 0:
                                send = False
                    #If not check if we are in a consecutive region
                        else:
                            for gene in new_genes:
                                if gene in current_genes:
                                    send = False
                        if send:
                            # If there is an intergenetic region we do not look at the compounds.
                            # The tasks are tuples like (variant_list, bool(if compounds))
                            self.tasks_queue.put(batch)
                            current_genes = new_genes
                            batch = self.add_variant({}, variant)
                        else:
                            current_genes = list(set(current_genes) | set(new_genes))
                            batch = self.add_variant(batch, variant) # Add variant batch
        tasks_queue.put(batch)
        
    def add_variant(self, batch, variant):
        """Adds the variant to the proper gene(s) in the batch."""
        for gene in variant.genes:
            if gene in batch:
                batch[gene][variant.variant_id] = variant
            else:
                batch[gene] = {variant.variant_id:variant}
        return batch
    
    def cmms_variant(self, variant_line):
        """Returns a variant object in the cmms format."""
        
        variant_info = OrderedDict()
        individual_genotypes = {} #DICT with {ind_id:{<genotype>}}
        counter = 0
        variant_line = variant_line.split('\t')
        ensemble_entry = variant_line[5]
        hgnc_entry = variant_line[6]
        
        # These must be parsed separately
        hgnc_genes = get_genes.get_genes(hgnc_entry, 'HGNC')
        ensemble_genes = get_genes.get_genes(ensemble_entry, 'Ensemble')
    
        for entry in range(len(variant_line)):
            
            if 'IDN' in self.header_line[entry]:
                # Looks like IDN:11-1-2A
                individual = self.header_line[entry].split(':')[-1]
                if individual not in self.individuals:
                    self.individuals.append(individual)
            variant_info[self.header_line[entry]] = variant_line[entry]
                
        chrom = variant_info['Chromosome']
        start = variant_info['Variant_start']
        stop = variant_info['Variant_stop']
        alternative = variant_info['Alternative_allele']
        reference = variant_info['Reference_allele']
        identity = variant_info['Dbsnp_rs_nr']
        my_variant = genetic_variant.Variant(chrom, start, stop, reference, alternative, identity, genes=hgnc_genes, all_info=variant_info)
    
        # Add the genotypes to variant:
        
        for individual in self.individuals:
            genotype_arguments = {} # args for genotype class
            key = 'IDN:' + individual
            # gt_info looks like 11-1-2A:GT=0/1:PL=32,3,2:...
            for gt_info in variant_info[key].split(':')[1:]:
                value_pair = gt_info.split('=')
                genotype_arguments[value_pair[0]] = value_pair[-1]
            my_genotype = genotype.Genotype(GT=genotype_arguments.get('GT','./.'), AD=genotype_arguments.get('AD','.,.'), DP=genotype_arguments.get('DP','0'), GQ=genotype_arguments.get('GQ','0'))
            my_variant.genotypes[individual] = my_genotype
   
        return my_variant


def main():
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    args = parser.parse_args()
    infile = args.variant_file[0]
    file_type = 'cmms'
    my_parser = VariantParser(infile, file_type)
    number_of_variants = 0
    for chrom in my_parser.chrom_shelves:
        current_db = shelve.open(my_parser.chrom_shelves[chrom])
        print chrom, len(current_db)
        for variant in current_db:
            print current_db[variant]
        number_of_variants += len(current_db)
        current_db.close()
        os.remove(my_parser.chrom_shelves[chrom])
    print 'Number of shelved:', number_of_variants


if __name__ == '__main__':
    main()
