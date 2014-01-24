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
from datetime import datetime

from genmod.variants import genetic_variant, genotype
from genmod.utils import get_genes
from collections import OrderedDict



class VariantParser(object):
    """docstring for VariantParser"""
    def __init__(self, variant_file, tasks_queue, individuals, verbosity):
        super(VariantParser, self).__init__()
        self.variant_file = variant_file
        self.tasks_queue = tasks_queue
        self.individuals = individuals
        self.header_line = []
        self.metadata = []
        self.verbosity = verbosity    
        variants = []
        current_chrom = '0'
        new_chrom = '0'
        chrom_time = datetime.now()
    
        # We only parse the variants here:
        with open(variant_file, 'rb') as f:
            
            beginning = True
            batch = {} # This is a dictionary to store the variant lines of a batch
            current_genes = []  # These are lists to keep track of the regions that we look at
            new_region = []
            for line in f:
                line = line.rstrip()
                if line[0] != '##':
                    #These are variant lines                    
                    variant = self.vcf_variant(line)
                    if self.verbosity:
                        new_chrom = variant.chr
                        if new_chrom != current_chrom:
                            print 'Now processing chromosome', new_chrom
                            print 'Time to process chromosome', current_chrom, ':', datetime.now() - chrom_time
                            chrom_time = datetime.now()
                            current_chrom = new_chrom
                    
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
    
    def vcf_variant(self, variant_line):
        """Read a VCF-variant line and returns a variant object and a 
        dictionary with the genotypes like {variant-id:[genotype1, genotype2,...]}"""
        variant_line = variant_line.split()
        format_info = [] # Information about the genotype format
        variant_info = OrderedDict()

        for entry in range(len(variant_line)):
            variant_info[self.header_lines[entry]] = variant_line[entry]
            
            if entry > 8:
                individual = self.header_lines[entry]
                if individual not in self.individuals:
                    self.individuals[individual] = {}
        
        chrom = variant_info['CHROM']
        position = int(variant_info['POS'])
        identity = variant_info['ID']
        reference = variant_info['REF']
        alternatives = variant_info['ALT']
        quality = variant_info['QUAL']
        filt = variant_info['FILTER']
        info = variant_info['INFO']
        format_info = variant_info['FORMAT'].split(':')
        number_of_individuals = len(variant_line[9:])
        
        if self.split_alternatives:
    
            for alternative in alternatives:
                # This is a SNV:
                if len(reference) == len(alternative):
                    start = position
                    stop = position
                # This is a deletion:
                elif len(reference) > len(alternative):
                    start = str(int(position) + 1)
                    stop = str(int(position) + len(reference) - 1)
                # This is a insertion
                elif len(reference) < len(alternative):
                    start = position
                    stop = position
                    alternative = alternative[1:]
                # Put the variants in a list:
                my_variant = genetic_variant.Variant(chrom, start, stop, reference, alternative, identity, variant_info)
                self.variants[my_variant.variant_id] = my_variant
                # Check the genotypes
        else:
            
            start = position
            stop = position
            my_variant = genetic_variant.Variant(chrom, start, stop, reference, alternatives, identity, variant_info)
            self.variants[my_variant.variant_id] = my_variant
        
        # Collect the genotypes:
        
        for individual in self.individuals:
                genotype_arguments = {} # args for Genotype class
                individual_info = variant_info[individual].split(':')
                for i in range(len(format_info)):
                # Fill the dictionary like {GT:'0/1', DP:'10'} and so on
                    if not i > len(individual_info)-1:
                        genotype_arguments[format_info[i]] = individual_info[i]
                my_genotype = genotype.Genotype(GT=genotype_arguments.get('GT','./.'), AD=genotype_arguments.get('AD','.,.'), DP=genotype_arguments.get('DP','0'), GQ=genotype_arguments.get('GQ','0'))
                
                # Add the genotypes to each individual:
                
                self.individuals[individual][my_variant.variant_id] = my_genotype


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
