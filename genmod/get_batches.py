#!/usr/bin/env python
# encoding: utf-8
"""
get_batches.py

Walks through variants and put batches of variants in a queue object.
A batch is a dictionary with variants that resides in the same feature(gene).
If genes are overlapping in a region then they will be considered as one batch.
An intergenic regions will become a batch. If the number of variants in an 
intergenic region is higher than 10000 it will be splitted in two batches.

Batch = {variant_id: variant_info, variant_2_id:variant_2_info, ... }   


This will be three batches, one for each gene 
and one for the space between them:

---------------------------------------------------------
|       Gene1                           Gene2           |
| --------------------        ----------------------    |
|                                                       |
|                                                       |
---------------------------------------------------------


This will be 4 bathces one for Gene1+Gene2, 
one for Gene3 and one for the space between them:

---------------------------------------------------------
|       Gene1                          Gene3            |
| --------------------             ---------------      |
|                Gene2                                  |
|        ----------------------                         |
---------------------------------------------------------


The variants must have region annotation otherwise the batches can not be
created. If they are annotated with VEP, this annotation will be used.
If that is not the case VariantParser will use an annotation built by genmod.



Created by Måns Magnusson on 2014-03-17.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

from __future__ import (print_function, unicode_literals)

import sys
import os
import argparse
from datetime import datetime
from codecs import open
  
from pprint import pprint as pp

from interval_tree import interval_tree
from genmod.errors import warning


# These are the SO-terms for genetic variants used by VEP:
INTERESTING_SO_TERMS = set(
            [
                'transcript_ablation',
                'splice_donor_variant',
                'splice_acceptor_variant',
                'stop_gained',
                'frameshift_variant',
                'stop_lost',
                'initiator_codon_variant',
                'transcript_amplification',
                'inframe_insertion',
                'inframe_deletion',
                'missense_variant',
                'splice_region_variant',
                'incomplete_terminal_codon_variant',
                'stop_retained_variant',
                'synonymous_variant',
                'coding_sequence_variant',
                'mature_miRNA_variant',
                '5_prime_UTR_variant',
                '3_prime_UTR_variant',
                'non_coding_exon_variant',
                'non_coding_transcript_exon_variant',
                'non_coding_transcript_variant',
                'nc_transcript_variant',
                'intron_variant',
                'NMD_transcript_variant',
                'non_coding_transcript_variant',
            ]
)

def get_batches(variant_parser, batch_queue, individuals, gene_trees={}, 
                exon_trees={}, phased=False, vep=False, whole_genes=False, 
                verbosity=False):
    """
    Create batches and put them into the queue.
    Annotate the variants with regions, either from the annotation built by
    genmod or check the VEP terms. The variants in one feature will be a 
    batch(default feature is a gene), if intergenic the batch sixe is 
    maximun 10000 variants long. After one batch is filled it is sent 
    to the variant queue.
    """
    beginning = True
    # A batch is a dictionary with variants
    batch = {}
    new_chrom = None
    current_chrom = None
    current_features = []
    haploblock_id = 1
    # Haploblocks is a dictionary with list of lists like {ind_id:[[start, stop, id],[start, stop,id],...], ...}
    haploblocks = {ind_id:[] for ind_id in individuals}
    nr_of_batches = 0
    chromosomes = []
    # Parse the vcf file:
    if verbosity:
        start_parsing_time = datetime.now()
        start_chrom_time = start_parsing_time
        start_twenty_time = start_parsing_time
        if batch_queue.full():
            warning('Queue full!!')
    
    nr_of_variants = 0
    for variant in variant_parser:
        
        variant_id = variant['variant_id']
        nr_of_variants += 1
        new_chrom = variant['CHROM']
        if new_chrom.startswith('chr'):
            new_chrom = new_chrom[3:]
        
        # Annotate which features the variant belongs to:
        annotate_variant(
                            variant, 
                            gene_trees, 
                            exon_trees, 
                            vep, 
                            whole_genes, 
                            verbosity
                        )
        
        new_features = variant['annotation']
        
        if verbosity:
            if nr_of_variants % 20000 == 0:
                print('%s variants parsed!' % nr_of_variants, file=sys.stderr)
                print('Last 20.000 took %s to parse.\n' % 
                        str(datetime.now() - start_twenty_time), file=sys.stderr)
                start_twenty_time = datetime.now()
        
        # If we look at the first variant, setup boundary conditions:
        if beginning:
            current_features = new_features
            # Add the variant to each of its features in a batch
            batch[variant_id] = variant
            current_chrom = new_chrom
            batch['haploblocks'] = {}
            if phased:
                # We collect the starts of the haploblocks
                haploblock_starts = {ind_id:int(variant['POS']) for ind_id in individuals}
            beginning = False
        else:
            # If we should put the batch in the queue:
            send = True
            
            if phased:
                for ind_id in individuals:
                    #A new haploblock is indicated by '/' if the data is phased
                    if '/' in variant.get(ind_id, './.'):
                    #If call is not passed we consider it to be on same haploblock(GATK recommendations)
                        if variant.get('FILTER', '.') == 'PASS':
                            haploblocks[ind_id].append(
                                                [   
                                                    haploblock_starts[ind_id], 
                                                    int(variant['POS']) - 1,
                                                    str(haploblock_id)
                                                ]
                                            )
                            haploblock_id += 1
                            haploblock_starts[ind_id] = int(variant['POS'])
            
        # Check if we are in a space between features:
            if len(new_features) == 0:
                if len(current_features) == 0:
                    # If the intergeneic region is bigger than 10000 we send it as a batch
                    if len(batch) < 10000:
                        send = False
        #If not check if we are in a region with overlapping features
            elif new_features.intersection(current_features):
                send = False
            
        # If we are at a new chromosome we finish the current batch:
            if new_chrom != current_chrom:
                chromosomes.append(current_chrom)
                # New chromosome means new batch
                send = True
                
                if verbosity:
                    print('Chromosome %s parsed!' % current_chrom, file=sys.stderr)
                    print('Time to parse chromosome %s' % 
                            str(datetime.now()-start_chrom_time), 
                            file=sys.stderr
                            )
                    start_chrom_time = datetime.now()
            
                current_chrom = new_chrom
            
            if send:
                if phased:
                # Create an interval tree for each individual with the phaing intervals 
                    for ind_id in individuals:
                        #Check if we have just finished an interval
                        if haploblock_starts[ind_id] != int(variant['POS']):                                        
                            haploblocks[ind_id].append(
                                                [
                                                    haploblock_starts[ind_id], 
                                                    int(variant['POS']), 
                                                    str(haploblock_id)
                                                ]
                                            )
                            haploblock_id += 1
                        # Create interval trees of the haploblocks
                        batch['haploblocks'][ind_id] = interval_tree.IntervalTree(
                                                                haploblocks[ind_id], 
                                                                haploblocks[ind_id][0][0]-1, 
                                                                haploblocks[ind_id][-1][1]+1
                                                        )
                    haploblocks = {ind_id:[] for ind_id in individuals}
                
                # Put the job in the queue
                batch_queue.put(batch)
                nr_of_batches += 1
                #Reset the variables
                current_features = new_features
                batch = {}
                batch[variant_id] = variant
                batch['haploblocks'] = {}
            else:
                current_features = current_features.union(new_features)
                batch[variant_id] = variant
    
    chromosomes.append(current_chrom)
    nr_of_batches += 1
    
    if verbosity:
        print('Chromosome %s parsed!' % current_chrom, file=sys.stderr)
        print('Time to parse chromosome %s \n' % str(datetime.now()-start_chrom_time), file=sys.stderr)
        print('Variants parsed!', file=sys.stderr)
        print('Time to parse variants:%s' % str(datetime.now() - start_parsing_time), file=sys.stderr)
        print('Number of variants in variant file:%s\n' % nr_of_variants, file=sys.stderr)
        print('Number of batches created:%s\n' % nr_of_batches, file=sys.stderr)
    
    if phased:
    # Create an interval tree for each individual with the phasing intervals
        for ind_id in individuals:
            #check if we have just finished an interval
            if haploblock_starts[ind_id] != int(variant['POS']):
                haploblocks[ind_id].append(
                                    [
                                        haploblock_starts[ind_id], 
                                        int(variant['POS']), 
                                        str(haploblock_id)
                                    ]
                                )
                haploblock_id += 1
            try:
                batch['haploblocks'][ind_id] = interval_tree.IntervalTree(
                                                        haploblocks[ind_id], 
                                                        haploblocks[ind_id][0][0]-1, 
                                                        haploblocks[ind_id][-1][1]+1
                                                )
            except IndexError:
                pass
    
    batch_queue.put(batch)
    
    return chromosomes

def check_vep_annotation(variant):
    """
    Return a set with the genes that vep has annotated this variant with.
    
    Input: A variant
    
    Returns: A set with genes"""
    
    annotation = set()
    # vep_info is a dictionary with genes as key and annotation as values
    
    for allele in variant.get('vep_info',{}):
        if allele != 'gene_ids':
            for vep_annotation in variant['vep_info'][allele]:
                for consequence in vep_annotation.get('Consequence', '').split('&'):
                    if consequence in INTERESTING_SO_TERMS:
                        annotation.add(vep_annotation.get('SYMBOL', ''))
    return annotation
        
    
def annotate_variant(variant, gene_trees, exon_trees, vep, whole_genes, verbosity):
    """
    Annotate variants with what regions the belong.
    Adds 'annotation' = set(set, of, genes) and 
    'compound_candidate' = Boolean to variant dictionary.
    Variants are compound candidates is the are exonic
    and in the same gene. 
    If 'while_gene' is used intronic variants are also
    compound candidates.
    
    Input: variant_dictionary
    
    Returns: variant_dictionary with annotation added    
    """
        
    variant['comp_candidate'] = False
    variant['annotation'] = set()
    
    # Internally we never use 'chr' in the chromosome names:        
    chrom = variant['CHROM']
    
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    
    alternatives = variant['ALT'].split(',')
    # When checking what features that are overlapped we use the longest alternative
    longest_alt = max([len(alternative) for alternative in alternatives])
    variant_position = int(variant['POS'])
    variant_interval = [variant_position, (variant_position + longest_alt-1)]
    
    #If annotated with vep we do not need to check interval trees
    if vep:
        variant['annotation'] = check_vep_annotation(variant)
        if len(variant['annotation']) > 0:
            variant['comp_candidate'] = True
    else:
        try:
            variant['annotation'] = set(gene_trees[chrom].find_range(variant_interval))
        except KeyError:
            if verbosity:
                warning(''.join(['Chromosome ', chrom, ' is not in annotation file!']))
        
        if whole_genes:
        # If compounds are to be checked in whole genes (including introns):
            if len(variant['annotation']) > 0:
                variant['comp_candidate'] = True
        else:
        #Check if exonic:
            try:
                if len(exon_trees[chrom].find_range(variant_interval)):
                    variant['comp_candidate'] = True
            except KeyError:
                if verbosity:
                    warning(''.join(['Chromosome ', chrom, ' is not in annotation file!']))
        return

def main():
    pass

if __name__ == '__main__':
    main()
