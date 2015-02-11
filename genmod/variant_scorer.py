#!/usr/bin/env python
# encoding: utf-8
"""
variant_scorer.py

Score a file with variant info in vcf format.

Creates batches and put them into a queue object.
The batches are dictionary objects with overlapping features where the feature
id:s are keys and the values are dictionarys with variants.

Batch = {feature_1_id:{variant_1_id:variant_1_info, variant_2_id:
variant_2_info}, feature_2_id:... }

Created by Måns Magnusson on 2014-03-17.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function, unicode_literals


import sys
import os
import argparse
from datetime import datetime
from tempfile import NamedTemporaryFile
from multiprocessing import Process

from pprint import pprint as pp

from genmod.models import score_variants

# Import third party library
# https://github.com/mitsuhiko/logbook
from logbook import Logger, StderrHandler
log = Logger('Logbook')
log_handler = StderrHandler()


class VariantScorer(Process):
    """Creates parser objects for scoring variants in vcf files"""
    def __init__(self, variant_queue, results, header,
                models_of_inheritance, family_id, alt_dict, score_dict, 
                value_dict, operation_dict, verbose):
        super(VariantScorer, self).__init__()
        self.variant_queue = variant_queue
        self.results = results
        self.header = header
        self.prefered_models = models_of_inheritance
        self.family_id = family_id
        self.alt_dict = alt_dict
        self.score_dict = score_dict
        self.value_dict = value_dict
        self.operation_dict = operation_dict
        self.verbose = verbose
        # self.verbosity = args.verbose
        # self.phased = args.phased
    
    def score_compounds(self, batch, family_id):
        """
        Score the compounds in a batch.
        Takes av input a batch with variants and scores all the compound pairs.
        If any of the compounds have a score<10 then the total score will be
        penalized with -6.
        
        Arguments:
            batch  : A dictionary with variants
        
        Returns:
            Nothing
            Only updates the variants in the bath.
        """
        for variant_id in batch:
            comp_list = []
            variant = batch[variant_id]
            only_compound = True
            for model in variant['genetic_models'].get(family_id,[]):
                if model not in ['AR_comp', 'AR_comp_dn']:
                    only_compound = False
            
            rank_score = int(variant.get('Individual_rank_score', 0))
            # If the variant is only AR_compound and one of the variants in the 
            # pair have a low compound score we want to downprioritize that
            # variant
            
            highest_scoring_partner = 0
            for compound in variant.get('compound_variants', {}).get(family_id, []):
                compound_id = compound['variant_id']
                
                pair_2_score = int(batch.get(compound_id, {}).get(
                                                'Individual_rank_score', 0))
                
                compound['compound_score'] = rank_score + pair_2_score
                
                if pair_2_score > highest_scoring_partner:
                    highest_scoring_partner = pair_2_score
            
            # If the variant is only 'AR_comp' and all compound partners are
            # low scoring variants we want to downprioritize the variant.
            if rank_score > 10:
                if only_compound:
                    if (highest_scoring_partner < 10):
                        rank_score -= 6
                        variant['Individual_rank_score'] = rank_score
                        for compound in variant.get(
                                        'compound_variants', {}).get(family_id, []):
                            compound['compound_score'] -= 6
    
    
    def make_scored_compound_string(self, variant, family_id):
        """
        Make a new compound string with scores in the proper format.
        
        Args:
            variant     : A dictionary with variant information
            family_id   : The id of the family that we want to change for.
        
        Returns:
            compound_string : A string on the proper compound format
        
        """
        compound_list = []
        for compound in variant.get('compound_variants', {}).get(family_id, []):
            compound_id = compound['variant_id']
            compound_score = str(compound['compound_score'])
            compound_list.append(compound_id + '>' + compound_score)
        
        return '|'.join(compound_list)
        
    def make_print_version(self, batch, family_id):
        """
        Add the new informatin in the proper vcf format.
        
        Args:
            batch  : A dictionary with variants
            family_id   : The id of the family that we want to change for.
        
        """
        for variant_id in batch:
            # Modify the INFO field:
            variant = batch[variant_id]
            
            # We turn the info field into a list
            info_field = variant['INFO'].split(';')
            # We need to replace the compound field with the annotated variants
            for position in range(len(info_field)):
                
                entry_info = info_field[position].split('=')
                if entry_info[0] == 'Compounds':
                    splitted_compunds = entry_info[1].split(',')
                    for family_number in range(len(splitted_compunds)):
                        if splitted_compunds[family_number].split(':')[0] == family_id:
                            scored_compounds = self.make_scored_compound_string(
                                    variant,
                                    family_id
                                )
                            splitted_compunds[family_number] = ':'.join([family_id, scored_compounds])
                    
                    info_field[position] = 'Compounds=' + ','.join(splitted_compunds)
            
            info_field.append('IndividualRankScore=' + family_id + ':' +
                              str(variant['original_rank_score']))
            
            # Add the rank score to the info field 
            info_field.append('RankScore=' + family_id + ':' +
                              str(variant['Individual_rank_score']))
            
            variant['INFO'] = ';'.join(info_field)
        
        return

    def run(self):
        """
        Score all variants in the batches.
        
        Take the batches from the queue and score all variants.
        Put the scored variants in the results queue for further processing.
        
        """
        proc_name = self.name
        
        if self.verbose:
            log.info('%s: Starting!' % proc_name)
        
        while True:
            # A batch is a dictionary on the form {variant_id:variant_dict}
            variant_batch = self.variant_queue.get()
            
            if variant_batch is None:
                self.variant_queue.task_done()
                if self.verbose:
                    log.info('%s: Exiting' % proc_name)
                break
            
            # We can now free som space by removing the haploblocks
            variant_batch.pop('haploblocks', None)
            
            score_variants(
                        variant_batch, 
                        self.prefered_models,
                        self.alt_dict, 
                        self.score_dict,
                        self.value_dict,
                        self.operation_dict,
                        self.verbose
                    )
            
            # We want to store the original rank score for traceability since
            # the rank score can change during the process:
            for variant_id in variant_batch:
                variant = variant_batch[variant_id]
                variant['original_rank_score'] = variant.get('Individual_rank_score', '0')
            
            self.score_compounds(variant_batch, self.family_id)
            self.make_print_version(variant_batch, self.family_id)
            
            self.results.put(variant_batch)
            self.variant_queue.task_done()


def main():
    from tempfile import NamedTemporaryFile
    from vcf_parser import parser
    from codecs import open

    parser = argparse.ArgumentParser(description=
                                     "Parse different kind of pedigree files.")
    parser.add_argument('variant_file',
                        type=str, nargs=1,
                        help='A file with variant information.'
                        )
    parser.add_argument('outfile',
                        type=str, nargs=1,
                        help='Specify the path to output.'
                        )
    parser.add_argument('-phased', '--phased',
                        action="store_true",
                        help='If variant file is phased.'
                        )
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        help='Increase output verbosity.'
                        )
    parser.add_argument('-vep', '--vep',
                        action="store_true",
                        help='If variants are annotated with vep.'
                        )

    args = parser.parse_args()
    var_file = args.variant_file[0]

    my_vcf_parser = vcf_parser.VCFParser(var_file)
    temporary_variant_file = NamedTemporaryFile(mode='w')

    my_parser = VariantScorer(my_vcf_parser, temporary_variant_file)
    my_parser.parse()
    temporary_variant_file.seek(0)

    with open(temporary_variant_file.name, 'r', encoding='utf-8') as f:
        for variant_line in f:
            print(variant_line.rstrip().split('\t')[7].split(';')[-1].split('=')[-1])

    # outFile=args.outfile[0]
    var_sorter = variant_sorter.FileSort(temporary_variant_file,
                                         outFile=args.outfile[0])
    var_sorter.sort()

    temporary_variant_file.close()

if __name__ == '__main__':
    main()
