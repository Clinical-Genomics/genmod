from __future__ import (print_function, absolute_import)

import sys
import os
import logging

from datetime import datetime
from collections import OrderedDict

import click

from . import get_annotation

def get_batches(variants, batch_queue, vep=False):
    """
    Create variant batches based on their annotation and put them into the 
    batch queue.
    
    Variants are given a new 'annotation' field in the variant dictionary and 
    also a 'exonic' field.
    get_batches will then use the annotation to search for sequences of variants
    with overlapping annotations. These are collected into one batch and gets
    put into a queue.
    Variants that are in between features will be in their own batch.
    
    Arguments:
         variants (Iterator):An iterator that returns variant dictionaries
         batch_queue (Queue): A queue where the batches will be putted
         vep (bool): If variant is annotated with vep
    
    Returns:
         Does not return but put the results in a queue
    """
    logger = logging.getLogger(__name__)
    #For testing only:
    logger = logging.getLogger("genmod.utils.get_batches")
    
    logger.debug("Set beginning to True")
    beginning = True
    logger.debug("Create first empty batch")
    # A batch is a ordered dictionary with variants
    batch = OrderedDict()
    new_chrom = None
    current_chrom = None
    current_features = []
    chromosomes = []
    
    start_parsing_time = datetime.now()
    start_chrom_time = start_parsing_time
    start_twenty_time = start_parsing_time
    
    nr_of_variants = 0
    nr_of_batches = 0
    
    logger.info("Start parsing the variants")
    for variant in variants:
        variant_id = variant['variant_id']
        logger.debug("Checking variant {0}".format(variant_id))
        
        nr_of_variants += 1
        new_chrom = variant['CHROM']
        if new_chrom.startswith('chr'):
            new_chrom = new_chrom[3:]

        logger.debug("Update new chrom to {0}".format(new_chrom))

        new_features = get_annotation(variant, vep)
        logger.debug("Adding {0} to variant {1}".format(
            ', '.join(new_features), variant_id
        ))
        variant['annotation'] = new_features
        
        
        if nr_of_variants % 20000 == 0:
            logger.info("{0} variants parsed".format(nr_of_variants))
            logger.info("Last 20.000 took {0} to parse.".format(
                str(datetime.now() - start_twenty_time)))
            start_twenty_time = datetime.now()
        
        if beginning:
            logger.debug("First variant.")
            current_features = new_features
            
            logger.debug("Adding {0} to variant batch".format(variant_id))
            batch[variant_id] = variant
            logger.debug("Updating current chrom to {0}".format(new_chrom))
            current_chrom = new_chrom
            
            chromosomes.append(current_chrom)
            logger.debug("Adding chr {0} to chromosomes".format(new_chrom)) 
            
            beginning = False
            logger.debug("Updating beginning to False")
        
        else:
            # If we should put the batch in the queue:
            logger.debug("Updating send to True") 
            send = True
            
            # Check if the variant ovelapps any features
            if len(new_features) != 0:
                # Check if the features overlap the previous variants features
                if new_features.intersection(current_features):
                    logger.debug("Set send to False since variant features overlap") 
                    send = False
            
            # If we are at a new chromosome we finish the current batch:
            if new_chrom != current_chrom:
                if current_chrom not in chromosomes:
                    chromosomes.append(current_chrom)
                logger.debug("Adding chr {0} to chromosomes".format(new_chrom)) 
                # New chromosome means new batch
                send = True
                logger.info("Chromosome {0} parsed. Time to parse"\
                            " chromosome: {0}".format(
                            current_chrom, datetime.now()-start_chrom_time))
                start_chrom_time = datetime.now()
                current_chrom = new_chrom
            
            if send:
                # Put the job in the queue
                logger.debug("Adding batch in queue") 
                batch_queue.put(batch)
                nr_of_batches += 1
                #Reset the variables
                current_features = new_features
                logger.debug("Initializing empty batch") 
                batch = {}
            else:
                current_features = current_features.union(new_features)
            
            logger.debug("Adding variant {0} to batch".format(variant_id)) 
            batch[variant_id] = variant
    
    if current_chrom not in chromosomes:
        logger.debug("Adding chr {0} to chromosomes".format(current_chrom))
        chromosomes.append(current_chrom)
    nr_of_batches += 1
    
    logger.info("Chromosome {0} parsed. Time to parse"\
                " chromosome: {0}".format(
                current_chrom, datetime.now()-start_chrom_time))
    
    logger.info("Variants parsed. Time to parse variants: {0}".format(
        str(datetime.now() - start_parsing_time)
    ))
    
    logger.info("Number of variants in variant file: {0}".format(nr_of_variants))
    logger.info("Number of batches created:{0}".format(nr_of_batches))
    
    logger.debug("Adding batch in queue") 
    
    batch_queue.put(batch)
    
    return chromosomes
        
@click.command()
@click.argument('variant_file', 
                    nargs=1, 
                    type=click.Path(),
                    metavar='<vcf_file> or -'
)
def cli(variant_file):
    """docstring for cli"""
    from Queue import Queue
    from vcf_parser import VCFParser
    
    if variant_file == '-':
        variant_parser = VCFParser(fsock = sys.stdin, split_variants=True)
    else:
        variant_parser = VCFParser(infile = variant_file, split_variants=True)
    
    batch_queue = Queue()
    chromosomes = get_batches(variant_parser, batch_queue)
    print(chromosomes)

if __name__ == '__main__':
    from genmod import logger
    from genmod.log import init_log
    init_log(logger, loglevel="DEBUG")
    cli()






