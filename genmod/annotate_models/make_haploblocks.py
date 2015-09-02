from __future__ import (print_function)
import logging
from interval_tree import IntervalTree

def get_haploblocks(variant_batch, individuals):
    """
    Take a variant batch and return the haploblocks for each of the individuals.
    The haploblocks are dictionaries with individual trees as keys.
    
    If a variant is phased it is denoted in the genotype call with a pipe 
    instead of a backslash.
    Unphased call: '0/1'
    Phased call: '0|1'
    A collection of censequtive phased variants makes a haploblock.
    The haploblocks are broken when a unphased call is seen.
    
    Arguments:
        variant_batch (dict): variant_batch is a dictionary with variant_id:s 
                              as keys and variant dictionaries as values.
        individuals (list): A list with strings that represents the
                            individual ids
    """
    logger = logging.getLogger(__name__)
    logger.debug("Init haploblocks")
    haploblocks = {ind_id:[] for ind_id in individuals}
    logger.debug("Haploblocks: {0}".format(haploblocks))
    haploblock_starts = {ind_id:None for ind_id in individuals}
    logger.debug("Set beginning to True")
    # This variable indicates if we are in a haploblock
    in_haploblock = {ind_id:False for ind_id in individuals}
    haploblock_id = 1
    interval_trees = {}
    for variant_id in variant_batch:
        logger.debug("Variant: {0}".format(variant_id))
        
        variant = variant_batch[variant_id]
        
        for ind_id in individuals:
            raw_gt_call = variant.get(ind_id, './.')
            # If the variant is phased we must check if it is the start of
            # a haploblock or in the middle of one
            if '|' in raw_gt_call:
                # Check if we are already in a haploblock
                logger.debug("Variant {0} is phased for individual {1}.".format(
                    variant_id, ind_id))
                
                if not in_haploblock[ind_id]:
                    logger.debug("Setting haploblock start to: {0}"\
                    " for individual {1}".format(
                        variant['POS'], ind_id))
                    
                    haploblock_starts[ind_id] = int(variant['POS'])
                    logger.debug("Setting is_haploblock to True for individual {0}".format(
                        ind_id
                    ))
                    in_haploblock[ind_id] = True
            
            # Here the variant is not phased
            else:
            #If call is not passed we consider it to be on same 
            # haploblock(GATK recommendations)
                if variant.get('FILTER', '.') == 'PASS':
                    # The intervals is a list wiht start, stop and id
                    if in_haploblock[ind_id]:
                        logger.debug("Creating a haploblock for individual {0}"\
                        " with start:{1}, stop:{2} and id:{3}".format(
                            ind_id,
                            haploblock_starts[ind_id],
                            int(variant['POS']) - 1,
                            str(haploblock_id)
                            ))
                        haploblocks[ind_id].append(
                            [
                                haploblock_starts[ind_id], 
                                int(variant['POS']) - 1,
                                str(haploblock_id)
                            ]
                        )
                        haploblock_id += 1
                        logger.debug("Setting haploblock id to {0}".format(
                            haploblock_id
                        ))
                        in_haploblock[ind_id] = False
                        logger.debug("Setting is_haploblock to False for individual {0}".format(
                            ind_id
                        ))
    
    for ind_id in individuals:
        #Check if we have just finished an interval
        
        if in_haploblock[ind_id]:
            logger.debug("Creating a haploblock for individual {0}"\
            " with start:{1}, stop:{2} and id:{3}".format(
                ind_id,
                haploblock_starts[ind_id],
                int(variant['POS']) - 1,
                str(haploblock_id)
                ))
            haploblocks[ind_id].append(
                                [
                                    haploblock_starts[ind_id], 
                                    int(variant['POS']), 
                                    str(haploblock_id)
                                ]
                            )
            
            haploblock_id += 1
            logger.debug("Setting haploblock id to {0}".format(
                haploblock_id
            ))
            
        # Create interval trees of the haploblocks
        if haploblocks[ind_id]:
            logger.debug("Creating IntervalTree for individual {0} with haploblocks:{1}, "\
            "start:{2}, stop:{3}".format(
                ind_id, haploblocks[ind_id], 
                haploblocks[ind_id][0][0]-1, haploblocks[ind_id][-1][1]+1
                ))
            interval_trees[ind_id] = IntervalTree(
                                            haploblocks[ind_id], 
                                            haploblocks[ind_id][0][0]-1, 
                                            haploblocks[ind_id][-1][1]+1
                                        )
        logger.debug("Interval tree created")
    
    return interval_trees
