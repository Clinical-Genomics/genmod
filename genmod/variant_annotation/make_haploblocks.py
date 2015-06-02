
from interval_tree import IntervalTree

def get_haploblocks(variant_batch, individuals):
    """
    Take a variant batch and return the haploblocks for each of the individuals.
    The haploblocks are dictionaries with individual trees as keys.
    
    Arguments:
        variant_batch (dict): variant_batch is a dictionary with variant_id:s 
                              as keys and variant dictionaries as values.
        individuals (list): A list with strings that represents the
                            individual ids
    """
    haploblocks = {ind_id:[] for ind_id in individuals}
    haploblock_id = 1
    interval_trees = {}
    
    beginning = True
    # We collect the starts of the haploblocks
    for variant_id in variant_batch:
        variant = variant_batch[variant_id]
        if beginning:
            # We need the start positions for the first haplotype block
            haploblock_starts = {ind_id:int(variant['POS']) 
                                    for ind_id in individuals}
            beginning = False
        
        else:
            for ind_id in individuals:
                raw_gt_call = variant.get(ind_id, './.')
                #A new haploblock is indicated by '/' if the data is phased
                if '/' in raw_gt_call:
                #If call is not passed we consider it to be on same 
                # haploblock(GATK recommendations)
                    if variant.get('FILTER', '.') == 'PASS':
                        # The intervals is a list wiht start, stop and id
                        haploblocks[ind_id].append(
                                            [   
                                                haploblock_starts[ind_id], 
                                                int(variant['POS']) - 1,
                                                str(haploblock_id)
                                            ]
                                        )
                        haploblock_id += 1
                        haploblock_starts[ind_id] = int(variant['POS'])
    
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
        interval_trees[ind_id] = IntervalTree(
                                                haploblocks[ind_id], 
                                                haploblocks[ind_id][0][0]-1, 
                                                haploblocks[ind_id][-1][1]+1
                                        )
    return interval_trees
    
    