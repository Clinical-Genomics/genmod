from genmod.vcf_tools import Genotype

def get_genotypes(variant, individuals):
    """Create genotype objects
        
        Create Genotype objects for all individuals and return them in a 
        dictionary
    
        Args:
            variant (dict): A variant dict
            individuals (list): A list with strings that are individual id:s
        
        Returns:
            genotype_dict (dict): A dictionary with individual id:s as strings
            and Genptype objects as keys
    """
    gt_format = variant.get('FORMAT', '').split(':')
    
    genotype_dict = {}
    
    for individual in individuals:
        gt_info = variant[individual].split(':')
        gt_call = dict(zip(gt_format, gt_info))
        
        #Create a genotype object for this individual
        genotype_dict[individual] = Genotype(**gt_call)
    
    return genotype_dict
    