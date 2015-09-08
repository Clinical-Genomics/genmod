import logging

def check_individuals(ped_individuals, vcf_individuals):
    """
    Check if the individuals from ped file is in vcf file
    
    Arguments:
        ped_individuals (iterator): An iterator with strings
        vcf_individuals (iterator): An iterator with strings
    
    Returns:
        bool: if the individuals exists
    """
    logger = logging.getLogger(__name__)
    logger.info("Checking if individuals in ped file exists in vcf file")
    for individual in ped_individuals:
        if individual not in vcf_individuals:
            logger.error("Individuals in ped:{0}, individuals in vcf: {0}".format(
                ' '.join(ped_individuals), ' '.join(vcf_individuals)
            ))
            raise IOError("Individuals in PED file must exist in VCF file") # Raise proper exception here
    
    return True