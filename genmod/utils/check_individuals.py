import logging

def check_individuals(family_parser, vcf_individuals):
    """
    Check if the individuals from ped file is in vcf file
    
    Arguments:
        family_parser (FamilyParser): A FamilyParser object
        vcf_individuals (list): the individuals found in the vcf file 
    """
    logger = logging.getLogger(__name__)
    logger.info("Checking if individuals in ped file exists in vcf file")
    for individual in family_parser.individuals:
        if individual not in vcf_individuals:
            raise IOError("Individuals in PED file must exist in VCF file") # Raise proper exception here
    return