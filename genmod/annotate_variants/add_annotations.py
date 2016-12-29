import logging

from genmod.vcf_tools import (add_metadata)

logger = logging.getLogger(__name__)

def add_regions(header):
    """Add region annotations to header"""
    logger.info("Adding 'Annotation' to vcf header")
    add_metadata(
        header,
        'info',
        'Annotation',
        annotation_number='.',
        entry_type='String',
        description='Annotates what feature(s) this variant belongs to.'
    )
    return

def add_exac(header):
    """Add exac annotation to vcf header"""
    logger.info("Adding 'EXACAF' to vcf header")    
    add_metadata(
        header,
        'info',
        'EXACAF',
        annotation_number='1',
        entry_type='Float',
        description="Frequency in the ExAC database."
    )
    return

def add_exac_max(header):
    """Add exac annotation to vcf header"""
    logger.info("Adding 'EXAC_MAX_AF' to vcf header")    
    add_metadata(
        header,
        'info',
        'EXAC_MAX_AF',
        annotation_number='1',
        entry_type='Float',
        description="The max af for ExAC populations."
    )
    return

def add_thousandg(header):
    """Add thousand genomes annotation to vcf header"""
    logger.info("Adding '1000GAF' to vcf header")    
    add_metadata(
        header,
        'info',
        '1000GAF',
        annotation_number='1',
        entry_type='Float',
        description="Frequency in the 1000G database."
    )
    return

def add_thousandg_max(header):
    """Add thousand genomes max annotation to vcf header"""
    logger.info("Adding '1000G_MAX_AF' to vcf header")    
    add_metadata(
        header,
        'info',
        '1000G_MAX_AF',
        annotation_number='1',
        entry_type='Float',
        description="The max af for thousand genomes populations."
    )
    return

def add_spidex(header):
    """Add spidex annotation to vcf header"""
    logger.info("Adding 'SPIDEX' to vcf header")    
    add_metadata(
        header,
        'info',
        'SPIDEX',
        annotation_number='1',
        entry_type='Float',
        description="Z score from the spidex database."
    )
    return

def add_cadd(header):
    """Add cadd annotation to vcf header"""
    logger.info("Adding 'CADD' to vcf header")    
    add_metadata(
        header,
        'info',
        'CADD',
        annotation_number='1',
        entry_type='Integer',
        description="The CADD relative score for this alternative."
    )
    return

def add_cadd_raw(header):
    """Add cadd annotation to vcf header"""
    logger.info("Adding 'CADD' to vcf header")    
    add_metadata(
        header,
        'info',
        'CADD_raw',
        annotation_number='1',
        entry_type='Float',
        description="The CADD raw score(s) for this alternative(s)."
    )
    return

def add_cosmic(header):
    """Add cosmic annotation to vcf header"""
    logger.info("Adding 'COSMIC' to vcf header")    
    add_metadata(
        header,
        'info',
        'COSMIC',
        annotation_number='0',
        entry_type='Flag',
        description="If variant is in COSMIC database."
    )
    return