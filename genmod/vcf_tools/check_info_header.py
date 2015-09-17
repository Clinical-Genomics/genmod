import logging

def check_info(info_key, head):
    """Check if a info key exists in the vcf header
        
        Args:
            info_key (str): The name of the info key
            head (HeaderParser): A header parser object
    """
    logger = logging.getLogger(__name__)
    logger.debug("Checking if {0} is in vcf header.".format(info_key))
    if info_key in head.info_dict:
        logger.debug("{0} is specified in vcf header.".format(info_key))
        return True
    logger.warning("{0} is not specified in vcf header.".format(info_key))
    
    return False
