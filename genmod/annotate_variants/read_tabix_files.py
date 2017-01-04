import os
import logging

import tabix

from tabix import TabixError

logger = logging.getLogger(__name__)

def get_tabixhandle(path):
    """Check if a file is zipped and that the index exists
        If something looks wierd raise a TabixError
    """
    if not path.endswith('.gz'):
        raise TabixError("File {0} does not end with '.gz'".format(path))
    index_file = path + '.tbi'
    if not os.path.isfile(index_file):
        raise TabixError("No index could be found for {0}".format(path))
    
    return tabix.open(path)

def get_tabix_records(tabix_reader, chrom, start):
    """Get the tabix records for some given coordinates
    
    Args:
        tabix_reader (Tabix.reader): A Tabix object
        chrom (str): The chromosome of the position
        start (str): The start position of the variant
        alt (str): The alternative sequence
    
    Returns:
        records (Iterable): The overlapping records found
    """
    records = []
    logger.debug("Looking for record in %s" % (tabix_reader))
    logger.debug("Looking for records with chr:%s, pos:%s" % (chrom, start))
    tabix_key = int(start)
    try:
        records = tabix_reader.query(chrom, tabix_key-1, tabix_key)
    except TypeError:            
        records = tabix_reader.query(str(chrom), tabix_key-1, tabix_key)
    except TabixError:
        logger.info("Chromosome {0} does not seem to exist in {1}".format(
            chrom, tabix_reader))
    except:
        pass
    
    return records

def get_frequencies(tabix_reader, chrom, start, alt):
    """
    Return the frequencies from a tabix indexed vcf file.
    
    Arguments:
        tabix_reader (Tabix.reader): A Tabix object
        chrom (str): The chromosome of the position
        start (str): The start position of the variant
        alt (str): The alternative sequence
    
    Returns:
        frequencies (dict): A dictionary with relevant frequencies
    """
    freq = None
    records = get_tabix_records(tabix_reader, chrom, start)

    frequencies = {}
    for record in records:
        logger.debug("Found record: %s" % '\t'.join(record))
        #We can get multiple rows so need to check each one
        #We also need to check each one of the alternatives per row
        for i,alternative in enumerate(record[4].split(',')):
            if alternative == alt:
                for info in record[7].split(';'):
                    info = info.split('=')
                    if info[0] in ['AF', 'MAX_AF']:
                        freqs = info[-1].split(',')
                        frequencies[info[0]] = freqs[i]

    return frequencies

def get_spidex_score(tabix_reader, chrom, start, alt):
    """
    Return the record from a spidex file.
    
    Arguments:
        tabix_reader (Tabix.reader): A Tabix object
        chrom (str): The chromosome of the position
        start (str): The start position of the variant
        alt (str): The alternative sequence
    
    Returns:
        spidex_score float: The spidex z scores for this position
    
    """
    records = get_tabix_records(tabix_reader, chrom, start)
    spidex_score = None
    
    for record in records:
        if record[3] == alt:
            #We need to send both cadd values
            spidex_score = float(record[5])
    
    logger.debug("Found spidex score: %s" % str(spidex_score))

    return spidex_score

def get_cosmic(tabix_reader, chrom, start, alt):
    """
    Return if record exists in cosmic database.
    
    Arguments:
        tabix_reader (Tabix.reader): A Tabix object
        chrom (str): The chromosome of the position
        start (str): The start position of the variant
        alt (str): The alternative sequence
    
    Returns:
        in_cosmic (bool): If variant is in COSMIC
    
    """
    records = get_tabix_records(tabix_reader, chrom, start)
    in_cosmic = False
    
    for record in records:
        if record[4] == alt:
            #We need to send both cadd values
            in_cosmic = True
    if in_cosmic:
        logger.debug("Variant was found in COSMIC")
    else:
        logger.debug("Variant was not found in COSMIC")

    return in_cosmic

def get_cadd_scores(tabix_reader, chrom, start, alt):
    """
    Return the record from a cadd file.
    
    Arguments:
        tabix_reader (Tabix.reader): A Tabix object
        chrom (str): The chromosome of the position
        start (str): The start position of the variant
        alternatives (str): The alternative sequence
    
    Returns:
        cadd_scores (dict): The cadd scores for this position
    
    """
    cadd_scores = {}
    records = get_tabix_records(tabix_reader, chrom, start)    
    # CADD values are only for snps:
    for record in records:
        if record[3] == alt:
            #We need to send both cadd values
            cadd_scores['cadd_raw'] = record[-2]
            cadd_scores['cadd_phred'] = record[-1]
    
    return cadd_scores

