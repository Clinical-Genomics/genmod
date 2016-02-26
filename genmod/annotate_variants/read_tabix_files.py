import logging

from tabix import TabixError

logger = logging.getLogger(__name__)


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
        logger.warning("Chromosome {0} does not seem to exist in {1}".format(
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
        logger.debug("Found record: {0}".format('\t'.join(record)))
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
    
    logger.debug("Found spidex score: {0}".format(spidex_score))
    
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
    cadd_scores = {
            'cadd_raw': None,
            'cadd_phred': None,
        }
    records = get_tabix_records(tabix_reader, chrom, start)    
    # CADD values are only for snps:
    for record in records:
        if record[3] == alt:
            #We need to send both cadd values
            cadd_scores['cadd_raw'] = record[-2]
            cadd_scores['cadd_phred'] = record[-1]
    
    return cadd_scores

def annotate_cadd_score(variant, cadd_raw, cadd_file, cadd_1000g, cadd_exac, 
                   cadd_ESP, cadd_InDels):
    """
    Add CADD scores to this variant.
    
    Take a variant and different file handles to tabix files with CADD info and
    tries to find the cadd scores for the variant.
    Default is to annotate the relative CADD score but if the flag cadd_raw is 
    used the raw cadd scores will be annotated to.
    
    Arguments:
        variant (dict): A variant dictionary
        cadd_raw (bool): If the raw CADD scores should be annotated
        cadd_file (tabix file handle)
        cadd_1000g (tabix file handle)
        cadd_exac (tabix file handle)
        cadd_ESP (tabix file handle)
        cadd_InDels (tabix file handle)
    
    """
    cadd_score = None
    cadd_relative = None
    cadd_absolute = None
    cadd_relative_scores = []
    cadd_absolute_scores = []
    #Check CADD file(s):
    for alt in variant['ALT'].split(','):
        if cadd_file:
            cadd_scores = get_cadd_scores(
                                    cadd_file, 
                                    variant['CHROM'], 
                                    variant['POS'], 
                                    alt
                                )
            cadd_relative = cadd_scores['cadd_phred']
            cadd_absolute = cadd_scores['cadd_raw']
        # If variant not found in big CADD file check the 1000G file:
        if not (cadd_relative and cadd_absolute) and cadd_1000g:
            cadd_scores = get_cadd_scores(
                                    cadd_1000g, 
                                    variant['CHROM'], 
                                    variant['POS'], 
                                    alt
                                )
            cadd_relative = cadd_scores['cadd_phred']
            cadd_absolute = cadd_scores['cadd_raw']
        
        if not (cadd_relative and cadd_absolute) and cadd_exac:
            cadd_scores = get_cadd_scores(
                                    cadd_exac, 
                                    variant['CHROM'], 
                                    variant['POS'], 
                                    alt
                                )
            cadd_relative = cadd_scores['cadd_phred']
            cadd_absolute = cadd_scores['cadd_raw']
        
        if not (cadd_relative and cadd_absolute) and cadd_ESP:
            cadd_scores = get_cadd_scores(
                                    cadd_ESP, 
                                    variant['CHROM'], 
                                    variant['POS'], 
                                    alt
                                )
            cadd_relative = cadd_scores['cadd_phred']
            cadd_absolute = cadd_scores['cadd_raw']
        
        if not (cadd_relative and cadd_absolute) and cadd_InDels:
            cadd_scores = get_cadd_scores(
                                    cadd_ESP, 
                                    variant['CHROM'], 
                                    variant['POS'], 
                                    alt
                                )
            cadd_relative = cadd_scores['cadd_phred']
            cadd_absolute = cadd_scores['cadd_raw']
        
        # If there are several alternatives we want to annotate with multiple
        # values
        if cadd_relative:
            cadd_relative_scores.append(str(cadd_relative))
        else:
            cadd_relative_scores.append('-')
        if cadd_absolute:
            cadd_absolute_scores.append(str(cadd_absolute))
        else:
            cadd_absolute_scores.append('-')
    
    any_annotation = False
    for annotation in cadd_relative_scores:
        if annotation != '-':
            any_annotation = True
    
    if any_annotation:
        variant['CADD'] = ','.join(cadd_relative_scores)
        if cadd_raw:
            variant['CADD_raw'] = ','.join(cadd_absolute_scores)
    
    return variant

