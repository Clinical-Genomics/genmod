
import logging

def get_frequency(tabix_reader, chrom, start, alt):
    """
    Return the frequency from a tabix indexed vcf file.
    
    Arguments:
        tabix_reader (Tabix.reader): A Tabix object
        chrom (str): The chromosome of the position
        start (str): The start position of the variant
        alt (str): The alternative sequence
    
    Returns:
        freq (str): The frequency for this position
    """
    logger = logging.getLogger(__name__)
    freq = None
    # CADD values are only for snps:
    cadd_key = int(start)
    try:
        for record in tabix_reader.query(chrom, cadd_key-1, cadd_key):
            i = 0
            #We can get multiple rows so need to check each one
            #We also need to check each one of the alternatives per row
            for alternative in record[4].split(','):
                if alternative == alt:
                    for info in record[7].split(';'):
                        info = info.split('=')
                        if info[0] == 'AF':
                            frequencies = info[-1].split(',')
                            return frequencies[i]
                i += 1
    except TypeError:            
        for record in tabix_reader.query(str(chrom), cadd_key-1, cadd_key):
            i = 0
            #We can get multiple rows so need to check each one
            #We also need to check each one of the alternatives per row
            for alternative in record[4].split(','):
                if alternative == alt:
                    for info in record[7].split(';'):
                        info = info.split('=')
                        if info[0] == 'AF':
                            frequencies = info[-1].split(',')
                            return frequencies[i]
                i += 1
    except:
        pass
    
    return freq

def get_cadd_scores(tabix_reader, chrom, start, alt=None):
    """
    Return the record from a cadd file.
    
    Arguments:
        tabix_reader (Tabix.reader): A Tabix object
        chrom (str): The chromosome of the position
        start (str): The start position of the variant
        alt (str): The alternative sequence
    
    Returns:
        cadd_scores (dict): The cadd scores for this position
    
    """
    cadd_scores = {
            'cadd_raw': None,
            'cadd_phred': None,
        }
    # CADD values are only for snps:
    cadd_key = int(start)
    try:
        for record in tabix_reader.query(chrom, cadd_key-1, cadd_key):
            if record[3] == alt:
                #We need to send both cadd values
                cadd_scores['cadd_raw'] = record[-2]
                cadd_scores['cadd_phred'] = record[-1]
    
    except TypeError:
        for record in tabix_reader.query(str(chrom), cadd_key-1, cadd_key):
            record = record.split('\t')
            if record[3] == alt:
                #We need to send both cadd values
                cadd_scores['cadd_raw'] = record[-2]
                cadd_scores['cadd_phred'] = record[-1]
    except:
        pass
    
    return cadd_scores

