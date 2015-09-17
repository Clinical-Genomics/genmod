import logging

# def check_tabix_index(compressed_file, file_type='cadd', verbose=False):
#     """
#     Check if a compressed file have a tabix index, if not build one.
#
#     Args:
#         compressed_file (str): Path to a file that is assumed to be compressed.
#         file_type (str): The type of the file. ('cadd' or 'vcf')
#         verbose (bool): Increase output verbosity
#
#     Returns:
#         0 if everything went ok.
#
#     """
#     if file_type == 'cadd':
#         try:
#             tabix_index(compressed_file, seq_col=0, start_col=1, end_col=1, meta_char='#')
#         except IOError as e:
#             pass
#     elif file_type == 'vcf':
#         try:
#             tabix_index(compressed_file, preset='vcf')
#         except IOError as e:
#             pass
#     return 0


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
    tabix_key = int(start)
    try:
        for record in tabix_reader.query(chrom, tabix_key-1, tabix_key):
            #We can get multiple rows so need to check each one
            #We also need to check each one of the alternatives per row
            for i,alternative in enumerate(record[4].split(',')):
                if alternative == alt:
                    for info in record[7].split(';'):
                        info = info.split('=')
                        if info[0] == 'AF':
                            frequencies = info[-1].split(',')
                            return frequencies[i]
    except TypeError:            
        for record in tabix_reader.query(str(chrom), tabix_key-1, tabix_key):
            #We can get multiple rows so need to check each one
            #We also need to check each one of the alternatives per row
            for i, alternative in enumerate(record[4].split(',')):
                if alternative == alt:
                    for info in record[7].split(';'):
                        info = info.split('=')
                        if info[0] == 'AF':
                            frequencies = info[-1].split(',')
                            return frequencies[i]
    except:
        pass
    
    return freq


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

