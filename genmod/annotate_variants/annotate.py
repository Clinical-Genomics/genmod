import logging

from genmod.annotate_regions.get_features import get_region
from genmod.annotate_variants.read_tabix_files import (get_frequencies, 
     get_spidex_score, get_cadd_scores)

logger = logging.getLogger(__name__)

def annotate_variant(variant, annotation_arguments):
    """Annotate a variant based on what arguments that are passed"""
    variant_info = variant.rstrip().split('\t')
    chrom = variant_info[0]
    if chrom.startswith(('chr', 'CHR', 'Chr')):
        chrom = chrom[3:]
    pos = int(variant_info[1])
    ref = variant_info[3]
    alt = variant_info[4]
    
    info = variant_info[7]
    if info == '.':
        info = []
    else:
        info = info.split(';')
    
    ## TODO this needs to be handeled different for SV:s
    start = pos
    end = pos + max(len(ref), len(alt))
    
    #Check which annotations that are available
    regions = None
    if 'region_trees' in annotation_arguments:
        regions = get_region(chrom, start, end, annotation_arguments['region_trees'])
        if regions:
            info.append("Annotation={0}".format(','.join(regions)))
    
    if 'exac' in annotation_arguments:
        reader = annotation_arguments['exac']
        frequencies = get_frequencies(reader, chrom, start, alt)
        if 'AF' in frequencies:
            info.append("EXACAF={0}".format(frequencies['AF']))
        if annotation_arguments.get('max_af'):
            if 'MAX_AF' in frequencies:
                info.append("EXAC_MAX_AF={0}".format(frequencies['MAX_AF']))

    if 'thousand_g' in annotation_arguments:
        reader = annotation_arguments['thousand_g']
        frequencies = get_frequencies(reader, chrom, start, alt)
        if 'AF' in frequencies:
            info.append("1000GAF={0}".format(frequencies['AF']))
        if annotation_arguments.get('max_af'):
            if 'MAX_AF' in frequencies:
                info.append("1000G_MAX_AF={0}".format(frequencies['MAX_AF']))

    if 'spidex' in annotation_arguments:
        reader = annotation_arguments['spidex']
        spidex_score = get_spidex_score(reader, chrom, start, alt)
        if spidex_score:
            info.append("SPIDEX={0}".format(spidex_score))

    if 'cadd_files' in annotation_arguments:
        readers = annotation_arguments['cadd_files']
        cadd_scores = {}
        for reader in readers:
            if not cadd_scores:
                cadd_scores = get_cadd_scores(reader, chrom, start, alt)

        if cadd_scores:
            info.append("CADD={0}".format(cadd_scores['cadd_phred']))
            if annotation_arguments.get('cadd_raw'):
                info.append("CADD_raw={0}".format(cadd_scores['cadd_raw']))

    # Rebuild the info string
    if len(info) > 0:
        info_string = ';'.join(info)
    else:
        info_string = '.'
    
    variant_info[7] = info_string
    
    return '\t'.join(variant_info)