from genmod.utils import get_rank_score
from genmod.vcf_tools import get_variant_dict, get_info_dict
from genmod.score_variants.rank_score_variant_definitions import RANK_SCORE_TYPE_NAMES

def get_variant_line():
    """Return a vcf formatted variant line."""
    pass

def test_get_rank_score():
    """docstring for test_get_rank_score"""
    variant_line = "1\t879537\t.\tT\tC\t100\tPASS\tMQ=1;GeneticModels=1:AR_hom;"\
    "ModelScore=1:55;Annotation=SAMD11;CADD=1.248;Exonic;RankScore=1:23;RankScoreNormalized=1:0.2\t"\
    "GT:AD:GQ\t0/1:10,10:60"
    
    assert float(get_rank_score(variant_line = variant_line)) == float('23')

def test_get_rank_score_no_score():
    """docstring for test_get_rank_score"""
    variant_line = "1\t879537\t.\tT\tC\t100\tPASS\tMQ=1;GeneticModels=1:AR_hom;"\
    "ModelScore=1:55;Annotation=SAMD11;CADD=1.248;Exonic\t"\
    "GT:AD:GQ\t0/1:10,10:60"

    for rank_score_type in RANK_SCORE_TYPE_NAMES:
        assert float(get_rank_score(variant_line=variant_line, rank_score_type=rank_score_type)) == float('-100')

def test_get_rank_score_multiple_families():
    """docstring for test_get_rank_score"""
    variant_line = "1\t879537\t.\tT\tC\t100\tPASS\tMQ=1;GeneticModels=1:AR_hom;"\
    "ModelScore=1:55;Annotation=SAMD11;CADD=1.248;Exonic;RankScore=1:23,2:12;RankScoreNormalized=1:0.2\t"\
    "GT:AD:GQ\t0/1:10,10:60"
    
    assert float(get_rank_score(variant_line = variant_line, family_id='2')) == float('12')

def test_get_rank_score_dict():
    """docstring for test_get_rank_score"""
    header_line = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t1"
    variant_info = "1\t879537\t.\tT\tC\t100\tPASS\tMQ=1;GeneticModels=1:AR_hom;"\
    "ModelScore=1:55;Annotation=SAMD11;CADD=1.248;Exonic;RankScore=1:23;RankScoreNormalized=1:0.2\t"\
    "GT:AD:GQ\t0/1:10,10:60"
    variant_dict = get_variant_dict(
        variant_line=variant_info, 
        header_line = header_line.split()
    )
    
    variant_dict['info_dict'] = get_info_dict(variant_dict['INFO'])
    
    assert float(get_rank_score(variant_dict = variant_dict)) == float('23')


def test_get_rank_score_normalized():
    """ Test fetching normalized rank score form VCF INFO field """
    # GIVEN an INFO field string with both raw and normalized rank score
    variant_line = "1\t879537\t.\tT\tC\t100\tPASS\tMQ=1;GeneticModels=1:AR_hom;" \
                   "ModelScore=1:55;Annotation=SAMD11;CADD=1.248;Exonic;RankScore=1:23;RankScoreNormalized=1:0.2\t" \
                   "GT:AD:GQ\t0/1:10,10:60"
    # WHEN fetching normalized rank score
    # THEN expect correct value fetched
    assert float(get_rank_score(variant_line=variant_line, rank_score_type='RankScoreNormalized')) == 0.2
