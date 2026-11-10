from genmod.plugin import Plugin
import pytest

vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'ADM1003A3']

def get_variant_line(chrom='2', pos='191964633', db_id='rs7574865', ref='T', 
alt='G', qual='1258.77', filt="PASS", info="", form="GT:AD:GQ:PL", 
gt_calls={'ADM1003A3': "1/1:0,35:99:1287,105,0"}):
    """
    Return a vcf formated variant line
    """
    if not info:
        info = "1000GAF=0.744609;AC=2;AF=1.00;AN=2;CLNSIG=255,10;"\
        "CSQ=G|ENSG00000138378|ENST00000358470|Transcript|intron_variant"\
        "||||||||-1|STAT4|HGNC|11365|protein_coding||CCDS2310.1|ENSP00000351255"\
        "|STAT4_HUMAN|Q53S87_HUMAN&Q53RU2_HUMAN&E9PBE2_HUMAN&C9JM11_HUMAN&C9JFG0_HUMAN"\
        "&B7SIX5_HUMAN|UPI00000015F2||||3/23||ENST00000358470.4:c.274-23582A>C|||||"\
        ",G|ENSG00000138378|ENST00000392320|Transcript|intron_variant|||||||"\
        "|-1|STAT4|HGNC|11365|protein_coding||CCDS2310.1|ENSP00000376134"\
        "|STAT4_HUMAN|Q53S87_HUMAN&Q53RU2_HUMAN&E9PBE2_HUMAN&C9JM11_HUMAN&"\
        "C9JFG0_HUMAN&B7SIX5_HUMAN|UPI00000015F2||||3/23||ENST00000392320.2:c.274-23582A>C|||||"\
        ";Clinical_db_gene_annotation=OMIM-141201,FullList;DB;DP=35;Dbsnp129LCAF=0.2493"

    variant_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(
        chrom, pos, db_id, ref, alt, qual, filt, info, form
    )
    for ind_id in gt_calls:
        variant_line = "{0}\t{1}".format(variant_line, gt_calls[ind_id])
    
    return variant_line

def get_variant_dict(chrom='2', pos='191964633', db_id='rs7574865', ref='T', 
alt='G', qual='1258.77', filt="PASS", info="", form="GT:AD:GQ:PL", 
gt_calls={'ADM1003A3': "1/1:0,35:99:1287,105,0"}):
    """
    Return a vcf formated variant line
    """
    if not info:
        info = "1000GAF=0.744609;AC=2;AF=1.00;AN=2;CLNSIG=255,10;"\
        "CSQ=G|ENSG00000138378|ENST00000358470|Transcript|intron_variant"\
        "||||||||-1|STAT4|HGNC|11365|protein_coding||CCDS2310.1|ENSP00000351255"\
        "|STAT4_HUMAN|Q53S87_HUMAN&Q53RU2_HUMAN&E9PBE2_HUMAN&C9JM11_HUMAN&C9JFG0_HUMAN"\
        "&B7SIX5_HUMAN|UPI00000015F2||||3/23||ENST00000358470.4:c.274-23582A>C|||||"\
        ",G|ENSG00000138378|ENST00000392320|Transcript|intron_variant|||||||"\
        "|-1|STAT4|HGNC|11365|protein_coding||CCDS2310.1|ENSP00000376134"\
        "|STAT4_HUMAN|Q53S87_HUMAN&Q53RU2_HUMAN&E9PBE2_HUMAN&C9JM11_HUMAN&"\
        "C9JFG0_HUMAN&B7SIX5_HUMAN|UPI00000015F2||||3/23||ENST00000392320.2:c.274-23582A>C|||||"\
        ";Clinical_db_gene_annotation=OMIM-141201,FullList;DB;DP=35;Dbsnp129LCAF=0.2493"
    
    info_dict = {}
    
    for info_entry in info.split(';'):
        splitted_info = info_entry.split('=')
        if len(splitted_info) > 1:
            info_dict[splitted_info[0]] = splitted_info[1]
        else:
            info_dict[splitted_info[0]] = []
    
    
    variant_dict = {
        'CHROM': chrom,
        'POS': pos,
        'ID': db_id,
        'REF': ref,
        'ALT': alt,
        'QUAL': qual,
        'FILTER': filt,
        'INFO': info,
        'FORMAT': form,
        'info_dict': info_dict
    }
    
    for ind_id in gt_calls:
        variant_dict[ind_id] = gt_calls[ind_id]
    
    return variant_dict


def test_chrom():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='Chrom', field='CHROM')
    test_chrom = '10'
    variant_line = get_variant_line(chrom=test_chrom)
    variant_dict = get_variant_dict(chrom=test_chrom)
    
    assert plugin.get_raw_entry(variant_line=variant_line) == test_chrom
    assert plugin.get_raw_entry(variant_dict=variant_dict) == test_chrom

def test_pos():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='Pos', field='POS')
    test_pos = '1000'
    variant_line = get_variant_line(pos=test_pos)
    variant_dict = get_variant_dict(pos=test_pos)
    
    assert plugin.get_raw_entry(variant_line=variant_line) == test_pos
    assert plugin.get_raw_entry(variant_dict=variant_dict) == test_pos

def test_id():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='ID', field='ID')
    db_id = 'rs001'
    variant_line = get_variant_line(db_id=db_id)
    variant_dict = get_variant_dict(db_id=db_id)
    
    assert plugin.get_raw_entry(variant_line=variant_line) == db_id
    assert plugin.get_raw_entry(variant_dict=variant_dict) == db_id

def test_ref():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='ref', field='REF')
    ref = 'rs001'
    variant_line = get_variant_line(ref=ref)
    variant_dict = get_variant_dict(ref=ref)
    
    assert plugin.get_raw_entry(variant_line=variant_line) == ref
    assert plugin.get_raw_entry(variant_dict=variant_dict) == ref

def test_alt():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='alt', field='ALT')
    alt = 'A'
    variant_line = get_variant_line(alt=alt)
    variant_dict = get_variant_dict(alt=alt)
    
    assert plugin.get_raw_entry(variant_line=variant_line) == alt
    assert plugin.get_raw_entry(variant_dict=variant_dict) == alt

def test_qual():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='QUAL', field='QUAL')
    qual = '1234.5'
    variant_line = get_variant_line(qual=qual)
    variant_dict = get_variant_dict(qual=qual)
    
    assert plugin.get_raw_entry(variant_line=variant_line) == qual
    assert plugin.get_raw_entry(variant_dict=variant_dict) == qual

def test_filter():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='Filter', field='FILTER')
    filt = 'PASS'
    variant_line = get_variant_line(filt=filt)
    variant_dict = get_variant_dict(filt=filt)
    
    assert plugin.get_raw_entry(variant_line=variant_line) == filt
    assert plugin.get_raw_entry(variant_dict=variant_dict) == filt

def test_1000G():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='thousand_g', field='INFO', info_key="1000GAF")
    variant_line = get_variant_line()
    variant_dict = get_variant_dict()
    
    test_value = '0.744609'
    
    assert plugin.get_raw_entry(variant_line=variant_line) == test_value
    assert plugin.get_raw_entry(variant_dict=variant_dict) == test_value
    

def test_CLNSIG():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='Clnsig', field='INFO', info_key="CLNSIG")
    test_value = '255,10'
    variant_line = get_variant_line()
    variant_dict = get_variant_dict()
    
    assert plugin.get_raw_entry(variant_line=variant_line) == test_value
    assert plugin.get_raw_entry(variant_dict=variant_dict) == test_value


def test_format():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='Format', field='FORMAT')
    form = 'GT:AD:GQ:PL'
    variant_line = get_variant_line(form=form)
    variant_dict = get_variant_dict(form=form)

    assert plugin.get_raw_entry(variant_line=variant_line) == form
    assert plugin.get_raw_entry(variant_dict=variant_dict) == form


def test_dict_two_separators_no_key():
    """Test to get the raw chromosome"""
    variant_dict = get_variant_dict(info="RankScore=1:12,2:11")
    variant_line = get_variant_line(info="RankScore=1:12,2:11")
    
    plugin = Plugin(
        name='rank_score', 
        field='INFO', 
        info_key="RankScore", 
        separators=[',', ':'],
        dict_entry=True
    )
    dict_entry = plugin.get_raw_entry(variant_dict=variant_dict)
    line_entry = plugin.get_raw_entry(variant_line=variant_line)
    assert dict_entry == '11'
    assert line_entry == '11'

def test_dict_two_separators_with_key():
    """Test to get the raw chromosome"""
    variant_dict = get_variant_dict(info="RankScore=1:12,2:11")
    variant_line = get_variant_line(info="RankScore=1:12,2:11")
    
    plugin = Plugin(
        name='rank_score', 
        field='INFO', 
        info_key="RankScore", 
        separators=[',', ':'],
        dict_entry=True
    )
    dict_entry = plugin.get_raw_entry(
        variant_dict=variant_dict,
        dict_key='1'
    )
    line_entry = plugin.get_raw_entry(
        variant_line=variant_line,
        dict_key='1'
    )
    
    assert dict_entry == '12'
    assert line_entry == '12'

def test_dict_three_separators_no_key():
    """Test to get the raw chromosome"""
    variant_dict = get_variant_dict(info="GeneticModels=1:AD|AD_dn,2:AR_hom")
    variant_line = get_variant_line(info="GeneticModels=1:AD|AD_dn,2:AR_hom")
    
    plugin = Plugin(
        name='genetic_models',
        field='INFO',
        info_key="GeneticModels",
        separators=[',', ':', '|'],
        dict_entry=True
    )

    dict_entry = plugin.get_raw_entry(variant_dict=variant_dict)
    line_entry = plugin.get_raw_entry(variant_line=variant_line)

    assert dict_entry == 'AR_hom'
    assert line_entry == 'AR_hom'

def test_dict_three_separators_with_key():
    """Test to get the raw chromosome"""
    variant_dict = get_variant_dict(info="GeneticModels=1:AD|AD_dn,2:AR_hom")
    variant_line = get_variant_line(info="GeneticModels=1:AD|AD_dn,2:AR_hom")
    
    plugin = Plugin(
        name='genetic_models',
        field='INFO',
        info_key="GeneticModels",
        separators=[',', ':', '|'],
        dict_entry=True
    )

    dict_entry = plugin.get_raw_entry(
        variant_dict=variant_dict,
        dict_key='1'
    )
    line_entry = plugin.get_raw_entry(
        variant_line=variant_line,
        dict_key='1'
    )

    assert dict_entry == 'AD|AD_dn'
    assert line_entry == 'AD|AD_dn'


def test_gt_call():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='gt_call', field='sample_id', gt_key='GT')
    
    vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'ADM1003A3']
    
    variant_line = get_variant_line()
    variant_dict = get_variant_dict()
    
    assert plugin.get_raw_entry(
        variant_line=variant_line, 
        vcf_header=vcf_header, 
        individual_id='ADM1003A3') == '1/1'

    assert plugin.get_raw_entry(
        variant_dict=variant_dict, 
        vcf_header=vcf_header, 
        individual_id='ADM1003A3') == '1/1'

def test_gt_call_no_header():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='gt_call', field='sample_id', gt_key='GT')
    variant_line = get_variant_line()
    
    with pytest.raises(IOError):
        plugin.get_raw_entry(variant_line=variant_line, individual_id='ADM1003A3')

def test_gt_call_no_gt_key():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='gt_call', field='sample_id')
    variant_line = get_variant_line()
    variant_dict = get_variant_dict()
    vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'ADM1003A3']

    with pytest.raises(IOError):
        plugin.get_raw_entry(variant_line=variant_line, vcf_header=vcf_header, individual_id='ADM1003A3')

    with pytest.raises(IOError):
        plugin.get_raw_entry(variant_dict=variant_dict, vcf_header=vcf_header, individual_id='ADM1003A3')
    

def test_gt_call_no_header_dict():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='gt_call', field='sample_id', gt_key='GT')
    variant_dict = get_variant_dict()
    
    assert plugin.get_raw_entry(variant_dict=variant_dict, individual_id='ADM1003A3') == '1/1'

def test_gt_call_no_individual():
    """Test to get the raw chromosome"""
    plugin = Plugin(name='gt_call', field='sample_id', gt_key='GT')
    vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'ADM1003A3']
    variant_line = get_variant_line()
    variant_dict = get_variant_dict()

    with pytest.raises(IOError):
        plugin.get_raw_entry(variant_line=variant_line, vcf_header=vcf_header)

    with pytest.raises(IOError):
        plugin.get_raw_entry(variant_dict=variant_dict, vcf_header=vcf_header)
