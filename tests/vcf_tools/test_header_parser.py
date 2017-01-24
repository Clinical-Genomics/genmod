from genmod.vcf_tools.header_parser import HeaderParser

def test_parse_info():
    ## GIVEN a header object
    head = HeaderParser()
    assert 'MQ' not in head.info_dict
    info_line = '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">'
    
    ## WHEN parsing a correct info line
    head.parse_meta_data(info_line)
    
    ## THEN assert it is added to the parser
    assert 'MQ' in head.info_dict


def test_parse_contig():
    ## GIVEN a header object
    head = HeaderParser()
    assert '1' not in head.contig_dict
    contig_line = '##contig=<ID=1,length=249250621,assembly=b37>'
    
    ## WHEN parsing a correct info line
    head.parse_meta_data(contig_line)
    
    ## THEN assert it is added to the parser
    assert '1' in head.contig_dict

def test_parse_contig_no_length():
    ## GIVEN a header object
    head = HeaderParser()
    assert '1' not in head.contig_dict
    contig_line = '##contig=<ID=1,assembly=b37>'
    
    ## WHEN parsing a correct info line
    head.parse_meta_data(contig_line)
    
    ## THEN assert it is added to the parser
    assert '1' in head.contig_dict

def test_parse_minimal_contig():
    ## GIVEN a header object
    head = HeaderParser()
    assert '1' not in head.contig_dict
    contig_line = '##contig=<ID=1>'
    
    ## WHEN parsing a correct info line
    head.parse_meta_data(contig_line)
    
    ## THEN assert it is added to the parser
    assert '1' in head.contig_dict