from genmod.vcf_tools.header_parser import HeaderParser


def test_parse_info():
    ## GIVEN a header object
    head = HeaderParser()
    assert "MQ" not in head.info_dict
    info_line = '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">'

    ## WHEN parsing a correct info line
    head.parse_meta_data(info_line)

    ## THEN assert it is added to the parser
    assert "MQ" in head.info_dict


def test_parse_info_with_source():
    ## GIVEN a header object
    head = HeaderParser()
    assert "MQ" not in head.info_dict
    info_line = '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality",Source="bwa">'

    ## WHEN parsing a correct info line
    head.parse_meta_data(info_line)

    ## THEN assert it is added to the parser
    assert "MQ" in head.info_dict


def test_parse_info_with_version():
    ## GIVEN a header object
    head = HeaderParser()
    assert "MQ" not in head.info_dict
    info_line = (
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality",Version="1.0.0">'
    )

    ## WHEN parsing a correct info line
    head.parse_meta_data(info_line)

    ## THEN assert it is added to the parser
    assert "MQ" in head.info_dict


def test_parse_info_with_source_and_version():
    ## GIVEN a header object
    head = HeaderParser()
    assert "MQ" not in head.info_dict
    info_line = '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality",Source="bwa",Version="1.0.0">'

    ## WHEN parsing a correct info line
    head.parse_meta_data(info_line)

    ## THEN assert it is added to the parser
    assert "MQ" in head.info_dict


def test_parse_contig():
    ## GIVEN a header object
    head = HeaderParser()
    assert "1" not in head.contig_dict
    contig_line = "##contig=<ID=1,length=249250621,assembly=b37>"

    ## WHEN parsing a correct info line
    head.parse_meta_data(contig_line)

    ## THEN assert it is added to the parser
    assert "1" in head.contig_dict


def test_parse_contig_no_length():
    ## GIVEN a header object
    head = HeaderParser()
    assert "1" not in head.contig_dict
    contig_line = "##contig=<ID=1,assembly=b37>"

    ## WHEN parsing a correct info line
    head.parse_meta_data(contig_line)

    ## THEN assert it is added to the parser
    assert "1" in head.contig_dict


def test_parse_minimal_contig():
    ## GIVEN a header object
    head = HeaderParser()
    assert "1" not in head.contig_dict
    contig_line = "##contig=<ID=1>"

    ## WHEN parsing a correct info line
    head.parse_meta_data(contig_line)

    ## THEN assert it is added to the parser
    assert "1" in head.contig_dict


def test_print_header_prints_fileformat_once_and_keeps_filedate() -> None:
    # GIVEN a header with a fileformat and a fileDate line
    head = HeaderParser()
    head.parse_meta_data("##fileformat=VCFv4.2")
    head.parse_meta_data("##fileDate=20260101")
    head.parse_header_line("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    # WHEN printing the header
    lines = head.print_header()

    # THEN the fileformat is printed once, first, and the fileDate line is kept as it was
    assert lines.count("##fileformat=VCFv4.2") == 1
    assert lines[0] == "##fileformat=VCFv4.2"
    assert "##fileDate=20260101" in lines
    assert lines[-1].startswith("#CHROM")
