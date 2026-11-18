from genmod.vcf_tools.header_parser import ANN_SPEC_COLUMNS, HeaderParser

# The ANN header line as written by SnpEff
SNPEFF_ANN_LINE = (
    '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
    "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | "
    "Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | "
    "CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' \">"
)


def test_parse_snpeff_ann_columns():
    ## GIVEN a header object
    head = HeaderParser()

    ## WHEN parsing the SnpEff ANN info line
    head.parse_meta_data(SNPEFF_ANN_LINE)

    ## THEN the ANN fields are stored without quotes or spaces, in the spec order
    assert head.ann_columns == ANN_SPEC_COLUMNS
    assert head.annotation_columns("ANN") == ANN_SPEC_COLUMNS
    assert head.vep_columns == []


def test_parse_ann_without_field_list_uses_spec_columns():
    ## GIVEN an ANN info line that does not list the fields
    head = HeaderParser()
    info_line = '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">'

    ## WHEN parsing it
    head.parse_meta_data(info_line)

    ## THEN the fields of the SnpEff spec are used
    assert head.ann_columns == ANN_SPEC_COLUMNS


def test_parse_vep_csq_columns():
    ## GIVEN a VEP CSQ info line
    head = HeaderParser()
    info_line = (
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from '
        'Ensembl VEP. Format: Allele|Consequence|Gene">'
    )

    ## WHEN parsing it
    head.parse_meta_data(info_line)

    ## THEN the CSQ columns are stored and there are no ANN columns
    assert head.annotation_columns("CSQ") == ["Allele", "Consequence", "Gene"]
    assert head.ann_columns == []


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
