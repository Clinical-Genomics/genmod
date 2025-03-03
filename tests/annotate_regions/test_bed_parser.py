from genmod.annotate_regions.parse_annotations import bed_parser


def test_bed_parser(bed_lines):
    # GIVEN an iterable with bed lines
    nr_regions = 0
    symbols = []
    for line in bed_lines:
        if not line.startswith("#"):
            nr_regions += 1
            symbols.append(line.split()[3])
    # WHEN parsing the bedlines
    for index, region in enumerate(bed_parser(bed_lines)):
        # THEN assert the symbols are found in the bed lines
        assert region["symbol"] in symbols

    assert index + 1 == nr_regions

def test_bed_parser_with_mt_chromosome():
    # Given a bed line with the mitochondrial chromosome encoded as MT
    bed_lines = [
        "MT\t1000\t5000\tgene_MT"
    ]
    # Then the output region should not be padded
    expected_output = [
        {'chrom': 'M', 'start': 1000, 'stop': 5000, 'symbol': 'gene_MT'}
    ]
    result = list(bed_parser(bed_lines, padding=4000))
    assert result == expected_output

def test_bed_parser_with_chrm_chromosome():
    # Given a bed line with the mitochondrial chromosome encoded as chrM instead of MT
    bed_lines = [
        "chrM\t1000\t5000\tgene_chrM"
    ]
    # Then the output region should not be padded and the chrom encoded as M
    expected_output = [
        {'chrom': 'M', 'start': 1000, 'stop': 5000, 'symbol': 'gene_chrM'}
    ]
    result = list(bed_parser(bed_lines, padding=4000))
    assert result == expected_output

