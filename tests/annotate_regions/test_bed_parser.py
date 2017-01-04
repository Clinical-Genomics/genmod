from genmod.annotate_regions.parse_annotations import bed_parser

def test_bed_parser(bed_lines):
    # GIVEN an iterable with bed lines
    nr_regions = 0
    symbols = []
    for line in bed_lines:
        if not line.startswith('#'):
            nr_regions += 1
            symbols.append(line.split()[3])
    # WHEN parsing the bedlines
    for index, region in enumerate(bed_parser(bed_lines)):
    # THEN assert the symbols are found in the bed lines
        assert region['symbol'] in symbols
    
    assert index+1 == nr_regions
        