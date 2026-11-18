from genmod.annotate_regions.parse_annotations import BedRegion, bed_parser


def test_bed_parser(bed_lines: list[str]) -> None:
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
        assert region.symbol in symbols

    assert index + 1 == nr_regions


def test_bed_parser_with_mt_chromosome() -> None:
    # Given a bed line with the mitochondrial chromosome encoded as MT
    bed_lines = ["MT\t1000\t5000\tgene_MT"]
    # Then the output region should not be padded
    expected_output = [BedRegion("M", 1000, 5000, "gene_MT")]
    result = list(bed_parser(bed_lines, padding=4000))
    assert result == expected_output


def test_bed_parser_with_chrm_chromosome() -> None:
    # Given a bed line with the mitochondrial chromosome encoded as chrM instead of MT
    bed_lines = ["chrM\t1000\t5000\tgene_chrM"]
    # Then the output region should not be padded and the chrom encoded as M
    expected_output = [BedRegion("M", 1000, 5000, "gene_chrM")]
    result = list(bed_parser(bed_lines, padding=4000))
    assert result == expected_output


def test_bed_parser_skips_comments_and_empty_lines() -> None:
    # GIVEN bed lines with a comment, an empty line and one region
    bed_lines = ["#header\n", "\n", "1\t100\t200\tGENE\n"]

    # WHEN parsing
    result = list(bed_parser(bed_lines, padding=0))

    # THEN only the region is returned
    assert result == [BedRegion("1", 100, 200, "GENE")]


def test_bed_parser_pads_start_but_not_below_zero() -> None:
    # GIVEN a gene close to the chromosome start and one further in
    bed_lines = ["1\t1000\t2000\tNEAR_START", "1\t10000\t11000\tINSIDE"]

    # WHEN parsing with padding 4000
    result = list(bed_parser(bed_lines, padding=4000))

    # THEN the start is padded but never negative, and the end is not padded
    assert result == [
        BedRegion("1", 0, 2000, "NEAR_START"),
        BedRegion("1", 6000, 11000, "INSIDE"),
    ]


def test_bed_parser_uses_line_index_as_symbol_without_name_column() -> None:
    # GIVEN a bed line with only chrom, start and end (it is the second line)
    bed_lines = ["#header", "2\t10\t20"]

    # WHEN parsing
    result = list(bed_parser(bed_lines, padding=0))

    # THEN the index of the line is used as the symbol
    assert result == [BedRegion("2", 10, 20, "1")]
