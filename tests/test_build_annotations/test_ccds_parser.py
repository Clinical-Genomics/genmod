from genmod.annotate_regions.parse_annotations import ccds_parser

def test_build_ccds_genes(ccds_lines):
    """
    Test if the parser builds the ccds gene dict correct
    """
    gene_trees, exon_trees = ccds_parser(ccds_lines)
    assert set(gene_trees['1'].keys()) == set(["LINC00115"])
    assert gene_trees['1']["LINC00115"] == [801942, 802433, "LINC00115"]
    assert exon_trees['1'] == {'801942802433':[801940, 802435, '801942802433']}
