from genmod.annotate_regions.parse_annotations import gff_parser


def test_build_gff_genes(gff_lines):
    """
    Test if the parser builds the gff gene dict correct
    """
    gene_trees, exon_trees = gff_parser(gff_lines)
    
    assert set(gene_trees['121CI_contig1'].keys()) == set(["121CI_00001", "121CI_00002"])
    assert gene_trees['121CI_contig1']["121CI_00001"] == [608, 748, "121CI_00001"]
    
    # assert set(exon_trees['13'].keys()) == set(["5057114250571899", "5058607050592603"])
    #
    # assert exon_trees['13']["5057114250571899"] == [50571140, 50571901, "5057114250571899"]
