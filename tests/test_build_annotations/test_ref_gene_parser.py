from genmod.annotate_regions.parse_annotations import gff_parser

def get_gff_lines():
    """
    Create a small ccds file
    """
    gff_lines = [
        "##gff-version 3",
        "##sequence-region 121CI_contig1 1 8876",
        "121CI_contig1\tProdigal:2.6\tCDS\t608\t748\t.\t-\t0\tID=121CI_00001;"\
        "inference=ab initio prediction:Prodigal:2.6;locus_tag=121CI_00001;"\
        "product=hypothetical protein",
        "121CI_contig1\tProdigal:2.6\tCDS\t832\t1593\t.\t-\t0\tID=121CI_00002;"\
        "inference=ab initio prediction:Prodigal:2.6,similar to AA sequence"\
        ":Protein_fasta_for_annot.faa:gi|15646195;locus_tag=121CI_00002"\
        ";note=protein of unknown function UPF0174 ;product=hypothetical protein"
    ]
    
    return gff_lines

def test_build_gff_genes():
    """
    Test if the parser builds the gff gene dict correct
    """
    gff_lines = get_gff_lines()
    gene_trees, exon_trees = gff_parser(gff_lines)
    
    assert set(gene_trees['121CI_contig1'].keys()) == set(["121CI_00001", "121CI_00002"])
    assert gene_trees['121CI_contig1']["121CI_00001"] == [608, 748, "121CI_00001"]
    
    # assert set(exon_trees['13'].keys()) == set(["5057114250571899", "5058607050592603"])
    #
    # assert exon_trees['13']["5057114250571899"] == [50571140, 50571901, "5057114250571899"]
